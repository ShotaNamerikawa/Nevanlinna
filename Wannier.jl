module Hamiltonian


module HR
using LinearAlgebra
using Distributed
using TensorOperations

mutable struct Hr{T}
    rpoints::Array{Int64,2}
    ir0::Int64
    hopping::Array{Complex{T},3}
    weight::Vector{Int64}
    num_wann::Int64
end

function make_hr_from_file(winfile::String;T=Float64)
    fp = open(winfile)
    weight = Int[]
    readline(fp)
    num_wann = parse(Int64,readline(fp))
    rpointsnum = parse(Int64,readline(fp))
    for i in 1:Int64(ceil(rpointsnum/15))
        line = readline(fp)
        append!(weight,parse.(Int64,split(line)))
    end
    lines = readlines(fp)
    close(fp)
    rpoints = zeros(Int64,(rpointsnum,3))
    hopping = zeros(Complex{T},(num_wann,num_wann,rpointsnum))
    ir0 = 0
    for x in 1:rpointsnum
        line = split.(lines[(x-1)*num_wann^2+1:x*num_wann^2])
        rcoord = parse.(Int64,line[1][1:3])
        #println(typeof(rcoord))
        if rcoord == [0, 0, 0]
            ir0 = x
        end
        rpoints[x,:] = rcoord 
        for i in 1:num_wann
            for j in 1:num_wann
                hopsite = parse.(Float64,popfirst!(line)[6:7])
                hopping[j,i,x] = hopsite[1] + 1.0im*hopsite[2]
            end
        end
    end
    Hr(rpoints,ir0,hopping,weight,num_wann)
end

function ftransform(hr::Hr{T},kpoints::Array) where{T} 
    phase = hr.rpoints * kpoints
    phasefactor=exp.((2.0im*pi).*phase)
    hkw = zeros(Complex{T},(size(hr.hopping,1),size(hr.hopping,2)))
    for j in 1:size(hr.hopping,2)
        for i in 1:size(hr.hopping,1)
            hkw[i,j] = dot(phasefactor,hr.hopping[i,j,:] ./ hr.weight)
        end
    end
    hkw
end

function generate_hkwa(hr::Hr{T},kpoints::Array) where {T}
    phase = hr.rpoints * kpoints
    phasefactor=exp.((2.0im*pi).*phase)
    hkwa = zeros(Complex{T},(size(hr.hopping,1),size(hr.hopping,2),3))
    for l in 1:3
        weight = 2.0im*pi*hr.rpoints[:,l].*phasefactor
        for j in 1:size(hr.hopping,2)
            for i in 1:size(hr.hopping,1)
                hkwa[i,j,l] = dot(weight,hr.hopping[i,j,:] ./ hr.weight)
            end
        end
    end
    hkwa
end

function generate_hkha(hkwa::Array{Complex{T},2},U::Array{Complex{T},2}) where {T}
    @tensoropt hkha[i,j,k] = conj(U)[l,i]*hkwa[l,m,k]*U[m,j] 
end

function generate_va(hkwa::Array{Complex{T},2},U::Array{Complex{T},2}) where {T}
    @tensoropt hkha[i,k] = conj(U)[l,i]*hkwa[l,m,k]*U[m,i] 
end

end

module HamKmesh
using ..HR
using Distributed
using LinearAlgebra

mutable struct Hamkmesh{T}
    kmesh::Array
    kp_list::Array
    kp_shift::Array
    nk::Int64
    hr::HR.Hr
    #hk_list::Array
    #hkwa_list::Array
    function Hamkmesh(kmesh,hr::HR.Hr{T}) where T
        hkm = new{T}()
        if size(kmesh)[1] == 6
            hkm.kmesh = kmesh[1:3]
            hkm.kp_shift = kmesh[4:6]
        elseif size(kmesh)[1] == 3
            hkm.kmesh = kmesh[1:3]
            hkm.kp_shift = zeros(Float64,3)
        else
            error("the length of kmesh has to be 3 or 6")
        end
        hkm.nk = prod(hkm.kmesh)
        hkm.kp_list = [zeros(Float64,3) for i in 1:hkm.nk]
        generate_kp_list!(hkm,hkm.kmesh)
        hkm.hr = hr
        hkm
    end
end

function generate_kp_list!(hkm::Hamkmesh{T},kmesh) where{T}
    for i in 0:hkm.kmesh[1]-1
    for j in 0:hkm.kmesh[2]-1
    for k in 0:hkm.kmesh[3]-1
        (hkm.kp_list[i*hkm.kmesh[2]*hkm.kmesh[3]+j*hkm.kmesh[3]+k+1][:] 
        = [(i+0.5*hkm.kp_shift[1]),(j+0.5*hkm.kp_shift[2]),
            (k+0.5*hkm.kp_shift[3])]./hkm.kmesh)
    end
    end
    end
end

function generate_hk(hkm::Hamkmesh{T}) where{T}
    hk_list = [zeros(Complex,hkm.hr.num_wann,hkm.hr.num_wann) for x in 1:hkm.nk]
    for x in 1:hkm.nk
        hk_list[x][:,:] = HR.ftransform(hkm.hr,hkm.kp_list[x][:])
    end
    hk_list
end

function diagonalize_hk(hk) 
    eigvals(hk)
end

#=
function apply_func_ksum(func,hkm::Hamkmesh)
    for i in 1:size(hkm.kp_list)[1]
        func(hkm.kp_list[i])
    end
end

function apply_func_ksum!(func,result,hkm::Hamkmesh)
    for i in 1:size(hkm.kp_list)[1]
        result[:,:] += func(hkm.kp_list[i])
    end
end
=#
function apply_func_ksum!(func,result::Array,hkm::Hamkmesh)
    #println("apply func ksum")
    for x in 1:nworkers()
        @spawnat x+1 func = func
    end
    pool = CachingPool(workers())
    value = sum(pmap(func,pool,hkm.kp_list))
    result .= value
end

end

end

#using DoubleFloats

#wan = Hr.make_hr_from_file("/home/shota/Nevanlinna/nev_test/Fe_hr.dat")
#d64 = Hamiltonian.ftrans(wan,convert.(Double64,[0.0,0.0,0.0]))
#f64 = ftransform(wan,convert.(Float64,[0.1,0.1,0.1]))
#=
for i in 1:wan.num_wann
    println(f64[i,i])
end
=#
#println("d64")
#println(inv(d64)[2,2])
#println("f64")
#println(inv(f64)[2,2])
