include("./Wannier.jl")

module CPA
using LinearAlgebra
using ..Hamiltonian
using Roots
using Distributed
using PyCall
sys = pyimport("sys")
sys.path = push!(sys.path,"/home/shota/wannier_utils/src/wannier_utils")
temperature_spir = pyimport("temperature_spir")

export loop

na = [CartesianIndex()]

mutable struct Cpa{T}
    temp::PyObject
    hr_list::Array
    ir0::Int64
    num_wann::Int64
    cpa_hr::Hamiltonian.HR.Hr{T}
    weight_list::Array{Number}
    cpa_hkm::Hamiltonian.HamKmesh.Hamkmesh
    pot_list::Vector{Vector}
    onsite_pot_list::Array{T}
    ne::Real
    mu::Real
    green::Vector{Array{Complex{T},2}}
    pot_iwn::Vector{Array{Complex{T},2}}
    tmat_iwn::Vector{Array{Complex{T},2}}

    function Cpa(temperature,kmesh,ne::Real,hr_list,weight_list::Array
        ,onsite_pot_list::Array{T};wmax=100) where{T}
        cpa = new{T}()
        cpa.temp = spir_instance_make(temperature,wmax)
        cpa.ne = ne
        cpa.hr_list = hr_list
        rpoints = cpa.hr_list[1].rpoints
        cpa.ir0 = cpa.hr_list[1].ir0
        rweight = cpa.hr_list[1].weight
        cpa.num_wann = cpa.hr_list[1].num_wann
        cpa.cpa_hr = Hamiltonian.HR.Hr(rpoints,cpa.ir0,zeros(Complex{T},size(cpa.hr_list[1].hopping))
                    ,rweight,cpa.num_wann)
        cpa.weight_list = weight_list
        cpa.onsite_pot_list = onsite_pot_list
        cpa.mu = 0.0
        cpa.green = [zeros(Complex{T},(cpa.num_wann,cpa.num_wann)) for
                    i in 1:size(cpa.temp.omegaF_sp)[1]]
        cpa.pot_iwn = [zeros(Complex{T},(cpa.num_wann,cpa.num_wann)) for
                    i in 1:size(cpa.temp.omegaF_sp)[1]]
        cpa.tmat_iwn = [zeros(Complex{T},(cpa.num_wann,cpa.num_wann)) for
                    i in 1:size(cpa.temp.omegaF_sp)[1]]
        println(cpa.hr_list[1].ir0)
        println(cpa.hr_list[2].ir0)
        init_cpa_hamiltonian!(cpa)
        cpa.cpa_hkm = Hamiltonian.HamKmesh.Hamkmesh(kmesh,cpa.cpa_hr)
        cpa
    end
end

function spir_instance_make(temperature,wmax)
    temperature_spir.Temperature(temperature,wmax)
end


function loop(cpa::Cpa;max_niter=100,cpa_thr=10^(-8))
    for i in 1:max_niter
        @time calc_green_function!(cpa)
        tmat_list = [calc_tmat(cpa,cpa.pot_list[i]) for i in 1:size(cpa.pot_list)[1]]
        cpa.tmat_iwn[:][:,:] = sum(cpa.weight_list.*tmat_list)
        sum_abs = sum([sum(abs.(cpa.tmat_iwn[i])) for i in 1: size(cpa.temp.omegaF_sp)[1]])
        if sum_abs < cpa_thr
            return i
        end
        print("iter = ")
        print(i)
        print(", delta = ")
        println(sum_abs)
        calc_new_pot(cpa)
    end
end


function init_cpa_hamiltonian!(cpa::Cpa)
    for (hr, pot) in zip(cpa.hr_list, cpa.onsite_pot_list)
        #@assert hr.ir0 == cpa.ir0
        #@assert hr.num_wann == cpa.num_wann
        for nw in 1:cpa.num_wann
            cpa.cpa_hr.hopping[nw,nw,cpa.ir0] -= pot
        end
    end
    cpa.pot_list = [zeros(Complex,(cpa.num_wann)) 
                    for i in 1:size(cpa.weight_list)[1]]
    for i in 1:size(cpa.weight_list)[1]
        cpa.pot_list[i][:] = diag(cpa.hr_list[i].hopping[:,:,cpa.ir0],0)
    end
    pot = sum(cpa.weight_list.*cpa.pot_list)
    pot_d = Diagonal(pot)
    for i in 1:size(cpa.temp.omegaF_sp)[1]
        cpa.pot_iwn[i][:,:] = pot_d
    end
    for i in 1:size(cpa.weight_list)[1]
        cpa.cpa_hr.hopping[:,:,:] .+= cpa.weight_list[i].*cpa.hr_list[i].hopping
    end

    for nw in 1:cpa.num_wann
        cpa.cpa_hr.hopping[nw,nw,cpa.ir0] = 0
    end

end

function calc_green_func_k(cpa,k::Array)
    hk = Hamiltonian.HR.ftransform(cpa.cpa_hr,k)
    E = Matrix{Complex}(I,cpa.num_wann,cpa.num_wann)
    hk_mu = hk - cpa.mu.*E
    @. inv(1.0im*cpa.temp.omegaF_sp*[E] - [hk_mu] - cpa.pot_iwn)
end

function calc_tr_gkiwn(cpa,k::Array)
    tr.(calc_green_func_k(cpa,k))
end

function calc_ne(cpa,mu)
    cpa.mu = mu
    trg = zeros(Complex,size(cpa.temp.omegaF_sp)[1])
    applyfunc(k) = calc_tr_gkiwn(cpa,k)/(cpa.cpa_hkm.nk)
    println("trg")
    @time Hamiltonian.HamKmesh.apply_func_ksum!(applyfunc,trg,cpa.cpa_hkm)
    println("gl")
    @time gl = cpa.temp.smpl_F.fit(trg)
    println("Ultau")
    @time Ultau = -cpa.temp.basis_F.u(cpa.temp.beta)
    ne = real(dot(gl,Ultau))
end

function calc_green_function!(cpa::Cpa{T}) where T
    dmu = 30
    f = (x -> (value = calc_ne(cpa,x) - cpa.ne;println(value);value))
    g = (x -> calc_green_func_k(cpa,x))
    cpa.mu = find_zero(f,(-dmu+cpa.mu,cpa.mu+dmu),Bisection())
    result = [zeros(Complex{T},cpa.num_wann,cpa.num_wann) for i in 1:size(cpa.temp.omegaF_sp)[1]]
    @time cpa.green[:] = Hamiltonian.HamKmesh.apply_func_ksum!(g,result,cpa.cpa_hkm)/cpa.cpa_hkm.nk
end

function calc_tmat(cpa::Cpa,pot)
    E = Matrix{Complex}(I,cpa.num_wann,cpa.num_wann)
    potd = [Diagonal(pot)]
    v = potd .- cpa.pot_iwn
    invtmp = inv.([E] .- cpa.green.*v)
    value = v.*invtmp
    value
end

function calc_new_pot(cpa::Cpa)
    tmp = cpa.green .* cpa.tmat_iwn
    tmp = inv.([E] .+ tmp)
    cpa.pot_iwn[:] += cpa.tmat_iwn.*tmp 
end


end



#=
using .hr
using  .hamkmesh
using .cpa
using LinearAlgebra
using DoubleFloats
using Printf

function write_trG(cpa,upfile::String,downfile::String)
trGup = zeros(Complex{Float64},size(cpa.temp.omegaF_sp)[1])
trGdn = zeros(Complex{Float64},size(cpa.temp.omegaF_sp)[1])

for i in 1:cpa.num_wann
    if i%2 != 0
        trGup += getindex.(cpa.green,LinearIndices(cpa.green[1])[i,i])
    elseif i%2 == 0
        trGdn += getindex.(cpa.green,LinearIndices(cpa.green[1])[i,i])
    end
end

fpup = open(upfile,"w")
fpdn = open(downfile,"w")

for i in 1:size(cpa.temp.omegaF_sp)[1]
    write(fpup,@sprintf("%20.16e %20.16e %20.16e\n",cpa.temp.omegaF_sp[i]
                        ,real(trGup[i]),imag(trGup[i])))
    write(fpdn,@sprintf("%20.16e %20.16e %20.16e\n",cpa.temp.omegaF_sp[i]
                        ,real(trGdn[i]),imag(trGdn[i])))
end

close(fpup)
close(fpdn)

end
=#


#hr_Fe = hr.readfile("/home/shota/Nevanlinna/nev_test/Fe_hr.dat")
#hr_Co = hr.readfile("/home/shota/Nevanlinna/nev_test/Co_hr.dat")
#println(hr_Fe.rpoints)
#=
for x in 1:597
    print("weight_inv")
    @printf("%d ",hr_Fe.weight[x])
    print(" real: ")
    @printf("%.30e ",real(hr_Fe.hopping[1,2,x]))
    print("imag: ")
    @printf("%.30e\n",imag(hr_Fe.hopping[1,2,x]))
end
#weight_inv = 1 ./ hr_Fe.weight[1:20]
#println(hr_Fe.hopping[1,2,1:20])
#println(sum(hr_Fe.hopping[1,2,1:20]./weight_inv))
for x in 1:597
    @printf("%.30e ",real(sum(hr_Fe.hopping[1,2,x] / hr_Fe.weight[x])))
    @printf("%.30e\n",imag(sum(hr_Fe.hopping[1,2,x] / hr_Fe.weight[x])))
end
println("sum")
for x in 1:597
    @printf("%.30e ",real(sum(hr_Fe.hopping[1,2,1:x] ./ hr_Fe.weight[1:x])))
    @printf("%.30e\n",imag(sum(hr_Fe.hopping[1,2,1:x] ./ hr_Fe.weight[1:x])))
end
=#

#=
hr_list = [hr_Fe hr_Co]
weightlist = [1.0 ,0.0]
onsite_pot_list = [29.84,24.95]
T=300
kmesh = [8,8,8]
test = cpa.Cpa(T,kmesh,8.0,hr_list,weightlist,onsite_pot_list)
cpa.loop(test)

write_trG(test,"jl_trG_up_Fe_10_Co_0.dat","jl_trG_dn_Fe_10_Co_0.dat")
=#
