using PyCall
sys = pyimport("sys")
sys.path = push!(sys.path,"/home/shota/wannier_utils/src/wannier_utils")
push!(LOAD_PATH,"/home/shota/Nevanlinna")
push!(LOAD_PATH,"/home/shota/Nevanlinna/wannier.jl")

module hr
    using LinearAlgebra
    
    export HR,ftransform,generate_hkwa

    mutable struct HR{T}
        rpoints::Array{Int64,2}
        ir0::Int64
        hopping::Array{Complex{T},3}
        weight::Vector{Int64}
        num_wann::Int64
    end

    function readfile(winfile::String;T=Float64)
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
        HR(rpoints,ir0,hopping,weight,num_wann)
    end

    function ftransform(hr::HR,kpoints::Array{T}) where{T}
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

    function generate_hkwa(hr::HR,kpoints::Array{T}) where {T}
        phase = hr.rpoints * kpoints
        phasefactor=exp.((2.0im*pi).*phase)
        hkwa = zeros(Complex{T},(size(hr.hopping,1),size(hr.hopping,2),3))
        for l in 1:3
            weight = 2.0im*pi*hk.hr.rpoints[:,l]*phasefactor
            for j in 1:size(hr.hopping,2)
                for i in 1:size(hr.hopping,1)
                    hkwa[i,j,l] = dot(weight,hr.hopping[i,j,:] ./ hr.weight)
                end
            end
        end
        hkwa
    end

end

module hamkmesh
    using ..hr
    export Hamkmesh,apply_func_ksum

    mutable struct Hamkmesh{T}
        kmesh::Array
        kp_list::Array
        kp_shift::Array
        nk::Int64
        hr::HR
        hk_list::Array
        #hkwa_list::Array
        function Hamkmesh(kmesh,hr::HR{T}) where T
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
            hkm.kp_list = zeros(Float64,hkm.nk,3) 
            generate_kp_list!(hkm,hkm.kmesh)
            hkm.hr = hr
            hkm.hk_list = zeros(Complex{T},(hkm.nk,hkm.hr.num_wann,hkm.hr.num_wann))
            generate_hk!(hkm)
            hkm
        end
    end

    function generate_kp_list!(hkm::Hamkmesh{T},kmesh) where{T}
        for i in 0:hkm.kmesh[1]-1
        for j in 0:hkm.kmesh[2]-1
        for k in 0:hkm.kmesh[3]-1
            (hkm.kp_list[i*hkm.kmesh[2]*hkm.kmesh[3]+j*hkm.kmesh[3]+k+1,:] 
            = [(i+0.5*hkm.kp_shift[1])/kmesh[1],(j+0.5*hkm.kp_shift[2])/kmesh[2],
              (k+0.5*hkm.kp_shift[3])/kmesh[3]])
        end
        end
        end
    end

    function generate_hk!(hkm::Hamkmesh{T}) where{T}
        for x in 1:hkm.nk
            hkm.hk_list[x,:,:] = ftransform(hkm.hr,hkm.kp_list[x,:])
        end
    end

    function apply_func_ksum(hkm::Hamkmesh,func)
        result = zeros(Complex,(hkm.hr.num_wann,hkm.hr.num_wann))
        for i in 1:hkm.nk
            result += func(hkm.kp_list[i,:])
        end
        result
    end
end

module cpa
    using LinearAlgebra
    using ..hr
    using ..hamkmesh
    using TensorOperations
    using Roots
    using PyCall
    temperature_spir = pyimport("temperature_spir")

    export loop

    na = [CartesianIndex()]

    mutable struct CPA{T}
        temp::PyObject
        hr_list::Array{HR}
        ir0::Int64
        num_wann::Int64
        cpa_hr::HR
        weight_list::Array{Number}
        cpa_hkm::Hamkmesh
        pot_list::Vector{Vector}
        onsite_pot_list::Array{T}
        ne::Real
        mu::Real
        green::Vector{Array{Complex{T},2}}
        pot_iwn::Vector{Array{Complex{T},2}}
        tmat_iwn::Vector{Array{Complex{T},2}}
        function CPA(temperature,kmesh,ne::Real,hr_list,weight_list::Array
            ,onsite_pot_list::Array{T};wmax=100) where{T}
            cpa = new{T}()
            cpa.temp = spir_instance_make(temperature,wmax)
            cpa.ne = ne
            cpa.hr_list = hr_list
            rpoints = cpa.hr_list[1].rpoints
            cpa.ir0 = cpa.hr_list[1].ir0
            rweight = cpa.hr_list[1].weight
            cpa.num_wann = cpa.hr_list[1].num_wann
            cpa.cpa_hr = HR(rpoints,cpa.ir0,zeros(Complex{T},size(cpa.hr_list[1].hopping))
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
            init_cpa_hamiltonian!(cpa)
            cpa.cpa_hkm = Hamkmesh(kmesh,cpa.cpa_hr)
            cpa
        end
    end
    
    function spir_instance_make(temperature,wmax)
        temperature_spir.Temperature(temperature,wmax)
    end

    function loop(cpa::CPA;max_niter=100,cpa_thr=10^(-8))
        for i in 1:max_niter
            calc_green_function(cpa)
            tmat_list = [calc_tmat(cpa.pot_list[i]) for i in 1:size(cpa.pot_list)[1]]
            cpa.tmat_iwn[:] = sum(cpa.weight_list.*tmat_list)
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


    function init_cpa_hamiltonian!(cpa::CPA)
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

    function calc_green_func_k(cpa::CPA,k::Array)
        hk = ftransform(cpa.cpa_hr,k)
        E = Matrix{Complex}(I,cpa.num_wann,cpa.num_wann)
        hk_mu = hk - cpa.mu.*E
        @. inv(1.0im*cpa.temp.omegaF_sp*[E] - [hk_mu] - cpa.pot_iwn)
    end

    function calc_tr_gkiwn(cpa::CPA,k::Array)
        tr.(calc_green_func_k(cpa,k))
    end

    function calc_ne(cpa::CPA,mu)
        cpa.mu = mu
        trg = apply_func_ksum(cpa.cpa_hkm,k -> calc_tr_gkiwn(cpa,k))/cpa.cpa_hkm.nk
        gl = cpa.temp.smpl_F.fit(trg)
        Ultau = -cpa.temp.basis_F.u(cpa.temp.beta)
        ne = real(dot(gl,Ultau))
    end

    function calc_green_function(cpa::CPA)
        dmu = 30
        f = (x -> calc_ne(cpa,x) - cpa.ne)
        cpa.mu = find_zero(f,(-dmu+cpa.mu,cpa.mu+dmu),Bisection())
        cpa.green[:] = apply_func_ksum(cpa.cpa_hkm,calc_green_func_k)/cpa.cpa_hkm.nk
    end
    
    function calc_tmat(cpa::CPA,pot)
        E = Matrix{Complex}(I,cpa.num_wann,cpa.num_wann)
        potd = [Diagonal(pot)]
        v = potd .- cpa.pot_iwn
        invtmp = inv.([E] .- cpa.green.*v)
        v.*invtmp
    end

    function calc_new_pot(cpa::CPA)
        tmp = cpa.green .* cpa.tmat_iwn
        tmp = inv.([E] .+ tmp)
        cpa.pot_iwn[:] += cpa.tmat_iwn.*tmp 
    end

end



import .hr
import .hamkmesh
import .cpa
using LinearAlgebra
using DoubleFloats

hr_Fe = hr.readfile("/home/shota/Nevanlinna/nev_test/Fe_hr.dat")
hr_Co = hr.readfile("/home/shota/Nevanlinna/nev_test/Co_hr.dat")
hr.ftransform(hr_Fe,[0.0,0.0,0.0])
kmesh = [8,8,8]
E= Matrix{ComplexF64}(I,16,16)
hkm =  hamkmesh.Hamkmesh(kmesh,hr_Fe)
f = x -> zeros(typeof(x[1]),16,16)
@time hamkmesh.apply_func_ksum(hkm,f)
@time hamkmesh.apply_func_ksum(hkm,f)
#=
hr_list = [hr_Fe hr_Co]
weightlist = [1.0 ,0.0]
onsite_pot_list = [29.84,24.95]

T=300
kmesh = [8,8,8]
test = cpa.CPA(T,kmesh,8.0,hr_list,weightlist,onsite_pot_list)
cpa.loop(test)
=#