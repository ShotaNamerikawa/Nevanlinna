include("./Hamiltonian.jl")
include("./temperature_spir.jl")

module TransPort
using ..Hamiltonian
using Roots
using LinearAlgebra
using QuadGK
using PhysicalConstants.CODATA2014
using Printf
using TensorOperations
using ..TemperatureSpir

ele = ElementaryCharge.val
Kelvin2eV = BoltzmannConstant.val / ElementaryCharge.val 
hbar = PlanckConstantOver2pi #in [J*s]
hbar_in_eV = PlanckConstantOver2pi.val / ElementaryCharge.val #in [eV*s]

nontransform = Matrix{Float64}(I,2,2)

mutable struct Transport{T}
    ir::TemperatureSpir.Temperature_spir
    u_sample::Array
    u_minus::Array
    hr::Hamiltonian.HR.Hr{T}
    hkm::Hamiltonian.HamKmesh.Hamkmesh
    transform_mat::Array #column lattice vectors represented by cart
    temperature::Float64
    energy::Array{Float64}
    ne::Real
    mu::Real
    delta::Real
    tau::Real
    constscat::Real
    E::Matrix{Complex{T}}
    function Transport(hr::Hamiltonian.HR.Hr{T},kmesh::Vector{Int64},temperature,energy,ne
        ;delta=0.1,mu=0.0,tau=100*10^(-15),transform_mat=nontransform) where T
        ir = TemperatureSpir.Temperature_spir(temperature,100.0)
        u_sample = ir.basis_F.u(ir.smpl_tau_F.sampling_points)
        u_minus = -ir.basis_F.u(ir.beta .- ir.smpl_tau_F.sampling_points)
        hkm = Hamiltonian.HamKmesh.Hamkmesh(kmesh,hr)
        mu = mu
        constscat = hbar_in_eV/(2*tau)
        E = Matrix{Complex}(I,hr.num_wann,hr.num_wann)
        new{T}(ir,u_sample,u_minus,hr,hkm,transform_mat,
        temperature,energy,ne,mu,delta,tau,constscat,E)
    end
end

function Transport(hr::Hamiltonian.HR.Hr{T},kmesh::Vector{Int64},temperature
    ;mu=0.0,tau=100*10^(-15),transform_mat=nontransform) where T
    Transport(hr,kmesh,temperature,undef,undef;
    mu=mu,tau=tau,transform_mat=transform_mat)
end

function calc_green_func_k(tp::Transport{T},k::Array) where{T}
    hk = Hamiltonian.HR.ftransform(tp.hr,k)
    E = Matrix{Complex}(I,tp.hr.num_wann,tp.hr.num_wann)
    hk_mu = [hk - tp.mu.*E]
    @. inv(1.0im*(tp.ir.omegaF_sp+sign(tp.ir.omegaF_sp)*tp.constscat )*[E] - hk_mu)
end

function calc_green_func_k!(tp::Transport{T},k::Array,gien) where T
    hk = Hamiltonian.HR.ftransform(tp.hr,k)
    E = Matrix{Complex}(I,tp.hr.num_wann,tp.hr.num_wann)
    for i in 1:size(tp.ir.omegaF_sp,1)
        gien[i,:,:] = inv((1.0im*(tp.ir.omegaF_sp[i] + sign(tp.ir.omegaF_sp[i])*tp.constscat)
                      + tp.mu)*E - hk)
    end
end

function calc_diagonalized_green_func_k!(tp::Transport{T},k::Array,gien) where T
    hk = Hamiltonian.HR.ftransform(tp.hr,k)
    ek,U = eigen(hk)
    v = Hamiltonian.HR.generate_hkwa(tp.hr,k)
    diagv = zeros(Complex,tp.hr.num_wann,3)
    for l in 1:3
        v[:,:,l] = U'*v[:,:,l]*U
        diagv[:,l] = diag(v[:,:,l])
    end
    for i in 1:size(tp.ir.omegaF_sp,1)
        gien[i,:] = 1 ./ ((1.0im*(tp.ir.omegaF_sp[i] + sign(tp.ir.omegaF_sp[i])*tp.constscat)
                      + tp.mu) .- ek)
    end
    diagv
end

function calc_green_func_k_2!(tp::Transport{T},k::Array,omega_lambda,gien,giwl) where T
    hk = Hamiltonian.HR.ftransform(tp.hr,k)
    E = Matrix{Complex}(I,tp.hr.num_wann,tp.hr.num_wann)
    hk_mu = [hk - tp.mu.*E]
    all_points = zeros(2*size(tp.ir.omegaF_sp,1))
    all_points[1:Int(size(tp.ir.omegaF_sp,1))] = tp.ir.omegaF_sp
    all_points[Int(size(tp.ir.omegaF_sp,1))+1:end] = tp.ir.omegaF_sp .+ omega_lambda
    gien .= @. inv(1.0im*(all_points + sign(all_points)*tp.constscat)*[E] - hk_mu)
    giwl .= @. inv(1.0im*(all_points - omega_lambda + 
       sign(all_points - omega_lambda)*tp.constscat)*[E] - hk_mu)   
end

function calc_green_k(tp::Transport,k::Array,omega)
    hk = Hamiltonian.HR.ftransform(tp.hr,k)
    E = Matrix{Complex}(I,tp.hr.num_wann,tp.hr.num_wann)
    hk_mu = hk - tp.mu*E
    green = inv((omega+1.0im*tp.delta)*E - hk_mu)
end

function calc_green_k!(green,tp::Transport,k::Array,omega)
    green[:,:] += inv((omega+1.0im*tp.delta+ tp.mu)*tp.E - Hamiltonian.HR.ftransform(tp.hr,k))
end

function calc_green(tp::Transport,omega)
    result = zeros(Complex,tp.hr.num_wann,tp.hr.num_wann)
    apply_func(k) = calc_green_k(tp,k,omega)
    Hamiltonian.HamKmesh.apply_func_ksum!(apply_func,result,tp.hkm)    
    result
end
#=
function calc_green(tp::Transport,omega)
    green = zeros(Complex,tp.hr.num_wann,tp.hr.num_wann)
    apply_func!(k) = calc_green_k!(green,tp,k,omega)
    Hamiltonian.HamKmesh.apply_func_ksum(apply_func!,tp.hkm)    
    green
end
=#
function calc_dos(tp,omega)
    -1/pi/tp.hkm.nk*imag(tr(calc_green(tp,omega)))
end

function calc_spin_dos(tp::Transport,omega)
    @time green = calc_green(tp,omega)
    dos_up = 0
    dos_down = 0
    for i in 1:tp.hr.num_wann
        if i % 2 != 0
            dos_up += -1/pi/tp.hkm.nk*imag(green[i,i])
        elseif i % 2 == 0
            dos_down += -1/pi/tp.hkm.nk*imag(green[i,i])
        end
    end
    dos_up,dos_down
end

function calc_ne(tp::Transport,mu)
    tp.mu = mu
    value = quadgk((omega -> 1/(exp(omega/(Kelvin2eV*tp.temperature))+1)
    *calc_dos(tp,omega)),-100.0,100.0)
    println("dos")
    println(value)
    value[1]
end

function calc_mu(tp::Transport)
    dmu = 30
    tp.mu = find_zero(mu->calc_ne(tp,mu)-tp.ne,(tp.mu-dmu,tp.mu+dmu),Bisection())
end

function write_dos(tp::Transport{T},up_file::String,dn_file::String) where T
    up_fp = open(up_file,"w")
    dn_fp = open(dn_file,"w")
    println(size(tp.energy)[1])
    for x in 1:size(tp.energy)[1]
        dos_up,dos_down = calc_spin_dos(tp,tp.energy[x])
        write(up_fp,@sprintf("%.16f %.16f\n",tp.energy[x],dos_up))
        write(dn_fp,@sprintf("%.16f %.16f\n",tp.energy[x],dos_down))
    end
end

function write_total_dos(tp::Transport{T},file::String) where T
    fp = open(file,"w")
    for x in 1:size(tp.energy)[1]
        dos = calc_dos(tp,tp.energy[x])
        write(fp,@sprintf("%.16f %.16f\n",tp.energy[x],dos))
    end
end


function cal_correlation_k(tp::Transport{T},k::Array,v,omega_lambda) where T
    result = zeros(Complex{T},3)
    gien = [zeros(Complex{T},tp.hr.num_wann,tp.hr.num_wann) for i in 1:2*size(tp.ir.omegaF_sp,1)]
    giwl = [zeros(Complex{T},tp.hr.num_wann,tp.hr.num_wann) for i in 1:2*size(tp.ir.omegaF_sp,1)]
    calc_green_func_k_2!(tp,k,omega_lambda,gien,giwl) 
    for i in 1:3
        result[i] = sum(tr.([v[:,:,i]].*gien
                    .*[v[:,:,i]].*giwl))
    end
    result
end

function cal_correlation_k(tp::Transport{T},k::Array,v) where T
    result = zeros(Complex{T},size(tp.ir.smpl_tau_B.sampling_points,1),3)
    gien = zeros(Complex{T},size(tp.ir.omegaF_sp,1),tp.hr.num_wann,tp.hr.num_wann) #gien=G(k,ien)
    gtau = zeros(Complex{T},size(tp.ir.smpl_tau_F.sampling_points,1),tp.hr.num_wann,tp.hr.num_wann) 
    #gtau=G(k,tau)
    gtau_minus = zeros(Complex{T},size(tp.ir.smpl_tau_F.sampling_points,1)
                       ,tp.hr.num_wann,tp.hr.num_wann) #gtau_minus=G(k,-tau)
    calc_green_func_k!(tp,k,gien) #G(k,ien) on sampling Matsubara points is calculated.
    gl = TemperatureSpir.fit(tp.ir.smpl_F,gien)#G(k,ien) -> G(k)_l
    #@tensoropt gtau[m,i,j] = tp.u_sample[l,m]*gl[l,i,j] #G(k)_l -> G(k,tau) (tau > 0)
    gtau = TemperatureSpir.evaluate(tp.ir.smpl_tau_F,gl) #G(k)_l -> G(k,tau) (tau > 0)
    #@tensoropt gtau_minus[m,i,j] = tp.u_minus[l,m]*gl[l,i,j] #G(k)_l -> G(k,-tau)
    gtau_minus = -1 .* TemperatureSpir.evaluate(tp.ir.smpl_tau_F,gl)[end:-1:1,:,:] #G(k)_l -> G(k,-tau)
    for x in 1:3
        for tau in 1:size(tp.ir.smpl_tau_B.sampling_points,1)
            sumval = 0
            @tensoropt sumval = v[:,:,x][j,n]*gtau[tau,:,:][n,l]*v[:,:,x][l,m]*gtau_minus[
                                 tau,:,:][m,j]
            result[tau,x] = sumval
        end
    end 
    result
end

function cal_diagonalized_correlation_k(tp::Transport{T},k::Array) where T
    result = zeros(Complex{T},size(tp.ir.smpl_tau_B.sampling_points,1),3)
    gien = zeros(Complex{T},size(tp.ir.omegaF_sp,1),tp.hr.num_wann) #gien=G(k,ien)
    gtau = zeros(Complex{T},size(tp.ir.smpl_tau_F.sampling_points,1),tp.hr.num_wann) 
    #gtau=G(k,tau)
    gtau_minus = zeros(Complex{T},size(tp.ir.smpl_tau_F.sampling_points,1)
                       ,tp.hr.num_wann) #gtau_minus=G(k,-tau)
    v = calc_diagonalized_green_func_k!(tp,k,gien) #G(k,ien) on sampling Matsubara points is calculated.
    println("v value")
    println(v[:,1])
    gl = TemperatureSpir.fit(tp.ir.smpl_F,gien)#G(k,ien) -> G(k)_l
    #@tensoropt gtau[m,i,j] = tp.u_sample[l,m]*gl[l,i,j] #G(k)_l -> G(k,tau) (tau > 0)
    gtau = TemperatureSpir.evaluate(tp.ir.smpl_tau_F,gl) #G(k)_l -> G(k,tau) (tau > 0)
    #@tensoropt gtau_minus[m,i,j] = tp.u_minus[l,m]*gl[l,i,j] #G(k)_l -> G(k,-tau)
    gtau_minus = -1 .* TemperatureSpir.evaluate(tp.ir.smpl_tau_F,gl)[end:-1:1,:] #G(k)_l -> G(k,-tau)
    for x in 1:3
        for tau in 1:size(tp.ir.smpl_tau_B.sampling_points,1)
            sumval = sum(v[:,x].*gtau[tau,:].*v[:,x].*gtau_minus[
                                 tau,:])
            result[tau,x] = sumval
        end
    end 
    result
end

function cal_correlation_sum(tp::Transport{T},omega_lambda::Real) where T
    func_k(k) = begin 
                v = Hamiltonian.HR.generate_hkwa(tp.hr,k)
                cal_correlation_k(tp,k,v,omega_lambda)
                end
    result = zeros(Complex,3)
    Hamiltonian.HamKmesh.apply_func_ksum!(func_k,result,tp.hkm)
    -Kelvin2eV*tp.temperature*result/tp.hkm.nk
end

function cal_correlation_sum(tp::Transport{T}) where T
    func_k(k) = begin 
                v = Hamiltonian.HR.generate_hkwa(tp.hr,k)
                for i in 1:tp.hr.num_wann
                    for j in 1:tp.hr.num_wann
                        v[i,j,:] = tp.transform_mat*v[i,j,:]/hbar_in_eV
                    end
                end
                cal_correlation_k(tp,k,v)
                end
    result = zeros(Complex{T},size(tp.ir.smpl_tau_B.sampling_points,1),3)
    @time Hamiltonian.HamKmesh.apply_func_ksum!(func_k,result,tp.hkm)
    result ./ tp.hkm.nk ./ det(tp.transform_mat) # Correlation time (-1)
end

function cal_diagonalized_correlation_sum(tp::Transport{T}) where T
    func_k(k) = begin 
                cal_diagonalized_correlation_k(tp,k)
                end
    result = zeros(Complex{T},size(tp.ir.smpl_tau_B.sampling_points,1),3)
    @time Hamiltonian.HamKmesh.apply_func_ksum!(func_k,result,tp.hkm)
    result ./ tp.hkm.nk ./ det(tp.transform_mat) # do not forget multiply by hbar_in_Js
                                                    #when derive conductivity 
end

function cal_correlation_func_all(tp::Transport{T}) where T
    pos_omegaB_sp = tp.ir.omegaB_sp[Int((size(tp.ir.omegaB_sp)[1]-1)/2)+2:end]
    correlate_func = zeros(Complex,(size(pos_omegaB_sp,1),3))
    for (i,omega) in enumerate(pos_omegaB_sp)
        @time correlate_func[i,:] = cal_correlation_sum(tp,omega)
    end
    correlate_func
end

function cal_correlation_func_all_2(tp::Transport{T}) where T
    correlate_func = cal_diagonalized_correlation_sum(tp)
    gl = TemperatureSpir.fit(tp.ir.smpl_tau_B,correlate_func)
    #(transpose(tp.ir.basis_B.uhat(tp.ir.smpl_B.sampling_points))*gl)[Int((size(tp.ir.omegaB_sp,1)-1)/2)+2:end,:]
    TemperatureSpir.evaluate(tp.ir.smpl_B,gl)[Int((size(tp.ir.omegaB_sp,1)-1)/2)+2:end,:]
end

function write_correlation(tp::Transport{T},out_file::String) where T
    println("export correlation function")
    correlate_func = cal_correlation_func_all_2(tp)
    println(size(correlate_func))
    println(correlate_func[:,1])
    println(correlate_func[:,2])
    println(correlate_func[:,3])
    fp = open(out_file,"w")
    pos_omegaB_sp = tp.ir.omegaB_sp[Int((size(tp.ir.omegaB_sp)[1]-1)/2)+2:end]
    for (i,omega) in enumerate(pos_omegaB_sp)
        #print(fp,@sprintf("%.16f",omega))
        print(fp,omega)
        for l in 1:3
            @printf(fp," %.15e %.15e",real(correlate_func[i,l]),
            imag(correlate_func[i,l]))
        end
        @printf(fp,"\n")
    end
    close(fp)
end
end

module SpCond
using ..Hamiltonian
using Roots
using LinearAlgebra
using QuadGK
using PhysicalConstants.CODATA2014
using Printf
using TensorOperations
using ..TemperatureSpir

ele = ElementaryCharge.val
Kelvin2eV = BoltzmannConstant.val / ElementaryCharge.val 
hbar = PlanckConstantOver2pi.val #in [J*s]
hbar_in_eV = PlanckConstantOver2pi.val / ElementaryCharge.val #in [eV*s]

nontransform = Matrix{Float64}(I,2,2)

mutable struct Spcond{T} 
    hr::Hamiltonian.HR.Hr{T}
    hkm::Hamiltonian.HamKmesh.Hamkmesh
    transform_mat::Array #column lattice vectors represented by cart
    temperature::Float64
    energy::Array{Float64}
    ne::Real
    mu::Real
    delta::Real
    tau::Real
    constscat::Real
    E::Matrix{Complex{T}}
    function Spcond(hr::Hamiltonian.HR.Hr{T},kmesh::Vector{Int64},temperature,energy,ne
        ;delta=0.1,mu=0.0,tau=100*10^(-15),transform_mat=nontransform) where T
        hkm = Hamiltonian.HamKmesh.Hamkmesh(kmesh,hr)
        mu = mu
        constscat = hbar_in_eV/(2*tau)
        E = Matrix{Complex}(I,hr.num_wann,hr.num_wann)
        new{T}(hr,hkm,transform_mat,
        temperature,energy,ne,mu,delta,tau,constscat,E)
    end
end

function cal_spectra(sc::Spcond,k::Array)
    hk = Hamiltonian.HR.ftransform(sc.hr,k)
    vals,vecs = eigen(Hermitian(hk))
    vals,vecs
end

function cal_spectralcond_k(sc::Spcond,k::Array)
    result = zeros(Float64,size(sc.energy,1),3,3)
    vel = zeros(Complex,sc.hr.num_wann,3)
    vals,vecs = cal_spectra(sc,k)
    vel_mat = Hamiltonian.HR.generate_hkwa(sc.hr,k)
    vel = Hamiltonian.HR.generate_va(vel_mat,vecs)

    for i in 1:sc.hr.num_wann
        vel[i,:] = sc.transform_mat*vel[i,:] ./ hbar_in_eV
    end
    
    vel = real.(vel)


    for i in 1:size(sc.energy,1)
        imag_gR = -sc.constscat ./ ((sc.energy[i] .- vals) .^ 2 .+ sc.constscat^2)
        for l in 1:3
            for m in 1:3
                result[i,l,m] = sum(vel[:,l] .* vel[:,m] .* imag_gR .^ 2)
            end
        end
    end

    result
end

function cal_spectralcond_sum(sc::Spcond)
    function func(k::Array)
        cal_spectralcond_k(sc,k)
    end
    result = zeros(Float64,size(sc.energy,1),3,3)
    Hamiltonian.HamKmesh.apply_func_ksum!(func,result,sc.hkm)
    (hbar/(pi*sc.hkm.nk*det(sc.transform_mat)))*result*10^10 # [S/m]
end

function write_spectralcond(sc::Spcond,outfile::String)
    spcond = cal_spectralcond_sum(sc)
    fp = open(outfile,"w")
    for i in 1:size(sc.energy,1)
        @printf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e \n",
        sc.energy[i],spcond[i,1,1],spcond[i,1,2],spcond[i,1,3]
        ,spcond[i,2,2],spcond[i,2,3],spcond[i,3,3])
    end 
    close(fp)
end

end

module TransportfromSC

using PhysicalConstants.CODATA2014
using Printf

ele = ElementaryCharge.val
Kelvin2eV = BoltzmannConstant.val / ElementaryCharge.val 
hbar = PlanckConstantOver2pi.val #in [J*s]
hbar_in_eV = PlanckConstantOver2pi.val / ElementaryCharge.val #in [eV*s]

mutable struct TS{S,T}
    energy::Array{S}
    sc::Array{T,3}
    temperature::Float64
end

function read_file(sc_file::String,temperature;T=Float64)
    fp = open(sc_file)
    s_data = readlines(fp)
    energy = zeros(T,size(s_data,1))
    sc = zeros(T,size(s_data,1),3,3)
    for i in 1:size(s_data,1)
        line = parse.(T,split(s_data[i]))
        energy[i] = line[1]
        sc[i,1,:] = line[2:4]
        sc[i,2,2:3] = line[5:6]
        sc[i,3,3] = line[7]
        sc[i,2,1] = line[3]
        sc[i,3,1] = line[4]
        sc[i,3,2] = line[6]
    end
    close(fp)
    TS(energy,sc,temperature)    
end

"""
   return derivative of -f_FD(epsilon) = -1/(e^beta(epsilon-mu) + 1)
"""

function simpson(sample_points::Array{S},func_value::Array{T}) where {S,T}
    if size(sample_points,1)%2 != 1
        error("the number of sample_points should be odd")
    elseif size(sample_points,1) != size(func_value,1)
        error("the number of sample points and func value is different!")
    end
    n = size(sample_points,1) - 1
    h = sample_points[2] - sample_points[1]
    result = h/3*sum(func_value[1:2:n,:,:] + 4*func_value[2:2:n,:,:] + func_value[3:2:n+1,:,:],dims=1)
    result[1,:,:]
end

function minus_fermi_diff(ts::TS,energy::T,mu::S) where {T,S}
    1/(4*Kelvin2eV*ts.temperature)*(sech.((energy .- mu)/(2*Kelvin2eV*ts.temperature))).^2
end

function cal_response_func(ts::TS,mu,n::Int64)
    simpson(ts.energy,minus_fermi_diff(ts,ts.energy,mu).*(ts.energy .- mu).^n.*ts.sc) 
end

function cal_all_thermo_electric(ts::TS,mu)
    L11 = cal_response_func(ts::TS,mu,0)
    L12 = 1/ele*cal_response_func(ts::TS,mu,1)
    Seebeck = -1/ts.temperature*L12*inv(L11)
    Powerfactor = 1/ts.temperature^2*L12^2*inv(L11)
    L11,L12,Seebeck,Powerfactor
end

end
 