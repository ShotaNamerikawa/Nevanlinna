module TemperatureSpir
using Printf
using SparseIR
using PhysicalConstants.CODATA2014

Kelvin2eV = BoltzmannConstant.val / ElementaryCharge.val

mutable struct Temperature_spir
    temerature_in_eV::Real
    beta::Real
    basis_F
    basis_B
    smpl_F
    smpl_B
    omegaF_sp::Vector
    omegaB_sp::Vector
    smpl_tau_F
    smpl_tau_B
    smpl_tau_F_on_B

    function Temperature_spir(temperature, wmax;epsilon=nothing, IR_eps=1e-8, use_fermionic_basis=true)
    
        temperature_in_eV = temperature * Kelvin2eV
        beta = 1/temperature_in_eV

        work_dtype = Float64
        basis_set = SparseIR.FiniteTempBasisSet(beta,wmax,epsilon)
        basis_F = basis_set.basis_f
        @printf("# for sparse_ir\n")
        @printf("# T = %f K\n",temperature)
        @printf("# wmax = %f eV\n",wmax)
        @printf("# irFdim = %d\n",size(basis_F.s)[1])

        # calculate sampling matsubara points
        smpl_F = basis_set.smpl_wn_f
        omegaF_sp = smpl_F.sampling_points * pi / beta
        @printf("# sampling Matsubara points = %d\n",size(smpl_F.sampling_points)[1])

        # calculate sampling tau points
        smpl_tau_F = basis_set.smpl_tau_f

        basis_B = basis_set.basis_b

        smpl_B = basis_set.smpl_wn_b
        omegaB_sp = smpl_B.sampling_points * pi / beta

        # calculate sampling tau points
        smpl_tau_B = basis_set.smpl_tau_b

        smpl_tau_F_on_B = SparseIR.TauSampling(
            basis_F,
            smpl_tau_B.sampling_points)

        new(temperature_in_eV,beta,basis_F,basis_B,smpl_F,smpl_B,omegaF_sp,omegaB_sp,
            smpl_tau_F,smpl_tau_B,smpl_tau_F_on_B)
    end

end

end