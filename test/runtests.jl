using Test
using QuadGK
using SparseIR
using PhysicalConstants.CODATA2014

k_B_in_eV = 1/ElementaryCharge.val*BoltzmannConstant.val

mutable struct test_param{T}
    temperature::Float64
    omega_range::Array{T}
end

function get_inv_temp(temperature)
    1/(k_B_in_eV*temperature)
end

function generate_spir_data(temperature,wmax)
    β = get_inv_temp(temperature)
    basis = FiniteTempBasis{Fermionic}(β,wmax,nothing)
    siw = MatsubaraSampling(basis; positive_only=true)
    matsu_freq = Array{ComplexF64}(undef,size(siw.sampling_points,1))
    
    for i in eachindex(matsu_freq)
        matsu_freq[i] = SparseIR.valueim(siw.sampling_points[i],β)
    end
    
    return (basis,siw,matsu_freq)
end

function fermi_kernel(x,freq)
    1/(freq-x)
end

function convert_spectral_to_matsu_green(spectral_func::Function, 
                                         matsu_freq,
                                         min_freq,
                                         max_freq
                                         ;
                                         kernel::Function=fermi_kernel)
    integral,error = quadgk(x -> kernel(x,matsu_freq)*spectral_func(x),min_freq,max_freq)    
end

function generate_green_from_spectral(spectral_func,
                                      matsu_freqs::Array{T}
                                      ;
                                      freq_window = (Big(-100.0),Big(100.0))  
                                      ) where T
    
    green = Array{T}(undef,size(matsu_freqs,1))
    
    for i in eachindex(green) 
        green[i] = convert_spectral_to_matsu_green(spectral_func,
                                                   matsu_freqs[i],
                                                   freq_window[1],
                                                   freq_window[2]
                                                   )[1]        
    end

    return green
end

function check_causality_green(green::Array{Complex{T}}) where T
    causality = true
    for i in eachindex(green)
        if imag(green[i]) > 0
            causality = false
        end        
    end
    return causality
end

function check_Schur_param_abs(Schur_param::Array{T};verbose = false) where T
    
    islessone = true

    for i in eachindex(Schur_param)
        abs_sp = abs(Schur_param[i])
        print("index:")
        print(i)
        print(" abs_value:")
        println(abs_sp)
        if abs_sp >= 1 
            islessone = false
            if verbose == true
                println("Absolute value is larger than 1. Something is wrong!")
            end
        end
    end
    
    return islessone
end

@testset "Nevanlinna.jl" begin
    """
    test for gaussian continuation by Nevanlinna AC
    """
    function gauss_test(mean,var,matsu_freq,wmax)
        function gaussian(arg)
            return 1/sqrt(2*pi*var^2)*exp(-(arg-mean)^2/(2*var^2))
        end
        green = generate_green_from_spectral(gaussian,
                                             matsu_freq
                                             ;
                                             freq_window = [-wmax,wmax]
                                            )                                           

        @testset "causality" begin 
            @test check_causality_green(green)
        end        
        
        
    end

    #test parameter
    mean = 0.0
    var = 1.0
    wmax = 10.0
    temperature = 300.0

    matsu_freq = generate_spir_data(temperature,wmax)[3]
    
    gauss_test(mean,var,matsu_freq,wmax)

end

nothing
