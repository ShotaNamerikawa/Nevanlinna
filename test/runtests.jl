using Nevanlinna
using Test
using QuadGK
using SparseIR
using PhysicalConstants.CODATA2014
using Plots
using Printf

k_B_in_eV = 1/ElementaryCharge.val*BoltzmannConstant.val

mutable struct test_param{T}
    temperature::Float64
    omega_range::Array{T}
end

function get_inv_temp(temperature)
    1/(k_B_in_eV*temperature)
end

function generate_spir_data(β,wmax)
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

"""
Back continue spectral function to its Green function
"""
function generate_green_from_spectral(spectral_func,
                                      matsu_freqs::Array{T}
                                      ;
                                      freq_window = (Big(-100.0),Big(100.0))  
                                      ) where T
    
    green = Array{Complex{eltype(freq_window)}}(undef,size(matsu_freqs,1))
    
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

function check_unity_thetas(thetas::Array{Complex{T}}) where T
    unity = true
    for i in eachindex(thetas)
        if abs(thetas[i]) > 1
            unity = false
        end
    end 
    unity
end

function check_Schur_param_abs(Schur_param::Array{T};verbose = false) where T
    
    islessone = true

    for i in eachindex(Schur_param)
        abs_sp = abs(Schur_param[i])
        if verbose == true
            print("index:")
            print(i)
            print(" abs_value:")
            @printf "%15.15f\n" abs_sp
        end
        if abs_sp >= 1 
            islessone = false
            if verbose == true
                println("Absolute value is larger than 1. Something is wrong!")
            end
        end
    end
    
    return islessone
end

function interpolate_plot(nevdata,original_rho,real_freq)
    interpolated_green = Array{ComplexF64}(undef,size(real_freq,1))

    for i in eachindex(real_freq)
        interpolated_green[i] = interpolate_nev_func_last_0(nevdata,real_freq[i])
    end

    p = plot(real_freq,original_rho.(real_freq);label="original",xlabel="energy [eV]",ylabel = "spectral [1/eV]")
    display(plot!(p,real_freq,1/pi.*imag.(interpolated_green);label="interpolation"))
end



@testset "Nevanlinna.jl" begin
    """
    test for spectral function continuation by Nevanlinna AC
    """
    function test_generate_green_from_lorentz(matsu_freq,
                                              peak,
                                              gamma
                                              ;
                                              wmax = 10.0
                                             )
        function lorentz(x)
            1/pi*gamma/((x - peak)^2 + gamma^2)
        end

        original_green_data = 1 ./ (matsu_freq .+ 1.0im*gamma)
        back_continued_green = generate_green_from_spectral(lorentz,
                                                            matsu_freq
                                                            ;
                                                            freq_window=[-big(wmax), big(wmax)]
                                                            )
        p = scatter(imag.(matsu_freq), imag.(back_continued_green); label="back continued")
        display(scatter!(p,imag.(matsu_freq), imag.(original_green_data);label = "original",markersize=2))
    end                

    function spectral_test(spectral_func,matsu_freq,wmax,real_freq;input_order="dsc")
        #get green function of gauss spectral function
        green = generate_green_from_spectral(spectral_func,
                                             matsu_freq
                                             ;
                                             freq_window = [-big(wmax),big(wmax)]
                                            )                                           

        println("eltype of green")
        println(eltype(green))
        @testset "causality" begin 
            @test check_causality_green(green)
        end                        

        nevdata = Nevanlinna.generate_NevanlinnaData(matsu_freq,green;input_order = input_order)
        @testset "unity of thetas" begin
            @test check_unity_thetas(nevdata.thetas)
        end
        
        set_N_imag!(nevdata)

        print("pick_num:")
        println(nevdata.Pick_num)


        print("N_imag:")
        println(nevdata.N_imag)

        @testset "existence of Nevanlinna func value" begin
            @test nevdata.Pick_num >= 1
        end

        @testset "phis value in disk" begin
            @test check_Schur_param_abs(nevdata.phis;verbose=true)
        end

        interpolate_plot(nevdata,gaussian,real_freq)
        
    end

    #test parameter
    mean = 0.0
    var = 1.0
    wmax = 10.0
    temperature = 300.0
    β = 100.0
    real_freq = collect(range(-10,10,1000))


    function gaussian(arg)
        return 1 / sqrt(2 * pi * var^2) * exp(-(arg - mean)^2 / (2 * var^2))
    end

    test_spectral(arg) = gaussian(arg)


    matsu_freq = generate_spir_data(β,wmax)[3]
    
    setprecision(128)
    spectral_test(test_spectral,matsu_freq,wmax,real_freq;input_order="dsc")

    #peak = 0.0
    #gamma = 1
    #test_generate_green_from_lorentz(matsu_freq,peak,gamma)
end

nothing
