include("temperature_spir.jl")

module GreeN
using ..TemperatureSpir
using QuadGK


function kernel(z,omega)
    1/(z - omega)
end

mutable struct Green{T}
    spectrum
    temperature
    wmin::T
    wmax::T
    temp_spir
    function Green(spectrum,temperature,wmin::T,wmax::T) where T
        temp_spir = TemperatureSpir.Temperature_spir(temperature,wmax)
        new{T}(spectrum,temperature,wmin,wmax,temp_spir)
    end
end

function temp_green_derive(green::Green;statistics="Fermi")
    if statistics=="Fermi"
        pos_omega = green.temp_spir.omegaF_sp[Int(size(green.temp_spir.omegaF_sp,1)/2)+1:end]
    elseif statistics=="Bose"
        pos_omega = green.temp_spir.omegaF_sp[(Int(size(green.temp_spir.omegaF_sp,1)-1)/2)+2:end]
    end
    integral = zeros(Complex{typeof(green.wmax)},size(pos_omega,1))
    for x in 1:size(pos_omega,1)
        integral[x] = QuadGK.quadgk((y->kernel(1.0im*pos_omega[x],y)*green.spectrum(y)),green.wmin,green.wmax)[1]
    end
    pos_omega,integral
end

end