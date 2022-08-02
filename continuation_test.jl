using Plots
using DoubleFloats
using MultiFloats
include("temperature_spir.jl")
include("nevanlinna.jl")
include("spectrum_to_green.jl")
function green(z,ek::AbstractFloat,gamma::AbstractFloat)
    1/(z - ek + 1.0im*gamma)
end

function nev_gen(T,gfunc,tempspir)
    pos_omega = Complex{T}.(1.0im.*tempspir.omegaF_sp[Int(size(tempspir.omegaF_sp,1)/2)+1:end])
    nev_value = zeros(Complex{T},size(pos_omega,1))
    for x in 1:size(pos_omega,1)
        nev_value[x] = -Complex{T}(gfunc(pos_omega[x]))
    end
    nevanlinna.Nevfunc(pos_omega,nev_value,order="asc")
end    

function compare(func,nev,omega_range,delta;title::String="")
    Plots.plot((x->-1/pi*imag(func(x+1.0im*delta))),-omega_range,omega_range;label="true",title=title)
    Plots.plot!((x->nevanlinna.spectral_derive(nev,x)),-omega_range,omega_range;label="interpolate")
end

function compare_spectrum(true_func,nev,omega_range;title::String="")
    Plots.plot(true_func,-omega_range,omega_range;label="true",title=title)
    Plots.plot!((x->nevanlinna.spectral_derive(nev,x)),-omega_range,omega_range;label="interpolate")
end

function compare_spectrum_boson(true_func,nev,omega_range;title::String="",ylims=nothing,
    optimize_cal=false)
    Plots.plot(true_func;xlims=(-omega_range,omega_range),ylims=ylims,label="true",title=title)
    Plots.plot!((x->nevanlinna.boson_spectral_derive(nev,x));xlims=(-omega_range,omega_range),
    ylims=ylims,label="interpolate")
    if optimize_cal == true
        Plots.plot!((x->nevanlinna.boson_spectral_hardy_derive(nev,x));xlims=(-omega_range,omega_range),
        ylims = ylims,label="optimized")
    end
end

function compare_spectrum_boson_delta(true_func,nev,omega_range;title::String="",ylims=nothing,
    optimize_cal=false,start=0.1,stop=1,step=0.1)
    Plots.plot(true_func;xlims=(-omega_range,omega_range),ylims=ylims,label="true",title=title)
    for delta in range(start,stop,step=step)
        Plots.plot!((x->nevanlinna.boson_spectral_derive(nev,x,delta));xlims=(-omega_range,omega_range),
        ylims=ylims,label="delta = $(delta)")
        if optimize_cal == true
            Plots.plot!((x->nevanlinna.boson_spectral_hardy_derive(nev,x));xlims=(-omega_range,omega_range),
            ylims = ylims,label="optimized")
        end
    end
end

function compare_spectrum_boson_delta(nev,omega_range;title::String="",ylims=nothing,
    optimize_cal=false,start=0.1,stop=1,step=0.1)
    flag = true
    for delta in range(start,stop,step=step)
        println("delta")
        println(delta)
        println("flag")
        println(flag)
        if flag == true
            Plots.plot((x->nevanlinna.boson_spectral_derive(nev,x,delta));xlims=(-omega_range,omega_range),
        ylims=ylims,label="delta = $(delta)")
            flag = false
        elseif flag == false
            Plots.plot!((x->nevanlinna.boson_spectral_derive(nev,x,delta));xlims=(-omega_range,omega_range),
            ylims=ylims,label="delta = $(delta)")
            if optimize_cal == true
                Plots.plot!((x->nevanlinna.boson_spectral_hardy_derive(nev,x));xlims=(-omega_range,omega_range),
                ylims = ylims,label="optimized")
            end
        end
    end
end


tempspir = TemperatureSpir.Temperature_spir(300.0,10.0)
real_points = collect(range(-10,10,6000))
delta = 0.001

#=
#simple Lorentzian
gamma = 1.0
lorentz(z) = green(z,0.0,gamma)
nev_func = nev_gen(Float64,lorentz,tempspir)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,delta,25)
compare(lorentz,nevan,10,0.001;title="Lorentz Float64")
png("Lorentz_gamma_1_Float64.png")
nev_func = nev_gen(BigFloat,lorentz,tempspir)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,delta,25)
compare(lorentz,nevan,10,0.001;title="Lorentz 128 bits")
png("Lorentz_gamma_1_128.png")
#Plots.plot((x->nevanlinna.spectral_derive(nevan,x)),-10,10)
=#

#=
#simple Lorentzian
lorentz(z) = 1/3*green(z,0.0,0.1) + 1/3*green(z,2.0,0.5) + 1/3*green(z,-4.0,0.2)
nev_func = nev_gen(Float64,lorentz,tempspir)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,delta,25)
compare(lorentz,nevan,10,0.001;title="3peek Lorentz Float64")
png("3peek_Lorentz_gamma_1_Float64.png")
nev_func = nev_gen(BigFloat,lorentz,tempspir)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,delta,25)
compare(lorentz,nevan,10,0.001;title="3peek Lorentz 128 bits")
png("3peek_Lorentz_gamma_1_128.png")
#Plots.plot((x->nevanlinna.spectral_derive(nevan,x)),-10,10)
=#

#=
#Co_up
nevan = nevanlinna.Nevanlinna("trG_up_Fe_0_Co_10.dat",real_points,0.1,25)
dos(x) = nevanlinna.spectral_derive(nevan,x)
Plots.plot(dos,-10,10;label="Double64")
nevan = nevanlinna.Nevanlinna("trG_up_Fe_0_Co_10.dat",real_points,0.1,25,T=BigFloat)
dos(x) = nevanlinna.spectral_derive(nevan,x)
Plots.plot!(dos,-10,10;label="128bits")
nevan = nevanlinna.Nevanlinna("trG_up_Fe_0_Co_10.dat",real_points,0.1,25,T=Float64x4)
dos(x) = nevanlinna.spectral_derive(nevan,x)
Plots.plot!(dos,-10,10;label="Float64x4")
=#

function gauss(x,mu,sigma)
    1/(sqrt(2*pi)*sigma)*exp(-1/2*(x-mu)^2/sigma^2)
end

#=
wmin = -100.0
wmax = 100.0
gauss_green = GreeN.Green((x->gauss(x,0.0,1.0)),300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(gauss_green)
nev_func= nevanlinna.Nevfunc(1.0im*pos_omega,-temp_green)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum((x->gauss(x,0.0,1.0)),nevan,10;title="gaussian sigma = 1.0 Float64")
png("gaussian_mu_0_simga_1_Float64.png")

setprecision(256)

wmin = -BigFloat(500.0)
wmax = BigFloat(500.0)
gauss_green = GreeN.Green((x->gauss(x,0.0,1.0)),300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(gauss_green)
nev_func= nevanlinna.Nevfunc(1.0im*BigFloat.(pos_omega),-temp_green)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum((x->gauss(x,0.0,1.0)),nevan,10;title="gaussian sigma = 1.0 256 bits")
png("gaussian_mu_0_simga_1_256_wmax_500.png")
=#

#=
wmin = -Float64x8(10.0)
wmax = Float64x8(10.0)
gauss_green = GreeN.Green((x->gauss(x,0.0,1.0)),300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(gauss_green)
nev_func= nevanlinna.Nevfunc(1.0im*Float64x8.(pos_omega),-temp_green)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum((x->gauss(x,0.0,1.0)),nevan,10;title="gaussian sigma = 1.0 Float64x8")
png("gaussian_mu_0_simga_1_Float64x8.png")
=#

function lorentzian(x,mu,gamma)
    1/pi*gamma/((x-mu)^2 + gamma^2)
end

#=
wmin = -100.0
wmax = 100.0
lorentz_green = GreeN.Green((x->lorentzian(x,0.0,1.0)),300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(lorentz_green)
nev_func= nevanlinna.Nevfunc(1.0im*pos_omega,-temp_green)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum((x->lorentzian(x,0.0,1.0)),nevan,10;title="lorentz gamma = 1.0 Float64")
png("lorentzian_mu_0_gamma_1_Float64.png")

wmin = -Double64(100.0)
wmax = Double64(100.0)
lorentz_green = GreeN.Green((x->lorentzian(x,0.0,1.0)),300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(lorentz_green)
nev_func= nevanlinna.Nevfunc(1.0im*Double64.(pos_omega),-temp_green)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum((x->lorentzian(x,0.0,1.0)),nevan,10;title="lorentz gamma = 1.0 Double64")
png("lorentzian_mu_0_gamma_1_Double64.png")
=#

#=
function mtanh2(x)
    2*sinh(x)/cosh(x)^3
end

wmin = -Double64(50.0)
wmax = Double64(50.0)
mtanh2_green = GreeN.Green(mtanh2,300.0,wmin,wmax)
pos_omega,temp_green = GreeN.temp_green_derive(mtanh2_green)
nev_value = -1.0im*pos_omega.*temp_green
nev_func= nevanlinna.Nevfunc(1.0im*Double64.(pos_omega),nev_value)
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
compare_spectrum_boson_delta(mtanh2,nevan,0.1;title="mtanh2 Double64",ylims=(-0.3,0.3),start=0.01,stop=0.1,step=0.01)
png("mtanh2_Double64_omega_0.1_delta_compare_2.png")
=#

#=
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,0.1;title="mtanh2 Double64",ylims=(-0.3,0.3),optimize_cal=true)
png("mtanh2_Double64_delta_0.001_omega_0.1.png")

nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.01,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,0.1;title="mtanh2 Double64",ylims=(-0.3,0.3),optimize_cal=true)
png("mtanh2_Double64_delta_0.01_omega_0.1.png")

nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.1,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,0.1;title="mtanh2 Double64",ylims=(-0.3,0.3),optimize_cal=true)
png("mtanh2_Double64_delta_0.1_omega_0.1.png")
=#



#=
nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.001,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,10;title="mtanh2 Double64",ylims=(-1,1),optimize_cal=true)
png("mtanh2_Double64_delta_0.001_omega_10.png")

nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.01,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,10;title="mtanh2 Double64",ylims=(-1,1),optimize_cal=true)
png("mtanh2_Double64_delta_0.01_omega_10.png")

nevan = nevanlinna.Nevanlinna(nev_func,real_points,0.1,25)
iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson(mtanh2,nevan,10;title="mtanh2 Double64",ylims=(-1,1),optimize_cal=true)
png("mtanh2_Double64_delta_0.1_omega_10.png")
=#

min = -1
max = 1
num = 1000
wmin = -100.0
wmax = 100.0
real_points = collect(range(min,max;length=num))
nevan = nevanlinna.Nevanlinna("Si_correlate_green_tau_100.dat",real_points,0.01,25;statistics="B")
#iinfo = nevanlinna.Integralinfo(nevan,-50.0,50.0,10^(-4))
#nevanlinna.optimize_spectrum!(iinfo,statistics="Bose",iterations=30)
compare_spectrum_boson_delta(nevan,10;title="current-current correlation Double64",start=0.1
,stop=1,step=0.1)
png("Si_current_current_Double64_delta_compare_omega_10.png")
compare_spectrum_boson_delta(nevan,0.1;title="current-current correlation Double64",start=0.1
,stop=1,step=0.1)
png("Si_current_current_Double64_delta_compare_omega_0.1.png")
compare_spectrum_boson_delta(nevan,0.001;title="current-current correlation Double64",start=0.001
,stop=0.01,step=0.001)
png("Si_current_current_Double64_delta_0.01_compare_omega_0.001.png")
