module nevanlinna

    using LinearAlgebra
    using QuadGK
    using ForwardDiff
    using Zygote
    using Optim
    using LineSearches
    using Random
    using SparseArrays
    using DoubleFloats
    using TOML
 
    linesearch = LineSearches.BackTracking(order=2)
    initialalpha = LineSearches.InitialStatic(scaled=true)
    #linesearch = LineSearches.Static()

    jacobian = ForwardDiff.jacobian
    derivative = ForwardDiff.derivative
    seconddr(f,x) = jacobian(y -> jacobian(f,y),x)

    const prec_type = Double64        

    E = Matrix{Complex}(I,2,2)

    mutable struct Nevfunc{T<:Complex}
        imag_points::Vector{T}
        nev_value::Vector{T}
        function Nevfunc(imag_points::Vector{S},nev_value::Vector{S}
            ;order::String = "asc") where {S<:Complex}
            if order == "asc"
                imag_points = reverse(imag_points)
                nev_value = reverse(nev_value)
            elseif order == "dsc"
                imag_points = imag_points
                nev_value = nev_value    
            end
            new{typeof(imag_points[1])}(imag_points,nev_value)
        end
    end

    function readgreenfile(greenfile::String;T=Double64,order::String = "asc",statistics::String="F")
        fp = open(greenfile,"r") 
        lines = readlines(fp)
        close(fp)
        imag_points = Complex[]
        nev_value = Complex[]
        for x in lines
            line = split(x," ")
            matsubara = T(parse(Float64,line[1]))
            G_real = T(parse(Float64,line[2]))
            G_imag = T(parse(Float64,line[3]))
            if matsubara > 0
                push!(imag_points,1im*matsubara)
                if statistics == "F"
                    push!(nev_value,-G_real-1im*G_imag)
                elseif statistics == "B"
                    push!(nev_value,-1.0im*matsubara*(G_real+1im*G_imag))
                    #push!(nev_value,-(1/-1.0im*matsubara)*(G_real+1im*G_imag))
                end
            end            
        end
        Nevfunc(imag_points,nev_value;order = "asc")
    end

    mutable struct Nevanlinna{T<:Number}
        imag_points::Vector{T}
        nev_value::Vector{T}
        p_num::Int64
        real_points::Vector{Float64}
        delta::Float64
        k_num::Int64
        hardy_coeff::Vector{Float64}
        thetas::Vector{T}
        phis::Vector{T}
        function Nevanlinna(nev_func::Nevfunc,real_points::Vector{Float64}
            ,delta::Float64,k_num::Int64)
            imag_points = nev_func.imag_points
            nev_value = nev_func.nev_value 
            p_num = size(imag_points)[1]
            real_points = real_points
            delta = delta
            k_num = k_num
            hardy_coeff = zeros(Float64,4*k_num)
            thetas = zeros(Complex,size(nev_func.imag_points)[1]) 
            phis = zeros(Complex,size(nev_func.imag_points)[1])
            nev = new{typeof(imag_points[1])}(imag_points,nev_value,p_num,real_points,delta,k_num,hardy_coeff,thetas,phis)
            thetas_derive(nev)
            check_psness(nev)
            phis_derive(nev)
            nev
        end
    end

    function Nevanlinna(greenfile::String,real_points::Vector{Float64},delta::Float64
        ,k_num::Int64;T=Double64,order::String="asc",statistics="F")
        nev_func = readgreenfile(greenfile;T=T,order=order,statistics=statistics)
        Nevanlinna(nev_func,real_points,delta,k_num)
    end

    function readparameterfile(file::String)
        prinln("not implemented")
    end

    mutable struct Integralinfo
        nev::Nevanlinna
        wmin::Union{Float64,Double64}
        wmax::Union{Float64,Double64}
        lamb::Float64
        spectrum::AbstractFloat
        spectrumsum::AbstractFloat
        function Integralinfo(nev,wmin,wmax,lamb::Float64;spectrumsum=1.0)
            nev = nev
            wmin = wmin
            wmax = wmax
            lamb = lamb
            spectrum = Double64(0)
            spectrumsum=spectrumsum
            println("spectum_sum")
            println(spectrumsum)
            new(nev,wmin,wmax,lamb,spectrum,spectrumsum)
        end
    end

    function mebius(arg::T) where{T<:Complex}
        (arg - T(1.0im))/(arg + T(1.0im))
    end

    function thetas_derive(nev::Nevanlinna)
        for x in 1:nev.p_num
            nev.thetas[x] = mebius(nev.nev_value[x])
        end
    end

    function pick_derive(nev::Nevanlinna{T}) where T 
        pick = zeros(T,nev.p_num,nev.p_num);
        for i in nev.p_num
            for j in nev.p_num
                pick[i,j] = (1 - nev.thetas[i]*conj(nev.thetas[j]))/(1 - mebius(nev.imag_points[i])*conj(mebius(nev.imag_points[j])));
            end
        end
        return pick
    end

    function check_psness(nev::Nevanlinna)
        pick = pick_derive(nev::Nevanlinna)
        try 
            cholesky(pick)
            println("Pick matrix is positive semi definite.")
        catch
            println("Pick matrix is NOT positive semi definite!")
        end
    end
        

    function phis_derive(nev::Nevanlinna)
        nev.phis[1] = nev.thetas[1]
        for x in 2:nev.p_num
            abcds = E
            for y in 1:x-1
                a = (nev.imag_points[x] - nev.imag_points[y])/(nev.imag_points[x] - conj(nev.imag_points[y]))
                b = nev.phis[y]
                c = a*conj(b)
                abcds = abcds*([a b; c 1.0])
            end
            nev.phis[x] = (-abcds[2,2]*nev.thetas[x] + abcds[1,2])/(abcds[2,1]*nev.thetas[x] - abcds[1,1])
        end
    end

    function abcds_derive(nev::Nevanlinna,z::Complex) 
        abcds = E
        for x in 1:nev.p_num
            a = (z - nev.imag_points[x])/(z - conj(nev.imag_points[x]))
            b = nev.phis[x]
            c = a*conj(b)
            abcds*= ([a b; c 1.0])
        end
        return abcds
    end

    function thetalast(nev::Nevanlinna,z::T) where{T<:Complex} 
        return 0
    end

    function hardyderive(z::Complex,k::Int64)
        value = 1/(sqrt(pi)*(z + 1.0im))*((z - 1.0im)/(z + 1.0im))^k
    end

    function hardyderive(z::T,k::Int64) where{T<:Number}
        value = 1/(sqrt(pi)*(z + 1.0im))*((z - 1.0im)/(z + 1.0im))^k
    end

    function hardysumderive(nev::Nevanlinna,z::T) where{T<:Number}
        hardyvalues = 0
        for k in 0:nev.k_num-1
            hardyvalues += (nev.hardy_coeff[k+1]+1.0im*nev.hardy_coeff[k+1+nev.k_num])*hardyderive(z,k) 
            hardyvalues += (nev.hardy_coeff[k+1+2*nev.k_num]
                           +1.0im*nev.hardy_coeff[k+1+3*nev.k_num])*conj(hardyderive(z,k))
        end
        hardyvalues
    end

    function interpolate(nev::Nevanlinna,z::T) where{T<:Number}
        abcds = abcds_derive(nev,z)
        value = (abcds[1,1]*thetalast(nev,z)+abcds[1,2])/(abcds[2,1]*thetalast(nev,z)+abcds[2,2])
        last_value = 1.0im*(1 + value)/(1 - value)
    end

    function interpolate_hardy(nev::Nevanlinna,z::T) where{T<:Number}
        abcds = abcds_derive(nev,z)
        value = (abcds[1,1]*hardysumderive(nev,z)+abcds[1,2])/(abcds[2,1]
        *hardysumderive(nev,z)+abcds[2,2])
        last_value = 1.0im*(1 + value)/(1 - value)
    end

    function interpolate_hardy(nev::Nevanlinna,z::T
        ,abcds::Vector{Complex}) where{T<:Complex}
        value = (abcds[1,1]*hardysumderive(nev,z)+abcds[1,2])/(abcds[2,1]
        *hardysumderive(nev,z)+abcds[2,2])
        last_value = 1.0im*(1 + value)/(1 - value)
    end

    function interpolate_univ(nev::Nevanlinna,z::T
        ,func) where{T<:Complex}
        abcds = abcds_derive(nev,z)
        value = (abcds[1,1]*func(z)+abcds[1,2])/(abcds[2,1]
        *func(z)+abcds[2,2])
        last_value = 1.0im*(1 + value)/(1 - value)
    end
        
    function spectral_derive(nev::Nevanlinna,x::T) where{T<:Number}
        1/pi*imag(interpolate(nev,x + nev.delta*1im))  
    end

    function spectral_hardy_derive(nev::Nevanlinna,x::T) where{T<:Number}
        1/pi*imag(interpolate_hardy(nev,x + nev.delta*1.0im))  
    end

    function spectral_derive_univ(nev::Nevanlinna,x::T,func) where{T<:Number}
        1/pi*imag(interpolate_univ(nev,x+nev.delta*1.0im,func))
    end

    #=
    function spectral_derive_para(nev::Nevanlinna,x::T) where{T<:Number}
        abcds = abcds_derive(nev,x+1.0im*nev.delta)
        coeff = 1/pi*imag(2.0im*(abcds[1,1]*abcds[2,2] - abcds[1,2]*abcds[2,1])
        /((abcds[2,1] - abcds[1,1])*hardysumderive(nev,x+1.0im*nev.delta)
        +abcds[2,2] -abcds[1,2] )^2)
        vector = [hardyderive(x+1.0im*nev.delta,k) for k in 0:nev.k_num - 1]
        vector = vcat(vector,[conj(hardyderive(x+1.0im*nev.delta,k)) for k in 0:nev.k_num - 1])
        coeff*vector
    end=#
    

    function coeff_derive(nev::Nevanlinna,x::T) where{T<:Number}
        abcds = abcds_derive(nev,x+1.0im*nev.delta)
        coeff = 1/pi*(2.0im*(abcds[1,1]*abcds[2,2] - abcds[1,2]*abcds[2,1])
        /((abcds[2,1] - abcds[1,1])*hardysumderive(nev,x+1.0im*nev.delta)
        +abcds[2,2] -abcds[1,2] )^2)
    end

    function spectral_derive_para(nev::Nevanlinna,x::T) where{T<:AbstractArray}
        coeff = coeff_derive(nev,x[1])
        vector = [hardyderive(x[1]+1.0im*nev.delta,k) for k in 0:nev.k_num - 1]
        vector = vcat(vector,[1.0im*hardyderive(x[1]+1.0im*nev.delta,k) for k in 0:nev.k_num - 1])
        vector = vcat(vector,[conj(hardyderive(x[1]+1.0im*nev.delta,k)) for k in 0:nev.k_num - 1])
        vector = vcat(vector,[1.0im*conj(hardyderive(x[1]+1.0im*nev.delta,k)) for k in 0:nev.k_num - 1])
        imag(coeff*vector)
    end

    function boson_spectral_derive(nev::Nevanlinna,x::Float64)
        1/pi*imag(1/(x+nev.delta*1im)*interpolate(nev,x + nev.delta*1im))
    end

    function boson_spectral_derive(nev::Nevanlinna,x::Float64,delta::Number)
        1/pi*imag(1/(x+delta*1im)*interpolate(nev,x + delta*1im))
    end

    function boson_spectral_hardy_derive(nev::Nevanlinna,x::Float64)
        1/pi*imag(1/(x+nev.delta*1im)*interpolate_hardy(nev,x + nev.delta*1im))
    end

    function spectral_derive_second(nev::Nevanlinna,x::T) where{T<:Number}
        derivative(y -> derivative(z -> 1/pi*imag(interpolate_hardy(nev,z+1.0im*nev.delta)),y),x)
    end

    function spectral_derive_para_second(nev::Nevanlinna,x::T) where{T<:Number}
        seconddr(y->spectral_derive_para(nev,y),[x])
    end

    function write_fermion_spectrum(nev::Nevanlinna,outfile)
        println("export Fermion spectrum")
        spectral_value = zeros(Float64,size(nev.real_points)[1])
        for x in 1:size(nev.real_points)[1]
            spectral_value[x] = spectral_derive(nev,nev.real_points[x])
        end
        open(outfile,"w") do fp
            for x in 1:size(nev.real_points)[1]
                println(fp,string(nev.real_points[x])*" "*string(spectral_value[x]))
            end
        end
    end

    function write_fermion_spectrum_hardy(nev::Nevanlinna,outfile)
        println("export Fermion spectrum")
        spectral_value = zeros(Float64,size(nev.real_points)[1])
        for x in 1:size(nev.real_points)[1]
            spectral_value[x] = spectral_derive_univ(nev,nev.real_points[x],
            y -> hardysumderive(nev,y))
        end
        open(outfile,"w") do fp
            for x in 1:size(nev.real_points)[1]
                println(fp,string(nev.real_points[x])*" "*string(spectral_value[x]))
            end
        end
    end


    function write_boson_spectrum(nev::Nevanlinna,outfile::String;hardy=true)
        println("export Boson spectrum")
        spectral_value = zeros(Float64,size(nev.real_points)[1])
        if hardy==false
            for x in 1:size(nev.real_points)[1]
                spectral_value[x] = boson_spectral_derive(nev,nev.real_points[x])
            end
        elseif hardy==true
            for x in 1:size(nev.real_points)[1]
                spectral_value[x] = boson_spectral_hardy_derive(nev,nev.real_points[x])
            end
        end
        open(outfile,"w") do fp
            for x in 1:size(nev.real_points)[1]
                println(fp,string(nev.real_points[x])*" "*string(spectral_value[x]))
            end
        end
    end

    function set_norm_conserve(iinfo::Integralinfo)
        f = x -> spectral_hardy_derive(iinfo.nev,x)
        (quadgk(f
        ,iinfo.wmin,iinfo.wmax)[1] - iinfo.spectrumsum)^2
    end

    function set_norm_conserve(iinfo::Integralinfo,f)
        iinfo.spectrum = quadgk(f,iinfo.wmin,iinfo.wmax)[1]
        (iinfo.spectrum - iinfo.spectrumsum)^2
    end

    function set_norm_conserve_vec(iinfo::Integralinfo,f)
        2*(iinfo.spectrum - iinfo.spectrumsum)*quadgk(f,iinfo.wmin,iinfo.wmax)[1]
    end

    function set_norm_conserve_vec(iinfo::Integralinfo,f,g)
        2*(quadgk(f,iinfo.wmin,iinfo.wmax)[1] - iinfo.spectrumsum)*quadgk(g,iinfo.wmax,iinfo.wmax)[1]
    end

    function refrain_curveture(iinfo::Integralinfo)
        quadgk(x -> spectral_derive_second(iinfo.nev,x)^2
        ,iinfo.wmin,iinfo.wmax)[1]
    end

    function refrain_curveture(iinfo::Integralinfo,f)
        quadgk(f,iinfo.wmin,iinfo.wmax)[1]
    end
    
    function spectralintegral(iinfo::Integralinfo)
        set_norm_conserve(iinfo) + iinfo.lamb*refrain_curveture(iinfo)
    end

    function spectralintegral(iinfo::Integralinfo,f,g)
        println("spectral integral")
        @time value = set_norm_conserve(iinfo,f) + iinfo.lamb*refrain_curveture(iinfo,g)
        value
    end

    function spectralintegralcurvonly(iinfo::Integralinfo,f)
        println("spectral integral")
        @time value = refrain_curveture(iinfo,f)
        value
    end

    function spectralgrad!(iinfo::Integralinfo,G::Vector,f,g)
        println("spectralgrad!") 
        @time G[:] = (set_norm_conserve_vec(iinfo,f) 
               + iinfo.lamb*refrain_curveture(iinfo,g))
    end

    function spectralgrad!(iinfo::Integralinfo,G::Vector,f,g,h) 
        G[:] = (set_norm_conserve_vec(iinfo,f,g) 
               + iinfo.lamb*refrain_curveture(iinfo,h))
    end

    function spectralgradcurvonly!(iinfo::Integralinfo,G::Vector,f)
        println("spectralgrad!") 
        @time G[:] = refrain_curveture(iinfo,f)
    end

    function optimize_spectrum!(iinfo::Integralinfo;statistics="Fermi",iterations=1000,
        linesearch = LineSearches.BackTracking(order=2),
        initialalpha = LineSearches.InitialStatic(scaled=true))
        println(statistics)
        start = time()
        norm = x -> spectral_hardy_derive(iinfo.nev,x)
        curv = x -> spectral_derive_second(iinfo.nev,x)^2     
        norm_grad = x -> spectral_derive_para(iinfo.nev,[x])
        curv_grad = x -> (y = spectral_derive_para_second(iinfo.nev,x);
                          z = spectral_derive_second(iinfo.nev,x);
                          2*z.*y)
        function fermi_functional(coeff)
            iinfo.nev.hardy_coeff = coeff
            spectralintegral(iinfo,norm,curv)
        end

        function fermi_functional_grad!(G,coeff) 
            iinfo.nev.hardy_coeff =coeff
            spectralgrad!(iinfo,G,norm_grad,curv_grad)
        end

        function bose_functional(coeff)
            iinfo.nev.hardy_coeff = coeff
            spectralintegralcurvonly(iinfo,curv)
        end

        function bose_functional_grad!(G,coeff) 
            iinfo.nev.hardy_coeff =coeff
            spectralgradcurvonly!(iinfo,G,curv_grad)
        end
            
        if statistics == "Fermi"
            functional = fermi_functional
            funktional_grad! = fermi_functional_grad!
        elseif statistics == "Bose"            
            functional = bose_functional
            funktional_grad! = bose_functional_grad!
        end
        
        function cb!(os)
            over = time()
            println("callback")
            println("iteration number is")
            println(os.iteration)
            println("Current value")
            println(os.metadata["x"])
            println("time to calculate is")
            println(over - start)
            println("callback/")
            return false
        end
        
        function funktional(iinfo::Integralinfo,coeff)
            iinfo.nev.hardy_coeff = coeff
            spectralintegralcurvonly(iinfo,curv)
        end
        function funktional_grad!(iinfo::Integralinfo,G,coeff) 
            iinfo.nev.hardy_coeff =coeff
            spectralgradcurvonly!(iinfo,G,curv_grad)
        end
        
        a = funktional
        b = funktional_grad!
        #precond(n) = spdiagm(-1=> -ones(n-1),0=>2*ones(n),1=>-ones(n-1))*(n+1)
        rng = MersenneTwister(100)
        random = rand(rng,iinfo.nev.k_num*4).*0.001
        println("optimize start")
        optimize_cal!(a,b,iinfo,linesearch
        ,initialalpha,iterations)
    end

    function optimize_cal!(functional,functional_grad!,iinfo::Integralinfo
        ,linesearch,initialalpha,iterations)
        res = optimize(coeff->functional(iinfo,coeff)
        ,(G,coeff)->functional_grad!(iinfo,G,coeff),zeros(Double64,iinfo.nev.k_num*4)
        ,LBFGS(alphaguess=initialalpha,linesearch=linesearch),Optim.Options(
        iterations = iterations,time_limit=32400,show_trace=true))
        #zeros(Float64,iinfo.nev.k_num*4)
        #=(G,coeff) -> (iinfo.nev.hardy_coeff 
        =coeff;spectralgrad!(iinfo,G,norm_grad,curv_grad)),=#
        iinfo.nev.hardy_coeff = Optim.minimizer(res)
        print(res)
    end

end
#=
module operation
    #using .nevanlinna

    struct Data
        greenfile::String
    end


    function read_file(data::Data,tomlfile::String)
        conf = TOML.parsefile(tomlfile)
        greenfile = conf["greenfile"] 
        start = conf["realpoints"]["start"]
    end

    function try_read(passindex::String)
        try 
            index
        catch
        end
    end
end
=#

#=
using .nevanlinna
using Zygote
using QuadGK
using LinearAlgebra
using DoubleFloats
using MultiFloats
num = 10000
min = -1
max = 1
wmin = -100.0
wmax = 100.0
real_points = collect(range(min,max;length=num))
nev = nevanlinna.Nevanlinna("correlate_green_tau_100.dat"
,real_points,0.01,25;T=Float64x4,order="asc",statistics="B")

nevanlinna.write_boson_spectrum(nev,"current_current_correlation_tau_100_delta_001.dat";hardy=false)
=#

#iinfo = nevanlinna.Integralinfo(nev,wmin,wmax,0.0001;spectrumsum=8)

#nev.hardy_coeff[:] = [0.000 for x in 1:4*nev.k_num]
#norm = x -> nevanlinna.spectral_hardy_derive(iinfo.nev,x)
#curv = x -> nevanlinna.spectral_derive_second(iinfo.nev,x)^2     
#norm_grad = x -> nevanlinna.spectral_derive_para(iinfo.nev,[x])
#curv_grad = x -> (y = nevanlinna.spectral_derive_para_second(iinfo.nev,x);y.*y)


#nevanlinna.optimize_spectrum!(iinfo)
#nevanlinna.write_fermion_spectrum_hardy(nev,"spectrum_optimize_Fe_100_Co_0.dat")
