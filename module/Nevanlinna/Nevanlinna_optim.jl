mutable struct Integralinfo
    nev::NevanlinnaData
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
#=
function sum_rule_cal(nev::Nevanlinna_data)
    nev.
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
=#