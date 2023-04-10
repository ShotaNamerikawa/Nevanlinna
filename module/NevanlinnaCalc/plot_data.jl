function cal_spec_fermion(ndata::NevanlinnaRealData{T}) where T
    imag.(ndata.real_data.val)./pi
end

function cal_spec_boson(ndata::NevanlinnaRealData{T}) where T
    imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq)) ./pi
end

function cal_opt_conductivity(ndata::NevanlinnaRealData{T}) where T
    hbar_in_eV.*imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq))./real.(ndata.real_data.freq)
end

function cal_opt_resistance(ndata::NevanlinnaRealData{T}) where T
    10^8 ./(hbar_in_eV.*imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq))./real.(ndata.real_data.freq))
end

function min_ind_cal(ndata::NevanlinnaRealData)
    min_ind = argmin(abs.(real.(ndata.real_data.freq)))
    if in(0.0,ndata.real_data.freq) == true
        min_ind += 1
    end
    min_ind
end

function cal_static_conductivity(ndata::NevanlinnaRealData{T}) where T
    min_ind = min_ind_cal(ndata)
    return cal_opt_conductivity(ndata)[min_ind]
end

function cal_static_resistance(ndata::NevanlinnaRealData{T}) where T
    #min_ind = min_ind_cal(ndata)
    min_ind = Int64((size(ndata.real_data.freq,1)+mod(size(ndata.real_data.freq,1),2))/2)
    return cal_opt_resistance(ndata)[min_ind]
end

function cal_spectral(ndata::NevanlinnaRealData{T};type="fermion") where T
    if type == "fermion"
        return cal_spec_fermion(ndata)
    elseif type == "boson"
        return cal_spec_boson(ndata)
    elseif type == "conductivity"
        return cal_opt_conductivty(ndata)
    elseif type == "resistance"
        return cal_opt_resistance(ndata)
    end
end

function plot_spec(ndata::NevanlinnaRealData{T};type="fermion",kwargs...) where T
    plot(real.(ndata.real_data.freq),cal_spectral(ndata,type=type);kwargs...)
end

function plot_spec(p,ndata::NevanlinnaRealData{T};type="fermion",kwargs...) where T
    plot!(p,real.(ndata.real_data.freq),cal_spectral(ndata,type=type);kwargs...)
end

function plot_spec(ndatas::Vector{NevanlinnaRealData},individual_args;type="fermion",kwargs...) 
    p=plot_spec(ndatas[1];type=type,individual_args[1]...,Dict(kwargs)...)
    for i in axes(ndatas,1)
        if i != 1
            plot_spec(p,ndatas[i];type=type,individual_args[i]...)
        end
    end
    p
end

function plot_spec(ndatas::Vector{NevanlinnaRealData};type="fermion",kwargs...) 
    p=plot_spec(ndatas[1];type=type,Dict(kwargs)...)
    for i in axes(ndatas,1)
        if i != 1
            plot_spec(p,ndatas[i];type=type,Dict(kwargs)...)
        end
    end
    p
end

function plot_spec_fermion(ndata::NevanlinnaRealData{T};kwargs...) where T
    plot(real.(ndata.real_data.freq), imag.(ndata.real_data.val)./pi;kwargs...)
end

function plot_spec_fermion(ndatas::Vector{NevanlinnaRealData};kwargs...) 
    p=plot(real.(ndatas[1].real_data.freq), imag.(ndatas[1].real_data.val)./pi;Dict(kwargs)...)
    for i in axes(ndatas,1)
        if i != 1
            plot!(p,real.(ndatas[i].real_data.freq), imag.(ndatas[i].real_data.val)./pi)
        end
    end
    p
end

function plot_spec_boson(ndata::NevanlinnaRealData{T};kwargs...) where T 
    plot(real.(ndata.real_data.freq), imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq));kwargs...)
end

function plot_spec_boson_conductivity(ndata::NevanlinnaRealData{T};kwargs...) where T
    plot(real.(ndata.real_data.freq), hbar_in_eV.*imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq))./real.(ndata.real_data.freq);kwargs...)
end

function plot_spec_boson_resistance(ndata::NevanlinnaRealData{T};kwargs...) where T
    plot(real.(ndata.real_data.freq), 10^8 ./(hbar_in_eV.*imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq))./real.(ndata.real_data.freq));kwargs...)
end

function plot_spec_resistance(ndata::NevanlinnaRealData{T};kwargs...) where T
    plot(real.(ndata.real_data.freq), 10^8 ./(hbar_in_eV.*imag.(ndata.real_data.val).*tanh.(0.5*ndata.beta*real.(ndata.real_data.freq))./real.(ndata.real_data.freq));kwargs...)
end

function get_sample_and_interpolate(ndata::NevanlinnaComplexData{T}) where T
    return (ndata.matsu_omega,ndata.green_iwn,ndata.complex_data.freq,ndata.complex_data.val)
end

function plot_interpolate_check(ndata::NevanlinnaComplexData{T};scaling_kbt = false) where T
    data = get_sample_and_interpolate(ndata)
    p = scatter(data[1],imag.(data[2]),label="Sampling points")
    return plot(p,imag.(ndata.complex_data.freq),imag.(ndata.complex_data.val),label="Interpolated value")
end