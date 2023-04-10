function cal_N_imags!(cdata::CalcData)
    cdata.N_imag_list = Int64[]
    for ndata in cdata.ndatas
        append!(cdata.N_imag_list,ndata.N_imag)
    end
end

function cal_static_conductivities!(cdata::CalcData)
    cdata.elcond_list = Float64[]
    for ndata in cdata.ndatas
        elcond = cal_static_conductivity(ndata)
        append!(cdata.elcond_list,elcond)
    end
end

function generate_plot_attributes(cdata::CalcData;ydata="elcond")
    plot_attribute = Dict()
    if cdata.xdata == "mu"
        plot_attribute[:xticks] = (
            range(minimum(cdata.xdata_list),maximum(cdata.xdata_list);step=1)
            )
        plot_attribute[:xlabel] = L"\mu\ \mathrm{[eV]}"
    elseif cdata.xdata == "concentration"
        plot_attribute[:xticks] = (
            range(minimum(cdata.xdata_list),maximum(cdata.xdata_list);step=0.1)
            )
        plot_attribute[:xlabel] = latexstring(@sprintf("\\mathrm{%s}_x",cdata.element))
    end
    if ydata == "elcond"
        plot_attribute[:ylabel] = L"\sigma_{xx}\ \mathrm{[S/m]}"
    elseif ydata == "rho"
        plot_attribute[:ylabel] = L"\rho_{xx} \ \mathrm{[\mu \Omega\cdot cm]}"
    elseif ydata == "N_imag"
        plot_attribute[:ylabel] = L"\mathrm{\mathrm{The number of Sampling points}}"
        plot_attribute[:yticks] = (
        vcat(0,collect(range(minimum(cdata.N_imag_list),maximum(cdata.N_imag_list))))
        )
        plot_attribute[:ylim] = (0,maximum(cdata.N_imag_list)+1)

    elseif ydata == "Pick"
        plot_attribute[:ylabel] = L"\mathrm{\mathrm{Max number matching Pick criterion }}"
        plot_attribute[:yticks] = (
        vcat(0,collect(range(minimum(cdata.N_imag_list .+ cdata.elcond_kwargs[:N_imag_reduce])
                            ,maximum(cdata.N_imag_list .+ cdata.elcond_kwargs[:N_imag_reduce]))))
        )
        plot_attribute[:ylim] = (0,maximum(cdata.N_imag_list .+ cdata.elcond_kwargs[:N_imag_reduce])+1)
        
    end
    plot_attribute
end

function generate_ydata(cdata::CalcData;ydata)
    if ydata=="elcond"
        return cdata.elcond_list
    elseif ydata == "rho"
        return 10^8 ./cdata.elcond_list
    elseif ydata == "N_imag"
        return cdata.N_imag_list
    elseif ydata == "Pick"
        return cdata.N_imag_list .+ cdata.elcond_kwargs[:N_imag_reduce]
    else
        throw(ErrorException("No ydata"))
    end
end

function plot_data(cdata::CalcData
                   ;
                   ydata="elcond",
                   plot_kwargs = Dict())
    plot_attributes = generate_plot_attributes(cdata;ydata=ydata)    
    ydata_list = generate_ydata(cdata;ydata=ydata)
    plot(cdata.xdata_list,
        ydata_list
        ;
        seriestype=:scatter,
        legend=false,
        plot_attributes...,
        plot_kwargs...)
end

function plot_data(p::Plots.Plot,cdata::CalcData
    ;
    ydata="elcond",
    plot_kwargs::Dict = Dict())
    ydata_list = generate_ydata(cdata;ydata=ydata)
    scatter(p,cdata.xdata_list,
        ydata_list
        ;
        plot_kwargs...)
end

function plot_elcond(cdata::CalcData;plot_kwargs::Dict=Dict())
    plot_data(cdata;ydata="elcond",plot_kwargs=plot_kwargs)
end

function plot_elcond(p,cdata::CalcData;plot_kwargs::Dict=Dict())
    plot_data(p,cdata;ydata="elcond",plot_kwargs=plot_kwargs)
end

function plot_rho(cdata::CalcData;plot_kwargs = Dict())
    plot_data(cdata;ydata="rho",plot_kwargs=plot_kwargs)
end

function plot_rho(p,cdata::CalcData;plot_kwargs = Dict())
    plot_data(p,cdata;ydata="rho",plot_kwargs=plot_kwargs)
end

function plot_N_imag(cdata::CalcData,plot_kwargs=Dict())
    plot_data(cdata;ydata="N_imag",plot_kwargs=plot_kwargs)
end

function plot_N_imag(p,cdata::CalcData,plot_kwargs=Dict())
    plot_data(p,cdata;ydata="N_imag",plot_kwargs=plot_kwargs)
end

function plot_Pick_criterion(cdata::CalcData,plot_kwargs=Dict())
    plot_data(cdata;ydata="Pick",plot_kwargs=plot_kwargs)
end

function plot_Pick_criterion(p,cdata::CalcData,plot_kwargs=Dict())
    plot_data(p,cdata;ydata="Pick",plot_kwargs=plot_kwargs)
end


