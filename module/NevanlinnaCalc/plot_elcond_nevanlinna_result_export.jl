"""
export datas calculated from cdata to files
"""
function export_data(cdata::CalcData,prefix::String)
    f = open(@sprintf("%s_info.toml",prefix),"w")
    data = Dict(
        "xdata" => cdata.xdata,
        "temperature" => cdata.temperature,
        "N_imag_reduce" => cdata.elcond_kwargs[:N_imag_reduce],
        "N_imag_min" => minimum(cdata.N_imag_list),
        "N_imag_max" => maximum(cdata.N_imag_list),
        "eta" => cdata.elcond_kwargs[:eta],
        "wmax" => cdata.ndatas[1].omega_max
    )
    TOML.print(f,data)
    close(f)
    g = open(@sprintf("%s_elcond.dat",prefix),"w")
    writedlm(g,[cdata.xdata_list cdata.elcond_list])
    close(g)
    h = open(@sprintf("%s_N_imag.dat",prefix),"w")
    writedlm(h,[cdata.xdata_list cdata.N_imag_list])
    close(h)
end

function export_all(cdata::CalcData,prefix::String;plot_type,data_output=true,plot_kwargs=Dict())
    if plot_type == "elcond"
        p = plot_elcond(cdata;plot_kwargs=plot_kwargs)
        savefig(p,@sprintf("%s_plot_elcond.pdf",prefix))
    elseif plot_type == "rho" 
        p = plot_rho(cdata;plot_kwargs=plot_kwargs)
        savefig(p,@sprintf("%s_plot_rho.pdf",prefix))
    #elseif plot_type == "cpa"
    #    p = plot_cpa_rho()
    else
        throw(ErrorException("plot_type is elcond or rho"))
    end
    if data_output == true
        export_data(cdata,prefix)
    end
    return p
end