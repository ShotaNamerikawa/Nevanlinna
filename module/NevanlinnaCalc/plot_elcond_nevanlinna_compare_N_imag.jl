"""
check the difference of data according to N_imag
cdata::CalcData
N_imag_reduce_list::Array{Int64}
calc_type: "elcond" , "rho", "cpa"
plot_kwargs: 
"""
function compare_N_imag(cdata::CalcData,
                        N_imag_reduce_list::Array{Int64}
                        ;
                        calc_type = "elcond",
                        plot_kwargs=Dict())
    reset_N_imag!(cdata;N_imag_reduce=N_imag_reduce_list[1])
    if N_imag_reduce_list[1] < 0
        label = @sprintf("Max Pick criterion num + %d",- N_imag_reduce_list[1])
    else
        label = @sprintf("Max Pick criterion num - %d",N_imag_reduce_list[1])
    end
    default_kwargs = Dict(:legend=>:outertopright,
                          :label=>label,
                          :legendtitle=>"Sampling points num")
    plot_kwargs = merge(default_kwargs,plot_kwargs)
    println(plot_kwargs)
    if calc_type == "elcond" 
        plot_func = plot_elcond
    elseif (calc_type == "rho") || (calc_type == "cpa")
        plot_func = plot_rho
    else
        throw(ErrorException("calc_type is wrong"))
    end
    println(cdata.xdata_list)
    p = plot_func(cdata;plot_kwargs=plot_kwargs)
    for N_imag_reduce in N_imag_reduce_list[2:end]
        reset_N_imag!(cdata;N_imag_reduce=N_imag_reduce)
        if N_imag_reduce < 0
            label = @sprintf("Max Pick criterion num + %d",-N_imag_reduce)
        else
            label = @sprintf("Max Pick criterion num - %d",N_imag_reduce)
        end
        p = plot_func(p,cdata;plot_kwargs=Dict(:label=>label))        
    end
    p
end

function data_N_imag(cdata::CalcData,
                     N_imag_reduce::Int64
                     ;
                     data_type = "elcond")
    reset_N_imag!(cdata;N_imag_reduce=N_imag_reduce)
    generate_ydata(cdata;ydata = data_type) 
end