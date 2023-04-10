mutable struct CalcData
    dir::String
    files::Array{String}
    file_expression::Regex
    xdata::String
    temperature::Float64
    element::Union{String,Nothing}
    elcond_kwargs::Dict
    ndatas::Array{NevanlinnaData}
    xdata_list::Array{Float64}
    elcond_list::Array{Float64}
    N_imag_list::Array{Int64}
    fermi_energy::Float64
    function CalcData(dir,
                      expression,
                      temperature
                      ;
                     elcond_kwargs=Dict(),
                     xdata="mu", 
                     element = nothing,
                     fermi_energy::Float64 = 0.0
                     )
        cdata = new()
        cdata.dir = dir
        cdata.file_expression = expression
        cdata.xdata = xdata
        cdata.temperature = temperature
        cdata.element = element
        cdata.fermi_energy = fermi_energy 
        default_kwargs=Dict(:N_real => 2,
                            :omega_max => 0.00001,
                            :eta => 0.000001,
                            :N_imag_reduce => 1)
        cdata.elcond_kwargs = merge(default_kwargs,elcond_kwargs)
        cdata.files = get_file(cdata)
        get_information_from_files!(cdata)
        cal_static_conductivities!(cdata)
        cal_N_imags!(cdata)
        cdata
    end
end

"""
reset N_imag or N_imag_reduce and recalculate conductivities
"""
function reset_N_imag!(cdata::CalcData;N_imag::Union{Int64,Nothing}=nothing,N_imag_reduce::Union{Int64,Nothing}=nothing)
    if isnothing(N_imag) == false
        cdata.elcond_kwargs[:N_imag] = N_imag
    elseif isnothing(N_imag_reduce) == false
        cdata.elcond_kwargs[:N_imag] = 0
        cdata.elcond_kwargs[:N_imag_reduce] = N_imag_reduce
    else
        return 1
    end
    get_information_from_files!(cdata)
    cal_static_conductivities!(cdata)
    cal_N_imags!(cdata)
    return 0
end