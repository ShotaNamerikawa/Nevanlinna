@with_kw mutable struct PostData
    xdata::String
    temperature::Float64
    prefix::String
    N_imag_reduce::Int64
    N_imag_min::Int64
    N_imag_max::Int64
    eta::Float64
    wmax::Float64
end


mutable struct PostDataList
    dir::String
    regex::Regex
    files::Array{String}
    function PostDataList(dir::String,regex::Regex)
        pdl = new()
        pdl.dir = dir
        pdl.regex = regex
        get_file!(pdl)
        pdl
    end
end

function convert_key_to_symbol(dict::Dict)
    Dict(Symbol(k) => v for (k,v) in pairs(dict))
end

function get_file!(p::PostDataList)
    p.files = get_file(p.regex,p.dir)
end

function plot_data(file::String)
    data = readdlm(file)
    plot(data[:,1],data[:,2])
end

"""
    map toml data to struct
    caution! struct has to be defined with @with_kw
"""
function map_toml_to_struct(file::String,T::Type;additional_attribute=Dict())
    toml_data = TOML.parsefile(file)
    toml_data = convert_key_to_symbol(toml_data)
    toml_data = merge(toml_data,additional_attribute)
    T(;toml_data...)
end

#function plot_data_type()

function plot_elconds(p::PostDataList)
    for (i,file) in enumerate(p.files)
        postdata = map_toml_to_struct(joinpath(p.dir,file),
                                      PostData
                                      ;
                                      additional_attribute=Dict(:prefix=>replace(file,"_info.toml" => "")))
        println("postdata")
        println(postdata)
        elcond_file = joinpath(p.dir,postdata.prefix*"_elcond.dat")
        @printf("elconf file is %s\n",elcond_file)
        elcond_data = readdlm(elcond_file)
        println("plot data is")
        println(plot_data)
        if i != 1
            p = plot(p,)
        else
            p = scatter(elcond_data[:,1],elcond_data[:,2];label =@sprintf("N_imag:%s",postdata.N_imag_reduce))
        end
    end
    return p
end