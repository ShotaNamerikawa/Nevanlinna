@with_kw struct PostData
    temperature::Float64
    prefix::String
    N_imag_reduce::Int64
    N_imag_min::Int64
    N_imag_max::Int64
    eta::Float64
    wmax::Float64
end

struct PostData_all
    dir::String
    regex
    files::Array{String}
end

function test3(a)
   println(a)
end

function convert_key_to_symbol!(dict::Dict)
    println("aaa")
    dict = Dict(Symbol(k) => v for (k,v) in pairs(dict))
end

function get_file!(p::PostData_all)
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
function map_toml_to_struct(file::String,T::Type)
    toml_data = TOML.parsefile(file)
    convert_key_to_symbol!(toml_data)
    println("aaa")
    println(toml_data)
    T(;toml_data...)
end

function test()
    println("sxxt fxxk")
end
#=
function plot_elconds(p::PostData_all,before_plot = nothing)

end
=#  