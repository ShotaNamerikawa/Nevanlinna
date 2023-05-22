module Nevanlinna

using LinearAlgebra
using QuadGK
using ForwardDiff
using Zygote
using Optim
using LineSearches
using Random
using SparseArrays
using DoubleFloats
using MultiFloats
using TOML

export readgreenfile
export NevanlinnaData
export calc_thetas!,set_N_imag!
export calc_pick_num!,derive_pick,derive_phis!

"""
Calculate mebius transformation.
"""
function mebius(arg::T) where{T<:Complex}
    (arg - T(1.0im))/(arg + T(1.0im))
end


include("convert_green.jl")
include("Nevanlinna_data.jl")


function readparameterfile(file::String)
    println("not implemented")
end


include("Nevanlinna_pre.jl")
include("Nevanlinna_impl.jl")
include("Nevanlinna_optim.jl")

end
