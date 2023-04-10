module Pade
export PadeData, interpolate
using LinearAlgebra
mutable struct PadeData{T}
    sample_x::Vector{Complex{T}}
    sample_func::Vector{Complex{T}}
    interpolate_x::Vector{Complex{T}}
    function PadeData(sample_x,sample_func,interpolate_x;T::Type = BigFloat) 
        println("type is")
        println(T)
        pd = new{T}()
        pd.sample_x = Complex{T}.(sample_x)
        pd.sample_func = Complex{T}.(sample_func)
        pd.interpolate_x = Complex{T}.(interpolate_x)
        pd
    end
end

function interpolate(pd::PadeData{T}) where T
    ndim = size(pd.sample_func,1)
    g = Array{Complex{T}}(undef,ndim,ndim)
    an = Array{Complex{T}}(undef,ndim)
    g[1:end,1] .= pd.sample_func
    for n in range(1,ndim-1)
        g[n+1:end,n+1] .= (g[n,n] .- g[n+1:end,n])./(pd.sample_x[n+1:end] .- pd.sample_x[n]) ./g[n+1:end,n]
    end
    an .= diag(g)
    p = zeros(Complex{T},size(an,1)+1,size(pd.interpolate_x,1))
    q = zeros(Complex{T},size(an,1)+1,size(pd.interpolate_x,1))
    q[1,1:end] .= 1
    q[2,1:end] .= 1
    p[1,1:end] .= 0
    p[2,2:end] .= an[1]
    for n in range(2,size(an,1))
        p[n+1,1:end] .= p[n,1:end] .+ an[n] * (pd.interpolate_x[1:end] .- pd.sample_x[n-1]) .* p[n-1,1:end]
        q[n+1,1:end] .= q[n,1:end] .+ an[n] * (pd.interpolate_x[1:end] .- pd.sample_x[n-1]) .* q[n-1,1:end]
    end
    return p[end,1:end] ./ q[end,1:end]
end

end

