mutable struct NevanlinnaData{T}
    matsu_freq::Array{T}
    matsu_green::Array{T}
    thetas::Array{T}
    Schur_param::Array{T}
    N_imag::Int64
end

