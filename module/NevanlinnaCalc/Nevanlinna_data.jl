abstract type NevanlinnaData{T}
end

mutable struct NevanlinnaRealData{T} <: NevanlinnaData{T}
    positive_omega::Bool
    nhw::Int64
    matsu_omega::Complex{T}
    green_iwn::Complex{T}
    N_imag::Int64
    N_imag_reduce::Int64
    imag_data::Nevanlinna.ImagDomainData
    real_data::Nevanlinna.RealDomainData
    phis::Array{Complex{T}}
    abcd
end

mutable struct NevanlinnaComplexData{T} <: NevanlinnaData{T}
    positive_omega::Bool
    nhw::Int64
    matsu_omega::Complex{T}
    green_iwn::Complex{T}
    N_imag::Int64
    N_imag_reduce::Int64
    imag_data::Nevanlinna.ImagDomainData
    complex_data::Nevanlinna.ComplexDomainData
end