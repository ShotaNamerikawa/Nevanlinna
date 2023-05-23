mutable struct NevanlinnaData{T}
    imag_points::Array{T}
    matsu_green::Array{T}
    thetas::Array{T}
    Schur_param::Array{T}
    Pick_num::Int64
    N_imag::Int64
    phis::Array{T}
    function NevanlinnaData(imag_points::Array{S}, 
                            matsu_green::Array{U}; 
                            type::Type = Complex{BigFloat}) where {S<:Complex,U<:Complex}
        nevdata = new{type}()
        nevdata.imag_points = imag_points
        nevdata.matsu_green = matsu_green
        calc_thetas!(nevdata)
        set_Pick_num!(nevdata)
        nevdata
    end
end

