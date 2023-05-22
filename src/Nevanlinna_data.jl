mutable struct NevanlinnaData{T}
    imag_points::Array{T}
    matsu_green::Array{T}
    thetas::Array{T}
    Schur_param::Array{T}
    Pick_num::Int64
    N_imag::Int64
    function NevanlinnaData(imag_points::Array{T}, matsu_green::Array{T}) where T
        nevdata = new{T}()
        nevdata.imag_points = imag_points
        nevdata.matsu_green = matsu_green
        calc_thetas!(nevdata)
        nevdata
    end
end

