mutable struct NevanlinnaData{T}
    imag_points::Array{T}
    matsu_green::Array{T}
    thetas::Array{T}
    Pick_num::Int64
    N_imag::Int64
    phis::Array{T}
    function NevanlinnaData(imag_points::Array{S}, 
                            matsu_green::Array{U}; 
                            type::Type = Complex{BigFloat}) where {S<:Complex,U<:Complex}
        nevdata = new{type}()
        nevdata.imag_points = imag_points
        nevdata.matsu_green = matsu_green
        return nevdata
    end
end

"""
make NevanlinnaData and calculate phis
"""
function generate_NevanlinnaData(imag_points::Array{S},
                                 matsu_green::Array{U};
                                 type::Type=Complex{BigFloat},
                                 N_imag=nothing,
                                 input_order = "dsc") where {S<:Complex,U<:Complex}

    nevdata = NevanlinnaData(imag_points,
                             matsu_green,
                             ;
                             type = type)
    calc_thetas!(nevdata)
    set_Pick_num!(nevdata)
    if isnothing(N_imag) == true
        set_N_imag!(nevdata)
    else
        set_N_imag!(nevdata,N_imag)
    end
    fit_input_data_to_Pick_num!(nevdata;order = input_order)
    derive_phis!(nevdata) 
    nevdata
end


