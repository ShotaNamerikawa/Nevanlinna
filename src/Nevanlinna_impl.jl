function set_phis!(nev::NevanlinnaData{T}) where T
    nev.phis = Array{T}(undef,nev.N_imag)
end


"""
Calculate phis' value from sampling points value.
"""
function derive_phis!(nev::NevanlinnaData{T}) where T
    set_phis!(nev)
    nev.phis[1] = nev.thetas[1]
    tmp_mat = ones(T,2,2)
    for x in 2:nev.N_imag
        abcds = Array{T}(I,2,2)
        for y in 1:x-1
            @views tmp_mat[1,1] = (nev.imag_points[x] - nev.imag_points[y])/(nev.imag_points[x] - conj(nev.imag_points[y]))
            @views tmp_mat[1,2] = nev.phis[y]
            @views tmp_mat[2,1] = tmp_mat[1,1]*conj(tmp_mat[1,2])
            @views abcds *= tmp_mat
        end
        nev.phis[x] = (-abcds[2,2]*nev.thetas[x] + abcds[1,2])/(abcds[2,1]*nev.thetas[x] - abcds[1,1])
    end
end

function derive_abcds(nev::NevanlinnaData{T},z) where{T} 
    abcds = Array{T}(I,2,2)
    tmp_mat = ones(T,2,2)
    for x in 1:nev.N_imag
        tmp_mat[1,1] = (z - nev.imag_points[x])/(z - conj(nev.imag_points[x]))
        tmp_mat[1,2] = nev.phis[x]
        tmp_mat[2,1] = tmp_mat[1,1]*conj(tmp_mat[1,2])
        abcds*= tmp_mat
    end
    return abcds
end

function thetalast(nev::NevanlinnaData,z) 
    return 0
end

function interpolate_nev_func_last_0(nev::NevanlinnaData{T},z) where T
    abcds = derive_abcds(nev,z)
    value = abcds[1, 2] / abcds[2, 2]
    return 1.0im*(1 + value)/(1 - value)
end

function interpolate_nev_func_theta_last(thetalast::Function,nev::NevanlinnaData,z::T) where{T<:Number}
    abcds = derive_abcds(nev,z)
    value = (abcds[1,1]*thetalast(nev,z)+abcds[1,2])/(abcds[2,1]*thetalast(nev,z)+abcds[2,2])
    last_value = 1.0im*(1 + value)/(1 - value)
end


function interpolate_nev_func_hardy(nev::NevanlinnaData,z::T) where{T<:Number}
    abcds = derive_abcds(nev,z)
    value = (abcds[1,1]*hardysumderive(nev,z)+abcds[1,2])/(abcds[2,1]
    *hardysumderive(nev,z)+abcds[2,2])
    last_value = 1.0im*(1 + value)/(1 - value)
end

function interpolate_nev_func_hardy(nev::NevanlinnaData,z::T
    ,abcds::Vector{Complex}) where{T<:Complex}
    value = (abcds[1,1]*hardysumderive(nev,z)+abcds[1,2])/(abcds[2,1]
    *hardysumderive(nev,z)+abcds[2,2])
    last_value = 1.0im*(1 + value)/(1 - value)
end

function interpolate_nev_func_gen_func(nev::NevanlinnaData{T},z
    ,func) where {T}
    abcds = derive_abcds(nev,z)
    value = (abcds[1,1]*func(z)+abcds[1,2])/(abcds[2,1]
    *func(z)+abcds[2,2])
    last_value = 1.0im*(1 + value)/(1 - value)
end

    
function derive_spectral(nev::NevanlinnaData{T},x::S) where{T,S<:AbstractFloat}
    1/pi*imag(interpolate_nev_func_theta_last(nev,x + nev.delta*1im))  
end


function derive_spectral_hardy(nev::NevanlinnaData,x::T) where{T<:Number}
    1/pi*imag(interpolate_nev_func_hardy(nev,x + nev.delta*1.0im))  
end

function derive_spectral_univ(nev::NevanlinnaData,x::T,func) where{T<:Number}
    1/pi*imag(interpolate_nev_func_gen_func(nev,x+nev.delta*1.0im,func))
end

function derive_coeff(nev::NevanlinnaData,x::T) where{T<:Number}
    abcds = derive_abcds(nev,x+1.0im*nev.delta)
    coeff = 1/pi*(2.0im*(abcds[1,1]*abcds[2,2] - abcds[1,2]*abcds[2,1])
    /((abcds[2,1] - abcds[1,1])*hardysumderive(nev,x+1.0im*nev.delta)
    +abcds[2,2] -abcds[1,2] )^2)
end

function derive_spectral_para(nev::NevanlinnaData,x::T) where{T<:AbstractArray}
    coeff = derive_coeff(nev,x[1])
    vector = [hardyderive(x[1]+1.0im*nev.delta,k) for k in 0:nev.k_num - 1]
    vector = vcat(vector,[1.0im*hardyderive(x[1]+1.0im*nev.delta,k) for k in 0:nev.k_num - 1])
    vector = vcat(vector,[conj(hardyderive(x[1]+1.0im*nev.delta,k)) for k in 0:nev.k_num - 1])
    vector = vcat(vector,[1.0im*conj(hardyderive(x[1]+1.0im*nev.delta,k)) for k in 0:nev.k_num - 1])
    imag(coeff*vector)
end

