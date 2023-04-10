function calc_pre_calculation(nev::NevanlinnaData{T}) where T
    calc_pick_denominator(nev)
end

function calc_pick_num!(nev::NevanlinnaData{T}) where T
    for i in 1:nev.p_num
        @views nev.pick_mat[1:i,i] .= (1.0 .- nev.thetas[1:i].*conj(nev.thetas[i])) ./ nev.pick_denominator[1:i,i]
        if issuccess(cholesky!(Hermitian(nev.pick_mat[1:i,1:i]),check=false)) == false
            nev.pick_num = i-1
            break
        end
    end
end

function set_N_imag!(nev::NevanlinnaData{T}) where{T}
    nev.N_imag = nev.pick_num
end

function set_N_imag!(nev::NevanlinnaData{T},N_imag::Int64) where{T}
    nev.N_imag = N_imag
end

"""
Caution! This function has to be called for a given green function sample
only at once!
"""
function fit_paramater_to_pick_num!(nev::NevanlinnaData{T}) where T
    nev.imag_points = nev.smpl_imag_points[nev.N_imag:-1:1]
    nev.thetas[1:nev.N_imag] .= @view(nev.thetas[nev.N_imag:-1:1])
end

function derive_pick(nev::NevanlinnaData{T}) where T 
    pick = zeros(Complex{T},nev.p_num,nev.p_num)
    for i in 1:nev.p_num
        for j in 1:nev.p_num
            pick[i,j] = (1 - nev.thetas[i]*conj(nev.thetas[j]))/(1 - mebius(nev.imag_points[i])*conj(mebius(nev.imag_points[j])));
        end
    end
    return pick
end

function derive_pick!(nev::NevanlinnaData{T}) where T 
    for i in 1:nev.p_num
        for j in 1:nev.p_num
            nev.pick[i,j] = (1 - nev.thetas[i]*conj(nev.thetas[j]))/(1 - mebius(nev.imag_points[i])*conj(mebius(nev.imag_points[j])));
        end
    end
    return pick
end

"""
Check Pick criterion for Nevanlinna function.
"""
function check_psness(nev::NevanlinnaData)
    pick = derive_pick(nev::NevanlinnaData)
    try 
        cholesky(pick)
        println("Pick matrix is positive semi definite.")
    catch
        println("Pick matrix is NOT positive semi definite!")
    end
end


function set_phis!(nev::NevanlinnaData{T}) where T
    nev.phis = Array{Complex{T}}(undef,nev.N_imag)
end


"""
Calculate phis' value from sampling points value.
"""
function derive_phis!(nev::NevanlinnaData{T}) where T
    nev.phis[1] = nev.thetas[1]
    tmp_mat = ones(Complex{T},2,2)
    for x in 2:nev.N_imag
        abcds = nev.E
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
    abcds = nev.E
    tmp_mat = ones(Complex{T},2,2)
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

function interpolate_nev_func_theta_last(nev::NevanlinnaData,z::T) where{T<:Number}
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

