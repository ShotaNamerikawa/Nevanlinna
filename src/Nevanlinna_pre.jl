"""
Transform Nevanlinna function into contractive function.
"""
function calc_thetas!(nev::NevanlinnaData{T}) where T
    nev.thetas = Array{T}(undef,size(nev.matsu_green,1))
    for i in eachindex(nev.matsu_green)
        nev.thetas[i] = mebius(- nev.matsu_green[i])       
    end
end

function calc_pick_mat(nev::NevanlinnaData{T}) where T
    pick_mat = Array{eltype(nev.imag_points)}(undef,
                                              size(nev.imag_points,1),
                                              size(nev.imag_points,1))
    for i in eachindex(nev.imag_points)
        for j in eachindex(nev.imag_points)
            pick_mat[i,j] = (1.0 - nev.thetas[i].*conj(nev.thetas[j])) / (1.0 - mebius(nev.imag_points[i])*conj(mebius(nev.imag_points[j])) )
        end
    end     
    pick_mat
end

function calc_pick_num!(nev::NevanlinnaData{T}) where T
    for i in eachindex(nev.imag_points)
        for j in 1:nev.p_num
            @views nev.pick_mat[1:i,i] .= (1.0 .- nev.thetas[1:i].*conj(nev.thetas[i])) ./ nev.pick_denominator[1:i,i]
            if issuccess(cholesky!(Hermitian(nev.pick_mat[1:i,1:i]),check=false)) == false
                nev.pick_num = i-1
                break
            end
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
function fit_input_points_to_pick_num!(nev::NevanlinnaData{T}) where T
    nev.imag_points = nev.smpl_imag_points[nev.N_imag:-1:1]
    nev.thetas[1:nev.N_imag] .= nev.thetas[nev.N_imag:-1:1]
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

