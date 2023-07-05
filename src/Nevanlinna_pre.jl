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
    for i in axes(nev.imag_points,1)
        for j in axes(nev.imag_points,1)
            pick_mat[i,j] = (1.0 - nev.thetas[i].*conj(nev.thetas[j])) / (1.0 - mebius(nev.imag_points[i])*conj(mebius(nev.imag_points[j])) )
        end
    end     
    Hermitian(pick_mat)
end

function calc_Pick_num(nev::NevanlinnaData{T}) where T
    pick_mat = calc_pick_mat(nev)
    Pick_num = 0
    for i in axes(nev.imag_points,1)
        if issuccess(cholesky(pick_mat[1:i,1:i],check=false)) == false
            break
        end
        Pick_num += 1
    end    
    return Pick_num         
end

function set_Pick_num!(nev::NevanlinnaData{T}) where{T}
    nev.Pick_num = calc_Pick_num(nev)
end

function set_N_imag!(nev::NevanlinnaData{T}) where{T}
    set_Pick_num!(nev)
    nev.N_imag = nev.Pick_num
end

function set_N_imag!(nev::NevanlinnaData{T},N_imag::Int64) where{T}
    nev.N_imag = N_imag
end

"""
Caution! This function has to be called for a given green function sample
only at once!
"""
function fit_input_data_to_Pick_num!(nev::NevanlinnaData{T};order = "dsc") where T
    if order == "dsc" #Fixme! change also green function value!
        nev.imag_points = nev.imag_points[nev.N_imag:-1:1]
        nev.matsu_green = nev.matsu_green[nev.N_imag:-1:1]
        nev.thetas = nev.thetas[nev.N_imag:-1:1] 
    elseif order == "asc"
        nev.imag_points = nev.imag_points[1:nev.N_imag]
        nev.matsu_green = nev.matsu_green[nev.N_imag:-1:1]
        nev.thetas = nev.thetas[1:nev.N_imag]
    else
        throw(ErrorException("choose order value from asc or dsc!"))
    end
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

