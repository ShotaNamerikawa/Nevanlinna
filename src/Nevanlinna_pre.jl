"""
Transform Nevanlinna function into contractive function.
"""
function derive_thetas!(nev::Nevanlinna_data{T}) where T
    for x in 1:nev.p_num
        nev.thetas[x] = mebius(nev.nev_value[x])
    end
end

"""
Calculate Pick matrix elements' denominator
"""
function calc_pick_denominator!(nev::Nevanlinna_data{T}) where T
    nev.pick_denominator = zeros(Complex{T},nev.p_num,nev.p_num)
    for i in 1:nev.p_num
        @views nev.pick_denominator[1:i,i] .= - mebius.(nev.smpl_imag_points[1:i])*conj(mebius(nev.smpl_imag_points[i]))
    end
end
