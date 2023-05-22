"""
Transform Nevanlinna function into contractive function.
"""
function calc_thetas!(nev::NevanlinnaData{T}) where T
    nev.thetas = Array{T}(undef,size(nev.matsu_green,1))
    for i in eachindex(nev.matsu_green)
        nev.thetas[i] = mebius(- nev.matsu_green[i])       
    end
end
