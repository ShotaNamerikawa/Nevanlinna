function calc_inv_temp(temp::T) where T
    1/(k_B_in_eV*temp)
end

function change_imag_data(ndata::S,
    matsu_omega::Array{Complex{T},1},
    green_iwn::Array{Complex{T},1},
    ;
    N_imag = 0,
    N_imag_reduce = 0) where {T,S <: NevanlinnaData{T}}
    ndata.matsu_omega = matsu_omega
    ndata.green_iwn = green_iwn
    ndata.N_imag = N_imag
    ndata.hnw = size(matsu_omega,1)
    set_calculation!(ndata)
    interpolate!(ndata)
end


function set_calculation!(ndata::NevanlinnaRealData{T}) where T
    if ndata.N_imag == 0
        ndata.N_imag  =  Nevanlinna.calc_opt_N_imag(ndata.hnw, ndata.matsu_omega, ndata.green_iwn) #- ndata.N_imag_reduce
        println("N_imag $(ndata.N_imag)")
        ndata.N_imag -= ndata.N_imag_reduce
    end
    ndata.imag_data = Nevanlinna.ImagDomainData(ndata.matsu_omega, ndata.green_iwn, ndata.N_imag)
    ndata.phis = Nevanlinna.calc_phis(ndata.imag_data)
end

function set_calculation!(ndata::Union{NevanlinnaRealData{T},NevanlinnaComplexData{T}},
                          matsu_omega,
                          green_iwn) where T
    if ndata.N_imag == 0
        ndata.N_imag  =  Nevanlinna.calc_opt_N_imag(ndata.hnw, matsu_omega, green_iwn) #- ndata.N_imag_reduce
        println("N_imag $(ndata.N_imag)")
        ndata.N_imag -= ndata.N_imag_reduce
    end
    ndata.imag_data = Nevanlinna.ImagDomainData(matsu_omega, green_iwn, ndata.N_imag)
    ndata.phis = Nevanlinna.calc_phis(ndata.imag_data)
end

function interpolate!(ndata::NevanlinnaRealData{T}) where T
    ndata.real_data = Nevanlinna.RealDomainData(ndata.N_real, ndata.omega_max, ndata.eta, ndata.sum, T=T)
    ndata.abcd = Nevanlinna.calc_abcd(ndata.imag_data, ndata.real_data, ndata.phis)
    hardy_matrix = Nevanlinna.calc_hardy_matrix(ndata.real_data, ndata.H_max) 
    Nevanlinna.evaluation!(ndata.real_data, ndata.abcd, ndata.H_max, ndata.ab_coeff, hardy_matrix)
end

function delete_data_rel_N_imag!(ndata::Union{NevanlinnaRealData,NevanlinnaComplexData}; rel_pos::Int64 = 0)
    if rel_pos > 0
        println("enter positive value for rel_pos")
        throw(Exception)
    end
    N_imag = Nevanlinna.calc_opt_N_imag(size(ndata.matsu_omega,1),ndata.matsu_omega,ndata.green_iwn)
    index = collect(range(1,N_imag)) 
    deleteat!(index,N_imag - rel_pos) #delete one green function data at N_imag - rel_pos point  
    ndata.N_imag = N_imag - 1
    ndata.imag_data = Nevanlinna.ImagDomainData(ndata.matsu_omega, ndata.green_iwn,index)
    ndata.phis = Nevanlinna.calc_phis(ndata.imag_data)
end