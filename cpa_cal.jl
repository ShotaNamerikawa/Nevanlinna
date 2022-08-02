@everywhere include("/home/shota/Nevanlinna/cpa.jl")
include("/home/shota/Nevanlinna/Wannier.jl")
@everywhere using .CPA
using .Hamiltonian.HR
using LinearAlgebra
using DoubleFloats
using Printf

function get_diag(A::Matrix,i::Int64)
    A[i,i]
end

function write_trG(cpa,upfile::String,downfile::String)
    trGup = zeros(Complex{Float64},size(cpa.temp.omegaF_sp)[1])
    trGdn = zeros(Complex{Float64},size(cpa.temp.omegaF_sp)[1])

    for i in 1:cpa.num_wann
        if i%2 != 0
            trGup += get_diag.(cpa.green,i)
        elseif i%2 == 0
            trGdn += get_diag.(cpa.green,i)
        end
    end

    fpup = open(upfile,"w")
    fpdn = open(downfile,"w")

    for i in 1:size(cpa.temp.omegaF_sp)[1]
        write(fpup,@sprintf("%20.16e %20.16e %20.16e\n",cpa.temp.omegaF_sp[i]
                           ,real(trGup[i]),imag(trGup[i])))
        write(fpdn,@sprintf("%20.16e %20.16e %20.16e\n",cpa.temp.omegaF_sp[i]
                           ,real(trGdn[i]),imag(trGdn[i])))
    end

    close(fpup)
    close(fpdn)

end

hr_Fe = HR.make_hr_from_file("./nev_test/Fe_hr.dat")
hr_Co = HR.make_hr_from_file("./nev_test/Co_hr.dat")

hr_list = [hr_Fe hr_Co]
weightlist = [1.0 ,0.0]
onsite_pot_list = [29.84,24.95]
T=300
kmesh = [8,8,8]
test = CPA.Cpa(T,kmesh,8.0,hr_list,weightlist,onsite_pot_list)
CPA.loop(test)

write_trG(test,"jl_trG_up_Fe_10_Co_0.dat","jl_trG_dn_Fe_10_Co_0.dat")
