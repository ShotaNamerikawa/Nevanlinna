using Distributed
addprocs(4)
using Plots
include("/home/shota/Nevanlinna/Hamiltonian.jl")
using .Hamiltonian

begin
using LinearAlgebra
include("/home/shota/Nevanlinna/Hamiltonian.jl")
using .Hamiltonian
using SharedArrays
end

begin
kmesh = [1,1,1]
hr = Hamiltonian.HR.make_hr_from_file("/home/shota/Ag/pwscf_hr.dat")
hkm = Hamiltonian.HamKmesh.Hamkmesh(kmesh,hr)
end

function calc_dos(hkm,energy;eta=0.01)
    dos = SharedArray{Float64}(size(energy,1))
    for i in 1:size(energy,1)
        for kp in hkm.kp_list
        e = eigvals(Hamiltonian.HR.ftransform(hkm.hr,kp))
        green = 1 ./ ((energy[i]+1.0im*eta) .- e)
        dos[i] += -1/pi*imag(sum(green))/hkm.nk
        end
    end     
    dos
end
GC.gc()
energy = collect(range(8,18;length=100))
println(energy)
@time dos = calc_dos(hkm,energy)
print(dos)
plot(energy,dos)

hkm.kp_list