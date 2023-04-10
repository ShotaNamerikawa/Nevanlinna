include("./Hamiltonian.jl")
using LinearAlgebra


hr = Hamiltonian.HR.make_hr_from_file("/home/shota/Silicon/si_hr.dat")
k = [0.5,0.5,0.5]
hk = Hamiltonian.HR.ftransform(hr,k)
hkwa = Hamiltonian.HR.generate_hkwa(hr,k)

ek,vecs = eigen(hk)

hkha = zeros(Complex,hr.num_wann,hr.num_wann,3)

for l in 1:3
    hkha[:,:,l] = vecs'* hkwa[:,:,l] * vecs
end

println(hkha[2,1,1])
