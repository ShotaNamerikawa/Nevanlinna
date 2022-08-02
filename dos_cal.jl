using Distributed
@everywhere include("./transport.jl")
@everywhere include("./Hamiltonian.jl")
@everywhere using .TransPort
#@everywhere using .Hamiltonian.HR
@everywhere using DoubleFloats

Fe_hr = TransPort.Hamiltonian.HR.make_hr_from_file("/home/shota/Nevanlinna/nev_test/Fe_hr.dat")
energy = collect(range(10,40;step=0.1))
kmesh = [8,8,8]
tp = TransPort.Transport(Fe_hr,kmesh,300,energy,8,delta=1)

TransPort.write_dos(tp,"Fe_k8_up_dos_sp.dat","Fe_k8_dn_dos_sp.dat")
