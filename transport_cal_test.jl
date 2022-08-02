@everywhere include("/home/shota/Nevanlinna/transport.jl")
#include("./transport.jl")
@everywhere using MultiFloats
@everywhere using DoubleFloats
using .TransPort
using .SpCond
using PhysicalConstants.CODATA2014
using LinearAlgebra

function read_from_winfile(win_file::String)
    bohr = BohrRadius.val*10^10 # bohr radius in Angstrom
    vals = Dict()
    lattice_vector = zeros(3,3)
    kmesh = zeros(3)
    fp = open(win_file)
    while true
        if eof(fp)
            break
        end

        line = readline(fp)
        
        if isnothing(match(r"num_wann",line)) == false
            vals["num_wann"] = parse(Int64,split(line)[end])
        end

        if isnothing(match(r"begin unit_cell_cart",line)) == false
            line = readline(fp)
            if isnothing(match(r"bohr",line)) == false
                scale = bohr
                line = readline(fp)
            elseif isnothing(match(r"ang",line)) == false
                scale = 1 # for Angstrom
                line = readline(fp)
            end
            lattice_vector[1,:] = scale*parse.(Float64,split(line))
            for i in 2:3
                line = readline(fp)
                lattice_vector[i,:] = scale*parse.(Float64,split(line))
            end
            vals["lattice_vector"] = lattice_vector
        end

        if isnothing(match(r"mp_grid",line)) == false
            mp_grid = parse.(Int64,split(line)[end-2:end])  
            vals["mp_grid"] = mp_grid
        end
    end
    close(fp)
    vals
end

Si_hr = TransPort.Hamiltonian.HR.make_hr_from_file("/home/shota/Silicon/si_hr.dat")
energy = collect(range(-10,10;step=0.001))
vals = read_from_winfile("/home/shota/Silicon/si.win")
println("vals")
println(vals)
kmesh = vals["mp_grid"]
num_wann = vals["num_wann"]
lattice_vector = vals["lattice_vector"]
println("lattice vector is ")
println(lattice_vector)

ts = TransportfromSC.read_file("Si_tau_100_kmesh_8_spcond.dat",300.0)
tc = TransportfromSC.cal_all_thermo_electric(ts,6.2449)
println(tc[1])
#sc = SpCond.Spcond(Si_hr,[8,8,8],300,energy,8
#;mu=6.2478,delta=0.01,transform_mat = transpose(lattice_vector))

#=
k =[0,0,0]

hk = TransPort.Hamiltonian.HR.ftransform(sc.hr,k)

ek,U = eigen(hk)

hkwa = TransPort.Hamiltonian.HR.generate_hkwa(sc.hr,k)

println("hkwa")
println(hkwa[1,:,2])

va = TransPort.Hamiltonian.HR.generate_va(hkwa,U)

cal = SpCond.cal_spectralcond_k(sc,[0,0,0])[1101,:,:]

print("gamma")
println(sc.constscat)

println("omega")
println(sc.energy[1101])

println("va")
println(va[1,:])

println("vel")
println((sc.transform_mat*va[1,:])./TransPort.hbar_in_eV)

println("sp_k")
println(cal)
=#
#SpCond.write_spectralcond(sc,"Si_tau_100_kmesh_8_spcond.dat")


#=
v = TransPort.Hamiltonian.HR.generate_hkwa(tp.hr,[0,0,0])
println("v value")
for x in 1:16
    println(v[x,x,1])
end
println("not diagonalized")
println(TransPort.cal_correlation_k(tp,[0,0,0],v)[1:5,1])

println("diagonalized ")
println(TransPort.cal_diagonalized_correlation_k(tp,[0,0,0])[1:5,1])
=#

#=
println("iomega is ")
println(tp.ir.omegaF_sp[1])

hk = TransPort.Hamiltonian.HR.ftransform(tp.hr,[0,0,0])

eigval1,eigvec1 = eigen(hk)

mat = ((1.0im*(tp.ir.omegaF_sp[1] + tp.constscat*sign(tp.ir.omegaF_sp[1])) + tp.mu)*tp.E
       -hk)

eigval2,eigvec2 = eigen(mat)

println("eigval hk ")
println(eigval1[:])
println("eigval plus")
println(eigval1[:])
println("eigval mu")
println(eigval2)
println("eigval green")
println(1 ./ eigval2[:])

gien = zeros(Complex,size(tp.ir.omegaF_sp,1),tp.hr.num_wann,tp.hr.num_wann)
TransPort.calc_green_func_k!(tp,[0,0,0],gien)
println("hk -> green")
println(eigen(gien[1,:,:]).values[1:16])

gien = zeros(Complex,size(tp.ir.omegaF_sp,1),tp.hr.num_wann)
TransPort.calc_diagonalized_green_func_k!(tp,[0,0,0],gien)
println("ek -> green")
println(gien[1,1:16])
=#

#TransPort.write_correlation(tp,"/home/shota/Nevanlinna/Si_correlate_green_tau_100.dat")
#=
pg = TransPort.calc_green_func_k(tp,[0.0,0.0,0.0])
mg = TransPort.calc_green_func_k(tp,[0.0,0.0,0.0],0.0)
println(size(pg))
println(size(mg))
println(size(pg[1]))
println(size(mg[1]))
println(pg[2][:,3])
println(mg[2][3,:])
=#
#TransPort.write_correlation(tp,"current_current_green_2.dat")
