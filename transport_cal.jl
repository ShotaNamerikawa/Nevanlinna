@everywhere include("/home/shota/Nevanlinna/transport.jl")
#include("./transport.jl")
@everywhere using MultiFloats
@everywhere using DoubleFloats
using .TransPort
using PhysicalConstants.CODATA2014

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
energy = [0]
vals = read_from_winfile("/home/shota/Silicon/si.win")
println("vals")
println(vals)
kmesh = vals["mp_grid"]
num_wann = vals["num_wann"]
lattice_vector = vals["lattice_vector"]
println(lattice_vector)

tp = TransPort.Transport(Si_hr,[1,1,1],300,energy,8
;mu=6.2478,delta=0.01,transform_mat = transpose(lattice_vector))
TransPort.write_correlation(tp,"/home/shota/Nevanlinna/Si_correlate_green_tau_100.dat")
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
