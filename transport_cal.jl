@everywhere include("./transport.jl")
@everywhere include("./Hamiltonian.jl")
#include("./transport.jl")
@everywhere using MultiFloats
@everywhere using DoubleFloats
using .TransPort
using PhysicalConstants.CODATA2014
using LinearAlgebra
using TensorOperations

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

Si_hr = TransPort.Hamiltonian.HR.make_hr_from_file("/home/shota/Silicon/nonrel/si_hr.dat")
hk = TransPort.Hamiltonian.HR.ftransform(Si_hr,[0,0,0])
e,v = eigen(hk)
hkwa = TransPort.Hamiltonian.HR.generate_hkwa(Si_hr,[0,0,0])
va = TransPort.Hamiltonian.HR.generate_va(hkwa,v)
energy = [0]
vals = read_from_winfile("/home/shota/Silicon/nonrel/si.win")
println("vals")
println(vals)
kmesh = vals["mp_grid"]
num_wann = vals["num_wann"]
lattice_vector = vals["lattice_vector"]
println(lattice_vector)
println("start")
tp = TransPort.Transport(Si_hr,[40,40,40],300,energy,8
;mu=6.2478,delta=0.01,transform_mat = transpose(lattice_vector))
println("end")
gien = zeros(Complex,size(tp.ir.omegaF_sp,1),size(e,1))


TransPort.calc_diagonalized_green_func_k!(tp,[0,0,0],gien)

vgvg = zeros(Complex,size(tp.ir.omegaF_sp,1),3,3)
vgvg_band = zeros(Complex,size(tp.ir.omegaF_sp,1),3,3,size(e,1))


for x in 1:size(tp.ir.omegaF_sp,1)
    for y in 1:3
        for z in 1:3
            vgvg[x,y,z] = sum(va[:,y].*gien[x,:].*va[:,z].*gien[x,:])
            vgvg_band[x,y,z,:] = va[:,y].*gien[x,:].*va[:,z].*gien[x,:]
        end
    end
end

gien_wan = zeros(ComplexF64,size(tp.ir.omegaF_sp,1),size(e,1),size(e,1))

vgvg_wan = zeros(Complex,size(tp.ir.omegaF_sp,1),3,3)
vgvg_wan_band_ham = zeros(Complex,size(tp.ir.omegaF_sp,1),3,3,size(e,1),size(e,1))

TransPort.calc_green_func_k!(tp,[0,0,0],gien_wan)

for x in 1:size(tp.ir.omegaF_sp,1)
    for y in 1:3
        for z in 1:3
            tmp = hkwa[:,:,y]*gien_wan[x,:,:]*hkwa[:,:,z]*gien_wan[x,:,:]
            vgvg_wan_band_ham[x,y,z,:,:] = v'*tmp*v 
            for i in 1:size(e,1)
                vgvg_wan[x,y,z] += tmp[i,i]
            end
        end
    end 
end

println("hamiltonain")
println(vgvg[1])

println("wannier")
println(vgvg_wan[1])

println("vgvg_ham")
println(vgvg_band[1,1,1,1])

println("vgvg_wan_ham")
println(vgvg_wan_band_ham[1,1,1,1,1])

hkha = TransPort.Hamiltonian.HR.generate_hkha(hkwa,v)
hkha_diag = zeros(Complex,size(hkha,1),3)
for x in 1:3
    hkha_diag[:,x] = diag(hkha[:,:,x])
end
println("hkha diag")
println(hkha_diag[:,1])
println("va")
println(va[:,1])
println("diag green")
println(gien[1,:])

gien_diag = zeros(ComplexF64,size(tp.ir.omegaF_sp,1),size(e,1),size(e,1))
#=
println("diag matrix type")
println(size(v))

println("v type")
println(typeof(v))
println("gien type")
println(typeof(gien_wan))
=#
@tensoropt gien_diag[i,j,k] = v'[j,l]*gien_wan[i,l,m]*v[m,k] 

for x in 1:size(e,1)
    println(gien_diag[1,x,x])
end


TransPort.write_correlation(tp,"/home/shota/Silicon/nonrel/Si_correlate_green_non_diagonalized_tau_100.dat")
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
