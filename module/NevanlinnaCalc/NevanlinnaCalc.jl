module NevanlinnaCalc
using LinearAlgebra
using MultiFloats
using PhysicalConstants.CODATA2014
using Plots
hbar_in_eV = PlanckConstantOver2pi.val/ElementaryCharge.val
k_B_in_eV = BoltzmannConstant.val/ElementaryCharge.val

export NevanlinnaRealData,NevanlinnaComplexData,init_calculation!,change_real_freq
export delete_data_rel_N_imag
export plot_spec_fermion,plot_spec_boson,plot_spec_boson_conductivity,plot_spec_boson_resistance,plot_spec,cal_opt_conductivity,cal_opt_resistance
export cal_static_conductivity,cal_static_resistance,get_sample_and_interpolate

include("/home/shota/Nevanlinna_namerikawa/src/Nevanlinna.jl")
using .Nevanlinna

include("Nevanlinna_data.jl")
include("read_green_file.jl")
include("calculation_from_data.jl")
include("plot_data.jl")

end
