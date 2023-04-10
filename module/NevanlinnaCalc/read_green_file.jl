function read_green_data_from_file(ndata::NevanlinnaRealData{T},green_file::String,freq_col=1) where T
    fp = open(green_file)
    lines = readlines(fp)
    matsu_omega = Complex{T}[]
    green_iwn = Complex{T}[]
    for line in lines
        data = parse.(T,split(line))
        if ndata.positive_omega == true && data[freq_col] <= 0
            continue
        end
        push!(matsu_omega,1.0im*data[freq_col])
        push!(green_iwn,data[ndata.component+1] + 1.0im*data[ndata.component+2])
    end
    ndata.hnw = size(matsu_omega,1)
    matsu_omega,green_iwn
end

function read_green_data_from_file(ndata::NevanlinnaComplexData{T},green_file::String) where T
    fp = open(green_file)
    lines = readlines(fp)
    matsu_omega = Complex{T}[]
    green_iwn = Complex{T}[]
    for line in lines
        data = parse.(T,split(line))
        if ndata.positive_omega == true && data[1] <= 0
            continue
        end
        push!(matsu_omega,1.0im*data[1])
        push!(green_iwn,data[ndata.component+1] + 1.0im*data[ndata.component+2])
    end
    ndata.hnw = size(matsu_omega,1)
    matsu_omega,green_iwn
end

function impose_particle_hole_2_green_iwn!(green_iwn;statistics="fermion")
    if statistics == "fermion"
        green_iwn .= 1.0im*imag.(green_iwn)
    elseif statistics == "boson"
        green_iwn .= real.(green_iwn)
    end
end

function derive_max_Pick_criterion_num(matsu_omega,green_iwn)
    N = size(matsu_omega,1) 
    println("smpl num")
    println(N)
    N_imag = Nevanlinna.calc_opt_N_imag(N,matsu_omega,green_iwn)
end

function cut_l_component(green_iwn::Array{Complex{T}},cut_num,smpl) where T
    gl = fit(smpl,green_iwn)
    gl[end-cut_num:end] .= T(0)
    green_iwn_new = evaluate(smpl,gl)
end