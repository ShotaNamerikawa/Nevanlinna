"""
    ConvertedGreenData(imag_points,nev_value)

Contain Nevanlinna function's values (nev_value) at sampling points (imag_points) 

"""
mutable struct ConvertedGreenData{T<:Complex}
    imag_points::Vector{T}
    nev_value::Vector{T}
    function ConvertedGreenData(imag_points::Vector{S},nev_value::Vector{S}
        ;order::String = "dsc") where {S<:Complex}
        if order == "asc"
            println("asc")
            imag_points = reverse(imag_points)
            nev_value = reverse(nev_value)
        elseif order == "dsc"
            println("dsc")
            imag_points = imag_points
            nev_value = nev_value    
        end
        new{typeof(imag_points[1])}(imag_points,nev_value)
    end
end

"""
    readfile(greenfile,[T,order,statistics])

Read samping points and Nevanlinna values on them from a file.
"""
function readgreenfile(greenfile::String;T=Float64,order::String = "dsc",statistics::String="F")
    fp = open(greenfile,"r") 
    lines = readlines(fp)
    close(fp)
    imag_points = Complex[]
    nev_value = Complex[]
    for x in lines
        line = split(x)
        matsubara = parse(T,line[1])
        G_real = parse(T,line[2])
        G_imag = parse(T,line[3])
        if matsubara > 0
            push!(imag_points,1.0im*matsubara)
            if statistics == "F"
                push!(nev_value,-G_real-1im*G_imag)
            elseif statistics == "B"
                push!(nev_value,-1.0im*matsubara*(G_real+1im*G_imag))
                #push!(nev_value,-(1/-1.0im*matsubara)*(G_real+1im*G_imag))
            end
        end            
    end
    ConvertedGreenData(imag_points,nev_value;order = order)
end