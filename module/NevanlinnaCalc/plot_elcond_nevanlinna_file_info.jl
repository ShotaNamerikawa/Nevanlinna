function get_file(file_expression::Regex,dir::String)
    matchfiles = filter(file -> false == isnothing(match(file_expression,file)),readdir(dir))
    return matchfiles
end

function get_file(cdata::CalcData)
    get_file(cdata.file_expression,cdata.dir)
end

function get_proper_regex(cdata::CalcData)
    if cdata.xdata == "mu"
        regex = r"[0-9]*\.[0-9]*"
    elseif cdata.xdata == "concentration"
        if isnothing(cdata.element) == true
            throw(ErrorException("Enter element for concentration plot"))
        else
            regex = Regex(@sprintf("%s_[0-9\\.]*",cdata.element))
        end
    end
    regex
end

function get_xdata_from_file!(cdata::CalcData,regex,file)    
    m = match(regex,file)
    #println("m.match")
    #println(m.match)
    m = match(r"[0-9\.][0-9\.]*",m.match)
    #println("m.match")
    #println(m.match)
    if isnothing(m) == true
        throw(ErrorException("file name does not contain mu or concentration value"))
    else
        push!(cdata.xdata_list,parse(Float64,m.match) - cdata.fermi_energy)
    end
end

function get_nevanlinna_data_from_file!(cdata::CalcData,file::String)
    file = joinpath([cdata.dir,file])
    nevan = NevanlinnaData(cdata.temperature,
                           file
                           ;
                           cdata.elcond_kwargs...)
    push!(cdata.ndatas,nevan)
end

function sort_data_as_to_xdata!(cdata::CalcData)
    index = sortperm(cdata.xdata_list)
    cdata.xdata_list = cdata.xdata_list[index]
    cdata.ndatas = cdata.ndatas[index]
end

function get_information_from_files!(cdata::CalcData)
    cdata.xdata_list = Float64[]
    cdata.ndatas = NevanlinnaData[]
    regex = get_proper_regex(cdata)
    for file in cdata.files
        get_xdata_from_file!(cdata,regex,file)
        get_nevanlinna_data_from_file!(cdata,file)
    end
    sort_data_as_to_xdata!(cdata)
end