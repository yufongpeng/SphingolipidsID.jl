quantification_mrm(project::Project, quantifier::Tuple; mz_tol = last(project.data).config[:mz_tol], rt_tol = last(project.data).config[:rt_tol], kwargs...) = quantification_mrm(project.analyte, quantifier; kwargs..., mz_tol, rt_tol, signal = project.appendix[:signal], anion = project.appendix[:anion])
quantification_mrm(aquery::AbstractQuery, quantifier::Tuple; mz_tol = last(aquery.project.data).config[:mz_tol], rt_tol = last(aquery.project.data).config[:rt_tol], kwargs...) = quantification_mrm(aquery.result, quantifier; kwargs..., mz_tol, rt_tol, signal = aquery.project.appendix[:signal], anion = aquery.project.appendix[:anion])
quantification_mrm(quantifier::Tuple, pq::Union{Project, AbstractQuery}; kwargs...) = quantification_mrm(pq, quantifier; kwargs...)

set_quantification_mrm!(project::Project, quantifier::Tuple; mz_tol = last(project.data).config[:mz_tol], rt_tol = last(project.data).config[:rt_tol], kwargs...) = 
    (project.quantification = quantification_mrm(project.analyte, quantifier; kwargs..., mz_tol, rt_tol, signal = project.appendix[:signal], anion = project.appendix[:anion]))
set_quantification_mrm!(aquery::AbstractQuery, quantifier::Tuple; mz_tol = last(aquery.project.data).config[:mz_tol], rt_tol = last(aquery.project.data).config[:rt_tol], kwargs...) = 
    (aquery.project.quantification = quantification_mrm(aquery.result, quantifier; kwargs..., mz_tol, rt_tol, signal = aquery.project.appendix[:signal], anion = aquery.project.appendix[:anion]))
set_quantification_mrm!(quantifier::Tuple, pq::Union{Project, AbstractQuery}; kwargs...) = set_quantification_mrm!(pq, quantifier; kwargs...)

function quantification_mrm(analyte::AbstractVector{<: AbstractAnalyteID}, quantifier::Tuple;
                            qualifier = Tuple[],
                            isd = nothing,
                            mz_tol = 0.35, rt_tol = 0.1,
                            signal = :area,
                            anion = :acetate,
                            default_mz2 = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1))),
                            coelution_fn = x -> lastindex(x.analytes),
                            isd_fn = x -> 0)
    analytetable = analytetable_mrm(analyte, quantifier; qualifier, mz_tol, rt_tol, anion, default_mz2, isd_fn)
    isnothing(isd) && throw(ArgumentError("Keyword argument isd can not be nothing; external calibration has not been implemented."))
    if isa(isd, NamedTuple)
        isdtable = analytetable_mrm(isd.analyte, 
                        :quantifier in propertynames(isd) ? getproperty(isd, :quantifier) : quantifier;
                        qualifier = :qualifier in propertynames(isd) ? getproperty(isd, :qualifier) : qualifier,
                        coelution_fn = x -> lastindex(x.analytes),
                        isd_fn = x -> -1,
                        mz_tol, rt_tol, anion, default_mz2
                        )
        conctable = ColumnDataTable(isdtable.analyte, :Level, Table(; Level = [1], (Symbol.(isdtable.analyte) .=> vectorize.(isd.concentration))...))
    else
        isdtable = analytetable_mrm(isd, quantifier; qualifier, mz_tol, rt_tol, anion, default_mz2, coelution_fn = x -> lastindex(x.analytes), isd_fn = x -> -1)
        conctable = ColumnDataTable(isdtable.analyte, :Level, Table(; Level = [1], (Symbol.(isdtable.analyte) .=> vectorize.(ones(Float64, length(isdtable.analyte))))...))
    end
    analytetable.id .+= length(isdtable)
    Quantification(Batch(MethodTable{Table}(vcat(isdtable, analytetable), signal, Int[], conctable, nothing)), dictionary(pairs((; signal, anion, rt_tol, mz_tol, quantifier, qualifier, default_mz2, isd_fn, coelution_fn))))
end

function transition_matching!(data::MRM, project = nothing; rt_tol = 0.1, mz_tol = 0.35, rt_correction = nothing)
    rawdata = data.table
    isnothing(project) && (rawdata.match_id .= rawdata.id; return rawdata)
    method =  project.quantification.batch.method
    if !isnothing(rt_correction) 
        method = MethodTable{Table}(Table(method.analytetable; rt = map(method.analytetable) do q
        rt_correction(class(q.analyte), rt(q))
    end), method.signal, method.pointlevel, method.conctable, method.signaltable)
    end
    rawdata.match_id .= 
        mapmany(groupview(getproperty(:id), rawdata)) do raw
            mz1 = mean(getproperty.(raw, :mz1))
            mz2 = data.mz2[first(raw).mz2_id]
            rt = mean(getproperty.(raw, :rt))
            id = findfirst(i -> between(mz1, i.mz1, mz_tol) && between(mz2, i.mz2, mz_tol) && between(rt, i.rt, rt_tol), method.analytetable)
            isnothing(id) && return repeat([0.0], length(raw))
            repeat([method.analytetable.id[id]], length(raw))
        end
    filter!(x -> x.match_id > 0, rawdata)
end

function analysistable(rawdata::Table, project = nothing; data = isnothing(project) ? :area : project.appendix[:signal], default = 0.0, maxncol = 300)
    method = isnothing(project) ? nothing : project.quantification.batch.method
    grawdata = groupview(getproperty(:match_id), rawdata)
    datafile = unique(rawdata.datafile)
    et = eltype(getproperty(rawdata, data))
    if !isa(default, et)
        et = promote_type(et, typeof(default))
    end
    nt = Vector{Pair{Symbol, Vector{et}}}(undef, length(grawdata))
    if isnothing(method)
        analyte = Vector{Symbol}(undef, length(grawdata))
        for (i, (match_id, raw)) in zip(eachindex(nt), pairs(grawdata))
            dt = getproperty(raw, data)
            nt[i] = Symbol(match_id) => map(datafile) do file
                id = findfirst(==(file), raw.datafile)
                isnothing(id) && return default
                dt[id]
            end
            analyte[i] = Symbol(match_id)
        end 
    else
        analyte = Vector{eltype(method.analytetable.analyte)}(undef, length(grawdata))
        for (i, (match_id, raw)) in enumerate(pairs(grawdata))
            dt = getproperty(raw, data)
            nt[i] = Symbol(method.analytetable.analyte[match_id]) => map(datafile) do file
                id = findfirst(==(file), raw.datafile)
                isnothing(id) && return default
                dt[id]
            end
            analyte[i] = method.analytetable.analyte[match_id]
        end
    end
    if length(nt) > maxncol && length(last(first(nt))) < length(nt)
        datafile = Symbol.(datafile)
        AnalysisTable([data], [
            RowDataTable(analyte, :analyte, datafile, Table(; analyte, (datafile .=> map(eachindex(datafile)) do i
                    map(x -> last(x)[i], nt)
                end)...))])
    else
        AnalysisTable([data], [ColumnDataTable(analyte, :Sample, Table(; Sample = datafile, nt...))])
    end
end

function set_qcdata_mrm!(project::Project, featuretable::Table; 
                    rt_correction = nothing,
                    name = r"pooledqc.*_\d*.*",
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    push!(project.data, qcdata_mrm!(featuretable, project; rt_correction, name, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    set!(project.quantification.config, :qcdata, last(project.data))
    last(project.data)
end

qcdata_mrm(project::Project, featuretable::Table; kwargs...) = qcdata_mrm!(project, deepcopy(featuretable); kwargs...)
qcdata_mrm(featuretable::Table, project = nothing; kwargs...) = qcdata_mrm!(deepcopy(featuretable), project; kwargs...)
qcdata_mrm!(project::Project, featuretable::Table; kwargs...) = qcdata_mrm!(featuretable, project; kwargs...)
function qcdata_mrm!(featuretable::Table, project = nothing;
                rt_correction = nothing,  
                name = r"pooledqc.*_\d*.*",
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :match_id], [only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    rawdata = filter!(x -> occursin(name, x.datafile), featuretable)
    n = length(unique(rawdata.datafile))
    mrm = MRM!(Table(rawdata; match_id = zeros(Int, length(rawdata))); combine = false, rt_tol, mz_tol, n = 1, signal = nothing, est_fn, err_fn, err_tol, other_fn)
    rawdata = transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    at = analysistable(mrm.table, project)
    if !isnothing(project)
        signal = :estimated_concentration
        set_quantification!(at, project.quantification.batch; estimated_concentration = signal)
    end
    config = dictionary(pairs((; rt_correction, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    QCData(mrm, at, config)
end

function qualitytable(data::QCData, signal::Symbol = data.config[:signal])
    rsds = map(data.config[:err_fn], eachanalyte(getproperty(data.table, signal)))
    Table(; analyte = data.analyte, mean = map(mean, eachanalyte(getproperty(data.table, signal))), rsd = rsds, pass = rsds .<= data.config[:err_tol])
end

function set_serialdilution_mrm!(project::Project, featuretable::Table, concentration::Vector; 
                    pointlevel = nothing, 
                    rt_correction = nothing,  
                    name = r"cal.*_(\d*)-r(\d*).*",
                    r2_threshold = 0.8,
                    nlevel = 5,
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    push!(project.data, serialdilution_mrm!(featuretable, concentration, project; pointlevel, rt_correction, name, r2_threshold, nlevel, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    set!(project.quantification.config, :serialdilution, last(project.data))
    last(project.data)
end

serialdilution_mrm(project::Project, featuretable::Table, concentration::Vector; kwargs...) = serialdilution_mrm!(project, deepcopy(featuretable), concentration; kwargs...)
serialdilution_mrm(featuretable::Table, concentration::Vector, project = nothing; kwargs...) = serialdilution_mrm!(deepcopy(featuretable), concentration, project; kwargs...)
serialdilution_mrm!(project::Project, featuretable::Table, concentration::Vector; kwargs...) = serialdilution_mrm!(featuretable, concentration, project; kwargs...)
function serialdilution_mrm!(featuretable::Table, concentration::Vector, project = nothing; 
                pointlevel = nothing, 
                rt_correction = nothing, 
                name = r"cal.*_(\d*)-r\d*.*",
                r2_threshold = 0.8,
                nlevel = 5,
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = std, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :match_id], [only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    rawdata = filter!(x -> occursin(name, x.datafile), featuretable)
    mrm = MRM!(Table(rawdata; match_id = zeros(Int, length(rawdata))); combine = false, rt_tol, mz_tol, n = 1, signal = nothing, est_fn, err_fn, err_tol, other_fn)
    rawdata = transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    config = dictionary(pairs((; concentration, rt_correction, r2_threshold, nlevel, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn)))
    at = analysistable(rawdata, project)
    signaltable = getproperty(at, signal)
    pointlevel = isnothing(pointlevel) ? map(x -> parse(Int, first(match(name, string(x)))), signaltable.sample) : pointlevel
    if signaltable isa ColumnDataTable
        conctable = ColumnDataTable(signaltable.analyte, :Level, Table(; Level = collect(eachindex(concentration)), (signaltable.analytename .=> repeat([concentration], length(signaltable.analyte)))...))
    else
        level = Symbol.(collect(eachindex(concentration)))
        conctable = RowDataTable(signaltable.analyte, :analyte, level, Table(; analyte = signaltable.analyte, (level .=> map(eachindex(level)) do i 
            repeat([concentration[i]], length(signaltable.analyte))
        end)...))
    end
    method = MethodTable{Table}(conctable, signaltable, signal, pointlevel)
    batch = Batch(method)
    Threads.@threads for cal in batch.calibration
        _, id, _ = rec_r2(collect(eachindex(cal.table.x)), cal.table.x, cal.table.y, nlevel, r2_threshold, 0)
        cal.table.include .= false
        cal.table.include[id] .= true
        update_calibration!(cal, method)
    end
    SerialDilution(mrm, batch, config)
end

function qualitytable(data::SerialDilution)
    r2s = map(data.batch.calibration) do c
        r2(c.model)
    end
    Table(; analyte = data.analyte, r² = r2s, pass = r2s .>= data.config[:r2_threshold], var"signal range" = signal_range.(data.batch.calibration))
end

function rec_r2(id, conc, value, nlevel, r2_threshold, r2_current)
    length(unique(conc)) < nlevel && return (false, Int[], -Inf)
    r2 = cor(conc, value) ^ 2
    r2 < r2_current && return (false, id, r2_current)
    r2 >= r2_threshold && return (true, id, r2)
    s = findfirst(!=(first(conc)), conc)
    e = findfirst(==(last(conc)), conc) - 1
    r2_current = r2
    id_current = id
    status = false
    for i in s:e
        idx = setdiff(eachindex(conc), i)
        status, id_new, r2_new = rec_r2(id[idx], conc[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && return (status, id_current, r2_current)
    for i in firstindex(id):(s - 1)
        idx = setdiff(eachindex(conc), i)
        status, id_new, r2_new = rec_r2(id[idx], conc[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && (status, id_current, r2_current)
    for i in (e + 1):lastindex(id)
        idx = setdiff(eachindex(conc), i)
        status, id_new, r2_new = rec_r2(id[idx], conc[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    (status, id_current, r2_current)
end

function set_quantdata_mrm!(project::Project, featuretable::Table; 
                    rt_correction = nothing,
                    name = r"S\d*.*",
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    push!(project.data, quantdata_mrm!(featuretable, project; rt_correction, name, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    set!(project.quantification.config, :quantdata, last(project.data))
    last(project.data)
end

quantdata_mrm(project::Project, featuretable::Table; kwargs...) = quantdata_mrm!(project, deepcopy(featuretable); kwargs...)
quantdata_mrm(featuretable::Table, project = nothing; kwargs...) = quantdata_mrm!(deepcopy(featuretable), project; kwargs...)
quantdata_mrm!(project::Project, featuretable::Table; kwargs...) = quantdata_mrm!(featuretable, project; kwargs...)
function quantdata_mrm!(featuretable::Table, project = nothing;
                rt_correction = nothing,  
                name = r"S\d*.*",
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :match_id], [only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    rawdata = filter!(x -> occursin(name, x.datafile), featuretable)
    n = length(unique(rawdata.datafile))
    mrm = MRM!(Table(rawdata; match_id = zeros(Int, length(rawdata))); combine = false, rt_tol, mz_tol, n = 1, signal = false, est_fn, err_fn, err_tol, other_fn)
    rawdata = transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    at = analysistable(mrm.table, project)
    if !isnothing(project)
        signal = :estimated_concentration
        update_quantification!(project.quantification.batch, at; estimated_concentration = signal)
    end
    config = dictionary(pairs((; rt_correction, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    QuantData(mrm, at, config)
end