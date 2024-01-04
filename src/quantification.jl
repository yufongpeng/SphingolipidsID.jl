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
    featuretable = data.table
    sort!(featuretable, [:id, :datafile])
    set!(data.config, :method, nothing)
    isnothing(project) && (featuretable.match_id .= featuretable.id; return sort!(featuretable, [:match_id, :match_score]))
    method = project.quantification.batch.method
    data.config[:method] = method
    i = 0
    for ft in groupview(getproperty(:id), featuretable)
        mz1 = mean(getproperty.(ft, :mz1))
        mz2 = data.mz2[first(ft).mz2_id]
        id = findfirst(i -> between(mz1, i.mz1, mz_tol) && between(mz2, i.mz2, mz_tol) && between(mean(rt_correction(class(i.analyte), Table(ft; rtx = repeat([i.rt], length(ft)))) .- ft.rt), 0, rt_tol), method.analytetable)
        featuretable.match_id[i + 1:i + length(ft)] .= isnothing(id) ? repeat([0.0], length(ft)) : repeat([method.analytetable.id[id]], length(ft))
        featuretable.match_score[i + 1:i + length(ft)] .= isnothing(id) ? repeat([0.0], length(ft)) : 
            repeat(
                [
                    ((method.analytetable.mz1[id] - mz1) / mz1) ^ 2 + 
                    ((method.analytetable.mz2[id] - mz2) / mz2) ^ 2 + 
                    (mean(rt_correction(class(method.analytetable.analyte[id]), Table(ft; rtx = repeat([method.analytetable.rt[id]], length(ft)))) .- ft.rt) / mean(ft.rt)) ^ 2
                ], length(ft)
                )
        i += length(ft)
    end
    filter!(x -> x.match_id > 0, featuretable)
    sort!(featuretable, [:match_id, :match_score])
end

function analysistable(featuretable::Table; data = :area, method = nothing, default = 0.0, maxncol = 300)
    gft = groupview(getproperty(:match_id), featuretable)
    datafile = unique(featuretable.datafile)
    et = eltype(getproperty(featuretable, data))
    if !isa(default, et)
        et = promote_type(et, typeof(default))
    end
    nt = Vector{Pair{Symbol, Vector{et}}}(undef, length(gft))
    if isnothing(method)
        analyte = Vector{Symbol}(undef, length(gft))
        for (i, (match_id, ft)) in zip(eachindex(nt), pairs(gft))
            dt = getproperty(ft, data)
            nt[i] = Symbol(match_id) => map(datafile) do file
                id = findfirst(==(file), ft.datafile)
                isnothing(id) ? default : dt[id]
            end
            analyte[i] = Symbol(match_id)
        end 
    else
        analyte = Vector{eltype(method.analytetable.analyte)}(undef, length(gft))
        for (i, (match_id, ft)) in enumerate(pairs(gft))
            dt = getproperty(ft, data)
            nt[i] = Symbol(method.analytetable.analyte[match_id]) => map(datafile) do file
                id = findfirst(==(file), ft.datafile)
                isnothing(id) ? default : dt[id]
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
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    push!(project.data, qcdata_mrm!(featuretable, project; rt_correction, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    set!(project.quantification.config, :qcdata, last(project.data))
    last(project.data)
end

qcdata_mrm(project::Project, featuretable::Table; kwargs...) = qcdata_mrm!(project, deepcopy(featuretable); kwargs...)
qcdata_mrm(featuretable::Table, project = nothing; kwargs...) = qcdata_mrm!(deepcopy(featuretable), project; kwargs...)
qcdata_mrm!(project::Project, featuretable::Table; kwargs...) = qcdata_mrm!(featuretable, project; kwargs...)
function qcdata_mrm!(featuretable::Table, project = nothing;
                rt_correction = nothing,  
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :injection_order, :match_id], [only ∘ unique, only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    n = length(unique(featuretable.datafile))
    mrm = MRM!(Table(featuretable; match_id = zeros(Int, length(featuretable)), match_score = -getproperty(featuretable, signal)); combine = false, rt_tol, mz_tol, n = 1, signal = nothing, est_fn, err_fn, err_tol, other_fn)
    transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    at = analysistable(mrm.table; method = mrm.config[:method], data = signal)
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
                    r2_threshold = 0.8,
                    nlevel = 5,
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    push!(project.data, serialdilution_mrm!(featuretable, concentration, project; pointlevel, rt_correction, r2_threshold, nlevel, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    set!(project.quantification.config, :serialdilution, last(project.data))
    last(project.data)
end

serialdilution_mrm(project::Project, featuretable::Table, concentration::Vector; kwargs...) = serialdilution_mrm!(project, deepcopy(featuretable), concentration; kwargs...)
serialdilution_mrm(featuretable::Table, concentration::Vector, project = nothing; kwargs...) = serialdilution_mrm!(deepcopy(featuretable), concentration, project; kwargs...)
serialdilution_mrm!(project::Project, featuretable::Table, concentration::Vector; kwargs...) = serialdilution_mrm!(featuretable, concentration, project; kwargs...)
function serialdilution_mrm!(featuretable::Table, concentration::Vector, project = nothing; 
                pointlevel = nothing, 
                name = r"cal.*_(\d*)-r\d*.*",
                rt_correction = nothing,
                r2_threshold = 0.8,
                nlevel = 5,
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = std, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :injection_order, :match_id], [only ∘ unique, only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    mrm = MRM!(Table(featuretable; match_id = zeros(Int, length(featuretable)), match_score = -getproperty(featuretable, signal)); combine = false, rt_tol, mz_tol, n = 1, signal = nothing, est_fn, err_fn, err_tol, other_fn)
    transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    config = dictionary(pairs((; concentration, rt_correction, r2_threshold, nlevel, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn)))
    at = analysistable(mrm.table; method = mrm.config[:method], data = signal)
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
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile, :injection_order, :match_id], [only ∘ unique, only ∘ unique, only ∘ unique])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    n = length(unique(featuretable.datafile))
    mrm = MRM!(Table(featuretable; match_id = zeros(Int, length(featuretable)), match_score = -getproperty(featuretable, signal)); combine = false, rt_tol, mz_tol, n = 1, signal = false, est_fn, err_fn, err_tol, other_fn)
    transition_matching!(mrm, project; rt_correction, rt_tol, mz_tol)
    at = analysistable(mrm.table; method = mrm.config[:method], data = signal)
    if !isnothing(project)
        signal = :estimated_concentration
        update_quantification!(project.quantification.batch, at; estimated_concentration = signal)
    end
    config = dictionary(pairs((; rt_correction, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    QuantData(mrm, at, config)
end