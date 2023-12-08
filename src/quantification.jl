quantification_mrm(project::Project, quantifier::Tuple; mz_tol = last(project.data).config[:mz_tol], rt_tol = last(project.data).config[:rt_tol], kwargs...) = quantification_mrm(project.analytes, quantifier; kwargs..., mz_tol, rt_tol, anion = project.appendix[:anion])
quantification_mrm(aquery::AbstractQuery, quantifier::Tuple; mz_tol = last(aquery.project.data).config[:mz_tol], rt_tol = last(aquery.project.data).config[:rt_tol], kwargs...) = quantification_mrm(aquery.result, quantifier; kwargs..., mz_tol, rt_tol, anion = aquery.project.appendix[:anion])
quantification_mrm(quantifier::Tuple, pq::Union{Project, AbstractQuery}; kwargs...) = quantification_mrm(pq, quantifier; kwargs...)

quantification_mrm!(project::Project, quantifier::Tuple; mz_tol = last(project.data).config[:mz_tol], rt_tol = last(project.data).config[:rt_tol], kwargs...) = 
    (project.quantification = quantification_mrm(project.analytes, quantifier; kwargs..., mz_tol, rt_tol, anion = project.appendix[:anion]))
quantification_mrm!(aquery::AbstractQuery, quantifier::Tuple; mz_tol = last(aquery.project.data).config[:mz_tol], rt_tol = last(aquery.project.data).config[:rt_tol], kwargs...) = 
    (aquery.project.quantification = quantification_mrm(aquery.result, quantifier; kwargs..., mz_tol, rt_tol, anion = aquery.project.appendix[:anion]))
quantification_mrm!(quantifier::Tuple, pq::Union{Project, AbstractQuery}; kwargs...) = quantification_mrm!(pq, quantifier; kwargs...)

function quantification_mrm(analytes::AbstractVector{<: AbstractAnalyteID}, quantifier::Tuple;
                            qualifier = Tuple[],
                            isd = nothing,
                            mz_tol = 0.35, rt_tol = 0.1,
                            anion = :acetate,
                            default_mz2 = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1))),
                            coelution_fn = x -> lastindex(x.analytes),
                            isd_fn = x -> 0)
    quanttable = mrmquanttable(analytes, quantifier; qualifier, mz_tol, rt_tol, anion, default_mz2, isd_fn)
    isdtable = isnothing(isd) ? empty(quanttable) : isa(isd, NamedTuple) ? 
        mrmquanttable(isd.analytes, 
                        :quantifier in propertynames(isd) ? getproperty(isd, :quantifier) : quantifier;
                        qualifier = :qualifier in propertynames(isd) ? getproperty(isd, :qualifier) : qualifier,
                        coelution_fn = x -> lastindex(x.analyte),
                        isd_fn = x -> 0,
                        mz_tol, rt_tol, anion, default_mz2
                        ) : 
                    mrmquanttable(isd, quantifier; qualifier, mz_tol, rt_tol, anion, default_mz2, coelution_fn = x -> lastindex(x.analyte), isd_fn = x -> 0)
    Quantification(quanttable, isdtable, dictionary(pairs((; rt_tol, mz_tol, quantifier, qualifier, default_mz2, isd_fn, coelution_fn))))
end

function transition_matching!(rawdata::Table, project = nothing; file = :repeat, rt_tol = 0.1, mz_tol = 0.35, rt_correction = nothing)
    quanttable, isdtable, config = isnothing(project) ? ([], [], []) : (project.quantification.quanttable, project.quantification.isdtable, project.quantification.config)
    if !isempty(quanttable)
        quanttable = isnothing(rt_correction) ? quanttable : Table(quanttable; rt = map(quanttable) do q
            rt_correction(class(q.analyte[config[:coelution_fn](q)]), rt(q))
        end)
        rawdata.match_id .= if isempty(isdtable)
            map(rawdata) do r
                id = findfirst(i -> between(r.mz1, i.mz1, mz_tol) && between(r.mz2, i.mz2, mz_tol) && between(r.rt, i.rt, rt_tol), quanttable)
                isnothing(id) ? 0 : quanttable.id[id]
            end
        else
            map(rawdata) do r
                id = findfirst(i -> between(r.mz1, i.mz1, mz_tol) && between(r.mz2, i.mz2, mz_tol) && between(r.rt, i.rt, rt_tol), isdtable)
                if isnothing(id)
                    id = findfirst(i -> between(r.mz1, i.mz1, mz_tol) && between(r.mz2, i.mz2, mz_tol) && between(r.rt, i.rt, rt_tol), quanttable)
                    isnothing(id) ? 0 : quanttable.id[id]
                else
                    -isdtable.id[id]
                end
            end
        end       
        if !isempty(isdtable)
            rawdata = Table(rawdata; ratio = map(eachindex(rawdata)) do i 
                id = rawdata.match_id[i]
                id_isd = >(id, 0) ? (quanttable.isd_id[id] == 0 ? (return getproperty(rawdata, project.appendix[:signal])[i]) : findfirst(x -> ==(x.match_id, -quanttable.isd_id[id]) && ==(getproperty(x, file), getproperty(rawdata, file)[i]), rawdata)) : return 0
                isnothing(id_isd) ? 0 : getproperty(rawdata, project.appendix[:signal])[i] / getproperty(rawdata, project.appendix[:signal])[id_isd]
            end)
        end
    end
    rawdata
end

function qcdata_mrm!(project::Project, featuretable::Table; 
                    rt_correction = nothing,
                    name = r"pooledqc.*_(\d*).*",
                    rt_tol = last(project.data).config[:rt_tol], 
                    mz_tol = last(project.data).config[:mz_tol], 
                    signal = project.appendix[:signal],
                    est_fn = mean, 
                    err_fn = rsd, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:repeat, :polarity, :match_id], [nothing, nothing, first])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    push!(project.data, qcdata_mrm(featuretable, project; rt_correction, name, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    last(project.data)
end

qcdata_mrm(project::Project, featuretable::Table; kwargs...) = qcdata_mrm(featiretable, project; kwargs...)
function qcdata_mrm(featuretable::Table, project = nothing;
                rt_correction = nothing,  
                name = r"pooledqc.*_(\d*).*",
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:repeat, :polarity, :match_id], [nothing, nothing, first])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    rawdata = filter(x -> occursin(name, x.datafile), featuretable)
    rep = map(x -> parse(Int, match(name, x)[1]), rawdata.datafile)
    n = length(unique(rep))
    rawdata = transition_matching!(Table(rawdata; datafile = nothing, repeat = rep, match_id = zeros(Int, length(rawdata))), project; rt_correction, rt_tol, mz_tol)
    signal = :ratio in propertynames(rawdata) ? :ratio : signal
    mrm = MRM!(rawdata; raw = true, rt_tol, mz_tol, n = 1, signal, est_fn, err_fn, err_tol, other_fn = Dictionary{Symbol, Any}())
    featuretable = filter_combine_features!(isnothing(project) ? deepcopy(rawdata) : filter(r -> r.match_id > 0, rawdata); use_match_id = !isnothing(project), raw_id = true, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
    config = dictionary(pairs((; rt_correction, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    QCData(mrm, featuretable, config)
end

function serialdilution_mrm!(project::Project, featuretable::Table, level::Vector; 
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
    default_other_fn = Dictionary{Symbol, Any}([:mz2_id, :collision_energy, :FWHM, :symmetry, :level, :polarity, :height, :area, :match_id], [nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, first])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    push!(project.data, serialdilution_mrm(featuretable, level, project; rt_correction, name, r2_threshold, nlevel, rt_tol, mz_tol, signal, est_fn, err_fn, err_tol, other_fn))
    last(project.data)
end

function serialdilution_mrm(featuretable::Table, level::Vector, project = nothing; 
                rt_correction = nothing, 
                name = r"cal.*_(\d*)-r(\d*).*",
                r2_threshold = 0.8,
                nlevel = 5,
                rt_tol = isnothing(project) ? 0.1 : last(project.data).config[:rt_tol], 
                mz_tol = isnothing(project) ? 0.35 : last(project.data).config[:mz_tol], 
                signal = isnothing(project) ? :area : project.appendix[:signal],
                est_fn = mean, 
                err_fn = std, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:mz2_id, :collision_energy, :FWHM, :symmetry, :level, :polarity, :height, :area, :ratio, :match_id], [nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, first])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    rawdata = filter(x -> occursin(name, x.datafile), featuretable)
    lv = map(x -> parse.(Int, match(name, x)), rawdata.datafile)
    rawdata = transition_matching!(Table(rawdata; datafile = nothing, level = lv, match_id = zeros(Int, length(rawdata))), project; rt_correction, rt_tol, mz_tol, file = :level)
    signal = :ratio in propertynames(rawdata) ? :ratio : signal
    mrm = MRM!(rawdata; raw = true, rt_tol, mz_tol, n = 1, signal, est_fn, err_fn, err_tol, other_fn = Dictionary{Symbol, Any}())
    featuretable = filter_combine_features!(isnothing(project) ? deepcopy(rawdata) : filter(r -> r.match_id > 0, rawdata); use_match_id = !isnothing(project), raw_id = true, rt_tol, mz_tol, n = 1, signal = false, est_fn, err_fn, err_tol, other_fn)
    config = dictionary(pairs((; rt_correction, r2_threshold, nlevel, rt_tol, mz_tol, n = 1, signal, est_fn, err_fn, err_tol, other_fn)))
    # r2, model, data_id
    r2s = zeros(Float64, length(featuretable))
    data_ids = Vector{Vector{Int}}(undef, length(featuretable))
    Threads.@threads for i in eachindex(featuretable)
        rawid = findall(x -> in(x, featuretable.raw_id[i]), mrm.table.id)
        v = getproperty(mrm.table, mrm.config[:signal])[rawid]
        l = getindex.(mrm.table.level[rawid], 1)
        ord = sortperm(l)
        l = level[l[ord]]
        v = v[ord]
        id = rawid[ord]
        _, data_ids[i], r2s[i] = rec_r2(id, l, v, nlevel, r2_threshold, 0)
    end
    SerialDilution(mrm, Table(featuretable; r2 = r2s, data_id = data_ids), level, config)
end

function rec_r2(id, level, value, nlevel, r2_threshold, r2_current)
    length(unique(level)) < nlevel && return (false, Int[], -Inf)
    r2 = cor(level, value) ^ 2
    r2 < r2_current && return (false, id, r2_current)
    r2 >= r2_threshold && return (true, id, r2)
    s = findfirst(!=(first(level)), level)
    e = findfirst(==(last(level)), level) - 1
    r2_current = r2
    id_current = id
    status = false
    for i in s:e
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && return (status, id_current, r2_current)
    for i in firstindex(id):(s - 1)
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && (status, id_current, r2_current)
    for i in (e + 1):lastindex(id)
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    (status, id_current, r2_current)
end

function set_quanttable!(project::Project, mrmquanttable::Table)
    set!(project.appendix, :quanttable, mrmquanttable)
    project
end

function append_quanttable!(project::Project, mrmquanttable::Table)
    old = get(project.appendix, :quanttable, nothing)
    isnothing(old) && return set_quanttable!(project, mrmquanttable)
    id = maximum(old.id)
    mrmquanttable = Table(mrmquanttable; id = collect(eachindex(mrmquanttable)) .+ id)
    append!(old, mrmquanttable)
    project
end

function set_isd!(project::Project, isdtable::Table)
    set!(project.appendix, :isd, isdtable)
    project
end

function append_isd!(project::Project, isdtable::Table)
    old = get(project.appendix, :isd, nothing)
    isnothing(old) && return set_isd!(project, isdtable)
    id = maximum(old.id)
    isdtable = Table(isdtable; id = collect(eachindex(isdtable)) .+ id)
    append!(old, isdtable)
    project
end