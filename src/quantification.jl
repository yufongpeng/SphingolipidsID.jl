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
        conctable = ColumnDataTable(isdtable.analyte, :Level, Table(; Level = [1], (Symbol.(isdtable.analyte) .=> vectorize.(convert.(Float64, isd.concentration)))...))
    else
        isdtable = analytetable_mrm(isd, quantifier; qualifier, mz_tol, rt_tol, anion, default_mz2, coelution_fn = x -> lastindex(x.analytes), isd_fn = x -> -1)
        conctable = ColumnDataTable(isdtable.analyte, :Level, Table(; Level = [1], (Symbol.(isdtable.analyte) .=> vectorize.(ones(Float64, length(isdtable.analyte))))...))
    end
    analytetable.id .+= length(isdtable)
    Quantification(Batch(MethodTable{Table}(vcat(isdtable, analytetable), signal, Int[], conctable, nothing)), dictionary(pairs((; signal, anion, rt_tol, mz_tol, quantifier, qualifier, default_mz2, isd_fn, coelution_fn))))
end

function transition_matching!(data::MRM, method = nothing; rt_tol = 0.1, mz_tol = 0.35, rt_correction = nothing)
    featuretable = data.table
    sort!(featuretable, [:id, :injection_order])
    isnothing(method) && return sort!(featuretable, [:match_id, :id, :injection_order, :match_score])
    i = 0
    for ft in groupview(getproperty(:id), featuretable)
        mz1 = mean(ft.mz1)
        mz2 = data.mz2[first(ft).mz2_id]
        id = findfirst(i -> between(mz1, i.mz1, mz_tol) && between(mz2, i.mz2, mz_tol) && between(mean(rt_correction(class(i.analyte), Table(ft; rtx = repeat([i.rt], length(ft)))) .- ft.rt), 0, rt_tol), method.analytetable)
        featuretable.match_id[i + 1:i + length(ft)] .= isnothing(id) ? repeat([0.0], length(ft)) : repeat([method.analytetable.id[id]], length(ft))
        featuretable.match_score[i + 1:i + length(ft)] .= isnothing(id) ? repeat([0.0], length(ft)) : 
        begin
            lmz1 = method.analytetable.mz1[id]
            lmz2 = method.analytetable.mz2[id]
            lrt = rt_correction(class(method.analytetable.analyte[id]), Table(ft; rtx = repeat([method.analytetable.rt[id]], length(ft))))
            @. (ft.mz1 / lmz1 - 1) ^ 2 + (mz2 / lmz2 - 1) ^ 2 + (ft.rt / lrt - 1) ^ 2
        end
        i += length(ft)
    end
    filter!(x -> x.match_id > 0, featuretable)
    sort!(featuretable, [:match_id, :id, :injection_order, :match_score])
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
            gdf = groupview(getproperty(:datafile), ft)
            nt[i] = Symbol(match_id) => map(datafile) do file
                num = get(gdf, file, nothing)
                isnothing(num) ? default : first(getproperty(num, data))
            end
            analyte[i] = Symbol(match_id)
        end 
    else
        analyte = Vector{eltype(method.analytetable.analyte)}(undef, length(gft))
        for (i, (match_id, ft)) in enumerate(pairs(gft))
            if length(unique(ft.id)) > 1
                dts = map(unique(ft.id)) do idx
                    subft = ft[ft.id .== idx]
                    gdf = groupview(getproperty(:datafile), subft)
                    map(datafile) do file
                        num = get(gdf, file, nothing)
                        isnothing(num) ? default : first(getproperty(num, data))
                    end
                end
                _, idx = findmax(dts) do dt
                    mean(dt)
                end
                dt = dts[idx]
            else
                gdf = groupview(getproperty(:datafile), ft)
                dt = map(datafile) do file
                    num = get(gdf, file, nothing)
                    isnothing(num) ? default : first(getproperty(num, data))
                end
            end
            nt[i] = Symbol(method.analytetable.analyte[match_id]) => dt
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

function match(mrm::MRM, project = nothing; kwargs...)
    # check columns
    hasproperty(mrm.table, :match_id) || throw(ArgumentError("Table in MRM dose not have column :match_id; please use `qtMRM` or `qtMRM!`."))
    hasproperty(mrm.table, :match_score) || throw(ArgumentError("Table in MRM dose not have column :match_score; please use `qtMRM` or `qtMRM!`."))
    @info "MRM | Transition matching"
    _match!(deepcopy(mrm), project, isnothing(project) ? nothing : project.quantification; kwargs...)
end

function match!(mrm::MRM, project = nothing; kwargs...)
    # check columns
    hasproperty(mrm.table, :match_id) || throw(ArgumentError("Table in MRM dose not have column :match_id; please use `qtMRM` or `qtMRM!`."))
    hasproperty(mrm.table, :match_score) || throw(ArgumentError("Table in MRM dose not have column :match_score; please use `qtMRM` or `qtMRM!`."))
    @info "MRM | Transition matching"
    _match!(mrm, project, isnothing(project) ? nothing : project.quantification; kwargs...)
end

function _match!(mrm::MRM, project::Project, qt::Quantification; 
                rt_correction = nothing, 
                rt_tol = qt.config[:rt_tol],
                mz_tol = qt.config[:mz_tol])
    transition_matching!(mrm, qt.batch.method; rt_correction, rt_tol, mz_tol)
    set!(mrm.config, :match_config, dictionary(pairs((; project, rt_correction, rt_tol, mz_tol))))
    mrm
end

function _match!(mrm::MRM, project::Union{Project, Nothing}, qt::Nothing; 
                rt_correction = nothing,  
                rt_tol = mrm.config[:rt_tol], 
                mz_tol = mrm.config[:mz_tol], 
    )
    transition_matching!(mrm, nothing; rt_correction, rt_tol, mz_tol)
    set!(mrm.config, :match_config, dictionary(pairs((; project, rt_correction, rt_tol, mz_tol))))
    mrm
end

function quantify(mrm::MRM, project = get(mrm.config, :match_config, Dictionary([:project], [nothing]))[:project]; kwargs...)
    hasproperty(mrm.table, :match_id) || throw(ArgumentError("Table in MRM dose not have column :match_id; please use `qtMRM` or `qtMRM!`."))
    hasproperty(mrm.table, :match_score) || throw(ArgumentError("Table in MRM dose not have column :match_score; please use `qtMRM` or `qtMRM!`."))
    @info "MRM | Quantifying analytes"
    _quantify(mrm, project, isnothing(project) ? nothing : project.quantification; kwargs...)
end

function _quantify(mrm::MRM, project::Project, qt::Quantification; 
                    signal = project.appendix[:signal],
                    rel_sig = :relative_signal,
                    est_conc = :estimated_concentration)
    at = analysistable(mrm.table; method = qt.batch.method, data = signal)
    set_quantification!(at, project.quantification.batch; rel_sig, est_conc)
end

function _quantify(mrm::MRM, project::Union{Project, Nothing}, qt::Nothing; 
                    signal = isnothing(project) ? :area : project.appendix[:signal],
                    rel_sig = :relative_signal,
                    est_conc = :estimated_concentration)
    analysistable(mrm.table; method = qt, data = signal)
end

function quantify!(mrm::MRM, project = get!(mrm.config, :match_config, Dictionary([:project], [nothing]))[:project]; kwargs...)
    hasproperty(mrm.table, :match_id) || throw(ArgumentError("Table in MRM dose not have column :match_id; please use `qtMRM` or `qtMRM!`."))
    hasproperty(mrm.table, :match_score) || throw(ArgumentError("Table in MRM dose not have column :match_score; please use `qtMRM` or `qtMRM!`."))
    @info "MRM | Quantifying analytes"
    _quantify!(mrm, project, isnothing(project) ? nothing : project.quantification; kwargs...)
end

function _quantify!(mrm::MRM, project::Project, qt::Quantification; 
                    signal = project.appendix[:signal],
                    rel_sig = :relative_signal,
                    est_conc = :estimated_concentration)
    id = findfirst(x -> ===(x, mrm), project.data)
    id = isnothing(id) ? (push!(project.data, mrm); lastindex(project.data)) : id
    at = analysistable(mrm.table; method = qt.batch.method, data = signal)
    update_quantification!(project.quantification.batch, at; rel_sig, est_conc)
    set!(project.quantification.config, :data_id, id)
    set!(project.quantification.config, :signal, signal)
    set!(project.quantification.config, :rel_sig, rel_sig)
    set!(project.quantification.config, :est_conc, est_conc)
    at
end

function _quantify!(mrm::MRM, project::Project, qt::Nothing; 
                    signal = project.appendix[:signal],
                    rel_sig = :relative_signal,
                    est_conc = :estimated_concentration)
    id = findfirst(x -> ===(x, mrm), project.data)
    id = isnothing(id) ? (push!(project.data, mrm); lastindex(project.data)) : id
    at = analysistable(mrm.table; method = qt, data = signal)
    signaltable = getproperty(at, signal)
    @info "MRM | Constructing new Quantification"
    conctable = if isa(signaltable, ColumnDataTable)
        tbl = Table((; signaltable.samplecol => [1], (Symbol.(signaltable.analyte) .=> [[1.0] for i in eachindex(signaltable.analyte)])...))
        ColumnDataTable(signaltable.analyte, signaltable.samplecol, tbl)
    else
        tbl = Table((; signaltable.analytecol => signaltable.analyte, Symbol(1) => [[1.0] for i in eachindex(signaltable.analyte)]))
        RowDataTable(signaltable.analyte, signaltable.analytecol, [Symbol(1)], tbl)
    end
    batch = Batch(MethodTable{Table}(
            Table(; analyte = signaltable.analyte, isd = repeat([0], length(signaltable.analyte)), calibration = collect(eachindex(signaltable.analyte))), 
            signal, 
            Int[], 
            conctable, 
            nothing), 
        at
    )
    project.quantification = Quantification(batch, Dictionary{Symbol, Any}(
                                [:data_id, :signal, :anion, :rt_tol, :mz_tol, :rel_sig, :est_conc], 
                                [id, project.appendix[:signal], project.appendix[:anion], mrm.config[:rt_ol], mrm.config[:mz_tol], rel_sig, est_conc])
                            )
    update_quantification!(project.quantification.batch, at; rel_sig, est_conc)
    at
end

function _quantify!(mrm::MRM, project::Nothing, qt::Nothing; 
            signal = :area,
            rel_sig = :relative_signal,
            est_conc = :estimated_concentration)
    analysistable(mrm.table; method = qt, data = signal)
end

set_qc_id!(project::Project, crit::T) where {T <: Function} = set_qc_id!(project.quantification, crit)
function set_qc_id!(qt::Quantification, crit::T) where {T <: Function}
    s = getproperty(qt.batch.data, qt.config[:signal]).sample
    id = findall(crit, s)
    set!(qt.config, :qc_id, length(id) == length(s) ? nothing : id)
end

qctable(project::Project; signal = :estimated_concentration, est_fn = mean, err_fn = default_error, err_tol = 0.25, id = project.quantification.config[:qc_id]) = qctable(project.quantification; signal, est_fn, err_fn, err_tol, id)
qctable(qt::Quantification; signal = :estimated_concentration, est_fn = mean, err_fn = default_error, err_tol = 0.25, id = qt.config[:qc_id]) = qctable(qt.batch; signal, est_fn, err_fn, err_tol, id)
qctable(batch::Batch; signal = :estimated_concentration, est_fn = mean, err_fn = default_error, err_tol = 0.25, id = nothing) = qctable(batch.data; signal, est_fn, err_fn, err_tol, id)
qctable(at::AnalysisTable; signal = :estimated_concentration, est_fn = mean, err_fn = default_error, err_tol = 0.25, id = nothing) = qctable(getproperty(at, signal); est_fn, err_fn, err_tol, id)
function qctable(dt::ChemistryQuantitativeAnalysis.AbstractDataTable; est_fn = mean, err_fn = default_error, err_tol = 0.25, id = nothing)
    if isnothing(id)
        rsds = map(err_fn, eachanalyte(dt))
        Table(; analyte = dt.analyte, mean = map(est_fn, eachanalyte(dt)), rsd = rsds, pass = rsds .<= err_tol)
    else
        rsds = map(x -> err_fn(x[id]), eachanalyte(dt))
        Table(; analyte = dt.analyte, mean = map(x -> est_fn(x[id]), eachanalyte(dt)), rsd = rsds, pass = rsds .<= err_tol)
    end
end

set_serialdilution_id!(project::Project, crit::T) where {T <: Function} = set_serialdilution_id!(project.quantification, crit)
function set_serialdilution_id!(qt::Quantification, crit::T) where {T <: Function}
    s = getproperty(qt.batch.data, qt.config[:signal]).sample
    id = findall(crit, s)
    set!(qt.config, :serialdilution_id, length(id) == length(s) ? nothing : id)
end

apply_or_return_vector(fn, v) = map(fn, v)
apply_or_return_vector(fn::Vector, v) = fn

run_serialdilution!(project::Project, pointlevel, concentration; signal = project.appendix[:signal], id = project.quantification.config[:serialdilution_id], r2_threshold = 0.8) = run_serialdilution!(project.quantification, pointlevel, concentration; signal, id, r2_threshold)
function run_serialdilution!(qt::Quantification, pointlevel, concentration; signal = qt.config[:signal], id = qt.config[:serialdilution_id], r2_threshold = 0.8)
    batch = run_serialdilution(qt.batch, pointlevel, concentration; signal, id, r2_threshold)
    set!(qt.config, :serialdilution, batch)
    batch
end

run_serialdilution(project::Project, pointlevel, concentration; signal = project.appendix[:signal], id = project.quantification.config[:serialdilution_id], r2_threshold = 0.8) = run_serialdilution(project.quantification, pointlevel, concentration; signal, id, r2_threshold)
run_serialdilution(qt::Quantification, pointlevel, concentration; signal = qt.config[:signal], id = qt.config[:serialdilution_id], r2_threshold = 0.8) = run_serialdilution(qt.batch, pointlevel, concentration; signal, id, r2_threshold)
run_serialdilution(batch::Batch, pointlevel, concentration; signal = :area, id = nothing, r2_threshold = 0.8) = run_serialdilution(batch.data, pointlevel, concentration; signal, id, r2_threshold)
run_serialdilution(at::AnalysisTable, pointlevel, concentration; signal = :area, id = nothing, r2_threshold = 0.8) = run_serialdilution(getproperty(at, signal), pointlevel, concentration; signal, id, r2_threshold)
function run_serialdilution(signaltable::ChemistryQuantitativeAnalysis.AbstractDataTable, pointlevel, concentration; signal = :area, id = nothing, r2_threshold = 0.8)
    if !isnothing(id)
        if signaltable isa ColumnDataTable
            signaltable = ColumnDataTable(signaltable.analyte, signaltable.samplecol, signaltable.table[id])
        else
            p = (signaltable.analytecol, signaltable.samplename[id]...)
            signaltable = RowDataTable(signaltable.analyte, signaltable.analytecol, signaltable.samplename[id], getproperties(signaltable.table, p))
        end
    end
    pointlevel = convert(Vector{Int}, apply_or_return_vector(pointlevel, signaltable.sample))
    concentration = convert(Vector{Float64}, apply_or_return_vector(concentration, signaltable.sample))
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
    nlevel = length(conctable)
    Threads.@threads for cal in batch.calibration
        _, id, _ = rec_r2(collect(eachindex(cal.table.x)), cal.table.x, cal.table.y, nlevel, r2_threshold, 0)
        cal.table.include .= false
        cal.table.include[id] .= true
        update_calibration!(cal, method)
    end
    batch
end

serialdilutiontable(project::Project; r2_threshold = 0.8) = serialdilutiontable(project.quantification; r2_threshold)
serialdilutiontable(qt::Quantification; r2_threshold = 0.8) = serialdilutiontable(qt.config[:serialdilution]; r2_threshold)
function serialdilutiontable(batch::Batch; r2_threshold = 0.8)
    r2s = map(batch.calibration) do c
        r2(c.model)
    end
    Table(; analyte = batch.analyte, r² = r2s, pass = r2s .>= r2_threshold, var"signal range" = signal_range.(batch.calibration))
end