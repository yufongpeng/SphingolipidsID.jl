"""
    q!(qc, project::Project)
    q!(project::Project, qc)
    q!(qc, aquery::Query)
    q!(aquery::Query, qc)
    q!(qc, reusable::ReusableQuery)
    q!(reusable::ReusableQuery, qc)

Create a query retrieving analytes from a project, previous query or reusable query.

# Query commands
1. Unary function returning bolean value
2. `:topscore => topN` or `:topsc => topN`
3. `attribute => criteria`\n
    Attribute: A unary function applying to `analyte`, e.g. `rt`, `mz`, `mw`, `ncb`, etc.\n
    Criteria:\n
    1. `ri\"interval\"`: a real interval, e.g. `ri\"[5, 5.5)\"`, `ri\"(-inf, 7]\"`, etc.
    2. A unary function returning bolean value
    3. The above objects wrapped in `FunctionalFunction`, e.g. `allow_unknown(ri\"[5, 10]" ∘ abs)`, `only_known(<(10))`, etc.
    3. `(lowerbound, upperbound)` equivalent to `ri\"(lowerbound, upperbound)\"`
    4. `[lowerbound, upperbound]` equivalent to `ri\"[lowerbound, upperbound]\"`
    5. `lowerbound` equivalent to `ri\"[lowerbound, lowerbound + 1)\"`
4. Specific structure
    * Subtype or instance of `ClassSP`
    * Instance of `LCB`, `ACYL`, `ChainSP`, or `SPID`
5. `Symbol`\n
    Query based on analyte status.\n
    Format: `attribute + logic`, e.g., `:class!`.\n
    Attribute:\n
    * `class`
    * `chain`
    * `rt`
    * `error`
    * `isf`
    * `total`
    Logic:\n
    * nothing: positive(`1`)
    * `_`: not determined(`0`) or positive(`1`)
    * `!`: negative(`-1`)
6. `MRM`
7. `QueryCommands`\n
    Any object of `QueryCommands` created by
    * `qnot`
    * `qor`
    * `qand`
    Any valid queries can be wrapped into `QueryCommands`.
# Keyword Arguments
* `view`: whether copy the `result`.
"""
function q!(project::Project; view = true)
    v = @view project[:]
    v = view ? v : Vector(v)
    Query(project, v, QueryAnd(QueryCommands[]), view)
end

q!(qc, project::Project; kwargs...) = q!(project::Project, qc; kwargs...)

function q!(project::Project, qc; view = true)
    v = @view project[:]
    v = view ? v : Vector(v)
    q!(Query(project, v, QueryAnd(QueryCommands[]), view), qc)
end

q!(qc, aquery::Query; kwargs...) = q!(aquery, qc; kwargs...)
reuse(aquery::Query) = ReusableQuery(aquery)
q!(qc, reusable::ReusableQuery; kwargs...) = q!(reusable, qc; kwargs...)
q!(reusable::ReusableQuery, qc; kwargs...) = q!(reuse_copy(reusable.query), qc; kwargs...)

q!(aquery::Query, qf::T; view = aquery.view, neg = false, kwargs...) where {T <: Function} = finish_query!(aquery, qcmd(qf, neg), qf; view)

function q!(aquery::Query, qc::Pair; view = aquery.view, neg = false, kwargs...)
    @match qc.first begin
        :topscore   => return finish_query!(aquery, qcmd(qc, neg), query_score(aquery, qc.second; kwargs...); view)
        :topsc      => return finish_query!(aquery, qcmd(qc, neg), query_score(aquery, qc.second; kwargs...); view)
        _           => nothing
    end
    fn2 = @match qc.second begin
        (lb, ub)            => in(RealInterval(lb, ub, <, <))
        [lb, ub]            => in(RealInterval(lb, ub, <=, <=))
        v::Number           => in(RealInterval(v, v + 1, <=, <))
        f::RealIntervals    => in(f)
        _                   => qc.second
    end
    qf = @match qc.first begin
        &mz => mz(fn2)
        fn1 => fn2 ∘ fn1
    end
    finish_query!(aquery, qcmd(qc.first => fn2, neg), qf; view)
end

query_fn(qc::Type{<: ClassSP}) = analyte -> isa(class(analyte), qc) ||
    any(isa(isomer, qc) for isomer in isomer_tuple(class(analyte)))
query_fn(qc::ClassSP) = analyte -> ≡(class(analyte), qc) ||
    any(≡(ia, iq) for (ia, iq) in Iterators.product(isomer_tuple(class(analyte)), isomer_tuple(qc)))

query_fn(qc::LCB) = analyte -> iscompatible(lcb(analyte), qc)
query_fn(qc::ACYL) = analyte -> iscompatible(acyl(analyte), qc)
query_fn(qc::ChainSP) = analyte -> iscompatible(chain(analyte), qc)
query_fn(qc::SPID) = analyte -> query_fn(qc.class)(analyte) && query_fn(qc.chain)(analyte)

q!(aquery::Query, qc::Type{<: ClassSP}; view = aquery.view, neg = false) = isempty(methods(qc, ())) ?
    finish_query!(aquery, qcmd(qc, neg), query_fn(qc); view) : q!(aquery, qc(); view, neg)
q!(aquery::Query, qc; view = aquery.view, neg = false) =
    finish_query!(aquery, qcmd(qc, neg), query_fn(qc); view)

function q!(aquery::Query, qc::Symbol; view = aquery.view, neg = false)
    qf = @match qc begin
        :class      => analyte -> ≡(analyte.state[state_id(:class)], 1)
        :chain      => analyte -> ≡(analyte.state[state_id(:chain)], 1)
        :rt         => analyte -> ≡(analyte.state[state_id(:rt)], 1)
        :error      => analyte -> ≡(analyte.state[state_id(:error)], 1)
        :isf        => analyte -> ≡(analyte.state[state_id(:isf)], 1)
        :total      => analyte -> ≡(analyte.state[state_id(:total)], 1)
        :manual     => analyte -> ≡(analyte.state[state_id(:manual)], 1)
        :class_     => analyte -> ≡(analyte.state[state_id(:class)], 0)
        :chain_     => analyte -> ≡(analyte.state[state_id(:chain)], 0)
        :rt_        => analyte -> ≡(analyte.state[state_id(:rt)], 0)
        :error_     => analyte -> ≡(analyte.state[state_id(:error)], 0)
        :isf_       => analyte -> ≡(analyte.state[state_id(:isf)], 0)
        :total_     => analyte -> ≡(analyte.state[state_id(:total)], 0)
        :manual_    => analyte -> ≡(analyte.state[state_id(:manual)], 0)
        :class!     => analyte -> ≡(analyte.state[state_id(:class)], -1)
        :chain!     => analyte -> ≡(analyte.state[state_id(:chain)], -1)
        :rt!        => analyte -> ≡(analyte.state[state_id(:rt)], -1)
        :error!     => analyte -> ≡(analyte.state[state_id(:error)], -1)
        :isf!       => analyte -> ≡(analyte.state[state_id(:isf)], -1)
        :total!     => analyte -> ≡(analyte.state[state_id(:total)], -1)
        :manual!    => analyte -> ≡(analyte.state[state_id(:manual)], -1)
        #:all    => analyte -> all(==(1), analyte.state)
    end
    finish_query!(aquery, qcmd(qc, neg), qf; view)
end

function q!(aquery::Query, mrm::MRM;
                view = aquery.view, neg = false, compatible = false,
                mode::Symbol = :default)
    products = @p mrm.mz2 map(id_product(_, mrm.polarity; mz_tol = mrm.config[:mz_tol]))
    if compatible
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.appedix[:anion])
            any(!isempty(products[i]) && any(iscompatible(product, cpd) for product in products[i]) &&
                any(any(between(ms, mz1, mrm.config[:mz_tol]) for ms in ms1) for mz1 in filterview(x -> ≡(x.mz2_id, i), mrm.table).mz1) for i in eachindex(products))
        end
    else
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.appendix[:anion])
            any(!isempty(products[i]) && any(iscomponent(product, cpd) for product in products[i]) &&
                any(any(between(ms, mz1, mrm.config[:mz_tol]) for ms in ms1) for mz1 in filterview(x -> ≡(x.mz2_id, i), mrm.table).mz1) for i in eachindex(products))
        end
    end
    finish_query!(aquery, qcmd(:validated_by => mrm, neg), qf; view, objects = last.(aquery))
end

function coelution(project::Project; analyte = project.analyte, polarity::Bool = true, mz_tol = 0.35, rt_tol = 0.1)
    result = NamedTuple{(:mz1, :rt, :analytes), Tuple{Vector{Float64}, Vector{Float64}, SubArray{AnalyteSP}}}[]
    for ana in analyte
        frags = last(ana).fragment
        id = findfirst(i -> isa(frags.mz1[i].adduct, Pos) ≡ polarity, reverse(eachindex(frags)))
        isnothing(id) && continue
        mz = query_data(project, frags.source[id], frags.id[id], :mz1)
        rt = query_data(project, frags.source[id], frags.id[id], :rt)
        pushed = false
        for group in result
            abs(mean(group.mz1) - mz) > mz_tol && continue
            abs(mean(group.rt - rt)) > rt_tol && continue
            push!(group.mz1, mz)
            push!(group.rt, rt)
            push!(group.analytes, ana)
            pushed = true
            break
        end
        pushed || push!(result, (mz1 = [mz], rt = [rt], analytes = [ana]))
    end
    result
end
"""
    normalized_sig_diff(Δ)

Create a function that takes `global_scores, local_scores, prev_score, current_score` as arguments and return if normalized significant difference is within Δ.
"""
normalized_sig_diff(Δ) =
    (global_scores, local_scores, prev_score, current_score) -> _normalized_sig_diff(global_scores, local_scores, prev_score, current_score, Δ)

_normalized_sig_diff(global_scores, local_scores, prev_score, current_score, Δ) = (prev_score - current_score)/current_score < Δ
"""
    abs_sig_diff(Δ)

Create a function that takes `global_scores, local_scores, prev_score, current_score` as arguments and return if absolute significant difference is within Δ.
"""
abs_sig_diff(Δ) =
    (global_scores, local_scores, prev_score, current_score) -> _abs_sig_diff(global_scores, local_scores, prev_score, current_score, Δ)

_abs_sig_diff(global_scores, local_scores, prev_score, current_score, Δ) = (prev_score - current_score) < Δ

function query_score(aquery::Query, topN = 0.5; is_wild_card = normalized_sig_diff(0.1))
    analyte = aquery.result
    # no isomer => has isomer
    id1 = Int[]
    id2 = Int[]
    for (i, ana) in enumerate(analyte)
        (hasisomer(class(ana)) || hasisomer(chain(ana))) ? push!(id2, i) : push!(id1, i)
    end
    dict = Dict{SPID, Vector{Tuple{Float64, Int}}}()
    foreach(id -> push!(get!(dict, SPID(last(analyte[id])), Tuple{Float64, Int}[]), (analyte[id].score.second, id)), id1)
    for id in id2
        sc = analyte[id].score.second
        pushed = false
        for (ky, vl) in dict
            iscompatible(analyte[id], ky) &&
                (pushed = true; push!(vl, (sc, id)))
        end
        pushed || push!(dict, SPID(last(analyte[id])) => [(sc, id)])
    end
    final = Int[]
    len = topN >= 1 ? (vl -> topN) : (vl -> max(1, round(Int, length(vl) * topN * (1 + eps(Float64)))))
    global_scores = first.(vcat(values(dict)...))
    for vl in values(dict)
        sort!(vl, by = first, lt = isless_nan_min)
        local_scores = copy(vl)
        prev_score = 0.0
        for _ in 1:len(vl)
            isempty(vl) && break
            prev_score, id = pop!(vl)
            push!(final, id)
        end
        while !isempty(vl)
            current_score, id = pop!(vl)
            is_wild_card(global_scores, local_scores, prev_score, current_score) || break
            push!(final, id)
        end
    end
    final
end

finish_query!(aquery::Query, qc::QueryCmd, qf::T; objects = aquery, view = aquery.view) where {T <: Function} =
    finish_query!(aquery, qc, findall(qf, objects); view)
finish_query!(aquery::Query, qc::QueryNot, qf::T; objects = aquery, view = aquery.view) where {T <: Function} =
    finish_query!(aquery, qc, findall(!qf, objects); view)
function finish_query!(aquery::Query, qc::QueryCommands, id; view = aquery.view)
    push!(aquery.query.qcmd, qc)
    aquery.view = view
    aquery.result = view ? Base.view(aquery, id) : aquery[id]
    aquery
end

flatten_query(qc::QueryCmd) = qc
flatten_query(qc::QueryAnd) = QueryAnd(map(flatten_query, qc.qcmd))
flatten_query(qc::QueryOr) = QueryOr(map(flatten_query, qc.qcmd))
flatten_query(qc::QueryNot{QueryCmd}) = qc
flatten_query(qc::QueryNot{QueryOr}) = QueryAnd(map(flatten_query ∘ qnot, qc.qcmd))
flatten_query(qc::QueryNot{QueryAnd}) = QueryOr(map(flatten_query ∘ qnot, qc.qcmd))

q!(aquery::Query, qc::QueryCommands; view = true) = _q!(aquery, flatten_query(qc); view)
_q!(aquery::Query, qc::QueryNot{QueryCmd}; view = true) = q!(aquery, qc.qcmd.query; view, neg = true)
_q!(aquery::Query, qc::QueryCmd; view = true) = q!(aquery, qc.query; view)
_q!(aquery::Query, qcs::QueryAnd; view = true) = length(qcs.qcmd) == 1 ? _q!(aquery, qcs.qcmd[1]; view) : _q!(_q!(aquery, popfirst!(qcs.qcmd); view), QueryAnd(qcs.qcmd); view)
function _q!(aquery::Query, qcs::QueryOr; view = true)
    pid = parentindices(aquery.result)[1]
    qs = map(qcs.qcmd) do qc
        _q!(Query(aquery.project, view ? Base.view(aquery.project, deepcopy(pid)) : copy_wo_project.(aquery.result), QueryAnd(QueryCommands[]), view), qc; view)
    end
    qs = union!(qs...)
    push!(aquery.query.qcmd, qs.query)
    aquery.view = view
    aquery.result = view ? qs.result : Vector(qs.result)
    aquery
end

"""
    query_data(project::Project, source::Int, id::Int)
    query_data(project::Project, source::Int, id::Int, col)

Make query to data in a `project`. The data is specified by `source`, `id`, and optionally `col`.
"""
query_data(project::Project, source::Int, id::Int) = project.data[source].table[findfirst(==(id), project.data[source].table.id)]
query_data(project::Project, source::Int, id::Int, col) = getproperty(project.data[source].table, col)[findfirst(==(id), project.data[source].table.id)]
"""
    set_data!(project::Project, source::Int, id::Int, col, val)

Replace query value with `val`.
"""
set_data!(project::Project, source::Int, id::Int, col, val) = (getproperty(project.data[source].table, col)[findfirst(==(id), project.data[source].table.id)] = val)
"""
    set_data_if!(project::Project, source::Int, id::Int, col, val, crit)

Replace query value with `val` if `crit(old_value) == true`.
"""
function set_data_if!(project::Project, source::Int, id::Int, col, val, crit)
    i = findfirst(==(id), project.data[source].table.id)
    crit(getproperty(project.data[source].table, col)[i]) || return
    getproperty(project.data[source].table, col)[i] = val
end
mz(fn) = analyte -> any(fn(query_data(last(analyte).project, frag.source, frag.id).mz1) for frag in last(analyte).fragment)