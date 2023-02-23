function q!(project::Project; view = true)
    v = @view project[:]
    v = view ? v : Vector(v)
    Query(project, v, QueryAnd(QueryCommands[]), view)
end

q!(qq, project::Project; kwargs...) = q!(project::Project, qq; kwargs...)

function q!(project::Project, qq; view = true)
    v = @view project[:]
    v = view ? v : Vector(v)
    q!(Query(project, v, QueryAnd(QueryCommands[]), view), qq)
end

q!(qq, aquery::Query; kwargs...) = q!(aquery, qq; kwargs...)
reuse(aquery::Query) = ReUseable(aquery)
q!(qq, reuseable::ReUseable; kwargs...) = q!(reuseable, qq; kwargs...)
q!(reuseable::ReUseable, qq; kwargs...) = q!(reuse_copy(reuseable.query), qq; kwargs...)

function q!(aquery::Query, qq::Pair; view = aquery.view, neg = false, kwargs...)
    @match qq.first begin
        :topscore   => return finish_query!(aquery, qcmd(qq, neg), query_score(aquery, qq.second; kwargs...); view)
        :topsc      => return finish_query!(aquery, qcmd(qq, neg), query_score(aquery, qq.second; kwargs...); view)
        _           => nothing
    end
    qf = @match qq.first begin
        :rt => analyte -> between(analyte.rt, qq.second)
        :mw => analyte -> between(mw(last(analyte)), qq.second)
        :mz => analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1,
                                                qq.second) for frag in last(analyte).fragments)
    end
    finish_query!(aquery, qcmd(qq, neg), qf; view)
end

query_fn(qq::Type{<: ClassSP}) = analyte -> isa(class(analyte), qq) ||
    any(isa(isomer, qq) for isomer in isomer_tuple(class(analyte)))
query_fn(qq::ClassSP) = analyte -> ≡(class(analyte), qq) ||
    any(≡(ia, iq) for (ia, iq) in Iterators.product(isomer_tuple(class(analyte)), isomer_tuple(qq)))

query_fn(qq::LCB) = analyte -> iscompatible(lcb(analyte), qq)
query_fn(qq::ACYL) = analyte -> iscompatible(acyl(analyte), qq)
query_fn(qq::ChainSP) = analyte -> iscompatible(chain(analyte), qq)
query_fn(qq::SPID) = analyte -> query_fn(qq.class)(analyte) && query_fn(qq.chain)(analyte)

q!(aquery::Query, qq::Type{<: ClassSP}; view = aquery.view, neg = false) = isempty(methods(qq, ())) ?
    finish_query!(aquery, qcmd(qq, neg), query_fn(qq); view) : q!(aquery, qq(); view, neg)
q!(aquery::Query, qq; view = aquery.view, neg = false) =
    finish_query!(aquery, qcmd(qq, neg), query_fn(qq); view)

function q!(aquery::Query, qq::Symbol; view = aquery.view, neg = false)
    qf = @match qq begin
        :class  => analyte -> ≡(analyte.states[states_id(:class)], 1)
        :chain  => analyte -> ≡(analyte.states[states_id(:chain)], 1)
        :rt     => analyte -> ≡(analyte.states[states_id(:rt)], 1)
        :error  => analyte -> ≡(analyte.states[states_id(:error)], 1)
        :isf    => analyte -> ≡(analyte.states[states_id(:isf)], 1)
        :class_ => analyte -> ≡(analyte.states[states_id(:class)], 0)
        :chain_ => analyte -> ≡(analyte.states[states_id(:chain)], 0)
        :rt_    => analyte -> ≡(analyte.states[states_id(:rt)], 0)
        :error_ => analyte -> ≡(analyte.states[states_id(:error)], 0)
        :isf_   => analyte -> ≡(analyte.states[states_id(:isf)], 0)
        :class! => analyte -> ≡(analyte.states[states_id(:class)], -1)
        :chain! => analyte -> ≡(analyte.states[states_id(:chain)], -1)
        :rt!    => analyte -> ≡(analyte.states[states_id(:rt)], -1)
        :error! => analyte -> ≡(analyte.states[states_id(:error)], -1)
        :isf!   => analyte -> ≡(analyte.states[states_id(:isf)], -1)
        :all    => analyte -> all(==(1), analyte.states)
    end
    finish_query!(aquery, qcmd(qq, neg), qf; view)
end

function q!(aquery::Query, mrm::MRM;
                view = aquery.view, neg = false, compatible = false, db_product = SPDB[mrm.polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                mode::Symbol = :default)
    products = @p mrm.mz2 map(id_product(_, mrm.polarity; db = db_product, mz_tol = mrm.mz_tol))
    if compatible
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscompatible(cpd, product) for product in products[i]) &&
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ≡(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    else
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscomponent(cpd, product) for product in products[i]) &&
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ≡(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    end
    finish_query!(aquery, qcmd(:validated_by => "mrm", neg), qf; view, objects = last.(aquery))
end

function coelution(project::Project; analytes = project.analytes, polarity::Bool = true, mz_tol = 0.35, rt_tol = 0.1)
    result = NamedTuple{(:mz1, :rt, :analytes), Tuple{Vector{Float64}, Vector{Float64}, SubArray{AnalyteSP}}}[]
    for analyte in analytes
        frags = last(analyte).fragments
        id = findfirst(i -> isa(frags.mz1[i].adduct, Pos) ≡ polarity, reverse(eachindex(frags)))
        isnothing(id) && continue
        mz = query_raw(project, frags.source[id], frags.id[id], :mz1)
        rt = query_raw(project, frags.source[id], frags.id[id], :rt)
        pushed = false
        for group in result
            abs(mean(group.mz1) - mz) > mz_tol && continue
            abs(mean(group.rt - rt)) > rt_tol && continue
            push!(group.mz1, mz)
            push!(group.rt, rt)
            push!(group.analytes, analyte)
            pushed = true
            break
        end
        pushed || push!(result, (mz1 = [mz], rt = [rt], analytes = [analyte]))
    end
    result
end

normalized_sig_diff(Δ) =
    (global_scores, local_scores, prev_score, current_score) -> _normalized_sig_diff(global_scores, local_scores, prev_score, current_score, Δ)

_normalized_sig_diff(global_scores, local_scores, prev_score, current_score, Δ) = (prev_score - current_score) < Δ

abs_sig_diff(Δ) =
    (global_scores, local_scores, prev_score, current_score) -> _abs_sig_diff(global_scores, local_scores, prev_score, current_score, Δ)

_abs_sig_diff(global_scores, local_scores, prev_score, current_score, Δ) = (prev_score - current_score) < Δ

function query_score(aquery::Query, topN = 0.5; is_wild_card = normalized_sig_diff(0.1), score_fn = x -> mean(filter(!isnan, x)))
    analytes = aquery.result
    # no isomer => has isomer
    id1 = Int[]
    id2 = Int[]
    for (i, analyte) in enumerate(analytes)
        (hasisomer(class(analyte)) || hasisomer(chain(analyte))) ? push!(id2, i) : push!(id1, i)
    end
    dict = Dict{SPID, Vector{Tuple{Float64, Int}}}()
    foreach(id -> push!(get!(dict, SPID(last(analytes[id])), Tuple{Float64, Int}[]), (score_fn(analytes[id].scores), id)), id1)
    for id in id2
        sc = score_fn(analytes[id].scores)
        pushed = false
        for (ky, vl) in dict
            iscompatible(analytes[id], ky) &&
                (pushed = true; push!(vl, (sc, id)))
        end
        pushed || push!(dict, SPID(last(analytes[id])) => [(sc, id)])
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

finish_query!(aquery::Query, qq::QueryCmd, qf::T; objects = aquery, view = aquery.view) where {T <: Function} =
    finish_query!(aquery, qq, findall(qf, objects); view)
finish_query!(aquery::Query, qq::QueryNot, qf::T; objects = aquery, view = aquery.view) where {T <: Function} =
    finish_query!(aquery, qq, findall(!qf, objects); view)
function finish_query!(aquery::Query, qq::QueryCommands, id; view = aquery.view)
    push!(aquery.query.qcmd, qq)
    aquery.view = view
    aquery.result = view ? Base.view(aquery, id) : aquery[id]
    aquery
end

flatten_query(qq::QueryCmd) = qq
flatten_query(qq::QueryAnd) = QueryAnd(map(flatten_query, qq.qcmd))
flatten_query(qq::QueryOr) = QueryOr(map(flatten_query, qq.qcmd))
flatten_query(qq::QueryNot{QueryCmd}) = qq
flatten_query(qq::QueryNot{QueryOr}) = QueryAnd(map(flatten_query ∘ qnot, qq.qcmd))
flatten_query(qq::QueryNot{QueryAnd}) = QueryOr(map(flatten_query ∘ qnot, qq.qcmd))

q!(aquery::Query, qq::QueryCommands; view = true) = _q!(aquery, flatten_query(qq); view)
_q!(aquery::Query, qq::QueryNot{QueryCmd}; view = true) = q!(aquery, qq.qcmd.query; view, neg = true)
_q!(aquery::Query, qq::QueryCmd; view = true) = q!(aquery, qq.query; view)
_q!(aquery::Query, qqs::QueryAnd; view = true) = _q!(_q!(aquery, popfirst!(qqs.qcmd); view), qqs.qcmd; view)
function _q!(aquery::Query, qqs::QueryOr; view = true)
    pid = parentindices(aquery.result)[1]
    qs = map(qqs.qcmd) do qq
        _q!(Query(aquery.project, view ? Base.view(aquery.project, deepcopy(pid)) : copy_wo_project.(aquery.result), QueryAnd(QueryCommands[]), view), qq; view)
    end
    qs = union!(qs...)
    push!(aquery.query.qcmd, qs.query)
    aquery.view = view
    aquery.result = view ? qs.result : Vector(qs.result)
    aquery
end

query_raw(project::Project, source::Int, id::Int) = project.data[source].raw[findfirst(==(id), project.data[source].raw.id)]
query_raw(project::Project, source::Int, id::Int, col) = getproperty(project.data[source].raw, col)[findfirst(==(id), project.data[source].raw.id)]
set_raw!(project::Project, source::Int, id::Int, col, val) = (getproperty(project.data[source].raw, col)[findfirst(==(id), project.data[source].raw.id)] = val)
function set_raw_if!(project::Project, source::Int, id::Int, col, val, crit)
    i = findfirst(==(id), project.data[source].raw.id)
    crit(getproperty(project.data[source].raw, col)[i]) || return
    getproperty(project.data[source].raw, col)[i] = val
end

# filter/comparator