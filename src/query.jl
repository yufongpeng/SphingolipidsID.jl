function query(project::Project; view = true) 
    v = @view project[:]
    v = view ? v : Vector(v)
    Query(project, v, [], view)
end

query(qq, project::Project; kwargs...) = query(project::Project, qq; kwargs...)

function query(project::Project, qq; view = true) 
    v = @view project[:]
    v = view ? v : Vector(v)
    query(Query(project, v, [], view), qq)
end

query(qq, aquery::Query; kwargs...) = query(aquery, qq; kwargs...)
reuse(aquery::Query) = ReUseable(aquery)
query(qq, reuseable::ReUseable; kwargs...) = query(reuseable, qq; kwargs...)
query(reuseable::ReUseable, qq; kwargs...) = query(reuse_copy(reuseable.query), qq; kwargs...)

function query(aquery::Query, qq::Pair{Symbol, <: Any}; view = aquery.view, fn = identity, kwargs...)
    @match qq.first begin
        :topscore   => return finish_query!(aquery, qq.first, qq.second, fn, query_score(aquery, qq.second; kwargs...), view)
        :topsc      => return finish_query!(aquery, qq.first, qq.second, fn, query_score(aquery, qq.second; kwargs...), view)
        _           => nothing
    end
    qf = @match qq.first begin
        :rt => analyte -> between(analyte.rt, qq.second)
        :mw => analyte -> between(mw(last(analyte)), qq.second)
        :mz => analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1, 
                                                qq.second) for frag in last(analyte).fragments)
    end
    finish_query!(aquery, qq.first, qq.second, fn, findall(fn ∘ qf, aquery), view)
end

function query(aquery::Query, qq::Type{<: ClassSP}; view = aquery.view, fn = identity)
    qf(analyte) = fn(
        isa(class(analyte), qq) || 
        hasisomer(class(analyte)) && 
        any(isa(isomer, qq) for isomer in class(analyte).isomer)
    )
    finish_query!(aquery, :class, qq, fn, findall(qf, aquery), view)
end

function query(aquery::Query, qq::Lcb; view = aquery.view, fn = identity)
    match = convert_type(LCB, qq)
    qf(analyte) = fn(!isnothing(_chain(analyte)) && isa(_lcb(analyte), match))
    finish_query!(aquery, :lcb, qq, fn, findall(qf, aquery), view)
end

function query(aquery::Query, qq::NACYL; view = aquery.view, fn = identity)
    match = convert_type(ACYL, qq)
    qf(analyte) = fn(
        !isnothing(_chain(analyte)) && 
        isa(_acyl(analyte), match) && 
        all((sumcomp(analyte) .- sumcomp(_lcb(analyte))) .== (qq.cb, qq.db, qq.ox))
    )
    finish_query!(aquery, :acyl, qq, fn, findall(qf, aquery), view)
end

function query(aquery::Query, qq::Type{<: LCB}; view = aquery.view, fn = identity)
    qf(analyte) = fn(!isnothing(_chain(analyte)) && isa(_lcb(analyte), qq))
    finish_query!(aquery, :lcb, qq, fn, findall(qf, aquery), view)
end

function query(aquery::Query, qq::CompoundID; view = aquery.view, fn = identity)
    match = convert_internal(qq)
    qf, sym = @match match begin
        (::Type{Nothing}, sum, ::Type{Nothing}, ::Type{Nothing})    => (analyte -> (==(sumcomp(analyte), sum)), :sum)
        (::Type{Nothing}, sum, lcb, acyl)                           => (analyte -> (==(sumcomp(analyte), sum) && 
                                                                                !isnothing(_chain(analyte)) && 
                                                                                isa(_lcb(analyte), lcb)) && 
                                                                                isa(_acyl(analyte), acyl), 
                                                                        :chain)
        (cls, sum, ::Type{Nothing}, ::Type{Nothing})                => (analyte -> (isa(class(analyte), cls) && ==(sumcomp(analyte), sum)), :cpd)
        (cls, sum, lcb, acyl)                                       => (analyte -> (isa(class(analyte), cls) &&
                                                                                ==(sumcomp(analyte), sum) && 
                                                                                !isnothing(_chain(analyte)) && 
                                                                                isa(_lcb(analyte), lcb)) && 
                                                                                isa(_acyl(analyte), acyl),
                                                                        :cpd)
    end
    finish_query!(aquery, sym, qq, fn, findall(fn ∘ qf, aquery), view)
end

function query(aquery::Query, qq::Symbol; view = aquery.view, fn = identity)
    qf = @match qq begin
        :class  => analyte -> ==(analyte.states[states_id(:class)], 1)
        :chain  => analyte -> ==(analyte.states[states_id(:chain)], 1)
        :rt     => analyte -> ==(analyte.states[states_id(:rt)], 1)
        :error  => analyte -> ==(analyte.states[states_id(:error)], 1)
        :isf    => analyte -> ==(analyte.states[states_id(:isf)], 1)
        :class_ => analyte -> ==(analyte.states[states_id(:class)], 0)
        :chain_ => analyte -> ==(analyte.states[states_id(:chain)], 0)
        :rt_    => analyte -> ==(analyte.states[states_id(:rt)], 0)
        :error_ => analyte -> ==(analyte.states[states_id(:error)], 0)
        :isf_   => analyte -> ==(analyte.states[states_id(:isf)], 0)
        :class! => analyte -> ==(analyte.states[states_id(:class)], -1)
        :chain! => analyte -> ==(analyte.states[states_id(:chain)], -1)
        :rt!    => analyte -> ==(analyte.states[states_id(:rt)], -1)
        :error! => analyte -> ==(analyte.states[states_id(:error)], -1)
        :isf!   => analyte -> ==(analyte.states[states_id(:isf)], -1)
        :all    => analyte -> all(==(1), analyte.states)
    end
    finish_query!(aquery, :id, qq, fn, findall(fn ∘ qf, aquery.result), view)
end

function query(aquery::Query, mrm::MRM; 
                view = aquery.view, fn = identity, compatible = false, db_product = SPDB[mrm.polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                mode::Symbol = :default)
    products = @p mrm.mz2 map(id_product(_, mrm.polarity; db = db_product, mz_tol = mrm.mz_tol))
    if compatible
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscompatible(cpd, product) for product in products[i]) && 
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ==(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    else
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscomponent(cpd, product) for product in products[i]) && 
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ==(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    end
    finish_query!(aquery, :validated_by, "mrm", fn, findall(fn ∘ qf, last.(aquery)), view)
end

function coelution(project::Project; analytes = project.analytes, polarity::Bool = true, mz_tol = 0.35, rt_tol = 0.1)
    result = NamedTuple{(:mz1, :rt, :analytes), Tuple{Vector{Float64}, Vector{Float64}, SubArray{AnalyteSP}}}[]
    for analyte in analytes
        frags = last(analyte).fragments
        id = findfirst(i -> isa(frags.mz1[i].adduct, Pos) == polarity, reverse(eachindex(frags)))
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
        (hasisomer(class(analyte)) || hasisomer(_chain(analyte))) ? push!(id2, i) : push!(id1, i)
    end
    dict = Dict{Tuple{<: ClassSP, NTuple{3, Int}, <: Union{Chain, Nothing}}, Vector{Tuple{Float64, Int}}}()

    for id in id1
        push!(get!(dict, (class(analytes[id]), sumcomp(analytes[id]), _chain(analytes[id])), Tuple{Float64, Int}[]), (score_fn(analytes[id].scores), id))
    end

    for id in id2
        sc = score_fn(analytes[id].scores)
        pushed = false
        for (ky, vl) in dict
            sumcomp(analytes[id]) == ky[2] || continue
            (hasisomer(class(analytes[id])) ? in(ky[1], class(analytes[id]).isomer) : ==(ky[1], class(analytes[id]))) || continue
            (isnothing(_chain(analytes[id])) || nhydroxyl(_acyl(analytes[id])) == nhydroxyl(ky[3].acyl)) && 
                (pushed = true; push!(vl, (sc, id)))
        end
        pushed || push!(dict, (class(analytes[id]), sumcomp(analytes[id]), _chain(analytes[id])) => [(sc, id)])
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

function finish_query!(aquery::Query, ql, qr, fn, id, view)
    fn === identity ? push!(aquery.query, ql => qr) : push!(aquery.query, ql => Inv(qr))
    aquery.view = view
    aquery.result = view ? Base.view(aquery, id) : aquery[id]
    aquery
end

query(aquery::Query, qq::Inv; view = true) = query(aquery, qq.arg; view, fn = !)
function query(aquery::Query, qqs::Tuple; view = true) 
    pid = parentindices(aquery.result)[1]
    qs = map(qqs) do qq
        query(Query(aquery.project, view ? Base.view(aquery.project, deepcopy(pid)) : copy_wo_project.(aquery.result), [], view), qq; view)
    end
    qs = union!(qs...)
    push!(aquery.query, tuple(qs.query...))
    aquery.view = view
    aquery.result = view ? qs.result : Vector(qs.result)
    aquery
end

query_raw(project::Project, source, id) = project.data[source].raw[findfirst(==(id), project.data[source].raw.id)]
query_raw(project::Project, source, id, col) = getproperty(project.data[source].raw, col)[findfirst(==(id), project.data[source].raw.id)]
set_raw!(project::Project, source, id, col, val) = (getproperty(project.data[source].raw, col)[findfirst(==(id), project.data[source].raw.id)] = val)
function set_raw_if!(project::Project, source, id, col, val, crit) 
    i = findfirst(==(id), project.data[source].raw.id)
    crit(getproperty(project.data[source].raw, col)[i]) || return
    getproperty(project.data[source].raw, col)[i] = val
end

# filter/comparator