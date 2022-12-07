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

function query(aquery::Query, qq::Pair{Symbol, <: Any}; view = aquery.view, fn = identity)
    score = @match qq.first begin
        :topscore => true
        :topsc    => true
        _         => false
    end
    score && return finish_query!(aquery, qq.first, qq.second, fn, query_score(aquery, qq.second), view)
    qf = @match qq.first begin
        :rt => analyte -> between(analyte.rt, qq.second)
        :mw => analyte -> between(mw(last(analyte)), qq.second)
        :mz => analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1, 
                                                qq.second) for frag in last(analyte).fragments)
    end
    finish_query!(aquery, qq.first, qq.second, fn, findall(fn ∘ qf, aquery), view)
end

function query(aquery::Query, qq::Type{<: ClassSP}; view = aquery.view, fn = identity)
    qf(cpd) = fn(
        isa(cpd.class, qq) || 
        hasisomer(cpd.class) && 
        any(isa(isomer, qq) for isomer in cpd.class.isomer)
    )
    finish_query!(aquery, :class, qq, fn, findall(qf, last.(aquery)), view)
end

function query(aquery::Query, qq::Lcb; view = aquery.view, fn = identity)
    match = convert_type(LCB, qq)
    qf(cpd) = fn(!isnothing(cpd.chain) && isa(cpd.chain.lcb, match))
    finish_query!(aquery, :lcb, qq, fn, findall(qf, last.(aquery)), view)
end

function query(aquery::Query, qq::NACYL; view = aquery.view, fn = identity)
    match = convert_type(ACYL, qq)
    qf(cpd) = fn(
        !isnothing(cpd.chain) && 
        isa(cpd.chain.acyl, match) && 
        all((cpd.sum .- sumcomp(cpd.chain.lcb)) .== (qq.cb, qq.db, qq.ox))
    )
    finish_query!(aquery, :acyl, qq, fn, findall(qf, last.(aquery)), view)
end

function query(aquery::Query, qq::Type{<: LCB}; view = aquery.view, fn = identity)
    qf(cpd) = fn(!isnothing(cpd.chain) && isa(cpd.chain.lcb, qq))
    finish_query!(aquery, :lcb, qq, fn, findall(qf, last.(aquery)), view)
end

function query(aquery::Query, qq::CompoundID; view = aquery.view, fn = identity)
    match = convert_internal(qq)
    qf, sym = @match match begin
        (::Type{Nothing}, sum, ::Type{Nothing}, ::Type{Nothing})    => (cpd -> (==(cpd.sum, sum)), :sum)
        (::Type{Nothing}, sum, lcb, acyl)                           => (cpd -> (==(cpd.sum, sum) && 
                                                                                !isnothing(cpd.chain) && 
                                                                                isa(cpd.chain.lcb, lcb)) && 
                                                                                isa(cpd.chain.acyl, acyl), 
                                                                        :chain)
        (class, sum, ::Type{Nothing}, ::Type{Nothing})              => (cpd -> (isa(cpd.class, class) && ==(cpd.sum, sum)), :cpd)
        (class, sum, lcb, acyl)                                     => (cpd -> (isa(cpd.class, class) &&
                                                                                ==(cpd.sum, sum) && 
                                                                                !isnothing(cpd.chain) && 
                                                                                isa(cpd.chain.lcb, lcb)) && 
                                                                                isa(cpd.chain.acyl, acyl),
                                                                        :cpd)
    end
    finish_query!(aquery, sym, qq, fn, findall(fn ∘ qf, last.(aquery)), view)
end

function query(aquery::Query, qq::Symbol; view = aquery.view, fn = identity)
    qf = @match qq begin
        :class  => analyte -> ==(analyte.states[1], 1)
        :chain  => analyte -> ==(analyte.states[2], 1)
        :class_ => analyte -> ==(analyte.states[1], 0)
        :chain_ => analyte -> ==(analyte.states[2], 0)
        :class! => analyte -> ==(analyte.states[1], -1)
        :chain! => analyte -> ==(analyte.states[2], -1)
        :both   => analyte -> all(==(1), analyte.states)
    end
    finish_query!(aquery, :id, qq, fn, findall(fn ∘ qf, aquery.result), view)
end

function query_score(aquery::Query, topN = 0.5)
    analytes = aquery.result
    # no isomer => has isomer
    id1 = Int[]
    id2 = Int[]
    for (i, analyte) in enumerate(analytes)
        (hasisomer(last(analyte).class) || hasisomer(last(analyte).chain)) ? push!(id2, i) : push!(id1, i)
    end
    dict = Dict{Tuple{<: ClassSP, NTuple{3, Int}, <: Union{Chain, Nothing}}, Vector{Tuple{Float64, Int}}}()

    for id in id1
        cpd = last(analytes[id])
        push!(get!(dict, (cpd.class, cpd.sum, cpd.chain), Tuple{Float64, Int}[]), (mean(filter(!isnan, analytes[id].scores)), id))
    end

    for id in id2
        cpd = last(analytes[id])
        pushed = false
        for (ky, vl) in dict
            cpd.sum == ky[2] || continue
            (hasisomer(cpd.class) ? in(ky[1], cpd.class.isomer) : ==(ky[1], cpd.class)) || continue
            (isnothing(cpd.chain) || nhydroxyl(cpd.chain.acyl) == nhydroxyl(ky[3].acyl)) && 
                (pushed = true; push!(vl, (mean(filter(!isnan, analytes[id].scores)), id)))
        end
        pushed || push!(dict, (cpd.class, cpd.sum, cpd.chain) => [(mean(filter(!isnan, analytes[id].scores)), id)])
    end
    final = Int[] 
    len = topN >= 1 ? (vl -> topN) : (vl -> max(1, round(Int, length(vl) * topN * (1 + eps(Float64)))))
    for vl in values(dict)
        sort!(vl, by = first)
        for _ in 1:len(vl)
            isempty(vl) && break
            push!(final, pop!(vl)[2])
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

# filter/comparator