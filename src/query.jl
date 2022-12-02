function query(project::Project; view = true) 
    v = @view project[:]
    v = view ? v : Vector(v)
    Query(project, v, [], view)
end

function query(project::Project, args...; view = true) 
    v = @view project[:]
    v = view ? v : Vector(v)
    query(Query(project, v, [], view), args...)
end

query(aquery::Query, id::Tuple; view = aquery.view, fn = identity) = query(aquery, id...; view, fn)

function query(aquery::Query, quantity::Symbol, low, up; view = aquery.view, fn = identity)
    qf = @match quantity begin
        :rt => analyte -> between(analyte.rt; low, up)
        :mw => analyte -> between(mw(last(analyte)); low, up)
        :mz => analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1; 
                                                low, up) for frag in eachrow(last(analyte).fragments))
    end
    finish_query!(aquery, quantity, (low, up), fn, findall(fn ∘ qf, aquery), view)
end

function query(aquery::Query, quantity::Symbol, param; view = aquery.view, fn = identity)
    id = @match quantity begin
        :topscore => query_score(aquery, param)
        :topsc    => query_score(aquery, param)
    end
    finish_query!(aquery, quantity, param, fn, id, view)
end

function query(aquery::Query, id::Type{<: ClassSP}; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
        last(analyte)
    end
    qf(cpd) = fn(
        isa(cpd.class, id) || 
        hasisomer(cpd.class) && 
        any(isa(isomer, id) for isomer in cpd.class.isomer)
    )
    finish_query!(aquery, :class, id, fn, findall(qf, cpds), view)
end

function query(aquery::Query, id::Lcb; view = aquery.view, fn = identity)
    match = convert_type(LCB, id)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    qf(cpd) = fn(!isnothing(cpd.chain) && isa(cpd.chain.lcb, match))
    finish_query!(aquery, :lcb, id, fn, findall(qf, cpds), view)
end

function query(aquery::Query, id::NACYL; view = aquery.view, fn = identity)
    match = convert_type(ACYL, id)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    qf(cpd) = fn(
        !isnothing(cpd.chain) && 
        isa(cpd.chain.acyl, match) && 
        all((cpd.sum .- sumcomp(cpd.chain.lcb)) .== (id.cb, id.db, id.ox))
    )
    finish_query!(aquery, :acyl, id, fn, findall(qf, cpds), view)
end

function query(aquery::Query, id::Type{<: LCB}; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    qf(cpd) = fn(!isnothing(cpd.chain) && isa(cpd.chain.lcb, id))
    finish_query!(aquery, :lcb, id, fn, findall(qf, cpds), view)
end

function query(aquery::Query, id::CompoundID; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    match = convert_internal(id)
    qf, sym = @match match begin
        (::Type{Nothing}, sum, ::Type{Nothing}, ::Type{Nothing})    => (cpd -> (==(cpd.sum, sum)), :sum)
        (::Type{Nothing}, sum, lcb, acyl)                           => (cpd -> (==(cpd.sum, sum) && 
                                                                                isa(cpd.chain.lcb, lcb)) && 
                                                                                isa(cpd.chain.acyl, acyl), 
                                                                        :chain)
        (class, sum, ::Type{Nothing}, ::Type{Nothing})              => (cpd -> (isa(cpd.class, class) && ==(cpd.sum, sum)), :cpd)
        (class, sum, lcb, acyl)                                     => (cpd -> (isa(cpd.class, class) &&
                                                                                ==(cpd.sum, sum) && 
                                                                                isa(cpd.chain.lcb, lcb)) && 
                                                                                isa(cpd.chain.acyl, acyl),
                                                                        :cpd)
    end
    finish_query!(aquery, sym, id, fn, findall(fn ∘ qf, cpds), view)
end

function query(aquery::Query, id::Symbol; view = aquery.view, fn = identity)
    qf = @match id begin
        :class  => analyte -> ==(analyte.states[1], 1)
        :chain  => analyte -> ==(analyte.states[2], 1)
        :class_ => analyte -> ==(analyte.states[1], 0)
        :chain_ => analyte -> ==(analyte.states[2], 0)
        :class! => analyte -> ==(analyte.states[1], -1)
        :chain! => analyte -> ==(analyte.states[2], -1)
        :both   => analyte -> all(==(1), analyte.states)
    end
    finish_query!(aquery, :id, id, fn, findall(fn ∘ qf, aquery.result), view)
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
        push!(get!(dict, (cpd.class, cpd.sum, cpd.chain), Tuple{Float64, Int}[]), (mean(analytes[id].scores), id))
    end

    for id in id2
        cpd = last(analytes[id])
        pushed = false
        for (ky, vl) in dict
            cpd.sum == ky[2] || continue
            (hasisomer(cpd.class) ? in(ky[1], cpd.class.isomer) : ==(ky[1], cpd.class)) || continue
            (isnothing(cpd.chain) || nhydroxyl(cpd.chain.acyl) == nhydroxyl(ky[3].acyl)) && 
                (pushed = true; push!(vl, (mean(analytes[id].scores), id)))
        end
        pushed || push!(dict, (cpd.class, cpd.sum, cpd.chain) => [(mean(analytes[id].scores), id)])
    end
    final = Int[] 
    len = topN >= 1 ? (vl -> topN) : (vl -> round(Int, length(vl) * topN * (1 + eps(Float64))))
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

query(aquery::Query, id::Inv; view = true) = query(aquery, id.args...; view, fn = !)
function query(aquery::Query, ids::Vararg{Union{Tuple, Inv, Type, CompoundID}, N}; view = true) where N 
    pid = parentindices(aquery.result)[1]
    qs = map(ids) do id
        query(Query(aquery.project, view ? Base.view(aquery.project, deepcopy(pid)) : copy_wo_project.(aquery.result), [], view), id; view)
    end
    qs = union!(qs...)
    push!(aquery.query, tuple(qs.query...))
    aquery.view = view
    aquery.result = view ? qs.result : Vector(qs.result)
    aquery
end

query_raw(project::Project, source, id) = @views project.data[source].raw[findfirst(==(id), project.data[source].raw.id), :]

# filter/comparator