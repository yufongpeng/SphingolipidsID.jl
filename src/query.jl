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
    v = @match quantity begin
        :rt => @view aquery[findall(analyte -> between(analyte.rt; low, up)|>fn, aquery)]
        :mw => @view aquery[findall(analyte -> between(mw(last(analyte)); low, up)|>fn, aquery)]
        :mz => @view aquery[findall(analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1; low, up) for frag in eachrow(last(analyte).fragments))|>fn, aquery)]
    end
    finish_query!(aquery, quantity, (low, up), fn, v, view)
end

function query(aquery::Query, id::Type{<: ClassSP}; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
        last(analyte)
    end
    v = @view aquery[findall(cpd -> isa(cpd.class, id)|>fn, cpds)]
    finish_query!(aquery, :class, id, fn, v, view)
end

function query(aquery::Query, id::Lcb; view = aquery.view, fn = identity)
    match = convert_type(LCB, id)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    v = @view aquery[findall(cpd -> (!isnothing(cpd.chain) && isa(cpd.chain.lcb, match))|>fn, cpds)]
    finish_query!(aquery, :lcb, id, fn, v, view)
end

function query(aquery::Query, id::NACYL; view = aquery.view, fn = identity)
    match = convert_type(ACYL, id)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    v = @view aquery[findall(cpd -> (!isnothing(cpd.chain) && isa(cpd.chain.acyl, match) && (cpd.sum .- sumcomp(cpd.chain.lcb)) .== (id.cb, id.db, id.ox))|>fn, cpds)]
    finish_query!(aquery, :acyl, id, fn, v, view)
end

function query(aquery::Query, id::Type{<: LCB}; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    v = @view aquery[findall(cpd -> (!isnothing(cpd.chain) && isa(cpd.chain.lcb, id))|>fn, cpds)]
    finish_query!(aquery, :lcb, id, fn, v, view)
end

function query(aquery::Query, id::CompoundID; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    match = convert_internal(id)
    v, sym = @match match begin
        (::Type{Nothing}, sum, ::Type{Nothing}, ::Type{Nothing})    => (Base.view(aquery, findall(cpd -> (==(cpd.sum, sum))|>fn, cpds)), :sum)
        (::Type{Nothing}, sum, lcb, acyl)                           => (Base.view(aquery, findall(cpd -> (==(cpd.sum, sum) && isa(cpd.chain.lcb, lcb)) && isa(cpd.chain.acyl, acyl)|>fn, cpds)), :chain)
        (class, sum, ::Type{Nothing}, ::Type{Nothing})              => (Base.view(aquery, findall(cpd -> (isa(cpd.class, class) && ==(cpd.sum, sum))|>fn, cpds)), :cpd)
        (class, sum, lcb, acyl)                                     => (Base.view(aquery, findall(cpd -> (isa(cpd.class, class) && ==(cpd.sum, sum) && isa(cpd.chain.lcb, lcb)) && isa(cpd.chain.acyl, acyl)|>fn, cpds)), :cpd)
    end
    finish_query!(aquery, sym, id, fn, v, view)
end

function query(aquery::Query, id::Symbol; view = aquery.view, fn = identity)
    v = @match id begin
        :class  => @view aquery[findall(analyte -> ==(analyte.states[1], 1)|>fn, aquery.result)]
        :chain  => @view aquery[findall(analyte -> ==(analyte.states[2], 1)|>fn, aquery.result)]
        :class_  => @view aquery[findall(analyte -> ==(analyte.states[1], 0)|>fn, aquery.result)]
        :chain_  => @view aquery[findall(analyte -> ==(analyte.states[2], 0)|>fn, aquery.result)]
        :class!  => @view aquery[findall(analyte -> ==(analyte.states[1], -1)|>fn, aquery.result)]
        :chain!  => @view aquery[findall(analyte -> ==(analyte.states[2], -1)|>fn, aquery.result)]
        :both   => @view aquery[findall(analyte -> all(==(1), analyte.states)|>fn, aquery.result)]
    end
    finish_query!(aquery, :id, id, fn, v, view)
end

function finish_query!(aquery::Query, ql, qr, fn, v, view)
    if fn === identity
        push!(aquery.query, ql => qr)
    else
        push!(aquery.query, ql => Inv(qr))
    end
    aquery.view = view
    aquery.result = view ? v : Vector(v)
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