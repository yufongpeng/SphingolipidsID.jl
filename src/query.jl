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

function query(aquery::Query, id::Symbol, type; view = aquery.view, fn = identity)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    v = @match id begin
        :class  => @view aquery[findall(cpd -> isa(cpd.class, type)|>fn, cpds)]
        :sum    => @view aquery[findall(cpd -> ==(cpd.sum, type)|>fn, cpds)]
        :lcb    => @view aquery[findall(cpd -> (!isnothing(cpd.chain) && ==(cpd.chain.lcb, type()))|>fn, cpds)]
    end
    finish_query!(aquery, id, type, fn, v, view)
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

function finish_query!(aquery::Query, ql::Symbol, qr, fn, v, view)
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
function query(aquery::Query, ids::Vararg{Union{Tuple, Inv}, N}; view = true) where N 
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