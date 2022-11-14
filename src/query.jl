function query(project::Project, args...; view = true) 
    v = @view project[:]
    v = view ? v : Vector(v)
    query(Query(project, v, [], view), args...)
end


function query(aquery::Query, quantity::Symbol, low, up; view = aquery.view)
    v = @views @match quantity begin
        :rt => aquery[findall(analyte -> between(analyte.rt; low, up), aquery)]
        :mw => aquery[findall(analyte -> between(mw(last(analyte)); low, up), aquery)]
        :mz => aquery[findall(analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz1; low, up) for frag in eachrow(last(analyte).fragments)), aquery)]
    end
    push!(aquery.query, quantity => (low, up))
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

function query(aquery::Query, id::Symbol, type; view = true)
    cpds = map(aquery) do analyte
            last(analyte)
    end
    v = @views @match id begin
        :class  => aquery[findall(cpd -> cpd.class isa type, cpds)]
        :sum    => aquery[findall(cpd -> cpd.sum == type, cpds)]
        :lcb    => aquery[findall(cpd -> !isnothing(cpd.chain) && cpd.chain.lcb == type(), cpds)]
    end
    push!(aquery.query, id => type)
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

function query(aquery::Query, id::Symbol; view = true)
    cpds = map(aquery) do analyte
        last(analyte)
    end
    v = @views @match id begin
        :class  => aquery[findall(cpd -> isa(cpd.states[1], Symbol), cpds)]
        :chain  => aquery[findall(cpd -> isa(cpd.states[2], Symbol), cpds)]
        :class_fail  => aquery[findall(cpd -> isnothing(cpd.states[1]), cpds)]
        :chain_fail  => aquery[findall(cpd -> isnothing(cpd.states[2]), cpds)]
        :class_wait  => aquery[findall(cpd -> !isnothing(cpd.states[1]) && !isa(cpd.states[1], Symbol), cpds)]
        :chain_wait  => aquery[findall(cpd -> !isnothing(cpd.states[2]) && !isa(cpd.states[2], Symbol), cpds)]
        :both   => aquery[findall(cpd -> all(r -> isa(r, Symbol), cpd.states), cpds)]
    end
    push!(aquery.query, :id => id)
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

query_raw(project::Project, source, id) = @views project.data[source].raw[findfirst(==(id), project.data[source].raw.id), :]