function query(project::Project, args...; view = true) 
    v = @view project.analytes[:]
    v = view ? v : Vector(v)
    query(Query(project, v, [], view), args...)
end


function query(aquery::Query, quantity::Symbol, low, up; view = aquery.view)
    v = @views @match quantity begin
        :rt => aquery.result[findall(analyte -> between(analyte.rt; low, up), aquery.result)]
        :mw => aquery.result[findall(analyte -> between(mw(last(analyte.identification)); low, up), aquery.result)]
        :mz => aquery.result[findall(analyte -> any(between(query_raw(aquery.project, frag.source, frag.id).mz; low, up) for frag in eachrow(last(analyte.identification).fragments)), aquery.result)]
    end
    push!(aquery.query, quantity => (low, up))
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

function query(aquery::Query, type::Symbol, object; view = true)
    cpds = map(aquery.result) do analyte
            last(analyte.identification)
    end
    v = @views @match type begin
        :class  => aquery.result[findall(cpd -> cpd.class isa object, cpds)]
        :sum    => aquery.result[findall(cpd -> cpd.sum == collect(object), cpds)]
        :lcb    => aquery.result[findall(cpd -> !isnothing(cpd.chain) && cpd.chain.lcb == object(), cpds)]
    end
    push!(aquery.query, type => object)
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

function query(aquery::Query, type::Symbol; view = true)
    cpds = map(aquery.result) do analyte
        last(analyte.identification)
    end
    v = @views @match type begin
        :class  => aquery.result[findall(cpd -> isnothing(cpd.sub_rule[1]), cpds)]
        :chain  => aquery.result[findall(cpd -> isnothing(cpd.sub_rule[1]), cpds)]
        :both   => aquery.result[findall(cpd -> all(isnothing, cpd.sub_rule), cpds)]
    end
    push!(aquery.query, :id => type)
    aquery.view = view
    aquery.result = view ? v : Vector(v)
    aquery
end

query_raw(project::Project, source, id) = @views project.data[source].raw[findfirst(==(id), project.data[source].raw.id), :]