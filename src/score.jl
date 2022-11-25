macro score(targetconvert, weight, objective, threshold)
    target, converter = @match targetconvert begin
        :class                                  => (:class, :first)   
        :chain                                  => (:chain, :last)   
        Expr(:call, :(=>), target, converter)   => (target, converter)
    end # any 1-arg fn: Vector -> Int 
    parameters = (target = target, converter = converter, weight = weight, objective = objective, threshold = threshold)
    weight = @match weight begin
        1 => :w1
        _ => weight
    end # any 2-arg fn: (analyte, cpd) -> Number
    replace_int = @λ begin
        x::Int  => Expr(:call, :get!, :dict, x, 0)
        :all    => Expr(:call, :sum, Expr(:call, :values, :dict))
        x       => x
    end
    objective = @match objective begin
        ::Expr  => Expr(objective.head, map(replace_int, objective.args)...)
        ::Int   => replace_int(objective)
    end
    threshold = Expr(threshold.head, threshold.args[1], :x, threshold.args[3])
    return quote
        transform_score(dict) = $objective
        apply_threshold(x) = $threshold
        function apply_score(analyte::AnalyteSP)
            dict = Dict{Int, Float64}()
            for cpd in analyte
                uw = $converter(cpd.states)
                dict[uw] = get!(dict, uw, 0) + $weight(analyte, cpd)
            end
            transform_score(dict)
        end
        Score(0.0, $parameters, apply_score, apply_threshold)
    end
end
filter_score!(aquery::Query, score = nothing) = (filter_score!(aquery.project, score; analytes = aquery.result); aquery)
function filter_score!(project::Project, score = nothing; analytes = project.analytes)
    printstyled("Score> \n", color = :green, bold = true)
    isnothing(score) || for analyte in analytes
        push!(analyte.scores, deepcopy(score))
    end 
    for analyte in analytes
        apply_score!(analyte)
        apply_threshold!(analyte)
    end
    project
end

apply_score(analyte::AnalyteSP) = last(analyte.scores).apply_score(analyte)
apply_score!(analyte::AnalyteSP) = (last(analyte.scores).score = apply_score(analyte)) 
apply_threshold(analyte::AnalyteSP) = last(analyte.scores).apply_threshold(analyte)
function apply_threshold!(analyte::AnalyteSP)
    sc = last(analyte.scores)
    id = @match sc.parameters.target begin
            :class => 1:1
            :chain => 2:2
            :both  => 1:2
    end
    r = sc.apply_threshold(sc.score) 
    r || (analyte.states[id] .= -1)
    r
end

nfrags(analyte::AnalyteSP, cpd::CompoundSP) = nfrags(cpd)
w1(analyte::AnalyteSP, cpd::CompoundSP) = 1

default_score = @score chain 1 -1 s <= 1

select_score!(aquery::Query, score = nothing; topN = 0.5, by = :cpd) = select_score!(aquery.project, score; aquery = aquery, analytes = aquery.result, view = aquery.view, topN = topN, by = by)
function select_score!(project::Project, score = nothing; topN = 0.5, by = :cpd, analytes = project.analytes, aquery = nothing, view = true)
    isnothing(score) || begin
        for analyte in analytes
            push!(analyte.scores, deepcopy(score))
        end 
        apply_score!(analytes)
    end
    # no isomer => has isomer
    id1 = Int[]
    id2 = Int[]
    for (i, analyte) in enumerate(analytes)
        (hasisomer(last(analyte).class) || hasisomer(last(analyte).chain)) ? push!(id2, i) : push!(id1, i)
    end
    dict = Dict{Tuple{<: ClassSP, NTuple{3, Int}, <: Union{Chain, Nothing}}, Vector{Tuple{Float64, Int}}}()

    for id in id1
        cpd = last(analytes[id])
        push!(get!(dict, (cpd.class, cpd.sum, cpd.chain), Tuple{Float64, Int}[]), (last(analytes[id].scores).score, id))
    end

    for id in id2
        cpd = last(analytes[id])
        pushed = false
        for (ky, vl) in dict
            cpd.sum == ky[2] || continue
            (hasisomer(cpd.class) ? in(ky[1], cpd.class.isomer) : ==(ky[1], cpd.class)) || continue
            (isnothing(cpd.chain) || nhydroxyl(cpd.chain.acyl) == nhydroxyl(ky[3].acyl)) && (pushed = true; push!(vl, (last(analytes[id].scores).score, id)))
        end
        pushed || push!(dict, (cpd.class, cpd.sum, cpd.chain) => [(last(analytes[id].scores).score, id)])
    end
    final = Int[] 
    len = @match topN begin
        ::Int       => (vl -> topN)
        ::Float64   => (vl -> round(Int, length(vl) * topN * (1 + eps(Float64))))
    end
    for vl in values(dict)
        sort!(vl, by = first)
        for _ in 1:len(vl)
            isempty(vl) && break
            push!(final, pop!(vl)[2])
        end
    end
    isnothing(aquery) && (aquery = query(project; view))
    aquery.result = view ? (@view project[final]) : project[final]
    push!(aquery.query, :select_score => Symbol("top$(topN)"))
    aquery
end

