SPDB[:SCORE] = Dict{NamedTuple, Function}()

macro score(targetconvert, weight, objective)
    target, converter = @match targetconvert begin
        :class                                  => (:class, :first)   
        :chain                                  => (:chain, :last)   
        Expr(:call, :(=>), target, converter)   => (target, converter)
    end # any 1-arg fn: Vector -> Int 
    parameters = (target = target, converter = converter, weight = weight, objective = objective)
    haskey(SPDB[:SCORE], parameters) && return quote SPDB[:SCORE][$parameters] end
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
    #threshold = Expr(threshold.head, threshold.args[1], :x, threshold.args[3])
    return quote
        transform_score(dict) = $objective
        transform_score_v(dict) = @match $parameters.target begin
            :class => (transform_score(dict), NaN) 
            :chain => (NaN, transform_score(dict)) 
            :both  => (transform_score(dict), transform_score(dict)) 
        end
        function score_fn(analyte::AnalyteSP)
            dict = Dict{Int, Float64}()
            for cpd in analyte
                uw = $converter(cpd.states)
                dict[uw] = get!(dict, uw, 0) + $weight(analyte, cpd)
            end
            transform_score_v(dict)
        end
        push!(SPDB[:SCORE], $parameters => score_fn)
        SPDB[:SCORE][$parameters]
    end
end
apply_score!(aquery::Query, score_fn) = (apply_score!(aquery.project, score_fn; analytes = aquery.result); aquery)
function apply_score!(project::Project, score_fn; analytes = project.analytes)
    printstyled("Score> \n", score_fn, "\n", color = :green, bold = true)
    for analyte in analytes
        apply_score!(analyte, score_fn)
    end
    project
end
apply_threshold!(aquery::Query, thresh_fn) = (apply_threshold!(aquery.project, thresh_fn; analytes = aquery.result); aquery)
function apply_threshold!(project::Project, thresh_fn; analytes = project.analytes)
    printstyled("Threshold> ", thresh_fn, "\n", color = :green, bold = true)
    for analyte in analytes
        apply_threshold!(analyte, thresh_fn)
    end
    project
end

apply_score(analyte::AnalyteSP, score_fn) = score_fn(analyte)
apply_score!(analyte::AnalyteSP, score_fn) = (analyte.scores = score_fn(analyte)) 
apply_threshold(analyte::AnalyteSP, thresh_fn) = thresh_fn(analyte)
function apply_threshold!(analyte::AnalyteSP, thresh_fn)
    id = @match analyte.scores begin
            (x, x) && if isnan(x) end   => 1:0     
            (_, x) && if isnan(x) end   => 1:1       
            (x, _) && if isnan(x) end   => 2:2
            (_, _)          => 1:2
    end
    map(id) do i
        r = thresh_fn(analyte.scores[i])
        r || (analyte.states[i] = -1)
        r
    end
end

nfrags(analyte::AnalyteSP, cpd::CompoundSP) = nfrags(cpd)
w1(analyte::AnalyteSP, cpd::CompoundSP) = 1