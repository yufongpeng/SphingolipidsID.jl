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
        x::Expr => Expr(x.head, map(replace_int, x.args)...)
        x       => x
    end
    objective = replace_int(objective)
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

apply_score!(object, fn::T; kwargs...) where T <: Function = apply_score!(fn, object; kwargs...)
apply_score!(score_fn::T, aquery::AbstractQuery) where T <: Function = (apply_score!(score_fn, aquery.project; analytes = aquery.result); aquery)
function apply_score!(score_fn::T, project::Project; analytes = project.analytes) where T <: Function
    printstyled("Score> ", color = :green, bold = true)
    print(score_fn, "\n")
    for analyte in analytes
        apply_score!(score_fn, analyte)
    end
    project
end
apply_threshold!(object, fn::T; kwargs...) where T <: Function = apply_threshold!(fn, object; kwargs...)
apply_threshold!(thresh_fn::T, aquery::AbstractQuery) where T <: Function = (apply_threshold!(thresh_fn, aquery.project; analytes = aquery.result); aquery)
function apply_threshold!(thresh_fn::T, project::Project; analytes = project.analytes) where T <: Function
    printstyled("Threshold> ", color = :green, bold = true)
    print(thresh_fn, "\n")
    for analyte in analytes
        apply_threshold!(thresh_fn, analyte)
    end
    project
end

apply_score(score_fn::T, analyte::AnalyteSP) where T <: Function = score_fn(analyte)
apply_score!(score_fn::T, analyte::AnalyteSP) where T <: Function = (analyte.scores = score_fn(analyte)) 
apply_threshold(thresh_fn, analyte::AnalyteSP) = thresh_fn(analyte)
function apply_threshold!(thresh_fn, analyte::AnalyteSP)
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