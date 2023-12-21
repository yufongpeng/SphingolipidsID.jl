SPDB[:SCORE] = (param = NamedTuple[], fn = Function[])
"""
    @cpdsc(targetconvert, weight, objective)

Create a score function based on compound.

* `targetconvert`: the target of score.
    * `:class`
    * `:chain`
    * `:both => converter`: converter is a 1-arg function where the parameter is the state vector of compounds.
* `weight`: the weight function for each compound. `1` indicates identity; otherwise, a full function should be given. The function should accept analyte and compounds as inputs and a number as output.
* `objective`: a combination of `1`, `0`, `-1` and other operations which generates a function for calculating score.
"""
macro cpdsc(targetconvert, weight, objective)
    target, converter = @match targetconvert begin
        :class                                  => (:class, :first)
        :chain                                  => (:chain, :last)
        Expr(:call, :(=>), target, converter)   => (target, converter)
    end # any 1-arg fn: Vector -> Int
    parameters = (target = target, converter = converter, weight = weight, objective = objective)
    id = findfirst(==(parameters), SPDB[:SCORE].param) 
    isnothing(id) || return quote SPDB[:SCORE].fn[$id] end
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
        push!(SPDB[:SCORE].param, $parameters)
        i = lastindex(SPDB[:SCORE].param)
        fn = (analyte::AnalyteSP, mutate::Bool) -> begin
            dict = Dict{Int, Float64}()
            for cpd in analyte
                uw = $converter(cpd.state)
                dict[uw] = get!(dict, uw, 0) + $weight(analyte, cpd)
            end
            result = i => transform_score(dict)
            mutate ? (analyte.cpdsc = result) : result
        end
        push!(SPDB[:SCORE].fn, fn)
        SPDB[:SCORE].fn[i]
    end
end

"""
    @score(targetobjective)

Create a score function based on analyte.

* `targetobjective`: the target of score.
    * `:all`: sum up all states except `:total` and `:manual`.
    * `target::Symbol`: use the state of `target`.
    * `target::Vector`: sum up selected states. 
    * `anything => objective`: use objective which is a 1-arg function taking analyte as input and returning a number.
"""
macro score(targetobjective)
    target, objective = @match targetobjective begin
        :all                                    => (:all, :(x -> sum(x.state[1:end - 2])))
        target::Symbol                          => (target, :(x -> x.state[state_id[target]]))
        target::Vector                          => (target, :(x -> sum(x.state[state_id[i]] for i in target)))
        Expr(:call, :(=>), target, objective)   => (target, objective)
    end # any 1-arg fn: AnalyteSP -> Number
    parameters = (target = target, objective = objective)
    return quote
        id = findfirst(==($parameters), SPDB[:SCORE].param)
        isnothing(id) && begin
            push!(SPDB[:SCORE].param, $parameters)
            id = lastindex(SPDB[:SCORE].param)
            push!(SPDB[:SCORE].fn, (analyte, mutate) -> mutate ? (analyte.score = id => $objective(analyte)) : (id => $objective(analyte)))
        end
        SPDB[:SCORE].fn[id]
    end
end

"""
    apply_score!(fn::Function, object::AbstractQuery; analyte)
    apply_score!(fn::Function, object::Project; analyte)
    apply_score!(fn::Function, object::AnalyteSP; analyte)
    apply_score!(object, fn::Function; analyte)

Apply score function `fn` to `object` and store the result. 
"""
apply_score!(object, fn::T; kwargs...) where T <: Function = apply_score!(fn, object; kwargs...)
apply_score!(score_fn!::T, aquery::AbstractQuery) where T <: Function = (apply_score!(score_fn!, aquery.project; analyte = aquery.result); aquery)
function apply_score!(score_fn!::T, project::Project; analyte = project.analyte) where T <: Function
    @p analyte foreach(score_fn!(_, true))
    project
end

"""
    calc_score(fn::Function, object::AbstractQuery; analyte)
    calc_score(fn::Function, object::Project; analyte)
    calc_score(fn::Function, object::AnalyteSP; analyte)
    calc_score(object, fn::Function; analyte)

Apply score function `fn` to `object` and return the results. 
"""
calc_score(object, fn::T; kwargs...) where T <: Function = calc_score(fn, object; kwargs...)
calc_score(score_fn::T, aquery::AbstractQuery) where T <: Function = calc_score(score_fn, aquery.project; analyte = aquery.result)
calc_score(score_fn::T, project::Project; analyte = project.analyte) where T <: Function = @p analyte map(score_fn(_, false))
calc_score(score_fn::T, analyte::AnalyteSP) where T <: Function = score_fn(analyte, false)

"""
    apply_threshold!(fn::Function, object::AbstractQuery; cpd = false, analyte)
    apply_threshold!(fn::Function, object::Project; cpd = false, analyte)
    apply_threshold!(fn::Function, object::AnalyteSP; cpd = false, analyte)
    apply_threshold!(object, fn::Function; cpd = false, analyte)

Apply threshold function `fn` to `object` and store the result. 
    
It uses `analyte.cpdsc` and mutates `analyte.state[state_id(:class)]`, `analyte.state[state_id(:chain)]`, or both if `cpd` is true; otherwise, it uses `analyte.score` and mutates `analyte.state[state_id(:total)]`.
"""
apply_threshold!(object, fn::T; kwargs...) where T <: Function = apply_threshold!(fn, object; kwargs...)
apply_threshold!(thresh_fn::T, aquery::AbstractQuery; cpd = false) where T <: Function = (apply_threshold!(thresh_fn, aquery.project; cpd, analyte = aquery.result); aquery)
function apply_threshold!(thresh_fn::T, project::Project; cpd = false, analyte = project.analyte) where T <: Function
    @p analyte foreach(apply_threshold!(thresh_fn, _; cpd))
    project
end

function apply_threshold!(thresh_fn::T, analyte::AnalyteSP; cpd = false) where T <: Function
    cpd || return (analyte.state[state_id(:total)] = thresh_fn(analyte.score.second) ? 1 : -1; analyte)
    result = thresh_fn(analyte.cpdsc.second) ? 1 : -1
    @match SPDB[:SCORE].param[analyte.cpdsc.first].target begin
        :class => (analyte.state[state_id(:class)] = min(analyte.state[state_id(:class)], result))
        :chain => (analyte.state[state_id(:chain)] = min(analyte.state[state_id(:chain)], result))
        :both  => (analyte.state[state_id.([:class, :chain])] .= min.(analyte.state[state_id.([:class, :chain])], result))
    end
    analyte
end

"""
    calc_threshold(fn::Function, object::AbstractQuery; cpd = false, analyte)
    calc_threshold(fn::Function, object::Project; cpd = false, analyte)
    calc_threshold(fn::Function, object::AnalyteSP; cpd = false, analyte)
    calc_threshold(object, fn::Function; cpd = false, analyte)

Apply threshold function `fn` to `object` and return the result.

It uses `analyte.cpdsc` and compares with `analyte.state[state_id(:class)]`, `analyte.state[state_id(:chain)]`, or both if `cpd` is true; otherwise, it uses `analyte.score`.
"""
calc_threshold(object, fn::T; kwargs...) where T <: Function = calc_threshold(fn, object; kwargs...)
calc_threshold(thresh_fn::T, aquery::AbstractQuery; cpd = false) where T <: Function = calc_threshold(thresh_fn, aquery.project; cpd, analyte = aquery.result)
calc_threshold(thresh_fn::T, project::Project; cpd = false, analyte = project.analyte) where T <: Function = @p analyte map(calc_threshold(thresh_fn, _; cpd))

function calc_threshold(thresh_fn::T, analyte::AnalyteSP; cpd = false) where T <: Function
    cpd || return (thresh_fn(analyte.score.second) ? 1 : -1)
    result = thresh_fn(analyte.cpdsc.second) ? 1 : -1
    @match SPDB[:SCORE].param[analyte.cpdsc.first].target begin
        :class => min(analyte.state[state_id(:class)], result)
        :chain => min(analyte.state[state_id(:chain)], result)
        :both  => min(analyte.state[state_id.([:class, :chain])]..., result)
    end
    analyte
end
"""
    calc_threshold(fn::Function, scores::Vector{<: Pair})

Apply threshold function to the output of `calc_score`.
"""
calc_threshold(thresh_fn::T, scores::Vector{<: Pair}) where T <: Function = @p scores map(thresh_fn(_.second))

# weights
"""
    w_nfrags(analyte::AnalyteSP, cpd::AbstractCompoundSP)

Number of fragments of `cpd`, i.e., `size(cpd.fragment, 1)`.
"""
w_nfrags(analyte::AnalyteSP, cpd::CompoundSP) = nfrags(cpd)
"""
    w_rank(analyte::AnalyteSP, cpd::CompoundSP)

The index of `cpd`.
"""
w_rank(analyte::AnalyteSP, cpd::CompoundSP) = findfirst(==(cpd), analyte)
"""
    w_rank_log(analyte::AnalyteSP, cpd::CompoundSP)

The log of index of `cpd`.
"""
w_rank_log(analyte::AnalyteSP, cpd::CompoundSP) = log(findfirst(==(cpd), analyte))
"""
    w_rank_exp(analyte::AnalyteSP, cpd::CompoundSP)

The exponential of index of `cpd`.
"""
w_rank_exp(analyte::AnalyteSP, cpd::CompoundSP) = exp(findfirst(==(cpd), analyte))
w1(analyte::AnalyteSP, cpd::CompoundSP) = 1