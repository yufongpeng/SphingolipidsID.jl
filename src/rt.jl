"""
    group_analyte(analytes)

Group `analytes` based on `class`.
"""
group_analyte(analytes) = @p analytes groupview(deisomerized(class(_))) map((first ∘ parentindices)(_))
"""
    initialize_clusters!(aquery::AbstractQuery; kwargs...)
    initialize_clusters!(project::Project; analytes = project.analytes)

Initialize `project.appendix[:clusters]` with `group_analyte(analytes)`.
"""
initialize_clusters!(aquery::AbstractQuery; kwargs...) = (initialize_clusters!(aquery.project; analytes = aquery.result, kwargs...); aquery)
initialize_clusters!(project::Project; analytes = project.analytes) = replace_clusters!(project, group_analyte(analytes))
"""
    analytes2clusters!(aquery::AbstractQuery; kwargs...)
    analytes2clusters!(project::Project; new = false, scale = 0.0, radius = 0.0, analytes = project.analytes, kwargs...)

Run dbscan algorithm on `analytes`, and put all results in `project.appendix[:clusters_possible]`.

# Keyword Arguments
* `new`: whether clear all existing clusters in `project.appendix[:clusters_possible]` or not.
* `scale`: the normalizing factor for mw.
* `radius`: radius of dbscan.

See `dbscan` for more detailed settings.
"""
analytes2clusters!(aquery::AbstractQuery; kwargs...) = (analytes2clusters!(aquery.project; analytes = aquery.result, kwargs...); aquery)
function analytes2clusters!(project::Project; new = false, scale = 0.0, radius = 0.0, analytes = project.analytes, kwargs...)
    valid_id = collect((first ∘ parentindices)(analytes))
    groups = @p project.appendix[:clusters] map(filter(x -> in(x, valid_id), _)) filter(length(_) > 2)
    ret = @p groups map(map(rt, @views project.analytes[_]))
    maxrt = maximum(maximum(r) for r in ret)
    mass = @p groups map(map(mw, @views project.analytes[_]))
    scale = scale > 0 ? scale : median(map(m -> quantile(m, 0.9) - quantile(m, 0.1), mass)) * 2 / maxrt
    mass = @p mass map(_ / scale)
    radius = radius > 0 ? radius : maxrt / 20 * sqrt(2)
    clusters = @p zip(ret, mass) map(dbscan(hcat(_...)', radius; kwargs...)) map(map(getproperty(:core_indices), _.clusters))
    clusters = @p zip(clusters, groups) map(map(x -> _[2][x], _[1]))
    clusters_possible = get!(project.appendix, :clusters_possible, Dictionary{ClassSP, Vector{Vector{Int}}}())
    new && empty!(clusters_possible)
    @p zip(keys(groups), clusters) foreach(union!(empty!(get!(clusters_possible, _[1], Vector{Int}[])), _[2]))
    printstyled("Possible clusters:\n", color = :green)
    display(project.appendix[:clusters_possible])
    println()
    project
end

@as_record UnitRange
@as_record StepRange
convert_negidx(x, i) = i > 0 ? i : lastindex(x) + i
index_fn_v(i::Vector, n::Int) = map(index_fn, i)
index_fn_v(i, n::Int) = repeat([index_fn(i)], n)
index_fn(i) = @match i begin
    ::Function              => i
    nothing                 => (x -> Int[])
    i::Int && if i > 0 end  => (x -> getindex(x, i))
    i::Int                  => (x -> getindex(x, lastindex(x) + i))
    UnitRange(s, e)         => (x -> getindex(x, UnitRange(convert_negidx(x, s), convert_negidx(x, e))))
    StepRange(s, p, e)      => (x -> getindex(x, StepRange(convert_negidx(x, s), p, convert_negidx(x, e))))
    ::Vector{Int}           => (x -> getindex(x, i))
end
"""
    select_clusters!(project::Project; by = identity, new = false, kwargs...)

Select clusters and store the results in `project.appendix[:clusters_possible]`.

# Keyword Arguments
* `by`: Either
    1. `fn::Function`: takes 'project.appendix[:clusters_possible][key]::Vector{Vector{Int}}`, and return mutiple clusters(`Vector{Vector{Int}}`) or a cluster(`Vector{Int}`).
    2. `nothing`: delete all clusters.
    3. `i::Int`: `i`th cluster. If `i` ≤ 0, it becomes `i+1`th cluster counting from the end.
    4. `r::Union{UnitRange, StepRange, Vector{Int}}`: multiple clusters.
* `new`: whether clear all existing clusters in `project.appendix[:clusters_candidate]` or not.
* Any keys in `project.appendix[:clusters_possible]`; the value can be set as descibed in `by`.
"""
function select_clusters!(project::Project; by = identity, new = false, kwargs...)
    clusters_possible = project.appendix[:clusters_possible]
    targets = collect(keys(clusters_possible))
    fn = convert(Vector{Any}, index_fn_v(by, length(targets)))
    for (key, val) in kwargs
        id = findfirst(==(eval(key)()), targets)
        isnothing(id) && continue
        fn[id] = index_fn(val)
    end
    clusters_candidate = get!(project.appendix, :clusters_candidate, Dictionary{ClassSP, Vector{Int}}())
    new && empty!(clusters_candidate)
    for (key, val) in zip(targets, fn)
        union!(get!(clusters_candidate, key, Int[]), val(clusters_possible[key])...)
    end
    filter!(!isempty, clusters_candidate)
    printstyled("Selected candidate clusters:\n", color = :green)
    display(clusters_candidate)
    println()
    project
end

find_tbl(expr) = @match expr begin
    Expr(:kw, args...)  => false
    Expr(head, args...) => any(find_tbl, args)
    :tbl                => true
    _                   => false
end

find_terms(expr) = @match expr begin
    Expr(:macrocall, arg, args...) && if arg == Symbol("@formula") end => unique(reduce_terms.(union(getterms(eval(expr))...)))
    Expr(:(.), :tbl, QuoteNode(x))                  => [x]
    Expr(head, args...)                             => begin
        terms = find_terms.(args)
        terms = union(terms...)
        isempty(terms) ? [] : terms
    end
    x                                               => []
end

reduce_terms(expr) = @match expr begin
    Expr(:call, :(^), arg, arg2)    => arg
    x                               => x
end

replace_formula(expr, fnterms, nmterms) = @match expr begin
    Expr(:macrocall, arg, args...) && if arg == Symbol("@formula") end => replace_terms(expr, fnterms, nmterms)
    Expr(head, args...) => Expr(head, map(x -> replace_formula(x, fnterms, nmterms), args)...)
    x                   => x
end
replace_terms(expr, fnterms, nmterms) = @match expr begin
    Expr(:call, :(∘), args...)  => nmterms[findfirst(==(expr), fnterms)]
    Expr(:(->), args...)        => nmterms[findfirst(==(expr), fnterms)]
    Expr(head, args...)         => Expr(head, map(x -> replace_terms(x, fnterms, nmterms), args)...)
    if expr in fnterms end      => nmterms[findfirst(==(expr), fnterms)]
    x                           => x
end

fnterm2symbol(e) = occursin("∘", string(e)) ? Symbol(replace(string(e), " ∘ " => "_")) : e

model_fn() = Expr(:call, :CurrentModelCall)
function model_fn(expr)
    expr == Expr(:tuple) && return Expr(:call, :CurrentModelCall)
    fnterms = find_terms(expr)
    cl = :cluster in fnterms
    filter!(!=(:cluster), fnterms)
    nmterms = fnterm2symbol.(fnterms)
    expr = find_tbl(expr) ? expr : Expr(expr.head, push!(expr.args, :tbl)...)
    expr = replace_formula(expr, fnterms, nmterms)
    Expr(:call, Expr(:curly, :ModelCall, cl), QuoteNode(nmterms), QuoteNode(fnterms), Expr(:(->), :tbl, Expr(:block, LineNumberNode(@__LINE__, @__FILE__), expr)))
end
"""
    @model()
    @model(expr)
    @model(expr...)

Create function(s) taking data as input and return a regression model.
    
# Arguuments
* `()` or no input represent current model, i.e., `preject.appendix[:model_formula]`.
* An expression as if fitting a regression model. `tbl` represents the data which can be omited or used in other arguments, e.g., `lm(@fomula(y~x); wts=tbl.wts)`.
"""
macro model()
    return quote $(model_fn()) end
end

macro model(expr)
    return quote $(model_fn(expr)) end
end

macro model(expr...)
    expr = model_fn.(expr)
    return quote map(eval, $expr) end
end
"""
    model_clusters!(project::Project)
    model_clusters!(model::RetentionModelCall, project::Project)

Fit `model` with `project.appendix[:clusters_candidate]`. The result is stored in `project.appendix[:clusters_model]`.
"""
function model_clusters!(modelcall::RetentionModelCall, project::Project)
    haskey(project.appendix, :clusters_model) && delete!(project.appendix, :clusters_model)
    @p pairs(project.appendix[:clusters_candidate]) |>
        map(Table(
                mw = map(mw, @views project.analytes[_[2]]),
                rt = map(rt, @views project.analytes[_[2]]),
                cluster = repeat([_[1]], length(_[2]))
            )) |>
        reduce(vcat) |>
        modelcall |> 
        RetentionModel(modelcall) |> 
        insert!(project.appendix, :clusters_model)
    display(project.appendix[:clusters_model].model)
    println()
    project
end
#model_clusters!(project::Project) = model_clusters!(ModelCall{true}([:rt, :mw, :cluster], [:rt, :mw, :cluster], tbl -> lm(@formula(rt ~ mw + cluster), tbl)), project)
"""
    model_rt!(modelcall::RetentionModelCall, aquery::AbstractQuery)
    model_rt!(modelcall::RetentionModelCall, project::Project; analytes = project.analytes)

Fit `model` with `analytes`. The result is stored in `project.appendix[:rt_model]`.
"""
model_rt!(modelcall::RetentionModelCall, aquery::AbstractQuery) = (model_rt!(modelcall, aquery.project; analytes = aquery.result); aquery)
function model_rt!(modelcall::RetentionModelCall, project::Project; analytes = project.analytes)
    haskey(project.appendix, :rt_model) && delete!(project.appendix, :rt_model)
    @p Table(; (modelcall.nmterms .=> map(x -> eval(x).(analytes), modelcall.fnterms))...) |>
        modelcall |> 
        RetentionModel(modelcall) |> 
        insert!(project.appendix, :rt_model)
    display(project.appendix[:rt_model].model)
    println()
    project
end

"""
    compare_models(project::Project, models::Tuple; analytes = project.analytes)
    compare_models(aquery::AbstractQuery, models::Tuple)
    compare_models(project::Project, models...; analytes = project.analytes)
    compare_models(aquery::AbstractQuery, models...)

Compare models by ANOVA.
"""
compare_models(aquery::AbstractQuery, models::Tuple) = compare_models(aquery.project, models...; analytes = aquery.result)
compare_models(aquery::AbstractQuery, models...) = compare_models(aquery.project, models...; analytes = aquery.result)
compare_models(project::Project, models::Tuple; analytes = project.analytes) = compare_models(project, models...; analytes = analytes)
function compare_models(project::Project, models...; analytes = project.analytes)
    if any(model -> isa(model, ModelCall{true}), models)
        models = replace(collect(models), CurrentModelCall() => project.appendix[:clusters_model].fn)
        tbl = @p pairs(project.appendix[:clusters_candidate]) |>
            map(Table(
                    mw = map(mw, @views project.analytes[_[2]]),
                    rt = map(rt, @views project.analytes[_[2]]),
                    cluster = repeat([_[1]], length(_[2]))
                )) |>
            reduce(vcat)
    else
        models = replace(collect(models), CurrentModelCall() => project.appendix[:rt_model].fn)
        nmterms = mapreduce(x -> x.nmterms, union, models)
        fnterms = mapreduce(x -> x.fnterms, union, models)
        tbl = Table(; (nmterms .=> map(x -> eval(x).(analytes), fnterms))...)
    end
    anova(map(f -> f(tbl), models)...)
end
"""
    predfn_clusters!(project::Project; replaces...)

Create predition function for clusters.

# Keyword argument
Any subtypes of `ClassSP`. The key represents class that is going to be assigned to other class. The value is the assigned class.
"""
function predfn_clusters!(project::Project; replaces...)
    replaces = @p pairs(Dictionary(replaces)) map(eval(first(_))() => last(_)()) filter(in(last(_), project.appendix[:clusters_model].model.mf.schema.schema[Term(:cluster)].contrasts.levels))
    function fn(model, analyte::AnalyteSP, rt_tol)
        cluster = replace!(ClassSP[deisomerized(class(analyte))], replaces...)[1]
        in(cluster, model.model.mf.schema.schema[Term(:cluster)].contrasts.levels) ?
            abs(rt(analyte) - predict(model.model, [(mw = mw(analyte), cluster = cluster)])[1]) <= rt_tol : false
    end
    haskey(project.appendix, :clusters_predict) && delete!(project.appendix, :clusters_predict)
    insert!(project.appendix, :clusters_predict, fn)
    project
end
"""
    expand_clusters!(project::Project; rt_tol = 1)

Expand each cluster by `project.appendix[:clusters_predict]`.
"""
function expand_clusters!(project::Project; rt_tol = 1)
    pred = project.appendix[:clusters_predict]
    for (i, analyte) in enumerate(project)
        pred(project.appendix[:clusters_model], analyte, rt_tol) &&
            push!(get!(project.appendix[:clusters_candidate], deisomerized(class(analyte)), Int[]), i)
    end
    foreach(unique!, project.appendix[:clusters_candidate])
    printstyled("Expanded candidate clusters:\n", color = :green)
    display(project.appendix[:clusters_candidate])
    println()
    project
end
"""
    show_clusters(project::Project)

Display clusters.
"""
function show_clusters(project::Project)
    printstyled("Current clusters:\n", color = :green)
    display(get(project.appendix, :clusters, Dictionary{ClassSP, Vector{Int}}()))
    println()
    printstyled("Possible clusters:\n", color = :green)
    display(get(project.appendix, :clusters_possible, Dictionary{ClassSP, Vector{Vector{Int}}}()))
    println()
    printstyled("Candidate clusters:\n", color = :green)
    display(get(project.appendix, :clusters_candidate, Dictionary{ClassSP, Vector{Int}}()))
    return
end
"""
    replace_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate])

Replace old clusters with `new_clusters`.
"""
replace_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate]) = (empty!(get!(project.appendix, :clusters, Dictionary{ClassSP, Vector{Int}}())); update_clusters!(project, new_clusters))
"""
    update_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate])

Update `project.appendix[:clusters]` with `new_clusters`.
"""
function update_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate])
    clusters = get!(project.appendix, :clusters, Dictionary{ClassSP, Vector{Int}}())
    @p pairs(new_clusters) foreach(union!(get!(clusters, _[1], Int[]), _[2]))
    printstyled("Updated clusters:\n", color = :green)
    display(clusters)
    println()
    project
end

apply_clusters!(aquery::AbstractQuery; kwargs...) = (apply_clusters!(aquery.project; analytes = aquery.result, kwargs...); aquery)
function apply_clusters!(project::Project; analytes = project.analytes)
    for (i, analyte) in enumerate(analytes)
        cls = deisomerized(class(analyte))
        analyte.states[states_id(:rt)] = in(cls, keys(project.appendix[:clusters])) ? (in((first ∘ parentindices)(analytes)[i], project.appendix[:clusters][cls]) ? 1 : -1) : 0
    end
    project
end
"""
    apply_rt!(aquery::AbstractQuery; kwargs...)
    apply_rt!(project::Project; analytes = project.analytes, atol = 0.5, rtol = 0.05)

Apply rt prediction to each analytes. `analyte.states[states_id(:rt)]` is set as `1` if the prediction error is under certain value; `-1` if the prediction error is too large; `0` if there is no prediction can be made. Prediction is based on `project.appendix[:rt_model]`.
"""
apply_rt!(aquery::AbstractQuery; kwargs...) = (apply_rt!(aquery.project; analytes = aquery.result, kwargs...); aquery)
function apply_rt!(project::Project; analytes = project.analytes, atol = 0.5, rtol = 0.05)
    tbl = Table(; (project.appendix[:rt_model].fn.nmterms .=> map(x -> eval(x).(analytes), project.appendix[:rt_model].fn.fnterms))...)
    include_id = find_predictable(project.appendix[:rt_model].model, tbl)
    ŷ = predict(project.appendix[:rt_model].model, tbl[include_id])
    y = @view getproperty(tbl, first(propertynames(tbl)))[include_id]
    for (v, v̂, analyte) in zip(y, ŷ, view(analytes, include_id))
        analyte.states[states_id(:rt)] = (abs(v - v̂) < atol || abs(v - v̂) / v < rtol) ? 1 : -1
    end
    project
end

function find_predictable(model, tbl)
    sch = model.mf.schema
    toval = [name for name in propertynames(tbl) if isa(get(sch, Term(name), nothing), CategoricalTerm)]
    todel = Int[]
    for name in toval
        level = sch[Term(name)].contrasts.levels
        union!(todel, findall(!in(level), getproperty(tbl, name)))
    end
    setdiff(eachindex(tbl), todel)
end
"""
    predict_rt(analyte::AnalyteSP)
    predict_rt(analytes::AbstractVector{AnalyteSP})

Predict rt of analyte(s) based on `project.appendix[:rt_model]`. If prediction is not available, it returns `nothing`.
"""
function predict_rt(analyte::AnalyteSP)
    project = last(analyte).project
    tbl = Table(; (project.appendix[:rt_model].fn.nmterms .=> map(x -> eval(x).([analyte]), project.appendix[:rt_model].fn.fnterms))...)
    sch = project.appendix[:rt_model].model.mf.schema
    toval = [name for name in propertynames(tbl) if isa(get(sch, Term(name), nothing), CategoricalTerm)]
    for name in toval
        level = sch[Term(name)].contrasts.levels
        any(!in(level), getproperty(tbl, name)) && return nothing
    end
    predict(project.appendix[:rt_model].model, tbl)
end
function predict_rt(analytes::AbstractVector{AnalyteSP})
    project = last(first(analytes)).project
    tbl = Table(; (project.appendix[:rt_model].fn.nmterms .=> map(x -> eval(x).(analytes), project.appendix[:rt_model].fn.fnterms))...)
    include_id = find_predictable(project.appendix[:rt_model].model, tbl)
    ŷ = Vector{Union{Float64, Nothing}}(undef, length(tbl))
    fill!(ŷ, nothing)
    ŷ[include_id] .= predict(project.appendix[:rt_model].model, tbl[include_id])
end
"""
    err_rt(analyte::AnalyteSP)
    err_rt(analytes::AbstractVector{AnalyteSP})

Compute the error of rt predition based on `project.appendix[:rt_model]`. If prediction is not available, it returns `nothing`.
"""
function err_rt(analyte::AnalyteSP)
    prt = predict_rt(analyte)
    isnothing(prt) ? prt : prt[1] - analyte.rt
end
function err_rt(analytes::AbstractVector{AnalyteSP})
    ŷs = predict_rt(analytes::AbstractVector{AnalyteSP})
    map(analytes, ŷs) do analyte, ŷ
        isnothing(ŷ) ? nothing : ŷ - analyte.rt
    end
end
"""
    abs_err_rt(analyte::AnalyteSP)
    abs_err_rt(analytes::AbstractVector{AnalyteSP})

Compute the absolute error of rt predition based on `project.appendix[:rt_model]`. If prediction is not available, it returns `nothing`.
"""
function abs_err_rt(analyte::AnalyteSP)
    v = err_rt(analyte)
    isnothing(v) ? v : abs(v)
end
function abs_err_rt(analytes::AbstractVector{AnalyteSP})
    vs = err_rt(analytes)
    map(vs) do v
        isnothing(v) ? v : abs(v)
    end
end