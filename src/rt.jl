"""
    group_analyte(analytes)

Group `analytes` based on `class`.
"""
group_analyte(analytes) = @p analytes groupview(deisomerized(class(_))) map((first ∘ parentindices)(_))
"""
    initialize_clusters!(aquery::AbstractQuery; kwargs...)
    initialize_clusters!(project::Project; analytes = project.analytes)

Initialize `project.clusters` with `group_analyte(analytes)`.
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
    groups = @p project.clusters map(filter(x -> in(x, valid_id), _)) filter(length(_) > 2)
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

model_fn(expr) = expr == Expr(:tuple) ? Expr(:tuple) :
    Expr(:(->), :tbl, Expr(:block, LineNumberNode(@__LINE__, @__FILE__), find_tbl(expr) ? expr : Expr(expr.head, push!(expr.args, :tbl)...)))
"""
    @model(expr)
    @model(expr...)

Create function(s) taking data as input and return `GLM` model.
    
# Arguuments
* `_` or no input represent current model, i.e., `preject.appendix[:model_fn]`.
* An expression as if fitting a `GLM` model. `tbl` represents the data which can be omited or used in other arguments, e.g., `lm(@fomula(y~x); wts=tbl.wts)`.
"""
macro model(expr)
    return quote $(model_fn(expr)) end
end

macro model(expr...)
    expr = model_fn.(expr)
    return quote map(eval, $expr) end
end
"""
    model_clusters!(project::Project)
    model_clusters!(model, project::Project)

Apply `model` to the `project`. `model` must be a 1-arg function taking data as input and return `GLM` model.
"""
model_clusters!(project::Project) = model_clusters!(tbl -> lm(@formula(rt ~ mw + cluster), tbl), project)
function model_clusters!(model, project::Project)
    haskey(project.appendix, :clusters_model) && delete!(project.appendix, :clusters_model)
    haskey(project.appendix, :clusters_formula) && delete!(project.appendix, :clusters_formula)
    insert!(project.appendix, :clusters_formula, model)
    @p pairs(project.appendix[:clusters_candidate]) |>
        map(Table(
                mw = map(mw, @views project.analytes[_[2]]),
                rt = map(rt, @views project.analytes[_[2]]),
                cluster = repeat([_[1]], length(_[2]))
            )) |>
        reduce(vcat) |>
        model |>
        insert!(project.appendix, :clusters_model)
    display(project.appendix[:clusters_model])
    println()
    project
end
"""
    compare_models(project::Project, models::Tuple)
    compare_models(project::Project, models...)

Compare models by ANOVA.
"""
compare_models(project::Project, models::Tuple) = compare_models(project, models...)
function compare_models(project::Project, models...)
    models = replace(collect(models), () => project.appendix[:clusters_formula])
    tbl = @p pairs(project.appendix[:clusters_candidate]) |>
        map(Table(
                mw = map(mw, @views project.analytes[_[2]]),
                rt = map(rt, @views project.analytes[_[2]]),
                cluster = repeat([_[1]], length(_[2]))
            )) |>
        reduce(vcat)
    anova(map(f -> f(tbl), models)...)
end
"""
    predfn_clusters!(project::Project; replaces...)

Create predition function for clusters.

# Keyword argument
Any suntypes of `ClassSP`. The key represents class that is going to be assigned to other class. The value is the assigned class.
"""
function predfn_clusters!(project::Project; replaces...)
    replaces = @p pairs(Dictionary(replaces)) map(eval(first(_))() => last(_)()) filter(in(last(_), project.appendix[:clusters_model].mf.schema.schema[Term(:cluster)].contrasts.levels))
    function fn(model, analyte::AnalyteSP, rt_tol)
        cluster = replace!(ClassSP[deisomerized(class(analyte))], replaces...)[1]
        in(cluster, model.mf.schema.schema[Term(:cluster)].contrasts.levels) ?
            abs(rt(analyte) - predict(model, [(mw = mw(analyte), cluster = cluster)])[1]) <= rt_tol : false
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
    display(project.clusters)
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
replace_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate]) = (empty!(project.clusters); update_clusters!(project, new_clusters))
"""
    update_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate])

Update `project.clusters` with `new_clusters`.
"""
function update_clusters!(project::Project, new_clusters = project.appendix[:clusters_candidate])
    @p pairs(new_clusters) foreach(union!(get!(project.clusters, _[1], Int[]), _[2]))
    printstyled("Updated clusters:\n", color = :green)
    display(project.clusters)
    println()
    project
end
"""
    apply_clusters!(aquery::AbstractQuery; kwargs...)
    apply_clusters!(project::Project; analytes = project.analytes)

Apply the clusters to each analytes. `analyte.states[states_id(:rt)]` is set as `1` if the analyte is in the correct cluster; `-1` if the analyte is in a wrong cluster; `0` if there is no cluster for the class. 
"""
apply_clusters!(aquery::AbstractQuery; kwargs...) = (apply_clusters!(aquery.project; analytes = aquery.result, kwargs...); aquery)
function apply_clusters!(project::Project; analytes = project.analytes)
    for (i, analyte) in enumerate(analytes)
        cls = deisomerized(class(analyte))
        analyte.states[states_id(:rt)] = in(cls, keys(project.clusters)) ? (in((first ∘ parentindices)(analytes)[i], project.clusters[cls]) ? 1 : -1) : 0
    end
    project
end