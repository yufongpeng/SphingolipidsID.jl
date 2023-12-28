"""
    group_analyte(analyte)

Group `analyte` based on `class`.
"""
group_analyte(analyte) = @p analyte groupview(deisomerized(class(_))) map((first ∘ parentindices)(_))
"""
    initialize_cluster!(aquery::AbstractQuery; kwargs...)
    initialize_cluster!(project::Project; analyte = project.analyte)

Initialize `project.appendix[:cluster]` with `group_analyte(analyte)`.
"""
initialize_cluster!(aquery::AbstractQuery; kwargs...) = (initialize_cluster!(aquery.project; analyte = aquery.result, kwargs...); aquery)
initialize_cluster!(project::Project; analyte = project.analyte) = replace_cluster!(project, group_analyte(analyte))
"""
    analyte2cluster!(aquery::AbstractQuery; kwargs...)
    analyte2cluster!(project::Project; new = false, scale = 0.0, radius = 0.0, analyte = project.analyte, kwargs...)

Run dbscan algorithm on `analyte`, and put all results in `project.appendix[:cluster_possible]`.

# Keyword Arguments
* `new`: whether clear all existing cluster in `project.appendix[:cluster_possible]` or not.
* `scale`: the normalizing factor for mw.
* `radius`: radius of dbscan.

See `dbscan` for more detailed settings.
"""
analyte2cluster!(aquery::AbstractQuery; kwargs...) = (analyte2cluster!(aquery.project; analyte = aquery.result, kwargs...); aquery)
function analyte2cluster!(project::Project; new = false, scale = 0.0, radius = 0.0, analyte = project.analyte, kwargs...)
    valid_id = collect((first ∘ parentindices)(analyte))
    groups = @p project.appendix[:cluster] map(filter(x -> in(x, valid_id), _)) filter(length(_) > 2)
    ret = @p groups map(map(rt, @views project.analyte[_]))
    maxrt = maximum(maximum(r) for r in ret)
    mass = @p groups map(map(mw, @views project.analyte[_]))
    scale = scale > 0 ? scale : median(map(m -> quantile(m, 0.9) - quantile(m, 0.1), mass)) * 2 / maxrt
    mass = @p mass map(_ / scale)
    radius = radius > 0 ? radius : maxrt / 20 * sqrt(2)
    cluster_ = @p zip(ret, mass) map(dbscan(hcat(_...)', radius; kwargs...)) map(map(getproperty(:core_indices), _.clusters))
    cluster_ = @p zip(cluster_, groups) map(map(x -> _[2][x], _[1]))
    cluster_possible = get!(project.appendix, :cluster_possible, Dictionary{ClassSP, Vector{Vector{Int}}}())
    new && empty!(cluster_possible)
    @p zip(keys(groups), cluster_) foreach(union!(empty!(get!(cluster_possible, _[1], Vector{Int}[])), _[2]))
    printstyled("Possible clusters:\n", color = :green)
    display(project.appendix[:cluster_possible])
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
    select_cluster!(project::Project; by = identity, new = false, kwargs...)

Select clusters and store the results in `project.appendix[:cluster_possible]`.

# Keyword Arguments
* `by`: Either
    1. `fn::Function`: takes 'project.appendix[:cluster_possible][key]::Vector{Vector{Int}}`, and return mutiple clusters(`Vector{Vector{Int}}`) or a cluster(`Vector{Int}`).
    2. `nothing`: delete all clusters.
    3. `i::Int`: `i`th cluster. If `i` ≤ 0, it becomes `i+1`th cluster counting from the end.
    4. `r::Union{UnitRange, StepRange, Vector{Int}}`: multiple clusters.
* `new`: whether clear all existing clusters in `project.appendix[:cluster_candidate]` or not.
* Any keys in `project.appendix[:cluster_possible]`; the value can be set as descibed in `by`.
"""
function select_cluster!(project::Project; by = identity, new = false, kwargs...)
    cluster_possible = project.appendix[:cluster_possible]
    targets = collect(keys(cluster_possible))
    fn = convert(Vector{Any}, index_fn_v(by, length(targets)))
    for (key, val) in kwargs
        id = findfirst(==(eval(key)()), targets)
        isnothing(id) && continue
        fn[id] = index_fn(val)
    end
    cluster_candidate = get!(project.appendix, :cluster_candidate, Dictionary{ClassSP, Vector{Int}}())
    new && empty!(cluster_candidate)
    for (key, val) in zip(targets, fn)
        union!(get!(cluster_candidate, key, Int[]), val(cluster_possible[key])...)
    end
    filter!(!isempty, cluster_candidate)
    printstyled("Selected candidate clusters:\n", color = :green)
    display(cluster_candidate)
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

Create function(s) taking data as input and returning a regression model.
    
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
    model_cluster!(project::Project)
    model_cluster!(model::RetentionModelCall, project::Project)

Fit `model` with `project.appendix[:cluster_candidate]`. The result is stored in `project.appendix[:cluster_model]`.
"""
function model_cluster!(modelcall::RetentionModelCall, project::Project)
    haskey(project.appendix, :cluster_model) && delete!(project.appendix, :cluster_model)
    @p pairs(project.appendix[:cluster_candidate]) |>
        map(Table(
                mw = map(mw, @views project.analyte[_[2]]),
                rt = map(rt, @views project.analyte[_[2]]),
                cluster = repeat([_[1]], length(_[2]))
            )) |>
        reduce(vcat) |>
        modelcall |> 
        RetentionModel(modelcall) |> 
        insert!(project.appendix, :cluster_model)
    display(project.appendix[:cluster_model].model)
    println()
    project
end
#model_cluster!(project::Project) = model_cluster!(ModelCall{true}([:rt, :mw, :cluster], [:rt, :mw, :cluster], tbl -> lm(@formula(rt ~ mw + cluster), tbl)), project)
"""
    model_rt!(modelcall::RetentionModelCall, aquery::AbstractQuery)
    model_rt!(modelcall::RetentionModelCall, project::Project; analyte = project.analyte)

Fit `model` with `analyte`. The result is stored in `project.appendix[:rt_model]`.
"""
model_rt!(modelcall::RetentionModelCall, aquery::AbstractQuery) = (model_rt!(modelcall, aquery.project; analyte = aquery.result); aquery)
function model_rt!(modelcall::RetentionModelCall, project::Project; analyte = project.analyte)
    haskey(project.appendix, :rt_model) && delete!(project.appendix, :rt_model)
    @p Table(; (modelcall.nmterms .=> map(x -> eval(x).(analyte), modelcall.fnterms))...) |>
        modelcall |> 
        RetentionModel(modelcall) |> 
        insert!(project.appendix, :rt_model)
    display(project.appendix[:rt_model].model)
    println()
    project
end

"""
    compare_models(project::Project, models::Tuple; analyte = project.analyte)
    compare_models(aquery::AbstractQuery, models::Tuple)
    compare_models(project::Project, models...; analyte = project.analyte)
    compare_models(aquery::AbstractQuery, models...)

Compare models by ANOVA.
"""
compare_models(aquery::AbstractQuery, models::Tuple) = compare_models(aquery.project, models...; analyte = aquery.result)
compare_models(aquery::AbstractQuery, models...) = compare_models(aquery.project, models...; analyte = aquery.result)
compare_models(project::Project, models::Tuple; analyte = project.analyte) = compare_models(project, models...; analyte)
function compare_models(project::Project, models...; analyte = project.analyte)
    if any(model -> isa(model, ModelCall{true}), models)
        models = replace(collect(models), CurrentModelCall() => project.appendix[:cluster_model].fn)
        tbl = @p pairs(project.appendix[:cluster_candidate]) |>
            map(Table(
                    mw = map(mw, @views project.analyte[_[2]]),
                    rt = map(rt, @views project.analyte[_[2]]),
                    cluster = repeat([_[1]], length(_[2]))
                )) |>
            reduce(vcat)
    else
        models = replace(collect(models), CurrentModelCall() => project.appendix[:rt_model].fn)
        nmterms = mapreduce(x -> x.nmterms, union, models)
        fnterms = mapreduce(x -> x.fnterms, union, models)
        tbl = Table(; (nmterms .=> map(x -> eval(x).(analyte), fnterms))...)
    end
    anova(map(f -> f(tbl), models)...)
end
"""
    predfn_cluster!(project::Project; replaces...)

Create predition function for clusters.

# Keyword argument
Any subtypes of `ClassSP`. The key represents class that is going to be assigned to other class. The value is the assigned class.
"""
function predfn_cluster!(project::Project; replaces...)
    replaces = @p pairs(Dictionary(replaces)) map(eval(first(_))() => last(_)()) filter(in(last(_), project.appendix[:cluster_model].model.mf.schema.schema[Term(:cluster)].contrasts.levels))
    function fn(model, analyte::AnalyteSP, rt_tol)
        cluster = replace!(ClassSP[deisomerized(class(analyte))], replaces...)[1]
        in(cluster, model.model.mf.schema.schema[Term(:cluster)].contrasts.levels) ?
            abs(rt(analyte) - predict(model.model, [(mw = mw(analyte), cluster = cluster)])[1]) <= rt_tol : false
    end
    haskey(project.appendix, :cluster_predict) && delete!(project.appendix, :cluster_predict)
    insert!(project.appendix, :cluster_predict, fn)
    project
end
"""
    expand_cluster!(project::Project; rt_tol = 1)

Expand each cluster by `project.appendix[:cluster_predict]`.
"""
function expand_cluster!(project::Project; rt_tol = 1)
    pred = project.appendix[:cluster_predict]
    for (i, analyte) in enumerate(project)
        pred(project.appendix[:cluster_model], analyte, rt_tol) &&
            push!(get!(project.appendix[:cluster_candidate], deisomerized(class(analyte)), Int[]), i)
    end
    foreach(unique!, project.appendix[:cluster_candidate])
    printstyled("Expanded candidate clusters:\n", color = :green)
    display(project.appendix[:cluster_candidate])
    println()
    project
end
"""
    show_cluster(project::Project)

Display clusters.
"""
function show_cluster(project::Project)
    printstyled("Current clusters:\n", color = :green)
    display(get(project.appendix, :cluster, Dictionary{ClassSP, Vector{Int}}()))
    println()
    printstyled("Possible clusters:\n", color = :green)
    display(get(project.appendix, :cluster_possible, Dictionary{ClassSP, Vector{Vector{Int}}}()))
    println()
    printstyled("Candidate clusters:\n", color = :green)
    display(get(project.appendix, :cluster_candidate, Dictionary{ClassSP, Vector{Int}}()))
    return
end
"""
    replace_cluster!(project::Project, new_cluster = project.appendix[:cluster_candidate])

Replace old clusters with `new_cluster`.
"""
replace_cluster!(project::Project, new_cluster = project.appendix[:cluster_candidate]) = (empty!(get!(project.appendix, :cluster, Dictionary{ClassSP, Vector{Int}}())); update_cluster!(project, new_cluster))
"""
    update_cluster!(project::Project, new_cluster = project.appendix[:cluster_candidate])

Update `project.appendix[:cluster]` with `new_cluster`.
"""
function update_cluster!(project::Project, new_cluster = project.appendix[:cluster_candidate])
    cluster_ = get!(project.appendix, :cluster, Dictionary{ClassSP, Vector{Int}}())
    @p pairs(new_cluster) foreach(union!(get!(cluster_, _[1], Int[]), _[2]))
    printstyled("Updated clusters:\n", color = :green)
    display(cluster_)
    println()
    project
end

apply_cluster!(aquery::AbstractQuery; kwargs...) = (apply_cluster!(aquery.project; analyte = aquery.result, kwargs...); aquery)
function apply_cluster!(project::Project; analyte = project.analyte)
    for (i, ana) in enumerate(analyte)
        cls = deisomerized(class(ana))
        ana.state[state_id(:rt)] = in(cls, keys(project.appendix[:cluster])) ? (in((first ∘ parentindices)(analyte)[i], project.appendix[:cluster][cls]) ? 1 : -1) : 0
    end
    project
end
"""
    apply_rt!(aquery::AbstractQuery; kwargs...)
    apply_rt!(project::Project; analyte = project.analyte, atol = 0.25, rtol = 0.05)

Apply rt prediction to each analyte. `analyte.state[state_id(:rt)]` is set as `1` if the prediction error is under certain value; `-1` if the prediction error is too large; `0` if there is no prediction can be made. 
Prediction is based on `project.appendix[:rt_model]`. `atol` and `rtol` can be a number or a function taking `analyte` as input and returning a number.
"""
apply_rt!(aquery::AbstractQuery; kwargs...) = (apply_rt!(aquery.project; analyte = aquery.result, kwargs...); aquery)
function apply_rt!(project::Project; analyte = project.analyte, atol = 0.5, rtol = 0.05)
    atol_fn = atol isa Number ? (x -> atol) : atol
    rtol_fn = rtol isa Number ? (x -> rtol) : rtol
    tbl = Table(; (project.appendix[:rt_model].fn.nmterms .=> map(x -> eval(x).(analyte), project.appendix[:rt_model].fn.fnterms))...)
    include_id = find_predictable(project.appendix[:rt_model].model, tbl)
    ŷ = predict(project.appendix[:rt_model].model, tbl[include_id])
    y = @view getproperty(tbl, first(propertynames(tbl)))[include_id]
    for (v, v̂, ana) in zip(y, ŷ, view(analyte, include_id))
        ana.state[state_id(:rt)] = (abs(v - v̂) < atol_fn(ana) || abs(v - v̂) / v < rtol_fn(ana)) ? 1 : -1
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
    for name in setdiff(propertynames(tbl), toval)
        union!(todel, findall(isnothing, getproperty(tbl, name)))
    end
    setdiff(eachindex(tbl), todel)
end
"""
    predict_rt(analyte::AnalyteSP)
    predict_rt(analyte::AbstractVector{AnalyteSP})

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
function predict_rt(analyte::AbstractVector{AnalyteSP})
    project = last(first(analyte)).project
    tbl = Table(; (project.appendix[:rt_model].fn.nmterms .=> map(x -> eval(x).(analyte), project.appendix[:rt_model].fn.fnterms))...)
    include_id = find_predictable(project.appendix[:rt_model].model, tbl)
    ŷ = Vector{Union{Float64, Nothing}}(undef, length(tbl))
    fill!(ŷ, nothing)
    ŷ[include_id] .= predict(project.appendix[:rt_model].model, tbl[include_id])
end
"""
    err_rt(analyte::AnalyteSP)
    err_rt(analyte::AbstractVector{AnalyteSP})

Compute the error of rt predition based on `project.appendix[:rt_model]`. If prediction is not available, it returns `nothing`.
"""
function err_rt(analyte::AnalyteSP)
    prt = predict_rt(analyte)
    isnothing(prt) ? prt : prt[1] - analyte.rt
end
function err_rt(analyte::AbstractVector{AnalyteSP})
    ŷs = predict_rt(analyte::AbstractVector{AnalyteSP})
    map(analyte, ŷs) do ana, ŷ
        isnothing(ŷ) ? nothing : ŷ - ana.rt
    end
end
"""
    abs_err_rt(analyte::AnalyteSP)
    abs_err_rt(analyte::AbstractVector{AnalyteSP})

Compute the absolute error of rt predition based on `project.appendix[:rt_model]`. If prediction is not available, it returns `nothing`.
"""
function abs_err_rt(analyte::AnalyteSP)
    v = err_rt(analyte)
    isnothing(v) ? v : abs(v)
end
function abs_err_rt(analyte::AbstractVector{AnalyteSP})
    vs = err_rt(analyte)
    map(vs) do v
        isnothing(v) ? v : abs(v)
    end
end

rt_correction(data::AbstractRawData, analytetable::Table; mz_tol = 0.35, rt_tol = 0.3, wfn = log2, wcol = :area) = rt_correction(analytetable, data; mz_tol, rt_tol, wfn, wcol)
function rt_correction(analytetable::Table, data::AbstractRawData; mz_tol = 0.35, rt_tol = 0.3, wfn = log2, wcol = :area)
    RTCorrection(data, analytetable, map(groupview(x -> class(x.analyte), analytetable)) do tbl
        rty = [Float64[] for i in eachindex(tbl)]
        wts = [Float64[] for i in eachindex(tbl)]
        rtx = tbl.rt
        for (i, t) in enumerate(tbl)
            for r in data.table
                (!between(t.rt, r.rt, rt_tol) || 
                !between(t.mz1, r.mz1, mz_tol) ||
                !between(t.mz2, data.mz2[r.mz2_id], mz_tol)) && continue
                push!(rty[i], r.rt)
                push!(wts[i], wfn(getproperty(r, wcol)))
            end
        end
        rtx = mapmany(zip(rtx, rty)) do (x, y)
            repeat([x], length(y))
        end
        rty = reduce(vcat, rty)
        wts = reduce(vcat, wts)
        sqrtw = diagm(sqrt.(max.(0, wts)))
        if length(unique(rtx)) > 1
            m = hcat(ones(Float64, length(rtx)), rtx)
            (sqrtw * m) \ (sqrtw * rty)
        elseif isempty(rtx)
            [0, 1]
        else
            [0, mean(rty, weights(wts)) / first(rtx)]
        end
    end)
end