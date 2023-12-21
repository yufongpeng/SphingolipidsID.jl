"""
    plot_rt_mw(aquery::AbstractQuery; kwargs...)
    plot_rt_mw(project::Project;
                analyte = project.analyte,
                all = true,
                cluster = nothing,
                groupby = deisomerized ∘ class,
                cluster_model = false,
                rt_prediction = false,
                xlabel = "Retention Time (min)",
                ylabel = "Molecular weight",
                title = "Analytes",
                kwargs...)

Plot retention time vs molecular weight on specified analyte.

# Keyword Arguments
* `analyte`: `Vector{AnalyteSP}`, analytes to be plotted.
* `all`: determine whether include all other analytes or not.
* `cluster`: determine which kind of clusters to be used.
    * `:possible`: create a plot for all clusters in `project.appendix[:cluster_possible][cls]` for each `cls` occured in `analyte`.
    * `:candidate`: take the union of `project.appendix[:cluster_candidate]` and `analyte`.
    * `:cluster`: take the union of `project.appendix[:cluster]` and `analyte`.
* `groupby`: function for grouping analytes. It applies to an analyte and return an unique property.
* `cluster_model`: determine whether plot prediction line for clusters.
* `rt_prediction`: determine whether plot rt prediction results.
* Any other keyword arguments for plot.
"""
plot_rt_mw(aquery::AbstractQuery; kwargs...) = plot_rt_mw(aquery.project; analyte = aquery.result, kwargs...)
function plot_rt_mw(project::Project;
                    analyte = project.analyte,
                    cluster = nothing,
                    groupby = deisomerized ∘ class,
                    cluster_model = false,
                    rt_prediction = false,
                    all = !rt_prediction,
                    xlabel = "Retention Time (min)",
                    ylabel = "Molecular weight",
                    title = "Analytes",
                    kwargs...)
    #≡(cluster, :possible) && return map(collect(keys(project.appendix[:cluster_possible]))) do cls
    #    plot_rt_mw(project, cls; all, xlabel, ylabel, deepcopy(kwargs)...)
    #end
    if ≡(cluster, :possible)
        cluster = project.appendix[:cluster_possible]
        kys = similar(analyte, String)
        Threads.@threads for id in eachindex(analyte)
            oid  = (first ∘ parentindices)(analyte)[id]
            ky = "Nothing"
            for (cls, clt) in pairs(cluster)
                nid = findfirst(x -> in(oid, x), clt)
                isnothing(nid) && continue
                ky = string(cls, "_", nid)
                break
            end
            kys[id] = ky
        end
        gid = group(kys, eachindex(analyte))
        gana = dictionary(ky => @view analyte[gid[ky]] for ky in unique!(kys))
    else
        gana = isnothing(groupby) ? Dictionary([:Analytes], [analyte]) : groupview(groupby, analyte)
        if !isnothing(cluster)
            used_cluster = @match cluster begin
                :cluster  => project.appendix[:cluster]
                :candidate => project.appendix[:cluster_candidate]
            end
            gana = @p pairs(gana) |>
                        map(filter(x -> in(x, get(used_cluster, _[1], Int[])), (first ∘ parentindices)(_[2]))) |>
                        filter(!isempty(_)) |>
                        map(@view project.analyte[_])
        end
    end
    data = Dictionary{keytype(gana), Table}()
    if rt_prediction
        gana_new = Dictionary{keytype(gana), eltype(gana)}()
        for (key, analyte) in pairs(gana)
            tbl = Table(; (union(project.appendix[:rt_model].fn.nmterms, [:mw]) .=> map(x -> eval(x).(analyte), union(project.appendix[:rt_model].fn.fnterms, [:mw])))...)
            include_id = find_predictable(project.appendix[:rt_model].model, tbl)
            if !isempty(include_id)
                insert!(data, key, Table(tbl[include_id]; predicted_rt = predict(project.appendix[:rt_model].model, tbl[include_id])))
                insert!(gana_new, key, @view analyte[include_id])
            end
        end
        gana = gana_new
    else
        for (key, analyte) in pairs(gana)
            insert!(data, key, Table(; rt = rt.(analyte), mw = mw.(analyte)))
        end
    end
    #mass = @p gana map(map(mw, _))
    #ret = @p gana map(map(rt, _))
    p = scatter()
    if all
        id = setdiff(eachindex(project.analyte), map(first ∘ parentindices, gana)...)
        others = project.analyte[id]
        isempty(others) || scatter!((@p others map((rt = rt(_), mw = mw(_))));
                    label = "Others", alpha = 0.2, hover = hover_repr(project, id), get_attributes!(kwargs)...)
    end
    if rt_prediction
        for (key, analyte) in pairs(gana)
            for row in data[key]
                plot!([row.rt, row.predicted_rt], [row.mw, row.mw]; alpha = 0.5, linestyle = :dot, color = :grey, label = false, get_attributes!(kwargs)...)
            end
            scatter!(data[key].predicted_rt, data[key].mw; alpha = 0.5, label = false, marker = :utriangle, hover = hover_repr(analyte, data[key].predicted_rt), get_attributes!(kwargs)...)
        end
    end
    for (key, analyte) in pairs(gana)
        scatter!(data[key].rt, data[key].mw; label = "$key", hover = hover_repr(analyte), get_attributes!(kwargs)...)
    end
    (cluster_model && haskey(project.appendix, :cluster_model)) || return scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    lim_mass = all ? extrema(mw.(project.analyte)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(mass))))
    lim_rt = all ? extrema(rt.(project.analyte)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(ret))))
    cls = project.appendix[:cluster_model].model.mf.data.cluster |> unique
    mwrange = range(lim_mass..., length = 100)
    tbl = Table(mw = repeat(mwrange, length(cls)), cluster = repeat(cls, inner = length(mwrange)))
    tbl = Table(tbl, rt = predict(project.appendix[:cluster_model].model, tbl))
    @p tbl filter!(between(_.rt, lim_rt))
    gtbl = groupview(getproperty(:cluster), tbl)
    for tbl in gtbl
        plot!(tbl.rt, tbl.mw; label = "Model_$(first(tbl.cluster))", get_attributes!(kwargs)...)
    end
    scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    p
end

function plot_rt_mw(project::Project, cls; all = true, xlabel = "Retention Time (min)", ylabel = "Molecular weight", title = cls, kwargs...)
    cluster_ = project.appendix[:cluster_possible][isa(cls, Type) ? cls() : cls]
    mass = @p cluster_ map(map(mw, @views project.analyte[_]))
    ret = @p cluster_ map(map(rt, @views project.analyte[_]))
    p = scatter()
    if all
        id = setdiff(eachindex(project.analyte), cluster_...)
        others = project.analyte[id]
        isempty(others) || scatter!((@p others map((rt = rt(_), mw = mw(_))));
                    label = "Others", alpha = 0.2, hover = hover_repr(project, id), get_attributes!(kwargs)...)
    end
    for key in keys(cluster_)
        scatter!(ret[key], mass[key]; label = "$key", hover = hover_repr(project, cluster_[key]), get_attributes!(kwargs)...)
    end
    scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    p
end

get_attributes!(x::Base.Pairs) = @p keys(x) map(_ => get_attributes!(getproperty(values(x), _))) filter(!isnothing(_[2]))
get_attributes!(v::Vector) = isempty(v) ? nothing : popfirst!(v)
get_attributes!(v) = v

"""
histogram_transition(tbl::Table; 
                        end_time = nothing,                     
                        xlabel = "Retention Time (min)",
                        ylabel = "Counts",
                        title = "Transitions", 
                        color = :steelblue, 
                        label = nothing,
                        linecolor = :steelblue,
                        kwargs...
                    )

Plot histogram of transitions. The input `tbl` can be a transition table or table from `concurrent_transition`.

# Keyword Arguments
* `end_time` is the end of aquisition.
* Any other keyword arguments for plot.
"""
function histogram_transition(tbl::Table; 
                                end_time = nothing,                     
                                xlabel = "Retention Time (min)",
                                ylabel = "Counts",
                                title = "Transitions", 
                                color = :steelblue, 
                                label = nothing,
                                linecolor = :steelblue,
                                kwargs...
                            )
    if !in(:count, propertynames(tbl)) 
        tbl = isnothing(end_time) ? concurrent_transition(tbl) : concurrent_transition(tbl; end_time)
    end
    x = (tbl.rt[1:end - 1] .+ tbl.rt[2:end]) ./ 2
    w = tbl.rt[2:end] .- tbl.rt[1:end - 1]
    bar(x, tbl.count[1:end - 1]; xlabel, ylabel, title, color, label, linecolor, bar_width = w, kwargs...)
end

function plot_rt_rt()
end

function hover_repr(analyte)
    ind = string.((first ∘ parentindices)(analyte), ". ")
    ss = @. split(repr(analyte), " st")
    map(ind, ss) do id, s
        st = string("st", s[2])
        string(id, s[1], "<br>", "\n" ^ (2 * length(id) - 1), rpad(st, 2 * length(s[1]) - length(st) - 12))
    end
end
hover_repr(project::Project, id) = map(id) do i
    ind = "$i. "
    s = split(repr(project[i]), " st")
    st = string("st", s[2])
    string(ind, s[1], "<br>", "\n" ^ (2 * length(ind) - 1), rpad(st, 2 * length(s[1]) - length(st) - 12))
end
function hover_repr(analyte, rts)
    ind = string.((first ∘ parentindices)(analyte), ". ")
    ss = @. split(repr(analyte), " st")
    Δrts = @. rts - getproperty(analyte, :rt)
    map(ind, ss, Δrts) do id, s, Δrt
        string(id, s[1], " Δrt=", round(Δrt; digits = 2))
    end
end