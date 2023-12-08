"""
    plot_rt_mw(aquery::AbstractQuery; kwargs...)
    plot_rt_mw(project::Project;
                analytes = project.analytes,
                all = true,
                clusters = nothing,
                groupby = deisomerized ∘ class,
                clusters_model = false,
                rt_prediction = false,
                xlabel = "Retention Time (min)",
                ylabel = "Molecular weight",
                title = "Analytes",
                kwargs...)

Plot retention time vs molecular weight on specified analytes.

# Keyword Arguments
* `analytes`: `Vector{AnalyteSP}`, analytes to be plotted.
* `all`: determine whether include all other analytes or not.
* `clusters`: determine which kind of clusters to be used.
    * `:possible`: create a plot for all clusters in `project.appendix[:clusters_possible][cls]` for each `cls` occured in `analytes`.
    * `:candidate`: take the union of `project.appendix[:clusters_candidate]` and `analytes`.
    * `:clusters`: take the union of `project.appendix[:clusters]` and `analytes`.
* `groupby`: function for grouping analytes. It applies to an analyte and return an unique property.
* `clusters_model`: determine whether plot prediction line for clusters.
* `rt_prediction`: determine whether plot rt prediction results.
* Any other keyword arguments for plot.
"""
plot_rt_mw(aquery::AbstractQuery; kwargs...) = plot_rt_mw(aquery.project; analytes = aquery.result, kwargs...)
function plot_rt_mw(project::Project;
                    analytes = project.analytes,
                    clusters = nothing,
                    groupby = deisomerized ∘ class,
                    clusters_model = false,
                    rt_prediction = false,
                    all = !rt_prediction,
                    xlabel = "Retention Time (min)",
                    ylabel = "Molecular weight",
                    title = "Analytes",
                    kwargs...)
    #≡(clusters, :possible) && return map(collect(keys(project.appendix[:clusters_possible]))) do cls
    #    plot_rt_mw(project, cls; all, xlabel, ylabel, deepcopy(kwargs)...)
    #end
    if ≡(clusters, :possible)
        clusters = project.appendix[:clusters_possible]
        kys = similar(analytes, String)
        Threads.@threads for id in eachindex(analytes)
            oid  = (first ∘ parentindices)(analytes)[id]
            ky = "Nothing"
            for (cls, cluster) in pairs(clusters)
                nid = findfirst(x -> in(oid, x), cluster)
                isnothing(nid) && continue
                ky = string(cls, "_", nid)
                break
            end
            kys[id] = ky
        end
        gid = group(kys, eachindex(analytes))
        gana = dictionary(ky => @view analytes[gid[ky]] for ky in unique!(kys))
    else
        gana = isnothing(groupby) ? Dictionary([:Analytes], [analytes]) : groupview(groupby, analytes)
        if !isnothing(clusters)
            used_clusters = @match clusters begin
                :clusters  => project.appendix[:clusters]
                :candidate => project.appendix[:clusters_candidate]
            end
            gana = @p pairs(gana) |>
                        map(filter(x -> in(x, get(used_clusters, _[1], Int[])), (first ∘ parentindices)(_[2]))) |>
                        filter(!isempty(_)) |>
                        map(@view project.analytes[_])
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
        id = setdiff(eachindex(project.analytes), map(first ∘ parentindices, gana)...)
        others = project.analytes[id]
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
    (clusters_model && haskey(project.appendix, :clusters_model)) || return scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    lim_mass = all ? extrema(mw.(project.analytes)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(mass))))
    lim_rt = all ? extrema(rt.(project.analytes)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(ret))))
    cls = project.appendix[:clusters_model].model.mf.data.cluster |> unique
    mwrange = range(lim_mass..., length = 100)
    tbl = Table(mw = repeat(mwrange, length(cls)), cluster = repeat(cls, inner = length(mwrange)))
    tbl = Table(tbl, rt = predict(project.appendix[:clusters_model].model, tbl))
    @p tbl filter!(between(_.rt, lim_rt))
    gtbl = groupview(getproperty(:cluster), tbl)
    for tbl in gtbl
        plot!(tbl.rt, tbl.mw; label = "Model_$(first(tbl.cluster))", get_attributes!(kwargs)...)
    end
    scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    p
end

function plot_rt_mw(project::Project, cls; all = true, xlabel = "Retention Time (min)", ylabel = "Molecular weight", title = cls, kwargs...)
    clusters = project.appendix[:clusters_possible][isa(cls, Type) ? cls() : cls]
    mass = @p clusters map(map(mw, @views project.analytes[_]))
    ret = @p clusters map(map(rt, @views project.analytes[_]))
    p = scatter()
    if all
        id = setdiff(eachindex(project.analytes), clusters...)
        others = project.analytes[id]
        isempty(others) || scatter!((@p others map((rt = rt(_), mw = mw(_))));
                    label = "Others", alpha = 0.2, hover = hover_repr(project, id), get_attributes!(kwargs)...)
    end
    for key in keys(clusters)
        scatter!(ret[key], mass[key]; label = "$key", hover = hover_repr(project, clusters[key]), get_attributes!(kwargs)...)
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

function hover_repr(analytes)
    ind = string.((first ∘ parentindices)(analytes), ". ")
    ss = @. split(repr(analytes), " st")
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
function hover_repr(analytes, rts)
    ind = string.((first ∘ parentindices)(analytes), ". ")
    ss = @. split(repr(analytes), " st")
    Δrts = @. rts - getproperty(analytes, :rt)
    map(ind, ss, Δrts) do id, s, Δrt
        string(id, s[1], " Δrt=", round(Δrt; digits = 2))
    end
end