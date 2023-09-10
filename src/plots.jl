"""
    plot_rt_mw(aquery::AbstractQuery; kwargs...)
    plot_rt_mw(project::Project;
                analytes = project.analytes,
                all = true,
                clusters = nothing,
                groupby = deisomerized ∘ class,
                model = false,
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
    * `:clusters`: take the union of `project.clusters` and `analytes`.
* `groupby`: function for grouping analytes. It applies to an analyte and return an unique property.
* `model`: determine whether plot prediction line or not.
* Any other keyword arguments for plot.
"""
plot_rt_mw(aquery::AbstractQuery; kwargs...) = plot_rt_mw(aquery.project; analytes = aquery.result, kwargs...)
function plot_rt_mw(project::Project;
                    analytes = project.analytes,
                    all = true,
                    clusters = nothing,
                    groupby = deisomerized ∘ class,
                    model = false,
                    xlabel = "Retention Time (min)",
                    ylabel = "Molecular weight",
                    title = "Analytes",
                    kwargs...)
    ≡(clusters, :possible) && return map(collect(keys(project.appendix[:clusters_possible]))) do cls
        plot_rt_mw(project, cls; all, xlabel, ylabel, deepcopy(kwargs)...)
    end
    gana = isnothing(groupby) ? Dictionary([:Analytes], [analytes]) : groupview(groupby, analytes)
    if !isnothing(clusters)
        used_clusters = @match clusters begin
            :clusters  => project.clusters
            :candidate => project.appendix[:clusters_candidate]
        end
        gana = @p pairs(gana) |>
                    map(filter(x -> in(x, get(used_clusters, _[1], Int[])), (first ∘ parentindices)(_[2]))) |>
                    filter(!isempty(_)) |>
                    map(@view project.analytes[_])
    end
    mass = @p gana map(map(mw, _))
    ret = @p gana map(map(rt, _))
    p = scatter()
    if all
        id = setdiff(eachindex(project.analytes), map(first ∘ parentindices, gana)...)
        others = project.analytes[id]
        isempty(others) || scatter!((@p others map((rt = rt(_), mw = mw(_))));
                    label = "Others", alpha = 0.2, hover = string.(id, ". ", repr.(others)), get_attributes!(kwargs)...)
    end
    for (key, analyte) in pairs(gana)
        scatter!(ret[key], mass[key]; label = "$key", hover = string.((first ∘ parentindices)(analyte), ": ", repr.(analyte)), get_attributes!(kwargs)...)
    end
    (model && haskey(project.appendix, :clusters_model)) || return scatter!(; xlabel, ylabel, title, legend = :outertopright, get_attributes!(kwargs)...)
    lim_mass = all ? extrema(mw.(project.analytes)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(mass))))
    lim_rt = all ? extrema(rt.(project.analytes)) : map(((f, x), ) -> f(x), zip([first, last], extrema(extrema.(ret))))
    cls = project.appendix[:clusters_model].mf.data.cluster |> unique
    mwrange = range(lim_mass..., length = 100)
    tbl = Table(mw = repeat(mwrange, length(cls)), cluster = repeat(cls, inner = length(mwrange)))
    tbl = Table(tbl, rt = predict(project.appendix[:clusters_model], tbl))
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
    all && scatter!((@p project.analytes[setdiff(eachindex(project.analytes), clusters...)] map((rt = rt(_), mw = mw(_))));
                    label = "Others", alpha = 0.2, get_attributes!(kwargs)...)
    for key in keys(clusters)
        scatter!(ret[key], mass[key]; label = "$key", get_attributes!(kwargs)...)
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