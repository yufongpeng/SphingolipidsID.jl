plot_rt_mz1(tbl; legend = false, xlabel = "Retention Time (min)", ylabel = "m/z", title = "Features", kwargs...) = 
    scatter(tbl.rt, tbl.mz1; legend, xlabel, ylabel, title, kwargs...)

function plot_rt_mz1(project, data_id; cluster = true, xlabel = "Retention Time (min)", ylabel = "m/z", title = "Features of Data_$data_id", kwargs...)
    if haskey(project.data[data_id].additional, :clusters) && cluster
        kwargs = Dict{Symbol, Any}(kwargs)
        label = pop!(kwargs, :label, cluster_label(project.data[data_id].additional[:clusters]))
        color = pop!(kwargs, :color, repeat([:auto], length(label) + 1))
        scatter()
        id = setdiff(project.data[data_id].raw.id, union(project.data[data_id].additional[:clusters].id...))
        scatter!((@p project.data[data_id].raw filterview(in(_.id, id)) getproperties(__, (:rt, :mz1)))
                    ; label = "None", color = popfirst!(color), alpha = 0.2)
        for id in project.data[data_id].additional[:clusters].id
            scatter!((@p project.data[data_id].raw filterview(in(_.id, id)) getproperties(__, (:rt, :mz1)))
                    ; label = popfirst!(label), color = popfirst!(color))
        end
        scatter!(; xlabel, ylabel, title, NamedTuple(kwargs)...) |> display
    else
        scatter(project.data[data_id].raw.rt, project.data[data_id].raw.mz1; xlabel, ylabel, title, kwargs...) |> display
    end
end

cluster_label(clusters) = 
    allequal(clusters.class) ? collect(eachindex(clusters)) : union_repr.(clusters.class)

union_repr(x) = repr(x)
union_repr(x::Union) = join(map(propertynames(x)) do i
    union_repr(getproperty(x, i))
end, ", ")

