function cluster_ion!(project, data_id; formula = @formula(rt ~ mz1), dev = 1, radius = 0.0, error_tol = 0.5, extend = true, pval_cluster = 0.001, kwargs...)
    haskey(project.data[data_id].additional, :clusters) && delete!(project.data[data_id].additional, :clusters)
    tbl = @p project.data[data_id].raw filterview(_.isf < 0 && _.error <= error_tol) deepcopy
    maxrt = round(maximum(tbl.rt))
    radius = radius > 0 ? radius : maxrt / 20
    scale = - foldl(-, extrema(tbl.mz1)) / maxrt
    tbl.mz1 ./= scale
    clusters = map(dbscan(hcat(tbl.rt, tbl.mz1)', radius; kwargs...)) do cluster
        cluster.core_indices
    end
   
    # test adjacent lm(rt ~ mz1 ^2 * cluster) vs lm(rt ~ mz1 ^2) if not significant, merge clusters
    uswts = (1 .- 1 ./ (1 .+ abs.(tbl.mz1[clusters[1]] .- median(tbl.mz1[clusters[1]]))))
    wts =  uswts .* length(tbl.mz1[clusters[1]]) ./ sum(uswts)
    models = [lm(formula, tbl[clusters[1]]; wts)]
    if length(clusters) > 1
        f = @match formula.rhs begin
            ::Tuple => FormulaTerm(formula.lhs, (Term(:cluster), formula.rhs..., map(formula.rhs) do r
                    @match r begin
                        ::Tuple => InteractionTerm((r..., Term(:cluster)))
                        _       => InteractionTerm((r, Term(:cluster)))
                    end
                end...))
            _       => FormulaTerm(formula.lhs, (Term(:cluster), formula.rhs, InteractionTerm((formula.rhs, Term(:cluster)))))
        end
        i = 1
        while i < length(clusters)
            j = i + 1
            subt = Table(tbl[union(clusters[i], clusters[j])], 
                        cluster = vcat(repeat([0], length(clusters[i])), repeat([1], length(clusters[j])))
                    )
            uswts = (1 .- 1 ./ (1 .+ abs.(subt.mz1 .- median(subt.mz1)))) 
            wts = uswts .* length(subt.mz1) ./ sum(uswts)
            l1 = lm(formula, subt; wts)
            l2 = lm(f, subt; wts)
            if last(pval(anova(l1, l2))) <= pval_cluster
                uswts = (1 .- 1 ./ (1 .+ abs.(tbl.mz1[clusters[j]] .- median(tbl.mz1[clusters[j]]))))
                wts =  uswts .* length(tbl.mz1[clusters[j]]) ./ sum(uswts)
                push!(models, lm(formula, tbl[clusters[j]]; wts))
                i += 1
            else
                union!(clusters[i], popat!(clusters, j))
                models[i] = l1
            end
        end
    end
    extend || return project.data[data_id].additional[:clusters] = Table(
        class = repeat(Any[ClassSP], length(clusters)),
        id = map(clusters) do cluster
            tbl.id[cluster]
        end
    )
    assign_parent!(project)
    ids = map(clusters) do cluster
        tbl.id[cluster]
    end
    tbl = @p project.data[data_id].raw filterview(_.isf < 0) deepcopy
    tbl.mz1 ./= scale
    clusters = map(ids) do id
        findall(x -> in(x, id), tbl.id)
    end
    to_clustered = setdiff(eachindex(tbl.id), union(clusters...))
    distances = [abs(tbl.rt[i] - predict(model, tbl[i:i])[1]) for (i, model) in Iterators.product(to_clustered, models)]
    while !isempty(to_clustered)
        dist, id = findmin(distances)
        i, j = Tuple(id)
        dist > dev && break
        push!(clusters[j], to_clustered[i])
        distances = distances[setdiff(eachindex(to_clustered), [i]), :]
        deleteat!(to_clustered, i)
        uswts = (1 .- 1 ./ (1 .+ abs.(tbl.mz1[clusters[j]] .- median(tbl.mz1[clusters[j]]))))
        wts =  uswts .* length(tbl.mz1[clusters[j]]) ./ sum(uswts)
        models[j] = lm(formula, tbl[clusters[j]]; wts)
        distances[:, j] = map(to_clustered) do point
            abs.(tbl.rt[point] .- predict(models[j], tbl[point:point])[1])
        end
    end
    project.data[data_id].additional[:clusters] = Table(
            class = repeat(Any[ClassSP], length(clusters)),
            id = map(clusters) do cluster
                tbl.id[cluster]
            end
        )
end

set_cluster_class!(project, data_id, iterated) = 
    haskey(project.data[data_id].additional, :clusters) ? (project.data[data_id].additional[:clusters].class .= iterated) : throw(ErrorException("Data_$data_id does not have clustering results. Try run `cluster_ion!(project, $data_id)`"))

function apply_cluster!(project; mode = ≥(0.5), data_id = 1, analytes = project.analytes)
    haskey(project.data[data_id].additional, :clusters) || throw(ErrorException("Data_$data_id does not have clustering results. Try run `cluster_ion!(project, $data_id)` and `set_cluster_class!(project, $data_id, class"))
    classes = map(analytes) do analyte
        cpd = last(analyte)
        frag = filterview(x -> x.source == data_id, cpd.fragments)
        classes = map(Iterators.reverse(frag.id)) do id
            id_cluster = findfirst(x -> in(id, x), project.data[data_id].additional[:clusters].id)
            isnothing(id_cluster) ? nothing : project.data[data_id].additional[:clusters].class[id_cluster]
        end
        filter!(!isnothing, classes)
    end
    vote!(analytes, classes, mode)
    project
end

vote!(analytes::Vector{AnalyteSP}, classes::Vector{Vector}, mode::T) where {S, T <: Base.Fix2{S, Float64}} = 
    for (analyte, class) in zip(analytes, classes)
        isempty(class) && (analyte.states[3] = 0; continue)
        analyte.states[3] = mode(count(x -> isa(last(analyte).class, x), class) / length(class)) ? 1 : -1
    end

vote!(analytes::Vector{AnalyteSP}, classes::Vector{Vector}, mode::T) where {S, T <: Base.Fix2{S, Int}} = 
    for (analyte, class) in zip(analytes, classes)
        isempty(class) && (analyte.states[3] = 0; continue)
        analyte.states[3] = mode(count(x -> isa(last(analyte).class, x), class)) ? 1 : -1
    end
#=

function align!(project::Project, align_to::id, compatible = false, db_product = SPDB[mrm.polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                mode::Symbol = :default)
    project.alignment = align_to
    products = @p mrm.mz2 map(id_product(_, mrm.polarity; db = db_product, mz_tol = mrm.mz_tol))
    if compatible
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscompatible(cpd, product) for product in products[i]) && 
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ==(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    else
        qf = cpd -> begin
            _, ms1 = generate_ms(cpd, mode, aquery.project.anion)
            any(!isempty(products[i]) && any(iscomponent(cpd, product) for product in products[i]) && 
                any(any(between(ms, mz1, mrm.mz_tol) for ms in ms1) for mz1 in filterview(x -> ==(x.mz2_id, i), mrm.raw).mz1) for i in eachindex(products))
        end
    end
=#