preis(anion = :acetate) = Project(AnalyteSP[], Data[], anion, Dictionary{ClassSP, Vector{Int}}(), 1, Dictionary{Symbol, Any}())
preis(
        featuretable,
        mz_range,
        polarity::Bool;
        db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
        db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
        mz_tol = 0.35,
        rt_tol = 0.1,
        anion = :acetate,
        data = -1,
        additional = Dict()
    ) = preis!(Project(AnalyteSP[], Data[], anion, Dictionary{ClassSP, Vector{Int}}(), 1, Dictionary{Symbol, Any}()), deepcopy(featuretable), mz_range, polarity; db, db_product, mz_tol, rt_tol, data, additional)

preis(
        project::Project,
        featuretable,
        mz_range,
        polarity::Bool;
        db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
        db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
        mz_tol = 0.35,
        rt_tol = 0.1,
        data = -1,
        additional = Dict()
    ) = preis!(deepcopy(project), deepcopy(featuretable), mz_range, polarity; db, db_product, mz_tol, rt_tol, data, additional)

function preprocess_preis!(data, data_id, featuretable, mz_range, mz_tol, polarity, additional)
    mzr = vectorize(mz_range)
    mz2v = Vector{Float64}[]
    mz2_loc = Vector{Int}[]
    for (id, mz2) in enumerate(featuretable.mz2)
        i = findfirst(x -> isapprox(mean(x), mz2; atol = 1e-4), mz2v)
        if isnothing(i)
            push!(mz2v, [mz2])
            push!(mz2_loc, [id])
        else
            push!(mz2v[i], mz2)
            push!(mz2_loc[i], id)
        end
    end
    mz2s = mean.(mz2v)
    length(mz2s) != length(mzr) && length(mzr) != 1 && throw(ArgumentError(
        "#mz_range does not match #mz2 (either equal to #mz2 or only 1)"
    ))
    length(mzr) ≡ 1 && (mzr = repeat(mzr, length(mz2s)))

    source, id_lb = length(data) ≡ 0 || data_id ≡ 0 ?
        (push!(data, PreIS(new_ft(featuretable, mz2_loc), mzr, mz2s, mz_tol, polarity, additional)); (length(data), 0)) :
            data_id < 0 ?
                merge_preis!(data, featuretable, mz2_loc, mzr, length(data) + 1 + data_id, mz_tol) :
                    merge_preis!(data, featuretable, mz2_loc, mzr, data_id, mz_tol)
    dt = data[source]
    sort!(dt.raw, [:mz2_id, :mz1, :rt])
    source, (@p dt.raw |> filterview(>(_.id, id_lb)))
end

function new_ft(featuretable, mz2_loc)
    Table(copy(featuretable);
        mz2_id = (@p eachindex(featuretable) |> map(findfirst(x -> in(_, x), mz2_loc))),
        mz2 = nothing,
        alignment = zeros(Float64, size(featuretable, 1)),
        isf = zeros(Int, size(featuretable, 1))
    )
end

function merge_preis!(data, featuretable, mz2_loc, mz_range, source, mz_tol)
    dt = data[source]
    id_lb = maximum(dt.raw.id)
    for loc in mz2_loc
        subft = @view featuretable[loc]
        mz2_id = findfirst(mz2 -> abs(subft.mz2[1] - mz2) < mz_tol, dt.mz2)
        mzb = popfirst!(mz_range)
        isnothing(mz2_id) ? push_ft!(dt, subft, mzb) :
            any(between(bd, dt.range[mz2_id]) for bd in mzb) ? merge_ft!(dt, subft, mzb, mz2_id) : push_ft!(dt, subft, mzb)
    end
    source, id_lb
end

function push_ft!(dt, subft, mzb)
    push!(dt.mz2, subft.mz2[1])
    push!(dt.range, mzb)
    subft = Table(copy(subft),
        id = maximum(dt.raw.id) .+ eachindex(subft),
        mz2_id = repeat([lastindex(dt.mz2)], size(subft, 1)),
        mz2 = nothing,
        alignment = zeros(Float64, size(subft, 1)),
        isf = zeros(Int, size(subft, 1))
    )
    append!(dt.raw, subft)
end

function merge_ft!(dt, subft, mzb, mz2_id)
    lb, ub = mzb
    dt.range[mz2_id] = (min(dt.range[mz2_id][1], lb), max(dt.range[mz2_id][2], ub))
    last_id = maximum(dt.raw.id)
    origin = @p dt.raw filter(≡(_.mz2_id, mz2_id))
    setdiff!(dt.raw, origin)
    subft = Table(copy(subft);
        mz2_id = repeat([mz2_id], size(subft, 1)),
        mz2 = nothing,
        alignment = zeros(Float64, size(subft, 1)),
        isf = zeros(Int, size(subft, 1))
    )
    for ft_id in eachindex(subft)
        rt = subft.rt[ft_id]
        mz1 = subft.mz1[ft_id]
        id = findfirst(id -> abs(origin.rt[id] - rt) <= rt_tol && abs(origin.mz1[id] - mz1) <= mz_tol, eachindex(origin))
        if isnothing(id)
            push!(origin, subft[ft_id])
            last_id += 1
            origin.id[end] = last_id
        else
            # record # to recalculate
            ft = (; subft[ft_id]..., id = origin.id[id])
            origin[id] = (; (x => mean(v) for (x, v) in zip(propertynames(origin), zip(origin[id], ft)))...)
        end
    end
    append!(dt.raw, origin)
end

function generate_cpd(db, db_id, products)
    cpds = compoundspvanilla.(Ref((Species = db.Species[db_id],
                            Abbreviation = db.Abbreviation[db_id],
                            Adduct = db.Adduct[db_id])), products)
    filter!(!isnothing, cpds)
    n = length(cpds)
    s = 1
    component = [[i] for i in eachindex(cpds)]
    while n > 0
        for i in s:n
            (!iscompatible(cpds[s], cpds[i]) || issubset(component[s], component[i])) && continue
            push!(cpds, union(cpds[s], cpds[i]))
            push!(component, union(component[s], component[i]))
            n += 1
        end
        n -= 1
        s += 1
    end
    cpds
end

function agg_cpd!(cpd, subanalytes, rt)
    existids = @views [findfirst(connected_cpd -> iscompatible(cpd, connected_cpd), a) for a in subanalytes]
    ids = findall(!isnothing, existids)
    if isempty(ids)
        connected_ids = [find_connected(cpd, analyte) for analyte in subanalytes]
        all(isempty, connected_ids) && (push!(cpd.project, AnalyteSP([cpd], rt)); return)
        n = length(cpd.project)
        for id in findall(!isempty, connected_ids)
            connected_id = connected_ids[id]
            analyte = subanalytes[id]
            id_new = findfirst(new -> iscompatible(last(new), cpd), @view cpd.project[n + 1:end])
            if !isnothing(id_new)
                new = cpd.project[n + id_new]
                union!(new, copy_wo_project.(analyte[connected_id]))
                new.rt = calc_rt(new)
                continue
            end
            length(connected_id) < length(analyte) && (analyte = push!(cpd.project, AnalyteSP(copy_wo_project.(analyte[connected_id]), analyte.rt)) |> last)
            push_cpd!(analyte, cpd)
        end
    else
        existids = @view existids[ids]
        existanalytes = @view subanalytes[ids]
        for (id, analyte) in zip(existids, existanalytes)
            isnothing(id) ? push_cpd!(analyte, cpd) : union!(analyte, id, cpd) # copy?
        end
    end
end

function preis!(
                    project::Project,
                    featuretable,
                    mz_range,
                    polarity::Bool;
                    db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
                    db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                    mz_tol = 0.35,
                    rt_tol = 0.1,
                    data_id = -1,
                    additional = Dict()
                )

    # check data compatibility ?
    source, ft = preprocess_preis!(project.data, data_id, featuretable, mz_range, mz_tol, polarity, additional)
    dt = project.data[source]
    for analyte in project
        analyte.rt = calc_rt(analyte)
    end

    for subft in groupview(getproperty(:mz2_id), ft)
        products = id_product(dt.mz2[subft.mz2_id[1]], polarity; mz_tol, db = db_product)
        printstyled("PreIS> ", color = :green, bold = true)
        println(products)
        for ft_id in eachindex(subft)
            cpdsvanilla = mapreduce(vcat, eachindex(db)) do db_id
                abs(db.mz[db_id] - subft.mz1[ft_id]) > mz_tol ? CompoundSPVanilla[] :
                    generate_cpd(db, db_id, products)
            end
            isempty(cpdsvanilla) && continue
            cpds = convert.(CompoundSP, cpdsvanilla)
            for cpd in cpds
                cpd.project = project
                cpd.fragments.id .= subft.id[ft_id]
                cpd.fragments.source .= source
                cpd.area = (subft.area[ft_id], subft.error[ft_id])
            end
            rt = subft.rt[ft_id]
            # Unidentified
            subanalytes = @views [analyte for analyte in project if abs(analyte.rt - rt) <= rt_tol]
            isempty(subanalytes) && (foreach(cpd -> push!(project, AnalyteSP([cpd], rt)), cpds); continue)
            # Aggregates
            foreach(cpd -> agg_cpd!(cpd, subanalytes, rt), cpds)
        end
    end
    project
end

function finish_profile!(project::Project; rt_tol = 0.1, err_tol = 0.3)
    printstyled("PreIS> ", color = :green, bold = true)
    println("Sorting compounds")
    for analyte in project
        sort!(analyte, lt = isless_class)
        for cpd in analyte
            sort!(cpd.fragments, :ion1; lt = isless_ion)
        end
    end
    printstyled("PreIS> ", color = :green, bold = true)
    println("Merging, splitting and deleting analytes")
    del = Int[]
    for (i, analyte) in enumerate(project)
        area_error = last(analyte).area
        todel = false
        for a in @view project[setdiff(eachindex(project), del, i)]
            abs(a.rt - analyte.rt) > rt_tol && continue
            id = findfirst(cpd -> iscompatible(cpd, last(analyte)), a)
            isnothing(id) && continue
            any(cpd.area[1] >= area_error[1] && cpd.area[2] <= err_tol for cpd in @view a[id:end]) || continue
            todel = true
            union!(a, analyte)
            a.rt = calc_rt(a)
        end

        todel && (push!(del, i); continue)

        for (i, cpd) in Iterators.reverse(enumerate(analyte))
            cpd.area[1] < area_error[1] && continue
            cpd.area[2] > err_tol && continue
            any(iscompatible(last(a), cpd) for a in project if abs(a.rt - analyte.rt) <= rt_tol) && continue
            push!(project, AnalyteSP(copy_wo_project.(analyte[1:i]), analyte.rt))
            last(project).rt = calc_rt(last(project))
            break
        end
        analyte.states[states_id(:error)] = area_error[2] > err_tol ? -1 : 1
        analyte.states[states_id(:isf)] = any(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), last(analyte).fragments.ion1) ? 1 : -1
    end
    unique!(del)
    deleteat!(project, del)
    del = Int[]
    for (i, analyte) in enumerate(project)
        analyte.states[states_id(:isf)] < 0 || continue
        todel = any(@view project[setdiff(eachindex(project), i)]) do a
            iscompatible(last(analyte), last(a))
        end
        todel && push!(del, i)
    end
    unique!(del)
    deleteat!(project, del)
    assign_isf_parent!(project)
    for analyte in project
        for cpd in analyte
            cpd.project = project
        end
    end
    project
end

assign_isf_parent!(project::Project) =
    for analyte in project
        for cpd in @view analyte[1:end - 1]
            set_raw!.(Ref(project), cpd.fragments.source, cpd.fragments.id, :isf, 1)
        end
        cpd = last(analyte)
        id = findall(ion -> !in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragments.ion1)
        set_raw_if!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, -1, <(1))
        id = setdiff(eachindex(cpd.fragments), id)
        set_raw!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, 1)
    end

assign_parent!(project::Project) =
    for analyte in project
        cpd = last(analyte)
        id = findall(ion -> !in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragments.ion1)
        set_raw!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, -1)
    end