preis(anion = :acetate) = Project(AnalyteSP[], Data[], anion, 1)
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
    ) = preis!(Project(AnalyteSP[], Data[], anion, 1), deepcopy(featuretable), mz_range, polarity; db, db_product, mz_tol, rt_tol, data, additional)

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

function preis!(
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
                )

    # check data compatibility ?
    mz_range = vectorize(mz_range)
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
    length(mz2s) != length(mz_range) && length(mz_range) != 1 && throw(ArgumentError(
        "#mz_range does not match #mz2 (either equal to #mz2 or only 1)"
    ))
    length(mz_range) == 1 && (mz_range = repeat(mz_range, length(mz2s)))
    source = if length(project.data) == 0; 1
    elseif data == 0; length(project.data) + 1
    elseif data < 0; length(project.data) + 1 + data
    else; data end
    
    if source > length(project.data)
        featuretable = Table(copy(featuretable); 
            mz2_id = (@p eachindex(featuretable) |> map(findfirst(x -> in(_, x), mz2_loc))), 
            mz2 = nothing, 
            alignment = zeros(Float64, size(featuretable, 1)), 
            isf = zeros(Int, size(featuretable, 1))
        )
        push!(project.data, PreIS(featuretable, mz_range, mz2s, mz_tol, polarity, additional))
        dt = last(project.data)
        newid = 0
    else 
        dt = project.data[source]
        newid = maximum(dt.raw.id)
        for loc in mz2_loc
            subft = @view featuretable[loc]
            mz2_id = findfirst(mz2 -> abs(subft.mz2[1] - mz2) < mz_tol, dt.mz2)
            if isnothing(mz2_id)
                push!(dt.mz2, subft.mz2[1])
                push!(dt.range, popfirst!(mz_range))
                subft = Table(copy(subft), 
                    id = maximum(dt.raw.id) .+ eachindex(subft), 
                    mz2_id = repeat([lastindex(dt.mz2)], size(subft, 1)), 
                    mz2 = nothing, 
                    alignment = zeros(Float64, size(subft, 1)),
                    isf = zeros(Int, size(subft, 1))
                )
                append!(dt.raw, subft)
            else
                prev_id = maximum(dt.raw.id)
                last_id = prev_id
                lb, ub = popfirst!(mz_range)
                if !between(lb, dt.range[mz2_id]) && !between(ub, dt.range[mz2_id])
                    push!(dt.mz2, subft.mz2[1])
                    push!(dt.range, (lb, ub))
                    subft = Table(copy(subft); 
                        mz2_id = repeat([lastindex(dt.mz2)], size(subft, 1)), 
                        mz2 = nothing, 
                        alignment = zeros(Float64, size(subft, 1)),
                        isf = zeros(Int, size(subft, 1))
                    )
                    continue
                end
                dt.range[mz2_id] = (min(dt.range[mz2_id][1], lb), max(dt.range[mz2_id][2], ub))
                subg = @p dt.raw |> group(==(_.mz2_id, mz2_id))
                setdiff!(dt.raw, subg[true])
                origin = subg[true]
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
        end
    end
    sort!(dt.raw, [:mz2_id, :mz1, :rt])
    featuretable = @p dt.raw |> filterview(>(_.id, newid))
    for analyte in project
        analyte.rt = calc_rt(analyte)
    end

    for subft in groupview(getproperty(:mz2_id), featuretable)
        products = id_product(dt.mz2[subft.mz2_id[1]], polarity; mz_tol, db = db_product)
        printstyled("PreIS> ", color = :green, bold = true)
        println(products)
        for ft_id in eachindex(subft)
            current_cpd = CompoundSP[]
            area_error = (subft.area[ft_id], subft.error[ft_id])
            mz1 = subft.mz1[ft_id]
            id = subft.id[ft_id]
            rt = subft.rt[ft_id]
            for db_id in eachindex(db)
                abs(db.mz[db_id] - mz1) > mz_tol && continue
                cpds = CompoundSP.(Ref(project), Ref((Species = db.Species[db_id], 
                                                    Abbreviation = db.Abbreviation[db_id], 
                                                    Adduct = db.Adduct[db_id])), 
                                    products, source, id, Ref(area_error))
                filter!(!isnothing, cpds)
                n = length(cpds)
                s = 1
                component = [[i] for i in eachindex(cpds)]
                while n > 0
                    cpd1 = cpds[s]
                    for i in s:n
                        (!iscompatible(cpd1, cpds[i]) || issubset(component[s], component[i])) && continue
                        push!(cpds, union(cpd1, cpds[i]))
                        push!(component, union(component[s], component[i]))
                        n += 1
                    end
                    n -= 1
                    s += 1
                end
                current_cpd = vcat(current_cpd, cpds)
            end
            # Unidentified
            isempty(current_cpd) && continue
            subanalytes = @views [analyte for analyte in project if abs(analyte.rt - rt) <= rt_tol]
            if isempty(subanalytes) 
                for cpd in current_cpd
                    push!(project, AnalyteSP([cpd], rt, [0, 0, 0, 0, 0]))
                end
                continue
            end
            # Aggregates
            for cpd in current_cpd
                agg = false
                existids = @views [findfirst(connected_cpd -> iscompatible(cpd, connected_cpd), a) for a in subanalytes]
                ids = findall(!isnothing, existids)
                if isempty(ids)
                    n = length(project)
                    for analyte in subanalytes
                        connected_id = find_connected(cpd, analyte)
                        isempty(connected_id) && continue
                        agg = true
                        pushed = false
                        for new in @view project[n + 1:end]
                            iscompatible(last(new), cpd) || continue
                            union!(new, copy_wo_project.(analyte[connected_id]))
                            new.rt = calc_rt(new)
                            pushed = true
                            break
                        end
                        pushed && continue
                        if length(connected_id) < length(analyte)
                            push!(project, AnalyteSP(copy_wo_project.(analyte[connected_id]), analyte.rt, [0, 0, 0, 0, 0]))
                            analyte = last(project)
                        end
                        push!(analyte, copy_wo_project(cpd))
                        sort!(analyte, lt = isless_class)
                        analyte.rt = calc_rt(analyte)
                    end
                else
                    agg = true
                    existids = @view existids[ids]
                    existanalytes = @view subanalytes[ids]
                    for (id, analyte) in zip(existids, existanalytes)
                        if isnothing(id)
                            push!(analyte, copy_wo_project(cpd))
                            sort!(analyte, lt = isless_class)
                        else
                            union!(project, analyte, id, copy_wo_project(cpd))
                        end
                        analyte.rt = calc_rt(analyte)
                    end
                end
                agg || push!(project, AnalyteSP([cpd], rt, [0, 0, 0, 0, 0]))
            end
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
            push!(project, AnalyteSP(copy_wo_project.(analyte[1:i]), analyte.rt, [0, 0, 0, 0, 0]))
            last(project).rt = calc_rt(last(project))
            break
        end
        analyte.states[4] = area_error[2] > err_tol ? 1 : -1        
        analyte.states[5] = any(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), last(analyte).fragments.ion1) ? -1 : 1
    end
    unique!(del)
    deleteat!(project, del)
    del = Int[]
    for (i, analyte) in enumerate(project)
        analyte.states[5] > 0 || continue
        cpd1 = last(analyte)
        todel = any(@view project[setdiff(eachindex(project), i)]) do a
            cpd2 = last(a)
            isclasscompatible(cpd1.class, cpd2.class) && ischaincompatible(cpd1.chain, cpd2.chain)
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
        id = findall(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragments.ion1)
        set_raw_if!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, -1, <(1))
        id = setdiff(eachindex(cpd.fragments), id)
        set_raw!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, 1)
    end

assign_parent!(project::Project) = 
    for analyte in project
        cpd = last(analyte)
        id = findall(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragments.ion1)
        set_raw!.(Ref(project), cpd.fragments.source[id], cpd.fragments.id[id], :isf, -1)
    end