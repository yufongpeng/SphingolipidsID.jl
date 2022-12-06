preis(anion = :acetate) = Project(AnalyteSP[], Data[], anion)
preis(
        featuretable, 
        mz_range,
        ms2, 
        polarity::Bool;
        db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
        db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
        mz_tol = 0.35, 
        rt_tol = 0.1,
        anion = :acetate,
        data = -1,
        additional = Dict()
    ) = preis!(Project(AnalyteSP[], Data[], anion), deepcopy(featuretable), mz_range, ms2, polarity; db, db_product, mz_tol, rt_tol, data, additional)

preis(
        project::Project,
        featuretable, 
        mz_range,
        ms2, 
        polarity::Bool;
        db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
        db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
        mz_tol = 0.35, 
        rt_tol = 0.1,
        data = -1,
        additional = Dict()
    ) = preis!(deepcopy(project), deepcopy(featuretable), mz_range, ms2, polarity; db, db_product, mz_tol, rt_tol, data, additional)

function preis!(
                    project::Project,
                    featuretable, 
                    mz_range,
                    ms2, 
                    polarity::Bool;
                    db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG],
                    db_product = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                    mz_tol = 0.35, 
                    rt_tol = 0.1,
                    data = -1,
                    additional = Dict()
                )
    products = id_product(ms2, polarity; mz_tol, db = db_product)
    featuretable = deepcopy(featuretable)
    printstyled("PreIS> ", color = :green, bold = true)
    println(products)
    # check data compatibility ?
    source = if length(project.data) == 0; 1
            elseif data == 0; length(project.data) + 1
            elseif data < 0; length(project.data) + 1 + data
            else; data end
    
    if source > length(project.data)
        featuretable.scan .= 1
        sort!(featuretable, [:mz1, :rt])
        push!(project.data, PreIS(featuretable, [mz_range], [ms2], mz_tol, polarity, additional))
    else 
        dt = project.data[source]
        scan_id = findfirst(mz2 -> abs(ms2 - mz2) < mz_tol, dt.mz2)
        if isnothing(scan_id)
            push!(dt.mz2, ms2)
            push!(dt.range, mz_range)
            featuretable.scan .= lastindex(dt.mz2)
            featuretable.id .= maximum(dt.raw.id) .+ (1:size(featuretable, 1))
            append!(dt.raw, featuretable)
            sort!(dt.raw, [:mz1, :rt])
            featuretable = filterview(x -> ==(x.scan, lastindex(dt.mz2)), dt.raw)
        else
            origin = filterview(x -> ==(x.scan, lastindex(dt.mz2)), dt.raw)
            prev_id = maximum(origin.id)
            last_id = prev_id
            featuretable.scan .= scan_id
            for ft_id in eachindex(featuretable)
                rt = featuretable.rt[ft_id]
                mz1 = featuretable.mz1[ft_id]
                id = findfirst(id -> abs(origin.rt[id] - rt) <= rt_tol && abs(origin.mz1[id] - mz1) <= mz_tol, eachindex(origin))
                if isnothing(id)
                    push!(project.data[source].raw, featuretable[ft_id])
                    last_id += 1
                    origin.id[end] = last_id
                else
                    ft = (; featuretable[ft_id]..., id = origin.id[id])
                    origin[id] = (; (x => mean(v) for (x, v) in zip(propertynames(origin), zip(origin[id], ft)))...)
                end
            end
            sort!(dt.raw, [:mz1, :rt])
            featuretable = filterview(x -> >(x.id, prev_id), dt.raw)
        end
    end

    for analyte in project
        analyte.rt = calc_rt(analyte)
    end

    for ft_id in eachindex(featuretable)
        current_cpd = CompoundSP[]
        area_error = (featuretable.area[ft_id], featuretable.error[ft_id])
        mz1 = featuretable.mz1[ft_id]
        id = featuretable.id[ft_id]
        rt = featuretable.rt[ft_id]
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
                push!(project, AnalyteSP([cpd], rt, [0, 0]))
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
                        push!(project, AnalyteSP(copy_wo_project.(analyte[connected_id]), analyte.rt, [0, 0]))
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
            agg || push!(project, AnalyteSP([cpd], rt, [0, 0]))
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
        todel = true
        if area_error[2] <= err_tol
            todel = false
            for a in @view project[setdiff(eachindex(project), del, i)]
                abs(a.rt - analyte.rt) > rt_tol && continue
                id = findfirst(cpd -> iscompatible(cpd, last(analyte)), a)
                isnothing(id) && continue
                any(cpd.area[1] >= area_error[1] for cpd in @view a[id + 1:end]) || continue
                push!(del, i)
                todel = true
                union!(a, analyte)
                a.rt = calc_rt(a)
            end
            todel && continue
        end

        for (i, cpd) in Iterators.reverse(enumerate(analyte))
            cpd.area[1] < area_error[1] && continue
            cpd.area[2] > err_tol && continue
            any(iscompatible(last(a), cpd) for a in project if abs(a.rt - analyte.rt) <= rt_tol) && continue            
            push!(project, AnalyteSP(copy_wo_project.(analyte[1:i]), analyte.rt, [0, 0]))
            last(project).rt = calc_rt(last(project))
            break
        end
        todel && continue
        any(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), last(analyte).fragments.ion1) || push!(del, i)
    end
    unique!(del)
    deleteat!(project, del)
    for analyte in project
        for cpd in analyte
            cpd.project = project
        end
    end
    project
end