preis(anion = :acetate) = Project(AnalyteGSL[], Data[], anion)
preis(
        featuretable, 
        mz_range,
        ms2, 
        polarity::Bool,
        eV = -1;
        db = polarity ? DEFAULT_POS_LIBRARY : DEFAULT_NEG_LIBRARY,
        db_product = polarity ? DEFAULT_POS_FRAGMENT : DEFAULT_NEG_FRAGMENT,
        mz_tol = 0.35, 
        rt_tol = 0.1,
        anion = :acetate,
        data = -1,
        additional = Dict()
    ) = preis!(Project(AnalyteGSL[], Data[], anion), deepcopy(featuretable), mz_range, ms2, polarity, eV; db, mz_tol, rt_tol, data, additional)

preis(
        project::Project,
        featuretable, 
        mz_range,
        ms2, 
        polarity::Bool,
        eV = -1;
        db = polarity ? DEFAULT_POS_LIBRARY : DEFAULT_NEG_LIBRARY,
        db_product = polarity ? DEFAULT_POS_FRAGMENT : DEFAULT_NEG_FRAGMENT,
        mz_tol = 0.35, 
        rt_tol = 0.1,
        data = -1,
        additional = Dict()
    ) = preis!(deepcopy(project), deepcopy(featuretable), mz_range, ms2, polarity, eV; db, mz_tol, rt_tol, data, additional)

function preis!(
                    project::Project,
                    featuretable, 
                    mz_range,
                    ms2, 
                    polarity::Bool,
                    eV = -1;
                    db = polarity ? DEFAULT_POS_LIBRARY : DEFAULT_NEG_LIBRARY,
                    db_product = polarity ? DEFAULT_POS_FRAGMENT : DEFAULT_NEG_FRAGMENT,
                    mz_tol = 0.35, 
                    rt_tol = 0.1,
                    data = -1,
                    additional = Dict()
                )
    products = id_product(ms2, polarity; mz_tol, db = db_product)
    featuretable = deepcopy(featuretable)
    println("PreIS> ", products)
    analytes = project.analytes
    # check data compatibility ?
    if length(project.data) == 0
        source = 1
    elseif data == 0
        source = length(project.data) + 1
    elseif data < 0
        source = length(project.data) + 1 + data
    else
        source = data
    end
    # Custom eV
    featuretable.eV .= eV
    if source > length(project.data)
        featuretable.mz2 .= 1
        push!(project.data, PreIS(featuretable, [mz_range], [ms2], mz_tol, polarity, additional))
    else 
        dt = project.data[source]
        mz2_id = findfirst(mz2 -> abs(ms2 - mz2) < mz_tol, dt.mz2)
        if isnothing(mz2_id)
            push!(dt.mz2, ms2)
            push!(dt.range, mz_range)
            featuretable.mz2 .= lastindex(dt.mz2)
            featuretable.id .= maximum(dt.raw.id) .+ (1:size(featuretable, 1))
            append!(dt.raw, featuretable)
            sort!(dt.raw, [:mz, :rt])
            featuretable = filter(:mz2 => ==(lastindex(dt.mz2)), dt.raw, view = true)
        else
            origin = filter(:mz2 => ==(lastindex(dt.mz2)), dt.raw)
            prev_id = maximum(origin.id)
            last_id = prev_id
            for ft in featuretable
                id = findfirst(row -> abs(row.rt - ft.rt) < rt_tol && abs(row.mz - ft.mz) < mz_tol, eachrow(origin))
                if isnothing(id)
                    push!(project.data.raw[mz2_id], ft)
                    last_id += 1
                    origin.id[end] = last_id
                else
                    ft.id = origin.id[id]
                    origin[id, :] = [mean(x) for x in zip(origin[id, :], ft)]
                end
            end
            sort!(dt.raw, [:mz, :rt])
            featuretable = filter(:id => >(prev_id), dt.raw, view = true)
        end
    end

    for analyte in analytes
        analyte.rt = mean(mean(query_raw(project, source, row.id).rt for row in eachrow(cpd.fragments)) for cpd in analyte.identification)
    end
    
    for ft in eachrow(featuretable)
        current_cpd = CompoundGSL[]
        for row in eachrow(db)
            if abs(row.var"m/z" - ft.mz) < mz_tol
                cpd = CompoundGSL(row, deepcopy(products), source, ft.id)
                isnothing(cpd) || push!(current_cpd, cpd)
            end
        end
        
        # Unidentified
        isempty(current_cpd) && continue

        subanalytes = @views [analyte for analyte in analytes if abs(analyte.rt - ft.rt) < rt_tol]

        if isempty(subanalytes) 
            for cpd in current_cpd
                push!(analytes, AnalyteGSL([cpd], ft.rt))
            end
        else
            # Aggregates
            for cpd in current_cpd
                agg = false
                for analyte in subanalytes
                    connected_id = find_connected(cpd, analyte)
                    isempty(connected_id) && continue
                    agg = true
                    if length(connected_id) < length(analyte.identification)
                        push!(analytes, AnalyteGSL(deepcopy(analyte.identification[connected_id]), analyte.rt))
                        analyte = last(analytes)
                    end
                    id = findfirst(connected_cpd -> iscompatible(cpd, connected_cpd), analyte.identification)
                    if isnothing(id)
                        push!(analyte.identification, deepcopy(cpd))
                    else
                        union!(analyte.identification[id], deepcopy(cpd))
                    end
                    analyte.rt = mean(mean(query_raw(project, source, row.id).rt for row in eachrow(cpd.fragments)) for cpd in analyte.identification)
                end
                agg || push!(analytes, AnalyteGSL([cpd], ft.rt))
            end
        end
        #=
        del = Int[]
        for (i, analyte) in enumerate(analytes)
            any(x.identification == analyte.identification for x in @view analytes[1:i - 1]) && push!(del, i)
        end
        deleteat!(analytes, del)
        =#
    end
    project
end

function finish_profile!(project::Project)
    println("PreIS> Sorting and deleting analytes")
    del = Int[]
    for (i, analyte) in enumerate(project.analytes)
        sort!(analyte.identification, lt = isless_class)
        any(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), @views last(analyte.identification).fragments[!, :ion1]) || (push!(del, i); continue)
        for cpd in analyte.identification
            sort!(cpd.fragments, :ion1, lt = isless_ion)
        end
    end
    deleteat!(project.analytes, del)
    project
end