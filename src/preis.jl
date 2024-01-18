function PreIS(featuretable::Table;
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 3, 
                signal = :area,
                anion = :acetate,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:mz_range, :datafile], [splat(union), nothing])
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    PreIS!(deepcopy(Table(featuretable; polarity = nothing)), polarity;
            rt_tol, 
            mz_tol, 
            n, 
            signal,
            anion,
            est_fn, 
            err_fn, 
            err_tol,
            other_fn
        )
end

function PreIS!(featuretable::Table;
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 3, 
                signal = :area,
                anion = :acetate,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:mz_range, :datafile], [splat(union), nothing])
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    PreIS!(Table(featuretable; polarity = nothing), polarity;
            rt_tol, 
            mz_tol, 
            n, 
            signal,
            anion,
            est_fn, 
            err_fn, 
            err_tol,
            other_fn
        )
end

function PreIS!(featuretable::Table, polarity::Bool;
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 3, 
                signal = :area,
                anion = :acetate,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:mz_range, :datafile], [splat(union), nothing])
    )
    featuretable = filter_combine_features!(featuretable; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
    config = dictionary(pairs((; rt_tol, mz_tol, n, signal, anion, est_fn, err_fn, err_tol, other_fn)))
    #=
    mz2v = Vector{Float64}[]
    mz2_loc = Vector{Int}[]
    mzr = RealIntervals[]
    for (id, row) in enumerate(featuretable)
        i = findfirst(x -> isapprox(mean(x), row.mz2; atol = 1e-4), mz2v)
        if isnothing(i)
            push!(mz2v, [row.mz2])
            push!(mz2_loc, [id])
            push!(mzr, row.mz_range)
        else
            push!(mz2v[i], row.mz2)
            push!(mz2_loc[i], id)
            mzr[i] = union(mzr[i], row.mz_range)
        end
    end
    mz2s = mean.(mz2v)
    =#
    mz2_loc = @p featuretable groupfind(x -> round.(getproperty(x, :mz2); digits = -floor(Int, log10(mz_tol)))) collect
    mzr = map(i -> union(featuretable.mz_range[i]...), mz2_loc)
    mz2s = map(i -> mean(featuretable.mz2[i]), mz2_loc)
    PreIS(new_ft(featuretable, collect(mz2_loc)), mzr, mz2s, polarity, config)
end

function Project(preis::PreIS; data_id = -1, kwargs...)
    appendix = Dictionary{Symbol, Any}(kwargs)
    project = Project(AnalyteSP[], AbstractData[], nothing, appendix)
    preis!(project, preis; data_id)
end

function preis!(project::Project, preis::PreIS; data_id = -1)
    # check config
    get(project.appendix, :anion, preis.config[:anion]) != preis.config[:anion] && throw(ArgumentError("Anions are different; can not add this data to the project"))
    get(project.appendix, :signal, preis.config[:signal]) != preis.config[:signal] && throw(ArgumentError("Signals are different; can not add this data to the project"))
    data = project.data
    if length(data) ≡ 0 || data_id ≡ 0
        push!(data, preis)
        dt = data[length(data)]
        idend = 0
        ft = (@p dt.table |> filterview(>(_.id, 0)))
        source = length(data)
    else
        source = data_id < 0 ? length(data) + 1 + data_id : data_id
        dt = data[source]
        for mz2 in dt.mz2
            any(x -> isapprox(x, mz2; atol = preis.config[:mz_tol]), preis.mz2) && throw(ArgumentError("The data contains mz2 value too close to existing mz2 value.\n current: $(dt.mz2)\n new: $(preis.mz2)"))
        end
        idend = maximum(dt.id)
        preis.table.id .+= idend
        preis.table.mz2_id .+= maximum(dt.mz2_id)
        append!(dt.mz2, preis.mz2)
        append!(dt.range, preis.range)
        append!(dt.table, preis.table)
    end
    set!(project.appendix, :anion, preis.config[:anion])
    set!(project.appendix, :signal, preis.config[:signal])
    sort!(dt.table, [:mz2_id, :mz1, :rt])
    ft = (@p dt.table |> filterview(>(_.id, idend)))
    for analyte in project
        analyte.rt = calc_rt(analyte)
    end

    db = SPDB[preis.polarity ? :LIBRARY_POS : :LIBRARY_NEG]
    for subft in groupview(getproperty(:mz2_id), ft)
        products = id_product(dt.mz2[subft.mz2_id[1]], preis.polarity; mz_tol = preis.config[:mz_tol])
        @info string("PreIS | ", products)
        for data in subft
            cpdsvanilla = mapreduce(vcat, eachindex(db)) do db_id
                abs(db.mz[db_id] - data.mz1) > preis.config[:mz_tol] ? CompoundSPVanilla[] :
                    generate_cpd(db, db_id, products)
            end
            isempty(cpdsvanilla) && continue
            cpds = convert.(CompoundSP, cpdsvanilla)
            for cpd in cpds
                cpd.project = project
                cpd.fragment.id .= data.id
                cpd.fragment.source .= source
                cpd.signal = (getproperty(data, project.appendix[:signal]), data.error)
            end
            # Unidentified
            subanalytes = @views [analyte for analyte in project if abs(analyte.rt - data.rt) <= preis.config[:rt_tol]]
            isempty(subanalytes) && (foreach(cpd -> push!(project, AnalyteSP([cpd], data.rt)), cpds); continue)
            # Aggregates
            foreach(cpd -> agg_cpd!(cpd, subanalytes, data.rt), cpds)
        end
    end
    project
end
"""
    preis(
            featuretable::Table,
            mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
            polarity::Bool;
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
            data_id = -1,
            appedix = Dictionary{Symbol, Any}()
        )

Create an empty `Project` and add precursor ion scan features. See `preis!` for details.
"""
preis(
        featuretable::Table,
        mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
        rt_tol = 0.1, 
        mz_tol = 0.35, 
        n = 3, 
        signal = :area,
        est_fn = mean, 
        err_fn = default_error, 
        err_tol = 0.5,
        other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
        data_id = -1,
        appedix = Dictionary{Symbol, Any}()
    ) = preis!(Project(; appedix), deepcopy(featuretable), mz_range, polarity; mz_tol, rt_tol, n, signal, est_fn, err_fn, err_tol, other_fn, data_id, appedix)

"""
    preis(
            project::Project,
            featuretable::Table,
            mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
            polarity::Bool;
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
            data_id = -1,
            appedix = Dictionary{Symbol, Any}()
        )

Copy the `project` and add precursor ion scan features. See `preis!` for details.
"""
preis(
        project::Project,
        featuretable::Table,
        mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
        polarity::Bool;
        rt_tol = 0.1, 
        mz_tol = 0.35, 
        n = 3, 
        signal = :area,
        est_fn = mean, 
        err_fn = default_error, 
        err_tol = 0.5,
        other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
        data_id = -1,
        appedix = Dictionary{Symbol, Any}()
    ) = preis!(deepcopy(project), deepcopy(featuretable), mz_range, polarity; mz_tol, rt_tol, n, signal, est_fn, err_fn, err_tol, other_fn, data_id, appedix)

function preprocess_preis!(data, data_id, featuretable, mz_range, polarity, config)
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
        (push!(data, PreIS(new_ft(featuretable, mz2_loc), mzr, mz2s, polarity, config)); (length(data), 0)) :
            data_id < 0 ?
                merge_preis!(data, featuretable, mz2_loc, mzr, length(data) + 1 + data_id, config[:mz_tol]) :
                    merge_preis!(data, featuretable, mz2_loc, mzr, data_id, config[:mz_tol])
    dt = data[source]
    sort!(dt.table, [:mz2_id, :mz1, :rt])
    source, (@p dt.table |> filterview(>(_.id, id_lb)))
end

function new_ft(featuretable, mz2_loc)
    Table(copy(featuretable);
        mz2_id = (@p eachindex(featuretable) |> map(findfirst(x -> in(_, x), mz2_loc))),
        mz2 = nothing,
        isf = zeros(Int, size(featuretable, 1)),
        n = ones(Int, size(featuretable, 1))
    )
end

function merge_preis!(data, featuretable, mz2_loc, mz_range, source, mz_tol)
    dt = data[source]
    id_lb = maximum(dt.table.id)
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
        id = maximum(dt.table.id) .+ eachindex(subft),
        mz2_id = repeat([lastindex(dt.mz2)], size(subft, 1)),
        mz2 = nothing,
        isf = zeros(Int, size(subft, 1)),
        n = ones(Int, size(subft, 1))
    )
    append!(dt.table, subft)
end

function merge_ft!(dt, subft, mzb, mz2_id)
    lb, ub = mzb
    dt.range[mz2_id] = (min(dt.range[mz2_id][1], lb), max(dt.range[mz2_id][2], ub))
    last_id = maximum(dt.table.id)
    origin = @p dt.table filter(≡(_.mz2_id, mz2_id))
    setdiff!(dt.table, origin)
    subft = Table(copy(subft);
        mz2_id = repeat([mz2_id], size(subft, 1)),
        mz2 = nothing,
        isf = zeros(Int, size(subft, 1)),
        n = ones(Int, size(subft, 1))
    )
    for data in subft
        id = findfirst(id -> abs(origin.rt[id] - data.rt) <= rt_tol && abs(origin.mz1[id] - data.mz1) <= mz_tol, eachindex(origin))
        if isnothing(id)
            push!(origin, data)
            last_id += 1
            origin.id[end] = last_id
        else
            # record # to recalculate # use config?
            #ft = (; subft[ft_id]..., id = origin.id[id], )
            origin[id] = merge_nt(origin[id], data)
        end
    end
    append!(dt.table, origin)
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
"""
    preis!(
            project::Project,
            featuretable::Table,
            mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
            polarity::Bool;
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
            data_id = -1,
            appedix = Dictionary{Symbol, Any}()
        )

Add precursor ion scan features.

# Arguments
* `project`: `Project`.
* `featuretable`: `Table`.
* `mz_range`: `Tuple{Float64, Float64}` or `Vector{Tuple{Float64, Float64}}`. The m/z range of PreIS is represented as `(low, up)`. When mutiple ranges are provided, the number should match the number of `mz2` in `featuretable`.
* `polarity`: `Bool`; `true` for positive ion mode; `false` for negative ion mode.

# Keyword Arguments
* `rt_tol`: Maximum allowable difference in retention time.
* `mz_tol`: Tolerance of m/z deviation or maximum allowable difference in m/z.
* `err_tol`: Tolerance of signal deviation when contructing this object.
* `est_fn`: Function for calculating signal estimate.
* `err_fn`: Function for calculating error.
* `n`: Minimal number of detection.
* `other_fn`: `Dictionary`. Stores functions for calculating other signal information.
* `data_id`: `0` indicates adding a new `PreIS`; `-n` indicates merging to last `n` data; `n` indicates merging to `n`th data.
* `appendix`: `Dictionary` containg additional information.
    * `:anion`: the anion used in mobile phase (`:acetate`, `:formate`).
    * `:signal`: `:area` or `:height` for concentration estimation.
"""
function preis!(
                    project::Project,
                    featuretable::Table,
                    mz_range::Union{Tuple{<: Number, <: Number}, Vector{<: Tuple{<: Number, <: Number}}},
                    polarity::Bool;
                    rt_tol = 0.1, 
                    mz_tol = 0.35, 
                    n = 3, 
                    signal = :area,
                    est_fn = mean, 
                    err_fn = default_error, 
                    err_tol = 0.5,
                    other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]),
                    data_id = -1,
                    appedix = Dictionary{Symbol, Any}()
                )

    db = SPDB[polarity ? :LIBRARY_POS : :LIBRARY_NEG]
    # check, update appendix and config
    featuretable = filter_combine_features!(featuretable; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
    config = dictionay(pairs((; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    source, ft = preprocess_preis!(project.data, data_id, featuretable, mz_range, polarity, config)
    dt = project.data[source]
    for analyte in project
        analyte.rt = calc_rt(analyte)
    end

    for subft in groupview(getproperty(:mz2_id), ft)
        products = id_product(dt.mz2[subft.mz2_id[1]], polarity; mz_tol)
        @info string("PreIS | ", products)
        for data in subft
            cpdsvanilla = mapreduce(vcat, eachindex(db)) do db_id
                abs(db.mz[db_id] - data.mz1) > mz_tol ? CompoundSPVanilla[] :
                    generate_cpd(db, db_id, products)
            end
            isempty(cpdsvanilla) && continue
            cpds = convert.(CompoundSP, cpdsvanilla)
            for cpd in cpds
                cpd.project = project
                cpd.fragment.id .= data.id
                cpd.fragment.source .= source
                cpd.signal = (getproperty(data, project.appendix[:signal]), data.error)
            end
            # Unidentified
            subanalytes = @views [analyte for analyte in project if abs(analyte.rt - data.rt) <= rt_tol]
            isempty(subanalytes) && (foreach(cpd -> push!(project, AnalyteSP([cpd], data.rt)), cpds); continue)
            # Aggregates
            foreach(cpd -> agg_cpd!(cpd, subanalytes, data.rt), cpds)
        end
    end
    project
end
"""
    finish_profile!(project::Project; rt_tol = 0.1, err_tol = 0.3)

Sort, merge, split, and delete analytes after all PreIS data are added.

* `rt_tol`: maximum allowable difference in retention time for two compounds to consider them as the same.
* `err_tol`: maximum allowable signal error.
"""
function finish_profile!(project::Project; rt_tol = 0.1, err_tol = 0.3)
    set!(project.appendix, :rt_tol, rt_tol)
    set!(project.appendix, :err_tol, err_tol)
    @info "PreIS | Sorting compounds"
    for analyte in project
        sort!(analyte, lt = isless_class)
        for cpd in analyte
            sort!(cpd.fragment, :ion1; lt = isless_ion)
        end
    end
    @info "PreIS | Merging, splitting and deleting analytes"
    del = Int[]
    for (i, analyte) in enumerate(project)
        signal = last(analyte).signal
        todel = false
        for a in @view project[setdiff(eachindex(project), del, i)]
            abs(a.rt - analyte.rt) > rt_tol && continue
            id = findfirst(cpd -> iscompatible(cpd, last(analyte)), a)
            isnothing(id) && continue
            any(cpd.signal[1] > signal[1] && cpd.signal[2] <= err_tol for cpd in @view a[id:end]) || continue
            todel = true
            union!(a, analyte)
            a.rt = calc_rt(a)
        end

        todel && (push!(del, i); continue)

        for (i, cpd) in Iterators.reverse(enumerate(analyte))
            cpd.signal[1] < signal[1] && continue
            cpd.signal[2] > err_tol && continue
            any(iscompatible(last(a), cpd) for a in project if abs(a.rt - analyte.rt) <= rt_tol) && continue
            push!(project, AnalyteSP(copy_wo_project.(analyte[1:i]), analyte.rt))
            last(project).rt = calc_rt(last(project))
            break
        end
        analyte.state[state_id(:error)] = signal[2] > err_tol ? -1 : 1
        analyte.state[state_id(:isf)] = any(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), last(analyte).fragment.ion1) ? 1 : -1
    end
    unique!(del)
    deleteat!(project, del)
    del = Int[]
    for (i, analyte) in enumerate(project)
        analyte.state[state_id(:isf)] < 0 || continue
        todel = any(@view project[setdiff(eachindex(project), i)]) do a
            iscompatible(last(analyte), last(a))
        end
        todel && push!(del, i)
    end
    unique!(del)
    deleteat!(project, del)
    @info "PreIS | Assigning ISF"
    assign_isf_parent!(project)
    for analyte in project
        cpd = last(analyte)
        id = findall(ion -> in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragment.ion1)
        v = query_data.(Ref(project), cpd.fragment.source[id], cpd.fragment.id[id], :isf)
        analyte.state[state_id(:isf)] = analyte.state[state_id(:isf)] == -1 ? -1 : any(==(1), v) ? 1 : any(==(-1), v) ? -1 : 0
        for cpd in analyte
            cpd.project = project
        end
    end
    @info "PreIS | Sorting analytes"
    sort!(project; by = x -> (mw(x), rt(x)))
    project
end

function assign_isf_parent!(project::Project)
    aq = q!(project, qnot(:isf!))
    for analyte in aq
        for cpd in @view analyte[1:end - 1]
            set_data!.(Ref(project), cpd.fragment.source, cpd.fragment.id, :isf, -1)
        end
    end
    ids = Dict{Int, Vector{Int}}()
    for analyte in aq
        cpd = last(analyte)
        id = findall(ion -> !in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragment.ion1)
        set_data!.(Ref(project), cpd.fragment.source[id], cpd.fragment.id[id], :isf, -1)
        id = setdiff(eachindex(cpd.fragment), id)
        for i in id
            push!(get!(ids, cpd.fragment.source[i], Int[]), cpd.fragment.id[i])
        end
    end
    for (k, v) in pairs(ids)
        replace_data!.(Ref(project), k, unique!(v), :isf, 0 => 1, -1 => 0)
    end
end

assign_parent!(project::Project) =
    for analyte in project
        cpd = last(analyte)
        id = findall(ion -> !in(ion.adduct, class_db_index(ion.molecule).parent_ion), cpd.fragment.ion1)
        set_data!.(Ref(project), cpd.fragment.source[id], cpd.fragment.id[id], :isf, -1)
    end
