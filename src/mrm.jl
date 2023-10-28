"""
    MRM(featuretable::Table; 
            raw = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]))

Create `MRM` from a feature table.

* `raw`: whether use the raw feature table or feature table processed by `filter_duplicate!`.
* `rt_tol`: maximum allowable difference in retention time for two compounds to consider them as the same.
* `mz_tol`: maximum allowable difference in m/z for two compounds to consider them as the same.
* `n`: minimal number of detection. The default is 3.
* `signal`: `:area` or `:height` for concentration estimation. The default is `:area`.
* `est_fn`: estimation function for calculating estimated value. The default is `mean`.
* `err_fn`: error function. If the length is three or above, the default is `rsd`; otherwise, it is `re`. 
* `err_tol`: maximum allowable signal error.
* `other_fn`: a `Dictionary` specifying functions for calculating other signal information.
"""
function MRM(featuretable::Table; 
            raw = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing])
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    MRM!(deepcopy(Table(featuretable; polarity = nothing)), polarity; 
        raw,
        rt_tol, 
        mz_tol, 
        n, 
        signal,
        est_fn, 
        err_fn, 
        err_tol,
        other_fn)
end

function MRM!(featuretable::Table; 
            raw = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing])
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    MRM!(Table(featuretable; polarity = nothing), polarity; 
        raw,
        rt_tol, 
        mz_tol, 
        n, 
        signal,
        est_fn, 
        err_fn, 
        err_tol,
        other_fn)
end

function MRM!(featuretable::Table, polarity::Bool; 
                raw = false,
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 3, 
                signal = :area,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]))
    featuretable = raw ? featuretable : filter_duplicate!(featuretable; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
    config = dictionary(pairs((; raw, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    #=
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
    tables = Table[]
    for (i, loc) in enumerate(mz2_loc)
        subft = @view featuretable[loc]
        push!(tables, Table(copy(subft), mz2_id = repeat([i], size(subft, 1)), mz2 = nothing))
    end
    =#
    mz2_loc = @p featuretable groupfind(x -> round.(getproperty(x, :mz2); digits = -floor(Int, log10(mz_tol)))) collect
    mz2s = map(i -> mean(featuretable.mz2[i]), mz2_loc)
    MRM(Table(copy(featuretable);
                mz2_id = (@p eachindex(featuretable) |> map(findfirst(x -> in(_, x), mz2_loc))),
                mz2 = nothing), mz2s, polarity, config)
end

function qcdata_mrm(featuretable::Table; name = r"pooledqc.*_(\d*).*",
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                signal = :area,
                est_fn = mean, 
                err_fn = rsd, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:repeat, :polarity], [nothing, nothing]))
    rawdata = filter(x -> occursin(name, x.datafile), featuretable)
    rep = map(x -> parse(Int, match(name, x)[1]), rawdata.datafile)
    n = length(unique(rep))
    rawdata = Table(rawdata; datafile = nothing, repeat = rep)
    mrm = MRM!(rawdata; raw = true, rt_tol, mz_tol, n = 1, signal, est_fn, err_fn, err_tol, other_fn = Dictionary{Symbol, Any}())
    # match project for IS
    featuretable = filter_duplicate!(deepcopy(rawdata); raw_id = true, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
    # add matched analyte
    config = dictionary(pairs((; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    QCData(mrm, featuretable, config)
end

function serialdilution_mrm(featuretable::Table, level::Vector; name = r"cal.*_(\d*).*",
                r2_threshold = 0.8,
                nlevel = 5,
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                signal = :area,
                est_fn = mean, 
                err_fn = std, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}([:mz2_id, :collision_energy, :FWHM, :symmetry, :level, :polarity, :height, :area], [nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing]))
    rawdata = filter(x -> occursin(name, x.datafile), featuretable)
    lv = map(x -> parse(Int, match(name, x)[1]), rawdata.datafile)
    rawdata = Table(rawdata; datafile = nothing, level = lv)
    mrm = MRM!(rawdata; raw = true, rt_tol, mz_tol, n = 1, signal, est_fn, err_fn, err_tol, other_fn = Dictionary{Symbol, Any}())
    # match project for IS
    featuretable = filter_duplicate!(deepcopy(rawdata); raw_id = true, rt_tol, mz_tol, n = 1, signal = false, est_fn, err_fn, err_tol, other_fn)
    # add matched analyte
    config = dictionary(pairs((; r2_threshold, nlevel, rt_tol, mz_tol, n = 1, signal = false, est_fn, err_fn, err_tol, other_fn)))
    # r2, model, data_id
    r2s = zeros(Float64, length(featuretable))
    data_ids = Vector{Vector{Int}}(undef, length(featuretable))
    Threads.@threads for i in eachindex(featuretable)
        v = getproperty(mrm.table, mrm.config[:signal])[featuretable.raw_id[i]]
        l = mrm.table.level[featuretable.raw_id[i]]
        ord = sortperm(l)
        l = level[l[ord]]
        v = v[ord]
        id = featuretable.raw_id[i][ord]
        _, data_ids[i], r2s[i] = rec_r2(id, l, v, nlevel, r2_threshold, 0)
    end
    SerialDilution(mrm, Table(featuretable; r2 = r2s, data_id = data_ids), level, config)
end

function rec_r2(id, level, value, nlevel, r2_threshold, r2_current)
    length(unique(level)) < nlevel && return (false, Int[], -Inf)
    r2 = cor(level, value) ^ 2
    r2 < r2_current && return (false, id, r2_current)
    r2 >= r2_threshold && return (true, id, r2)
    s = findfirst(!=(first(level)), level)
    e = findfirst(==(last(level)), level) - 1
    r2_current = r2
    id_current = id
    status = false
    for i in s:e
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && return (status, id_current, r2_current)
    for i in firstindex(id):(s - 1)
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    status && (status, id_current, r2_current)
    for i in (e + 1):lastindex(id)
        idx = setdiff(eachindex(level), i)
        status, id_new, r2_new = rec_r2(id[idx], level[idx], value[idx], nlevel, r2_threshold, r2_current)
        if r2_new > r2_current
            r2_current = r2_new
            id_current = id_new
            status && break
        end
    end
    (status, id_current, r2_current)
end