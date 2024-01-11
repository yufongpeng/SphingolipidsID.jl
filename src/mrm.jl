"""
    MRM(featuretable::Table; 
            raw = false,
            filter = true,
            combine = false,
            use_match_id = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 1, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}())

Create `MRM` from a feature table.

* `raw`: whether use the raw feature table or feature table processed by `filter_combine_features!`.
* `filter`: whether filter data based on number of detection and signal errors.
* `combine`: whether combine data after grouping based on `mz1`, `mz2`, `rt`, and `collision_energy`. `datafile` is ignored.
* `use_match_id`: whether using `match_id` instead of `mz1`, `mz2`, `rt` to group data.
* `raw_id`: a `Bool` determining whether preserving the original `tbl.id` as `tbl.raw_id` for staying connection to the raw featuretable.
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
            filter = true,
            combine = false,
            use_match_id = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 1, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}()
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    MRM!(deepcopy(Table(deepcopy(featuretable); polarity = nothing)), polarity; 
        raw,                            
        filter,
        combine,
        use_match_id,
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
            filter = true,
            combine = false,
            use_match_id = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 1, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}()
    )
    polarity = all(featuretable.polarity) ? true : any(featuretable.polarity) ? throw(ArgumentError(
        "featuretable should contain only positive or negative data."
    )) : false
    MRM!(Table(featuretable; polarity = nothing), polarity; 
        raw,
        filter,
        combine,
        use_match_id,
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
                filter = true,
                combine = false,
                use_match_id = false,
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 1, 
                signal = :area,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    featuretable = raw ? featuretable : filter_combine_features!(featuretable; 
        filter,
        combine,
        use_match_id, 
        rt_tol, 
        mz_tol, 
        n, 
        signal, 
        est_fn, 
        err_fn, 
        err_tol, 
        other_fn
    )
    config = dictionary(pairs((; raw, filter, combine, use_match_id, rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)))
    mz2_loc = @p featuretable groupfind(x -> round.(getproperty(x, :mz2); digits = -floor(Int, log10(mz_tol)))) collect
    mz2s = map(i -> mean(featuretable.mz2[i]), mz2_loc)
    MRM(Table(copy(featuretable);
                mz2_id = (@p eachindex(featuretable) |> map(findfirst(x -> in(_, x), mz2_loc))),
                mz2 = nothing), mz2s, polarity, config)
end

function qtMRM!(featuretable::Table; rt_tol = 0.1, mz_tol = 0.35, signal = :area)
    mrm = MRM!(Table(featuretable; 
            match_id = hasproperty(featuretable, :match_id) ? featuretable.match_id : zeros(Int, length(featuretable)), 
            match_score = hasproperty(featuretable, :match_score) ? featuretable.match_score : -getproperty(featuretable, signal)
            ); 
        rt_tol, mz_tol
    )
    mrm.table.match_id .= mrm.table.id
    mrm
end

function qtMRM(featuretable::Table; rt_tol = 0.1, mz_tol = 0.35, signal = :area)
    mrm = MRM(Table(featuretable; 
            match_id = hasproperty(featuretable, :match_id) ? featuretable.match_id : zeros(Int, length(featuretable)), 
            match_score = hasproperty(featuretable, :match_score) ? featuretable.match_score : -getproperty(featuretable, signal)
            ); 
        rt_tol, mz_tol
    )
    mrm.table.match_id .= mrm.table.id
    mrm
end