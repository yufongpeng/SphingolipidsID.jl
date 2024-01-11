"""
    MRM(featuretable::Table; 
            raw = false,
            filter = true,
            combine = true,
            use_match_id = false,
            rt_tol = 0.1, 
            mz_tol = 0.35, 
            n = 3, 
            signal = :area,
            est_fn = mean, 
            err_fn = default_error, 
            err_tol = 0.5,
            other_fn = Dictionary{Symbol, Any}([:datafile], [nothing]))

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
            combine = true,
            use_match_id = false,
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
            combine = true,
            use_match_id = false,
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
                combine = true,
                use_match_id = false,
                rt_tol = 0.1, 
                mz_tol = 0.35, 
                n = 3, 
                signal = :area,
                est_fn = mean, 
                err_fn = default_error, 
                err_tol = 0.5,
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile], [nothing])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
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