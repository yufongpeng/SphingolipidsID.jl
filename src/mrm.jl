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

* `raw`: whether use the raw feature table or feature table processed by `filter_combine_features!`.
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
                other_fn = Dictionary{Symbol, Any}())
    default_other_fn = Dictionary{Symbol, Any}([:datafile], [nothing])
    if !isempty(other_fn)
        for (k, v) in pairs(other_fn)
            set!(default_other_fn, k, v)
        end
    end
    other_fn = default_other_fn
    featuretable = raw ? featuretable : filter_combine_features!(featuretable; rt_tol, mz_tol, n, signal, est_fn, err_fn, err_tol, other_fn)
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