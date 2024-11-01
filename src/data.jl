"""
    sort_data(tbl::Table)

Create a new `Table` and sort it by `:injection_order`.
"""
sort_data(tbl::Table) = sort_data!(deepcopy(tbl))
"""
    sort_data!(tbl::Table)

Sort the table by `:injection_order` in place.
"""
function sort_data!(tbl::Table)
    sort!(tbl, :injection_order)
    tbl
end

"""
    file_order(tbl::Table)

Return unique `tbl.datafile` for filling CE, mz2 or m/z range by the `injection_order`.
"""
file_order(tbl::Table) = !in(:datafile, propertynames(tbl)) ? throw(ArgumentError("This table doesn't contain `datafile`.")) : 
    !in(:datafile, propertynames(tbl)) ? throw(ArgumentError("This table doesn't contain `injection_order`.")) : 
    map(x -> tbl.datafile[findfirst(==(x), tbl.injection_order)], sort!(unique(tbl.injection_order)))

"""
    fill_polarity!(tbl::Table, polarity::Union{<: AbstractVector, <: Tuple})
    fill_polarity!(tbl::Table, polarity::Union{<: AbstractDictionary, <: AbstractDict})
    fill_polarity!(tbl::Table, polarity::Bool)

Fill column `polarity` with `polarity`.

When the input is a number, it fills all values with the number.
When the input is a `Vector` or `Tuple`, the values should represent polarity of each `tbl.injection_order`.
When the input is a Dictionary, the keys should contain each `tbl.datafile`.
"""
function fill_polarity!(tbl::Table, polarity)
    tbl.polarity .= fill_data(tbl, polarity)
    tbl
end

function fill_polarity!(tbl::Table, polarity::Bool)
    fill!(tbl.polarity, polarity)
    tbl
end

"""
    fill_polarity(tbl::Table, polarity::Union{<: AbstractVector, <: Tuple})
    fill_polarity(tbl::Table, polarity::Union{<: AbstractDictionary, <: AbstractDict})
    fill_polarity(tbl::Table, polarity::Bool)

Create a table with column `polarity` and fill it with `polarity`.

See `fill_polarity!`.
"""
fill_polarity(tbl::Table, polarity) = Table(tbl; polarity = fill_data(tbl, polarity))

"""
    fill_mz2!(tbl::Table, mz2::Union{<: AbstractVector, <: Tuple})
    fill_mz2!(tbl::Table, mz2::Union{<: AbstractDictionary, <: AbstractDict})
    fill_mz2!(tbl::Table, mz2::Float64)

Fill column `mz2` with `mz2`.

When the input is a number, it fills all values with the number.
When the input is a `Vector` or `Tuple`, the values should represent mz2 of each `tbl.injection_order`.
When the input is a Dictionary, the keys should contain each `tbl.datafile`.
"""
function fill_mz2!(tbl::Table, mz2)
    tbl.mz2 .= fill_data(tbl, mz2)
    tbl
end

function fill_mz2!(tbl::Table, mz2::Float64)
    fill!(tbl.mz2, mz2)
    tbl
end

"""
    fill_mz2(tbl::Table, mz2::Union{<: AbstractVector, <: Tuple})
    fill_mz2(tbl::Table, mz2::Union{<: AbstractDictionary, <: AbstractDict})
    fill_mz2(tbl::Table, mz2::Float64)
    
Create a table with column `mz2` and fill it with `mz2`.

See `fill_mz2!`.
"""
fill_mz2(tbl::Table, mz2) = Table(tbl; mz2 = fill_data(tbl, mz2))

"""
    fill_ce!(tbl::Table, ce::Union{<: AbstractVector, <: Tuple})
    fill_ce!(tbl::Table, ce::Union{<: AbstractDictionary, <: AbstractDict})
    fill_ce!(tbl::Table, ce::Float64)

Fill column `collision_energy` with `ce`.

When the input is a number, it fills all values with the number.
When the input is a `Vector` or `Tuple`, the values should represent collision energy of each `tbl.injection_order`.
When the input is a Dictionary, the keys should contain each `tbl.datafile`.
"""
function fill_ce!(tbl::Table, ce)
    tbl.collision_energy .= fill_data(tbl, ce)
    tbl
end

function fill_ce!(tbl::Table, ce::Float64)
    fill!(tbl.collision_energy, ce)
    tbl
end

"""
    fill_ce(tbl::Table, ce::Union{<: AbstractVector, <: Tuple})
    fill_ce(tbl::Table, ce::Union{<: AbstractDictionary, <: AbstractDict})
    fill_ce(tbl::Table, ce::Float64)

Create a table with column `collision_energy` and fill it with `ce`.

See `fill_ce!`.
"""
fill_ce(tbl::Table, ce) = Table(tbl; collision_energy = fill_data(tbl, ce))

"""
    fill_range!(tbl::Table, mzr::Union{<: AbstractVector, <: Tuple})
    fill_range!(tbl::Table, mzr::Union{<: AbstractDictionary, <: AbstractDict})
    fill_range!(tbl::Table, mzr::RealIntervals)

Fill column `mz_range` with `mzr`.

When the input is a interval, it fills all values with the interval.
When the input is a `Vector` or `Tuple`, the values should represent m/z range of each `tbl.injection_order`.
When the input is a Dictionary, the keys should contain each `tbl.datafile`.
"""
function fill_range!(tbl::Table, mzr)
    tbl.mz_range .= fill_data(tbl, mzr)
    tbl
end

function fill_range!(tbl::Table, mzr::RealIntervals)
    fill!(tbl.mz_range, mzr)
    tbl
end

"""
    fill_range(tbl::Table, mzr::Union{<: AbstractVector, <: Tuple})
    fill_range(tbl::Table, mzr::Union{<: AbstractDictionary, <: AbstractDict})
    fill_range(tbl::Table, mzr::RealIntervals)

Create a table with column `mz_range` and fill it with `mzr`.

See `fill_range!`.
"""
fill_range(tbl::Table, mzr) = Table(tbl; mz_range = fill_data(tbl, mzr))

fill_data(tbl::Table, data::Union{<: AbstractVector, <: Tuple}) = getindex.(Ref(data), tbl.injection_order)
fill_data(tbl::Table, data::Union{<: AbstractDictionary, <: AbstractDict}) = getindex.(Ref(data), tbl.datafile)
fill_data(tbl::Table, data) = repeat([data], length(tbl))

"""
    rsd(v)

Relative standard deviation, i.e., `std(v) / mean(v)`.
"""
rsd(v) = std(v) / mean(v)
"""
    rmae(v)

Relative mean absolute error.
"""
function rmae(v) 
    m = mean(v)
    mean(abs, v .- m) / m
end
default_error(v) = length(v) > 2 ? rsd(v) : rmae(v)

"""
    split_datafile(tbl::Table; name = r"(.*)-r\\d*\\.(?:d|mzML)", rep = r".*-r(\\d*)\\.(?:d|mzML)")
Split the feature table accoding to data file name. `name` specifies the name pattern, and the resulting name will be the keys of the dictionary. `rep` specifies the repeating pattern, and will be furhter parsed as integer. If `rep` is nothing, the column `:datafile` will be removed; otherwise, the resulting number will be stored in this column.
"""
function split_datafile(tbl::Table; name = r"(.*)-r\d*\.(?:d|mzML)", rep = r".*-r(\d*)\.(?:d|mzML)")
    datafile = (@p tbl.datafile map(match(name, _)[1]))
    id = findall(!isnothing, datafile)
    tbl = tbl[id]
    tbl = Table(tbl; datafile = isnothing(rep) ? nothing : (@p tbl.datafile map(parse(Int, match(rep, _)[1]))))
    group(datafile[id], tbl)
end

"""
    filter_combine_features!(tbl::Table; 
                            filter = true,
                            combine = true,
                            use_match_id = false,
                            raw_id = false,
                            rt_tol = 0.1, 
                            mz_tol = 0.35, 
                            n = 3, 
                            err_fn = default_error, 
                            err_tol = 0.5,
                            other_fn = Dictionary{Symbol, Any}([:mz_range, :polarity, :datafile, :injection_order], [splat(union), only ∘ unique, nothing, nothing]))

Filter out and combine duplicated or unstable features. Each `tbl.id` will map to a unique feature.

* `filter`: whether filter data based on number of detection and signal errors. Only works when `signal` is not `nothing`.
* `combine`: whether combine data after grouping based on `mz1`, `mz2`, `rt`, and `collision_energy`.
* `use_match_id`: whether using `match_id` instead of `mz1`, `mz2`, `rt` to group data.
* `raw_id`: a `Bool` determining whether preserving the original `tbl.id` as `tbl.raw_id` for staying connection to the raw featuretable. Only works when `combine` is true.
* `rt_tol`: maximum allowable difference in retention time for two compounds to consider them as the same.
* `mz_tol`: maximum allowable difference in m/z for two compounds to consider them as the same.
* `n`: minimal number of detection. The default is 3.
* `signal`: `:area` or `:height` for concentration estimation. The default is `:area`. If `nothing` is given, this function will not calculate errors.
* `est_fn`: estimation function for calculating estimated value. The default is `mean`.
* `err_fn`: error function. If the length is three or above, the default is `rsd`; otherwise, it is `rmae`. 
* `err_tol`: maximum allowable signal error.
* `other_fn`: a `Dictionary` specifying functions for combining other signal information. The keys are column names, and the values are functions. If the value is `nothing`, the column will be removed.
For columns not mentioned, `est_fn` will be used for numeric value, and `first` will be used for other types. 
"""
function filter_combine_features!(tbl::Table; 
                            filter = true, 
                            combine = true,
                            use_match_id = false,
                            raw_id = false,
                            rt_tol = 0.1, 
                            mz_tol = 0.35, 
                            n = 3, 
                            signal = :area,
                            est_fn = mean, 
                            err_fn = default_error, 
                            err_tol = 0.5,
                            other_fn = Dictionary{Symbol, Any}([:mz_range, :polarity, :datafile, :injection_order], [splat(union), only ∘ unique, nothing, nothing]))
    sort!(tbl, [:mz2, :mz1, :rt])
    calc_signal = !isnothing(signal)
    filter = calc_signal && filter
    if use_match_id 
        locs = collect(groupfind(getproperties((:match_id, :collision_energy)), tbl)) 
    else
        gids = groupfind(x -> round.(getproperty(x, :mz2); digits = -floor(Int, log10(mz_tol))), tbl)
        locs = mapreduce(∪, collect(gids); init = Vector{Int}[]) do ids
            locs = Vector{Int}[]
            for i in ids
                rrt = tbl.rt[i]
                rmz1 = tbl.mz1[i]
                rmz2 = tbl.mz2[i]
                rcollision_energy = tbl.collision_energy[i]
                new = true
                for loc in locs
                    (!between(mean(tbl.rt[loc]), rrt, rt_tol) || 
                    !between(mean(tbl.mz1[loc]), rmz1, mz_tol) ||
                    tbl.collision_energy[loc[1]] != rcollision_energy ||
                    !between(mean(tbl.mz2[loc]), rmz2, mz_tol)) && continue
                    push!(loc, i)
                    new = false
                    break
                end
                new && push!(locs, [i])
            end
            locs
        end
    end
    filter && n > 1 && filter!(loc -> (length(loc) >= n && err_fn(getproperty(tbl, signal)[loc]) <= err_tol), locs)
    if raw_id && combine
        ids = map(loc -> tbl.id[loc], locs)
    end
    fill!(tbl.id, 0)
    for (i, loc) in enumerate(locs)
        tbl.id[loc] .= i
    end
    filter!(x -> >(x.id, 0), tbl)
    sort!(tbl, :id)
    combine || return tbl
    del = [k for (k, v) in pairs(other_fn) if isnothing(v)]
    tbl = Table(tbl; (del .=> nothing)...)
    calc_signal && set!(other_fn, signal, est_fn)
    get!(other_fn, :id, only ∘ unique)
    for k in propertynames(tbl)
        eltype(getproperty(tbl, k)) isa Number ? get!(other_fn, k, est_fn) : get!(other_fn, k, first)
    end
    gtbl = @p tbl groupview(getproperty(:id))
    calc_signal && (errors = @p gtbl map(err_fn(getproperty(_, signal))) collect)
    nt = @p gtbl map(map((k, v) -> k => other_fn[k](v), propertynames(_), columns(_))) map(NamedTuple) 
    if calc_signal && raw_id
        Table(nt; error = errors, raw_id = ids)
    elseif calc_signal
        Table(nt; error = errors)
    elseif raw_id
        Table(nt; raw_id = ids)
    else
        Table(nt)
    end
end

"""
    compoundspvanilla(cpd::NamedTuple, product::Ion)

Create a `CompoundSPVanilla` by a row of database(matched by mz1) and product. 

If there are conflicts, it returns `nothing`.
"""
function compoundspvanilla(cpd::NamedTuple, product::Ion)
    cls = (cpd.Abbreviation)()
    sc = sumcomp(cpd.Species)
    adduct_class = object_adduct(cpd.Adduct)
    pr = process_product(product, cls, sc)
    isnothing(pr) && return pr
    ion2, sc = pr
    CompoundSPVanilla(cls, sc,
        Table(ion1 = Ion[Ion(adduct_class, cls)], ion2 = Ion[ion2], source = [1], id = [1])
    )
end

"""
    process_product(product::Ion, cls::ClassSP, sc::SumChain)

Infer the structure of side chain and adjust `product` if neccessary.

This function returns `(product, sidechain)` if there are no conflicts; otherwise, it returns `nothing`.
"""
process_product(product::Ion, cls::ClassSP, sc::Union{SumChain, SumChainIS}) =
    (in(product, NANA) && !hasnana(cls)) ? nothing : (product, sc)
function process_product(product::Ion{<: Adduct, <: LCB}, cls::ClassSP, sc::SumChain)
    lcb_o = nox(product.molecule)
    lcb_db = ndb(product.molecule)
    sum_db = ndb(sc)
    sum_o = nox(sc)
    #=
    1. sum_db + sum_o < lcb_db + lcb_o (total unsaturation deficiency)
    2. sum_o < lcb_o
        0 < Δ = lcb_o - sum_o ≤ sum_db - lcb_db
        lcb_o' = lcb_o - Δ = sum_o => sum_o - lcb_o' = 0
        lcb_db' = lcb_db + Δ ≤ sum_db
    3. sum_db < lcb_db
        0 < Δ = lcb_db - sum_db ≤ sum_o - lcb_o
        lcb_db' = lcb_db - Δ = sum_db
        lcb_o' = lcb_o + Δ ≤ sum_o => sum_o - lcb_o' = sum_o - lcb_o - Δ ≥ 0
    =#
    (sum_db + sum_o) < (lcb_o + lcb_db) && return nothing
    Δ = sum_o - lcb_o
    Δ = Δ < 0 ? Δ : max(0, lcb_db - sum_db)
    product_new = @match Δ begin
        0 => product
        _ => begin
            lcb_new = hydroxyl_shift(product.molecule, Δ)
            add_new = hydroxyl_shift(product.adduct, Δ)
            (isnothing(lcb_new) || isnothing(add_new)) && return nothing
            Ion(add_new, lcb_new)
        end
    end
    lcb_new = product_new.molecule
    product_new, DiChain(lcb_new, Acyl(ncb(sc) - ncb(lcb_new), sum_db - ndb(lcb_new), sum_o - nox(lcb_new)))
end

function process_product(product::Ion{<: Adduct, <: LCB}, cls::ClassSP, sc::SumChainIS)
    lcb_o = nox(product.molecule)
    lcb_db = ndb(product.molecule)
    lcb_13C = n13C(product.molecule)
    lcb_D = nD(product.molecule)
    sum_db = ndb(sc)
    sum_o = nox(sc)
    sum_13C = n13C(sc)
    sum_D = nD(sc)
    #=
    1. sum_db + sum_o < lcb_db + lcb_o (total unsaturation deficiency)
    2. sum_o < lcb_o
        0 < Δ = lcb_o - sum_o ≤ sum_db - lcb_db
        lcb_o' = lcb_o - Δ = sum_o => sum_o - lcb_o' = 0
        lcb_db' = lcb_db + Δ ≤ sum_db
    3. sum_db < lcb_db
        0 < Δ = lcb_db - sum_db ≤ sum_o - lcb_o
        lcb_db' = lcb_db - Δ = sum_db
        lcb_o' = lcb_o + Δ ≤ sum_o => sum_o - lcb_o' = sum_o - lcb_o - Δ ≥ 0
    =#
    (sum_db + sum_o) < (lcb_o + lcb_db) && return nothing
    sum_13C < lcb_13C && return nothing
    sum_D < lcb_D && return nothing
    Δ = sum_o - lcb_o
    Δ = Δ < 0 ? Δ : max(0, lcb_db - sum_db)
    product_new = @match Δ begin
        0 => product
        _ => begin
            lcb_new = hydroxyl_shift(product.molecule, Δ)
            add_new = hydroxyl_shift(product.adduct, Δ)
            (isnothing(lcb_new) || isnothing(add_new)) && return nothing
            Ion(add_new, lcb_new)
        end
    end
    lcb_new = product_new.molecule
    acyl_13C = sum_13C - n13C(lcb_new)
    acyl_D = sum_D - nD(lcb_new)
    (product_new, acyl_13C == 0 && acyl_D == 0 ? DiChain(lcb_new, Acyl(ncb(sc) - ncb(lcb_new), sum_db - ndb(lcb_new), sum_o - nox(lcb_new))) : 
        DiChain(lcb_new, AcylIS(ncb(sc) - ncb(lcb_new), sum_db - ndb(lcb_new), sum_o - nox(lcb_new), (n13C = acyl_13C, nD = acyl_D))))
end

"""
    id_product(ms2::Float64, polarity::Bool; mz_tol = 0.35)

Identify possible products based on given ms2 value.

Return `Vector{Ion}`.
"""
function id_product(ms2::Float64, polarity::Bool; mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG])
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end

"""
    union_transition(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)

Take the union of two MRM transition tables and create a new `Table`. MRM data will be merged into one data if they are considered close enough.

The tables should containing the following columns:
* `analyte`: `String`; possible analyte name separated by "|".
* `mz1`: m/z of parent ion.
* `mz2`: m/z of fragment ion.
* `rt`: retention time.
* `Δrt`: the difference between the maximum and minimum retention time.
* `collision_energy`: collision energy.
* `polarity`: `Bool`; `true` for positive ion,; `false` for negative ion.
"""
union_transition(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = union_transition!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
"""
    union_transition!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)

Take the union of two MRM transition tables and mutate `tbl1` in place. See `union_transition`.
"""
function union_transition!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    for data2 in tbl2
        pushed = false
        for (id1, data1) in enumerate(tbl1)
            (!between(data1.mz1, data2.mz1, mz_tol) || !between(data1.mz2, data2.mz2, mz_tol)) && continue
            data1.collision_energy == data2.collision_energy || continue
            rt_l = min(data1.rt - data1.Δrt / 2, data2.rt - data2.Δrt / 2)
            rt_r = max(data1.rt + data1.Δrt / 2, data2.rt + data2.Δrt / 2)
            analytes1 = split(data1.analyte, "|")
            analytes2 = split(data2.analyte, "|")
            n = length(analytes1)
            tbl1[id1] = (;
                analyte = join(union!(analytes1, analytes2), "|"),
                mz1 = (data1.mz1 * n + data2.mz1) / (n + 1),
                mz2 = (data1.mz2 * n + data2.mz2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = data1.collision_energy,
                polarity = data1.polarity
            )
            pushed = true
        end
        pushed || push!(tbl1, (analyte = data2.analyte, mz1 = data2.mz1, mz2 = data2.mz2, rt = data2.rt, Δrt = rt_tol * 2, collision_energy = data2.collision_energy, polarity = data2.polarity))
    end
    sort!(tbl1, [:mz1, :rt])
end

"""
    diff_transition(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)

Take the set difference of two MRM transition tables. MRM data will be merged into one data if they are considered close enough.

This function return a `NamedTuple`:
* `extended`: data from `tbl1` that has been modified with data from `tbl2`.
* `new`: data from `tbl1` that do not exist in `tbl2`.

The tables should containing the following columns:
* `analyte`: `String`; possible analyte name separated by "|".
* `mz1`: m/z of parent ion.
* `mz2`: m/z of fragment ion.
* `rt`: retention time.
* `Δrt`: the difference between the maximum and minimum retention time.
* `collision_energy`: collision energy.
* `polarity`: `Bool`; `true` for positive ion,; `false` for negative ion.
"""
diff_transition(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = diff_transition!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
"""
    diff_transition!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)

Take the set difference of two MRM transition tables and mutate `tbl1` in place. See `diff_transition`.
"""
function diff_transition!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    extended = Int[]
    l = size(tbl1, 1)
    for data2 in tbl2
        pushed = false
        for (id1, data1) in enumerate(tbl1)
            (!between(data1.mz1, data2.mz1, mz_tol) || !between(data1.mz2, data2.mz2, mz_tol)) && continue
            data1.collision_energy == data2.collision_energy || continue
            rt_l = min(data1.rt - data1.Δrt / 2, data2.rt - data2.Δrt / 2)
            rt_r = max(data1.rt + data1.Δrt / 2, data2.rt + data2.Δrt / 2)
            analytes1 = split(data1.analyte, "|")
            analytes2 = split(data2.analyte, "|")
            n = length(analytes1)
            isempty(setdiff(analytes2, analytes1)) || push!(extended, id1)
            tbl1[id1] = (;
                analyte = join(union!(analytes1, analytes2), "|"),
                mz1 = (data1.mz1 * n + data2.mz1) / (n + 1),
                mz2 = (data1.mz2 * n + data2.mz2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = data1.collision_energy,
                polarity = data1.polarity
            )
            pushed = true
        end
        pushed && continue
        #push!(new, id2)
        push!(tbl1, (analyte = data2.analyte, mz1 = data2.mz1, mz2 = data2.mz2, rt = data2.rt, Δrt = rt_tol * 2, collision_energy = data2.collision_energy, polarity = data2.polarity))
    end
    (extended = tbl1[extended], new = size(tbl1, 1) > l ? tbl1[l + 1:end] : nothing)
end