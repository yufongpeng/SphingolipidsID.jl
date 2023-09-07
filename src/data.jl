"""
    sort_data(tbl::Table)

Create a new `Table` and sort it by `:datafile`.
"""
sort_data(tbl::Table) = sort_data!(deepcopy(tbl))
"""
    sort_data!(tbl::Table)

Sort the table by `:datafile` in place.
"""
function sort_data!(tbl::Table)
    sort!(tbl, :datafile)
    tbl
end

"""
    file_order(tbl::Table)

Return unique `tbl.datafile` for filling CE or mz2 by the order.
"""
file_order(tbl::Table) = in(:datafile, propertynames(tbl)) ? unique(tbl.datafile) : throw(ArgumentError("This table doesn't contain datafile."))

"""
    fill_mz2!(tbl::Table, mz2::Union{<: Vector, <: Tuple})
    fill_mz2!(tbl::Table, mz2::Dict)
    fill_mz2!(tbl::Table, mz2::Float64)

Fill column `mz2` with `mz2`.

If the input is not a number, it requires that column `datafile` exists and uses the order or values as keys to assign mz2.

When the input is a `Vector` or `Tuple`, the values should represent mz2 of each data file in the order as same as `file_order(tbl)`.
"""
function fill_mz2!(tbl::Table, mz2::Union{<: Vector, <: Tuple})
    in(:datafile, propertynames(tbl)) || throw(ArgumentError("Column `datafile` is not in the table"))
    mapping = Dict(unique(tbl.datafile) .=> mz2)
    tbl.mz2 .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

function fill_mz2!(tbl::Table, mz2::Dict)
    in(:datafile, propertynames(tbl)) || throw(ArgumentError("Column `datafile` is not in the table"))
    tbl.mz2 .= getindex.(Ref(mz2), tbl.datafile)
    tbl
end

function fill_mz2!(tbl::Table, mz2::Float64)
    fill!(tbl.mz2, mz2)
    tbl
end

"""
    fill_ce!(tbl::Table, ce::Union{<: Vector, <: Tuple})
    fill_ce!(tbl::Table, ce::Dict)
    fill_ce!(tbl::Table, ce::Float64)

Fill column `ce` with `ce`.

If the input is not a number, it requires that column `datafile` exists and uses the order or values as keys to assign ce.

When the input is a `Vector` or `Tuple`, the values should represent ce of each data file in the order as same as `file_order(tbl)`.
"""
function fill_ce!(tbl::Table, ce::Union{<: Vector, <: Tuple})
    in(:datafile, propertynames(tbl)) || throw(ArgumentError("Column `datafile` is not in the table"))
    mapping = Dict(unique(tbl.datafile) .=> ce)
    tbl.collision_energy .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

function fill_ce!(tbl::Table, ce::Dict)
    in(:datafile, propertynames(tbl)) || throw(ArgumentError("Column `datafile` is not in the table"))
    tbl.collision_energy .= getindex.(Ref(ce), tbl.datafile)
    tbl
end

function fill_ce!(tbl::Table, ce::Float64)
    fill!(tbl.collision_energy, ce)
    tbl
end

"""
    rsd(v)

Relative standard deviation, i.e., `std(v) / mean(v)`.
"""
rsd(v) = std(v) / mean(v)
"""
    re(v)

Relative error, i.e., `(maximum(v) - minimum(v)) / 2 / mean(v)`.
"""
re(v) =  - foldl(-, extrema(v)) / mean(v) / 2
default_error(v) = length(v) > 2 ? rsd(v) : re(v)

"""
    filter_duplicate!(tbl::Table; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5)

Filter out duplicated or unstable data.

* `rt_tol`: maximum allowable difference in retention time for two compounds to consider them as the same.
* `mz_tol`: maximum allowable difference in m/z for two compounds to consider them as the same.
* `err`: error function. If the length is three or above, the default is `rsd`; otherwise, it is `re`. 
* `err_tol`: maximum allowable signal error.
"""
function filter_duplicate!(tbl::Table; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5)
    tbl = Table(tbl; datafile = nothing)
    sort!(tbl, [:mz2, :mz1, :rt])
    fill!(tbl.id, 0)
    locs = Vector{Int}[]
    for (i, ft) in enumerate(tbl)
        new = true
        for loc in locs
            abs(mean(tbl.rt[loc]) - ft.rt) > rt_tol && continue
            abs(mean(tbl.mz1[loc]) - ft.mz1) > mz_tol && continue
            abs(mean(tbl.mz2[loc]) - ft.mz2) > mz_tol && continue
            tbl.collision_energy[loc[1]] == ft.collision_energy || continue
            push!(loc, i)
            new = false
            break
        end
        new && push!(locs, [i])
    end
    n > 1 && filter!(loc -> (length(loc) >= n && err(tbl.area[loc]) <= err_tol), locs)
    for (i, loc) in enumerate(locs)
        tbl.id[loc] .= i
    end
    @p tbl filter!(>(_.id, 0))
    #@p tbl |> DataFrame |> groupby(__, :id) |> combine(__, All() .=> mean, :area => err => :error, renamecols = false) |> Table
    gtbl = @p tbl groupview(getproperty(:id))
    errors = @p gtbl map(err(_.area)) collect
    @p gtbl map(map(mean, columns(Table(_; id = nothing)))) Table(Table(id = unique!(tbl.id)); error = errors)
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
process_product(product::Ion, cls::ClassSP, sc::SumChain) =
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
* `compound`: `Vector{CompoundSP}`; possible compounds.
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
            n = length(vectorize(data1.compound))
            tbl1[id1] = (;
                compound = vectorize(data1.compound),
                mz1 = (data1.mz1 * n + data2.mz1) / (n + 1),
                mz2 = (data1.mz2 * n + data2.mz2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = data1.collision_energy,
                polarity = data1.polarity
            )
            union!(data1.compound, vectorize(data2.compound))
            pushed = true
        end
        pushed || push!(tbl1, (compound = data2.compound, mz1 = data2.mz1, mz2 = data2.mz2, rt = data2.rt, Δrt = rt_tol * 2, collision_energy = data2.collision_energy, polarity = data2.polarity))
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
* `compound`: `Vector{CompoundSP}`; possible compounds.
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
            n = length(vectorize(data1.compound))
            isempty(setdiff(vectorize(data2.compound), vectorize(data1.compound))) || push!(extended, id1)
            tbl1[id1] = (;
                compound = vectorize(data1.compound),
                mz1 = (data1.mz1 * n + data2.mz1) / (n + 1),
                mz2 = (data1.mz2 * n + data2.mz2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = data1.collision_energy,
                polarity = data1.polarity
            )
            union!(data1.compound, vectorize(data2.compound))
            pushed = true
        end
        pushed && continue
        #push!(new, id2)
        push!(tbl1, (compound = data2.compound, mz1 = data2.mz1, mz2 = data2.mz2, rt = data2.rt, Δrt = rt_tol * 2, collision_energy = data2.collision_energy, polarity = data2.polarity))
    end
    (extended = tbl1[extended], new = size(tbl1, 1) > l ? tbl1[l + 1:end] : nothing)
end