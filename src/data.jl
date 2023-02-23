sort_data(tbl::Table) = sort_data!(tbl)
function sort_data!(tbl::Table)
    sort!(tbl, :datafile)
    tbl
end

file_order(tbl::Table) = in(:datafile, propertynames(tbl)) ? unique(tbl.datafile) : throw(ArgumentError("This table doesn't contain datafile."))
function fill_mz2!(tbl::Table, mz2::Float64)
    fill!(tbl.mz2, mz2)
    tbl
end

function fill_mz2!(tbl::Table{
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile),
        Tuple{Int64, Float64, Float64, Float64, Float64, Float64, Int64, Float64, Float64, SubString{String}}}, 1,
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile),
        Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{SubString{String}}}}
    },
    mz2::Union{<: Vector, <: Tuple})
    mapping = Dict(unique(tbl.datafile) .=> mz2)
    tbl.mz2 .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

function fill_ce!(tbl::Table, eV::Float64)
    fill!(tbl.collision_energy, eV)
    tbl
end

function fill_ce!(tbl::Table{
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile),
        Tuple{Int64, Float64, Float64, Float64, Float64, Float64, Int64, Float64, Float64, SubString{String}}}, 1,
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile),
        Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{SubString{String}}}}
    },
    eV::Union{<: Vector, <: Tuple})
    mapping = Dict(unique(tbl.datafile) .=> eV)
    tbl.collision_energy .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

fill_mz2!(tbl, mz2) = tbl
fill_ce!(tbl, eV) = tbl

rsd(v) = std(v) / mean(v)
re(v) =  - foldl(-, extrema(v)) / mean(v) / 2
default_error(v) = length(v) > 2 ? rsd(v) : re(v)

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

function compoundspvanilla(cpd, product)
    cls = (cpd.Abbreviation)()
    sc = sumcomp(cpd.Species)
    adduct_class = object_adduct(cpd.Adduct)
    pr = process_product(product, cls, sc)
    isnothing(pr) && return pr
    ion2, sc = pr
    CompoundSPVanilla(cls, sc,
        Table(ion1 = Ion[Ion(adduct_class, cls)], ion2 = ion2, source = [1], id = [1])
    )
end

process_product(product::Ion, cls::ClassSP, sc::SumChain) =
    (in(product, NANA) && !hasnana(cls)) ? nothing : (Ion[product], sc)
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
    Ion[product_new], DiChain(lcb_new, Acyl(ncb(sc) - ncb(lcb_new), sum_db - ndb(lcb_new), sum_o - nox(lcb_new)))
end

AnalyteSP(compounds::Vector{CompoundSP}, rt::Float64) = AnalyteSP(compounds, rt, repeat([0], 5), repeat([0.0], 5), 0)

function id_product(ms2, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(db)
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end

union_mrm(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = union_mrm!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
function union_mrm!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    for id2 in eachindex(tbl2)
        pushed = false
        ms1 = tbl2.mz1[id2]
        ms2 = tbl2.mz2[id2]
        ce = tbl2.collision_energy[id2]
        for id1 in eachindex(tbl1)
            (!between(tbl1.mz1[id1], ms1, mz_tol) || !between(tbl1.mz2[id1], ms2, mz_tol)) && continue
            tbl1.collision_energy[id1] == ce || continue
            rt_l = min(tbl1.rt[id1] - tbl1.Δrt[id1] / 2, tbl2.rt[id2] - tbl2.Δrt[id2] / 2)
            rt_r = max(tbl1.rt[id1] + tbl1.Δrt[id1] / 2, tbl2.rt[id2] + tbl2.Δrt[id2] / 2)
            n = length(vectorize(tbl1.compound[id1]))
            tbl1[id1] = (;
                compound = vectorize(tbl1.compound[id1]),
                mz1 = (tbl1.mz1[id1] * n + ms1) / (n + 1),
                mz2 = (tbl1.mz2[id1] * n + ms2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = tbl1.collision_energy[id1],
                polarity = tbl1.polarity[id1]
            )
            union!(tbl1.compound[id1], vectorize(tbl2.compound[id2]))
            pushed = true
        end
        pushed || push!(tbl1, (compound = tbl2.compound[id2], mz1 = ms1, mz2 = ms2, rt = tbl2.rt[id2], Δrt = rt_tol * 2, collision_energy = ce, polarity = tbl2.polarity[id2]))
    end
    sort!(tbl1, [:mz1, :rt])
end

diff_mrm(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = diff_mrm!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
function diff_mrm!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    extended = Int[]
    l = size(tbl1, 1)
    for id2 in eachindex(tbl2)
        pushed = false
        ms1 = tbl2.mz1[id2]
        ms2 = tbl2.mz2[id2]
        ce = tbl2.collision_energy[id2]
        for id1 in eachindex(tbl1)
            (!between(tbl1.mz1[id1], ms1, mz_tol) || !between(tbl1.mz2[id1], ms2, mz_tol)) && continue
            tbl1.collision_energy[id1] == ce || continue
            rt_l = min(tbl1.rt[id1] - tbl1.Δrt[id1] / 2, tbl2.rt[id2] - tbl2.Δrt[id2] / 2)
            rt_r = max(tbl1.rt[id1] + tbl1.Δrt[id1] / 2, tbl2.rt[id2] + tbl2.Δrt[id2] / 2)
            n = length(vectorize(tbl1.compound[id1]))
            isempty(setdiff(vectorize(tbl2.compound[id2]), vectorize(tbl1.compound[id1]))) || push!(extended, id1)
            tbl1[id1] = (;
                compound = vectorize(tbl1.compound[id1]),
                mz1 = (tbl1.mz1[id1] * n + ms1) / (n + 1),
                mz2 = (tbl1.mz2[id1] * n + ms2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = tbl1.collision_energy[id1],
                polarity = tbl1.polarity[id1]
            )
            union!(tbl1.compound[id1], vectorize(tbl2.compound[id2]))
            pushed = true
        end
        pushed && continue
        push!(new, id2)
        push!(tbl1, (compound = tbl2.compound[id2], mz1 = ms1, mz2 = ms2, rt = tbl2.rt[id2], Δrt = rt_tol * 2, collision_energy = ce, polarity = tbl2.polarity[id2]))
    end
    (extended = tbl1[extended], new = size(tbl1, 1) > l ? tbl1[l + 1:end] : nothing)
end