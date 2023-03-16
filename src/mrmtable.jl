"""
    mrmtable(project::Project, adduct, product, polarity::Bool; kwargs...)
    mrmtable(aquery::AbstractQuery, adduct, product, polarity::Bool; kwargs...)
    mrmtable(adduct, product, polarity::Bool, pq::Union{Project, AbstractQuery}; kwargs...)
    mrmtable(analytes::AbstractVector{AnalyteSP}, adduct, product, polarity::Bool;
            mz_tol = 0.35, rt_tol = 0.5,
            db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
            anion = :acetate,
            default = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1))))

Create a table of MRM transition.

# Arguuments
* `analytes`: analytes to be included.
* `adduct`: an adduct or a vector of adducts for parent ions.
* `product`: `Ion` or `LCB`; product ion. There are default adducts for `LCB`. 
* `polarity`: `true` for positive ion mode; `false` for negative ion mode.

# Keyword Arguuments
* `mz_tol`: Maximum allowable difference in m/z.
* `anion`: `:formate` or `acetate`.
* `default`: default m/z value.
"""
mrmtable(project::Project, adduct, product, polarity::Bool; kwargs...) = mrmtable(project.analytes, adduct, product, polarity; kwargs..., anion = project.anion)
mrmtable(aquery::AbstractQuery, adduct, product, polarity::Bool; kwargs...) = mrmtable(aquery.result, adduct, product, polarity; kwargs..., anion = aquery.project.anion)
mrmtable(adduct, product, polarity::Bool, pq::Union{Project, AbstractQuery}; kwargs...) = mrmtable(pq, adduct, product, polarity; kwargs...)

function mrmtable(analytes::AbstractVector{AnalyteSP}, adduct, product, polarity::Bool;
                    mz_tol = 0.35, rt_tol = 0.5,
                    anion = :acetate,
                    default = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1))))
    adduct = proc_adduct(vectorize(adduct))
    cpd_list = isempty(adduct) ? cpdlist(analytes, polarity, anion) :
        cpdlist(analytes, filter_adduct(adduct, polarity, anion), anion)
    filter_cpdlist!(cpd_list, product)
    product_list = productlist(cpd_list, product, polarity)
    ce_list = celist(cpd_list, product)
    tbl = NamedTuple{(:compound, :mz1, :mz2, :rt, :Δrt, :collision_energy), Tuple{Any, Float64, Float64, Float64, Float64, Int64}}[]
    for (rcpd, ms2, ce) in zip(cpd_list, product_list, ce_list)
        if ms2 == 0
            isnothing(default) && continue
            ms2 = default
        end
        pushed = false
        ms1 = mz(rcpd.cpd, rcpd.add)
        for (id, rtbl) in enumerate(tbl)
            (!between(rtbl.mz1, ms1, mz_tol) || !between(rtbl.mz2, ms2, mz_tol)) && continue
            #!between(rcpd.rt - rt_tol, rtbl.rt, rtbl.Δrt / 2) && !between(rcpd.rt + rt_tol, rtbl.rt, rtbl.Δrt / 2) && continue
            rtbl.collision_energy == ce || continue
            rt_l = min(rtbl.rt - rtbl.Δrt / 2, rcpd.rt - rt_tol)
            rt_r = max(rtbl.rt + rtbl.Δrt / 2, rcpd.rt + rt_tol)
            tbl[id] = (;
                compound = vectorize(rtbl.compound),
                mz1 = (rtbl.mz1 * length(rtbl.compound) + ms1) / (length(rtbl.compound) + 1),
                mz2 = (rtbl.mz2 * length(rtbl.compound) + ms2) / (length(rtbl.compound) + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = rtbl.collision_energy
            )
            push!(tbl[id].compound, repr(rcpd.cpd))
            pushed = true
        end
        pushed || (push!(tbl, (compound = repr(rcpd.cpd), mz1 = ms1, mz2 = ms2, rt = rcpd.rt, Δrt = rt_tol * 2, collision_energy = ce)))
    end
    Table(tbl, polarity = repeat([polarity ? "Positive" : "Negative"], size(tbl, 1)))
    #=
    tbl = DataFrame("Compound Name"     => map(repr, cpd_list[:, 1]),
                    "Precursor Ion"     => map(mz, cpd_list[:, 1], cpd_list[:, 2]),
                    "Product Ion"       => product_list,
                    "Ret Time (min)"    => cpd_list[:, 3],
                    "Delta Ret Time"    => repeat([rt_tol], n),
                    "Collision Energy"  => ce_list,
                    "Polarity"          => repeat([polarity ? "Positive" : "Negative"], n)
    )=#
end

proc_adduct(adduct::Vector{Symbol}) = @match adduct begin
    [:default] => Ion[]
end
proc_adduct(adduct::Vector{Int}) = SPDB[:ADDUCTCODE].object[adduct]
proc_adduct(adduct::Vector{<: AbstractString}) = object_adduct.(adduct)
proc_adduct(adduct::Vector{<: Type{T}}) where {T <: Ion} = map(t -> t(), adduct)
function filter_adduct(adduct::Vector{<: Adduct}, polarity::Bool, anion::Symbol)
    if polarity
        filter(t -> isa(t, Pos), adduct)
    elseif anion ≡ :acetate
        filter(t -> isa(t, Neg) && !=(t, AddHCOO()), adduct)
    elseif anion ≡ :formate
        filter(t -> isa(t, Neg) && !=(t, AddOAc()), adduct)
    end
end

cpd_add_rt(analyte, add) = (cpd = last(analyte), add = add, rt = analyte.rt)
cpd_add_rt(analyte, polarity, anion) = Iterators.map(add -> cpd_add_rt(analyte, add), filter_adduct(class_db_index(class(analyte)).default_ion, polarity, anion))

cpdlist(analytes::AbstractVector{AnalyteSP}, adduct::Vector, anion) =
    productview(cpd_add_rt, adduct, analytes) |> splitdimsview |> flatten
    # Transducers
    # Iterators.product(adduct, analytes) |> MapSplat(cpd_add_rt) |> collect
    # DataPipes and BangBang:
    # @p Iterators.product(adduct, analytes) |> Iterators.map(cpd_add_rt(_...)) |> reduce(push!!; init = [])
    # Much more allocation for 1st
    # SplitApplyCombine
    # productview(cpd_add_rt, adduct, analytes) |> splitdimsview |> flatten
    # A little more allocation for 1st, less allocation for 2nd

cpdlist(analytes::AbstractVector{AnalyteSP}, polarity::Bool, anion) =
    @p analytes |> mapmany(cpd_add_rt(_, polarity, anion))
    # Transducers
    # analytes |> MapCat(analyte -> cpd_add_rt(analyte, polarity, anion)) |> collect
    # DataPipes and BangBang:
    # @p analytes |> Iterators.map(cpd_add_rt(_, polarity, anion)) |> reduce(vcat)
    # Much less allocation for 1st
    # DataPipes and SplitApplyCombine
    # Even less allocation for 1st

filter_cpdlist!(cpd_list, product) = cpd_list
filter_cpdlist!(cpd_list, ::Ion{<: Pos, NeuAc}) = filter!(cpd -> isa(class(cpd.cpd), CLS.fg.nana), cpd_list)
filter_cpdlist!(cpd_list, ::LCB) = filter!(cpd -> !isspceieslevel(cpd.cpd), cpd_list)

productlist(cpd_list::Vector, ::Type{LCB}, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG]) =
    map(cpd_list) do row
        isspceieslevel(row.cpd) && return 0
        id = findfirst(x -> ≡(lcb(row.cpd), x.molecule), db[:, 1])
        isnothing(id) ? mz(default_adduct(lcb(row.cpd))) : db[id, 2]
    end
productlist(cpd_list::Vector, ion::Ion{<: Pos, NeuAc}, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG]) =
    repeat([db[findfirst(==(ion), db[:, 1]), 2]], size(cpd_list, 1))

celist(cpd_list::Vector, ::Type{LCB}) =
    map(cpd_list) do row
        id = findfirst(x -> ≡(class(row.cpd), SPDB[:CE].ms1[x]) && ≡(row.add, SPDB[:CE].adduct1[x]) && ==(SPDB[:CE].ms2[x], "LCB"), eachindex(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE].eV[id]
    end

celist(cpd_list::Vector, ion::Ion{<: Pos, NeuAc}) =
    map(cpd_list) do row
        id = findfirst(x -> ≡(class(row.cpd), SPDB[:CE].ms1[x]) && ≡(row.add, SPDB[:CE].adduct1[x]) &&
                            ==(SPDB[:CE].ms2[x], "NeuAc") && ≡(ion.adduct, SPDB[:CE].adduct2[x]), eachindex(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE].eV[id]
    end