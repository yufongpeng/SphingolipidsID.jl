generate_mrm(project::Project, adduct, product, polarity::Bool; kwargs...) = generate_mrm(project.analytes, adduct, product, polarity; kwargs..., anion = project.anion)
generate_mrm(aquery::Query, adduct, product, polarity::Bool; kwargs...) = generate_mrm(aquery.result, adduct, product, polarity; kwargs..., anion = aquery.project.anion)
generate_mrm(adduct, product, polarity::Bool, pq::Union{Project, Query}; kwargs...) = generate_mrm(pq, adduct, product, polarity; kwargs...)

function generate_mrm(analytes::AbstractVector{AnalyteSP}, adduct, product, polarity::Bool; 
                    mz_tol = 0.35, rt_tol = 0.5, 
                    db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], 
                    anion = :acetate,
                    default = mz(Ion(ProtonationNL2H2O(), SPB3{2, 18}())))
    adduct = generate_adduct(vectorize(adduct))
    cpdlist = isempty(adduct) ? generate_cpdlist(analytes, polarity, anion) : 
        generate_cpdlist(analytes, filter_adduct(adduct, polarity, anion), anion)
    filter_cpdlist!(cpdlist, product)
    productlist = generate_productlist(cpdlist, product, polarity)
    celist = generate_celist(cpdlist, product)
    tbl = NamedTuple{(:compound, :mz1, :mz2, :rt, :Δrt, :collision_energy), Tuple{Any, Float64, Float64, Float64, Float64, Int64}}[]
    for (rcpd, ms2, ce) in zip(cpdlist, productlist, celist)
        if ms2 == 0 
            isnothing(default) && continue
            ms2 = default
        end
        pushed = false
        ms1 = mz(rcpd.cpd, rcpd.add)
        for (id, rtbl) in enumerate(tbl)
            (!between(rtbl.mz1, ms1, mz_tol) || !between(rtbl.mz2, ms2, mz_tol)) && continue
            !between(rcpd.rt - rt_tol, rtbl.rt, rtbl.Δrt / 2) && !between(rcpd.rt + rt_tol, rtbl.rt, rtbl.Δrt / 2) && continue
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
    tbl = DataFrame("Compound Name"     => map(repr, cpdlist[:, 1]), 
                    "Precursor Ion"     => map(mz, cpdlist[:, 1], cpdlist[:, 2]),
                    "Product Ion"       => productlist,
                    "Ret Time (min)"    => cpdlist[:, 3],
                    "Delta Ret Time"    => repeat([rt_tol], n),
                    "Collision Energy"  => celist, 
                    "Polarity"          => repeat([polarity ? "Positive" : "Negative"], n)
    )=#
end

generate_adduct(adduct::Vector{Symbol}) = @match adduct begin
    [:default] => Ion[]
end
generate_adduct(adduct::Vector{Int}) = SPDB[:ADDUCTCODE].object[adduct]
generate_adduct(adduct::Vector{<: AbstractString}) = object_adduct.(adduct)
generate_adduct(adduct::Vector{<: Type{T}}) where {T <: Ion} = map(t -> t(), adduct)
function filter_adduct(adduct::Vector{<: Adduct}, polarity::Bool, anion::Symbol)
    if polarity 
        filter(t -> isa(t, Pos), adduct)
    elseif anion == :acetate
        filter(t -> isa(t, Neg) && !=(t, AddHCOO()), adduct)
    elseif anion == :formate
        filter(t -> isa(t, Neg) && !=(t, AddOAc()), adduct)
    end
end

cpd_add_rt(analyte, add) = (cpd = last(analyte), add = add, rt = analyte.rt)
cpd_add_rt(analyte, polarity, anion) = Iterators.map(add -> cpd_add_rt(analyte, add), filter_adduct(class_db_index(last(analyte).class).default_ion, polarity, anion))

generate_cpdlist(analytes::AbstractVector{AnalyteSP}, adduct::Vector, anion) = 
    productview(cpd_add_rt, adduct, analytes) |> splitdimsview |> flatten
    # Transducers
    # Iterators.product(adduct, analytes) |> MapSplat(cpd_add_rt) |> collect
    # DataPipes and BangBang:
    # @p Iterators.product(adduct, analytes) |> Iterators.map(cpd_add_rt(_...)) |> reduce(push!!; init = [])
    # Much more allocation for 1st
    # SplitApplyCombine
    # productview(cpd_add_rt, adduct, analytes) |> splitdimsview |> flatten
    # A little more allocation for 1st, less allocation for 2nd

generate_cpdlist(analytes::AbstractVector{AnalyteSP}, polarity::Bool, anion) = 
    @p analytes |> mapmany(cpd_add_rt(_, polarity, anion))
    # Transducers
    # analytes |> MapCat(analyte -> cpd_add_rt(analyte, polarity, anion)) |> collect
    # DataPipes and BangBang:
    # @p analytes |> Iterators.map(cpd_add_rt(_, polarity, anion)) |> reduce(vcat)
    # Much less allocation for 1st
    # DataPipes and SplitApplyCombine
    
    # Even less allocation for 1st

filter_cpdlist!(cpdlist, product) = cpdlist
filter_cpdlist!(cpdlist, ::LCB) = filter!(cpd -> !isnothing(cpd.cpd.chain), cpdlist)

generate_productlist(cpdlist::Vector, ::Type{LCB}, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG]) = 
    map(cpdlist) do row
        isnothing(row.cpd.chain) && return 0
        id = findfirst(x -> ==(row.cpd.chain.lcb, x.molecule), db[:, 1])
        isnothing(id) ? mz(default_adduct(row.cpd.chain.lcb)) : db[id, 2]
    end

generate_celist(cpdlist::Vector, ::Type{LCB}) = 
    map(cpdlist) do row
        id = findfirst(x -> ==(row.cpd.class, SPDB[:CE].ms1[x]) && ==(row.add, SPDB[:CE].adduct1[x]) && ==(SPDB[:CE].ms2[x], "LCB"), eachindex(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE].eV[id]
    end

write_mrm(io, tbl::Table) = CSV.write(io, tbl; 
    header = ["Compound Name", "Precursor Ion", "Product Ion", "Ret Time (min)", "Delta Ret Time", "Collision Energy", "Polarity"])
