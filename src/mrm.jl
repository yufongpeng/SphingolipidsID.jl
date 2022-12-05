# get rid of dataframe => tables.rowtables
generate_mrm(project::Project, args...; kwargs...) = generate_mrm(project.analytes, args...; kwargs..., anion = project.anion)
generate_mrm(aquery::Query, args...; kwargs...) = generate_mrm(aquery.result, args...; kwargs..., anion = aquery.project.anion)
function generate_mrm(analytes::AbstractVector{AnalyteSP}, adduct, product, polarity::Bool = true; 
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
    tbl = DataFrame("Compound Name"     => Any[],
                    "Precursor Ion"     => Float64[],
                    "Product Ion"       => Float64[],
                    "Ret Time (min)"    => Float64[],
                    "Delta Ret Time"    => Float64[],
                    "Collision Energy"  => Int[],
                    )

    for (rcpd, mz2, ce) in zip(cpdlist, productlist, celist)
        if mz2 == 0 
            isnothing(default) && continue
            mz2 = default
        end
        pushed = false
        ms1 = mz(rcpd.cpd, rcpd.add)
        for row in eachrow(tbl)
            (!between(row[2], ms1, mz_tol) || !between(row[3], mz2, mz_tol)) && continue
            !between(rcpd.rt - rt_tol, row[4], row[5] / 2) && !between(rcpd.rt + rt_tol, row[4], row[5] / 2) && continue
            row[6] == ce || continue
            rt_l = min(row[4] - row[5] / 2, rcpd.rt - rt_tol)
            rt_r = max(row[4] + row[5] / 2, rcpd.rt + rt_tol)
            row[1] = vectorize(row[1])
            row[2] = (row[2] * length(row[1]) + ms1) / (length(row[1]) + 1)
            row[3] = (row[3] * length(row[1]) + mz2) / (length(row[1]) + 1)
            row[4] = (rt_l + rt_r) / 2
            row[5] = rt_r - rt_l
            push!(row[1], repr(rcpd.cpd))
            pushed = true
        end
        pushed || (push!(tbl, (repr(rcpd.cpd), ms1, mz2, rcpd.rt, rt_tol * 2, ce)))
    end
    insertcols!(tbl, "Polarity" => repeat([polarity ? "Positive" : "Negative"], nrow(tbl)))
    #=
    tbl = DataFrame("Compound Name"     => map(repr, cpdlist[:, 1]), 
                    "Precursor Ion"     => map(mz, cpdlist[:, 1], cpdlist[:, 2]),
                    "Product Ion"       => productlist,
                    "Ret Time (min)"    => cpdlist[:, 3],
                    "Delta Ret Time"    => repeat([rt_tol], n),
                    "Collision Energy"  => celist, 
                    "Polarity"          => repeat([polarity ? "Positive" : "Negative"], n)
    )=#
    tbl
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
filter_cpdlist!(cpdlist, product::LCB) = filter!(cpd -> !isnothing(cpd.cpd.chain), cpdlist)

generate_productlist(cpdlist::Vector, product::Type{LCB}, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG]) = 
    map(cpdlist) do row
        isnothing(row.cpd.chain) && return 0
        id = findfirst(x -> ==(row.cpd.chain.lcb, x.molecule), @view db[:, 1])
        isnothing(id) ? mz(default_adduct(row.cpd.chain.lcb)) : db[id, 2]
    end

generate_celist(cpdlist::Vector, product::Type{LCB}) = 
    map(cpdlist) do row
        id = findfirst(x -> ==(row.cpd.class, x[1]) && ==(row.add, x[2]) && ==(x[3], "LCB"), eachrow(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE].eV[id]
    end