function featuretable_mzmine(paths)
    tbl = CSV.read(paths, DataFrame; select = (i, name) -> any(startswith(String(name), text) for text in ["id", "rt", "mz", "height", "area", "intensity"]))
    sort!(tbl, [:mz, :rt])
end
function filter_duplicate(tbl::DataFrame; rt_tol = 0.1, mz_tol = 0.35)
    ids = Vector{Int}[]
    for (i, ft) in enumerate(eachrow(tbl))
        new = true
        for id in ids
            if abs(mean(tbl[id, :rt]) - ft.rt) < rt_tol && abs(mean(tbl[id, :mz]) - ft.mz) < mz_tol
                push!(id, i)
                new = false
                break
            end
        end
        new && push!(ids, [i])
    end
    for (i, id) in enumerate(ids)
        tbl[id, :id] .= i
    end
    combine(groupby(tbl, :id), All() .=> mean, renamecols = false)
end

function CompoundGSL(cpd, products, source, id)
    sum_cb, sum_db, sum_o = match(r".+ (\d+):(\d+);(\d*)O", cpd.Species).captures
    sum_cb = parse(Int, sum_cb)
    sum_db = parse(Int, sum_db)
    sum_o = isempty(sum_o) ? 1 : parse(Int, sum_o)
    class = (cpd.Abbreviation)()
    adduct_class = object_adduct(cpd.Adduct)
    chain = nothing
    del = Int[]
    for (i, product) in enumerate(products)
        @match product begin
            Ion(adduct, lcb::LCB) => begin 
                lcb_o = nhydroxyl(lcb)
                lcb_db = ndb(lcb)
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
                if (sum_db + sum_o) < (lcb_o + lcb_db)
                    push!(del, i)
                    continue
                end
                Δ = sum_o - lcb_o 
                Δ = Δ < 0 ? Δ : max(0, lcb_db - sum_db)
                acyl_o = sum_o - lcb_o - Δ
                if Δ != 0
                    lcb_new = hydroxyl_shift(lcb, Δ)
                    add_new = hydroxyl_shift(product.adduct, Δ)
                    if isnothing(lcb_new) || isnothing(add_new)
                        push!(del, i)
                        continue
                    else
                        products[i] = Ion(add_new, lcb_new) 
                    end
                end
                chain = Chain(products[i].molecule, Acyl{acyl_o}())
            end
            x && if x in NANA end       => begin
                hasnana(cpd) || push!(del, i)
            end
        end
    end
    fragments = DataFrame(ion1 = Ion[], ion2 = Ion[], source = Int[], id = Int[])
    deleteat!(products, del)
    isempty(products) && return nothing
    for product in products
        push!(fragments, (ion1 = Ion(adduct_class, class), ion2 = product, source = source, id = id))
    end
    CompoundGSL(class, [sum_cb, sum_db, sum_o], chain, 
        fragments, 
        Union{Missing, Pair, IonMode, IonUnion, Nothing}[missing, missing]
    )
end

function id_product(ms2, polarity; db = polarity ? DEFAULT_POS_FRAGMENT : DEFAULT_NEG_FRAGMENT, mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(db)
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end