

generate_mrm(project::Project, args...; kwargs...) = generate_mrm(project.analytes, args...; kwargs..., anion = project.anion)
generate_mrm(aquery::Query, args...; kwargs...) = generate_mrm(aquery.result, args...; kwargs..., anion = aquery.project.anion)
function generate_mrm(analytes::AbstractVector{AnalyteSP}, adduct, product, polarity::Bool = true; 
                    mz_tol = 0.35, rt_tol = 0.5, 
                    db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], 
                    anion = :acetate,
                    default = mz(Ion(ProtonationNL2H2O(), SPB3{2, 18}())))
    adduct = @match vectorize(adduct) begin
        [:default]                      => polarity
        add::Vector{Int}                => SPDB[:ADDUCTCODE].object[add]
        add::Vector{<: AbstractString}  => object_adduct.(add)
        add::Vector{<: Type}            => map(t -> t(), add)
        _                               => adduct
    end

    if isa(adduct, Bool)
        cpdlist = generate_cpdlist(analytes, adduct, anion)
    else
        if polarity 
            adduct = filter(t -> isa(t, Pos), adduct)
        elseif anion == :acetate
            adduct = filter(t -> isa(t, Neg) && !=(t, AddHCOO()), adduct)
        elseif anion == :formate
            adduct = filter(t -> isa(t, Neg) && !=(t, AddOAc()), adduct)
        end
    
        cpdlist = generate_cpdlist(analytes, adduct, anion)
    end
    @match product begin
        ::LCB   => (cpdlist = cpdlist[findall(cpd -> !isnothing(cpd.chain), analytes)])
        _       => product
    end
    productlist = generate_productlist(cpdlist, product, polarity)
    celist = generate_celist(cpdlist, product)
    tbl = DataFrame("Compound Name"     => Any[],
                    "Precursor Ion"     => Float64[],
                    "Product Ion"       => Float64[],
                    "Ret Time (min)"    => Float64[],
                    "Delta Ret Time"    => Float64[],
                    "Collision Energy"  => Int[],
                    )

    for (rcpd, product, ce) in zip(eachrow(cpdlist), productlist, celist)
        if product == 0 
            isnothing(default) && continue
            product = default
        end
        pushed = false
        ms1 = mz(rcpd[1], rcpd[2])
        for row in eachrow(tbl)
            (!between(row[2], ms1, mz_tol) || !between(row[3], product, mz_tol)) && continue
            !between(rcpd[3] - rt_tol, row[4], row[5] / 2) && !between(rcpd[3] + rt_tol, row[4], row[5] / 2) && continue
            row[6] == ce || continue
            rt_l = min(row[4] - row[5] / 2, rcpd[3] - rt_tol)
            rt_r = max(row[4] + row[5] / 2, rcpd[3] + rt_tol)
            row[1] = vectorize(row[1])
            row[2] = (row[2] * length(row[1]) + ms1) / (length(row[1]) + 1)
            row[3] = (row[3] * length(row[1]) + product) / (length(row[1]) + 1)
            row[4] = (rt_l + rt_r) / 2
            row[5] = rt_r - rt_l
            push!(row[1], repr(rcpd[1]))
            pushed = true
        end
        pushed || (push!(tbl, (repr(rcpd[1]), ms1, product, rcpd[3], rt_tol * 2, ce)))
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

function generate_cpdlist(analytes::AbstractVector{AnalyteSP}, adduct::Vector, anion)
    mapreduce(vcat, analytes) do analyte
        cpd = last(analyte)
        rt = analyte.rt
        mapreduce(vcat, adduct) do add
            [cpd add rt]
        end
    end
end

function generate_cpdlist(analytes::AbstractVector{AnalyteSP}, polarity::Bool, anion)
    mapreduce(vcat, analytes) do analyte
        cpd = last(analyte)
        rt = analyte.rt
        adduct = class_db_index(cpd.class).default_ion
        if polarity 
            adduct = filter(t -> isa(t, Pos), adduct)
        elseif anion == :acetate
            adduct = filter(t -> isa(t, Neg) && !=(t, AddHCOO()), adduct)
        elseif anion == :formate
            adduct = filter(t -> isa(t, Neg) && !=(t, AddOAc()), adduct)
        end
        mapreduce(vcat, adduct) do add
            [cpd add rt]
        end
    end
end

function generate_productlist(cpdlist::Matrix, product::Type{LCB}, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG])
    map(eachrow(cpdlist)) do row
        isnothing(row[1].chain) && return 0
        id = findfirst(x -> ==(row[1].chain.lcb, x.molecule), @view db[:, 1])
        isnothing(id) ? mz(default_adduct(row[1].chain.lcb)) : db[id, 2]
    end
end

function generate_celist(cpdlist::Matrix, product::Type{LCB})
    map(eachrow(cpdlist)) do row
        id = findfirst(x -> ==(row[1].class, x[1]) && ==(row[2], x[2]) && ==(x[3], "LCB"), eachrow(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE][id, :eV]
    end
end