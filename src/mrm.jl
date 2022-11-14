

generate_mrm(project::Project, args...; kwargs...) = generate_mrm(project.analytes, args...; kwargs..., anion = project.anion)
generate_mrm(aquery::Query, args...; kwargs...) = generate_mrm(aquery.result, args...; kwargs..., anion = aquery.project.anion)
function generate_mrm(analytes::Vector{AnalyteGSL}, adduct, product, polarity::Bool = true; mz_tol = 0.35, rt_tol = 1, db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], anion = :acetate)
    adduct = @match vectorize(adduct) begin
        [:default]                      => polarity
        add::Vector{Int}                => SPDB[:ADDUCTCODE].object[add]
        add::Vector{<: AbstractString}  => object_adduct.(add)
        add::Vector{<: Type}            => map(t -> t(), add)
        _                               => adduct
    end

    if polarity 
        adduct = filter(t -> isa(t, Pos), adduct)
    elseif anion == :acetate
        adduct = filter(t -> isa(t, Neg) && !=(t, AddHCOO()), adduct)
    elseif anion == :formate
        adduct = filter(t -> isa(t, Neg) && !=(t, AddOAc()), adduct)
    end

    cpdlist = generate_cpdlist(analytes, adduct, anion)
    productlist = generate_productlist(cpdlist, product, polarity)
    celist = generate_celist(cpdlist, product)
    n = length(celist)
    tbl = DataFrame("Compound Name"     => map(repr, cpdlist[:, 1]), 
                    "Precursor Ion"     => map(mz, cpdlist[:, 1], cpdlist[:, 2]),
                    "Product Ion"       => productlist,
                    "Ret Time (min)"    => cpdlist[:, 3],
                    "Delta Ret Time"    => repeat([rt_tol], n),
                    "Collision Energy"  => repeat([40], n), 
                    "Polarity"          => repeat([polarity ? "Positive" : "Negative"], n)
    )
    tbl
end

function generate_cpdlist(analytes::Vector{AnalyteGSL}, adduct::Vector, anion)
    mapreduce(vcat, analytes) do analyte
        cpd = last(analyte)
        rt = analyte.rt
        mapreduce(vcat, adduct) do add
            [cpd add rt]
        end
    end
end

function generate_cpdlist(analytes::Vector{AnalyteGSL}, polarity::Bool, anion)
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
        db[findfirst(x -> ==(row[1].chain.spb, x.molecule), db), 2]
    end
end

function generate_celist(cpdlist::Matrix, product::Type{LCB})
    map(eachrow(cpdlist)) do row
        SPDB[:CE][findfirst(x -> ==(row[1].class, x[1]) && ==(row[2], x[2]) && ==(row[3], "LCB"), eachrow(SPDB[:CE])), :eV]
    end
end