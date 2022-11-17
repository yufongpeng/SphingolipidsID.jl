

generate_mrm(project::Project, args...; kwargs...) = generate_mrm(project.analytes, args...; kwargs..., anion = project.anion)
generate_mrm(aquery::Query, args...; kwargs...) = generate_mrm(aquery.result, args...; kwargs..., anion = aquery.project.anion)
function generate_mrm(analytes::AbstractVector{AnalyteGSL}, adduct, product, polarity::Bool = true; mz_tol = 0.35, rt_tol = 1, db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], anion = :acetate)
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
    n = length(celist)
    tbl = DataFrame("Compound Name"     => map(repr, cpdlist[:, 1]), 
                    "Precursor Ion"     => map(mz, cpdlist[:, 1], cpdlist[:, 2]),
                    "Product Ion"       => productlist,
                    "Ret Time (min)"    => cpdlist[:, 3],
                    "Delta Ret Time"    => repeat([rt_tol], n),
                    "Collision Energy"  => celist, 
                    "Polarity"          => repeat([polarity ? "Positive" : "Negative"], n)
    )
    tbl
end

function generate_cpdlist(analytes::AbstractVector{AnalyteGSL}, adduct::Vector, anion)
    mapreduce(vcat, analytes) do analyte
        cpd = last(analyte)
        rt = analyte.rt
        mapreduce(vcat, adduct) do add
            [cpd add rt]
        end
    end
end

function generate_cpdlist(analytes::AbstractVector{AnalyteGSL}, polarity::Bool, anion)
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
        isnothing(row[1].chain) && return nothing
        id = findfirst(x -> ==(row[1].chain.lcb, x.molecule), @view db[:, 1])
        isnothing(id) ? mz(default_adduct(row[1].chain.lcb)) : db[id , 2]
    end
end

function generate_celist(cpdlist::Matrix, product::Type{LCB})
    map(eachrow(cpdlist)) do row
        id = findfirst(x -> ==(row[1].class, x[1]) && ==(row[2], x[2]) && ==(x[3], "LCB"), eachrow(SPDB[:CE]))
        isnothing(id) ? 40 : SPDB[:CE][id, :eV]
    end
end