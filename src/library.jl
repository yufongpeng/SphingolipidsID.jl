function read_adduct_code(file::String) 
    code = CSV.read(file, DataFrame, header = [:repr, :object])
    transform!(code, :object => ByRow(eval ∘ Meta.parse), renamecols = false)
end

const ADDUCTCODE = read_adduct_code(joinpath(@__DIR__, "..", "config", "ADDUCTCODE.csv"))

SPDB[:ADDUCTCODE] = ADDUCTCODE
SPDB[:NLH2O] = @view SPDB[:ADDUCTCODE].object[1:4]

repr_adduct(x::Adduct) = SPDB[:ADDUCTCODE].repr[findfirst(==(x), SPDB[:ADDUCTCODE].object)]
object_adduct(x::AbstractString) = SPDB[:ADDUCTCODE].object[findfirst(==(x), SPDB[:ADDUCTCODE].repr)]

function read_class_db(file::String, 
                        col_name = :Abbreviation, 
                        col_regex = :regex, 
                        col_formula = [:formula, :u1, :u2], 
                        col_unit = ["#us", "#u2"], 
                        col_ions = [:default_ion, :parent_ion, :adduct_ion, :isf_ion]
                        )

    class_db = CSV.read(file, DataFrame)
    #regexs = filter(x -> !isnothing(match(r"^regex", x)), names(class_db))
    transform!(class_db, col_name => ByRow(eval ∘ Meta.parse), renamecols = false)
    transform!(class_db, col_regex .=> ByRow(eval ∘ Meta.parse), renamecols = false)
    transform!(class_db, col_formula .=> ByRow(parse_compound ∘ String), renamecols = false)
    transform!(class_db, col_ions .=> ByRow(x -> SPDB[:ADDUCTCODE].object[(eval ∘ Meta.parse)(x)]), renamecols = false)
    rename!(class_db, [col_name => :Abbreviation, col_regex => :regex, (col_formula .=> [:formula, :u1, :u2])..., (col_ions .=> [:default_ion, :parent_ion, :adduct_ion, :isf_ion])...])
    rename!(class_db, col_unit .=> ["#us", "#u2"])
    class_db
end

const CLASSDB = read_class_db(joinpath(@__DIR__, "..", "config", "CLASSDB.csv"))
SPDB[:CLASSDB] = CLASSDB

class_db_index(::LCB) = @views SPDB[:CLASSDB][findfirst(==(SPB), SPDB[:CLASSDB][!, :Abbreviation]), :]
class_db_index(cls::T) where {T <: ClassSP} = @views SPDB[:CLASSDB][findfirst(==(deisomerized(T)), SPDB[:CLASSDB][!, :Abbreviation]), :]

library(class::Vector, adduct::Vector{<: Adduct}, range::Vector{<: Tuple}) = library(class, map(repr_adduct, adduct), range)
#library(class::Vector{Type{<: ClassSP}}, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple}) = library([cls() for cls in class], adduct, range)

function library(class::Vector, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple})
    range = map(range) do r
        vectorize.(r)
    end
    n = mapreduce(rng -> mapreduce(length, *, rng), +, range) * length(adduct) * length(class)
    df = DataFrame(
                    :Abbreviation => Vector{Type{<: ClassSP}}(undef, n),
                    :Species => Vector{String}(undef, n),
                    :Formula => Vector{String}(undef, n),
                    :Adduct => repeat(adduct, Int(n / length(adduct))),
                    Symbol("m/z") => zeros(Float64, n)
                    )

    i = 1
    #init_cndb == [init_cn, init_db] || throw(ArgumentError("Unmatched carbon and double bonds number of initial compounds."))
    map!(add -> replace(add, "FA-H" => "HCOO", "OAc" => "CH3COO"), adduct, adduct)
    add_fn = parse_adduct.(adduct)

    # n_within? # per adduct
    for (cls, rng) in Iterators.product(class, range)            
        # specific for species
        # match cls with Abbreviation first
        id_compound = findfirst(==(cls), SPDB[:CLASSDB][!, :Abbreviation])
        # isnothing(id_compound) && (id_compound = findfirst(abbr -> match(Regex(cls * "-.*"), abbr), CLASSDB[!, :Abbreviation]))
        init_elements = SPDB[:CLASSDB][id_compound, :formula]
        unit = SPDB[:CLASSDB][id_compound, [:u1, :u2]]
        init_us = SPDB[:CLASSDB][id_compound, ["#u1", "#u2"]]

        if length(rng) == 3
            posts = vectorize(last(rng))
            rng = rng[1:end - 1]
        else
            posts = ""
        end

        if length(rng) == 0
            rng = (0, 0)
        elseif length(rng) == 1
            rng = (first(rng), 0)
        end

        if posts == ""
            post_sep = [split(SPDB[:CLASSDB][id_compound, :regex].pattern, ";")[2:end]]
            if any(x -> occursin("\\d", x), post_sep[1])
                throw(ArgumentError("Custom modification, i.e. hydroxylation, glycosylation, etc, must be provided"))
            end
            posts = [join(post_sep, ";")]
        else
            post_sep = map(x -> split(x, ";"), posts)
        end
        init_elements = map(post_sep) do post
            add_addition(init_elements, post)
        end
        pre = ""


        for (post, elements) in zip(posts, init_elements), u2 in last(rng), u1 in first(rng)
            Δu = [u1 - init_us[1], u2 - init_us[2]]
            newf = merge_formula(elements, unit, Δu; sign = (:+, :-))
            for fn in add_fn
                df[i, :Abbreviation] = cls
                df[i, :Species] = merge_species(u1, u2, repr(cls), pre, post)
                df[i, :Formula] = newf
                df[i, Symbol("m/z")] = fn(mw(newf))
                i += 1
            end
        end
    end
    df
end

library(class, adduct, range) = library(vectorize(class), vectorize(adduct), vectorize(range))

const LIBRARY_POS = reduce(append!, (
    library(Cer, ["[M+H]+"], [(32:46, 1:2, "O"), (32:46, 0:4, "2O"), (32:46, 0:3, "3O")]),
    library([HexCer, Hex2Cer, Hex3Cer, GM3], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "2O"), (32:46, 0:3, "3O")]),
    library(HexNAcHex2Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "2O")),
    library(HexNAcHex3Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "2O")),
    library([GM2, GM1, GD3, GD2, GD1, GT3, GT2, GT1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, "2O")),
    library(GQ1, "[M+2H]2+", (32:46, 0:4, "2O")),
    library(GP1, ["[M+2H]2+", "[M+3H]3+"], (34:2:38, 0:4, "2O")),
))
SPDB[:LIBRARY_POS] = LIBRARY_POS

#LIBRARY_NEG

const FRAGMENT_POS = let
    frags = [SPB2{1, 18}(), SPB2{2, 18}(), 
            SPB3{2, 16}(), SPB3{2, 17}(), SPB3{2, 18}(), SPB3{3, 18}(), SPB3{2, 19}(), SPB3{2, 20}(), 
            SPB4{2, 16}(), SPB4{2, 18}(), SPB4{3, 18}(), SPB4{2, 20}()]
    ion = map(default_adduct, frags)
    append!(ion, [Ion(ProtonationNL2H2O(), NeuAc()), 
                    Ion(ProtonationNLH2O(), NeuAc()), 
                    HexNAcHex, 
                    HexNAcHexNANA, 
                    HexNAcHexNANA2, 
                    Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())), 
                    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc())), 
                    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc()))
                    ])
    [ion map(mz, ion)]
end

SPDB[:FRAGMENT_POS] = FRAGMENT_POS

#FRAGMENT_NEG


const CONNECTION = Dict{ClassSP, Any}(
    GP1c()                  => (GQ1c(), GQ1b()),
    GP1cα()                 => (GQ1bα(), GQ1c()),
    GQ1c()                  => (GT1c(), GT1b()),
    GQ1b()                  => (GT1a(), GT1b()),
    GQ1bα()                 => (GT1aα(), GT1b()),
    GT1c()                  => GT2(),
    GT1b()                  => (GD1b(), GD1a()),
    GT1a()                  => (GD1a(), GD1c()),
    GT1aα()                 => (GD1a(), GD1α()),
    GD1b()                  => GD2(),
    GD1a()                  => (GM1a(), GM1b()),
    GD1c()                  => GM1b(),
    GD1α()                  => GM1b(),
    GM1a()                  => GM2(),
    GM1b()                  => Hex_HexNAc_Hex2Cer(),
    GT2()                   => (GT3(), GD2()),
    GD2()                   => (GD3(), GM2()),
    GM2()                   => (GM3(), HexNAc_Hex2Cer()),
    GT3()                   => GD3(),
    GD3()                   => GM3(),
#    GP1()                   => GQ1(),
#    GQ1()                   => GT1(),
#    GT1()                   => GD1(),
#    GD1()                   => GM1(),
    GM3()                   => Hex2Cer(),
    Hex2Cer()               => HexCer(),
    Cer()                   => Cer(),
    SHexCer()               => HexCer(),
    SHexHexCer()            => Hex2Cer(),
    GM4()                   => HexCer(),
    Hex3Cer()               => Hex2Cer(),
    HexNAc_Hex2Cer()        => Hex2Cer(),
    Hex_HexNAc_Hex2Cer()    => HexNAc_Hex2Cer(),
    HexNAc_Hex3Cer()        => Hex3Cer(),
    HexCer()                => Cer()
)

SPDB[:CONNECTION] = CONNECTION


function read_ce(file::String) 
    ce = CSV.read(file, DataFrame; select = 1:5)
    transform!(ce, :ms1 => ByRow(t -> (eval ∘ Meta.parse)(t)()), :adduct1 => ByRow(x -> SPDB[:ADDUCTCODE].object[x]), renamecols = false)
    ce.addduct2 = map(ce.adduct2) do adduct
        @match adduct begin
            0 => :default
            _ => SPDB[:ADDUCTCODE].object[adduct]
        end
    end
    ce
end

const CE = read_ce(joinpath(@__DIR__, "..", "config", "CE.csv"))
SPDB[:CE] = CE