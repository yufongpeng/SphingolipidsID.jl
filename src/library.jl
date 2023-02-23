function read_adduct_code(file::String)
    code = CSV.read(file, Table, header = [:repr, :object])
    Table(
        repr = code.repr,
        object = (eval ∘ Meta.parse).(code.object)
    )
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
                        col_unit = [:nu1, :nu2],
                        col_ions = [:default_ion, :parent_ion, :adduct_ion, :isf_ion]
                        )

    class_db = CSV.read(file, Table)
    parse_ion = x -> SPDB[:ADDUCTCODE].object[(eval ∘ Meta.parse)(x)]
    Table(
        #Abbreviation = map(x -> (eval ∘ Meta.parse)(x)(), getproperty(class_db, col_name)),
        Abbreviation = (eval ∘ Meta.parse).(getproperty(class_db, col_name)),
        regex = (eval ∘ Meta.parse).(getproperty(class_db, col_regex)),
        formula = (parse_compound ∘ String).(getproperty(class_db, col_formula[1])),
        nu1 = getproperty(class_db, col_unit[1]),
        u1 = (parse_compound ∘ String).(getproperty(class_db, col_formula[2])),
        nu2 = getproperty(class_db, col_unit[2]),
        u2 = (parse_compound ∘ String).(getproperty(class_db, col_formula[3])),
        default_ion = parse_ion.(getproperty(class_db, col_ions[1])),
        parent_ion = parse_ion.(getproperty(class_db, col_ions[2])),
        adduct_ion = parse_ion.(getproperty(class_db, col_ions[3])),
        isf_ion = parse_ion.(getproperty(class_db, col_ions[4]))
    )
end

const CLASSDB = read_class_db(joinpath(@__DIR__, "..", "config", "CLASSDB.csv"))
SPDB[:CLASSDB] = CLASSDB

class_db_index(spb::LCB) = SPDB[:CLASSDB][findfirst(==(spb), SPDB[:CLASSDB].Abbreviation)]
#class_db_index(cls::T) where {T <: ClassSP} = SPDB[:CLASSDB][findfirst(==(deisomerized(cls)), SPDB[:CLASSDB].Abbreviation)]
class_db_index(cls::T) where {T <: ClassSP} = SPDB[:CLASSDB][findfirst(==(deisomerized(T)), SPDB[:CLASSDB].Abbreviation)]

library(class::Vector, adduct::Vector{<: Adduct}, range::Vector{<: Tuple}) = library(class, map(repr_adduct, adduct), range)
#library(class::Vector{Type{<: ClassSP}}, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple}) = library([cls() for cls in class], adduct, range)

function library(class::Vector, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple})
    range = map(range) do r
        vectorize.(r)
    end
    n = mapreduce(rng -> mapreduce(length, *, rng), +, range) * length(adduct) * length(class)
    tbl = Table(
                    Abbreviation = Vector{Type{<: ClassSP}}(undef, n),
                    Species = Vector{SPID}(undef, n),
                    Formula = Vector{String}(undef, n),
                    Adduct = repeat(adduct, Int(n / length(adduct))),
                    mz = zeros(Float64, n)
                )

    i = 1
    #init_cndb == [init_cn, init_db] || throw(ArgumentError("Unmatched carbon and double bonds number of initial compounds."))
    map!(add -> replace(add, "FA-H" => "HCOO", "OAc" => "CH3COO"), adduct, adduct)
    add_fn = parse_adduct.(adduct)

    # n_within? # per adduct
    for (cls, rng) in Iterators.product(class, range)
        # specific for species
        # match cls with Abbreviation first
        #id_compound = findfirst(==(cls()), SPDB[:CLASSDB].Abbreviation)
        id_compound = findfirst(==(cls), SPDB[:CLASSDB].Abbreviation)
        # isnothing(id_compound) && (id_compound = findfirst(abbr -> match(Regex(cls * "-.*"), abbr), CLASSDB[!, :Abbreviation]))
        init_elements = SPDB[:CLASSDB].formula[id_compound]
        unit = getproperties(SPDB[:CLASSDB], (:u1, :u2))[id_compound]
        init_us = getproperties(SPDB[:CLASSDB], (:nu1, :nu2))[id_compound]

        if length(rng) ≡ 3
            posts = vectorize(last(rng))
            rng = rng[1:end - 1]
        else
            posts = ""
        end

        if length(rng) ≡ 0
            rng = (0, 0)
        elseif length(rng) ≡ 1
            rng = (first(rng), 0)
        end

        if posts == ""
            post_sep = [split(SPDB[:CLASSDB].regex[id_compound].pattern, ";")[2:end]]
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

        sp_posts = map(posts) do post
            o = split(post, "O")[2]
            isempty(o) ? 1 : parse(Int, o)
        end

        for (post, elements) in zip(sp_posts, init_elements), u2 in last(rng), u1 in first(rng)
            Δu = [u1 - init_us[1], u2 - init_us[2]]
            newf = merge_formula(elements, unit, Δu; sign = (:+, :-))
            for fn in add_fn
                tbl.Abbreviation[i] = cls
                tbl.Species[i] = spid(cls, u1, u2, post)
                tbl.Formula[i] = newf
                tbl.mz[i] = fn(mw(newf))
                i += 1
            end
        end
    end
    tbl
end

library(class, adduct, range) = library(vectorize(class), vectorize(adduct), vectorize(range))

const LIBRARY_POS = reduce(append!, (
    library(Cer, ["[M+H]+"], [(32:46, 1:2, "O"), (32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
    library([HexCer, Hex2Cer, Hex3Cer, GM3], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
    library(HexNAcHex2Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "O2")),
    library(HexNAcHex3Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "O2")),
    library([GM2, GM1, GD3, GD2, GD1, GT3, GT2, GT1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, "O2")),
    library(GQ1, "[M+2H]2+", (32:46, 0:4, "O2")),
    library(GP1, ["[M+2H]2+", "[M+3H]3+"], (34:2:38, 0:4, "O2")),
))
SPDB[:LIBRARY_POS] = LIBRARY_POS

#LIBRARY_NEG

const FRAGMENT_POS = let
    frags = [lcb(18, 1, 1), lcb(18, 0, 2),
            lcb(16, 1, 2), lcb(17, 1, 2), lcb(18, 1, 2), lcb(18, 0, 3), lcb(19, 1, 2), lcb(20, 1, 2),
            lcb(16, 2, 2), lcb(18, 2, 2), lcb(18, 1, 3), lcb(20, 2, 2)]
    ion = map(default_adduct, frags)
    append!(ion, [Ion(ProtonationNL2H2O(), NeuAc()),
                    Ion(ProtonationNLH2O(), NeuAc()),
                    HexNAcHex,
                    HexNAcHexNANA,
                    HexNAcHexNANA2,
                    Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())),
                    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc())),
                    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())),
                    Ion(Protonation(), PhosphoCholine())
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
    ce = CSV.read(file, Table; select = 1:5)
    Table(
        ms1 = [t() for t in (eval ∘ Meta.parse).(ce.ms1)],
        adduct1 = [SPDB[:ADDUCTCODE].object[x] for x in ce.adduct1],
        ms2 = ce.ms2,
        adduct2 = map(ce.adduct2) do adduct
            @match adduct begin
                0 => :default
                _ => SPDB[:ADDUCTCODE].object[adduct]
            end
        end,
        eV = ce.eV
    )
end

const CE = read_ce(joinpath(@__DIR__, "..", "config", "CE.csv"))
SPDB[:CE] = CE