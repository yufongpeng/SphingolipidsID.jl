
"""
    set_db!(key::Symbol, db = eval(key); kwargs...)

Set `SPDB[key]` to `db`. If customized db is used for `:ADDUCTCODE`, keyword argument `NLH2O` must provided.
"""
function set_db!(key::Symbol, db = eval(key); kwargs...)
    SPDB[key] = db
    if key == :ADDUCTCODE
        v = db === ADDUCTCODE ? (1:4) : begin
            haskey(kwargs, :NLH2O) || throw(ArgumentError("Keyword argument NLH2O is required."))
            getindex(kwargs, :NLH2O)
        end
        SPDB[:NLH2O] = @view SPDB[:ADDUCTCODE].object[v]
    end
    db
end

"""
    read_adduct_code(file::String; 
                        silencewarnings::Bool = false,
                        maxwarnings::Int = 10, 
                        debug::Bool = false
                    )

Read config file for `SPDB[:ADDUCTCODE]`. The csv file should not contain header and the columns are the string representation and type of the adduct, respectively.
"""
function read_adduct_code(file::String; 
                            silencewarnings::Bool = false,
                            maxwarnings::Int = 10, 
                            debug::Bool = false
                        )
    code = CSV.read(file, Table; header = [:repr, :object], silencewarnings, maxwarnings, debug)
    Table(
        repr = code.repr,
        object = (@p code.object map((eval ∘ Meta.parse)(_)()))
    )
end
"""
    const ADDUCTCODE

Default `SPDB[:ADDUCTCODE]`.
"""
const ADDUCTCODE = read_adduct_code(joinpath(@__DIR__, "..", "config", "ADDUCTCODE.csv"))
set_db!(:ADDUCTCODE)

repr_adduct(x::Adduct) = SPDB[:ADDUCTCODE].repr[findfirst(==(x), SPDB[:ADDUCTCODE].object)]
object_adduct(x::AbstractString) = SPDB[:ADDUCTCODE].object[findfirst(==(x), SPDB[:ADDUCTCODE].repr)]

"""
    read_class_db(file::String;
                    col_name = :Abbreviation,
                    col_regex = :regex,
                    col_formula = [:formula, :u1, :u2],
                    col_unit = [:nu1, :nu2],
                    col_ions = [:default_ion, :parent_ion, :adduct_ion, :isf_ion],
                    silencewarnings::Bool = false,
                    maxwarnings::Int = 10, 
                    debug::Bool = false
                )

Read config file for `SPDB[:CLASSDB]`. The csv file should contain columns controlled by the following keyword arguments:
* `col_name`: the abbreviation of the class.
* `col_regex`: regular expression for mathching this class.
* `col_formula`: the first column is the formula of the base compound; the second and third one are the difference in formula when adding one to each units. For sphingolipids, this is always "CH2" and "-H2".
* `col_unit`: number of units(chain and double bonds) of base compound.
* `col_ions`: ions that will be generated. The first column is the default ions for generating MRM; the second one is parent ions which a compound must contain at least one; 
the third one is adduct ion which are all possible ions; the fourth one is isf ions generating from in-source fragmentation. Ions are represented by the index in `SPDB[:ADDUCTCODE]`.
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_class_db(file::String;
                        col_name = :Abbreviation,
                        col_regex = :regex,
                        col_formula = [:formula, :u1, :u2],
                        col_unit = [:nu1, :nu2],
                        col_ions = [:default_ion, :parent_ion, :adduct_ion, :isf_ion],
                        silencewarnings::Bool = false,
                        maxwarnings::Int = 10, 
                        debug::Bool = false
                    )

    class_db = CSV.read(file, Table; silencewarnings, maxwarnings, debug)
    parse_ion = x -> SPDB[:ADDUCTCODE].object[(eval ∘ Meta.parse)(x)]
    Table(
        #Abbreviation = map(x -> (eval ∘ Meta.parse)(x)(), getproperty(class_db, col_name)),
        Abbreviation = (eval ∘ Meta.parse).(getproperty(class_db, col_name)),
        regex = (eval ∘ Meta.parse).(getproperty(class_db, col_regex)),
        formula = (parse_compound ∘ String).(getproperty(class_db, col_formula[1])),
        nu1 = getproperty(class_db, col_unit[1]),
        u1 = parse_unit.(getproperty(class_db, col_formula[2])),
        nu2 = getproperty(class_db, col_unit[2]),
        u2 = parse_unit.(getproperty(class_db, col_formula[3])),
        default_ion = parse_ion.(getproperty(class_db, col_ions[1])),
        parent_ion = parse_ion.(getproperty(class_db, col_ions[2])),
        adduct_ion = parse_ion.(getproperty(class_db, col_ions[3])),
        isf_ion = parse_ion.(getproperty(class_db, col_ions[4]))
    )
end

"""
    const CLASSDB

Default `SPDB[:CLASSDB]`.
"""
const CLASSDB = read_class_db(joinpath(@__DIR__, "..", "config", "CLASSDB.csv"))
set_db!(:CLASSDB)

class_db_index(spb::LCB) = SPDB[:CLASSDB][findfirst(==(spb), SPDB[:CLASSDB].Abbreviation)]
#class_db_index(cls::T) where {T <: ClassSP} = SPDB[:CLASSDB][findfirst(==(deisomerized(cls)), SPDB[:CLASSDB].Abbreviation)]
class_db_index(cls::T) where {T <: ClassSP} = SPDB[:CLASSDB][findfirst(==(deisomerized(T)), SPDB[:CLASSDB].Abbreviation)]

"""
    library(class, adduct, range)

Create a library. 

* `class`: a type of `ClassSP` or a `Vector{<: ClassSP}`.
* `adduct`: `Adduct` or `Vector{<: Adduct}`.
* `range`: `Tuple` or `Vector{<: Tuple}`. Each tuple contains the range of carbon chain, double bonds, and addtional elements as `String`.
"""
library(class::Vector, adduct::Vector{<: Adduct}, range::Vector{<: Tuple}) = library(class, map(repr_adduct, adduct), range)
#library(class::Vector{DataTyple}, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple}) = library([cls() for cls in class], adduct, range)
function library(class::Vector, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple})
    range = map(range) do r
        vectorize.(r)
    end
    n = mapreduce(rng -> mapreduce(length, *, rng), +, range) * length(adduct) * length(class)
    tbl = Table(
                    Abbreviation = Vector{DataType}(undef, n),
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
            newf = merge_formula(elements, unit, Δu)
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

"""
    const LIBRARY_POS

Default `SPDB[:LIBRARY_POS]`.
"""
const LIBRARY_POS = reduce(append!, (
    library(Cer, ["[M+H]+"], [(32:46, 1:2, "O"), (32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
    library([HexCer, Hex2Cer, Hex3Cer, GM3], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
    library(HexNAcHex2Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "O2")),
    library(HexNAcHex3Cer, ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, "O2")),
    library([GM2, GM1, GD3, GD2, GD1, GT3, GT2, GT1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, "O2")),
    library(GQ1, "[M+2H]2+", (32:46, 0:4, "O2")),
    library(GP1, ["[M+2H]2+", "[M+3H]3+"], (34:2:38, 0:4, "O2")),
))

#LIBRARY_NEG
"""
    const FRAGMENT_POS

Default `SPDB[:FRAGMENT_POS]`.
"""
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

#FRAGMENT_NEG

"""
    const CONNECTION

Default `SPDB[:CONNECTION]`.
"""
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

"""
    read_ce(file::String; 
            col_ms = [:ms1, :ms2], 
            col_adduct = [:adduct1, :adduct2], 
            col_ev = :eV,
            silencewarnings::Bool = false,
            maxwarnings::Int = 10, 
            debug::Bool = false
            )

Read config file for `SPDB[:CE]`. The csv file should contain columns controlled by the following keyword arguments:
* `col_ms`: 2-element vector containing the type of ms1 and ms2, repectively.
* `col_adduct`: 2-element vector containing the index of adducts for ms1 and ms2, repectively. `0` indicates default value.
* `col_ev`: collision aenergy.
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_ce(file::String; 
                col_ms = [:ms1, :ms2],
                col_adduct = [:adduct1, :adduct2],
                col_ev = :eV,
                silencewarnings::Bool = false,
                maxwarnings::Int = 10, 
                debug::Bool = false
                )
    ce = CSV.read(file, Table; select = 1:5, silencewarnings, maxwarnings, debug)
    Table(
        ms1 = [t() for t in (eval ∘ Meta.parse).(getproperty(ce, col_ms[1]))],
        adduct1 = [SPDB[:ADDUCTCODE].object[x] for x in getproperty(ce, col_adduct[1])],
        ms2 = getproperty(ce, col_ms[2]),
        adduct2 = map(getproperty(ce, col_adduct[2])) do adduct
            @match adduct begin
                0 => :default
                _ => SPDB[:ADDUCTCODE].object[adduct]
            end
        end,
        eV = getproperty(ce, col_ev)
    )
end
"""
    const CE

Default `SPDB[:CE]`.
"""
const CE = read_ce(joinpath(@__DIR__, "..", "config", "CE.csv"))

set_db!.([:LIBRARY_POS, :FRAGMENT_POS, :CONNECTION, :CE])