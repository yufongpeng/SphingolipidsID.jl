
function read_adduct_code(file::String) 
    code = CSV.read(file, DataFrame, header = [:repr, :object])
    transform!(code, :object => ByRow(eval ∘ Meta.parse), renamecols = false)
end
const ADDUCTCODE = read_adduct_code(joinpath(@__DIR__, "..", "config", "ADDUCTCODE.csv"))
const ADDUCTS = ADDUCTCODE.object
const NLH2O = @view ADDUCTS[1:4]

repr_adduct(x::Adduct) = ADDUCTCODE.repr[findfirst(==(x), ADDUCTS)]
object_adduct(x::AbstractString) = ADDUCTCODE.object[findfirst(==(x), ADDUCTCODE.repr)]

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
    transform!(class_db, col_ions .=> ByRow(x -> ADDUCTS[(eval ∘ Meta.parse)(x)]), renamecols = false)
    rename!(class_db, [col_name => :Abbreviation, col_regex => :regex, (col_formula .=> [:formula, :u1, :u2])..., (col_ions .=> [:default_ion, :parent_ion, :adduct_ion, :isf_ion])...])
    rename!(class_db, col_unit .=> ["#us", "#u2"])
    class_db
end

const CLASSDB = read_class_db(joinpath(@__DIR__, "..", "config", "config\\CLASSDB.csv"))

class_db_index(::LCB) = @views CLASSDB[findfirst(==(SPB), CLASSDB[!, :Abbreviation]), :]
class_db_index(cls::T) where {T <: ClassGSL} = @views CLASSDB[findfirst(==(deisomerized(T)), CLASSDB[!, :Abbreviation]), :]

library(class::Vector, adduct::Vector{<: Adduct}, range::Vector{<: Tuple}) = library(class, map(repr_adduct, adduct), range)
#library(class::Vector{Type{<: ClassGSL}}, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple}) = library([cls() for cls in class], adduct, range)

function library(class::Vector, adduct::Vector{<: AbstractString}, range::Vector{<: Tuple})

    length(range) != 1 && length(range != length(class)) &&
        throw(ArgumentError("Number of ranges should be one or as same as number of classes"))
    length(range) == 1 && (range = repeat(range, length(class)))  

    n_within = isempty(length(range)) ? repeat([length(adduct)], length(class)) : map(range) do rng
        mapreduce(length, *, rng) * length(adduct)
    end

    n = reduce(+, n_within)
    df = DataFrame(
                    :Abbreviation => Vector{Type{<: ClassGSL}}(undef, n),
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
    for (cls, rng) in zip(class, range)            
        # specific for species
        # match cls with Abbreviation first
        id_compound = findfirst(==(cls), CLASSDB[!, :Abbreviation])
        # isnothing(id_compound) && (id_compound = findfirst(abbr -> match(Regex(cls * "-.*"), abbr), CLASSDB[!, :Abbreviation]))
        init_elements = CLASSDB[id_compound, :formula]
        unit = CLASSDB[id_compound, [:u1, :u2]]
        init_us = CLASSDB[id_compound, ["#u1", "#u2"]]

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
            post_sep = [split(CLASSDB[id_compound, :regex].pattern, ";")[2:end]]
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

const DEFAULT_POS_LIBRARY = reduce(append!, (
    library([Cer], ["[M+H]+"], (32:46, 1:2, ["O"])),
    library([Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
    library([HexCer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([HexCer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
    library([Hex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([Hex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
    library([Hex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([Hex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
    library([HexNAcHex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([HexNAcHex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([GM3], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
    library([GM3], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
    library([GM2], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GM1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GD3], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GD2], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GD1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GT3], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GT2], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GT1], ["[M+H]+", "[M+2H]2+"], (32:46, 0:4, ["2O"])),
    library([GQ1], ["[M+2H]2+", "[M+3H]3+"], (32:46, 0:4, ["2O"])),
    library([GP1], ["[M+2H]2+", "[M+3H]3+"], (34:2:38, 0:4, ["2O"])),
))

const DEFAULT_POS_LIBRARY_NANA = @view DEFAULT_POS_LIBRARY[findfirst(==(GM3), DEFAULT_POS_LIBRARY[!, 1]):end, :]

#DEFAULT_NEG_LIBRARY

const DEFAULT_POS_FRAGMENT = let
    frags = [SPB2{18, 1}(), SPB2{18, 2}(), 
            SPB3{16, 2}(), SPB3{17, 2}(), SPB3{18, 2}(), SPB3{18, 3}(), SPB3{19, 2}(), SPB3{20, 2}(), 
            SPB4{16, 2}(), SPB4{18, 2}(), SPB4{18, 3}(), SPB4{20, 2}()]
    add = map(default_adduct, frags)
    append!(add, [Ion(ProtonationNL2H2O(), NeuAc()), Ion(ProtonationNLH2O(), NeuAc()), HexNAcHex, HexNAcHexNANA, HexNAcHexNANA2])
    [add map(mz, add)]
end

#DEFAULT_NEG_FRAGMENT