# molecular weight to be coordinated with https://www.sisweb.com/referenc/tools/exactmass.htm? and LipidMaps
const MW = Dict(
    "C"  => 12u"g",
    "H"  => 1.007825u"g",
    "O"  => 15.994915u"g",
    "N"  => 14.003074u"g",
    "P"  => 30.973763u"g",
    "S"  => 31.972072u"g",
    "Li" => 7.016005u"g",
    "Na" => 22.989770u"g",
    "K"  => 38.963708u"g",
    "Cl" => 34.968853u"g",
    "Ag" => 106.905095u"g"
)

mw(formula::AbstractString) = mw(string(formula))
mw(formula::String) = mw(parse_compound(formula))

function mw(elements::Vector)
    # Vector of el => #el
    weight = 0.0u"g"
    for (el, n) in elements
        weight += MW[string(el)] * n
    end

    ustrip(weight)
end

function nmw(formula::AbstractString)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        mw(formula)
    else
        formula = replace(formula, n.match => "")
        mw(formula) * parse(Int, n.match)
    end
end

mz(cpd::CompoundGSL, ion::ISF) = parse_adduct(ion.adduct)(mw(CompoundGSL(ion.molecule, cpd.sum, cpd.chain, DataFrame(), 0, [], nothing)))
function mz(cpd::CompoundGSL, ions::AcylIon{T}) where T
    acyl_cb, acyl_db, acyl_o = cpd.sum .- sumcomp(T())
    map(ions.ions) do ion
        @match ion begin
            Ion(_, ::LCB)   => mz(ion)
            Ion(_, ::ACYL)  => mz(ion, acyl_cb, acyl_db, acyl_o)
        end
    end
end
mz(cpd::CompoundGSL, ion::Ion) = mz(ion)
mz(cpd::CompoundGSL, add::Adduct) = parse_adduct(add)(mw(cpd))
mz(cpd::CompoundGSL, ion::Ion{<: Adduct, <: ClassGSL}) = parse_adduct(ion.adduct)(mw(CompoundGSL(ion.molecule, cpd.sum, nothing, DataFrame(), 0, [], nothing)))
mz(ion::Union{Ion, ISF}) = parse_adduct(ion.adduct)(mw(ion.molecule))
mz(ion::Ion{S, <: ACYL{N}}, cb, db, o = N) where {S, N} = parse_adduct(ion.adduct)(mw(ion.molecule, cb, db, o))
mw(::T, cb, db, o = N) where {N, T <: ACYL{N}} = mw("C") * cb + mw("H") * (2 * cb - 2 * db - 1) + mw("O") * (o + 1)
#hydrosyn(a, b) = a + b - nmw("H2O")
mw(::Hex) = mw("C6H12O6") 
mw(::HexNAc) = mw("C8H15NO6") 
mw(::NeuAc) = mw("C11H19NO9") 

mw(glycan::Glycan) = mapreduce(mw, +, glycan.chain) - (length(glycan.chain) - 1) * mw("H2O")

function mw(lcb::LCB{N, C}) where {N, C}
    cls = class_db_index(lcb)
    unit = cls[[:u1, :u2]]
    init_us = cls[["#u1", "#u2"]]
    Δu = [C - init_us[1], nunsa(lcb) - N - init_us[2]]
    ms = mw(merge_formula(cls.formula, unit, Δu; sign = (:+, :-)))
    ms + mw("O") * N
end

function mw(cpd::CompoundGSL)
    cls = class_db_index(cpd.class)
    unit = cls[[:u1, :u2]]
    init_us = cls[["#u1", "#u2"]]
    Δu = [cpd.sum[1] - init_us[1], cpd.sum[2] - init_us[2]]
    ms = mw(merge_formula(cls.formula, unit, Δu; sign = (:+, :-)))
    ms + mw("O") * cpd.sum[3]
end

merge_species(cn::Int, db::Int, class::AbstractString, pre::AbstractString, post::AbstractString) = 
    *(class, " ", isempty(pre) ? pre : (pre * "-"), "$cn:$db", isempty(post) ? post : (";" * post))


function merge_formula(elements::Vector, Δcn::Int, Δdb::Int) 
    Δdb -= Δcn
    mapreduce(*, elements) do element
        el, num = element
        if el == "C"
            num = num + Δcn
        elseif el == "H"
            num = num - Δdb * 2
        end
        num == 1 ? "$el" : "$el$num"
    end
end

function merge_formula(elements::Vector, units, Δu; sign = ntuple(i -> :+, length(init_unit)))
    adds = Pair[]
    for (id, (unit, nu)) in enumerate(zip(units, Δu))
        if isnothing(unit)
            break
        else
            append!(adds, map(unit) do (el, num)
                el => eval(sign[id])(1) * num  * nu
            end)
        end
    end

    if isempty(adds) 
        mapreduce(*, elements) do element
            el, num = element
            num == 1 ? "$el" : "$el$num"
        end
    else
        mapreduce(*, elements) do element
            el, num = element
            for (el_add, num_add) in adds
                el_add == el && (num += num_add)
            end
            num == 1 ? "$el" : "$el$num"
        end
    end
end

function add_addition(init_elements, posts)
    init_elements = deepcopy(init_elements)
    for post in posts
        m = match(r"(\d*)O(\d*)", post)
        if !isnothing(match)
            if all(==(""), m.captures)
                add = 1
            elseif first(m.captures) == ""
                add = parse(Int, last(m.captures))
            else
                add = parse(Int, first(m.captures))
            end
            for (id, (el, num)) in enumerate(init_elements)
                if el == "O"
                    init_elements[id] = el => num + add
                end
            end
        end
    end
    init_elements
end

function parse_adduct(adduct::AbstractString)
    adduct = replace(adduct, "FA-H" => "HCOO")
    adduct = replace(adduct, "OAc" => "CH3COO")
    ion, charge = split(adduct, "]")
    charge = split(charge, "+", keepempty = false)
    charge = isempty(charge) ? 1 : begin
        charge = split(only(charge), "-", keepempty = false)
        isempty(charge) ? 1 : parse(Int, only(charge))
    end
    ion = replace(ion, "[" => "")
    # #compounds
    nm, ion = split(ion, "M")
    nm = isempty(nm) ? 1 : parse(Int, nm)

    pos_adds = split(ion, "+", keepempty = false)
    madd = mapreduce(+, pos_adds) do pos_add
        neg_adds = split(pos_add, "-")
        # M-H / M+FA-H
        isempty(first(neg_adds)) ? (- mapreduce(nmw, +, neg_adds[2:end])) : mapreduce(nmw, -, neg_adds)
    end

    x -> (x * nm + madd) / charge 
end

parse_adduct(::Protonation) = (x -> (x + nmw("H")))
parse_adduct(::ProtonationNLH2O) = (x -> (x - nmw("OH")))
parse_adduct(::ProtonationNL2H2O) = (x -> (x - nmw("H3O2")))
parse_adduct(::ProtonationNL3H2O) = (x -> (x - nmw("H5O3")))
parse_adduct(::DiProtonation) = (x -> (x / 2 + nmw("H")))
parse_adduct(::TriProtonation) = (x -> (x / 3 + nmw("H")))
parse_adduct(::Deprotonation) = (x -> (x - nmw("H")))
parse_adduct(::DeprotonationNLH2O) = (x -> (x - nmw("H3O")))
parse_adduct(::DiDeprotonation) = (x -> (x / 2 - nmw("H")))
parse_adduct(::TriDeprotonation) = (x -> (x / 3 - nmw("H")))
parse_adduct(::AddOAc) = (x -> (x + nmw("C2H3O2")))
parse_adduct(::AddHCOO) = (x -> (x + nmw("CHO2")))
parse_adduct(::LossCH2O) = (x -> (x - nmw("CH2O")))
parse_adduct(::AddO) = (x -> (x + nmw("O")))
parse_adduct(::AddC2H2O) = (x -> (x + nmw("C2H2O")))
parse_adduct(::LossCH8NO) = (x -> (x - nmw("CH8NO")))
parse_adduct(::LossC2H8NO) = (x -> (x - nmw("C2H8NO")))
parse_adduct(::AddC3H5NO) = (x -> (x + nmw("C3H5NO")))
parse_adduct(::AddC2H5NO) = (x -> (x + nmw("C2H5NO")))
