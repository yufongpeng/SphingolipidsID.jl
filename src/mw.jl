# molecular weight to be coordinated with https://www.sisweb.com/referenc/tools/exactmass.htm? and LipidMaps
const MW = Dict(
    ""   => 0u"g",
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

mz(cpd::CompoundSP, ion::ISF) = parse_adduct(ion.adduct)(mw(CompoundSP(ion.molecule, cpd.sum, cpd.chain)))
mz(cpd::CompoundSP, ions::AcylIon{T}) where T = mz.(ions.ions, (cpd.sum .- sumcomp(T()))...)
mz(cpd::CompoundSP, ion::Ion) = mz(ion)
mz(cpd::CompoundSP, add::Adduct) = parse_adduct(add)(mw(cpd))
mz(cpd::CompoundSP, ion::Ion{<: Adduct, <: ClassSP}) = parse_adduct(ion.adduct)(mw(CompoundSP(ion.molecule, cpd.sum, nothing)))
mz(ion::Union{Ion, ISF}) = parse_adduct(ion.adduct)(mw(ion.molecule))
mz(ion::Ion, cb, db, o) = parse_adduct(ion.adduct)(mw(ion.molecule, cb, db, o))
mw(::T, cb, db, o = N) where {N, T <: ACYL{N}} = mw("C") * cb + mw("H") * (2 * cb - 2 * db - 1) + mw("O") * (o + 1)
mw(molecule, cb, db, o) = mw(molecule)
#hydrosyn(a, b) = a + b - nmw("H2O")
mw(::Hex) = mw("C6H12O6") 
mw(::HexNAc) = mw("C8H15NO6") 
mw(::NeuAc) = mw("C11H19NO9") 

mw(glycan::Glycan) = mapreduce(mw, +, glycan.chain) - (length(glycan.chain) - 1) * mw("H2O")

mw(lcb::LCB{N, C}) where {N, C} = C * mw("C") + mw("O") * N + (2 * (C - nunsa(lcb) + N) + 3) * mw("H") + mw("N")

function mw(cpd::CompoundSP)
    cls = class_db_index(cpd.class)
    unit = cls[[:u1, :u2]]
    init_us = cls[[:nu1, :nu2]]
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
        num = @match el begin
            "C" => num + Δcn
            "H" => num - Δdb * 2
        end
        num == 1 ? "$el" : "$el$num"
    end
end

function merge_formula(elements::Vector, units, Δu; sign = ntuple(i -> :+, length(init_unit)))
    adds = mapmany(units, Δu, sign) do unit, nu, s
        map(unit) do (el, num)
            el => eval(s)(1) * num  * nu
        end
    end

    mapreduce(*, elements) do element
        el, num = element
        num += sum(num_add for (el_add, num_add) in adds if el_add == el; init = 0)
        num == 1 ? "$el" : "$el$num"
    end
end

function add_addition(init_elements, posts)
    init_elements = deepcopy(init_elements)
    ido = findfirst(x -> x.first == "O", init_elements)
    for post in posts
        m = match(r"(\d*)O(\d*)", post)
        isnothing(m) && continue
        add = @match m.captures begin
            ["", ""] => 1
            ["", s]  => parse(Int, s)
            [s, _]   => parse(Int, s)
        end
        init_elements[ido] = init_elements[ido].first => init_elements[ido].second + add
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
        # M-H / M+FA-H
        mapreduce(nmw, -, split(pos_add, "-"))
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
