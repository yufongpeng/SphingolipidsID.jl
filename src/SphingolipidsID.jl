module SphingolipidsID

using DataFrames, CSV, MLStyle, Statistics, PrettyTables
using UnitfulMoles: parse_compound, ustrip, @u_str
export library, rule, featuretable_mzmine, add_ce_mzmine!, featuretable_masshunter_mrm, filter_duplicate, 
        preis, preis!, finish_profile!, 
        apply_rules!, id_product, query,
        generate_mrm,

        SPDB, 
        
        LIBRARY_POS, FRAGMENT_POS, ADDUCTCODE, CLASSDB, 

        read_adduct_code, read_class_db, read_ce, class_db_index, 

        mw, mz, 

        ClassGSL, SPB, Cer, CerP, HexCer, SHexCer, Hex2Cer, SHexHexCer, Hex3Cer, 
        ClassHexNAcHex2Cer, HexNAcHex2Cer, ClassHexNAcHex3Cer, HexNAcHex3Cer, HexNAc_Hex3Cer, Hex_HexNAc_Hex2Cer, ClasswoNANA,
        GM4, GM3, GM2, ClassGM1, GM1, GM1a, GM1b,
        GD3, GD2, ClassGD1, GD1, GD1a, GD1b, GD1c, GD1α, 
        GT3, GT2, GT1, ClassGT1, GT1, GT1a, GT1b, GT1c, GT1aα,
        ClassGQ1, GQ1, GQ1b, GQ1c, GQ1bα,
        ClassGP1, GP1, GP1c, GP1cα,

        LCB, LCB2, LCB3, LCB4, SPB2, SPB3, SPB4, PhytoSPB, PhytoSPB3, PhytoSPB4,
        ACYL, Acyl, Acylα, Acylβ, 
        Chain, 

        Sugar, Glycan, Hex, HexNAc, NeuAc,

        Adduct, Pos, Neg, 
        Ion, ISF,
        Protonation, ProtonationNLH2O, ProtonationNL2H2O, ProtonationNL3H2O, DiProtonation, TriProtonation, 
        AddNH4, AddHNH4, Add2NH4, Sodization, SodizationProtonation, DiSodization, 
        Deprotonation, DeprotonationNLH2O, DiDeprotonation, TriDeprotonation, AddOAc, AddHCOO, 
        LossCH2O, AddO, AddC2H2O, LossCH8NO, LossC2H8NO, AddC3H5NO, AddC2H5NO,

        CompoundGSL, AnalyteGSL, 
        IonPlus, IonComparison, IonMode, IonUnion, AcylIon,

        Data, PreIS, MRM, Project, Query, not

import Base: show, print, isless, isempty, keys, length, union, union!, deleteat!, 
        iterate, getindex, view, firstindex, lastindex, sort, sort!, push!, pop!, popat!, popfirst!, reverse, reverse!

abstract type ClassGSL end

struct SPB <: ClassGSL end
struct Cer <: ClassGSL end
struct CerP <: ClassGSL end
struct HexCer <: ClassGSL end
struct SHexCer <: ClassGSL end
struct SHexHexCer <: ClassGSL end
struct Hex2Cer <: ClassGSL end
struct Hex3Cer <: ClassGSL end

abstract type ClassHexNAcHex2Cer <: ClassGSL end
struct HexNAcHex2Cer <: ClassHexNAcHex2Cer end

abstract type ClassHexNAcHex3Cer <: ClassGSL end
struct HexNAcHex3Cer <: ClassHexNAcHex3Cer 
    isomer
end
struct HexNAc_Hex3Cer <: ClassHexNAcHex3Cer end
struct Hex_HexNAc_Hex2Cer <: ClassHexNAcHex3Cer end

const ClasswoNANA = Union{SPB, Cer, CerP, HexCer, SHexCer, Hex2Cer, SHexHexCer, Hex3Cer, <: ClassHexNAcHex2Cer, <: ClassHexNAcHex3Cer}

struct GM4 <: ClassGSL end
struct GM3 <: ClassGSL end
struct GM2 <: ClassGSL end

abstract type ClassGM1 <: ClassGSL end
struct GM1 <: ClassGM1 
    isomer
end
struct GM1a <: ClassGM1 end
struct GM1b <: ClassGM1 end

struct GD3 <: ClassGSL end
struct GD2 <: ClassGSL end

abstract type ClassGD1 <: ClassGSL end
struct GD1 <: ClassGD1
    isomer
end
struct GD1a <: ClassGD1 end
struct GD1b <: ClassGD1 end
struct GD1c <: ClassGD1 end
struct GD1α <: ClassGD1 end

struct GT3 <: ClassGSL end
struct GT2 <: ClassGSL end

abstract type ClassGT1 <: ClassGSL end
struct GT1 <: ClassGT1
    isomer
end
struct GT1a <: ClassGT1 end
struct GT1aα <: ClassGT1 end
struct GT1b <: ClassGT1 end
struct GT1c <: ClassGT1 end

abstract type ClassGQ1 <: ClassGSL end
struct GQ1 <: ClassGQ1
    isomer
end
struct GQ1b <: ClassGQ1 end
struct GQ1bα <: ClassGQ1 end
struct GQ1c <: ClassGQ1 end

abstract type ClassGP1 <: ClassGSL end
struct GP1 <: ClassGP1    
    isomer
end
struct GP1c <: ClassGP1 end
struct GP1cα <: ClassGP1 end

for (class, super) in zip((:HexNAcHex3Cer, :GM1, :GD1, :GT1, :GQ1, :GP1), (:ClassHexNAcHex3Cer, :ClassGM1, :ClassGD1, :ClassGT1, :ClassGQ1, :ClassGP1))
    @eval begin 
        $class(cls::Vararg{<: $super}) = $class(cls)
    end
end

# N = # OH
abstract type LCB{N, C} end
abstract type LCB4{N, C} <: LCB{N, C} end
abstract type LCB3{N, C} <: LCB{N, C} end
abstract type LCB2{N, C} <: LCB{N, C} end
struct SPB4{N, C} <: LCB4{N, C} end
struct PhytoSPB4{N, C} <: LCB4{N, C} end
struct SPB3{N, C} <: LCB3{N, C} end
struct PhytoSPB3{N, C} <: LCB3{N, C} end
struct SPB2{N, C} <: LCB2{N, C} end
const PhytoSPB{N, C} = Union{PhytoSPB3{N, C}, PhytoSPB4{N, C}}

abstract type ACYL{N} end
struct Acyl{N} <: ACYL{N} end
struct Acylα{N} <: ACYL{N} end
struct Acylβ{N} <: ACYL{N} end

struct Chain{S <: LCB, T <: ACYL}
    lcb::S
    acyl::T
end

@as_record Chain
@as_record LCB
@as_record ACYL

abstract type Sugar end

struct Glycan{T} 
    chain::T
end

Glycan(sugar::Sugar...) = Glycan(sugar)

struct NeuAc <: Sugar end
struct Hex <: Sugar end
struct HexNAc <: Sugar end

abstract type Adduct end
abstract type Pos <: Adduct end
abstract type Neg <: Adduct end

struct Ion{S <: Adduct, T <: Union{<: Sugar, <: Glycan, <: ClassGSL, <: LCB, <: ACYL}}
    adduct::S
    molecule::T
end

struct ISF{S <: Adduct, T <: ClassGSL}
    adduct::S
    molecule::T
end

struct Protonation <: Pos end
struct ProtonationNLH2O <: Pos end
struct ProtonationNL2H2O <: Pos end
struct ProtonationNL3H2O <: Pos end
struct DiProtonation <: Pos end
struct TriProtonation <: Pos end
struct AddNH4 <: Pos end
struct AddHNH4 <: Pos end
struct Add2NH4 <: Pos end
struct Sodization <: Pos end
struct SodizationProtonation <: Pos end
struct DiSodization <: Pos end
struct Deprotonation <: Neg end
struct DeprotonationNLH2O <: Neg end
struct DiDeprotonation <: Neg end
struct TriDeprotonation <: Neg end
struct AddOAc <: Neg end
struct AddHCOO <: Neg end


struct LossCH2O <: Neg end
struct AddO <: Neg end
struct AddC2H2O <: Neg end
struct LossCH8NO <: Neg end
struct LossC2H8NO <: Neg end
struct AddC3H5NO <: Neg end
struct AddC2H5NO <: Neg end

const NANA = (Ion(ProtonationNL2H2O(), NeuAc()), Ion(DeprotonationNLH2O(), NeuAc()), Ion(ProtonationNLH2O(), NeuAc()))
const NANA1 = Ion(DeprotonationNLH2O(), NeuAc())
const NANA2 = Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc()))
const NANA3 = Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc(), NeuAc()))
const NANA23 = (Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc())), Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc(), NeuAc())))
const HexNAcHex = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex()))
const HexNAcHexNANA = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex(), NeuAc()))
const HexNAcHexNANA2 = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex(), NeuAc(), NeuAc()))

@as_record Ion
@as_record Adduct

mutable struct CompoundGSL
    class::ClassGSL
    sum::NTuple{3, Int}   # (#C, #db, #OH)
    chain::Union{Chain, Nothing}
    fragments::DataFrame # ion1, ion2, source, id
    area::Float64
    states::Vector
    project
end

mutable struct AnalyteGSL
    compounds::Vector{CompoundGSL}
    rt::Float64
    states::Vector{Int}
end
# ==()

struct IonPlus{T}
    ions::T
end

struct IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}}
    ions::Tuple{S, T}
end

IonComparison(ion1, ion2) = IonComparison((ion1, ion2))

struct IonUnion
    rules
end

struct IonMode
    rules
end

struct Hydroxyl 
    spb::LCB
end

struct AcylIon{T <: LCB}
    ions
end

AcylIon{T}(ions...) where {T<: LCB} = AcylIon{T}(ions)

for fn in (:IonPlus, :IonUnion, :IonMode)
    @eval begin
        $fn(x...) = $fn(x)
    end
end

@as_record IonPlus

### connect db
abstract type Data end

# precursor/product, ion => name, ms/mz => m/z, mw => molecular weight
struct PreIS <: Data
    raw::DataFrame #id, mz1, scan(index of mz2), area, ev, rt
    range::Vector{Tuple{Float64, Float64}} 
    mz2::Vector{Float64}
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end

struct MRM <: Data
    raw::DataFrame #id, mz1, mz2, area, ev, rt
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end


struct Project
    analytes::Vector{AnalyteGSL}
    data::Vector{Data}
    anion
end

mutable struct Query
    project::Project
    result::AbstractVector
    query::Vector
    view::Bool
end

struct Inv
    args
end

not(id::Symbol, args...) = Inv((id, args...))

include("io.jl")
include("utils.jl")
include("mw.jl")
include("library.jl")
include("rule.jl")
include("data.jl")
include("preis.jl")
include("ruleID.jl")
include("query.jl")

# S: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-237.2224
# DHS: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-239.2224
# Phyto(4-OH): Acyl + C3H5NO-422.4004, SPB;3O - CH8NO-267.2329, Acyl + C2H5NO-410.4004, SPB;3O - C2H8NO-255.233
# α-hydroxyl: Acyl - CH2O-337.3476, Acyl + O-383.3531
# β-hydroxyl: SPB;2O + C2H2O - 340.246


end
