module SphingolipidsID

using CSV, MLStyle, Statistics, PrettyTables, DataPipes, SplitApplyCombine, TypedTables
using DataFrames: DataFrame, groupby, combine, All
using UnitfulMoles: parse_compound, ustrip, @u_str
export SPDB, LIBRARY_POS, FRAGMENT_POS, ADDUCTCODE, CLASSDB, 

        read_adduct_code, read_class_db, read_ce, class_db_index, 

        mw, mz, 

        ClassSP, SPB, Cer, CerP, HexCer, SHexCer, Hex2Cer, SHexHexCer, Hex3Cer, 
        HexNAcHex2Cer, HexNAcHex2Cer_, HexNAc_Hex2Cer, HexNAcHex3Cer, HexNAcHex3Cer_, HexNAc_Hex3Cer, Hex_HexNAc_Hex2Cer, ClasswoNANA,
        GM4, GM3, GM2, GM1, GM1_, GM1a, GM1b,
        GD3, GD2, GD1, GD1_, GD1a, GD1b, GD1c, GD1α, 
        GT3, GT2, GT1, GT1_, GT1a, GT1b, GT1c, GT1aα,
        GQ1, GQ1_, GQ1b, GQ1c, GQ1bα,
        GP1, GP1_, GP1c, GP1cα,
        CLS,

        LCB, LCB2, LCB3, LCB4, SPB2, SPB3, SPB4, NotPhyto, NotPhyto3, NotPhyto4,
        ACYL, Acyl, Acylα, Acylβ, 
        Lcb, Nacyl, Nacylα, Nacylβ, NACYL,
        Chain, 

        Sugar, Glycan, Hex, HexNAc, NeuAc,

        Adduct, Pos, Neg, 
        Ion, ISF,
        Protonation, ProtonationNLH2O, ProtonationNL2H2O, ProtonationNL3H2O, DiProtonation, TriProtonation, 
        AddNH4, AddHNH4, Add2NH4, Sodization, SodizationProtonation, DiSodization, 
        Deprotonation, DeprotonationNLH2O, DiDeprotonation, TriDeprotonation, AddOAc, AddHCOO, 
        LossCH2O, AddO, AddC2H2O, LossCH8NO, LossC2H8NO, AddC3H5NO, AddC2H5NO,

        CompoundSP, AnalyteSP, 
        #IonPlus, IonComparison, RuleMode, RuleUnion, AcylIon,

        Data, PreIS, MRM, Project, Query, 

        library, rule, 
        
        featuretable_mzmine, fill_ce_mzmine!, featuretable_masshunter_mrm, filter_duplicate!, rsd, re,

        preis, preis!, finish_profile!, 

        apply_rules!, id_product, 

        query, not, cpd, lcb, acyl, acylα, acylβ,

        nfrags, @score, filter_score!, apply_score, apply_score!, apply_threshold, apply_threshold!, select_score!,

        generate_mrm, nMRM, write_mrm



import Base: show, print, isless, isempty, keys, length, union, union!, deleteat!, 
        iterate, getindex, view, firstindex, lastindex, sort, sort!, push!, pop!, popat!, popfirst!, reverse, reverse!

const SPDB = Dict{Symbol, Any}()
const ISOMER = Dict{Type, Tuple}()
abstract type ClassSP end

struct SPB <: ClassSP end
struct Cer <: ClassSP end
struct CerP <: ClassSP end
struct HexCer <: ClassSP end
struct SHexCer <: ClassSP end
struct SHexHexCer <: ClassSP end
struct Hex2Cer <: ClassSP end
struct Hex3Cer <: ClassSP end

abstract type HexNAcHex2Cer <: ClassSP end
struct HexNAcHex2Cer_ <: HexNAcHex2Cer 
    isomer
end
struct HexNAc_Hex2Cer <: HexNAcHex2Cer end
push!(ISOMER, HexNAcHex2Cer_ => (HexNAc_Hex2Cer(),))

abstract type HexNAcHex3Cer <: ClassSP end
struct HexNAcHex3Cer_ <: HexNAcHex3Cer 
    isomer
end
struct HexNAc_Hex3Cer <: HexNAcHex3Cer end
struct Hex_HexNAc_Hex2Cer <: HexNAcHex3Cer end
push!(ISOMER, HexNAcHex3Cer_ => (HexNAc_Hex3Cer(), Hex_HexNAc_Hex2Cer()))

struct GM4 <: ClassSP end
struct GM3 <: ClassSP end
struct GM2 <: ClassSP end

abstract type GM1 <: ClassSP end
struct GM1_ <: GM1 
    isomer
end
struct GM1a <: GM1 end
struct GM1b <: GM1 end
push!(ISOMER, GM1_ => (GM1a(), GM1b()))

struct GD3 <: ClassSP end
struct GD2 <: ClassSP end

abstract type GD1 <: ClassSP end
struct GD1_ <: GD1
    isomer
end
struct GD1a <: GD1 end
struct GD1b <: GD1 end
struct GD1c <: GD1 end
struct GD1α <: GD1 end
push!(ISOMER, GD1_ => (GD1a(), GD1b(), GD1c(), GD1α()))


struct GT3 <: ClassSP end
struct GT2 <: ClassSP end

abstract type GT1 <: ClassSP end
struct GT1_ <: GT1
    isomer
end
struct GT1a <: GT1 end
struct GT1aα <: GT1 end
struct GT1b <: GT1 end
struct GT1c <: GT1 end
push!(ISOMER, GT1_ => (GT1a(), GT1b(), GT1c(), GT1aα()))

abstract type GQ1 <: ClassSP end
struct GQ1_ <: GQ1
    isomer
end
struct GQ1b <: GQ1 end
struct GQ1bα <: GQ1 end
struct GQ1c <: GQ1 end
push!(ISOMER, GQ1_ => (GQ1b(), GQ1c(), GQ1bα()))

abstract type GP1 <: ClassSP end
struct GP1_ <: GP1    
    isomer
end
struct GP1c <: GP1 end
struct GP1cα <: GP1 end
push!(ISOMER, GP1_ => (GP1c(), GP1cα()))
SPDB[:ISOMER] = ISOMER

for (class, super) in zip((:HexNAcHex2Cer_, :HexNAcHex3Cer_, :GM1_, :GD1_, :GT1_, :GQ1_, :GP1_), (:HexNAcHex2Cer, :HexNAcHex3Cer, :GM1, :GD1, :GT1, :GQ1, :GP1))
    @eval begin 
        $class(cls::Vararg{<: $super}) = $class(cls)
        $class() = $class(SPDB[:ISOMER][$class])
        $super(cls::Vararg{<: $super}) = $class(cls)
        $super() = $class(SPDB[:ISOMER][$class])

    end
end

const CLS = (
    nana = Dict{Int, Type}(
        0 => Union{SPB, Cer, CerP, HexCer, SHexCer, Hex2Cer, SHexHexCer, Hex3Cer, HexNAcHex2Cer, HexNAcHex3Cer},
        1 => Union{GM4, GM3, GM2, GM1},
        2 => Union{GD3, GD2, GD1},
        3 => Union{GT3, GT2, GT1},
        4 => GQ1,
        5 => GP1
    ),
    series = (
        as = Union{Hex2Cer, HexNAc_Hex2Cer, Hex_HexNAc_Hex2Cer, GM1b, GD1c, GD1α},
        a  = Union{GM3, GM2, GM1a, GD1a, GT1a, GT1aα},
        b  = Union{GD3, GD2, GD1b, GT1b, GQ1b, GQ1bα},
        c  = Union{GT3, GT2, GT1c, GQ1c, GP1c, GP1cα}
    ),
    fg = (
        nana = Union{GM4, GM3, GM2, GM1, GD3, GD2, GD1, GT3, GT2, GT1, GQ1, GP1},
        sulfate = Union{SHexCer, SHexHexCer}
    ) 
)

# N = # OH
abstract type LCB{N, C} end
abstract type LCB4{N, C} <: LCB{N, C} end
abstract type LCB3{N, C} <: LCB{N, C} end
abstract type LCB2{N, C} <: LCB{N, C} end
struct SPB4{N, C} <: LCB4{N, C} end
struct NotPhyto4{N, C} <: LCB4{N, C} end
struct SPB3{N, C} <: LCB3{N, C} end
struct NotPhyto3{N, C} <: LCB3{N, C} end
struct SPB2{N, C} <: LCB2{N, C} end
const NotPhyto{N, C} = Union{NotPhyto3{N, C}, NotPhyto4{N, C}}

abstract type ACYL{N} end
struct Acyl{N} <: ACYL{N} end
struct Acylα{N} <: ACYL{N} end
struct Acylβ{N} <: ACYL{N} end

# Used for query
struct Lcb 
    cb::Int
    db::Int
    ox::Int
end
abstract type NACYL end
struct Nacyl <: NACYL
    cb::Int
    db::Int
    ox::Int
end
struct Nacylα <: NACYL
    cb::Int
    db::Int
    ox::Int
end
struct Nacylβ <: NACYL
    cb::Int
    db::Int
    ox::Int
end

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

struct Ion{S <: Adduct, T <: Union{Sugar, Glycan, ClassSP, LCB, ACYL}}
    adduct::S
    molecule::T
end

struct ISF{S <: Adduct, T <: ClassSP}
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

mutable struct Score
    score::Float64
    parameters::NamedTuple
    apply_score
    apply_threshold
end

struct CompoundID{C <: Union{ClassSP, Nothing}}
    sum::NTuple{3, Int}
    lcb::Union{Lcb, Nothing}
    acyl::Union{NACYL, Nothing}
end

mutable struct CompoundSP
    class::ClassSP
    sum::NTuple{3, Int}   # (#C, #db, #OH)
    chain::Union{Chain, Nothing}
    fragments::Table # ion1, ion2, source, id
    area::Tuple{Float64, Float64}
    states::Vector{Int}
    results::Vector
    project
    CompoundSP(class, sum, chain, fragments, area, states, results, project) = new(class, sum, chain, fragments, area, states, results, project)
    CompoundSP(class, sum, chain) = new(class, sum, chain)
end

mutable struct AnalyteSP
    compounds::Vector{CompoundSP}
    rt::Float64
    states::Vector{Int}
    scores::Tuple{Float64, Float64}
end
# ==()

struct IonPlus{T}
    ions::T
end

struct IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}}
    ions::Tuple{S, T}
end

IonComparison(ion1, ion2) = IonComparison((ion1, ion2))

abstract type AbstractRule end

struct Rule{T} <: AbstractRule
    criteria::T
    rule
end

struct EmptyRule <: AbstractRule end

@as_record Rule

struct RuleSet
    mode::Symbol
    rule::Union{AbstractRule, ClassSP, Chain}
end

struct RuleUnion <: AbstractRule
    rules
    exception
    RuleUnion(x...; exception = EmptyRule()) = new(x, exception)
end

struct RuleMode <: AbstractRule
    rules
    exception
    RuleMode(x...; exception = EmptyRule()) = new(x, exception)
end

abstract type AbstractResult{T} end

struct Result{T} <: AbstractResult{T}
    matched::Bool
    rule::T
end
@as_record Result

struct PartialResult{T, C <: Union{Chain, ClassSP}} <: AbstractResult{T}
    matched::Bool
    rule::T
    result::C
end

@as_record PartialResult

struct Hydroxyl 
    spb::LCB
end

struct AcylIon{T <: LCB}
    ions
end

AcylIon{T}(ions...) where {T <: LCB} = AcylIon{T}(ions)

for fn in (:IonPlus, )
    @eval begin
        $fn(x...) = $fn(x)
    end
end

@as_record IonPlus

### connect db
abstract type Data end

# precursor/product, ion => name, ms/mz => m/z, mw => molecular weight
struct PreIS <: Data
    raw::Table #id, mz1, scan(index of mz2), area, ev, rt 
    range::Vector{Tuple{Float64, Float64}} 
    mz2::Vector{Float64}
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end

struct MRM <: Data
    raw::Table #id, mz1, mz2, area, ev, rt
    mz2::Vector{Float64}
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end


struct Project
    analytes::Vector{AnalyteSP}
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
    arg
end
not(arg) = Inv(arg)

include("io.jl")
include("interface.jl")
include("utils.jl")
include("mw.jl")
include("library.jl")
include("rule.jl")
include("data.jl")
include("preis.jl")
include("ruleID.jl")
include("score.jl")
include("query.jl")
include("mrm.jl")

# S: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-237.2224
# DHS: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-239.2224
# Phyto(4-OH): Acyl + C3H5NO-422.4004, SPB;3O - CH8NO-267.2329, Acyl + C2H5NO-410.4004, SPB;3O - C2H8NO-255.233
# α-hydroxyl: Acyl - CH2O-337.3476, Acyl + O-383.3531
# β-hydroxyl: SPB;2O + C2H2O - 340.246


end
