"""
    HeadGroup

Abstract type for head groups.
"""
abstract type HeadGroup end
"""
    Sugar <: HeadGroup

Abstract type for all monosaccharides.
"""
abstract type Sugar <: HeadGroup end
"""
    NeuAc <: Sugar

N-acetyl neuraminic acid(nana).
"""
struct NeuAc <: Sugar end
"""
    NeuAc <: Sugar

hexose.
"""
struct Hex <: Sugar end
"""
    HexNAc <: Sugar

N-acetyl hexosamine.
"""
struct HexNAc <: Sugar end
"""
    Glycan{T} <: HeadGroup

Polysaccharides.
"""
struct Glycan{T} <: HeadGroup
    chain::T
end
Glycan(sugar::Sugar...) = Glycan(sugar)
"""
    PhosphoCholine <: HeadGroup

Phosphocholine.
"""
struct PhosphoCholine <: HeadGroup end
"""
    Adduct

Abstract type for all kinds of adducts.
"""
abstract type Adduct end
"""
    Pos <: Adduct

Abstract type for adducts with positive charges.
"""
abstract type Pos <: Adduct end
"""
    Neg <: Adduct

Abstract type for adducts with negative charges.
"""
abstract type Neg <: Adduct end
"""
    AbstractIon{S, T}

Abstract type for Ions.
"""
abstract type AbstractIon{S, T} end
"""
    Ion{S <: Adduct, T <: Union{HeadGroup, ClassSP, ChainSP}} <: AbstractIon{S, T}

Ion forming in the collisional cell.
"""
struct Ion{S <: Adduct, T <: Union{HeadGroup, ClassSP, ChainSP}} <: AbstractIon{S, T}
    adduct::S
    molecule::T
end
"""
    ISF{S <: Adduct, T <: ClassSP} <: AbstractIon{S, T}

Ion forming in the ion source(in-source fragmentation).
"""
struct ISF{S <: Adduct, T <: ClassSP} <: AbstractIon{S, T}
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
struct LossCH3 <: Neg end

const NANA = (Ion(ProtonationNL2H2O(), NeuAc()), Ion(DeprotonationNLH2O(), NeuAc()), Ion(ProtonationNLH2O(), NeuAc()))
const NANA1 = Ion(DeprotonationNLH2O(), NeuAc())
const NANA2 = Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc()))
const NANA3 = Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc(), NeuAc()))
const NANA23 = (Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc())), Ion(DeprotonationNLH2O(), Glycan(NeuAc(), NeuAc(), NeuAc())))
const HexNAcHex = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex()))
const HexNAcHexNANA = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex(), NeuAc()))
const HexNAcHexNANA2 = Ion(ProtonationNLH2O(), Glycan(HexNAc(), Hex(), NeuAc(), NeuAc()))

@as_record AbstractIon
@as_record Adduct