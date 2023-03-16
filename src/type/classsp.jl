const ISOMER = Dict{Type, Tuple}()

"""
    ClassSP

Abstract type for all classes of sphingolipids.
"""
abstract type ClassSP end

"""
    SPB <: ClassSP

Sphingoid bases.
"""
struct SPB <: ClassSP end
"""
    Cer <: ClassSP

Ceramides.
"""
struct Cer <: ClassSP end
"""
    CerP <: ClassSP

Ceramide-1-phosphates.
"""
struct CerP <: ClassSP end
"""
    SM <: ClassSP

Sphingomyelins.
"""
struct SM <: ClassSP end
"""
    HexCer <: ClassSP

Monohexosylceramides.
"""
struct HexCer <: ClassSP end
"""
    SHexCer <: ClassSP

Monohexosylsulfatides.
"""
struct SHexCer <: ClassSP end
"""
    SHexHexCer <: ClassSP

Dihexosylsulfatides.
"""
struct SHexHexCer <: ClassSP end
"""
    Hex2Cer <: ClassSP

Dihexosylceramides.
"""
struct Hex2Cer <: ClassSP end
"""
    Hex3Cer <: ClassSP

Trihexosylceramides, including Gb3 and iGb3.
"""
struct Hex3Cer <: ClassSP end

"""
    HexNAcHex2Cer <: ClassSP

Abstract type for all isomers with 1 HexNAc and 2 Hex.
"""
abstract type HexNAcHex2Cer <: ClassSP end
"""
    HexNAcHex2Cer_ <: HexNAcHex2Cer

A type representing all isomers with 1 HexNAc and 2 Hex.
"""
struct HexNAcHex2Cer_ <: HexNAcHex2Cer
    isomer
end
"""
    HexNAc_Hex2Cer <: HexNAcHex2Cer

HexNAc-Hex-Hex-Cer, including GA2.
"""
struct HexNAc_Hex2Cer <: HexNAcHex2Cer end
push!(ISOMER, HexNAcHex2Cer_ => (HexNAc_Hex2Cer(),))
"""
    HexNAcHex3Cer <: ClassSP

Abstract type for all isomers with 1 HexNAc and 3 Hex.
"""
abstract type HexNAcHex3Cer <: ClassSP end
"""
    HexNAcHex3Cer_ <: HexNAcHex3Cer

A type representing all isomers with 1 HexNAc and 3 Hex.
"""
struct HexNAcHex3Cer_ <: HexNAcHex3Cer
    isomer
end
"""
    HexNAc_Hex3Cer <: HexNAcHex3Cer

HexNAc-Hex-Hex-Hex-Cer, including Gb4 and iGb4.
"""
struct HexNAc_Hex3Cer <: HexNAcHex3Cer end
"""
    Hex_HexNAc_Hex2Cer <: HexNAcHex3Cer

Hex-HexNAc-Hex-Hex-Cer, including GA1.
"""
struct Hex_HexNAc_Hex2Cer <: HexNAcHex3Cer end
push!(ISOMER, HexNAcHex3Cer_ => (HexNAc_Hex3Cer(), Hex_HexNAc_Hex2Cer()))
"""
    GM4 <: ClassSP

GM4, i.e., NeuAc-Gal-Cer.
"""
struct GM4 <: ClassSP end
"""
    GM3 <: ClassSP

GM3, i.e., NeuAc-Gal-Glc-Cer.
"""
struct GM3 <: ClassSP end
"""
    GM2 <: ClassSP

GM2, i.e., GalNAc-(NeuAc-)Gal-Glc-Cer.
"""
struct GM2 <: ClassSP end
"""
    GM1 <: ClassSP

Abstract type for all GM1, i.e., [NeuAc]Gal-GalNAc-Gal-Glc-Cer.
"""
abstract type GM1 <: ClassSP end
"""
    GM1_ <: GM1

A type representing multiple isomers of GM1, i.e., [NeuAc]Gal-GalNAc-Gal-Glc-Cer.
"""
struct GM1_ <: GM1
    isomer
end
"""
    GM1a <: GM1

GM1a, i.e., Gal-GalNAc-(NeuAc-)Gal-Glc-Cer.
"""
struct GM1a <: GM1 end
"""
    GM1b <: GM1

GM1b, i.e., NeuAc-Gal-GalNAc-Gal-Glc-Cer.
"""
struct GM1b <: GM1 end
push!(ISOMER, GM1_ => (GM1a(), GM1b()))
"""
    GD3 <: ClassSP

GD3, i.e., NeuAc-NeuAc-Gal-Glc-Cer.
"""
struct GD3 <: ClassSP end
"""
    GD2 <: ClassSP

GD2, i.e., GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GD2 <: ClassSP end
"""
    GD1 <: ClassSP

Abstract type for all GD1, i.e., [NeuAc₂]Gal-GalNAc-Gal-Glc-Cer.
"""
abstract type GD1 <: ClassSP end
"""
    GD1_ <: GD1

A type representing multiple isomers of GD1, i.e., [NeuAc₂]Gal-GalNAc-Gal-Glc-Cer.
"""
struct GD1_ <: GD1
    isomer
end
"""
    GD1a <: GD1

GD1a, i.e., NeuAc-Gal-GalNAc-(NeuAc-)Gal-Glc-Cer.
"""
struct GD1a <: GD1 end
"""
    GD1b <: GD1

GD1b, i.e., Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GD1b <: GD1 end
"""
    GD1c <: GD1

GD1c, i.e., NeuAc-NeuAc-Gal-GalNAc-Gal-Glc-Cer.
"""
struct GD1c <: GD1 end
"""
    GD1α <: GD1

GD1α, i.e., NeuAc-Gal-(NeuAc-)GalNAc-Gal-Glc-Cer.
"""
struct GD1α <: GD1 end
push!(ISOMER, GD1_ => (GD1a(), GD1b(), GD1c(), GD1α()))
"""
    GT3 <: ClassSP

GT3, i.e., NeuAc-NeuAc-NeuAc-Gal-Glc-Cer.
"""
struct GT3 <: ClassSP end
"""
    GT2 <: ClassSP

GT2, i.e., GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GT2 <: ClassSP end
"""
    GT1 <: ClassSP

Abstract type for all GT1, i.e., [NeuAc₃]Gal-GalNAc-Gal-Glc-Cer.
"""
abstract type GT1 <: ClassSP end
"""
    GT1_ <: GT1

A type representing multiple isomers of GT1, i.e., [NeuAc₃]Gal-GalNAc-Gal-Glc-Cer.
"""
struct GT1_ <: GT1
    isomer
end
"""
    GT1a <: GT1

GT1a, i.e., NeuAc-NeuAc-Gal-GalNAc-(NeuAc-)Gal-Glc-Cer.
"""
struct GT1a <: GT1 end
"""
    GT1aα <: GT1
    
GT1aα, i.e., NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-)Gal-Glc-Cer.
"""
struct GT1aα <: GT1 end
"""
    GT1b <: GT1

GT1b, i.e., NeuAc-Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GT1b <: GT1 end
"""
    GT1c <: GT1

GT1c, i.e., NeuAc-NeuAc-NeuAc-Gal-GalNAc-Gal-Glc-Cer.
"""
struct GT1c <: GT1 end
push!(ISOMER, GT1_ => (GT1a(), GT1b(), GT1c(), GT1aα()))
"""
    GQ1 <: ClassSP

Abstract type for all GQ1, i.e., [NeuAc₄]Gal-GalNAc-Gal-Glc-Cer.
"""
abstract type GQ1 <: ClassSP end
"""
    GQ1_ <: GQ1

multiple isomers of GQ1, i.e., [NeuAc₄]Gal-GalNAc-Gal-Glc-Cer.
"""
struct GQ1_ <: GQ1
    isomer
end
"""
    GQ1b <: GQ1

GQ1b, i.e., NeuAc-NeuAc-Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GQ1b <: GQ1 end
"""
    GQ1bα <: GQ1

GQ1bα, i.e., NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GQ1bα <: GQ1 end
"""
    GQ1c <: GQ1

GQ1c, i.e., NeuAc-Gal-GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GQ1c <: GQ1 end
push!(ISOMER, GQ1_ => (GQ1b(), GQ1c(), GQ1bα()))
"""
    GP1 <: ClassSP

Abstract type for all GP1, i.e., [NeuAc₅]Gal-GalNAc-Gal-Glc-Cer.
"""
abstract type GP1 <: ClassSP end
"""
    GP1_ <: GP1

A type representing multiple isomers of GP1, i.e., [NeuAc₅]Gal-GalNAc-Gal-Glc-Cer.
"""
struct GP1_ <: GP1
    isomer
end
"""
    GP1c <: GP1

GP1c, i.e., NeuAc-NeuAc-Gal-GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GP1c <: GP1 end
"""
    GP1cα <: GP1

GP1cα, i.e., NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer.
"""
struct GP1cα <: GP1 end
push!(ISOMER, GP1_ => (GP1c(), GP1cα()))

SPDB[:ISOMER] = ISOMER

for (class, super) in zip((:HexNAcHex2Cer_, :HexNAcHex3Cer_, :GM1_, :GD1_, :GT1_, :GQ1_, :GP1_), (:HexNAcHex2Cer, :HexNAcHex3Cer, :GM1, :GD1, :GT1, :GQ1, :GP1))
    @eval begin
        $class(cls::Vararg{<: $super}) = $class(cls)
        $class() = $class(SPDB[:ISOMER][$class])
        $super(cls::Vararg{<: $super}) = $class(cls)
        $super() = $class(SPDB[:ISOMER][$class])
        Isomer(cls::Vararg{<: Type{T}}) where {T <: $super} = $class((c() for c in cls))
    end
end

"""
    const CLS

A `NamedTuple` contains specific groups of sphingolipids.

* `nana`: a `Dict` where the keys are `Int` and the values are types with corresponding number of `NeuAc`.
* `series`: a `NamedTuple` specifying the series of gangliosides
    * `as`: asialo-seris 
    * `a`: a-series
    * `b`: b-series
    * `c`: c-series
* `fg`: a `NamedTuple` specifying types with specific functional group
    * `nana`: types with `NeuAc`(s)
    * `sulfate`: types with sulfate(s)
"""
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