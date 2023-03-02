const ISOMER = Dict{Type, Tuple}()

"""
    ClassSP

Abstract type for all classes of sphingolipids.
"""
abstract type ClassSP end

"""
    SPB <: ClassSP

The type repressenting sphingoid bases.
"""
struct SPB <: ClassSP end
"""
    Cer <: ClassSP

The type repressenting ceramides.
"""
struct Cer <: ClassSP end
"""
    CerP <: ClassSP

The type repressenting ceramide-1-phosphates.
"""
struct CerP <: ClassSP end
"""
    SM <: ClassSP

The type repressenting sphingomyelin.
"""
struct SM <: ClassSP end
"""
    HexCer <: ClassSP

The type repressenting monohexosylceramides.
"""
struct HexCer <: ClassSP end
"""
    SHexCer <: ClassSP

The type repressenting monohexosylsulfatides.
"""
struct SHexCer <: ClassSP end
"""
    SHexHexCer <: ClassSP

The type repressenting dihexosylsulfatides.
"""
struct SHexHexCer <: ClassSP end
"""
    Hex2Cer <: ClassSP

The type repressenting dihexosylceramides.
"""
struct Hex2Cer <: ClassSP end
"""
    Hex3Cer <: ClassSP

The type repressenting trihexosylceramides, including Gb3, iGb3.
"""
struct Hex3Cer <: ClassSP end

"""
    HexNAcHex2Cer <: ClassSP

Abstract type for all iomers with 1 HexNAc and 2 Hex.
"""
abstract type HexNAcHex2Cer <: ClassSP end
"""
    HexNAcHex2Cer_ <: HexNAcHex2Cer

The type repressenting all iomers with 1 HexNAc and 2 Hex.
"""
struct HexNAcHex2Cer_ <: HexNAcHex2Cer
    isomer
end
"""
    HexNAc_Hex2Cer <: HexNAcHex2Cer

The type repressenting HexNAc-Hex-Hex-Cer, including GA2.
"""
struct HexNAc_Hex2Cer <: HexNAcHex2Cer end
push!(ISOMER, HexNAcHex2Cer_ => (HexNAc_Hex2Cer(),))
"""
    HexNAcHex3Cer <: ClassSP

Abstract type for all iomers with 1 HexNAc and 3 Hex.
"""
abstract type HexNAcHex3Cer <: ClassSP end
"""
    HexNAcHex3Cer_ <: HexNAcHex3Cer

The type repressenting all iomers with 1 HexNAc and 3 Hex.
"""
struct HexNAcHex3Cer_ <: HexNAcHex3Cer
    isomer
end
"""
    HexNAc_Hex3Cer <: HexNAcHex3Cer

The type repressenting HexNAc-Hex-Hex-Hex-Cer, including Gb4, iGb4.
"""
struct HexNAc_Hex3Cer <: HexNAcHex3Cer end
"""
    Hex_HexNAc_Hex2Cer <: HexNAcHex3Cer

The type repressenting Hex-HexNAc-Hex-Hex-Cer, including GA1.
"""
struct Hex_HexNAc_Hex2Cer <: HexNAcHex3Cer end
push!(ISOMER, HexNAcHex3Cer_ => (HexNAc_Hex3Cer(), Hex_HexNAc_Hex2Cer()))
"""
    GM4 <: ClassSP

The type repressenting GM4 (NeuAc-Gal-Cer).
"""
struct GM4 <: ClassSP end
"""
    GM3 <: ClassSP

The type repressenting GM3 (NeuAc-Gal-Glc-Cer).
"""
struct GM3 <: ClassSP end
"""
    GM2 <: ClassSP

The type repressenting GM2 (GalNAc-(NeuAc-)Gal-Glc-Cer).
"""
struct GM2 <: ClassSP end
"""
    GM1 <: ClassSP

Abstract type for all GM1 ((NeuAc-)Gal-GalNAc-Gal-Glc-Cer).
"""
abstract type GM1 <: ClassSP end
"""
    GM1_ <: GM1

The type repressenting iosmers of GM1 ((NeuAc-)Gal-GalNAc-Gal-Glc-Cer).
"""
struct GM1_ <: GM1
    isomer
end
"""
    GM1a <: GM1

The type repressenting GM1a (Gal-GalNAc-(NeuAc-)Gal-Glc-Cer).
"""
struct GM1a <: GM1 end
"""
    GM1b <: GM1

The type repressenting GM1b (NeuAc-Gal-GalNAc-Gal-Glc-Cer).
"""
struct GM1b <: GM1 end
push!(ISOMER, GM1_ => (GM1a(), GM1b()))
"""
    GD3 <: ClassSP

The type repressenting GD3 (NeuAc-NeuAc-Gal-Glc-Cer).
"""
struct GD3 <: ClassSP end
"""
    GD2 <: ClassSP

The type repressenting GD2 (GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GD2 <: ClassSP end
"""
    GD1 <: ClassSP

Abstract type for all GD1 ((NeuAc2-)Gal-GalNAc-Gal-Glc-Cer).
"""
abstract type GD1 <: ClassSP end
"""
    GD1_ <: GD1

The type repressenting isomers of GD1 ((NeuAc2-)Gal-GalNAc-Gal-Glc-Cer).
"""
struct GD1_ <: GD1
    isomer
end
"""
    GD1a <: GD1

The type repressenting GD1a (NeuAc-Gal-GalNAc-(NeuAc-)Gal-Glc-Cer).
"""
struct GD1a <: GD1 end
"""
    GD1b <: GD1

The type repressenting GD1b (Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GD1b <: GD1 end
"""
    GD1c <: GD1

The type repressenting GD1c (NeuAc-NeuAc-Gal-GalNAc-Gal-Glc-Cer).
"""
struct GD1c <: GD1 end
"""
    GD1α <: GD1

The type repressenting GD1α (NeuAc-Gal-(NeuAc-)GalNAc-Gal-Glc-Cer).
"""
struct GD1α <: GD1 end
push!(ISOMER, GD1_ => (GD1a(), GD1b(), GD1c(), GD1α()))
"""
    GT3 <: ClassSP

The type repressenting GT3 (NeuAc-NeuAc-NeuAc-Gal-Glc-Cer).
"""
struct GT3 <: ClassSP end
"""
    GT2 <: ClassSP

The type repressenting GT2 (GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GT2 <: ClassSP end
"""
    GT1 <: ClassSP

Abstract type for all GT1 ((NeuAc3-)Gal-GalNAc-Gal-Glc-Cer).
"""
abstract type GT1 <: ClassSP end
"""
    GT1_ <: GT1

The type repressenting isomers of GT1 ((NeuAc3-)Gal-GalNAc-Gal-Glc-Cer).
"""
struct GT1_ <: GT1
    isomer
end
"""
    GT1a <: GT1

The type repressenting GT1a (NeuAc-NeuAc-Gal-GalNAc-(NeuAc-)Gal-Glc-Cer).
"""
struct GT1a <: GT1 end
"""
    GT1aα <: GT1
    
The type repressenting GT1aα (NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-)Gal-Glc-Cer).
"""
struct GT1aα <: GT1 end
"""
    GT1b <: GT1

The type repressenting GT1b (NeuAc-Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GT1b <: GT1 end
"""
    GT1c <: GT1

The type repressenting GT1c (NeuAc-NeuAc-NeuAc-Gal-GalNAc-Gal-Glc-Cer).
"""
struct GT1c <: GT1 end
push!(ISOMER, GT1_ => (GT1a(), GT1b(), GT1c(), GT1aα()))
"""
    GQ1 <: ClassSP

Abstract type for all GQ1 ((NeuAc4-)Gal-GalNAc-Gal-Glc-Cer).
"""
abstract type GQ1 <: ClassSP end
"""
    GQ1_ <: GQ1

The type repressenting isomers of GQ1 ((NeuAc4-)Gal-GalNAc-Gal-Glc-Cer).
"""
struct GQ1_ <: GQ1
    isomer
end
"""
    GQ1b <: GQ1

The type repressenting GQ1b (NeuAc-NeuAc-Gal-GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GQ1b <: GQ1 end
"""
    GQ1bα <: GQ1

The type repressenting GQ1bα (NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GQ1bα <: GQ1 end
"""
    GQ1c <: GQ1

The type repressenting GQ1c (NeuAc-Gal-GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GQ1c <: GQ1 end
push!(ISOMER, GQ1_ => (GQ1b(), GQ1c(), GQ1bα()))
"""
    GP1 <: ClassSP

Abstract type for all GP1 ((NeuAc5-)Gal-GalNAc-Gal-Glc-Cer).
"""
abstract type GP1 <: ClassSP end
"""
    GP1_ <: GP1

The type repressenting isomers of GP1 ((NeuAc5-)Gal-GalNAc-Gal-Glc-Cer).
"""
struct GP1_ <: GP1
    isomer
end
"""
    GP1c <: GP1

The type repressenting GP1c (NeuAc-NeuAc-Gal-GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer).
"""
struct GP1c <: GP1 end
"""
    GP1cα <: GP1

The type repressenting GP1cα (NeuAc-Gal-(NeuAc-)GalNAc-(NeuAc-NeuAc-NeuAc-)Gal-Glc-Cer).
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

A `NamedTuple` contains specifc groups of sphingolipids.

* `nana`: a `Dict` where the keys are `Int` and the values are types with corresponding number of `NeuAc`.
* `series`: a `NamedTuple` specfying series of gangliosides
    * `as`: asialo-seris 
    * `a`: a-series
    * `b`: b-series
    * `c`: c-series
* `fg`: a `NamedTuple` specfying types with specific functional group
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