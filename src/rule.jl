
rule(c::ClassGSL, x...) = rule(c)
rule(::Cer) = ISF(ProtonationNLH2O(), Cer()) => (Cer(), ())
rule(::HexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(ProtonationNLH2O(), Cer())) => (HexCer(), ())
rule(::SHexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(ProtonationNLH2O(), Cer())) => (SHexCer(), ())
rule(::Hex2Cer) =  ISF(ProtonationNLH2O(), Cer()) => (Hex2Cer(), ())
rule(::Hex3Cer) = (ISF(Protonation(), Hex2Cer()), ISF(ProtonationNLH2O(), Cer())) => (Hex3Cer(), Hex3Cer())
rule(::HexNAcHex2Cer) = HexNAcHex2Cer()
rule(::HexNAcHex3Cer) = HexNAcHex => (Hex_HexNAc_Hex2Cer(), HexNAc_Hex3Cer())

rule(::GM4) = NANA => (GM4(), ())
rule(::GM3) = NANA => (GM3(), ())
rule(::GM2) = NANA => (GM2(), ())
rule(::GM1) = NANA => (
    IonUnion(
        (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) => (GM1a(), ()),
        IonComparison(HexNAcHexNANA, HexNAcHex) => (
            0           => GM1a(),
            (0, 5)      => GM1(GM1a(), GM1b()),
            (5, Inf64)  => GM1b(),
            NaN         => ()
        )
    ),
    ()
)

rule_ab(::GM1) = NANA => (
    HexNAcHex => (
        GM1a(),
        ()
    ),
    ()
)

rule(::GD3) = NANA => (
    NANA2 => (GD3(), ()),
    ()
)

rule(::GD2) = NANA => (
    NANA2 => (GD2(), ()),
    ()
)

rule(::ClassGD1, mode::Symbol) = 
    mode == :ab ? rule_ab(GD1()) : rule(GD1())

rule_ab(::ClassGD1) = NANA => (
    IonUnion(
        NANA2           => (GD1b(), ()),
        HexNAcHexNANA   => (GD1a(), ())
    ),
    ()
)    
#=
rule(::GD1) = NANA => (
    HexNAcHexNANA2 => (
        IonUnion(
            IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
                (0, 5)      => (GD1c(), GD1a()),
                (5, Inf64)  => (GD1c(),),
                NaN         => ()
            ),
            IonComparison(IonPlus(HexNAcHexNANA2, HexNAcHexNANA), HexNAcHex) => (
                (0, 5)      => (GD1b()),
                (5, Inf64)  => (),
                NaN         => ()
            )
        ),
        rule_ab(GD1())
    ),
    ()
)
=#
rule(::ClassGD1) = NANA => (
    NANA2 => (
        IonUnion(
            HexNAcHexNANA2 => (
                GD1c(),
                ()
            ),
            IonComparison(HexNAcHexNANA2, HexNAcHex) => (
                (0, 5) => GD1b(),
                (5, Inf) => ()
            )
        ),
        IonUnion(
            HexNAcHexNANA2 => (
                GD1α(),
                ()
            ),
            (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) => (
                GD1a(),
                ()
            )
        )
        
    ),
    ()
)


rule(::GT3) = NANA => (
    NANA23 => (GT3(), ()),
    ()
)

rule(::GT2) = NANA => (
    NANA23 => (GT2(), ()),
    ()
)

rule(::ClassGT1, mode::Symbol) = 
    mode == :ab ? rule_ab(GT1()) : rule(GT1())

rule_α(::GT1a) = IonMode(
    Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) => (
        GT1aα(),
        GT1a()
    ),
    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) => (
        GT1a(),
        GT1aα()
    )
)

rule_ab(::GT1) = NANA => (
    NANA2 => (
        HexNAcHexNANA2 => (
            IonUnion(
                rule_α(GT1a()),
                IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
                    (0, 5) => GT1b(),
                    (5, Inf) => ()
                )
            ), 
            GT1b()
        ),
        HexNAcHexNANA2 => (
            GT1aα(),
            HexNAcHexNANA => (
                GT1b(),
                ()
            )
        )
    ),
    ()
)

rule_ab_neg(::GT1) = IonComparison(NANA1, NANA2) => (
    (0, 0.4)    => GT1a(),
    (0.4, 7)    => GT1(GT1a(), GT1b()),
    (7, Inf64)  => GT1b()
)

#=
rule_ab(::GT1) = NANA => (
    IonUnion(
        (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) => (
            rule_α(GT1()),
            ()
        ),
        IonMode(
            IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
                0           => (GT1b(),),
                (0, 5)      => IonUnion(rule_α(GT1()), GT1b()),
                (5, Inf64)  => (rule_α(GT1()),),
                NaN         => ()
            ),
            IonComparison(NANA1, NANA2) => (
                (0, 0.4)    => (rule_α(GT1()),),
                (0.4, 7)    => IonUnion(rule_α(GT1()), GT1b()),
                (7, Inf64)  => (GT1b(),)
            )
        )
    )
    ,
    ()
)

rule_α(::GT1) = Ion(Protonation(), Glycan(HexNAc(), NeuAc())) => (
    GT1aα(),
    ()
)

rule(::GT1) = NANA => (
    IonComparison(IonPlus(HexNAcHexNANA2, HexNAcHexNANA), HexNAcHex) => (
        0           => GT1c(),
        (0, 5)      => (
            IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
                0           => (GT1b(), GT1c()),
                (0, 5)      => IonUnion(rule_α(GT1()), GT1b(), GT1c()),
                (5, Inf64)  => IonUnion(rule_α(GT1()), GT1c())
            )
        ),
        (5, Inf64)  => rule_ab(GT1()),
        NaN         => ()
    ),
    ()
)
=#
rule(::GT1) = NANA => (
    NANA2 => (
        HexNAcHexNANA2 => (
            IonUnion(
                rule_α(GT1a()),
                IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
                    (0, 5) => GT1b(),
                    (5, Inf) => ()
                ),
                IonComparison(IonPlus(HexNAcHexNANA2, HexNAcHexNANA), HexNAcHex) => (
                    (0, 5) => GT1c(),
                    (5, Inf) => ()
                )
            ), 
            HexNAcHexNANA => (
                IonComparison(HexNAcHexNANA, HexNAcHex) => (
                    (0, 5) => GT1(GT1c(), GT1b()),
                    (5, Inf) => GT1b(),
                ), 
                HexNAcHex => (
                    GT1c(),
                    rule_ab_neg(GT1())
                )
            )
        ),
        HexNAcHexNANA2 => (
            GT1aα(),
            HexNAcHexNANA => (
                GT1b(),
                ()
            )
        )
    )
)

rule(::GQ1, mode::Symbol) = 
    mode == :ab ? rule_ab(GQ1()) : rule(GQ1())

rule_ab(::GQ1) = NANA => (
    HexNAcHexNANA2 => (
        IonMode(
            Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) => (
                GQ1bα(),
                GQ1b()
            ),
            Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) => (
                GQ1b(),
                GQ1bα()
            )
        )
    ),
    ()
)

rule(::GQ1) = NANA => (
    IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
        0           => GQ1c(),
        (0, 5)      => IonUnion(GQ1c(), rule_ab(GQ1())),
        (5, Inf64)  => rule_ab(GQ1()),
        NaN         => ()
    ),
    ()
)

rule_ab(::GP1) = ()

rule(::GP1) = NANA => (
    HexNAcHexNANA2 => (
        IonMode(
            Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) => (
                GP1cα(),
                GP1c()
            ),
            Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) => (
                GP1c(),
                GP1cα()
            )
        ), 
        ()
    ),
    ()
)

rule_acyl(spb::T, ::Acyl{N}) where {T <: LCB, N} = AcylIon{T}(Ion(LossCH2O(), Acyl{N}()), Ion(AddO(), Acyl{N}())) => (
    Chain(spb, Acylα{N}()),
    Ion(AddC2H2O(), spb) => (
        Chain(spb, Acylβ{N}()),
        Chain(spb, Acyl{N}())
    )
)

to_phyto(::SPB3{C, N}) where {C, N} = PhytoSPB3{C, N}()
to_phyto(::SPB4{C, N}) where {C, N} = PhytoSPB4{C, N}()

rule_phytoacyl(spb::T, ::Acyl{0}) where {T <: LCB} = AcylIon{T}(Ion(LossCH8NO(), spb), Ion(AddC3H5NO(), Acyl{0}()), Ion(LossC2H8NO(), spb), Ion(AddC2H5NO(), Acyl{0}())) => (
    Chain(to_phyto(spb), Acyl{0}()),
    Chain(spb, Acyl{0}())
)

rule_phytoacyl(spb::T, acyl::Acyl) where {T <: LCB} =  AcylIon{T}((Ion(AddC3H5NO(), acyl),)) => (
    rule_acyl(to_phyto(spb), acyl),
    rule_acyl(spb, acyl)
)

rule(::Chain{<: SPB2{C}, <: ACYL}) where C = IonComparison(Ion(ProtonationNL2H2O(), SPB2{C, 2}()), Ion(ProtonationNLH2O(), SPB2{C, 2}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB2{C, 2}()) => (
            0 => Chain(SPB2{C, 2}(), Acyl{0}()),
            1 => rule_acyl(SPB2{C, 2}(), Acyl{1}()),
            2 => rule_acyl(SPB2{C, 2}(), Acyl{2}())
        )

    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB2{C, 1}()) => (
            0 => Chain(SPB2{C, 1}(), Acyl{0}()),
            1 => rule_acyl(SPB2{C, 1}(), Acyl{1}()),
            2 => rule_acyl(SPB2{C, 1}(), Acyl{2}())
        )
    )
)

rule(::Chain{<: SPB3{C}, <: ACYL}) where C = IonComparison(Ion(ProtonationNL3H2O(), SPB3{C, 3}()), Ion(ProtonationNL2H2O(), SPB3{C, 3}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB3{C, 3}()) => (
            0 => rule_phytoacyl(SPB3{C, 3}(), Acyl{0}()),
            1 => rule_phytoacyl(SPB3{C, 3}(), Acyl{1}()),
            2 => rule_phytoacyl(SPB3{C, 3}(), Acyl{2}())
        )
    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB3{C, 2}()) => (
            0 => Chain(SPB3{C, 2}(), Acyl{0}()),
            1 => rule_acyl(SPB3{C, 2}(), Acyl{1}()),
            2 => rule_acyl(SPB3{C, 2}(), Acyl{2}())
        )
    )
)

rule(::Chain{<: SPB4{C}, <: ACYL}) where C = IonComparison(Ion(ProtonationNL3H2O(), SPB4{C, 3}()), Ion(ProtonationNL2H2O(), SPB4{C, 3}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB4{C, 3}()) => (
            0 => rule_phytoacyl(SPB4{C, 3}(), Acyl{0}()),
            1 => rule_phytoacyl(SPB4{C, 3}(), Acyl{1}()),
            2 => rule_phytoacyl(SPB4{C, 3}(), Acyl{2}())
        )
    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB4{C, 2}()) => (
            0 => Chain(SPB4{C, 2}(), Acyl{0}()),
            1 => rule_acyl(SPB4{C, 2}(), Acyl{1}()),
            2 => rule_acyl(SPB4{C, 2}(), Acyl{2}())
        )
    )
)

const CONNECTION = Dict{ClassGSL, ClassGSL}(
    GP1()                   => GQ1(),
    GQ1()                   => GT1(),
    GT1()                   => GD1(),
    GD1()                   => GM1(),
    GM1()                   => GM3(),
    GM3()                   => Hex2Cer(),
    Hex2Cer()               => HexCer(),
    Cer()                   => Cer(),
    GD3()                   => GM3(),
    GT3()                   => GD3(),
    GM2()                   => GM3(),
    GD2()                   => GD3(),
    GT2()                   => GT3(),
    SHexCer()               => HexCer(),
    GM4()                   => HexCer(),
    Hex3Cer()               => Hex2Cer(),
    HexNAcHex2Cer()         => Hex2Cer(),
    HexNAcHex3Cer()         => Hex2Cer(),
    Hex_HexNAc_Hex2Cer()    => HexNAcHex2Cer(),
    HexNAc_Hex3Cer()        => Hex3Cer(),
    HexCer()                => Cer()
)

find_connected(cpd::CompoundGSL, analyte::AnalyteGSL) = findall(id -> connected(cpd, id), analyte.identification)
connected(cpd::CompoundGSL, id::CompoundGSL) = ischaincompatible(cpd, id) && begin
    class1 = haskey(CONNECTION, cpd.class) ? cpd.class : deisomerized(cpd.class)
    class2 = haskey(CONNECTION, id.class) ? id.class : deisomerized(id.class)
    connected(class1, class2) || connected(class1, class2)
end

connected(cls1::ClassGSL, cls2::ClassGSL) = cls2 == cls1 || connected(CONNECTION[cls1], cls2)
connected(cls1::ClassGSL, cls2::Cer) = true
connected(cls1::Cer, cls2::ClassGSL) = false
connected(cls1::Cer, cls2::Cer) = true

isf(cls::ClassGSL) = push!(isf(CONNECTION[cls]), cls)
isf(cls::Cer) = ClassGSL[cls]
