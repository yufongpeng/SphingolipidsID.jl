
rule(c::ClassGSL, x...) = rule(c)
rule(::Cer) = ISF(ProtonationNLH2O(), Cer()) => (Cer(), ())
rule(::HexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(ProtonationNLH2O(), Cer())) => (HexCer(), ())
rule(::SHexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) => (SHexCer(), ())
rule(::Hex2Cer) =  (ISF(ProtonationNLH2O(), Cer()), ISF(ProtonationNLH2O(), Hex2Cer())) => (
    Hex2Cer(), 
    Ion(ProtonationNLH2O(), Cer()) => (Hex2Cer(), ())
)
rule(::SHexHexCer) = (ISF(ProtonationNLH2O(), Hex2Cer()), ISF(Protonation(), Hex2Cer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) => (SHexHexCer(), ())
rule(::Hex3Cer) = (ISF(Protonation(), Hex2Cer()), ISF(ProtonationNLH2O(), Hex2Cer()), ISF(ProtonationNLH2O(), Cer())) => (
    Hex3Cer(),    
    Ion(ProtonationNLH2O(), Cer()) => (Hex3Cer(), ())
)
rule(::HexNAcHex2Cer) = Ion(ProtonationNLH2O(), HexNAc()) => (HexNAcHex2Cer(), ())
rule(::HexNAcHex3Cer) = Ion(ProtonationNLH2O(), HexNAc()) => (
    Ion(Protonation(), Hex3Cer()) => (Hex_HexNAc_Hex2Cer(), HexNAc_Hex3Cer()),
    ()
)

rule(::GM4) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) => (
    GM4(),    
    Ion(ProtonationNLH2O(), Cer()) => (GM4(), ())
)
rule(::GM3) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) => (
    GM3(),
    Ion(ProtonationNLH2O(), Cer()) => (GM3(), ())
)
rule(::GM2) = Ion(ProtonationNLH2O(), HexNAc()) => (GM2(), ())
rule(::GM1) = IonUnion(
    (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) => (GM1a(), ()),
    IonComparison(HexNAcHexNANA, HexNAcHex) => (
        0           => GM1a(),
        (0, 5)      => GM1(GM1a(), GM1b()),
        (5, Inf64)  => GM1b(),
        NaN         => ()
    )
)

rule_ab(::GM1) = HexNAcHex => (
    GM1a(),
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
        IonComparison(HexNAcHexNANA, HexNAcHex) => (
            0 => GD1b(),
            (0, 1) => GD1(GD1b(), GD1a()),
            (1, 10) => GD1(GD1a(), GD1b()),
            (10, Inf) => GD1a()
        ),
        NANA2           => (GD1b(), ())
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
        IonComparison(HexNAcHexNANA2, HexNAcHex) => (
            0 => GD1b(),
            (0, 1) => GD1(GD1b(), GD1c()),
            (1, 10) => GD1(GD1c(), GD1b()),
            (10, Inf) => GD1c()
        ),
        IonUnion(
            (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) => (
                GD1a(),
                ()
            ),
            HexNAcHexNANA2 => (
                GD1α(),
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
        IonComparison(HexNAcHexNANA2, HexNAcHexNANA) => (
            0 => GT1b(),
            (0, 1) => IonUnion(GT1b(), rule_α(GT1a())),
            (1, 10) => IonUnion(rule_α(GT1a()), GT1b()),
            (10, Inf) => rule_α(GT1a())
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
                    (0, 10) => GT1b(),
                    (10, Inf) => ()
                ),
                IonComparison(IonPlus(HexNAcHexNANA2, HexNAcHexNANA), HexNAcHex) => (
                    (0, 10) => GT1c(),
                    (10, Inf) => ()
                )
            ), 
            HexNAcHexNANA => (
                IonComparison(HexNAcHexNANA, HexNAcHex) => (
                    (0, 10) => GT1(GT1c(), GT1b()),
                    (10, Inf) => GT1b(),
                ), 
                HexNAcHex => (
                    GT1c(),
                    rule_ab_neg(GT1()) # Assume a/b only
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

to_phyto(::LCB3{N, C}) where {N, C} = PhytoSPB3{N, C}()
to_phyto(::LCB4{N, C}) where {N, C} = PhytoSPB4{N, C}()

rule_phytoacyl(spb::T, ::Acyl{0}) where {T <: LCB} = AcylIon{T}(Ion(AddC2H5NO(), Acyl{0}()), Ion(AddC3H5NO(), Acyl{0}()), Ion(LossCH8NO(), spb), Ion(LossC2H8NO(), spb)) => (
    Chain(to_phyto(spb), Acyl{0}()),
    Chain(spb, Acyl{0}())
)

rule_phytoacyl(spb::T, acyl::Acyl) where {T <: LCB} =  AcylIon{T}((Ion(AddC3H5NO(), acyl),)) => (
    rule_acyl(to_phyto(spb), acyl),
    rule_acyl(spb, acyl)
)

rule(::Chain{<: LCB2{N, C}, <: ACYL}) where {N, C} = IonComparison(Ion(ProtonationNL2H2O(), SPB2{2, C}()), Ion(ProtonationNLH2O(), SPB2{2, C}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB2{2, C}()) => (
            0 => Chain(SPB2{2, C}(), Acyl{0}()),
            1 => rule_acyl(SPB2{2, C}(), Acyl{1}()),
            2 => rule_acyl(SPB2{2, C}(), Acyl{2}())
        )

    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB2{1, C}()) => (
            0 => Chain(SPB2{1, C}(), Acyl{0}()),
            1 => rule_acyl(SPB2{1, C}(), Acyl{1}()),
            2 => rule_acyl(SPB2{1, C}(), Acyl{2}())
        )
    )
)

rule(::Chain{<: LCB3{N, C}, <: ACYL}) where {N, C} = IonComparison(Ion(ProtonationNL3H2O(), SPB3{3, C}()), Ion(ProtonationNL2H2O(), SPB3{3, C}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB3{3, C}()) => (
            0 => rule_phytoacyl(SPB3{3, C}(), Acyl{0}()),
            1 => rule_phytoacyl(SPB3{3, C}(), Acyl{1}()),
            2 => rule_phytoacyl(SPB3{3, C}(), Acyl{2}())
        )
    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB3{2, C}()) => (
            0 => Chain(SPB3{2, C}(), Acyl{0}()),
            1 => rule_acyl(SPB3{2, C}(), Acyl{1}()),
            2 => rule_acyl(SPB3{2, C}(), Acyl{2}())
        )
    )
)

rule(::Chain{<: LCB4{N, C}, <: ACYL}) where {N, C} = IonComparison(Ion(ProtonationNL3H2O(), SPB4{3, C}()), Ion(ProtonationNL2H2O(), SPB4{3, C}())) => (
    (0.0, 1.0) => (
        Hydroxyl(SPB4{3, C}()) => (
            0 => rule_phytoacyl(SPB4{3, C}(), Acyl{0}()),
            1 => rule_phytoacyl(SPB4{3, C}(), Acyl{1}()),
            2 => rule_phytoacyl(SPB4{3, C}(), Acyl{2}())
        )
    ),
    (1.0, Inf64) => (
        Hydroxyl(SPB4{2, C}()) => (
            0 => Chain(SPB4{2, C}(), Acyl{0}()),
            1 => rule_acyl(SPB4{2, C}(), Acyl{1}()),
            2 => rule_acyl(SPB4{2, C}(), Acyl{2}())
        )
    )
)

find_connected(cpd::CompoundGSL, analyte::AnalyteGSL) = findall(id -> connected(cpd, id), analyte)
connected(cpd::CompoundGSL, id::CompoundGSL) = ischaincompatible(cpd, id) && begin
    class1 = haskey(SPDB[:CONNECTION], cpd.class) ? cpd.class : deisomerized(cpd.class)
    class2 = haskey(SPDB[:CONNECTION], id.class) ? id.class : deisomerized(id.class)
    connected(class1, class2) || connected(class2, class1)
end

connected(cls1::ClassGSL, cls2::ClassGSL) = cls2 == cls1 || connected(SPDB[:CONNECTION][cls1], cls2)
connected(cls1::ClassGSL, cls2::Cer) = true
connected(cls1::Cer, cls2::ClassGSL) = false
connected(cls1::Cer, cls2::Cer) = true

isf(cls::ClassGSL) = push!(isf(SPDB[:CONNECTION][cls]), cls)
isf(cls::Cer) = ClassGSL[cls]
