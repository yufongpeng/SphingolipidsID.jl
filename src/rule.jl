"""
    rule(c::ClassSP, x...)

Identification rules.
"""
rule(c::ClassSP, x...) = rule(c)
rule(::Missing) = EmptyRule()
rule(::Nothing) = EmptyRule()
rule(::Cer, sc::SumChain) = nox(sc) < 2 ? Cer() : rule(Cer())
rule(::Cer, sc::DiChain) = ((is4e(sc.lcb) && nox(sc.lcb) > 1) || nox(sc.acyl) > 0) ? rule(Cer()) : Cer()
@rule rule(::Cer) = ISF(ProtonationNLH2O(), Cer()) |> (Cer(), ())
@rule rule(::HexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(ProtonationNLH2O(), Cer())) |> (HexCer(), ())
@rule rule(::SHexCer) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) |> (SHexCer(), ())
@rule rule(::Hex2Cer) = (ISF(ProtonationNLH2O(), Cer()), ISF(ProtonationNLH2O(), Hex2Cer())) |> (
    Hex2Cer(),
    Ion(ProtonationNLH2O(), Cer()) |> (Hex2Cer(), ())
)
@rule rule(::SHexHexCer) = (ISF(ProtonationNLH2O(), Hex2Cer()), ISF(Protonation(), Hex2Cer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) |> (SHexHexCer(), ())
@rule rule(::Hex3Cer) = (ISF(Protonation(), Hex2Cer()), ISF(ProtonationNLH2O(), Hex2Cer()), ISF(ProtonationNLH2O(), Cer())) |> (
    Hex3Cer(),
    Ion(ProtonationNLH2O(), Cer()) |> (Hex3Cer(), ())
)
@rule rule(::HexNAcHex2Cer) = ISF(Protonation(), Hex2Cer()) |> (
    HexNAc_Hex2Cer(),
    Ion(ProtonationNLH2O(), HexNAc()) |> (
        Ion(Protonation(), Hex2Cer()) |> (HexNAc_Hex2Cer(), ()),
        ()
    )
)
@rule rule(::HexNAcHex3Cer_) = RuleUnion(
    ISF(Protonation(), Hex3Cer()) |> (HexNAc_Hex3Cer(), ()),
    ISF(Protonation(), HexNAcHex2Cer()) |> (Hex_HexNAc_Hex2Cer(), ());
    exception = Ion(ProtonationNLH2O(), HexNAc()) |> (
        RuleUnion(
            Ion(Protonation(), Hex3Cer()) |> (HexNAc_Hex3Cer(), ()),
            ISF(Protonation(), HexNAcHex2Cer()) |> (Hex_HexNAc_Hex2Cer(), ())
        ),
        ()
    )
)
@rule rule(::GM4) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) |> (
    GM4(),
    Ion(ProtonationNLH2O(), Cer()) |> (GM4(), ())
)
@rule rule(::GM3) = (ISF(ProtonationNLH2O(), HexCer()), ISF(Protonation(), HexCer()), ISF(ProtonationNLH2O(), Cer()), ISF(Protonation(), Cer())) |> (
    GM3(),
    Ion(ProtonationNLH2O(), Cer()) |> (GM3(), ())
)
@rule rule(::GM2) = Ion(ProtonationNLH2O(), HexNAc()) |> (GM2(), ())

@rule rule(::GM1, mode::Symbol) =
    mode ≡ :ab ? rule_ab(GM1_()) : rule(GM1_())

@rule rule_ab(::GM1) = HexNAcHex |> (GM1a(), ())

@rule rule(::GM1) = RuleUnion(
    (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) |> (GM1a(), ()),
    IonComparison(HexNAcHexNANA, HexNAcHex) |> (
        0           => GM1a(),
        ri"(0, 5)"  => GM1_(GM1a(), GM1b()),
        ri"[5, ∞]"  => GM1b(),
        NaN         => ()
    )
)
@rule rule(::GD3) = NANA |> (
    NANA2 |> (GD3(), ()),
    ()
)
@rule rule(::GD2) = NANA |> (
    NANA2 |> (GD2(), ()),
    ()
)
@rule rule(::GD1, mode::Symbol) =
    mode ≡ :ab ? rule_ab(GD1_()) : rule(GD1_())

rule_ab(::GD1) = NANA |> (
    RuleUnion(
        IonComparison(HexNAcHexNANA, HexNAcHex) |> (
            0           => GD1b(),
            ri"(0, 1)"  => GD1_(GD1b(), GD1a()),
            ri"[1, 10)" => GD1_(GD1a(), GD1b()),
            ri"[10, ∞]" => GD1a()
        ),
        NANA2           |> (GD1b(), ())
    ),
    ()
)
@rule rule(::GD1) = NANA |> (
    NANA2 |> (
        IonComparison(HexNAcHexNANA2, HexNAcHex) |> (
            0           => GD1b(),
            ri"(0, 1)"  => GD1_(GD1b(), GD1c()),
            ri"[1, 10)" => GD1_(GD1c(), GD1b()),
            ri"[10, ∞]" => GD1c()
        ),
        RuleUnion(
            (ISF(Protonation(), GM3()), ISF(ProtonationNLH2O(), GM3())) |> (
                GD1a(),
                ()
            ),
            HexNAcHexNANA2 |> (
                GD1α(),
                ()
            )
        )
    ),
    ()
)
@rule rule(::GT3) = NANA |> (
    NANA23 |> (GT3(), ()),
    ()
)
@rule rule(::GT2) = NANA |> (
    NANA23 |> (GT2(), ()),
    ()
)
@rule rule(::GT1, mode::Symbol) =
    mode ≡ :ab ? rule_ab(GT1_()) : rule(GT1_())

@rule rule_α(::GT1a) = RuleMode(
    Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) |> (
        GT1aα(),
        GT1a()
    ),
    Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) |> (
        GT1a(),
        GT1aα()
    )
)
@rule rule_ab(::GT1) = NANA |> (
    NANA2 |> (
        IonComparison(HexNAcHexNANA2, HexNAcHexNANA) |> (
            0           => GT1b(),
            ri"(0, 1)"  => RuleUnion(GT1b(), rule_α(GT1a())),
            ri"[1, 10)" => RuleUnion(rule_α(GT1a()), GT1b()),
            ri"[10, ∞]" => rule_α(GT1a())
        ),
        HexNAcHexNANA2 |> (
            GT1aα(),
            HexNAcHexNANA |> (
                GT1b(),
                ()
            )
        )
    ),
    ()
)
@rule rule_ab_neg(::GT1) = IonComparison(NANA1, NANA2) |> (
    ri"[0, 0.4)"    => GT1a(),
    ri"[0.4, 7)"    => GT1_(GT1a(), GT1b()),
    ri"[7, ∞]"      => GT1b()
)
@rule rule(::GT1) = NANA |> (
    NANA2 |> (
        HexNAcHexNANA2 |> (
            RuleUnion(
                rule_α(GT1a()),
                IonComparison(HexNAcHexNANA2, HexNAcHexNANA) |> (
                    ri"[0, 10)" => GT1b(),
                    ri"[10, ∞]" => ()
                ),
                IonComparison(IonPlus(HexNAcHexNANA2, HexNAcHexNANA), HexNAcHex) |> (
                    ri"[0, 10)" => GT1c(),
                    ri"[10, ∞]" => ()
                )
            ),
            HexNAcHexNANA |> (
                IonComparison(HexNAcHexNANA, HexNAcHex) |> (
                    ri"[0, 10)" => GT1_(GT1c(), GT1b()),
                    ri"[10, ∞]" => GT1b(),
                ),
                HexNAcHex |> (
                    GT1c(),
                    rule_ab_neg(GT1_()) # Assume a/b only
                )
            )
        ),
        HexNAcHexNANA2 |> (
            GT1aα(),
            HexNAcHexNANA |> (
                GT1b(),
                ()
            )
        )
    )
)
@rule rule(::GQ1, mode::Symbol) =
    mode ≡ :ab ? rule_ab(GQ1_()) : rule(GQ1_())

@rule rule_ab(::GQ1) = NANA |> (
    HexNAcHexNANA2 |> (
        RuleMode(
            Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) |> (
                GQ1bα(),
                GQ1b()
            ),
            Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) |> (
                GQ1b(),
                GQ1bα()
            )
        )
    ),
    ()
)
@rule rule(::GQ1) = NANA |> (
    IonComparison(HexNAcHexNANA2, HexNAcHexNANA) |> (
        0           => GQ1c(),
        ri"(0, 5)"  => RuleUnion(GQ1c(), rule_ab(GQ1_())),
        ri"[5, ∞]"  => rule_ab(GQ1_()),
        NaN         => ()
    ),
    ()
)
rule_ab(::GP1) = ()
@rule rule(::GP1) = NANA |> (
    HexNAcHexNANA2 |> (
        RuleMode(
            Ion(ProtonationNLH2O(), Glycan(HexNAc(), NeuAc())) |> (
                GP1cα(),
                GP1c()
            ),
            Ion(ProtonationNLH2O(), Glycan(Hex(), NeuAc(), NeuAc())) |> (
                GP1c(),
                GP1cα()
            )
        ),
        ()
    ),
    ()
)

@rule rule_acyl(sc::DiChain) = (Ion(LossCH2O(), sc.acyl), Ion(AddO(), sc.acyl)) |> (
    DiChain(sc.lcb, alpha_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc.acyl), nox(sc.acyl))),
    Ion(AddC2H2O(), sc.lcb) |> (
        DiChain(sc.lcb, beta_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc.acyl), nox(sc.acyl))),
        DiChain(sc.lcb, default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc.acyl), nox(sc.acyl)))
    )
)
to_phyto(spb::LCB3) = default_constructor(spb)(ncb(spb), ndb(spb))
to_phyto(spb::LCB4) = default_constructor(spb)(ncb(spb), ndb(spb))
de_phyto(spb::LCB3) = notphyto_constructor(spb)(ncb(spb), ndb(spb))
de_phyto(spb::LCB4) = notphyto_constructor(spb)(ncb(spb), ndb(spb))

@rule rule_phytoacyl0(sc::DiChain) = (Ion(AddC2H5NO(), sc.acyl), Ion(AddC3H5NO(), sc.acyl), Ion(LossCH8NO(), sc.lcb), Ion(LossC2H8NO(), sc.lcb)) |> (
    DiChain(to_phyto(sc.lcb), sc.acyl),
    DiChain(de_phyto(sc.lcb), sc.acyl)
)
@rule rule_phytoacyl(sc::DiChain) = Ion(AddC3H5NO(), sc.acyl) |> (
    rule_acyl(DiChain(to_phyto(sc.lcb), sc.acyl)),
    rule_acyl(DiChain(de_phyto(sc.lcb), sc.acyl))
)
@rule rule(sc::DiChain{<: LCB2}) = IonComparison(Ion(ProtonationNL2H2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 0)), Ion(ProtonationNLH2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 0))) |> (
    ri"[0.0, 1.5)"  => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 0)) |> (
            0 => DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 0)),
            1 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 1))),
            2 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 2)))
        )
    ),
    ri"[1.5, ∞]"    => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 1)) |> (
            0 => DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 0)),
            1 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 1))),
            2 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 2)))
        )
    )
)
@rule rule(sc::DiChain{<: LCB3}) = IonComparison(Ion(ProtonationNL3H2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 0)), Ion(ProtonationNL2H2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 0))) |> (
    ri"[0.0, 1.5)"  => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 0)) |> (
            0 => rule_phytoacyl0(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 0))),
            1 => rule_phytoacyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 1))),
            2 => rule_phytoacyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 0), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc), 2)))
        )
    ),
    ri"[1.5, ∞]"    => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 1)) |> (
            0 => DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 0)),
            1 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 1))),
            2 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 2)))
        )
    )
)
@rule rule(sc::DiChain{<: LCB4}) = IonComparison(Ion(ProtonationNL3H2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 1)), Ion(ProtonationNL2H2O(), default_constructor(sc.lcb)(ncb(sc.lcb), 1))) |> (
    ri"[0.0, 1.5)" => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 1)) |> (
            0 => rule_phytoacyl0(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 0))),
            1 => rule_phytoacyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 1))),
            2 => rule_phytoacyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 1), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 1, 2)))
        )
    ),
    ri"[1.5, ∞]"    => (
        CalcOx(default_constructor(sc.lcb)(ncb(sc.lcb), 2)) |> (
            0 => DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 2), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 2, 0)),
            1 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 2), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 2, 1))),
            2 => rule_acyl(DiChain(default_constructor(sc.lcb)(ncb(sc.lcb), 2), default_constructor(sc.acyl)(ncb(sc.acyl), ndb(sc) - 2, 2)))
        )
    )
)
