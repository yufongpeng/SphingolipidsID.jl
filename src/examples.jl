@time using SphingolipidsID, Chain
files = joinpath.(".\\test\\data\\mzmine", readdir(".\\test\\data\\mzmine"))
fts = featuretable_mzmine.(files)

ces = repeat([[30, 30, 60, 60, 60, 60]], 8)
insert!(ces, 6, [30, 30, 60 ,60 , 60])
insert!(ces, 9, [45, 45, 45, 45])
ms = [236.238, 250.253, 284.295, 264.269, 262.253, 278.285, 292.3, 266.285, 274.093, 282.28]
rang = repeat([(400, 1500)], 9)
insert!(rang, 9, (900, 1650))

# Custom DB
db = reduce(append!, (
            library(Cer, ["[M+H]+", "[M+H-H2O]+"], [(32:46, 1:2, "O"), (32:46, 0:4, "2O"), (32:46, 0:3, "3O")]),
            library([HexCer, Hex2Cer, Hex3Cer, GM3], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "2O"), (32:46, 0:3, "3O")]),
            library([SHexCer, SHexHexCer, HexNAcHex2Cer, HexNAcHex3Cer], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "2O")]),
            library(GM1, ["[M+H]+", "[M+2H]2+"], (34:2:42, 1:3, "2O"))
           )
        )

# PreIS 
SPDB[:LIBRARY_POS] = db

pj = preis()
for (m, r, ce, ft) in zip(ms, rang, ces, deepcopy(fts))
    ft = filter_duplicate(add_ce_mzmine!(ft, ce), n = 2)
    preis!(pj, ft, r, m, true; rt_tol = 0.1)
end
finish_profile!(pj)
# ID: Class, 1361
apply_rules!(pj, :class)
# ID: Chain
apply_rules!(pj, :chain)
# Score
@chain pj begin
    apply_score!(@score chain 1 -1)
    apply_threshold!(<=(1))
    apply_score!(@score chain 1 (0 + 1) / all)
    apply_threshold!(>=(0.5))
    #filter_score!(@score chain (analyte, cpd) -> nrow(cpd.fragments) (0 + 1)  x >= 2)
end

# Queries
@chain pj begin
    query(GM1)
end

@chain pj begin
    query(HexNAcHex3Cer)
    query(cpd(lcb(18, 1, 2), acyl(16, 0, 0)))
    _[1]
end

@chain pj begin
    query(acyl(24, 0, 1))
end

@chain pj begin
    query(cpd(Cer, (36, 1, 2)))
end

@chain pj begin
    query(cpd(36, 1, 2))
end

@chain pj begin
    query(:class)
    query(not(:chain!))
    query(HexCer, cpd(36, 1, 2), cpd(34, 1, 2))
    query(lcb(18, 1, 2))
end

@chain pj begin
    query(cpd(lcb(18, 1, 2), acyl(24, 0, 1)))
end

@chain pj begin
    query(cpd(42, 1, 3))
end

# Partial id
@chain pj begin
    query(HexCer)
    query(:chain!)
    apply_rules!
end

# MRM
@time @chain pj begin
    query(not(:class!))
    query(not(:chain!))
    generate_mrm(:default, LCB)
end

@chain pj begin
    query(:topsc, 0.5)
    query(not(:class!))
    query(not(:chain!))
    generate_mrm(:default, LCB)
end

@chain pj begin
    query(:class)
    query(not(:chain!))
    generate_mrm(:default, LCB)
end

@chain pj begin
    query(:topsc, 0.5)
    query(:class)
    query(not(:chain!))
    generate_mrm(:default, LCB)
end
