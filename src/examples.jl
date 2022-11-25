@time using SphingolipidsID, Chain
files = joinpath.(".\\test\\data\\mzmine", readdir(".\\test\\data\\mzmine"))
fts = map(files) do file
    filter_duplicate(add_ce_mzmine!(featuretable_mzmine(file), -1))
end;

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
for (m, r, ft) in zip(ms, rang, fts)
    preis!(pj, ft, r, m, true; rt_tol = 0.1)
end
finish_profile!(pj)
# ID: Class, 1614
apply_rules!(pj, :class)
# ID: Chain
apply_rules!(pj, :chain)
# Score
@chain pj begin
    filter_score!
    filter_score!(@score chain 1 (0 + 1) / all  x >= 0.5)
    filter_score!(@score chain (analyte, cpd) -> nrow(cpd.fragments) (0 + 1)  x >= 2)
end

# Queries
@chain pj begin
    query(not(:class!))
end

@chain pj begin
    query(:class)
    query(HexNAcHex3Cer)
end

@chain pj begin
    query(:class)
    query(HexNAcHex2Cer)
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
    query(:class, HexCer)
    query(:chain!)
    apply_rules!
end

# MRM
@chain pj begin
    query(not(:class!))
    query(not(:chain!))
    generate_mrm(:default, LCB)
end

@chain pj begin
    select_score!
    query(not(:class!))
    query(not(:chain!))
    generate_mrm(:default, LCB)
    nMRM
end

@chain pj begin
    query(:class)
    query(not(:chain!))
    generate_mrm(:default, LCB)
end

@chain pj begin
    select_score!()
    query(:class)
    query(not(:chain!))
    generate_mrm(:default, LCB)
end
