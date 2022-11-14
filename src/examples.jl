using SphingolipidsID, Chain
files = joinpath.(".\\test\\data\\mzmine", readdir(".\\test\\data\\mzmine"))
fts = map(files) do file
    filter_duplicate(add_ce_mzmine!(featuretable_mzmine(file), -1))
end;

ms = [236.238, 250.253, 284.295, 264.269, 262.253, 278.285, 292.3, 266.285, 274.093, 282.28]
range = repeat([(400, 1500)], 9)
insert!(range, 9, (900, 1650))

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
for (m, r, ft) in zip(ms, range, fts)
    preis!(pj, ft, r, m, true; rt_tol = 0.1)
end
finish_profile!(pj)

# ID: Class
apply_rules!(pj, :class)

# ID: Chain
apply_rules!(pj, :chain)

# Queries
@chain pj begin
    query(:class, GM3)
    query(:sum, (36, 1, 2))
end

@chain pj begin
    query(:class)
    query(:chain_wait)
end

@chain pj begin
    query(:class, HexCer)
    query(:chain_fail)
end

# Partial id
@chain pj begin
    query(:class, HexCer)
    query(:chain_fail)
    apply_rules!(;chain_fail = :del)
end