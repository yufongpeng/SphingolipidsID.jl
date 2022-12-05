using SphingolipidsID, DataPipes
files = joinpath.(".\\test\\data\\mzmine", readdir(".\\test\\data\\mzmine"))
fts = featuretable_mzmine.(files);

ces = repeat([[30, 30, 60, 60, 60, 60]], 8)
insert!(ces, 6, [30, 30, 60 ,60 , 60])
insert!(ces, 9, [45, 45, 45, 45])

fts = fill_ce_mzmine!.(fts, ces);
fts = filter_duplicate!.(fts; n = 2);
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
# ID: Class, 1361
apply_rules!(pj; match_mode = :class)
# ID: Chain
apply_rules!(pj; match_mode = :chain)
# Score
@p pj |>
    apply_score!(@score chain 1 -1) |> apply_threshold!(<=(1)) |>
    apply_score!(@score chain 1 (0 + 1) / all) |> apply_threshold!(>=(0.5))
    # apply_score!(@score chain => (analyte, cpd) -> size(cpd.fragments, 1) (0 + 1)) |> apply_threshold!(>=(2))

# Queries
@p pj |> query(GM1) |> query(:rt => (3, 4))
@p pj |> query(HexNAcHex3Cer) |> query(cpd(lcb(18, 1, 2), acyl(16, 0, 0))) |> __[1]
@p pj |> query(acyl(24, 0, 1)) |> query(:mz => (810, 813))
query(pj, cpd(Cer, (36, 1, 2)))

@p pj |>
    query(:class) |>
    query(not(:chain!)) |>
    query((HexCer, cpd(36, 1, 2), cpd(34, 1, 2))) |>
    query(lcb(18, 1, 2)) 

query(pj, (cpd(lcb(18, 1, 2), acyl(24, 0, 1))))

# Partial id
@p pj |> query(HexCer) |> query(:chain!) |> apply_rules!

# MRM
@p pj |> query(not(:class!)) |> query(not(:chain!)) |> generate_mrm(:default, LCB, true)
@p pj |> query(:topsc => 0.5) |> query(not(:class!)) |> query(not(:chain!)) |> generate_mrm(:default, LCB, true) |> write_mrm("lcb_mrm.csv")
@p pj |> query(:class!) |> query(not(:chain!)) |> generate_mrm(:default, LCB, true)
@p pj |> query(:topsc => 0.5) |> query(:class) |> query(not(:chain!)) |> generate_mrm(:default, LCB, true)

# Chain.jl 
# @chain project begin
#     query
#     query
# end