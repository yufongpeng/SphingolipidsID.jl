using SphingolipidsID, DataPipes
file1 = joinpath.(".\\test\\data\\mzmine3.3", readdir(".\\test\\data\\mzmine3.3"))
file2 = joinpath.(".\\test\\data\\mzmine3.2", readdir(".\\test\\data\\mzmine3.2"))
fts = @p [file1, file2] map(map(featuretable_mzmine, _)) zip(__...) map(append!(_...)) map(sort_data!);

ces = repeat([[60, 60, 30, 30, 60, 60]], 9)
insert!(ces, 9, [45, 45, 45, 45])
fts = fill_ce!.(fts, ces);

ms = [236.238, 250.253, 284.295, 264.269, 262.253, 278.285, 292.3, 266.285, 274.093, 282.28]
rang = repeat([(400, 1500)], 9)
insert!(rang, 9, (900, 1650))
fts = fill_mz2!.(fts, ms);

fts = filter_duplicate!.(fts; n = 3, err_tol = 0.7);

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
for (r, ft) in zip(rang, fts)
    preis!(pj, ft, r, true; rt_tol = 0.1)
end
finish_profile!(pj; err_tol = 0.5)
# ID: Class, 1047/1112
apply_rules!(pj)
apply_rules!(pj; match_mode = :class)
# ID: Chain
apply_rules!(pj; match_mode = :chain)
# Score
@p pj |>
    apply_score!(@score chain 1 -1) |> apply_threshold!(<=(1)) |>
    apply_score!(@score chain 1 (0 + 1) / all) |> apply_threshold!(>=(0.5))
    #apply_score!(@score chain 1 1 / all)
    # apply_score!(@score chain => (analyte, cpd) -> size(cpd.fragments, 1) (0 + 1)) |> apply_threshold!(>=(2))

# Queries
query(pj, cpd(GM3, (36, 1, 2)))
@p pj  query(GM1)  query(:rt => (3, 4))
@p pj  query(HexNAcHex3Cer)  query(cpd(lcb(18, 1, 2), acyl(16, 0, 0))) 
@p pj  query(acyl(24, 0, 1))  query(:mz => (810, 813))
@p pj  query(cpd(GM3, 42, 1, 2))
@p pj  query(cpd(Cer, 42, 0, 3))
@p pj  query(:mz => (780, 781))  query(SHexCer)
@p pj  query(cpd(Hex2Cer, 36, 1, 2))

@p pj |>
    query(:class) |>
    query(not(:chain!)) |>
    query((HexCer, cpd(36, 1, 2), cpd(34, 1, 2))) |>
    query(lcb(18, 1, 2)) 

@p pj  query(lcb(18, 1, 2))  query(CLS.nana[1])
@p pj  query(lcb(18, 1, 2))  query(CLS.series.as)
@p pj  query(not(:chain!))  query(not(:class!))  query(CLS.fg.sulfate)
# Partial id
@p pj  query(CLS.fg.nana)  query(:chain!)  apply_rules!

# MRM
@p pj  query(not(:class!))  query(not(:chain!))  generate_mrm(:default, LCB, true)
t1 = @p pj  query(not(:class!))  query(not(:chain!))  query(:topsc => 0.5)  generate_mrm(:default, LCB, true)
t2 = @p pj  query((CLS.fg.nana, CLS.series.as)) query(not(Hex2Cer)) query(not(:class!))  query(not(:chain!))  query(:topsc => 0.5) generate_mrm(:default, LCB, true)
t3 = @p pj  query(not(:class!))  query(not(:chain!))  query(:topsc => 0.5)  generate_mrm(:default, Ion(ProtonationNL2H2O(), NeuAc()), true)
append!(t1, t2)
append!(t1, t3) |> nMRM
@p t1 write_mrm("mrm_list.csv")
@p pj  query(:class!)  query(not(:chain!))  generate_mrm(:default, LCB, true)
@p pj  query(:topsc => 0.5)  query(:class)  query(not(:chain!))  generate_mrm(:default, LCB, true)

# Chain.jl 
# @chain project begin
#      query
#      query
# end

ft = featuretable_masshunter_mrm(".\\test\\data\\agilent\\qc0.csv");
ft1 = filter_duplicate!(deepcopy(ft), n = 3, err_tol = 0.25)
ft2 = filter_duplicate!(deepcopy(ft), n = 6, err_tol = 0.25)

mrm1 = MRM(ft1);
mrm2 = MRM(ft2);

@p pj query(mrm1) query(not(:class!))  query(not(:chain!)) generate_mrm(:default, LCB, true) write_mrm("mrm_list.csv")
@p pj query(mrm2) query(not(:class!))  query(not(:chain!)) generate_mrm(:default, LCB, true)