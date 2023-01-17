@time using SphingolipidsID, DataPipes, SplitApplyCombine, TypedTables
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

@time fts = filter_duplicate!.(fts; n = 3, err_tol = 0.7);

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
pj # 1950
finish_profile!(pj; err_tol = 0.5) # 1170/1652

apply_rules!(pj) # 1053
# ID: Class
# apply_rules!(pj; match_mode = :class)
# ID: Chain
# apply_rules!(pj; match_mode = :chain)
# Score
@p pj |>
    apply_score!(@cpdscore chain 1 -1) |> apply_threshold!(:chain, <=(1)) |>
    apply_score!(@cpdscore chain 1 (0 + 1) / all) |> apply_threshold!(:chain, >=(0.5))
    # apply_score!(@cpdscore chain 1 1 / all)
    # apply_score!(@cpdscore chain => (analyte, cpd) -> size(cpd.fragments, 1) (0 + 1)) |> apply_threshold!(>=(2))

# Queries
plotlyjs()
q!(pj, spid(GM3, (36, 1, 2)))
aq = @p pj q!(not(:class!)) q!(not(:chain!)) q!(not(:isf!)) reuse
@p aq q!(GM3) q!(:rt => (3, 4))
@p pj q!(GM3) q!(:rt => (3, 4))
@p aq q!(GM3) plot_rt_mw
@p aq q!(Hex2Cer) plot_rt_mw(; groupby = sumcomp)
@p aq q!(Hex3Cer) plot_rt_mw(; groupby = acyl)
@p aq q!(SHexCer) plot_rt_mw(; groupby = sidechain)
@p aq q!(HexNAcHex3Cer) plot_rt_mw(; groupby = lcb)
@p pj q!(HexNAcHex3Cer) q!(DiChain(lcb(18, 1, 2), acyl(16, 0, 0))) 
@p pj q!(acyl(24, 0, 1)) q!(:mz => (810, 813)) plot_rt_mw
@p pj q!(cpd(GM3, 42, 1, 2))
@p pj q!(cpd(Cer, 42, 0, 3))
@p pj q!(:mz => (780, 781)) q!(SHexCer)
@p pj q!(cpd(Hex2Cer, 36, 1, 2))
@p aq |>
    q!((HexCer, cpd(36, 1, 2), cpd(34, 1, 2))) |>
    q!(lcb(18, 1, 2)) 

@p pj q!(lcb(18, 1, 2)) q!(CLS.nana[1]) plot_rt_mw
@p pj q!(lcb(18, 1, 2)) q!(CLS.series.as) plot_rt_mw
@p aq q!(CLS.fg.sulfate) plot_rt_mw
# Partial id
@p pj q!(CLS.fg.nana) q!(:chain!) apply_rules!

# Clustering
generate_clusters!(pj)
plot_rt_mw(pj; clusters = :clusters)
plot_rt_mw(aq)
@p aq |>
    q!(not(GM1)) |> 
    q!(not(SHexHexCer)) |> 
    analytes2clusters!(; min_cluster_size = 10, new = true) |>
    plot_rt_mw(; clusters = :possible)

@p aq q!(HexNAcHex2Cer) analytes2clusters!(; min_cluster_size = 2, radius = 2) plot_rt_mw(; clusters = :possible)
@p aq q!(SHexCer) analytes2clusters!(; min_cluster_size = 2, radius = 2) plot_rt_mw(; clusters = :possible)
@p aq q!(HexCer) analytes2clusters!(; min_cluster_size = 10, radius = 1) plot_rt_mw(; clusters = :possible)
show_clusters(pj)
@p pj select_clusters!(; by = first, new = true, GM3 = 1:2, Hex2Cer = 1:2, Hex3Cer = 1:2, SHexCer = [1, 3], HexNAcHex2Cer = nothing) plot_rt_mw(; clusters = :candidate)
@p pj model_clusters!(:default) plot_rt_mw(; model = true, clusters = :candidate, linewidth = 2)
compare_models(pj, @model(lm(@formula(rt ~ mw))), @model(), @model(lm(@formula(rt ~ mw * mw + cluster))))
compare_models(pj, @model(lm(@formula(rt ~ mw)), :default, lm(@formula(rt ~ mw * cluster))))
@p pj model_clusters!(@model(lm(@formula(rt ~ mw * cluster)))) plot_rt_mw(; model = true, clusters = :candidate, linewidth = 2)
@p pj generate_clusters_prediction!(; HexNAcHex2Cer = Hex3Cer,  SHexHexCer = Hex3Cer) expand_clusters! replace_clusters! plot_rt_mw(; clusters = :clusters)
apply_clusters!(pj)
@p pj q!(not(:class!)) q!(not(:chain!)) q!(:rt) q!(:topsc => 0.5) plot_rt_mw
aq = @p pj q!(not(:class!)) q!(not(:chain!)) q!(:rt) q!(not(:isf!)) 
@time generate_mrm(:default, LCB, true, aq)
# MRM
@p pj q!(not(:class!)) q!(not(:chain!)) q!(:rt) q!(not(:isf!)) q!(:topsc => 0.5) generate_mrm(:default, LCB, true)
@p pj q!(not(:class!)) q!(not(:chain!)) generate_mrm(:default, LCB, true)
t1 = @p pj q!(not(:class!)) q!(not(:chain!)) q!(:topsc => 0.5) generate_mrm(:default, LCB, true)
t2 = @p pj q!((CLS.fg.nana, CLS.series.as)) q!(not(Hex2Cer)) q!(not(:class!)) q!(not(:chain!)) q!(:topsc => 0.5) generate_mrm(:default, LCB, true)
t3 = @p pj q!(not(:class!)) q!(not(:chain!)) q!(:topsc => 0.5) generate_mrm(:default, Ion(ProtonationNL2H2O(), NeuAc()), true)
append!(t1, t2)
append!(t1, t3) |> nMRM
@p t1 write_mrm("mrm_list.csv")
@p pj q!(:class!) q!(not(:chain!)) generate_mrm(:default, LCB, true)
@p pj q!(:topsc => 0.5) q!(:class) q!(not(:chain!)) generate_mrm(:default, LCB, true)

# Chain.jl 
# @chain project begin
#      q!
#      q!
# end

ft = featuretable_masshunter_mrm(".\\test\\data\\agilent\\qc0.csv");
ft1 = filter_duplicate!(deepcopy(ft), n = 3, err_tol = 0.25)
ft2 = filter_duplicate!(deepcopy(ft), n = 6, err_tol = 0.25)

mrm1 = MRM(ft1);
mrm2 = MRM(ft2);

@p pj q!(mrm1) q!(not(:class!)) q!(not(:chain!)) generate_mrm(:default, LCB, true) write_mrm("mrm_list.csv")
@p pj q!(mrm2) q!(not(:class!)) q!(not(:chain!)) generate_mrm(:default, LCB, true)