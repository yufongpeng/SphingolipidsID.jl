using SphingolipidsID, DataPipes, SplitApplyCombine, TypedTables; import PlotlyJS
using Test
import Base: show

test_show(x) = show(IOBuffer(), x)
macro test_error(x)
    return quote
        try 
            $x
            true
        catch e
            false
        end
    end
end

global dir = joinpath(dirname(@__FILE__), "data")

@testset "SphingolipidsID.jl" begin
    @testset "Input" begin
        global fts = @p joinpath(dir, "preis") readdir map(joinpath(dir, "preis", _)) map(read_featuretable_mzmine3) map(sort_data!)
        @test (@p fts map(length)) == [330, 270, 354, 1437, 482, 127, 157, 334, 85, 584]
        @test propertynames(fts[1]) == (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile)
    end
    @testset "Library" begin
        set_db!(:LIBRARY_POS, 
            reduce(append!, 
                (
                    library(Cer, ["[M+H]+", "[M+H-H2O]+"], [(32:46, 1:2, "O"), (32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
                    library([HexCer, Hex2Cer, Hex3Cer, GM3], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "O2"), (32:46, 0:3, "O3")]),
                    library([SHexCer, SHexHexCer, HexNAcHex2Cer, HexNAcHex3Cer], ["[M+H]+", "[M+H-H2O]+"], [(32:46, 0:4, "O2")]),
                    library(GM1, ["[M+H]+", "[M+2H]2+"], (34:2:42, 1:3, "O2"))
                )
            )
        )
        @test length(SPDB[:LIBRARY_POS]) == 2040
    end
    @testset "Preprocessing" begin
        file_order.(fts)
        ces = repeat([[60, 60, 30, 30, 60, 60]], 9)
        insert!(ces, 9, [45, 45, 45, 45])
        global fts_ce = fill_ce!.(fts, ces)
        global fts_mz2 = fill_mz2!.(fts_ce, [236.238, 250.253, 284.295, 264.269, 262.253, 278.285, 292.3, 266.285, 274.093, 282.28])
        global fts_unique = filter_duplicate!.(fts_mz2; n = 2, err_tol = 0.7)
        @test @p fts_unique map(all(getproperty(_, :error) .< 0.7)) all
    end
    @testset "PreIS" begin
        rang = repeat([(400, 1500)], 9)
        insert!(rang, 9, (900, 1650))
        global pj = Project()
        for (r, ft) in zip(rang, fts_unique)
            preis!(pj, ft, r, true; rt_tol = 0.1)
        end
        @test length(pj) == 1763
        @test @test_error test_show(pj)
        finish_profile!(pj; err_tol = 0.5)
        @test length(pj) == 1509
    end
    @testset "ID" begin
        apply_rules!(pj)
        @test length(pj) == 1509
    end
    @testset "Criteria" begin
        @p pj |>
            apply_score!(@cpdsc chain 1 -1) |> apply_threshold!(<=(1); cpd = true) |>
            apply_score!(@cpdsc chain 1 (0 + 1) / all) |> apply_threshold!(>=(0.5); cpd = true)
        @test length(pj) == 1509
    end
    @testset "Query" begin
        q!(pj, spid(GM3, (36, 1, 2)))
        global aq = @p pj q!(qnot(:class!)) q!(qnot(:chain!)) q!(qnot(:isf!)) reuse
        @test length(aq) == 574
        @p aq q!(GM3) q!(:rt => (3.5, 4))
        @p aq q!(Cer) q!(:rt => (7, 7.5))
        @p pj q!(HexNAcHex3Cer) q!(chain(lcb(18, 1, 2), acyl(16, 0, 0)))
        @p pj q!(acyl(24, 0, 1)) q!(:mz => (810, 813))
        @p aq |>
            q!(qor(GM1, chain(36, 1, 2), chain(34, 1, 2))) |>
            q!(qnot(:class!))
        @p pj q!(lcb(18, 1, 2)) q!(CLS.nana[1])
    end

    @testset "Clustering" begin
        @p pj initialize_clusters! 
        @test length(pj.clusters) == 10
        @p aq |>
            q!(qnot(GM1)) |>
            q!(qnot(SHexHexCer)) |>
            analytes2clusters!(; min_cluster_size = 10, new = true) |>
            
        @test isempty(pj.appendix[:clusters_possible][SHexCer()])
        @test length(pj.appendix[:clusters_possible]) == 8
        @p aq q!(SHexCer) analytes2clusters!(; min_cluster_size = 2, radius = 2) 
        @p aq q!(HexCer) analytes2clusters!(; min_cluster_size = 10, radius = 1) 
        @test length(pj.appendix[:clusters_possible][SHexCer()]) == 3
        @p pj select_clusters!(; by = first, new = true, Hex2Cer = 1:2, Hex3Cer = 1:2, SHexCer = [1, 3], HexNAcHex2Cer = nothing)
        @test length(pj.appendix[:clusters_candidate]) == 7
        @p pj model_clusters! 
        compare_models(pj, @model(lm(@formula(rt ~ mw)), (), lm(@formula(rt ~ mw * cluster))))
        @p pj predfn_clusters!(; HexNAcHex2Cer = Hex3Cer,  SHexHexCer = Hex3Cer) expand_clusters! replace_clusters! 
        @test length(pj.clusters) == 9
        apply_clusters!(pj)
        @p pj apply_score!(@score all) apply_threshold!(>=(3))
        global valid = @p pj q!(qnot(:class!)) q!(qnot(:chain!)) q!(:rt) q!(qnot(:isf!)) q!(:total) reuse
        @test length(valid) == 455
    end

    @testset "Plots" begin
        plotlyjs()
	    @test @test_error @p valid q!(GM3) plot_rt_mw
	    @test @test_error @p valid q!(Hex2Cer) plot_rt_mw(; groupby = sumcomp)
	    @test @test_error @p valid q!(Hex3Cer) plot_rt_mw(; groupby = acyl)
	    @test @test_error @p valid q!(SHexCer) plot_rt_mw(; groupby = chain)
	    @test @test_error @p valid q!(HexNAcHex3Cer) plot_rt_mw(; groupby = lcb)
        @test @test_error plot_rt_mw(pj; clusters = :clusters)
        @test @test_error plot_rt_mw(pj; clusters = :possible)
        @test @test_error plot_rt_mw(pj; clusters = :candidate)
        @test @test_error plot_rt_mw(pj; model = true, clusters = :candidate, linewidth = 2)
    end
    @testset "MRM" begin
        global translist1 = transitiontable(valid, :default, LCB, true) 
        @test length(translist1) == 312
        append!(translist1, transitiontable(valid, :default, Ion(ProtonationNL2H2O(), NeuAc()), true)) 
        @test length(translist1) == 340
        @test @test_error histogram_transition(translist1)
        global ft_mrm = read_featuretable(joinpath(dir, "mrm", "qc0.csv"), :mh_mrm; silencewarnings = true)
        @test length(ft_mrm) == 3959
        ft_mrm = filter_duplicate!(ft_mrm, n = 3, err_tol = 0.25)
        @test length(ft_mrm) == 469
        global mrm = MRM(ft_mrm)
    end
    @testset "Utils and Interfaces"  begin
        
    end
    @testset "Output" begin
        @test @test_error write_transition(joinpath(dir, "processed", "translist1.csv"), translist1)
    end
end
