module SphingolipidsID

using CSV, PrettyTables,
        MLStyle, DataPipes, SplitApplyCombine, TypedTables, Dictionaries,
        Statistics, Clustering, Plots, GLM, ThreadsX
using AnovaBase: getterms, anova, @formula
using UnitfulMoles: parse_compound, ustrip, @u_str
export SPDB, LIBRARY_POS, FRAGMENT_POS, ADDUCTCODE, CLASSDB,
        # config input
        read_adduct_code, read_class_db, read_ce, class_db_index, set_db!,
        # Type Class
        ClassSP, SPB, Cer, CerP, SM, HexCer, SHexCer, Hex2Cer, SHexHexCer, Hex3Cer,
        HexNAcHex2Cer, HexNAcHex2Cer_, HexNAc_Hex2Cer, HexNAcHex3Cer, HexNAcHex3Cer_, HexNAc_Hex3Cer, Hex_HexNAc_Hex2Cer, ClasswoNANA,
        GM4, GM3, GM2, GM1, GM1_, GM1a, GM1b,
        GD3, GD2, GD1, GD1_, GD1a, GD1b, GD1c, GD1α,
        GT3, GT2, GT1, GT1_, GT1a, GT1b, GT1c, GT1aα,
        GQ1, GQ1_, GQ1b, GQ1c, GQ1bα,
        GP1, GP1_, GP1c, GP1cα,
        CLS,
        # Type LCB/ACYL
        LCB, LCB2, LCB3, LCB4, SPB2, SPB3, SPB4, NotPhyto, NotPhyto3, NotPhyto4,
        ACYL, Acyl, Acylα, Acylβ,
        Lcb, Nacyl, Nacylα, Nacylβ, NACYL,
        ChainSP, SumChain, DiChain,
        # Type Glycan
        Sugar, Glycan, Hex, HexNAc, NeuAc,
        # Type Ion/Adduct
        Ion, ISF,
        Adduct, Pos, Neg,
        Protonation, ProtonationNLH2O, ProtonationNL2H2O, ProtonationNL3H2O, DiProtonation, TriProtonation,
        AddNH4, AddHNH4, Add2NH4, Sodization, SodizationProtonation, DiSodization,
        Deprotonation, DeprotonationNLH2O, DiDeprotonation, TriDeprotonation, AddOAc, AddHCOO,
        LossCH2O, AddO, AddC2H2O, LossCH8NO, LossC2H8NO, AddC3H5NO, AddC2H5NO, LossCH3,

        # Core Data structures
        CompoundSP, AnalyteSP,
        #IonPlus, IonComparison, RuleMode, RuleUnion, AcylIon,
        Data, PreIS, MRM, Project, Query,
        # Create library/rule
        library, rule,
        # Data inputs
        read_featuretable, read_featuretable_mzmine3, sort_data!, sort_data, file_order, 
        fill_polarity!, fill_polarity, fill_mz2!, fill_ce!, fill_range!, fill_mz2, fill_ce, fill_range, 
        split_datafile, read_featuretable_masshunter_mrm, filter_duplicate!, rsd, re, default_error,
        # PreIS
        preis, preis!, finish_profile!,
        # ID
        apply_rules!, id_product,
        # Score
        nfrags, @cpdsc, @score, calc_score, apply_score!, calc_threshold, apply_threshold!, w_rank, w_nfrags, w_rank_log, w_rank_exp, 
        normalized_sig_diff, abs_sig_diff,
        # Transition
        transitiontable, concurrent_transition, write_transition, read_transition, union_transition!, union_transition, diff_transition!, diff_transition,
        # Utils/Query
        q!, qand, qor, qnot, spid, lcb, acyl, acylα, acylβ, reuse,
        new_project, allow_unknown, only_known, ncb, ndb, nox, ndbox, mw, mz, class, chain, sumcomp, rt, cluster, incluster, @ri_str,
        assign_parent!, assign_isf_parent!,
        # RT prediction
        initialize_clusters!, analytes2clusters!, select_clusters!,
        model_clusters!, compare_models, @model, predfn_clusters!,
        expand_clusters!, update_clusters!, replace_clusters!,
        show_clusters, model_rt!, apply_rt!, predict_rt, err_rt, abs_err_rt,
        # Plots
        plot_rt_mw, plotlyjs, gr, histogram_transition

import Base: show, print, isless, isequal, isempty, keys, length, iterate, getindex, view, firstindex, lastindex, sort, sort!, in,
            union, union!, intersect, setdiff, deleteat!, delete!, push!, pop!, popat!, popfirst!, reverse, reverse!, getproperty, copy, convert

"""
    const SPDB

A `Dict` contains package configurations.

* `:IOSOMER`: `Dict` where keys are deisomerized types and values are a tuple of isomers.
* `:LIBRARY_POS`: `Table`; library for compounds in positive ion mode.
* `:LIBRARY_NEG`: `Table`; library for compounds in negative ion mode.
* `:FRAGMENT_POS`: `Matrix`; library for fragments in positive ion mode.
* `:FRAGMENT_NEG`: `Matrix`; library for fragments in negative ion mode.
* `:CLASSDB`: `Table`; information of class for generating library.
* `:ADDUCTCODE`: `Table`; code used in `:CLASSDB` for each adduct.
* `:NLH2O`: adducts that involve neutral loss of H2O.
* `:CONNECTION`: `Dict` represents structure relationship where the keys are more complicated compounds and the values are compounds losing one monosaccharide.
* `:CE`: `Table`; collision energy for generating MRM transition list.
* `:SCORE`: `Dict` where keys are description of score function and values are the function.
"""
const SPDB = Dict{Symbol, Any}()

include(joinpath("type", "classsp.jl"))
include(joinpath("type", "chainsp.jl"))
include(joinpath("type", "ion.jl"))
include(joinpath("type", "utils.jl"))

"""
    AbstractProject

Abstarct type for projects.
"""
abstract type AbstractProject end

"""
    AbstractCompoundSP

Abstract type for compound identification.
"""
abstract type AbstractCompoundSP end
"""
    SPID{C <: ClassSP, S <: ChainSP}

A `NamedTuple` contains only the minimal information of a compound.

* `class`: `ClassSP`.
* `chain`: `ChainSP`.
"""
const SPID{C <: ClassSP, S <: ChainSP} = @NamedTuple begin
    class::C
    chain::S
end
"""
    SPID(cls::ClassSP, sc::ChainSP)
    SPID(cpd::AbstractCompoundSP)

Function for creating `SPID`.
"""
SPID(cls::ClassSP, sc::ChainSP) = (class = cls, chain = sc)
SPID(cpd::AbstractCompoundSP) = SPID(class(cpd), chain(cpd))
"""
    CompoundSPVanilla <: AbstractCompoundSP

A simpler version of `CompoundSP` with only `class`, `chain`, and `fragments`.
"""
mutable struct CompoundSPVanilla <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragments::Table
    CompoundSPVanilla(class::ClassSP, chain::ChainSP, fragments = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragments)
end
"""
    CompoundSP <: AbstractCompoundSP

A compound which has a specific structure after ionization and in-source fragmentation. For example, [M+H]+ and [M+H-H2O]+ of Cer 18:1;O2/18:0 are considered as the same compound. 

* `class`: `ClassSP`.
* `chain`: `ChainSP`.
* `fragments`: `Table` with 4 columns.
    * `ion1`: parent ion.
    * `ion2`: fragments.
    * `source`: the index of the original data in `project.data`.
    * `id`: the index of the original data in the table `project.data[source]`.
* `signal`: `Tuple{Float64, Float64}`. First element is the estimate; second element is the error.
* `states`: `Vector{Int}`. First element represents class; second element represents chain. `1` means this compound was matched; `0` means this compounds was uncertain; `-1` means there was contradiction.
* `results`: `Vector{RuleSet}`. First element represents class; second element represents chain.
* `project`: `AbstractProject`; project associated with this compound.
"""
mutable struct CompoundSP <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragments::Table # ion1, ion2, source, id
    signal::Tuple{Float64, Float64}
    states::Vector{Int}
    results::Vector{RuleSet}
    project::AbstractProject
    CompoundSP(class::ClassSP, chain::ChainSP, fragments, signal, states, results, project) = new(class, chain, fragments, signal, states, results, project)
    CompoundSP(class::ClassSP, chain::ChainSP, fragments = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragments, (0.0, 0.0), [0, 0], RuleSet[RuleSet(:missing, EmptyRule()), RuleSet(:missing, EmptyRule())])
end
"""
    const CompoundID = Union{<: SPID, <: AbstractCompoundSP}

A Union type for compounds identification.
"""
const CompoundID = Union{<: SPID, <: AbstractCompoundSP}
#CompoundSP(cpd::CompoundSP, chain::ChainSP) = CompoundSP(cpd.class, chain, cpd.fragments, cpd.signal, cpd.states, cpd.results, cpd.project)

"""
    AnalyteSP

An analyte which has specific structure in sample but it may become several compounds in MS due to in-source fragmentation.

* `compounds`: `Vector{CompoundSP}`. The last element is the true identification of this analyte; the preceeding elements represent in-source fragments.
* `rt`: retention time in minutes.
* `states`: `Vector{Int}`
    1. class
    2. chain
    3. rt
    4. error
    5. isf
    6. total
    7. manual
    `1` indicates that it is qualified; `0` indicates that it is uncertain or not tested; `1` indicates that it is disqualified.
* `cpdsc`: `Pair{Int, Float64}`; score from `analyte.compounds` in the form of `index of function in SPDB[:SCORE] => score`
* `score`: `Pair{Int, Float64}`; score from analyte itself in the form of `index of function in SPDB[:SCORE] => score`

This type is considered as `Vector{CompoundSP}` and it supports iteration and most functions for `Vector`. 
"""
mutable struct AnalyteSP
    compounds::Vector{CompoundSP}
    rt::Float64
    states::Vector{Int}
    cpdsc::Pair{Int, Float64}
    score::Pair{Int, Float64}
    #manual_check::Int
end
# ==()
"""
    AnalyteSP(compounds::Vector{CompoundSP}, rt::Float64) 

Create a `AnalyteSP` with a vector of `CompoundSP` and their retention time.
"""
AnalyteSP(compounds::Vector{CompoundSP}, rt::Float64) = AnalyteSP(compounds, rt, repeat([0], 7), 0 => 0.0, 0 => 0.0)

### connect db
"""
    Data

Abstract type for LC-MS data.
"""
abstract type Data end
# precursor/product, ion => name, ms/mz => m/z, mw => molecular weight
"""
    PreIS <: Data

Precursor ion scan data.

* `table`: processed feature table.
    * `id`: a unique `Int`.
    * `mz1`: m/z of parent ion.
    * `mz2_id`: index of fragments in the attribute `mz2`.
    * `height`: peak height.
    * `area`: peak area.
    * `collision_energy`: collision energy (eV).
    * `FWHM`: Full width at half maximum.
    * `symmetry`: peak symmetry.
    * `error`: signal error.
    * `isf`: `1` indicates that it is not a isf; `0` indicates that it is uncertain; `-1` indicates that it is a isf.
* `range`: `Vector{RealIntervals}` representing the scan range of the corresponding `mz2`.
* `mz2`: `Vector{Float64}`, m/z of fragments.
* `polarity`: `Bool`; `true` for positive ion mode; `false` for negative ion mode.
* `config`: `Dictionary` containing config information.
    * `signal`: `:area` or `:height` for concentration estimation.
    * `:mz_tol`: `Float64`. Tolerance of m/z deviation when contructing this object.
    * `:rt_tol`: `Float64`. Tolerance of rt deviation when contructing this object.
    * `:err_tol`: `Float64`. Tolerance of signal deviation when contructing this object.
    * `:est_fn`: `Function`. Function for calculating signal estimate.
    * `:err_fn`: `Function`. Function for calculating error.
    * `:n`: `Int`. Minimal number of detection.
    * `:other_fn`: `Dictionary`. Stores functions for calculating other signal information.
"""
struct PreIS <: Data
    table::Table # id, mz1, mz2_id, height, area, collision_energy, FWHM, symmetry, error, isf
    range::Vector{RealIntervals}
    mz2::Vector{Float64}
    polarity::Bool
    config::Dictionary
end
"""
    MRM <: Data

Multiple reaction monitoring data.

* `table`: processed feature table.
    * `id`: a unique `Int`.
    * `mz1`: m/z of parent ion.
    * `mz2_id`: index of fragments in the attribute `mz2`.
    * `collision_energy`: collision energy (eV).    
    * `height`: peak height.
    * `area`: peak area.
    * `FWHM`: Full width at half maximum.
    * `symmetry`: peak symmetry.
    * `error`: signal error.
* `mz2`: `Vector{Float64}`, m/z of fragments.
* `polarity`: `Bool`; `true` for positive ion mode; `false` for negative ion mode.
* `config`: `Dictionary` containing config information.
    * `signal`: `:area` or `:height` for signal estimation.
    * `:mz_tol`: `Float64`. Tolerance of m/z deviation when contructing this object.
    * `:rt_tol`: `Float64`. Tolerance of rt deviation when contructing this object.
    * `:err_tol`: `Float64`. Tolerance of signal deviation when contructing this object.
    * `:est_fn`: `Function`. Function for calculating signal estimate.
    * `:err_fn`: `Function`. Function for calculating error.
    * `:n`: `Int`. Minimal number of detection.
    * `:other_fn`: `Dictionary`. Stores functions for calculating other signal information.
"""
struct MRM <: Data
    table::Table # id, mz1, mz2_id, height, area, collision_energy, FWHM, symmetry, error
    mz2::Vector{Float64}
    polarity::Bool
    config::Dictionary
end

struct QCData{T} <: Data
    raw::T
    table::Table # id, analyte, FWHM, symmetry, estimate, error, raw_id
    config::Dictionary
end
#QCData(raw::T, table::Table, config::Dictionary) where {T <: Data} = QCData{T}(raw, table, config)

struct SerialDilution{T} <: Data
    raw::T # id, mz1, mz2, height, area, collision_energy, FWHM, symmetry, level
    table::Table # id, analyte, raw_id, r2, model, data_id
    level::Vector{Float64}
    config::Dictionary
end
#SerialDilution(raw::T, table::Table, level::Vector, config::Dictionary) where {T <: Data} = SerialDilution{T}(raw, table, level, config)

# Calibration
# 

"""
    Project <: AbstractProject

An LC-MS analysis project.

* `analytes`: `Vector{AnalyteSP}`.
* `data`: `Vector{Data}`. Each data is considered as different experiments, and they will not be compared during identification.
* `appendix`: `Dictionary` containg additional information.
    * `:anion`: the anion used in mobile phase (`:acetate`, `:formate`).
    * `:signal`: `:area` or `:height` for concentration estimation.

This type is considered as `Vector{AnalyteSP}` and it supports iteration and most functions for `Vector`.
"""
struct Project <: AbstractProject
    analytes::Vector{AnalyteSP}
    data::Vector{Data}
    #anion::Symbol
    #clusters::Dictionary
    #alignment::Int
    appendix::Dictionary
end
"""
    Project(; kwargs...)

Create an empty `Project`. All keyword arguments will be put into `project.appendix`. By default, `project.appendix[:anion]` will be `:acetate`; `project.appendix[:signal]` will be `:area`
"""
function Project(; kwargs...)
    appendix = Dictionary{Symbol, Any}(kwargs)
    get!(appendix, :anion, :acetate)
    get!(appendix, :signal, :area)
    Project(AnalyteSP[], Data[], appendix)
end
"""
    AbstractQuery

Abstract type for query.
"""
abstract type AbstractQuery end
"""
    Query <: AbstractQuery

Query object.

* `project`: queried `Project`.
* `result`: the result of query.
* `query`: the query commands.
* `view`: `Bool` determining whether `result` is a view or a new vector.
"""
mutable struct Query <: AbstractQuery
    project::Project
    result::AbstractVector
    query::QueryCommands
    view::Bool
end
"""
    ReusableQuery <: AbstractQuery

Reusable query object.
"""
struct ReusableQuery <: AbstractQuery
    query::Query
end

include("io.jl")
include("interface.jl")
include("utils.jl")
include("mw.jl")
include("library.jl")
include("rule.jl")
include("data.jl")
include("mrm.jl")
include("preis.jl")
include("rt.jl")
include("plots.jl")
include("ruleID.jl")
include("score.jl")
include("query.jl")
include("transition.jl")

# S: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-237.2224
# DHS: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-239.2224
# Phyto(4-OH): Acyl + C3H5NO-422.4004, SPB;3O - CH8NO-267.2329, Acyl + C2H5NO-410.4004, SPB;3O - C2H8NO-255.233
# α-hydroxyl: Acyl - CH2O-337.3476, Acyl + O-383.3531
# β-hydroxyl: SPB;2O + C2H2O - 340.246

end
