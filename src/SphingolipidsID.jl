module SphingolipidsID

using CSV, PrettyTables,
        MLStyle, DataPipes, SplitApplyCombine, TypedTables, Dictionaries,
        Statistics, StatsBase, Clustering, Plots, GLM, ThreadsX, LinearAlgebra, ChemistryQuantitativeAnalysis
using AnovaBase: getterms, anova, @formula
using UnitfulMoles: parse_compound, ustrip, @u_str
export SPDB, LIBRARY_POS, FRAGMENT_POS, ADDUCTCODE, CLASSDB,
        # config input
        read_adduct_code, read_class_db, read_ce, class_db_index, set_db!, push_library!,
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
        LCBIS,SPB2IS, SPB3IS, SPB4IS, NotPhyto3IS, NotPhyto4IS, 
        ACYL, Acyl, Acylα, Acylβ,
        ACYLIS, AcylIS, AcylαIS, AcylβIS,
        ChainSP, SumChain, DiChain,
        ChainSPIS, SumChainIS, DiChainIS,
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
        PreIS, MRM, QuantData, QCData, SerialDilution, Quantification, Project, Query,
        # Create library/rule
        library, rule,
        # Data inputs
        read_featuretable, read_featuretable_mzmine3, sort_data!, sort_data, file_order, 
        fill_polarity!, fill_polarity, fill_mz2!, fill_ce!, fill_range!, fill_mz2, fill_ce, fill_range, 
        split_datafile, read_featuretable_masshunter_mrm, filter_combine_features!, rsd, rmae, default_error,
        # PreIS
        preis, preis!, finish_profile!,
        # ID
        apply_rules!, id_product,
        # Score
        nfrags, @cpdsc, @score, calc_score, apply_score!, calc_threshold, apply_threshold!, w_rank, w_nfrags, w_rank_log, w_rank_exp, 
        normalized_sig_diff, abs_sig_diff,
        # Transition
        analytetable_mrm, transitiontable, concurrent_transition, write_transition, 
        read_transition, union_transition!, union_transition, diff_transition!, diff_transition,
        # Quantification
        quantification_mrm, set_quantification_mrm!, set_qcdata_mrm!, qcdata_mrm, qcdata_mrm!, 
        set_serialdilution_mrm!, serialdilution_mrm, serialdilution_mrm!, qualitytable, 
        set_quantdata_mrm!, quantdata_mrm, quantdata_mrm!,
        # Utils/Query
        q!, qand, qor, qnot, spid, analyteid, transitionid, lcb, acyl, acylα, acylβ, reuse,
        new_project, allow_unknown, only_known, ncb, ndb, nox, ndbox, compound_formula, 
        mw, mz, class, chain, sumcomp, rt, cluster, incluster, @ri_str,
        assign_parent!, assign_isf_parent!, mode,
        # RT prediction
        initialize_cluster!, analyte2cluster!, select_cluster!,
        model_cluster!, compare_models, @model, predfn_cluster!,
        expand_cluster!, update_cluster!, replace_cluster!,
        show_cluster, model_rt!, apply_rt!, predict_rt, err_rt, abs_err_rt, rt_correction,
        # Plots
        plot_rt_mw, plotlyjs, gr, histogram_transition

import Base: getproperty, propertynames, show, print, isless, isequal, isempty, keys, length, iterate, getindex, view, firstindex, lastindex, sort, sort!, in,
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
    AbstractCompoundID

Abstract type for compound identification.
"""
abstract type AbstractCompoundID end
"""
    AbstractCompoundSP <: AbstractCompoundID

Abstract type for compound identification with fragments and other information.
"""
abstract type AbstractCompoundSP <: AbstractCompoundID end
"""
    SPID{C <: ClassSP, S <: ChainSP} <: AbstractCompoundID

A type contains only the minimal information of a compound.

* `class`: `ClassSP`.
* `chain`: `ChainSP`.
"""
struct SPID{C <: ClassSP, S <: ChainSP} <: AbstractCompoundID
    class::C
    chain::S
end
"""
    CompoundSPVanilla <: AbstractCompoundSP

A simpler version of `CompoundSP` with only `class`, `chain`, and `fragment`.
"""
mutable struct CompoundSPVanilla <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragment::Table
    CompoundSPVanilla(class::ClassSP, chain::ChainSP, fragment = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragment)
end
"""
    CompoundSP <: AbstractCompoundSP

A compound which has a specific structure after ionization and in-source fragmentation. For example, [M+H]+ and [M+H-H2O]+ of Cer 18:1;O2/18:0 are considered as the same compound. 

* `class`: `ClassSP`.
* `chain`: `ChainSP`.
* `fragment`: `Table` with 4 columns.
    * `ion1`: parent ion.
    * `ion2`: fragment.
    * `source`: the index of the original data in `project.data`.
    * `id`: the index of the original data in the table `project.data[source]`.
* `signal`: `Tuple{Float64, Float64}`. First element is the estimate; second element is the error.
* `state`: `Vector{Int}`. First element represents class; second element represents chain. `1` means this compound was matched; `0` means this compound was uncertain; `-1` means there was contradiction.
* `result`: `Vector{RuleSet}`. First element represents class; second element represents chain.
* `project`: `AbstractProject`; project associated with this compound.
"""
mutable struct CompoundSP <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragment::Table # ion1, ion2, source, id
    signal::Tuple{Float64, Float64}
    state::Vector{Int}
    result::Vector{RuleSet}
    project::AbstractProject
    CompoundSP(class::ClassSP, chain::ChainSP, fragment, signal, state, result, project) = new(class, chain, fragment, signal, state, result, project)
    CompoundSP(class::ClassSP, chain::ChainSP, fragment = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragment, (0.0, 0.0), [0, 0], RuleSet[RuleSet(:missing, EmptyRule()), RuleSet(:missing, EmptyRule())])
end
#CompoundSP(cpd::CompoundSP, chain::ChainSP) = CompoundSP(cpd.class, chain, cpd.fragment, cpd.signal, cpd.state, cpd.result, cpd.project)
"""
    AbstractAnalyteID

Abstract type for analyte identification.
"""
abstract type AbstractAnalyteID end
"""
    AnalyteSP <: AbstractAnalyteID

An analyte which has specific structure in sample but it may become several compounds in MS due to in-source fragmentation.

* `compound`: `Vector{CompoundSP}`. The last element is the true identification of this analyte; the preceeding elements represent in-source fragments.
* `rt`: retention time in minutes.
* `state`: `Vector{Int}`
    1. class
    2. chain
    3. rt
    4. error
    5. isf
    6. total
    7. manual
    `1` indicates that it is qualified; `0` indicates that it is uncertain or not tested; `1` indicates that it is disqualified.
* `cpdsc`: `Pair{Int, Float64}`; score from `analyte.compound` in the form of `index of function in SPDB[:SCORE] => score`
* `score`: `Pair{Int, Float64}`; score from analyte itself in the form of `index of function in SPDB[:SCORE] => score`

This type is considered as `Vector{CompoundSP}` and it supports iteration and most functions for `Vector`. 
"""
mutable struct AnalyteSP <: AbstractAnalyteID
    compound::Vector{CompoundSP}
    rt::Float64
    state::Vector{Int}
    cpdsc::Pair{Int, Float64}
    score::Pair{Int, Float64}
    #manual_check::Int
end
"""
    AnalyteSP(compound::Vector{CompoundSP}, rt::Float64) 

Create a `AnalyteSP` with a vector of `CompoundSP` and their retention time.
"""
AnalyteSP(compound::Vector{CompoundSP}, rt::Float64) = AnalyteSP(compound, rt, repeat([0], 7), 0 => 0.0, 0 => 0.0)
"""
    AnalyteID <: AbstractAnalyteID

A type contains only compounds (in-source fragments) of an analyte and retention time.

* `compound`: `Vector{SPID}`.
* `rt`: `Float64`.
"""
struct AnalyteID <: AbstractAnalyteID
    compound::Vector{SPID}
    rt::Float64
end
"""
    TransitionID

A type represents a transition of an analyte and whether it is a quantifier. 
"""
struct TransitionID
    compound::SPID
    quantifier::Bool
end
"""
    AbstractData
    AbstractRawData <: AbstractData
    AbstractQuantData{T} <: AbstractData

Abstract type for LC-MS data.
"""
abstract type AbstractData end
abstract type AbstractRawData <: AbstractData end
abstract type AbstractQuantData{T} <: AbstractData end
# precursor/product, ion => name, ms/mz => m/z, mw => molecular weight
"""
    PreIS <: AbstractRawData

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
struct PreIS <: AbstractRawData
    table::Table # id, mz1, mz2_id, height, area, collision_energy, FWHM, symmetry, error, isf
    range::Vector{RealIntervals}
    mz2::Vector{Float64}
    polarity::Bool
    config::Dictionary
end
"""
    MRM <: AbstractRawData

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
    * `error`: signal error, optinal.
    * `raw_id`: `id` of original table, optional.
    * `match_id`: `id` of matched transition from `project.quantification`.
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
struct MRM <: AbstractRawData
    table::Table
    mz2::Vector{Float64}
    polarity::Bool
    config::Dictionary
end

"""
    QuantData{T} <: AbstractQuantData{T}

Type containing raw data and tables related to quantification. 

* `raw`: raw data.
* `table`: `AnalysisTable`, data for quantification (`table.area` or `table.height`), and quantification result (`table.estimated_concentration`).
* `config`: `Dictionary` containing config information.
"""
struct QuantData{T} <: AbstractQuantData{T}
    raw::T
    table::AnalysisTable
    config::Dictionary
end

function getproperty(dt::QuantData, p::Symbol)
    if !in(p, fieldnames(QuantData)) && in(p, propertynames(getfield(dt, :table))) 
        getproperty(getfield(dt, :table), p)
    else
        getfield(dt, p)
    end
end
propertynames(dt::QuantData) = Tuple(unique((:raw, :table, :config, propertynames(getfield(dt, :table))...)))

"""
    QCData{T} <: AbstractQuantData{T}

Type for QC samples containing raw data and tables related to quantification. 

* `raw`: raw data.
* `table`: `AnalysisTable`, data for quantification (`table.area` or `table.height`), and quantification result (`table.estimated_concentration`).
* `config`: `Dictionary` containing config information.
"""
struct QCData{T} <: AbstractQuantData{T}
    raw::T
    table::AnalysisTable # id, analyte, FWHM, symmetry, estimate, error, raw_id
    config::Dictionary
end

function getproperty(qc::QCData, p::Symbol)
    if !in(p, fieldnames(QCData)) && in(p, propertynames(getfield(qc, :table))) 
        getproperty(getfield(qc, :table), p)
    else
        getfield(qc, p)
    end
end
propertynames(qc::QCData) = Tuple(unique((:raw, :table, :config, propertynames(getfield(qc, :table))...)))
#QCData(raw::T, table::Table, config::Dictionary) where {T <: Data} = QCData{T}(raw, table, config)
"""
    SerialDilution{T} <: AbstractQuantData{T}

Type for serial dilution samples containing raw data and calibration data. 

* `raw`: raw data.
* `batch`: `Batch` containg analyte settings and calibration curves.
* `config`: `Dictionary` containing config information.
"""
struct SerialDilution{T} <: AbstractQuantData{T}
    raw::T
    batch::Batch
    config::Dictionary
end

function getproperty(sd::SerialDilution, p::Symbol)
    if p in propertynames(getfield(sd, :batch)) 
        getproperty(getfield(sd, :batch), p)
    else
        getfield(sd, p)
    end
end
propertynames(sd::SerialDilution) = Tuple(unique((:raw, :batch, :config, propertynames(getfield(sd, :batch))...)))
#SerialDilution(raw::T, table::Table, level::Vector, config::Dictionary) where {T <: Data} = SerialDilution{T}(raw, table, level, config)

"""
    Quantification

Type containg information of analyte settings, calibration data, qc data, serial dilution data, and sample data.

* `batch`: `Batch` containg analyte settings, calibration data, and sample data.
* `config`: `Dictionary` containing config information. `config[:qcdata]` stores a `QCData`, `config[:serialdilution]` stores a `SerialDilution`, and `config[:quantdata]` stores a `QuantData` (sample data).
"""
struct Quantification
    batch::Batch
    config::Dictionary
    Quantification() = new()
    Quantification(batch::Batch, config::Dictionary) = new(batch, config)
end # Wrap MethodTable

function getproperty(quant::Quantification, p::Symbol)
    if p in propertynames(getfield(quant, :batch)) 
        getproperty(getfield(quant, :batch), p)
    else
        getfield(quant, p)
    end
end
propertynames(quant::Quantification) = (:batch, :config, propertynames(getfield(quant, :batch))...)

"""
    RTCorrection

Type holding rt data and regression line for correcting rt from old batch (data for identification) to new batch (data for quantification).

* `data`: raw data from new batch.
* `table`: old transition table generated from `project.analye` and function `analytetable_mrm`.
* `fn`: `Dictionary` containg regression coefficents for each class.

This is a callable object taking a class and rt or vector of class and a vector of rt as input, and returning new rt or a vector of new rt. If the class is not in `fn`, it will return the original rt. 
"""
struct RTCorrection
    data::AbstractRawData
    table::Table
    fn::Dictionary
end
(fn::RTCorrection)(cls::ClassSP, rt) = [1, rt]'get(fn.fn, cls, [0, 1])
(fn::RTCorrection)(cls::Vector{<: ClassSP}, rt::Vector) = [repeat([1], length(rt)) rt]'get.(Ref(fn.fn), cls, Ref([0, 1]))
function getproperty(rtcor::RTCorrection, p::Symbol)
    if !in(p, fieldnames(RTCorrection)) && in(p, propertynames(getfield(rtcor, :table))) 
        getproperty(getfield(rtcor, :table), p)
    else
        getfield(rtcor, p)
    end
end
propertynames(rtcor::RTCorrection) = Tuple(unique((:data, :table, :fn, propertynames(getfield(rtcor, :table))...)))

"""
    Project <: AbstractProject

An LC-MS analysis project.

* `analyte`: `Vector{AnalyteSP}`.
* `data`: `Vector{<: AbstractData}`. Each data is considered as different experiments, and they will not be compared during identification.
* `appendix`: `Dictionary` containg additional information.
    * `:anion`: the anion used in mobile phase (`:acetate`, `:formate`).
    * `:signal`: `:area` or `:height` for concentration estimation.

This type is considered as `Vector{AnalyteSP}` and it supports iteration and most functions for `Vector`.
"""
mutable struct Project <: AbstractProject
    analyte::Vector{AnalyteSP}
    data::Vector{AbstractData}
    quantification::Quantification
    appendix::Dictionary
end
"""
    Project(; kwargs...)

Create an empty `Project`. All keyword arguments will be put into `project.appendix`. By default, `project.appendix[:anion]` will be `:acetate`; `project.appendix[:signal]` will be `:area`.
"""
function Project(; kwargs...)
    appendix = Dictionary{Symbol, Any}(kwargs)
    get!(appendix, :anion, :acetate)
    get!(appendix, :signal, :area)
    Project(AnalyteSP[], AbstractData[], Quantification(), appendix)
end
"""
    Project(qt::Quantification; kwargs...)

Create a `Project` from `Quantification`. All keyword arguments will be put into `project.appendix`. By default, `project.appendix[:anion]` will be `:acetate`; `project.appendix[:signal]` will be `:area`.
"""
function Project(qt::Quantification; kwargs...)
    appendix = Dictionary{Symbol, Any}(kwargs)
    get!(appendix, :anion, get(qt.appendix, :anion, :acetate))
    get!(appendix, :signal, get(qt.appendix, :signal, :area))
    analyte = isa(qt.analyte, Vector{AnalyteSP}) ? qt.analyte : map(qt.analyte, qt.batch.analytetable.rt) do ana, rt
        isa(ana, AnalyteSP) ? ana : AnalyteSP([CompoundSP(class(ana), chain(ana))], rt)
    end
    Project(analyte, AbstractData[], qt, appendix)
end
"""
    AbstractQuery

Abstract type for query.
"""
abstract type AbstractQuery end
"""
    Query <: AbstractQuery

Query of a project.

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
include("quantification.jl")

# S: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-237.2224
# DHS: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-239.2224
# Phyto(4-OH): Acyl + C3H5NO-422.4004, SPB;3O - CH8NO-267.2329, Acyl + C2H5NO-410.4004, SPB;3O - C2H8NO-255.233
# α-hydroxyl: Acyl - CH2O-337.3476, Acyl + O-383.3531
# β-hydroxyl: SPB;2O + C2H2O - 340.246

end
