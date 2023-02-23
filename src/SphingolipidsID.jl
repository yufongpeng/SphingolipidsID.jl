module SphingolipidsID

using CSV, PrettyTables,
        MLStyle, DataPipes, SplitApplyCombine, TypedTables, Dictionaries,
        Statistics, Clustering, Plots, GLM, AnovaGLM
using UnitfulMoles: parse_compound, ustrip, @u_str
export SPDB, LIBRARY_POS, FRAGMENT_POS, ADDUCTCODE, CLASSDB,
        # config input
        read_adduct_code, read_class_db, read_ce, class_db_index,
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
        featuretable_mzmine, sort_data!, sort_data, file_order, fill_mz2!, fill_ce!, featuretable_masshunter_mrm, filter_duplicate!, rsd, re,
        # PreIS
        preis, preis!, finish_profile!,
        # ID
        apply_rules!, id_product,
        # Score
        nfrags, @cpdscore, filter_score!, apply_score, apply_score!, apply_threshold, apply_threshold!, select_score!,
        normalized_sig_diff, abs_sig_diff,
        # Data output: MRM
        generate_mrm, nMRM, write_mrm, read_mrm, union_mrm!, union_mrm, diff_mrm!, diff_mrm,
        # Utils/Query
        q!, qand, qor, qnot, spid, lcb, acyl, acylα, acylβ, reuse,
        new_project, mw, mz, class, chain, sumcomp, rt,
        assign_parent!, assign_isf_parent!,
        generate_clusters!, analytes2clusters!, select_clusters!,
        model_clusters!, compare_models, @model, generate_clusters_prediction!,
        expand_clusters!, update_clusters!, replace_clusters!,
        show_clusters, apply_clusters!,
        # Plots
        plot_rt_mw, plotlyjs, gr

import Base: show, print, isless, isempty, keys, length, iterate, getindex, view, firstindex, lastindex, sort, sort!,
            union, union!, deleteat!, delete!, push!, pop!, popat!, popfirst!, reverse, reverse!, getproperty, copy, convert

"""
    const SPDB

A `Dict` contains package configurations.

* `:IOSOMER`: a `Dict` where keys are deisomerized types and values are a tuple of isomers.

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

A `NamedTuple` contains only the minimal information of a compounds.

* `class`: a `ClassSP`.
* `chain`: a `ChainSP`.
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

A simpler version of `CompoundSP`.
"""
mutable struct CompoundSPVanilla <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragments::Table
    CompoundSPVanilla(class::ClassSP, chain::ChainSP, fragments = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragments)
end
"""
    CompoundSP <: AbstractCompoundSP

A type for storing the identification information.

* `class`: a `ClassSP`.
* `chain`: a `ChainSP`.
* `fragments`: a `Table` with 4 columns.
    * `ion1`: parent ion.
    * `ion2`: fragments.
    * `source`: the index of the original data in `project.data`.
    * `id`: the index of the original data in the table `project.data[source]`.
* `area`: a `Tuple`. First is mean area; second is relative error.
* `states`: a `Vector{Int}`. First value repressents class; second value repressents chain. 1 means this compound was matched; 0 means this compounds was uncertain; -1 means there was contradiction.
* `results`: a `Vector{RuleSet}`. First element repressents class; second element repressents chain.
* `project`: a `AbstractProject`.
"""
mutable struct CompoundSP <: AbstractCompoundSP
    class::ClassSP
    chain::ChainSP
    fragments::Table # ion1, ion2, source, id
    area::Tuple{Float64, Float64}
    states::Vector{Int}
    results::Vector{RuleSet}
    project::AbstractProject
    CompoundSP(class::ClassSP, chain::ChainSP, fragments, area, states, results, project) = new(class, chain, fragments, area, states, results, project)
    CompoundSP(class::ClassSP, chain::ChainSP, fragments = Table(ion1 = Ion[], ion2 = ion2, source = Int[], id = Int[])) = new(class, chain, fragments, (0.0, 0.0), [0, 0], RuleSet[RuleSet(:missing, EmptyRule()), RuleSet(:missing, EmptyRule())])
end
"""
    const CompoundID = Union{<: SPID, <: AbstractCompoundSP}

A Union type for compounds identification.
"""
const CompoundID = Union{<: SPID, <: AbstractCompoundSP}
#CompoundSP(cpd::CompoundSP, chain::ChainSP) = CompoundSP(cpd.class, chain, cpd.fragments, cpd.area, cpd.states, cpd.results, cpd.project)

"""
    AnalyteSP

A type repressenting an analyte. 

* `compounds`: a `Vector{CompoundSP}`. The last element is the true compound identification of this analyte; the preceeding elements repressent in-source fragmentation.
* `rt`: retention time in minutes.
* `states`: a `Vector{Int}`. 
    1. class
    2. chain
    3. rt
    4. error
    5. isf
* `scores` a `Vector{Float64}`
    1. class
    2. chain
* `manual_check`: 1 means checked and approved; 0 means unchecked, -1 means checked but disapproved.
"""
mutable struct AnalyteSP
    compounds::Vector{CompoundSP}
    rt::Float64
    states::Vector{Int}
    scores::Vector{Float64}
    manual_check::Int
end
# ==()

### connect db
"""
    Data

Abstract type for LC-MS data.
"""
abstract type Data end
# precursor/product, ion => name, ms/mz => m/z, mw => molecular weight
"""
    PreIS <: Data

A type for precursor ion scan data.

* `raw`: a `Table`; raw data processed by MZMINE.
    * `id`: a unique `Int`.
    * `mz1`: m/z of parent ion.
    * `mz2_id`: index of fragments in the attribute `mz2`.
    * `height`: peak height.
    * `area`: peak area.
    * `collision_energy`: collision energy (eV).
    * `FWHM`: Full width at half maximum.
    * `symmetry`: peak symmetry.
    * `error`: relative error of `area`.
    * `alignment`: index of `analyte` for alignment.
    * `isf`: 1 means not a isf; 0 means uncertain; -1 means a isf.
* `range`: a `Vector{Tuple{Float64, Float64}}`. Each tuple repressents the scan range of the correspond `mz2`.
* `mz2`: a `Vector{Float64}`. m/z of fragments.
* `mz_tol`: a `Float64`. Tolerance of m/z deviation when contructing this object.
* `polarity`: a `Bool`; `true` for positive ion mode; `false` for negative ion mode.
* `additional`: a `Dict` containing additional information.
"""
struct PreIS <: Data
    raw::Table #id, mz1, mz2_id, area, ev, rt, alignment, isf
    range::Vector{Tuple{Float64, Float64}}
    mz2::Vector{Float64}
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end
"""
    MRM <: Data

A type for multiple reaction monitoring data.

* `raw`: a `Table`; raw data processed by MZMINE.
    * `id`: a unique `Int`.
    * `mz1`: m/z of parent ion.
    * `mz2_id`: index of fragments in the attribute `mz2`.
    * `height`: peak height.
    * `area`: peak area.
    * `collision_energy`: collision energy (eV).
    * `FWHM`: Full width at half maximum.
    * `symmetry`: peak symmetry.
    * `error`: relative error of `area`.
    * `alignment`: index of `analyte` for alignment.
    * `isf`: 1 means not a isf; 0 means uncertain; -1 means a isf.
* `mz2`: a `Vector{Float64}`. m/z of fragments.
* `mz_tol`: a `Float64`. Tolerance of m/z deviation when contructing this object.
* `polarity`: a `Bool`; `true` for positive ion mode; `false` for negative ion mode.
* `additional`: a `Dict` containing additional information.
"""
struct MRM <: Data
    raw::Table #id, mz1, mz2_id, area, ev, rt, alignment, isf
    mz2::Vector{Float64}
    mz_tol::Float64
    polarity::Bool
    additional::Dict
end

"""
    Project <: AbstractProject

The type containg all information of a analysis project.

* `analytes`: a `Vector{AnalyteSP}`.
* `data`: a `Vector{Data}`.
* `anion`: the anion used in mobile phase (`:acetate`, `:formate`).
* `clusters`: a `Dictionary` containg clusters of analytes.
* `alignment`: the index of data which all other data was aligned to.
* `appendix`: a `Dictionary` containg additional information.
"""
struct Project <: AbstractProject
    analytes::Vector{AnalyteSP}
    data::Vector{Data}
    anion::Symbol
    clusters::Dictionary
    alignment::Int
    appendix::Dictionary
end

"""
    AbstractQuery

Abstract type for query.
"""
abstract type AbstractQuery end
"""
    Query <: AbstractQuery

The type for query.

* `project`: queried `Project`.
* `result`: the result of query.
* `query`: the query commands.
* `view`: a `Bool` determining whether `result` is a view or a new vector.
"""
mutable struct Query <: AbstractQuery
    project::Project
    result::AbstractVector
    query::QueryCommands
    view::Bool
end
"""
    ReUseable <: AbstractQuery

A type of query object which is reusable.
"""
struct ReUseable <: AbstractQuery
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
include("mrmtable.jl")

# S: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-237.2224
# DHS: Acyl + C2H3N-392.3898, SPB;2O - C2H7NO-239.2224
# Phyto(4-OH): Acyl + C3H5NO-422.4004, SPB;3O - CH8NO-267.2329, Acyl + C2H5NO-410.4004, SPB;3O - C2H8NO-255.233
# α-hydroxyl: Acyl - CH2O-337.3476, Acyl + O-383.3531
# β-hydroxyl: SPB;2O + C2H2O - 340.246

end
