"""
    AbstractCriteria

Abstarct type for criteria in rules.
"""
abstract type AbstractCriteria end
"""
    IonPlus{T} <: AbstractCriteria

Adding all signals in the attribute `ions`.
"""
struct IonPlus{T} <: AbstractCriteria
    ions::T
end
for fn in (:IonPlus, )
    @eval begin
        $fn(x...) = $fn(x)
    end
end
"""
    IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}} <: AbstractCriteria

Comparing two signals by computing the ratio.
"""
struct IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}} <: AbstractCriteria
    ions::Tuple{S, T}
end
IonComparison(ion1, ion2) = IonComparison((ion1, ion2))
"""
    CalcOx <: AbstractCriteria

Calculating the number of additional oxygen of the attribute `lcb`.
"""
struct CalcOx <: AbstractCriteria
    lcb::LCB
end
@as_record AbstractCriteria
"""
    AbstractRule

Abstarct type for all rules.
"""
abstract type AbstractRule end
"""
    Rule{T} <: AbstractRule

Typical rule for identification.

* `criteria`: the condition to be tested
* `rule`: next rule
"""
struct Rule{T} <: AbstractRule
    criteria::T
    rule
end
"""
    EmptyRule <: AbstractRule

Empty rule.
"""
struct EmptyRule <: AbstractRule end
"""
    RuleUnion <: AbstractRule

Taking the union of the attribute `rules` as the final result.

When all results are empty, the attribute `exception` will be returned.
"""
struct RuleUnion <: AbstractRule
    rules
    exception
    RuleUnion(x...; exception = EmptyRule()) = new(x, exception)
end
"""
    RuleMode <: AbstractRule

Taking the mode of the attribute `rules` as the final result.

When all results are empty, the attribute `exception` will be returned.
"""
struct RuleMode <: AbstractRule
    rules
    exception
    RuleMode(x...; exception = EmptyRule()) = new(x, exception)
end
@as_record AbstractRule
"""
    RuleSet

A type for storing identification mode and previous matched result or on-going matching state.
"""
struct RuleSet
    mode::Symbol
    rule::Union{AbstractRule, ClassSP, ChainSP}
end
"""
    AbstractResult{T}

Abstract type for identification result.
"""
abstract type AbstractResult{T} end
"""
    Result{T} <: AbstractResult{T}

Typical identification result.

* `matched`: a `Bool`; whether there's a match or not.
* `rule`: on-going matching state or final identification.
"""
struct Result{T} <: AbstractResult{T}
    matched::Bool
    rule::T
end
"""
    PartialResult{T, C <: Union{ChainSP, ClassSP}} <: AbstractResult{T}

Identification result including multiple possible structures.

* `matched`: a `Bool`; whether there's a match or not.
* `rule`: on-going matching state.
* `result`: identification which the exact structure is not settled.
"""
struct PartialResult{T, C <: Union{ChainSP, ClassSP}} <: AbstractResult{T}
    matched::Bool
    rule::T
    result::C
end
@as_record AbstractResult
"""
    QueryCommands

Abstract type for query logics.
"""
abstract type QueryCommands end
"""
    QueryCmd <: QueryCommands

Ordinary query commands.
"""
struct QueryCmd <: QueryCommands
    query
end
"""
    QueryAnd <: QueryCommands

Taking the intersection of attribute `qcmd`. 
"""
struct QueryAnd <: QueryCommands
    qcmd::Vector{QueryCommands}
end
"""
    QueryOr <: QueryCommands

Taking the union of attribute `qcmd`. 
"""
struct QueryOr <: QueryCommands
    qcmd::Vector{QueryCommands}
end
"""
    QueryNot{T} <: QueryCommands

Taking the negation of attribute `qcmd`. 
"""
struct QueryNot{T} <: QueryCommands
    qcmd::T
end

abstract type RetentionModelCall end
"""
    ModelCall{T} <: RetentionModelCall

Model for retention time prediction. `T` is `true` if `:cluster` is in formula; otherwise `false`. 
"""
struct ModelCall{T} <: RetentionModelCall
    nmterms
    fnterms
    expr
    fn
end
"""
    CurrentModelCall <: RetentionModelCall

Null struct indicates `project.appendix[:rt_model].fn` or `project.appendix[:cluster_model].fn`.
"""
struct CurrentModelCall <: RetentionModelCall end
(rm::ModelCall)(x) = rm.fn(x)
"""
    RetentionModel{P <: RetentionModelCall}

A wrapper for `RetentionModelCall` and fitted model.
"""
mutable struct RetentionModel{P <: RetentionModelCall}
    fn::P
    model
    RetentionModel(fn::P) where {P <: RetentionModelCall} = new{P}(fn)
    RetentionModel(fn::P, model) where {P <: RetentionModelCall} = new{P}(fn, model)
end
"""
    FunctionalFunction{F, G} <: Function

A type for storing original functions for a functional applying to a function.
* `fl`: a functional.
* `fn`: a function.
* `fln`: a function equivalent to `fl(fn)`
"""
struct FunctionalFunction{F, G} <: Function
    fl::F
    fn::G
    fln
end
(fln::FunctionalFunction)(x...) = fln.fln(x...)

abstract type RealIntervals end
struct EmptyInterval <: RealIntervals end
"""
    RealInterval <: RealIntervals

A type representing a real interval.
"""
struct RealInterval <: RealIntervals
    lowerbound
    upperbound
    leftoperator
    rightoperator
end
#(ri::RealInterval)(x) = between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator)
"""
    UnionInterval{N} <: RealIntervals

A type representing an union of multiple disjoint real intervals.
"""
struct UnionInterval{N} <: RealIntervals
    intervals::NTuple{N, RealInterval}
end
#(ui::UnionInterval)(x) = any(ri -> between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator), ui.intervals)

struct DictionarySerialization{C, T}
    d::Dictionary{Symbol, T}
end

JLD2.writeas(::Type{Dictionary{ClassSP, T}}) where T = DictionarySerialization{ClassSP, T}
JLD2.wconvert(::Type{DictionarySerialization{ClassSP, T}}, v::Dictionary{ClassSP, T}) where T = DictionarySerialization{ClassSP, T}(dictionary(Symbol.(typeof.(keys(v))) .=> values(v)))
JLD2.rconvert(::Type{Dictionary{ClassSP, T}}, v::DictionarySerialization{ClassSP, T}) where T = Dictionary{ClassSP, T}(collect(eval.(Expr.(:call, keys(v.d)))), collect(values(v.d)))

struct ModelCallSerialization{T} <: RetentionModelCall
    nmterms
    fnterms
    expr
end

JLD2.writeas(::Type{ModelCall{T}}) where T = ModelCallSerialization{T}
JLD2.wconvert(::Type{ModelCallSerialization{T}}, v::ModelCall{T}) where T = ModelCallSerialization{T}(v.nmterms, v.fnterms, v.expr)
JLD2.rconvert(::Type{ModelCall{T}}, v::ModelCallSerialization{T}) where T = ModelCall{T}(v.nmterms, v.fnterms, v.expr, eval(Expr(:(->), :tbl, Expr(:block, LineNumberNode(@__LINE__, @__FILE__), v.expr))))