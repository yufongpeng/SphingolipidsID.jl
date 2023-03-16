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