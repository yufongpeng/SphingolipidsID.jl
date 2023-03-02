"""
    IonPlus{T}

The type indicates adding all signals of in the attribute `ions`.
"""
struct IonPlus{T}
    ions::T
end
for fn in (:IonPlus, )
    @eval begin
        $fn(x...) = $fn(x)
    end
end
@as_record IonPlus
"""
    IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}}

This type indicates comparing two signals by computing the ratio.
"""
struct IonComparison{S <: Union{Ion, IonPlus}, T <: Union{Ion, IonPlus}}
    ions::Tuple{S, T}
end
IonComparison(ion1, ion2) = IonComparison((ion1, ion2))
"""
    AbstractRule

Abstarct type for all rules.
"""
abstract type AbstractRule end
"""
    Rule{T} <: AbstractRule

Typical rules for identification.

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

This type indicates taking the union of results of the attribute `rules` as final result.

When all results are empty, the attribute `exception` will be returned.
"""
struct RuleUnion <: AbstractRule
    rules
    exception
    RuleUnion(x...; exception = EmptyRule()) = new(x, exception)
end
"""
    RuleMode <: AbstractRule

This type indicates taking the mode of results of the attribute `rules` as final result.

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

This type is for storing identification mode and previous matched result or ongoing matching state.
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
* `rule`: ongoing matching state or final identification.
"""
struct Result{T} <: AbstractResult{T}
    matched::Bool
    rule::T
end
"""
    PartialResult{T, C <: Union{ChainSP, ClassSP}} <: AbstractResult{T}

Identification result including multiple possible structures.

* `matched`: a `Bool`; whether there's a match or not.
* `rule`: ongoing matching state.
* `result`: identification which some structures might not be specific.
"""
struct PartialResult{T, C <: Union{ChainSP, ClassSP}} <: AbstractResult{T}
    matched::Bool
    rule::T
    result::C
end
@as_record AbstractResult
"""
    CalcOx

A type indicates calculating the number of additional oxygen of the attribute `lcb`.
"""
struct CalcOx
    lcb::LCB
end
"""
    QueryCommands

Abstract type for query logics.
"""
abstract type QueryCommands end
"""
    QueryCmd <: QueryCommands

A type for ordinary query commands.
"""
struct QueryCmd <: QueryCommands
    query
end
"""
    QueryAnd <: QueryCommands

A type for taking intersection of the queries in the attribute `qcmd`. 
"""
struct QueryAnd <: QueryCommands
    qcmd::Vector{QueryCommands}
end
"""
    QueryOr <: QueryCommands

A type for taking union of the queries in the attribute `qcmd`. 
"""
struct QueryOr <: QueryCommands
    qcmd::Vector{QueryCommands}
end
"""
    QueryNot{T} <: QueryCommands

A type for taking negation of the query in the attribute `qcmd`. 
"""
struct QueryNot{T} <: QueryCommands
    qcmd::T
end