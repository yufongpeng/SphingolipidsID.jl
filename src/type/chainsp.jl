"""
    ChainSP

Abstract type for all kinds of side chains for sphingolipids.
"""
abstract type ChainSP end
"""
    LCB <: ChainSP

Abstract type for long chain bases.
"""
abstract type LCB <: ChainSP end
"""
    LCB4 <: LCB

Abstract type for long chain bases with number of double bonds and additional oxygen adding up to 4.
"""
abstract type LCB4 <: LCB end
"""
    LCB3 <: LCB

Abstract type for long chain bases with number of double bonds and additional oxygen adding up to 3.
"""
abstract type LCB3 <: LCB end
"""
    LCB2 <: LCB

Abstract type for long chain bases with number of double bonds and additional oxygen adding up to 2.
"""
abstract type LCB2 <: LCB end
"""
    SPB4 <: LCB4

The type repressenting sphingosine or sphinganine with number of double bonds and additional oxygen adding up to 4.

When the number of addtional oxygen is 3 and above, it's assumed that the long chain base is a phytospingosine.
"""
struct SPB4 <: LCB4
    carbon::Int
    doublebond::Int
end
"""
    NotPhyto4 <: LCB4

The type repressenting sphingosine or sphinganine with number of double bonds and additional oxygen adding up to 4.

When the number of addtional oxygen is 3 and above, it's assumed that the long chain base is not a phytospingosine.
"""
struct NotPhyto4 <: LCB4
    carbon::Int
    doublebond::Int
end
"""
    SPB3 <: LCB3

The type repressenting sphingosine or sphinganine with number of double bonds and additional oxygen adding up to 3.

When the number of addtional oxygen is 3 and above, it's assumed that the long chain base is a phytospingosine.
"""
struct SPB3 <: LCB3
    carbon::Int
    doublebond::Int
end
"""
    NotPhyto3 <: LCB3

The type repressenting sphingosine or sphinganine with number of double bonds and additional oxygen adding up to 3.

When the number of addtional oxygen is 3 and above, it's assumed that the long chain base is not a phytospingosine.
"""
struct NotPhyto3 <: LCB3
    carbon::Int
    doublebond::Int
end
"""
    SPB2 <: LCB2

The type repressenting sphingosine or sphinganine with number of double bonds and additional oxygen adding up to 2.
"""
struct SPB2 <: LCB2
    carbon::Int
    doublebond::Int
end
const NotPhyto = Union{NotPhyto3, NotPhyto4}
"""
    ACYL <: ChainSP

Abstract type for N-acyl chain.
"""
abstract type ACYL <: ChainSP end
"""
    Acyl <: ACYL

The type repressenting N-acyl chain with additional oxygen at unknown position or no additional oxygen.
"""
struct Acyl <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    Acylα <: ACYL

The type repressenting N-acyl chain with additional oxygen at α position.
"""
struct Acylα <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    Acylβ <: ACYL

The type repressenting N-acyl chain with additional oxygen at β position.
"""
struct Acylβ <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    SumChain <: ChainSP

The type repressenting chains with only species-level information(sum composition).
"""
struct SumChain <: ChainSP
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    DiChain{S <: LCB, T <: ACYL} <: ChainSP

The type repressenting two side chains, one long chain base and one N-acyl chain.
"""
struct DiChain{S <: LCB, T <: ACYL} <: ChainSP
    lcb::S
    acyl::T
end

@as_record ChainSP