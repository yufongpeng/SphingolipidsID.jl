"""
    ChainSP

Abstract type for all kinds of side chains for sphingolipids.
"""
abstract type ChainSP end
"""
    LCB <: ChainSP

Abstract type for long-chain bases.
"""
abstract type LCB <: ChainSP end
"""
    LCB4 <: LCB

Abstract type for long-chain bases with double bonds and additional oxygen adding up to four.
"""
abstract type LCB4 <: LCB end
"""
    LCB3 <: LCB

Abstract type for long-chain bases with double bonds and additional oxygen adding up to three.
"""
abstract type LCB3 <: LCB end
"""
    LCB2 <: LCB

Abstract type for long-chain bases with double bonds and additional oxygen adding up to two.
"""
abstract type LCB2 <: LCB end
"""
    SPB4 <: LCB4

Sphingosine or sphinganine with double bonds and additional oxygen adding up to four.

When the additional oxygen is three and above, it's assumed that the long-chain base is a phytosphingosine.
"""
struct SPB4 <: LCB4
    carbon::Int
    doublebond::Int
end
"""
    NotPhyto4 <: LCB4

Sphingosine or sphinganine with double bonds and additional oxygen adding up to four.

When the additional oxygen is three and above, it's assumed that the long-chain base is not a phytosphingosine.
"""
struct NotPhyto4 <: LCB4
    carbon::Int
    doublebond::Int
end
"""
    SPB3 <: LCB3

Sphingosine or sphinganine with double bonds and additional oxygen adding up to three.

When the additional oxygen is three and above, it's assumed that the long-chain base is a phytosphingosine.
"""
struct SPB3 <: LCB3
    carbon::Int
    doublebond::Int
end
"""
    NotPhyto3 <: LCB3

Sphingosine or sphinganine with double bonds and additional oxygen adding up to three.

When the additional oxygen is three and above, it's assumed that the long-chain base is not a phytosphingosine.
"""
struct NotPhyto3 <: LCB3
    carbon::Int
    doublebond::Int
end
"""
    SPB2 <: LCB2

Sphingosine or sphinganine with double bonds and additional oxygen adding up to two.
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

N-acyl chain with additional oxygen at an unknown position or no additional oxygen.
"""
struct Acyl <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    Acylα <: ACYL

N-acyl chain with additional oxygen at the α position.
"""
struct Acylα <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    Acylβ <: ACYL

N-acyl chain with additional oxygen at the β position.
"""
struct Acylβ <: ACYL
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    SumChain <: ChainSP

Side chain with only species-level information(sum composition).
"""
struct SumChain <: ChainSP
    carbon::Int
    doublebond::Int
    oxygen::Int
end
"""
    DiChain{S <: LCB, T <: ACYL} <: ChainSP

Two known side chains, one long-chain base, and one N-acyl chain.
"""
struct DiChain{S <: LCB, T <: ACYL} <: ChainSP
    lcb::S
    acyl::T
end

@as_record ChainSP