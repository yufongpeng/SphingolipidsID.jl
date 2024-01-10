# public
"""
    nfrags(cpd::AbstractCompoundSP)

Number of fragments of a compund, i.e., `size(cpd.fragment, 1)`.
"""
nfrags(cpd::AbstractCompoundSP) = size(cpd.fragment, 1)
"""
    concurrent_transition(tbl::Table; end_time = maximum(tbl.rt .+ tbl.Δrt))

Number of transitions at each time bins. The `rt` represents the start of bin. `end_time` is the end of aquisition.
"""
function concurrent_transition(tbl::Table; end_time = nothing)
    tbl_s = Table((rt = tbl.rt .- tbl.Δrt ./ 2, count = repeat([1], length(tbl.rt))))
    tbl_e = Table((rt = tbl.rt .+ tbl.Δrt ./ 2, count = repeat([-1], length(tbl.rt))))
    tbl_t = vcat(tbl_s, tbl_e)
    tbl_t = @p tbl_t group(getproperty(:rt)) map(getproperty(:count)) map(sum) pairs map((rt = _.first, count = _.second)) Table
    sort!(tbl_t, :rt)
    tbl_t.count .= accumulate(+, tbl_t.count)
    end_time = isnothing(end_time) ? last(tbl_t.rt) * 1.1 : max(end_time, last(tbl_t.rt))
    push!(tbl_t, (rt = end_time, count = 0))
    tbl_t
end

"""
    new_project(aquery::AbstractQuery; copy_data = false)
    new_project(project::Project; copy_data = false, analyte = copy_wo_project.(project.analyte))

Create a new project from a query or project with given analyte. `copy_data` determines whether copying the data or not.
"""
new_project(aquery::AbstractQuery; copy_data = false) = new_project(aquery.project; copy_data, analyte = aquery.view ? copy_wo_project.(aquery.result) : aquery.result)
function new_project(project::Project; copy_data = false, analyte = copy_wo_project.(project.analyte))
    appendix = empty(project.appendix)
    insert!(appendix, :anion, project.appendix[:anion])
    insert!(appendix, :signal, project.appendix[:signal])
    project_new = Project(analyte, copy_data ? deepcopy(project.data) : project.data, copy_data ? deepcopy(project.quantification) : project.quantification, appendix)
    for analyte in project_new
        for cpd in analyte
            cpd.project = project_new
        end
    end
    project_new
end
"""
    between(num::Number, range; lop = <=, rop = <=)
    between(num::Number; low, up, lop = <=, rop = <=)
    between(num::Number, value, tol; lop = <=, rop = <=)

Determine whether `num` lies in given range. 

The lower bound is `first(range)`, `low`, or `value - tol`; the upper bound is `last(range)`, `up`, or `value + tol`. 

`lop` and `rop` determine the operators for the lower bound and upper bound respectively.
"""
between(num::Number, range; lop = <=, rop = <=) = lop(first(range), num) && rop(num, last(range))
between(num::Number; low, up, lop = <=, rop = <=) = lop(low, num) && rop(num, up)
between(num::Number, value, tol; lop = <=, rop = <=) = lop(value - tol, num) && rop(num, value + tol)
"""
    @ri_str -> RealInterval

Create a `RealInterval` by mathematical real interval notation.

# Examples
```julia
julia> f = ri"[4, 7)"
(::SphingolipidsID.RealInterval) (generic function with 1 method)

julia> f(4)
true

julia> f(7)
false
```
"""
macro ri_str(expr)
    return real_interval(expr)
end
function real_interval(expr)
    lc, lv, rv, rc = match(r" *([\(\[]) *([+-]*[\d∞Inf]*\.*\d*) *, *([+-]*[\d∞Inf]*\.*\d*) *([\)\]]) *", expr)
    lop = @match lc begin
        "[" => <=
        "(" => <
    end
    rop = @match rc begin
        "]" => <=
        ")" => <
    end
    lv = @match lv begin
        "+∞"    => Inf
        "+Inf"  => Inf
        "∞"     => Inf
        "Inf"   => Inf
        "-∞"    => -Inf
        "-Inf"  => -Inf
        if occursin(".", lv) end   => parse(Float64, lv)
        _       => parse(Int, lv)
    end
    rv = @match rv begin
        "+∞"    => Inf
        "+Inf"  => Inf
        "∞"     => Inf
        "Inf"   => Inf
        "-∞"    => -Inf
        "-Inf"  => -Inf
        if occursin(".", rv) end   => parse(Float64, rv)
        _       => parse(Int, rv)
    end
    real_interval(lv, rv, lop, rop)
end
function real_interval(lb, ub, lop, rop)
    lop(lb, ub) || return EmptyInterval()
    rop(lb, ub) || return EmptyInterval()
    RealInterval(lb, ub, lop, rop)
end

"""
    intersection(range...)
    intersection(range::AbstractVector)

Return the intersection of all ranges in `ranges`. The range is a tuple `(low, up)`.
"""
intersection(ranges...) = (maximum(first(r) for r in ranges), minimum(last(r) for r in ranges))
intersection(ranges::AbstractVector) = (maximum(first(r) for r in ranges), minimum(last(r) for r in ranges))
"""
    calc_rt(analyte::AnalyteSP)

Calculate retention time.
"""
calc_rt(analyte::AnalyteSP) = mean(mean(query_data(cpd.project, fragment.source, fragment.id).rt for fragment in cpd.fragment) for cpd in analyte)
"""
    state_id

Return the index for certain evaluation parameters in the states of an analyte.
"""
state_id = @λ begin
    :class  => 1
    :chain  => 2
    :rt     => 3
    :error  => 4
    :isf    => 5
    :total  => 6
    :manual => 7
end
# user interface for query
"""
    qand(qcmds...)

Create `QueryAnd`; taking the intersection of `qcmds`.
"""
qand(qcmds...) = QueryAnd(collect(map(qcmd, qcmds)))
"""
    qor(qcmds...)

Create `QueryOr`; taking the union of `qcmds`.
"""
qor(qcmds...) = QueryOr(collect(map(qcmd, qcmds)))
"""
    qnot(qcmd::QueryCommands)
    qnot(qcmd::QueryNot)
    qnot(qcmd)

Create `QueryNot`; taking the negation of `qcmd`.
"""
qnot(qcmd::QueryCommands) = QueryNot(qcmd)
qnot(qcmd::QueryNot) = qcmd.qcmd
qnot(qcmd) = QueryNot(QueryCmd(qcmd))
"""
    qcmd(qcmd::QueryCommands)
    qcmd(qcmd)
    qcmd(query, neg)

Wrap qcmd into `QueryCmd` or `QueryNot`.
"""
qcmd(qcmd::QueryCommands) = qcmd
qcmd(qcmd) = QueryCmd(qcmd)
qcmd(query, neg) = neg ? qnot(query) : qcmd(query)

"""
    lcb(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0)

Create `LCB` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
function lcb(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0)
    if n13C == 0 && nD == 0 
        Cons = @match d + n begin
            2 => SPB2
            3 => SPB3
            4 => SPB4
        end
        Cons(c, d)
    else
        Cons = @match d + n begin
            2 => SPB2IS
            3 => SPB3IS
            4 => SPB4IS
        end
        Cons(c, d, (; n13C, nD))
    end
end
"""
    acyl(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0)

Create `Acyl` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acyl(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0) = n13C == 0 && nD == 0 ? Acyl(c, d, n) : AcylIS(c, d, n, (; n13C, nD))
"""
    acylα(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0)

Create `Acylα` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acylα(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0) = n13C == 0 && nD == 0 ? Acylα(c, d, n) : AcylαIS(c, d, n, (; n13C, nD))
"""
    acylβ(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0)

Create `Acylβ` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acylβ(c::Int, d::Int, n::Int; n13C::Int = 0, nD::Int = 0) = n13C == 0 && nD == 0 ? Acylβ(c, d, n) : AcylβIS(c, d, n, (; n13C, nD))
"""
    chain(spb::LCB, nacyl::ACYL)
    chain(c::Int, n::Int, o::Int; n13C::Int = 0, nD::Int = 0)

Create `ChainSP` with long-chain base and N-acyl or carbons, double bonds, and addtional oxygens.
"""
chain(spb::LCB, nacyl::ACYL) = DiChain(spb, nacyl)
chain(c::Int, n::Int, o::Int; n13C::Int = 0, nD::Int = 0) = n13C == 0 && nD == 0 ? SumChain(c, n, o) : SumChainIS(c, d, n, (; n13C, nD))
"""
    spid(cpd::AbstractAnalyteID)
    spid(cpd::AbstractCompoundID)
    spid(transition::TransitionID) 
    spid(cls::Type{<: ClassSP}, sc::ChainSP)
    spid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL)
    spid(cls::Type{<: ClassSP}, sum::NTuple{3, Int}; n13C::Int = 0, nD::Int = 0)
    spid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int; n13C::Int = 0, nD::Int = 0)

Create `SPID` with given class and chain information.
"""
spid(cpd::AbstractAnalyteID) = spid(last(cpd))
spid(cpd::AbstractCompoundID) = SPID(class(cpd), chain(cpd))
spid(transition::TransitionID) = SPID(class(transition), chain(transition))
spid(cls::Type{<: ClassSP}, sc::ChainSP) = SPID(cls(), sc)
spid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL) = SPID(cls(), DiChain(spb, nacyl))
spid(cls::Type{<: ClassSP}, sum::NTuple{3, Int}; n13C::Int = 0, nD::Int = 0) = SPID(cls(), chain(sum...; n13C, nD))
spid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int; n13C::Int = 0, nD::Int = 0) = SPID(cls(), chain(c, n, o; n13C, nD))

analyteid(cpd::AbstractCompoundID, rt::Float64) = AnalyteID([SPID(class(cpd), chain(cpd))], rt)
analyteid(cls::Type{<: ClassSP}, sc::ChainSP, rt::Float64) = AnalyteID([SPID(cls(), sc)], rt)
analyteid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL, rt::Float64) = AnalyteID([SPID(cls(), DiChain(spb, nacyl))], rt)
analyteid(cls::Type{<: ClassSP}, sum::NTuple{3, Int}, rt::Float64; n13C::Int = 0, nD::Int = 0) = AnalyteID([SPID(cls(), chain(sum...; n13C, nD))], rt)
analyteid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int, rt::Float64; n13C::Int = 0, nD::Int = 0) = AnalyteID([SPID(cls(), chain(c, n, o; n13C, nD))], rt)

transitionid(cpd::SPID, quantifier::Bool) = TransitionID(cpd, quantifier)
# isomer
"""
    deisomerized(::ClassSP)
    deisomerized(::Type{<: ClassSP})

Return the object or type representing all isomers.
"""
deisomerized(::GM1) = GM1()
deisomerized(::GD1) = GD1()
deisomerized(::GT1) = GT1()
deisomerized(::GQ1) = GQ1()
deisomerized(::GP1) = GP1()
deisomerized(::HexNAcHex2Cer) = HexNAcHex2Cer()
deisomerized(::HexNAcHex3Cer) = HexNAcHex3Cer()
deisomerized(x::ClassSP) = x
deisomerized(::Type{<: GM1}) = GM1
deisomerized(::Type{<: GD1}) = GD1
deisomerized(::Type{<: GT1}) = GT1
deisomerized(::Type{<: GQ1}) = GQ1
deisomerized(::Type{<: GP1}) = GP1
deisomerized(::Type{<: HexNAcHex2Cer}) = HexNAcHex2Cer
deisomerized(::Type{<: HexNAcHex3Cer}) = HexNAcHex3Cer
deisomerized(x::Type{<: ClassSP}) = x
"""
    hasisomer(::ClassSP)
    hasisomer(::Type)

Determine whether the object or type has isomers.
"""
hasisomer(::GM1_) = true
hasisomer(::GD1_) = true
hasisomer(::GT1_) = true
hasisomer(::GQ1_) = true
hasisomer(::GP1_) = true
hasisomer(::HexNAcHex2Cer_) = true
hasisomer(::HexNAcHex3Cer_) = true
hasisomer(::ClassSP) = false
hasisomer(sc::DiChain) = nox(sc.acyl) > 0
hasisomer(::SumChain) = true
hasisomer(::LCB) = false
hasisomer(sc::ACYL) = nox(sc) > 0
hasisomer(::Nothing) = true
hasisomer(::Type{T}) where T = hasisomer(T())
"""
    hasnana

Determine whether the object has NeuAc.
"""
hasnana(::CLS.nana[0]) = false
hasnana(cls) = true

isomer_tuple(cls::ClassSP) = hasisomer(cls) ? cls.isomer : (cls, )

# instance api
"""
    cluster(analyte::AnalyteSP)

Return the cluster `analyte` belongs to. If none, it returns an empty array.
"""
function cluster(analyte::AnalyteSP)
    cls = deisomerized(class(analyte))
    project = last(analyte).project
    cluster_ = project.appendix[:cluster]
    in(cls, keys(cluster_)) && in(findfirst(==(analyte), project), cluster_[cls]) ? cluster_[cls] : Int[]
end
"""
    incluster(analyte::AnalyteSP)

Whether `analyte` belongs to a cluster or not.
"""
incluster(analyte::AnalyteSP) = !isempty(cluster(analyte))
"""
    class(analyte::AbstractAnalyteID)
    class(cpd::AbstractCompoundID)
    class(transition::TransitionID)

Return the class of an analyte or compound.
"""
class(analyte::AbstractAnalyteID) = class(last(analyte))
class(cpd::AbstractCompoundID) = cpd.class
class(transition::TransitionID) = class(transition.compound)
"""
    chain(analyte::AbstractAnalyteID)
    chain(cpd::AbstractCompoundID)
    chain(transition::TransitionID)

Return the chain of an analyte or compound.
"""
chain(analyte::AbstractAnalyteID) = chain(last(analyte))
chain(cpd::AbstractCompoundID) = cpd.chain
chain(transition::TransitionID) = chain(transition.compound)
"""
    sumcomp(analyte::AbstractAnalyteID)
    sumcomp(cpd::AbstractCompoundID)
    sumcomp(transition::TransitionID)
    sumcomp(sc::ChainSP)

Return sumcomposition information, i.e., a `SumChain` or `SumChainIS`.
"""
sumcomp(analyte::AbstractAnalyteID) = sumcomp(chain(analyte))
sumcomp(cpd::AbstractCompoundID) = sumcomp(chain(cpd))
sumcomp(transition::TransitionID) = sumcomp(transition.compound)
sumcomp(sc::SumChain) = sc
sumcomp(sc::SumChainIS) = sc
sumcomp(sc::ChainSP) = SumChain(ncb(sc), ndb(sc), nox(sc))
sumcomp(sc::ChainSPIS) = SumChainIS(ncb(sc), ndb(sc), nox(sc), (n13C = n13C(sc), nD = nD(sc)))
"""
    chainconstituent(analyte::AbstractAnalyteID)
    chainconstituent(cpd::AbstractCompoundID)
    chainconstituent(transition::TransitionID)
    chainconstituent(sc::ChainSP)

Return molecular-level chain composition. If there are only species-level information(sumcomposition), it returns `nothing`.
"""
chainconstituent(analyte::AbstractAnalyteID) = chainconstituent(chain(analyte))
chainconstituent(cpd::AbstractCompoundID) = chainconstituent(chain(cpd))
chainconstituent(transition::TransitionID) = chainconstituent(transition.compound)
chainconstituent(sc::ChainSP) = sc
chainconstituent(sc::SumChain) = nothing
chainconstituent(sc::SumChainIS) = nothing
chainconstituent(sc::ACYL) = nothing
chainconstituent(sc::ACYLIS) = nothing
"""
    isspceieslevel

Determine whether the object has only species-level information.
"""
const isspceieslevel = isnothing ∘ chainconstituent
"""
    lcb(analyte::AbstractAnalyteID)
    lcb(cpd::AbstractCompoundID)
    lcb(transition::TransitionID)
    lcb(sc::ChainSP)

Return long-chain base. If there is no long-chain base, it returns `nothing`.
"""
lcb(analyte::AbstractAnalyteID) = lcb(chain(analyte))
lcb(cpd::AbstractCompoundID) = lcb(chain(cpd))
lcb(transition::TransitionID) = lcb(transition.compound)
lcb(sc::ChainSP) = nothing
lcb(sc::DiChain) = sc.lcb
lcb(sc::LCB) = sc
"""
    acyl(analyte::AbstractAnalyteID)
    acyl(cpd::AbstractCompoundID)
    acyl(transition::TransitionID)
    acyl(sc::ChainSP)

Return N-acyl chain. If there is no N-acyl chain, it returns `nothing`.
"""
acyl(analyte::AbstractAnalyteID) = acyl(chain(analyte))
acyl(cpd::AbstractCompoundID) = acyl(chain(cpd))
acyl(transition::TransitionID) = acyl(transition.compound)
acyl(sc::ChainSP) = nothing
acyl(sc::DiChain) = sc.acyl
acyl(sc::ACYL) = sc
"""
    rt(analyte::AbstractAnalyteID)

Return retention time.
"""
rt(analyte::AbstractAnalyteID) = analyte.rt
rt(analyte) = analyte.rt
"""
    ncb(analyte::AbstractAnalyteID)
    ncb(cpd::AbstractCompoundID)
    ncb(transition::TransitionID)
    ncb(sc::ChainSP)
    ncb(ion::Ion)

Return number of carbons in chains. For ions, adducts are not taken into account.
"""
ncb(analyte::AbstractAnalyteID) = ncb(chain(analyte))
ncb(cpd::AbstractCompoundID) = ncb(chain(cpd))
ncb(transition::TransitionID) = ncb(transition.compound)
ncb(sc::ChainSP) = sc.carbon
ncb(sc::DiChain) = ncb(sc.lcb) + ncb(sc.acyl)
"""
    ndb(analyte::AbstractAnalyteID)
    ndb(cpd::AbstractCompoundID)
    ndb(transition::TransitionID)
    ndb(sc::ChainSP)
    ndb(ion::Ion)

Return number of double bonds in chains. For ions, adducts are not taken into account except in cases of neutral loss of H2O.
"""
ndb(analyte::AbstractAnalyteID) = ndb(chain(analyte))
ndb(cpd::AbstractCompoundID) = ndb(chain(cpd))
ndb(transition::TransitionID) = ndb(transition.compound)
ndb(sc::ChainSP) = sc.doublebond
ndb(sc::DiChain) = ndb(sc.lcb) + ndb(sc.acyl)
"""
    nox(analyte::AbstractAnalyteID)
    nox(cpd::AbstractCompoundID)
    nox(transition::TransitionID)
    nox(sc::ChainSP)
    nox(ion::Ion)

Return number of additional oxygens. For ions, adducts are not taken into account except in cases of neutral loss of H2O.
"""
nox(analyte::AbstractAnalyteID) = nox(chain(analyte))
nox(cpd::AbstractCompoundID) = nox(chain(cpd))
nox(transition::TransitionID) = nox(transition.compound)
nox(sc::ChainSP) = sc.oxygen
nox(sc::DiChain) = nox(sc.lcb) + nox(sc.acyl)
nox(sc::LCB) = ndbox(sc) - ndb(sc)
ndbox(analyte::AbstractAnalyteID) = ndbox(chain(analyte))
ndbox(cpd::AbstractCompoundID) = ndbox(chain(cpd))
ndbox(transition::TransitionID) = ndbox(transition.compound)
ndbox(sc::SumChain) = sc.doublebond + sc.oxygen
ndbox(sc::DiChain) = ndbox(sc.lcb) + ndbox(sc.acyl)
ndbox(::LCB2) = 2
ndbox(::LCB3) = 3
ndbox(::LCB4) = 4
ndbox(sc::ACYL) = ndb(sc) + nox(sc)

ncb(ion::Ion) = ncb(ion.molecule)
ndbox(ion::Ion) = ndbox(ion.molecule)
ndb(ion::Ion) = ndb(ion.molecule)
ndb(ion::Ion{Protonation}) = ndb(ion.molecule)
ndb(ion::Ion{ProtonationNLH2O}) = ndb(ion.molecule) + 1
ndb(ion::Ion{ProtonationNL2H2O}) = ndb(ion.molecule) + 2
ndb(ion::Ion{ProtonationNL3H2O}) = ndb(ion.molecule) + 3
nox(ion::Ion) = nox(ion.molecule)
nox(ion::Ion{Protonation}) = nox(ion.molecule)
nox(ion::Ion{ProtonationNLH2O}) = nox(ion.molecule) - 1
nox(ion::Ion{ProtonationNL2H2O}) = nox(ion.molecule) - 2
nox(ion::Ion{ProtonationNL3H2O}) = nox(ion.molecule) - 3

n13C(ion::Ion) = n13C(ion.molecule)
n13C(analyte::AbstractAnalyteID) = n13C(chain(analyte))
n13C(cpd::AbstractCompoundID) = n13C(chain(cpd))
n13C(transition::TransitionID) = n13C(transition.compound)
n13C(sumchain::SumChainIS) = sumchain.isotope.n13C
n13C(dichain::DiChainIS) = n13C(lcb(dichain)) + n13C(acyl(dichain))
n13C(lcb::LCBIS) = lcb.isotope.n13C
n13C(acyl::ACYLIS) = acyl.isotope.n13C
n13C(x) = 0

nD(ion::Ion) = nD(ion.molecule)
nD(analyte::AbstractAnalyteID) = nD(chain(analyte))
nD(cpd::AbstractCompoundID) = nD(chain(cpd))
nD(transition::TransitionID) = nD(transition.compound)
nD(sumchain::SumChainIS) = sumchain.isotope.nD
nD(dichain::DiChainIS) = nD(lcb(dichain)) + nD(acyl(dichain))
nD(lcb::LCBIS) = lcb.isotope.nD
nD(acyl::ACYLIS) = acyl.isotope.nD
nD(x) = 0
const n2H = nD

for fn in [:ncb, :ndb, :nox, :ndbox, :n13C, :nD]
    @eval $fn(::Nothing) = nothing
end

"""
    allow_unknown(f::Function)

Return a `FunctionalFunction` which returns true if the input `x` is `nothing` or `f(x)` is true.
"""
allow_unknown(f::T) where {T <: Function} = FunctionalFunction(allow_unknown, f, x -> isnothing(x) || f(x))
"""
    only_known(f::Function)

Return a `FunctionalFunction` which returns true if the input `x` is not `nothing` and `f(x)` is true.
"""
only_known(f::T) where {T <: Function} = FunctionalFunction(only_known, f, x -> !isnothing(x) && f(x))

# internal
is4e(sc::LCB) = nox(sc) < 3
is4e(sc::LCB2) = nox(sc) < 2

function default_adduct(sc::LCB, polarity::Bool)
    if polarity 
        is4e(sc) ? Ion(SPDB[:NLH2O][nox(sc) + 1], sc) : Ion(SPDB[:NLH2O][nox(sc)], sc)
    else
        throw(ErrorException("Not yet implement"))
    end
end

function hydroxyl_shift(adduct::Adduct, Δ::Int)
    id = findfirst(==(adduct), SPDB[:NLH2O]) + Δ
    0 < id <= length(SPDB[:NLH2O]) ? SPDB[:NLH2O][id] : nothing
end

hydroxyl_shift(sc::T, Δ::Int) where {T <: LCB} = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? T(ncb(sc), ndb(sc) - Δ) : nothing
hydroxyl_shift(sc::Union{NotPhyto3, NotPhyto3IS}, Δ::Int) = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? default_constructor(sc)(ncb(sc), ndb(sc) - Δ) : nothing
hydroxyl_shift(sc::Union{NotPhyto4, NotPhyto4IS}, Δ::Int) = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? default_constructor(sc)(ncb(sc), ndb(sc) - Δ) : nothing

default_constructor(::T) where T = T
default_constructor(lcb::T) where {T <: LCBIS} = (x...) -> T(x..., lcb.isotope)
default_constructor(acyl::T) where {T <: ACYLIS} = (x...) -> AcylIS(x..., acyl.isotope)
alpha_constructor(acyl::T) where {T <: ACYLIS} = (x...) -> AcylαIS(x..., acyl.isotope)
beta_constructor(acyl::T) where {T <: ACYLIS} = (x...) -> AcylβIS(x..., acyl.isotope)
default_constructor(::T) where {T <: ACYL} = (x...) -> Acyl(x...)
alpha_constructor(::T) where {T <: ACYL} = (x...) -> Acylα(x...)
beta_constructor(::T) where {T <: ACYL} = (x...) -> Acylβ(x...)
default_constructor(::NotPhyto3) = SPB3
notphyto_constructor(::LCB3) = NotPhyto3
default_constructor(::NotPhyto4) = SPB4
notphyto_constructor(::LCB4) = NotPhyto4
default_constructor(lcb::NotPhyto3IS) = (x...) -> SPB3IS(x..., lcb.isotope)
notphyto_constructor(lcb::LCBIS) = lcb isa LCB3 ? (x...) -> NotPhyto3IS(x..., lcb.isotope) : (x...) -> NotPhyto4IS(x..., lcb.isotope)
default_constructor(lcb::NotPhyto4IS) = (x...) -> SPB4IS(x..., lcb.isotope)
default_constructor(::Acylα) = Acyl
default_constructor(::Acylβ) = Acyl

# connection
find_connected(cpd::AbstractCompoundID, analyte::AbstractAnalyteID) = findall(id -> connected(cpd, id), analyte)

connected(cpd::AbstractCompoundID, id::AbstractCompoundID) = iscompatible(chain(cpd), chain(id)) && begin
    class1 = hasisomer(cpd.class) ? cpd.class.isomer : (cpd.class, )
    class2 = hasisomer(id.class) ? id.class.isomer : (id.class, )
    any(connected(cls1, cls2) || connected(cls2, cls1) for (cls1, cls2) in Iterators.product(class1, class2))
end
connected(cls1::ClassSP, cls2::ClassSP) = cls1 ≡ cls2 || any(connected(cls1, cls2) for cls1 in tuplize(SPDB[:CONNECTION][cls1]))
connected(cls1::ClassSP, cls2::Cer) = true
connected(cls1::Cer, cls2::ClassSP) = false
connected(cls1::Cer, cls2::Cer) = true

isf(cls::ClassSP) = hasisomer(cls) ? isf(cls.isomer) : push!(isf(SPDB[:CONNECTION][cls]), cls)
isf(cls::Tuple) = union!((push!(isf(SPDB[:CONNECTION][cls1]), cls1) for cls1 in cls)...)
isf(cls::Cer) = ClassSP[cls]

macro rule(expr)
    expr = pipe2rule(expr)
    return quote
        $(esc(expr))
    end
end

pipe2rule(expr) =
    @match expr begin
        Expr(:call, :(|>), arg1, arg2)  => Expr(:call, :(Rule), pipe2rule(arg1), pipe2rule(arg2))
        Expr(head, args...)             => Expr(head, map(pipe2rule, args)...)
        e                               => e
    end

# miscellaneous
vectorize(x::AbstractVector) = x
vectorize(x::UnitRange) = x
vectorize(x) = [x]

tuplize(x::Tuple) = x
tuplize(x) = (x,)

function merge_nt(old, new; tosum = (:n, ), toold = (:id, ), tonew = ())
    pn = propertynames(old)
    (;(p => _merge_nt(old, new, p; tosum, toold, tonew) for p in pn)...)
end

function _merge_nt(old, new, p; tosum = (:n, ), toold = (:id, ), tonew = ())
    if p in tosum
        +(getproperty(old, p), getproperty(new, p))
    elseif p in toold
        getproperty(old, p)
    elseif p in tonew
        getproperty(new, p)
    else
        (old.n * getproperty(old, p) + new.n * getproperty(new, p)) / (old.n + new.n)
    end
end

function mode_ion(ions)
    dict = Dict{Union{<: ClassSP, <: ChainSP}, Int64}()
    for is in ions
        for i in is
            dict[i] = get(dict, i, 0) + 1
        end
    end
    maxk, maxi = Union{<: ClassSP, <: ChainSP}[], 0
    for (k, v) in dict
        if maxi < v
            maxk = [k]
        elseif maxi ≡ v
            push!(maxk, k)
        end
        maxi = max(maxi, v)
    end
    maxk
end

function mode(v::Base.AbstractVecOrTuple{T}) where T
    X = try
        T
    catch _
        Any
    end
    dict = Dict{X, Int64}()
    for i in v
        dict[i] = get(dict, i, 0) + 1
    end
    maxk, maxi = X[], 0
    for (k, v) in dict
        if maxi < v
            maxk = X[k]
        elseif maxi ≡ v
            push!(maxk, k)
        end
        maxi = max(maxi, v)
    end
    maxk
end