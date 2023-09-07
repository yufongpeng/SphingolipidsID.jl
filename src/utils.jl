# public
"""
    nfrags(cpd::AbstractCompoundSP)

Number of fragments of a compund, i.e., `size(cpd.fragments, 1)`.
"""
nfrags(cpd::AbstractCompoundSP) = size(cpd.fragments, 1)
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
    new_project(project::Project; copy_data = false, analytes = copy_wo_project.(project.analytes))

Create a new project from a query or project with given analytes. `copy_data` determines whether copying the raw data or not.
"""
new_project(aquery::AbstractQuery; copy_data = false) = new_project(aquery.project; copy_data, analytes = aquery.view ? copy_wo_project.(aquery.result) : aquery.result)
function new_project(project::Project; copy_data = false, analytes = copy_wo_project.(project.analytes))
    project_new = Project(analytes, copy_data ? deepcopy(project.data) : project.data, project.anion, empty(project.clusters), empty(project.appendix))
    for analyte in project_new
        for cpd in analyte
            cpd.project = project_new
        end
    end
    project_new
end
"""
    between(num::Number, range)
    between(num::Number; low, up)
    between(num::Number, value, tol)

Determine whether `num` lies in given range. 

The lower bound is `first(range)`, `low`, or `value - tol`; the upper bound is `last(range)`, `up`, or `value + tol`. 

It is inclusive on both sides.
"""
between(num::Number, range) = first(range) <= num <= last(range)
between(num::Number; low, up) = low <= num <= up
between(num::Number, value, tol) = value - tol <= num <= value + tol
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
calc_rt(analyte::AnalyteSP) = mean(mean(query_raw(cpd.project, fragments.source, fragments.id).rt for fragments in cpd.fragments) for cpd in analyte)
"""
    states_id

Return the index for certain evaluation parameters in the states of an analyte.
"""
states_id = @λ begin
    :class  => 1
    :chain  => 2
    :rt     => 3
    :error  => 4
    :isf    => 5
    :total  => 6
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
    lcb(c::Int, d::Int, n::Int)

Create `LCB` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
function lcb(c::Int, d::Int, n::Int)
    Cons = @match d + n begin
        2 => SPB2
        3 => SPB3
        4 => SPB4
    end
    Cons(c, d)
end
"""
    acyl(c::Int, d::Int, n::Int)

Create `Acyl` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acyl(c::Int, d::Int, n::Int) = Acyl(c, d, n)
"""
    acylα(c::Int, d::Int, n::Int)

Create `Acylα` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acylα(c::Int, d::Int, n::Int) = Acylα(c, d, n)
"""
    acylβ(c::Int, d::Int, n::Int)

Create `Acylβ` with given carbons(`c`), double bonds(`d`), and addtional oxygens(`n`).
"""
acylβ(c::Int, d::Int, n::Int) = Acylβ(c, d, n)
"""
    chain(spb::LCB, nacyl::ACYL)
    chain(c::Int, n::Int, o::Int)

Create `ChainSP` with long-chain base and N-acyl or carbons, double bonds, and addtional oxygens.
"""
chain(spb::LCB, nacyl::ACYL) = DiChain(spb, nacyl)
chain(c::Int, n::Int, o::Int) = SumChain(c, n, o)
"""
    spid(cls::Type{<: ClassSP}, sc::ChainSP)
    spid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL)
    spid(cls::Type{<: ClassSP}, sum::NTuple{3, Int})
    spid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int)

Create `SPID` with given class and chain information.
"""
spid(cls::Type{<: ClassSP}, sc::ChainSP) = SPID(cls(), sc)
spid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL) = SPID(cls(), DiChain(spb, nacyl))
spid(cls::Type{<: ClassSP}, sum::NTuple{3, Int}) = SPID(cls(), SumChain(sum...))
spid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int) = SPID(cls(), SumChain(c, n, o))

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
    class(analyte::AnalyteSP)
    class(cpd::CompoundID)

Return the class of an analyte or compound.
"""
class(analyte::AnalyteSP) = class(last(analyte))
class(cpd::CompoundID) = cpd.class
"""
    chain(analyte::AnalyteSP)
    chain(cpd::CompoundID)

Return the chain of an analyte or compound.
"""
chain(analyte::AnalyteSP) = chain(last(analyte))
chain(cpd::CompoundID) = cpd.chain
"""
    sumcomp(analyte::AnalyteSP)
    sumcomp(cpd::CompoundID)
    sumcomp(sc::ChainSP)

Return sumcomposition information, i.e., a `SumChain`.
"""
sumcomp(analyte::AnalyteSP) = sumcomp(chain(analyte))
sumcomp(cpd::CompoundID) = sumcomp(chain(cpd))
sumcomp(sc::SumChain) = sc
sumcomp(sc::ChainSP) = SumChain(ncb(sc), ndb(sc), nox(sc))
"""
    chainconstituent(analyte::AnalyteSP)
    chainconstituent(cpd::CompoundID)
    chainconstituent(sc::ChainSP)

Return molecular-level chain composition. If there are only species-level information(sumcomposition), it returns `nothing`.
"""
chainconstituent(analyte::AnalyteSP) = chainconstituent(chain(analyte))
chainconstituent(cpd::CompoundID) = chainconstituent(chain(cpd))
chainconstituent(sc::ChainSP) = sc
chainconstituent(sc::SumChain) = nothing
chainconstituent(sc::ACYL) = nothing
"""
    isspceieslevel

Determine whether the object has only species-level information.
"""
const isspceieslevel = isnothing ∘ chainconstituent
"""
    lcb(analyte::AnalyteSP)
    lcb(cpd::CompoundID)
    lcb(sc::ChainSP)

Return long-chain base. If there is no long-chain base, it returns `nothing`.
"""
lcb(analyte::AnalyteSP) = lcb(chain(analyte))
lcb(cpd::CompoundID) = lcb(chain(cpd))
lcb(sc::ChainSP) = nothing
lcb(sc::DiChain) = sc.lcb
lcb(sc::LCB) = sc
"""
    acyl(analyte::AnalyteSP)
    acyl(cpd::CompoundID)
    acyl(sc::ChainSP)

Return N-acyl chain. If there is no N-acyl chain, it returns `nothing`.
"""
acyl(analyte::AnalyteSP) = acyl(chain(analyte))
acyl(cpd::CompoundID) = acyl(chain(cpd))
acyl(sc::ChainSP) = nothing
acyl(sc::DiChain) = sc.acyl
acyl(sc::ACYL) = sc
"""
    rt(analyte::AnalyteSP)

Return retention time.
"""
rt(analyte::AnalyteSP) = analyte.rt
"""
    ncb(analyte::AnalyteSP)
    ncb(cpd::CompoundID)
    ncb(sc::ChainSP)
    ncb(ion::Ion)

Return number of carbons. For ions, adducts are not taken into account.
"""
ncb(analyte::AnalyteSP) = ncb(chain(analyte))
ncb(cpd::CompoundID) = ncb(chain(cpd))
ncb(sc::ChainSP) = sc.carbon
ncb(sc::DiChain) = ncb(sc.lcb) + ncb(sc.acyl)
"""
    ndb(analyte::AnalyteSP)
    ndb(cpd::CompoundID)
    ndb(sc::ChainSP)
    ndb(ion::Ion)

Return number of double bonds. For ions, adducts are not taken into account except in cases of neutral loss of H2O.
"""
ndb(analyte::AnalyteSP) = ndb(chain(analyte))
ndb(cpd::CompoundID) = ndb(chain(cpd))
ndb(sc::ChainSP) = sc.doublebond
ndb(sc::DiChain) = ndb(sc.lcb) + ndb(sc.acyl)
"""
    nox(analyte::AnalyteSP)
    nox(cpd::CompoundID)
    nox(sc::ChainSP)
    nox(ion::Ion)

Return number of additional oxygens. For ions, adducts are not taken into account except in cases of neutral loss of H2O.
"""
nox(analyte::AnalyteSP) = nox(chain(analyte))
nox(cpd::CompoundID) = nox(chain(cpd))
nox(sc::ChainSP) = sc.oxygen
nox(sc::DiChain) = nox(sc.lcb) + nox(sc.acyl)
nox(sc::LCB) = ndbox(sc) - ndb(sc)
ndbox(analyte::AnalyteSP) = ndbox(chain(analyte))
ndbox(cpd::CompoundID) = ndbox(chain(cpd))
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

# internal
is4e(sc::LCB) = nox(sc) < 3
is4e(sc::LCB2) = nox(sc) < 2

default_adduct(sc::LCB) =
    is4e(sc) ? Ion(SPDB[:NLH2O][nox(sc) + 1], sc) : Ion(SPDB[:NLH2O][nox(sc)], sc)

function hydroxyl_shift(adduct::Adduct, Δ::Int)
    id = findfirst(==(adduct), SPDB[:NLH2O]) + Δ
    0 < id <= length(SPDB[:NLH2O]) ? SPDB[:NLH2O][id] : nothing
end

hydroxyl_shift(sc::T, Δ::Int) where {T <: LCB} = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? T(ncb(sc), ndb(sc) - Δ) : nothing
hydroxyl_shift(sc::NotPhyto3, Δ::Int) = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? SPB3(ncb(sc), ndb(sc) - Δ) : nothing
hydroxyl_shift(sc::NotPhyto4, Δ::Int) = (0 <= ndb(sc) - Δ <= ndbox(sc)) ? SPB4(ncb(sc), ndb(sc) - Δ) : nothing

# connection
find_connected(cpd::CompoundID, analyte::AnalyteSP) = findall(id -> connected(cpd, id), analyte)

connected(cpd::CompoundID, id::CompoundID) = iscompatible(chain(cpd), chain(id)) && begin
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
        (old.n * getproperty(old, p) + new.n) / (old.n + new.n)
    end
end

function mode(ions)
    dict = Dict{Union{ClassSP, ChainSP}}()
    for is in ions
        for i in is
            dict[i] = get(dict, i, 0) + 1
        end
    end
    maxk, maxi = Union{ClassSP, ChainSP}[], 0
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