# public
nfrags(cpd::AbstractCompoundSP) = size(cpd.fragments, 1)
function nMRM(tbl::Table, precision::Float64 = 0.1; end_time = maximum(tbl.rt .+ tbl.Δrt))
    trans = Dict(0:precision:end_time .=> 0)
    for id in eachindex(tbl)
        for k in keys(trans)
            if k <= tbl.rt[id] + tbl.Δrt[id] / 2 && k >= tbl.rt[id] - tbl.Δrt[id] / 2
                trans[k] += 1
            end
        end
    end
    sort!(collect(trans), by = (x -> x.second))
end

new_project(aquery::AbstractQuery) = new_project(aquery.project; analytes = aquery.view ? copy_wo_project.(aquery.result) : aquery.result)
function new_project(project::Project; analytes = copy_wo_project.(project.analytes))
    project_new = Project(analytes, project.data, project.anion)
    for analyte in project_new
        for cpd in analyte
            cpd.project = project_new
        end
    end
    project_new
end

# user interface for query
qand(qcmds...) = QueryAnd(collect(map(qcmd, qcmds)))
qor(qcmds...) = QueryOr(collect(map(qcmd, qcmds)))
qnot(qcmd::QueryCommands) = QueryNot(qcmd)
qnot(qcmd::QueryNot) = qcmd.qcmd
qnot(qcmd) = QueryNot(QueryCmd(qcmd))

qcmd(qcmd::QueryCommands) = qcmd
qcmd(qcmd) = QueryCmd(qcmd)
qcmd(query, neg) = neg ? qnot(query) : qcmd(query)

function lcb(c::Int, d::Int, n::Int)
    Cons = @match d + n begin
        2 => SPB2
        3 => SPB3
        4 => SPB4
    end
    Cons(c, d)
end
acyl(c::Int, d::Int, n::Int) = Acyl(c, d, n)
acylα(c::Int, d::Int, n::Int) = Acylα(c, d, n)
acylβ(c::Int, d::Int, n::Int) = Acylβ(c, d, n)
chain(spb::LCB, nacyl::ACYL) = DiChain(spb, nacyl)
chain(c::Int, n::Int, o::Int) = SumChain(c, n, o)
spid(cls::Type{<: ClassSP}, sc::ChainSP) = SPID(cls(), sc)
spid(cls::Type{<: ClassSP}, spb::LCB, nacyl::ACYL) = SPID(cls(), DiChain(spb, nacyl))
spid(cls::Type{<: ClassSP}, sum::NTuple{3, Int}) = SPID(cls(), SumChain(sum...))
spid(cls::Type{<: ClassSP}, c::Int, n::Int, o::Int) = SPID(cls(), SumChain(c, n, o))
#=
convert_type(::Type{LCB}, lcb::Lcb) = @match (lcb.db + lcb.ox, lcb.ox) begin
    (2, N) => LCB2{N, lcb.cb}
    (3, N) => LCB3{N, lcb.cb}
    (4, N) => LCB4{N, lcb.cb}
end
convert_type(::Type{LCB}, lcb) = Nothing

convert_type(::Type{ACYL}, acyl::Nacylα) = Acylα{acyl.ox}
convert_type(::Type{ACYL}, acyl::Nacylβ) = Acylβ{acyl.ox}
convert_type(::Type{ACYL}, acyl::Nacyl) = Acyl{acyl.ox}
convert_type(::Type{ACYL}, acyl) = Nothing
convert_internal(cpd::CompoundID{C}) where C = (C, cpd.sum, convert_type(LCB, cpd.lcb), convert_type(ACYL, cpd.acyl))
=#
# isomer
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

hasnana(::CLS.nana[0]) = false
hasnana(cls) = true

isomer_tuple(cls::ClassSP) = hasisomer(cls) ? cls.isomer : (cls, )

# instance api
class(analyte::AnalyteSP) = class(last(analyte))
class(cpd::CompoundID) = cpd.class
chain(analyte::AnalyteSP) = chain(last(analyte))
chain(cpd::CompoundID) = cpd.chain

sumcomp(analyte::AnalyteSP) = sumcomp(chain(analyte))
sumcomp(cpd::CompoundID) = sumcomp(chain(cpd))
sumcomp(sc::SumChain) = sc
sumcomp(sc::ChainSP) = SumChain(ncb(sc), ndb(sc), nox(sc))
chainconstituent(analyte::AnalyteSP) = chainconstituent(chain(analyte))
chainconstituent(cpd::CompoundID) = chainconstituent(chain(cpd))
chainconstituent(sc::ChainSP) = sc
chainconstituent(sc::SumChain) = nothing
chainconstituent(sc::ACYL) = nothing
const isspceieslevel = isnothing ∘ chainconstituent

lcb(analyte::AnalyteSP) = lcb(chain(analyte))
lcb(cpd::CompoundID) = lcb(chain(cpd))
lcb(sc::ChainSP) = nothing
lcb(sc::DiChain) = sc.lcb
lcb(sc::LCB) = sc
acyl(analyte::AnalyteSP) = acyl(chain(analyte))
acyl(cpd::CompoundID) = acyl(chain(cpd))
acyl(sc::ChainSP) = nothing
acyl(sc::DiChain) = sc.acyl
acyl(sc::ACYL) = sc
rt(analyte::AnalyteSP) = analyte.rt

ndbox(analyte::AnalyteSP) = ndbox(chain(analyte))
ndbox(cpd::CompoundID) = ndbox(chain(cpd))
ndbox(sc::SumChain) = sc.doublebond + sc.oxygen
ndbox(sc::DiChain) = ndbox(sc.lcb) + ndbox(sc.acyl)
ndbox(::LCB2) = 2
ndbox(::LCB3) = 3
ndbox(::LCB4) = 4
ndbox(sc::ACYL) = ndb(sc) + nox(sc)

ncb(analyte::AnalyteSP) = ncb(chain(analyte))
ncb(cpd::CompoundID) = ncb(chain(cpd))
ncb(sc::ChainSP) = sc.carbon
ncb(sc::DiChain) = ncb(sc.lcb) + ncb(sc.acyl)
nox(analyte::AnalyteSP) = nox(chain(analyte))
nox(cpd::CompoundID) = nox(chain(cpd))
nox(sc::ChainSP) = sc.oxygen
nox(sc::DiChain) = nox(sc.lcb) + nox(sc.acyl)
nox(sc::LCB) = ndbox(sc) - ndb(sc)
ndb(analyte::AnalyteSP) = ndb(chain(analyte))
ndb(cpd::CompoundID) = ndb(chain(cpd))
ndb(sc::ChainSP) = sc.doublebond
ndb(sc::DiChain) = ndb(sc.lcb) + ndb(sc.acyl)

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

between(num::Number, range) = first(range) <= num <= last(range)
between(num::Number; up, low) = low <= num <= up # REVERSE~~~~
between(num::Number, value, tol) = value - tol <= num <= value + tol
intersection(range...) = (maximum(first(r) for r in range), minimum(last(r) for r in range))
intersection(range::AbstractVector) = (maximum(first(r) for r in range), minimum(last(r) for r in range))

function equivalent_in_ion2(cpd::AbstractCompoundSP, criteria, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    equivalent_in(criteria, filterview(eqin, cpd.fragments).ion2)
end

function equivalent_in_ion2(cpd::AbstractCompoundSP, criteria::Tuple, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    any(equivalent_in(frag, criteria) for frag in filterview(eqin, cpd.fragments).ion2)
end

equivalent_in_ion2(analyte::AnalyteSP, criteria, prec) = any(equivalent_in_ion2(cpd, criteria, prec) for cpd in analyte)

equivalent_in_ion1(cpd::AbstractCompoundSP, criteria) = equivalent_in(criteria, cpd.fragments.ion1)
equivalent_in_ion1(cpd::AbstractCompoundSP, criteria::Tuple) = any(equivalent_in(frag, criteria) for frag in cpd.fragments.ion1)
equivalent_in_ion1(analyte::AnalyteSP, criteria) = any(equivalent_in_ion1(cpd, criteria) for cpd in analyte)

calc_rt(analyte::AnalyteSP) = mean(mean(query_raw(cpd.project, cpd.fragments.source[id], cpd.fragments.id[id]).rt for id in eachindex(cpd.fragments)) for cpd in analyte)
#=
lcb(analyte::AnalyteSP) = _lcb(last(analyte))
lcb(cpd::CompoundSP) = _lcb(cpd.chain)
lcb(chain::Chain) = chain.lcb
acyl(analyte::AnalyteSP) = acyl(last(analyte))
acyl(::Nothing) = nothing
function acyl(cpd::CompoundSP)
    spb_c, spb_db, spb_o = sumcomp(cpd.chain.lcb)
    string("Acyl ", cpd.sum[1] - spb_c, ":", cpd.sum[2] - spb_db, repr_hydroxl(cpd.chain.acyl))
end
chain(analyte::AnalyteSP) = chain(last(analyte))
chain(::Nothing) = nothing
function chain(cpd::CompoundSP)
    spb_c, spb_db, spb_o = sumcomp(cpd.chain.lcb)
    string(spb_c, ":", spb_db, ";", "O", spb_o > 1 ? spb_o : "", "/", cpd.sum[1] - spb_c, ":", cpd.sum[2] - spb_db, repr_hydroxl(cpd.chain.acyl))
end
sumcomps(analyte::AnalyteSP) = sumcomps(last(analyte))
sumcomps(cpd::CompoundSP) = string(cpd.sum[1], ":", cpd.sum[2], ";", "O", cpd.sum[3] > 1 ? cpd.sum[3] : "")
=#
states_id = @λ begin
    :class  => 1
    :chain  => 2
    :rt     => 3
    :error  => 4
    :isf    => 5
end