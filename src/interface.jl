isempty(::ClassSP) = false

length(project::Project) = length(project.analyte)
length(aquery::AbstractQuery) = length(aquery.result)
length(analyte::AbstractAnalyteID) = length(analyte.compound)

iterate(project::Project) = iterate(project.analyte)
iterate(project::Project, i) = iterate(project.analyte, i)
iterate(aquery::AbstractQuery) = iterate(aquery.result)
iterate(aquery::AbstractQuery, i) = iterate(aquery.result, i)
iterate(analyte::AbstractAnalyteID) = iterate(analyte.compound)
iterate(analyte::AbstractAnalyteID, i) = iterate(analyte.compound, i)

getindex(project::Project, i) = getindex(project.analyte, i)
getindex(aquery::AbstractQuery, i) = getindex(aquery.result, i)
getindex(analyte::AbstractAnalyteID, i) = getindex(analyte.compound, i)

view(project::Project, i) = view(project.analyte, i)
view(aquery::AbstractQuery, i) = view(aquery.result, i)
view(analyte::AbstractAnalyteID, i) = view(analyte.compound, i)

firstindex(project::Project) = firstindex(project.analyte)
firstindex(aquery::AbstractQuery) = firstindex(aquery.result)
firstindex(analyte::AbstractAnalyteID) = firstindex(analyte.compound)

lastindex(project::Project) = lastindex(project.analyte)
lastindex(aquery::AbstractQuery) = lastindex(aquery.result)
lastindex(analyte::AbstractAnalyteID) = lastindex(analyte.compound)

sort(project::Project; kwargs...) = sort!(deepcopy(project); kwargs...)
sort!(project::Project; kwargs...) = sort!(project.analyte; kwargs...)
sort(aquery::AbstractQuery; kwargs...) = sort!(copy_wo_project(aquery); kwargs...)
sort!(aquery::AbstractQuery; kwargs...) = sort!(aquery.result; kwargs...)
sort(analyte::AbstractAnalyteID; kwargs...) = sort!(copy_wo_project(analyte); kwargs...)
sort!(analyte::AbstractAnalyteID; kwargs...) = sort!(analyte.compound; kwargs...)

push!(project::Project, analyte::AnalyteSP) = push!(project.analyte, analyte)
push!(analyte::AnalyteSP, cpd::CompoundSP) = push!(analyte.compound, cpd)
push!(analyte::AnalyteID, cpd::SPID) = push!(analyte.compound, cpd)

pop!(project::Project) = pop!(project.analyte)
pop!(analyte::AbstractAnalyteID) = pop!(analyte.compound)

popfirst!(project::Project) = popfirst!(project.analyte)
popfirst!(analyte::AbstractAnalyteID) = popfirst!(analyte.compound)

popat!(project::Project, del::Vector{Int}) = popat!(project.analyte, del)
popat!(analyte::AbstractAnalyteID, del::Vector{Int}) = popat!(analyte.compound, del)
popat!(project::Project, del::Int) = popat!(project.analyte, del)
popat!(analyte::AbstractAnalyteID, del::Int) = popat!(analyte.compound, del)

deleteat!(project::Project, del::Vector{Int}) = deleteat!(project.analyte, del)
deleteat!(analyte::AbstractAnalyteID, del::Vector{Int}) = deleteat!(analyte.compound, del)
deleteat!(project::Project, del::Int) = deleteat!(project.analyte, del)
deleteat!(analyte::AbstractAnalyteID, del::Int) = deleteat!(analyte.compound, del)

deleteat!(analyte::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = deleteat!(parent(analyte), parentindices(analyte)[1][del])
deleteat!(analyte::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Int) = deleteat!(parent(analyte), parentindices(analyte)[1][del])
popat!(analyte::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = popat!(parent(analyte), parentindices(analyte)[1][del])
popat!(analyte::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Int) = popat!(parent(analyte), parentindices(analyte)[1][del])

deleteat!(analyte::SubArray{AnalyteID, 1, Vector{AnalyteID}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = deleteat!(parent(analyte), parentindices(analyte)[1][del])
deleteat!(analyte::SubArray{AnalyteID, 1, Vector{AnalyteID}, Tuple{Vector{Int64}}, false}, del::Int) = deleteat!(parent(analyte), parentindices(analyte)[1][del])
popat!(analyte::SubArray{AnalyteID, 1, Vector{AnalyteID}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = popat!(parent(analyte), parentindices(analyte)[1][del])
popat!(analyte::SubArray{AnalyteID, 1, Vector{AnalyteID}, Tuple{Vector{Int64}}, false}, del::Int) = popat!(parent(analyte), parentindices(analyte)[1][del])

function delete!(aquery::AbstractQuery, target::Symbol)
    delete!(aquery.project, target; analyte = query.result)
    aquery.view ? (printstyled("Please re-query to get the correct result\n"; bold = true, color = :red); aquery.project) : aquery
end

function delete!(project::Project, target::Symbol; analyte = project.analyte)
    printstyled("Delete> ", uppercasefirst(string(target)), "\n"; color = :green, bold = true)
    deleteat!(analyte, findall(ana -> ana.state[state_id(target)] < 0, analyte))
    project
end

keys(project::Project) = LinearIndices(project.analyte)
keys(aquery::AbstractQuery) = LinearIndices(aquery.result)
keys(analyte::AbstractAnalyteID) = LinearIndices(analyte.compound)

reverse(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(deepcopy(project), start, stop)
reverse!(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(project.analyte, start, stop)
reverse(aquery::AbstractQuery, start::Int = 1, stop::Int = length(aquery)) = reverse!(copy_wo_project(aquery), start, stop)
reverse!(aquery::AbstractQuery, start::Int = 1, stop::Int = length(aquery)) = reverse!(aquery.result, start, stop)
reverse(analyte::AbstractAnalyteID, start::Int = 1, stop::Int = length(analyte)) = reverse!(copy_wo_project(analyte), start, stop)
reverse!(analyte::AbstractAnalyteID, start::Int = 1, stop::Int = length(analyte)) = reverse!(analyte.compound, start, stop)

in(x::Number, ri::EmptyInterval) = false
in(x::Number, ri::RealInterval) = between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator)
in(x::Number, ui::UnionInterval) = any(ri -> between(x; low = ri.lowerbound, up = ri.upperbound, lop = ri.leftoperator, rop = ri.rightoperator), ui.intervals)

function union(ris::Vararg{<: RealIntervals, N}) where N
    vri = RealInterval[]
    for ri in ris
        @match ri begin
            ::UnionInterval => union!(vri, ri.intervals)
            ::RealInterval  => push!(vri, ri)
            ::EmptyInterval => nothing
        end
    end
    unique!(vri)
    sort!(vri)
    nvri = RealInterval[]
    i = 1
    push!(nvri, vri[1])
    while i < length(vri)
        i += 1
        if nvri[end].upperbound in vri[i]
            nvri[end] = RealInterval(nvri[end].lowerbound, vri[i].upperbound, nvri[end].leftoperator, vri[i].rightoperator)
        elseif vri[i].lowerbound in nvri[end]
            if nvri[end].upperbound == vri[i].lowerbound
                nvri[end] = RealInterval(nvri[end].lowerbound, vri[i].upperbound, nvri[end].leftoperator, vri[i].rightoperator)
            end
        else
            push!(nvri, vri[i])
        end
    end
    length(nvri) == 1 ? nvri[1] : UnionInterval(Tuple(nvri))
end

intersect(ris::Vararg{<: Union{RealInterval, UnionInterval}, N}) where N = reduce(intersect, ris)
intersect(ri1::EmptyInterval, ri2::RealIntervals) = EmptyInterval()
intersect(ri1::RealIntervals, ri2::EmptyInterval) = EmptyInterval()
intersect(ri1::EmptyInterval, ri2::EmptyInterval) = EmptyInterval()

function intersect(ri1::RealInterval, ri2::RealInterval)
    ri1 === ri2 && return ri1
    if ri1 > ri2
        ri1, ri2 = ri2, ri1
    end
    if !xor((ri2.rightoperator)(ri1.upperbound, ri2.upperbound), (ri1.rightoperator)(ri2.upperbound, ri1.upperbound))
        return ri2
    elseif ri1.upperbound == ri2.upperbound
        return (ri2.rightoperator)(ri1.upperbound, ri2.upperbound) ? RealInterval(ri2.lowerbound, ri1.upperbound, ri2.leftoperator, ri1.rightoperator) : ri2
    elseif (ri1.rightoperator)(ri2.upperbound, ri1.upperbound)
        return ri2
    elseif (ri2.leftoperator)(ri2.lowerbound, ri1.upperbound) && (ri1.rightoperator)(ri2.lowerbound, ri1.upperbound)
        return RealInterval(ri2.lowerbound, ri1.upperbound, ri2.leftoperator, ri1.rightoperator)
    else 
        return EmptyInterval()
    end
end

intersect(ri1::UnionInterval, ri2::RealInterval) = union((intersect(ri, ri2) for ri in ri1)...)
intersect(ri1::RealInterval, ri2::UnionInterval) = union((intersect(ri1, ri) for ri in ri2)...)
intersect(ri1::UnionInterval, ri2::UnionInterval) = union((intersect(ria, rib) for ria in ri1, rib in ri2)...)
    
function union!(qs::Vararg{Query, N}) where N
    q = qs[1]
    length(qs) <= 1 && return q
    q.query = QueryOr([q.query])
    if q.view
        ids = parentindices(q.result)[1]
        for qo in qs[2:end]
            push!(q.query.qcmd, qo.query)
            union!(ids, parentindices(qo.result)[1])
        end
        q.result = @views parent(q.result)[ids]
    end
    q
end

function push_cpd!(analyte::AnalyteSP, cpd::CompoundSP)
    @p analyte push!(__, copy_wo_project(cpd)) sort!(; lt = isless_class)
    analyte.rt = calc_rt(analyte)
    analyte
end

union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP) = union!(analyte, id, cpd2, chain(cpd2))
function union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP, ::Union{SumChain, SumChainIS})
    cpd1 = analyte[id]
    append!(cpd1.fragment, cpd2.fragment)
    unique!(cpd1.fragment)
    _, i = findmax((cpd1.signal[1], cpd2.signal[1]))
    cpd1.signal = i ≡ 1 ? cpd1.signal : cpd2.signal
    analyte
end

function _union!(analyte::AnalyteSP, id::Int, ::Union{SumChain, SumChainIS}, sc::ChainSP)
    push!(last(analyte).project, copy_wo_project(analyte))
    for a in analyte
        a.chain = sc
    end
    analyte[id]
end

function _union!(analyte::AnalyteSP, id::Int, sc1::ChainSP, sc2::ChainSP)
    nox(lcb(sc1)) <= nox(lcb(sc2)) && (return analyte[id])
    for a in analyte
        a.chain = sc2
    end
    analyte[id]
end

function union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP, ::ChainSP)
    cpd1 = _union!(analyte, id, chain(analyte[id]), chain(cpd2))
    append!(cpd1.fragment, cpd2.fragment)
    unique!(cpd1.fragment)
    _, i = findmax((cpd1.signal[1], cpd2.signal[1]))
    cpd1.signal = i ≡ 1 ? cpd1.signal : cpd2.signal
    @p analyte sort!(; lt = isless_class)
    analyte.rt = calc_rt(analyte)
    analyte
end

union(cpd1::CompoundSP, cpd2::CompoundSP) = union!(copy_wo_project(cpd1), cpd2)
union(cpd1::CompoundSPVanilla, cpd2::CompoundSPVanilla) = union!(deepcopy(cpd1), cpd2)
function union!(cpd1::CompoundSPVanilla, cpd2::CompoundSPVanilla)
    append!(cpd1.fragment, cpd2.fragment)
    unique!(cpd1.fragment)
    _union!(cpd1, cpd2, cpd1.chain, cpd2.chain)
end

_union!(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID, ::Union{SumChain, SumChainIS}, ::Union{SumChain, SumChainIS}) = cpd1
_union!(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID, ::ChainSP, ::Union{SumChain, SumChainIS}) = cpd1
_union!(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID, ::Union{SumChain, SumChainIS}, sc::ChainSP) = (cpd1.chain = sc; cpd1)
_union!(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID, ::ChainSP, sc::ChainSP) =
    nox(lcb(cpd1)) > nox(lcb(cpd2)) ? (cpd1.chain = sc; cpd1) : cpd1

function union!(cpd1::CompoundSP, cpd2::CompoundSP)
    append!(cpd1.fragment, cpd2.fragment)
    unique!(cpd1.fragment)
    _, i = findmax((cpd1.signal[1], cpd2.signal[1]))
    cpd1.signal = i ≡ 1 ? cpd1.signal : cpd2.signal
    _union!(cpd1, cpd2, cpd1.chain, cpd2.chain)
end

union(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(copy_wo_project(analyte1), analyte2.compound, analyte2.state)
union(analyte1::AnalyteSP, cpds::Vector{CompoundSP}) = union!(copy_wo_project(analyte1), cpds)
union!(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(analyte1, analyte2.compound, analyte2.state)

function union!(analyte1::AnalyteSP, cpds::Vector{CompoundSP}, states2 = [0, 0, 0, 0, 0, 0])
    for cpd2 in cpds
        id = findfirst(cpd1 -> iscompatible(cpd1, cpd2), analyte1)
        isnothing(id) ? push!(analyte1, cpd2) : union!(analyte1[id], cpd2)
    end
    sort!(analyte1, lt = isless_class)
    analyte1.state = min.(analyte1.state, states2)
    analyte1
end

getproperty(reusable::ReusableQuery, sym::Symbol) = sym ≡ :query ? getfield(reusable, :query) : getfield(getfield(reusable, :query), sym)

convert(::Type{CompoundSP}, cpd::CompoundSP) = cpd
convert(::Type{CompoundSPVanilla}, cpd::CompoundSPVanilla) = cpd
convert(::Type{SPID}, cpd::SPID) = cpd
convert(::Type{CompoundSPVanilla}, cpd::CompoundSP) = CompoundSPVanilla(cpd.class, cpd.chain, cpd.fragment)
convert(::Type{CompoundSP}, cpd::CompoundSPVanilla) = CompoundSP(cpd.class, cpd.chain, cpd.fragment)
convert(::Type{SPID}, cpd::CompoundSP) = SPID(cpd.calss, cpd.chain)
convert(::Type{CompoundSP}, cpd::SPID) = CompoundSP(cpd.calss, cpd.chain)
convert(::Type{SPID}, cpd::CompoundSPVanilla) = SPID(cpd.calss, cpd.chain)
convert(::Type{CompoundSPVanilla}, cpd::SPID) = CompoundSPVanilla(cpd.calss, cpd.chain)
convert(::Type{AnalyteSP}, analyte::AnalyteSP) = analyte
convert(::Type{AnalyteID}, analyte::AnalyteID) = analyte
convert(::Type{TransitionID}, analyte::TransitionID) = analyte
convert(::Type{AnalyteSP}, analyte::AnalyteID) = AnalyteSP(analyte.compound, analyte.rt)
convert(::Type{AnalyteID}, analyte::AnalyteSP) = AnalyteID(analyte.compound, analyte.rt)

# variant of interface
# iscomponent: whether ion is a component of cpd
iscomponent(ion::Ion{<: Adduct, <: ClassSP}, cpd::AbstractCompoundID) = iscompatible(ion.molecule, cpd.class)
iscomponent(ion::Ion{<: Adduct, <: ChainSP}, cpd::AbstractCompoundID) = iscomponent(ion, cpd.chain)
iscomponent(ion::Ion{<: Adduct, <: ChainSP}, sc::SumChain) =
    ncb(ion) <= ncb(sc) && ndbox(ion) <= ndbox(sc) &&
        any(ox <= nox(sc) && db <= ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))
iscomponent(ion::Ion{<: Adduct, <: ChainSP}, sc::SumChainIS) =
    n13C(ion) <= n13C(sc) && nD(ion) <= nD(sc) && ncb(ion) <= ncb(sc) && ndbox(ion) <= ndbox(sc) &&
        any(ox <= nox(sc) && db <= ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))
    
iscomponent(ion::Ion{<: Adduct, <: LCB}, sc::DiChain) = iscomponent(ion, sc.lcb)
iscomponent(ion::Ion{<: Adduct, <: ACYL}, sc::DiChain) = iscomponent(ion, sc.acyl)
iscomponent(ion::Ion{<: Adduct, <: LCB}, sc::LCB) =
    ncb(ion) ≡ ncb(sc) && ndbox(ion) ≡ ndbox(sc) &&
        any(ox ≡ nox(sc) && db ≡ ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))
iscomponent(ion::Ion{<: Adduct, <: ACYL}, sc::ACYL) =
    ncb(ion) ≡ ncb(sc) && ndb(ion) ≡ ndb(sc) && nox(ion) ≡ ndb(sc)
iscomponent(ion::Ion{<: Adduct, <: LCB}, sc::LCBIS) =
    n13C(ion) ≡ n13C(sc) && nD(ion) ≡ nD(sc) && ncb(ion) ≡ ncb(sc) && ndbox(ion) ≡ ndbox(sc) &&
        any(ox ≡ nox(sc) && db ≡ ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))
iscomponent(ion::Ion{<: Adduct, <: ACYL}, sc::ACYLIS) =
    n13C(ion) ≡ n13C(sc) && nD(ion) ≡ nD(sc) && ncb(ion) ≡ ncb(sc) && ndb(ion) ≡ ndb(sc) && nox(ion) ≡ ndb(sc)

iscomponent(::Ion{<: Adduct, NeuAc}, cpd::AbstractCompoundID) = isa(cpd.class, CLS.fg.nana)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{NeuAc, NeuAc}}}, cpd::AbstractCompoundID) =
    isa(cpd.class, CLS.series.b) || isa(cpd.class, CLS.series.c) || isa(cpd.class, GT1a) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1c)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{NeuAc, NeuAc, NeuAc}}}, cpd::AbstractCompoundID) = isa(cpd.class, CLS.series.c)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, Hex}}}, cpd::AbstractCompoundID) =
    (isa(cpd.class, CLS.series.as) && !isa(cpd.class, Hex2Cer)) ||
    (isa(cpd.class, CLS.series.a) && !isa(cpd.class, GM3)) ||
    (isa(cpd.class, CLS.series.b) && !isa(cpd.class, GD3)) ||
    (isa(cpd.class, CLS.series.c) && !isa(cpd.class, GT3)) ||
    isa(cpd.class, HexNAcHex2Cer) ||
    isa(cpd.class, HexNAcHex3Cer)

iscomponent(::Ion{S, Glycan{Tuple{HexNAc, Hex, NeuAc}}}, cpd::AbstractCompoundID) where {S <: Adduct} =
    isa(cpd.class, CLS.fg.nana) && iscomponent(cpd, Ion(S(), Glycan(HexNAc(), Hex())))

iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, Hex, NeuAc, NeuAc}}}, cpd::AbstractCompoundID) =
    (isa(cpd.class, CLS.series.b) && !isa(cpd.class, GD3)) ||
    (isa(cpd.class, CLS.series.c) && !isa(cpd.class, GT3)) ||
    isa(cpd.class, GT1a) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1c) || isa(cpd.class, GD1α)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, NeuAc}}}, cpd::AbstractCompoundID) =
    isa(cpd.class, GP1cα) || isa(cpd.class, GQ1bα) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1α)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{Hex, NeuAc}}}, cpd::AbstractCompoundID) =
    isa(cpd.class, CLS.fg.nana)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{Hex, NeuAc, NeuAc}}}, cpd::AbstractCompoundID) =
    isa(cpd.class, CLS.series.b) ||
    isa(cpd.class, CLS.series.c) ||
    isa(cpd.class, GT1a) ||
    isa(cpd.class, GD1c)

# iscompatible: not exact, allow isomer
iscompatible(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID) =
    iscompatible(cpd1.class, cpd2.class) && iscompatible(cpd1.chain, cpd2.chain)

iscompatible(cls1::ClassSP, cls2::ClassSP) = cls1 ≡ cls2 ||
    any(cls1 ≡ cls2 for (cls1, cls2) in Iterators.product(isomer_tuple(cls1), isomer_tuple(cls2)))

iscompatible(x1, x2) = x1 ≡ x2
iscompatible(sc1::SumChain, sc2::SumChain) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::SumChain, sc2::ChainSP) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::ChainSP, sc2::SumChain) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::DiChain, sc2::DiChain) = iscompatible(sc1.lcb, sc2.lcb) && iscompatible(sc1.acyl, sc2.acyl)
iscompatible(sc1::LCB, sc2::LCB) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::ACYL, sc2::ACYL) = sc1 ≡ sc2
iscompatible(sc1::Acyl, sc2::Acyl) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::Acyl, sc2::ACYL) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::ACYL, sc2::Acyl) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(ion::Ion, cpd::AbstractCompoundID) = iscomponent(ion, cpd)

copy_wo_project(cpd::CompoundSP) = CompoundSP(cpd.class, cpd.chain, deepcopy(cpd.fragment), cpd.signal, deepcopy(cpd.state), deepcopy(cpd.result), cpd.project)
copy_wo_project(cpd::SPID) = SPID(cpd.class, cpd.chain)
copy_wo_project(analyte::AnalyteSP) = AnalyteSP(copy_wo_project.(analyte.compound), analyte.rt, deepcopy(analyte.state), analyte.cpdsc, analyte.score)
copy_wo_project(analyte::AnalyteID) = AnalyteID(copy_wo_project.(analyte.compound), analyte.rt)
copy_wo_project(project::Project) = Project(copy_wo_project.(project.analyte), deepcopy(project.data), deepcopy(project.quantification), deepcopy(project.appendix))
copy_wo_project(aquery::Query) = Query(aquery.project, copy_wo_project.(aquery.result), deepcopy(aquery.query), false)
reuse_copy(aquery::Query) = Query(aquery.project, reuse_copy(aquery.result), deepcopy(aquery.query), true)
reuse_copy(v::Vector) = copy(v)
reuse_copy(v::T) where {T <: SubArray} = T(parent(v), v.indices, v.offset1, v.stride1)

# equivalent: ion1 and ion2 are the same
equivalent_in(ion, collection) = any(equivalent(ion, x) for x in collection)
equivalent(ion1::Ion{<: Pos, <: LCB}, ion2::Ion{<: Pos, <: LCB}) =
    ncb(ion1) ≡ ncb(ion2) && ndb(ion1) ≡ ndb(ion2) && nox(ion1) ≡ nox(ion2)
    #any(ox1 ≡ ox2 && db1 ≡ db2 for ((ox1, db1), (ox2, db2)) in Iterators.product(zip(nox(ion1):nox(ion1.molecule), ndb(ion1):-1:ndb(ion1.molecule)), zip(nox(ion2):nox(ion2.molecule), ndb(ion2):-1:ndb(ion2.molecule))))

equivalent(ion1::AbstractIon{S}, ion2::AbstractIon{S}) where {S <: Adduct} = iscompatible(ion1.molecule, ion2.molecule)
equivalent(ion1::AbstractIon, ion2::AbstractIon) = false

function equivalent_in_ion2(cpd::AbstractCompoundSP, criteria, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    equivalent_in(criteria, filterview(eqin, cpd.fragment).ion2)
end

function equivalent_in_ion2(cpd::AbstractCompoundSP, criteria::Tuple, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    any(equivalent_in(frag, criteria) for frag in filterview(eqin, cpd.fragment).ion2)
end

equivalent_in_ion2(analyte::AnalyteSP, criteria, prec) = any(equivalent_in_ion2(cpd, criteria, prec) for cpd in analyte)

equivalent_in_ion1(cpd::AbstractCompoundSP, criteria) = equivalent_in(criteria, cpd.fragment.ion1)
equivalent_in_ion1(cpd::AbstractCompoundSP, criteria::Tuple) = any(equivalent_in(frag, criteria) for frag in cpd.fragment.ion1)
equivalent_in_ion1(analyte::AnalyteSP, criteria) = any(equivalent_in_ion1(cpd, criteria) for cpd in analyte)

# order
isless_class(cpd1::AbstractCompoundID) = true
isless_class(cpd1::AbstractCompoundID, cpd2::AbstractCompoundID) =
    any(connected(cls2, cls1) for (cls1, cls2) in Iterators.product(isomer_tuple(cpd1.class), isomer_tuple(cpd2.class)))

isless_ion(ion1) = true
function isless_ion(ion1::Ion{<: Adduct, <: ClassSP}, ion2::Ion{<: Adduct, <: ClassSP})
    id = class_db_index(ion1.molecule)
    level = (:default_ion, :parent_ion, :adduct_ion)
    id1 = findfirst(p -> in(ion1.adduct, getproperty(id, p)), level)
    id2 = findfirst(p -> in(ion2.adduct, getproperty(id, p)), level)
    id1 > id2 && return true
    id1 < id2 && return false
    id1 = findfirst(==(ion1.adduct), getproperty(id, level[id1]))
    id1 = findfirst(==(ion2.adduct), getproperty(id, level[id2]))
    id1 > id2 && return true
    false
end

isless_ion(a, b) = isless(a, b)

isless_nan_min(a, b) = isnan(a) || !isnan(b) && isless(a, b)
isless_nan_max(a, b) = isnan(b) || !isnan(a) && isless(a, b)

function isless(x::RealInterval, y::RealInterval)
    x.lowerbound < y.lowerbound && return true
    x.lowerbound > y.lowerbound && return false
    xor((x.leftoperator)(x.lowerbound, y.lowerbound), (y.leftoperator)(x.lowerbound, y.lowerbound)) && return (x.leftoperator)(x.lowerbound, y.lowerbound)
    x.upperbound < y.upperbound && return true
    x.upperbound > y.upperbound && return false
    xor((x.rightoperator)(x.upperbound, y.upperbound), (y.rightoperator)(x.upperbound, y.upperbound)) && return (y.rightoperator)(x.upperbound, y.upperbound)
    false
end

isequal(x::RealInterval, y::RealInterval) = 
    x.lowerbound == y.lowerbound && x.upperbound == y.upperbound && !xor((x.rightoperator)(x.upperbound, y.upperbound), (y.rightoperator)(x.upperbound, y.upperbound))

sort(tbl::Table, i::Symbol; kwargs...) = tbl[sortperm(getproperty(tbl, i); kwargs...)]
function sort!(tbl::Table, i::Symbol; kwargs...)
    tbl[:] = tbl[sortperm(getproperty(tbl, i); kwargs...)]
end
sort(tbl::Table, v::Vector{Symbol}; kwargs...) = tbl[sortperm(collect(zip(getproperty.(Ref(tbl), v)...)); kwargs...)]
sort!(tbl::Table, v::Vector{Symbol}; kwargs...) = (tbl[:] .= tbl[sortperm(collect(zip(getproperty.(Ref(tbl), v)...)); kwargs...)])
sort(tbl::Table, i::AbstractString; kwargs...) = sort(tbl, Symbol(i); kwargs...)
sort!(tbl::Table, i::AbstractString; kwargs...) = sort!(tbl, Symbol(i); kwargs...)
sort(tbl::Table, v::Vector{<: AbstractString}; kwargs...) = sort(tbl, Symbol.(v); kwargs...)
sort!(tbl::Table, v::Vector{<: AbstractString}; kwargs...) = sort!(tbl, Symbol.(v); kwargs...)
sort(tbl::Table, i::Int; kwargs...) = sort(tbl, propertynames(tbl)[i]; kwargs...)
sort!(tbl::Table, i::Int; kwargs...) = sort!(tbl, propertynames(tbl)[i]; kwargs...)
sort(tbl::Table, v::Vector{<: Int}; kwargs...) = sort(tbl, propertynames(tbl)[v]; kwargs...)
sort!(tbl::Table, v::Vector{<: Int}; kwargs...) = sort!(tbl, propertynames(tbl)[v]; kwargs...)