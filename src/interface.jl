isempty(::ClassSP) = false

length(project::Project) = length(project.analytes)
length(aquery::AbstractQuery) = length(aquery.result)
length(analyte::AnalyteSP) = length(analyte.compounds)

iterate(project::Project) = iterate(project.analytes)
iterate(project::Project, i) = iterate(project.analytes, i)
iterate(aquery::AbstractQuery) = iterate(aquery.result)
iterate(aquery::AbstractQuery, i) = iterate(aquery.result, i)
iterate(analyte::AnalyteSP) = iterate(analyte.compounds)
iterate(analyte::AnalyteSP, i) = iterate(analyte.compounds, i)

getindex(project::Project, i) = getindex(project.analytes, i)
getindex(aquery::AbstractQuery, i) = getindex(aquery.result, i)
getindex(analyte::AnalyteSP, i) = getindex(analyte.compounds, i)

view(project::Project, i) = view(project.analytes, i)
view(aquery::AbstractQuery, i) = view(aquery.result, i)
view(analyte::AnalyteSP, i) = view(analyte.compounds, i)

firstindex(project::Project) = firstindex(project.analytes)
firstindex(aquery::AbstractQuery) = firstindex(aquery.result)
firstindex(analyte::AnalyteSP) = firstindex(analyte.compounds)

lastindex(project::Project) = lastindex(project.analytes)
lastindex(aquery::AbstractQuery) = lastindex(aquery.result)
lastindex(analyte::AnalyteSP) = lastindex(analyte.compounds)

sort(project::Project; kwargs...) = sort!(deepcopy(project); kwargs...)
sort!(project::Project; kwargs...) = sort!(project.analytes; kwargs...)
sort(aquery::AbstractQuery; kwargs...) = sort!(copy_wo_project(aquery); kwargs...)
sort!(aquery::AbstractQuery; kwargs...) = sort!(aquery.result; kwargs...)
sort(analyte::AnalyteSP; kwargs...) = sort!(copy_wo_project(analyte); kwargs...)
sort!(analyte::AnalyteSP; kwargs...) = sort!(analyte.compounds; kwargs...)

push!(project::Project, analyte::AnalyteSP) = push!(project.analytes, analyte)
push!(analyte::AnalyteSP, cpd::CompoundSP) = push!(analyte.compounds, cpd)

pop!(project::Project) = pop!(project.analytes)
pop!(analyte::AnalyteSP) = pop!(analyte.compounds)

popfirst!(project::Project) = popfirst!(project.analytes)
popfirst!(analyte::AnalyteSP) = popfirst!(analyte.compounds)

popat!(project::Project, del::Vector{Int}) = popat!(project.analytes, del)
popat!(analyte::AnalyteSP, del::Vector{Int}) = popat!(analyte.compounds, del)
popat!(project::Project, del::Int) = popat!(project.analytes, del)
popat!(analyte::AnalyteSP, del::Int) = popat!(analyte.compounds, del)

deleteat!(project::Project, del::Vector{Int}) = deleteat!(project.analytes, del)
deleteat!(analyte::AnalyteSP, del::Vector{Int}) = deleteat!(analyte.compounds, del)
deleteat!(project::Project, del::Int) = deleteat!(project.analytes, del)
deleteat!(analyte::AnalyteSP, del::Int) = deleteat!(analyte.compounds, del)

deleteat!(analytes::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = deleteat!(parent(analytes), parentindices(analytes)[1][del])
deleteat!(analytes::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Int) = deleteat!(parent(analytes), parentindices(analytes)[1][del])
popat!(analytes::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = popat!(parent(analytes), parentindices(analytes)[1][del])
popat!(analytes::SubArray{AnalyteSP, 1, Vector{AnalyteSP}, Tuple{Vector{Int64}}, false}, del::Int) = popat!(parent(analytes), parentindices(analytes)[1][del])

function delete!(aquery::AbstractQuery, target::Symbol) 
    delete!(aquery.project, target; analytes = query.result)
    aquery.view ? (printstyled("Please re-query to get the correct result\n"; bold = true, color = :red); aquery.project) : aquery
end

function delete!(project::Project, target::Symbol; analytes = project.analytes)
    printstyled("Delete> ", uppercasefirst(string(target)), "\n"; color = :green, bold = true)
    deleteat!(analytes, findall(analyte -> analyte.states[states_id(target)] < 0, analytes))
    project
end

keys(project::Project) = LinearIndices(project.analytes)
keys(aquery::AbstractQuery) = LinearIndices(aquery.result)
keys(analyte::AnalyteSP) = LinearIndices(analyte.compounds)

reverse(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(deepcopy(project), start, stop)
reverse!(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(project.analytes, start, stop)
reverse(aquery::AbstractQuery, start::Int = 1, stop::Int = length(aquery)) = reverse!(copy_wo_project(aquery), start, stop)
reverse!(aquery::AbstractQuery, start::Int = 1, stop::Int = length(aquery)) = reverse!(aquery.result, start, stop)
reverse(analyte::AnalyteSP, start::Int = 1, stop::Int = length(analyte)) = reverse!(copy_wo_project(analyte), start, stop)
reverse!(analyte::AnalyteSP, start::Int = 1, stop::Int = length(analyte)) = reverse!(analyte.compounds, start, stop)

function union!(qs::Vararg{Query, N}) where N
    q = qs[1]
    length(qs) <= 1 && return q
    q.query = [q.query]
    if q.view
        ids = parentindices(q.result)[1]
        for qo in qs[2:end]
            push!(q.query, qo.query)
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

union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP) = union!(analyte, id, cpd2, sidechain(cpd2))
function union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP, ::SumChain)
    cpd1 = analyte[id]
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _, i = findmax((cpd1.area[1], cpd2.area[1]))
    cpd1.area = i ≡ 1 ? cpd1.area : cpd2.area
    analyte
end

function _union!(analyte::AnalyteSP, id::Int, ::SumChain, sc::SideChain)
    push!(last(analyte).project, copy_wo_project(analyte))
    for a in analyte
        a.sidechain = sc
    end
    analyte[id]
end

function _union!(analyte::AnalyteSP, id::Int, sc1::SideChain, sc2::SideChain)
    nox(lcb(sc1)) <= nox(lcb(sc2)) && (return analyte[id])
    for a in analyte
        a.sidechain = sc2
    end
    analyte[id]
end

function union!(analyte::AnalyteSP, id::Int, cpd2::CompoundSP, ::SideChain)
    cpd1 = _union!(analyte, id, sidechain(analyte[id]), sidechain(cpd2))
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _, i = findmax((cpd1.area[1], cpd2.area[1]))
    cpd1.area = i ≡ 1 ? cpd1.area : cpd2.area
    @p analyte sort!(; lt = isless_class)
    analyte.rt = calc_rt(analyte)
    analyte
end

union(cpd1::CompoundSP, cpd2::CompoundSP) = union!(copy_wo_project(cpd1), cpd2)
union(cpd1::CompoundSPVanilla, cpd2::CompoundSPVanilla) = union!(deepcopy(cpd1), cpd2)
function union!(cpd1::CompoundSPVanilla, cpd2::CompoundSPVanilla)
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _union!(cpd1, cpd2, cpd1.sidechain, cpd2.sidechain)
end

_union!(cpd1::CompoundID, cpd2::CompoundID, ::SumChain, ::SumChain) = cpd1
_union!(cpd1::CompoundID, cpd2::CompoundID, ::SideChain, ::SumChain) = cpd1
_union!(cpd1::CompoundID, cpd2::CompoundID, ::SumChain, sc::SideChain) = (cpd1.sidechain = sc; cpd1)
_union!(cpd1::CompoundID, cpd2::CompoundID, ::SideChain, sc::SideChain) = 
    nox(lcb(cpd1)) > nox(lcb(cpd2)) ? (cpd1.sidechain = sc; cpd1) : cpd1

function union!(cpd1::CompoundSP, cpd2::CompoundSP)
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _, i = findmax((cpd1.area[1], cpd2.area[1]))
    cpd1.area = i ≡ 1 ? cpd1.area : cpd2.area
    _union!(cpd1, cpd2, cpd1.sidechain, cpd2.sidechain)
end

union(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(copy_wo_project(analyte1), analyte2.compounds, analyte2.states)
union(analyte1::AnalyteSP, cpds::Vector{CompoundSP}) = union!(copy_wo_project(analyte1), cpds)
union!(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(analyte1, analyte2.compounds, analyte2.states)

function union!(analyte1::AnalyteSP, cpds::Vector{CompoundSP}, states2 = [0, 0, 0, 0, 0])
    for cpd2 in cpds
        id = findfirst(cpd1 -> iscompatible(cpd1, cpd2), analyte1)
        isnothing(id) ? push!(analyte1, cpd2) : union!(analyte1[id], cpd2)
    end
    sort!(analyte1, lt = isless_class)
    analyte1.states = min.(analyte1.states, states2)
    analyte1.scores .= 0
    analyte1
end

getproperty(reuseable::ReUseable, sym::Symbol) = sym ≡ :query ? getfield(reuseable, :query) : getfield(getfield(reuseable, :query), sym)

convert(::Type{CompoundSPVanilla}, cpd::CompoundSP) = CompoundSPVanilla(cpd.class, cpd.sidechain, cpd.fragments)
convert(::Type{CompoundSP}, cpd::CompoundSPVanilla) = CompoundSP(cpd.class, cpd.sidechain, cpd.fragments)
convert(::Type{SPID}, cpd::SPID) = cpd
# variant of interface
# iscomponent: whether ion is a component of cpd
iscomponent(ion::Ion{<: Adduct, <: ClassSP}, cpd::CompoundID) = iscompatible(ion.molecule, cpd.class)
iscomponent(ion::Ion{<: Adduct, <: SideChain}, cpd::CompoundID) = iscomponent(ion, cpd.sidechain)
iscomponent(ion::Ion{<: Adduct, <: SideChain}, sc::SumChain) = 
    ncb(ion) <= ncb(sc) && ndbox(ion) <= ndbox(sc) && 
        any(ox <= nox(sc) && db <= ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))

iscomponent(ion::Ion{<: Adduct, <: LCB}, sc::DiChain) = iscomponent(ion, sc.lcb)
iscomponent(ion::Ion{<: Adduct, <: ACYL}, sc::DiChain) = iscomponent(ion, sc.acyl)
iscomponent(ion::Ion{<: Adduct, <: LCB}, sc::LCB) =
    ncb(ion) ≡ ncb(sc) && ndbox(ion) ≡ ndbox(sc) && 
        any(ox ≡ nox(sc) && db ≡ ndb(sc) for (ox, db) in zip(nox(ion):nox(ion.molecule), ndb(ion):-1:ndb(ion.molecule)))
iscomponent(ion::Ion{<: Adduct, <: ACYL}, sc::ACYL) =
    ncb(ion) ≡ ncb(sc) && ndb(ion) ≡ ndb(sc) && nox(ion) ≡ ndb(sc)
    
iscomponent(::Ion{<: Adduct, NeuAc}, cpd::CompoundID) = isa(cpd.class, CLS.fg.nana)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{NeuAc, NeuAc}}}, cpd::CompoundID) = 
    isa(cpd.class, CLS.series.b) || isa(cpd.class, CLS.series.c) || isa(cpd.class, GT1a) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1c)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{NeuAc, NeuAc, NeuAc}}}, cpd::CompoundID) = isa(cpd.class, CLS.series.c)
iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, Hex}}}, cpd::CompoundID) = 
    (isa(cpd.class, CLS.series.as) && !isa(cpd.class, Hex2Cer)) ||
    (isa(cpd.class, CLS.series.a) && !isa(cpd.class, GM3)) ||
    (isa(cpd.class, CLS.series.b) && !isa(cpd.class, GD3)) ||
    (isa(cpd.class, CLS.series.c) && !isa(cpd.class, GT3)) ||
    isa(cpd.class, HexNAcHex2Cer) ||
    isa(cpd.class, HexNAcHex3Cer)

iscomponent(::Ion{S, Glycan{Tuple{HexNAc, Hex, NeuAc}}}, cpd::CompoundID) where {S <: Adduct} = 
    isa(cpd.class, CLS.fg.nana) && iscomponent(cpd, Ion(S(), Glycan(HexNAc(), Hex())))

iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, Hex, NeuAc, NeuAc}}}, cpd::CompoundID) = 
    (isa(cpd.class, CLS.series.b) && !isa(cpd.class, GD3)) ||
    (isa(cpd.class, CLS.series.c) && !isa(cpd.class, GT3)) ||
    isa(cpd.class, GT1a) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1c) || isa(cpd.class, GD1α)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{HexNAc, NeuAc}}}, cpd::CompoundID) = 
    isa(cpd.class, GP1cα) || isa(cpd.class, GQ1bα) || isa(cpd.class, GT1aα) || isa(cpd.class, GD1α)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{Hex, NeuAc}}}, cpd::CompoundID) = 
    isa(cpd.class, CLS.fg.nana)

iscomponent(::Ion{<: Adduct, Glycan{Tuple{Hex, NeuAc, NeuAc}}}, cpd::CompoundID) = 
    isa(cpd.class, CLS.series.b) ||
    isa(cpd.class, CLS.series.c) ||
    isa(cpd.class, GT1a) || 
    isa(cpd.class, GD1c)

# iscompatible: not exact, allow isomer
iscompatible(cpd1::CompoundID, cpd2::CompoundID) = 
    iscompatible(cpd1.class, cpd2.class) && iscompatible(cpd1.sidechain, cpd2.sidechain)

iscompatible(cls1::ClassSP, cls2::ClassSP) = cls1 ≡ cls2 || 
    any(cls1 ≡ cls2 for (cls1, cls2) in Iterators.product(isomer_tuple(cls1), isomer_tuple(cls2)))

iscompatible(x1, x2) = x1 ≡ x2
iscompatible(sc1::SumChain, sc2::SumChain) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::SumChain, sc2::SideChain) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::SideChain, sc2::SumChain) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::DiChain, sc2::DiChain) = iscompatible(sc1.lcb, sc2.lcb) && iscompatible(sc1.acyl, sc2.acyl)
iscompatible(sc1::LCB, sc2::LCB) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::ACYL, sc2::ACYL) = sc1 ≡ sc2
iscompatible(sc1::Acyl, sc2::Acyl) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::Acyl, sc2::ACYL) = sumcomp(sc1) ≡ sumcomp(sc2)
iscompatible(sc1::ACYL, sc2::Acyl) = sumcomp(sc1) ≡ sumcomp(sc2)

iscompatible(ion::Ion, cpd::CompoundID) = iscomponent(ion, cpd)
#=
isclasscompatible(cls1::ClassSP, cls2::ClassSP) = cls1 ≡ cls2 || begin
    class1 = hasisomer(cls1) ? cls1.isomer : (cls1, )
    class2 = hasisomer(cls2) ? cls2.isomer : (cls2, )
    any(cls1 ≡ cls2 for (cls1, cls2) in Iterators.product(class1, class2))
end

ischaincompatible(cpd1::CompoundID, cpd2::CompoundID) = 
    cpd1.sum ≡ cpd2.sum && ischaincompatible(cpd1.chain, cpd2.chain)

ischaincompatible(chain1::Nothing, chain2::Chain) = true
ischaincompatible(chain1::Chain, chain2::Nothing) = true
ischaincompatible(chain1::Nothing, chain2::Nothing) = true

ischaincompatible(chain1::Chain, chain2::Chain) = ischaincompatible(chain1.lcb, chain2.lcb) && ischaincompatible(chain1.acyl, chain2.acyl)
ischaincompatible(lcb1::LCB, lcb2::LCB) = 
    @match (lcb1, lcb2) begin
        ::Tuple{<: LCB2{N1, C}, <: LCB2{N2, C}} where {C, N1, N2}   => true
        ::Tuple{<: LCB3{N1, C}, <: LCB3{N2, C}} where {C, N1, N2}   => true
        ::Tuple{<: LCB4{N1, C}, <: LCB4{N2, C}} where {C, N1, N2}   => true
        _                                                           => false
    end

ischaincompatible(acyl1::ACYL, acyl2::ACYL) = 
    @match (acyl1, acyl2) begin
        ::Tuple{<: Acyl{N}, <: ACYL{N}} where N => true
        ::Tuple{<: ACYL{N}, <: Acyl{N}} where N => true
        _                                       => acyl1 ≡ acyl2
    end

ischainequal(chain1::Chain, chain2::Chain) = !isnothing(chain1) && !isnothing(chain2) && sumcomp(chain1.lcb) ≡ sumcomp(chain2.lcb) 
=#
copy_wo_project(cpd::CompoundSP) = CompoundSP(cpd.class, cpd.sidechain, deepcopy(cpd.fragments), cpd.area, deepcopy(cpd.states), deepcopy(cpd.results), cpd.project)
copy_wo_project(analyte::AnalyteSP) = AnalyteSP(copy_wo_project.(analyte.compounds), analyte.rt, deepcopy(analyte.states), deepcopy(analyte.scores), analyte.manual_check)
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

isless_class(cpd1::CompoundID) = true
isless_class(cpd1::CompoundID, cpd2::CompoundID) = 
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

sort(tbl::Table, i::Symbol; kwargs...) = tbl[sortperm(getproperty(tbl, i); kwargs...)]
function sort!(tbl::Table, i::Symbol; kwargs...)
    tbl[:] = tbl[sortperm(getproperty(tbl, i); kwargs...)]
end
sort(tbl::Table, v::Vector{Symbol}; kwargs...) = tbl[sortperm(collect(zip(getproperty.(Ref(tbl), v)...)); kwargs...)]
sort!(tbl::Table, v::Vector{Symbol}; kwargs...) = (tbl[:] = tbl[sortperm(collect(zip(getproperty.(Ref(tbl), v)...)); kwargs...)])
sort(tbl::Table, i::AbstractString; kwargs...) = sort(tbl, Symbol(i); kwargs...)
sort!(tbl::Table, i::AbstractString; kwargs...) = sort!(tbl, Symbol(i); kwargs...)
sort(tbl::Table, v::Vector{<: AbstractString}; kwargs...) = sort(tbl, Symbol.(v); kwargs...)
sort!(tbl::Table, v::Vector{<: AbstractString}; kwargs...) = sort!(tbl, Symbol.(v); kwargs...)
sort(tbl::Table, i::Int; kwargs...) = sort(tbl, propertynames(tbl)[i]; kwargs...)
sort!(tbl::Table, i::Int; kwargs...) = sort!(tbl, propertynames(tbl)[i]; kwargs...)
sort(tbl::Table, v::Vector{<: Int}; kwargs...) = sort(tbl, propertynames(tbl)[v]; kwargs...)
sort!(tbl::Table, v::Vector{<: Int}; kwargs...) = sort!(tbl, propertynames(tbl)[v]; kwargs...)