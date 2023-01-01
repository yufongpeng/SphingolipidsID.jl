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

function union!(project::Project, analyte::AnalyteSP, id::Int, cpd2::CompoundSP)
    cpd1 = analyte[id]
    if isnothing(cpd2.chain) 
        append!(cpd1.fragments, cpd2.fragments)
        unique!(cpd1.fragments)
        _, i = findmax((cpd1.area[1], cpd2.area[1]))
        cpd1.area = i == 1 ? cpd1.area : cpd2.area
        return cpd1
    end
    if isnothing(cpd1.chain)
        push!(project, copy_wo_project(analyte))
        analyte = last(project)
        for a in analyte
            a.chain = cpd2.chain
        end
        cpd1 = analyte[id]
    else
        chain = @match (_lcb(cpd1), _lcb(cpd2)) begin
            (::Tuple{<: LCB{N1, C}, <: LCB{N2, C}} where {C, N1, N2}) && if N1 > N2 end => cpd2.chain
            _                                                                           => cpd1.chain
        end
        for a in analyte
            a.chain = chain
        end
    end
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _, i = findmax((cpd1.area[1], cpd2.area[1]))
    cpd1.area = i == 1 ? cpd1.area : cpd2.area
    sort!(analyte, lt = isless_class)
    cpd1
end

union(cpd1::CompoundSP, cpd2::CompoundSP) = union!(copy_wo_project(cpd1), cpd2)

function union!(cpd1::CompoundSP, cpd2::CompoundSP)
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    _, i = findmax((cpd1.area[1], cpd2.area[1]))
    cpd1.area = i == 1 ? cpd1.area : cpd2.area
    if isnothing(cpd2.chain) 
        return cpd1
    end
    if isnothing(cpd1.chain)
        cpd1.chain = cpd2.chain
    else
        cpd1.chain = @match (_lcb(cpd1), _lcb(cpd2)) begin
            (::Tuple{<: LCB{N1, C}, <: LCB{N2, C}} where {C, N1, N2}) && if N1 > N2 end => cpd2.chain
            _                                                                           => cpd1.chain
        end
    end
    cpd1
end

union(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(copy_wo_project(analyte1), analyte2.compounds, analyte2.states)
union(analyte1::AnalyteSP, cpds::Vector{CompoundSP}) = union!(copy_wo_project(analyte1), cpds)
union!(analyte1::AnalyteSP, analyte2::AnalyteSP) = union!(analyte1, analyte2.compounds, analyte2.states)

function union!(analyte1::AnalyteSP, cpds::Vector{CompoundSP}, states2 = [0, 0, 0, 0, 0])
    for cpd2 in cpds
        id = findfirst(cpd1 -> iscompatible(cpd1, cpd2), analyte1)
        if isnothing(id)
            push!(analyte1, cpd2)
        else
            union!(analyte1[id], cpd2)
        end
    end
    sort!(analyte1, lt = isless_class)
    analyte1.states = min.(analyte1.states, states2)
    analyte1.scores .= 0
    analyte1
end

getproperty(reuseable::ReUseable, sym::Symbol) = sym == :query ? getfield(reuseable, :query) : getfield(getfield(reuseable, :query), sym)

# variant of interface
function iscomponent(cpd1::CompoundSP, ion::Ion{S, <: LCB}) where S 
    no = nhydroxyl(ion)
    cpd1.sum[3] >= no && (isnothing(cpd1.chain) || _lcb(cpd1) == ion.molecule)
end
iscomponent(cpd1::CompoundSP, ion::Ion{S, <: ACYL}) where S = iscompatible(cpd1, ion)
iscomponent(cpd1::CompoundSP, ::Ion{S, NeuAc}) where S = isa(cpd1.class, CLS.fg.nana)
iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{NeuAc, NeuAc}}}) where S = 
    isa(cpd1.class, CLS.series.b) || isa(cpd1.class, CLS.series.c) || isa(cpd1.class, GT1a) || isa(cpd1.class, GT1aα) || isa(cpd1.class, GD1c)
iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{NeuAc, NeuAc, NeuAc}}}) where S = isa(cpd1.class, CLS.series.c)
iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{HexNAc, Hex}}}) where S = 
    (isa(cpd1.class, CLS.series.as) && !isa(cpd1.class, Hex2Cer)) ||
    (isa(cpd1.class, CLS.series.a) && !isa(cpd1.class, GM3)) ||
    (isa(cpd1.class, CLS.series.b) && !isa(cpd1.class, GD3)) ||
    (isa(cpd1.class, CLS.series.c) && !isa(cpd1.class, GT3)) ||
    isa(cpd1.class, HexNAcHex2Cer) ||
    isa(cpd1.class, HexNAcHex3Cer)

iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{HexNAc, Hex, NeuAc}}}) where S = 
    isa(cpd1.class, CLS.fg.nana) && iscomponent(cpd1, Ion(S(), Glycan(HexNAc(), Hex())))

iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{HexNAc, Hex, NeuAc, NeuAc}}}) where S = 
    (isa(cpd1.class, CLS.series.b) && !isa(cpd1.class, GD3)) ||
    (isa(cpd1.class, CLS.series.c) && !isa(cpd1.class, GT3)) ||
    isa(cpd1.class, GT1a) || isa(cpd1.class, GT1aα) || isa(cpd1.class, GD1c) || isa(cpd1.class, GD1α)

iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{HexNAc, NeuAc}}}) where S = 
    isa(cpd1.class, GP1cα) || isa(cpd1.class, GQ1bα) || isa(cpd1.class, GT1aα) || isa(cpd1.class, GD1α)

iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{Hex, NeuAc}}}) where S = 
    isa(cpd1.class, CLS.fg.nana)

iscomponent(cpd1::CompoundSP, ::Ion{S, Glycan{Tuple{Hex, NeuAc, NeuAc}}}) where S = 
    isa(cpd1.class, CLS.series.b) ||
    isa(cpd1.class, CLS.series.c) ||
    isa(cpd1.class, GT1a) || 
    isa(cpd1.class, GD1c)

iscompatible(cpd1::CompoundSP, cpd2::CompoundSP) = 
    isclasscompatible(cpd1.class, cpd2.class) && ischaincompatible(cpd1, cpd2)

iscompatible(cpd1::CompoundSP, ion::Ion) = iscomponent(cpd1, ion)
iscompatible(cpd1::CompoundSP, ion::Ion{S, <: ClassSP}) where S = 
    isclasscompatible(cpd1.class, ion.molecule)

function iscompatible(cpd1::CompoundSP, ion::Ion{S, <: LCB}) where S 
    no = nhydroxyl(ion)
    cpd1.sum[3] >= no && (isnothing(cpd1.chain) || 
        nhydroxyl(_lcb(cpd1)) >= no && ischaincompatible(_lcb(cpd1), ion.molecule))
end

function iscompatible(cpd1::CompoundSP, ion::Ion{S, <: ACYL}) where S 
    no = nhydroxyl(ion)
    cpd1.sum[3] >= no && (isnothing(cpd1.chain) || 
        nhydroxyl(_acyl(cpd1)) == no)
end

isclasscompatible(cls1::ClassSP, cls2::ClassSP) = cls1 == cls2 || begin
    class1 = hasisomer(cls1) ? cls1.isomer : (cls1, )
    class2 = hasisomer(cls2) ? cls2.isomer : (cls2, )
    any(cls1 == cls2 for (cls1, cls2) in Iterators.product(class1, class2))
end

ischaincompatible(cpd1::CompoundSP, cpd2::CompoundSP) = 
    cpd1.sum == cpd2.sum && ischaincompatible(cpd1.chain, cpd2.chain)

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
        _                                       => acyl1 == acyl2
    end

ischainequal(chain1::Chain, chain2::Chain) = !isnothing(chain1) && !isnothing(chain2) && sumcomp(chain1.lcb) == sumcomp(chain2.lcb) 

copy_wo_project(cpd::CompoundSP) = CompoundSP(cpd.class, cpd.sum, cpd.chain, deepcopy(cpd.fragments), cpd.area, deepcopy(cpd.states), deepcopy(cpd.results), cpd.project)
copy_wo_project(analyte::AnalyteSP) = AnalyteSP(copy_wo_project.(analyte.compounds), analyte.rt, deepcopy(analyte.states), deepcopy(analyte.scores), analyte.manual_check)
copy_wo_project(aquery::Query) = Query(aquery.project, copy_wo_project.(aquery.result), deepcopy(aquery.query), false)
reuse_copy(aquery::Query) = Query(aquery.project, reuse_copy(aquery.result), deepcopy(aquery.query), true)
reuse_copy(v::Vector) = copy(v)
reuse_copy(v::T) where {T <: SubArray} = T(parent(v), v.indices, v.offset1, v.stride1)


equivalent_in(ion, collection) = any(equivalent(ion, x) for x in collection)
function equivalent(ion1::Ion{<: Pos, <: LCB}, ion2::Ion{<: Pos, <: LCB})
    nunsa(ion1.molecule) == nunsa(ion2.molecule) || return false
    nhydroxyl(ion1.molecule) - findfirst(==(ion1.adduct), SPDB[:NLH2O]) == nhydroxyl(ion2.molecule) - findfirst(==(ion2.adduct), SPDB[:NLH2O])
end

equivalent(ion1::Ion{S, T1}, ion2::Ion{S, T2}) where {S, T1 <: ClassSP, T2 <: ClassSP} = 
    isclasscompatible(ion1.molecule, ion2.molecule)

equivalent(ion1::Ion, ion2::Ion) = ==(ion1, ion2)
equivalent(ion1::ISF, ion2::Ion) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)
equivalent(ion1::Ion, ion2::ISF) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)
equivalent(ion1::ISF, ion2::ISF) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)

isless_class(cpd1::CompoundSP) = true
function isless_class(cpd1::CompoundSP, cpd2::CompoundSP)
    class1 = hasisomer(cpd1.class) ? cpd1.class.isomer : (cpd1.class, )
    class2 = hasisomer(cpd2.class) ? cpd2.class.isomer : (cpd2.class, )
    any(connected(cls2, cls1) for (cls1, cls2) in Iterators.product(class1, class2))
end
isless_ion(ion1) = true
function isless_ion(ion1::Ion{S, <: ClassSP}, ion2::Ion{T, <: ClassSP}) where {S, T}
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