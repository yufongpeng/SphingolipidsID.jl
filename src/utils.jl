deisomerized(::ClassGM1) = GM1()
deisomerized(::ClassGD1) = GD1()
deisomerized(::ClassGT1) = GT1()
deisomerized(::ClassGQ1) = GQ1()
deisomerized(::ClassGP1) = GP1()
deisomerized(::ClassHexNAcHex3Cer) = HexNAcHex3Cer()
deisomerized(x::ClassGSL) = x
deisomerized(::Type{<: ClassGM1}) = GM1
deisomerized(::Type{<: ClassGD1}) = GD1
deisomerized(::Type{<: ClassGT1}) = GT1
deisomerized(::Type{<: ClassGQ1}) = GQ1
deisomerized(::Type{<: ClassGP1}) = GP1
deisomerized(::Type{<: ClassHexNAcHex3Cer}) = HexNAcHex3Cer
deisomerized(x::Type{<: ClassGSL}) = x

hasisomer(::Type{GM1}) = true
hasisomer(::Type{GD1}) = true
hasisomer(::Type{GT1}) = true
hasisomer(::Type{GQ1}) = true
hasisomer(::Type{GP1}) = true
hasisomer(::Type{HexNAcHex3Cer}) = true
hasisomer(x::Type{<: ClassGSL}) = false

hasnana(::ClasswoNANA) = false
hasnana(class) = true

nunsa(::LCB2) = 2
nunsa(::LCB3) = 3
nunsa(::LCB4) = 4

ncb(::LCB{N, C}) where {N, C} = C
nhydroxyl(::LCB{N, C}) where {N, C} = N
ndb(lcb) = nunsa(lcb) - nhydroxyl(lcb)
sumcomp(lcb::LCB) = (ncb(lcb), ndb(lcb), nhydroxyl(lcb))

isphyto(::PhytoSPB) = true
isphyto(lcb::LCB{N, C}) where {N, C} = (nunsa(lcb) - N) == 0

default_adduct(lcb::LCB{N, C}) where {N, C} = 
    isphyto(lcb) ? Ion(SPDB[:NLH2O][N], lcb) : Ion(SPDB[:NLH2O][N + 1], lcb)

function hydroxyl_shift(adduct::Adduct, Δ::Int)
    id = findfirst(==(adduct), SPDB[:NLH2O]) + Δ
    0 < id <= length(SPDB[:NLH2O]) ? SPDB[:NLH2O][id] : nothing
end

hydroxyl_shift(lcb::LCB{N, C}, Δ::Int) where {N, C} = (0 <= N + Δ <= nunsa(lcb)) ? _hydroxyl_shift(lcb, Δ) : nothing 
_hydroxyl_shift(lcb::SPB2{N, C}, Δ::Int) where {N, C} = SPB2{N + Δ, C}()
_hydroxyl_shift(lcb::SPB3{N, C}, Δ::Int) where {N, C} = SPB3{N + Δ, C}()
_hydroxyl_shift(lcb::SPB4{N, C}, Δ::Int) where {N, C} = SPB4{N + Δ, C}()
_hydroxyl_shift(lcb::PhytoSPB3{N, C}, Δ::Int) where {N, C} = PhytoSPB3{N + Δ, C}()
_hydroxyl_shift(lcb::PhytoSPB4{N, C}, Δ::Int) where {N, C} = PhytoSPB4{N + Δ, C}()

function mode(ions)
    dict = Dict{Union{ClassGSL, Chain}}()
    for is in ions
        for i in is 
            dict[i] = get(dict, i, 0) + 1
        end
    end
    maxk, maxi = Union{ClassGSL, Chain}[], 0
    for (k, v) in dict
        if maxi < v
            maxk = [k]
        elseif maxi == v
            push!(maxk, k)
        end
    end
    length(maxk) > 1 ? typeof(deisomerized(maxk[1]))(maxk) : first(maxk) 
end

iscompatible(cpd1::CompoundGSL, cpd2::CompoundGSL) = 
    cpd1.class == cpd2.class && ischaincompatible(cpd1, cpd2)

function ischaincompatible(cpd1::CompoundGSL, cpd2::CompoundGSL)
    cpd1.sum == cpd2.sum && begin 
        isnothing(cpd1.chain) || isnothing(cpd2.chain) || begin
            @match (cpd1.chain.lcb, cpd2.chain.lcb) begin
                ::Tuple{<: LCB2{N1, C}, <: LCB2{N2, C}} where {C, N1, N2}   => true
                ::Tuple{<: LCB3{N1, C}, <: LCB3{N2, C}} where {C, N1, N2}   => true
                ::Tuple{<: LCB4{N1, C}, <: LCB4{N2, C}} where {C, N1, N2}   => true
                _                                                           => false
            end
        end
    end
end

length(project::Project) = length(project.analytes)
length(aquery::Query) = length(aquery.result)
length(analyte::AnalyteGSL) = length(analyte.compounds)

iterate(project::Project) = iterate(project.analytes)
iterate(project::Project, i) = iterate(project.analytes, i)
iterate(aquery::Query) = iterate(aquery.result)
iterate(aquery::Query, i) = iterate(aquery.result, i)
iterate(analyte::AnalyteGSL) = iterate(analyte.compounds)
iterate(analyte::AnalyteGSL, i) = iterate(analyte.compounds, i)

getindex(project::Project, i) = getindex(project.analytes, i)
getindex(aquery::Query, i) = getindex(aquery.result, i)
getindex(analyte::AnalyteGSL, i) = getindex(analyte.compounds, i)

view(project::Project, i) = view(project.analytes, i)
view(aquery::Query, i) = view(aquery.result, i)
view(analyte::AnalyteGSL, i) = view(analyte.compounds, i)

firstindex(project::Project) = firstindex(project.analytes)
firstindex(aquery::Query) = firstindex(aquery.result)
firstindex(analyte::AnalyteGSL) = firstindex(analyte.compounds)

lastindex(project::Project) = lastindex(project.analytes)
lastindex(aquery::Query) = lastindex(aquery.result)
lastindex(analyte::AnalyteGSL) = lastindex(analyte.compounds)

sort(project::Project; kwargs...) = sort!(deepcopy(project); kwargs...)
sort!(project::Project; kwargs...) = sort!(project.analytes; kwargs...)
sort(aquery::Query; kwargs...) = sort!(deepcopy(aquery); kwargs...)
sort!(aquery::Query; kwargs...) = sort!(aquery.result; kwargs...)
sort(analyte::AnalyteGSL; kwargs...) = sort!(deepcopy(analyte); kwargs...)
sort!(analyte::AnalyteGSL; kwargs...) = sort!(analyte.compounds; kwargs...)

push!(project::Project, analyte::AnalyteGSL) = push!(project.analytes, analyte)
push!(analyte::AnalyteGSL, cpd::CompoundGSL) = push!(analyte.compounds, cpd)

pop!(project::Project) = pop!(project.analytes)
pop!(analyte::AnalyteGSL) = pop!(analyte.compounds)

popfirst!(project::Project) = popfirst!(project.analytes)
popfirst!(analyte::AnalyteGSL) = popfirst!(analyte.compounds)

popat!(project::Project, del::Vector{Int}) = popat!(project.analytes, del)
popat!(analyte::AnalyteGSL, del::Vector{Int}) = popat!(analyte.compounds, del)
popat!(project::Project, del::Int) = popat!(project.analytes, del)
popat!(analyte::AnalyteGSL, del::Int) = popat!(analyte.compounds, del)

deleteat!(project::Project, del::Vector{Int}) = deleteat!(project.analytes, del)
deleteat!(analyte::AnalyteGSL, del::Vector{Int}) = deleteat!(analyte.compounds, del)
deleteat!(project::Project, del::Int) = deleteat!(project.analytes, del)
deleteat!(analyte::AnalyteGSL, del::Int) = deleteat!(analyte.compounds, del)

deleteat!(analytes::SubArray{AnalyteGSL, 1, Vector{AnalyteGSL}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = deleteat!(parent(analytes), parentindices(analytes)[1][del])
deleteat!(analytes::SubArray{AnalyteGSL, 1, Vector{AnalyteGSL}, Tuple{Vector{Int64}}, false}, del::Int) = deleteat!(parent(analytes), parentindices(analytes)[1][del])
popat!(analytes::SubArray{AnalyteGSL, 1, Vector{AnalyteGSL}, Tuple{Vector{Int64}}, false}, del::Vector{Int}) = popat!(parent(analytes), parentindices(analytes)[1][del])
popat!(analytes::SubArray{AnalyteGSL, 1, Vector{AnalyteGSL}, Tuple{Vector{Int64}}, false}, del::Int) = popat!(parent(analytes), parentindices(analytes)[1][del])

keys(project::Project) = LinearIndices(project.analytes)
keys(aquery::Query) = LinearIndices(aquery.result)
keys(analyte::AnalyteGSL) = LinearIndices(analyte.compounds)

reverse(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(deepcopy(project), start, stop)
reverse!(project::Project, start::Int = 1, stop::Int = length(project)) = reverse!(project.analytes, start, stop)
reverse(aquery::Query, start::Int = 1, stop::Int = length(aquery)) = reverse!(deepcopy(aquery), start, stop)
reverse!(aquery::Query, start::Int = 1, stop::Int = length(aquery)) = reverse!(aquery.result, start, stop)
reverse(analyte::AnalyteGSL, start::Int = 1, stop::Int = length(analyte)) = reverse!(deepcopy(analyte), start, stop)
reverse!(analyte::AnalyteGSL, start::Int = 1, stop::Int = length(analyte)) = reverse!(analyte.compounds, start, stop)

function union!(project::Project, analyte::AnalyteGSL, id::Int, cpd2::CompoundGSL)
    cpd1 = analyte[id]
    if isnothing(cpd2.chain) 
        append!(cpd1.fragments, cpd2.fragments)
        unique!(cpd1.fragments)
        cpd1.area = max(cpd1.area, cpd2.area)
        return cpd1
    end
    if isnothing(cpd1.chain)
        push!(project, deepcopy(analyte))
        analyte = last(project)
        for a in analyte
            a.chain = cpd2.chain
        end
        cpd1 = analyte[id]
    else
        chain = @match (cpd1.chain.lcb, cpd2.chain.lcb) begin
            (::Tuple{<: LCB{N1, C}, <: LCB{N2, C}} where {C, N1, N2}) && if N1 > N2 end => cpd2.chain
            _                                                                           => cpd1.chain
        end
        for a in analyte
            a.chain = chain
        end
    end
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    cpd1.area = max(cpd1.area, cpd2.area)
    sort!(analyte, lt = isless_class)
    cpd1
end

union(cpd1::CompoundGSL, cpd2::CompoundGSL) = union!(deepcopy(cpd1), cpd2)

function union!(cpd1::CompoundGSL, cpd2::CompoundGSL)
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    cpd1.area = max(cpd1.area, cpd2.area)
    if isnothing(cpd2.chain) 
        return cpd1
    end
    if isnothing(cpd1.chain)
        cpd1.chain = cpd2.chain
    else
        cpd1.chain = @match (cpd1.chain.lcb, cpd2.chain.lcb) begin
            (::Tuple{<: LCB{N1, C}, <: LCB{N2, C}} where {C, N1, N2}) && if N1 > N2 end => cpd2.chain
            _                                                                           => cpd1.chain
        end
    end
    cpd1
end

union(analyte1::AnalyteGSL, analyte2::AnalyteGSL) = union!(deepcopy(analyte1), analyte2.identi)
union(analyte1::AnalyteGSL, cpds::Vector{CompoundGSL}) = union!(deepcopy(analyte1), cpds)
union!(analyte1::AnalyteGSL, analyte2::AnalyteGSL) = union!(analyte1, analyte2.compounds)

function union!(analyte1::AnalyteGSL, cpds::Vector{CompoundGSL})
    for cpd2 in cpds
        id = findfirst(cpd1 -> iscompatible(cpd1, cpd2), analyte1)
        if isnothing(id)
            push!(analyte1, cpd2)
        else
            union!(analyte1[id], cpd2)
        end
    end
    sort!(analyte1, lt = isless_class)
    analyte1
end

equivalent_in(ion, collection) = any(equivalent(ion, x) for x in collection)
function equivalent(ion1::Ion{<: Pos, <: LCB}, ion2::Ion{<: Pos, <: LCB})
    nunsa(ion1.molecule) == nunsa(ion2.molecule) || return false
    nhydroxyl(ion1.molecule) - findfirst(==(ion1.adduct), SPDB[:NLH2O]) == nhydroxyl(ion2.molecule) - findfirst(==(ion2.adduct), SPDB[:NLH2O])
end

equivalent(ion1::Ion, ion2::Ion) = ==(ion1, ion2)
equivalent(ion1::ISF, ion2::Ion) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)
equivalent(ion1::Ion, ion2::ISF) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)
equivalent(ion1::ISF, ion2::ISF) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)

isless_class(cpd1::CompoundGSL) = true
function isless_class(cpd1::CompoundGSL, cpd2::CompoundGSL)
    class1 = haskey(CONNECTION, cpd1.class) ? cpd1.class : deisomerized(cpd1.class)
    class2 = haskey(CONNECTION, cpd2.class) ? cpd2.class : deisomerized(cpd2.class)
    class1 == class2 ? false : connected(class2, class1)
end
isless_ion(ion1) = true
function isless_ion(ion1, ion2)
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

between(num::Number, range) = first(range) <= num <= last(range) 
between(num::Number; up, low) = low <= num <= up
between(num::Number, value, tol) = value - tol <= num <= value + tol
intersection(range...) = (maximum(first(r) for r in range), minimum(last(r) for r in range))
intersection(range::AbstractVector) = (maximum(first(r) for r in range), minimum(last(r) for r in range))
isempty(::ClassGSL) = false

vectorize(x::AbstractVector) = x
vectorize(x::UnitRange) = x
vectorize(x) = [x]

tuplize(x::Tuple) = x
tuplize(x) = (x,)