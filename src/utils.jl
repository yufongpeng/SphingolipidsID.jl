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

ncb(::LCB{C, N}) where {C, N} = C
nhydroxyl(::LCB{C, N}) where {C, N} = N
ndb(lcb) = nunsa(lcb) - nhydroxyl(lcb)
sumcomp(lcb::LCB) = (ncb(lcb), ndb(lcb), nhydroxyl(lcb))

isphyto(::PhytoSPB) = true
isphyto(lcb::LCB{C, N}) where {C, N} = (nunsa(lcb) - N) == 0

default_adduct(lcb::LCB{C, N}) where {C, N} = 
    isphyto(lcb) ? Ion(NLH2O[N], lcb) : Ion(NLH2O[N + 1], lcb)

function hydroxyl_shift(adduct::Adduct, Δ::Int)
    id = findfirst(==(adduct), NLH2O) + Δ
    0 < id <= length(NLH2O) ? NLH2O[id] : nothing
end

hydroxyl_shift(lcb::LCB{C, N}, Δ::Int) where {C, N} = (0 <= N + Δ <= nunsa(lcb)) ? _hydroxyl_shift(lcb, Δ) : nothing 
_hydroxyl_shift(lcb::SPB2{C, N}, Δ::Int) where {C, N} = SPB2{C, N + Δ}()
_hydroxyl_shift(lcb::SPB3{C, N}, Δ::Int) where {C, N} = SPB3{C, N + Δ}()
_hydroxyl_shift(lcb::SPB4{C, N}, Δ::Int) where {C, N} = SPB4{C, N + Δ}()
_hydroxyl_shift(lcb::PhytoSPB3{C, N}, Δ::Int) where {C, N} = PhytoSPB3{C, N + Δ}()
_hydroxyl_shift(lcb::PhytoSPB4{C, N}, Δ::Int) where {C, N} = PhytoSPB4{C, N + Δ}()

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
                ::Tuple{<: LCB2{C, N1}, <: LCB2{C, N2}} where {C, N1, N2}   => true
                ::Tuple{<: LCB3{C, N1}, <: LCB3{C, N2}} where {C, N1, N2}   => true
                ::Tuple{<: LCB4{C, N1}, <: LCB4{C, N2}} where {C, N1, N2}   => true
                _                                                           => false
            end
        end
    end
end

function union!(cpd1::CompoundGSL, cpd2::CompoundGSL)
    append!(cpd1.fragments, cpd2.fragments)
    unique!(cpd1.fragments)
    isnothing(cpd2.chain) && return cpd1
    if isnothing(cpd1.chain)
        cpd1.chain = cpd2.chain
    else
        cpd1.chain = @match (cpd1.chain.lcb, cpd2.chain.lcb) begin
            (::Tuple{<: LCB{C, N1}, <: LCB{C, N2}} where {C, N1, N2}) && if N1 > N2 end => cpd2.chain
            _                                                                           => cpd1.chain
        end
    end
    cpd1
end

equivalent_in(ion, collection) = in(ion, collection)
equivalent_in(ion::Union{ISF, Ion{<: Pos, <: LCB}}, collection) = any(equivalent(ion, x) for x in collection)
function equivalent(ion1::Ion{<: Pos, <: LCB}, ion2::Ion{<: Pos, <: LCB})
    nunsa(ion1.molecule) == nunsa(ion2.molecule) || return false
    nhydroxyl(ion1.molecule) - findfirst(==(ion1.adduct), NLH2O) == nhydroxyl(ion2.molecule) - findfirst(==(ion2.adduct), NLH2O)
end

equivalent(ion1::Ion, ion2::Ion) = ==(ion1, ion2)
equivalent(ion1::ISF, ion2::Ion) = ==(ion1.adduct, ion2.adduct) && ==(ion1.molecule, ion2.molecule)

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
vectorize(x) = [x]

tuplize(x::Tuple) = x
tuplize(x) = (x,)