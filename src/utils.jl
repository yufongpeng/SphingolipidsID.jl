# public
nfrags(cpd::CompoundSP) = size(cpd.fragments, 1)
function nMRM(tbl::DataFrame, precision::Float64 = 0.1; end_time = maximum(tbl[!, "Ret Time (min)"] .+ tbl[!, "Delta Ret Time"]))
    trans = Dict(0:precision:end_time .=> 0)
    for r in eachrow(tbl)
        for k in keys(trans)
            if k <= r["Ret Time (min)"] + r["Delta Ret Time"] / 2 && k >= r["Ret Time (min)"] - r["Delta Ret Time"] / 2
                trans[k] += 1
            end
        end
    end
    sort!(collect(trans), by = (x -> x.second))
end

# private
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
hasisomer(::Chain{S, Acyl{0}}) where S = false
hasisomer(::Chain{S, T}) where {S, T} = true
hasisomer(::Nothing) = true
hasisomer(::Type{T}) where T = hasisomer(T())

hasnana(::CLS.nana[0]) = false
hasnana(class) = true

# chain compositopn
nunsa(::LCB2) = 2
nunsa(::LCB3) = 3
nunsa(::LCB4) = 4

ncb(::LCB{N, C}) where {N, C} = C
nhydroxyl(::LCB{N, C}) where {N, C} = N
nhydroxyl(::ACYL{N}) where N = N
ndb(chain) = nunsa(chain) - nhydroxyl(chain)
sumcomp(chain) = (ncb(chain), ndb(chain), nhydroxyl(chain))

is4e(lcb::LCB{3}) = false
is4e(lcb::LCB2{2}) = false
is4e(lcb::LCB) = true

default_adduct(lcb::LCB{N, C}) where {N, C} = 
    is4e(lcb) ? Ion(SPDB[:NLH2O][N + 1], lcb) : Ion(SPDB[:NLH2O][N], lcb)

function hydroxyl_shift(adduct::Adduct, Δ::Int)
    id = findfirst(==(adduct), SPDB[:NLH2O]) + Δ
    0 < id <= length(SPDB[:NLH2O]) ? SPDB[:NLH2O][id] : nothing
end

hydroxyl_shift(lcb::LCB{N, C}, Δ::Int) where {N, C} = (0 <= N + Δ <= nunsa(lcb)) ? _hydroxyl_shift(lcb, Δ) : nothing 
_hydroxyl_shift(lcb::SPB2{N, C}, Δ::Int) where {N, C} = SPB2{N + Δ, C}()
_hydroxyl_shift(lcb::SPB3{N, C}, Δ::Int) where {N, C} = SPB3{N + Δ, C}()
_hydroxyl_shift(lcb::SPB4{N, C}, Δ::Int) where {N, C} = SPB4{N + Δ, C}()
_hydroxyl_shift(lcb::NotPhyto3{N, C}, Δ::Int) where {N, C} = SPB3{N + Δ, C}()
_hydroxyl_shift(lcb::NotPhyto4{N, C}, Δ::Int) where {N, C} = SPB4{N + Δ, C}()

# connection
find_connected(cpd::CompoundSP, analyte::AnalyteSP) = findall(id -> connected(cpd, id), analyte)

connected(cpd::CompoundSP, id::CompoundSP) = ischaincompatible(cpd, id) && begin
    class1 = hasisomer(cpd.class) ? cpd.class.isomer : (cpd.class, )
    class2 = hasisomer(id.class) ? id.class.isomer : (id.class, )
    any(connected(cls1, cls2) || connected(cls2, cls1) for (cls1, cls2) in Iterators.product(class1, class2))
end
connected(cls1::ClassSP, cls2::ClassSP) = cls1 == cls2 || any(connected(cls1, cls2) for cls1 in tuplize(SPDB[:CONNECTION][cls1]))
connected(cls1::ClassSP, cls2::Cer) = true
connected(cls1::Cer, cls2::ClassSP) = false
connected(cls1::Cer, cls2::Cer) = true

isf(cls::ClassSP) = hasisomer(cls) ? isf(cls.isomer) : push!(isf(SPDB[:CONNECTION][cls]), cls)
isf(cls::Tuple) = union!((push!(isf(SPDB[:CONNECTION][cls1]), cls1) for cls1 in cls)...)
isf(cls::Cer) = ClassSP[cls]

# user interface for query
lcb(c::Int, d::Int, n::Int) = Lcb(c, d, n)
acyl(c::Int, d::Int, n::Int) = Nacyl(c, d, n)
acylα(c::Int, d::Int, n::Int) = Nacylα(c, d, n)
acylβ(c::Int, d::Int, n::Int) = Nacylβ(c, d, n)
cpd(class::Type{<: ClassSP}, sum::NTuple{3, Int}, lcb::Lcb, acyl::NACYL) = CompoundID{class}(sum, lcb, acyl)
cpd(class::Type{<: ClassSP}, c::Int, n::Int, o::Int, lcb::Lcb, acyl::NACYL) = CompoundID{class}((c, n, o), lcb, acyl)
cpd(class::Type{<: ClassSP}, sum::NTuple{3, Int}) = CompoundID{class}(sum, nothing, nothing)
cpd(class::Type{<: ClassSP}, c::Int, n::Int, o::Int) = CompoundID{class}((c, n, o), nothing, nothing)
cpd(sum::NTuple{3, Int}, lcb::Lcb, acyl::NACYL) = CompoundID{Nothing}(sum, lcb, acyl)
cpd(c::Int, n::Int, o::Int, lcb::Lcb, acyl::NACYL) = CompoundID{Nothing}((c, n, o), lcb, acyl)
cpd(sum::NTuple{3, Int}) = CompoundID{Nothing}(sum, nothing, nothing)
cpd(c::Int, n::Int, o::Int) = CompoundID{Nothing}((c, n, o), nothing, nothing)
cpd(lcb::Lcb, acyl::NACYL) = cpd(Nothing, lcb, acyl)
cpd(class, lcb::Lcb, acyl::NACYL) = CompoundID{class}((lcb.cb + acyl.cb, lcb.db + acyl.db, lcb.ox + acyl.ox), lcb, acyl)

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
    dict = Dict{Union{ClassSP, Chain}}()
    for is in ions
        for i in is 
            dict[i] = get(dict, i, 0) + 1
        end
    end
    maxk, maxi = Union{ClassSP, Chain}[], 0
    for (k, v) in dict
        if maxi < v
            maxk = [k]
        elseif maxi == v
            push!(maxk, k)
        end
        maxi = max(maxi, v)
    end
    maxk
end

between(num::Number, range) = first(range) <= num <= last(range) 
between(num::Number; up, low) = low <= num <= up
between(num::Number, value, tol) = value - tol <= num <= value + tol
intersection(range...) = (maximum(first(r) for r in range), minimum(last(r) for r in range))
intersection(range::AbstractVector) = (maximum(first(r) for r in range), minimum(last(r) for r in range))

function equivalent_in_ion2(cpd::CompoundSP, criteria, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    equivalent_in(criteria, filterview(eqin, cpd.fragments).ion2)
end

function equivalent_in_ion2(cpd::CompoundSP, criteria::Tuple, prec)
    eqin = x -> equivalent_in(x.ion1, prec)
    any(equivalent_in(frag, criteria) for frag in filterview(eqin, cpd.fragments).ion2)
end

equivalent_in_ion2(analyte::AnalyteSP, criteria, prec) = any(equivalent_in_ion2(cpd, criteria, prec) for cpd in analyte)

equivalent_in_ion1(cpd::CompoundSP, criteria) = equivalent_in(criteria, cpd.fragments.ion1)
equivalent_in_ion1(cpd::CompoundSP, criteria::Tuple) = any(equivalent_in(frag, criteria) for frag in cpd.fragments.ion1)
equivalent_in_ion1(analyte::AnalyteSP, criteria) = any(equivalent_in_ion1(cpd, criteria) for cpd in analyte)

calc_rt(analyte::AnalyteSP) = mean(mean(query_raw(cpd.project, cpd.fragments.source[id], cpd.fragments.id[id]).rt for id in eachindex(cpd.fragments)) for cpd in analyte)
