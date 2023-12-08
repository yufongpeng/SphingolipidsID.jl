"""
    transitiontable(project::Project, adduct, product, polarity::Bool; kwargs...)
    transitiontable(aquery::AbstractQuery, adduct, product, polarity::Bool; kwargs...)
    transitiontable(adduct, product, polarity::Bool, pq::Union{Project, AbstractQuery}; kwargs...)
    transitiontable(analytes::AbstractVector{<: AbstractAnalyteID}, adduct, product, polarity::Bool;
                    mz_tol = 0.35, rt_tol = 0.5,
                    db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG],
                    anion = :acetate,
                    default = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1)))
                    )

Create a table of MRM transitions.

# Arguuments
* `analytes`: analytes to be included.
* `adduct`: an adduct or a vector of adducts for parent ions.
* `product`: `Ion` or `LCB`; product ion. There are default adducts for `LCB`. 
* `polarity`: `true` for positive ion mode; `false` for negative ion mode.

# Keyword Arguuments
* `mz_tol`: Maximum allowable difference in m/z.
* `anion`: `:formate` or `acetate`.
* `default`: default product m/z value.
"""
transitiontable(project::Project, transitions; kwargs...) = transitiontable(project.analytes, transitions; kwargs..., anion = project.appendix[:anion])
transitiontable(aquery::AbstractQuery, transitions; kwargs...) = transitiontable(aquery.result, transitions; kwargs..., anion = aquery.project.appendix[:anion])
transitiontable(transitions, pq::Union{Project, AbstractQuery}; kwargs...) = transitiontable(pq, transitions; kwargs...)

function transitiontable(analytes::AbstractVector{<: AbstractAnalyteID}, transitions::Vector;
                        mz_tol = 0.35, rt_tol = 0.5,
                        anion = :acetate,
                        default_mz2 = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1)))
                        )
    polarity, transitions = zip(proc_transition.(transitions)...)
    polarity = all(polarity) ? true : any(polarity) ? throw(ArgumentError("Polairties of adducts must be the same")) : false
    analyte_list = analytelist(analytes, collect(transitions), polarity, anion)
    #filter_analytelist!(analyte_list, last(quantifier))
    mz2_list = calc_mz2(analyte_list, polarity)
    ce_list = celist(analyte_list, polarity)
    tbl = NamedTuple{(:analyte, :mz1, :mz2, :rt, :Δrt, :collision_energy), Tuple{Any, Float64, Float64, Float64, Float64, Int64}}[]
    for (ranalyte, ms2, ce) in zip(analyte_list, mz2_list, ce_list)
        if ms2 == 0
            isnothing(default_mz2) && continue
            ranalyte.prod == LCB || continue
            ms2 = default_mz2
        end
        pushed = false
        ms1 = mz(last(ranalyte.analyte), ranalyte.add)
        for (id, rtbl) in enumerate(tbl)
            (!between(rtbl.mz1, ms1, mz_tol) || !between(rtbl.mz2, ms2, mz_tol)) && continue
            #!between(ranalyte.rt - rt_tol, rtbl.rt, rtbl.Δrt / 2) && !between(ranalyte.rt + rt_tol, rtbl.rt, rtbl.Δrt / 2) && continue
            rtbl.collision_energy == ce || continue
            rt_l = min(rtbl.rt - rtbl.Δrt / 2, ranalyte.rt - rt_tol)
            rt_r = max(rtbl.rt + rtbl.Δrt / 2, ranalyte.rt + rt_tol)
            tbl[id] = (;
                analyte = rtbl.analyte,
                mz1 = (rtbl.mz1 * length(rtbl.analyte) + ms1) / (length(rtbl.analyte) + 1),
                mz2 = (rtbl.mz2 * length(rtbl.analyte) + ms2) / (length(rtbl.analyte) + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = rtbl.collision_energy
            )
            push!(tbl[id].analyte, repr(last(ranalyte.analyte)))
            pushed = true
            break
        end
        pushed || (push!(tbl, (analyte = [repr(last(ranalyte.analyte))], mz1 = ms1, mz2 = ms2, rt = ranalyte.rt, Δrt = rt_tol * 2, collision_energy = ce)))
    end
    Table((; id = collect(eachindex(tbl)), ), tbl; analyte = join.(getproperty.(tbl, :analyte), " | "), polarity = repeat([polarity ? "Positive" : "Negative"], size(tbl, 1)))
end

mrmquanttable(project::Project, quantifier::Tuple; kwargs...) = mrmquanttable(project.analytes, quantifier; kwargs..., mz_tol = last(project.data).config[:mz_tol], rt_tol = last(project.data).config[:rt_tol], anion = project.appendix[:anion])
mrmquanttable(aquery::AbstractQuery, quantifier::Tuple; kwargs...) = mrmquanttable(aquery.result, quantifier; kwargs..., mz_tol = last(aquery.project.data).config[:mz_tol], rt_tol = last(aquery.project.data).config[:rt_tol], anion = aquery.project.appendix[:anion])
mrmquanttable(quantifier::Tuple, pq::Union{Project, AbstractQuery}; kwargs...) = mrmquanttable(pq, quantifier; kwargs...)

function mrmquanttable(analytes::AbstractVector{<: AbstractAnalyteID}, quantifier::Tuple;
                        qualifier = Tuple[],
                        mz_tol = 0.35, rt_tol = 0.1,
                        anion = :acetate,
                        default_mz2 = mz(Ion(ProtonationNL2H2O(), SPB3(18, 1))),
                        coelution_fn = x -> lastindex(x.analyte),
                        isd_fn = x -> 0
                        )
    polarity, transitions = zip(proc_transition.([quantifier, qualifier...])...)
    polarity = all(polarity) ? true : any(polarity) ? throw(ArgumentError("Polairties of adducts must be the same")) : false
    analyte_list = analytelist(analytes, collect(transitions), polarity, anion)
    filter_analytelist!(analyte_list, last(quantifier))
    mz2_list = calc_mz2(analyte_list, polarity)
    tbl = NamedTuple{(:mz1, :mz2, :rt, :quantifier, :analyte), Tuple{Float64, Float64, Float64, Bool, Any}}[]
    prev = nothing
    sid = 1
    for (ranalyte, ms2) in zip(analyte_list, mz2_list)
        if ms2 == 0
            isnothing(default_mz2) && continue
            ranalyte.prod == LCB || continue
            ms2 = default_mz2
        end
        pushed = false
        prev, quant = ranalyte.analyte == prev ? (prev, false) : (ranalyte.analyte, true) 
        ms1 = mz(last(ranalyte.analyte), ranalyte.add)
        test = sid:length(tbl)
        if length(test) == 0 
            push!(tbl, (mz1 = ms1, mz2 = ms2, rt = ranalyte.rt, quantifier = quant, analyte = [ranalyte.analyte]))
            continue
        elseif !between(mean(getproperty.(tbl[test], :mz1)), ms1, mz_tol * 2) 
            push!(tbl, (mz1 = ms1, mz2 = ms2, rt = ranalyte.rt, quantifier = quant, analyte = [ranalyte.analyte]))
            sid = length(tbl)
            continue
        end
        for (id, rtbl) in zip(test, @view tbl[sid:end])
            (quant != rtbl.quantifier || !between(rtbl.mz1, ms1, mz_tol) || !between(rtbl.mz2, ms2, mz_tol) || !between(rtbl.rt, ranalyte.rt, rt_tol)) && continue
            tbl[id] = (;
                mz1 = (rtbl.mz1 * length(rtbl.analyte) + ms1) / (length(rtbl.analyte) + 1),
                mz2 = (rtbl.mz2 * length(rtbl.analyte) + ms2) / (length(rtbl.analyte) + 1),
                rt = (rtbl.rt * length(rtbl.analyte) + ranalyte.rt) / (length(rtbl.analyte) + 1),
                quantifier = quant,
                analyte = vectorize(rtbl.analyte)
            )
            push!(tbl[id].analyte, ranalyte.analyte)
            pushed = true
            break
        end
        pushed || (push!(tbl, (mz1 = ms1, mz2 = ms2, rt = ranalyte.rt, quantifier = quant, analyte = [ranalyte.analyte])))
    end
    tbl = Table((; id = collect(eachindex(tbl)), name = join.(map(analytes -> repr.(last.(analytes)), getproperty.(tbl, :analyte)), " | "),), tbl; 
            polarity = repeat([polarity ? "Positive" : "Negative"], size(tbl, 1)), 
            coelution_id = coelution_fn.(tbl),
            analyte = getproperty.(tbl, :analyte)
        )
    Table(tbl; isd_id = isd_fn.(tbl))
end

proc_transition(transition::Tuple{Symbol, T}) where T = (first(transition) == :default_pos, transition)
function proc_transition(transition::Tuple{Int, T}) where T
    add = SPDB[:ADDUCTCODE].object[first(transition)]
    (add isa Pos, (add, last(transition)))
end
function proc_transition(transition::Tuple{AbstractString, T}) where T 
    add = object_adduct(first(transition))
    (add isa Pos, (add, last(transition)))
end
proc_transition(transition::Tuple{Type{<: Adduct}, T}) where T  = (first(transition)() isa Pos, (first(transition)(), last(transition)))
proc_transition(transition::Tuple{Adduct, T}) where T  = (first(transition) isa Pos, transition)
proc_transition(transition::Tuple{Ion, T}) where T  = (first(transition).adduct isa Pos, transition)
function filter_adduct(adduct::Vector, anion::Symbol, polarity = nothing)
    if !isnothing(polarity) && polarity
        return filter(t -> isa(t, Pos), adduct)
    elseif !isnothing(polarity) && !polarity
        adduct = filter(t -> isa(t, Neg), adduct)
    end
    if anion ≡ :acetate
        filter(t -> !=(t, AddHCOO()), adduct)
    elseif anion ≡ :formate
        filter(t -> !=(t, AddOAc()), adduct)
    end
end

function analyte_transition_rt(analyte, transition, polarity, anion)
    add = first(transition) isa Symbol ? only(filter_adduct(class_db_index(class(analyte)).default_ion, anion, polarity)) : first(transition)
    (analyte = analyte, add = add, prod = last(transition), rt = analyte.rt)
end
analytelist(analytes::AbstractVector{<: AbstractAnalyteID}, transitions::Vector, polarity::Bool, anion) =
    productview((x, y) -> analyte_transition_rt(y, x, polarity, anion), transitions, analytes) |> splitdimsview |> flatten
    # Transducers
    # Iterators.product(adduct, analytes) |> MapSplat(analyte_add_rt) |> collect
    # DataPipes and BangBang:
    # @p Iterators.product(adduct, analytes) |> Iterators.map(analyte_add_rt(_...)) |> reduce(push!!; init = [])
    # Much more allocation for 1st
    # SplitApplyCombine
    # productview(analyte_add_rt, adduct, analytes) |> splitdimsview |> flatten
    # A little more allocation for 1st, less allocation for 2nd

#analytelist(analytes::AbstractVector{AbstractAnalyteID}, anion) =
#    @p analytes |> mapmany(analyte_add_rt(_, anion))
    # Transducers
    # analytes |> MapCat(analyte -> analyte_add_rt(analyte, polarity, anion)) |> collect
    # DataPipes and BangBang:
    # @p analytes |> Iterators.map(analyte_add_rt(_, polarity, anion)) |> reduce(vcat)
    # Much less allocation for 1st
    # DataPipes and SplitApplyCombine
    # Even less allocation for 1st

filter_analytelist!(analyte_list, transition) = analyte_list
filter_analytelist!(analyte_list, ::Ion{<: Adduct, NeuAc}) = filter!(analyte -> isa(class(analyte.analyte), CLS.fg.nana), analyte_list)
filter_analytelist!(analyte_list, ::LCB) = filter!(analyte -> !isspceieslevel(analyte.analyte), analyte_list)

calc_mz2(analyte_list::Vector, polarity::Bool) = 
    map(x -> calc_mz2(x, x.prod, polarity), analyte_list)
function calc_mz2(analyte::NamedTuple, products::Type{LCB}, polarity::Bool; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG])
    isspceieslevel(analyte.analyte) && return 0
    id = findfirst(x -> ≡(lcb(analyte.analyte), x.molecule), db[:, 1])
    isnothing(id) ? mz(default_adduct(lcb(analyte.analyte), polarity)) : db[id, 2]
end
calc_mz2(analyte::NamedTuple, ion::Ion{<: T}, polarity::Bool; db = SPDB[T <: Pos ? :FRAGMENT_POS : :FRAGMENT_NEG]) where T =
    isa(class(analyte.analyte), CLS.fg.nana) ? db[findfirst(==(ion), db[:, 1]), 2] : 0

celist(analyte_list::Vector, polarity::Bool) = 
    map(x -> celist(x, x.prod, polarity), analyte_list)
function celist(analyte::NamedTuple, ::Type{LCB}, polarity::Bool)
    polarity || throw(ArgumentError("Not yet implement this fragments"))
    id = findfirst(x -> ≡(class(analyte.analyte), SPDB[:CE].ms1[x]) && ≡(analyte.add, SPDB[:CE].adduct1[x]) && ==(SPDB[:CE].ms2[x], "LCB"), eachindex(SPDB[:CE]))
    isnothing(id) ? 40 : SPDB[:CE].eV[id]
end

function celist(analyte::NamedTuple, ion::Ion{<: Pos, NeuAc}, polarity::Bool)
    id = findfirst(x -> ≡(class(analyte.analyte), SPDB[:CE].ms1[x]) && ≡(analyte.add, SPDB[:CE].adduct1[x]) &&
                        ==(SPDB[:CE].ms2[x], "NeuAc") && ≡(ion.adduct, SPDB[:CE].adduct2[x]), eachindex(SPDB[:CE]))
    isnothing(id) ? 40 : SPDB[:CE].eV[id]
end