# ==================================================================================
# Rule-based ID
"""
    apply_rules!(aquery::AbstractQuery; match_mode::Symbol = :both, kwargs...)
    apply_rules!(project::Project; analytes = project.analytes, 
                    match_mode::Symbol = :both,
                    class_mode::Symbol = :default,
                    chain_mode::Symbol = :isf,
                    class_rule::Symbol = :ab,
                    class_fail::Symbol = :pop_cpd,
                    chain_fail::Symbol = :ignore
                )

Apply identification rules to the analytes in a project or query.

* `match_mode`: `:both`, `:class`, or `:chain`.
* `class_mode`, `chain_mode`: determine which kind of ions is used. `:default` for default ions; `parent` for parent ions; `adduct` for adduct ions; `isf` for isf ions.
* `class_rule`: addtional rule for matching complex ganglioside. `:ab` for matching only a- and b-series gangliosides; `:all` for all gangliosides.
* `class_fail`: the next step after the class id is failed. `:pop_cpd` for deleting the current identification of the analyte and rematching with the previous compounds; `:ignore` for doing nothing.
* `chain_fail`: the next step after the chain id is failed. `:ignore` for doing nothing.
"""
apply_rules!(aquery::AbstractQuery; match_mode::Symbol = :both, kwargs...) = (apply_rules!(aquery.project; match_mode, analytes = aquery.result, kwargs...); aquery)
function apply_rules!(project::Project; analytes = project.analytes, 
                        match_mode::Symbol = :both,
                        class_mode::Symbol = :default,
                        chain_mode::Symbol = :isf,
                        class_rule::Symbol = :ab,
                        class_fail::Symbol = :pop_cpd,
                        chain_fail::Symbol = :ignore
                    )
    printstyled("ID> ", color = :green, bold = true)
    tomatch = @match match_mode begin
        :both   => (println("Class and Chain"); (true, true))
        :class  => (println("Class only"); (true, false))
        :chain  => (println("Chain only"); (false, true))
    end
    Threads.@threads for analyte in analytes
        _apply_rules!(analyte, tomatch, class_mode, chain_mode, class_rule, class_fail, chain_fail)
    end
    project
end

function _apply_rules!(analyte::AnalyteSP,
                        tomatch::Tuple{Bool, Bool},
                        class_mode::Symbol,
                        chain_mode::Symbol,
                        class_rule::Symbol,
                        class_fail::Symbol,
                        chain_fail::Symbol
                    )
    isempty(analyte) && return
    cpd = last(analyte)
    anion = cpd.project.appendix[:anion]
    # Class
    if tomatch[1]
        r = if class_mode ≡ cpd.results[1].mode
            cpd.results[1].rule
        elseif ≡(cpd.class, Cer())
            rule(cpd.class, cpd.chain)
        else
            rule(cpd.class, class_rule)
        end
        id_result_class = match_rules(analyte, generate_ms(cpd, class_mode, anion)..., class_mode, r)
        # Put in result
        cpd.results[1] = RuleSet(class_mode, id_result_class.rule)
        cpd.states[1] = analyte.states[states_id(:class)] = if !id_result_class.matched; 0
                                            elseif id_result_class.rule ≡ EmptyRule(); -1
                                            else cpd.class = id_result_class.rule; 1 end

        if analyte.states[states_id(:class)] < 0
            @match class_fail begin
                :pop_cpd    => ((length(analyte) ≡ 1 ? (return) : pop!(analyte)); _apply_rules!(analyte, tomatch, class_mode, chain_mode, class_rule, class_fail, chain_fail))
                :ignore     => nothing
            end
        end
    end
    # Chain
    if tomatch[2]
        id_result_chain = if chain_mode ≡ :isf
            id = findfirst(!isspceieslevel, @view analyte[end:-1:1])
            isnothing(id) ? Result(false, EmptyRule()) : match_rules(analyte,
                    generate_ms(cpd, chain_mode, anion)..., chain_mode, rule(analyte[end - id + 1].chain))
        elseif isspceieslevel(cpd)
            Result(false, EmptyRule())
        elseif chain_mode ≡ cpd.results[2].mode
            match_rules(analyte, generate_ms(cpd, chain_mode, anion)..., chain_mode, cpd.results[2].rule)
        else
            match_rules(analyte, generate_ms(cpd, chain_mode, anion)..., chain_mode, rule(cpd.chain))
        end
        # Put in result
        cpd.results[2] = RuleSet(chain_mode, id_result_chain.rule)
        cpd.states[2] = analyte.states[states_id(:chain)] = if !id_result_chain.matched
                                                @when PartialResult(_, _, r) = id_result_chain @inline cpd.chain = id_result_chain.result; 0
                                            elseif id_result_chain.rule ≡ EmptyRule(); -1
                                            else cpd.chain = id_result_chain.rule; 1 end
    # check isf
        prec, ms1 = generate_ms(cpd, :isf_sep, anion)
        for (id, cf) in enumerate(@view analyte[1:end - 1])
            isspceieslevel(cf) && continue
            r = @match cf.results[2].mode begin
                :isf_sep    => cf.results[2].rule
                _           => rule(cf.chain)
            end
            p = findfirst(ions -> iscompatible(ions.first, cf.class), prec)
            result = match_rules(analyte, prec[p].second, ms1[p].second, id, r)
            cf.results[2] = RuleSet(:isf_sep, id_result_chain.rule)
            cf.states[2] = if !result.matched
                                @when PartialResult(_, _, r) = result @inline cf.chain = result.result
                                iscompatible(cf.chain, cpd.chain) ? 0 : -1
                            elseif result.rule ≡ EmptyRule(); -1
                            else cf.chain = result.rule; iscompatible(cf.chain, cpd.chain) ? 1 : -1 end
        end
        if analyte.states[states_id(:chain)] < 0
            @match chain_fail begin
                :ignore => nothing
            end
        end
    end
    return
end

# One line match
#match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule) = match_rules(analyte, prec, ms1, mmode, rule, rule.criteria)

# Mutiline match
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::T) where {T <: Union{RuleUnion, RuleMode}}
    ions = Union{Tuple, AbstractRule}[]
    done = true
    for r in rule.rules
        result = match_rules(analyte, prec, ms1, mmode, r)
        done = done || result.matched
        @match result begin
            Result(false, x::AbstractRule)                                  => push!(ions, result.rule)
            Result(true, cls::T where T <: ClassSP) && if hasisomer(T) end  => push!(ions, result.rule.isomer)
            Result(true, ::EmptyRule)                                       => continue
            _                                                               => push!(ions, tuplize(result.rule))
        end
    end
    ions = @match (rule, done) begin
        (::RuleUnion, _)    => union((), ions...)
        (_, true)           => mode_ion(ions)
        (_, false)          => ions
    end
    done || return Result(false, T(ions; exception = rule.exception))
    @match ions begin
        []                          => match_rules(analyte, prec, ms1, mmode, rule.exception)
        [x]                         => Result(true, x)
        [x::S where S, xs...]       => Result(true, deisomerized(S)(ions))
    end
end

# match result
#match_rules(::AnalyteSP, ::Vector, ::Vector, mmode::Symbol, rule::Nothing) = (true, rule)
match_rules(::AnalyteSP, ::Vector, ::Vector, mmode, rule) = Result(true, rule)
match_rules(::AnalyteSP, ::Vector, ::Vector, mmode, rule::Tuple) = Result(true, rule)
match_rules(::AnalyteSP, ::Vector, ::Vector, mmode, rule::Tuple{}) = Result(true, EmptyRule())
match_rules(::AnalyteSP, ::Vector, ::Vector, mmode, rule::EmptyRule) = Result(true, rule)

# Y/N ==
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule{<: AbstractIon})
    cpd = get_cpd(analyte, mmode)
    found = isa(rule.criteria, ISF) ? equivalent_in_ion1(analyte, rule.criteria) :
                mmode ≡ :isf ? equivalent_in_ion2(analyte, rule.criteria, prec) : equivalent_in_ion2(cpd, rule.criteria, prec)
    if found
        match_rules(analyte, prec, ms1, mmode, first(rule.rule))
    elseif isa(rule.criteria, ISF) ? isf_checker(cpd, (rule.criteria,)) : ion_checker(cpd, prec, ms1, rule.criteria)
        match_rules(analyte, prec, ms1, mmode, last(rule.rule))
    else
        Result(false, rule)
    end
end

# Y/N in
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule{<: Tuple})
    cpd = get_cpd(analyte, mmode)
    found = isa(first(rule.criteria), ISF) ? equivalent_in_ion1(analyte, rule.criteria) :
        mmode ≡ :isf ? equivalent_in_ion2(analyte, rule.criteria, prec) : equivalent_in_ion2(cpd, rule.criteria, prec)
    if found
        match_rules(analyte, prec, ms1, mmode, first(rule.rule))
    elseif isa(first(rule.criteria), ISF) ? isf_checker(cpd, rule.criteria) : ion_checker(cpd, prec, ms1, rule.criteria)
        match_rules(analyte, prec, ms1, mmode, last(rule.rule))
    else
        Result(false, rule)
    end
end
#=
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule{<: AcylIon})
    cpd = get_cpd(analyte, mmode)
    found = mmode ≡ :isf ? equivalent_in_ion2(analyte, rule.criteria.ions, prec) : equivalent_in_ion2(cpd, rule.criteria.ions, prec)
    if found
        match_rules(analyte, prec, ms1, mmode, first(rule.rule))
    elseif ion_checker(cpd, prec, ms1, rule.criteria)
        match_rules(analyte, prec, ms1, mmode, last(rule.rule))
    else
        Result(false, rule)
    end
end
=#
# Comparison
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule{<: IonComparison})
    ratios = ion_comparison(analyte, prec, ms1, mmode, rule.criteria)
    isempty(ratios) && return Result(false, rule)
    for ratio in Iterators.reverse(ratios)
        for level in rule.rule
            result = @match level.first begin
                n::RealIntervals        => in(ratio, n)
                x && if isnan(x) end    => isnan(ratio)
                n::Number               => ==(ratio. n)
            end
            result && return match_rules(analyte, prec, ms1, mmode, level.second)
        end
    end
    Result(true, EmptyRule())
end

# oxygen
function match_rules(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::Rule{<: CalcOx})
    cpd = get_cpd(analyte, mmode)
    ox = nox(cpd) - nox(rule.criteria.lcb)
    for level in rule.rule
        if ox ≡ level.first
            result = match_rules(analyte, prec, ms1, mmode, level.second)
            lcb_new = rule.criteria.lcb
            return result.matched ? result : ndb(cpd) < ndb(lcb_new) ? Result(true, EmptyRule()) :
                PartialResult(false, result.rule, DiChain(lcb_new, Acyl(ncb(cpd) - ncb(lcb_new), ndb(cpd) - ndb(lcb_new), nox(cpd) - nox(lcb_new))))
        end
    end
    return Result(true, EmptyRule())
end

# ===========================================================================
get_cpd(analyte, n::Int) = analyte[n]
get_cpd(analyte, _) = last(analyte)

function generate_ms(cpd::CompoundSP, mmode::Symbol, anion::Symbol)
    add = @match mmode begin
        :isf        => return generate_ms_isf(cpd, anion)
        :isf_sep    => return generate_ms_isf_sep(cpd, anion)
        :adduct     => class_db_index(cpd.class).adduct_ion
        :parent     => class_db_index(cpd.class).parent_ion
        :default    => class_db_index(cpd.class).default_ion
    end
    if anion ≡ :acetate
        add = filter(!=(AddHCOO()), add)
    elseif anion ≡ :formate
        add = filter(!=(AddOAc()), add)
    end
    Ion.(add, Ref(cpd.class)), mz.(Ref(cpd), add)
end

function generate_ms_isf(cpd::CompoundSP, anion::Symbol)
    all_isf = isf(cpd.class)
    if anion ≡ :acetate
        fn = x -> filter(!=(AddHCOO()), x)
    elseif anion ≡ :formate
        fn = x -> filter(!=(AddOAc()), x)
    end
    ions = mapmany(all_isf) do isf_class
        Ion.(fn(class_db_index(isf_class).isf_ion), Ref(isf_class))
    end # transducer

    ions, mz.(Ref(cpd), ions)
end

function generate_ms_isf_sep(cpd::CompoundSP, anion::Symbol)
    all_isf = isf(cpd.class)
    if anion ≡ :acetate
        fn = x -> filter(!=(AddHCOO()), x)
    elseif anion ≡ :formate
        fn = x -> filter(!=(AddOAc()), x)
    end
    ions = map(all_isf) do isf_class
        isf_class => Ion.(fn(class_db_index(isf_class).isf_ion), Ref(isf_class))
    end

    ions, map(subions -> subions.first => mz.(Ref(cpd), subions.second), ions)
end
# ===================================================================================
#isf_checker(analyte::AnalyteSP, left::Tuple) = true

function isf_checker(cpd::CompoundSP, criteria)
    polarity = isa(first(criteria).adduct, Pos)
    ms1 = [mz(cpd, ion) for ion in criteria]
    any(cpd.project.data) do data
        ⊻(polarity, data.polarity) && return false
        table_id = isspceieslevel(cpd) ? eachindex(data.mz2) :
            findall(mz2 -> abs(rem(mz(Ion(Protonation(), lcb(cpd))) - mz2, mw("H2O"), RoundNearest)) < data.config[:mz_tol], data.mz2)
        isempty(table_id) && return false
        @match data begin
            ::PreIS => any(any(m1 in data.range[id] for m1 in ms1) for id in table_id)
            ::MRM   => any(any(between(m1, mz1, data.config[:mz_tol]) for m1 in ms1) for mz1 in (@p data.table filterview(in(_.id, table_id)) getproperty(:mz1)))
        end
    end
end

#ion_checker(precursor::Vector, ms1::Vector, rule::ISF) = true # Assume covered by PreIS

function ion_checker(cpd::CompoundSP, prec::Vector, ms1::Vector, rule::Ion{T}) where T
    P = @match T begin
        ::Type{<: Pos} => Pos
        ::Type{<: Neg} => Neg
    end
    id = findall(ion -> isa(ion.adduct, P), prec)
    ms1 = @view ms1[id]
    ms2 = mz(cpd, rule)
    polarity = ≡(P, Pos)
    any(cpd.project.data) do data
        ⊻(polarity, data.polarity) || return false
        table_id = findfirst(mz2 -> between(ms2, mz2, data.config[:mz_tol]), data.mz2)
        isnothing(table_id) && return false
        @match data begin
            ::PreIS => any(any(m1 in range for m1 in ms1) for range in @view data.range[table_id])
            ::MRM   => any(any(between(m1, mz1, data.config[:mz_tol]) for m1 in ms1) for mz1 in @views data.table.mz1[data.table.mz2_id .≡ table_id])
        end

    end
end

function separate_charge(ms1, prec)
    id1 = findall(ion -> isa(ion.adduct, Pos), prec)
    id2 = setdiff(eachindex(ms1), id1)
    (isempty(id1) ? empty(ms1) : ms1[id1], isempty(id2) ? empty(ms1) : ms1[id2])
end

function ion_checker(cpd::CompoundSP, prec::Vector, ms1::Vector, rule::Tuple)
    length(rule) ≡ 1 && return ion_checker(cpd, prec, ms1, rule)
    pos, neg = zip(separate_charge(ms1, prec), separate_charge(mz.(Ref(cpd), rule), rule))
    any(cpd.project.data) do data
        ms1, ms2 = data.polarity ? pos : neg
        table_id = findfirst(mz2 -> any(between(m2, mz2, data.config[:mz_tol]) for m2 in ms2), data.mz2)
        isnothing(table_id) && return false
        @match data begin
            ::PreIS => any(any(m1 in range for m1 in ms1) for range in @view data.range[table_id])
            ::MRM   => any(any(between(m1, mz1, data.config[:mz_tol]) for m1 in ms1) for mz1 in @views data.table.mz1[data.table.mz2_id .≡ table_id])
        end
    end
end
#=
function ion_checker(cpd::CompoundSP, prec::Vector, ms1::Vector, rule::AcylIon)
    id = findall(ion -> isa(ion.adduct, Neg), prec)
    ms1 = isempty(id) ? empty(ms1) : ms1[id]
    ms2 = mz(cpd, rule)
    any(cpd.project.data) do data
        data.polarity && return false
        table_id = findfirst(mz2 -> any(between(m2, mz2, data.config[:mz_tol]) for m2 in ms2), data.mz2)
        isnothing(table_id) && return false
        @match data begin
            ::PreIS => any(any(between(m1, range) for m1 in ms1) for range in @view data.range[table_id])
            ::MRM   => any(any(between(m1, mz1, data.config[:mz_tol]) for m1 in ms1) for mz1 in @views data.table.mz1[data.table.mz2_id .≡ rwa_id])
        end
    end
end
=#
function ion_comparison(analyte::AnalyteSP, prec::Vector, ms1::Vector, mmode, rule::IonComparison{S, T}) where {S, T}
    all_ions = map(rule.ions) do ions
        @match ions begin
            IonPlus(_)  => map(i -> (ion = i, mz = mz(i)), ions.ions)
            _           => [(ion = ions, mz = mz(ions))]
        end
    end
    ratios = Float64[]
    @match mmode begin
        :isf    => for cpd in analyte
                        _ion_comparison!(ratios, cpd, prec, ms1, all_ions)
                    end
        n::Int  => _ion_comparison!(ratios, analyte[n], prec, ms1, all_ions)
        _       => _ion_comparison!(ratios, last(analyte), prec, ms1, all_ions)
    end
    ratios
end

function _ion_comparison!(ratios::Vector{Float64}, cpd::CompoundSP, prec::Vector, ms1::Vector, all_ions)
    gdf = @p cpd.fragments |> filterview(equivalent_in(_.ion1, prec)) |> groupview(getproperties((:source, :ion1)))
    for subdf in gdf
        done = true
        ion1 = subdf.ion1[1]
        ms1_filter = ms1[findfirst(x -> equivalent(x, ion1), prec)]
        source = subdf.source[1]
        ratio = map(all_ions) do ions
            map(ions) do ion
                done || return 0
                id = findfirst(x -> equivalent(ion.ion, x), subdf.ion2)
                isnothing(id) || return query_data(cpd.project, source, subdf.id[id], cpd.project.appendix[:signal])
                ion_checker(cpd.project.data[source], ms1_filter, ion.mz) || (done = false)
                return 0
            end
        end
        done && push!(ratios, mapreduce(sum, /, ratio))
    end
    ratios
end

ion_checker(data::PreIS, ms1::Float64, ms2::Float64) =
    any(between(mz2, ms2, data.config[:mz_tol]) && ms1 in range for (mz2, range) in zip(data.mz2, data.range))
ion_checker(data::MRM, ms1::Float64, ms2::Float64) =
    any(between(mz2, ms2, data.config[:mz_tol]) && any(between(mz1, ms1, data.config[:mz_tol]) for mz1 in @view data.table.mz1[data.table.mz2 .≡ i]) for (i, mz2) in enumerate(data.mz2))