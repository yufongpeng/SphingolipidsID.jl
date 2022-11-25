# ==================================================================================
# Rule-based ID
# start point
apply_rules!(aquery::Query, match_mode::Symbol = :both; kwargs...) = (apply_rules!(aquery.project, part, aquery.result; kwargs...); aquery)

function apply_rules!(project::Project, match_mode::Symbol = :both, analytes = project.analytes; 
                        class_mode::Symbol = :default, 
                        chain_mode::Symbol = :isf, 
                        class_rule::Symbol = :ab,
                        class_fail::Symbol = :pop_cpd, 
                        chain_fail::Symbol = :ignore
                    )
    printstyled("ID> ", color = :green, bold = true)
    match_mode = @match match_mode begin
        :both   => (println("Class and Chain"); (true, true))
        :class  => (println("Class only"); (true, false))
        :chain  => (println("Chain only"); (false, true))
    end
    i = 1
    n = length(analytes)
    while i <= n
        analyte = analytes[i]
        isempty(analyte) && (i += 1; continue)
        cpd = last(analyte)
        # Class
        if match_mode[1]
            r = if class_mode == cpd.results[1].mode 
                cpd.results[1].rule 
            elseif ==(cpd.class, Cer()) && ((cpd.sum[3] < 1) || (!isnothing(cpd.chain) && (isa(cpd.chain.lcb, SPB2))))
                rule(cpd.class, false)
            else
                rule(cpd.class, class_rule)
            end
            id_result_class = match_rules(project, analyte, generate_ms(cpd, class_mode, project.anion)..., class_mode, r)               
            # Put in result
            cpd.results[1] = RuleSet(class_mode, id_result_class.rule)
            cpd.states[1] = analyte.states[1] = if !id_result_class.matched; 0
                                                elseif id_result_class.rule === EmptyRule(); -1
                                                else cpd.class = id_result_class.rule; 1 end

            if analyte.states[1] < 0
                @match class_fail begin
                    :pop_cpd    => begin
                        if length(analyte) == 1
                            i += 1
                        else
                            pop!(analyte)
                        end
                        continue
                    end 
                    :ignore     => nothing
                end
            end
        end                                    
        # Chain
        if match_mode[2]
            id_result_chain = if chain_mode == cpd.results[2].mode
                match_rules(project, analyte, generate_ms(cpd, chain_mode, project.anion)..., chain_mode, cpd.results[2].rule)
            elseif chain_mode == :isf  
                id = findfirst(cpd -> !isnothing(cpd.chain), reverse(analyte))
                isnothing(id) ? Result(false, EmptyRule()) : match_rules(project, analyte, 
                        generate_ms(cpd, chain_mode, project.anion)..., chain_mode, rule(analyte[end - id + 1].chain))
            elseif isnothing(cpd.chain)
                Result(false, EmptyRule())
            else
                match_rules(project, analyte, generate_ms(cpd, chain_mode, project.anion)..., chain_mode, rule(cpd.chain))
            end
            # Put in result
            cpd.results[2] = RuleSet(chain_mode, id_result_chain.rule)
            cpd.states[2] = analyte.states[2] = if !id_result_chain.matched
                                                    @when PartialResult(_, _, r) = id_result_chain @inline cpd.chain = id_result_chain.result; 0
                                                elseif id_result_chain.rule === EmptyRule(); -1
                                                else cpd.chain = id_result_chain.rule; 1 end
        # check isf
            prec, ms1 = generate_ms(cpd, :isf_sep, project.anion)
            chain = cpd.chain
            for (id, cpd) in enumerate(@view analyte[1:end - 1])
                isnothing(cpd.chain) && continue                
                r = @match cpd.results[2].mode begin
                    :isf_sep => cpd.results[2].mode
                    _       => rule(cpd.chain)
                end
                p = findfirst(ions -> isclasscompatible(ions.first, cpd.class), prec)
                result = match_rules(project, analyte, prec[p].second, ms1[p].second, id, r)
                cpd.results[2] = RuleSet(:isf_sep, id_result_chain.rule)
                cpd.states[2] = if !result.matched
                                    @when PartialResult(_, _, r) = result @inline cpd.chain = result.result
                                    ischaincompatible(cpd.chain, chain) ? 0 : -1
                                elseif result.rule === EmptyRule(); -1
                                else cpd.chain = result.rule; ischainequal(cpd.chain, chain) ? 1 : -1 end
            end
            if analyte.states[2] < 0
                @match chain_fail begin
                    :ignore => nothing
                end
            end
        end
        i += 1
    end
    project
end
#=
function filter_score!(project::Project, score)
end


nq = -(mapreduce(isf_cpd -> !isnothing(isf_cpd.states[2]) && ischaincompatible(isf_cpd.chain, cpd.chain), +, analyte), 
        mapreduce(isf_cpd -> isa(isf_cpd.states[2], Symbol) && !ischainequal(isf_cpd.chain, cpd.chain), +, analyte))
ndq = length(analyte) - nq
analyte.states[2] = (ndq < first(chain!_isf) && ndq / (ndq + nq) < last(chain!_isf)) ? 1 : -1

if fail[2]
    @match chain_fail begin
        :del    => push!(del, i)
        :ignore => nothing
    end
end

if fail[1]
    @match class_fail begin
        :del        => (push!(del, i); i+= 1; continue)
        :pop_cpd    => begin
            pop!(analyte)
            isempty(analyte) && (push!(del, i); i += 1)
            continue
        end 
        :ignore     => nothing
    end
end

                        chain!_isf::Tuple = (1, 0.5),
                        class_fail::Symbol = :pop_cpd, 
                        chain_fail::Symbol = :ignore
=#


# One line match
#match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule) = match_rules(project, analyte, prec, ms1, mode, rule, rule.criteria)

# Mutiline match
function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Union{RuleUnion, RuleMode})
    ions = Union{Tuple, AbstractRule}[]
    done = true
    for r in rule.rules
        result = match_rules(project, analyte, prec, ms1, mode, r)
        done = done || result.matched
        @match result begin
            Result(false, x::AbstractRule)                                  => push!(ions, result.rule)
            Result(true, cls::T where T <: ClassSP) && if hasisomer(T) end  => push!(ions, result.rule.isomer)
            Result(true, ::EmptyRule)                                       => continue
            _                                                               => push!(ions, tuplize(result.rule))
        end
    end
    @match rule begin
        ::RuleUnion  => begin
            ions = union((), ions...)
            if done 
                @match ions begin
                    []                          => match_rules(project, analyte, prec, ms1, mode, rule.exception)
                    [x::T where T <: ClassSP]   => Result(true, T())
                    [x::T where T, xs...]       => Result(true, deisomerized(T)(ions))
                end
            else
                Result(false, RuleUnion(ions; exception = rule.exception))
            end
        end
        ::RuleMode   => done ? Result(false, RuleMode(ions)) : Result(true, mode(ions))
    end
end

# match result
#match_rules(project::Project, analyte::AnalyteSP, precursor::Vector, ms1::Vector, mode::Symbol, rule::Nothing) = (true, rule)
match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule) = Result(true, rule)
match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Tuple) = Result(true, rule)
match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Tuple{}) = Result(true, EmptyRule())
match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::EmptyRule) = Result(true, rule)

# Y/N ==
function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule{<: Union{Ion, ISF}})
    found = @match (rule.criteria, mode) begin
        (criteria::ISF, _)  => any(equivalent_in(criteria, @view cpd.fragments[!, :ion1]) for cpd in analyte)
        (criteria, :isf)    => any(equivalent_in(criteria, @views filter(:ion1 => (x -> equivalent_in(x, prec)), cpd.fragments, view = true)[!, :ion2]) for cpd in analyte)
        (criteria, n::Int)  => equivalent_in(criteria, @views filter(:ion1 => (x -> equivalent_in(x, prec)), analyte[n].fragments, view = true)[!, :ion2])
        (criteria, _)       => equivalent_in(criteria, @views filter(:ion1 => (x -> equivalent_in(x, prec)), last(analyte).fragments, view = true)[!, :ion2])
    end

    cpd = @match mode begin
        n::Int  => analyte[n]
        _       => last(analyte)
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.rule))
    elseif isa(rule.criteria, ISF) ? isf_checker(project, cpd, (rule.criteria,)) : ion_checker(project, cpd, prec, ms1, rule.criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.rule))
    else
        Result(false, rule)
    end
end

# Y/N in
function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule{<: Tuple})
    found = @match (first(rule.criteria), mode) begin
        (::ISF, _)  => any(any(equivalent_in(frag, rule.criteria) for frag in @view cpd.fragments[!, :ion1]) for cpd in analyte)
        (_, :isf)   => any(any(equivalent_in(frag, rule.criteria) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)),  cpd.fragments, view = true)[!, :ion2]) for cpd in analyte)
        (_, n::Int) => any(equivalent_in(frag, rule.criteria) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)), analyte[n].fragments, view = true)[!, :ion2])
        _           => any(equivalent_in(frag, rule.criteria) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)), last(analyte).fragments, view = true)[!, :ion2])
    end
    
    cpd = @match mode begin
        n::Int  => analyte[n]
        _       => last(analyte)
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.rule))
    elseif isa(first(rule.criteria), ISF) ? isf_checker(project, cpd, rule.criteria) : ion_checker(project, cpd, prec, ms1, rule.criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.rule))
    else
        Result(false, rule)
    end
end

function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule{<: AcylIon})    
    found = @match mode begin
        :isf    => any(any(equivalent_in(frag, rule.criteria.ions) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)), cpd.fragments, view = true)[!, :ion2]) for cpd in analyte)
        n::Int  => any(equivalent_in(frag, rule.criteria.ions) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)), analyte[n].fragments, view = true)[!, :ion2])
        _       => any(equivalent_in(frag, rule.criteria.ions) for frag in @views filter(:ion1 => (x -> equivalent_in(x, prec)), last(analyte).fragments, view = true)[!, :ion2])
    end

    cpd = @match mode begin
        n::Int  => analyte[n]
        _       => last(analyte)
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.rule))
    elseif ion_checker(project, cpd, prec, ms1, rule.criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.rule))
    else
        Result(false, rule)
    end
end
# Comparison
function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule{<: IonComparison})
    ratios = ion_comparison(project, analyte, prec, ms1, mode, rule.criteria)
    if !isempty(ratios)
        for ratio in Iterators.reverse(ratios)
            for level in rule.rule
                result = @match level.first begin
                    n::Tuple                => between(ratio, n)
                    x && if isnan(x) end    => isnan(ratio)
                    n::Number               => ==(ratio. n)
                end 
                result && return match_rules(project, analyte, prec, ms1, mode, level.second)
            end
        end
        Result(true, EmptyRule())
    else
        Result(false, rule)
    end
end

# hydroxl
function match_rules(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::Rule{<: Hydroxyl})
    cpd = @match mode begin
        n::Int  => analyte[n]
        _       => last(analyte)
    end
    hydroxyl = cpd.sum[3] - nhydroxyl(rule.criteria.spb)
    for level in rule.rule
        if hydroxyl == level.first
            result = match_rules(project, analyte, prec, ms1, mode, level.second)
            return result.matched ? result : cpd.sum[2] < ndb(rule.criteria.spb) ? Result(true, EmptyRule()) : PartialResult(false, result.rule, Chain(rule.criteria.spb, Acyl{hydroxyl}()))
        end
    end
    return Result(true, EmptyRule())
end

# ===========================================================================
function generate_ms(cpd::CompoundSP, mode::Symbol, anion::Symbol)
    add = @match mode begin
        :isf        => return generate_ms_isf(cpd, anion)
        :isf_sep    => return generate_ms_isf_sep(cpd, anion)
        :adduct     => class_db_index(cpd.class).adduct_ion
        :parent     => class_db_index(cpd.class).parent_ion
        :default    => class_db_index(cpd.class).default_ion
    end
    if anion == :acetate
        add = filter(!=(AddHCOO()), add)
    elseif anion == :formate
        add = filter(!=(AddOAc()), add)
    end
    Ion.(add, Ref(cpd.class)), mz.(Ref(cpd), add)
end

function generate_ms_isf(cpd::CompoundSP, anion::Symbol)
    all_isf = isf(cpd.class)
    if anion == :acetate
        fn = x -> filter(!=(AddHCOO()), x)
    elseif anion == :formate
        fn = x -> filter(!=(AddOAc()), x)
    end
    ions = mapreduce(vcat, all_isf) do isf_class
        map(fn(class_db_index(isf_class).isf_ion)) do add
            Ion(add, isf_class)
        end
    end

    ions, mz.(Ref(cpd), ions)
end

function generate_ms_isf_sep(cpd::CompoundSP, anion::Symbol)
    all_isf = isf(cpd.class)
    if anion == :acetate
        fn = x -> filter(!=(AddHCOO()), x)
    elseif anion == :formate
        fn = x -> filter(!=(AddOAc()), x)
    end
    ions = map(all_isf) do isf_class
        isf_class => map(fn(class_db_index(isf_class).isf_ion)) do add
            Ion(add, isf_class)
        end
    end

    ions, map(subions -> subions.first => mz.(Ref(cpd), subions.second), ions)
end
#isf_checker(project::Project, analyte::AnalyteSP, left::Tuple) = true

function isf_checker(project::Project, cpd::CompoundSP, criteria)
    polarity = isa(first(criteria).adduct, Pos)
    ms1 = [mz(cpd, ion) for ion in criteria]
    any(project.data) do data
        polarity == data.polarity || return false
        if isnothing(cpd.chain)
            raw_id = eachindex(data.mz2)
        else
            ###
            raw_id = findall(mz2 -> abs(rem(mz(Ion(Protonation(), cpd.chain.lcb)) - mz2, mw("H2O"), RoundNearest)) < data.mz_tol, data.mz2)
            isempty(raw_id) && return false
        end
        @match data begin
            ::PreIS => any(any(between(m1, data.range[id]) for m1 in ms1) for id in raw_id)
            ::MRM   => any(any(between(m1, mz1, data.mz_tol) for m1 in ms1) for mz1 in filter(:id => in(raw_id), data.raw, view = true).mz1)
        end
    end
end

#ion_checker(project::Project, precursor::Vector, ms1::Vector, rule::ISF) = true # Assume covered by PreIS

function ion_checker(project::Project, cpd::CompoundSP, prec::Vector, ms1::Vector, rule::Ion{T}) where T
    polarity = @match T begin
        ::Type{<: Pos} => Pos
        ::Type{<: Neg} => Neg
    end
    id = findall(ion -> isa(ion.adduct, polarity), prec)
    ms1 = @view ms1[id]
    ms2 = mz(cpd, rule)
    polarity = ==(polarity, Pos)
    any(project.data) do data
        polarity == data.polarity && begin
            raw_id = findfirst(mz2 -> between(ms2, mz2, data.mz_tol), data.mz2)
            isnothing(raw_id) && return false
            @match data begin
                ::PreIS => any(any(between(m1, range) for m1 in ms1) for range in @view data.range[raw_id])
                ::MRM   => any(any(between(m1, mz1, data.mz_tol) for m1 in ms1) for mz1 in @views data.raw[data.raw.mz2 .== raw_id, :mz])
            end
        end
    end
end

function ion_checker(project::Project, cpd::CompoundSP, prec::Vector, ms1::Vector, rule::Tuple)
    length(rule) == 1 && return ion_checker(project, cpd, prec, ms1, rule)
    id1 = findall(ion -> isa(ion.adduct, Pos), prec)
    id2 = setdiff(eachindex(ms1), id1)
    ms1_pos = isempty(id1) ? nothing : @view ms1[id1]
    ms1_neg = isempty(id2) ? nothing : @view ms1[id2]
    ms2 = mz.(Ref(cpd), rule)
    id1 = findall(ion -> isa(ion.adduct, Pos), rule)
    id2 = setdiff(eachindex(ms2), id1)
    ms2_pos = isempty(id1) ? nothing : ms2[id1]
    ms2_neg = isempty(id2) ? nothing : ms2[id2]
    any(project.data) do data
        ms1, ms2 = data.polarity ? (ms1_pos, ms2_pos) : (ms1_neg, ms2_neg)
        raw_id = findfirst(mz2 -> any(between(m2, mz2, data.mz_tol) for m2 in ms2), data.mz2)
        isnothing(raw_id) && return false
        @match data begin
            ::PreIS => any(any(between(m1, range) for m1 in ms1) for range in @view data.range[raw_id])
            ::MRM   => any(any(between(m1, mz1, data.mz_tol) for m1 in ms1) for mz1 in @views data.raw[data.raw.mz2 .== raw_id, :mz])
        end
    end
end

function ion_checker(project::Project, cpd::CompoundSP, prec::Vector, ms1::Vector, rule::AcylIon)
    id = findall(ion -> isa(ion.adduct, Neg), prec)
    ms1 = isempty(id) ? nothing : @view ms1[id]
    ms2 = mz(cpd, rule)
    any(project.data) do data
        data.polarity && return false
        raw_id = findfirst(mz2 -> any(between(m2, mz2, data.mz_tol) for m2 in ms2), data.mz2)
        isnothing(raw_id) && return false
        @match data begin
            ::PreIS => any(any(between(m1, range) for m1 in ms1) for range in @view data.range[raw_id])
            ::MRM   => any(any(between(m1, mz1, data.mz_tol) for m1 in ms1) for mz1 in @views data.raw[data.raw.mz2 .== rwa_id, :mz])
        end
    end
end

function ion_comparison(project::Project, analyte::AnalyteSP, prec::Vector, ms1::Vector, mode, rule::IonComparison{S, T}) where {S, T}
    all_ions = map(rule.ions) do ions
        @match ions begin
            IonPlus(_)  => map(i -> (i, mz(i)), ions.ions)
            _           => [(ions, mz(ions))]
        end
    end
    if mode == :isf
        ratios = Float64[]
        for cpd in analyte
            gdf = groupby(filter(:ion1 => (x -> equivalent_in(x, prec)), cpd.fragments, view = true), [:source, :ion1])
            for subdf in gdf
                done = true
                ion1 = subdf[1, :ion1]
                ms1_filter = ms1[findfirst(==(ion1), prec)]
                source = subdf[1, :source]
                ratio = map(all_ions) do ions
                    map(ions) do ion
                        id = findfirst(x -> equivalent(ion[1], x), @view subdf[!, :ion2])
                        if isnothing(id)
                            ion_checker(project.data[source], ms1_filter, ion[2]) || (done = false)
                            0
                        else
                            query_raw(project, source, subdf[id, :id]).area
                        end
                    end
                end
                done && push!(ratios, mapreduce(sum, /, ratio))
            end
        end
    else
        cpd = @match mode begin
            n::Int  => analyte[n]
            _       => last(analyte)
        end
        gdf = groupby(filter(:ion1 => (x -> equivalent_in(x, prec)), cpd.fragments, view = true), [:source, :ion1])
        ratios = Float64[]
        for subdf in gdf
            done = true
            ion1 = subdf[1, :ion1]
            ms1_filter = ms1[findfirst(==(ion1), prec)]
            source = subdf[1, :source]
            ratio = map(all_ions) do ions
                map(ions) do ion
                    id = findfirst(x -> equivalent(ion[1], x), @view subdf[!, :ion2])
                    if isnothing(id)
                        ion_checker(project.data[source], ms1_filter, ion[2]) || (done = false)
                        0
                    else
                        query_raw(project, source, subdf[id, :id]).area
                    end
                end
            end
            done && push!(ratios, mapreduce(sum, /, ratio))
        end
    end
    ratios
end

function ion_checker(data::Data, ms1::Float64, ms2::Float64)
    @match data begin
        ::PreIS => any(between(mz2, ms2, data.mz_tol) && between(ms1, range) for (mz2, range) in zip(data.mz2, data.range))
        ::MRM   => any(between(mz2, ms2, data.mz_tol) && any(between(mz1, ms1, data.mz_tol) for mz1 in filter(:mz2 => ==(i), data.raw, view = true)[!, :mz]) for (i, mz2) in enumerate(data.mz2))
    end
end

#=
function ion_checker(project::Project, precursor::Vector, ms1::Vector, rule::IonComparison)
    polarity = @match rule.ions[1] begin
        IonPlus((Ion(::Pos, _), _)) => true
        IonPlus                     => false
        Ion(::Pos, _)               => true
        Ion(::Neg, _)               => false
    end
    id = findall(ion -> isa(ion.adduct, Pos) == polarity, precursor)
    ms1 = ms1[id]
    ions = map(rule.ions) do ion
        @match ion begin
            ::IonPlus   => collect(mw.(ion.ions))
            _           => mz(ion)
        end
    end
    ms2 = vcat(ions...)
    any(project.data) do data
        polarity == data.polarity && begin
            ids_prod = map(ions) do ms2
                findfirst(mz2 -> between(ms2, mz2, data.mz_tol), data.products)
            end
            any(isnothing, ids_prod) && return false
            @match data begin
                ::PreIS => any(between(ms, intersection(@view data.range[ids_prod]) for ms in ms1))
                ::MRM   => !isempty(
                    mapreduce(intersect, ids_prod) do id_prod
                        [ms for mz1 in @views data.precursors[data.precursors.mz2 .== id_prod, :mz] for ms in ms1 if between(ms, mz1, data.mz_tol)]
                    end
                )
            end
        end
    end
end =#