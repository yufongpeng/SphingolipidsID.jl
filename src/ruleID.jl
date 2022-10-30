# ==================================================================================
# Rule-based ID
# start point
function apply_rules!(project::Project, part::Symbol = :both; class_mode::Symbol = :default, chain_mode::Symbol = :isf, class_rule_mode::Symbol = :ab)
    del = Int[]
    @match part begin
        :both   => println("ID> Class and Chain")
        :class  => println("ID> Class only")
        :chain  => println("ID> Chain only")
    end

    for (i, analyte) in enumerate(project.analytes)
        cpd = last(analyte.identification)
        if part == :both || part == :class 
            id_result_class = match_rules(project, analyte, generate_ms(cpd, class_mode, project.anion)..., class_mode,
            ismissing(first(cpd.sub_rule)) ? rule(cpd.class, class_rule_mode) : first(cpd.sub_rule)
            )
        else
            id_result_class = nothing
        end

        if part == :both || part == :chain
            r = @match (last(last(analyte.identification).sub_rule), last(analyte.identification).chain, chain_mode) begin
                (missing, nothing, :isf)    => begin
                    id = findfirst(cpd -> !isnothing(cpd.chain), reverse(analyte.identification))
                    isnothing(id) ? id : rule(analyte.identification[end - id + 1].chain)
                end
                (missing, nothing, _)       => nothing
                (missing, chain, _)         => rule(chain)
                (r, _, _)                   => r
            end
            isnothing(r) && continue
            id_result_chain = match_rules(project, analyte, generate_ms(last(analyte.identification), chain_mode, project.anion)..., chain_mode, r)
        else
            id_result_chain = nothing
        end

        if !isnothing(id_result_class) 
            if first(id_result_class)
                if last(id_result_class) === ()
                    push!(del, i)
                else
                    last(analyte.identification).sub_rule[1] = nothing
                    last(analyte.identification).class = last(id_result_class)
                end
            else
                last(analyte.identification).sub_rule[1] = last(id_result_class) 
            end
        end

        if !isnothing(id_result_chain)
            if first(id_result_chain)
                if last(id_result_chain) === ()
                    push!(del, i)
                else
                    last(analyte.identification).sub_rule[2] = nothing
                    last(analyte.identification).chain = last(id_result_chain)
                end
            else
                last(analyte.identification).sub_rule[2] = last(id_result_chain) 
            end
        end
    end
    unique!(del)
    deleteat!(project.analytes, del)
    project
end

# One line match
match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule::Pair) = match_rules(project, analyte, prec, ms1, mode, rule, rule.first)

function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule::Union{IonUnion, IonMode})
    ions = Tuple[]
    done = true
    for r in rule.rules
        result = match_rules(project, analyte, prec, ms1, mode, r)
        done = done || first(result)
        push!(ions, tuplize(last(result)))
    end
    @match rule begin
        ::IonUnion  => begin
            ions = union(ions)
            done ? (false, IonUnion(ions)) : (true, typeof(deisomerized(ions[1]))(ions))
        end
        ::IonMode   => done ? (false, IonMode(ions)) : (true, mode(ions))
    end
end

# match result
#match_rules(project::Project, analyte::AnalyteGSL, precursor::Vector, ms1::Vector, mode::Symbol, rule::Nothing) = (true, rule)
match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule) = true, rule
match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule::Tuple) = true, rule

# Y/N ==
function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule, criteria::Union{Ion, ISF})
    found = @match (criteria, mode) begin
        (::ISF, _)  => any(equivalent_in(criteria, @view cpd.fragments[!, :ion1]) for cpd in analyte.identification)
        (_, :isf)   => any(equivalent_in(criteria, @views filter(:ion1 => in(prec), cpd.fragments, view = true)[!, :ion2]) for cpd in analyte.identification)
        _           => equivalent_in(criteria, @views filter(:ion1 => in(prec), last(analyte.identification).fragments, view = true)[!, :ion2])
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.second))
    elseif isa(criteria, ISF) ? isf_checker(project, analyte, (criteria,)) : ion_checker(project, analyte, prec, ms1, criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.second))
    else
        false, rule
    end
end

# Y/N in
function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule, criteria::Tuple)
    found = @match (first(criteria), mode) begin
        (::ISF, _)  => any(any(equivalent_in(frag, criteria) for frag in @view cpd.fragments[!, :ion1]) for cpd in analyte.identification)
        (_, :isf)   => any(any(equivalent_in(frag, criteria) for frag in @views filter(:ion1 => in(prec),  cpd.fragments, view = true)[!, :ion2]) for cpd in analyte.identification)
        _           => any(equivalent_in(frag, criteria) for frag in @views filter(:ion1 => in(prec), last(analyte.identification).fragments, view = true)[!, :ion2])
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.second))
    elseif isa(first(criteria), ISF) ? isf_checker(project, analyte, criteria) : ion_checker(project, analyte, prec, ms1, criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.second))
    else
        false, rule
    end
end

function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule, criteria::AcylIon)    
    found = @match mode begin
        :isf    => any(any(equivalent_in(frag, criteria.ions) for frag in  @views filter(:ion1 => in(prec), cpd.fragments, view = true)[!, :ion2]) for cpd in analyte.identification)
        _       => any(equivalent_in(frag, criteria.ions) for frag in @views filter(:ion1 => in(prec), last(analyte.identification).fragments, view = true)[!, :ion2])
    end

    if found
        match_rules(project, analyte, prec, ms1, mode, first(rule.second))
    elseif ion_checker(project, analyte, prec, ms1, criteria)
        match_rules(project, analyte, prec, ms1, mode, last(rule.second))
    else
        false, rule
    end
end
# Comparison
function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule, criteria::IonComparison)
    ratios = ion_comparison(project, analyte, prec, ms1, mode, criteria)
    if !isempty(ratios)
        for ratio in Iterators.reverse(ratios)
            for level in rule.second
                result = @match level.first begin
                    n::Tuple                => between(ratio, n)
                    x && if isnan(x) end    => isnan(ratio)
                    n::Number               => ==(ratio. n)
                end 
                result && return match_rules(project, analyte, prec, ms1, mode, level.second)
            end
        end
        true, ()
    else
        false, rule
    end
end

# hydroxl
function match_rules(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule, criteria::Hydroxyl)
    hydroxyl = last(analyte.identification).sum[3] - nhydroxyl(criteria.spb)
    for level in rule.second
        hydroxyl == level.first && return match_rules(project, analyte, prec, ms1, mode, level.second)
    end
    return true, ()
end

# ===========================================================================
function generate_ms(cpd::CompoundGSL, mode::Symbol, anion::Symbol)
    add = @match mode begin
        :isf        => return generate_ms_isf(cpd, anion)
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

function generate_ms_isf(cpd::CompoundGSL, anion::Symbol)
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
#isf_checker(project::Project, analyte::AnalyteGSL, left::Tuple) = true

function isf_checker(project::Project, analyte::AnalyteGSL, criteria)
    polarity = isa(first(criteria).adduct, Pos)
    cpd = last(analyte.identification)
    ms1 = [mz(cpd, ion) for ion in criteria]
    any(project.data) do data
        polarity == data.polarity || return false
        if isnothing(cpd.chain)
            raw_id = eachindex(data.mz2)
        else
            ###
            raw_id = findall(mz2 -> abs(rem(mz(Ion(Protonation(), cpd.chain.lcb)) - mz2, mw("H2O"))) < data.mz_tol, data.mz2)
            isnothing(raw_id) && return false
        end
        @match data begin
            ::PreIS => any(any(between(m1, data.range[id]) for m1 in ms1) for id in raw_id)
            ::MRM   => any(any(between(m1, mz1, data.mz_tol) for m1 in ms1) for mz1 in filter(:id => in(raw_id), data.raw, view = true).mz1)
        end
    end
end

#ion_checker(project::Project, precursor::Vector, ms1::Vector, rule::ISF) = true # Assume covered by PreIS

function ion_checker(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, rule::Ion{T}) where T
    #cpd = last(analyte.identification)
    polarity = @match T begin
        ::Type{<: Pos} => Pos
        ::Type{<: Neg} => Neg
    end
    id = findall(ion -> isa(ion.adduct, polarity), prec)
    ms1 = @view ms1[id]
    ms2 = mz(rule)
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

function ion_checker(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, rule::Tuple)
    #cpd = last(analyte.identification)
    length(rule) == 1 && return ion_checker(project, analyte, prec, ms1, rule)
    id1 = findall(ion -> isa(ion.adduct, Pos), prec)
    id2 = setdiff(eachindex(ms1), id1)
    ms1_pos = isempty(id1) ? nothing : @view ms1[id1]
    ms1_neg = isempty(id2) ? nothing : @view ms1[id2]
    ms2 = mz.(rule)
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

function ion_checker(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, rule::AcylIon)
    cpd = last(analyte.identification)
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

function ion_comparison(project::Project, analyte::AnalyteGSL, prec::Vector, ms1::Vector, mode::Symbol, rule::IonComparison{S, T}) where {S, T}
    all_ions = map(rule.ions) do ions
        @match ions begin
            IonPlus(_)  => map(i -> (i, mz(i)), ions.ions)
            _           => [(ions, mz(ions))]
        end
    end
    if mode == :isf
        ratios = Float64[]
        for cpd in analyte.identification
            gdf = groupby(filter(:ion1 => in(prec), cpd.fragments, view = true), [:source, :ion1])
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
        cpd = last(analyte.identification)
        gdf = groupby(filter(:ion1 => in(prec), cpd.fragments, view = true), [:source, :ion1])
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