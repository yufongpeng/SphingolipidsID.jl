
function Base.show(io::IO, ::MIME"text/plain", class::T) where {T <: ClassSP}
    hasisomer(T) ? begin
        str = join(map(repr, class.isomer), ", ")
        str = isempty(str) ? "" : str 
        print(io, replace(repr(T), "_" => "?"), isempty(str) ? "" : "($str)") 
    end : print(io, T)
end

Base.show(io::IO, class::T) where {T <: ClassSP} = print(io, replace(repr(T), r"_$" => ""))

Base.show(io::IO, ion::Ion) = print(io, repr_adduct(ion.adduct), " of ", ion.molecule)
Base.show(io::IO, sugar::T) where {T <: Sugar} = print(io, T)
Base.show(io::IO, glycan::Glycan{T}) where {T <: Tuple} = print(io, join(T.parameters, "_"))
Base.show(io::IO, lcb::LCB{N}) where N = print(io, "SPB ", ncb(lcb), ":", ndb(lcb), ";", "O", N > 1 ? nhydroxyl(lcb) : "")
Base.show(io::IO, lcb::Lcb) = print(io, "SPB ", lcb.cb, ":", lcb.db, ";", "O", lcb.ox > 1 ? lcb.ox : "")

Base.show(io::IO, acyl::NACYL) = print(io, "Acyl ", acyl.cb, ":", acyl.db, repr_hydroxl(acyl))

Base.show(io::IO, acyl::Acyl{0}) = print(io, "Acyl")
Base.show(io::IO, acyl::Acyl{N}) where N = print(io, "Acyl (", N, " OH" , ")")
Base.show(io::IO, acyl::Acylα{N}) where N = print(io, "Acyl (", 1, " α-OH" , N - 1, " other OH", ")")
Base.show(io::IO, acyl::Acylβ{N}) where N = print(io, "Acyl (", 1, " β-OH" , N - 1, " other OH", ")")

repr_hydroxl(acyl::Nacyl) = acyl.ox == 0 ? "" : acyl.ox == 1 ? ";O" : ";O$(acyl.ox)"
repr_hydroxl(acyl::Nacylα) = acyl.ox == 0 ? "" : acyl.ox == 1 ? ";(2OH)" : ";(2OH);O$(acyl.ox - 1)"
repr_hydroxl(acyl::Nacylβ) = acyl.ox == 0 ? "" : acyl.ox == 1 ? ";(3OH)" : ";(3OH);O$(acyl.ox - 1)"

repr_hydroxl(acyl::Acyl{0}) = ""
repr_hydroxl(acyl::Acyl{1}) = ";O"
repr_hydroxl(acyl::Acyl{N}) where N = ";O$N"
repr_hydroxl(acyl::Acylα{1}) = ";(2OH)"
repr_hydroxl(acyl::Acylα{N}) where N = ";(2OH);O$(N - 1)"
repr_hydroxl(acyl::Acylβ{1}) = ";(3OH)"
repr_hydroxl(acyl::Acylβ{N}) where N = ";(3OH);O$(N - 1)"

function Base.show(io::IO, score::Score)
    println(io, score.score)
end

function Base.show(io::IO, ::MIME"text/plain", score::Score)
    println(io, "Score:")
    println(io, "∘ Current score: ", score.score)
    println(io, "∘ Target: ", score.parameters.target)
    println(io, "∘ Converter: ", score.parameters.converter)
    println(io, "∘ Weight: ", score.parameters.weight)
    println(io, "∘ Objective: ", score.parameters.objective)
    println(io, "∘ Threshold: ", score.parameters.threshold)
end

fragment_table(cpd) = map(cpd.fragments) do row
    s = query_raw(cpd.project, row.source, row.id)
    mz2 = in(:mz2, propertynames(s)) ? cpd.project.data[row.source].mz2[s.mz2] : cpd.project.data[row.source].mz2[s.scan]
    mode = isa(cpd.project.data[row.source], PreIS) ? "PreIS" : "MRM"
    (ion1 = row.ion1, mz1 = s.mz1, ion2 = row.ion2, mz2 = mz2, area = s.area, error = s.error, CE = s.collision_energy, rt = s.rt, source = mode)
end

show_cpdid_class(c::Type{T}) where {T <: ClassSP}= "$T "
show_cpdid_class(c::Type{Nothing}) = ""

function Base.show(io::IO, cpd::CompoundID{C}) where C
    if isnothing(cpd.lcb)
        print(io, show_cpdid_class(C), cpd.sum[1], ":", cpd.sum[2], ";", "O", cpd.sum[3] > 1 ? cpd.sum[3] : "")
    else
        print(io, show_cpdid_class(C), cpd.lcb.cb, ":", cpd.lcb.db, ";", "O", cpd.lcb.ox > 1 ? cpd.lcb.ox : "", "/")
        print(io, cpd.acyl.cb, ":", cpd.acyl.db, repr_hydroxl(cpd.acyl))
    end
end

function Base.show(io::IO, ::MIME"text/plain", cpd::CompoundSP)
    class = cpd.states[1] == 1 ? "🟢" : cpd.states[1] == -1 ? "🔴" : "🟡" 
    chain = cpd.states[2] == 1 ? "🟢"  : cpd.states[2] == -1 ? "🔴" : "🟡"
    print(io, "Compound with ", size(cpd.fragments, 1), " fragments ($class,$chain):")
    print(io, "\n∘ ID: ")
    show(io, MIME"text/plain"(), cpd.class)
    if isnothing(cpd.chain)
        print(io, " ", cpd.sum[1], ":", cpd.sum[2], ";", "O", cpd.sum[3] > 1 ? cpd.sum[3] : "")
    else
        spb_c, spb_db, spb_o = sumcomp(cpd.chain.lcb)
        print(io, " ", spb_c, ":", spb_db, ";", "O", spb_o > 1 ? spb_o : "", "/")
        print(io, cpd.sum[1] - spb_c, ":", cpd.sum[2] - spb_db, repr_hydroxl(cpd.chain.acyl))
    end
    print(io, "\n∘ Area: ", cpd.area[1])
    print(io, "\n∘ Error: ", cpd.area[2])
    dt = fragment_table(cpd)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r])
end

function Base.show(io::IO, cpd::CompoundSP)
    if isnothing(cpd.chain)
        print(io, cpd.class, " ", cpd.sum[1], ":", cpd.sum[2], ";", "O", cpd.sum[3] > 1 ? cpd.sum[3] : "")
    else
        spb_c, spb_db, spb_o = sumcomp(cpd.chain.lcb)
        print(io, cpd.class, " ", spb_c, ":", spb_db, ";", "O", spb_o > 1 ? spb_o : "", "/")
        print(io, cpd.sum[1] - spb_c, ":", cpd.sum[2] - spb_db, repr_hydroxl(cpd.chain.acyl))
    end
end

function Base.show(io::IO, ::MIME"text/plain", analyte::AnalyteSP)
    class = analyte.states[1] == 1 ? "🟢" : analyte.states[1] == -1 ? "🔴" : "🟡" 
    chain = analyte.states[2] == 1 ? "🟢"  : analyte.states[2] == -1 ? "🔴" : "🟡" 
    print(io, "Analytes with ", length(analyte), " compounds @", round(analyte.rt, digits = 2), " ($class,$chain):")
    print(io, "\n∘ Score: ", analyte.scores)
    print(io, "\n∘ Compounds:")
    for cpd in analyte
        print(io, "\n ", cpd)
    end
    dt = mapreduce(fragment_table, vcat, analyte)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r])
end

function Base.show(io::IO, analyte::AnalyteSP)
    class = analyte.states[1] == 1 ? "🟢" : analyte.states[1] == -1 ? "🔴" : "🟡" 
    chain = analyte.states[2] == 1 ? "🟢"  : analyte.states[2] == -1 ? "🔴" : "🟡" 
    print(io, isempty(analyte.compounds) ? "?" : last(analyte), " @", round(analyte.rt, digits = 2), " ($class,$chain)")
end

function Base.show(io::IO, data::PreIS)
    println(io, "PreIS: ")
    for dt in zip(data.range, data.mz2)
        println(io, " ", dt[1][1], " ~ ", dt[1][2], " -> ", dt[2])
    end
end

function Base.show(io::IO, pj::Project)
    println(io, "Project with ", length(pj), " analytes:")
    println(io, "∘ Salt: ", pj.anion)
    println(io, "∘ Data: ")
    for data in pj.data
        show(io, data)
    end
    print(io, "∘ Analytes: ")
    if length(pj) > 10
        for analyte in @view pj[1:5]
            print(io, "\n ", analyte)
        end
        println(io, "\n\t", ".")
        println(io, "\t", ".")
        print(io, "\t", ".")
        for analyte in @view pj[end - 4:end]
            print(io, "\n ", analyte)
        end
    else
        for analyte in pj
            print(io, "\n ", analyte)
        end
    end
end

function Base.show(io::IO, aquery::Query)
    print(io, "Query with ", length(aquery), " analytes: \n")
    print(io, "∘ Queries: ")
    print(io, join(map(x -> replace(repr(x), "Any" => ""), aquery.query), ", "))
    println(io)
    print(io, "∘ Result")
    aquery.view ? print(io, " (view):") : print(io, ":")
    if length(aquery) > 10
        for r in @view aquery[1:5]
            print(io, "\n ", r)
        end
        println(io, "\n\t", ".")
        println(io, "\t", ".")
        print(io, "\t", ".")
        for r in @view aquery[end - 4:end]
            print(io, "\n ", r)
        end
    else
        for r in aquery
            print(io, "\n ", r)
        end
    end
end

function Base.show(io::IO, pred::Inv)
    print(io, "Not(", pred.arg, ")")
end