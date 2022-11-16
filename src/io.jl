
function Base.show(io::IO, class::T) where {T <: ClassGSL}
    hasisomer(T) ? begin
        str = join(map(x -> repr, class.isomer), ", ")
        str = isempty(str) ? "" : str 
        print(io, repr(T), isempty(str) ? "" : "($str)") 
    end : print(io, T)
end

Base.show(io::IO, ion::Ion) = print(io, repr_adduct(ion.adduct), " of ", ion.molecule)
Base.show(io::IO, sugar::T) where {T <: Sugar} = print(io, T)
Base.show(io::IO, glycan::Glycan{T}) where {T <: Tuple} = print(io, join(T.parameters, "_"))
Base.show(io::IO, lcb::LCB{N}) where N = print(io, "SPB ", ncb(lcb), ":", ndb(lcb), ";", N > 1 ? nhydroxyl(lcb) : "", "O")
Base.show(io::IO, acyl::Acyl{0}) = print(io, "Acyl")
Base.show(io::IO, acyl::Acyl{N}) where N = print(io, "Acyl (", N, " OH" , ")")
Base.show(io::IO, acyl::Acylα{N}) where N = print(io, "Acyl (", 1, " α-OH" , N - 1, " other OH", ")")
Base.show(io::IO, acyl::Acylβ{N}) where N = print(io, "Acyl (", 1, " β-OH" , N - 1, " other OH", ")")
repr_hydroxl(acyl::Acyl{0}) = ""
repr_hydroxl(acyl::Acyl{1}) = ";O"
repr_hydroxl(acyl::Acyl{N}) where N = ";$(N)O"
repr_hydroxl(acyl::Acylα{1}) = ";(2OH)"
repr_hydroxl(acyl::Acylα{N}) where N = ";(2OH);$(N - 1)O"
repr_hydroxl(acyl::Acylβ{1}) = ";(3OH)"
repr_hydroxl(acyl::Acylβ{N}) where N = ";(3OH);$(N - 1)O"

fragment_table(cpd) = map(eachrow(cpd.fragments)) do row
    s = query_raw(cpd.project, row.source, row.id)
    mz2 = in(:mz2, propertynames(s)) ? s.mz2 : cpd.project.data[row.source].mz2[s.scan]
    (ion1 = row.ion1, mz1 = s.mz1, ion2 = row.ion2, mz2 = mz2, area = s.area, CE = s.collision_energy, rt = s.rt)
end

function Base.show(io::IO, ::MIME"text/plain", cpd::CompoundGSL)
    class = cpd.states[1] == 1 ? "🟢" : cpd.states[1] == -1 ? "🔴" : "🟡" 
    chain = cpd.states[2] == 1 ? "🟢"  : cpd.states[2] == -1 ? "🔴" : "🟡"
    print(io, "Compound with ", size(cpd.fragments, 1), " fragments ($class,$chain):")
    print(io, "\n∘ ID: ", cpd)
    print(io, "\n∘ Area: ", cpd.area)
    dt = fragment_table(cpd)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "CE (eV)", "RT (min)"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r])
end

function Base.show(io::IO, cpd::CompoundGSL)
    if isnothing(cpd.chain)
        print(io, cpd.class, " ", cpd.sum[1], ":", cpd.sum[2], ";", cpd.sum[3], "O")
    else
        spb_c, spb_db, spb_o = sumcomp(cpd.chain.lcb)
        print(io, cpd.class, " ", spb_c, ":", spb_db, ";", spb_o > 1 ? spb_o : "", "O/")
        print(io, cpd.sum[1] - spb_c, ":", cpd.sum[2] - spb_db, repr_hydroxl(cpd.chain.acyl))
    end
end

function Base.show(io::IO, ::MIME"text/plain", analyte::AnalyteGSL)
    class = analyte.states[1] == 1 ? "🟢" : analyte.states[1] == -1 ? "🔴" : "🟡" 
    chain = analyte.states[2] == 1 ? "🟢"  : analyte.states[2] == -1 ? "🔴" : "🟡" 
    print(io, "Analytes with ", length(analyte), "compounds @", round(analyte.rt, digits = 2), " ($class,$chain):")
    print(io, "\n∘ Compounds:")
    for cpd in analyte
        print(io, "\n ", cpd)
    end
    dt = mapreduce(fragment_table, vcat, analyte)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "CE (eV)", "RT (min)"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r])
end

function Base.show(io::IO, analyte::AnalyteGSL)
    class = analyte.states[1] == 1 ? "🟢" : analyte.states[1] == -1 ? "🔴" : "🟡" 
    chain = analyte.states[2] == 1 ? "🟢"  : analyte.states[2] == -1 ? "🔴" : "🟡" 
    print(io, last(analyte), " @", round(analyte.rt, digits = 2), " ($class,$chain)")
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
    print(io, "Not(", pred.args, ")")
end