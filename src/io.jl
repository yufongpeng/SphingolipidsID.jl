
function Base.show(io::IO, class::T) where {T <: ClassGSL}
    hasisomer(T) ? begin
        str = join(map(repr, class.isomer), ", ")
        str = isempty(str) ? "all isomers" : str 
        print(io, repr(T), ": ", str) 
    end : print(io, T)
end

Base.show(io::IO, ion::Ion) = print(io, repr_adduct(ion.adduct), " of ", ion.molecule)
Base.show(io::IO, lcb::LCB) = print(io, "SPB ", ncb(lcb), ":", ndb(lcb), ";", nhydroxyl(lcb), "O")
Base.show(io::IO, acyl::Acyl{N}) where N = print(io, "Acyl (", N, " OH" , ")")
Base.show(io::IO, acyl::Acylα{N}) where N = print(io, "Acyl (", 1, " α-OH" , N - 1, " other OH", ")")
Base.show(io::IO, acyl::Acylβ{N}) where N = print(io, "Acyl (", 1, " β-OH" , N - 1, " other OH", ")")

function Base.show(io::IO, cpd::CompoundGSL)
    println(io, cpd.class, " ", cpd.sum[1], ":", cpd.sum[2], ";", cpd.sum[3], "O")
    isnothing(cpd.chain) && return
    println(io, "  ", cpd.chain.lcb)
    print(io, "  ", cpd.chain.acyl)
end

function Base.show(io::IO, analyte::AnalyteGSL)
    println(io, "Compounds at ", round(analyte.rt, digits = 2), ":")
    for cpd in analyte.identification
        println(io, " ", cpd)
    end
end

function Base.show(io::IO, data::PreIS)
    println(io, "PreIS: ")
    for dt in zip(data.range, data.mz2)
        println(io, " ", dt[1][1], " - ", dt[1][2], " -> ", dt[2])
    end
end

function Base.show(io::IO, pj::Project)
    println(io, "Project with ", length(pj.analytes), " analytes:")
    println(io, "∘ Salt: ", pj.anion)
    println(io, "∘ Data: ")
    for data in pj.data
        show(io, data)
    end
    println(io, "∘ Analytes: ")
    if length(pj.analytes) > 10
        for analyte in @view pj.analytes[1:5]
            println(io, analyte)
        end
        println(io, "\t", ".")
        println(io, "\t", ".")
        println(io, "\t", ".")
        for analyte in @view pj.analytes[end - 4:end]
            println(io, analyte)
        end
    else
        for analyte in pj.analytes
            println(io, analyte)
        end
    end
end

function Base.show(io::IO, aquery::Query)
    print(io, "Query: ")
    for q in aquery.query
        print(io, q, ", ")
    end
    println(io)
    print(io, "Result")
    aquery.view ? println(io, "(view):") : println(io, ":")
    if length(aquery.result) > 10
        for r in @view aquery.result[1:5]
            println(io, r)
        end
        println(io, "\t", ".")
        println(io, "\t", ".")
        println(io, "\t", ".")
        for r in @view aquery.result[end - 4:end]
            println(io, r)
        end
    else
        for r in aquery.result
            println(io, r)
        end
    end
end