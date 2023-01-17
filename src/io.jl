function featuretable_mzmine(path)
    head = split(readline(path), ",")
    idmain = findall(x -> any(==(x, text) for text in ["id", "rt", "mz", "height", "area"]), head)
    idfwhm = findall(x-> endswith(x, "fwhm"), head)
    idsym = findall(x-> endswith(x, "asymmetry_factor"), head)
    tbl = CSV.read(path, Table; select = idmain)
    n = size(tbl, 1)
    fwhm = CSV.read(path, Table; select = idfwhm)
    sym = CSV.read(path, Table; select = idsym)
    datafile = Dict(propertynames(fwhm) .=> map(col -> match(r".*:(.*):.*", string(col))[1], propertynames(fwhm)))
    id = findfirst.(!ismissing, fwhm)
    tbl = Table(
        id = zeros(Int, n), 
        mz1 = tbl.mz, 
        mz2 = zeros(Float64, n), 
        rt = tbl.rt,
        height = tbl.height, 
        area = tbl.area,
        collision_energy = zeros(Int, n), 
        FWHM = getindex.(fwhm, id), 
        symmetry = get.(sym, replace!(findfirst.(!ismissing, sym), nothing => :_symmetry), 1.0), 
        datafile = getindex.(Ref(datafile), id)
    )
    tbl
end

function featuretable_masshunter_mrm(path)
    strs = readlines(path)
    starts = Int[]
    data = Table(ms1 = Float64[], ms2 = Float64[], eV = Int[])
    ends = Int[]
    status = false
    for (i, l) in enumerate(strs)
        if !status
            transition = match(r"\((\d*.\d*) -> (\d*.\d*)\)", l)
            isnothing(transition) && continue
            ms1, ms2 = parse.(Float64, transition)
            eV = parse(Float64, match(r"CID@(\d*.\d*)", l)[1])
            push!(data, (; ms1, ms2, eV))
            push!(starts, i + 1)
            status = true
        elseif l == ""
            push!(ends, i - 1)
            status = false
        elseif occursin("No integration results available", l)
            pop!(starts)
            pop!(data)
            status = false
        end
    end
    str = map(starts, ends) do st, ed
        IOBuffer(join(strs[st:ed], "\n"))
    end
    rep = ends .- starts
    txt = ["RT", "Height", "Area", "Symmetry", "FWHM"]
    tbl = CSV.read(str, Table; select = (i, name) -> any(==(text, String(name)) for text in txt))
    n = size(tbl, 1)
    Table(
        id = zeros(Int, n), 
        mz1 = (@p zip(data.ms1, rep) |> mapmany(repeat([_[1]], _[2]))),
        mz2 = (@p zip(data.ms2, rep) |> mapmany(repeat([_[1]], _[2]))),
        rt = tbl.RT,
        height = tbl.Height,
        area = tbl.Area,
        collision_energy = (@p zip(data.eV, rep) |> mapmany(repeat([_[1]], _[2]))),
        FWHM = tbl.FWHM,
        symmetry = tbl.Symmetry
    )
end

function read_mrm(path::String; vendor = :agilent)
    tbl = CSV.read(path, Table)
    if vendor ≡ :agilent
        Table(
            compound = [occursin("[", x) ? eval(Meta.parse(x)) : x for x in getproperty(tbl, Symbol("Compound Name"))],
            mz1 = getproperty(tbl, Symbol("Precursor Ion")),
            mz2 = getproperty(tbl, Symbol("Product Ion")),
            rt = getproperty(tbl, Symbol("Ret Time (min)")),
            Δrt = getproperty(tbl, Symbol("Delta Ret Time")),
            collision_energy = getproperty(tbl, Symbol("Collision Energy")),
            polarity = getproperty(tbl, Symbol("Polarity"))
        )
    end
end

function write_mrm(io, tbl::Table; vendor = :agilent) 
    if vendor ≡ :agilent
        CSV.write(io, tbl; 
            header = ["Compound Name", "Precursor Ion", "Product Ion", "Ret Time (min)", "Delta Ret Time", "Collision Energy", "Polarity"])
    end
end

function Base.show(io::IO, ::MIME"text/plain", class::T) where {T <: ClassSP}
    hasisomer(T) ? begin
        str = join(map(repr, class.isomer), ", ")
        str = isempty(str) ? "" : str 
        print(io, replace(repr(T), "_" => "?"), isempty(str) ? "" : "($str)") 
    end : print(io, T)
end

Base.show(io::IO, ::T) where {T <: ClassSP} = print(io, replace(repr(T), r"_$" => ""))

Base.show(io::IO, ion::Ion) = print(io, repr_adduct(ion.adduct), " of ", ion.molecule)
Base.show(io::IO, ::T) where {T <: Sugar} = print(io, T)
Base.show(io::IO, ::Glycan{T}) where {T <: Tuple} = print(io, join(T.parameters, "_"))

repr_ox = @λ begin
    0 => ""
    1 => ";O"
    n => ";O" * string(n)
end

Base.show(io::IO, sc::LCB) = print(io, ncb(sc), ":", ndb(sc), repr_ox(nox(sc)))
Base.show(io::IO, ::MIME"text/plain", sc::LCB) = print(io, "SPB ", sc)
Base.show(io::IO, ::MIME"text/plain", sc::ACYL) = print(io, "Acyl ", sc)
Base.show(io::IO, sc::Acyl) = print(io, ncb(sc), ":", ndb(sc), repr_ox(nox(sc)))
Base.show(io::IO, sc::Acylα) = print(io, ncb(sc), ":", ndb(sc), 
                                        @match nox(sc) begin
                                            0 => ""
                                            1 => ";(2OH)"
                                            2 => ";(2OH);O"
                                            o => ";(2OH);O" * string(o - 1)
                                        end
                                    )
Base.show(io::IO, sc::Acylβ) = print(io, ncb(sc), ":", ndb(sc), 
                                        @match nox(sc) begin
                                            0 => ""
                                            1 => ";(3OH)"
                                            2 => ";(3OH);O"
                                            o => ";(3OH);O" * string(o - 1)
                                        end
                                    )

Base.show(io::IO, score::Score) = println(io, score.score)

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
    mz2 = cpd.project.data[row.source].mz2[s.mz2_id]
    mode = isa(cpd.project.data[row.source], PreIS) ? "PreIS" : "MRM"
    (ion1 = row.ion1, mz1 = s.mz1, ion2 = row.ion2, mz2 = mz2, area = s.area, error = s.error, CE = s.collision_energy, rt = s.rt, source = mode)
end

Base.show(io::IO, cpd::SPID) = print(io, cpd.class, " ", cpd.sidechain)
Base.show(io::IO, cpd::CompoundSPVanilla) = print(io, cpd.class, " ", cpd.sidechain)

states_color = @λ begin
    1  => "🟢"
    0  => "🟡"
    -1 => "🔴"
end

Base.show(io::IO, sc::SideChain) = print(io, sc.lcb, "/", sc.acyl)
Base.show(io::IO, sc::SumChain) = print(io, ncb(sc), ":", ndb(sc), repr_ox(nox(sc)))

function Base.show(io::IO, ::MIME"text/plain", cpd::CompoundSP)
    class, chain = states_color.(cpd.states)
    print(io, "Compound with ", size(cpd.fragments, 1), " fragments ($class$chain):")
    print(io, "\n∘ ID: ")
    show(io, MIME"text/plain"(), cpd.class)
    print(io, " ", cpd.sidechain)
    print(io, "\n∘ Area: ", cpd.area[1])
    print(io, "\n∘ Error: ", cpd.area[2])
    dt = fragment_table(cpd)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r])
end

Base.show(io::IO, cpd::CompoundSP) = print(io, cpd.class, " ", cpd.sidechain)

function Base.show(io::IO, ::MIME"text/plain", analyte::AnalyteSP)
    sc = states_color.(analyte.states)
    class = sc[states_id(:class)]
    chain = sc[states_id(:chain)]
    rt = sc[states_id(:rt)]
    diq = sc[states_id(:error)]
    isf = sc[states_id(:isf)]
    print(io, "Analytes with ", length(analyte), " compounds @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4), " (id$class$chain,rt$rt,signal$diq$isf):")
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
    sc = states_color.(analyte.states)
    class = sc[states_id(:class)]
    chain = sc[states_id(:chain)]
    rt = sc[states_id(:rt)]
    diq = sc[states_id(:error)]
    isf = sc[states_id(:isf)]
    print(io, isempty(analyte.compounds) ? "?" : last(analyte), " @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4), " (id$class$chain,rt$rt,signal$diq$isf):")
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

Base.show(io::IO, reuseable::ReUseable) = print(io, "ReUsable ", reuseable.query)
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

Base.show(io::IO, pred::Inv) = print(io, "Not(", pred.arg, ")")