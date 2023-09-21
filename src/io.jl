"""
    read_featuretable(path, type = :mzmine3; kwargs...)

Read feature data. 

# Arguments
* `path`: csv file path.
* `type`: the type of csv file. `:mzmine3` for data from MZMINE3; `:masshunter_mrm` or `:mh_mrm` for MRM feature data from MassHunter. 

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable(path, type = :mzmine3; kwargs...)
    if type == :mzmine3
        read_featuretable_mzmine3(path; kwargs...)
    elseif type == :masshunter_mrm || type == :mh_mrm
        read_featuretable_masshunter_mrm(path; kwargs...)
    else
        nothing
    end
end 

"""
    read_featuretable_mzmine3(path;
                                silencewarnings::Bool = false,
                                maxwarnings::Int = 10, 
                                debug::Bool = false
                            )

Read feature data from MZMINE3.

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable_mzmine3(path;
                                silencewarnings::Bool = false,
                                maxwarnings::Int = 10, 
                                debug::Bool = false
                            )
    head = split(readline(path), ",")
    idmain = findall(x -> any(==(x, text) for text in ["id", "rt", "mz", "height", "area"]), head)
    idfwhm = findall(x-> endswith(x, "fwhm"), head)
    idsym = findall(x-> endswith(x, "asymmetry_factor"), head)
    tbl = CSV.read(path, Table; select = idmain, silencewarnings, maxwarnings, debug)
    n = size(tbl, 1)
    fwhm = CSV.read(path, Table; select = idfwhm, silencewarnings, maxwarnings, debug)
    sym = CSV.read(path, Table; select = idsym, silencewarnings, maxwarnings, debug)
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

"""
    read_featuretable_masshunter_mrm(path;
                                        silencewarnings::Bool = false,
                                        maxwarnings::Int = 10, 
                                        debug::Bool = false
                                    )

Read MRM feature data from MassHunter.

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable_masshunter_mrm(path;
                                        silencewarnings::Bool = false,
                                        maxwarnings::Int = 10, 
                                        debug::Bool = false
                                    )
    strs = readlines(path)
    starts = Int[]
    data = Table(ms1 = Float64[], ms2 = Float64[], eV = Int[], datafile = String[])
    ends = Int[]
    status = false
    for (i, l) in enumerate(strs)
        if !status
            transition = match(r"\((\d*.\d*) -> (\d*.\d*)\)", l)
            isnothing(transition) && continue
            ms1, ms2 = parse.(Float64, transition)
            eV = parse(Float64, match(r"CID@(\d*.\d*)", l)[1])
            datafile = match(r"(.*\.d).*\.d", l)[1]
            push!(data, (; ms1, ms2, eV, datafile))
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
    tbl = CSV.read(str, Table; select = (i, name) -> any(==(text, String(name)) for text in txt), silencewarnings, maxwarnings, debug)
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
        symmetry = tbl.Symmetry,
        datafile = (@p zip(data.datafile, rep) |> mapmany(repeat([_[1]], _[2]))),
    )
end

"""
    read_transition(path::String, vendor = :agilent;
                        silencewarnings::Bool = false,
                        maxwarnings::Int = 10, 
                        debug::Bool = false
                    )
                    
Read MRM transition table.

# Arguments
* `path`: csv file path.
* `vendor`: the vendor of data source. `:agilent` for data from Agilent. 

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_transition(path::String, vendor = :agilent;
                            silencewarnings::Bool = false,
                            maxwarnings::Int = 10, 
                            debug::Bool = false
                        )
    tbl = CSV.read(path, Table; silencewarnings, maxwarnings, debug)
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

"""
    write_transition(io, tbl::Table, vendor = :agilent)

Write MRM transition table into csv file.

# Arguments
* `io`: csv file path.
* `tbl`: transition table.
* `vendor`: desired vendor for transition table. `:agilent` for Agilent. 
"""
function write_transition(io, tbl::Table, vendor = :agilent)
    if vendor ≡ :agilent
        CSV.write(io, tbl;
            header = ["Compound Name", "Precursor Ion", "Product Ion", "Ret Time (min)", "Delta Ret Time", "Collision Energy", "Polarity"])
    end
end

function Base.show(io::IO, ::MIME"text/plain", class::T) where {T <: ClassSP}
    hasisomer(T) ? begin
        str = join(map(repr, class.isomer), ", ")
        str = isempty(str) ? "" : str
        print(io, repr(class), isempty(str) ? "" : "($str)")
    end : print(io, repr(class))
end

Base.show(io::IO, ::T) where {T <: ClassSP} = print(io, replace(repr(T), r"_$" => "", "SphingolipidsID." => ""))

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

fragmenttable(cpd) = map(cpd.fragments) do row
    s = query_raw(cpd.project, row.source, row.id)
    mz2 = cpd.project.data[row.source].mz2[s.mz2_id]
    mode = isa(cpd.project.data[row.source], PreIS) ? "PreIS" : "MRM"
    (ion1 = row.ion1, mz1 = s.mz1, ion2 = row.ion2, mz2 = mz2, area = s.area, error = s.error, CE = s.collision_energy, rt = s.rt, source = mode)
end

Base.show(io::IO, cpd::SPID) = print(io, cpd.class, " ", cpd.chain)
Base.show(io::IO, ::MIME"text/html", cpd::SPID) = print(io, cpd.class, " ", cpd.chain)
Base.show(io::IO, cpd::CompoundSPVanilla) = print(io, cpd.class, " ", cpd.chain)

states_color = @λ begin
    1  => "🟢"
    0  => "🟡"
    -1 => "🔴"
end

Base.show(io::IO, sc::ChainSP) = print(io, sc.lcb, "/", sc.acyl)
Base.show(io::IO, sc::SumChain) = print(io, ncb(sc), ":", ndb(sc), repr_ox(nox(sc)))

function Base.show(io::IO, ::MIME"text/plain", cpd::CompoundSP)
    class, chain = states_color.(cpd.states)
    print(io, "Compound with ", size(cpd.fragments, 1), " fragments ($class$chain):")
    print(io, "\n∘ ID: ")
    show(io, MIME"text/plain"(), cpd.class)
    print(io, " ", cpd.chain)
    print(io, "\n∘ Area: ", cpd.area[1])
    print(io, "\n∘ Error: ", cpd.area[2])
    dt = fragmenttable(cpd)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r])
end

Base.show(io::IO, cpd::CompoundSP) = print(io, cpd.class, " ", cpd.chain)

function Base.show(io::IO, ::MIME"text/plain", analyte::AnalyteSP)
    sc = states_color.(analyte.states)
    class = sc[states_id(:class)]
    chain = sc[states_id(:chain)]
    rt = sc[states_id(:rt)]
    diq = sc[states_id(:error)]
    isf = sc[states_id(:isf)]
    total = sc[states_id(:total)]
    manual = states_color(analyte.manual_check)
    sc = map([analyte.cpdsc, analyte.score]) do sc
        sc.first > 0 ? *(string(SPDB[:SCORE].param[sc.first].target), " => ", string(sc.second)) : ""
    end
    scs = join(filter!(!isempty, sc), ", ")
    print(io, "Analytes with ", length(analyte), " compounds @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4), " $(manual)$(total)id$(class)$(chain)rt$(rt)sig$(diq)$(isf):")
    print(io, "\n∘ Score: ", scs)
    print(io, "\n∘ Compounds:")
    for cpd in analyte
        print(io, "\n ", cpd)
    end
    dt = mapreduce(fragmenttable, vcat, analyte)
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
    total = sc[states_id(:total)]
    manual = states_color(analyte.manual_check)
    print(io, isempty(analyte.compounds) ? "?" : last(analyte), " @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4), " $(manual)$(total)id$(class)$(chain)rt$(rt)sig$(diq)$(isf):")
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

function Base.show(io::IO, qcmd::QueryCommands)
    length(qcmd.qcmd) > 1 && print(io, "(")
    print_init(io, qcmd)
    length(qcmd.qcmd) > 1 && print(io, ")")
end
function Base.show(io::IO, qcmd::QueryNot)
    print(io, "¬(")
    print_init(io, qcmd.qcmd)
    print(io, ")")
end
Base.show(io::IO, qcmd::QueryCmd) = print(io, qcmd.query)
print_init(io::IO, qcmd::QueryAnd) = print(io, join(repr.(qcmd.qcmd), " ∧ "))
print_init(io::IO, qcmd::QueryOr) = print(io, join(repr.(qcmd.qcmd), " ∨ "))
print_init(io::IO, qcmd::QueryCmd) = print(io, qcmd)
print_init(io::IO, qcmd::QueryNot) = print(io, qcmd)
Base.show(io::IO, reuseable::ReUseable) = print(io, "ReUsable ", reuseable.query)
function Base.show(io::IO, aquery::Query)
    print(io, "Query with ", length(aquery), " analytes: \n")
    print(io, "∘ Queries: ")
    print_init(io, aquery.query)
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