"""
    read_featuretable(file, type = :mzmine3; kwargs...)

Read feature data. 

# Arguments
* `file`: csv file path.
* `type`: the type of csv file. `:mzmine3` for data from MZMINE3; `:masshunter_mrm` or `:mh_mrm` for MRM feature data from MassHunter. 

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable(file::String, type = :mzmine3; kwargs...)
    if type == :mzmine3
        read_featuretable_mzmine3(file; kwargs...)
    elseif type == :masshunter_mrm || type == :mh_mrm
        read_featuretable_masshunter_mrm(file; kwargs...)
    else
        nothing
    end
end 

"""
    read_featuretable_mzmine3(file;
                                silencewarnings::Bool = false,
                                maxwarnings::Int = 10, 
                                debug::Bool = false, 
                                delim = ',', kwargs...
                            )

Read feature data from MZMINE3.

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable_mzmine3(file::String;
                                silencewarnings::Bool = false,
                                maxwarnings::Int = 10, 
                                debug::Bool = false, 
                                delim = ',', kwargs...
                            )
    head = split(readline(file), ",")
    idmain = findall(x -> any(==(x, text) for text in ["id", "rt", "mz", "height", "area"]), head)
    idfwhm = findall(x-> endswith(x, "fwhm"), head)
    idsym = findall(x-> endswith(x, "asymmetry_factor"), head)
    tbl = CSV.read(file, Table; select = idmain, silencewarnings, maxwarnings, debug, delim, kwargs...)
    n = size(tbl, 1)
    fwhm = CSV.read(file, Table; select = idfwhm, silencewarnings, maxwarnings, debug, delim, kwargs...)
    sym = CSV.read(file, Table; select = idsym, silencewarnings, maxwarnings, debug, delim, kwargs...)
    datafile = Dict(propertynames(fwhm) .=> map(col -> match(r".*:(.*):.*", string(col))[1], propertynames(fwhm)))
    id = findfirst.(!ismissing, fwhm)
    uid = unique(id)
    tbl = Table(
        id = collect(1:n),
        mz1 = tbl.mz,
        mz2 = zeros(Float64, n),
        mz_range = repeat([RealInterval(0, 0, <=, <=)], n),
        rt = tbl.rt,
        height = tbl.height,
        area = tbl.area,
        collision_energy = zeros(Int, n),
        FWHM = getindex.(fwhm, id),
        symmetry = get.(sym, replace!(findfirst.(!ismissing, sym), nothing => :_symmetry), 1.0),
        polarity = trues(n),
        datafile = getindex.(Ref(datafile), id),
        injection_order = map(x -> findfirst(==(x), uid), id)
    )
    tbl
end

"""
    read_featuretable_masshunter_mrm(file;
                                        silencewarnings::Bool = false,
                                        maxwarnings::Int = 10, 
                                        debug::Bool = false, 
                                        delim = ',', kwargs...
                                    )

Read MRM feature data from MassHunter.

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_featuretable_masshunter_mrm(file::String;
                                        silencewarnings::Bool = false,
                                        maxwarnings::Int = 10, 
                                        debug::Bool = false, 
                                        delim = ',', kwargs...
                                    )
    strs = readlines(file)
    starts = Int[]
    data = Table(ms1 = Float64[], ms2 = Float64[], eV = Int[], polarity = Bool[], datafile = String[])
    ends = Int[]
    status = false
    for (i, l) in enumerate(strs)
        if !status
            pol = match(r"([+-])ESI MRM Frag", l)
            isnothing(pol) && continue
            polarity = pol[1] == "+"
            transition = match(r"\((\d*.\d*) -> (\d*.\d*)\)", l)
            isnothing(transition) && continue
            ms1, ms2 = parse.(Float64, transition)
            eV = parse(Float64, match(r"CID@(\d*.\d*)", l)[1])
            datafile = match(r"(.*\.d).*\.d", l)[1]
            push!(data, (; ms1, ms2, eV, polarity, datafile))
            push!(starts, i + 2)
            status = true
        elseif l == "" || unique(l) == [',']
            push!(ends, i - 1)
            status = false
        elseif occursin("No integration results available", l)
            pop!(starts)
            pop!(data)
            status = false
        end
    end
    if length(ends) == length(starts) - 1
        push!(ends, length(strs))
    end
    rep = @. ends - starts + 1
    starts[1] -= 1
    str = IOBuffer(join(mapreduce(vcat, starts, ends) do st, ed
        strs[st:ed]
    end, "\n"))
    txt = ["RT", "Height", "Area", "Symmetry", "FWHM"]
    tbl = CSV.read(str, Table; select = (i, name) -> any(==(text, String(name)) for text in txt), silencewarnings, maxwarnings, debug, delim, kwargs...)
    n = size(tbl, 1)
    tbl = Table(
        id = collect(1:n),
        mz1 = (@p zip(data.ms1, rep) |> mapmany(repeat([_[1]], _[2]))),
        mz2 = (@p zip(data.ms2, rep) |> mapmany(repeat([_[1]], _[2]))),
        rt = tbl.RT,
        height = tbl.Height,
        area = tbl.Area,
        collision_energy = (@p zip(data.eV, rep) |> mapmany(repeat([_[1]], _[2]))),
        FWHM = tbl.FWHM,
        symmetry = tbl.Symmetry,
        polarity = (@p zip(data.polarity, rep) |> mapmany(repeat([_[1]], _[2]))),
        datafile = (@p zip(data.datafile, rep) |> mapmany(repeat([_[1]], _[2])))
    )
    tbl = tbl[tbl.symmetry .!= "MM"]
    udf = unique(tbl.datafile)
    Table(tbl; symmetry = parse.(Float64, string.(tbl.symmetry)), injection_order = map(x -> findfirst(==(x), udf), tbl.datafile))
end

"""
    read_transition(file::String, vendor = :agilent;
                        silencewarnings::Bool = false,
                        maxwarnings::Int = 10, 
                        debug::Bool = false, 
                        delim = ',', kwargs...
                    )
                    
Read MRM transition table.

# Arguments
* `file`: csv file path.
* `vendor`: the vendor of data source. `:agilent` for data from Agilent. 

# Keyword Arguments
* `silencewarnings`: whether invalid value warnings should be silenced.
* `maxwarnings`: if more than `maxwarnings` number of warnings are printed while parsing, further warnings will be
silenced by default; for multithreaded parsing, each parsing task will print up to `maxwarnings`.
* `debug`: passing `true` will result in many informational prints while a dataset is parsed; can be useful when
reporting issues or figuring out what is going on internally while a dataset is parsed.
"""
function read_transition(file::String, vendor = :agilent;
                            silencewarnings::Bool = false,
                            maxwarnings::Int = 10, 
                            debug::Bool = false, 
                            delim = ',', kwargs...
                        )
    tbl = CSV.read(file, Table; silencewarnings, maxwarnings, debug, delim, kwargs...)
    if vendor ≡ :agilent
        Table(
            id = collect(eachindex(tbl)),
            analyte = [occursin("|", x) ? split(x, " | ") : x for x in getproperty(tbl, Symbol("Compound Name"))],
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
    write_transition(file, tbl::Table, vendor = :agilent, delim = ',', kwargs...)

Write MRM transition table into csv file.

# Arguments
* `file`: csv file path.
* `tbl`: transition table.
* `vendor`: desired vendor for transition table. `:agilent` for Agilent. 
"""
function write_transition(file, tbl::Table, vendor = :agilent)
    if vendor ≡ :agilent
        CSV.write(file, tbl;
            header = ["Compound Name", "Precursor Ion", "Product Ion", "Ret Time (min)", "Delta Ret Time", "Collision Energy", "Polarity"], delim, kwargs...)
    end
end

function write(file::String, qt::Quantification{A}; delim = '\t') where A
    mkpath(file)
    ChemistryQuantitativeAnalysis.write(joinpath(file, "batch.batch"), qt.batch; delim)
    JLD2.save_object(joinpath(file, "config.jld2"), set!(copy(qt.config), :analytetype, A))
end

data_extension(::PreIS, files::Vararg{String}) = string(joinpath(files...), ".preis")
data_extension(::MRM, files::Vararg{String}) = string(joinpath(files...), ".mrm")
data_extension(::Type{PreIS}, files::Vararg{String}) = string(joinpath(files...), ".preis")
data_extension(::Type{MRM}, files::Vararg{String}) = string(joinpath(files...), ".mrm")
data_extension(::QCData, files::Vararg{String}) = string(joinpath(files...), ".qcd")
data_extension(::SerialDilution, files::Vararg{String}) = string(joinpath(files...), ".sdd")
data_extension(::QuantData, files::Vararg{String}) = string(joinpath(files...), ".qtd")

function write(file::String, data::PreIS; delim = '\t', kwargs...)
    mkpath(file)
    JLD2.save_object(joinpath(file, "config.jld2"), set!(data.config, :table, pairs((; kwargs..., delim))))
    delete!(data.config, :table)
    JLD2.save_object(joinpath(file, "mz2.jld2"), data.mz2)
    JLD2.save_object(joinpath(file, "range.jld2"), data.range)
    tbl = Table((id = data.table.id, mz1 = data.table.mz1, mz2 = map(i -> data.mz2[i], data.table.mz2_id), rt = data.table.rt), data.table; 
            mz2_id = nothing, mz_range = map(i -> data.range[i], data.table.mz2_id), polarity = repeat([data.polarity], length(data.table)))
    CSV.write(joinpath(file, "table.txt"), tbl; delim, kwargs...)
end

function write(file::String, data::MRM; delim = '\t', kwargs...)
    mkpath(file)
    JLD2.save_object(joinpath(file, "config.jld2"), set!(data.config, :table, pairs((; kwargs..., delim))))
    delete!(data.config, :table)
    JLD2.save_object(joinpath(file, "mz2.jld2"), data.mz2)
    tbl = Table((id = data.table.id, mz1 = data.table.mz1, mz2 = map(i -> data.mz2[i], data.table.mz2_id), rt = data.table.rt), data.table; 
            mz2_id = nothing, polarity = repeat([data.polarity], length(data.table)))
    CSV.write(joinpath(file, "table.txt"), tbl; delim, kwargs...)
end

function write(file::String, data::QCData{T}; delim = '\t', kwargs...) where T
    mkpath(file)
    JLD2.save_object(joinpath(file, "config.jld2"), data.config)
    write(data_extension(T, file, "raw"), data.raw; delim, kwargs...)
    ChemistryQuantitativeAnalysis.write(joinpath(file, "table.at"), data.table; delim)
end

function write(file::String, data::SerialDilution{T}; delim = '\t', kwargs...) where T
    mkpath(file)
    JLD2.save_object(joinpath(file, "config.jld2"), data.config)
    write(data_extension(T, file, "raw"), data.raw; delim, kwargs...)
    ChemistryQuantitativeAnalysis.write(joinpath(file, "batch.batch"), data.batch; delim)
end

function write(file::String, data::QuantData{T}; delim = '\t', kwargs...) where T
    mkpath(file)
    JLD2.save_object(joinpath(file, "config.jld2"), data.config)
    write(data_extension(T, file, "raw"), data.raw; delim, kwargs...)
    ChemistryQuantitativeAnalysis.write(joinpath(file, "table.at"), data.table; delim)
end

function write(file::String, project::Project; delim = '\t', kwargs...)
    mkpath(file)
    emptyp = Project()
    analyte = map(project.analyte) do ana
        ana = copy(ana)
        ana.compound = copy.(ana.compound)
        for c in ana
            c.project = emptyp
        end
        ana
    end
    JLD2.save_object(joinpath(file, "analyte.jld2"), analyte)
    methods = Vector{Union{Nothing, Missing, MethodTable}}(undef, length(project.data))
    rt_corrections = Vector{Union{Nothing, RTCorrection}}(undef, length(project.data))
    fill!(methods, nothing)
    fill!(rt_corrections, nothing)
    for (i, data) in enumerate(project.data)
        if !haskey(data.config, :method) && !haskey(data.config, :rt_correction)
            write(data_extension(data, file, "data", string(i)), data; delim, kwargs...)
        else
            if haskey(data.config, :method)
                methods[i] = data.config[:method]
                set!(data.config, :method, nothing)
            end
            if haskey(data.config, :rt_correction)
                rt_corrections[i] = data.config[:rt_correction]
                set!(data.config, :rt_correction, nothing)
            end
            write(data_extension(data, file, "data", string(i)), data; delim, kwargs...)
            haskey(data.config, :method) && set!(data.config, :method, methods[i])
            haskey(data.config, :rt_correction) && set!(data.config, :rt_correction, rt_corrections[i])
        end
    end
    if isdefined(project.quantification, :batch)
        methods = map(methods) do m
            m === project.quantification.batch.method ? missing : m
        end
    end
    JLD2.save_object(joinpath(file, "methods.jld2"), methods)
    JLD2.save_object(joinpath(file, "rt_corrections.jld2"), rt_corrections)
    qt = project.quantification
    if !isdefined(qt, :config) || (!haskey(qt.config, :quantdata) && !haskey(qt.config, :qcdata) && !haskey(qt.config, :serialdilution))
        write(joinpath(file, "quantification.qt"), qt; delim)
    else
        qtdata = Dictionary{Symbol, Int}()
        for type in [:qcdata, :serialdilution, :quantdata]
            data = get(qt.config, type, nothing)
            id = isnothing(data) ? 0 : findfirst(x -> x === data, project.data)
            isnothing(id) && throw(ErrorException("Data $data is not in `project.data`"))
            set!(qtdata, type, id)
            set!(qt.config, type, nothing)
        end
        JLD2.save_object(joinpath(file, "qtdata.jld2"), qtdata)
        write(joinpath(file, "quantification.qt"), qt; delim)
        for type in [:qcdata, :serialdilution, :quantdata]
            i = get(qtdata, type, nothing)
            isnothing(i) || (i == 0 ? delete!(qt.config, type) : set!(qt.config, type, project.data[i]))
        end
    end
    JLD2.save_object(joinpath(file, "appendix.jld2"), project.appendix)
end

function read_quantification(file::String)
    endswith(file, ".qt") || throw(ArgumentError("The file is not a valid Quantification directory"))
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    batch = ChemistryQuantitativeAnalysis.read(joinpath(file, "batch.batch"), Table; analytetype = config[:analytetype])
    Quantification(batch, delete!(config, :analytetype))
end

function read_mrm(file::String; analytetype = TransitionID)
    endswith(file, ".mrm") || throw(ArgumentError("The file is not a valid MRM directory"))
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    tbl = CSV.read(joinpath(file, "table.txt"), Table; config[:table]...)
    mz2 = JLD2.load_object(joinpath(file, "mz2.jld2"))
    polarity = only(unique(tbl.polarity))
    MRM(Table(tbl; mz2_id = map(x -> findfirst(==(x), mz2), tbl.mz2), mz2 = nothing, polarity = nothing), mz2, polarity, delete!(config, :table))
end

function read_preis(file::String; analytetype = TransitionID)
    endswith(file, ".preis") || throw(ArgumentError("The file is not a valid PreIS directory"))
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    tbl = CSV.read(joinpath(file, "table.txt"), Table; config[:table]...)
    mz2 = JLD2.load_object(joinpath(file, "mz2.jld2"))
    range = JLD2.load_object(joinpath(file, "range.jld2"))
    polarity = only(unique(tbl.polarity))
    PreIS(Table(tbl; mz2_id = map(x -> findfirst(==(x), mz2), tbl.mz2), mz2 = nothing, range = nothing, polarity = nothing), range, mz2, polarity,  delete!(config, :table))
end

function read_qcdata(file::String; analytetype = TransitionID)
    endswith(file, ".qcd") || throw(ArgumentError("The file is not a valid QCData directory"))
    rawfile = only(filter(startswith("raw"), readdir(file)))
    raw = read(joinpath(file, rawfile))
    table = ChemistryQuantitativeAnalysis.read(joinpath(file, "table.at"), Table; analytetype)
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    QCData(raw, table, config)
end

function read_serialdilution(file::String; analytetype = TransitionID)
    endswith(file, ".sdd") || throw(ArgumentError("The file is not a valid SerialDilution directory"))
    rawfile = only(filter(startswith("raw"), readdir(file)))
    raw = read(joinpath(file, rawfile))
    batch = ChemistryQuantitativeAnalysis.read(joinpath(file, "batch.batch"), Table; analytetype)
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    SerialDilution(raw, batch, config)
end

function read_quantdata(file::String; analytetype = TransitionID)
    endswith(file, ".qtd") || throw(ArgumentError("The file is not a valid QuantData directory"))
    rawfile = only(filter(startswith("raw"), readdir(file)))
    raw = read(joinpath(file, rawfile))
    table = ChemistryQuantitativeAnalysis.read(joinpath(file, "table.at"), Table; analytetype)
    config = JLD2.load_object(joinpath(file, "config.jld2"))
    QuantData(raw, table, config)
end

function read_project(file::String)
    endswith(file, ".project") || throw(ArgumentError("The file is not a valid Project directory"))
    appendix = JLD2.load_object(joinpath(file, "appendix.jld2"))
    quantification = read_quantification(joinpath(file, "quantification.qt"))
    analytetype = isdefined(quantification, :batch) ? first(typeof(quantification).parameters) : Symbol
    methods = map(JLD2.load_object(joinpath(file, "methods.jld2"))) do x
        ismissing(x) ? quantification.batch.method : x
    end
    rt_corrections = JLD2.load_object(joinpath(file, "rt_corrections.jld2"))
    data = map(x -> read(x; analytetype), readdir(joinpath(file, "data"); join = true))
    for (dt, m, r) in zip(data, methods, rt_corrections)
        set!(dt.config, :method, m)
        set!(dt.config, :rt_correction, r)
    end
    if in("qtdata.jld2", readdir(file)) && isdefined(quantification, :config)
        qcdata = JLD2.load_object(joinpath(file, "qtdata.jld2"))
        for (k, v) in pairs(qcdata)
            v > 0 ? set!(quantification.config, k, data[v]) : delete!(quantification.config, k)
        end
    end
    analyte = JLD2.load_object(joinpath(file, "analyte.jld2"))
    project = last(first(analyte)).project
    project.analyte = analyte
    project.data = data
    project.quantification = quantification
    project.appendix = appendix
    project
end

function read(file::String; analytetype = TransitionID)
    if endswith(file, ".project")
        read_project(file)
    elseif endswith(file, ".qtd")
        read_quantdata(file; analytetype)
    elseif endswith(file, ".sdd")
        read_serialdilution(file; analytetype)
    elseif endswith(file, ".qcd")
        read_qcdata(file; analytetype)
    elseif endswith(file, ".mrm")
        read_mrm(file)
    elseif endswith(file, ".preis")
        read_preis(file)
    elseif endswith(file, ".qt")
        read_quantification(file)
    end
end

function TransitionID(analyte::AbstractString)
    analyte == "nothing" && return nothing
    cls, chn = split(String(analyte), " ")
    cls = eval(Meta.parse(cls))
    if occursin("/", chn)
        lb, acy = split(chn, "/")
        cb_lcb, db_lcb = match(r"(\d*):(\d*)", lb)
        cb_acyl, db_acyl = match(r"(\d*):(\d*)", acy)
        ox_lcb = match(r";O(\d*)", lb)
        ox_acyl = match(r";O(\d*)", acy)
        ox_lcb = isnothing(ox_lcb) ? 0 : isempty(first(ox_lcb)) ? 1 : parse(Int, first(ox_lcb))
        ox_acyl = isnothing(ox_acyl) ? 0 : isempty(first(ox_acyl)) ? 1 : parse(Int, first(ox_acyl))
        cons_acyl = @match acy begin
            r";(2OH)" => (ox_acyl += 1; acylα)
            r";(3OH)" => (ox_acyl += 1; acylβ)
            _         => acyl
        end
        n13C_lcb = match(r"\(.*13C(\d).*\)", lb)
        n13C_acyl = match(r"\(.*13C(\d).*\)", acy)
        nD_lcb = match(r"\(.*D(\d).*\)", lb)
        nD_acyl = match(r"\(.*D(\d).*\)", acy)
        n13C_lcb = isnothing(n13C_lcb) ? 0 : isempty(first(n13C_lcb)) ? 1 : parse(Int, first(n13C_lcb))
        n13C_acyl = isnothing(n13C_acyl) ? 0 : isempty(first(n13C_acyl)) ? 1 : parse(Int, first(n13C_acyl))
        nD_lcb = isnothing(nD_lcb) ? 0 : isempty(first(nD_lcb)) ? 1 : parse(Int, first(nD_lcb))
        nD_acyl = isnothing(nD_acyl) ? 0 : isempty(first(nD_acyl)) ? 1 : parse(Int, first(nD_acyl))
        cpd = spid(cls, chain(
            lcb(parse(Int, cb_lcb), parse(Int, db_lcb), ox_lcb; n13C = n13C_lcb, nD = nD_lcb), 
            cons_acyl(parse(Int, cb_acyl), parse(Int, db_acyl), ox_acyl; n13C = n13C_acyl, nD = nD_acyl))
            )
    else
        cb, db, ox = match(r"(\d*):(\d*);O(\d*)", chn)
        n13C = match(r"\(.*13C(\d).*\)", chn)
        nD = match(r"\(.*D(\d).*\)", chn)
        n13C = isnothing(n13C) ? 0 : isempty(first(n13C)) ? 1 : parse(Int, first(n13C))
        nD = isnothing(nD) ? 0 : isempty(first(nD)) ? 1 : parse(Int, first(nD))
        ox = isnothing(ox) ? 0 : isempty(ox) ? 1 : parse(Int, ox)
        cpd = spid(cls, chain(parse(Int, cb), parse(Int, db), ox; n13C, nD))
    end
    if occursin("(qualifier)", analyte)
        TransitionID(cpd, false)
    else
        TransitionID(cpd, true)
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

repr_ox(sc::Union{Acylα, AcylαIS}) = 
    @match nox(sc) begin
        0 => ""
        1 => ";(2OH)"
        2 => ";(2OH);O"
        o => ";(2OH);O" * string(o - 1)
    end
repr_ox(sc::Union{Acylβ, AcylβIS}) = 
    @match nox(sc) begin
        0 => ""
        1 => ";(3OH)"
        2 => ";(3OH);O"
        o => ";(3OH);O" * string(o - 1)
    end
repr_ox(sc) = 
    @match nox(sc) begin
        0 => ""
        1 => ";O"
        n => ";O" * string(n)
    end

print_comp(io::IO, x) = print(io, x)
print_comp(io::IO, x...) = foreach(s -> print_comp(io, s), x)
Base.show(io::IO, sc::LCB) = (print(io, "SPB "); print_comp(io, sc))
print_comp(io::IO, sc::LCB) = print(io, ncb(sc), ":", ndb(sc), repr_ox(sc))
print_comp(io::IO, sc::LCBIS) = print(io, ncb(sc), ":", ndb(sc), repr_ox(sc), repr_is(sc.isotope))
#Base.show(io::IO, ::MIME"text/plain", sc::LCB) = print(io, "SPB ", sc)
#Base.show(io::IO, ::MIME"text/plain", sc::ACYL) = print(io, "Acyl ", sc)
Base.show(io::IO, sc::Acyl) = (print(io, "Acyl "); print_comp(io, sc))
print_comp(io::IO, sc::ACYL) = print(io, ncb(sc), ":", ndb(sc), repr_ox(sc))
print_comp(io::IO, sc::ACYLIS) = print(io, ncb(sc), ":", ndb(sc), repr_ox(sc), repr_is(sc.isotope))

fragmenttable(cpd) = map(cpd.fragment) do row
    s = query_data(cpd.project, row.source, row.id)
    mz2 = cpd.project.data[row.source].mz2[s.mz2_id]
    mode = isa(cpd.project.data[row.source], PreIS) ? "PreIS" : "MRM"
    (ion1 = row.ion1, mz1 = s.mz1, ion2 = row.ion2, mz2 = mz2, area = s.area, height = s.height, error = s.error, CE = s.collision_energy, rt = s.rt, source = mode)
end

Base.show(io::IO, cpd::SPID) = print(io, cpd.class, " ", cpd.chain)
Base.show(io::IO, ::MIME"text/html", cpd::SPID) = print(io, cpd.class, " ", cpd.chain)
Base.show(io::IO, cpd::CompoundSPVanilla) = print(io, cpd.class, " ", cpd.chain)

states_color = @λ begin
    1  => "🟢"
    0  => "🟡"
    -1 => "🔴"
end

Base.show(io::IO, sc::ChainSP) = print_comp(io, sc.lcb, "/", sc.acyl)
Base.show(io::IO, sc::SumChain) = print_comp(io, ncb(sc), ":", ndb(sc), repr_ox(sc))
Base.show(io::IO, sc::SumChainIS) = print_comp(io, ncb(sc), ":", ndb(sc), repr_ox(sc), repr_is(sc.isotope))
function repr_is(is::NamedTuple)
    ls = String[]
    is.n13C > 0 && push!(ls, string("13C", is.n13C))
    is.nD > 0 && push!(ls, string("D", is.nD))
    string("(", join(ls, ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", cpd::CompoundSP)
    class, chain = states_color.(cpd.state)
    print(io, "Compound with ", size(cpd.fragment, 1), " fragments|$class$chain:")
    print(io, "\n∘ ID: ")
    show(io, MIME"text/plain"(), cpd.class)
    print(io, " ", cpd.chain)
    print(io, "\n∘ Estimate: ", cpd.signal[1])
    print(io, "\n∘ Error: ", cpd.signal[2])
    dt = fragmenttable(cpd)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Height", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r, :r])
end

Base.show(io::IO, cpd::CompoundSP) = print(io, cpd.class, " ", cpd.chain)

function Base.show(io::IO, ::MIME"text/plain", analyte::AnalyteSP)
    sc = states_color.(analyte.state)
    class = sc[state_id(:class)]
    chain = sc[state_id(:chain)]
    rt = sc[state_id(:rt)]
    diq = sc[state_id(:error)]
    isf = sc[state_id(:isf)]
    total = sc[state_id(:total)]
    manual = sc[state_id(:manual)]
    sc = map([analyte.cpdsc, analyte.score]) do sc
        sc.first > 0 ? *(string(SPDB[:SCORE].param[sc.first].target), " => ", string(sc.second)) : ""
    end
    scs = join(filter!(!isempty, sc), ", ")
    print(io, "Analytes with ", length(analyte), " compounds @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 2), " st$(manual)$(total)id$(class)$(chain)rt$(rt)sig$(diq)$(isf):")
    print(io, "\n∘ Score: ", scs)
    print(io, "\n∘ Compounds:")
    for cpd in analyte
        print(io, "\n ", cpd)
    end
    dt = mapreduce(fragmenttable, vcat, analyte)
    println(io, "\n∘ Fragments: ")
    PrettyTables.pretty_table(io, dt; header = ["Ion1", "m/z", "Ion2", "m/z", "Area", "Height", "Error", "CE (eV)", "RT (min)", "Source"], header_alignment = :l, alignment = [:r, :l, :r, :l, :r, :r, :r, :r, :r, :r])
end

function Base.show(io::IO, analyte::AnalyteSP)
    sc = states_color.(analyte.state)
    class = sc[state_id(:class)]
    chain = sc[state_id(:chain)]
    rt = sc[state_id(:rt)]
    diq = sc[state_id(:error)]
    isf = sc[state_id(:isf)]
    total = sc[state_id(:total)]
    manual = sc[state_id(:manual)]
    print(io, isempty(analyte.compound) ? "?" : last(analyte), " @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4), " st$(manual)$(total)id$(class)$(chain)rt$(rt)sig$(diq)$(isf)")
end

Base.show(io::IO, analyte::AnalyteID) = print(io, isempty(analyte.compound) ? "?" : last(analyte), " @", round(analyte.rt, digits = 2), " MW=", round(mw(analyte), digits = 4))
Base.show(io::IO, transition::TransitionID) = print(io, transition.compound, transition.quantifier ? "" : " (qualifier)")

Base.show(io::IO, data::PreIS) = print(io, "PreIS in ", data.polarity ? "positive" : "negative", " ion mode with ", length(data.table), " features")
function Base.show(io::IO, ::MIME"text/plain", data::PreIS)
    print(io, data, ":\n")
    for dt in zip(data.range, data.mz2)
        println(io, " ", dt[1], " -> ", round(dt[2]; digits = 4))
    end
end

Base.show(io::IO, data::MRM) = print(io, "MRM in ", data.polarity ? "positive" : "negative", " ion mode with ", length(data.table), " features")
function Base.show(io::IO, ::MIME"text/plain", data::MRM)
    print(io, data, ":\n")
    for (i, m) in enumerate(data.mz2)
        v = @view data.table.mz1[data.table.mz2_id .== i]
        v = length(v) > 10 ? string(join(round.(v[1:5]; digits = 4), ", "), ", ..., ", join(round.(v[end - 4:end]; digits = 4))) : join(round.(v; digits = 4), ", ")
        println(io, " ", v, " -> ", round(m; digits = 4))
    end
end

Base.show(io::IO, data::QuantData{T}) where T = print(io, "QuantData{$T} in ", data.raw.polarity ? "positive" : "negative", " ion mode with ", length(data.analyte), " analytes")
function Base.show(io::IO, ::MIME"text/plain", data::QuantData{T}) where T 
    print(io, data, ":\n")
    show(io, MIME"text/plain"(), data.table)
    println(io)
end
Base.show(io::IO, data::QCData{T}) where T = print(io, "QCData{$T} in ", data.raw.polarity ? "positive" : "negative", " ion mode with ", length(data.analyte), " analytes")
function Base.show(io::IO, ::MIME"text/plain", data::QCData{T}) where T 
    print(io, data, ":\n")
    show(io, MIME"text/plain"(), qualitytable(data))
    println(io)
end
Base.show(io::IO, data::SerialDilution{T}) where T = print(io, "SerialDilution{$T} in ", data.raw.polarity ? "positive" : "negative", " ion mode with ", length(data.analyte), " analytes")
function Base.show(io::IO, ::MIME"text/plain", data::SerialDilution{T}) where T 
    print(io, data, ":\n")
    show(io, MIME"text/plain"(), qualitytable(data))
    println(io)
end
Base.show(io::IO, qt::Quantification) = print(io, "Quantification with ", length(qt.analyte), " analytes")
function Base.show(io::IO, ::MIME"text/plain", qt::Quantification)
    print(io, qt, ":\n")
    show(io, MIME"text/plain"(), qt.batch)
    println(io)
end
Base.show(io::IO, rtcor::RTCorrection) = print(io, "RTCorrection with ", length(rtcor.model), " correction functions")
function Base.show(io::IO, ::MIME"text/plain", rtcor::RTCorrection)
    print(io, rtcor, ":\n")
    show(io, MIME"text/plain"(), rtcor.model)
    println(io)
end

Base.show(io::IO, pj::Project) = print(io, "Project with ", length(pj), " analytes")
function Base.show(io::IO, ::MIME"text/plain", pj::Project)
    print(io, pj, ":\n")
    println(io, "∘ Salt: ", get(pj.appendix, :anion, "unknown"))
    println(io, "∘ Signal: ", get(pj.appendix, :signal, "unknown"))
    println(io, "∘ Data: ")
    for data in pj.data
        show(io, MIME"text/plain"(), data)
        println(io)
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
Base.show(io::IO, qcmd::QueryCmd) = print_qcmd(io, qcmd.query)
function print_init(io::IO, qcmd::QueryAnd)
    print_qcmd(io, first(qcmd.qcmd))
    for q in @view qcmd.qcmd[2:end]
        print(io, " ∧ ")
        print_qcmd(io, q)
    end
end
function print_init(io::IO, qcmd::QueryOr)
    print_qcmd(io, first(qcmd.qcmd))
    for q in @view qcmd.qcmd[2:end]
        print(io, " ∨ ")
        print_qcmd(io, q)
    end
end
print_init(io::IO, qcmd::QueryCmd) = print_qcmd(io, qcmd)
print_init(io::IO, qcmd::QueryNot) = print_qcmd(io, qcmd)
Base.show(io::IO, reusable::ReusableQuery) = print(io, "Reusable ", reusable.query)
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

function Base.show(io::IO, ri::UnionInterval)
    print(io, ri.intervals[1])
    i = 1
    while i < length(ri.intervals)
        i += 1
        print(io, " ∪ ", ri.intervals[i])
    end
end

function Base.show(io::IO, ::MIME"text/plain", ri::UnionInterval)
    print(io, typeof(ri), ":\n ")
    print(io, ri)
end

function Base.show(io::IO, ri::RealInterval)
    @match ri.leftoperator begin
        (&<)    => print(io, "(", ri.lowerbound, ", ")
        (&<=)   => print(io, "[", ri.lowerbound, ", ")
        x       => print(io, ri.lowerbound, " ", x ," x ")
    end
    @match ri.rightoperator begin
        (&<)    => print(io, ri.upperbound, ")")
        (&<=)   => print(io, ri.upperbound, "]")
        x       => print(io, x, " ", ri.upperbound)
    end
end

function Base.show(io::IO, ::MIME"text/plain", ri::RealInterval)
    print(io, typeof(ri), ":\n ")
    print(io, ri)
end

Base.show(io::IO, ri::EmptyInterval) = 
    print(io, "∅")

Base.show(io::IO, ::MIME"text/plain", ri::EmptyInterval) = 
    print(io, typeof(ri), "()")

Base.show(io::IO, fln::FunctionalFunction) = 
    print(io, fln.fl, "(", fln.fn, ")")

function print_qcmd(io::IO, c::ComposedFunction)
    c.outer isa ComposedFunction ? print(io, c.outer) : _printcomposed(io, c.outer)
    print(io, " ∘ ")
    _printcomposed(io, c.inner)
end
print_qcmd(io::IO, x...) = print(io, x...)
print_qcmd(io::IO, mrm::MRM) = print(io, "MRM")
print_qcmd(io::IO, p::Pair) = (print_qcmd(io, p.first); print(io, " => "); print_qcmd(io, p.second))
function print_qcmd(io::IO, p::Pair{T, F}) where {T, F <: FunctionalFunction}
    print_qcmd(io, p.first => p.second.fn)
    print(io, "[")
    print_qcmd(io, p.second.fl)
    print(io, "]")
end
function print_qcmd(io::IO, p::Pair{T, F}) where {T, G <: Function, F <: ComposedFunction{<: Union{Base.Fix1, Base.Fix2, RealInterval}, G}}
    print_qcmd(io, p.second.inner ∘ p.first => p.second.outer)
end
function print_qcmd(io::IO, p::Pair{T, F}) where {T, F <: Base.Fix2}
    if Base.isoperator(repr(p.second.f))
        print(io, "(")
        print_qcmd(io, p.first)
        print(io, " ")
        print_qcmd(io, p.second.f)
        print(io, " ")
        print_qcmd(io, p.second.x)
        print(io, ")")
    elseif repr(p.second.f) == "in"
        print(io, "(")
        print_qcmd(io, p.first)
        print(io, " ")
        print_qcmd(io, "∈")
        print(io, " ")
        print_qcmd(io, p.second.x)
        print(io, ")")
    else
        print_qcmd(io, p.second.f)
        print_qcmd(io, "(")
        print_qcmd(io, p.first)
        print_qcmd(io, ", ")
        print_qcmd(io, p.second.x)
        print_qcmd(io, ")")
    end
end

function print_qcmd(io::IO, p::Pair{T, F}) where {T, F <: Base.Fix1}
    if Base.isoperator(repr(p.second.f))
        print(io, "(")
        print_qcmd(io, p.second.x)
        print(io, " ")
        print_qcmd(io, p.second.f)
        print(io, " ")
        print_qcmd(io, p.first)
        print(io, ")")
    else
        print_qcmd(io, p.second.f)
        print_qcmd(io, "(")
        print_qcmd(io, p.second.x)
        print_qcmd(io, ", ")
        print_qcmd(io, p.first)
        print_qcmd(io, ")")
    end
end

function print_qcmd(io::IO, p::Pair{T, F}) where {T, F <: RealInterval}
    print(io, "(")
    print(io, p.second.lowerbound)
    print(io, " ")
    print(io, repr(p.second.leftoperator))
    print(io, " ")
    print_qcmd(io, p.first)
    print(io, " ")
    print(io, repr(p.second.rightoperator))
    print(io, " ")
    print(io, p.second.upperbound)
    print(io, ")")
end
#shows !f instead of (!) ∘ f when ! is the outermost function
function print_qcmd(io::IO, c::ComposedFunction{typeof(!)})
    print(io, '!')
    _printcomposed(io, c.inner)
end

_printcomposed(io::IO, x) = print(io, x)
#display operators like + and - inside parens
_printcomposed(io::IO, f::Function) = Base.isoperator(Symbol(f)) ? (print(io, '('); print(io, f); print(io, ')')) : print(io, f)
#nesting for chained composition
_printcomposed(io::IO, f::ComposedFunction) = (print(io, '('); print_qcmd(io, f); print(io, ')'))
#no nesting when ! is the outer function in a composition chain
_printcomposed(io::IO, f::ComposedFunction{typeof(!)}) = print(io, f)