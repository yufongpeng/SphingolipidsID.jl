function featuretable_mzmine(path)
    head = split(readline(path), ",")
    idmain = findall(x -> any(startswith(x, text) for text in ["id", "rt", "mz", "height", "area", "intensity"]), head)
    idfwhm = findall(x-> endswith(x, "fwhm"), head)
    tbl = CSV.read(path, Table; select = idmain)
    head = [:mz, :rt, :height, :area, Symbol("mz_range:min"), Symbol("mz_range:max"), Symbol("rt_range:min"), Symbol("rt_range:max"), Symbol("intensity_range:min"), Symbol("intensity_range:max")]
    tbl = Table(; (head .=> getproperty.(Ref(tbl), head))...)
    n = size(tbl, 1)
    fwhm = CSV.read(path, Table; select = idfwhm)
    id = findfirst.(!ismissing, fwhm)
    datafile = Dict(propertynames(fwhm) .=> map(col -> match(r".*:(.*):.*", string(col))[1], propertynames(fwhm)))
    tbl = Table(Table(id = zeros(Int, n), mz1 = tbl.mz, scan = ones(Int, n)), tbl, collision_energy = zeros(Int, n), FWHM = zeros(Float64, n), datafile = getindex.(Ref(datafile), id); mz = nothing)
    sort!(tbl, :datafile)
    println("DataFiles Order: ")
    for i in unique(tbl.datafile)
        println(" ", i)
    end
    tbl
end

function fill_ce_mzmine!(tbl, eV::Float64)
    tbl.collision_energy .= eV
    sort!(tbl, [:mz1, :rt])
    tbl.id .= 1:size(tbl, 1)
    Table(tbl; datafile = nothing)
end

function fill_ce_mzmine!(tbl, eV)
    mapping = Dict(unique(tbl.datafile) .=> eV)
    tbl.collision_energy .= getindex.(Ref(mapping), tbl.datafile)
    sort!(tbl, [:mz1, :rt])
    tbl.id .= 1:size(tbl, 1)
    Table(tbl; datafile = nothing)
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
        else
            #if occursin("No integration results available", l)
            pop!(starts)
            pop!(data)
            status = false
        end
    end
    str = map(starts, ends) do st, ed
        IOBuffer(join(strs[st:ed], "\n"))
    end
    rep = ends .- starts
    txt = ["Start", "End", "RT", "Height", "Area", "Symmetry", "Y", "FWHM"]
    tbl = CSV.read(str, Table; select = (i, name) -> any(occursin(text, String(name)) for text in txt) && all(!occursin(text, String(name)) for text in ["B1, B2"]))
    n = size(tbl, 1)
    del = Symbol.(["RT", "Height", "Area", "Start", "End", "Start Y", "End Y", "Symmetry", "Start BL Y", "End BL Y", "Max Y"])
    tbl = Table(Table(
                id = zeros(Int, n), 
                mz1 = (@p zip(data.ms1, rep) |> mapmany(repeat([_[1]], _[2]))),
                mz2 = (@p zip(data.ms2, rep) |> mapmany(repeat([_[1]], _[2]))),
                collision_energy = (@p zip(data.eV, rep) |> mapmany(repeat([_[1]], _[2])))
            ),
            (; (Symbol.(["rt", "height", "area", "rt_range:min", "rt_range:max", 
                        "intensity_range:start", "intensity_range:end", 
                        "symmetry", "baseline:start", "baseline:end", 
                        "intensity_range:max"]) .=> getproperty.(Ref(tbl), del))...
            ), 
            tbl; (del .=> nothing)...)
    sort!(tbl, [:mz1, :rt])
    tbl.id .= 1:n
    tbl
end

rsd(v) = std(v) / mean(v)
re(v) =  - foldl(-, extrema(v)) / mean(v) / 2
default_error(v) = length(v) > 2 ? rsd(v) : re(v)

function filter_duplicate!(tbl::Table; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5)
    ids = Vector{Int}[]
    for (i, ft) in enumerate(tbl)
        new = true
        for id in ids
            if abs(mean(tbl.rt[id]) - ft.rt) < rt_tol && abs(mean(tbl.mz1[id]) - ft.mz1) < mz_tol && tbl.collision_energy[id[1]] == ft.collision_energy 
                push!(id, i)
                new = false
                break
            end
        end
        new && push!(ids, [i])
    end
    n > 1 && filter!(id -> (length(id) >= n && err(tbl.area[id]) <= err_tol), ids)
    fill!(tbl.id, 0)
    for (i, id) in enumerate(ids)
        tbl.id[id] .= i
    end
    @p tbl |> filter!(>(_.id, 0))
    @p tbl |> DataFrame |> groupby(__, [:id, :scan]) |> combine(__, All() .=> mean, :area => err => :error, renamecols = false) |> Table
    #=
    gtbl = @p tbl |> groupview(getproperty(:id))
    tbl = @p gtbl |> map(map(mean, columns(_))) |> Table
    errors = @p gtbl |> map(err(_.area)) |> collect
    tbl = Table(tbl, id = Int.(tbl.id), collision_energy = Int.(tbl.collision_energy), error = errors)
    =#
end

function CompoundSP(project::Project, cpd, product, source, id, area)
    sum_cb, sum_db, sum_o = match(r".+ (\d+):(\d+);(\d*)O", cpd.Species).captures
    sum_cb = parse(Int, sum_cb)
    sum_db = parse(Int, sum_db)
    sum_o = isempty(sum_o) ? 1 : parse(Int, sum_o)
    class = (cpd.Abbreviation)()
    adduct_class = object_adduct(cpd.Adduct)
    chain = nothing
    ion2 = Ion[product]
    @match product begin
        Ion(adduct, lcb::LCB) => begin 
            lcb_o = nhydroxyl(lcb)
            lcb_db = ndb(lcb)
            #=
            1. sum_db + sum_o < lcb_db + lcb_o (total unsaturation deficiency)
            2. sum_o < lcb_o 
                0 < Δ = lcb_o - sum_o ≤ sum_db - lcb_db
                lcb_o' = lcb_o - Δ = sum_o => sum_o - lcb_o' = 0
                lcb_db' = lcb_db + Δ ≤ sum_db 
            3. sum_db < lcb_db
                0 < Δ = lcb_db - sum_db ≤ sum_o - lcb_o
                lcb_db' = lcb_db - Δ = sum_db 
                lcb_o' = lcb_o + Δ ≤ sum_o => sum_o - lcb_o' = sum_o - lcb_o - Δ ≥ 0
            =#
            if (sum_db + sum_o) < (lcb_o + lcb_db)
                return nothing
            end
            Δ = sum_o - lcb_o 
            Δ = Δ < 0 ? Δ : max(0, lcb_db - sum_db)
            acyl_o = sum_o - lcb_o - Δ
            if Δ != 0
                lcb_new = hydroxyl_shift(lcb, Δ)
                add_new = hydroxyl_shift(product.adduct, Δ)
                if isnothing(lcb_new) || isnothing(add_new)
                    return nothing
                else
                    ion2[1] = Ion(add_new, lcb_new) 
                end
            end
            chain = Chain(ion2[1].molecule, Acyl{acyl_o}())
        end
        x && if x in NANA end       => begin
            hasnana(class) || return nothing
        end
    end
    CompoundSP(class, (sum_cb, sum_db, sum_o), chain, 
        Table(ion1 = Ion[Ion(adduct_class, class)], ion2 = ion2, source = [source], id = [id]), area,
        [0, 0],
        RuleSet[RuleSet(:missing, EmptyRule()), RuleSet(:missing, EmptyRule())],
        project
    )
end

AnalyteSP(compounds::Vector{CompoundSP}, rt::Float64, states::Vector{Int}) = AnalyteSP(compounds, rt, states, (NaN, NaN))

function id_product(ms2, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(db)
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end