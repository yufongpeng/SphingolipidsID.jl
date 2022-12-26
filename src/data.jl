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

sort_data(tbl::Table) = sort_data!(tbl)
function sort_data!(tbl::Table)
    sort!(tbl, :datafile)
    for i in unique(tbl.datafile)
        println(" ", i)
    end
    println()
    tbl
end
function fill_mz2!(tbl::Table, mz2::Float64)
    fill!(tbl.mz2, mz2)
    tbl
end

function fill_mz2!(tbl::Table{
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile), 
        Tuple{Int64, Float64, Float64, Float64, Float64, Float64, Int64, Float64, Float64, SubString{String}}}, 1, 
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile), 
        Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{SubString{String}}}}
    }, 
    mz2::Union{<: Vector, <: Tuple})
    mapping = Dict(unique(tbl.datafile) .=> mz2)
    tbl.mz2 .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

function fill_ce!(tbl::Table, eV::Float64)
    fill!(tbl.collision_energy, eV)
    tbl
end

function fill_ce!(tbl::Table{
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile), 
        Tuple{Int64, Float64, Float64, Float64, Float64, Float64, Int64, Float64, Float64, SubString{String}}}, 1, 
    NamedTuple{
        (:id, :mz1, :mz2, :rt, :height, :area, :collision_energy, :FWHM, :symmetry, :datafile), 
        Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Float64}, Vector{Float64}, Vector{SubString{String}}}}
    }, 
    eV::Union{<: Vector, <: Tuple})
    mapping = Dict(unique(tbl.datafile) .=> eV)
    tbl.collision_energy .= getindex.(Ref(mapping), tbl.datafile)
    tbl
end

fill_mz2!(tbl, mz2) = tbl
fill_ce!(tbl, eV) = tbl

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

rsd(v) = std(v) / mean(v)
re(v) =  - foldl(-, extrema(v)) / mean(v) / 2
default_error(v) = length(v) > 2 ? rsd(v) : re(v)

function filter_duplicate!(tbl::Table; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5)
    tbl = Table(tbl; datafile = nothing)
    sort!(tbl, [:mz2, :mz1, :rt])
    fill!(tbl.id, 0)
    locs = Vector{Int}[]
    for i in eachindex(tbl)
        new = true
        rt = tbl.rt[i]
        mz1 = tbl.mz1[i]
        mz2 = tbl.mz2[i]
        collision_energy = tbl.collision_energy[i]
        for loc in locs
            abs(mean(tbl.rt[loc]) - rt) > rt_tol && continue
            abs(mean(tbl.mz1[loc]) - mz1) > mz_tol && continue
            abs(mean(tbl.mz2[loc]) - mz2) > mz_tol && continue
            tbl.collision_energy[loc[1]] == collision_energy || continue
            push!(loc, i)
            new = false
            break
        end
        new && push!(locs, [i])
    end
    n > 1 && filter!(loc -> (length(loc) >= n && err(tbl.area[loc]) <= err_tol), locs)
    for (i, loc) in enumerate(locs)
        tbl.id[loc] .= i
    end
    @p tbl |> filter!(>(_.id, 0))
    @p tbl |> DataFrame |> groupby(__, :id) |> combine(__, All() .=> mean, :area => err => :error, renamecols = false) |> Table
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

AnalyteSP(compounds::Vector{CompoundSP}, rt::Float64, states::Vector{Int}) = AnalyteSP(compounds, rt, states, (NaN, NaN), 0)

function id_product(ms2, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(db)
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end

function read_mrm(path::String; vendor = :agilent)
    tbl = CSV.read(path, Table)
    if vendor == :agilent
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

union_mrm(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = union_mrm!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
function union_mrm!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    for id2 in eachindex(tbl2)
        pushed = false
        ms1 = tbl2.mz1[id2]
        ms2 = tbl2.mz2[id2]
        ce = tbl2.collision_energy[id2]
        for id1 in eachindex(tbl1)
            (!between(tbl1.mz1[id1], ms1, mz_tol) || !between(tbl1.mz2[id1], ms2, mz_tol)) && continue
            tbl1.collision_energy[id1] == ce || continue
            rt_l = min(tbl1.rt[id1] - tbl1.Δrt[id1] / 2, tbl2.rt[id2] - tbl2.Δrt[id2] / 2)
            rt_r = max(tbl1.rt[id1] + tbl1.Δrt[id1] / 2, tbl2.rt[id2] + tbl2.Δrt[id2] / 2)
            n = length(vectorize(tbl1.compound[id1]))
            tbl1[id1] = (;
                compound = vectorize(tbl1.compound[id1]),
                mz1 = (tbl1.mz1[id1] * n + ms1) / (n + 1),
                mz2 = (tbl1.mz2[id1] * n + ms2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = tbl1.collision_energy[id1],
                polarity = tbl1.polarity[id1]
            )
            union!(tbl1.compound[id1], vectorize(tbl2.compound[id2]))
            pushed = true
        end
        pushed || push!(tbl1, (compound = tbl2.compound[id2], mz1 = ms1, mz2 = ms2, rt = tbl2.rt[id2], Δrt = rt_tol * 2, collision_energy = ce, polarity = tbl2.polarity[id2]))
    end
    sort!(tbl1, [:mz1, :rt])
end

diff_mrm(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5) = diff_mrm!(deepcopy(tbl1), tbl2; mz_tol, rt_tol)
function diff_mrm!(tbl1::Table, tbl2::Table; mz_tol = 0.35, rt_tol = 0.5)
    sort!(tbl1, [:mz1, :rt])
    extended = Int[]
    l = size(tbl1, 1)
    for id2 in eachindex(tbl2)
        pushed = false
        ms1 = tbl2.mz1[id2]
        ms2 = tbl2.mz2[id2]
        ce = tbl2.collision_energy[id2]
        for id1 in eachindex(tbl1)
            (!between(tbl1.mz1[id1], ms1, mz_tol) || !between(tbl1.mz2[id1], ms2, mz_tol)) && continue
            tbl1.collision_energy[id1] == ce || continue
            rt_l = min(tbl1.rt[id1] - tbl1.Δrt[id1] / 2, tbl2.rt[id2] - tbl2.Δrt[id2] / 2)
            rt_r = max(tbl1.rt[id1] + tbl1.Δrt[id1] / 2, tbl2.rt[id2] + tbl2.Δrt[id2] / 2)
            n = length(vectorize(tbl1.compound[id1]))
            isempty(setdiff(vectorize(tbl2.compound[id2]), vectorize(tbl1.compound[id1]))) || push!(extended, id1)
            tbl1[id1] = (;
                compound = vectorize(tbl1.compound[id1]),
                mz1 = (tbl1.mz1[id1] * n + ms1) / (n + 1),
                mz2 = (tbl1.mz2[id1] * n + ms2) / (n + 1),
                rt = (rt_l + rt_r) / 2,
                Δrt = rt_r - rt_l,
                collision_energy = tbl1.collision_energy[id1],
                polarity = tbl1.polarity[id1]
            )
            union!(tbl1.compound[id1], vectorize(tbl2.compound[id2]))
            pushed = true
        end
        pushed && continue
        push!(new, id2)
        push!(tbl1, (compound = tbl2.compound[id2], mz1 = ms1, mz2 = ms2, rt = tbl2.rt[id2], Δrt = rt_tol * 2, collision_energy = ce, polarity = tbl2.polarity[id2]))
    end
    (extended = tbl1[extended], new = size(tbl1, 1) > l ? tbl1[l + 1:end] : nothing)
end

function write_mrm(io, tbl::Table; vendor = :agilent) 
    if vendor == :agilent
        CSV.write(io, tbl; 
            header = ["Compound Name", "Precursor Ion", "Product Ion", "Ret Time (min)", "Delta Ret Time", "Collision Energy", "Polarity"])
    end
end
