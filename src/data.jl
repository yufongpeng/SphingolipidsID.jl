function featuretable_mzmine(paths)
    tbl = CSV.read(paths, DataFrame; select = (i, name) -> any(startswith(String(name), text) for text in ["id", "rt", "mz", "height", "area", "intensity"]))
    fwhm = CSV.read(paths, DataFrame; select = (i, name) -> endswith(String(name), "fwhm"))
    id = findfirst.(!ismissing, eachrow(fwhm))
    tbl.FWHM = getindex.(eachrow(fwhm), id)
    datafile = Dict(propertynames(fwhm) .=> map(col -> match(r".*:(.*):.*", col)[1], names(fwhm)))
    tbl.datafile = getindex.(Ref(datafile), id)
    rename!(tbl, :mz => :mz1)
    sort!(tbl, :datafile)
    println("DataFiles Order: ")
    for i in unique(tbl.datafile)
        println(" ", i)
    end
    tbl
end

function fill_ce_mzmine!(tbl, eV::Float64)
    tbl.collision_energy .= eV
    select!(tbl, Not(:datafile))
    sort!(tbl, [:mz1, :rt])
    tbl.id = 1:size(tbl, 1)
    tbl
end

function fill_ce_mzmine!(tbl, eV)
    mapping = Dict(unique(tbl.datafile) .=> eV)
    tbl.collision_energy = getindex.(Ref(mapping), tbl.datafile)
    select!(tbl, Not(:datafile))
    sort!(tbl, [:mz1, :rt])
    tbl.id = 1:size(tbl, 1)
    tbl
end

function featuretable_masshunter_mrm(path)
    strs = readlines(path)
    starts = Int[]
    data = Vector[]
    ends = Int[]
    status = false
    for (i, l) in enumerate(strs)
        if status
            if occursin("No integration results available", l)
                pop!(starts)
                pop!(data)
                status = false
            elseif l == ""
                push!(ends, i - 1)
                status = false
            end
        else
            transition = match(r"\((\d*.\d*) -> (\d*.\d*)\)", l)
            isnothing(transition) && continue
            ms1, ms2 = parse.(Float64, transition)
            eV = parse(Float64, match(r"CID@(\d*.\d*)", l)[1])
            push!(data, [ms1, ms2, eV])
            push!(starts, i + 1)
            status = true
        end
    end
    str = map(starts, ends) do st, ed
        IOBuffer(join(strs[st:ed], "\n"))
    end
    txt = ["Start", "End", "RT", "Height", "Area", "Symmetry", "Y", "FWHM"]
    tbl = CSV.read(str, DataFrame; select = (i, name) -> any(occursin(text, String(name)) for text in txt) && all(!occursin(text, String(name)) for text in ["B1, B2"]))
    rename!(tbl,
        "RT"        => "rt",
        "Start"     => "rt_range:min",
        "End"       => "rt_range:max",
        "Height"    => "height",
        "Area"      => "area",
        "Start Y"   => "intensity_range:start",
        "End Y"     => "intensity_range:end",
        "Symmetry"  => "symmetry",
        "Start BL Y"=> "baseline:start",
        "End BL Y"  => "baseline:end",
        "Max Y"     => "intensity_range:max"
        )
    datas = eachcol(mapreduce(vcat, zip(data, starts, ends)) do (dt, st, ed)
        repeat(dt', ed - st)
    end)
    insertcols!(tbl, ([:mz1, :mz2, :collision_energy] .=> datas)...)
    sort!(tbl, [:mz1, :rt])
    tbl.id = 1:size(tbl, 1)
    tbl
end

rsd(v) = std(v) / mean(v)
re(v) =  - foldl(-, extrema(v)) / mean(v) / 2
default_error(v) = length(v) > 2 ? rsd(v) : re(v)

filter_duplicate!(tbls::Vector{DataFrame}; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5) = 
    (tbls[:] = filter_duplicate.(tbls; rt_tol, mz_tol, n, err, err_tol))

function filter_duplicate(tbl::DataFrame; rt_tol = 0.1, mz_tol = 0.35, n = 3, err = default_error, err_tol = 0.5)
    ids = Vector{Int}[]
    for (i, ft) in enumerate(eachrow(tbl))
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
    filter!(:id => >(0), tbl)
    combine(groupby(tbl, :id), All() .=> mean, :area => err => :error, renamecols = false)
end

function CompoundSP(project::Project, cpd, product, source, id, area)
    sum_cb, sum_db, sum_o = match(r".+ (\d+):(\d+);(\d*)O", cpd.Species).captures
    sum_cb = parse(Int, sum_cb)
    sum_db = parse(Int, sum_db)
    sum_o = isempty(sum_o) ? 1 : parse(Int, sum_o)
    class = (cpd.Abbreviation)()
    adduct_class = object_adduct(cpd.Adduct)
    chain = nothing
    fragments = DataFrame(ion1 = Ion[Ion(adduct_class, class)], ion2 = Ion[product], source = [source], id = [id])
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
                    fragments.ion2[1] = Ion(add_new, lcb_new) 
                end
            end
            chain = Chain(fragments.ion2[1].molecule, Acyl{acyl_o}())
        end
        x && if x in NANA end       => begin
            hasnana(class) || return nothing
        end
    end
    CompoundSP(class, (sum_cb, sum_db, sum_o), chain, 
        fragments, area,
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