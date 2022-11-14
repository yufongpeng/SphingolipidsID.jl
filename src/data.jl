function featuretable_mzmine(paths)
    tbl = CSV.read(paths, DataFrame; select = (i, name) -> any(startswith(String(name), text) for text in ["id", "rt", "mz", "height", "area", "intensity"]))
    fwhm = CSV.read(paths, DataFrame; select = (i, name) -> endswith(String(name), "fwhm"))
    id = map(eachrow(fwhm)) do row
        findfirst(!ismissing, row)
    end
    tbl.FWHM = [row[i] for (row, i) in zip(eachrow(fwhm), id)]
    datafile = Dict(propertynames(fwhm) .=> map(names(fwhm)) do col
        match(r".*:(.*):.*", col)[1]
    end)
    tbl.datafile = [datafile[i] for i in id]
    rename!(tbl, :mz => :mz1)
    sort!(tbl, [:mz1, :rt])
    tbl.id = 1:size(tbl, 1)
    println("DataFiles Order: ", unique(tbl.datafile)...)
    tbl
end

function add_ce_mzmine!(tbl, eV::Float64)
    tbl.collision_energy .= eV
    delete!(tbl, :datafile)
    tbl
end

function add_ce_mzmine!(tbl, eV)
    mapping = Dict(unique(tbl.datafile) .=> eV)
    tbl.collision_energy = getindex.(Ref(mapping), tbl.datafile)
    select!(tbl, Not(:datafile))
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

function filter_duplicate(tbl::DataFrame; rt_tol = 0.1, mz_tol = 0.35)
    ids = Vector{Int}[]
    for (i, ft) in enumerate(eachrow(tbl))
        new = true
        for id in ids
            if abs(mean(tbl[id, :rt]) - ft.rt) < rt_tol && abs(mean(tbl[id, :mz1]) - ft.mz1) < mz_tol && tbl[id[1], :collision_energy] == ft.collision_energy 
                push!(id, i)
                new = false
                break
            end
        end
        new && push!(ids, [i])
    end
    for (i, id) in enumerate(ids)
        tbl[id, :id] .= i
    end
    combine(groupby(tbl, :id), All() .=> mean, renamecols = false)
end

function CompoundGSL(cpd, product, source, id, area)
    sum_cb, sum_db, sum_o = match(r".+ (\d+):(\d+);(\d*)O", cpd.Species).captures
    sum_cb = parse(Int, sum_cb)
    sum_db = parse(Int, sum_db)
    sum_o = isempty(sum_o) ? 1 : parse(Int, sum_o)
    class = (cpd.Abbreviation)()
    adduct_class = object_adduct(cpd.Adduct)
    chain = nothing
    fragments = DataFrame(ion1 = Ion[Ion(adduct_class, class)], ion2 = Ion[product], source = [source], id = [id])
    result = @match product begin
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
                    fragments[1, :ion2] = Ion(add_new, lcb_new) 
                end
            end
            chain = Chain(fragments[1, :ion2].molecule, Acyl{acyl_o}())
        end
        x && if x in NANA end       => begin
            hasnana(class) || return nothing
        end
    end
    CompoundGSL(class, (sum_cb, sum_db, sum_o), chain, 
        fragments, area,
        Union{Missing, Pair, IonMode, IonUnion, Symbol, Nothing}[missing, missing]
    )
end

function id_product(ms2, polarity; db = SPDB[polarity ? :FRAGMENT_POS : :FRAGMENT_NEG], mz_tol = 0.35)
    products = Ion[]
    for row in eachrow(db)
        between(ms2, row[2], mz_tol) && push!(products, row[1])
    end
    products
end