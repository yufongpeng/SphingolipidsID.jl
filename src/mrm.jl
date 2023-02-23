function MRM(featuretable::Table, polarity::Bool = true; mz_tol = 0.35, additional = Dict())
    mz2v = Vector{Float64}[]
    mz2_loc = Vector{Int}[]
    for (id, mz2) in enumerate(featuretable.mz2)
        i = findfirst(x -> isapprox(mean(x), mz2; atol = 1e-4), mz2v)
        if isnothing(i)
            push!(mz2v, [mz2])
            push!(mz2_loc, [id])
        else
            push!(mz2v[i], mz2)
            push!(mz2_loc[i], id)
        end
    end
    mz2s = mean.(mz2v)
    raws = Table[]
    for (i, loc) in enumerate(mz2_loc)
        subft = @view featuretable[loc]
        push!(raws, Table(copy(subft), mz2_id = repeat([i], size(subft, 1)), mz2 = nothing))
    end
    MRM(reduce(append!, raws), mz2s, mz_tol, polarity, additional)
end