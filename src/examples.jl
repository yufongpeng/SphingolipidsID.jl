using IDGSL
files = joinpath.(".\\data", readdir(".\\data"))
fts = map(files) do file
    filter_duplicate(featuretable_mzmine(file))
end
ms = [236.238, 250.253, 284.295, 264.269, 262.253, 278.285, 292.3, 266.285, 274.093, 282.28]
range = vcat(repeat([(400, 1500)], 8), [(900, 1650), (400, 1500)])

db = reduce(append!, (
           library([Cer], ["[M+H]+"], (32:46, 1:2, ["O"])),
           library([Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
           library([HexCer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([HexCer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
           library([Hex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([Hex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
           library([Hex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([Hex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
           library([HexNAcHex2Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([HexNAcHex3Cer], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([GM3], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:4, ["2O"])),
           library([GM3], ["[M+H]+", "[M+H-H2O]+"], (32:46, 0:3, ["3O"])),
           library([GM1], ["[M+H]+", "[M+2H]2+"], (34:2:42, 1:2, ["2O"]))
           )
        )

pj = preis()
for (m, r, ft) in zip(ms, range, fts)
    preis!(pj, ft, r, m, true; db = db, rt_tol = 0.1)
end

finish_profile!(pj)

apply_rules!(pj, :class)

pj2 = deepcopy(pj)
pj2 = Project(deepcopy(query(pj, :mw, 667, 668)[2:2]), deepcopy(pj.data), :anion)
analyte = @chain pj begin
    query(:class, GM3)
    query(:sum, (40, 2, 2))
    _.result[end-1]
end
pj2 = Project(deepcopy([analyte]), deepcopy(pj.data), :acetate)
pj2 = Project(deepcopy(query(pj, :rt, 9.59, 9.61).result), deepcopy(pj.data), :acetate)

apply_rules!(pj, :chain)

pj2 = Project(deepcopy(query(pj, :mz, 1180, 1181)[5:5]), deepcopy(pj.data), :anion)

@chain pj begin
    query(:class, GM3)
    query(:sum, (36,1,2))
end