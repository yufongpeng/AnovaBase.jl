using Pkg, Documenter
cd("..\\..")
Pkg.activate("AnovaBase")
using AnovaBase
DocMeta.setdocmeta!(AnovaBase, :DocTestSetup, :(using AnovaBase); recursive = true)

makedocs(
    modules=[AnovaBase],
    sitename = "AnovaBase",
    doctest = false,
    pages = [
        "index.md",
        "GLM.md",
        "MixedModels.md",
        "FixedEffectModels.md",
        "API.md"
    ],
)

deploydocs(
    repo = "github.com/yufongpeng/AnovaBase.jl.git",
    push_preview = true, 
    devbranch = "master"
)
