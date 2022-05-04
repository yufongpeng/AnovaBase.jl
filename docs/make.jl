using Pkg, Documenter, MixedAnova
DocMeta.setdocmeta!(MixedAnova, :DocTestSetup, :(using MixedAnova); recursive = true)
glm_init()
mm_init()
fem_init()

makedocs(
    modules=[MixedAnova],
    sitename = "MixedAnova",
    doctest = true,
    pages = [
        "index.md",
        "GLM.md",
        "MixedModels.md",
        "FixedEffectModels.md",
        "API.md"
    ],
)

deploydocs(
    repo = "github.com/yufongpeng/MixedAnova.jl.git",
    push_preview = true
)
