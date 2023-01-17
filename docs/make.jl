push!(LOAD_PATH,"../src/")
using Pkg, Documenter
using AnovaBase, GLM, AnovaGLM, AnovaMixedModels, AnovaFixedEffectModels, DataFrames
#DocMeta.setdocmeta!(AnovaBase, :DocTestSetup, :(using AnovaBase); recursive = true)

makedocs(
    modules=[AnovaBase, AnovaGLM, AnovaMixedModels, AnovaFixedEffectModels],
    sitename = "AnovaBase",
    doctest = false,
    pages = [
        "index.md",
        "Examples" => [
            "Examples_GLM.md",
            "Examples_MixedModels.md",
            "Examples_FixedEffectModels.md"
        ],
        "Algorithm" => [
            "Algorithm_AnovaGLM.md",
            "Algorithm_AnovaMixedModels.md",
            "Algorithm_AnovaFixedEffectModels.md"
        ],
        "API" => [  
            "AnovaBase.md", 
            "AnovaGLM.md",
            "AnovaMixedModels.md",
            "AnovaFixedEffectModels.md"
        ]
    ],
)

deploydocs(
    repo = "github.com/yufongpeng/AnovaBase.jl.git",
    push_preview = true, 
    devbranch = "dev"
)
