push!(LOAD_PATH,"../src/")
using Documenter
using AnovaBase, GLM, AnovaGLM, AnovaMixedModels, AnovaFixedEffectModels, DataFrames

makedocs(
    sitename = "AnovaBase",
    pages = [
        "index.md",
        "Examples" => [
            "Examples_GLM.md",
            "Examples_MixedModels.md",
            "Examples_FixedEffectModels.md"
        ],
        "Interfacing AnovaBase.jl" => [
            "Interface.md"
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
    ]
)

deploydocs(
    repo = "github.com/yufongpeng/AnovaBase.jl.git",
    push_preview = true, 
    devbranch = "dev"
)
