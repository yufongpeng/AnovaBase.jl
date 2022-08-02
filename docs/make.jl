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
            "AnovaGLM.md",
            "AnovaMixedModels.md",
            "AnovaFixedEffectModels.md"
        ],
        "API" => [  "AnovaBase" => "API\\AnovaBase.md", 
                    "AnovaGLM" => "API\\AnovaGLM.md",
                    "AnovaMixedModels" => "API\\AnovaMixedModels.md",
                    "AnovaFixedEffectModels" => "API\\AnovaFixedEffectModels.md"]
    ],
)

deploydocs(
    repo = "github.com/yufongpeng/AnovaBase.jl.git",
    push_preview = true, 
    devbranch = "main"
)
