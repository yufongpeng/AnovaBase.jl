push!(LOAD_PATH,"../src/")
using Documenter
using AnovaBase, GLM, AnovaGLM, AnovaMixedModels, AnovaFixedEffectModels, DataFrames
import AnovaBase: _diff, teststat, AnovaTable, pval
import AnovaMixedModels: GLM_MODEL, anovatable, AnovaResult, LRT, NestedModels
#DocMeta.setdocmeta!(AnovaBase, :DocTestSetup, :(using AnovaBase); recursive = true)
teststat(aov::AnovaResult) = aov.teststat
pval(aov::AnovaResult) = aov.pval
function anovatable(aov::AnovaResult{NestedModels{M, N}, LRT}; 
                    rownames = string.(1:N)) where {M <: Union{GLM_MODEL, MixedModel}, N}
    if last(aov.anovamodel.model) isa GLM_MODEL
        AnovaTable([
                dof(aov), 
                [NaN, _diff(dof(aov))...], 
                dof_residual(aov), 
                deviance(aov), 
                teststat(aov), 
                pval(aov)
            ],
            ["DOF", "ΔDOF", "Res.DOF", "Deviance", "χ²", "Pr(>|χ²|)"],
            rownames, 6, 5)
    elseif last(aov.anovamodel.model).optsum.REML 
        AnovaTable([
                        dof(aov), 
                        [NaN, _diff(dof(aov))...], 
                        dof_residual(aov), 
                        deviance(aov), 
                        teststat(aov), 
                        pval(aov)
                    ],
                ["DOF", "ΔDOF", "Res.DOF", "-2 logLik", "χ²", "Pr(>|χ²|)"],
                rownames, 6, 5)
    else
        AnovaTable([
                        dof(aov), 
                        [NaN, _diff(dof(aov))...], 
                        dof_residual(aov), 
                        aic.(aov.anovamodel.model),
                        bic.(aov.anovamodel.model),
                        deviance(aov), 
                        teststat(aov), 
                        pval(aov)
                    ],
                ["DOF", "ΔDOF", "Res.DOF", "AIC", "BIC", "-2 logLik", "χ²", "Pr(>|χ²|)"],
                rownames, 8, 7)
    end
end

makedocs(
    sitename = "AnovaBase",
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
    ]
)

deploydocs(
    repo = "github.com/yufongpeng/AnovaBase.jl.git",
    push_preview = true, 
    devbranch = "dev"
)
