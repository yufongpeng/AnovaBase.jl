# AnovaFixedEffectModels
```@setup fem
using AnovaFixedEffectModels, DataFrames, CSV, CategoricalArrays
gpa = CSV.read("gpa.csv", DataFrame)
transform!(gpa,
        7 => x->replace(x, "yes" => true, "no" => false, "NA" => missing),
        4 => x->categorical(x, levels = ["1 hour", "2 hours", "3 hours"], ordered = true),
        renamecols = false)
transform!(gpa, [1, 2, 5, 7] .=> categorical, renamecols = false)
import AnovaBase: prednames
function prednames(aov::AnovaResult{T, FTest}; kwargs...) where {T <: StatsModels.TableRegressionModel{<: FixedEffectModel}}
    v = prednames(aov.model)
    push!(v, "(Residuals)")
    v
end
prednames(trm::StatsModels.TableRegressionModel{<: FixedEffectModel}) = AnovaFixedEffectModels.vectorize(prednames(formula(trm).rhs.terms[unique(trm.mm.assign)])
```
```@example fem
using AnovaFixedEffectModels
```
`AnovaFixedEffectModels.jl` supports [`FixedEffectModels`](https://github.com/FixedEffects/FixedEffectModels.jl); however, because `anova` relies on model schema, the output of `FixedEffectModels.reg` is not compatible. 

To solve this issue, fitting model using `lfe` instead of `reg`.
```@example fem
fem1 = lfe(@formula(gpa ~ fe(student) + occasion + job), gpa)
```
If a model is already fitted by `reg`, use `to_trm` to convert it into [`StatsModels.TableRegressionModel`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.TableRegressionModel).
```@example fem
model = reg(gpa, @formula(gpa ~ fe(student) + occasion + job))
fem1 = to_trm(model, gpa)
aovf = anova(fem1)
```
!!! note
    `lfe` is actually slower because it re-compiles every execution.
!!! note 
    Only F-test is available for `FixedEffectModel`.