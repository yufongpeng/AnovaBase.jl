# AnovaFixedEffectModels
```@setup fem
using AnovaFixedEffectModels, DataFrames, CSV, CategoricalArrays, AnovaGLM
gpa = CSV.read("gpa.csv", DataFrame)
transform!(gpa,                                                                                                                   
           7 => x->replace(x, "yes" => true, "no" => false, "NA" => missing),                                                            
           5 => x->replace(x, "male" => 1, "female" => 0),                                                            
           4 => x->replace(x, "1 hour" => 1, "2 hours" => 2, "3 hours" => 3),                                                            
           renamecols = false)
transform!(gpa, [1, 2, 5, 7] .=> categorical, 4 => ByRow(Int), renamecols = false)
```
```@example fem
using AnovaFixedEffectModels
```
`AnovaFixedEffectModels.jl` supports [`FixedEffectModels`](https://github.com/FixedEffects/FixedEffectModels.jl). 

`lfe` is as same as `reg`, but the order of arguments is closer to other modeling packages.
```@example fem
fem1 = lfe(@formula(gpa ~ fe(student) + occasion + job), gpa)
aovf = anova(fem1)
```
Comparing between `FixedEffectModels`s and [`LinearModel`](https://juliastats.org/GLM.jl/stable/api/#GLM.LinearModel) is also available.
```@example fem
fems = nestedmodels(FixedEffectModel, @formula(gpa ~ fe(student) + occasion + job), gpa)
aovf = anova(fems)
```
In this case, `LinearModel` has to be the simplest model.
```@example fem
aovf = anova(lm(@formula(gpa ~ occasion + job), gpa), lfe(@formula(gpa ~ fe(student) + occasion + job), gpa))
```
Likelihood-ratio test is available for nested models.
```@example fem
fems = nestedmodels(FixedEffectModel, @formula(gpa ~ fe(student) + occasion + job), gpa)
anova(LRT, fems)
```