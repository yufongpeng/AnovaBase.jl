# AnovaMixedModels
```@setup mm
using AnovaMixedModels, DataFrames, GLM, CSV, CategoricalArrays
anxiety = CSV.read("anxiety.csv", DataFrame)
transform!(anxiety, :id => categorical, renamecols = false)
toenail = CSV.read("toenail.csv", DataFrame)
transform!(toenail, [1, 2, 3] .=> categorical, renamecols = false)
```
The implementation of ANOVA for [mixed-effects models](https://en.wikipedia.org/wiki/Mixed_model) is primarily based on [`MixedModels`](https://juliastats.org/MixedModels.jl/stable/). The syntax is similar to anova for `GLM`.   
```@example mm
using AnovaMixedModels
```
## Linear mixed-effects model
We get a dataset from `R` directly by [`RCall`](https://juliainterop.github.io/RCall.jl/stable/).
```julia
R"""data("anxiety", package = "datarium")"""
anxiety = stack(rcopy(R"anxiety"), [:t1, :t2, :t3], [:id, :group], variable_name = :time, value_name = :score)
anxiety = combine(anxiety, Not(:time), :time => ByRow(x->parse(Int, replace(String(x), "t"=>""))) => :time)
```
We can fit a linear mixed-effects model first. 
```@example mm
lmm1 = lmm(@formula(score ~ group * time + (1|id)), anxiety)
anova(lmm1)
```
Alternatively, we can use `anova_lmm`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` fitted by [Restricted maximum likelihood](https://en.wikipedia.org/wiki/Restricted_maximum_likelihood).
```@example mm
aov = anova_lmm(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
```
```@example mm
aov.anovamodel.model.optsum.REML
```
For [likeihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test), all submodels are fitted. The model should be fitted by [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation).
```@example mm
anova(LRT, lmm1)
```
When comparing multiple mixed models, likelihood-ratio test is used by default. 
It's also identical to [`StatsModels.lrtest`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.lrtest) and [`MixedModels.likelihoodratiotest`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LikelihoodRatioTest).
```@example mm
lmms = nestedmodels(lmm1)
anova(lmms) # == anova(LRT, lmm1)
```
```@example mm
MixedModels.likelihoodratiotest(lmms.model[2:end]...)
``` 
Comparing between [`LinearModel`](https://juliastats.org/GLM.jl/stable/api/#GLM.LinearModel) and [`LinearMixedModel`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LinearMixedModel) is also available.
```@example mm
lm1 = lm(@formula(score ~ group * time), anxiety)
lmm2 = lmm(@formula(score ~ group * time + (group|id)), anxiety)
anova(lm1, lmm1, lmm2)
```
## Generalized linear mixed-effects model
The following is an example of generalized mixed model.
```julia
R"""data("toenail", package = "HSAUR2")"""
toenail = rcopy(R"toenail")
```
```@example mm
glmm1 = glmm(@formula(outcome ~ visit + treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
glmm2 = glmm(@formula(outcome ~ visit * treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
anova(glmm1, glmm2)
```
!!! note
    Only likelihood-ratio test is available now for `GeneralizedLinearMixedModel`.