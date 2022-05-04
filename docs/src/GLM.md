# GLM.jl
```@setup glm
using MixedAnova, Pkg
glm_init()
Pkg.activate(joinpath("..", "..", "test"))
using RDatasets, DataFrames
iris = dataset("datasets", "iris")
quine = dataset("MASS", "quine")
mtcars = dataset("datasets", "mtcars")
Pkg.activate(joinpath("..", ".."))
```
To use `anova` on [`GLM objects`](https://juliastats.org/GLM.jl/stable/) , we need to initiate some files first.
```julia
glm_init()
```
This function will export all functions from [`GLM`](https://juliastats.org/GLM.jl/stable/) and related function in this package, including `anova`, `anova_lm`, `anova_glm`.
## Ordinary linear model
We first import the well-known [`iris` dataset](https://en.wikipedia.org/wiki/Iris_flower_data_set) from [`RDatasets`](https://github.com/JuliaStats/RDatasets.jl).
```julia
iris = dataset("datasets", "iris")
```
There's two way to perform ANOVA. `anova_lm` accepts a formula and data like [`GLM.lm`](https://juliastats.org/GLM.jl/stable/api/#GLM.lm).
```@example glm
anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris)
```
We can specify the type of sum of squares by keyword argument `type`. Let's use type II SS.
```@example glm
anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris, type = 2)
```
A [`StatsModels.TableRegressionModel`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.TableRegressionModel) object is fitted and stored in the output of `anova_lm`.  

We can fit a model first and call `anova` instead. `anova` store the model as well.
!!! warn
    It doesn't create a copy, so any in-place change of the original model should be noticed. 
```@example glm
lm1 = lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris)
anova(lm1, type = 3)
```
Multiple models can be compared through the same function.  
!!! note
    The checker for nested models is not implemented now, so it should be ensured that the later model is more saturated than the previous one.  
```@example glm
lms = nestedmodels(LinearModel, @formula(SepalLength ~ SepalWidth * Species), iris, dropcollinear = false)
anova(lms...)
```
This result is a little bit different from [`GLM.ftest`](https://juliastats.org/GLM.jl/stable/api/#GLM.ftest):
```@example glm
ftest(getproperty.(lms[2:end], :model)...)
```
In `anova`, the F value is calculated by dividing [MSR](https://en.wikipedia.org/wiki/Mean_squared_error) (mean of ΔDeviance) with mean of [RSS](https://en.wikipedia.org/wiki/Residual_sum_of_squares) of saturated model just like [`anova` in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova), while in [`GLM.ftest`](https://juliastats.org/GLM.jl/stable/api/#GLM.ftest), the denominator is replaced by [RSS](https://en.wikipedia.org/wiki/Residual_sum_of_squares) of subsequant model.
## Generalized linear models 
```julia
quine = dataset("MASS", "quine")
```
We fit a [negative binomial regression](https://en.wikipedia.org/wiki/Generalized_linear_model) on [`quine` dataset](https://www.rdocumentation.org/packages/MASS/versions/7.3-57/topics/quine) from [`MASS`](https://www.rdocumentation.org/packages/MASS/versions/7.3-57).
```@example glm
nbm = glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink())
anova(nbm)
```
There's also `anova_glm` similar to `anova_lm`.  

`anova` will automatically select test from [F-test](https://en.wikipedia.org/wiki/F-test) or [likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test) depending on the type of [distribution](https://juliastats.org/GLM.jl/stable/#Fitting-GLM-models). For distribution of `Bernoulli()`, `Binomial()`, `Poisson()` that have fixed dispersion, likelihood-ratio test is preferred. For other distribution, F-test is preferred.  

The next one is an axample of [logistic regression](https://en.wikipedia.org/wiki/Logistic_regression).
```julia
mtcars = dataset("datasets", "mtcars")
```
We want to predict if the `AM` is 0 or 1. Let's use logistic regression with and without interaction terms, and compare this two models by likelihood-ratio test. 
```@example glm
glm1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink())
glm2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), LogitLink())
anova(glm1, glm2)
```
```@example glm
lrtest(glm1, glm2)
```
This function works identically as [`StatsModels.lrtest`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.lrtest).
!!! note
    We can also specify test by keword arguments `test` or putting test in the first argument.