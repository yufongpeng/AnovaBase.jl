# AnovaBase.jl Documentation

*AnovaBase.jl* is a Julia package providing a simple framework for Analysis of Varaincae (ANOVA) on various types of julia statistical models.
It is similar to function [anova in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova).

## Functionality overview
### anova
```
anova(<model>; <type>, <test>)
anova(<test>, <model>; <type>)
anova(<models>; <test>)
anova(<test>, <models>)

anova_lm(<formula>, <data>; <type>, <test>)
anova_lm(<test>, <formula>, <data>; <type>)
anova_glm(<formula>, <data>, <distr>, <link>; <type>, <test>)
anova_glm(<test>, <formula>, <data>, <distr>, <link>; <type>)
anova_lme(<formula>, <data>; <type>, <test>)
anova_lme(<test>, <formula>, <data>; <type>)
anova_lfe(<formula>, <data>, <vcov>; <type>, <test>)
anova_lfe(<test>, <formula>, <data>, <vcov>; <type>)
```
### nestedmodels
```
nestedmodels(<modeltype>, <formula>, <data>)
nestedmodels(<model>)
```
## Usage
The usage is integrated with [`GLM.jl`](https://juliastats.org/GLM.jl/stable/), [`MixedModels.jl`](https://juliastats.org/MixedModels.jl/stable/) and [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl). 

This package is not intentded to be used directly.
Use the following packages for different models:
1. `AnovaGLM.jl` for `GLM.jl`.
2. `AnovaMixedModels.jl` for `MixedModels.JL`
3. `AnovaFixedEffectModels.jl` for `FixedEffectModels.jl`.

### Statistical Models
1. `TableRegressionModel{<: LinearModel, T}` fitted by `GLM.lm`.
2. `TableRegressionModel{<: GeneralizedLinearModel, T}` fitted by `GLM.glm`.
3. `LinearMixedModel` fitted by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`.
4. `GeneralizedLinearMixedModel` fitted by `MixedAnova.glme` or `fit(GeneralizedLinearMixedModel, ...)`
5. `TableRegressionModel{<: FixedEffectModel, T}` fitted by `MixedAnova.lfe`.

### Tests for Goodness of Fit
1. `FTest`: [F-test](https://en.wikipedia.org/wiki/F-test)
2. `LikelihoodRatioTest`, `LRT`: [likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test)

### Types of Estimable Functions
[Type I, II, III SS](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introglmest_sect001.htm)  are supported. 

## Table of Contents
```@contents
Pages = [
    "GLM.md",
    "MixedModels.md",
    "FixedEffectModels.md",
    "API.md"
]
Depth = 2
```