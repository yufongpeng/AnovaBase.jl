# AnovaBase.jl Documentation

*AnovaBase.jl* is a Julia package providing a simple framework for Analysis of Varaincae (ANOVA) on various types of julia statistical models.
It is similar to function [anova in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova).

## Functionality overview
### ANOVA
```
anova(model; type, test)
anova(test, model; type)
anova(models; test)
anova(test, models)
```
### Model-specific ANOVA
#### AnovaGLM
```
anova_lm(formula, data; type, test)
anova_lm(test, formula, data; type)
anova_glm(formula, data, distr, link; type, test)
anova_glm(test, formula, data, distr, link; type)
```
#### AnovaMixedModels
```
anova_lme(formula, data; type, test)
anova_lme(test, formula, data; type)
```
#### AnovaFixedEffecModels
```
anova_lfe(formula, data, vcov; type, test)
anova_lfe(test, formula, data, vcov; type)
```
### nestedmodels
```
nestedmodels(modeltype, formula, data)
nestedmodels(model)
```
## Usage
This package is not intentded to be used directly.
Use the following packages for different models:

|Packages for models|Packages for ANOVA|Models|Fitted by|
|-------------------|------------------|------|---------|
|[GLM.jl](https://juliastats.org/GLM.jl/stable/)|[AnovaGLM.jl](https://github.com/yufongpeng/AnovaGLM.jl)|`TableRegressionModel{<: LinearModel}`|`GLM.lm` or `fit(LinearModel, _)`|
|||`TableRegressionModel{<: GeneralizedLinearModel}`|`GLM.glm` or `fit(GeneralizedLinearModel, _)`|
|[MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/)|[AnovaMixedModels.jl](https://github.com/yufongpeng/AnovaMixedModels.jl)|`LinearMixedModel`|`AnovaMixedModels.lme` or `fit(LinearMixedModel, _)`|
|||`GeneralizedLinearMixedModel`|`AnovaGLM.glme` or `fit(GeneralizedLinearMixedModel, _)`|
|[FixedEffectModels.jl](https://github.com/FixedEffects/FixedEffectModels.jl)|[AnovaFixedEffectModels.jl](https://github.com/yufongpeng/AnovaFixedEffectModels.jl)|`FixedEffectModel`|`AnovaFixedEffectModels.lfe` or `FixedEffectModels.reg`|

### Tests for Goodness of Fit
1. `FTest`: [F-test](https://en.wikipedia.org/wiki/F-test)
2. `LikelihoodRatioTest`, `LRT`: [likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test)

### Types of Estimable Functions
[Type I, II, III SS](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_introglmest_sect001.htm)  are supported. 

## Table of Contents
### Examples
```@contents
Pages = [
    "Examples_GLM.md",
    "Examples_MixedModels.md",
    "Examples_FixedEffectModels.md"
]
Depth = 2
```
### Algorithm
```@contents
Pages = [
    "Algorithm_AnovaGLM.md",
    "Algorithm_AnovaMixedModels.md",
    "Algorithm_AnovaFixedEffectModels.md"
]
Depth = 2
```
### Interfacing AnovaBase.jl
```@contents
Pages = ["Interface.md"]
Depth = 2
```
### API
```@contents
Pages = [
    "AnovaBase.md",
    "AnovaGLM.md",
    "AnovaMixedModels.md",
    "AnovaFixedEffectModels.md"
]
```