# MixedAnova.jl Documentation

*MixedAnova.jl* is a Julia package providing functionality of Analysis of Varaincae (ANOVA) for various types of julia statistical models.
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
The usage is integrated with [`GLM.jl`](https://juliastats.org/GLM.jl/stable/), [`MixedModels.jl`](https://juliastats.org/MixedModels.jl/stable/) and [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl). The examples below are just for dementration of ANOVA. To understand details of models, please refer to each documentations.

```@contents
Pages = [
    "GLM.md",
    "MixedModels.md",
    "FixedEffectModels.md",
    "API.md"
]
Depth = 2
```