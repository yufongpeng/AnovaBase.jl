# MixedAnova
|Documentation|CI status|Coverage|
|:-----------:|:-------:|:------:|
|[![Stable Docs][docs-stable-img]][docs-stable-url] [![Dev Docs][docs-dev-img]][docs-dev-url]| [![TravisCI][travis-img]][travis-url] [![][ci-img]][ci-url]| [![][codecov-img]][codecov-url]|

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://yufongpeng.github.io/MixedAnova.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://yufongpeng.github.io/MixedAnova.jl/stable
[travis-img]: https://travis-ci.com/yufongpeng/MixedAnova.jl.svg?branch=master
[travis-url]: https://travis-ci.com/github/yufongpeng/MixedAnova.jl
[ci-img]: https://github.com/yufongpeng/MixedAnova.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/yufongpeng/MixedAnova.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/yufongpeng/MixedAnova.jl/coveage.svg
[codecov-url]: https://codecov.io/gh/yufongpeng/MixedAnova.jl

Implement one-way and multi-way anova, including type 1, 2 and 3 sum of squares. The syntax and output resemble package `GLM`. 
The types of models supported:
1. `TableRegressionModel{<: LinearModel, T}` fit by `GLM.lm`.
2. `TableRegressionModel{<: GeneralizedLinearModel, T}` fit by `GLM.glm`.
3. `LinearMixedModel` fit by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`.
4. `TableRegressionModel{<: FixedEffectModel, T}` fit by `MixedAnova.lfe`.

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
## TO DO
1. Likelihood ratio test for `FixedEffectModels`.
2. Implementation of `Rao` and `Mallow's Cp`.




