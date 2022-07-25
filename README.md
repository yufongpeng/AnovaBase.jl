# AnovaBase
|Documentation|CI status|Coverage|
|:-----------:|:-------:|:------:|
|[![Stable Docs][docs-stable-img]][docs-stable-url] [![Dev Docs][docs-dev-img]][docs-dev-url]| [![TravisCI][travis-img]][travis-url] [![][ci-img]][ci-url]| [![][codecov-img]][codecov-url]|

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://yufongpeng.github.io/AnovaBase.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://yufongpeng.github.io/AnovaBase.jl/stable
[travis-img]: https://travis-ci.com/yufongpeng/AnovaBase.jl.svg?branch=master
[travis-url]: https://travis-ci.com/github/yufongpeng/AnovaBase.jl
[ci-img]: https://github.com/yufongpeng/AnovaBase.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/yufongpeng/AnovaBase.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/yufongpeng/AnovaBase.jl/coveage.svg
[codecov-url]: https://codecov.io/gh/yufongpeng/AnovaBase.jl

*AnovaBase.jl* is a Julia package providing a simple framework for Analysis of Variance (ANOVA) on various types of Julia statistical models.
It is similar to function [anova in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova).

## Usage
The usage is integrated with [`GLM.jl`](https://juliastats.org/GLM.jl/stable/), [`MixedModels.jl`](https://juliastats.org/MixedModels.jl/stable/) and [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl). 

This package is not intentded to be used directly.
Use the following packages for different models:
1. `AnovaGLM.jl` for `TableRegressionModel{<: LinearModel, T}` fitted by `GLM.lm` and `TableRegressionModel{<: GeneralizedLinearModel, T}` fitted by `GLM.glm`.
2. `AnovaMixedModels.jl` for `LinearMixedModel` fitted by `AnovaMixedModels.lme` or `fit(LinearMixedModel, ...)` and `GeneralizedLinearMixedModel` fitted by `AnovaGLM.glme` or `fit(GeneralizedLinearMixedModel, ...)`
3. `AnovaFixedEffectModels.jl` for `TableRegressionModel{<: FixedEffectModel, T}` fitted by `AnovaFixedEffectModels.lfe`.

## TO DO
1. Likelihood ratio test for `FixedEffectModels`.
2. `nestedmodels` for `FixedEffectModels`.
3. Implementation of Rao and Mallow's Cp.
