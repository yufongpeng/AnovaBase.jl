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

*MixedAnova.jl* is a Julia package providing functionality of Analysis of Varaincae (ANOVA) for various types of julia statistical models.
It is similar to function [anova in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova).

## Usage
The usage is integrated with [`GLM.jl`](https://juliastats.org/GLM.jl/stable/), [`MixedModels.jl`](https://juliastats.org/MixedModels.jl/stable/) and [`FixedEffectModels.jl`](https://github.com/FixedEffects/FixedEffectModels.jl). 

Supported models:
1. `TableRegressionModel{<: LinearModel, T}` fitted by `GLM.lm`.
2. `TableRegressionModel{<: GeneralizedLinearModel, T}` fitted by `GLM.glm`.
3. `LinearMixedModel` fitted by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`.
4. `GeneralizedLinearMixedModel` fitted by `MixedAnova.glme` or `fit(GeneralizedLinearMixedModel, ...)`
5. `TableRegressionModel{<: FixedEffectModel, T}` fitted by `MixedAnova.lfe`.

## TO DO
1. Likelihood ratio test for `FixedEffectModels`.
2. Implementation of `Rao` and `Mallow's Cp`.




