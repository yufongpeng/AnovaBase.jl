# AnovaBase
|Documentation|CI status|Coverage|
|:-----------:|:-------:|:------:|
|[![Stable Docs][docs-stable-img]][docs-stable-url] [![Dev Docs][docs-dev-img]][docs-dev-url]| [![][ci-img]][ci-url]| [![][codecov-img]][codecov-url]|

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://yufongpeng.github.io/AnovaBase.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://yufongpeng.github.io/AnovaBase.jl/stable
[ci-img]: https://github.com/yufongpeng/AnovaBase.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/yufongpeng/AnovaBase.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/yufongpeng/AnovaBase.jl/branch/main/graph/badge.svg?token=XXS9TJB8IU
[codecov-url]: https://codecov.io/gh/yufongpeng/AnovaBase.jl

*AnovaBase.jl* is a Julia package providing a simple framework for Analysis of Variance (ANOVA) on various types of Julia statistical models.
It is similar to function [anova in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova).

## Usage
This package is not intentded to be used directly.
Use the following packages for different models:
|Packages for models|Packages for ANOVA|Models|Fitted by|
|:-----------------:|:----------------:|:----:|:-------:|
|[GLM.jl](https://juliastats.org/GLM.jl/stable/)|[AnovaGLM.jl](https://github.com/yufongpeng/AnovaGLM.jl)|`TableRegressionModel{<: LinearModel}`|`GLM.lm` or `fit(LinearModel, _)`|
| | |`TableRegressionModel{<: GeneralizedLinearModel}`|`GLM.glm` or `fit(GeneralizedLinearModel, _)`|
|[MixedModels.jl](https://juliastats.org/MixedModels.jl/stable/)|[AnovaMixedModels.jl](https://github.com/yufongpeng/AnovaMixedModels.jl)|`LinearMixedModel`|`MixedModels.lmm` or `fit(LinearMixedModel, _)`|
| | |`GeneralizedLinearMixedModel`|`MixedModels.glmm` or `fit(GeneralizedLinearMixedModel, _)`|
|[FixedEffectModels.jl](https://github.com/FixedEffects/FixedEffectModels.jl)|[AnovaFixedEffectModels.jl](https://github.com/yufongpeng/AnovaFixedEffectModels.jl) |`FixedEffectModel`|`AnovaFixedEffectModels.lfe` or `FixedEffectModels.reg`|

## TO DO
1. Implementation of Rao and Mallow's Cp.