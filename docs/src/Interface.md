# Interfacing AnovaBase.jl
Given model types `SomeModel` and `OtherModel`, the following functions have to be defined or used.

## [anova](./AnovaBase.md#AnovaBase.anova-Tuple{Type{<:GoodnessOfFit},%20AnovaModel})
`anova` can be overloaded in the following arguments signature.
1. `anova(::SomeModel; test, type, kwargs...)` 
2. `anova(::GoodnessOfFit, ::SomeModel; type, kwargs...)`
3. `anova(::Vararg{SomeModel}; test, kwargs...)`
4. `anova(::GoodnessOfFit, ::Vararg{SomeModel}; kwargs...)`

Other arguments signatures are acceptable, but there will be `Type piracy` warnings.

Models should be wrapped in the following types before storing in the returned `AnovaResult`.

|Model|Wrapper|
|:---:|:-----:|
|Single model|`FullModel`|
|Nested models|`NestedModels`|
|Nested models with different types|`MixedAovModels`|

It is recommended that `anova` additionally dispatched on wrapper types and the main algorithm is implemented in it.

## [predictors](./AnovaBase.md#AnovaBase.predictors-Tuple{RegressionModel})
This function returns a tuple of terms which can be used in ANOVA (some terms may not be used because of ANOVA type or model itself).

There is a default method for `RegressionModel`, i.e., `formula(model).rhs.terms`. If the formula for `SomeModel` has special structure like `MixedModel`, this function should be overloaded. 

## AnovaModels
### [NestedModels](./AnovaBase.md#AnovaBase.NestedModels)
Wrap nested models in the following way
1. `NestedModels{SomeModel}((model1, model2, ...))`
2. `NestedModels{SomeModel}(model1, model2, ...)`

### [MixedAovModels](./AnovaBase.md#AnovaBase.MixedAovModels)
This is for comparing models with different types.
1. `MixedAovModels{Union{SomeModel, OtherModel}}((model1, model2, ...))`
2. `MixedAovModels{Union{SomeModel, OtherModel}}(model1, model2, ...)`

Both `NestedModels` and `MixedAovModels` do not check if the order is correct, so user should be careful to put simpler models former and more complex models later.
ANOVA can be computed along with some basic statistics. See [`ftest_nested`](./AnovaBase.md#AnovaBase.ftest_nested) and [`lrt_nested`](./AnovaBase.md#AnovaBase.lrt_nested).

### [FullModel](./AnovaBase.md#AnovaBase.FullModel)
`FullModel` wraps the model along with the index of predictors(`pred_id`) that is actually used in ANOVA and ANOVA type.

`AnovaBase` provides a method which automatically determines `pred_id` based on ANOVA type, whether empty model is allowed, and whether intercept is tested. See [`FullModel`](./AnovaBase.md#AnovaBase.FullModel).

To customize `FullModel`, there are many helper functions for manipulating terms. See [Developer utility](./AnovaBase.md#Developer-utility)

## [anovatable](./AnovaBase.md#AnovaBase.anovatable-Tuple{AnovaResult{<:FullModel}})
This function returns a `AnovaTable` for showing `AnovaResult`. For `NestedModels` and `MixedAovModels`, `AnovaBase` provides a default interface. When dealing with `FullModel` or customizing methods on `SomeModel`, the following methods should be defined.
1. `anovatable(::AnovaResult{<: FullModel{SomeModel}}; rownames)`
2. `anovatable(::AnovaResult{<: NestedModels{SomeModel}}; rownames)`
3. `anovatable(::AnovaResult{<: MixedAovModels{SomeModel}}; rownames)`

See [`AnovaTable`](./AnovaBase.md#AnovaBase.AnovaTable) for the argument signature.

## [nestedmodels](./AnovaBase.md#AnovaBase.nestedmodels-Tuple{RegressionModel})
This function is not essential for ANOVA; it is just for convenience to create nested models(`NestedModels`). It should be defined as follow
1. `nestedmodels(SomeModel, ::FormulaTerm, data; kwargs)`
2. `nestedmodels(::SomeModel; kwargs)`

`AnovaBase` provides a lot of functions to work on formula, terms and contrasts. See [Developer utility](./AnovaBase.md#Developer-utility)

## Other function
* [`dof_residual`](./AnovaBase.md#StatsAPI.dof_residual-Tuple{AnovaResult}) applies `dof_residual` to all models by default. If `dof_residual(::SomeModel)` is not valid for ANOVA, customize `dof_residual(::AnovaResult{<: AnovaModel{SomeModel}})` alternatively.