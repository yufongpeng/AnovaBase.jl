# Interface
Assume there is a package defining a model type `Some_Model`.

## [anova](./AnovaBase.md#AnovaBase.anova)
`anova` can be overloaded in the following arguments signature.
1. `anova(::Some_Model; test, type, kwargs...)` 
2. `anova(::GoodnessOfFit, ::Some_Model; type, kwargs...)`
3. `anova(::Vararg{Some_Model}; test, kwargs...)`
4. `anova(::GoodnessOfFit, ::Vararg{Some_Model}; kwargs...)`

If `nestedmodels` is defined, the method for `NestedModels` should be defined in order to use `anova` on the output directly.

Models should be wrapped in the following types

|Model|Wrapper|
|:---:|:-----:|
|Single model|`FullModel`|
|Nested models|`NestedModels`|
|Nested models with different types|`MixedAovModels`|

The main algorithm is recommended to be implemented in the method dispatched on wrapper type and should eventually return `AnovaResult`.

## [predictors](./AnovaBase.md#AnovaBase.predictors)
This function returns a tuple of terms which can be used in ANOVA (some terms may not be used because of ANOVA type or model itself).
There is a default method for `RegressionModel`, i.e., `formula(model).rhs.terms`. If the formula for `Some_Model` has special structure like `MixedModel`, this function should be overloaded. 

## AnovaModels
### [NestedModels](./AnovaBase.md#AnovaBase.NestedModels)
Wrap nested models in the following way
1. `NestedModels{Some_Model}((model1, model2, ...))`
2. `NestedModels{Some_Model}(model1, model2, ...)`

### [MixedAovModels](./AnovaBase.md#AnovaBase.MixedAovModels)
This is for comparing with models with different types.
1. `MixedAovModels{Union{Some_Model, Other_Model}}((model1, model2, ...))`
2. `MixedAovModels{Union{Some_Model, Other_Model}}(model1, model2, ...)`

Both `NestedModels` and `MixedAovModels` do not check if the models' order is correct, so user should be careful to put simpler models former and more complex models later.
ANOVA can be computed along with some basic statistics. See [`ftest_nested`](./AnovaBase.md#AnovaBase.ftest_nested) and [`lrt_nested`](./AnovaBase.md#AnovaBase.lrt_nested).

### [FullModel](./AnovaBase.md#AnovaBase.FullModel)
`FullModel` wraps the model along with the index of predictors(`pred_id`) that is actually used in ANOVA and ANOVA type.
`AnovaBase` provides a method which automatically determines `pred_id` based on ANOVA type, whether empty model is allowed, and whether intercept is tested. See [`FullModel`](./AnovaBase.md#AnovaBase.FullModel).

To customize `FullModel`, there are many helper functions for manipulating terms. See [Developer utility](./AnovaBase.md#developer-utility)

## [anovatable](./AnovaBase.md#AnovaBase.anovatable)
This function returns a `AnovaTable` for showing `AnovaResult`. For `NestedModels` and `MixedAovModels`, `AnovaBase` provides a default interface. When dealing with `FullModel` or customizing methods on `Some_Model`, the following methods should be defined.
1. `anovatable(::AnovaResult{<: FullModel{Some_Model}}; rownames)`
2. `anovatable(::AnovaResult{<: NestedModels{Some_Model}}; rownames)`
3. `anovatable(::AnovaResult{<: MixedAovModels{Some_Model}}; rownames)`

See [`AnovaTable`](./AnovaBase.md#AnovaBase.AnovaTable) for the argument signature.

## [nestedmodels](./AnovaBase.md#AnovaBase.nestedmodels)
This function is not essential for ANOVA; it is just for convenience to create nested models(`NestedModels`). It should be defined as follow
1. `nestedmodels(Some_Model, ::FormulaTerm, data; kwargs)`
2. `nestedmodels(::Some_Model; kwargs)`

`AnovaBase` provides a lot of functions to work on formula, terms and contrasts. See [Developer utility](./AnovaBase.md#developer-utility)

## Other function
* [`dof_residual`](./AnovaBase.md#AnovaBase.dof_residual) applies `dof_residual` to all models by default. If `dof_residual(::Some_Model)` is not valid for ANOVA, customize `dof_residual(::AnovaResult{<: AnovaModel{Some_Model}})` alternatively.