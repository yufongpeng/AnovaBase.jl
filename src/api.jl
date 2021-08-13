# General api
"""
    formula(model::TableRegressionModel)
    formula(model::MixedModel)

Unify formula api.
"""
formula(model) = error("formula is not defined for $(typeof(model)).")

"""
    nestedmodels(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; null::Bool = false, <keyword arguments>)
    nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = false, <keyword arguments>)
    nestedmodels(model::LinearMixedModel; null::Bool = false, <keyword arguments>)

    nestedmodels(::Type{LinearModel}, formula, data; null::Bool = true, <keyword arguments>)
    nestedmodels(::Type{GeneralizedLinearModel}, formula, data, distr::UnivariateDistribution, link::Link = canonicallink(d); null::Bool = true, <keyword arguments>)
    nestedmodels(::Type{LinearMixedModel}, f::FormulaTerm, tbl; null::Bool = true, wts = [], contrasts = Dict{Symbol, Any}(), verbose::Bool = false, REML::Bool = false)

Generate nested models from a saturated model or formula and data. \n
The null model will be a model with at least one factor (including intercept) if the link function does not allow factors to be 0 (factors in denominators). \n
* `InverseLink` for `Gamma`
* `InverseSquareLink` for `InverseGaussian`
Otherwise, it will be a model with no factors.
"""
nestedmodels(model) = error("nestedmodels is not defined for $(typeof(model)).")

const FixDispDist = Union{Bernoulli, Binomial, Poisson}
"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return LRT if the distribution has fixed dispersion
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest

"""
    anova(<models>...; test::Type{T}) where {T <: GoodnessOfFit}

Analysis of variance.

* `models`: model objects
    1. `TableRegressionModel{<: LinearModel, <: AbstractArray}` fit by `GLM.lm`
    2. `TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}` fit by `GLM.glm`
    3. `LinearMixedModel` fit by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`
    If mutiple models are provided, they should be nested and the last one is the most saturated.
* `test`: test statistics for goodness of fit. Available tests are `LikelihoodRatioTest` (`LRT`) and `FTest`. \n
    If no test argument is provided, the function will automatically determine based on the model type:
    1. `TableRegressionModel{<: LinearModel, <: AbstractArray}`: `FTest`.
    2. `TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}`: based on distribution function, see `canonicalgoodnessoffit`.
    3. `LinearMixedModel`: `FTest` for one model, `LRT` for nested models.

For fitting new models and conducting anova at the same time,  
see `anova_lm` for `LinearModel`, `anova_lme` for `LinearMixedModel`, `anova_glm` for `GeneralizedLinearModel`.
"""
anova(model) = error("anova is not defined for $(typeof(model)).")

"""
    anova(::Type{FTest}, <model>; kwargs...)
    anova(::Type{FTest}, <models>...; kwargs...)

Analysis of Variance by F-test.

* `type` specifies type of anova. For one `LinearModel` `1, 2, 3` are valid; for one `LinearMixedModel` `1, 3` are valid. For others, only `1` is valid.
* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
anova(::Type{FTest}, model) = error("anova by F-test is not defined for $(typeof(model)).")

"""
    anova(::Type{LRT}, <model>; kwargs...)
    anova(::Type{LRT}, <models>...; kwargs...)

Analysis of Variance by likelihood-ratio test.

* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
anova(::Type{LRT}, model) = error("anova by likelihood-ratio test is not defined for $(typeof(model)).")
