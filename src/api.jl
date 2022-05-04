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

<<<<<<< Updated upstream
Analysis of Variance by likelihood-ratio test.
=======
Analysis of Variance by likelihood-ratio test. \n
Return `AnovaResult{M, LRT, N}`. See `AnovaResult` for details.
"""
anova(::Type{LRT}, model) = error("anova by likelihood-ratio test is not defined for $(typeof(model)).")

# across different kind of models
function _lrt_nested(models::NTuple{N, RegressionModel}, df, dev, σ²; nestedwarn::Bool = true) where N
    nestedwarn && @warn "Could not check whether models are nested: results may not be meaningful"
    Δdf = _diff(df)
    Δdev = _diffn(dev)
    lrstat = Δdev ./ σ²
    pval = map(zip(Δdf, lrstat)) do (dof, lr)
        lr > 0 ? ccdf(Chisq(dof), lr) : NaN
    end
    AnovaResult{LRT}(models, 1, df, dev, (NaN, lrstat...), (NaN, pval...), NamedTuple())
end

"""
    isnullable(::CholeskyPivoted)
    isnullable(::Cholesky
    isnullable(::InverseLink)
    isnullable(::InverseSquareLink)
    isnullable(::Link)
    isnullable(::LinearModel)
    isnullable(model::GeneralizedLinearModel)
    isnullable(::LinearMixedModel)

Return `true` if empty model can be fitted.
"""
isnullable(m) = error("isnullable is not defined for $(typeof(m)).")

"""
    nobs(aov::AnovaResult{<: Tuple})
    nobs(aov::AnovaResult)

Apply `nobs` to all models in `aov.model`
"""
StatsBase.nobs(aov::AnovaResult{<: Tuple}) = nobs.(aov.model)
StatsBase.nobs(aov::AnovaResult) = nobs(aov.model)

"""
    dof(aov::AnovaResult)

Degree of freedom of models or factors.
"""
StatsBase.dof(aov::AnovaResult) = aov.dof

"""
    dof_residual(aov::AnovaResult{<: Tuple})
    dof_residual(aov::AnovaResult)
    dof_residual(aov::AnovaResult{<: MixedModel, FTest})

Degree of freedom of residuals.
Default is applying `dof_residual` to models in `aov.model`.
For `MixedModels` applying `FTest`, it is calculated by between-within method. See `calcdof` for details.
"""
StatsBase.dof_residual(aov::AnovaResult{<: Tuple}) = dof_residual.(aov.model)
StatsBase.dof_residual(aov::AnovaResult) = dof_residual(aov.model)
>>>>>>> Stashed changes

* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
<<<<<<< Updated upstream
anova(::Type{LRT}, model) = error("anova by likelihood-ratio test is not defined for $(typeof(model)).")
=======
    deviance(aov::AnovaResult)

Return the stored devaince. The value repressents different statistics for different models and tests.
1. `LinearModel`: Sum of Squares.
2. `GeneralizedLinearModel`: `deviance(model)`
3. `LinearMixedModel`: `NaN` when applying `FTest`; `-2loglikelihood(model) == deviance(model)` when applying `LRT`.
When `LinearModel` is compared to `LinearMixedModel`, the deviance is alternatively `-2loglikelihood(model)`.
"""
StatsBase.deviance(aov::AnovaResult) = aov.deviance

"""
    teststat(aov::AnovaResult)

Values of test statiscics of `anova`.
"""
teststat(aov::AnovaResult) = aov.teststat

"""
    teststat(aov::AnovaResult)

P-values of test statiscics of `anova`.
"""
pval(aov::AnovaResult) = aov.pval

"""
    anova_test(::AnovaResult)

Test statiscics of `anova`.
"""
anova_test(::AnovaResult{M, T}) where {M, T <: GoodnessOfFit} = T

"""
    anova_type(aov::AnovaResult)

Type of `anova`.
"""
anova_type(aov::AnovaResult) = aov.type

# Calculate dof from assign
function StatsBase.dof(v::Vector{Int})
    dofv = zeros(Int, v[end])
    prev = 1
    ind = 1
    n = length(v)
    while ind <= n
        v[ind] == prev || (prev = v[ind])
        dofv[prev] += 1
        ind += 1
    end
    Int.(dofv)
end
>>>>>>> Stashed changes
