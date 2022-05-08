# General api
"""
    formula(trm::TableRegressionModel)
    formula(model::MixedModel)

Unified api for formula of each statistical model.
"""
formula(model) = error("formula is not defined for $(typeof(model)).")

"""
    nestedmodels(trm::TableRegressionModel{<: LinearModel}; null::Bool = false, <keyword arguments>)
    nestedmodels(trm::TableRegressionModel{<: GeneralizedLinearModel}; null::Bool = false, <keyword arguments>)
    nestedmodels(model::LinearMixedModel; null::Bool = false, <keyword arguments>)

    nestedmodels(::Type{LinearModel}, formula, data; null::Bool = true, <keyword arguments>)
    nestedmodels(::Type{GeneralizedLinearModel}, formula, data, distr::UnivariateDistribution, link::Link = canonicallink(d); null::Bool = true, <keyword arguments>)
    nestedmodels(::Type{LinearMixedModel}, f::FormulaTerm, tbl; null::Bool = true, wts = [], contrasts = Dict{Symbol, Any}(), verbose::Bool = false, REML::Bool = false)

Generate nested models from a model or formula and data.

The null model will be a model with at least one factor (including intercept) if the link function does not allow factors to be 0 (factors in denominators).
* `InverseLink` for `Gamma`
* `InverseSquareLink` for `InverseGaussian`
* `LinearModel` fitted with `CholeskyPivoted` when `dropcollinear = true`
Otherwise, it will be an empty model.
"""
nestedmodels(model) = error("nestedmodels is not defined for $(typeof(model)).")

# implement drop1/add1 in R?

const FixDispDist = Union{Bernoulli, Binomial, Poisson}
"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return LRT if the distribution has a fixed dispersion.
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest

"""
    anova(<models>...; test::Type{<: GoodnessOfFit})

Analysis of variance.

Return `AnovaResult{M, test, N}`. See [`AnovaResult`](@ref) for details.

* `models`: model objects
    1. `TableRegressionModel{<: LinearModel}` fitted by `GLM.lm`
    2. `TableRegressionModel{<: GeneralizedLinearModel}` fitted by `GLM.glm`
    3. `LinearMixedModel` fitted by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`
    4. `GeneralizedLinearMixedModel` fitted by `MixedAnova.glme` or `fit(GeneralizedLinearMixedModel, ...)`
    5. `TableRegressionModel{<: FixedEffectModel}` fitted by `MixedAnova.lfe`.
    If mutiple models are provided, they should be nested and the last one is the most saturated.
* `test`: test statistics for goodness of fit. Available tests are [`LikelihoodRatioTest`](@ref) ([`LRT`](@ref)) and [`FTest`](@ref). The default is based on the model type.
    1. `TableRegressionModel{<: LinearModel}`: `FTest`.
    2. `TableRegressionModel{<: GeneralizedLinearModel}`: based on distribution function, see `canonicalgoodnessoffit`.
    3. `LinearMixedModel`: `FTest` for one model fit; `LRT` for nested models.
    4. `GeneralizedLinearMixedModel`: `LRT` for nested models.
    5. `TableRegressionModel{<: FixedEffectModel}`: `FTest`.

When multiple models are provided:  
* `check`: allows to check if models are nested. Defalut value is true. Some checkers are not implemented now.
* `isnested`: true when models are checked as nested (manually or automatically). Defalut value is false. 

For fitting new models and conducting anova at the same time, see [`anova_lm`](@ref) for `LinearModel`, [`anova_glm`](@ref) for `GeneralizedLinearModel`, [`anova_lme`](@ref) for `LinearMixedModel`, and [`anova_lfe`](@ref) for `FixedEffectModel`.
"""
anova(model) = error("anova is not defined for $(typeof(model)).")

"""
    anova(::Type{FTest}, <model>; kwargs...)
    anova(::Type{FTest}, <models>...; kwargs...)

Analysis of Variance by F-test.

Return `AnovaResult{M, FTest, N}`. See [`AnovaResult`](@ref) for details.

* `type` specifies type of anova: 
    1. One `LinearModel` or `GeneralizedLinearModel`: 1, 2, 3 are valid
    2. One `LinearMixedModel`: 1, 3 are valid. 
    3. Others: only 1 is valid.  
* `adjust_sigma` determines if adjusting to REML when `LinearMixedModel` is fitted by maximum likelihood.
!!! note
    The result with adjustments will be slightly deviated from that of model fitted directly by REML.
"""
anova(::Type{FTest}, model) = error("anova by F-test is not defined for $(typeof(model)).")

"""
    anova(::Type{LRT}, <model>; kwargs...)
    anova(::Type{LRT}, <models>...; kwargs...)

Analysis of Variance by likelihood-ratio test.

Return `AnovaResult{M, LRT, N}`. See [`AnovaResult`](@ref) for details.
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

Apply `nobs` to all models in `aov.model`. See [`AnovaResult`](@ref) for details.
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
!!! note 
    For `MixedModels` applying [`FTest`](@ref), it is calculated by between-within method. See [`calcdof`](@ref) for details.
"""
StatsBase.dof_residual(aov::AnovaResult{<: Tuple}) = dof_residual.(aov.model)
StatsBase.dof_residual(aov::AnovaResult) = dof_residual(aov.model)

"""
    deviance(aov::AnovaResult)

Return the stored devaince. The value repressents different statistics for different models and tests.
1. `LinearModel`: Sum of Squares.
2. `GeneralizedLinearModel`: `deviance(model)`
3. `MixedModel`: `NaN` when applying [`FTest`](@ref); `-2loglikelihood(model) == deviance(model)` when applying [`LRT`](@ref).
When `LinearModel` is compared to `LinearMixedModel`, the deviance is alternatively `-2loglikelihood(model)`.
"""
StatsBase.deviance(aov::AnovaResult) = aov.deviance

"""
    teststat(aov::AnovaResult)

Values of test statiscics of [`anova`](@ref). See [`AnovaResult`](@ref) for details.
"""
teststat(aov::AnovaResult) = aov.teststat

"""
    teststat(aov::AnovaResult)

P-values of test statiscics of [`anova`](@ref). See [`AnovaResult`](@ref) for details.
"""
pval(aov::AnovaResult) = aov.pval

"""
    anova_test(::AnovaResult)

Test statiscics of [`anova`](@ref). See [`AnovaResult`](@ref) for details.
"""
anova_test(::AnovaResult{M, T}) where {M, T <: GoodnessOfFit} = T

"""
    anova_type(aov::AnovaResult)

Type of [`anova`](@ref). See [`AnovaResult`](@ref) for details.
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
