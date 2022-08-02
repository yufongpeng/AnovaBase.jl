# General api
"""
    formula(model)
    formula(trm::TableRegressionModel)

Unified api for formula of each statistical model.
"""
formula(trm::TableRegressionModel) = trm.mf.f

"""
    coefnames(aov::AnovaResult)
    coefnames(<model>, anova::Val{:anova})

Customize coefnames for each models.
"""
coefnames(aov::AnovaResult) = coefnames(aov.model, Val(:anova))
coefnames(model::RegressionModel, anova::Val{:anova}) = coefnames(formula(model).rhs, anova)

"""
    nestedmodels(<model>; <keyword arguments>)
    nestedmodels(<model type>, formula, data, <keyword arguments>)

Generate nested models from a model or formula and data.
"""
function nestedmodels end

# implement drop1/add1 in R?
"""
    anova(<models>...; test::Type{<: GoodnessOfFit})
    anova(test::Type{FTest}, <model>; <keyword arguments>)
    anova(test::Type{FTest}, <models>...; <keyword arguments>)
    anova(test::Type{LRT}, <model>; <keyword arguments>)
    anova(test::Type{LRT}, <models>...; <keyword arguments>)

Analysis of variance.

Return `AnovaResult{M, test, N}`. See [`AnovaResult`](@ref) for details.

* `models`: model objects. If mutiple models are provided, they should be nested and the last one is the most saturated.
* `test`: test statistics for goodness of fit. Available tests are [`LikelihoodRatioTest`](@ref) ([`LRT`](@ref)) and [`FTest`](@ref).
"""
function anova end

"""
    nobs(aov::AnovaResult{<: Tuple})
    nobs(aov::AnovaResult)

Apply `nobs` to all models in `aov.model`. See [`AnovaResult`](@ref) for details.
"""
nobs(aov::AnovaResult{<: Tuple}) = nobs.(aov.model)
nobs(aov::AnovaResult) = nobs(aov.model)

"""
    dof(aov::AnovaResult)

Degree of freedom of models or factors.
"""
dof(aov::AnovaResult) = aov.dof

"""
    dof_residual(aov::AnovaResult{<: Tuple})
    dof_residual(aov::AnovaResult)
    dof_residual(aov::AnovaResult{<: MixedModel, FTest})

Degree of freedom of residuals.

Default is applying `dof_residual` to models in `aov.model`.
"""
dof_residual(aov::AnovaResult{<: Tuple}) = dof_residual.(aov.model)
dof_residual(aov::AnovaResult) = dof_residual(aov.model)

"""
    deviance(aov::AnovaResult)

Return the stored devaince. The value repressents different statistics for different models and tests.
1. Ordinary linear regression: `sum of squares`
2. Generalized linear regression: `deviance`
3. Mixed-effect model: `NaN` when applying [`FTest`](@ref); `-2×loglikelihood` or `deviance` when applying [`LRT`](@ref)
When ordinary linear regression is compared to linear mixed-effect model, the deviance is alternatively `-2×loglikelihood`
"""
deviance(aov::AnovaResult) = aov.deviance

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
