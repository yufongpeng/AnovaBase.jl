# AnovaResult api
"""
    nobs(aov::AnovaResult)
    nobs(aov::AnovaResult{<: MultiAovModels})

Number of observations.
"""
nobs(aov::AnovaResult) = round(Int, nobs(aov.anovamodel.model))
nobs(aov::AnovaResult{<: MultiAovModels}) = round(Int, nobs(first(aov.anovamodel.model)))

"""
    dof(aov::AnovaResult)

Degrees of freedom of each models or predictors.
"""
dof(aov::AnovaResult) = aov.dof

"""
    deviance(aov::AnovaResult)

Return the stored devaince. The value repressents different statistics for different models and tests. 
It may be deviance, Δdeviance, -2loglikelihood or other measures of model performance.
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
    anova_type(model::MultiAovModels)
    anova_type(model::FullModel)

Type of [`anova`](@ref), either 1, 2 or 3.
"""
anova_type(aov::AnovaResult) = anova_type(aov.anovamodel)
anova_type(model::MultiAovModels) = 1
anova_type(model::FullModel) = model.type

