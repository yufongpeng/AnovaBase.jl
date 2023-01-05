# AnovaResult api
const doc_nestedmodels = """
    nestedmodels(<model>; <keyword arguments>)
    nestedmodels(<model type>, formula, data; <keyword arguments>)

Generate nested models from a model or formula and data.
"""
@doc doc_nestedmodels doc_nestedmodels

@doc doc_nestedmodels
function nestedmodels(::T; kwargs...) where {T <: RegressionModel} 
    throw(function_arg_error(nestedmodels, T))
end

@doc doc_nestedmodels
function nestedmodels(::Type{T}, f::FormulaTerm, tbl::S; kwargs...) where {T <: RegressionModel, S}
    throw(function_arg_error(nestedmodels, "::Type{$T}), ::FormulaTerm, ::$S"))
end

# implement drop1/add1 in R?
const doc_anova = """
    anova(<models>...; test::Type{<: GoodnessOfFit}, <keyword arguments>)
    anova(Test::Type{FTest}, <model>; <keyword arguments>)
    anova(Test::Type{FTest}, <models>...; <keyword arguments>)
    anova(Test::Type{LRT}, <model>; <keyword arguments>)
    anova(Test::Type{LRT}, <models>...; <keyword arguments>)

Analysis of variance.

Return `AnovaResult{M, Test, N}`. See [`AnovaResult`](@ref) for details.

* `models`: model objects. If mutiple models are provided, they should be nested, fitted with the same data and the last one is the most complex.
* `Test`: test statistics for goodness of fit. Available tests are [`LikelihoodRatioTest`](@ref) ([`LRT`](@ref)) and [`FTest`](@ref).
"""
@doc doc_anova doc_anova
@doc doc_anova
function anova(models::Vararg{T}; test::Type{S}, kwargs...) where {T <: RegressionModel, S <: GoodnessOfFit}
    throw(function_arg_error(anova, "::Vararg{$T}; test::Type{$S})"))
end
@doc doc_anova
function anova(Test::Type{FTest}, models::T; kwargs...) where {T <: RegressionModel}
    throw(function_arg_error(anova, "::Type{FTest}, ::$T"))
end
@doc doc_anova
function anova(Test::Type{FTest}, models::Vararg{T}; kwargs...) where {T <: RegressionModel}
    throw(function_arg_error(anova, "::Type{FTest}, ::Vararg{$T}"))
end
@doc doc_anova
function anova(Test::Type{LRT}, models::T; kwargs...) where {T <: RegressionModel}
    throw(function_arg_error(anova, "::Type{LRT}, ::$T"))
end
@doc doc_anova
function anova(Test::Type{LRT}, models::Vararg{T}; kwargs...) where {T <: RegressionModel}
    throw(function_arg_error(anova, "::Type{LRT}, ::Vararg{$T}"))
end

"""
    nobs(aov::AnovaResult)
    nobs(aov::AnovaResult{<: Tuple})

Number of observations.
"""
nobs(aov::AnovaResult) = nobs(aov.model)
nobs(aov::AnovaResult{<: Tuple}) = nobs(first(aov.model))

"""
    dof(aov::AnovaResult)

Degrees of freedom of each models or predictors.
"""
dof(aov::AnovaResult) = aov.dof

"""
    dof_residual(aov::AnovaResult{<: Tuple})
    dof_residual(aov::AnovaResult)

Degrees of freedom of residuals.

By default, it applies `dof_residual` to models in `aov.model`.
"""
dof_residual(aov::AnovaResult{<: Tuple}) = dof_residual.(aov.model)
dof_residual(aov::AnovaResult{M, T, N}) where {M, T, N} = ntuple(x -> dof_residual(aov.model), N)

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

Type of [`anova`](@ref). See [`AnovaResult`](@ref) for details.
"""
anova_type(aov::AnovaResult) = aov.type
