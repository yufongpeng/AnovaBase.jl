# Interfaces need defined in user-end packages
"""
    nestedmodels(<model>; <keyword arguments>)
    nestedmodels(<model type>, formula, data; <keyword arguments>)

Generate nested models `NestedModels` from a model or modeltype, formula and data.
"""
function nestedmodels(::T; kwargs...) where {T <: RegressionModel} 
    throw(function_arg_error(nestedmodels, T))
end
function nestedmodels(::Type{T}, f::FormulaTerm, tbl::S; kwargs...) where {T <: RegressionModel, S}
    throw(function_arg_error(nestedmodels, "::Type{$T}), ::FormulaTerm, ::$S"))
end

# implement drop1/add1 in R?
"""
    anova(<models>...; test::Type{<: GoodnessOfFit}, <keyword arguments>)
    anova(Test::Type{<: GoodnessOfFit}, <model>; <keyword arguments>)
    anova(Test::Type{<: GoodnessOfFit}, <models>...; <keyword arguments>)

Analysis of variance.

Return `AnovaResult{M, Test, N}`. See [`AnovaResult`](@ref) for details.

* `models`: model objects. If mutiple models are provided, they should be nested, fitted with the same data and the last one is the most complex.
* `Test`: test statistics for goodness of fit. Available tests are [`LikelihoodRatioTest`](@ref) (`LRT`) and [`FTest`](@ref).
"""
function anova(Test::Type{T}, model::S; kwargs...) where {T <: GoodnessOfFit, S <: RegressionModel}
    throw(function_arg_error(anova, "::Type{$T}, ::$S"))
end
function anova(Test::Type{T}, model::Vararg{S}; kwargs...) where {T <: GoodnessOfFit, S <: RegressionModel}
    throw(function_arg_error(anova, "::Type{$T}, ::Vararg{$S}"))
end
function anova(models::Vararg{T}; test::Type{S}, kwargs...) where {T <: RegressionModel, S <: GoodnessOfFit}
    throw(function_arg_error(anova, "::Vararg{$T}; test::Type{$S})"))
end

"""
    dof_residual(aov::AnovaResult)    
    dof_residual(aov::AnovaResult{<: NestedModels})

Degrees of freedom of residuals.

By default, it applies `dof_residual` to models in `aov.anovamodel`.
"""
dof_residual(aov::AnovaResult{M, T, N}) where {M, T, N} = ntuple(x -> dof_residual(aov.anovamodel.model), N)
dof_residual(aov::AnovaResult{<: NestedModels}) = dof_residual.(aov.anovamodel.model)

"""
    predictors(model::RegressionModel)

Return a tuple of `Terms` which are predictors of the model. 

By default, it returns `formula(model).rhs.terms`; if the formula has special structures, this function should be overloaded.
"""
predictors(model::RegressionModel) = formula(model).rhs.terms

"""
    anovatable(aov::AnovaResult{<: FullModel, Test}; rownames = prednames(aov))
    anovatable(aov::AnovaResult{<: NestedModels, Test}; rownames = string.(1:N))
    anovatable(aov::AnovaResult{<: NestedModels, FTest, N}; rownames = string.(1:N)) where N
    anovatable(aov::AnovaResult{<: NestedModels, LRT, N}; rownames = string.(1:N)) where N

Return a table with coefficients and related statistics of ANOVA.
When displaying `aov` in repl, `rownames` will be `prednames(aov)` for `FullModel` and `string.(1:N)` for `NestedModels`. 

For nested models, there are two default methods for `FTest` and `LRT`; one can also define new methods dispatching on `::NestedModels{M}` where `M` is a model type. 

For a single model, no default api is implemented.

The returned `AnovaTable` object implements the [`Tables.jl`](https://github.com/JuliaData/Tables.jl/) interface, and can be  
converted e.g. to a DataFrame via `using DataFrames; DataFrame(anovatable(aov))`.
"""
function anovatable(aov::AnovaResult{T}; rownames = prednames(aov)) where {T <: FullModel}
    throw(function_arg_error(anovatable, AnovaResult{T}))
end

function anovatable(::AnovaResult{T, S, N}; rownames = string.(1:N)) where {T <: NestedModels, S, N}
    throw(function_arg_error(anovatable, AnovaResult{T}))
end

# default anovatable api for comparing multiple models
function anovatable(aov::AnovaResult{<: NestedModels, FTest, N}; rownames = string.(1:N)) where N
    AnovaTable([
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    [NaN, _diffn(deviance(aov))...], 
                    teststat(aov), 
                    pval(aov)
                ],
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "ΔDeviance", "F value", "Pr(>|F|)"],
              rownames, 7, 6)
end 

function anovatable(aov::AnovaResult{<: NestedModels, LRT, N}; rownames = string.(1:N)) where N
    AnovaTable([
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    teststat(aov), 
                    pval(aov)
                ],
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "χ²", "Pr(>|χ²|)"],
              rownames, 6, 5)
end 