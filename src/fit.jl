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

# across different kind of models
"""
    ftest_nested(models::NestedModels{M, N}, df, dfr, dev, σ²) where {M <: RegressionModel, N}

Calculate F-statiscics and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedoms of each models
* `dfr`: degrees of freedom of residuals of each models
* `dev`: deviances of each models, i.e. [unit deviance](https://en.wikipedia.org/wiki/Deviance_(statistics))
* `σ²`: squared dispersion of each models

F-statiscic is `(devᵢ - devᵢ₋₁) / (dfᵢ₋₁ - dfᵢ) / σ²` for the ith predictor.
"""
function ftest_nested(models::NestedModels{M, N}, df, dfr, dev, σ²) where {M <: RegressionModel, N}
    length(df) == length(dfr) == length(dev) || throw(ArgumentError("`df`, `dfr` and `dev` must have the same length."))
    Δdf = _diff(df)
    msr = _diffn(dev) ./ Δdf
    fstat = msr ./ σ²
    pval = map(zip(Δdf, dfr[2:end], fstat)) do (dof, dofr, fs)
        fs > 0 ? ccdf(FDist(dof, dofr), fs) : NaN
    end
    AnovaResult{FTest}(models, df, dev, (NaN, fstat...), (NaN, pval...), NamedTuple())
end

"""
    lrt_nested(models::NestedModels{M, N}, df, dev, σ²) where {M <: RegressionModel, N}

Calculate likelihood ratio and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedom of each models
* `dev`: deviances of each models, i.e. [unit deviance](https://en.wikipedia.org/wiki/Deviance_(statistics))
* `σ²`: squared dispersion of each models

The likelihood ratio of the ith predictor is `LRᵢ = (devᵢ - devᵢ₋₁) / σ²`.

If `dev` is alternatively `-2loglikelihood`, `σ²` should be set to 1.
"""
function lrt_nested(models::NestedModels{M, N}, df, dev, σ²) where {M <: RegressionModel, N}
    length(df) == length(dev) || throw(ArgumentError("`df` and `dev` must have the same length."))
    Δdf = _diff(df)
    Δdev = _diffn(dev)
    lrstat = Δdev ./ σ²
    pval = map(zip(Δdf, lrstat)) do (dof, lr)
        lr > 0 ? ccdf(Chisq(dof), lr) : NaN
    end
    AnovaResult{LRT}(models, df, dev, (NaN, lrstat...), (NaN, pval...), NamedTuple())
end

# Calculate dof from assign
"""
    dof_asgn(v::Vector{Int})

Calculate degrees of freedom of each predictors. 'v' can be obtained by `StatsModels.asgn(f::FormulaTerm)`. For a given `trm::RegressionModel`, it is as same as `trm.mm.assign`.
"""
function dof_asgn(v::Vector{Int})
    dofv = zeros(Int, maximum(v))
    for i in v
        @inbounds dofv[i] += 1
    end
    dofv
end

@deprecate dof(v::Vector{Int}) dof_asgn(v::Vector{Int})

const FixDispDist = Union{Bernoulli, Binomial, Poisson}
"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return `LRT` if the distribution has a fixed dispersion; otherwise, `FTest`.
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest

"""
    _diff(t::NTuple)

Return a tuple of difference between adjacent elements of a tuple(later - former). 
"""
_diff(t::NTuple{N, T}) where {N, T} = ntuple(i -> t[i + 1] - t[i], N - 1)

"""
    _diff(t::NTuple)

Return a tuple of difference between adjacent elements of a tuple(former - later). 
"""
_diffn(t::NTuple{N, T}) where {N, T} = ntuple(i -> t[i] - t[i + 1], N - 1)