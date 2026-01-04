# across different kind of models
"""
    ftest_aov(models::MultiAovModels{M, N}, df, dfr, dev, σ²) where {M <: RegressionModel, N}

Calculate F-statiscics and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedoms of each models
* `dfr`: degrees of freedom of residuals of each models
* `dev`: deviances of each models, i.e. [unit deviance](https://en.wikipedia.org/wiki/Deviance_(statistics))
* `σ²`: squared dispersion of each models

F-statiscic is `(devᵢ - devᵢ₋₁) / (dfᵢ₋₁ - dfᵢ) / σ²` for the ith predictor.
"""
function ftest_aov(models::MultiAovModels{M, N}, df, dfr, dev, σ²) where {M <: RegressionModel, N}
    length(df) ≡ length(dfr) ≡ length(dev) || throw(ArgumentError("`df`, `dfr` and `dev` must have the same length."))
    Δdf = _diff(df)
    msr = _diffn(dev) ./ Δdf
    fstat = msr ./ σ²
    pval = map(zip(Δdf, dfr[2:end], fstat)) do (dof, dofr, fs)
        fs > 0 ? ccdf(FDist(dof, dofr), fs) : NaN
    end
    AnovaResult(models, FTest, df, dev, (NaN, fstat...), (NaN, pval...), NamedTuple())
end

"""
    lrt_aov(models::MultiAovModels{M, N}, df, loglik) where {M <: RegressionModel, N}

Calculate likelihood ratio and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedom of each models
* `loglik`: loglikelihood of each models

The likelihood ratio of the ith predictor is `LRᵢ = -2loglikᵢ + 2loglikᵢ₋₁`.
"""
function lrt_aov(models::MultiAovModels{M, N}, df, loglik) where {M <: RegressionModel, N}
    length(df) ≡ length(loglik) || throw(ArgumentError("`df` and `loglik` must have the same length."))
    Δdf = _diff(df)
    lrstat = -2 .* _diffn(loglik)
    pval = map(zip(Δdf, lrstat)) do (dof, lr)
        lr > 0 ? ccdf(Chisq(dof), lr) : NaN
    end
    AnovaResult(models, LRT, df, loglik, (NaN, lrstat...), (NaN, pval...), NamedTuple())
end

"""
    slrt_aov(models::MultiAovModels{M, N}, df, loglik, σ²) where {M <: RegressionModel, N}

Calculate scaled likelihood ratio and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedom of each models
* `loglik`: loglikelihood of each models
* `σ²`: squared dispersion of each models

The scaled likelihood ratio of the ith predictor is `SLRᵢ = (- 2loglikᵢ + 2loglikᵢ₋₁) / σ²`.
"""
function slrt_aov(models::MultiAovModels{M, N}, df, loglik, σ²) where {M <: RegressionModel, N}
    length(df) ≡ length(loglik) || throw(ArgumentError("`df` and `loglik` must have the same length."))
    Δdf = _diff(df)
    lrstat = _diffn(loglik) ./ (-σ² / 2)
    pval = map(zip(Δdf, lrstat)) do (dof, lr)
        lr > 0 ? ccdf(Chisq(dof), lr) : NaN
    end
    AnovaResult(models, SLRT, df, loglik, (NaN, lrstat...), (NaN, pval...), NamedTuple())
end

# Calculate dof from assign
"""
    dof_asgn(assign::Vector{Int})

Calculate degrees of freedom of each predictors. `assign` can be obtained by `StatsModels.asgn(f::FormulaTerm)`. For a given `trm::RegressionModel`, it is as same as `trm.mm.assign`.

The index of the output matches values in the orinal `assign`. If any index value is not in `assign`, the default is 0.

# Examples
```julia
julia> dof_asgn([1, 2, 2, 3, 3, 3])
3-element Vector{Int64}:
 1
 2
 3

julia> dof_asgn([2, 2, 3, 3, 3])
3-element Vector{Int64}:
 0
 2
 3

```
"""
function dof_asgn(v::Vector{Int})
    dofv = zeros(Int, maximum(v))
    for i in v
        @inbounds dofv[i] += 1
    end
    dofv
end

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