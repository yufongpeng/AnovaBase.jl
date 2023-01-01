"""
    _diff(t::NTuple)

Return a tuple of difference between adjacent elements of a tuple(later - former). 
"""
_diff(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i + 1] - t[i], N - 1)

"""
    _diff(t::NTuple)

Return a tuple of difference between adjacent elements of a tuple(former - later). 
"""
_diffn(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i] - t[i + 1], N - 1)

# across different kind of models
"""
    ftest_nested(models::NTuple{N, RegressionModel}, df, dfr, dev, σ²) where N

Calculate F-statiscics and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedoms of each models
* `dfr`: degrees of freedom of residuals of each models
* `dev`: deviances of each models
* `σ²`: squared dispersion of each models

F-statiscic is `(devᵢ - devᵢ₋₁) / (dfᵢ₋₁ - dfᵢ) / σ²` for the ith factor.
"""
function ftest_nested(models::NTuple{N, RegressionModel}, df, dfr, dev, σ²) where N
    length(df) == length(dfr) == length(dev) || throw(ArgumentError("`df`, `dfr` and `dev` must have the same length."))
    Δdf = _diff(df)
    msr = _diffn(dev) ./ Δdf
    fstat = msr ./ σ²
    pval = map(zip(Δdf, dfr[2:end], fstat)) do (dof, dofr, fs)
        fs > 0 ? ccdf(FDist(dof, dofr), fs) : NaN
    end
    AnovaResult{FTest}(models, 1, df, dev, (NaN, fstat...), (NaN, pval...), NamedTuple())
end

"""
    lrt_nested(models::NTuple{N, RegressionModel}, df, dev, σ²) where N

Calculate likelihood ratio and p-values based on given parameters.

* `models`: nested models 
* `df`: degrees of freedoms of each models
* `dev`: deviances of each models
* `σ²`: squared dispersion of each models

The likelihood ratio of the ith factor is `LRᵢ = (devᵢ - devᵢ₋₁) / σ²`.

If `dev` is alternatively `-2loglikelihood`, `σ²` should be 1.
"""
function lrt_nested(models::NTuple{N, RegressionModel}, df, dev, σ²) where N
    length(df) == length(dev) || throw(ArgumentError("`df` and `dev` must have the same length."))
    Δdf = _diff(df)
    Δdev = _diffn(dev)
    lrstat = Δdev ./ σ²
    pval = map(zip(Δdf, lrstat)) do (dof, lr)
        lr > 0 ? ccdf(Chisq(dof), lr) : NaN
    end
    AnovaResult{LRT}(models, 1, df, dev, (NaN, lrstat...), (NaN, pval...), NamedTuple())
end

# Calculate dof from assign
"""
    dof(v::Vector{Int})

Calculate degrees of freedom of each predictors. For a given `trm::RegressionModel`, `v` is from `trm.mm.assign` and must be a non-decreasing array of integers.
"""
function dof(v::Vector{Int})
    dofv = zeros(Int, v[end])
    prev = 1
    ind = 1
    n = length(v)
    while ind <= n
        v[ind] == prev || (prev = v[ind])
        dofv[prev] += 1
        ind += 1
    end
    dofv
end

const FixDispDist = Union{Bernoulli, Binomial, Poisson}
"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return LRT if the distribution has a fixed dispersion.
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest