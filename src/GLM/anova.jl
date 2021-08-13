# ===========================================================================================
# Main API

using GLM
@reexport using GLM
import GLM: glm, LinPredModel, LinearModel, LmResp, DensePred, DensePredChol, SparsePredChol, QRCompactWY, LinPred, installbeta!, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!,  AbstractGLM, FP, SparseMatrixCSC, Link, dispersion

"""
    glm(f, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...)

Automatically transform dependent variable into 0/1 for family `Binomial`
"""
glm(f::FormulaTerm, df::DataFrame, d::Binomial, l::Link, args...; kwargs...) = 
    fit(GeneralizedLinearModel, f, 
        combine(df, : , f.lhs.sym => ByRow(x -> x == unique(df[:, f.lhs.sym])[end]) => f.lhs.sym), 
        d, l, args...; kwargs...)


anova(models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N}; 
        test::Type{T} = FTest,
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

anova(models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        test::Type{T} = canonicalgoodnessoffit(models[1].model.rr.d),
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

# ==================================================================================================================
# ANOVA by F test 
# LinearModels

function anova(::Type{FTest}, 
                model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
                type::Int = 1)
    @assert (type in [1,2,3]) "Invalid type"

    assign = model.mm.assign
    ss = SS(model, type = type)
    df = dof(assign)
    push!(df, Int(nobs(model) - sum(df))) # res dof
    first(assign) == 1 || popfirst!(df)
    MSR = ss ./ df
    fstat = (MSR[1:(end-1)] / last(MSR)..., NaN)
    pvalue = (ccdf.(FDist.(df, last(df)), abs.(fstat))[1:(end-1)]..., NaN)
    AnovaResult(model, FixedAnovaStatsF{LinearModel, length(df)}(type, nobs(model), tuple(df...), ss, fstat, pvalue))
end

# ----------------------------------------------------------------------------------------
# ANOVA for genaralized linear models
# λ = -2ln(𝓛(̂θ₀)/𝓛(θ)) ~ χ²ₙ , n = difference of predictors

function anova(::Type{FTest}, 
            model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; 
            kwargs...)
    null = first(formula(model).rhs.terms) != InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova(FTest, models)
end

function anova(::Type{FTest}, 
        models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}}) where N
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    Δdev = _diffn(dev)
    mdev = Δdev ./Δdf
    σ² = dispersion(last(models).model) ^ 2
    fstat = mdev ./ σ²
    pval = ccdf.(FDist.(Δdf, dfr[2:end]), abs.(fstat))
    if first(formula(first(models)).rhs.terms) == InterceptTerm{false}()
        AnovaResult(models, FixedAnovaStatsF{GeneralizedLinearModel, length(Δdf)}(1, n, Δdf, Δdev , fstat, pval))
    else
        AnovaResult(models, FixedAnovaStatsF{GeneralizedLinearModel, 1 + length(Δdf)}(1, n, (1, Δdf...), (NaN, Δdev...) , (NaN, fstat...), (NaN, pval...)))
    end
end

# ==================================================================================================================
# ANOVA by Likehood-ratio test 
# LinearModels

function anova(::Type{LRT}, 
            model::TableRegressionModel{<: LinearModel, <: AbstractArray})
    ss = SS(model, type = 1)
    df = tuple(dof(model.mm.assign)...)
    den = last(ss) / (nobs(model) - dof(model) + 1)
    lrstat = ss[1:end - 1] ./ den
    σ² = dispersion(model.model, true)
    n = length(lrstat)
    dev = zeros(Float64, n)
    i = n - 1
    dev[end] = deviance(model)
    while i > 0
        dev[i] = σ² * lrstat[i + 1] + dev[i + 1]
        i -= 1
    end
    pval = ccdf.(Chisq.(df), abs.(lrstat))
    AnovaResult(model, FixedAnovaStatsLRT{LinearModel, n}(1, nobs(model), df, tuple(dev...), lrstat, pval))
end


# ------------------------------------------------------------------------------------------------------------
# ANOVA for GeneralizedLinearModel

function anova(::Type{LRT}, 
        model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; 
        kwargs...)
    @warn "fit all submodels"
    null = first(formula(model).rhs.terms) != InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova(LRT, models)
end

function anova(::Type{LRT}, 
        models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}}) where N
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    Δdev = _diffn(dev)
    σ² = dispersion(last(models).model, true)
    lrstat = Δdev ./ σ²
    pval = ccdf.(Chisq.(Δdf), abs.(lrstat))
    AnovaResult(models, FixedAnovaStatsLRT{GeneralizedLinearModel, length(Δdf)}(1, n, Δdf, dev[2:end], lrstat, pval))
end

# =================================================================================================================
# Nested models 

function anova(::Type{FTest}, 
        models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N}; 
        testnested::Bool = true) where N
    
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(last(models).model, true)
    fstat = (NaN, msr./σ²...)
    pval = (NaN, ccdf.(FDist.(Δdf, dfr[2:end]), abs.(fstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsF{length(df)}(n, df, dev, fstat, pval))
end

function anova(::Type{FTest}, 
        models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested::Bool = true) where N

    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(last(models).model, true)
    fstat = (NaN, msr./σ²...)
    pval = (NaN, ccdf.(FDist.(Δdf, dfr[2:end]), abs.(fstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsF{length(df)}(n, df, dev, fstat, pval))
end

function anova(::Type{LikelihoodRatioTest}, 
            models::Vararg{TableRegressionModel, N}; 
            testnested::Bool = true) where N
    # AIC and BIC
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    σ² = dispersion(last(models).model, true)
    dev = deviance.(models)
    Δdev = _diffn(dev)
    lrstat = (NaN, Δdev ./ σ² ...)
    pval = (NaN, ccdf.(Chisq.(Δdf), abs.(lrstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsLRT{length(df)}(n, df, dev, lrstat, pval))
end

# =================================================================================================================================
# Fit new models

"""
    anova_lm(X, y; test::Type{T} = FTest, <keyword arguments>) 

    anova_lm(test::Type{T}, X, y; <keyword arguments>)

    anova(test::Type{T}, ::Type{LinearModel}, X, y; 
        type::Int = 1, 
        <keyword arguments>)

ANOVA for simple linear regression.

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`. \n

* `type` specifies type of anova.
* `dropcollinear` controls whether or not lm accepts a model matrix which is less-than-full rank. If true (the default), only the first of each set of linearly-dependent columns  
is used. The coefficient for redundant linearly dependent columns is 0.0 and all associated statistics are set to NaN.

`anova_lm` generate a `TableRegressionModel` object, which is fitted by `lm`.
"""
anova_lm(X, y; 
        test::Type{T} = FTest, 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearModel, X, y; kwargs...)

anova_lm(test::Type{T}, X, y; kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearModel, X, y; kwargs...)

function anova(test::Type{T}, ::Type{LinearModel}, X, y; 
        type::Int = 1, 
        kwargs...) where {T <: GoodnessOfFit}
    model = lm(X, y; kwargs...)
    anova(test, model; type = type)
end

"""
    anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); 
            test::Type{T} = canonicalgoodnessoffit(d), <keyword arguments>)

    anova_glm(test::Type{T}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

    anova(test::Type{T}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

ANOVA for genaralized linear models.

* `d`: a `GLM.UnivariateDistribution`.
* `l`: a `GLM.Link`

For other keyword arguments, see `fit`.
"""
anova_glm(X, y, 
        d::UnivariateDistribution, l::Link = canonicallink(d); 
        test::Type{T} = canonicalgoodnessoffit(d), 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

anova_glm(test::Type{T}, X, y, 
        d::UnivariateDistribution, l::Link = canonicallink(d); 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

function anova(test::Type{T}, ::Type{GeneralizedLinearModel}, X, y, 
            d::UnivariateDistribution, l::Link = canonicallink(d);
            kwargs...) where {T <: GoodnessOfFit}

    @warn "fit all submodels"
    model = glm(X, y, d, l; kwargs...)
    null = first(formula(model).rhs.terms) != InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova(test, models)
end   
