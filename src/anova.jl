# ===========================================================================================
# Main API

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
anova(models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N}; 
        test::Type{T} = FTest,
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

anova(models::Vararg{LinearMixedModel, N}; 
        test::Type{T} = length(models) > 1 ? (LRT) : (FTest), 
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

anova(models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        test::Type{T} = canonicalgoodnessoffit(models[1].model.rr.d),
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

# ==================================================================================================================
# ANOVA by F test 

"""
    anova(::Type{FTest}, <model>; kwargs...)
    anova(::Type{FTest}, <models>...; kwargs...)

Analysis of Variance by F test.

* `type` specifies type of anova. For one `LinearModel` `1, 2, 3` are valid; for one `LinearMixedModel` `1, 3` are valid. For others, only `1` is valid.
* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
function anova(::Type{FTest}, 
                model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
                type::Int = 1, 
                pivot::Bool = false)
    @assert (type in [1,2,3]) "Invalid type"

    assign = model.mm.assign
    ss = SS(model, type = type, pivot = pivot)
    df = dof(assign)
    push!(df, Int(nobs(model) - sum(df))) # res dof
    first(assign) == 1 || popfirst!(df)
    MSR = ss ./ df
    fstat = (MSR[1:(end-1)] / last(MSR)..., NaN)
    pvalue = (ccdf.(FDist.(df, last(df)), abs.(fstat))[1:(end-1)]..., NaN)
    AnovaResult(model, FixedAnovaStatsF{LinearModel, length(df)}(type, nobs(model), tuple(df...), ss, fstat, pvalue))
end

# ----------------------------------------------------------------------------------------------------
# Linear mixed-effect models

function anova(::Type{FTest}, 
            model::LinearMixedModel; 
            type::Int = 1, 
            adjust_sigma::Bool = true)

    @assert (type in [1,2,3]) "Invalid type"
    @assert (type in [1,3]) "Type 2 anova is not supported now"

    varβ = vcov(model) 
    β = fixef(model)

    assign = asgn(first(model.formula.rhs))
    intercept = first(assign) == 1

    # calculate degree of freedom for factors and residuals
    df, resdf, dofinfo = calcdof(model)
    
    # use MMatrix/SizedMatrix ?
    if type == 1
        invvarfixchol = cholesky(inv(varβ)|> Hermitian).L 
        # adjust σ like linear regression
        model.optsum.REML || adjust_sigma && begin
            invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - length(β)))
        end 
        fs = invvarfixchol'β
        uniqas = unique(assign)
        fstat = ntuple(lastindex(uniqas)) do fac
            mapreduce(val->fs[val] ^ 2, +, findall(==(fac), assign)) / df[fac]
        end
    else 
        # calculate block by block
        adjust = 1.0
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) 
        offset = 0
        intercept || (offset = 1)
        fstat = ntuple(last(assign) - offset) do factor
            select = findall(==(factor + offset), assign)
            invvarfix = inv(varβ[select, select]) 
            view(β, select)' * invvarfix * view(β, select) / rank(invvarfix) * adjust
        end
    end

    pvalue = ntuple(lastindex(fstat)) do id
            ccdf(FDist(df[id], resdf[id]), abs(fstat[id]))
    end
    AnovaResult(model, MixedAnovaStatsF{LinearMixedModel, length(fstat)}(type, nobs(model), df, resdf, fstat, pvalue, dofinfo))
end

# ----------------------------------------------------------------------------------------
# ANOVA for genaralized linear models
# λ = -2ln(𝓛(̂θ₀)/𝓛(θ)) ~ χ²ₙ , n = difference of predictors

function anova(::Type{FTest}, 
            model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; 
            kwargs...)
    null = first(formula(model).rhs.terms) == InterceptTerm{false}()
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

"""
    anova(::Type{LRT}, <model>; kwargs...)
    anova(::Type{LRT}, <models>...; kwargs...)

Analysis of Variance by likelihood-ratio test.

* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
function anova(::Type{LRT}, 
            model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
            pivot::Bool = false)
    ss = SS(model, type = 1, pivot = pivot)
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

# ----------------------------------------------------------------------------------------------------------
# ANOVA for LinearMixedModel

function anova(::Type{LRT}, model::LinearMixedModel)
    @warn "fit all submodels"
    null = first(first(formula(model).rhs).terms) == InterceptTerm{false}()
    models = nestedmodels(model; null = null)
    anova(LRT, models)
end

function anova(::Type{LRT}, models::NTuple{N, LinearMixedModel}) where N
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dev = deviance.(models)
    lrstat = _diffn(dev)
    pval = ccdf.(Chisq.(abs.(Δdf)), abs.(lrstat))
    AnovaResult(models, MixedAnovaStatsLRT{LinearMixedModel, length(Δdf)}(1, n, Δdf, dev[2:end], lrstat, pval))
end

# ------------------------------------------------------------------------------------------------------------
# ANOVA for GeneralizedLinearModel

function anova(::Type{LRT}, 
        model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; 
        kwargs...)
    @warn "fit all submodels"
    null = first(formula(model).rhs.terms) == InterceptTerm{false}()
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
        models::Vararg{<: LinearMixedModel, N}; 
        testnested::Bool = true) where N
    
    @warn "F test can only be appplied to a single linear mixed-effects model"
    aov = anova(FTest, last(models))
    AnovaResult(models, NestedAnovaStatsF{length(models)}(aov.stats.nobs, dof.(models), deviance.(models), (NaN, aov.stats.fstat...), (NaN, aov.stats.pval...)))
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

function anova(::Type{LikelihoodRatioTest}, 
            models::Vararg{<: LinearMixedModel, N}; 
            testnested::Bool = true) where N
    # AIC and BIC
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dev = deviance.(models)
    Δdev = _diffn(dev)
    lrstat = (NaN, Δdev...)
    pval = (NaN, ccdf.(Chisq.(abs.(Δdf)), abs.(lrstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsLRT{length(df)}(n, df, dev, lrstat, pval))
end

# =================================================================================================================================
# Fit new models

"""
    anova_lm(X, y; test::Type{T} = FTest, <keyword arguments>) 

    anova_lm(test::Type{T}, X, y; <keyword arguments>)

    anova(test::Type{T}, ::Type{LinearModel}, X, y; 
        pivot::Bool = false, 
        type::Int = 1, 
        <keyword arguments>)

ANOVA for simple linear regression.

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`. \n

* `type` specifies type of anova.
* `pivot` determines if pivot is used, if modelmatrix is rank deficient.

`anova_lm` generate a `TableRegressionModel` object, which is fitted by `lm`.
"""
anova_lm(X, y; 
        test::Type{T} = FTest, 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearModel, X, y; kwargs...)

anova_lm(test::Type{T}, X, y; kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearModel, X, y; kwargs...)

function anova(test::Type{T}, ::Type{LinearModel}, X, y; 
        pivot::Bool = false, 
        type::Int = 1, 
        kwargs...) where {T <: GoodnessOfFit}
    model = lm(X, y, pivot; kwargs...)
    anova(test, model; pivot = pivot, type = type)
end

"""
    anova_lme(f::FormulaTerm, tbl; test::Type{T} = FTest, <keyword arguments>)

    anova_lme(test::Type{T}, f::FormulaTerm, tbl; <keyword arguments>)

    anova(test::Type{T}, ::Type{LinearMixedModel}, f::FormulaTerm, tbl;
            type::Int = 1, 
            adjust_sigma::Bool = true)

ANOVA for linear mixed-effect models.

The arguments `f` and `tbl` are `Formula` and `DataFrame`.

* `type` specifies type of anova. only `1, 3` are valid.
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fit by REML.

Other keyword arguments
* `wts = []`
* `contrasts = Dict{Symbol,Any}()`
* `verbose::Bool = false`
* `REML::Bool = true`

`anova_lme` generate a `LinearMixedModel` object through calling `anova`, which is fit by `lme` with REML.
"""
anova_lme(f::FormulaTerm, tbl; 
        test::Type{T} = FTest,
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearMixedModel, f, tbl; kwargs...)

anova_lme(test::Type{T}, f::FormulaTerm, tbl; 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearMixedModel, f, tbl; kwargs...)

function anova(test::Type{T}, ::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
        type::Int = 1, 
        adjust_sigma::Bool = true,
        wts = [], 
        contrasts = Dict{Symbol,Any}(), 
        verbose::Bool = false, 
        REML::Bool = true) where {T <: GoodnessOfFit}
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    anova(test, model; type = type, adjust_sigma = adjust_sigma)
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
    null = first(formula(model).rhs.terms) == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova(test, models)
end   