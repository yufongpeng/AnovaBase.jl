# ============================================================================================================
# Main API

using MixedModels
@reexport using MixedModels
<<<<<<< Updated upstream
import MixedModels: FeMat, createAL, reweight!, getθ
=======
import MixedModels: FeMat, createAL, reweight!, getθ,
                     _iscomparable, _criterion,
                     deviance, dof, dof_residual, nobs
import StatsModels: width
>>>>>>> Stashed changes

"""
    lme(f::FormulaTerm, tbl; wts, contrasts, verbose, REML)

An alias for `fit(LinearMixedModel, f, tbl; wts, contrasts, verbose, REML)`.

"""
lme(f::FormulaTerm, tbl; 
    wts = [], 
    contrasts = Dict{Symbol, Any}(), 
    verbose::Bool = false, 
    REML::Bool = false) = 
    fit(LinearMixedModel, f, tbl, 
    wts =  wts, contrasts = contrasts, verbose = verbose, REML = REML)

anova(models::Vararg{LinearMixedModel, N}; 
        test::Type{T} = length(models) > 1 ? (LRT) : (FTest), 
        kwargs...) where {N, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

# ==================================================================================================================
# ANOVA by F test
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
<<<<<<< Updated upstream
        invvarfixchol = cholesky(inv(varβ)|> Hermitian).L 
        # adjust σ like linear regression
        model.optsum.REML || adjust_sigma && begin
            invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - length(β)))
        end 
        fs = invvarfixchol'β
        uniqas = unique(assign)
        fstat = ntuple(lastindex(uniqas)) do fac
            mapreduce(val->fs[val] ^ 2, +, findall(==(fac), assign)) / df[fac]
=======
        fs = abs2.(cholesky(Hermitian(inv(varβ))).U  * β)
        # adjust σ like linear regression
        adjust = 1.0
        offset = first(assign) - 1
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) 
        fstat = ntuple(last(assign) - offset) do fix
            sum(fs[findall(==(fix + offset), assign)]) / df[fix] * adjust
>>>>>>> Stashed changes
        end
    else 
        # calculate block by block
        adjust = 1.0
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) 
<<<<<<< Updated upstream
        offset = 0
        intercept || (offset = 1)
        fstat = ntuple(last(assign) - offset) do factor
            select = findall(==(factor + offset), assign)
            invvarfix = inv(varβ[select, select]) 
            view(β, select)' * invvarfix * view(β, select) / rank(invvarfix) * adjust
=======
        offset = first(assign) - 1
        fstat = ntuple(last(assign) - offset) do fix
            select = findall(==(fix + offset), assign)
            β[select]' * (varβ[select, select] \ β[select]) / df[fix] * adjust
>>>>>>> Stashed changes
        end
    end

    pvalue = ntuple(lastindex(fstat)) do id
            ccdf(FDist(df[id], resdf[id]), abs(fstat[id]))
    end
    AnovaResult(model, MixedAnovaStatsF{LinearMixedModel, length(fstat)}(type, nobs(model), df, resdf, fstat, pvalue, dofinfo))
end

# ==================================================================================================================
# ANOVA by Likehood-ratio test 
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

# =================================================================================================================
# Nested models 

function anova(::Type{FTest}, 
    models::Vararg{<: LinearMixedModel, N}; 
    testnested::Bool = true) where N

@warn "F test can only be appplied to a single linear mixed-effects model"
aov = anova(FTest, last(models))
AnovaResult(models, NestedAnovaStatsF{length(models)}(aov.stats.nobs, dof.(models), deviance.(models), (NaN, aov.stats.fstat...), (NaN, aov.stats.pval...)))
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
    anova_lme(f::FormulaTerm, tbl; test::Type{T} = FTest, <keyword arguments>)

    anova_lme(test::Type{T}, f::FormulaTerm, tbl; <keyword arguments>)

    anova(test::Type{T}, ::Type{LinearMixedModel}, f::FormulaTerm, tbl;
            type::Int = 1, 
            adjust_sigma::Bool = true)

ANOVA for linear mixed-effect models.

The arguments `f` and `tbl` are `Formula` and `DataFrame`.

<<<<<<< Updated upstream
* `type` specifies type of anova. only `1, 3` are valid.
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fit by REML.
=======
* `test`: `GoodnessOfFit`. The default is `FTest`.
* `type`: type of anova. Only 1, 3 are valid.
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fit by REML.  
    The result will be slightly deviated from that of model fit by REML.
>>>>>>> Stashed changes

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
<<<<<<< Updated upstream
        verbose::Bool = false, 
        REML::Bool = true) where {T <: GoodnessOfFit}
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    anova(test, model; type = type, adjust_sigma = adjust_sigma)
end
=======
        progress::Bool = true, 
        REML::Bool = test == FTest ? true : false)
    model = lme(f, tbl; wts, contrasts, progress, REML)
    anova(test, model; type, adjust_sigma)
end

"""
    lme(f::FormulaTerm, tbl; wts, contrasts, progress, REML)

An alias for `fit(LinearMixedModel, f, tbl; wts, contrasts, progress, REML)`.
"""
lme(f::FormulaTerm, tbl; 
    wts = [], 
    contrasts = Dict{Symbol, Any}(), 
    progress::Bool = true, 
    REML::Bool = false) = 
    fit(LinearMixedModel, f, tbl; 
        wts, contrasts, progress, REML)

"""
    glme(f::FormulaTerm, tbl, d::Distribution, l::Link; kwargs...)

An alias for `fit(GeneralizedLinearMixedModel, f, tbl, d, l; kwargs...)`
"""
glme(f::FormulaTerm, tbl, d::Distribution = Normal(), l::Link = canonicallink(d);
    kwargs...) = 
    fit(GeneralizedLinearMixedModel, f, tbl, d, l; kwargs...)
>>>>>>> Stashed changes
