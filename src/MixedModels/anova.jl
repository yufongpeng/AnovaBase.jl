# ============================================================================================================
# Main API

using MixedModels, GLM
@reexport using MixedModels
import MixedModels: FeMat, createAL, reweight!, getθ,
                     _iscomparable, _criterion,
                     deviance, dof, dof_residual, nobs

anova(models::Vararg{<: LinearMixedModel}; 
        test::Type{<: GoodnessOfFit} = length(models) > 1 ? LRT : FTest, 
        kwargs...) = 
    anova(test, models...; kwargs...)

# ==================================================================================================================
# ANOVA by F test
# Linear mixed-effect models

function anova(::Type{FTest}, 
        model::LinearMixedModel; 
        type::Int = 1, 
        adjust_sigma::Bool = true,
        kwargs...)

    @assert (type in [1,2,3]) "Invalid type"
    @assert (type in [1,3]) "Type 2 anova is not supported now"

    varβ = vcov(model) 
    β = fixef(model)

    assign = asgn(first(model.formula.rhs))

    # calculate degree of freedom for factors and residuals
    df, resdf = calcdof(model)

    # use MMatrix/SizedMatrix ?
    if type == 1
        invvarfixchol = cholesky(Hermitian(inv(varβ))).L 
        fs = abs2.(invvarfixchol'β)
        # adjust σ like linear regression
        adjust = 1.0
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) 
        fstat = ntuple(lastindex(unique(assign))) do fix
            sum(fs[findall(==(fix), assign)]) / df[fix] * adjust
        end
    else 
        # calculate block by block
        adjust = 1.0
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) 
        offset = 0
        first(assign) == 1 || (offset = 1)
        fstat = ntuple(last(assign) - offset) do fix
            select = findall(==(fix + offset), assign)
            β[select]' * inv(varβ[select, select])  * β[select] / df[fix] * adjust
        end
    end

    pvalue = ntuple(lastindex(fstat)) do id
            ccdf(FDist(df[id], resdf[id]), abs(fstat[id]))
    end
    AnovaResult{FTest}(model, type, df, ntuple(x->NaN, length(fstat)), fstat, pvalue, (resdof = resdf,))
end

# ==================================================================================================================
# ANOVA by Likehood-ratio test 
# Linear mixed-effect models

function anova(::Type{LRT}, model::LinearMixedModel; kwargs...)
    # check if fitted by ML 
    # nested random effects for REML ?
    model.optsum.REML && throw(ArgumentError("Likelihood-ratio tests for REML-fitted models are only valid when the fixed-effects specifications are identical"))
    @warn "fit all submodels"
    models = nestedmodels(model; null = isnullable(model))
    anova(LRT, models...; check = false, isnested = true)
end

# =================================================================================================================
# Nested models 
#= 
function anova(::Type{FTest}, 
                models::Vararg{<: LinearMixedModel, N}; 
                testnested::Bool = true, 
                warn::Bool = true) where N
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dev = deviance.(models)
    Δdev = _diffn(dev)
    AnovaResult(models, NestedAnovaStatsF{length(models)}(aov.stats.nobs, df, dev, (NaN, aov.stats.fstat...), (NaN, aov.stats.pval...)))
end
=#
function anova(::Type{LikelihoodRatioTest}, 
                models::Vararg{<: LinearMixedModel}; 
                check::Bool = true,
                isnested::Bool = false,
                kwargs...)

    check && (_iscomparable(models...) || throw(
        ArgumentError("""Models are not comparable: are the objectives, data
                         and, where appropriate, the link and family the same?
        """)))
    # isnested is not part of _iscomparable:  
    # isnested = true  
    df = dof.(models)
    ord = sortperm(collect(df))
    df = df[ord]
    models = models[ord]
    _lrt_nested(models, df, deviance.(models), 1; nestedwarn = isnested)
end

# Compare to GLM
function anova(m0::Union{TableRegressionModel{<: Union{LinearModel, GeneralizedLinearModel}}, LinearModel, GeneralizedLinearModel},
                m::T,
                ms::Vararg{T};
                check::Bool = true,
                isnested::Bool = false,
                kwargs...) where {T <: MixedModel}

    check && (_iscomparable(m0, m) ||
        throw(ArgumentError("""Models are not comparable: are the objectives, data
                                and, where appropriate, the link and family the same?
                            """)))
    check && (_iscomparable(m, ms...) || throw(
        ArgumentError("""Models are not comparable: are the objectives, data
                        and, where appropriate, the link and family the same?
        """)))
    m = [m, ms...]
    df = dof.(m)
    ord = sortperm(df)
    df = (dof(m0), df[ord]...)
    models = (m0, m[ord]...)
    # isnested is not part of _iscomparable:  
    # isnested = true 
    dev = (_criterion(m0), deviance.(models[2:end])...)
    _lrt_nested(models, df, dev, 1; nestedwarn = isnested)
end

# =================================================================================================================================
# Fit new models

"""
    anova_lme(f::FormulaTerm, tbl; test::Type{<: GoodnessOfFit} = FTest, <keyword arguments>)

    anova_lme(test::Type{<: GoodnessOfFit}, f::FormulaTerm, tbl; <keyword arguments>)

    anova(test::Type{<: GoodnessOfFit}, ::Type{<: LinearMixedModel}, f::FormulaTerm, tbl;
            type::Int = 1, 
            adjust_sigma::Bool = true, <keyword arguments>)

ANOVA for linear mixed-effect models.

The arguments `f` and `tbl` are `Formula` and `DataFrame`.

* `test`: `GoodnessOfFit`. The default is `FTest`.
* `type`: type of anova. Only 1, 3 are valid.
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fit by REML. \n 
    The result will be slightly deviated from that of model fit by REML.

Other keyword arguments
* `wts = []`
* `contrasts = Dict{Symbol,Any}()`
* `progress::Bool = true`
* `REML::Bool = true`

`anova_lme` generate a `LinearMixedModel` object which is fit by `lme` with REML if applying `FTest`; otherwise with ML.
"""
anova_lme(f::FormulaTerm, tbl; 
        test::Type{<: GoodnessOfFit} = FTest,
        kwargs...) = 
    anova(test, LinearMixedModel, f, tbl; kwargs...)

anova_lme(test::Type{<: GoodnessOfFit}, f::FormulaTerm, tbl; 
        kwargs...) = 
    anova(test, LinearMixedModel, f, tbl; kwargs...)

function anova(test::Type{<: GoodnessOfFit}, ::Type{<: LinearMixedModel}, f::FormulaTerm, tbl; 
        type::Int = 1, 
        adjust_sigma::Bool = true,
        wts = [], 
        contrasts = Dict{Symbol,Any}(), 
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