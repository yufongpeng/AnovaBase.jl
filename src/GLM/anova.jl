# ===========================================================================================
# Main API

using GLM
@reexport using GLM
import GLM: glm, 
            # Model
            LinPredModel, AbstractGLM, GeneralizedLinearModel, LinearModel, 
            LmResp, GlmResp, 
            # Pred
            LinPred, DensePred, 
            DensePredChol, SparsePredChol, QRCompactWY, SparseMatrixCSC, 
            # prediction
            installbeta!, delbeta!, linpred, linpred!,
            updateμ!, cholfactors, 
            # other
            FP, BlasReal, Link, dispersion, deviance, dof, dof_residual, nobs

anova(trms::Vararg{TableRegressionModel{<: LinearModel}}; 
        test::Type{<: GoodnessOfFit} = FTest,
        kwargs...) = 
    anova(test, trms...; kwargs...)

anova(trms::Vararg{TableRegressionModel{<: GeneralizedLinearModel}}; 
        test::Type{<: GoodnessOfFit} = canonicalgoodnessoffit(trms[1].model.rr.d),
        kwargs...) = 
    anova(test, trms...; kwargs...)

# ==================================================================================================================
# ANOVA by F test 

anova(::Type{FTest}, 
    trm::TableRegressionModel{<: Union{LinearModel, GeneralizedLinearModel{<: GLM.GlmResp{T, <: Normal, IdentityLink}}}}; 
    type::Int = 1,
    kwargs...) where T = _anova_vcov(trm; type, kwargs...)

function _anova_vcov(trm::TableRegressionModel{<: LinPredModel}; 
                    type::Int = 1, kwargs...)
    @assert (type in [1,2,3]) "Invalid type"

    assign = trm.mm.assign
    df = dof(assign)
    filter!(>(0), df)
    push!(df, Int(dof_residual(trm)))
    df = tuple(df...)
    if type in [1, 3] 
        # vcov methods
        varβ = vcov(trm.model)
        β = trm.model.pp.beta0
        if type == 1
            fs = abs2.(cholesky(Hermitian(inv(varβ))).L' * β) 
            fstat = ntuple(lastindex(unique(assign))) do fix
                sum(fs[findall(==(fix), assign)]) / df[fix]
            end
        else
            # calculate block by block
            offset = first(assign) == 1 ? 0 : 1
            fstat = ntuple(last(assign) - offset) do fix
                select = findall(==(fix + offset), assign)
                β[select]' * inv(varβ[select, select]) * β[select] / df[fix]
            end
        end
        σ² = dispersion(trm.model, true)
        devs = (fstat .* σ²..., σ²) .* df
    else
        # refit methods
        devs = deviances(trm; type, kwargs...)
        MSR = devs ./ df
        fstat = MSR[1:end - 1] ./ dispersion(trm.model, true)
    end
    pvalue = (ccdf.(FDist.(df[1:end - 1], last(df)), abs.(fstat))..., NaN)
    AnovaResult{FTest}(trm, type, df, devs, (fstat..., NaN), pvalue, NamedTuple())
end


function anova(::Type{FTest}, 
                trm::TableRegressionModel{<: LinPredModel}; 
                type::Int = 1, kwargs...)
    @assert (type in [1,2,3]) "Invalid type"

    assign = trm.mm.assign
    devs = deviances(trm; type, kwargs...)
    df = dof(assign)
    filter!(>(0), df)
    push!(df, Int(dof_residual(trm)))
    length(df) == length(devs) + 1 && popfirst!(df)
    df = tuple(df...)
    MSR = devs ./ df
    fstat = MSR[1:end - 1] ./ dispersion(trm.model, true)
    pvalue = (ccdf.(FDist.(df[1:end - 1], last(df)), abs.(fstat))..., NaN)
    AnovaResult{FTest}(trm, type, df, devs, (fstat..., NaN), pvalue, NamedTuple())
end

# ==================================================================================================================
# ANOVA by Likehood-ratio test 
# λ = -2ln(𝓛(̂θ₀)/𝓛(θ)) ~ χ²ₙ , n = difference of predictors

function anova(::Type{LRT}, 
            trm::TableRegressionModel{<: LinPredModel})
    Δdev = deviances(trm, type = 1)
    df = dof(trm.mm.assign)
    filter!(>(0), df)
    isnullable(trm.model) || popfirst!(df)
    df = tuple(df...)
    # den = last(ss) / (nobs(trm) - dof(trm) + 1)
    # lrstat = ss[1:end - 1] ./ den
    σ² = dispersion(trm.model, true)
    lrstat = Δdev[1:end - 1] ./ σ²
    n = length(lrstat)
    dev = collect(Δdev)
    i = n
    while i > 0
        dev[i] += dev[i + 1]
        i -= 1
    end
    pval = ccdf.(Chisq.(df), abs.(lrstat))
    AnovaResult{LRT}(trm, 1, df, tuple(dev[2:end]...), lrstat, pval, NamedTuple())
end

# =================================================================================================================
# Nested models 

function anova(::Type{FTest}, 
        trms::Vararg{TableRegressionModel{<: LinPredModel}}; 
        check::Bool = true,
        isnested::Bool = false)
    df = dof.(trms)
    ord = sortperm(collect(df))
    df = df[ord]
    trms = trms[ord]

    # check comparable and nested
    check && @warn "Could not check whether models are nested: results may not be meaningful"

    n = Int(nobs(first(trms)))
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(trms))
    dev = deviance.(trms)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(last(trms).model, true)
    fstat = msr ./ σ²
    pval = map(zip(Δdf, dfr[2:end], fstat)) do (dof, dofr, fs)
        fs > 0 ? ccdf(FDist(dof, dofr), fs) : NaN
    end
    AnovaResult{FTest}(trms, 1, df, dev, (NaN, fstat...), (NaN, pval...), NamedTuple())
end

function anova(::Type{LRT}, 
        trms::Vararg{<: TableRegressionModel}; 
        check::Bool = true,
        isnested::Bool = false)
    df = dof.(trms)
    ord = sortperm(collect(df))
    trms = trms[ord]
    # check comparable and nested
    df = df[ord]
    _lrt_nested(trms, df, deviance.(trms), dispersion(last(trms).model, true); nestedwarn = true)
end


# =================================================================================================================================
# Fit new models

"""
    anova_lm(X, y; test::Type{<: GoodnessOfFit} = FTest, <keyword arguments>) 

    anova_lm(test::Type{<: GoodnessOfFit}, X, y; <keyword arguments>)

    anova(test::Type{<: GoodnessOfFit}, ::Type{LinearModel}, X, y; 
        type::Int = 1, 
        dropcollinear::Bool = true,
        <keyword arguments>)

ANOVA for simple linear regression.

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`. \n

* `type`: type of anova.
* `dropcollinear`: whether or not lm accepts a model matrix which is less-than-full rank. If true (default), only the first of each set of linearly-dependent columns  
is used. The coefficient for redundant linearly dependent columns is 0.0 and all associated statistics are set to NaN.

`anova_lm` generate a `TableRegressionModel` object, which is fitted by `lm`.
"""
anova_lm(X, y; 
        test::Type{<: GoodnessOfFit} = FTest, 
        kwargs...)= 
    anova(test, LinearModel, X, y; kwargs...)

anova_lm(test::Type{<: GoodnessOfFit}, X, y; kwargs...) = 
    anova(test, LinearModel, X, y; kwargs...)

function anova(test::Type{<: GoodnessOfFit}, ::Type{LinearModel}, X, y; 
        type::Int = 1, 
        dropcollinear::Bool = true,
        kwargs...)
    trm = lm(X, y; dropcollinear, kwargs...)
    anova(test, trm; type)
end

"""
    anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); 
            test::Type{<: GoodnessOfFit} = canonicalgoodnessoffit(d), <keyword arguments>)

    anova_glm(test::Type{<: GoodnessOfFit}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

    anova(test::Type{<: GoodnessOfFit}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

ANOVA for genaralized linear models.

* `d`: a `GLM.UnivariateDistribution`.
* `l`: a `GLM.Link`

For other keyword arguments, see `fit`.
"""
anova_glm(X, y, 
        d::UnivariateDistribution, l::Link = canonicallink(d); 
        test::Type{<: GoodnessOfFit} = canonicalgoodnessoffit(d), 
        kwargs...) = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

anova_glm(test::Type{<: GoodnessOfFit}, X, y, 
        d::UnivariateDistribution, l::Link = canonicallink(d); 
        kwargs...) = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

function anova(test::Type{<: GoodnessOfFit}, ::Type{GeneralizedLinearModel}, X, y, 
            d::UnivariateDistribution, l::Link = canonicallink(d);
            type::Int = 1,
            kwargs...)
    trm = glm(X, y, d, l; kwargs...)
    anova(test, trm; type, kwargs... )
end 


"""
    glm(f, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...)

Automatically transform dependent variable into 0/1 for family `Binomial`.
"""
glm(f::FormulaTerm, df::DataFrame, d::Binomial, l::Link, args...; kwargs...) = 
    fit(GeneralizedLinearModel, f, 
        combine(df, : , f.lhs.sym => ByRow(x -> x == unique(df[!, f.lhs.sym])[end]) => f.lhs.sym), 
        d, l, args...; kwargs...)
