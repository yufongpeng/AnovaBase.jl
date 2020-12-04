# ===========================================================================================
# Main algorithm 

@doc """
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
""" anova

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


"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return LRT if the distribution has fixed dispersion
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest


# ==================================================================================================================
# ANOVA by F test 

"""
    anova(::Type{FTest}, models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N};
        testnested::Bool = true,
        type::Int = 1, 
        pivot::Bool = false)

    anova(::Type{FTest}, models::Vararg{LinearMixedModel, N}; 
        testnested::Bool = true,
        type::Int = 1, 
        adjust_sigma::Bool = true)
    
    anova(::Type{FTest}, models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested ::Bool = true,
        kwargs...)

Analysis of Variance bt F test.

* `type` specifies type of anova. For one `LinearModel` `1, 2, 3` are valid; for one `LinearMixedModel` `1, 3` are valid. For others, only `1` is valid.
* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""
anova(::Type{FTest}, models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N};
        testnested::Bool = true,
        type::Int = 1, 
        pivot::Bool = false) where N = 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; type = type, pivot = pivot)

anova(::Type{FTest}, models::Vararg{LinearMixedModel, N}; 
        testnested::Bool = true,
        type::Int = 1, 
        adjust_sigma::Bool = true) where N = 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; type = type, adjust_sigma = adjust_sigma)

    
anova(::Type{FTest}, models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested ::Bool = true,
        kwargs...) where N = 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; kwargs...)

# ------------------------------------------------------------------------------------------------------------------
# Simple linear regression

function anova1(::Type{FTest}, model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
                type::Int = 1, 
                pivot::Bool = false) where T 
    @assert (type in [1,2,3]) "Invalid type"
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df, Int(size(mm.m, 1) - sum(df))) # res dof
    assign = mm.assign
    f = model.mf.f.rhs
    if type == 1
        exclude = Set(assign)
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        @inbounds for id in 1:last(assign)
            delete!(exclude, id)
            ss[id] = SS(model, exclude, pivot)
        end
        ss = _diff(ss)
    elseif type == 2
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        ss[1] = SS(model, Set(assign[2:end]), pivot)
        @inbounds for id in 2:last(assign)
            delcoef = selectcoef(f, id)
            ss[id] = SS(model, delcoef, pivot) - SS(model, delete!(delcoef, id), pivot)
        end
        ss[end] = SS(model, 0, pivot)
    else
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        sse = SS(model, 0, pivot)
        fill!(ss, -sse)
        @inbounds for id in 1:last(assign)
            ss[id] += SS(model, id, pivot)
        end
        ss[end] = sse
    end
    SS(model, 0, pivot) # ensure model unchanged
    first(assign) == 1 || (popfirst!(df); popfirst!(ss))
    ss = tuple(ss...)
    MSR = ss ./ df
    fstat = (MSR[1:(end-1)] / last(MSR)..., NaN)
    pvalue = (ccdf.(FDist.(df, last(df)), abs.(fstat))[1:(end-1)]..., NaN)
    AnovaResult(model, FixedAnovaStatsF{LinearModel, length(df)}(type, size(mm.m, 1), tuple(df...), ss, fstat, pvalue))
end

# --------------------------------------------------------------------------------------------       
# Calculate SS
# use MMatrix/SizedMatrix ?
function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}, exclude::Int, pivot::Bool)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, assign.!= exclude)
    p.delbeta = repeat([0], size(X, 2))
    # cholesky for sparse array in delbeta!
    isdensechol(p) && begin
        F = X'X
        p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    end
    delbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 3

function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}, exclude::Set{Int}, pivot::Bool)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, map(x->!in(x, exclude), assign))
    p.delbeta = repeat([0], size(X, 2))
    # cholesky for sparse array in delbeta!
    isdensechol(p) && begin
        F = X'X
        p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    end 
    delbeta!(p, X, model.model.rr.y) # use delbeta to skip beta0
    updateμ!(model.model.rr, linpred(p, X))
end # for type 1 and 2

isdensechol(p::DensePredChol) = true
isdensechol(p::SparsePredChol) = false

function delbeta!(p::DensePredChol{T, <: Cholesky}, X::SubArray, r::Vector{T}) where T <: BlasReal
    ldiv!(p.chol, mul!(p.delbeta, transpose(X), r))
    p
end
# β = (X'X)⁻¹X'y

function delbeta!(p::DensePredChol{T, <: CholeskyPivoted}, X::SubArray, r::Vector{T}) where T <: BlasReal
    ch = p.chol
    beta = mul!(p.delbeta, adjoint(X), r)
    rnk = rank(ch)
    if rnk == length(beta)
        ldiv!(ch, beta)
    else
        permute!(beta, ch.piv)
        for k = (rnk + 1):length(beta)
            beta[k] = -zero(T)
        end
        LAPACK.potrs!(ch.uplo, view(ch.factors, 1:rnk, 1:rnk), view(beta, 1:rnk))
        invpermute!(beta, ch.piv)
    end
    p
end

# weighted least squares
function delbeta!(p::DensePredChol{T, <: Cholesky}, X::SubArray, r::Vector{T}, wt::Vector{T}) where T <: BlasReal
    scr = mul!(p.scratchm1, Diagonal(wt), X)
    cholesky!(Hermitian(mul!(cholfactors(p.chol), transpose(scr), X), :U))
    mul!(p.delbeta, transpose(scr), r)
    ldiv!(p.chol, p.delbeta)
    p
end

function delbeta!(p::DensePredChol{T, <: CholeskyPivoted}, X::SubArray, r::Vector{T}, wt::Vector{T}) where T <: BlasReal
    cf = cholfactors(p.chol)
    piv = p.chol.piv
    cf .= mul!(p.scratchm2, adjoint(LinearAlgebra.mul!(p.scratchm1, Diagonal(wt), X)), X)[piv, piv]
    cholesky!(Hermitian(cf, Symbol(p.chol.uplo)))
    ldiv!(p.chol, mul!(p.delbeta, transpose(p.scratchm1), r))
    p
end

function delbeta!(p::SparsePredChol{T}, X::SubArray, r::Vector{T}, wt::Vector{T}) where T
    scr = mul!(p.scratch, Diagonal(wt), X)
    XtWX = X'*scr
    c = p.chol = cholesky(Symmetric{eltype(XtWX),typeof(XtWX)}(XtWX, 'L'))
    p.delbeta = c \ mul!(p.delbeta, adjoint(scr), r)
end

linpred(p::LinPred, X::SubArray) = linpred!(Vector{eltype(p.X)}(undef, size(p.X, 1)), p, X)

function linpred!(out, p::LinPred, X::SubArray)
    mul!(out, X, p.delbeta)
end

_diff(v::Vector{T}) where T = cat(v[1], -diff(v), dims = 1)

# ----------------------------------------------------------------------------------------------------
# Linear mixed-effect models

function anova1(::Type{FTest}, model::LinearMixedModel; type::Int = 1, 
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

"""
    calcdof(model::LinearMixedModel)

Calculate degree of freedom of factors and residuals for linear mixed effect models
DOF of residuals are estimated by between-within method:
    dofᵢ = nobsᵢ - dofᵢ₋₁ - nfixᵢ
"""
function calcdof(model::LinearMixedModel)
    randoms = collect(formula(model).rhs)
    fixs = popfirst!(randoms)
    nfac = length(fixs.terms) # number of factors
    isempty(randoms) && (return repeat([nobs(model) - nfac], nfac))
    remat = reverse(model.allterms[1:end - 2])
    femat = model.allterms[end - 1]

    # Determine affected fix effects for each random effects
    affectfixef = Dict{String, Set{String}}()
    for pair in randoms
        lhs = coefnames(pair.lhs)
        isa(lhs, String) ? (lhs == "(Intercept)" ? (continue) : (lhs = [lhs])) : filter!(x->x != "(Intercept)", lhs)
        rhs = coefnames(pair.rhs, Val(:anova))
        haskey(affectfixef, rhs) ? union!(affectfixef[rhs], Set(lhs)) : (affectfixef[rhs] = Set(lhs))
    end 

    fixef = femat.cnames
    nfix = length(fixef) # number of fix effects

    # Determines if fix effects vary at each random effects group
    within = zeros(Bool, nfix, size(remat, 1) + 2)
    for (ranid, ranef) in enumerate(remat)
        levels = unique(ranef.refs)
        within[:, ranid + 1] .= map(1:nfix) do fixid
            any(map(levels) do level
                vals = view(femat.wtx, findall(==(level), ranef.refs), fixid)
                val = first(vals)
                for i in vals
                    i != val && (return true)
                end
                false
            end)
        end
    end 
    within[:, 1] = repeat([true], nfix) # population level
    within[1, 1] = false
    within[:, end] = repeat([false], nfix) # obsevation level

    # Turn affected fix effects into false
    ranef = ["", coefnames.(getproperty.(remat, :trm), Val(:anova))..., ""]
    for (key, value) in affectfixef
        within[findfirst(x->x in value, fixef), findfirst(==(key), ranef)] = false
    end 
    # Find first false for each fix effect, aka level
    level = mapslices(within, dims = 2) do fixef
        findfirst(!, fixef)
    end 
    # Number of fix effects per level
    nfixperlv = zeros(Int, length(ranef))
    for i in level
        nfixperlv[i] += 1
    end 

    # dfᵢ = nobsᵢ - dfᵢ₋₁ - nfixᵢ
    n0 = 1
    n = nobs(model)
    nobsperlv = (n0, length.(getproperty.(remat, :levels))..., n) # number of random effects observation per level
    ndiff = _diff(nobsperlv)
    dfperlv = ndiff .- nfixperlv[2:end]
    dfperlv = (last(dfperlv), dfperlv...)

    # Assign df to each factors based on the level
    assign = asgn(fixs)
    offset = 0
    intercept = first(assign) == 1
    intercept || (offset = 1)
    dfr = ntuple(last(assign) - offset) do i
        dfperlv[level[findlast(x->x == i + offset, assign)]] # use the last one to avoid wrong assign of the first arguments
    end
    df = dof(assign)
    intercept || popfirst!(df)
    tuple(df...), dfr, Dict(:level => tuple(level...), :lvdof => dfperlv)
end

# ----------------------------------------------------------------------------------------
# ANOVA for genaralized linear models
# λ = -2ln(𝓛(̂θ₀)/𝓛(θ)) ~ χ²ₙ , n = difference of predictors

function anova1(::Type{FTest}, model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; kwargs...)
    null = model.mf.f.rhs.terms[1] == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova1(FTest, models)
end

function anova1(::Type{FTest}, models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}}) where N
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
    if first(models).mf.f.rhs.terms[1] == InterceptTerm{false}()
        AnovaResult(models, FixedAnovaStatsF{GeneralizedLinearModel, length(Δdf)}(1, n, Δdf, Δdev , fstat, pval))
    else
        AnovaResult(models, FixedAnovaStatsF{GeneralizedLinearModel, 1 + length(Δdf)}(1, n, (1, Δdf...), (NaN, Δdev...) , (NaN, fstat...), (NaN, pval...)))
    end
end

# generate nested models
"""
    nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = true, kwargs...)
    
Generate nested models from a saturated model. \n
The null model will be a model with at least one factor (including inyercept) if the link function does not allow factors to be 0 (factors in denominators). \n
* `InverseLink` for `Gamma`
* `InverseSquareLink` for `InverseGaussian`
Otherwise, it will be a model with no factors.
"""
function nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = true, kwargs...)
    f = model.mf.f
    # fit models
    l = typeof(model.model.rr).parameters[3]
    l <: InverseLink && (null = true)
    l = l()
    d = model.model.rr.d
    wts = model.model.rr.wts
    offset = model.model.rr.offset
    range = null ? (1:length(f.rhs.terms) - 1) : (0:length(f.rhs.terms) - 1)
    models = map(range) do id
        # create sub-formula, modify schema, create mf and mm
        subf = subformula(model.mf.f.lhs, model.mf.f.rhs, id)
        terms = setdiff(getterms(f.rhs), getterms(subf.rhs))
        schema = deepcopy(model.mf.schema)
        @inbounds for term in terms
            pop!(schema.schema, Term(term))
        end
        pair = collect(pairs(model.mf.data))
        filter!(x->!in(x.first, terms), pair)
        mf = ModelFrame(subf, schema, (; pair...), model.mf.model) 
        mm = ModelMatrix(mf)
        y = response(mf)
        TableRegressionModel(fit(model.mf.model, mm.m, y, d, l; wts = wts, offset = offset, kwargs...), mf, mm)
    end
    (models..., model)
end

subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int) = 
    id > 0 ? FormulaTerm(lhs, collect_matrix_terms(rhs.terms[1:id])) : FormulaTerm(lhs, collect_matrix_terms((InterceptTerm{false}(),)))


# ==================================================================================================================
# ANOVA by Likehood-ratio test 

"""
    anova(::Type{LRT}, models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N};
        testnested::Bool = true,
        pivot::Bool = false)

    anova(::Type{LRT}, models::Vararg{<: LinearMixedModel, N}; 
        testnested::Bool = true,
        adjust_sigma::Bool = true)
    
    anova(::Type{LRT}, models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested ::Bool = true,
        kwargs...)

Analysis of Variance bt F test.

* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
"""

anova(::Type{LRT}, models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N};
        testnested::Bool = true,
        pivot::Bool = false) where N = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; type = type, pivot = pivot)

anova(::Type{LRT}, models::Vararg{<: LinearMixedModel, N}; 
        testnested::Bool = true,
        adjust_sigma::Bool = true) where N = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; type = 1, adjust_sigma = adjust_sigma)

    
anova(::Type{LRT}, models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested ::Bool = true,
        kwargs...) where N = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; kwargs...)

# ------------------------------------------------------------------------------------------------------------
# ANOVA for GeneralizedLinearModel

function anova1(::Type{LRT}, model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; kwargs...)
    null = model.mf.f.rhs.terms[1] == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova1(LRT, models)
end

function anova1(::Type{LRT}, models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}}) where N
    n = Int(nobs(first(models)))
    df = Int.(dof.(models))
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    Δdev = _diffn(dev)
    σ² = dispersion(last(models).model) ^ 2
    lrstat = Δdev ./ σ²
    pval = ccdf.(Chisq.(Δdf), abs.(lrstat))
    AnovaResult(models, FixedAnovaStatsLRT{GeneralizedLinearModel, length(Δdf)}(1, n, Δdf, Δdev, lrstat, pval))
end


# --------------------------------------------------------------------------------------------------
# Nested models 

function anovaN(::Type{FTest}, models::Vararg{TableRegressionModel{<: LinearModel, <: AbstractArray}, N}; 
        testnested::Bool = true) where N
    
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(last(models).model) ^ 2
    fstat = (NaN, msr./σ²...)
    pval = (NaN, ccdf.(FDist.(Δdf, dfr[2:end]), abs.(fstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsF{length(df)}(n, df, dev, fstat, pval))
end

function anovaN(::Type{FTest}, models::Vararg{<: LinearMixedModel, N}; 
        testnested::Bool = true) where N
end

function anovaN(::Type{FTest}, models::Vararg{TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}, N}; 
        testnested::Bool = true) where N

    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(last(models).model) ^ 2
    fstat = (NaN, msr./σ²...)
    pval = (NaN, ccdf.(FDist.(Δdf, dfr[2:end]), abs.(fstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsF{length(df)}(n, df, dev, fstat, pval))
end

function anovaN(::Type{LikelihoodRatioTest}, models::Vararg{TableRegressionModel, N}; testnested::Bool = true) where N
    # AIC and BIC
    n = Int(nobs(first(models)))
    df = dof.(models)
    Δdf = _diff(df)
    σ² = dispersion(last(models).model) ^ 2
    dev = deviance.(models)
    Δdev = _diffn(dev)
    lrstat = (NaN, Δdev ./ σ² ...)
    pval = (NaN, ccdf.(Chisq.(Δdf), abs.(lrstat[2:end]))...)
    AnovaResult(models, NestedAnovaStatsLRT{length(df)}(n, df, dev, lrstat, pval))
end

function anovaN(::Type{LikelihoodRatioTest}, models::Vararg{<: LinearMixedModel, N}; testnested::Bool = true) where N
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
    anova1(test, model; pivot = pivot, type = type)
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

anova_lme(test::Type{T}, f::FormulaTerm, tbl; kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, LinearMixedModel, f, tbl; kwargs...)

function anova(test::Type{T}, ::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
        type::Int = 1, 
        adjust_sigma::Bool = true,
        wts = [], 
        contrasts = Dict{Symbol,Any}(), 
        verbose::Bool = false, 
        REML::Bool = true) where {T <: GoodnessOfFit}
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    anova1(test, model; type = type, adjust_sigma = adjust_sigma)
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
anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); 
        test::Type{T} = canonicalgoodnessoffit(d), 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

anova_glm(test::Type{T}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d); 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, GeneralizedLinearModel, X, y, d, l; kwargs...)

function anova(test::Type{T}, ::Type{GeneralizedLinearModel}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d);
            kwargs...) where {T <: GoodnessOfFit}
    model = glm(X, y, d, l; kwargs...)
    null = model.mf.f.rhs.terms[1] == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova1(test, models)
end    
