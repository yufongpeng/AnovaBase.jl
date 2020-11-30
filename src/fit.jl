# ===========================================================================================
# Main algorithm 

@doc """
    anova(<models>...; test::Type{T}) where {T <: GoodnessOfFit}

Analysis of variance.

* `models`: model objects
    1. `TableRegressionModel{<: LinearModel, T}` fit by `GLM.lm`
    2. `TableRegressionModel{<: GeneralizedLinearModel, T}` fit by `GLM.glm`
    3. `LinearMixedModel` fit by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`
    If mutiple models are provided, they should be nested and the last one is the most saturated.
* `test`: test statistics for goodness of fit. Available tests are `LikelihoodRatioTest` (`LRT`) and `FTest`. \n
    If no test argument is provided, the function will automatically determine based on the model type:
    1. `TableRegressionModel{<: LinearModel, T}`: `FTest`.
    2. `TableRegressionModel{<: GeneralizedLinearModel, T}`: based on distribution function, see `canonicalgoodnessoffit`.
    3. `LinearMixedModel`: `FTest` for one model, `LRT` for nested models.

For fitting new models and conducting anova at the same time,  
see `anova_lm` for `LinearModel`, `anova_lme` for `LinearMixedModel`, `anova_glm` for `GeneralizedLinearModel`.
""" anova

anova(models::TableRegressionModel{<: LinearModel, S}...; 
        test::Type{T} = FTest,
        kwargs...) where {S, T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

anova(models::LinearMixedModel...; 
        test::Type{T} = length(models) > 1 ? (LRT) : (FTest), 
        kwargs...) where {T <: GoodnessOfFit} = 
    anova(test, models...; kwargs...)

anova(models::TableRegressionModel{<: GeneralizedLinearModel, S}...; 
        test::Type{T} = canonicalgoodnessoffit(models[1].model.rr.d),
        kwargs...) where {S, T <: GoodnessOfFit} = 
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
    anova(::Type{FTest}, models::TableRegressionModel{<: LinearModel, T}...;
        testnested::Bool = true,
        type::Int = 1, 
        pivot::Bool = false)

    anova(::Type{FTest}, models::T...; 
        testnested::Bool = true,
        type::Int = 1, 
        adjust_sigma::Bool = true, 
        between::Union{Nothing,Array{Int64,1}} = nothing) where {T <: LinearMixedModel}
    
    anova(::Type{FTest}, models::TableRegressionModel{<: GeneralizedLinearModel, T}...; 
        testnested ::Bool = true,
        kwargs...)

Analysis of Variance bt F test.

* `type` specifies type of anova. For one `LinearModel` `1, 2, 3` are valid; for one `LinearMixedModel` `1, 3` are valid. For others, only `1` is valid.
* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
* `between` specifies the variable that manually assigned to between-subjects. 
"""
anova(::Type{FTest}, models::TableRegressionModel{<: LinearModel, T}...;
        testnested::Bool = true,
        type::Int = 1, 
        pivot::Bool = false) where T = 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; type = type, pivot = pivot)

anova(::Type{FTest}, models::T...; 
        testnested::Bool = true,
        type::Int = 1, 
        adjust_sigma::Bool = true, 
        between::Union{Nothing,Array{Int64,1}} = nothing) where {T <: LinearMixedModel}= 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; type = type, adjust_sigma = adjust_sigma, between = between)

    
anova(::Type{FTest}, models::TableRegressionModel{<: GeneralizedLinearModel, T}...; 
        testnested ::Bool = true,
        kwargs...) where T = 
    length(models) > 1 ? anovaN(FTest, models...; testnested = testnested) : anova1(FTest, models...; kwargs...)

# ------------------------------------------------------------------------------------------------------------------
# Simple linear regression

function anova1(::Type{FTest}, model::TableRegressionModel{<: LinearModel, T}; 
                type::Int = 1, 
                pivot::Bool = false) where T 
    @assert (type in [1,2,3]) "Invalid type"
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df, Int(size(mm.m, 1) - sum(df)))
    assign = mm.assign
    f = model.mf.f.rhs
    if type == 1
        exclude = Set(assign)
        ss = map(1:last(assign)) do id
            delete!(exclude, id)
            SS(model, exclude, pivot)
        end
        push!(ss,0)
        ss = _diff(ss)
    elseif type == 2
        sse = SS(model, 0, pivot)
        ss = map(1:last(assign)) do id
            ifelse(id  == 1, SS(model, Set(assign[2:end]), pivot),
            SS(model, selectcoef(f, id), pivot) - SS(model, delete!(selectcoef(f, id), id), pivot))
        end
        push!(ss,sse)
    else
        sse = SS(model, 0, pivot)
        ss = map(1:last(assign)) do id
            SS(model, id, pivot) - sse
        end
        push!(ss, sse)
    end
    ss = tuple(ss...)
    SS(model, 0, pivot) # ensure model unchanged
    width(first(f.terms)) == 0 && (popfirst!(df); popfirst!(ss))
    MSR = ss ./ df
    fstat = (MSR[1:(end-1)] / last(MSR)..., NaN)
    pvalue = (ccdf.(FDist.(df, last(df)), abs.(fstat))[1:(end-1)]..., NaN)
    AnovaResult(model, FixedAnovaStatsF{LinearModel, length(df)}(type, size(mm.m, 1), tuple(df...), ss, fstat, pvalue))
end

# --------------------------------------------------------------------------------------------       
# Calculate SS
# use MMatrix/SizedMatrix ?
function SS(model::TableRegressionModel{<: LinearModel,T}, exclude::Int, pivot::Bool) where T
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, assign.!= exclude)
    p.delbeta = repeat([0], size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    delbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 3

function SS(model::TableRegressionModel{<: LinearModel,T}, exclude::Set{Int}, pivot::Bool) where T
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, map(x->!in(x, exclude), assign))
    p.delbeta = repeat([0], size(X, 2))
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
        for k=(rnk+1):length(beta)
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
            adjust_sigma::Bool = true, 
            between::Union{Nothing,Array{Int64,1}} = nothing)

    @assert (type in [1,2,3]) "Invalid type"
    @assert (type in [1,3]) "Type 2 anova is not supported now"
    fet = first(model.formula.rhs)
    ret = model.formula.rhs[2:end]
    femat = first(model.feterms)
    remat = model.reterms
    @assert (length(ret) == 1) "Multiple random factor design is not implemented now"

    varβ = vcov(model) 
    β = fixef(model)

    assign = asgn(fet)
    # Determine between/within
    btw = isbtw(fet, assign, first(remat), model.X) # to be modify for multiple random effects 
    isnothing(between) || (btw[between] .= true)
    intercept = width(first(fet.terms)) == 1
    ngroups = map(x->size(x, 2), remat)
    nbetween = Int(prod(nlevels.(fet.terms[btw])))
    n = first(ngroups) / nbetween 
    btw = intercept ? (btw) : (btw[2:end])

    lst = last(assign) - first(assign) + 1
    df = zeros(Int64, dof(model) - length(β) + lst)

    df[lst + 1] = nbetween * (n-1) # to be modify for multiple random effects 
    df[end] = nobs(model) - sum(df) - length(β)

    df[1:lst] .= intercept ? dof(assign) : dof(assign)[2:end]
    # use MMatrix/SizedMatrix ?
    if type == 1
        invvarfixchol = cholesky(cholesky(varβ) \ Matrix(I, size(varβ)...) |> Hermitian).L # column factor contains between factor should × -1
        model.optsum.REML || adjust_sigma && (invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - length(β)))) # little offset, overestimated
        fs = invvarfixchol'β
        uniqas = unique(assign)
        fstat = ntuple(lastindex(uniqas)) do loc
                sum(fs[assign .== uniqas[loc]] .^ 2) / df[loc]
        end
    else 
        # calculate block by block
        adjust = 1.0
        model.optsum.REML || adjust_sigma && (adjust = (nobs(model) - length(β)) / nobs(model)) # little offset, overestimated
        if first(assign) > 1
            fstat = ntuple(last(assign) - 1) do factor
                select = assign .== (factor + 1)
                varfix = varβ[select, select]
                invvarfix = inv(varfix) 
                β[select]' * invvarfix * β[select] / length(select) * adjust
            end
        else
            fstat = ntuple(last(assign)) do factor
                select = assign .== factor
                varfix = varβ[select, select]
                invvarfix = inv(varfix) 
                β[select]' * invvarfix * β[select] / length(select)
            end
        end
    end

    pvalue = ntuple(lastindex(fstat)) do id
        if btw[id]
            ccdf(FDist(df[id], df[lst + 1]), abs(fstat[id]))
        else
            ccdf(FDist(df[id], last(df)), abs(fstat[id]))
        end
    end
    AnovaResult(model, MixedAnovaStatsF{LinearMixedModel, length(fstat)}(type, nobs(model), first(ngroups), tuple(Bool.(btw)...), tuple(df[1:lst]...), fstat, pvalue))
end

# Determine between subjects vaiable
function isbtw(fet::MatrixTerm, assign::Array{Int64,1}, remat::ReMat, X::Matrix)
    n = length(fet.terms)
    between = ones(Bool, n)
    select = 1:length(assign)
    for id in 1:size(remat, 2)
        loc = findall(==(1), view(remat, :, id))
        x = view(X, loc, :)
        for level in select
            if length(unique(x[:, level])) > 1
                factor = assign[level]
                select = select[assign[select] .!= factor]
                between[factor] = false
            end
        end
    end
    between[1] = false
    between
end

# ----------------------------------------------------------------------------------------
# ANOVA for genaralized linear models
# λ = -2ln(𝓛(̂θ₀)/𝓛(θ)) ~ χ²ₙ , n = difference of predictors

function anova1(::Type{FTest}, model::TableRegressionModel{<: GeneralizedLinearModel, T}; kwargs...) where T
    null = model.mf.f.rhs.terms[1] == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova1(FTest, models)
end

function anova1(::Type{FTest}, models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, T}}) where {N, T}
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
    nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, T}; null::Bool = true, kwargs...) where T 
    
Generate nested models from a saturated model. \n
The null model will be a model with at least one factor (including inyercept) if the link function does not allow factors to be 0 (factors in denominators). \n
* `InverseLink` for `Gamma`
* `InverseSquareLink` for `InverseGaussian`
Otherwise, it will be a model with no factors.
"""
function nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, T}; null::Bool = true, kwargs...) where T
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
        # modify mf.f.rhs ::Tuple, mf.schema key ::Dict,  mf.data key ::NamedTuple
        subf = subformula(model.mf.f.lhs, model.mf.f.rhs, id)
        terms = setdiff(getterms(f.rhs), getterms(subf.rhs))
        schema = deepcopy(model.mf.schema)
        @inbounds for term in terms
            pop!(schema.schema, Term(term))
        end
        pair = collect(pairs(model.mf.data))
        filter!(x->!(x.first in terms), pair)
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
    anova(::Type{LRT}, models::TableRegressionModel{<: LinearModel, T}...;
        testnested::Bool = true,
        pivot::Bool = false)

    anova(::Type{LRT}, models::T...; 
        testnested::Bool = true,
        adjust_sigma::Bool = true, 
        between::Union{Nothing,Array{Int64,1}} = nothing) where {T <: LinearMixedModel}
    
    anova(::Type{LRT}, models::TableRegressionModel{<: GeneralizedLinearModel, T}...; 
        testnested ::Bool = true,
        kwargs...)

Analysis of Variance bt F test.

* `testnested` checks if models are nested, when multiple models are provided. Not implemented now.
* `pivot` determinea if pivot is used, if modelmatrix is rank deficient, 
* `adjust_sigma` determines if adjusting to REML if `LinearMixedModel` is fit by maximum likelihood. The result is slightly different with that of model fit by REML. This problem is be fixed.
* `between` specifies the variable that manually assigned to between-subjects. 
"""

anova(::Type{LRT}, models::TableRegressionModel{<: LinearModel, T}...;
        testnested::Bool = true,
        pivot::Bool = false) where T = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; type = type, pivot = pivot)

anova(::Type{LRT}, models::T...; 
        testnested::Bool = true,
        adjust_sigma::Bool = true, 
        between::Union{Nothing,Array{Int64,1}} = nothing) where {T <: LinearMixedModel} = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; type = 1, adjust_sigma = adjust_sigma, between = between)

    
anova(::Type{LRT}, models::TableRegressionModel{<: GeneralizedLinearModel, T}...; 
        testnested ::Bool = true,
        kwargs...) where T = 
    length(models) > 1 ? anovaN(LRT, models...; testnested = testnested) : anova1(LRT, models...; kwargs...)

# ------------------------------------------------------------------------------------------------------------
# ANOVA for GeneralizedLinearModel

function anova1(::Type{LRT}, model::TableRegressionModel{<: GeneralizedLinearModel, T}; kwargs...) where T
    null = model.mf.f.rhs.terms[1] == InterceptTerm{false}()
    models = nestedmodels(model; null = null, kwargs...)
    anova1(LRT, models)
end

function anova1(::Type{LRT}, models::NTuple{N, TableRegressionModel{<: GeneralizedLinearModel, T}}) where {N, T}
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

function anovaN(::Type{FTest}, models::TableRegressionModel{<: LinearModel, T}...; 
        testnested::Bool = true) where T
    
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

function anovaN(::Type{FTest}, models::LinearMixedModel...; 
        testnested::Bool = true)
end

function anovaN(::Type{FTest}, models::TableRegressionModel{<: GeneralizedLinearModel, T}...; 
        testnested::Bool = true) where T

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

function anovaN(::Type{LikelihoodRatioTest}, models::TableRegressionModel...; testnested::Bool = true)
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

function anovaN(::Type{LikelihoodRatioTest}, models::MixedModel...; testnested::Bool = true)
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
            between::Union{Nothing,Array{Int64,1}} = nothing, 
            adjust_sigma::Bool = true)

ANOVA for linear mixed-effect models.

The arguments `f` and `tbl` are `Formula` and `DataFrame`.

* `type` specifies type of anova. only `1, 3` are valid.
* `between` specifies the variable that manually assigned to between-subjects. 
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
        between::Union{Nothing, Array{Int64,1}} = nothing,
        wts = [], 
        contrasts = Dict{Symbol,Any}(), 
        verbose::Bool = false, 
        REML::Bool = true) where {T <: GoodnessOfFit}
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    anova1(test, model; type = type, adjust_sigma = adjust_sigma, between = between)
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
