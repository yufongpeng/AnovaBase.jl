# ===========================================================================================
# Main algorithm 

@doc """
    anova(models::StatisticalModel...)
    
    anova(models::StatisticalModel...; test::Type{GoodnessOfFit} = GoodnessOfFit, testnested::Bool = true...)
    
    anova(::Type{GoodnessOfFit}, models::StatisticalModel...)

Analysis of variance of nested models.   \n

* models: objects fit by `GLM.lm`, `GLM.glm`, `lme` (`LinearMixedModel`). They should be nested and the last one is the most saturated.
* test: test statistics for goodness of fit. Available tests are `LikelihoodRatioTest` (`LRT`) and `FTest`.

See `anova_lm` for `LinearModel`, `anova_lme` for `LinearMixedModel`, `anova_glm` for `GeneralizedLinearModel`.
""" anova
#---------------------------------------------------------------------------------------------
# ANOVA for simple linear regression

"""
    anova_lm(X, y, allowrankdeficient::Bool = false; <keyword arguments>)

    anova(::Type{LinearModel}, X, y, 
        allowrankdeficient::Bool = false; type::Int = 1, wts = similar(y, 0))

    anova(model::TableRegressionModel{M, T}, 
        allowrankdeficient::Bool = false; type::Int = 1) where {M <: LinearModel, T <: AbstractArray}

ANOVA for simple linear regression.

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`.

The keyword argument `type` specifies type of anova.

`anova_lm` generate a `TableRegressionModel` object through calling `anova`, which is fitted by `lm`.
"""
anova_lm(X, y, allowrankdeficient::Bool = false; kwargs...) = 
        anova(LinearModel, X, y, allowrankdeficient; kwargs...)

function anova(::Type{LinearModel}, X, y, 
               allowrankdeficient::Bool = false; type::Int = 1, kwargs...)
    model = lm(X, y, allowrankdeficient; kwargs...)
    anova(model, allowrankdeficient; type = type)
end

function anova(model::TableRegressionModel{M,T}, 
               allowrankdeficient::Bool = false; type::Int = 1) where {M <: LinearModel, T <: AbstractArray}
    @assert (type in [1,2,3]) "Invalid type"
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    f = model.mf.f.rhs
    if type == 1
        exclude = Set(assign[1:end])
        ss = map(1:(assign[end])) do id
            delete!(exclude, id)
            SS(model, exclude, allowrankdeficient)
        end
        push!(ss,0)
        ss = _diff(ss)
    elseif type == 2
        sse = SS(model, 0, allowrankdeficient)
        ss = map(1:assign[end]) do id
            ifelse(id  == 1, SS(model, Set(assign[2:end]), allowrankdeficient),
            SS(model, selectcoef(f, id), allowrankdeficient) - SS(model, delete!(selectcoef(f, id), id), allowrankdeficient))
        end
        push!(ss,sse)
    else
        sse = SS(model, 0, allowrankdeficient)
        ss = map(1:assign[end]) do id
            SS(model, id, allowrankdeficient) - sse
        end
        push!(ss, sse)
    end
    SS(model, 0, allowrankdeficient) # ensure model unchanged
    width(f.terms[1]) == 0 && (popfirst!(df); popfirst!(ss))
    MSR = ss ./ df
    fstat = [MSR[1:(end-1)] ./ MSR[end]..., NaN]
    pvalue = [ccdf.(FDist.(df, df[end]), abs.(fstat))[1:(end-1)]..., NaN]
    AnovaResult(model, AnovaStats(type,size(mm.m, 1), df, ss, fstat, pvalue))
end

# --------------------------------------------------------------------------------------------       
# Calculate SS
function SS(model::TableRegressionModel{M,T}, exclude::Int, pivot::Bool) where {M <: LinearModel, T <: AbstractArray}
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, assign.!= exclude)
    p.delbeta = repeat([0], size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    delbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 3

function SS(model::TableRegressionModel{M,T}, exclude::Set{Int}, pivot::Bool) where {M <: LinearModel, T <: AbstractArray}
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
# ANOVA for linear mixed-effect models
"""
    anova_lme(X, y; <keyword arguments>)

    anova(::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
            type::Int = 1,
            between::Union{Nothing,Array{Int64,1}} = nothing,
            adjust_sigma::Bool = true,    
            <keyword arguments>)

    anova(model::LinearMixedModel; type::Int = 1, 
            between::Union{Nothing,Array{Int64,1}} = nothing, 
            adjust_sigma::Bool = true)
    
ANOVA for linear mixed-effect models.

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`.

* `type` specifies type of anova.
* `between` specifies the variable that manually assigned to between-subjects. 
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fitted by REML.

Other keyword arguments
* `wts = []`
* `contrasts = Dict{Symbol,Any}()`
* `verbose::Bool = false`
* `REML::Bool = true`

`anova_lme` generate a `LinearMixedModel` object through calling `anova`, which is fitted by `lme` with REML.
"""
anova_lme(X, y; kwargs...) = 
        anova(LinearMixedModel, X, y; kwargs...)


function anova(::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
                type::Int = 1, 
                between::Union{Nothing,Array{Int64,1}} = nothing, 
                adjust_sigma::Bool = true,
                wts = [], 
                contrasts = Dict{Symbol,Any}(), 
                verbose::Bool = false, 
                REML::Bool = true)
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    anova(model, type = type, between = between, adjust_sigma = adjust_sigma)
end

function anova(model::LinearMixedModel; type::Int = 1, 
                between::Union{Nothing,Array{Int64,1}} = nothing, 
                adjust_sigma::Bool = true)

    @assert (type in [1,2,3]) "Invalid type"
    @assert (type in [1,3]) "Type 2 anova is not supported now"
    fet = model.formula.rhs[1]
    ret = model.formula.rhs[2:end]
    femat = model.feterms[1]
    remat = model.reterms
    @assert (length(ret) == 1) "Multiple random factor design is not implemented now"
    
    varβ = vcov(model) 
    β = fixef(model)

    assign = asgn(fet)
    # Determine between/within
    btw = isbtw(fet, assign, remat[1], model.X) # to be modify for multiple random effects 
    isnothing(between) || (btw[between] .= true)
    intercept = width(fet.terms[1]) == 1
    ngroups = map(x->size(x, 2), remat)
    nbetween = Int(prod(nlevels.(fet.terms[btw])))
    n = ngroups[1] / nbetween 
    btw = intercept ? (btw) : (btw[2:end])

    last = assign[end] - assign[1] + 1
    fstat = zeros(Float64, dof(model) - length(β) + last) 
    ss = copy(fstat)
    df = zeros(Int64, dof(model) - length(β) + last)

    fstat[(last+1):end] .= NaN
    @inbounds for id in last+1:length(fstat)-1
        push!(btw, 1)
    end
    push!(btw, 0)
    df[last + 1] = nbetween * (n-1) # to be modify for multiple random effects 
    df[end] = nobs(model) - sum(df) - length(β)
    ss[last + 1] = sum(residuals(model).^2) # to be modify for multiple random effects 
    ss[end] = varest(model) * df[end]

    df[1:last] .= intercept ? dof(assign) : dof(assign)[2:end]
    if type == 1
        invvarfixchol = cholesky(cholesky(varβ) \ Matrix(I, size(varβ)...) |> Hermitian).L # column factor contains between factor should × -1
        model.optsum.REML || adjust_sigma && (invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - size(invvarfixchol, 1))))
        fs = invvarfixchol'β
        fstat[1:last] .= map(((loc, factor), )->sum(fs[assign .== factor] .^ 2) / df[loc],
                            enumerate(unique(assign)))
    else 
        fstat[1:last] = map(unique(assign)) do factor
            select = assign .== factor
            id1 = collect(1:size(assign, 1))[select]
            id2 = collect(1:size(assign, 1))[(!).(select)]
            id = cat(id2, id1, dims = 1)
            varfix = varβ[id, id]
            invvarfixchol = cholesky(cholesky(varfix) \ Matrix(I, size(varfix)...) |> Hermitian).L # column factor contains between factor should × -1
            model.optsum.REML || adjust_sigma && (invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - size(invvarfixchol, 1))))
            fs = invvarfixchol'*(β[id])
            sum(fs[end - length(id1) + 1:end].^2) / length(id1)
        end
    end
    ss[1:last] .= map(1:last) do id
        btw[id] ? fstat[id] * ss[last + 1] * df[id] / df[last+1] : fstat[id] * varest(model) * df[id]
    end
    pvalue = map(1:lastindex(fstat)) do id
        if id > last
            NaN
        elseif btw[id]
            ccdf(FDist(df[id], df[last + 1]), abs(fstat[id]))
        else
            ccdf(FDist(df[id], df[end]), abs(fstat[id]))
        end
    end
    AnovaResult(model, AnovaStatsGrouped(type, nobs(model), ngroups, Bool.(btw), df, ss, fstat, pvalue))
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

"""
    anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

    anova(::Type{M}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d);
        test::GoodnessOfFit = LikelihoodRatioTest(), <keyword arguments>)
    
    anova(model::AbstractGLM; test::GoodnessOfFit = LikelihoodRatioTest())

ANOVA for genaralized linear models.

"""
anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); kwargs...) = 
        anova(GeneralizedLinearModel, X, y, d, l; kwargs...)


function anova(::Type{M}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d);
        test::Type{GoodnessOfFit} = LRT,
        dofit::Bool = true,
        wts::AbstractVector{<:Real} = similar(y, 0),
        offset::AbstractVector{<:Real} = similar(y, 0),
        fitargs...) where {M <: AbstractGLM}
        model = glm(X, y, d, l, dofit = dofit, wts = wts, offset = offset, fitargs...)
        anova(model; test = test, fitargs...)
end

    
function anova(model::AbstractGLM; test::Type{GoodnessOfFit} = LRT, fitargs...)
    f = model.mf.f
    # fit models
    models = map(1:length(f.rhs)-1) do id
        mf = deepcopy(model.mf)
        # modify mf.f.rhs , mf.schema key ,  mf.data key
        mm = ModelMatrix(mf)
        y = response(mf)
        # d, l ,offset, wts in model.rr
        TableRegressionModel(fit(typeof(model), mm.m, y, d, l; wts = wts, offset = offset, fitargs...), mf, mm)
    end
    anova(models..., model, test = typeof(test), testnested = false)
end

# --------------------------------------------------------------------------------------------------
# ANOVA for nested models

# Auto-determination of test
function anova(models::TableRegressionModel ...; test::Type{GoodnessOfFit} = GoodnessOfFit, testnested::Bool = true)
    testnested && (print("")) # isnested
    (test == GoodnessOfFit) || (return anova(test, models...))
    typeof(models[1].model) <: LinearModel && (return anova(FTest, models...))
    # Bernoulli, Binomial, and Poisson fits: LRT
    # Other fits: FTest
    (typeof(models[1].model.rr.d) <: FixDispDist) ? anova(LRT, models...) : anova(FTest, models...)
end

function anova(models::MixedModel ...; test::Type{GoodnessOfFit} = GoodnessOfFit, testnested::Bool = true)
    testnested && (print("")) # isnested
    (test == GoodnessOfFit) || (return anova(test, models...))
    anova(LRT, models...)
end


function anova(::Type{FTest}, models::StatisticalModel...)
    n = Int(nobs(models[1]))
    df = dof.(models)
    Δdf = _diff(df)
    dfr = Int.(dof_residual.(models))
    dev = deviance.(models)
    msr = _diffn(dev) ./Δdf
    σ² = dispersion(models[end].model)^ 2
    fstat = (NaN, msr./σ²...)
    pval = (NaN, ccdf.(FDist.(abs.(Δdf), dfr[2:end]), abs.(fstat[2:end]))...)
    AnovaResult(models, AnovaStatsF(n, df, dev, fstat, pval))
end

function anova(::Type{LikelihoodRatioTest}, models::TableRegressionModel...)
    # AIC and BIC
    n = Int(nobs(models[1]))
    df = dof.(models)
    Δdf = _diff(df)
    σ² = dispersion(models[end].model)^ 2
    dev = deviance.(models)
    Δdev = _diffn(dev)
    lrstat = (NaN, Δdev ./ σ² ...)
    pval = (NaN, ccdf.(Chisq.(abs.(Δdf)), abs.(lrstat[2:end]))...)
    AnovaResult(models, AnovaStatsLRT(n, df, dev, lrstat, pval))
end

function anova(::Type{LikelihoodRatioTest}, models::MixedModel...)
    # AIC and BIC
    n = Int(nobs(models[1]))
    df = dof.(models)
    Δdf = _diff(df)
    dev = deviance.(models)
    Δdev = _diffn(dev)
    lrstat = (NaN, Δdev...)
    pval = (NaN, ccdf.(Chisq.(abs.(Δdf)), abs.(lrstat[2:end]))...)
    AnovaResult(models, AnovaStatsLRT(n, df, dev, lrstat, pval))
end


"""
function anova_QR()
    @assert (type in [1,2,3]) "Invalid type"

    model = TableRegressionModel(LmResp(y, similar(y, 0)),AnovaDensePredQR(mm.m))
    mf = model.mf
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    effect = (transpose(model.pp.qr.Q)*model.rr.y).^2
    whole_id = 1
    var_id = 1
    ss = zeros(Float64,assign[end])
    while whole_id < length(assign)
        whole_id += 1
        var_id == assign[whole_id] || (var_id += 1)
        ss[var_id] += effect[whole_id]
    end
    ss[end] = sum(effect[(whole_id+1):end])
    popfirst!(df)
end
"""