
# Old version functions
@deprecate function isbtw(fet::MatrixTerm, assign::Array{Int64,1}, remat::ReMat, X::Matrix)
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

# Calculate number of groups
@deprecate nlevels(term::CategoricalTerm) = length(term.contrasts.levels)
@deprecate nlevels(term::ContinuousTerm) = 1 
@deprecate nlevels(term::InterceptTerm) = 1 
@deprecate nlevels(term::InteractionTerm) = prod(nlevels.(term.terms))
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

# -------------------------------------------------------------------------------------------------------------
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
AnovaResult(model, FixedAnovaStats(type,size(mm.m, 1), df, ss, fstat, pvalue))
end

end

# ========================================================================================
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
* `adjust_sigma` determines whether adjust σ to match that of linear mixed-effect model fit by REML.

Other keyword arguments
* `wts = []`
* `contrasts = Dict{Symbol,Any}()`
* `verbose::Bool = false`
* `REML::Bool = true`

`anova_lme` generate a `LinearMixedModel` object through calling `anova`, which is fit by `lme` with REML.
# --------------------------------------------------------------------------------------------------------
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
    fstat = similar(β) 
    df = zeros(Int64, dof(model) - length(β) + last)

    df[last + 1] = nbetween * (n-1) # to be modify for multiple random effects 
    df[end] = nobs(model) - sum(df) - length(β)

    df[1:last] .= intercept ? dof(assign) : dof(assign)[2:end]
    # use MMatrix/SizedMatrix ?
    if type == 1
        invvarfixchol = cholesky(cholesky(varβ) \ Matrix(I, size(varβ)...) |> Hermitian).L # column factor contains between factor should × -1
        model.optsum.REML || adjust_sigma && (invvarfixchol = invvarfixchol / sqrt(nobs(model) / (nobs(model) - size(invvarfixchol, 1))))
        fs = invvarfixchol'β
        fstat[1:last] .= map(((loc, factor), )->sum(fs[assign .== factor] .^ 2) / df[loc],
                            enumerate(unique(assign)))
    else 
        # calculate block by block
        fstat[1:last] = map(unique(assign)) do factor
            select = assign .== factor
            varfix = varβ[select, select]
            invvarfix = inv(varfix) # column factor contains between factor should × -1
            fs = invvarfixchol'*(β[id])
            β[select]' invvarfix * β[select] / length(select)
        end
        model.optsum.REML || adjust_sigma && (fstat[1:last] ./= (nobs(model) / (nobs(model) - length(β))))
    end

    pvalue = map(1:lastindex(fstat)) do id
        if btw[id]
            ccdf(FDist(df[id], df[last + 1]), abs(fstat[id]))
        else
            ccdf(FDist(df[id], df[end]), abs(fstat[id]))
        end
    end
    AnovaResult(model, MixedAnovaStats(type, nobs(model), ngroups, Bool.(btw), df[1:last],  fstat, pvalue))
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

# ======================================================================================================================
    anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); <keyword arguments>)

    anova(::Type{M}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d);
        test::GoodnessOfFit = LikelihoodRatioTest(), <keyword arguments>)

ANOVA for genaralized linear models.

* `d`: a `GLM.UnivariateDistribution`.
* `l`: a `GLM.Link`

For other keyword arguments, see `fit`.
# -----------------------------------------------------------------------------------------------------------------
anova_glm(X, y, d::UnivariateDistribution, l::Link = canonicallink(d); kwargs...) = 
        anova(GeneralizedLinearModel, X, y, d, l; kwargs...)


function anova(::Type{GeneralizedLinearModel}, X, y, d::UnivariateDistribution, l::Link = canonicallink(d);
        test::Type{T} = GoodnessOfFit, kwargs...) where {T <: GoodnessOfFit}
        model = glm(X, y, d, l; kwargs...)
        models = nestedmodels(model; kwargs...)
        anova(models...; test = test, testnested = false)
end

    
function nestedmodels(model::TableRegressionModel{M, T}; kwargs...) where {M <: AbstractGLM, T}
    f = model.mf.f
    # fit models
    l = typeof(model.model.rr).parameters[3]()
    d = model.model.rr.d
    wts = model.model.rr.wts
    offset = model.model.rr.offset
    models = map(1:length(f.rhs.terms) - 1) do id
        # modify mf.f.rhs ::Tuple, mf.schema key ::Dict,  mf.data key ::NamedTuple
        subf = FormulaTerm(model.mf.f.lhs, collect_matrix_terms(model.mf.f.rhs.terms[1:id]))
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
    (models...,model)
end

# =================================================================================================================
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