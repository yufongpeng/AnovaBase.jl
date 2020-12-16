# ==========================================================================================================
# funcion utilized by fitting

"""
    canonicalgoodnessoffit(::FixDispDist) = LRT
    canonicalgoodnessoffit(::UnivariateDistribution) = FTest

    const FixDispDist = Union{Bernoulli, Binomial, Poisson}
    
Return LRT if the distribution has fixed dispersion
"""
canonicalgoodnessoffit(::FixDispDist) = LRT
canonicalgoodnessoffit(::UnivariateDistribution) = FTest

# Calculate dof from assign
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
    Int.(dofv)
end

# Calculate SS
function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
    type::Int = 1, 
    pivot::Bool = false)

    assign = model.mm.assign
    f = formula(model).rhs
    if type == 1
        exclude = Set(assign)
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 2)
        @inbounds for id in 0:last(assign)
            delete!(exclude, id)
            ss[id + 1] = SS(model, exclude, pivot)
        end
        ss = -diff(ss)
    elseif type == 2
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        ss[1] = SS(model, Set(assign), pivot) - SS(model, Set(assign[2:end]), pivot)
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
    first(assign) == 1 || popfirst!(ss)
    tuple(ss...)
end

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
    isempty(model.model.rr.wts) ? delbeta!(p, X, model.model.rr.y) : delbeta!(p, X, model.model.rr.y, model.model.rr.wts)
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
    isempty(model.model.rr.wts) ? delbeta!(p, X, model.model.rr.y) : delbeta!(p, X, model.model.rr.y, model.model.rr.wts) # use delbeta to skip beta0
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
    scr = mul!(similar(X), Diagonal(wt), X)
    cholesky!(Hermitian(mul!(cholfactors(p.chol), transpose(scr), X), :U))
    mul!(p.delbeta, transpose(scr), r)
    ldiv!(p.chol, p.delbeta)
    p
end

# experimental
function delbeta!(p::DensePredChol{T, <: CholeskyPivoted}, X::SubArray, r::Vector{T}, wt::Vector{T}) where T <: BlasReal
    cf = cholfactors(p.chol)
    piv = p.chol.piv
    scr = mul!(similar(X), Diagonal(wt), X)
    cf .= mul!(similar(p.chol.factors), adjoint(scr), X)[piv, piv]
    cholesky!(Hermitian(cf, Symbol(p.chol.uplo)))
    ldiv!(p.chol, mul!(p.delbeta, transpose(scr), r))
    p
end

# experimental
function delbeta!(p::SparsePredChol{T}, X::SubArray, r::Vector{T}, wt::Vector{T}) where T
    scr = mul!(similar(X), Diagonal(wt), X)
    XtWX = X'*scr
    c = p.chol = cholesky(Symmetric{eltype(XtWX),typeof(XtWX)}(XtWX, 'L'))
    p.delbeta = c \ mul!(p.delbeta, adjoint(scr), r)
end

linpred(p::LinPred, X::SubArray) = linpred!(Vector{eltype(p.X)}(undef, size(p.X, 1)), p, X)

function linpred!(out, p::LinPred, X::SubArray)
    mul!(out, X, p.delbeta)
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

# generate nested models
"""
    nestedmodels(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; null::Bool = false, kwargs...)
    nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = false, kwargs...)
    nestedmodels(model::LinearMixedModel; null::Bool = false, kwargs...)

Generate nested models from a saturated model. \n
The null model will be a model with at least one factor (including inyercept) if the link function does not allow factors to be 0 (factors in denominators). \n
* `InverseLink` for `Gamma`
* `InverseSquareLink` for `InverseGaussian`
Otherwise, it will be a model with no factors.
"""
function nestedmodels(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; null::Bool = false, kwargs...)
    f = formula(model)
    range = null ? (1:length(f.rhs.terms) - 1) : (0:length(f.rhs.terms) - 1)
    wts = model.model.rr.wts
    models = map(range) do id
        # create sub-formula, modify schema, create mf and mm
        subf = subformula(f.lhs, f.rhs, id)
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
        TableRegressionModel(fit(model.mf.model, mm.m, y; wts = wts, kwargs...), mf, mm)
    end
    (models..., model)
end

function nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = false, kwargs...)
    f = formula(model)
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
        subf = subformula(f.lhs, f.rhs, id)
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

function nestedmodels(model::LinearMixedModel; null::Bool = false, kwargs...)
    f = formula(model)
    range = null ? (1:length(f.rhs[1].terms) - 1) : (0:length(f.rhs[1].terms) - 1)
    assign = asgn(first(f.rhs))
    REML = model.optsum.REML
    pivot = first(model.feterms).piv
    X = first(model.feterms).x
    y = last(model.feterms)
    T = promote_type(Float64, eltype(y.x))
    size(X, 2) > 0 && (X = pivot == collect(1:size(X, 2)) ? X : X[:, pivot])
    reterms = model.reterms
    sqrtwts = model.sqrtwts
    models = map(range) do id
        subf = subformula(f.lhs, f.rhs, id)
        cnames = coefnames(first(subf.rhs))
        # modify X by assign
        select = findall(x->x <= id, assign)
        subX = X[:, select]
        feterms = FeMat{T}[]
        push!(feterms, FeMat(subX, isa(cnames, String) ? [cnames] : collect(cnames)))
        push!(feterms, y)
        allterms = convert(Vector{Union{AbstractReMat{T},FeMat{T}}}, vcat(reterms, feterms))
        reweight!.(allterms[end - 1:end], Ref(sqrtwts))
        A, L = createAL(allterms)
        lbd = foldl(vcat, lowerbd(c) for c in reterms)
        θ = foldl(vcat, getθ(c) for c in reterms)
        subX = first(feterms)
        optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
        fill!(optsum.xtol_abs, 1.0e-10)
        lmm = LinearMixedModel(
            subf,
            allterms,
            reterms,
            feterms,
            sqrtwts,
            model.parmap,
            (n = size(subX, 1), p = subX.rank, nretrms = length(reterms)),
            A,
            L,
            optsum,
            )
        fit!(lmm, REML = REML, verbose = false)
        lmm
    end
    (models..., model)
end

subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int) = 
    id > 0 ? FormulaTerm(lhs, collect_matrix_terms(rhs.terms[1:id])) : FormulaTerm(lhs, collect_matrix_terms((InterceptTerm{false}(),)))

subformula(lhs::AbstractTerm, rhs::NTuple{N, AbstractTerm}, id::Int) where N = 
    id > 0 ? FormulaTerm(lhs, (collect_matrix_terms(first(rhs).terms[1:id]), rhs[2:end]...)) : FormulaTerm(lhs, (collect_matrix_terms((InterceptTerm{false}(),)), rhs[2:end]...))


 
