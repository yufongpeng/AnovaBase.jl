# ==========================================================================================================
# Backend funcion

formula(model::TableRegressionModel) = model.mf.f

# Backend for LinearModel
function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; 
    type::Int = 1)

    assign = model.mm.assign
    f = formula(model).rhs
    if type == 1
        exclude = Set(assign)
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 2)
        @inbounds for id in 0:last(assign)
            delete!(exclude, id)
            ss[id + 1] = SS(model, exclude)
        end
        ss = -diff(ss)
    elseif type == 2
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        ss[1] = SS(model, Set(assign)) - SS(model, Set(assign[2:end]))
        @inbounds for id in 2:last(assign)
            delcoef = selectcoef(f, id)
            ss[id] = SS(model, delcoef) - SS(model, delete!(delcoef, id))
        end
        ss[end] = SS(model, 0)
    else
        # Preallocate without push!
        ss = zeros(Float64, last(assign) + 1)
        sse = SS(model, 0)
        fill!(ss, -sse)
        @inbounds for id in 1:last(assign)
            ss[id] += SS(model, id)
        end
        ss[end] = sse
    end
    SS(model, 0) # ensure model unchanged
    first(assign) == 1 || popfirst!(ss)
    tuple(ss...)
end

# use MMatrix/SizedMatrix ?
function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}, exclude::Int)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, assign.!= exclude)
    p.delbeta = repeat([0], size(X, 2))
    # cholesky for sparse array in delbeta!
    isdensechol(p) && begin
        F = Hermitian(float(X'X))
        p.chol = genchol(p.chol, F)
    end
    isempty(model.model.rr.wts) ? delbeta!(p, X, model.model.rr.y) : delbeta!(p, X, model.model.rr.y, model.model.rr.wts)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 3

function SS(model::TableRegressionModel{<: LinearModel, <: AbstractArray}, exclude::Set{Int})
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X, :, map(x->!in(x, exclude), assign))
    p.delbeta = repeat([0], size(X, 2))
    # cholesky for sparse array in delbeta!
    isdensechol(p) && begin
        F = Hermitian(float(X'X))
        p.chol = genchol(p.chol, F)
    end 
    isempty(model.model.rr.wts) ? delbeta!(p, X, model.model.rr.y) : delbeta!(p, X, model.model.rr.y, model.model.rr.wts) # use delbeta to skip beta0
    updateμ!(model.model.rr, linpred(p, X))
end # for type 1 and 2

isdensechol(p::DensePredChol) = true
isdensechol(p::SparsePredChol) = false

genchol(chol::CholeskyPivoted{<: Number, <: AbstractMatrix{<: Number}}, F::AbstractMatrix{<: Number}) = 
    cholesky!(F, Val(true), tol = -one(eltype(F)), check = false)

genchol(chol::Cholesky{<: Number, <: AbstractMatrix{<: Number}}, F::AbstractMatrix{<: Number}) = 
    cholesky!(F)

# Backend for GeneralizedLinearModel
#= 
function _fit!(m::AbstractGLM, verbose::Bool, maxiter::Integer, minstepfac::Real,
               atol::Real, rtol::Real, start)

    # Return early if model has the fit flag set
    m.fit && return m

    # Check arguments
    maxiter >= 1       || throw(ArgumentError("maxiter must be positive"))
    0 < minstepfac < 1 || throw(ArgumentError("minstepfac must be in (0, 1)"))

    # Extract fields and set convergence flag
    cvg, p, r = false, m.pp, m.rr
    lp = r.mu

    # Initialize β, μ, and compute deviance
    if start == nothing || isempty(start)
        # Compute beta update based on default response value
        # if no starting values have been passed
        delbeta!(p, wrkresp(r), r.wrkwt)
        linpred!(lp, p)
        updateμ!(r, lp)
        installbeta!(p)
    else
        # otherwise copy starting values for β
        copy!(p.beta0, start)
        fill!(p.delbeta, 0)
        linpred!(lp, p, 0)
        updateμ!(r, lp)
    end
    devold = deviance(m)

    for i = 1:maxiter
        f = 1.0 # line search factor
        local dev

        # Compute the change to β, update μ and compute deviance
        try
            delbeta!(p, r.wrkresid, r.wrkwt)
            linpred!(lp, p)
            updateμ!(r, lp)
            dev = deviance(m)
        catch e
            isa(e, DomainError) ? (dev = Inf) : rethrow(e)
        end

        # Line search
        ## If the deviance isn't declining then half the step size
        ## The rtol*dev term is to avoid failure when deviance
        ## is unchanged except for rouding errors.
        while dev > devold + rtol*dev
            f /= 2
            f > minstepfac || error("step-halving failed at beta0 = $(p.beta0)")
            try
                updateμ!(r, linpred(p, f))
                dev = deviance(m)
            catch e
                isa(e, DomainError) ? (dev = Inf) : rethrow(e)
            end
        end
        installbeta!(p, f)

        # Test for convergence
        verbose && println("Iteration: $i, deviance: $dev, diff.dev.:$(devold - dev)")
        if devold - dev < max(rtol*devold, atol)
            cvg = true
            break
        end
        @assert isfinite(dev)
        devold = dev
    end
    cvg || throw(ConvergenceException(maxiter))
    m.fit = true
    m
end
=#

# Linear prediction
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
    elseif size(ch.factors) == (0, 0) # Unstable when the matrix is 0x0
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

# Create nestedmodels
function nestedmodels(model::TableRegressionModel{<: LinearModel, <: AbstractArray}; null::Bool = true, kwargs...)
    f = formula(model)
    null && (isnullable(model.model.pp.chol) || (null = false))
    range = null ? (0:length(f.rhs.terms) - 1) : (1:length(f.rhs.terms) - 1)
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

function nestedmodels(model::TableRegressionModel{<: GeneralizedLinearModel, <: AbstractArray}; null::Bool = true, kwargs...)
    f = formula(model)
    # fit models
    link = typeof(model.model.rr).parameters[3]()
    null && (isnullable(l) || (null = false))
    dist = model.model.rr.d
    wts = model.model.rr.wts
    offset = model.model.rr.offset
    range = null ? (0:length(f.rhs.terms) - 1) : (1:length(f.rhs.terms) - 1)
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
        TableRegressionModel(fit(model.mf.model, mm.m, y, dist, link; wts = wts, offset = offset, kwargs...), mf, mm)
    end
    (models..., model)
end


nestedmodels(::Type{LinearModel}, formula, data; null::Bool = true, kwargs...) = 
    nestedmodels(lm(formula, data; kwargs...); null = null, kwargs...)
nestedmodels(::Type{GeneralizedLinearModel}, formula, data, distr::UnivariateDistribution, link::Link = canonicallink(d); null::Bool = true, kwargs...) = 
    nestedmodels(glm(formula, data, distr, link; kwargs...); null = null, kwargs...)
    
# Null model for CholeskyPivoted is unstable now
isnullable(chol::CholeskyPivoted{<: Number, <: AbstractMatrix{<: Number}}) = false
isnullable(chol::Cholesky{<: Number, <: AbstractMatrix{<: Number}}) = true
isnullable(::T) where {T <: InverseLink} = false
isnullable(::T) where {T <: Link} = true

