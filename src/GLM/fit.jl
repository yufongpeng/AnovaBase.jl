# ==========================================================================================================
# Backend funcion

formula(trm::TableRegressionModel) = trm.mf.f

function deviances(trm::TableRegressionModel{<: LinPredModel, <: AbstractArray{T}}; 
                    type::Int = 1, kwargs...) where {T <: BlasReal}
    @assert isa(trm.model.pp, DensePredChol) "Other PredChol types are not implemented"

    assign, f = trm.mm.assign, formula(trm).rhs
    # Determine null model by link function
    ## model with 0 factors is null if allowed
    ## if intercept is not included, 
    ### the first factor is null for type 1
    ### the first and only factor is null for type 2
    # start value
    start = isnullable(trm.model) ? 0 : first(assign)
    if type == 1
        todel = union(start, assign)
        # ~ 0 + A, ~ 1
        start > 0 && @assert length(todel) > 1 "Empty model is not valid"
        exclude = Set(assign)
        devs = zeros(T, length(todel) + 1)
        @inbounds for (id, del) in enumerate(todel)
            delete!(exclude, del)
            devs[id] = deviance(trm, exclude; kwargs...)
        end
        devs = -diff(devs)
    elseif type == 2
        todel = unique(assign)
        if start > 0 
            selectcoef(f, Val(first(todel))) == Set(todel) && popfirst!(todel)
            # ~ 0 + A, ~ 1
            @assert !isempty(todel) "Empty model is not valid"
        end
        devs = zeros(T, length(todel) + 1)
        # cache fitted
        dict_devs = Dict{Set{Int}, T}()
        @inbounds for (id, del) in enumerate(todel)
            delcoef = selectcoef(f, Val(del))
            dev1 = get(dict_devs, delcoef, (push!(dict_devs, delcoef => deviance(trm, delcoef; kwargs...)); dict_devs[delcoef]))
            delete!(delcoef, del)
            dev2 = get(dict_devs, delcoef, (push!(dict_devs, delcoef => deviance(trm, delcoef; kwargs...)); dict_devs[delcoef]))
            devs[id] = dev1 - dev2
        end
        devs[end] = deviance(trm, 0; kwargs...)
    else
        todel = unique(assign)
        # ~ 0 + A, ~ 1
        start > 0 && @assert length(todel) > 1 "Empty model is not valid"
        devs = zeros(T, length(todel) + 1)
        @inbounds for (id, del) in enumerate(todel)
            devs[id] = deviance(trm, del; kwargs...)
        end
        devr = deviance(trm, 0; kwargs...)
        devs .-= devr
        devs[end] = devr
    end
    # deviance(trm, 0; kwargs...) 
    installbeta!(trm.model.pp) # ensure model unchanged
    # first(assign) == 1 || popfirst!(devs)
    tuple(devs...)
end

# for type3
exclude_X(p, assign, exclude::Int) = p.X[:, assign .!= exclude]
# for type 1/2
exclude_X(p, assign, exclude::Set{Int}) = p.X[:, map(!in(exclude), assign)]

# Backend for LinearModel
# Only used by type 2 SS
function deviance(trm::TableRegressionModel{<: LinearModel}, exclude)
    p, assign = trm.model.pp, trm.mm.assign

    # subset X, reset beta
    X = exclude_X(p, assign, exclude)
    p.beta0 = p.delbeta = repeat([0], size(X, 2))
    p.scratchbeta = similar(p.beta0)

    # cholesky 
    if isdensechol(p)
        F = Hermitian(float(X'X))
        p.chol = updatechol(p.chol, F)
    end

    # reset scratch
    p.scratchm1 = similar(X)
    p.scratchm2 = similar(p.chol.factors)

    isempty(trm.model.rr.wts) ? delbeta!(p, X, trm.model.rr.y) : delbeta!(p, X, trm.model.rr.y, trm.model.rr.wts)
    updateμ!(trm.model.rr, linpred(p, X))
    # installbeta is ommited
end 

isdensechol(::DensePredChol) = true
isdensechol(::SparsePredChol) = false

updatechol(::CholeskyPivoted, F::AbstractMatrix{<: BlasReal}) = 
    cholesky!(F, Val(true), tol = -one(eltype(F)), check = false)

updatechol(::Cholesky, F::AbstractMatrix{<: BlasReal}) = cholesky!(F)

# Backend for GeneralizedLinearModel
function deviance(trm::TableRegressionModel{<: GeneralizedLinearModel}, exclude; 
                    verbose::Bool = false, 
                    maxiter::Integer = 30, 
                    minstepfac::Real = 0.001,
                    atol::Real = 1e-6, 
                    rtol::Real = 1e-6, 
                    kwargs...)

    # Check arguments
    maxiter >= 1       || throw(ArgumentError("maxiter must be positive"))
    0 < minstepfac < 1 || throw(ArgumentError("minstepfac must be in (0, 1)"))

    # Extract fields and set convergence flag
    cvg, p, r, assign = false, trm.model.pp, trm.model.rr, trm.mm.assign
    lp = r.mu

    # subset X, reset beta
    X = exclude_X(p, assign, exclude)
    p.beta0 = p.delbeta = repeat([0], size(X, 2))
    p.scratchbeta = similar(p.beta0)
    
    # cholesky 
    if isdensechol(p)
        F = Hermitian(float(X'X))
        p.chol = updatechol(p.chol, F)
    end

    # reset scratch
    p.scratchm1 = similar(X)
    p.scratchm2 = similar(p.chol.factors)

    # Initialize β, μ, and compute deviance
    delbeta!(p, X, GLM.wrkresp(r), r.wrkwt)
    linpred!(lp, p, X)
    updateμ!(r, lp)
    installbeta!(p)
    devold = deviance(trm.model)

    for i = 1:maxiter
        f = 1.0 # line search factor
        local dev

        # Compute the change to β, update μ and compute deviance
        try
            delbeta!(p, X, r.wrkresid, r.wrkwt)
            linpred!(lp, p, X)
            updateμ!(r, lp)
            dev = deviance(trm.model)
        catch e
            isa(e, DomainError) ? (dev = Inf) : rethrow(e)
        end

        # Line search
        ## If the deviance isn't declining then half the step size
        ## The rtol*dev term is to avoid failure when deviance
        ## is unchanged except for rouding errors.
        while dev > devold + rtol*dev
            println("ok")
            f /= 2
            f > minstepfac || error("step-halving failed at beta0 = $(p.beta0)")
            try
                updateμ!(r, linpred(p, X, f))
                dev = deviance(trm.model)
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
    devold
end

# Linear prediction
function delbeta!(p::DensePredChol{<: BlasReal, <: Cholesky}, 
                    X::Matrix{<: BlasReal}, 
                    r::Vector{<: BlasReal})
    ldiv!(p.chol, mul!(p.delbeta, transpose(X), r))
    p
end
# β = (X'X)⁻¹X'y

function delbeta!(p::DensePredChol{<: BlasReal, <: CholeskyPivoted}, 
                    X::Matrix{<: BlasReal}, 
                    r::Vector{<: BlasReal})
    ch = p.chol
    beta = mul!(p.delbeta, adjoint(X), r)
    rnk = rank(ch)
    if rnk == length(beta)
        ldiv!(ch, beta)
    elseif size(ch.factors) == (0, 0) # Rank undtermined when the matrix is 0x0
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

# exterior weights
function delbeta!(p::DensePredChol{<: BlasReal, <: Cholesky}, 
                    X::Matrix{<: BlasReal}, 
                    r::Vector{<: BlasReal}, 
                    wt::Vector{<: BlasReal})
    scr = mul!(p.scratchm1, Diagonal(wt), X)
    cholesky!(Hermitian(mul!(cholfactors(p.chol), transpose(scr), X), :U))
    mul!(p.delbeta, transpose(scr), r)
    ldiv!(p.chol, p.delbeta)
    p
end

function delbeta!(p::DensePredChol{<: BlasReal, <: CholeskyPivoted}, 
                    X::Matrix{<: BlasReal}, 
                    r::Vector{<: BlasReal}, 
                    wt::Vector{<: BlasReal})
    cf = cholfactors(p.chol)
    piv = p.chol.piv
    cf .= mul!(p.scratchm2, adjoint(mul!(p.scratchm1, Diagonal(wt), X)), X)[piv, piv]
    cholesky!(Hermitian(cf, Symbol(p.chol.uplo)))
    ldiv!(p.chol, mul!(p.delbeta, transpose(p.scratchm1), r))
    p
end

# experimental
function delbeta!(p::SparsePredChol{<: BlasReal}, 
                    X::Matrix{<: BlasReal}, 
                    r::Vector{<: BlasReal}, 
                    wt::Vector{<: BlasReal})
    scr = mul!(p.scratch, Diagonal(wt), X)
    XtWX = X'*scr
    c = p.chol = cholesky(Symmetric{eltype(XtWX),typeof(XtWX)}(XtWX, 'L'))
    p.delbeta = c \ mul!(p.delbeta, adjoint(scr), r)
end

# skip updating beta0
linpred(p::LinPred, X::Matrix, f::Real = 1.0) = 
    linpred!(Vector{eltype(p.X)}(undef, size(p.X, 1)), p, X, f)

function linpred!(out, p::LinPred, X::Matrix, f::Real = 1.0)
    mul!(out, X, iszero(f) ? p.beta0 : broadcast!(muladd, p.scratchbeta, f, p.delbeta, p.beta0))
end

# Create nestedmodels
function nestedmodels(trm::TableRegressionModel{<: LinearModel}; null::Bool = true, kwargs...)
    f = formula(trm)
    null && (isnullable(trm.model.pp.chol) || (null = false))
    assign = unique(trm.mm.assign)
    pop!(assign)
    dropcollinear, range = null ? (false, union(0, assign)) : (true, assign)
    wts = trm.model.rr.wts
    trms = map(range) do id
        # create sub-formula, modify schema, create mf and mm
        subf = subformula(f.lhs, f.rhs, id)
        terms = setdiff(getterms(f.rhs), getterms(subf.rhs))
        schema = deepcopy(trm.mf.schema)
        for term in terms
            pop!(schema.schema, Term(term))
        end
        pair = collect(pairs(trm.mf.data))
        filter!(x->!in(x.first, terms), pair)
        mf = ModelFrame(subf, schema, (; pair...), trm.mf.model) 
        mm = ModelMatrix(mf)
        y = response(mf)
        TableRegressionModel(fit(trm.mf.model, mm.m, y; wts, dropcollinear, kwargs...), mf, mm)
    end
    (trms..., trm)
end

function nestedmodels(trm::TableRegressionModel{<: GeneralizedLinearModel}; null::Bool = true, kwargs...)
    f = formula(trm)
    # fit models
    null && (isnullable(trm.model) || (null = false))
    distr = trm.model.rr.d
    link = typeof(trm.model.rr).parameters[3]()
    wts = trm.model.rr.wts
    offset = trm.model.rr.offset
    assign = unique(trm.mm.assign)
    pop!(assign)
    range = null ? union(0, assign) : assign
    trms = map(range) do id
        # create sub-formula, modify schema, create mf and mm
        subf = subformula(f.lhs, f.rhs, id)
        terms = setdiff(getterms(f.rhs), getterms(subf.rhs))
        schema = deepcopy(trm.mf.schema)
        for term in terms
            pop!(schema.schema, Term(term))
        end
        pair = collect(pairs(trm.mf.data))
        filter!(x->!in(x.first, terms), pair)
        mf = ModelFrame(subf, schema, (; pair...), trm.mf.model) 
        mm = ModelMatrix(mf)
        y = response(mf)
        TableRegressionModel(fit(trm.mf.model, mm.m, y, distr, link; wts, offset, kwargs...), mf, mm)
    end
    (trms..., trm)
end


nestedmodels(::Type{LinearModel}, formula, data; null::Bool = true, kwargs...) = 
    nestedmodels(lm(formula, data; kwargs...); null, kwargs...)
nestedmodels(::Type{GeneralizedLinearModel}, formula, data, 
                distr::UnivariateDistribution, link::Link = canonicallink(distr); 
                null::Bool = true, kwargs...) = 
    nestedmodels(glm(formula, data, distr, link; kwargs...); null, kwargs...)
    
# Null model for CholeskyPivoted is unstable
## For nestedmodels
isnullable(::CholeskyPivoted) = false
isnullable(::Cholesky) = true
isnullable(::InverseLink) = false
isnullable(::InverseSquareLink) = false
isnullable(::Link) = true

## For deviances
isnullable(::LinearModel) = true
isnullable(model::GeneralizedLinearModel) = isnullable(typeof(model.rr).parameters[3]())
