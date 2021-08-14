# ==========================================================================================================
# Backend funcion

formula(model::MixedModel) = model.formula

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

nestedmodels(::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
    null::Bool = true, 
    wts = [], 
    contrasts = Dict{Symbol, Any}(), 
    verbose::Bool = false, 
    REML::Bool = false) = 
    nestedmodels(fit(LinearMixedModel, f, tbl, 
        wts =  wts, contrasts = contrasts, verbose = verbose, REML = REML); null = null)

function nestedmodels(model::LinearMixedModel; null::Bool = true, kwargs...)
    f = formula(model)
    range = null ? (0:length(f.rhs[1].terms) - 1) : (1:length(f.rhs[1].terms) - 1)
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