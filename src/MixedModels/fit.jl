# ==========================================================================================================
# Backend funcion

formula(model::MixedModel) = model.formula

lhs_no_intercept(lhs::String) = Set([lhs])
lhs_no_intercept(lhs) = Set(filter!(!=("(Intercept)"), lhs))

"""
    calcdof(model::LinearMixedModel)

Calculate degree of freedom of factors and residuals for linear mixed effect models
DOF of residuals are estimated by between-within method:
    dofᵢ = nobsᵢ - dofᵢ₋₁ - nfixᵢ
"""
function calcdof(model::LinearMixedModel)
    randoms = collect(formula(model).rhs) # ranef formula terms
    fixs = popfirst!(randoms) # fixef formula terms
    nfac = length(fixs.terms) # number of factors
    isempty(randoms) && (return repeat([nobs(model) - nfac], nfac))
    reterms = reverse(model.reterms) # Vector of ReMat
    fixef = model.feterm.cnames # fixef name
    nfix = length(fixef) # number of fix effects

    # Determine affected fix effects for each random effects
    affectfixef = Dict{String, Set{String}}()
    for pair in randoms
        lhs = coefnames(pair.lhs)
        lhs = lhs_no_intercept(lhs)
        rhs = coefnames(pair.rhs, Val(:anova))
        haskey(affectfixef, rhs) ? union!(affectfixef[rhs], lhs) : (affectfixef[rhs] = lhs)
    end 

    # Determines if fix effects vary at each random effects group
    within = zeros(Bool, nfix, size(reterms, 1) + 2)
    for (reid, remat) in enumerate(reterms)
        levels = unique(remat.refs)
        within[:, reid + 1] .= map(1:nfix) do fixid
            mapreduce(*, levels) do level
                vals = view(model.Xymat.wtxy, findall(==(level), remat.refs), fixid)
                val = first(vals)
                for i in vals
                    i != val && (return false)
                end
                true
            end
        end
    end 
    within[:, 1] = repeat([false], nfix) # population level
    within[1, 1] = true
    within[:, end] = repeat([true], nfix) # observation level

    # Turn affected fix effects into true
    ranef = ["", coefnames.(getproperty.(reterms, :trm), Val(:anova))..., ""]
    for (key, value) in affectfixef
        within[findfirst(in(value), fixef), findfirst(==(key), ranef)] = true
    end 
    # Find first true for each fix effect, aka level
    level = mapslices(within, dims = 2) do fixef
        findfirst(fixef)
    end 
    # Number of fix effects per level
    nfixperlv = zeros(Int, length(ranef))
    for i in level
        nfixperlv[i] += 1
    end 

    # dfᵢ = nobsᵢ - dfᵢ₋₁ - nfixᵢ
    n0 = 1
    n = nobs(model)
    nobsperlv = (n0, length.(getproperty.(reterms, :levels))..., n) # number of random effects observation per level
    ndiff = _diff(nobsperlv)
    dfperlv = ndiff .- nfixperlv[2:end]
    dfperlv = (last(dfperlv), dfperlv...)

    # Assign df to each factors based on the level
    assign = asgn(fixs)
    offset = 0
    intercept = first(assign) == 1
    intercept || (offset = 1)
    dfr = ntuple(last(assign) - offset) do i
        dfperlv[level[findlast(==(i + offset), assign)]] # use the last one to avoid wrong assign of the first arguments
    end
    df = dof(assign)
    intercept || popfirst!(df)
    tuple(df...), dfr
end

nestedmodels(::Type{<: LinearMixedModel}, f::FormulaTerm, tbl; 
                null::Bool = true, 
                wts = [], 
                contrasts = Dict{Symbol, Any}(), 
                progress::Bool = true, 
                REML::Bool = false) = 
    nestedmodels(fit(LinearMixedModel, f, tbl; wts, contrasts, progress, REML); null)

function nestedmodels(model::LinearMixedModel; null::Bool = true, kwargs...)
    f = formula(model)
    range = null ? (0:length(f.rhs[1].terms) - 1) : (1:length(f.rhs[1].terms) - 1)
    assign = asgn(first(f.rhs))
    REML = model.optsum.REML
    pivot = model.feterm.piv
    X = copy(model.X)
    y = copy(model.y)
    T = promote_type(Float64, eltype(y))
    size(X, 2) > 0 && (X = pivot == collect(1:size(X, 2)) ? X : X[:, pivot])
    reterms = deepcopy(model.reterms)
    sqrtwts = copy(model.sqrtwts)
    models = map(range) do id
        subf = subformula(f.lhs, f.rhs, id)
        cnames = coefnames(first(subf.rhs))
        # modify X by assign
        select = findall(x->x <= id, assign)
        subX = X[:, select]
        feterms = MixedModels.FeTerm{T}[]
        push!(feterms, MixedModels.FeTerm(subX, isa(cnames, String) ? [cnames] : collect(cnames)))
        feterm = only(feterms)
        Xy = MixedModels.FeMat(feterm, vec(y))
        reweight!(Xy, sqrtwts)
        A, L = createAL(reterms, Xy)
        lbd = foldl(vcat, lowerbd(c) for c in reterms)
        θ = foldl(vcat, getθ(c) for c in reterms)
        optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
        fill!(optsum.xtol_abs, 1.0e-10)
        lmm = LinearMixedModel(
            subf,
            reterms,
            Xy,
            feterm,
            sqrtwts,
            model.parmap,
            (n = length(y), p = feterm.rank, nretrms = length(reterms)),
            A,
            L,
            optsum,
            )
        fit!(lmm; REML, progress = true)
        lmm
    end
    (models..., model)
end

# For nestedmodels
isnullable(::LinearMixedModel) = true

# Specialized dof_residual
dof_residual(aov::AnovaResult{<: MixedModel, FTest}) = aov.tests.resdof