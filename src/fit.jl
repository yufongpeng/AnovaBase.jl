

"""
    anova(X,y,allowrankdeficient::Bool=false;type::Int=1)
    anova(model::TableRegressionModel,allowrankdeficient::Bool=false;type::Int=1)

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`.

The keyword argument `type` specifies type of anova.

`TableRegressionModel` created by `lm` can also be analyzed.
"""
function anova(model::Type{T}, X, y, 
               allowrankdeficient::Bool = false; type::Int = 1) where {T <:StatisticalModel}
    model = fit(model,X,y,allowrankdeficient)
    return anova(model,allowrankdeficient,type=type)
end

function anova(model::TableRegressionModel, 
               allowrankdeficient::Bool = false; type::Int = 1)
    @assert (type in [1,2,3]) "Invalid type"
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    f = model.mf.f.rhs
    if type == 1
        exclude = Set(assign[1:end])
        ss = map(1:(assign[end])) do id
            delete!(exclude,id)
            SS(model,exclude,allowrankdeficient)
        end
        push!(ss,0)
        ss = _diff(ss)
    elseif type == 2
        sse = SS(model,0,allowrankdeficient)
        ss = map(1:assign[end]) do id
            ifelse(id  == 1,SS(model,Set(assign[2:end]),allowrankdeficient),
            SS(model,selectcoef(f,id),allowrankdeficient) - SS(model,delete!(selectcoef(f,id),id),allowrankdeficient))
        end
        push!(ss,sse)
    else
        sse = SS(model,0,allowrankdeficient)
        ss = map(1:assign[end]) do id
            SS(model,id,allowrankdeficient)-sse
        end
        push!(ss,sse)
    end
    width(f.terms[1]) == 0 && (popfirst!(df);popfirst!(ss))
    MSR = ss./df
    fstat = [MSR[1:(end-1)]./MSR[end]...,NaN]
    pvalue = [ccdf.(FDist.(df, df[end]), abs.(fstat))[1:(end-1)]...,NaN]
    AnovaResult(model,AnovaStats(type,size(mm.m, 1),ss,df,fstat,pvalue))
end

anova_lm(X, y, allowrankdeficient::Bool = false; type::Int=1) = 
        anova(LinearModel, X, y, allowrankdeficient, type = type)

            
# calculate SS
function SS(model::TableRegressionModel, exclude::Int, pivot::Bool)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X,:,assign.!=exclude)
    p.beta0 = repeat([0],size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    installbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 3

function SS(model::TableRegressionModel, exclude::Set{Int}, pivot::Bool)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X,:,map(x->!in(x,exclude),assign))
    p.beta0 = repeat([0],size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    installbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end # for type 1 and 2

function installbeta!(p::DensePredChol{T,<:Cholesky}, X::SubArray, r::Vector{T}) where T<:BlasReal
    ldiv!(p.chol, mul!(p.beta0, transpose(X), r))
    p
end
# β = (X'X)⁻¹X'y

function installbeta!(p::DensePredChol{T,<:CholeskyPivoted}, X::SubArray, r::Vector{T}) where T<:BlasReal
    ch = p.chol
    beta = mul!(p.beta0, adjoint(X), r)
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

linpred(p::DensePredChol, X::SubArray) = linpred!(Vector{eltype(p.X)}(undef, size(p.X, 1)), p, X)

function linpred!(out, p::DensePredChol, X::SubArray)
    mul!(out, X, p.beta0)
end

_diff(v::Vector{T}) where T = map(i->(i>1) ? (v[i-1]-v[i]) : (v[i]),1:(length(v)))

# calculate dof from model.mm.assign
function dof(v::Vector{Int})
    dofv = zeros(v[end])
    prev = 1
    ind = 1
    while ind <= length(v)
        v[ind] == prev || (prev = v[ind])
        dofv[prev] += 1
        ind += 1
    end
    dofv
end


"""
anova for LinearMixedModel
"""


function anova(model::Type{LinearMixedModel}, f::FormulaTerm, tbl; 
                wts = [], 
                contrasts = Dict{Symbol,Any}(), 
                verbose::Bool = false, 
                REML::Bool = true, 
                type::Int = 1, 
                between = nothing)
    model = lme(f, tbl, wts = wts, contrasts = contrasts, verbose = verbose, REML = REML)
    return anova(model, type = type, between = between)
end

function anova(model::MixedModel; type::Int = 1, 
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
    invvarfix = cholesky(varβ)\Matrix(I,size(varβ)...) |> Hermitian
    invvarfixchol = cholesky(invvarfix).L # column factor contains between factor should × -1
    model.optsum.REML || adjust_sigma && (invvarfixchol = invvarfixchol/sqrt(nobs(model)/(nobs(model) - size(invvarfixchol,1))))
    fs = invvarfixchol'β

    assign = asgn(fet)
    # Determine between/within
    btw = isbetween(fet,assign,remat[1],model.X) # to be modify for multiple random effects 
    isnothing(between) || (btw[between] .= true)
    intercept = width(fet.terms[1]) == 1
    ngroups = map(x->size(x)[2],remat)
    nbetween = Int(prod(nlevels.(fet.terms[btw])))
    n = ngroups[1]/nbetween 
    btw = intercept ? (btw) : (btw[2:end])

    last = assign[end]-assign[1]+1
    fstat = zeros(Float64,dof(model)-length(β)+last) 
    ss = copy(fstat)
    df = zeros(Int64,dof(model)-length(β)+last)

    fstat[(last+1):end] .= NaN
    df[last+1] = nbetween*(n-1) # to be modify for multiple random effects 
    df[end] = nobs(model) - sum(df) - length(β)
    ss[last+1] = sum(residuals(model).^2) # to be modify for multiple random effects 
    ss[end] = varest(model)*df[end]

    df[1:last] .= intercept ? dof(assign) : dof(assign)[2:end]
    if type == 1
        fstat[1:last] .= map((loc,factor)->sum(fs[assign.==factor].^2)/df[loc],
                            1:last,unique(assign))
    else 
        varfixchol = inv(invvarfixchol)
        ## modify this ugly code
        fstat[1:last] .= map((loc,factor)->sum((qr(varfixchol[:,assign.== factor]).Q'*fs)[1:df[loc]].^2)/df[loc],
            1:last,unique(assign))
    end
    ss[1:last] .= map(1:last) do id
        (btw[id] ? fstat[id]*ss[last+1]*df[id]/df[last+1] : fstat[id]*varest(model)*df[id])
    end
    pvalue = map(1:lastindex(fstat)) do id
        if id > last
            NaN
        elseif btw[id]
            ccdf(FDist(df[id], df[last+1]), abs(fstat[id]))
        else
            ccdf(FDist(df[id], df[end]), abs(fstat[id]))
        end
    end
    AnovaResult(model,AnovaStatsGrouped(type,nobs(model),ngroups,Bool.(btw),ss,df,fstat,pvalue))
end

anova_lme(X, y; type::Int = 1,
        wts = [],
        contrasts = Dict{Symbol,Any}(), 
        verbose::Bool = false, 
        REML::Bool = true) = 
        anova(LinearMixedModel, X, y, type = type,wts =  wts, contrasts = contrasts, verbose = verbose, REML = REML)

# assign
function assigncols(fet::MatrixTerm)
    start = 0
    assign = map(fet.terms) do term
        wid = width(term)
        start += wid
        (start-wid+1):start
    end
end

# Determine between subjects vaiable
function isbetween(fet::MatrixTerm,assign::Array{Int64,1},remat::ReMat,X::Matrix) where N
    n = length(fet.terms)
    between = ones(Bool,n)
    select = 1:length(assign)
    for id in 1:size(remat,2)
        loc = findall(==(1),view(remat,:,id))
        x = view(X,loc,:)
        for level in select
            if length(unique(x[:,level])) > 1
                factor = assign[level]
                select = select[assign[select].!= factor]
                between[factor] = false
            end
        end
    end
    between[1] = false
    between
end

"""
for each id there is a fixed number of replicate of each within variable
Ex:
pa    time    BP      drug
1     am      110     ACEI
1     pm      130     ACEI
2     am      120     ARB
2     pm      125     ARB
3     am      140     ACEI
3     pm      165     ACEI
4     am      130     ARB
4     pm      145     ARB
BP ~ time*drug + (1|pa), time is within, drug is between, pa is id.


I: intercept, Aᵢ: between, Bⱼ: within, Cₖ: continuous (Not considering C|S), S: random
Πaᵢ⋅Πbⱼ⋅n = N, n(S) = Πaᵢ⋅n
Factors: 
    dof(I) = 1
    dof(Aᵢ) = aᵢ-1
    dof(Bⱼ) = bⱼ-1
    dof(Cₖ) = 1

    In general: 
    dof(ΠₚAᵢ×ΠᵣBⱼ×ΠₛCₖ) = Πₚ(aᵢ-1)⋅Πᵣ(bⱼ-1) 

Error term:
    dof(S|A) = Πaᵢ⋅(n-1)  # for between subject factors
    dof(B×C×S|A) = N-dof(I)-ΣₚΣᵣΣₛdof(ΠₚAᵢ×ΠₗBⱼ×ΠₛCₖ)-dof(S|A) 
               = Πaᵢ⋅(Πbⱼ-1)(n-1) # for within subject factors

SS(S|A) = SSR = sum(residuals(model).^2)
SS(B×C×S|A) = varest(model)*dof(B×C×S|A)

S₁|S₂|E
S₁|A₁ = S₁+E
S₂|A₂ = S₂+E
B×C×(S₁|A₁+S₂|A₂) = E



function group(start::Int64,term::CategoricalTerm,X::Matrix)
    wid = width(term)
    obvs = Array{Array{Int64,1},1}()
    for i in 1:(wid-1)
        obvs[i] = findall(==(1),X[:,start+i])
    end
    obvs
end

function group(start::Int64,term::InteractionTerm,X::Matrix)
    wid = width(term)
    obvs = Array{Array{Int64,1},1}()
    for i in 1:(wid-1)
        obvs[i] = findall(==(1),X[:,start+i])
    end
    obvs
end

group(start::Int64,term::InterceptTerm,X::Matrix) = [collect(1:size(X)[1])]

function group(start::Int64,term::ContinuousTerm,X::Matrix)
    vals = unique(X[:,start+1])
    obvs = Array{Array{Int64,1},1}()
    for (i,val) in enumertae(vals)
        obvs[i] = findall(==(val),X[:,start+1])
    end
    obvs
end

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