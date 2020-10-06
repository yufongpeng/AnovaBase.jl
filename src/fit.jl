

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
    if type == 1
        ss = zeros(Float64,assign[end]-1)
        sse = SS(model,0,allowrankdeficient)
        exclude = Set(assign[1:end])
        ss = map(1:(assign[end]-1)) do id
            delete!(exclude,id)
            SS(model,exclude,allowrankdeficient)
        end
        push!(ss,sse,0)
        ss = _diff(ss)
        popfirst!(df)
    elseif type == 2
        f = model.mf.f.rhs
        sse = SS(model,0,allowrankdeficient)
        ss = map(2:assign[end]) do id
            SS(model,selectcoef(f,id),allowrankdeficient) - SS(model,delete!(selectcoef(f,id),id),allowrankdeficient)
        end
        push!(ss,sse)
        popfirst!(df)
    else
        sse = SS(model,0,allowrankdeficient)
        ss = map(1:assign[end]) do id
            SS(model,id,allowrankdeficient)-sse
        end
        push!(ss,sse)
    end
    MSR = ss./df
    fstat = [MSR[1:(end-1)]./MSR[end]...,NaN]
    pvalue = [ccdf.(FDist.(df, df[end]), abs.(fstat))[1:(end-1)]...,NaN]
    AnovaResult(model,AnovaStats(type,size(mm.m, 1),ss,df,fstat,pvalue))
end

anova_lm(X, y, allowrankdeficient::Bool = false; type::Int=1) = 
        anova(LinearModel, X, y, allowrankdeficient, type = type)

"""
anova for LinearMixedModel


function anova(model::type{LinearMixedModel}, f::FormulaTerm, tbl; 
               contrasts = Dict{Symbol,Any}(), wts = [], type::Int = 1, within = nothing)
    model = LinearMixedModel(f, tbl, contrasts, wts)
    return anova(model, type = type, within = within)
end

function anova(model::MixedModel; type::Int = 1, within = nothing)
    @assert (type in [1,2,3]) "Invalid type"
    fet = model.formula.rhs[1]
    ret = model.formula.rhs[2:end]
    @assert (length(ret) == 1) "Multi-id design is not implemented now"
    ## Determine between/within
    if isnothing(within)
        # for each id there is a fixed number of replicate of each within variable
        # Ex:
        # pa    time    BP      drug
        # 1     am      110     ACEI
        # 1     pm      130     ACEI
        # 2     am      120     ARB
        # 2     pm      125     ARB
        # 3     am      140     ACEI
        # 3     pm      165     ACEI
        # 4     am      130     ARB
        # 4     pm      145     ARB
        # BP ~ time*drug + (1|pa), time is within, drug is between, pa is id.
        #for fev in fet
        #end
    end 
    ## 
end

anova_lme(X, y; type::Int=1) = 
        anova(LinearMixedModel, X, y, type = type)

# Determine within subjects vaiable
function iswithin(fet::MixedModel.FeMat,ret::ReMat)

end
"""

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

_diff(v::Vector{T}) where T = map(i->(v[i]-v[i+1]),1:(length(v)-1))

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


function anova_QR()
    @assert (type in [1,2,3]) "Invalid type"
    """
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
    """
end