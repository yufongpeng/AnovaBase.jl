

"""
    anova(X,y,allowrankdeficient::Bool=false;type::Int=1)
    anova(model::TableRegressionModel,allowrankdeficient::Bool=false;type::Int=1)

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`.

The keyword argument `type` specifies type of anova.

`TableRegressionModel` created by `lm` can also be analyzed.
"""
function anova(X,y,allowrankdeficient::Bool=false;type::Int=1)

    @assert (type in [1,2,3]) "Invalid type"
    model = fit(LinearModel,X,y,allowrankdeficient)
    mm = model.mm
    mf = model.mf
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    if type == 1
        ss = zeros(Float64,assign[end]-1)
        exclude = Set(assign[1:end])
        var_id = 1
        while var_id < assign[end]
            delete!(exclude,var_id)
            ss[var_id] = SS(model,exclude,allowrankdeficient)
            var_id += 1
        end
        sse = SS(model,0,allowrankdeficient)
        push!(ss,sse,0)
        ss = _diff(ss)
        popfirst!(df)
    elseif type == 2
        f = mf.f.rhs
        sse = SS(model,0,allowrankdeficient)
        ss = map(2:assign[end]) do ref
            SS(model,selectcoef(f,ref),allowrankdeficient) - SS(model,delete!(selectcoef(f,ref),ref),allowrankdeficient)
        end
        SS(model,0,allowrankdeficient) # back to full model
        push!(ss,sse)
        popfirst!(df)
    else
        sse = SS(model,0,allowrankdeficient)
        ss = map(1:assign[end]) do ref
            SS(model,ref,allowrankdeficient)-sse
        end
        SS(model,0,allowrankdeficient) # back to full model
        push!(ss,sse)
    end
    MSR = ss./df
    fstat = [MSR[1:(end-1)]./MSR[end]...,NaN]
    pvalue = [ccdf.(FDist.(df, df[end]), abs.(fstat))[1:(end-1)]...,NaN]
    AnovaResult(model,AnovaStats(type,size(mm.m, 1),ss,df,fstat,pvalue))
end

function anova(model::TableRegressionModel,allowrankdeficient::Bool=false;type::Int=1)

    @assert (type in [1,2,3]) "Invalid type"
    mf = model.mf
    mm = model.mm
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    if type == 1
        ss = zeros(Float64,assign[end]-1)
        sse = SS(model,0,allowrankdeficient)
        exclude = Set(assign[1:end])
        var_id = 1
        while var_id < assign[end]
            delete!(exclude,var_id)
            ss[var_id] = SS(model,exclude,allowrankdeficient)
            var_id += 1
        end
        push!(ss,sse,0)
        ss = _diff(ss)
        popfirst!(df)
    elseif type == 2
        f = mf.f.rhs
        sse = SS(model,0,allowrankdeficient)
        ss = map(2:assign[end]) do ref
            SS(model,selectcoef(f,ref),allowrankdeficient) - SS(model,delete!(selectcoef(f,ref),ref),allowrankdeficient)
        end
        push!(ss,sse)
        popfirst!(df)
    else
        sse = SS(model,0,allowrankdeficient)
        ss = map(1:assign[end]) do ref
            SS(model,ref,allowrankdeficient)-sse
        end
        push!(ss,sse)
    end
    MSR = ss./df
    fstat = [MSR[1:(end-1)]./MSR[end]...,NaN]
    pvalue = [ccdf.(FDist.(df, df[end]), abs.(fstat))[1:(end-1)]...,NaN]
    AnovaResult(model,AnovaStats(type,size(mm.m, 1),ss,df,fstat,pvalue))
end

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
end

function SS(model::TableRegressionModel, exclude::Set{Int}, pivot::Bool)
    p = model.model.pp
    assign = model.mm.assign
    X = view(p.X,:,map(x->!in(x,exclude),assign))
    p.beta0 = repeat([0],size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    installbeta!(p, X, model.model.rr.y)
    updateμ!(model.model.rr, linpred(p, X))
end

function installbeta!(p::DensePredChol{T,<:Cholesky}, X::SubArray, r::Vector{T}) where T<:BlasReal
    ldiv!(p.chol, mul!(p.beta0, transpose(X), r))
    p
end


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


"""
function anova2lm(model,include::Vector{Symbol}=[],exclude=::Vector{Symbol}=[])
    include = Set(Term.(include))
    exclude = Set(Term.(exclude))
    # for term in model.mf.f.rhs
    #   select terms
    # => new FormulaTerm(model.mf.f.lhs,new MatrixTerm)
    # for (term,info) in model.mf.schema.schema
    #   add to new dict
    #  => new schema
    # for sym in keys(model.mf.data)
    #   select Sym
    # => new ntuple data
    # f = apply_schema(f, sch, LinearModel)
    # ModelFrame(f, sch, data, LinearModel)
    # mm = ModelMatrix(mf)
    # y = response(mf)
    # TableRegressionModel(fit(LinearModel, mm.m, y, args...; kwargs...), mf, mm)
end
"""

# customize coef name
const TableModels = Union{TableStatisticalModel, TableRegressionModel}
StatsBase.coefnames(model::TableModels,anova::Val{:anova}) = coefnames(model.mf,anova)
StatsBase.coefnames(mf::ModelFrame,anova::Val{:anova}) = vectorize(coefnames(mf.f.rhs,anova))
StatsBase.coefnames(t::MatrixTerm,anova::Val{:anova}) = mapreduce(coefnames, vcat, t.terms, repeat([anova],length(t.terms)))
StatsBase.coefnames(t::FormulaTerm,anova::Val{:anova}) = (coefnames(t.lhs), coefnames(t.rhs))
StatsBase.coefnames(::InterceptTerm{H},anova::Val{:anova}) where {H} = H ? "(Intercept)" : []
StatsBase.coefnames(t::ContinuousTerm,anova::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::CategoricalTerm,anova::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::FunctionTerm,anova::Val{:anova}) = string(t.exorig)
StatsBase.coefnames(ts::StatsModels.TupleTerm,anova::Val{:anova}) = reduce(vcat, coefnames.(ts))

StatsBase.coefnames(t::InteractionTerm,anova::Val{:anova}) =
    kron_insideout((args...) -> join(args, " & "), vectorize.(coefnames.(t.terms,anova))...)

Base.show(io::IO, t::FunctionTerm) = print(io, "$(t.exorig)")

# subsetting coef names for type 2 anova
getterms(term::AbstractTerm) = Union{Symbol,Expr}[term.sym]
getterms(term::InterceptTerm) = Union{Symbol,Expr}[Symbol(1)]
getterms(term::InteractionTerm) = map(i->getterm(i),term.terms)
getterms(term::FunctionTerm) = Union{Symbol,Expr}[term.exorig]
getterm(term::AbstractTerm) = term.sym
getterm(term::FunctionTerm) = term.exorig
getterm(term::InterceptTerm) = Symbol(1)

isinteract(f::MatrixTerm,ref::Int,comp::Int) = issubset(getterms(f.terms[ref]),getterms(f.terms[comp]))
selectcoef(f::MatrixTerm,ref::Int) = Set([comp for comp in 1:length(f.terms) if isinteract(f,ref,comp)])

# coeftable implementation
function StatsBase.coeftable(model::AnovaResult; kwargs...)
    ct = coeftable(model.stats, kwargs...)
    cfnames = coefnames(model.model.mf,Val(:anova))
    push!(cfnames,"(Residual)")
    model.stats.type == 3 ||(popfirst!(cfnames))
    if length(ct.rownms) == length(cfnames)
        ct.rownms = cfnames
    end
    ct
end

function StatsBase.coeftable(stats::AnovaStats; kwargs...)
    ct = CoefTable(hcat([stats.dof...],[stats.ss...],[(stats.ss)./stats.dof...],[stats.fstat...],[stats.pval...]),
              ["DOF","Sum of Squares","Mean of Squares","F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    ct
    
end  # function coeftable

# show function that delegates to coeftable
function Base.show(io::IO, model::AnovaResult)
    ct = coeftable(model)
    println(io, typeof(model))
    println(io)
    println(io,"Type $(model.stats.type) ANOVA")
    println(io)
    println(io, model.model.mf.f)
    println(io)
    println(io,"Coefficients:")
    show(io, ct)
end


#if QR
   # model = AnovaModel(LmResp(y, similar(y, 0)),AnovaDensePredQR(mm.m))
    #effect = (transpose(model.pp.qr.Q)*model.rr.y).^2
    #whole_id = 1
    #var_id = 1
    #ss = zeros(Float64,assign[end])
    #while whole_id < length(assign)
    #    whole_id += 1
    #    var_id == assign[whole_id] || (var_id += 1)
    #    ss[var_id] += effect[whole_id]
    #end
    #ss[end] = sum(effect[(whole_id+1):end])
    #popfirst!(df)
#else