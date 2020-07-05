__precompile__()

module Anova

using GLM, Statistics, StatsBase, StatsModels, LinearAlgebra, Distributions, Reexport
import GLM: LinPredModel, LinearModel, LmResp, DensePred,
            DensePredChol, QRCompactWY, LinPred, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!
import StatsBase: fit!, fit
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable
import LinearAlgebra.BlasReal
import Tables.istable

export
    # models
    AnovaModel, AnovaDensePredQR, AnovaDensePredChol, AnovaResult,

    # functions
    anova

@reexport using GLM

mutable struct AnovaModel{L<:LmResp,T<:LinPred} <: LinPredModel
    rr::L
    pp::T
end


mutable struct AnovaDensePredQR{T<:BlasReal} <: DensePred
    X::Matrix{T}                   # model matrix
    qr::QRCompactWY{T}
end

function AnovaDensePredQR(X::Matrix{T}) where T
    AnovaDensePredQR{T}(X,   qr(X))
end


mutable struct AnovaDensePredChol{T<:BlasReal,C} <: DensePred
    X::Matrix{T}                   # model matrix
    assign::Vector{Int}            # X' = view(X,:,assign.!=exclude)
    beta::Vector{T}                # base vector for coefficients
    chol::C
    scratchm1::Matrix{T}
    scratchm2::Matrix{T}
end

function AnovaDensePredChol(X::StridedMatrix,assign::Vector{Int}, pivot::Bool)
    F = Hermitian(float(X'X))
    T = eltype(F)
    F = pivot ? cholesky!(F, Val(true), tol = -one(T), check = false) : cholesky!(F)
    AnovaDensePredChol(AbstractMatrix{T}(X),
        assign,
        zeros(T,size(X, 2)),
        F,
        similar(X, T),
        similar(cholfactors(F)))
end

mutable struct AnovaResult{AnovaModel,T}
    model::AnovaModel
    type::Int
    nobs::Int
    ss::Vector{Float64}
    dof::Vector{Int}
    fstat::Vector{Float64}
    pval::Vector{Float64}
    mf::ModelFrame
    mm::ModelMatrix{T}
end

"""
Two-way anova:
Type I  : SS(A) -> SS(B | A) -> SS(AB | A,B)
Type II : SS(A | B) -> SS(B | A) -> SS(AB | A,B)
Type III: SS(A | A,AB), SS(B | AB,B), SS(AB | A,B) equivalent to linear regression

A, B: fixed effect -> denominator: SSR
A, B: random effect -> denominator: SS(AB | A,B)
A/ B: fixed effect/ random effect -> denominator: SS(AB | A,B)/ SSR
"""


"""
    anova(X,y,allowrankdeficient::Bool=false;type::Int=1)

The arguments `X` and `y` can be a `Matrix` and a `Vector` or a `Formula` and a `DataFrame`.

The keyword argument `type` specifies type of anova.
"""
function anova(X,y,allowrankdeficient::Bool=false;type::Int=1)

    @assert (type in [1,2,3]) "Invalid type"
    y, mf, mm = frame(AnovaModel,X,y,allowrankdeficient)
    df = Int.(dof(mm.assign))
    push!(df,Int(size(mm.m, 1)-sum(df)))
    assign = mm.assign
    model = AnovaModel(LmResp(y, similar(y, 0)),
                            AnovaDensePredChol(mm.m,assign,allowrankdeficient))
    if type == 1
        assign .-= 1
        ss = zeros(Float64,assign[end])
        sse = SS(model,-1,allowrankdeficient)
        exclude = Set(assign[2:end])
        var_id = 1
        while var_id <= assign[end]
            ss[var_id] = SS(model,exclude,allowrankdeficient)
            delete!(exclude,var_id)
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
    AnovaResult(model,type,size(mm.m, 1),ss,df,fstat,pvalue,mf,mm)
end

# customize fit! to prevent constructing model
function frame(::Type{T}, f::FormulaTerm, data, args...;
    contrasts::Dict{Symbol,<:Any} = Dict{Symbol,Any}(),kwargs...) where T<:StatisticalModel

    istable(data) || throw(ArgumentError("expected data in a Table, got $(typeof(data))"))
    cols = columntable(data)
    mf = ModelFrame(f, cols, model=T, contrasts=contrasts)
    mm = ModelMatrix(mf)
    y = response(mf)
    return y, mf, mm
end

# calculate SS
function SS(obj::AnovaModel, exclude::Int, pivot::Bool)
    p = obj.pp
    X = view(p.X,:,p.assign.!=exclude)
    p.beta = repeat([0],size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    delbeta!(p, X, obj.rr.y)
    updateμ!(obj.rr, linpred(p, X))
end

function SS(obj::AnovaModel, exclude::Set{Int}, pivot::Bool)
    p = obj.pp
    X = view(p.X,:,map(x->!in(x,exclude),p.assign))
    p.beta = repeat([0],size(X, 2))
    F = X'X
    p.chol = pivot ? cholesky!(F, Val(true), tol = -one(eltype(F)), check = false) : cholesky!(F)
    delbeta!(p, X, obj.rr.y)
    updateμ!(obj.rr, linpred(p, X))
end

function delbeta!(p::AnovaDensePredChol{T,<:Cholesky}, X::SubArray, r::Vector{T}) where T<:BlasReal
    ldiv!(p.chol, mul!(p.beta, transpose(X), r))
    p
end


function delbeta!(p::AnovaDensePredChol{T,<:CholeskyPivoted}, X::SubArray, r::Vector{T}) where T<:BlasReal
    ch = p.chol
    delbeta = mul!(p.beta, adjoint(X), r)
    rnk = rank(ch)
    if rnk == length(delbeta)
        ldiv!(ch, delbeta)
    else
        permute!(delbeta, ch.piv)
        for k=(rnk+1):length(delbeta)
            delbeta[k] = -zero(T)
        end
        LAPACK.potrs!(ch.uplo, view(ch.factors, 1:rnk, 1:rnk), view(delbeta, 1:rnk))
        invpermute!(delbeta, ch.piv)
    end
    p
end

linpred(p::AnovaDensePredChol, X::SubArray) = linpred!(Vector{eltype(p.X)}(undef, size(p.X, 1)), p, X)

function linpred!(out, p::AnovaDensePredChol, X::SubArray)
    mul!(out, X, p.beta)
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
const TableModels = Union{TableStatisticalModel, TableRegressionModel,AnovaResult}
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
    ct = CoefTable(hcat([model.dof...],[model.ss...],[(model.ss)./model.dof...],[model.fstat...],[model.pval...]),
              ["DOF","Sum of Squares","Mean of Squares","F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(model.dof)], 5, 4)
    cfnames = coefnames(model.mf,Val(:anova))
    push!(cfnames,"(Residual)")
    model.type == 3 ||(popfirst!(cfnames))
    if length(ct.rownms) == length(cfnames)
        ct.rownms = cfnames
    end
    ct
end

# show function that delegates to coeftable
function Base.show(io::IO, model::AnovaResult)
    ct = coeftable(model)
    println(io, typeof(model))
    println(io)
    println(io, model.mf.f)
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

end

