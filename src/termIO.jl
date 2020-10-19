
# customize coef name
const TableModels = Union{TableStatisticalModel, TableRegressionModel}
StatsBase.coefnames(model::TableModels,anova::Val{:anova}) = coefnames(model.mf,anova)
StatsBase.coefnames(mf::ModelFrame,anova::Val{:anova}) = begin
    v = vectorize(coefnames(mf.f.rhs,anova))
    push!(v,"Residual")
    v
end

StatsBase.coefnames(model::MixedModel,anova::Val{:anova}) = begin 
    v = vectorize(coefnames(model.formula.rhs[1],anova))
    push!(v,"Residual (between-subjects)","Residual (within-subjects)")
    v
end

StatsBase.coefnames(t::MatrixTerm,anova::Val{:anova}) = mapreduce(coefnames, vcat, t.terms, repeat([anova],length(t.terms)))
StatsBase.coefnames(t::FormulaTerm,anova::Val{:anova}) = (coefnames(t.lhs), coefnames(t.rhs))
StatsBase.coefnames(::InterceptTerm{H},anova::Val{:anova}) where {H} = H ? "(Intercept)" : []
StatsBase.coefnames(t::ContinuousTerm,anova::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::CategoricalTerm,anova::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::FunctionTerm,anova::Val{:anova}) = string(t.exorig)
StatsBase.coefnames(ts::StatsModels.TupleTerm,anova::Val{:anova}) = reduce(vcat, coefnames.(ts))

StatsBase.coefnames(t::InteractionTerm,anova::Val{:anova}) = begin
    join(coefnames.(t.terms,anova), " & ")
end
    
Base.show(io::IO, t::FunctionTerm) = print(io, "$(t.exorig)")

# subsetting coef names for type 2 anova
getterms(term::AbstractTerm) = Union{Symbol,Expr}[term.sym]
getterms(term::InterceptTerm) = Union{Symbol,Expr}[Symbol(1)]
getterms(term::InteractionTerm) = map(i->getterm(i),term.terms)
getterms(term::FunctionTerm) = Union{Symbol,Expr}[term.exorig]
getterm(term::AbstractTerm) = term.sym
getterm(term::FunctionTerm) = term.exorig
getterm(term::InterceptTerm) = Symbol(1)

isinteract(f::MatrixTerm,id1::Int,id2::Int) = issubset(getterms(f.terms[id1]),getterms(f.terms[id2]))
selectcoef(f::MatrixTerm,id::Int) = Set([comp for comp in 1:length(f.terms) if isinteract(f,id,comp)])

# unify formula api
formula(model::TableModels) = model.mf.f
formula(model::MixedModel) = model.formula

nlevels(term::CategoricalTerm) = length(term.contrasts.levels)
nlevels(term::ContinousTerm) = 1 
nlevels(term::InteractionTerm) = prod(nlevels.(term.terms))

# coeftable implementation
function StatsBase.coeftable(model::AnovaResult; kwargs...)
    ct = coeftable(model.stats, kwargs...)
    cfnames = coefnames(model.model,Val(:anova))
    if length(ct.rownms) == length(cfnames)
        ct.rownms = cfnames
    end
    ct
end

## customized for mixed model
function StatsBase.coeftable(stats::AnovaStats; kwargs...)
    ct = CoefTable(hcat([stats.dof...],[stats.ss...],[(stats.ss)./stats.dof...],[stats.fstat...],[stats.pval...]),
              ["DOF","Sum of Squares","Mean of Squares","F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    ct
    
end 

function StatsBase.coeftable(stats::AnovaStatsGrouped; kwargs...)
    ct = CoefTable(hcat([stats.dof...],[stats.betweensubjects...,repeat([NaN],length(stats.dof)-length(stats.betweensubjects))...],[stats.ss...],[(stats.ss)./stats.dof...],[stats.fstat...],[stats.pval...]),
              ["DOF","Between-subjects","Sum of Squares","Mean of Squares","F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 6, 4)
    ct

end

# show function that delegates to coeftable
function Base.show(io::IO, model::AnovaResult)
    ct = coeftable(model)
    println(io, typeof(model))
    println(io)
    println(io,"Type $(model.stats.type) ANOVA")
    println(io)
    println(io, formula(model.model))
    println(io)
    println(io,"Coefficients:")
    show(io, ct)
end