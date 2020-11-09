__precompile__()

module Anova

using GLM, MixedModels, Statistics, StatsBase, StatsModels, LinearAlgebra, Distributions, HypothesisTests, Reexport, DataFrames, Printf
import GLM: LinPredModel, LinearModel, LmResp, DensePred,
            DensePredChol, SparsePredChol, QRCompactWY, LinPred, installbeta!, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!, glm, AbstractGLM, FP, SparseMatrixCSC, Link
import StatsBase: fit!, fit, dof
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable, asgn
import LinearAlgebra.BlasReal
import Tables.istable
import Base: isbetween, show

export
    # models
    AnovaResult, AnovaStats,AnovaStatsGrouped,

    # functions
    anova, anova_lm, lme, anova_lme,

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT

@reexport using GLM
@reexport using MixedModels
@reexport using StatsModels

# Wrapper for ANOVA
abstract type AbstractAnovaStats end
abstract type SequentialAnovaStats end

"""
    AnovaResult{TableRegressionModel,AnovaStats}
    AnovaResult{LinearMixednModel,AnovaStatsGrouped}

Returned object of `anova`.

## Fields

* `model` is the full model or an array of  tested models.
* `stats` contains result of ANOVA, including sum of squares, fstats, etc.
"""
struct AnovaResult{T,S}
    model::T 
    stats::S

    AnovaResult(model::T, stats::S) where {T <: StatisticalModel, S <: AbstractAnovaStats} = new{T,S}(model, stats)
    AnovaResult(model::T, stats::S) where {T <: NTuple{N, M}, S <: SequentialAnovaStats} where {N, M <: StatisticalModel} = new{T,S}(model, stats)
end


"""
    AnovaStats <: AbstractAnovaStats

Object contains result of ANOVA

## Fields

* `type`: the type of ANOVA.
* `nobs`: the number of observations.
* `ss`: sum of squares for each predictors. 
* `dof`: degre of freedom.
* `fstat`: f statiscics for each predictor.
* `pval`: p-values for each predictor.
"""
mutable struct AnovaStats <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::Vector{Int}
    ss::Vector{Float64}
    fstat::Vector{Float64}
    pval::Vector{Float64}
end

"""
    AnovaStatsGrouped <: AbstractAnovaStats

Object contains result of ANOVA from mixed-effect models

## Fields

* `ngroups`: number of groups for each random effect
* `betweensubjects`: whether a variable is between-subjects.

For other fields, please see `AnovaStats`
"""
mutable struct AnovaStatsGrouped <: AbstractAnovaStats
    type::Int
    nobs::Int
    ngroups::Vector{Int}
    betweensubjects::Vector{Bool}
    dof::Vector{Int}
    ss::Vector{Float64}
    fstat::Vector{Float64}
    pval::Vector{Float64}
end




"""
"""
mutable struct AnovaStatsF{N} <: SequentialAnovaStats where N
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    fstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
"""
mutable struct AnovaStatsLRT{N} <: SequentialAnovaStats where N
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    lrstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
"""
mutable struct AnovaStatsRao <: SequentialAnovaStats
end

"""
"""
mutable struct AnovaStatsCp <: SequentialAnovaStats
end

abstract type GoodnessOfFit end

struct FTest <: GoodnessOfFit end

struct LikelihoodRatioTest <: GoodnessOfFit end

const LRT = LikelihoodRatioTest

const FixDispDist = Union{Bernoulli, Binomial, Poisson}

# Alias for LinearMixedModel
"""
    lme(f::FormulaTerm, tbl; wts, contrasts, verbose, REML)
An alias for `fit(LinearMixedModel,f, tbl; wts, contrasts, verbose, REML)`.

"""
lme(f::FormulaTerm, tbl; 
    wts = [], 
    contrasts = Dict{Symbol,Any}(), 
    verbose::Bool = false, 
    REML::Bool = false) = 
    fit(LinearMixedModel,f, tbl, 
    wts =  wts, contrasts = contrasts, verbose = verbose, REML = REML)

"""
    glm(f, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...)

Automated transform dependent variable into 0/1 for family `Binomial`
"""
glm(f::FormulaTerm, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...) =
    fit(GeneralizedLinearModel, f, 
        combine(df, : , f.lhs.sym => ByRow(x -> x == unique(df[:, f.lhs.sym])[end]) => f.lhs.sym), 
        d, l, args...; kwargs...)

glm(f::FormulaTerm, df::DataFrame, d::Poisson, l::GLM.Link, args...; kwargs...) =
    fit(GeneralizedLinearModel, f, 
        combine(df, : , f.lhs.sym => ByRow(x -> x == unique(df[:, f.lhs.sym])[end]) => f.lhs.sym), 
        d, l, args...; kwargs...)

_diff(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i+1]-t[i], N-1)
_diffn(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i]-t[i+1], N-1)

include("fit.jl")
include("termIO.jl")

# Appendix
"""
# Appendix I: Type of sum of squares

For exmaple, a two-way ANOVA with factor A and B \n  
Type I  : SS(A) -> SS(B | A) -> SS(AB | A,B)  \n
Type II : SS(A | B) -> SS(B | A) -> SS(AB | A,B) \n   
Type III: SS(A | A,AB) -> SS(B | AB,B) -> SS(AB | A,B) # equivalent to linear regression  \n


# Appendix II: Degree of freedoms of for linear mixed-effect model

I: intercept, Aᵢ: between-subjects, Bⱼ: within-subjects, Cₖ: continuous (Not considering C|S), S: random-effect (1|S)  \n
n: number of observations in each cells, n(S) = Πaᵢ⋅n = number of groups, Πaᵢ⋅Πbⱼ⋅n = N = total number of observations
## Factors 
``` 
    dof(I) = 1
    dof(Aᵢ) = aᵢ-1
    dof(Bⱼ) = bⱼ-1
    dof(Cₖ) = 1
``` 
In general:
```  
    dof(ΠₚAᵢ×ΠᵣBⱼ×ΠₛCₖ) = Πₚ(aᵢ-1)⋅Πᵣ(bⱼ-1) 
``` 
## Error term
For between-subjects factors:
```
    dof(S|A) = Πaᵢ⋅(n-1)
```  
For within-subjects factors: 
```     
    dof(B×C×S|A) = N-dof(I)-ΣₚΣᵣΣₛdof(ΠₚAᵢ×ΠₗBⱼ×ΠₛCₖ)-dof(S|A) 
                = Πaᵢ⋅(Πbⱼ-1)(n-1) 
```

# Appendix III: Sum of squares of error terms for linear mixed-effect model  

SS(S|A) = SSR = `sum(residuals(model).^2)`  \n  
SS(B×C×S|A) = `varest(model)`*dof(B×C×S|A) \n


# Appendix IV: Examples for linear mixed-effect model
Balanced design: for each level of random effect, there is a fixed number of observations of each within variable
```
16×4 DataFrame    
│ Row │ pa   │ time │ BP    │ drug │  
│     │ Cat… │ Cat… │ Int64 │ Cat… │  
├─────┼──────┼──────┼───────┼──────┤  
│ 1   │ 1    │ am   │ 110   │ ACEI │  
│ 2   │ 1    │ am   │ 112   │ ACEI │  
│ 3   │ 1    │ pm   │ 130   │ ACEI │
│ 4   │ 1    │ pm   │ 135   │ ACEI │
│ 5   │ 2    │ am   │ 120   │ ARB  │
│ 6   │ 2    │ am   │ 127   │ ARB  │
│ 7   │ 2    │ pm   │ 125   │ ARB  │
│ 8   │ 2    │ pm   │ 121   │ ARB  │
│ 9   │ 3    │ am   │ 140   │ ACEI │
│ 10  │ 3    │ am   │ 141   │ ACEI │
│ 11  │ 3    │ pm   │ 165   │ ACEI │
│ 12  │ 3    │ pm   │ 152   │ ACEI │
│ 13  │ 4    │ am   │ 130   │ ARB  │
│ 14  │ 4    │ am   │ 124   │ ARB  │
│ 15  │ 4    │ pm   │ 145   │ ARB  │
│ 16  │ 4    │ pm   │ 151   │ ARB  │ 
```
BP ~ time * drug + (1|pa), time is within-subjects, drug is between-subjects, pa is random effect.
     
"""
appendix
end

