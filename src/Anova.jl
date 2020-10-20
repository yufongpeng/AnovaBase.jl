__precompile__()

module Anova

using GLM, MixedModels, Statistics, StatsBase, StatsModels, LinearAlgebra, Distributions, Reexport
import GLM: LinPredModel, LinearModel, LmResp, DensePred,
            DensePredChol, QRCompactWY, LinPred, installbeta!, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!
import StatsBase: fit!, fit, dof
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable, asgn
import LinearAlgebra.BlasReal
import Tables.istable
import Base: levels, isbetween

export
    # models
    AnovaResult, AnovaStats,AnovaStatsGrouped,

    # functions
    anova, anova_lm, lme, anova_lme

@reexport using GLM
@reexport using MixedModels
@reexport using StatsModels

"""
    AnovaResult{TableRegressionModel,AnovaStats}
    AnovaResult{MixednModel,AnovaStatsGrouped}

Returned object of `anova`.

## Fields

* `model` is the full model.
* `stats` contains result of ANOVA, including sum of squares, fstats, etc.
"""
struct AnovaResult{T,S}
    model::T
    stats::S
end


"""
    AnovaStats

Object contains result of ANOVA

## Fields

* `type`: the type of ANOVA.
* `nobs`: the number of observations.
* `ss`: sum of squares for each predictors. 
      if `type` is 3, than the first predictor is the intercept; otherwise, it is the first variable.
* `dof`: degre of freedom.
* `fstat`: f statiscics for each predictor.
* `pval`: p-values for each predictor.
"""
mutable struct AnovaStats
    type::Int
    nobs::Int
    ss::Vector{Float64}
    dof::Vector{Int}
    fstat::Vector{Float64}
    pval::Vector{Float64}
end

"""
    AnovaStatsGrouped

Object contains result of ANOVA from mixed-effect models

## Fields

* `ngroups`: number of groups for each random effect

For other fields, please see `AnovaStats`
"""
mutable struct AnovaStatsGrouped
    type::Int
    nobs::Int
    ngroups::Vector{Int}
    betweensubjects::Vector{Bool}
    ss::Vector{Float64}
    dof::Vector{Int}
    fstat::Vector{Float64}
    pval::Vector{Float64}
end

include("fit.jl")
include("termIO.jl")
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
# Appendix I: Type of sum of squares
```
For exmaple, a two-way ANOVA with factor A and B   
Type I  : SS(A) -> SS(B | A) -> SS(AB | A,B)  
Type II : SS(A | B) -> SS(B | A) -> SS(AB | A,B)    
Type III: SS(A | A,AB) -> SS(B | AB,B) -> SS(AB | A,B) # equivalent to linear regression  
```

# Appendix II: Degree of freedoms of for linear mixed-effect model

I: intercept, Aᵢ: between-subjects, Bⱼ: within-subjects, Cₖ: continuous (Not considering C|S), S: random-effect (1|S)  

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
For between-subject factors:
```
    dof(S|A) = Πaᵢ⋅(n-1)
```  
For within-subject factors: 
```     
    dof(B×C×S|A) = N-dof(I)-ΣₚΣᵣΣₛdof(ΠₚAᵢ×ΠₗBⱼ×ΠₛCₖ)-dof(S|A) 
                = Πaᵢ⋅(Πbⱼ-1)(n-1) 
```

# Appendix III: Sum of squares of error terms for linear mixed-effect model  
``` 
    SS(S|A) = SSR = sum(residuals(model).^2)    
    SS(B×C×S|A) = varest(model)*dof(B×C×S|A) 
```  

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
BP ~ time * drug + (1|pa), time is within-subject, drug is between-subject, pa is random effect.
     
"""
appendix
end

