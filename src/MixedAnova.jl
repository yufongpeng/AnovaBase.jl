module MixedAnova

using GLM, MixedModels, FixedEffectModels, Statistics, StatsBase, StatsModels, LinearAlgebra, Distributions, Reexport, DataFrames, Printf
import GLM: LinPredModel, LinearModel, LmResp, DensePred,
            DensePredChol, SparsePredChol, QRCompactWY, LinPred, installbeta!, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!, glm, AbstractGLM, FP, SparseMatrixCSC, Link
import MixedModels: FeMat, createAL, reweight!, getθ
import StatsBase: fit!, fit, dof, coefnames
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable, asgn, collect_matrix_terms
import LinearAlgebra.BlasReal
import Tables.istable
import Base: show

export
    # wrappers
    AnovaResult, FixedAnovaStatsF, MixedAnovaStatsF, FixedAnovaStatsLRT, 
    MixedAnovaStatsLRT, NestedAnovaStatsF, NestedAnovaStatsLRT,

    # anova functions
    anova, anova_lm, anova_lme, anova_glm,

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, canonicalgoodnessoffit, 

    # Others
    lme, lfe, nestedmodels, calcdof, formula, getterms, getterm

@reexport using GLM
@reexport using MixedModels
@reexport using FixedEffectModels

# Wrapper for ANOVA
abstract type AbstractAnovaStats end
abstract type NestedAnovaStats end

"""
    AnovaResult{TableRegressionModel,FixedAnovaStats}
    AnovaResult{LinearMixednModel,MixedAnovaStats}

Returned object of `anova`.

## Fields

* `model` is the full model or an tuple of  tested models.
* `stats` contains result of ANOVA, including sum of squares, fstats, etc.
"""
struct AnovaResult{T,S}
    model::T 
    stats::S

    AnovaResult(model::T, stats::S) where {T <: StatisticalModel, S <: AbstractAnovaStats} = new{T,S}(model, stats)
    AnovaResult(model::T, stats::S) where {T <: NTuple{N, M}, S <: AbstractAnovaStats} where {N, M <: StatisticalModel} = new{T,S}(model, stats)
    AnovaResult(model::T, stats::S) where {T <: NTuple{N, M}, S <: NestedAnovaStats} where {N, M <: StatisticalModel} = new{T,S}(model, stats)
end


# Use tuple/svector instead vector?
"""
    FixedAnovaStatsF <: AbstractAnovaStats

Object contains result of ANOVA

## Fields

* `type`: the type of ANOVA.
* `nobs`: the number of observations.
* `deviance`: deviance (sum of squares) for each predictor. 
* `dof`: degree of freedom.
* `fstat`: f statiscics for each predictor.
* `pval`: p-values for each predictor.
"""
mutable struct FixedAnovaStatsF{M, N} <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    fstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
    MixedAnovaStatsF <: AbstractAnovaStats

Object contains result of ANOVA from linear mixed-effect models

## Fields

* `resdof`: degree of freedom for residuals.
* `dofinfo`: dictionary contains information for `redof`.
    * `level`: evaluating level for each variable.
    * `lvdof`: dgree of freedom for each level.

For other fields, please see `FixedAnovaStatsF`.
"""
mutable struct MixedAnovaStatsF{M, N} <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::NTuple{N, Int}
    resdof::NTuple{N, Int}
    fstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
    dofinfo::Dict{Symbol, Tuple}
end

"""
    FixedAnovaStatsF <: AbstractAnovaStats

Object contains result of ANOVA

## Fields

* `type`: the type of ANOVA.
* `nobs`: the number of observations.
* `deviance`: deviance (sum of squares) for each predictor. 
* `dof`: degree of freedom.
* `lrstat`: likelihoood-ratio for each predictor.
* `pval`: p-values for each predictor.
"""
mutable struct FixedAnovaStatsLRT{M, N} <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    lrstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
    MixedAnovaStatsF <: AbstractAnovaStats

Object contains result of ANOVA from linear mixed-effect models

See `FixedAnovaStatsLRT`

"""
mutable struct MixedAnovaStatsLRT{M, N} <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    lrstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
    NestedAnovaStatsF <: NestedAnovaStats

Object contains result of ANOVA for nested models.

## Fields

* `nobs`: the number of observations.
* `deviance`: deviance (sum of squares) for each model. 
* `dof`: degree of freedom.
* `fstat`: f statiscics for each model.
* `pval`: p-values for each model.
"""
mutable struct NestedAnovaStatsF{N} <: NestedAnovaStats
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    fstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
    NestedAnovaStatsLRT <: NestedAnovaStats

Object contains result of ANOVA

## Fields

* `nobs`: the number of observations.
* `deviance`: deviance (sum of squares) for each model. 
* `dof`: degree of freedom.
* `lrstat`: likelihoood-ratio for each model.
* `pval`: p-values for each model.
"""
mutable struct NestedAnovaStatsLRT{N} <: NestedAnovaStats
    nobs::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    lrstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
end

"""
"""
mutable struct NestedAnovaStatsRao <: NestedAnovaStats
end

"""
"""
mutable struct NestedAnovaStatsCp <: NestedAnovaStats
end

const AnovaStatsF = Union{FixedAnovaStatsF, MixedAnovaStatsF, NestedAnovaStatsF}
const AnovaStatsLRT = Union{FixedAnovaStatsLRT, MixedAnovaStatsLRT, NestedAnovaStatsLRT}
const AnovaStatsRao = Union{NestedAnovaStatsRao}
const AnovaStatsCp = Union{NestedAnovaStatsCp}

abstract type GoodnessOfFit end

struct FTest <: GoodnessOfFit end

struct LikelihoodRatioTest <: GoodnessOfFit end

const LRT = LikelihoodRatioTest

const FixDispDist = Union{Bernoulli, Binomial, Poisson}

# Alias for LinearMixedModel
"""
    lme(f::FormulaTerm, tbl; wts, contrasts, verbose, REML)

An alias for `fit(LinearMixedModel, f, tbl; wts, contrasts, verbose, REML)`.

"""
lme(f::FormulaTerm, tbl; 
    wts = [], 
    contrasts = Dict{Symbol,Any}(), 
    verbose::Bool = false, 
    REML::Bool = false) = 
    fit(LinearMixedModel, f, tbl, 
    wts =  wts, contrasts = contrasts, verbose = verbose, REML = REML)

"""
    glm(f, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...)

Automatically transform dependent variable into 0/1 for family `Binomial`
"""
glm(f::FormulaTerm, df::DataFrame, d::Binomial, l::GLM.Link, args...; kwargs...) = 
    fit(GeneralizedLinearModel, f, 
        combine(df, : , f.lhs.sym => ByRow(x -> x == unique(df[:, f.lhs.sym])[end]) => f.lhs.sym), 
        d, l, args...; kwargs...)

"""
    lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...)

An alias for `reg`
"""
lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...) = 
    reg(df, formula, vcov, kwargs...)

_diff(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i + 1] - t[i], N - 1)
_diffn(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i] - t[i + 1], N - 1)

include("fit.jl")
include("termIO.jl")
include("anova.jl")

# Appendix
"""
# Appendix I: Type of sum of squares

For exmaple, a two-way ANOVA with factor A and B \n  
Type I  : SS(A) -> SS(B | A) -> SS(AB | A,B)  \n
Type II : SS(A | B) -> SS(B | A) -> SS(AB | A,B) \n   
Type III: SS(A | A,AB) -> SS(B | AB,B) -> SS(AB | A,B) # equivalent to linear regression  \n

# Appendix II: Examples for linear mixed-effect model
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

