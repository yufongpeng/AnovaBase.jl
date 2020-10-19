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
Two-way anova:
Type I  : SS(A) -> SS(B | A) -> SS(AB | A,B)
Type II : SS(A | B) -> SS(B | A) -> SS(AB | A,B)
Type III: SS(A | A,AB), SS(B | AB,B), SS(AB | A,B) equivalent to linear regression

A, B: fixed effect -> denominator: SSR
A, B: random effect -> denominator: SS(AB | A,B)
A/ B: fixed effect/ random effect -> denominator: SS(AB | A,B)/ SSR
"""


end

