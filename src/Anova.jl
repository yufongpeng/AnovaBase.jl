__precompile__()

module Anova

using GLM, Statistics, StatsBase, StatsModels, LinearAlgebra, Distributions, Reexport
import GLM: LinPredModel, LinearModel, LmResp, DensePred,
            DensePredChol, QRCompactWY, LinPred, installbeta!, delbeta!,  linpred!,
            updateμ!, linpred, cholfactors, updateμ!
import StatsBase: fit!, fit
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable
import LinearAlgebra.BlasReal
import Tables.istable

export
    # models
    AnovaResult, AnovaStats,

    # functions
    anova

@reexport using GLM
@reexport using StatsModels

"""
    AnovaResult{TableRegressionModel,AnovaStats}

Returned object of `anova`.

`model` is the full model.
`stats` contains result of ANOVA, including sum of squares, fstats, etc.
"""
struct AnovaResult{TableRegressionModel,AnovaStats}
    model::TableRegressionModel
    stats::AnovaStats
end

"""
    AnovaStats

Object contains result of ANOVA

`type`: the type of ANOVA.
`nobs`: the number of observations.
`ss`: sum of squares for each predictors. 
    if `type` is 3, than the first predictor is the intercept; otherwise, it is the first variable.
`dof`: degre of freedom.
`fstat`: f statiscics for each predictor.
`pval`: p-values for each predictor.
"""
mutable struct AnovaStats
    type::Int
    nobs::Int
    ss::Vector{Float64}
    dof::Vector{Int}
    fstat::Vector{Float64}
    pval::Vector{Float64}
end

include("fixedeffect.jl")
include("mixedmodel.jl")


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

