module MixedAnova
using DataFrames, Statistics, StatsBase, LinearAlgebra, Distributions, Reexport, Printf
@reexport using StatsModels
import StatsBase: fit!, fit, dof, dof_residual, deviance, nobs, coefnames
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout ,
                    ModelFrame, ModelMatrix, response, columntable, asgn, collect_matrix_terms
import LinearAlgebra.BlasReal
import Tables.istable
import Base: show

export
    # wrappers
    AnovaResult, 

    # anova functions
    anova, anova_lm, anova_lme, anova_glm,

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, #canonicalgoodnessoffit, 

    # Others
    lme, lfe, nestedmodels, calcdof, formula, #getterms, getterm,
    teststat, pval, anova_test, anova_type,

    # init
    glm_init, mm_init, fem_init


# Test 
abstract type GoodnessOfFit end
struct FTest <: GoodnessOfFit end
struct LikelihoodRatioTest <: GoodnessOfFit end
const LRT = LikelihoodRatioTest

# Wrapper for ANOVA
"""
    AnovaResult{M, T, N}

Returned object of `anova`. \n
`M` is a subtype of `Tuple` if multiple models are provided; otherwise, a typeof model. \n
`T` is a subtype of `GoodnessOfFit`; either `FTest` or `LRT`. \n
`N` is the length of parameters.

## Fields

* `model`: full model or tuple of tested models.
* `type`: type of `anova`.
* `dof`: degree of freedoms of models or factors.
* `deviance`: deviance(s) for calculating test statistics. See `deviance` for more details.
* `teststat`: value(s) of test statiscics.
* `pval`: p-value(s) of test statiscics.
* `tests`: `NamedTuple` contained extra statistics.
"""
struct AnovaResult{M, T, N}
    model::M
    type::Int
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    teststat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
    tests::NamedTuple
end

AnovaResult{T}(model::M,
                type::Int,
                dof::NTuple{N, Int},
                deviance::NTuple{N, Float64},
                teststat::NTuple{N, Float64},
                pval::NTuple{N, Float64},
                tests::NamedTuple) where {M, N, T <: GoodnessOfFit} = 
    AnovaResult{M, T, N}(model, type, dof, deviance, teststat, pval, tests)

# Use tuple/svector instead vector?
#=
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

For other fields, please see `FixedAnovaStatsF`.
"""
mutable struct MixedAnovaStatsF{M, N} <: AbstractAnovaStats
    type::Int
    nobs::Int
    dof::NTuple{N, Int}
    resdof::NTuple{N, Int}
    fstat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
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


#mutable struct NestedAnovaStatsRao <: NestedAnovaStats end

#mutable struct NestedAnovaStatsCp <: NestedAnovaStats end

const AnovaStatsF = Union{FixedAnovaStatsF, MixedAnovaStatsF, NestedAnovaStatsF}
const AnovaStatsLRT = Union{FixedAnovaStatsLRT, MixedAnovaStatsLRT, NestedAnovaStatsLRT}
#const AnovaStatsRao = Union{NestedAnovaStatsRao}
#const AnovaStatsCp = Union{NestedAnovaStatsCp}
=#

# Calculate dof from assign
function dof(v::Vector{Int})
    dofv = zeros(Int, v[end])
    prev = 1
    ind = 1
    n = length(v)
    while ind <= n
        v[ind] == prev || (prev = v[ind])
        dofv[prev] += 1
        ind += 1
    end
    Int.(dofv)
end

_diff(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i + 1] - t[i], N - 1)
_diffn(t::NTuple{N, T}) where {N, T} = ntuple(i->t[i] - t[i + 1], N - 1)

include("api.jl")
include("termIO.jl")
const path = @__DIR__()

function glm_init() 
    include(joinpath(path, "GLM", "anova.jl"))
    include(joinpath(path, "GLM", "io.jl"))
    include(joinpath(path, "GLM", "fit.jl"))
    return
end

function mm_init()
    include(joinpath(path, "MixedModels", "anova.jl"))
    include(joinpath(path, "MixedModels", "io.jl"))
    include(joinpath(path, "MixedModels", "fit.jl"))
    return
end

function fem_init()
    include(joinpath(path, "FixedEffectModels", "anova.jl"))
    include(joinpath(path, "FixedEffectModels", "io.jl"))
    include(joinpath(path, "FixedEffectModels", "fit.jl"))
    return
end

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

