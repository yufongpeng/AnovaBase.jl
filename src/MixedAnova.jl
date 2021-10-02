module MixedAnova
using DataFrames, Statistics, StatsBase, LinearAlgebra, Distributions, Reexport, Printf
@reexport using StatsModels
import StatsBase: fit!, fit, dof, dof_residual, deviance, nobs, coefnames
import StatsModels: TableStatisticalModel, TableRegressionModel, vectorize, kron_insideout,
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
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, canonicalgoodnessoffit, 

    # Others
    lme, glme, lfe, nestedmodels, calcdof, formula, #getterms,
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

