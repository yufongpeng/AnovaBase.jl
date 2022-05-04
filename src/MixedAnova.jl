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
    anova, anova_lm, anova_lme, anova_glm, anova_lfe,

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, canonicalgoodnessoffit, 

    # Others
    lme, glme, lfe, nestedmodels, calcdof, formula, #getterms,
    teststat, pval, anova_test, anova_type, to_trm,

    # init
    glm_init, mm_init, fem_init

# Test 
"""
    abstract type GoodnessOfFit end
"""
abstract type GoodnessOfFit end
"""
    struct FTest <: GoodnessOfFit end

ANOVA by F-test. It can be the first argument or keyword argument `test`.
"""
struct FTest <: GoodnessOfFit end
"""
    struct LikelihoodRatioTest <: GoodnessOfFit end
    const LRT = LikelihoodRatioTest

ANOVA by Likelihood-ratio test. It can be the first argument or keyword argument `test`.
"""
struct LikelihoodRatioTest <: GoodnessOfFit end
"""
    const LRT = LikelihoodRatioTest

See `LikelihoodRatioTest`.
"""
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

end

