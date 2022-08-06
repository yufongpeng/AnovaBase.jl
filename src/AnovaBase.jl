module AnovaBase

using Statistics, Distributions, Reexport, Printf
@reexport using StatsModels
import StatsBase: fit!, fit, dof, dof_residual, deviance, nobs, vcov
import StatsModels: TableRegressionModel, vectorize, collect_matrix_terms, coefnames
import Base: show

export
    # wrappers
    AnovaResult, 

    # anova functions
    anova,

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, canonicalgoodnessoffit, 

    # Others
    nestedmodels, dof, dof_residual, deviance, nobs,
    teststat, pval, anova_test, anova_type

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

ANOVA by likelihood-ratio test. It can be the first argument or keyword argument `test`.
"""
struct LikelihoodRatioTest <: GoodnessOfFit end
"""
    const LRT = LikelihoodRatioTest

See [`LikelihoodRatioTest`](@ref).
"""
const LRT = LikelihoodRatioTest

# Wrapper for ANOVA
"""
    AnovaResult{M, T, N}

Returned object of `anova`.
* `M` is a subtype of `Tuple` if multiple models are provided; otherwise, a typeof model.
* `T` is a subtype of `GoodnessOfFit`; either `FTest` or `LRT`.
* `N` is the length of parameters.

## Fields

* `model`: full model or tuple of tested models.
* `type`: type of `anova`.
* `dof`: degree of freedoms of models or factors.
* `deviance`: deviance(s) for calculating test statistics. See [`deviance`](@ref) for more details.
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

include("fit.jl")
include("api.jl")
include("termIO.jl")

end
