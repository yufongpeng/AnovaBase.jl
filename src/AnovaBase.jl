module AnovaBase

using Statistics, Distributions, Reexport, Printf
@reexport using StatsModels
import StatsBase: fit!, fit, dof, dof_residual, deviance, nobs, vcov
import StatsModels: TableRegressionModel, vectorize, collect_matrix_terms, coefnames, formula
import Base: show

export
    # Wrappers
    AnovaResult, 

    # anova functions
    anova, nestedmodels, 

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, canonicalgoodnessoffit, 

    # Attributes
    dof, dof_residual, deviance, nobs, formula,
    teststat, pval, anova_test, anova_type, 

    # Utils
    ftest_nested, lrt_nested, dof_asgn, _diff, _diffn,
    getterms, isinteract, 
    select_super_interaction, select_sub_interaction, 
    select_not_super_interaction, select_not_sub_interaction, 
    subformula, extract_contrasts, clear_schema, 

    # IO
    prednames, testname, anovatable, add_prednames!

# Test 
"""
    abstract type GoodnessOfFit end

An abstract type as super type of goodness of fit.
"""
abstract type GoodnessOfFit end
"""
    struct FTest <: GoodnessOfFit end

Type indicates conducting ANOVA by F-test. It can be the first argument or keyword argument `test`.
"""
struct FTest <: GoodnessOfFit end

const doc_lrt = """
    struct LikelihoodRatioTest <: GoodnessOfFit end
    const LRT = LikelihoodRatioTest

Type indicates conducting ANOVA by likelihood-ratio test. It can be the first argument or keyword argument `test`.
"""
@doc doc_lrt
struct LikelihoodRatioTest <: GoodnessOfFit end
@doc doc_lrt
const LRT = LikelihoodRatioTest

# Wrapper for ANOVA
"""
    AnovaResult{M, T, N}

Returned object of `anova`.
* `M` is a subtype of `Tuple` if multiple models are provided; otherwise, a type of statistical model.
* `T` is a subtype of `GoodnessOfFit`; either `FTest` or `LRT`.
* `N` is the length of parameters.

## Fields

* `model`: full model or a tuple of tested models.
* `type`: type of `anova`.
* `dof`: degrees of freedom of models or predictors.
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

function_arg_error(fn, type) = ErrorException("Arguments are valid for $fn; however, no method match $fn(::$type)")
function_arg_error(fn, type::AbstractString) = ErrorException("Arguments are valid for $fn; however, no method match $fn($type)")

include("fit.jl")
include("api.jl")
include("term.jl")
include("io.jl")

end
