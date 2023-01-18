module AnovaBase

using Statistics, Distributions, Reexport, Printf
@reexport using StatsModels
import StatsBase: fit!, fit, dof, dof_residual, deviance, nobs, vcov, coeftable
import StatsModels: TableRegressionModel, vectorize, collect_matrix_terms, coefnames, formula, asgn, TupleTerm
import Base: show

export
    # Wrappers
    AnovaResult, AnovaModel, NestedModels, FullModel,

    # Main function
    anova, nestedmodels, 

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT,

    # Attributes
    dof, dof_residual, deviance, nobs, formula,
    teststat, pval, anova_test, anova_type,

    # IO
    prednames, anovatable, AnovaTable

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

include("term.jl")

"""
    abstract type AnovaModel{M, N} end

An abstract type as super type of any models for ANOVA.
"""
abstract type AnovaModel{M, N} <: StatisticalModel end
"""
    NestedModels{M, N} <: AnovaModel{M, N}

A wrapper of nested models for conducting ANOVA.
* `M` is a type of regression model.
* `N` is the number of models.

# Fields
* `models`: a tuple of models.
"""
struct NestedModels{M, N} <: AnovaModel{M, N}
    model::Tuple

    NestedModels{M}(model::T) where {M, N, T <: NTuple{N, M}} = new{M, N}(model)
end

NestedModels{M}(model...) where M = NestedModels{M}(model)
NestedModels{M}(model::T...) where {M, T <: Tuple} = throw(ArgumentError("Some models in $T are not subtype of $M"))

"""
    FullModel{M, N} <: AnovaModel{M, N}

A wrapper of full model for conducting ANOVA.
* `M` is a type of regression model.
* `N` is the number of predictors.

# Fields
* `model`: a regression model.
* `pred_id`: the index of terms included in ANOVA. The source iterable can be obtained by `predictors(model)`. This value may depend on `type` for certain model, e.g. type 1 ANOVA for a gamma regression model with inverse link.
* `type`: type of ANOVA, either 1, 2 or 3.
"""
struct FullModel{M, N} <: AnovaModel{M, N}
    model::M
    pred_id::NTuple{N, Int}
    type::Int
end

"""
    FullModel(model::RegressionModel, type::Int, null::Bool, test_intercept::Bool)

Create a `FullModel` with several model-specific parameters.

* `model`: a regression model.
* `type`: type of ANOVA, either 1, 2 or 3.
* `null`: whether `y ~ 0` is allowed.
* `test_intercept`: whether intercept is going to be tested.
"""
function FullModel(model::RegressionModel, type::Int, null::Bool, test_intercept::Bool)
    err1 = ArgumentError("Invalid set of model specification for ANOVA; not enough variables provided.")
    #err2 = ArgumentError("Invalid set of model specification for ANOVA; all coefficents are aliased with 1.")
    preds = predictors(model)
    pred_id = collect(eachindex(preds))
    has_intercept(preds) || popfirst!(pred_id)
    isempty(pred_id) && throw(err1) # ~ 0
    if type ≡ 1
        # ~ 0 + A + B..., ~ 1 + B..., ~ B as null
        null || popfirst!(pred_id)
    elseif type ≡ 2
        # ~ 0 + A + A & B + A & ..., all terms are related to A, ~ A as null
        null || isempty(select_not_super_interaction(preds, first(pred_id))) && popfirst!(pred_id)
        # ~ 1 + A..., all terms are aliased with 1, ~ 1 as null
        null || any_not_aliased_with_1(preds) || filter!(!=(1), pred_id)
    elseif type ≡ 3
        # ~ 1 + A + B + C... when !null and all aliased with 1, ~ 1 as null
        null || any_not_aliased_with_1(preds) || filter!(!=(1), pred_id)
        # ~ 1
        null || length(pred_id) ≡ 1 && throw(err1)
        
    else
        throw(ArgumentError("Invalid type of ANOVA"))
    end
    test_intercept || filter!(!=(1), pred_id) # ~ 1
    isempty(pred_id) && throw(err1)
    FullModel(model, tuple(pred_id...), type)
end

# Wrapper for ANOVA
"""
    AnovaResult{M, T, N}

Returned object of `anova`.
* `M` is `NestedModels` or `FullModel`.
* `T` is a subtype of `GoodnessOfFit`; either `FTest` or `LRT`.
* `N` is the length of parameters.

# Fields
* `anovamodel`: a `NestedModels` or a `FullModel`.
* `dof`: degrees of freedom of models or predictors.
* `deviance`: deviance(s) for calculating test statistics. See [`deviance`](@ref) for more details.
* `teststat`: value(s) of test statiscics.
* `pval`: p-value(s) of test statiscics.
* `otherstat`: `NamedTuple` contained extra statistics.
"""
struct AnovaResult{M, T, N}
    anovamodel::M
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    teststat::NTuple{N, Float64}
    pval::NTuple{N, Float64}
    otherstat::NamedTuple
end

AnovaResult{T}(anovamodel::M,
            dof::NTuple{N, Int},
            deviance::NTuple{N, Float64},
            teststat::NTuple{N, Float64},
            pval::NTuple{N, Float64},
            otherstat::NamedTuple) where {N, M <: AnovaModel{<: RegressionModel, N}, T <: GoodnessOfFit} = 
    AnovaResult{M, T, N}(anovamodel, dof, deviance, teststat, pval, otherstat)

function_arg_error(fn, type) = ErrorException("Arguments are valid for $fn; however, no method match $fn(::$type)")
function_arg_error(fn, type::AbstractString) = ErrorException("Arguments are valid for $fn; however, no method match $fn($type)")

include("fit.jl")
include("attr.jl")
include("io.jl")
include("interface.jl")

end
