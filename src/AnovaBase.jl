module AnovaBase

using Statistics, Distributions, Reexport, Printf
@reexport using StatsModels
using StatsBase: fit!, fit, dof, dof_residual, deviance, vcov, coeftable
import StatsAPI: dof_residual, deviance, nobs, dof
using StatsModels: TableRegressionModel, vectorize, collect_matrix_terms, coefnames, formula, asgn, TupleTerm, hasintercept
import Base: show

export
    # Wrappers
    AnovaResult, AnovaModel, NestedModels, MixedAovModels, FullModel, MultiAovModels,

    # Main function
    anova, nestedmodels, 

    # GoodnessOfFit
    GoodnessOfFit, FTest, LikelihoodRatioTest, LRT, ScaledLikelihoodRatioTest, SLRT,

    # Attributes
    dof, dof_residual, deviance, nobs, formula,
    teststat, pvalue, anova_test, anova_type,

    # IO
    prednames, anovatable, AnovaTable

# Test 
"""
    abstract type GoodnessOfFit end

An abstract type as super type of all types of goodness of fit.
"""
abstract type GoodnessOfFit end
"""
    struct FTest <: GoodnessOfFit end

Type indicates conducting ANOVA by F-test. It can be the first argument or keyword argument `test` in `anova` function.
"""
struct FTest <: GoodnessOfFit end

const doc_lrt = """
    struct LikelihoodRatioTest <: GoodnessOfFit end
    const LRT = LikelihoodRatioTest

Type indicates conducting ANOVA by likelihood-ratio test. It can be the first argument or keyword argument `test` in `anova` function.
"""
@doc doc_lrt
struct LikelihoodRatioTest <: GoodnessOfFit end
@doc doc_lrt
const LRT = LikelihoodRatioTest

const doc_slrt = """
    struct ScaledLikelihoodRatioTest <: GoodnessOfFit end
    const SLRT = ScaledLikelihoodRatioTest

Type indicates conducting ANOVA by scaled likelihood-ratio test. It can be the first argument or keyword argument `test` in `anova` function.
"""
@doc doc_slrt
struct ScaledLikelihoodRatioTest <: GoodnessOfFit end
@doc doc_slrt
const SLRT = ScaledLikelihoodRatioTest

include("term.jl")

"""
    abstract type AnovaModel{M, N} end

An abstract type as super type of any models for ANOVA.
"""
abstract type AnovaModel{M, N} <: StatisticalModel end
"""
    NestedModels{M, N} <: AnovaModel{M, N}

A wrapper of nested models of the same type for conducting ANOVA.
* `M` is a type of regression model.
* `N` is the number of models.

# Fields
* `model`: a tuple of models.

# Constructors
    NestedModels(model::Vararg{M, N}) where {M, N}
    NestedModels(model::NTuple{N, M}) where {M, N}
"""
struct NestedModels{M, N} <: AnovaModel{M, N}
    model::NTuple{N, M}
end
NestedModels(model::Vararg{M, N}) where {M, N} = NestedModels(model)
NestedModels(model::T) where {T <: Tuple} = throw(ArgumentError("`NestedModels` only accept models of the same type; use `MixedAovModels` instead."))
NestedModels(model...) = throw(ArgumentError("`NestedModels` only accept models of the same type; use `MixedAovModels` instead."))
"""
    MixedAovModels{M, N} <: AnovaModel{M, N}

A wrapper of nested models of multiple types for conducting ANOVA.
* `M` is a union type of regression models.
* `N` is the number of models.

# Fields
* `model`: a tuple of models.

# Constructors
    MixedAovModels{M}(model...) where M 
    MixedAovModels{M}(model::T) where {M, T <: Tuple}
"""
struct MixedAovModels{M, N} <: AnovaModel{M, N}
    model::Tuple
end
MixedAovModels(model...) = MixedAovModels{Union{typeof.(model)...}, length(model)}(model)
MixedAovModels(model::Vararg{M, N}) where {M, N} = throw(ArgumentError("`MixedAovModels` only accept models of different types; use `NestedModels` instead."))
MixedAovModels(model::T) where {T <: Tuple} = MixedAovModels{Union{T.parameters...}, length(model)}(model)
MixedAovModels(model::NTuple{N, M}) where {M, N} = throw(ArgumentError("`MixedAovModels` only accept models of different types; use `NestedModels` instead."))
"""
    const MultiAovModels{M, N} = Union{NestedModels{M, N}, MixedAovModels{M, N}} where {M, N}

Wrappers of mutiple models.
"""
const MultiAovModels{M, N} = Union{NestedModels{M, N}, MixedAovModels{M, N}} where {M, N}
"""
    MultiAovModels(model::NTuple{N, M}) where {M, N} -> NestedModels{M, N}
    MultiAovModels(model::Vararg{M, N}) where {M, N} -> NestedModels{M, N}
    MultiAovModels(model::T) where {T <: Tuple}      -> MixedAovModels
    MultiAovModels(model...)                         -> MixedAovModels

Construct `NestedModels` or `MixedAovModels` based on model types.
"""
MultiAovModels(model::NTuple{N, M}) where {M, N} = NestedModels(model)
MultiAovModels(model::Vararg{M, N}) where {M, N} = NestedModels(model)
MultiAovModels(model::T) where {T <: Tuple} = MixedAovModels(model)
MultiAovModels(model...) = MixedAovModels(model)

"""
    FullModel{M, N} <: AnovaModel{M, N}

A wrapper of a regression model for conducting ANOVA.
* `M` is a type of regression model.
* `N` is the number of predictors.

# Fields
* `model`: a regression model.
* `pred_id`: the index of terms included in ANOVA. The source iterable can be obtained by `predictors(model)`. This value may depend on `type` for certain model, e.g. type 1 ANOVA for a gamma regression model with inverse link.
* `type`: type of ANOVA, either 1, 2 or 3.

# Constructor
    FullModel(model::RegressionModel, type::Int, null::Bool, test_intercept::Bool)

* `model`: a regression model.
* `type`: type of ANOVA, either 1, 2 or 3.
* `null`: whether `y ~ 0` is allowed.
* `test_intercept`: whether intercept is going to be tested.
"""
struct FullModel{M, N} <: AnovaModel{M, N}
    model::M
    pred_id::NTuple{N, Int}
    type::Int
end
function FullModel(model::RegressionModel, type::Int, null::Bool, test_intercept::Bool)
    err1 = ArgumentError("Invalid set of model specification for ANOVA; not enough variables provided.")
    #err2 = ArgumentError("Invalid set of model specification for ANOVA; all coefficents are aliased with 1.")
    preds = predictors(model)::TupleTerm
    pred_id = collect(eachindex(preds))
    hasintercept(preds) || popfirst!(pred_id)
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
* `T` is a subtype of `GoodnessOfFit`; `FTest`, `LRT`, or `SLRT`.
* `N` is the length of parameters.

# Fields
* `anovamodel`: [`NestedModels`](@ref), [`MixedAovModels`](@ref), or [`FullModel`](@ref).
* `dof`: degrees of freedom of models or predictors.
* `deviance`: deviance(s) for calculating test statistics. See [`deviance`](@ref) for more details.
* `teststat`: value(s) of test statiscics.
* `pvalue`: p-value(s) of test statiscics.
* `otherstat`: `NamedTuple` contained extra statistics.

# Constructor
    AnovaResult(
            anovamodel::M,
            ::Type{T},
            dof::NTuple{N, Int},
            deviance::NTuple{N, Float64},
            teststat::NTuple{N, Float64},
            pvalue::NTuple{N, Float64},
            otherstat::NamedTuple
    ) where {N, M <: AnovaModel{<: RegressionModel, N}, T <: GoodnessOfFit}
"""
struct AnovaResult{M, T, N}
    anovamodel::M
    dof::NTuple{N, Int}
    deviance::NTuple{N, Float64}
    teststat::NTuple{N, Float64}
    pvalue::NTuple{N, Float64}
    otherstat::NamedTuple
end

AnovaResult(
    anovamodel::M,
    ::Type{T},
    dof::NTuple{N, Int},
    deviance::NTuple{N, Float64},
    teststat::NTuple{N, Float64},
    pvalue::NTuple{N, Float64},
    otherstat::NamedTuple
) where {N, M <: AnovaModel{<: RegressionModel, N}, T <: GoodnessOfFit} = 
    AnovaResult{M, T, N}(anovamodel, dof, deviance, teststat, pvalue, otherstat)

function_arg_error(fn, type) = ErrorException("Arguments are valid for $fn; however, no method match $fn(::$type)")
function_arg_error(fn, type::AbstractString) = ErrorException("Arguments are valid for $fn; however, no method match $fn($type)")

include("fit.jl")
include("attr.jl")
include("io.jl")
include("interface.jl")

end
