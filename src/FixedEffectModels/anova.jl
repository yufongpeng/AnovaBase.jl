# ================================================================================================
# Main API

using FixedEffectModels
@reexport using FixedEffectModels
import StatsModels: width, apply_schema

anova(trms::Vararg{TableRegressionModel{<: FixedEffectModel}}; 
        test::Type{<: GoodnessOfFit} = FTest,
        kwargs...) = 
    anova(test, trms...; kwargs...)

# ================================================================================================
# ANOVA by F test

function anova(::Type{FTest}, 
            trm::TableRegressionModel{<: FixedEffectModel};
            type::Int = 1, kwargs...)

    type == 2           && throw(ArgumentError("Type 2 anova is not implemented"))
    type in [1, 2, 3]   || throw(ArgumentError("Invalid type"))
    assign = trm.mm.assign
    df = dof(assign)
    filter!(>(0), df)
    # May exist some floating point error from dof_residual
    push!(df, round(Int, dof_residual(trm)))
    df = tuple(df...)
    if type in [1, 3] 
        # vcov methods
        varβ = vcov(trm)
        β = trm.model.coef
        if type == 1
            fs = abs2.(cholesky(Hermitian(inv(varβ))).U * β) 
            offset = first(assign) - 1
            fstat = ntuple(last(assign) - offset) do fix
                sum(fs[findall(==(fix + offset), assign)]) / df[fix]
            end
        else
            # calculate block by block
            offset = first(assign) - 1
            fstat = ntuple(last(assign) - offset) do fix
                select = findall(==(fix + offset), assign)
                β[select]' * (varβ[select, select] \ β[select]) / df[fix]
            end
        end
        σ² = rss(trm.model) / last(df)
        devs = (fstat .* σ²..., σ²) .* df
    end
    pvalue = (ccdf.(FDist.(df[1:end - 1], last(df)), abs.(fstat))..., NaN)
    AnovaResult{FTest}(trm, type, df, devs, (fstat..., NaN), pvalue, NamedTuple())
end

# =================================================================================================================
# Nested models 

function anova(::Type{FTest}, 
                trms::Vararg{TableRegressionModel{<: FixedEffectModel}}; 
                check::Bool = true,
                isnested::Bool = false)

    df = dof.(trms)
    ord = sortperm(collect(df))
    df = df[ord]
    trms = trms[ord]

    # check comparable and nested
    check && @warn "Could not check whether models are nested: results may not be meaningful"

    Δdf = _diff(df)
    # May exist some floating point error from dof_residual
    dfr = round.(Int, dof_residual.(trms))
    dev = ntuple(length(trms)) do i 
        trms[i].model.rss
    end
    msr = _diffn(dev) ./Δdf
    σ² = last(dev) / last(dfr)
    fstat = msr ./ σ²
    pval = map(zip(Δdf, dfr[2:end], fstat)) do (dof, dofr, fs)
        fs > 0 ? ccdf(FDist(dof, dofr), fs) : NaN
    end
    AnovaResult{FTest}(trms, 1, df, dev, (NaN, fstat...), (NaN, pval...), NamedTuple())
end

"""
    lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...)

Fit a `FixedEffectModel` and wrap it into `TableRegressionModel`. 
!!! warn
    This function currently does not perform well. It re-compiles everytime; may be due to `@nonspecialize` for parameters of `reg`.
"""
lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...) = 
    to_trm(reg(df, formula, vcov; kwargs...), df)

"""
    to_trm(model, df)

Wrap fitted `FixedEffectModel` into `TableRegressionModel`.
"""
function to_trm(model::FixedEffectModel, df)
    f = model.formula
    has_fe_intercept = any(fe_intercept(f))
    rhs = vectorize(f.rhs)
    f = isa(first(rhs), ConstantTerm) ? f : FormulaTerm(f.lhs, (ConstantTerm(1), rhs...))
    s = schema(f, df, model.contrasts)
    f = apply_schema(f, s, FixedEffectModel, has_fe_intercept)
    mf = ModelFrame(f, s, Tables.columntable(df[!, getproperty.(keys(s), :sym)]), FixedEffectModel)
    # Fake modelmatrix
    assign = asgn(f)
    has_fe_intercept && popfirst!(assign)
    mm = ModelMatrix(ones(Float64, 1, 1), assign)
    TableRegressionModel(model, mf, mm)
end

# =================================================================================================================================
# Fit new models
"""
    anova_lfe(f::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); 
            test::Type{<: GoodnessOfFit} = FTest, <keyword arguments>)
    anova_lfe(test::Type{<: GoodnessOfFit}, f::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); <keyword arguments>)

ANOVA for fixed-effect linear regression.

* `type`: type of anova.

`anova_lfe` generate a `TableRegressionModel{<: FixedEffectModel}` object, which is fitted by `lfe`.
"""
anova_lfe(f::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); 
        test::Type{<: GoodnessOfFit} = FTest, 
        kwargs...)= 
    anova(test, FixedEffectModel, f, df, vcov; kwargs...)

anova_lfe(test::Type{<: GoodnessOfFit}, f::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...) = 
    anova(test, FixedEffectModel, f, df, vcov; kwargs...)

function anova(test::Type{<: GoodnessOfFit}, ::Type{FixedEffectModel}, f::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); 
        type::Int = 1, 
        kwargs...)
    trm = to_trm(reg(df, f, vcov; kwargs...), df)
    anova(test, trm; type)
end


