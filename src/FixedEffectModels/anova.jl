# ================================================================================================
# Main API

using FixedEffectModels
@reexport using FixedEffectModels

"""
    lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...)

An alias for `reg`
"""
lfe(formula::FormulaTerm, df, vcov::CovarianceEstimator = Vcov.simple(); kwargs...) = 
    reg(df, formula, vcov, kwargs...)