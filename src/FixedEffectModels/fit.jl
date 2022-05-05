# ==========================================================================================================
# Backend funcion

formula(trm::TableRegressionModel) = trm.mf.f

fe_intercept(f::FormulaTerm) = fe_intercept(f.rhs)
fe_intercept(term::StatsModels.TupleTerm) = map(fe_intercept, term)
fe_intercept(term::FunctionTerm) = first(term.exorig.args) == :fe 
fe_intercept(term) = false

width(term::FunctionTerm) = first(term.exorig.args) == :fe ? 0 : 1

# Variable dispersion
dof(model::FixedEffectModel) = model.nobs - model.dof_residual + 1