# AnovaBase.jl
```@meta
CurrentModule = AnovaBase
```
```@index
Modules = [AnovaBase]
Order   = [:type, :function]
```
## Type
```@autodocs
Modules = [AnovaBase]
Order   = [:type]
Private = false
```

### ANOVA
```@docs
AnovaBase.anova
```
## Model fit
```@docs
AnovaBase.nestedmodels
```
## Attributes
```@docs
AnovaBase.formula
AnovaBase.anova_test
AnovaBase.anova_type
AnovaBase.pval
AnovaBase.teststat
AnovaBase.factornames
AnovaBase.deviance
AnovaBase.dof(::AnovaResult)
AnovaBase.dof_residual
AnovaBase.nobs
```

## Developer utility
```@docs
AnovaBase.ftest_nested
AnovaBase.lrt_nested
AnovaBase.dof_asgn(::Vector{Int})
AnovaBase.canonicalgoodnessoffit
AnovaBase.isinteract
AnovaBase.select_super_interaction
AnovaBase.select_not_super_interaction
AnovaBase.select_sub_interaction
AnovaBase.select_not_sub_interaction
AnovaBase.subformula
AnovaBase.getterms
AnovaBase.clearschema
AnovaBase.extract_contrasts
AnovaBase._diff
AnovaBase._diffn
```

## IO interface
```@docs
AnovaBase.anovatable
AnovaBase.TestStat
AnovaBase.PValue
AnovaBase.OtherStat
AnovaBase.NoQuote
```