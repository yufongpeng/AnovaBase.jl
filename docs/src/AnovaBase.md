# AnovaBase.jl
```@meta
CurrentModule = AnovaBase
```
```@index
Modules = [AnovaBase]
Order   = [:type, :function]
```
```@index
Modules = [StatsModels]
Order   = [:function]
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
StatsModels.formula
AnovaBase.anova_test
AnovaBase.anova_type
AnovaBase.pval
AnovaBase.teststat
AnovaBase.deviance
AnovaBase.dof(::AnovaResult)
AnovaBase.dof_residual
AnovaBase.nobs
```

## Goodness of fit
```@docs
AnovaBase.GoodnessOfFit
AnovaBase.FTest
AnovaBase.LikelihoodRatioTest
AnovaBase.canonicalgoodnessoffit
```

## Developer utility
```@docs
AnovaBase.ftest_nested
AnovaBase.lrt_nested
AnovaBase.dof_asgn(::Vector{Int})
AnovaBase.getterms
AnovaBase.isinteract
AnovaBase.select_super_interaction
AnovaBase.subformula
AnovaBase.clear_schema
AnovaBase.extract_contrasts
AnovaBase._diff
AnovaBase._diffn
```

## IO interface
```@docs
AnovaBase.anovatable
AnovaBase.add_prednames!
AnovaBase.prednames
AnovaBase.testname
```