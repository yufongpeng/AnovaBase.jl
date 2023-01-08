# AnovaBase.jl
```@meta
CurrentModule = AnovaBase
```

## Models
```@docs
AnovaBase.AnovaModel
AnovaBase.FullModel
AnovaBase.NestedModels
AnovaBase.nestedmodels(::RegressionModel)
```

## ANOVA
```@docs
AnovaBase.AnovaResult
AnovaBase.anova(::Type{<: GoodnessOfFit}, ::RegressionModel)
```

## Attributes
```@docs
AnovaBase.anova_test(aov::AnovaResult)
AnovaBase.anova_type(aov::AnovaResult)
AnovaBase.pval(aov::AnovaResult)
AnovaBase.teststat(aov::AnovaResult)
AnovaBase.deviance(aov::AnovaResult)
AnovaBase.dof(::AnovaResult)
AnovaBase.dof_residual(aov::AnovaResult)
AnovaBase.nobs(aov::AnovaResult)
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
AnovaBase.dof_asgn
AnovaBase.predictors(::RegressionModel)
AnovaBase.getterms
AnovaBase.isinteract
AnovaBase.select_super_interaction
AnovaBase.select_sub_interaction
AnovaBase.select_not_super_interaction
AnovaBase.select_not_sub_interaction
AnovaBase.subformula
AnovaBase.clear_schema
AnovaBase.extract_contrasts
AnovaBase._diff
AnovaBase._diffn
```

## IO interface
```@docs
AnovaBase.anovatable(::AnovaResult{<: RegressionModel})
AnovaBase.AnovaTable
AnovaBase.prednames
AnovaBase.testname
```