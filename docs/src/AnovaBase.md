# AnovaBase.jl
```@meta
CurrentModule = AnovaBase
```

## Models
```@docs
AnovaBase.AnovaModel
AnovaBase.FullModel
AnovaBase.NestedModels
AnovaBase.MixedAovModels
AnovaBase.MultiAovModels
AnovaBase.nestedmodels(::RegressionModel)
```

## ANOVA
```@docs
AnovaBase.AnovaResult
AnovaBase.anova(::Type{<: GoodnessOfFit}, ::AnovaModel)
```

## Attributes
```@docs
AnovaBase.anova_test(aov::AnovaResult)
AnovaBase.anova_type(aov::AnovaResult)
AnovaBase.pval(aov::AnovaResult)
AnovaBase.teststat(aov::AnovaResult)
AnovaBase.deviance(aov::AnovaResult)
AnovaBase.dof(::AnovaResult)
AnovaBase.nobs(aov::AnovaResult)
```

## Goodness of fit
```@docs
AnovaBase.GoodnessOfFit
AnovaBase.FTest
AnovaBase.LikelihoodRatioTest
AnovaBase.canonicalgoodnessoffit
```

## Other interface
```@docs
AnovaBase.ftest_nested
AnovaBase.lrt_nested
AnovaBase.dof_residual(aov::AnovaResult)
AnovaBase.predictors(::RegressionModel)
AnovaBase.anovatable(::AnovaResult{<: FullModel})
```

## Developer utility
```@docs
AnovaBase.dof_asgn
AnovaBase.prednames
AnovaBase.has_intercept
AnovaBase.any_not_aliased_with_1
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
AnovaBase.AnovaTable
AnovaBase.testname
```