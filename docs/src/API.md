```@meta
CurrentModule = AnovaBase
```
# API
```@index
Modules = [AnovaBase]
Order   = [:type, :function]
Private = false
```

## Type
```@autodocs
Modules = [AnovaBase]
Order   = [:type]
Private = false
```

## ANOVA
```@docs
anova
anova_lm
anova_glm
anova_lme
anova_lfe
```
## Model fit
```@docs
nestedmodels
GLM.glm(::FormulaTerm, ::DataFrame, ::Binomial, ::Link, ::Vararg{Any})
lme
glme
lfe
to_trm
```
## Attributes
```@docs
formula
anova_test
anova_type
pval
teststat
coefnames
deviance
dof
dof_residual
nobs
```
## Miscellaneous
```@docs
canonicalgoodnessoffit
calcdof
```
