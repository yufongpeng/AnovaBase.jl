# AnovaGLM.jl

```@meta
CurrentModule = AnovaGLM
```

```@index
Modules = [AnovaGLM]
```

## ANOVA
```@docs
AnovaGLM.anova(::Val{:AnovaGLM})
anova_lm
anova_glm
```

## Model fit
```@docs
AnovaGLM.nestedmodels(::Val{:AnovaGLM})
GLM.glm(::FormulaTerm, ::DataFrame, ::Binomial, ::Link, ::Vararg{Any})
```
