# AnovaGLM.jl
```@meta
CurrentModule = AnovaGLM
```

## ANOVA
```@docs
AnovaGLM.anova(::Val{:AnovaGLM})
anova_lm
anova_glm
```

## Models
```@docs
AnovaGLM.nestedmodels(::Val{:AnovaGLM})
GLM.glm(::FormulaTerm, ::DataFrame, ::Binomial, ::Link, ::Vararg{Any})
```
