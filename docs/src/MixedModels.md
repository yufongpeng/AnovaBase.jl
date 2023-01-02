# AnovaMixedModels
```@setup mm
using AnovaMixedModels, DataFrames, GLM, CSV, CategoricalArrays
anxiety = CSV.read("anxiety.csv", DataFrame)
transform!(anxiety, :id => categorical, renamecols = false)
toenail = CSV.read("toenail.csv", DataFrame)
transform!(toenail, [1, 2, 3] .=> categorical, renamecols = false)
```
The implementation of ANOVA for [mixed-effects models](https://en.wikipedia.org/wiki/Mixed_model) is primarily based on [`MixedModels`](https://juliastats.org/MixedModels.jl/stable/). The syntax is similar to anova for `GLM`.   
```@example mm
using AnovaMixedModels
```
## Linear mixed-effects model
We get a dataset from `R` directly by [`RCall`](https://juliainterop.github.io/RCall.jl/stable/).
```julia
R"""data("anxiety", package = "datarium")"""
anxiety = stack(rcopy(R"anxiety"), [:t1, :t2, :t3], [:id, :group], variable_name = :time, value_name = :score)
anxiety = combine(anxiety, Not(:time), :time => ByRow(x->parse(Int, replace(String(x), "t"=>""))) => :time)
```
We can fit a linear mixed-effects model first. `lme` is an alias for [`fit(LinearMixedModel, formula, data, args...)`](https://juliastats.org/MixedModels.jl/stable/constructors/#Examples-of-linear-mixed-effects-model-fits).
```@example mm
lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety)
anova(lmm1)
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` fit by [Restricted maximum likelihood](https://en.wikipedia.org/wiki/Restricted_maximum_likelihood).
```@example mm
aov = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
```
```@example mm
aov.model.optsum.REML
```
!!! note
    Type 2 sum of squares is not implemented now.  
For [likeihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test), all submodels are fitted. The model should be fitted by [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation).
```@example mm
anova(LRT, lmm1)
```
When comparing multiple mixed models, likelihood-ratio test is used by default. 
It's also identical to [`StatsModels.lrtest`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.lrtest) and [`MixedModels.likelihoodratiotest`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LikelihoodRatioTest).
```@example mm
lmms = nestedmodels(lmm1)
anova(lmms...) # as same as anova(LRT, lmm1)
MixedModels.likelihoodratiotest(lmms[2:end]...)
``` 
Comparing between [`LinearModel`](https://juliastats.org/GLM.jl/stable/api/#GLM.LinearModel) and [`LinearMixedModel`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LinearMixedModel) is also available.
```@example mm
lm1 = lm(@formula(score ~ group * time), anxiety)
lmm2 = lme(@formula(score ~ group * time + (group|id)), anxiety)
anova(lm1, lmm1, lmm2)
```
## Generalized linear mixed-effects model
The following is an example of generalized mixed model. `glme` is an alias for [`fit(GeneralizedLinearMixedModel, formula, data, args...)`](https://juliastats.org/MixedModels.jl/stable/constructors/#Fitting-generalized-linear-mixed-models).
```julia
R"""data("toenail", package = "HSAUR2")"""
toenail = rcopy(R"toenail")
```
```@example mm
glmm1 = glme(@formula(outcome ~ visit + treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
glmm2 = glme(@formula(outcome ~ visit * treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
glm1 = glm(@formula(outcome ~ visit + treatment), toenail, Binomial(), LogitLink());
anova(glm1, glmm1, glmm2)
```
!!! note
    Only likelihood-ratio test is available now for `GeneralizedLinearMixedModel`.

## Algorithm
### F-test
Given a model $M$, $n$ is the number of factors. Each factors are assigned a natural number $f_k$ sequentially.

Let the number of columns of model matrix of $M$ $m$.

A map $id_X: [1, m] \mapsto [1, n]$ maps the index of columns into the corresponding factors.

We can define a vector of index set for each factors:
```math
\mathbf{I} = (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) = f_k$.

The degrees of freedom (dof) is:
```math
\mathbf{df} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$

F-value is a vector:
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, \mathbf{df_r}}
```
where $\mathbf df_r$ is estimated by between-within method.

F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the model. 
#### Type I

Factors are sequentially added to the models, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$

Calculate the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiply with $\boldsymbol \beta$:
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \mathbf{f} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{f_k^2}}{df_j}
\end{aligned}
```

#### Type III

```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```

### LRT
Given a vector of models:
```math
\mathbf{M} = (M_1, ..., M_n)
``` 
The $\mathcal{D}$ is a vector of loglikelihoods $-2loglikelihood(M_i)$ for a linear mixed-effect model (linear model); unit deviance for a generalized linear mixed-effect model (generalized linear model).

The likelihood ratio is a vector:
```math
\begin{aligned}
    \mathbf{LR} &= \mathcal{D}_{[1, n - 1]} - \mathcal{D}_{[2, n]}\\\\
    \mathbf{LR} &\sim \chi^2_{\mathbf{df}}
\end{aligned}
```
where $df_i = dof(M_i) - dof(M_{i+1})$