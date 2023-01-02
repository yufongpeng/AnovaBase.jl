# AnovaGLM
```@setup glm
using AnovaGLM, RDatasets, DataFrames
iris = dataset("datasets", "iris")
quine = dataset("MASS", "quine")
mtcars = dataset("datasets", "mtcars")
```
To use `anova` on [`GLM objects`](https://juliastats.org/GLM.jl/stable/) , we need `AnovaGLM.jl`.
```@example glm
using AnovaGLM
```
This function will export all functions from [`GLM`](https://juliastats.org/GLM.jl/stable/) and related function in this package, including `anova`, `anova_lm`, `anova_glm`.
## Ordinary linear model
We first import the well-known [`iris` dataset](https://en.wikipedia.org/wiki/Iris_flower_data_set) from [`RDatasets`](https://github.com/JuliaStats/RDatasets.jl).
```julia
iris = dataset("datasets", "iris")
```
There's two way to perform ANOVA. `anova_lm` accepts a formula and data like [`GLM.lm`](https://juliastats.org/GLM.jl/stable/api/#GLM.lm).
```@example glm
anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris)
```
We can specify the type of sum of squares by keyword argument `type`. Let's use type II SS.
```@example glm
anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris, type = 2)
```
A [`StatsModels.TableRegressionModel`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.TableRegressionModel) object is fitted and stored in the output of `anova_lm`.  

We can fit a model first and call `anova` instead. `anova` store the model as well.
!!! warn
    It doesn't create a copy, so any in-place change of the original model should be noticed. 
```@example glm
lm1 = lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris)
anova(lm1, type = 3)
```
Multiple models can be compared through the same function.  
!!! note
    The checker for nested models is not implemented now, so it should be ensured that the later model is more complex than the previous one.  
```@example glm
lms = nestedmodels(LinearModel, @formula(SepalLength ~ SepalWidth * Species), iris, dropcollinear = false)
anova(lms...)
```
This result is a little bit different from [`GLM.ftest`](https://juliastats.org/GLM.jl/stable/api/#GLM.ftest):
```@example glm
ftest(getproperty.(lms[2:end], :model)...)
```
In `anova`, the F value is calculated by dividing [MSR](https://en.wikipedia.org/wiki/Mean_squared_error) (mean of ΔDeviance) with mean of [RSS](https://en.wikipedia.org/wiki/Residual_sum_of_squares) of the most complex model just like [`anova` in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/anova), while in [`GLM.ftest`](https://juliastats.org/GLM.jl/stable/api/#GLM.ftest), the denominator is replaced by [RSS](https://en.wikipedia.org/wiki/Residual_sum_of_squares) of subsequent model.
## Generalized linear models 
```julia
quine = dataset("MASS", "quine")
```
We fit a [negative binomial regression](https://en.wikipedia.org/wiki/Generalized_linear_model) on [`quine` dataset](https://www.rdocumentation.org/packages/MASS/versions/7.3-57/topics/quine) from [`MASS`](https://www.rdocumentation.org/packages/MASS/versions/7.3-57).
```@example glm
nbm = glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink())
anova(nbm)
```
There's also `anova_glm` similar to `anova_lm`.  

`anova` will automatically select test from [F-test](https://en.wikipedia.org/wiki/F-test) or [likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test) depending on the type of [distribution](https://juliastats.org/GLM.jl/stable/#Fitting-GLM-models). For distribution of `Bernoulli()`, `Binomial()`, `Poisson()` that have fixed dispersion, likelihood-ratio test is preferred. For other distribution, F-test is preferred.  

The next one is an axample of [logistic regression](https://en.wikipedia.org/wiki/Logistic_regression).
```julia
mtcars = dataset("datasets", "mtcars")
```
We want to predict if the `AM` is 0 or 1. Let's use logistic regression with and without interaction terms, and compare this two models by likelihood-ratio test. 
```@example glm
glm1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink())
glm2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), LogitLink())
anova(glm1, glm2)
```
```@example glm
lrtest(glm1, glm2)
```
This function works identically as [`StatsModels.lrtest`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.lrtest).
!!! note
    We can also specify test by keyword arguments `test` or putting test in the first argument.

## Algorithm
There are vectors of models and the corresponding base models:
```math
\begin{aligned}
    \mathbf{M} &= (M_1, ..., M_n)\\\\
    \mathbf{B} &= (B_1, ..., B_n)
\end{aligned}
```
When $m$ models are given, $\mathbf{M} = (M_2, ..., M_m)$, $\mathbf{B} = (M_1, ..., M_{m-1})$. 

When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $\mathbf M$ and $\mathbf B$ depends on the type of ANOVA.

For the most complex model $M_n$, each factors are assigned a natural number $f_k$ sequentially where $f_l$ is the last factor.

Let the number of columns of model matrix of $M_n$ $m$.

The included factors of $M_j$ and $B_j$ are
```math
\begin{aligned}
    MF_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } M_i\}\\\\
    BF_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } B_i\}
\end{aligned}
```
A map $id_X: [1, m] \mapsto [1, f_l]$ maps the index of columns into the corresponding factors.

We can define a vector of index sets for each model:
```math
\mathbf{I} = (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) \in MF_k\setminus BF_k$.

The deviances are:
```math
\begin{aligned}
    \mathcal{D} &= (\mathcal{D}_1, ..., \mathcal{D}_n)\\\\
    \mathcal{R} &= (\mathcal{R}_1, ..., \mathcal{R}_n)
\end{aligned}
```
where $\mathcal{D}_i$ and $\mathcal{R}_i$ are the sum of [squared deviance residuals (unit deviance)](https://en.wikipedia.org/wiki/Deviance_(statistics)) of $M_i$ and $B_i$. 
It is equivalent to the residual sum of squares for ordinary linear regression.

The difference of $\mathcal{D}$ and $\mathcal{R}$ is:
```math
\boldsymbol{\Delta} \mathcal{D} = \mathcal{D} - \mathcal{R}
```
The degree of freedom is:
```math
\mathbf{dof} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$

The $\sigma$ is the estimated dispersion (or scale) parameter for $M_n$'s distribution

For ordinary linear regression, 
```math
\sigma^2 =\frac{\text{ residual sum of squares}}{\text{dof(residuals)}}
```

### F-test
F-value is a vector:
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{dof}, dof_{res}}
```
where 
```math
F_i = \frac{\Delta \mathcal{D}_i}{\sigma^2 \times dof_i}
```
and $dof_{res}$ is $dof(\text{residuals of } B_n)$

For a single model, F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the most complex model; the deviance is calculated backward. Each $M_j$ corresponds to a factor $f_j$, i.e. $id_X[I_j] = \{f_j\}$.
#### Type I

Factors are sequentially added to the models, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$

Calculate the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiply with $\boldsymbol \beta$: 
    
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \mathbf{f} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{f_k^2}}{dof_j}
\end{aligned}
```

#### Type II
    
The included facrors are defined as follows:
```math
\begin{aligned}
    BF_j &= \{f_k \in [1, f_l]\, |\, f_k \text{ is not an interaction term of }f_j \text{ and other terms}\}\\\\
    MF_j &= BF_j \cup \{f_j\}
\end{aligned}
```
Define two vectors of index sets $\mathbf J$ and $\mathbf K$ where 
```math
\begin{aligned}
    J_j &= \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }f_j \text{ and other terms}\}\\\\
    K_j &= J_j \cup I_j
\end{aligned}
```
F-value is: 
```math
F_j = \frac{\boldsymbol{\beta}_{K_j}^T \boldsymbol{\Sigma}_{K_j; K_j}^{-1} \boldsymbol{\beta}_{K_j} - \boldsymbol{\beta}_{J_j}^T \boldsymbol{\Sigma}_{J_j; J_j}^{-1} \boldsymbol{\beta}_{J_j}}{dof_j}
```

#### Type III

The models are all $M_n$, the base models are models without each factors.  

F-value is:
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{dof_j}
```

### LRT
The likelihood ratio is a vector:
```math
\begin{aligned} 
    \mathbf{LR} &= \boldsymbol{\Delta} \mathcal{D}/\sigma^2
    \mathbf{LR} &\sim \chi^2_{\mathbf{dof}}
\end{aligned}
```
When a single model is provided, $\mathbf{LR}$ is computed directly by the variance-covariance matrix.

Calculate the the upper factor of Cholesky factorization of $\sigma^2 \boldsymbol{\Sigma}^{-1}$ and multiply with $\boldsymbol \beta$:
```math
\begin{aligned}
    \sigma^2 \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \mathbf{f} &= \mathbf{U}\boldsymbol{\beta}\\\\
    LR_j &= \sum_{k \in I_j}{f_k^2}
\end{aligned}
```