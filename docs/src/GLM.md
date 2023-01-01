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
    \bf{model} &= (model_1, ..., model_n)\\\\
    \bf{basemodel} &= (basemodel_1, ..., basemodel_n)
\end{aligned}
```
When $m$ models are given, $\bf model$ $= (model_2, ..., model_m)$, $\bf basemodel$ $ = (model_1, ..., model_{m-1})$. 

When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $\bf model$ and $\bf basemodel$ depends on the type of ANOVA.

For the most complex model $model_n$, each factors are assigned a natural number $f_k$ sequentially where $f_l$ is the last factor.

Let the number of columns of model matrix of $model_n$ $m$.

The included factors of $model_j$ and $basemodel_j$ are
```math
\begin{aligned}
    mf_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } model_i\}\\\\
    bf_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } basemodel_i\}
\end{aligned}
```
A map $id_X: [1, m] \mapsto [1, f_l]$ maps the index of columns into the corresponding factors.

We can define a vector of index sets for each model:
```math
\bf{I} = (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) \in mf_k\setminus bf_k$.

The deviances are:
```math
\begin{aligned}
    \bf{dev} &= (devres_1, ..., devres_n)\\\\
    \bf{basedev} &= (basedevres_1, ..., basedevres_n)
\end{aligned}
```
where $devres_i$ and $basedevres_i$ are the sum of [squared deviance residuals (unit deviance)](https://en.wikipedia.org/wiki/Deviance_(statistics)) of $model_i$ and $basemodel_i$. 
It is equivalent to the residual sum of squares for ordinary linear regression.

The difference of $\bf dev$ and $\bf basedev$ is:
```math
\bf{\Delta dev = basedev - dev}
```
The degree of freedom is:
```math
\bf{dof} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$

The $\sigma$ is the estimated dispersion (or scale) parameter for $model_n$'s distribution

For ordinary linear regression, 
```math
\sigma^2 =\frac{\text{ residual sum of squares}}{\text{degree of freedom of residuals}}
```

### F-test
F-value is a vector $\bf F$ where 
```math
\begin{aligned}
    F_i &= \frac{\Delta dev_i}{\sigma^2 \times dof_i}\\\\
    \bf{F} &\sim mathcal{F}_{\bf{dof}, dof_{res}}
\end{aligned}
```
where $dof_{res}$ is $dof(\text{residuals of } basemodel_n)$

For a single model, F-value is computed directly by the variance-covariance matrix ($\bf \Sigma$) and the coefficients ($\beta$) of the most complex model; the deviance is calculated backward. Each $model_j$ corresponds to a factor $f_j$, i.e. $id_X[I_j] = \{f_j\}$.
1. Type I:

    Factors are sequentially added to the models, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$

    Calculate the the upper factor of Cholesky factorization of $\bf \Sigma^{-1}$ and multiply with $\beta$:
    ```math
    \begin{aligned}
        \bf{\Sigma}^{-1} &= \bf{LU}\\\\
        \bf{f} &= \bf{U}\beta\\\\
        F_j &= \frac{\sum_{k \in I_j}{f_k^2}}{dof_j}
    \end{aligned}
    ```

2. Type II:
    
    For each $j$, $bf_j$ includes other factors irrevalent to $f_j$, i.e. 
    ```math
    bf_j = \{f_k in [1, f_l]\, |\, f_k \text{ is not an interaction term of }f_j \text{ and other terms}\}
    ```
    Define two vectors of index sets $J$ and $K$ where 
    ```math
    \begin{aligned}
        J_j &= \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }f_j \text{ and other terms}\}
        K_j &= \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }f_j \text{ and other terms or } id_X(i) = f_j\}
    \end{aligned}
    ```
    F-value is: 
    ```math
    F_j = \frac{(\beta_{K_j}^T \bf{\Sigma}_{K_j; K_j}^{-1} \beta_{K_j} - \beta_{J_j}^T \bf{\Sigma}_{J_j; J_j}^{-1} \beta_{J_j})}{dof_j}
    ```

3. Type III:

    The models are all $model_n$, the base models are models without each factors.  

    F-value is:
    ```math
    F_j = \frac{\beta_{I_j}^T \bf{\Sigma}_{I_j; I_j}^{-1} \beta_{I_j}}{dof_j}
    ```

## LRT
The likelihood ratio is defined as $\bf \Delta dev$ $/ \sigma^2$. 

When a single model is provided, lrt is computed directly by the variance-covariance matrix.

First, calculate the the upper factor of Cholesky factorization of $\sigma^2 \bf \Sigma^{-1}$ and multiply with $\beta$:
```math
\begin{aligned}
    \sigma^2 \bf{\Sigma}^{-1} &= \bf{LU}\\\\
    \bf{f} &= \bf{U}\beta
\end{aligned}
```
The likelihood ratio is a vector $\bf LR$ where 
```math
\begin{aligned}
    LR_j = &\sum_{k \in I_j}{f_k^2}\\\\
    \bf{LR} &\sim \chi^2_{\bf{dof}}
\end{aligned}
```
