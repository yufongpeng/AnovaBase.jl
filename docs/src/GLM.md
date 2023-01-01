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
Given vectors of models and the corresponding base models:
$$model = (model_1, ..., model_n)$$ 
$$basemodel = (basemodel_1, ..., basemodel_n)$$ 
When $m$ models are given, $model = (model_2, ..., model_m), basemodel = (model_1, ..., model_{m-1})$. 

When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $model$ and $basemodel$ depends on the type of ANOVA.

For the most complex model $model_n$, each factors are assigned a natural number sequentially.
$factor_i$ is the included factors of $model_i$; $basefactor_i$ is the included factors of $basemodel_i$. 

Let the number of columns of model matrix of $model_n$ $m$.

A map $id_X: [1, m] \mapsto factor_n$ maps the index of columns into the corresponding factors.

We can define a vector of index set for each model:
$$I = (I_1, ..., I_{n})$$
where $ \forall i \in I_k, id_X(i) \notin basefactor_k, id_X(i) \in factor_{k}$

The deviances are:
$$dev = (devres_1, ..., devres_n)$$ 
$$basedev = (basedevres_1, ..., basedevres_n)$$ 
where $devres_i$ and $basedevres_i$ are the sum of [squared deviance residuals (unit deviance)](https://en.wikipedia.org/wiki/Deviance_(statistics)) of $model_i$ and $basemodel_i$. 
It is equivalent to the residual sum of squares for ordinary linear regression.

The difference of $dev$ and $basedev$ is:
$$\Delta dev = basedev - dev$$
The degree of freedom is:
$$dof = (n(I_1), ..., n(I_{n}))$$ 
where $n(I)$ is the size of $I$

$dispersion$ is the estimated dispersion (or scale) parameter for $model_n$'s distribution.

For ordinary linear regression, 
$$dispersion^2 =\frac{\text{ residual sum of squares}}{\text{degree of freedom of residuals}}$$

### F-test
The attribute `deviance` of the returned object is $(\Delta dev_{1}, ..., \Delta dev_n, NaN)$.

F-value is a vector $F$ where 
$$ F_i = \frac{\Delta dev_i}{dispersion^2 \times dof_i}
$$

For a single model, F-value is computed directly by the variance-covariance matrix ($vcov$) and the coefficients ($\beta$) of the most complex model; the deviance is calculated backward. Each $model_l$ corresponds to a $factor_l$.
1. Type I:

    Factors are sequentially added to the models.
    $basemodel_l$ is $model_l$ without $factor_l$.

    First, calculate the the upper factor of Cholesky factorization of $vcov^{-1}$ and multiply with $\beta$:
    $$ vcov^{-1} = LU $$
    $$ f = U\beta $$
    F-value is a vector $F$ where 
    $$ F_l = \frac{\sum_{k \in I_l}{f_k^2}}{dof_l}
    $$
2. Type II:
    
    All included factors for $mode_l$ are other factors irrevalent to $factor_l$.
    The $basemodel_l$ is $model_l$ without $factor_l$. 

    A vector $J$ is a vector of index sets for each $basemodel_l$ where $J_l = \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }l \text{ and other terms}\}$;

    A vector $K$ is a vector of index sets for each $model_l$ where $K_l = \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }l \text{ and other terms or } id_X(i) = l\}$.

    F-value is a vector $F$ where 
    $$ F_l = \frac{(\beta_{K_l}^T  vcov_{K_l; K_l}^{-1} \beta_{K_l} - \beta_{J_l}^T  vcov_{J_l; J_l}^{-1} \beta_{J_l})}{dof_l}
    $$
3. Type III:

    The models are all the most comlex model $model_n$, the basemodels are models without each factors.  

    A vector of index sets $L$ where $L_l = \{i \in [1, m]\, |\, id_X(i) = l\}$.

    F-value is a vector $F$ where 
    $$ F_l = \frac{\beta_{L_l}^T  vcov_{L_l; L_l}^{-1} \beta_{L_l}}{dof_l}
    $$

## LRT
The attribute `deviance` of the returned object is is $dev$.

The likelihood ratio is defined as $Δdev / dispersion²$. 

When a single model is provided, lrt is computed directly by the variance-covariance matrix.

First, calculate the the upper factor of Cholesky factorization of $dispersion² \times vcov^{-1}$ and multiply with $\beta$:
$$ dispersion² \times vcov^{-1} = LU $$
$$ f = U\beta $$
The likelihood ratio is a vector $LR$ where $ LR_j = \sum_{k \in I_j}{f_k^2}$
