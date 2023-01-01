# AnovaFixedEffectModels
```@setup fem
using AnovaFixedEffectModels, DataFrames, CSV, CategoricalArrays
gpa = CSV.read("gpa.csv", DataFrame)
transform!(gpa,
        7 => x->replace(x, "yes" => true, "no" => false, "NA" => missing),
        4 => x->categorical(x, levels = ["1 hour", "2 hours", "3 hours"], ordered = true),
        renamecols = false)
transform!(gpa, [1, 2, 5, 7] .=> categorical, renamecols = false)
```
```@example fem
using AnovaFixedEffectModels
```
`AnovaFixedEffectModels.jl` supports [`FixedEffectModels`](https://github.com/FixedEffects/FixedEffectModels.jl); however, because `anova` relies on model schema, the output of `FixedEffectModels.reg` is not compatible. 

To solve this issue, fitting model using `lfe` instead of `reg`.
```@example fem
fem1 = lfe(@formula(gpa ~ fe(student) + occasion + job), gpa)
```
If a model is already fitted by `reg`, use `to_trm` to convert it into [`StatsModels.TableRegressionModel`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.TableRegressionModel).
```@example fem
model = reg(gpa, @formula(gpa ~ fe(student) + occasion + job))
fem1 = to_trm(model, gpa)
aovf = anova(fem1)
```
!!! note
    `lfe` is actually slower because it re-compiles every execution.
!!! note 
    Only F-test is available for `FixedEffectModel`.

## Algorithm
Given vectors of models and the corresponding base models:
$$model = (model_1, ..., model_n)$$ 
$$basemodel = (basemodel_1, ..., basemodel_n)$$ 
When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $model$ and $basemodel$ depends on the type of ANOVA.

For the most complex models $model_n$, each factors are assigned a natural number sequentially.
$factor_i$ is the included factors of $model_i$; $basefactor_i$ is the included factors of $basemodel_i$. 

Let the number of columns of model matrix of $model_n$ $m$.

A map $id_X: [1, m] \mapsto factor_n$ maps the index of columns into the corresponding factors.

We can define a vector of index set for each model:
$$I = (I_1, ..., I_{n})$$
where $ \forall i \in I_k, id_X(i) \notin basefactor_k, id_X(i) \in factor_{k}$

The degree of freedom is:
$$dof = (n(I_1), ..., n(I_{n}))$$ 
where $n(I)$ is the size of $I$

No deviance is computed. F-value is computed directly by the variance-covariance matrix ($vcov$) and the coefficients ($\beta$) of the model. Each $model_l$ corresponds to a $factor_l$.
1. Type I:

    Factors are sequentially added to the models.
    $basemodel_l$ is $model_l$ without $factor_l$.

    First, calculate the the upper factor of Cholesky factorization of $vcov^{-1}$ and multiply with $\beta$:
    $$ vcov^{-1} = LU $$
    $$ f = U\beta $$
    F-value is a vector $F$ where 
    $$ F_j = \frac{\sum_{k \in I_j}{f_k^2}}{dof_j}
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
