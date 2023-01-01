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
There are vectors of models and the corresponding base models:
```math
\begin{aligned}
    \bf{M} &= (M_1, ..., M_n)\\\\
    \bf{B} &= (B_1, ..., B_n)
\end{aligned}
```
When $m$ models are given, $\bf M$ $= (M_2, ..., M_m)$, $\bf B$ $= (M_1, ..., M_{m-1})$. 

When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $\bf M$ and $\bf B$ depends on the type of ANOVA.

For the most complex model $M_n$, each factors are assigned a natural number $f_k$ sequentially where $f_l$ is the last factor.

Let the number of columns of model matrix of $M_n$ $m$.

The included factors of $M_j$ and $B_j$ are
```math
\begin{aligned}
    mf_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } M_i\}\\\\
    bf_j &= \{i \in [1, m]\, |\, f_i \text{ is a factor of } B_i\}
\end{aligned}
```
A map $id_X: [1, m] \mapsto [1, f_l]$ maps the index of columns into the corresponding factors.

We can define a vector of index sets for each model:
```math
\bf{I} &= (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) \in mf_k\setminus bf_k$.

The deviances are:
```math
\begin{aligned}
    \mathcal{D} &= (\mathcal{D}_1, ..., \mathcal{D}_n)\\\\
    \mathcal{R} &= (\mathcal{R}_1, ..., \mathcal{R}_n)
\end{aligned}
```
which is equivalent to the residual sum of squares.

The difference of $\mathcal{D}$ and $\mathcal{R}$ is:
```math
\bf{\Delta \mathcol{D} = \mathcol{D} - \mathcol{R}}
```
The degree of freedom is:
```math
\bf{dof} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$

The $rss$ is the residual sum of squares of $B_n$; $dof_{res}$ is the degree of freedoms of the residuals.

F-value is a vector:
```math
\bf{F} \sim \mathcal{F}_{\bf{dof}, dof_{res}}
```
where 
```math
F_i = \frac{\Delta dev_i \times dof_{res}}{rss^2 \times dof_i}
```
For a single model, F-value is computed directly by the variance-covariance matrix ($\bf \Sigma$) and the coefficients ($\beta$) of the model; the deviance is calculated backward. Each $M_j$ corresponds to a factor $f_j$, i.e. $id_X[I_j] = \{f_j\}$.
1. Type I:

    Factors are sequentially added to the models, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$

    Calculate the the upper factor of Cholesky factorization of $\bf \Sigma^{-1}$ and multiply with $\beta$: 
    $\bf{\Sigma}^{-1} = \bf{LU}$, $\bf{f} = \bf{U}\beta$
```math
F_j = \frac{\sum_{k \in I_j}{f_k^2}}{dof_j}
```

3. Type III:

    The models are all $M_n$, the base models are models without each factors.  
```math
F_j = \frac{\beta_{I_j}^T \bf{\Sigma}_{I_j; I_j}^{-1} \beta_{I_j}}{dof_j}
```