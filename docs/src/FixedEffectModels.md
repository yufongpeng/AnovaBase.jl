# AnovaFixedEffectModels
```@setup fem
using AnovaFixedEffectModels, DataFrames, CSV, CategoricalArrays
gpa = CSV.read("gpa.csv", DataFrame)
transform!(gpa,
        7 => x->replace(x, "yes" => true, "no" => false, "NA" => missing),
        4 => x->categorical(x, levels = ["1 hour", "2 hours", "3 hours"], ordered = true),
        renamecols = false)
transform!(gpa, [1, 2, 5, 7] .=> categorical, renamecols = false)
import AnovaBase: factornames
function factornames(aov::AnovaResult{T, FTest}; kwargs...) where {T <: TableRegressionModel{<: FixedEffectModel}}
    v = factornames(aov.model)
    push!(v, "(Residuals)")
    v
end
factornames(trm::TableRegressionModel{<: FixedEffectModel}) =
    vectorize(factornames(formula(trm).rhs.terms[unique(trm.mm.assign)])
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
Let a vector of models $\mathbf{M}$ and the corresponding base models $\mathbf{B}$:
```math
\begin{aligned}
    \mathbf{M} &= (M_1, ..., M_n)\\\\
    \mathbf{B} &= (B_1, ..., B_n)
\end{aligned}
```
When $m$ models, $(M_1, ..., M_m)$, are given, $\mathbf{M} = (M_2, ..., M_m)$, $\mathbf{B} = (M_1, ..., M_{m-1})$. 

When one model is given, $n$ is the number of factors except for the factors used in the simplest model. The $\mathbf M$ and $\mathbf B$ depends on the type of ANOVA.

Let the number of columns of $M_n$'s model matrix, $m$ and the number of factors of $M_n$, $l$. 

A map $id_X: [1, m] \mapsto [1, l]$ maps the index of columns into the corresponding factor sequentially, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$ and $\forall i \in [1, m], id_X(i) = k \implies \text{column}_i \text{ is a component of } k\text{th factor}$.

The included factors of $M_j$ and $B_j$ are:
```math
\begin{aligned}
    \mathcal{M}_j &= \{f \in [1, l]\, |\, f \text{ is a factor of } M_j\}\\\\
    \mathcal{B}_j &= \{f \in [1, l]\, |\, f \text{ is a factor of } B_j\}
\end{aligned}
```
We can define a vector of index sets for each model:
```math
\mathbf{I} = (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) \in \mathcal{M}_k\setminus \mathcal{B}_k$.

The deviances for models and base models are:
```math
\begin{aligned}
    \mathcal{D} &= (\mathcal{D}_1, ..., \mathcal{D}_n)\\\\
    \mathcal{R} &= (\mathcal{R}_1, ..., \mathcal{R}_n)
\end{aligned}
```
which is equivalent to the residual sum of squares.

The difference of $\mathcal{D}$ and $\mathcal{R}$ is:
```math
\boldsymbol{\Delta} \mathcal{D} = \mathcal{D} - \mathcal{R}
```
The degrees of freedom (dof) is:
```math
\mathbf{df} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$.

F-value is a vector:
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, df_r}
```
where 
```math
F_i = \frac{\Delta \mathcal{D}_i \times df_r}{rss^2 \times df_i}
```
and $rss$ is the residual sum of squares of $B_n$; $df_r$ is the degrees of freedom of the residuals.

For a single model, F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the model; the deviance is calculated backward. Each $M_j$ corresponds to a factor $f_j$, i.e. $id_X[I_j] = \{f_j\}$.
### Type I
Factors are sequentially added to the models, i.e. $\forall i, j \in [1, n], i < j \implies (\mathcal{B}_i \subset \mathcal{B}_j) \cap (\mathcal{M}_i \subset \mathcal{M}_j)$.

Calculate the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiply with $\boldsymbol \beta$: 
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \boldsymbol{\eta} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{\eta_k^2}}{df_j}
\end{aligned}
```

### Type III:
The models are all $M_n$, the base models are models without each factors.  
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```