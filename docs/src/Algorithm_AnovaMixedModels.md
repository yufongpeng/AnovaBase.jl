# AnovaMixedModels
## F-test
Given a model $M$, $n$ is the number of factors, $m$ is the number of columns of $M$'s model matrix.

A map $id_X: [1, m] \mapsto [1, n]$ maps the index of columns into the corresponding factor sequentially, i.e. $\forall i, j \in [1, m], i \lt j \implies id_X(i) \leq id_X(j)$ and $\forall i \in [1, m], id_X(i) = k \implies \text{column}_i \text{ is a component of } k\text{th factor}$.

We can define a vector of index set for each factors:
```math
\mathbf{I} = (I_1, ..., I_n)
```
where $\forall i \in I_k, id_X(i) = k$.

The degrees of freedom (dof) is:
```math
\mathbf{df} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$.

F-value is a vector:
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, \mathbf{df_r}}
```
where $\mathbf{df_r}$ is estimated by between-within method.

F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the model. 
### Type I
Calculate the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiply with $\boldsymbol \beta$:
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \boldsymbol{\eta} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{\eta_k^2}}{df_j}
\end{aligned}
```

### Type III
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```

## LRT
Given a vector of models:
```math
\mathbf{M} = (M_1, ..., M_n)
``` 
The $\mathcal{D}$ is $-2loglikelihood(\mathbf{M})$ for linear mixed-effect models or ordinary linear models; unit deviance for generalized linear mixed-effect model or generalized linear models.

The likelihood ratio is a vector:
```math
\begin{aligned}
    \mathbf{L} &= \mathcal{D}_{[1, n - 1]} - \mathcal{D}_{[2, n]}\\\\
    \mathbf{L} &\sim \chi^2_{\mathbf{df}}
\end{aligned}
```
where $df_i = dof(M_i) - dof(M_{i+1})$