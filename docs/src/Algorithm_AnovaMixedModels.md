# AnovaMixedModels
## F-test
Given a model $M$, $n$ is the number of predictors, $m$ is the number of columns of $M$'s model matrix.

Define two sets, $\mathcal{C} = \{x \in \mathbb{N}\, |\, 1 \leq x \leq m\}$, the index of columns and $\mathcal{P} = \{x \in \mathbb{N}\, |\, 1 \leq x \leq n\}$, the index of predictors.

A map $id_X: \mathcal{C} \mapsto \mathcal{P}$ maps the index of columns into the corresponding predictor sequentially, i.e.,
```math
\begin{aligned}
    \forall i \in \mathcal{C}, id_X(i) = k &\implies i\text{th column} \text{ is a level of } k\text{th predictor}\\\\
    \forall i, j \in \mathcal{C}, i \lt j &\implies id_X(i) \leq id_X(j)
\end{aligned}
```
We can define a vector of index set for each predictors,
```math
\begin{aligned}
    \mathbf{I} &= (I_1, ..., I_n)\\\\
    \mathbf{df} &= (n(I_1), ..., n(I_n))
\end{aligned}
```
where $\forall i \in I_k, id_X(i) = k$, and $n(I)$ is the size of $I$.

F-value is a vector
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, \mathbf{df_r}}
```
where $\mathbf{df_r}$ is estimated by [between-within method](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#why-doesnt-lme4-display-denominator-degrees-of-freedomp-values-what-other-options-do-i-have).

F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the model. 
### Type I
Calculate F-value by the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiplying with $\boldsymbol \beta$
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \boldsymbol{\eta} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{\eta_k^2}}{df_j}
\end{aligned}
```
### Type II
Define two vectors of index sets $\mathbf J$ and $\mathbf K$ where
```math
\begin{aligned}
    J_j &= \{i \in \mathcal{C}\, |\, id_X(i) \text{ is an interaction term of }j\text{th predictor and other terms}\}\\\\
    K_j &= J_j \cup I_j
\end{aligned}
```
And F-value is
```math
F_j = \frac{\boldsymbol{\beta}_{K_j}^T \boldsymbol{\Sigma}_{K_j; K_j}^{-1} \boldsymbol{\beta}_{K_j} - \boldsymbol{\beta}_{J_j}^T \boldsymbol{\Sigma}_{J_j; J_j}^{-1} \boldsymbol{\beta}_{J_j}}{df_j}
```
### Type III
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```

## LRT
Given a vector of models
```math
\mathbf{M} = (M_1, ..., M_n)
``` 
Define the deviance $\mathbf{D}$ as -2loglikelihood of $\mathbf{M}$ and $\mathbf{dm}$ as the degrees of freedom of $\mathbf{M}$, the likelihood ratio is a vector
```math
\begin{aligned}
    \mathbf{L} &= \mathbf{D}_{[1, n - 1]} - \mathbf{D}_{[2, n]}\\\\
    \mathbf{df} &= \mathbf{dm}_{[1, n - 1]} - \mathbf{dm}_{[2, n]}\\\\
    \mathbf{L} &\sim \chi^2_{\mathbf{df}}
\end{aligned}
```