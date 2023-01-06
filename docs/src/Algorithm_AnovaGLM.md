
# AnovaGLM
Let a vector of models $\mathbf{M}$ and the corresponding base models $\mathbf{B}$:
```math
\begin{aligned}
    \mathbf{M} &= (M_1, ..., M_n)\\\\
    \mathbf{B} &= (B_1, ..., B_n)
\end{aligned}
```
where $M_1$ is the simplest model and $M_n$ is the most complex model.

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
where $\mathcal{D}_i$ and $\mathcal{R}_i$ are the sum of [squared deviance residuals (unit deviance)](https://en.wikipedia.org/wiki/Deviance_(statistics)) of $M_i$ and $B_i$. 
It is equivalent to the residual sum of squares for ordinary linear regression.

The difference of $\mathcal{D}$ and $\mathcal{R}$ is:
```math
\boldsymbol{\Delta} \mathcal{D} = \mathcal{D} - \mathcal{R}
```
The degrees of freedom (dof) are:
```math
\mathbf{df} = (n(I_1), ..., n(I_n))
```
where $n(I)$ is the size of $I$.

The $\sigma$ is the estimated dispersion (or scale) parameter for $M_n$'s distribution.

For ordinary linear regression, 
```math
\sigma^2 =\frac{rss}{df_r}
```
where $rss$ is the residual sum of squares of $B_n$; $df_r$ is the degrees of freedom of the residuals.

## F-test
F-value is a vector:
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, df_r}
```
where 
```math
F_i = \frac{\Delta \mathcal{D}_i}{\sigma^2 \times df_i}
```
For a single model, F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the most complex model; the deviance is calculated backward. Each $M_j$ corresponds to a factor $f_j$, i.e. $id_X[I_j] = \{f_j\}$.
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

### Type II 
The included facrors are defined as follows:
```math
\begin{aligned}
    \mathcal{B}_j &= \{k \in [1, l]\, |\, k \text{ is not an interaction term of }f_j \text{ and other terms}\}\\\\
    \mathcal{M}_j &= \mathcal{B}_j \cup \{f_j\}
\end{aligned}
```
Define two vectors of index sets $\mathbf J$ and $\mathbf K$ where 
```math
\begin{aligned}
    J_j &= \{i \in [1, m]\, |\, id_X(i) \text{ is an interaction term of }f_j \text{ and other terms}\}\\\\
    K_j &= J_j \cup I_j
\end{aligned}
```
And F-value is: 
```math
F_j = \frac{\boldsymbol{\beta}_{K_j}^T \boldsymbol{\Sigma}_{K_j; K_j}^{-1} \boldsymbol{\beta}_{K_j} - \boldsymbol{\beta}_{J_j}^T \boldsymbol{\Sigma}_{J_j; J_j}^{-1} \boldsymbol{\beta}_{J_j}}{df_j}
```

### Type III
The models are all $M_n$, the base models are models without each factors.  

F-value is:
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```

## LRT
The likelihood ratio is a vector:
```math
\begin{aligned} 
    \mathbf{L} &= \boldsymbol{\Delta} \mathcal{D}/\sigma^2\\\\
    \mathbf{L} &\sim \chi^2_{\mathbf{df}}
\end{aligned}
```
When a single model is provided, $\mathbf{L}$ is computed directly by the variance-covariance matrix.

Calculate the the upper factor of Cholesky factorization of $\sigma^2 \boldsymbol{\Sigma}^{-1}$ and multiply with $\boldsymbol \beta$:
```math
\begin{aligned}
    \sigma^2 \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \boldsymbol{\eta} &= \mathbf{U}\boldsymbol{\beta}\\\\
    L_j &= \sum_{k \in I_j}{\eta_k^2}
\end{aligned}
```