# AnovaFixedEffectModels
Define a vector of models $\mathbf{M}$ and the corresponding base models $\mathbf{B}$
```math
\begin{aligned}
    \mathbf{M} &= (M_1, ..., M_n)\\\\
    \mathbf{B} &= (B_1, ..., B_n)
\end{aligned}
```
where $M_1$ is the simplest model with fixed effects and $M_n$ is the most complex model.

When $m$ models, $(M_1, ..., M_m)$, are given, $\mathbf{M} = (M_2, ..., M_m)$, $\mathbf{B} = (M_1, ..., M_{m-1})$. 

When one model is given, $n$ is the number of predictors except for the predictors used in the simplest model. The $\mathbf M$ and $\mathbf B$ depends on the type of ANOVA.

Let $m$, the number of columns of $M_n$'s model matrix; $l$, the number of predictors of $M_n$. 

Define two sets, $\mathcal{C} = \{x \in \mathbb{N}\, |\, 1 \leq x \leq m\}$, the index of columns and $\mathcal{P} = \{x \in \mathbb{N}\, |\, 1 \leq x \leq l\}$, the index of predictors.

A map $id_X: \mathcal{C} \mapsto \mathcal{P}$ maps the index of columns into the corresponding predictor sequentially, i.e., 
```math
\begin{aligned}
    \forall i \in \mathcal{C}, id_X(i) = k &\implies i\text{th column} \text{ is a level of } k\text{th predictor}\\\\
    \forall i, j \in \mathcal{C}, i \lt j &\implies id_X(i) \leq id_X(j)
\end{aligned}
```
The included predictors of $M_j$ and $B_j$ are $\mathcal{M}_j \subset \mathcal{P}$,  $\mathcal{B}_j \subset \mathcal{P}$, respectively.

We can define a vector of index sets for each model, and calulate degrees of freedom (dof) of each predictor
```math
\begin{aligned}
    \mathbf{I} &= (I_1, ..., I_n)\\\\
    \mathbf{df} &= (n(I_1), ..., n(I_n))
\end{aligned}
```
where $\forall i \in I_k, id_X(i) \in \mathcal{M}_k\setminus \mathcal{B}_k$, and $n(I)$ is the size of $I$.

The explained deviance of each predictor is the difference of $\mathbf{D}$ and $\mathbf{S}$
```math
\mathbf{E} = \mathbf{D} - \mathbf{S}
```
The mean explained deviance $\epsilon_i^2$ is therefore
```math
\epsilon_i^2 = \frac{E_i}{df_i}
```
The mean residual deviance $\sigma^2$ 
```math
\sigma^2 =\frac{D_n}{df_r}
```
where $D_n$ is the residual sum of squares of $M_n$; $df_r$ is the degrees of freedom of the residuals, i.e. $df_r = nob - n(\mathcal{C})$, where $nob$ is number of observations.

## F-test
F-value is a vector
```math
\mathbf{F} \sim \mathcal{F}_{\mathbf{df}, df_r}
```
where 
```math
F_i = \frac{\epsilon_i^2}{\sigma^2}
```
For a single model, F-value is computed directly by the variance-covariance matrix ($\boldsymbol \Sigma$) and the coefficients ($\boldsymbol \beta$) of the model, the deviance is calculated backward; each $M_j$ corresponds to a predictor $p_j$, i.e. $id_X[I_j] = \{j\}$.
### Type I
Predictors are sequentially added to the null model with fixed effects $B_1$, i.e., 
```math
\begin{aligned}
    \forall i, j \in \{x \in \mathbb{N}\, |\, 1\leq x\leq n\}, i < j &\implies (\mathcal{B}_i \subset \mathcal{B}_j) \land (\mathcal{M}_i \subset \mathcal{M}_j)\\\\
    \mathcal{M}_i &= \mathcal{B}_i \cup \{p_i\}\\\\
    \mathcal{B}_{i+1} &= \mathcal{M}_i
\end{aligned}
```
Calculate F-value by the the upper factor of Cholesky factorization of $\boldsymbol \Sigma^{-1}$ and multiplying with $\boldsymbol \beta$: 
```math
\begin{aligned}
    \boldsymbol{\Sigma}^{-1} &= \mathbf{LU}\\\\
    \boldsymbol{\eta} &= \mathbf{U}\boldsymbol{\beta}\\\\
    F_j &= \frac{\sum_{k \in I_j}{\eta_k^2}}{df_j}
\end{aligned}
```
### Type II 
The included predictors are defined as follows,
```math
\begin{aligned}
    \mathcal{B}_j &= \{k \in \mathcal{P}\, |\, k \text{ is not an interaction term of }p_j \text{ and other terms}\}\\\\
    \mathcal{M}_j &= \mathcal{B}_j \cup \{p_j\}
\end{aligned}
```
Define two vectors of index sets $\mathbf J$ and $\mathbf K$ where 
```math
\begin{aligned}
    J_j &= \{i \in \mathcal{C}\, |\, id_X(i) \text{ is an interaction term of }p_j \text{ and other terms}\}\\\\
    K_j &= J_j \cup I_j
\end{aligned}
```
And F-value is
```math
F_j = \frac{\boldsymbol{\beta}_{K_j}^T \boldsymbol{\Sigma}_{K_j; K_j}^{-1} \boldsymbol{\beta}_{K_j} - \boldsymbol{\beta}_{J_j}^T \boldsymbol{\Sigma}_{J_j; J_j}^{-1} \boldsymbol{\beta}_{J_j}}{df_j}
```
### Type III
All elements of $\mathbf{M}$ are the most complex model, and the base models are models without each predictors, i.e.
```math
\begin{aligned}
    \mathcal{M}_j &= \mathcal{P}\\\\
    \mathcal{B}_j &= \mathcal{P} \setminus \{p_j\}
\end{aligned}
```
And F-value is
```math
F_j = \frac{\boldsymbol{\beta}_{I_j}^T \boldsymbol{\Sigma}_{I_j; I_j}^{-1} \boldsymbol{\beta}_{I_j}}{df_j}
```
## LRT
The likelihood ratio is a vector
```math
\begin{aligned} 
    \mathbf{L} &= \frac{\mathbf{E}}{\sigma^2}\\\\
    \mathbf{L} &\sim \chi^2_{\mathbf{df}}
\end{aligned}
```