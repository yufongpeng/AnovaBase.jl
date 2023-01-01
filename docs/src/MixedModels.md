# AnovaMixedModels
```@setup mm
using AnovaMixedModels, DataFrames, GLM, CSV, CategoricalArrays
anxiety = CSV.read("anxiety.csv", DataFrame)
transform!(anxiety, :id => categorical, renamecols = false)
toenail = CSV.read("toenail.csv", DataFrame)
transform!(toenail, [1, 2, 3] .=> categorical, renamecols = false)
```
The implementation of ANOVA for [mixed-effects models](https://en.wikipedia.org/wiki/Mixed_model) is primarily based on [`MixedModels`](https://juliastats.org/MixedModels.jl/stable/). The syntax is similar to anova for `GLM`.   
```@example mm
using AnovaMixedModels
```
## Linear mixed-effects model
We get a dataset from `R` directly by [`RCall`](https://juliainterop.github.io/RCall.jl/stable/).
```julia
R"""data("anxiety", package = "datarium")"""
anxiety = stack(rcopy(R"anxiety"), [:t1, :t2, :t3], [:id, :group], variable_name = :time, value_name = :score)
anxiety = combine(anxiety, Not(:time), :time => ByRow(x->parse(Int, replace(String(x), "t"=>""))) => :time)
```
We can fit a linear mixed-effects model first. `lme` is an alias for [`fit(LinearMixedModel, formula, data, args...)`](https://juliastats.org/MixedModels.jl/stable/constructors/#Examples-of-linear-mixed-effects-model-fits).
```@example mm
lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety)
anova(lmm1)
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` fit by [Restricted maximum likelihood](https://en.wikipedia.org/wiki/Restricted_maximum_likelihood).
```@example mm
aov = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
```
```@example mm
aov.model.optsum.REML
```
!!! note
    Type 2 sum of squares is not implemented now.  
For [likeihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test), all submodels are fitted. The model should be fitted by [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation).
```@example mm
anova(LRT, lmm1)
```
When comparing multiple mixed models, likelihood-ratio test is used by default. 
It's also identical to [`StatsModels.lrtest`](https://juliastats.org/StatsModels.jl/stable/api/#StatsModels.lrtest) and [`MixedModels.likelihoodratiotest`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LikelihoodRatioTest).
```@example mm
lmms = nestedmodels(lmm1)
anova(lmms...) # as same as anova(LRT, lmm1)
MixedModels.likelihoodratiotest(lmms[2:end]...)
``` 
Comparing between [`LinearModel`](https://juliastats.org/GLM.jl/stable/api/#GLM.LinearModel) and [`LinearMixedModel`](https://juliastats.org/MixedModels.jl/stable/api/#MixedModels.LinearMixedModel) is also available.
```@example mm
lm1 = lm(@formula(score ~ group * time), anxiety)
lmm2 = lme(@formula(score ~ group * time + (group|id)), anxiety)
anova(lm1, lmm1, lmm2)
```
## Generalized linear mixed-effects model
The following is an example of generalized mixed model. `glme` is an alias for [`fit(GeneralizedLinearMixedModel, formula, data, args...)`](https://juliastats.org/MixedModels.jl/stable/constructors/#Fitting-generalized-linear-mixed-models).
```julia
R"""data("toenail", package = "HSAUR2")"""
toenail = rcopy(R"toenail")
```
```@example mm
glmm1 = glme(@formula(outcome ~ visit + treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
glmm2 = glme(@formula(outcome ~ visit * treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));
glm1 = glm(@formula(outcome ~ visit + treatment), toenail, Binomial(), LogitLink());
anova(glm1, glmm1, glmm2)
```
!!! note
    Only likelihood-ratio test is available now for `GeneralizedLinearMixedModel`.

## Algorithm
### F-test
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

## LRT
Vectors of the given models:
$$model = (model_1, ..., model_n)$$ 
The attribute `deviance` of the returned object ($dev$) is a vector of loglikelihoods $-2loglikelihood(model_i)$ for a linear mixed-effect model (linear model); unit deviance for a generalized linear mixed-effect model (generalized linear model).

The likelihood ratio is defined as $dev_{[1, n - 1]} - dev_{[2, n]}$. 