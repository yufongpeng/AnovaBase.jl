# Anova
Implement one-way and multi-way anova, including type 1, 2 and 3 sum of squares. The syntax and output resemble package `GLM`. 
The types of models supported:
1. Objects return by `GLM.lm`
2. `MixedModels.LinearMixedModel`
3. Objects return by `GLM.glm` (Multiple nested models only)

## Examples
### Simple linear model
```
julia> using RCall, Anova, DataFrames

julia> df = rcopy(R"iris")
150×5 DataFrame
│ Row │ Sepal_Length │ Sepal_Width │ Petal_Length │ Petal_Width │ Species   │
│     │ Float64      │ Float64     │ Float64      │ Float64     │ Cat…      │
├─────┼──────────────┼─────────────┼──────────────┼─────────────┼───────────┤
│ 1   │ 5.1          │ 3.5         │ 1.4          │ 0.2         │ setosa    │
│ 2   │ 4.9          │ 3.0         │ 1.4          │ 0.2         │ setosa    │
│ 3   │ 4.7          │ 3.2         │ 1.3          │ 0.2         │ setosa    │
⋮
│ 147 │ 6.3          │ 2.5         │ 5.0          │ 1.9         │ virginica │
│ 148 │ 6.5          │ 3.0         │ 5.2          │ 2.0         │ virginica │
│ 149 │ 6.2          │ 3.4         │ 5.4          │ 2.3         │ virginica │
│ 150 │ 5.9          │ 3.0         │ 5.1          │ 1.8         │ virginica │

```
There's two way to perform a ANOVA. First, fit a model with `@formula` like `GLM.lm`.
```
julia> anova_lm(@formula(Sepal_Length ~ Sepal_Width * Petal_Length * Petal_Width * Species), df)
Analysis of Variance

Type 1 test / F Test

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Petal_Width + Species + Sepal_Width & Petal_Length + Sepal_Width & Petal_Width + Petal_Length & Petal_Width + Sepal_Width & Species + Petal_Length & Species + Petal_Width & Species + Sepal_Width & Petal_Length & Petal_Width + Sepal_Width & Petal_Length & Species + Sepal_Width & Petal_Width & Species + Petal_Length & Petal_Width & Species + Sepal_Width 
& Petal_Length & Petal_Width & Species

Table:
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                    DOF  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                                           1        102.17           102.17    1108.6328    <1e-63
Sepal_Width                                           1          1.4122           1.4122    15.3242    0.0001
Petal_Length                                          1         84.43            84.43     916.1245    <1e-58
Petal_Width                                           1          1.8834           1.8834    20.4364    <1e-4
Species                                               2          0.8889           0.4445     4.8229    0.0096
Sepal_Width & Petal_Length                            1          0.2754           0.2754     2.9887    0.0863
Sepal_Width & Petal_Width                             1          0.0097           0.0097     0.1053    0.7461
Petal_Length & Petal_Width                            1          0.2120           0.2120     2.3006    0.1318
Sepal_Width & Species                                 2          0.2195           0.1097     1.1907    0.3074
Petal_Length & Species                                2          0.5075           0.2537     2.7532    0.0676
Petal_Width & Species                                 2          0.2738           0.1369     1.4854    0.2303
Sepal_Width & Petal_Length & Petal_Width              1          0.0191           0.0191     0.2074    0.6496
Sepal_Width & Petal_Length & Species                  2          0.1780           0.0890     0.9658    0.3835
Sepal_Width & Petal_Width & Species                   2          0.0665           0.0333     0.3608    0.6978
Petal_Length & Petal_Width & Species                  2          0.0322           0.0161     0.1748    0.8399
Sepal_Width & Petal_Length & Petal_Width & Species    2          0.1510           0.0755     0.8192    0.4431
(Residual)                                          126         11.61             0.0922
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
We can specify type of sum of squares by keyword argument `type`. Let's discard all interaction terms and use type II SS.
```
julia> aov = anova_lm(@formula(Sepal_Length ~ Sepal_Width + Petal_Length + Petal_Width + Species), df, type = 2)
Analysis of Variance

Type 2 test / F Test

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Petal_Width + Species + Sepal_Width & Petal_Length

Table:
─────────────────────────────────────────────────────────────────────────────────────
                            DOF  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
─────────────────────────────────────────────────────────────────────────────────────
(Intercept)                   1        102.17           102.17    1100.0685    <1e-68
Sepal_Width                   1          3.1250           3.1250    33.6476    <1e-7
Petal_Length                  1         13.79            13.79     148.4298    <1e-23
Petal_Width                   1          0.2638           0.2638     2.8404    0.0941
Species                       2          1.0058           0.5029     5.4150    0.0054
Sepal_Width & Petal_Length    1          0.2754           0.2754     2.9656    0.0872
(Residual)                  143         13.28             0.0929
─────────────────────────────────────────────────────────────────────────────────────
```
`anoava_lm` fit and store a `StatsModels.TableRegressionModel`.
```
julia> aov.model
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Petal_Width + Species

Coefficients:
──────────────────────────────────────────────────────────────────────────────────
                         Coef.  Std. Error      t  Pr(>|t|)  Lower 95%   Upper 95%
──────────────────────────────────────────────────────────────────────────────────
(Intercept)           2.17127    0.279794    7.76    <1e-11   1.61823    2.7243
Sepal_Width           0.495889   0.0860699   5.76    <1e-7    0.325765   0.666013
Petal_Length          0.829244   0.0685276  12.10    <1e-22   0.693794   0.964694
Petal_Width          -0.315155   0.151196   -2.08    0.0389  -0.614005  -0.0163054
Species: versicolor  -0.723562   0.240169   -3.01    0.0031  -1.19827   -0.24885
Species: virginica   -1.0235     0.333726   -3.07    0.0026  -1.68313   -0.363863
──────────────────────────────────────────────────────────────────────────────────
```
We can fit a model first and call `anova` instead. `anova` store the model as well, but it doesn't create a copy, so any in-place change of the original model should be noticed. 
```
julia> lm1 = lm(@formula(Sepal_Length ~ Sepal_Width + Petal_Length + Petal_Width + Species), df);

julia> anova(lm1, type = 3)
Analysis of Variance

Type 3 test / F Test

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Petal_Width + Species

Table:
──────────────────────────────────────────────────────────────────────
              DOF  Sum of Squares  Mean of Squares   F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────
(Intercept)     1          5.6694           5.6694   60.2211    <1e-11
Sepal_Width     1          3.1250           3.1250   33.1945    <1e-7
Petal_Length    1         13.79            13.79    146.4310    <1e-22
Petal_Width     1          0.4090           0.4090    4.3448    0.0389
Species         2          0.8889           0.4445    4.7212    0.0103
(Residual)    144         13.56             0.0941
──────────────────────────────────────────────────────────────────────
```
### Linear mixed-effect model
The implementation of ANOVA for linear mixed-effect model is primarily based on `MixedModels`. The syntax is similar to above examples. Only one random factor on intercept is supported now.
```
julia> using RCall, Anova, DataFrames

julia> R"""
       data("anxiety", package = "datarium")
       """
       df = stack(rcopy(R"anxiety"),[:t1,:t2,:t3],[:id,:group],variable_name=:time,value_name=:score)
       df = combine(df, Not(:time), :time => ByRow(x->parse(Int, replace(String(x), "t"=>""))) => :time)
135×4 DataFrame
│ Row │ id   │ group │ score   │ time  │
│     │ Cat… │ Cat…  │ Float64 │ Int64 │
├─────┼──────┼───────┼─────────┼───────┤
│ 1   │ 1    │ grp1  │ 14.1    │ 1     │
│ 2   │ 2    │ grp1  │ 14.5    │ 1     │
│ 3   │ 3    │ grp1  │ 15.7    │ 1     │
⋮
│ 132 │ 42   │ grp3  │ 13.8    │ 3     │
│ 133 │ 43   │ grp3  │ 15.4    │ 3     │
│ 134 │ 44   │ grp3  │ 15.1    │ 3     │
│ 135 │ 45   │ grp3  │ 15.5    │ 3     │
```
Fit a linear mixed-effect model with a `LinearMixedModel` object. `lme` is an alias for `fit(LinearMixedModel, formula, data, args...)`.
```
julia> fm1 = lme(@formula(score ~ group * time + (1|id)), df);

julia> anova(fm1)
Analysis of Variance

Type 1 test / F Test

score ~ 1 + group + time + group & time + (1 | id)

Table:
─────────────────────────────────────────────────────────────────────────────────────────
              DOF  Between-subjects  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
─────────────────────────────────────────────────────────────────────────────────────────
(Intercept)     1                 0        526.53           526.53    5019.3950    <1e-77
group           2                 1          2.0188           1.0094     4.4554    0.0176
time            1                 0         63.44            63.44     604.7794    <1e-40
group & time    2                 0         33.43            16.71     159.3210    <1e-29
(Residual)     42                 1          9.5155           0.2266
(Residual)     87                 0          9.1263           0.1049
─────────────────────────────────────────────────────────────────────────────────────────
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` with REML.
```
julia> aov = anova_lme(@formula(score ~ group * time + (1|id)), df, type = 3)
Analysis of Variance

Type 3 test / F Test

score ~ 1 + group + time + group & time + (1 | id)

Table:
─────────────────────────────────────────────────────────────────────────────────────────
              DOF  Between-subjects  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
─────────────────────────────────────────────────────────────────────────────────────────
(Intercept)     1                 0        190.63           190.63    1756.6503    <1e-58
group           2                 1          1.4193           0.7097     3.1340    0.0539
time            1                 0          2.5230           2.5230    23.2498    <1e-5
group & time    2                 0         34.98            17.49     161.1736    <1e-29
(Residual)     42                 1          9.5104           0.2264
(Residual)     87                 0          9.4410           0.1085
─────────────────────────────────────────────────────────────────────────────────────────

julia> aov.model
Linear mixed model fit by REML
 score ~ 1 + group + time + group & time + (1 | id)
 REML criterion at convergence: 276.9530790162688

Variance components:
            Column   VarianceStd.Dev.
id       (Intercept)  2.33800 1.52905
Residual              0.10852 0.32942
 Number of obs: 135; levels of grouping factors: 45

  Fixed-effects parameters:
────────────────────────────────────────────────────────────
                         Coef.  Std. Error       z  Pr(>|z|)
────────────────────────────────────────────────────────────
(Intercept)         17.42        0.415629    41.91    <1e-99
group: grp2         -0.0866667   0.587788    -0.15    0.8828
group: grp3          1.22889     0.587788     2.09    0.0366
time                -0.29        0.0601435   -4.82    <1e-5
group: grp2 & time  -0.27        0.0850558   -3.17    0.0015
group: grp3 & time  -1.43667     0.0850558  -16.89    <1e-63
────────────────────────────────────────────────────────────
```
To be noticed, type 2 sum of squares is not implemented now.

### Comparison between nested models
Nested models can be compared through Likelihood Ratio Test, F test, or other goodness of fit tests. For now, the first two are implemented.

Let's load the data first.
```
julia> R"library(ISLR)"; df = rcopy(R"Smarket")
1250×9 DataFrame
│ Row  │ Year    │ Lag1    │ Lag2    │ Lag3    │ Lag4    │ Lag5    │ Volume  │ Today   │ Direction │
│      │ Float64 │ Float64 │ Float64 │ Float64 │ Float64 │ Float64 │ Float64 │ Float64 │ Cat…      │
├──────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼───────────┤
│ 1    │ 2001.0  │ 0.381   │ -0.192  │ -2.624  │ -1.055  │ 5.01    │ 1.1913  │ 0.959   │ Up        │
│ 2    │ 2001.0  │ 0.959   │ 0.381   │ -0.192  │ -2.624  │ -1.055  │ 1.2965  │ 1.032   │ Up        │
│ 3    │ 2001.0  │ 1.032   │ 0.959   │ 0.381   │ -0.192  │ -2.624  │ 1.4112  │ -0.623  │ Down      │
⋮
│ 1247 │ 2005.0  │ 0.043   │ 0.422   │ 0.252   │ -0.024  │ -0.584  │ 1.28581 │ -0.955  │ Down      │
│ 1248 │ 2005.0  │ -0.955  │ 0.043   │ 0.422   │ 0.252   │ -0.024  │ 1.54047 │ 0.13    │ Up        │
│ 1249 │ 2005.0  │ 0.13    │ -0.955  │ 0.043   │ 0.422   │ 0.252   │ 1.42236 │ -0.298  │ Down      │
│ 1250 │ 2005.0  │ -0.298  │ 0.13    │ -0.955  │ 0.043   │ 0.422   │ 1.38254 │ -0.489  │ Down      │
```
We want to predict if the `Direction` is `Up` or `Down`. Let's use logistic regression with and without interaction terms, and compare this two models by Likelihood Ratio Test. The checker for nested models is not implemented now, so it should be ensured that the order of models should be sorted by increasing saturation of models. 
```
julia> glm_full = glm(@formula(Direction ~ Lag1 * Lag2 * Lag3 * Lag4 * Lag5 * Volume), df, Binomial(), LogitLink());

julia> glm_nointer = glm(@formula(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume), df, Binomial(), LogitLink());

julia> anova(glm_nointer, glm_full)
Analysis of Variance

Likelihood Ratio Test

Model 1: Direction ~ 1 + Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume
Model 2: Direction ~ 1 + Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume + Lag1 & Lag2 + Lag1 & Lag3 + Lag2 & Lag3 + Lag1 & Lag4 + Lag2 & Lag4 + Lag3 & Lag4 + Lag1 & Lag5 + Lag2 & Lag5 + Lag3 & Lag5 + Lag4 & Lag5 + Lag1 & Volume + Lag2 & Volume + Lag3 & Volume + Lag4 & Volume + Lag5 & Volume + Lag1 & Lag2 & Lag3 + Lag1 & Lag2 & Lag4 + Lag1 & Lag3 & Lag4 + Lag2 & Lag3 & Lag4 + Lag1 & Lag2 & Lag5 + Lag1 & Lag3 & Lag5 + Lag2 & Lag3 & Lag5 + Lag1 & Lag4 & Lag5 + Lag2 & Lag4 & Lag5 + Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Volume + Lag1 & Lag3 & Volume + Lag2 & Lag3 & Volume + Lag1 & Lag4 & Volume + Lag2 & Lag4 
& Volume + Lag3 & Lag4 & Volume + Lag1 & Lag5 & Volume + Lag2 & Lag5 & Volume + Lag3 & Lag5 & Volume + Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 + Lag1 & Lag2 & Lag3 & Lag5 + Lag1 & Lag2 & Lag4 
& Lag5 + Lag1 & Lag3 & Lag4 & Lag5 + Lag2 & Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Lag3 & Volume + Lag1 & Lag2 & Lag4 & Volume + Lag1 & Lag3 & Lag4 & Volume + Lag2 & Lag3 & Lag4 & Volume + Lag1 & Lag2 & Lag5 & Volume + Lag1 & Lag3 & Lag5 & Volume + Lag2 & Lag3 & Lag5 & Volume + Lag1 & Lag4 & Lag5 & Volume + Lag2 & Lag4 & Lag5 & Volume + Lag3 & Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Lag3 & Lag4 & Volume + Lag1 & Lag2 & Lag3 & Lag5 & Volume + Lag1 & Lag2 & Lag4 & Lag5 & Volume + Lag1 & Lag3 & Lag4 & Lag5 & Volume + Lag2 & Lag3 & Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 & Lag5 & Volume

Table:
─────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res. DOF  Deviance  Likelihood Ratio  Pr(>|χ²|)
─────────────────────────────────────────────────────────────
1    7            1244   1727.58
2   64    57      1187   1664.57           63.0155     0.2720
─────────────────────────────────────────────────────────────
```
Noticed that `anova` choose Likelihood Ratio Test automatically by detecting the type of model. For families of `Bernoulli()`, `Binomial()`, `Poisson()` that have fixed dispersion and `LinearMixedModel`, Likelihood Ratio Test is preferred. For other models, F test is preferred. To specify the test, we can use `anova(models..., test = <test>)` or `anova(<test>, models...)`. Use `FTest` for F test, `LikelihoodRatioTest` or `LRT` for Likelihood Ratio Test. 
```
julia> anova(FTest, glm_nointer, glm_full)
Analysis of Variance

F Test

Model 1: Direction ~ 1 + Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume
Model 2: Direction ~ 1 + Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume + Lag1 & Lag2 + Lag1 & Lag3 + Lag2 & Lag3 + Lag1 & Lag4 + Lag2 & Lag4 + Lag3 & Lag4 + Lag1 & Lag5 + Lag2 & Lag5 + Lag3 & Lag5 + Lag4 & Lag5 + Lag1 & Volume + Lag2 & Volume + Lag3 & Volume + Lag4 & Volume + Lag5 & Volume + Lag1 & Lag2 & Lag3 + Lag1 & Lag2 & Lag4 + Lag1 & Lag3 & Lag4 + Lag2 & Lag3 & Lag4 + Lag1 & Lag2 & Lag5 + Lag1 & Lag3 & Lag5 + Lag2 & Lag3 & Lag5 + Lag1 & Lag4 & Lag5 + Lag2 & Lag4 & Lag5 + Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Volume + Lag1 & Lag3 & Volume + Lag2 & Lag3 & Volume + Lag1 & Lag4 & Volume + Lag2 & Lag4 
& Volume + Lag3 & Lag4 & Volume + Lag1 & Lag5 & Volume + Lag2 & Lag5 & Volume + Lag3 & Lag5 & Volume + Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 + Lag1 & Lag2 & Lag3 & Lag5 + Lag1 & Lag2 & Lag4 
& Lag5 + Lag1 & Lag3 & Lag4 & Lag5 + Lag2 & Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Lag3 & Volume + Lag1 & Lag2 & Lag4 & Volume + Lag1 & Lag3 & Lag4 & Volume + Lag2 & Lag3 & Lag4 & Volume + Lag1 & Lag2 & Lag5 & Volume + Lag1 & Lag3 & Lag5 & Volume + Lag2 & Lag3 & Lag5 & Volume + Lag1 & Lag4 & Lag5 & Volume + Lag2 & Lag4 & Lag5 & Volume + Lag3 & Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 & Lag5 + Lag1 & Lag2 & Lag3 & Lag4 & Volume + Lag1 & Lag2 & Lag3 & Lag5 & Volume + Lag1 & Lag2 & Lag4 & Lag5 & Volume + Lag1 & Lag3 & Lag4 & Lag5 & Volume + Lag2 & Lag3 & Lag4 & Lag5 & Volume + Lag1 & Lag2 & Lag3 & Lag4 & Lag5 & Volume

Table:
───────────────────────────────────────────────────
   DOF  ΔDOF  Res. DOF  Deviance  F value  Pr(>|F|)
───────────────────────────────────────────────────
1    7            1244   1727.58
2   64    57      1187   1664.57   1.1055    0.2780
───────────────────────────────────────────────────
```
### TO DO
1. More statitics will be printed to pressent more information. 
2. Ommit some terms if the formula contains more than 10 terms. 
3. `anova` for single `GeneralizedLinearModel` will be implemented in different way compared to other models. All submodels will be fitted and stored. 
4. Implementation of `Rao` and `Mallow's Cp`.




