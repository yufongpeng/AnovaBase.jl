# MixedAnova
|CI status| Coverage|
|:-------:|:-------:|
| [![][travis-img]][travis-url] [![][ci-img]][ci-url]| [![][codecov-img]][codecov-url]|

[travis-img]: https://travis-ci.com/Jejulia/MixedAnova.jl.svg?branch=master
[travis-url]: https://travis-ci.com/github/Jejulia/MixedAnova.jl
[ci-img]: https://github.com/Jejulia/MixedAnova.jl/worlflows/CI/badge.svg
[ci-url]: https://github.com/Jejulia/MixedAnova.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/Jejulia/MixedAnova.jl/coveage.svg
[codecov-url]: https://codecov.io/gh/Jejulia/MixedAnova.jl

Implement one-way and multi-way anova, including type 1, 2 and 3 sum of squares. The syntax and output resemble package `GLM`. 
The types of models supported:
1. `TableRegressionModel{<: LinearModel, T}` fit by `GLM.lm`
2. `LinearMixedModel` fit by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`
3. `TableRegressionModel{<: GeneralizedLinearModel, T}` fit by `GLM.glm`

## Examples
### Simple linear model
```
julia> using RDatasets, MixedAnova, DataFrames

julia> df = dataset("datasets", "iris")
150×5 DataFrame
│ Row │ SepalLength │ SepalWidth │ PetalLength │ PetalWidth │ Species   │
│     │ Float64     │ Float64    │ Float64     │ Float64    │ Cat…      │
├─────┼─────────────┼────────────┼─────────────┼────────────┼───────────┤
│ 1   │ 5.1         │ 3.5        │ 1.4         │ 0.2        │ setosa    │
⋮
│ 149 │ 6.2         │ 3.4        │ 5.4         │ 2.3        │ virginica │
│ 150 │ 5.9         │ 3.0        │ 5.1         │ 1.8        │ virginica │

```
There's two way to perform a ANOVA. First, fit a model with `@formula` like `GLM.lm`.
```
julia> anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), df)
Analysis of Variance

Type 1 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
──────────────────────────────────────────────────────────────────────
             DOF  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────
(Intercept)    1        102.17           102.17    1085.2548    <1e-68
SepalWidth     1          1.4122           1.4122    15.0011    0.0002
PetalLength    1         84.43            84.43     896.8059    <1e-62
PetalWidth     1          1.8834           1.8834    20.0055    <1e-4
Species        2          0.8889           0.4445     4.7212    0.0103
(Residual)   144         13.56             0.0941
──────────────────────────────────────────────────────────────────────
```
We can specify type of sum of squares by keyword argument `type`. Let's use type II SS.
```
julia> aov = anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), df, type = 2)
Analysis of Variance

Type 2 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
──────────────────────────────────────────────────────────────────────
             DOF  Sum of Squares  Mean of Squares    F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────
(Intercept)    1        102.17           102.17    1085.2548    <1e-68
SepalWidth     1          3.1250           3.1250    33.1945    <1e-7
PetalLength    1         13.79            13.79     146.4310    <1e-22
PetalWidth     1          0.4090           0.4090     4.3448    0.0389
Species        2          0.8889           0.4445     4.7212    0.0103
(Residual)   144         13.56             0.0941
──────────────────────────────────────────────────────────────────────
```
`anova_lm` fit and store a `StatsModels.TableRegressionModel`.
```
julia> aov.model
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Coefficients:
──────────────────────────────────────────────────────────────────────────────────
                         Coef.  Std. Error      t  Pr(>|t|)  Lower 95%   Upper 95%
──────────────────────────────────────────────────────────────────────────────────
(Intercept)           2.17127    0.279794    7.76    <1e-11   1.61823    2.7243
SepalWidth            0.495889   0.0860699   5.76    <1e-7    0.325765   0.666013
PetalLength           0.829244   0.0685276  12.10    <1e-22   0.693794   0.964694
PetalWidth           -0.315155   0.151196   -2.08    0.0389  -0.614005  -0.0163054
Species: versicolor  -0.723562   0.240169   -3.01    0.0031  -1.19827   -0.24885
Species: virginica   -1.0235     0.333726   -3.07    0.0026  -1.68313   -0.363863
──────────────────────────────────────────────────────────────────────────────────
```
We can fit a model first and call `anova` instead. `anova` store the model as well, but it doesn't create a copy, so any in-place change of the original model should be noticed. 
```
julia> lm1 = lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), df);

julia> anova(lm1, type = 3)
Analysis of Variance

Type 3 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
─────────────────────────────────────────────────────────────────────
             DOF  Sum of Squares  Mean of Squares   F value  Pr(>|F|)
─────────────────────────────────────────────────────────────────────
(Intercept)    1          5.6694           5.6694   60.2211    <1e-11
SepalWidth     1          3.1250           3.1250   33.1945    <1e-7
PetalLength    1         13.79            13.79    146.4310    <1e-22
PetalWidth     1          0.4090           0.4090    4.3448    0.0389
Species        2          0.8889           0.4445    4.7212    0.0103
(Residual)   144         13.56             0.0941
─────────────────────────────────────────────────────────────────────
```
### Linear mixed-effect model
#### Random intercept
The implementation of ANOVA for linear mixed-effect model is primarily based on `MixedModels`. The syntax is similar to above examples. 
```
julia> using RCall, MixedAnova, DataFrames

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

Type 1 test / F test

score ~ 1 + group + time + group & time + (1 | id)

Table:
───────────────────────────────────────────────
              DOF  Res.DOF    F value  Pr(>|F|)
───────────────────────────────────────────────
(Intercept)     1       87  5019.3950    <1e-77
group           2       42     4.4554    0.0176
time            1       87   604.7794    <1e-40
group & time    2       87   159.3210    <1e-29
───────────────────────────────────────────────
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` with REML.
```
julia> aov = anova_lme(@formula(score ~ group * time + (1|id)), df, type = 3)
Analysis of Variance

Type 3 test / F test

score ~ 1 + group + time + group & time + (1 | id)

Table:
───────────────────────────────────────────────
              DOF  Res.DOF    F value  Pr(>|F|)
───────────────────────────────────────────────
(Intercept)     1       87  1756.6503    <1e-58
group           2       42     3.1340    0.0539
time            1       87    23.2498    <1e-5
group & time    2       87   161.1736    <1e-29
───────────────────────────────────────────────

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
#### Random slope
Multiple random effects and random slopes are also available.
```
julia> lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety);

julia> lmm2 = lme(@formula(score ~ group * time + (group|id)), anxiety);

julia> anova(lmm1)
Analysis of Variance

Type 1 test / F test

score ~ 1 + group + time + group & time + (1 | id)

Table:
───────────────────────────────────────────────
              DOF  Res.DOF    F value  Pr(>|F|)
───────────────────────────────────────────────
(Intercept)     1       87  5019.3950    <1e-77
group           2       42     4.4554    0.0176
time            1       87   604.7794    <1e-40
group & time    2       87   159.3210    <1e-29
───────────────────────────────────────────────

julia> anova(lmm2)
Analysis of Variance

Type 1 test / F test

score ~ 1 + group + time + group & time + (1 + group | id)

Table:
───────────────────────────────────────────────
              DOF  Res.DOF    F value  Pr(>|F|)
───────────────────────────────────────────────
(Intercept)     1       87  5098.9600    <1e-78
group           2       42     4.9178    0.0121
time            1       87   604.7808    <1e-40
group & time    2       87   159.3214    <1e-29
───────────────────────────────────────────────
```
#### Nested random effects
```
julia> aov = anova_lme(@formula(extro ~ open + agree + social + (1|school) + (1|school&class)), school)
Analysis of Variance

Type 1 test / F test

extro ~ 1 + open + agree + social + (1 | school) + (1 | school & class)

Table:
─────────────────────────────────────────────
             DOF  Res.DOF   F value  Pr(>|F|)
─────────────────────────────────────────────
(Intercept)    1     1173  227.2490    <1e-46
open           1     1173    1.4797    0.2241
agree          1     1173    1.8027    0.1796
social         1     1173    0.0851    0.7705
─────────────────────────────────────────────
```
Nested random slopes are also available.
### Generalized linear models
Not like the above examples, `anova_glm` and `anova` for `GLM.GeneralizedLinearModel` take the input model as saturated model, fit all simpler model and store all of them. 

```
julia> df = dataset("MASS", "quine")
146×5 DataFrame
│ Row │ Eth  │ Sex  │ Age  │ Lrn  │ Days  │
│     │ Cat… │ Cat… │ Cat… │ Cat… │ Int32 │
├─────┼──────┼──────┼──────┼──────┼───────┤
│ 1   │ A    │ M    │ F0   │ SL   │ 2     │
⋮
│ 145 │ N    │ F    │ F3   │ AL   │ 22    │
│ 146 │ N    │ F    │ F3   │ AL   │ 37    │

julia> nbm = glm(@formula(Days ~ Eth + Sex + Age + Lrn), df, NegativeBinomial(2.0), LogLink())
StatsModels.TableRegressionModel{GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},NegativeBinomial{Float64},LogLink},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

Days ~ 1 + Eth + Sex + Age + Lrn

Coefficients:
────────────────────────────────────────────────────────────────────────────
                  Coef.  Std. Error      z  Pr(>|z|)   Lower 95%   Upper 95%
────────────────────────────────────────────────────────────────────────────
(Intercept)   2.88645      0.227144  12.71    <1e-36   2.44125     3.33164
Eth: N       -0.567515     0.152449  -3.72    0.0002  -0.86631    -0.26872
Sex: M        0.0870771    0.159025   0.55    0.5840  -0.224606    0.398761
Age: F1      -0.445076     0.239087  -1.86    0.0627  -0.913678    0.0235251
Age: F2       0.0927999    0.234502   0.40    0.6923  -0.366816    0.552416
Age: F3       0.359485     0.246586   1.46    0.1449  -0.123814    0.842784
Lrn: SL       0.296768     0.185934   1.60    0.1105  -0.0676559   0.661191
────────────────────────────────────────────────────────────────────────────

julia> aov = anova(nbm)
Analysis of Variance

Type 1 test / F test

Days ~ 1 + Eth + Sex + Age + Lrn

Table:
──────────────────────────────────────────────────────────────────
             DOF   Deviance  Mean of deviance    F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────
(Intercept)    1  3667.69           3667.69    2472.1194    <1e-92
Eth            1    19.92             19.92      13.4252    0.0003
Sex            1     2.8515            2.8515     1.9220    0.1678
Age            3    14.42              4.8073     3.2403    0.0241
Lrn            1     3.8782            3.8782     2.6140    0.1082
──────────────────────────────────────────────────────────────────

julia> typeof(aov.model)
NTuple{6,StatsModels.TableRegressionModel{GeneralizedLinearModel{GLM.GlmResp{Array{Float64,1},NegativeBinomial{Float64},LogLink},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}}
```
To refit with other test, use the fit models from previous call. Arguments for other tests are pressented in the next section.
```
julia> anova(aov.model..., test = LRT)
Analysis of Variance

Likelihood-ratio test

Model 1: Days ~ 0
Model 2: Days ~ 1
Model 3: Days ~ 1 + Eth
Model 4: Days ~ 1 + Eth + Sex
Model 5: Days ~ 1 + Eth + Sex + Age
Model 6: Days ~ 1 + Eth + Sex + Age + Lrn

Table:
─────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res. DOF  Deviance  Likelihood Ratio  Pr(>|χ²|)
─────────────────────────────────────────────────────────────
1    1             146   3947.87
2    2     1       145    280.18         2472.1194     <1e-99
3    3     1       144    260.26           13.4252     0.0002
4    4     1       143    257.41            1.9220     0.1656
5    7     3       140    242.99            9.7208     0.0211
6    8     1       139    239.11            2.6140     0.1059
─────────────────────────────────────────────────────────────
```

### Comparison between nested models
Nested models can be compared through likelihood-ratio test, F test, or other goodness of fit tests. For now, the first two are implemented.
```
julia> df = dataset("datasets", "mtcars")
32×12 DataFrame
│ Row │ Model         │ MPG     │ Cyl   │ Disp    │ HP    │ DRat    │ WT      │ QSec    │ VS    │ AM    │ Gear  │ Carb  │
│     │ String        │ Float64 │ Int64 │ Float64 │ Int64 │ Float64 │ Float64 │ Float64 │ Int64 │ Int64 │ Int64 │ Int64 │
├─────┼───────────────┼─────────┼───────┼─────────┼───────┼─────────┼─────────┼─────────┼───────┼───────┼───────┼───────┤
│ 1   │ Mazda RX4     │ 21.0    │ 6     │ 160.0   │ 110   │ 3.9     │ 2.62    │ 16.46   │ 0     │ 1     │ 4     │ 4     │
⋮
│ 31  │ Maserati Bora │ 15.0    │ 8     │ 301.0   │ 335   │ 3.54    │ 3.57    │ 14.6    │ 0     │ 1     │ 5     │ 8     │
│ 32  │ Volvo 142E    │ 21.4    │ 4     │ 121.0   │ 109   │ 4.11    │ 2.78    │ 18.6    │ 1     │ 1     │ 4     │ 2     │
```
We want to predict if the `AM` is 0 or 1. Let's use logistic regression with and without interaction terms, and compare this two models by likelihood-ratio test. The checker for nested models is not implemented now, so it should be ensured that the later model is more saturated than previous one. 
```
julia> glm1 = glm(@formula(AM ~ Cyl + HP + WT), df, Binomial(), LogitLink());

julia> glm2 = glm(@formula(AM ~ Cyl * HP * WT), df, Binomial(), LogitLink());

julia> anova(glm1, glm2)
Analysis of Variance

Likelihood-ratio test

Model 1: AM ~ 1 + Cyl + HP + WT
Model 2: AM ~ 1 + Cyl + HP + WT + Cyl & HP + Cyl & WT + HP & WT + Cyl & HP & WT

Table:
─────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res. DOF  Deviance  Likelihood Ratio  Pr(>|χ²|)
─────────────────────────────────────────────────────────────
1    4              29    9.8415
2    8     4        25   <1e-6              9.8415     0.0432
─────────────────────────────────────────────────────────────
```
Noticed that `anova` choose likelihood-ratio test automatically by detecting the type of model. For families of `Bernoulli()`, `Binomial()`, `Poisson()` that have fixed dispersion and `LinearMixedModel`, likelihood-ratio test is preferred. For other models, F test is preferred. To specify the test, we can use `anova(models..., test = <test>)` or `anova(<test>, models...)`. Use `FTest` for F test, `LikelihoodRatioTest` or `LRT` for likelihood-ratio test. 
```
julia> anova(FTest, glm1, glm2)
Analysis of Variance

F test

Model 1: AM ~ 1 + Cyl + HP + WT
Model 2: AM ~ 1 + Cyl + HP + WT + Cyl & HP + Cyl & WT + HP & WT + Cyl & HP & WT

Table:
──────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res. DOF  Deviance  ΔDeviance  F value  Pr(>|F|)
──────────────────────────────────────────────────────────────
1    4              29    9.8415
2    8     4        25   <1e-6       9.8415   2.4604    0.0715
──────────────────────────────────────────────────────────────
```
## TO DO
1. More statitics will be printed to pressent more information. 
2. Ommit some terms if the formula contains more than 10 terms. 
3. Implementation of `Rao` and `Mallow's Cp`.
4. `anova` for `GeneralizedLinearMixedModels` and `FixedEffectModels`.




