# Anova
Implement one-way and multi-way anova, including type 1, 2 and 3 sum of squares. The syntax and output resemble package GLM. 
Linear model for either fixed-effect only or mixed-effect is available. GLM will be further developed.

## Examples
### Simple linear model
```
julia> using RCall, Anova, DataFrame

julia> df = rcopy(R"iris")
150×5 DataFrame
│ Row │ Sepal_Length │ Sepal_Width │ Petal_Length │ Petal_Width │ Species   │
│     │ Float64      │ Float64     │ Float64      │ Float64     │ Cat…      │
├─────┼──────────────┼─────────────┼──────────────┼─────────────┼───────────┤
│ 1   │ 5.1          │ 3.5         │ 1.4          │ 0.2         │ setosa    │
│ 2   │ 4.9          │ 3.0         │ 1.4          │ 0.2         │ setosa    │
⋮

│ 148 │ 6.5          │ 3.0         │ 5.2          │ 2.0         │ virginica │
│ 149 │ 6.2          │ 3.4         │ 5.4          │ 2.3         │ virginica │
│ 150 │ 5.9          │ 3.0         │ 5.1          │ 1.8         │ virginica │

```
There's two way to perform a ANOVA. First, fit a model with `@formula` like `GLM.lm`.
```
julia> anova_lm(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df)
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 1 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────
                                        DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                             1.0     102.168          102.168      1117.01    <1e-67
Sepal_Width                             1.0       1.41224          1.41224      15.44    0.0001
Petal_Length                            1.0      84.4273          84.4273      923.05    <1e-62
Species                                 2.0       2.36325          1.18162      12.92    <1e-5
Sepal_Width & Petal_Length              1.0       0.42066          0.42066       4.60    0.0337
Sepal_Width & Species                   2.0       0.300328         0.150164      1.64    0.1974
Petal_Length & Species                  2.0       0.557498         0.278749      3.05    0.0507
Sepal_Width & Petal_Length & Species    2.0       0.0647334        0.0323667     0.35    0.7026
Residual                              138.0      12.6223           0.0914659   NaN       NaN
───────────────────────────────────────────────────────────────────────────────────────────────
```
We can specify type of sum of squares by keyword argument `type`.
```
julia> aov = anova_lm(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df,type=2)
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 2 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────
                                        DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                             1.0     102.168          102.168      1117.01    <1e-67
Sepal_Width                             1.0       2.91691          2.91691      31.89    <1e-7
Petal_Length                            1.0      14.445           14.445       157.93    <1e-23
Species                                 2.0       2.43369          1.21685      13.30    <1e-5
Sepal_Width & Petal_Length              1.0       0.0844784        0.0844784     0.92    0.3382
Sepal_Width & Species                   2.0       0.0672857        0.0336429     0.37    0.6929
Petal_Length & Species                  2.0       0.557498         0.278749      3.05    0.0507
Sepal_Width & Petal_Length & Species    2.0       0.0647334        0.0323667     0.35    0.7026
Residual                              138.0      12.6223           0.0914659   NaN       NaN
───────────────────────────────────────────────────────────────────────────────────────────────
```
`anoava_lm` fit and store a `StatsModels.TableRegressionModel`.
```
julia> aov.model
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                      Coef.  Std. Error      t  Pr(>|t|)   Lower 95%  Upper 95%
───────────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                                       -1.36862     3.9089    -0.35    0.7268   -9.0977     6.36046
Sepal_Width                                        1.723       1.12057    1.54    0.1264   -0.492706   3.9387
Petal_Length                                       2.87592     2.74902    1.05    0.2973   -2.55972    8.31156
Species: versicolor                                3.17458     5.49721    0.58    0.5646   -7.69508   14.0442
Species: virginica                                -0.994886    5.27248   -0.19    0.8506  -11.4202     9.4304
Sepal_Width & Petal_Length                        -0.743829    0.785394  -0.95    0.3453   -2.29679    0.809133
Sepal_Width & Species: versicolor                 -1.35643     1.85653   -0.73    0.4662   -5.02736    2.3145
Sepal_Width & Species: virginica                  -0.427356    1.65855   -0.26    0.7970   -3.70681    2.8521
Petal_Length & Species: versicolor                -2.06748     2.89505   -0.71    0.4763   -7.79188    3.65692
Petal_Length & Species: virginica                 -1.42331     2.81651   -0.51    0.6141   -6.99242    4.14579
Sepal_Width & Petal_Length & Species: versicolor   0.716114    0.856938   0.84    0.4048   -0.978313   2.41054
Sepal_Width & Petal_Length & Species: virginica    0.56492     0.812908   0.69    0.4883   -1.04245    2.17229
───────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
As we can see, no significant interaction occurs, so let's discard interaction term and fiting a model first.
```
julia> lm1 = lm(@formula(Sepal_Length~Sepal_Width+Petal_Length+Species),df)
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species

Coefficients:
─────────────────────────────────────────────────────────────────────────────────
                         Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────────────────
(Intercept)           2.39039    0.262268    9.11    <1e-15   1.87203    2.90875
Sepal_Width           0.432217   0.0813898   5.31    <1e-6    0.271354   0.593081
Petal_Length          0.775629   0.0642457  12.07    <1e-22   0.648651   0.902608
Species: versicolor  -0.955812   0.215199   -4.44    <1e-4   -1.38114   -0.530481
Species: virginica   -1.3941     0.285661   -4.88    <1e-5   -1.95869   -0.829501
─────────────────────────────────────────────────────────────────────────────────

julia> anova(lm1,type=3)
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 3 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species

Coefficients:
───────────────────────────────────────────────────────────────────────
                DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────
(Intercept)     1.0         8.00083        8.00083      83.07    <1e-15
Sepal_Width     1.0         2.71614        2.71614      28.20    <1e-6
Petal_Length    1.0        14.0382        14.0382      145.75    <1e-22
Species         2.0         2.36325        1.18162      12.27    <1e-4
Residual      145.0        13.9655         0.0963139   NaN       NaN
───────────────────────────────────────────────────────────────────────
```
### Linear mixed-effect model
The implementation of ANOVA for linear mixed-effect model is primarily based on `MixedModels` and R's nlme package. The syntax is similar to above examples. Only one random factor on intercept is supported now.
```
julia> using RCall, Anova, DataFrame

julia> R"""
       data("anxiety", package = "datarium")
       """
       df = stack(rcopy(R"anxiety"),[:t1,:t2,:t3],[:id,:group],variable_name=:time,value_name=:score)
135×4 DataFrame
│ Row │ id   │ group │ time │ score   │
│     │ Cat… │ Cat…  │ Cat… │ Float64 │
├─────┼──────┼───────┼──────┼─────────┤
│ 1   │ 1    │ grp1  │ t1   │ 14.1    │
│ 2   │ 2    │ grp1  │ t1   │ 14.5    │
⋮
│ 133 │ 43   │ grp3  │ t3   │ 15.4    │
│ 134 │ 44   │ grp3  │ t3   │ 15.1    │
│ 135 │ 45   │ grp3  │ t3   │ 15.5    │
```
Fit a linear mixed-effect model with a `LinearMixedModel` object. `lme` is an alias for `fit(LinearMixedModel,formula,data,args...)`.
```
julia> fm1 = lme(@formula(score ~ group * time + (1|id)),df)
Linear mixed model fit by maximum likelihood
 score ~ 1 + group + time + group & time + (1 | id)
   logLik   -2 logLik     AIC        BIC
 -119.76928  239.53857  261.53857  293.49659

Variance components:
            Column    Variance   Std.Dev.
id       (Intercept)  2.18966944 1.47975317
Residual              0.07867654 0.28049338
 Number of obs: 135; levels of grouping factors: 45

  Fixed-effects parameters:
────────────────────────────────────────────────────────────────
                             Coef.  Std. Error       z  Pr(>|z|)
────────────────────────────────────────────────────────────────
(Intercept)             17.0867       0.388874   43.94    <1e-99
group: grp2             -0.44         0.549951   -0.80    0.4237
group: grp3             -0.0733333    0.549951   -0.13    0.8939
time: t2                -0.16         0.102422   -1.56    0.1182
time: t3                -0.58         0.102422   -5.66    <1e-7
group: grp2 & time: t2  -0.02         0.144846   -0.14    0.8902
group: grp3 & time: t2  -1.84         0.144846  -12.70    <1e-36
group: grp2 & time: t3  -0.54         0.144846   -3.73    0.0002
group: grp3 & time: t3  -2.87333      0.144846  -19.84    <1e-86
────────────────────────────────────────────────────────────────

julia> anova(fm1)
AnovaResult{LinearMixedModel{Float64},AnovaStatsGrouped}

Type 1 ANOVA

score ~ 1 + group + time + group & time + (1 | id)

Coefficients:
──────────────────────────────────────────────────────────────────────────────────────────────────────────
                              DOF  Between-subjects  Sum of Squares  Mean of Squares     F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                   1.0               0.0       385.725             385.72  4902.67       <1e-75
group                         2.0               1.0         1.47605             0.74     4.35181    0.0192
time                          2.0               0.0        62.1402             31.07   394.91       <1e-42
group & time                  4.0               0.0        34.6767              8.67   110.188      <1e-31
Residual (between-subjects)  42.0             NaN           7.12279             0.17   NaN          NaN
Residual (within-subjects)   84.0             NaN           6.60883             0.08   NaN          NaN
──────────────────────────────────────────────────────────────────────────────────────────────────────────
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` with REML.
```
julia> aov = anova_lme(@formula(score ~ group * time + (1|id)),df,type=3)
AnovaResult{LinearMixedModel{Float64},AnovaStatsGrouped}

Type 3 ANOVA

score ~ 1 + group + time + group & time + (1 | id)

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────────────────
                              DOF  Between-subjects  Sum of Squares  Mean of Squares      F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                   1.0               0.0       151.895             151.89  1801.91        <1e-57
group                         2.0               1.0         0.11633             0.06     0.342975    0.7116
time                          2.0               0.0         2.692               1.35    15.9675      <1e-5
group & time                  4.0               0.0        37.1536              9.29   110.188       <1e-31
Residual (between-subjects)  42.0             NaN           7.12279             0.17   NaN           NaN
Residual (within-subjects)   84.0             NaN           7.08089             0.08   NaN           NaN
───────────────────────────────────────────────────────────────────────────────────────────────────────────

julia> aov.model
Linear mixed model fit by REML
 score ~ 1 + group + time + group & time + (1 | id)
 REML criterion at convergence: 256.63488307009317

Variance components:
            Column    Variance   Std.Dev.
id       (Intercept)  2.34607326 1.53168967
Residual              0.08429631 0.29033827
 Number of obs: 135; levels of grouping factors: 45

  Fixed-effects parameters:
────────────────────────────────────────────────────────────────
                             Coef.  Std. Error       z  Pr(>|z|)
────────────────────────────────────────────────────────────────
(Intercept)             17.0867       0.402523   42.45    <1e-99
group: grp2             -0.44         0.569253   -0.77    0.4396
group: grp3             -0.0733333    0.569253   -0.13    0.8975
time: t2                -0.16         0.106017   -1.51    0.1312
time: t3                -0.58         0.106017   -5.47    <1e-7
group: grp2 & time: t2  -0.02         0.14993    -0.13    0.8939
group: grp3 & time: t2  -1.84         0.14993   -12.27    <1e-33
group: grp2 & time: t3  -0.54         0.14993    -3.60    0.0003
group: grp3 & time: t3  -2.87333      0.14993   -19.16    <1e-81
────────────────────────────────────────────────────────────────
```
To be noticed, type 2 sum of squares is not implemented now.

We treat `time` as categorical variable; however time should be an ordered quantity. Let's turn time into continuous variable.
```
julia> using DataFrames

julia> df2 = combine(df,Not(:time),:time => ByRow(x->parse(Int,replace(String(x),"t"=>""))) => :time)
135×4 DataFrame
│ Row │ id   │ group │ score   │ time  │
│     │ Cat… │ Cat…  │ Float64 │ Int64 │
├─────┼──────┼───────┼─────────┼───────┤
│ 1   │ 1    │ grp1  │ 14.1    │ 1     │
│ 2   │ 2    │ grp1  │ 14.5    │ 1     │
⋮
│ 133 │ 43   │ grp3  │ 15.4    │ 3     │
│ 134 │ 44   │ grp3  │ 15.1    │ 3     │
│ 135 │ 45   │ grp3  │ 15.5    │ 3     │

julia> fm2 = lme(@formula(score ~ group * time + (1|id)),df2)
Linear mixed model fit by maximum likelihood
 score ~ 1 + group + time + group & time + (1 | id)
   logLik   -2 logLik     AIC        BIC    
 -132.71409  265.42819  281.42819  304.67039

Variance components:
            Column    Variance   Std.Dev.
id       (Intercept)  2.18092849 1.47679670
Residual              0.10489999 0.32388268
 Number of obs: 135; levels of grouping factors: 45

  Fixed-effects parameters:
────────────────────────────────────────────────────────────
                         Coef.  Std. Error       z  Pr(>|z|)
────────────────────────────────────────────────────────────
(Intercept)         17.42        0.402136    43.32    <1e-99
group: grp2         -0.0866667   0.568706    -0.15    0.8789
group: grp3          1.22889     0.568706     2.16    0.0307
time                -0.29        0.0591326   -4.90    <1e-6
group: grp2 & time  -0.27        0.0836261   -3.23    0.0012
group: grp3 & time  -1.43667     0.0836261  -17.18    <1e-65
────────────────────────────────────────────────────────────

julia> anova(fm2,type=3)
AnovaResult{LinearMixedModel{Float64},AnovaStatsGrouped}

Type 3 ANOVA

score ~ 1 + group + time + group & time + (1 | id)

Coefficients:
──────────────────────────────────────────────────────────────────────────────────────────────────────────
                              DOF  Between-subjects  Sum of Squares  Mean of Squares     F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                   1.0               0.0       188.097             188.10  1793.11       <1e-59
group                         2.0               1.0         1.44956             0.72     3.19908    0.0509
time                          1.0               0.0         2.41087             2.41    22.9825     <1e-5
group & time                  2.0               0.0        33.4255             16.71   159.321      <1e-29
Residual (between-subjects)  42.0             NaN           9.51549             0.23   NaN          NaN
Residual (within-subjects)   87.0             NaN           9.1263              0.10   NaN          NaN
──────────────────────────────────────────────────────────────────────────────────────────────────────────
```

