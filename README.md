# Anova
Implement one-way and multi-way anova, including type 1, 2 and 3. The syntax and output resemble GLM. 
Only fixed effect, independent ANOVA is supported now. Other variants will be further developed.

## Examples
```
julia> using RCall, Anova

julia> df = rcopy(R"iris")
150×5 DataFrames.DataFrame
│ Row │ Sepal_Length │ Sepal_Width │ Petal_Length │ Petal_Width │ Species   │
│     │ Float64      │ Float64     │ Float64      │ Float64     │ Cat…      │
├─────┼──────────────┼─────────────┼──────────────┼─────────────┼───────────┤
│ 1   │ 5.1          │ 3.5         │ 1.4          │ 0.2         │ setosa    │
│ 2   │ 4.9          │ 3.0         │ 1.4          │ 0.2         │ setosa    │
│ 3   │ 4.7          │ 3.2         │ 1.3          │ 0.2         │ setosa    │
│ 4   │ 4.6          │ 3.1         │ 1.5          │ 0.2         │ setosa    │
│ 5   │ 5.0          │ 3.6         │ 1.4          │ 0.2         │ setosa    │
│ 6   │ 5.4          │ 3.9         │ 1.7          │ 0.4         │ setosa    │
│ 7   │ 4.6          │ 3.4         │ 1.4          │ 0.3         │ setosa    │
│ 8   │ 5.0          │ 3.4         │ 1.5          │ 0.2         │ setosa    │
│ 9   │ 4.4          │ 2.9         │ 1.4          │ 0.2         │ setosa    │
│ 10  │ 4.9          │ 3.1         │ 1.5          │ 0.1         │ setosa    │
⋮
│ 140 │ 6.9          │ 3.1         │ 5.4          │ 2.1         │ virginica │
│ 141 │ 6.7          │ 3.1         │ 5.6          │ 2.4         │ virginica │
│ 142 │ 6.9          │ 3.1         │ 5.1          │ 2.3         │ virginica │
│ 143 │ 5.8          │ 2.7         │ 5.1          │ 1.9         │ virginica │
│ 144 │ 6.8          │ 3.2         │ 5.9          │ 2.3         │ virginica │
│ 145 │ 6.7          │ 3.3         │ 5.7          │ 2.5         │ virginica │
│ 146 │ 6.7          │ 3.0         │ 5.2          │ 2.3         │ virginica │
│ 147 │ 6.3          │ 2.5         │ 5.0          │ 1.9         │ virginica │
│ 148 │ 6.5          │ 3.0         │ 5.2          │ 2.0         │ virginica │
│ 149 │ 6.2          │ 3.4         │ 5.4          │ 2.3         │ virginica │
│ 150 │ 5.9          │ 3.0         │ 5.1          │ 1.8         │ virginica │

```
There's two way to conduct a ANOVA. First, fit a model with @formula like `GLM.lm`:
```
julia> model1 = anova(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df) # type 1
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 1 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────
                                        DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────
Sepal_Width                             1.0       1.41224          1.41224      15.44    0.0001
Petal_Length                            1.0      84.4273          84.4273      923.05    <1e-62
Species                                 2.0       2.36325          1.18162      12.92    <1e-5
Sepal_Width & Petal_Length              1.0       0.42066          0.42066       4.60    0.0337
Sepal_Width & Species                   2.0       0.300328         0.150164      1.64    0.1974
Petal_Length & Species                  2.0       0.557498         0.278749      3.05    0.0507
Sepal_Width & Petal_Length & Species    2.0       0.0647334        0.0323667     0.35    0.7026
(Residual)                            138.0      12.6223           0.0914659   NaN       NaN
───────────────────────────────────────────────────────────────────────────────────────────────

julia> model2 = anova(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df,type=2)
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 2 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────
                                        DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────
Sepal_Width                             1.0       2.91691          2.91691      31.89    <1e-7
Petal_Length                            1.0      14.445           14.445       157.93    <1e-23
Species                                 2.0       2.43369          1.21685      13.30    <1e-5
Sepal_Width & Petal_Length              1.0       0.0844784        0.0844784     0.92    0.3382
Sepal_Width & Species                   2.0       0.0672857        0.0336429     0.37    0.6929
Petal_Length & Species                  2.0       0.557498         0.278749      3.05    0.0507
Sepal_Width & Petal_Length & Species    2.0       0.0647334        0.0323667     0.35    0.7026
(Residual)                            138.0      12.6223           0.0914659   NaN       NaN
───────────────────────────────────────────────────────────────────────────────────────────────
```
Another one is conducting ANOVA with fitted model, only linear model is supported now: 
```
julia> lmmodel = anova(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df,type=3)
StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                   Estimate  Std. Error    t value  Pr(>|t|)   Lower 95%  Upper 95%
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                                       -1.36862     3.9089    -0.35013     0.7268   -9.0977     6.36046
Sepal_Width                                        1.723       1.12057    1.53761     0.1264   -0.492706   3.9387
Petal_Length                                       2.87592     2.74902    1.04616     0.2973   -2.55972    8.31156
Species: versicolor                                3.17458     5.49721    0.577488    0.5646   -7.69508   14.0442
Species: virginica                                -0.994886    5.27248   -0.188694    0.8506  -11.4202     9.4304
Sepal_Width & Petal_Length                        -0.743829    0.785394  -0.947077    0.3453   -2.29679    0.809133
Sepal_Width & Species: versicolor                 -1.35643     1.85653   -0.730627    0.4662   -5.02736    2.3145
Sepal_Width & Species: virginica                  -0.427356    1.65855   -0.257669    0.7970   -3.70681    2.8521
Petal_Length & Species: versicolor                -2.06748     2.89505   -0.714142    0.4763   -7.79188    3.65692
Petal_Length & Species: virginica                 -1.42331     2.81651   -0.505346    0.6141   -6.99242    4.14579
Sepal_Width & Petal_Length & Species: versicolor   0.716114    0.856938   0.835666    0.4048   -0.978313   2.41054
Sepal_Width & Petal_Length & Species: virginica    0.56492     0.812908   0.694938    0.4883   -1.04245    2.17229
───────────────────────────────────────────────────────────────────────────────────────────────────────────────────

julia> model3 = anova(lmmodel,type=3)
AnovaResult{StatsModels.TableRegressionModel{LinearModel{GLM.LmResp{Array{Float64,1}},GLM.DensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}},AnovaStats}

Type 3 ANOVA

Sepal_Length ~ 1 + Sepal_Width + Petal_Length + Species + Sepal_Width & Petal_Length + Sepal_Width & Species + Petal_Length & Species + Sepal_Width & Petal_Length & Species

Coefficients:
───────────────────────────────────────────────────────────────────────────────────────────────
                                        DOF  Sum of Squares  Mean of Squares  F value  Pr(>|F|)
───────────────────────────────────────────────────────────────────────────────────────────────
(Intercept)                             1.0       0.0112129        0.0112129     0.12    0.7268
Sepal_Width                             1.0       0.216248         0.216248      2.36    0.1264
Petal_Length                            1.0       0.100105         0.100105      1.09    0.2973
Species                                 2.0       0.061306         0.030653      0.34    0.7158
Sepal_Width & Petal_Length              1.0       0.0820408        0.0820408     0.90    0.3453
Sepal_Width & Species                   2.0       0.049004         0.024502      0.27    0.7654
Petal_Length & Species                  2.0       0.0625356        0.0312678     0.34    0.7111
Sepal_Width & Petal_Length & Species    2.0       0.0647334        0.0323667     0.35    0.7026
(Residual)                            138.0      12.6223           0.0914659   NaN       NaN
───────────────────────────────────────────────────────────────────────────────────────────────

```