# Anova
Implement one-way and multi-way anova, including type 1, 2 and 3. The syntax and output resemble GLM. Other variants of anova will be developed.

## Example
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

julia> model1 = anova(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df) # type 1
AnovaResult{AnovaModel{GLM.LmResp{Array{Float64,1}},AnovaDensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

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
AnovaResult{AnovaModel{GLM.LmResp{Array{Float64,1}},AnovaDensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

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

julia> model3 = anova(@formula(Sepal_Length~Sepal_Width*Petal_Length*Species),df,type=3)
AnovaResult{AnovaModel{GLM.LmResp{Array{Float64,1}},AnovaDensePredChol{Float64,LinearAlgebra.Cholesky{Float64,Array{Float64,2}}}},Array{Float64,2}}

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