# MixedAnova
|CI status|Coverage|
|:-------:|:------:|
| [![TravisCI][travis-img]][travis-url] [![][ci-img]][ci-url]| [![][codecov-img]][codecov-url]|

[travis-img]: https://travis-ci.com/Jejulia/MixedAnova.jl.svg?branch=master
[travis-url]: https://travis-ci.com/github/Jejulia/MixedAnova.jl
[ci-img]: https://github.com/Jejulia/MixedAnova.jl/workflows/CI/badge.svg
[ci-url]: https://github.com/Jejulia/MixedAnova.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/Jejulia/MixedAnova.jl/coveage.svg
[codecov-url]: https://codecov.io/gh/Jejulia/MixedAnova.jl

Implement one-way and multi-way anova, including type 1, 2 and 3 sum of squares. The syntax and output resemble package `GLM`. 
The types of models supported:
1. `TableRegressionModel{<: LinearModel, T}` fit by `GLM.lm`
2. `TableRegressionModel{<: GeneralizedLinearModel, T}` fit by `GLM.glm`
3. `LinearMixedModel` fit by `MixedAnova.lme` or `fit(LinearMixedModel, ...)`

## Functionality
### anova
```
anova(<model>; <type>, <test>)
anova(<test>, <model>; <type>)
anova(<models>; <test>)
anova(<test>, <models>)

anova_lm(<formula>, <data>; <type>, <test>)
anova_lm(<test>, <formula>, <data>; <type>)
anova_glm(<formula>, <data>, <distr>, <link>; <type>, <test>)
anova_glm(<test>, <formula>, <data>, <distr>, <link>; <type>)
anova_lme(<formula>, <data>; <type>, <test>)
anova_lme(<test>, <formula>, <data>; <type>)
```
### nestedmodels
```
nestedmodels(<modeltype>, <formula>, <data>)
nestedmodels(<model>)
```
## Examples
### Simple linear model
```
julia> using MixedAnova

julia> glm_init() # To enable GLM fuctionality
```
This function will export all functions from `GLM` and related function in this package, including `anova`, `anova_lm`, `anova_glm`.
```
julia> using RDatasets, DataFrames

julia> iris = dataset("datasets", "iris")
150×5 DataFrame
 Row │ SepalLength  SepalWidth  PetalLength  PetalWidth  Species   
     │ Float64      Float64     Float64      Float64     Cat…      
─────┼─────────────────────────────────────────────────────────────
   1 │         5.1         3.5          1.4         0.2  setosa
   2 │         4.9         3.0          1.4         0.2  setosa
   3 │         4.7         3.2          1.3         0.2  setosa
  ⋮  │      ⋮           ⋮            ⋮           ⋮           ⋮
 149 │         6.2         3.4          5.4         2.3  virginica
 150 │         5.9         3.0          5.1         1.8  virginica
                                                   145 rows omitted

```
There's two way to perform a ANOVA. First, fit a model with `@formula` like `GLM.lm`.
```
julia> anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris)
Analysis of Variance

Type 1 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
──────────────────────────────────────────────────────────────
             DOF     Exp.SS  Mean Square     F value  Pr(>|F|)
──────────────────────────────────────────────────────────────
(Intercept)    1  5121.68      5121.68    54403.6419    <1e-99
SepalWidth     1     1.4122       1.4122     15.0011    0.0002
PetalLength    1    84.43        84.43      896.8059    <1e-62
PetalWidth     1     1.8834       1.8834     20.0055    <1e-04
Species        2     0.8889       0.4445      4.7212    0.0103
(Residuals)  144    13.56         0.0941
──────────────────────────────────────────────────────────────
```
We can specify type of sum of squares by keyword argument `type`. Let's use type II SS.
```
julia> anova_lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris, type = 2)
Analysis of Variance

Type 2 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
──────────────────────────────────────────────────────────────
             DOF     Exp.SS  Mean Square     F value  Pr(>|F|)
──────────────────────────────────────────────────────────────
(Intercept)    1  5121.68         0.0002  54403.6419    <1e-99
SepalWidth     1     3.1250       0.3200     33.1945    <1e-07
PetalLength    1    13.79         0.0725    146.4310    <1e-22
PetalWidth     1     0.4090       2.4448      4.3448    0.0389
Species        2     0.8889       2.2499      4.7212    0.0103
(Residuals)  144    13.56        10.62
──────────────────────────────────────────────────────────────
```
`anova_lm` fit and store a `StatsModels.TableRegressionModel`.  

We can fit a model first and call `anova` instead. `anova` store the model as well, but it doesn't create a copy, so any in-place change of the original model should be noticed. 
```
julia> lm1 = lm(@formula(SepalLength ~ SepalWidth + PetalLength + PetalWidth + Species), iris);

julia> anova(lm1, type = 3)
Analysis of Variance

Type 3 test / F test

SepalLength ~ 1 + SepalWidth + PetalLength + PetalWidth + Species

Table:
──────────────────────────────────────────────────────────
             DOF   Exp.SS  Mean Square   F value  Pr(>|F|)
──────────────────────────────────────────────────────────
(Intercept)    1   5.6694       5.6694   60.2211    <1e-11
SepalWidth     1   3.1250       3.1250   33.1945    <1e-07
PetalLength    1  13.79        13.79    146.4310    <1e-22
PetalWidth     1   0.4090       0.4090    4.3448    0.0389
Species        2   0.8889       0.4445    4.7212    0.0103
(Residuals)  144  13.56         0.0941
──────────────────────────────────────────────────────────
```
Multiple models can be compared through the same function.  
The checker for nested models is not implemented now, so it should be ensured that the later model is more saturated than the previous one.  
```
julia> lms = nestedmodels(LinearModel, @formula(SepalLength ~ SepalWidth * Species), iris, dropcollinear = false);

julia> anova(lms...)
┌ Warning: Could not check whether models are nested: results may not be meaningful
└ @ MixedAnova ...
Analysis of Variance

Type 1 test / F test

Model 1: SepalLength ~ 0
Model 2: SepalLength ~ 1
Model 3: SepalLength ~ 1 + SepalWidth
Model 4: SepalLength ~ 1 + SepalWidth + Species
Model 5: SepalLength ~ 1 + SepalWidth + Species + SepalWidth & Species

Table:
──────────────────────────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF        R²      ΔR²   Res.SS     Exp.SS     F value  Pr(>|F|)
──────────────────────────────────────────────────────────────────────────────────
1    1            150  -50.13             5223.85
2    2     1      149   <1e-15   50.13     102.17  5121.68    26485.3010    <1e-99
3    3     1      148    0.0138   0.0138   100.76     1.4122      7.3030    0.0077
4    5     2      146    0.7259   0.7121    28.00    72.75      188.1091    <1e-40
5    7     2      144    0.7274   0.0015    27.85     0.1572      0.4064    0.6668
──────────────────────────────────────────────────────────────────────────────────
```
The result is a little bit different from `GLM.ftest`
```
julia> ftest(getproperty.(lms[2:end], :model)...)
F-test: 4 models fitted on 150 observations
────────────────────────────────────────────────────────────────────
     DOF  ΔDOF       SSR      ΔSSR      R²     ΔR²        F*   p(>F)
────────────────────────────────────────────────────────────────────
[1]    2        102.1683            0.0000
[2]    3     1  100.7561   -1.4122  0.0138  0.0138    2.0744  0.1519
[3]    5     2   28.0037  -72.7524  0.7259  0.7121  189.6512  <1e-40
[4]    7     2   27.8465   -0.1572  0.7274  0.0015    0.4064  0.6668
────────────────────────────────────────────────────────────────────
```
In `anova`, the F value is calculated by dividing MSR(mean of ΔDeviance) with mean of RSS of saturated model, like `anova` in R programming language, while in `GLM.ftest`, the denominator is replaced by RSS of subsequant model.
### Generalized linear models 
```
julia> quine = dataset("MASS", "quine")
146×5 DataFrame
 Row │ Eth   Sex   Age   Lrn   Days  
     │ Cat…  Cat…  Cat…  Cat…  Int32 
─────┼───────────────────────────────
   1 │ A     M     F0    SL        2
   2 │ A     M     F0    SL       11
   3 │ A     M     F0    SL       14
  ⋮  │  ⋮     ⋮     ⋮     ⋮      ⋮
 145 │ N     F     F3    AL       22
 146 │ N     F     F3    AL       37
                     141 rows omitted

julia> nbm = glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink());

julia> anova(nbm)
Analysis of Variance

Type 1 test / F test

Days ~ 1 + Eth + Sex + Age + Lrn

Table:
───────────────────────────────────────────────────────────
             DOF  ΔDeviance  Mean ΔDev    F value  Pr(>|F|)
───────────────────────────────────────────────────────────
(Intercept)    1  3667.69    3667.69    2472.0054    <1e-89
Eth            1    19.92      19.92      13.4245    0.0004
Sex            1     2.8515     2.8515     1.9219    0.1679
Age            3    14.42       4.8074     3.2401    0.0241
Lrn            1     3.8781     3.8781     2.6139    0.1082
(Residuals)  139   239.11       1.7202
───────────────────────────────────────────────────────────
```
There's also `anova_glm` similar to `anova_lm`.  

`anova` will automatically select test from f test or likelihood-ratio test depending on the type of distribution. For distribution of `Bernoulli()`, `Binomial()`, `Poisson()` that have fixed dispersion, likelihood-ratio test is preferred. For other distribution, F test is preferred.  
```
julia> mtcars = dataset("datasets", "mtcars")
32×12 DataFrame
 Row │ Model              MPG      Cyl    Disp     HP     DRat     WT       QSec     VS    AM     Gear   Carb  
     │ Cat…               Float64  Int64  Float64  Int64  Float64  Float64  Float64  Cat…  Int64  Int64  Int64 
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Mazda RX4             21.0      6    160.0    110     3.9     2.62     16.46  0         1      4      4
   2 │ Mazda RX4 Wag         21.0      6    160.0    110     3.9     2.875    17.02  0         1      4      4
   3 │ Datsun 710            22.8      4    108.0     93     3.85    2.32     18.61  1         1      4      1
  ⋮  │         ⋮             ⋮       ⋮       ⋮       ⋮       ⋮        ⋮        ⋮      ⋮      ⋮      ⋮      ⋮
  31 │ Maserati Bora         15.0      8    301.0    335     3.54    3.57     14.6   0         1      5      8
  32 │ Volvo 142E            21.4      4    121.0    109     4.11    2.78     18.6   1         1      4      2
                                                                                                27 rows omitted
```
We want to predict if the `AM` is 0 or 1. Let's use logistic regression with and without interaction terms, and compare this two models by likelihood-ratio test. 
```
julia> glm1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink());

julia> glm2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), LogitLink());

julia> anova(glm1, glm2)
┌ Warning: Could not check whether models are nested: results may not be meaningful
└ @ MixedAnova ...
Analysis of Variance

Type 1 test / Likelihood-ratio test

Model 1: AM ~ 1 + Cyl + HP + WT
Model 2: AM ~ 1 + Cyl + HP + WT + Cyl & HP + Cyl & WT + HP & WT + Cyl & HP & WT

Table:
──────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF  Deviance      χ²  Pr(>|χ²|)
──────────────────────────────────────────────────
1    4             29    9.8415
2    8     4       25   <1e-06   9.8415     0.0432
──────────────────────────────────────────────────

julia> lrtest(glm1, glm2)
┌ Warning: Could not check whether models are nested as model type TableRegressionModel does not implement isnested: results may not be meaningful
└ @ StatsModels ...
Likelihood-ratio test: 2 models fitted on 32 observations
──────────────────────────────────────────────
     DOF  ΔDOF  Deviance  ΔDeviance  p(>Chisq)
──────────────────────────────────────────────
[1]    4          9.8415
[2]    8     4    0.0000    -9.8415     0.0432
──────────────────────────────────────────────
```
This function works identically as `StatsModels.lrtest`

We can also specify test by keword arguments `test` or putting test in the first argument.

### Linear mixed-effect model
The implementation of ANOVA for linear mixed-effect model is primarily based on `MixedModels`. The syntax is similar to above examples.   
Likewise, to enable `MixedModels` fuctionality: 
```
julia> mm_init()

julia> using RCall, DataFrames

julia> R"""data("anxiety", package = "datarium")"""; anxiety = stack(rcopy(R"anxiety"), [:t1, :t2, :t3], [:id, :group], variable_name = :time, value_name = :score); anxiety = combine(anxiety, Not(:time), :time => ByRow(x->parse(Int, replace(String(x), "t"=>""))) => :time)
135×4 DataFrame
 Row │ id    group  score    time  
     │ Cat…  Cat…   Float64  Int64 
─────┼─────────────────────────────
   1 │ 1     grp1      14.1      1
   2 │ 2     grp1      14.5      1
   3 │ 3     grp1      15.7      1
  ⋮  │  ⋮      ⋮       ⋮       ⋮
 134 │ 44    grp3      15.1      3
 135 │ 45    grp3      15.5      3
                   130 rows omitted
```
Fit a linear mixed-effect model with a `LinearMixedModel` object. `lme` is an alias for `fit(LinearMixedModel, formula, data, args...)`.
```
julia> lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety);

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
```
Alternatively, we can use `anova_lme`. Like `anova_lm`, this function will fit and store a model; in this case, a `LinearMixedModel` fit by REML.
```
julia> aov = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
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

julia> aov.model.optsum.REML
true
```
To be noticed, type 2 sum of squares is not implemented now.  
For likeihood-ratio test, all submodels are fitted. The model should be fit by ML.
```
julia> anova(LRT, lmm1)
┌ Warning: fit all submodels
└ @ MixedAnova ...
┌ Warning: Could not check whether models are nested: results may not be meaningful
└ @ MixedAnova ...
Analysis of Variance

Type 1 test / Likelihood-ratio test

Model 1: score ~ 0 + (1 | id)
Model 2: score ~ 1 + (1 | id)
Model 3: score ~ 1 + group + (1 | id)
Model 4: score ~ 1 + group + time + (1 | id)
Model 5: score ~ 1 + group + time + group & time + (1 | id)

Table:
─────────────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF     AIC     BIC  -2 logLik        χ²  Pr(>|χ²|)
─────────────────────────────────────────────────────────────────────
1    2            133  705.74  711.55     701.74
2    3     1      132  501.55  510.27     495.55  206.1822     <1e-46
3    5     2      130  497.08  511.61     487.08    8.4747     0.0144
4    6     1      129  416.81  434.24     404.81   82.2717     <1e-18
5    8     2      127  281.43  304.67     265.43  139.3790     <1e-30
─────────────────────────────────────────────────────────────────────
```
When comparing multiple mixed models, likelihood-ratio test is used by default. 
It's also identical to `StatsModels.lrtest` and `MixedMOdels.likelihoodratiotest`.
```
julia> lmms = nestedmodels(lmm1);

julia> anova(lmms...)
Analysis of Variance

Type 1 test / Likelihood-ratio test

Model 1: score ~ 0 + (1 | id)
Model 2: score ~ 1 + (1 | id)
Model 3: score ~ 1 + group + (1 | id)
Model 4: score ~ 1 + group + time + (1 | id)
Model 5: score ~ 1 + group + time + group & time + (1 | id)

Table:
─────────────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF     AIC     BIC  -2 logLik        χ²  Pr(>|χ²|)
─────────────────────────────────────────────────────────────────────
1    2            133  705.74  711.55     701.74
2    3     1      132  501.55  510.27     495.55  206.1822     <1e-46
3    5     2      130  497.08  511.61     487.08    8.4747     0.0144
4    6     1      129  416.81  434.24     404.81   82.2717     <1e-18
5    8     2      127  281.43  304.67     265.43  139.3790     <1e-30
─────────────────────────────────────────────────────────────────────

julia> MixedModels.likelihoodratiotest(lmms[2:end]...)
Model Formulae
1: score ~ 1 + (1 | id)
2: score ~ 1 + group + (1 | id)
3: score ~ 1 + group + time + (1 | id)
4: score ~ 1 + group + time + group & time + (1 | id)
───────────────────────────────────────────────────
     model-dof  -2 logLik        χ²  χ²-dof  P(>χ²)
───────────────────────────────────────────────────
[1]          3   495.5536
[2]          5   487.0789    8.4747       2  0.0144
[3]          6   404.8072   82.2717       1  <1e-18
[4]          8   265.4282  139.3790       2  <1e-30
───────────────────────────────────────────────────
``` 
Comparing between `LinearModel` and `LinearMixedModels` is also available.
```
julia> lm1 = lm(@formula(score ~ group * time), anxiety); lmm2 = lme(@formula(score ~ group * time + (group|id)), anxiety);

julia> anova(lm1, lmm1, lmm2)
Analysis of Variance

Type 1 test / Likelihood-ratio test

Model 1: score ~ 1 + group + time + group & time
Model 2: score ~ 1 + group + time + group & time + (1 | id)
Model 3: score ~ 1 + group + time + group & time + (1 + group | id)

Table:
─────────────────────────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF     AIC     BIC  -2 logLik        χ²  Pr(>|χ²|)
─────────────────────────────────────────────────────────────────────
1    7            129  508.72  529.06     494.72
2    8     1      127  281.43  304.67     265.43  229.2935     <1e-51
3   13     5      122  290.79  328.56     264.79    0.6349     0.9864
─────────────────────────────────────────────────────────────────────
```
### Generalized linear mixed-effect model
Fit a linear mixed-effect model with a `GeneralizedLinearMixedModel` object. `glme` is an alias for `fit(GeneralizedLinearMixedModel, formula, data, args...)`.
```
julia> R"""data("toenail", package = "HSAUR2")"""; toenail = rcopy(R"toenail")
1908×5 DataFrame
  Row │ patientID  outcome             treatment    time       visit 
      │ Cat…       Cat…                Cat…         Float64    Int64 
──────┼──────────────────────────────────────────────────────────────
    1 │ 1          moderate or severe  terbinafine   0.0           1
    2 │ 1          moderate or severe  terbinafine   0.857143      2
    3 │ 1          moderate or severe  terbinafine   3.53571       3
  ⋮   │     ⋮              ⋮                ⋮           ⋮        ⋮
 1907 │ 383        none or mild        terbinafine   7.28571       5
 1908 │ 383        none or mild        terbinafine  10.7857        6
                                                    1903 rows omitted

julia> glmm1 = glme(@formula(outcome ~ visit + treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));

julia> glmm2 = glme(@formula(outcome ~ visit * treatment + (1|patientID)), toenail, Binomial(), LogitLink(), nAGQ = 20, wts = ones(Float64, size(toenail, 1)));

julia> glm1 = glm(@formula(outcome ~ visit + treatment), toenail, Binomial(), LogitLink());

julia> anova(glm1, glmm1, glmm2)
Analysis of Variance

Type 1 test / Likelihood-ratio test

Model 1: outcome ~ 1 + visit + treatment
Model 2: outcome ~ 1 + visit + treatment + (1 | patientID)
Model 3: outcome ~ 1 + visit + treatment + visit & treatment + (1 | patientID)

Table:
────────────────────────────────────────────────────
   DOF  ΔDOF  Res.DOF  Deviance        χ²  Pr(>|χ²|)
────────────────────────────────────────────────────
1    3           1906   1816.19
2    4     1     1904   1246.12  570.0715     <1e-99
3    5     1     1903   1242.39    3.7270     0.0535
────────────────────────────────────────────────────
```
Only likelihood-ratio test is available now.
## TO DO
1. More statitics will be printed to pressent more information. 
2. Implementation of `Rao` and `Mallow's Cp`.
3. `anova` for `GeneralizedLinearMixedModels` and `FixedEffectModels`.




