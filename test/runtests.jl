using MixedAnova, CSV, RDatasets, DataFrames, CategoricalArrays
using Test
import Base.isapprox, MixedAnova.dof_residual

test_show(x) = show(IOBuffer(), x)
macro test_error(x)
    return quote
        try 
            $x
            false
        catch e
            @error e
            true
        end
    end
end

glm_init()
mm_init()
fem_init()

const anova_datadir = joinpath(dirname(@__FILE__), "..", "data")

iris = dataset("datasets", "iris")

qog18 = CSV.read(joinpath(anova_datadir, "qog18.csv"), DataFrame)

iq = CSV.read(joinpath(anova_datadir, "iqsize.csv"), DataFrame)

quine = dataset("MASS", "quine")

mtcars = dataset("datasets", "mtcars")
transform!(mtcars, [:VS, :Model] .=> categorical, renamecols = false)

sim = CSV.read(joinpath(anova_datadir, "poisson_sim.csv"), DataFrame)
transform!(sim, [:id, :prog] .=> categorical, renamecols = false)

gpa = CSV.read(joinpath(anova_datadir, "gpa.csv"), DataFrame)
transform!(gpa, 
    7 => x->replace(x, "yes" => true, "no" => false, "NA" => missing), 
    4 => x->categorical(x, levels = ["1 hour", "2 hours", "3 hours"], ordered = true),
    renamecols = false)
transform!(gpa, [1, 2, 5, 7] .=> categorical, renamecols = false)

nurses = CSV.read(joinpath(anova_datadir, "nurses.csv"), DataFrame)
transform!(nurses, 
    [1, 2, 3, 4, 6, 9, 11] .=> categorical, 
    10 => x->categorical(x, levels = ["small", "medium", "large"], ordered = true), 
    renamecols = false)

anxiety = CSV.read(joinpath(anova_datadir, "anxiety.csv"), DataFrame)
transform!(anxiety, :id => categorical, renamecols = false)

# custimized approx
isapprox(x::NTuple{N, Float64}, y::NTuple{N, Float64}, atol::NTuple{N, Float64} = x ./ 1000) where N = 
    all(map((a, b, c)->isapprox(a, b, atol = c > eps(Float64) ? c : eps(Float64)), x, y, atol))

dof_residual(x::Int) = x
@testset "LinearModel" begin 
    @testset "Simple linear regression" begin
        lm0, lm1, lm2, lm3, lm4 = nestedmodels(LinearModel, @formula(SepalLength ~ SepalWidth * Species), iris, dropcollinear = false)
        global aov3 = anova(lm4, type = 3)
        global aov2 = anova_lm(@formula(SepalLength ~ SepalWidth * Species), iris, type = 2)
        global aov1 = anova(lm4)
        global aovf = anova(lm0, lm1, lm2, lm3, lm4)
        global aovlr = anova(LRT, lm0, lm1, lm2, lm3, lm4)
        global aov1lr = anova(LRT, lm4)
        global aovlf = anova_lm(FTest, @formula(wdi_lifexp ~ log(gle_rgdpc) * ti_cpi), qog18, type = 2)
        ft = ftest(lm1.model, lm2.model, lm3.model, lm4.model)
        @test !(@test_error(test_show(aov1)))
        @test !(@test_error(test_show(aovf)))
        @test !(@test_error(test_show(aov1lr)))
        @test !(@test_error(test_show(aovlr)))
        @test !(@test_error(test_show(aovlf)))
        @test anova_type(aov1) == 1
        @test nobs(aov2) == ft.nobs
        @test dof(aov3) == (1, 1, 2, 2, 144)
        @test dof_residual(aov1) == 144
        @test isapprox(deviance(aovf)[2:end], ft.ssr)
        @test isapprox(deviance(aov1)[1:end - 1], MixedAnova._diffn(deviance(aovf)))
        @test isapprox(filter(!isnan, teststat(aov1)), filter(!isnan, teststat(aovf)))
        @test isapprox(filter(!isnan, teststat(aov2)), (26485.300978452644, 56.637879295914324, 188.10911669921464, 0.4064420847481629))
        @test isapprox(filter(!isnan, teststat(aov1lr)), filter(!isnan, teststat(aovlr)))
        @test isapprox(coef(lm4), coef(aov3.model))
        @test coefnames(lm4, Val(:anova)) ==  ["(Intercept)", "SepalWidth", "Species", "SepalWidth & Species"]
        @test collect(coefnames(formula(lm1), Val(:anova))) == coefnames((formula(lm1).lhs, formula(lm1).rhs), Val(:anova))

    end

    @testset "Linear regression with frequency weights" begin
        wlm1, wlm2 = nestedmodels(LinearModel, @formula(PIQ ~ Brain), iq, wts = repeat([1/2], size(iq, 1)))
        global aov = anova(wlm2)
        global aovf = anova(wlm1, wlm2)
        @test !(@test_error test_show(aov))
        @test !(@test_error test_show(aovf))
        @test nobs(aov) == size(iq, 1) / 2
        @test dof(aov) == (1, 1, 36)
        @test isapprox(filter(!isnan, teststat(aov))[2:end], filter(!isnan, teststat(aovf)))
        @test wlm1.model.rr.wts == repeat([1/2], size(iq, 1))
    end
end

@testset "GeneralizedLinearModel" begin
    @testset "Gamma regression" begin
        global aov = anova_glm(FTest, @formula(SepalLength ~ 0 + SepalWidth * Species), iris, Gamma(), type = 2)
        @test !(@test_error test_show(aov))
        @test nobs(aov) == nobs(aov.model)
        @test dof(aov) == (1, 3, 2, 144)
        @test isapprox(filter(!isnan, deviance(aov)), MixedAnova.deviances(aov.model, type = anova_type(aov)))
        @test isapprox(filter(!isnan, teststat(aov)), (62.53219308804019, 290.56265173997156, 0.4945888881568289))
        @test isapprox(filter(!isnan, pval(aov)), (6.278059283406317e-13, 7.345501482060886e-61, 0.610853639884267))
    end

    @testset "NegativeBinomial regression" begin
        global aov = anova_glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink(), type = 3)
        @test !(@test_error test_show(aov))
        @test nobs(aov) == nobs(aov.model)
        @test dof(aov) == (1, 1, 1, 3, 1, 139)
        @test anova_test(aov) == FTest
        @test isapprox(deviance(aov), MixedAnova.deviances(aov.model, type = 3))
        @test isapprox(filter(!isnan, teststat(aov)), (227.97422313423752, 13.180680587112887, 0.2840882132754838, 4.037856229143672, 2.6138558930314595))
        @test isapprox(filter(!isnan, pval(aov)), (4.251334236308285e-31, 0.00039667906264309605, 0.5948852219093326, 0.008659894572621351, 0.10820118567663468))
        @test coefnames(aov.model, Val(:anova)) == ["(Intercept)", "Eth", "Sex", "Age", "Lrn"]
    end

    @testset "Poisson regression" begin
        gmp = nestedmodels(GeneralizedLinearModel, @formula(num_awards ~ prog * math), sim, Poisson())
        global aov = anova(gmp...)
        lr = lrtest(gmp...)
        @test !(@test_error test_show(aov))
        @test first(nobs(aov)) == lr.nobs
        @test dof(aov) == lr.dof
        @test anova_test(aov) == LRT
        @test isapprox(deviance(aov), lr.deviance)
        @test isapprox(filter(!isnan, pval(aov)), filter(!isnan, lr.pval))
    end

    @testset "Logit regression" begin
        gml = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink())
        global aov = anova(gml)
        lr = lrtest(nestedmodels(gml)...)
        @test !(@test_error test_show(aov))
        @test nobs(aov) == lr.nobs
        @test dof(aov) == MixedAnova._diff(lr.dof)
        @test isapprox(deviance(aov), lr.deviance[2:end])
        @test isapprox(MixedAnova.deviances(aov.model)[2:end], MixedAnova._diffn((deviance(aov)..., 0.0)))
        @test isapprox(filter(!isnan, pval(aov)), filter(!isnan, lr.pval))
    end

    @testset "Probit regression" begin
        gmp0 = glm(@formula(AM ~ 1), mtcars, Binomial(), ProbitLink())
        gmp1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), ProbitLink())
        gmp2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), ProbitLink())
        global aov = anova(gmp0, gmp1, gmp2)
        lr = lrtest(gmp0, gmp1, gmp2)
        @test !(@test_error test_show(aov))
        @test first(nobs(aov)) == lr.nobs
        @test dof(aov) == lr.dof
        @test isapprox(deviance(aov), lr.deviance)
        @test isapprox(filter(!isnan, pval(aov)), filter(!isnan, lr.pval))
    end

    @testset "InverseGaussian regression" begin
        gmi = glm(@formula(SepalLength ~ SepalWidth * Species), iris, InverseGaussian())
        global aov = anova(gmi)
        @test !(@test_error test_show(aov))
        @test nobs(aov) == nobs(gmi)
        @test dof(aov) == (1, 2, 2, 144)
        @test isapprox(filter(!isnan, teststat(aov)), (8.172100334461327, 217.34941014323272, 1.8933247444892272))
        @test isapprox(filter(!isnan, pval(aov)), (0.004885873691352542, 3.202701170052312e-44, 0.15429963193830074))
    end
end

# test contrast
@testset "LinearMixedModel" begin
    @testset "One random effect on intercept" begin
        lmm1 = lme(@formula(gpa ~ occasion + sex + job + (1|student)), gpa)
        lmm2 = lme(@formula(gpa ~ occasion * sex + job + (1|student)), gpa)
        lmm3 = lme(@formula(gpa ~ occasion * sex * job + (1|student)), gpa)
        lme3 = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
        lme1 = anova_lme(LRT, @formula(score ~ group * time + (1|id)), anxiety)
        global aovf = anova(lmm1)
        global aovlr = anova(lmm1, lmm2, lmm3)
        global aov3 = anova(lmm1, type = 3)
        lr = MixedModels.likelihoodratiotest(lmm1, lmm2, lmm3)
        @test !(@test_error test_show(aovf))
        @test !(@test_error test_show(aovlr))
        @test !(@test_error test_show(aov3))
        @test anova_type(aov3) == 3
        @test first(nobs(aovlr)) == nobs(lmm1)
        @test dof(aovf) == (1, 1, 1, 2)
        @test dof_residual(aovf) == (997, 997, 198, 997)
        @test isapprox(teststat(aovf), (28899.81231026455, 714.3647106859062, 19.19804098862722, 45.425989531937134))
        @test isapprox(deviance(aovlr), tuple(lr.deviance...))
        @test isapprox(filter(!isnan, pval(aovlr)), tuple(lr.pvalues...))
        @test lme3.model.optsum.REML 
        @test !last(lme1.model).optsum.REML
        @test coefnames(lmm2, Val(:anova)) == ["(Intercept)",  "occasion", "sex", "job", "occasion & sex"]
    end 
    @testset "Random effect on slope and intercept" begin
        lmm1 = lme(@formula(gpa ~ occasion + (1|student)), gpa, REML = true)
        lmm2 = lme(@formula(gpa ~ occasion + (occasion|student)), gpa, REML = true)
        lms = nestedmodels(LinearMixedModel, @formula(gpa ~ occasion + (1|student)), gpa, REML = true)
        global aovlr = anova(lmm1, lmm2)
        global aovf = anova(lmm2)
        lr = MixedModels.likelihoodratiotest(lmm1, lmm2)
        @test !(@test_error test_show(aovlr))
        @test !(@test_error test_show(aovf))
        @test first(nobs(aovlr)) == nobs(lmm1)
        @test dof(aovlr) == tuple(lr.dof...)
        @test dof_residual(aovf) == (1000, 198)
        @test isapprox(deviance(aovlr), tuple(lr.deviance...))
        @test isapprox(filter(!isnan, pval(aovlr)), tuple(lr.pvalues...))
    end
    @testset "Multiple random effects" begin
        lmm1 = lme(@formula(stress ~ age  + sex + experience + treatment + wardtype + hospsize + (1|hospital)), nurses)
        lmm2 = lme(@formula(stress ~ age  + sex + experience + treatment + wardtype + hospsize + (1|hospital/wardid)), nurses)
        global aov = anova(lmm1, lmm2)
        lr = MixedModels.likelihoodratiotest(lmm1, lmm2)
        @test !(@test_error test_show(aov))
        @test first(nobs(aov)) == nobs(lmm1)
        @test dof(aov) == tuple(lr.dof...)
        @test isapprox(deviance(aov), tuple(lr.deviance...))
        @test isapprox(filter(!isnan, pval(aov)), tuple(lr.pvalues...))
    end
    #=
    @testset "Random effects on slope and intercept" begin
        aov1 = anova_lme(@formula(extro ~ open + agree + social + (1|school) + (1|class)), school, REML = true)
        aov2 = anova_lme(@formula(extro ~ open + agree + social + (open|school) + (open|class)), school, REML = true)
        @test aov1.stats.type == 1
        @test aov2.stats.nobs == 1200
        @test all(aov1.stats.dof .== (1, 1, 1, 1))
        @test all(aov2.stats.resdof .== (1192, 2, 1192, 1192))
        @test isapprox(aov1.stats.fstat, (208.372624411485, 1.673863453540288, 0.3223490172035519, 0.32161077338130784))
        @test isapprox(aov2.stats.pval, (1.5264108749941817e-43, 0.32443726673109696, 0.5622508640865989, 0.5681019984771887))
    end

    @testset "Nested random effects on intercept" begin
        aov = anova_lme(@formula(extro ~ open + agree + social + (1|school) + (1|school&class)), school)
        @test aov.stats.type == 1
        @test aov.stats.nobs == 1200
        @test all(aov.stats.dof .== (1, 1, 1, 1))
        @test all(aov.stats.resdof .== (1173, 1173, 1173, 1173))
        @test isapprox(aov.stats.fstat, (227.24904612916632, 1.4797389900418527, 1.8026931823189165, 0.08511022388038622))
        @test isapprox(aov.stats.pval, (4.508653887500674e-47, 0.22406007872324163, 0.1796470040402073, 0.7705396333261272))
    end

    @testset "Nested random effects on slope" begin
        aov = anova_lme(@formula(extro ~ open + agree + social + (open|school) + (open|school&class)), school)
        @test aov.stats.type == 1
        @test aov.stats.nobs == 1200
        @test all(aov.stats.dof .== (1, 1, 1, 1))
        @test all(aov.stats.resdof .== (1174, 4, 1174, 1174))
        @test isapprox(aov.stats.fstat, (250.00542522864583, 1.2322678515772565, 2.1135395635863543, 0.10258998684862923))
        @test isapprox(aov.stats.pval, (3.364459604379112e-51, 0.3292014599294774, 0.1462685388135273, 0.748800408618393))
    end
    =#
end

@testset "LinearModel and LinearMixedModel" begin
    lm1 = lm(@formula(score ~ group * time), anxiety)
    lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety)
    lmm2 = lme(@formula(score ~ group * time + (group|id)), anxiety)
    global aov = anova(lm1, lmm1, lmm2)
    lr = MixedModels.likelihoodratiotest(lm1, lmm1, lmm2)
    @test !(@test_error test_show(aov))
    @test first(nobs(aov)) == nobs(lmm1)
    @test dof(aov) == tuple(lr.dof...)
    @test isapprox(deviance(aov), tuple(lr.deviance...))
    @test isapprox(pval(aov)[2:end], tuple(lr.pvalues...))
end

@testset "Miscellaneous" begin
    global ft = AnovaResult{FTest}(ntuple(identity, 7), 
                            1, 
                            ntuple(identity, 7), 
                            ntuple(one ∘ float, 7),
                            ntuple(one ∘ float, 7),
                            ntuple(zero ∘ float, 7),
                            NamedTuple())
    global lrt = AnovaResult{LRT}(ntuple(identity, 7), 
                            1, 
                            ntuple(identity, 7), 
                            ntuple(one ∘ float, 7),
                            ntuple(one ∘ float, 7),
                            ntuple(zero ∘ float, 7),
                            NamedTuple())
    @test @test_error test_show(ft)
    @test @test_error test_show(lrt)
    @test @test_error formula(1)
    @test @test_error nestedmodels(1)
    @test @test_error anova(1)
    @test @test_error anova(FTest, 1)
    @test @test_error anova(LRT, 1)
    @test @test_error MixedAnova.isnullable(1)
end



 