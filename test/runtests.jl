using MixedAnova, CSV, RDatasets, DataFrames
using Test

const anova_datadir = joinpath(dirname(@__FILE__), "..", "data")

iris = dataset("datasets", "iris")

anxiety = CSV.read(joinpath(anova_datadir, "anxiety.csv"), DataFrame)
categorical!(anxiety, :id)

quine = dataset("MASS", "quine")

mtcars = dataset("datasets", "mtcars")
categorical!(mtcars, :Model)
categorical!(mtcars, :VS)

sim = CSV.read(joinpath(anova_datadir, "poisson_sim.csv"), DataFrame)
categorical!(sim, :prog)
categorical!(sim, :id)

@testset "One model" begin
    @testset "LinearModel" begin
        lm1 = lm(@formula(SepalLength ~ SepalWidth * Species), iris)
        aov3 = anova(lm1, type = 3)
        aov2 = anova_lm(@formula(SepalLength ~ SepalWidth * Species), iris, type = 2)
        aov1 = anova(lm1)
        @test aov1.stats.type == 1
        @test aov2.stats.nobs == 150
        @test all(aov3.stats.dof .== [1, 1, 2, 2, 144])
        @test all(isapprox.(filter(!isnan, aov1.stats.fstat), (528.3340970625617, 7.30297949463606, 188.10911669921464, 0.4064420847481629)))
        @test all(isapprox.(filter(!isnan, aov2.stats.fstat), (528.3340970625617, 56.637879295914324, 188.10911669921464, 0.4064420847481629)))
        @test all(isapprox.(filter(!isnan, aov3.stats.pval), (8.526118154774506e-6, 5.311035138434277e-5, 0.27896344508083537, 0.6667772844756328)))
        @test all(isapprox.(coef(lm1), coef(aov3.model)))
    end

    # add test for weighted linear model and w/o intercept

    @testset "LinearMixedModel" begin
        fm1 = lme(@formula(score ~ group * time + (1|id)), anxiety)
        aov1 = anova(fm1)
        aov2 = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
        aov3 = anova_lme(@formula(score ~ 0 + time + (1|id)), anxiety, type = 3)
        @test aov1.stats.type == 1
        @test aov2.stats.nobs == 135
        @test all(aov2.stats.ngroups .== (45,))
        @test all(aov1.stats.betweensubjects .== (0, 1, 0, 0))
        @test all(aov1.stats.dof .== (1, 2, 1, 2))
        @test all(isapprox.(filter(!isnan, aov1.stats.fstat), (5019.394971941254, 4.455424160957743, 604.7793884113175, 159.32101757273566)))
        @test all(isapprox.(filter(!isnan, aov3.stats.pval), (1.1117230614428629e-18,)))
        @test_broken all(aov1.stats.fstat .== aov2.stats.fstat) # not passing
        @warn "Problem to be fixed"
        @test aov2.model.optsum.REML
    end

    @testset "GeneralizedLinearModel" begin
        @testset "NegativeBinomial regression" begin
            aov = anova_glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink())
            @test length(aov.model) == 6
            @test aov.stats.nobs == 146
            @test all(aov.stats.dof .== (1, 1, 1, 3, 1))
            @test all(MixedAnova._diffn(deviance.(aov.model)) .== aov.stats.deviance)
            @test all(isapprox.(aov.stats.fstat, (2472.119426154856, 13.4251751124643, 1.9219804119281718, 3.2402799596242686, 2.614021047835751)))
            @test all(isapprox.(aov.stats.pval, (5.491419784674088e-93, 0.0003481860061298712, 0.16779623562803148, 0.024071499409229473, 0.10819014518459427))) 
        end

        @testset "Gamma regression" begin
            gmm = glm(@formula(SepalLength ~ SepalWidth * Species), iris, Gamma())
            aov = anova(gmm)
            @test length(aov.model) == 4
            @test aov.stats.nobs == 150
            @test all(aov.stats.dof .== (1, 1, 2, 2))
            @test all(isapprox.(filter(!isnan, aov.stats.deviance), MixedAnova._diffn(deviance.(aov.model))))
            @test all(isapprox.(filter(!isnan, aov.stats.fstat), (7.955957684363591, 209.23937708735653, 0.49461923218008963)))
            @test all(isapprox.(filter(!isnan, aov.stats.pval), (0.005450780832897245, 1.3411693736689225e-43, 0.6108352308672702)))
        end

        @testset "Logit regression" begin
            gmb = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink())
            aov = anova(gmb)
            @test length(aov.model) == 5
            @test aov.stats.nobs == 32
            @test all(aov.stats.dof .== (1, 1, 1, 1))
            @test all(MixedAnova._diffn(deviance.(aov.model)) .== aov.stats.deviance)
            @test all(isapprox.(aov.stats.lrstat, (1.1316862789786413, 9.278394637080218, 5.325918053163626, 18.783928942376562)))
            @test all(isapprox.(aov.stats.pval, (0.28741593819608247, 0.0023187256005545906, 0.021010537653650668, 1.4639554391712643e-5)))
        end
    end
end

@testset "Nestedm models" begin
    @testset "LinearModels" begin
        lm0 = lm(@formula(SepalLength ~ 1), iris)
        lm1 = lm(@formula(SepalLength ~ SepalWidth), iris)
        lm2 = lm(@formula(SepalLength ~ SepalWidth + Species), iris)
        lm3 = lm(@formula(SepalLength ~ SepalWidth * Species), iris)
        aov = anova(lm0, lm1, lm2, lm3)
        @test aov.stats.nobs == 150
        @test all(aov.stats.dof .== (2, 3, 5, 7))
        @test all(deviance.(aov.model) .== aov.stats.deviance)
        @test all(isapprox.(filter(!isnan, aov.stats.fstat), (7.302979494635763, 188.10911669921464, 0.40644208474813515)))
        @test all(isapprox.(filter(!isnan, aov.stats.pval), (0.0076887684006947754, 3.9315172882928e-41, 0.6667772844756514)))
    end

    @testset "LinearMixedModels" begin
        fm0 = lme(@formula(score ~ time + (1|id)), anxiety)
        fm1 = lme(@formula(score ~ group + time + (1|id)), anxiety)
        fm2 = lme(@formula(score ~ group * time + (1|id)), anxiety)
        aov = anova(fm0, fm1, fm2)
        @test aov.stats.nobs == 135
        @test all(aov.stats.dof .== (4, 6, 8))
        @test all(isapprox.(filter(!isnan, aov.stats.lrstat), (8.474747188077288, 139.37898838211368)))
        @test all(isapprox.(filter(!isnan, aov.stats.pval), (0.014445481763593991, 5.42297030318914e-31)))
    end

    @testset "GeneralizedLinearModel" begin
        @testset "Poisson regression" begin
            gmp0 = glm(@formula(num_awards ~ 1), sim, Poisson())
            gmp1 = glm(@formula(num_awards ~ prog), sim, Poisson())
            gmp2 = glm(@formula(num_awards ~ prog + math), sim, Poisson())
            gmp3 = glm(@formula(num_awards ~ prog * math), sim, Poisson())
            aov = anova(gmp0, gmp1, gmp2, gmp3)
            @test aov.stats.nobs == 200
            @test all(aov.stats.dof .== (1, 3, 4, 6))
            @test all(isapprox.(filter(!isnan, aov.stats.lrstat), (53.21226085915595, 45.0103536834161, 0.3480013904814143)))
            @test all(isapprox.(filter(!isnan, aov.stats.pval), (2.7867908096413152e-12, 1.9599542632566133e-11, 0.8402963134500252)))
        end

        @testset "Probit regression" begin
            gmp0 = glm(@formula(AM ~ 1), mtcars, Binomial(), ProbitLink())
            gmp1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), ProbitLink())
            gmp2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), ProbitLink())
            aov = anova(gmp0, gmp1, gmp2)
            @test aov.stats.nobs == 32
            @test all(aov.stats.dof .== (1, 4, 8))
            @test all(isapprox.(filter(!isnan, aov.stats.lrstat), (33.63555518087605, 9.594177758919395)))
            @test all(isapprox.(filter(!isnan, aov.stats.pval), (2.3651087229710803e-7, 0.047847662863957405)))
        end

        @testset "InverseGaussian regression" begin
            gmi0 = glm(@formula(SepalLength ~ 1), iris, InverseGaussian())
            gmi1 = glm(@formula(SepalLength ~ SepalWidth), iris, InverseGaussian())
            gmi2 = glm(@formula(SepalLength ~ SepalWidth + Species), iris, InverseGaussian())
            gmi3 = glm(@formula(SepalLength ~ SepalWidth * Species), iris, InverseGaussian())
            aov = anova(gmi0, gmi1, gmi2, gmi3)
            @test aov.stats.nobs == 150
            @test all(aov.stats.dof .== (2, 3, 5, 7))
            @test all(isapprox.(filter(!isnan, aov.stats.fstat), (8.172100334461327, 217.34941014323272, 1.8933247444892272)))
            @test all(isapprox.(filter(!isnan, aov.stats.pval), (0.004868229756510584, 1.6956631363452993e-44, 0.1542995990843151)))
        end
    end
end
 