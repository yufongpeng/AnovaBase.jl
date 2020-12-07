using MixedAnova, CSV, RDatasets, DataFrames
using Test
import Base: isapprox

const anova_datadir = joinpath(dirname(@__FILE__), "..", "data")

iris = dataset("datasets", "iris")

iq = CSV.read(joinpath(anova_datadir, "iqsize.csv"), DataFrame)

anxiety = CSV.read(joinpath(anova_datadir, "anxiety.csv"), DataFrame)
transform!(anxiety, :id => categorical, renamecols = false)

school = CSV.read(joinpath(anova_datadir, "school.csv"), DataFrame)
transform!(school, :id => categorical, renamecols = false)

oxide = CSV.read(joinpath(anova_datadir, "oxide.csv"), DataFrame)
for i in 1:4
    transform!(oxide, names(oxide, i) => categorical, renamecols = false)
end

quine = dataset("MASS", "quine")

mtcars = dataset("datasets", "mtcars")
transform!(mtcars, :Model => categorical, renamecols = false)
transform!(mtcars, :VS => categorical, renamecols = false)

sim = CSV.read(joinpath(anova_datadir, "poisson_sim.csv"), DataFrame)
transform!(sim, :prog => categorical, renamecols = false)
transform!(sim, :id => categorical, renamecols = false)

# custimize approx
isapprox(x::NTuple{N, Number}, y::NTuple{N, Number}, atol::NTuple{N, Number} = x ./ 1000) where N = 
    all(map((a, b, c)->isapprox(a, b, atol = c), x, y, atol))


println()
@info "One model"
@testset "LinearModel" begin
    lm1 = lm(@formula(SepalLength ~ SepalWidth * Species), iris)
    aov3 = anova(lm1, type = 3)
    aov2 = anova_lm(@formula(SepalLength ~ SepalWidth * Species), iris, type = 2)
    aov1 = anova(lm1)
    @test aov1.stats.type == 1
    @test aov2.stats.nobs == 150
    @test all(aov3.stats.dof .== (1, 1, 2, 2, 144))
    @test isapprox(filter(!isnan, aov1.stats.fstat), (528.3340970625617, 7.30297949463606, 188.10911669921464, 0.4064420847481629))
    @test isapprox(filter(!isnan, aov2.stats.fstat), (528.3340970625617, 56.637879295914324, 188.10911669921464, 0.4064420847481629))
    @test isapprox(filter(!isnan, aov3.stats.pval), (8.526118154774506e-6, 5.311035138434277e-5, 0.27896344508083537, 0.6667772844756328))
    @test isapprox(coef(lm1), coef(aov3.model))
    @test all(coefnames(lm1, Val(:anova)) .==  ["(Intercept)", "SepalWidth", "Species", "SepalWidth & Species", "(Residual)"])
end

@testset "LinearModel with frequency weights" begin
    wlm = lm(@formula(PIQ ~ Brain + Height + Weight), iq, wts = repeat([1/2], size(iq, 1)))
    aov = anova(wlm, type = 3)
    @test aov.stats.nobs == 38
    @test all(aov.stats.dof .== (1, 1, 1, 1, 34))
    @test isapprox(filter(!isnan, aov.stats.fstat), (3.1269870172690166, 13.371594967250223, 4.937785526448175, 8.073364237764428e-6))
    @test isapprox(filter(!isnan, aov.stats.pval), (0.0859785407765702, 0.0008556321566144881, 0.03303381869828931, 0.9977495267695052))
end
@testset "LinearMixedModel" begin
    @testset "One random effect on intercept" begin
        lmm1 = lme(@formula(score ~ group * time + (1|id)), anxiety)
        aov1 = anova(lmm1)
        aov2 = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
        aov3 = anova_lme(@formula(score ~ 0 + time + (1|id)), anxiety, type = 3)
        @test aov1.stats.type == 1
        @test aov2.stats.nobs == 135
        @test all(aov1.stats.dof .== (1, 2, 1, 2))
        @test all(aov2.stats.resdof .== (87, 42, 87, 87))
        @test isapprox(aov1.stats.fstat, (5019.394971941254, 4.455424160957743, 604.7793884113175, 159.32101757273566))
        @test isapprox(aov3.stats.pval, (1.1117230614428629e-18,))
        @test aov2.model.optsum.REML
        @test all(coefnames(lmm1, Val(:anova)) .== ["(Intercept)", "group", "time", "group & time"])
    end 

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
end

@testset "GeneralizedLinearModel" begin
    @testset "NegativeBinomial regression" begin
        aov = anova_glm(@formula(Days ~ Eth + Sex + Age + Lrn), quine, NegativeBinomial(2.0), LogLink())
        @test length(aov.model) == 6
        @test aov.stats.nobs == 146
        @test all(aov.stats.dof .== (1, 1, 1, 3, 1))
        @test all(MixedAnova._diffn(deviance.(aov.model)) .== aov.stats.deviance)
        @test isapprox(aov.stats.fstat, (2472.119426154856, 13.4251751124643, 1.9219804119281718, 3.2402799596242686, 2.614021047835751))
        @test isapprox(aov.stats.pval, (5.491419784674088e-93, 0.0003481860061298712, 0.16779623562803148, 0.024071499409229473, 0.10819014518459427))
        @test all(coefnames(last(aov.model), Val(:anova)) .== ["(Intercept)", "Eth", "Sex", "Age", "Lrn"])
        @test isa(first(formula(first(aov.model)).rhs.terms), InterceptTerm{false})
    end

    @testset "Gamma regression" begin
        gmm = glm(@formula(SepalLength ~ SepalWidth * Species), iris, Gamma())
        aov = anova(gmm)
        @test length(aov.model) == 4
        @test aov.stats.nobs == 150
        @test all(aov.stats.dof .== (1, 1, 2, 2))
        @test isapprox(filter(!isnan, aov.stats.deviance), MixedAnova._diffn(deviance.(aov.model)))
        @test isapprox(filter(!isnan, aov.stats.fstat), (7.955957684363591, 209.23937708735653, 0.49461923218008963))
        @test isapprox(filter(!isnan, aov.stats.pval), (0.005450780832897245, 1.3411693736689225e-43, 0.6108352308672702))
        @test isa(first(formula(first(aov.model)).rhs.terms), InterceptTerm{true})
    end

    @testset "Logit regression" begin
        gmb = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), LogitLink())
        aov = anova(gmb)
        @test length(aov.model) == 5
        @test aov.stats.nobs == 32
        @test all(aov.stats.dof .== (1, 1, 1, 1))
        @test all(MixedAnova._diffn(deviance.(aov.model)) .== aov.stats.deviance)
        @test isapprox(aov.stats.lrstat, (1.1316862789786413, 9.278394637080218, 5.325918053163626, 18.783928942376562))
        @test isapprox(aov.stats.pval, (0.28741593819608247, 0.0023187256005545906, 0.021010537653650668, 1.4639554391712643e-5))
        @test isa(first(formula(first(aov.model)).rhs.terms), InterceptTerm{false})
    end
end

println()
@info "Nested models"

@testset "LinearModels" begin
    lm0 = lm(@formula(SepalLength ~ 1), iris)
    lm1 = lm(@formula(SepalLength ~ SepalWidth), iris)
    lm2 = lm(@formula(SepalLength ~ SepalWidth + Species), iris)
    lm3 = lm(@formula(SepalLength ~ SepalWidth * Species), iris)
    aov = anova(lm0, lm1, lm2, lm3)
    @test aov.stats.nobs == 150
    @test all(aov.stats.dof .== (2, 3, 5, 7))
    @test all(deviance.(aov.model) .== aov.stats.deviance)
    @test isapprox(filter(!isnan, aov.stats.fstat), (7.302979494635763, 188.10911669921464, 0.40644208474813515))
    @test isapprox(filter(!isnan, aov.stats.pval), (0.0076887684006947754, 3.9315172882928e-41, 0.6667772844756514))
end

@testset "LinearMixedModels" begin
    fm0 = lme(@formula(score ~ time + (1|id)), anxiety)
    fm1 = lme(@formula(score ~ group + time + (1|id)), anxiety)
    fm2 = lme(@formula(score ~ group * time + (1|id)), anxiety)
    aov = anova(fm0, fm1, fm2)
    @test aov.stats.nobs == 135
    @test all(aov.stats.dof .== (4, 6, 8))
    @test isapprox(filter(!isnan, aov.stats.lrstat), (8.474747188077288, 139.37898838211368))
    @test isapprox(filter(!isnan, aov.stats.pval), (0.014445481763593991, 5.42297030318914e-31))
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
        @test isapprox(filter(!isnan, aov.stats.lrstat), (53.21226085915595, 45.0103536834161, 0.3480013904814143))
        @test isapprox(filter(!isnan, aov.stats.pval), (2.7867908096413152e-12, 1.9599542632566133e-11, 0.8402963134500252))
    end

    @testset "Probit regression" begin
        gmp0 = glm(@formula(AM ~ 1), mtcars, Binomial(), ProbitLink())
        gmp1 = glm(@formula(AM ~ Cyl + HP + WT), mtcars, Binomial(), ProbitLink())
        gmp2 = glm(@formula(AM ~ Cyl * HP * WT), mtcars, Binomial(), ProbitLink())
        aov = anova(gmp0, gmp1, gmp2)
        @test aov.stats.nobs == 32
        @test all(aov.stats.dof .== (1, 4, 8))
        @test isapprox(filter(!isnan, aov.stats.lrstat), (33.63555518087605, 9.594177758919395))
        @test isapprox(filter(!isnan, aov.stats.pval), (2.3651087229710803e-7, 0.047847662863957405))
    end

    @testset "InverseGaussian regression" begin
        gmi0 = glm(@formula(SepalLength ~ 1), iris, InverseGaussian())
        gmi1 = glm(@formula(SepalLength ~ SepalWidth), iris, InverseGaussian())
        gmi2 = glm(@formula(SepalLength ~ SepalWidth + Species), iris, InverseGaussian())
        gmi3 = glm(@formula(SepalLength ~ SepalWidth * Species), iris, InverseGaussian())
        aov = anova(gmi0, gmi1, gmi2, gmi3)
        @test aov.stats.nobs == 150
        @test all(aov.stats.dof .== (2, 3, 5, 7))
        @test isapprox(filter(!isnan, aov.stats.fstat), (8.172100334461327, 217.34941014323272, 1.8933247444892272))
        @test isapprox(filter(!isnan, aov.stats.pval), (0.004868229756510584, 1.6956631363452993e-44, 0.1542995990843151))
    end
end

 