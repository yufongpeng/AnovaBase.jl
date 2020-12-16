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


@testset "LinearModel" begin
    lm0, lm1, lm2, lm3, lm4 = nestedmodels(lm(@formula(SepalLength ~ SepalWidth * Species), iris))
    aov3 = anova(lm4, type = 3)
    aov2 = anova_lm(@formula(SepalLength ~ SepalWidth * Species), iris, type = 2)
    aov1 = anova(lm4)
    aovf = anova(lm0, lm1, lm2, lm3, lm4)
    aovlr = anova(LRT, lm0, lm1, lm2, lm3, lm4)
    aov1lr = anova(LRT, lm4)
    @test aov1.stats.type == 1
    @test aov2.stats.nobs == 150
    @test all(aov3.stats.dof .== (1, 1, 2, 2, 144))
    @test isapprox(aov1.stats.deviance[1:end - 1], MixedAnova._diffn(aovf.stats.deviance))
    @test isapprox(filter(!isnan, aov1.stats.fstat), filter(!isnan, aovf.stats.fstat))
    @test isapprox(filter(!isnan, aov2.stats.fstat), (26485.300978452644, 56.637879295914324, 188.10911669921464, 0.4064420847481629))
    @test isapprox(filter(!isnan, aov1lr.stats.lrstat), filter(!isnan, aovlr.stats.lrstat))
    @test isapprox(coef(lm4), coef(aov3.model))
    @test all(coefnames(lm4, Val(:anova)) .==  ["(Intercept)", "SepalWidth", "Species", "SepalWidth & Species"])
end

@testset "LinearModel with frequency weights" begin
    wlm0, wlm1, wlm2 = nestedmodels(lm(@formula(PIQ ~ Brain), iq, wts = repeat([1/2], size(iq, 1))))
    aov = anova(wlm2)
    aovf = anova(wlm0, wlm1, wlm2)
    @test aov.stats.nobs == size(iq, 1) / 2
    @test all(aov.stats.dof .== (1, 1, 36))
    @test isapprox(filter(!isnan, aov.stats.fstat), filter(!isnan, aovf.stats.fstat))
    @test all(wlm0.model.rr.wts .== repeat([1/2], size(iq, 1)))
end
@testset "LinearMixedModel" begin
    @testset "One random effect on intercept" begin
        lmm0, lmm1, lmm2, lmm3, lmm4 = nestedmodels(lme(@formula(score ~ group * time + (1|id)), anxiety))
        aovf = anova(lmm4)
        aovlr = anova(lmm0, lmm1, lmm2, lmm3, lmm4)
        aovreml = anova_lme(@formula(score ~ group * time + (1|id)), anxiety, type = 3)
        @test aovf.stats.type == 1
        @test aovlr.stats.nobs == 135
        @test all(aovf.stats.dof .== (1, 2, 1, 2))
        @test all(aovf.stats.resdof .== (87, 42, 87, 87))
        @test isapprox(aovf.stats.fstat, (5019.394971941254, 4.455424160957743, 604.7793884113175, 159.32101757273566))
        @test isapprox(filter(!isnan, aovlr.stats.lrstat), (206.18217074629206, 8.474747188024821, 82.27167937164779, 139.37898838212664))
        @test aovreml.model.optsum.REML
        @test all(coefnames(lmm4, Val(:anova)) .== ["(Intercept)", "group", "time", "group & time"])
    end 
    """
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
    """
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
        @test all(deviance.(aov.model)[2:end] .== aov.stats.deviance)
        @test isapprox(aov.stats.lrstat, (1.1316862789786413, 9.278394637080218, 5.325918053163626, 18.783928942376562))
        @test isapprox(aov.stats.pval, (0.28741593819608247, 0.0023187256005545906, 0.021010537653650668, 1.4639554391712643e-5))
        @test isa(first(formula(first(aov.model)).rhs.terms), InterceptTerm{false})
    end

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

 