using AnovaBase
import AnovaBase: dof_residual, nobs, anovatable, factornames, formula
using Distributions: Gamma, Binomial
import Base: show
using Test

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

dof_residual(x::Int) = x
nobs(x) = ntuple(one, length(x))
nobs(x::Int) = one(x)
anovatable(::AnovaResult{StatsModels.TableRegressionModel{Int64, Matrix{Float64}}, LikelihoodRatioTest, 1}) = 
    AnovaBase.AnovaTable(hcat(AnovaBase.vectorize.((
        [1, 1, 1, 1, 1], 
        [NaN, 1, 1e-5, 7.7, 77.7], 
        [1, 1, 1, 1, 1], 
        [NaN, 1, 1e-5, 7.7, 77.7], 
        [NaN, 1, 1e-5, 7.7, 77.7], 
        [NaN, 1, 1e-5, 7.7, 77.7],  
        [NaN, 1.0, 1e-5, 1e-5, 1e-5], 
        ))...),
    ["DOF", "ΔDOF", "Res.DOF", "Deviance", "ΔDeviance", "Χ²", "Pr(>|Χ²|)"],
    ["$i" for i in eachindex([1, 1, 1, 1, 1])], 7, 6)

factornames(::StatsModels.TableRegressionModel{Int64, Matrix{Float64}}) = ["1", "2", "3", "4", "5"]
formula(::Int) = 1

@testset "AnovaBase.jl" begin
    global ft = AnovaResult{FTest}(ntuple(identity, 7), 
                                1, 
                                ntuple(identity, 7), 
                                ntuple(one ∘ float, 7),
                                ntuple(one ∘ float, 7),
                                ntuple(zero ∘ float, 7),
                                NamedTuple())
    global model = StatsModels.TableRegressionModel(1, ModelFrame(@formula(y~x), 1, 1, Int), ModelMatrix([1.0 1.0;], [1]))
    global lrt = AnovaResult{LRT}(model, 
                                1, 
                                ntuple(identity, 1), 
                                ntuple(one ∘ float, 1),
                                ntuple(one ∘ float, 1),
                                ntuple(zero ∘ float, 1),
                                NamedTuple())
    global lrt2 = AnovaResult{LRT}(ntuple(identity, 7), 
                                1, 
                                ntuple(identity, 7), 
                                ntuple(one ∘ float, 7),
                                ntuple(one ∘ float, 7),
                                ntuple(zero ∘ float, 7),
                                NamedTuple())
    @testset "api.jl" begin
        @test @test_error formula(nothing)
        @test @test_error nestedmodels(1)
        @test @test_error anova(1)
        @test @test_error anova(FTest, 1)
        @test @test_error anova(LRT, 1)
        @test nobs(ft) == (1, 1, 1, 1, 1, 1, 1)
        @test nobs(lrt) == 1
        @test deviance(ft) == ft.deviance
        @test teststat(lrt) == lrt.teststat
        @test pval(ft) == ft.pval
        @test anova_test(lrt) == LRT
        @test anova_type(ft) == ft.type

    end
    @testset "fit.jl" begin
        @test AnovaBase._diff((1,2,3)) == (1, 1)
        @test AnovaBase._diffn((1,2,3)) == (-1, -1)
        @test canonicalgoodnessoffit(Gamma()) == FTest
        @test canonicalgoodnessoffit(Binomial()) == LRT
        @test AnovaBase.lrt_nested((model, model), (1,2), (1.5, 1.5), 0.1).teststat[2] == 0.0
        @test AnovaBase.ftest_nested((model, model), (1,2), (10, 10), (1.5, 1.5), 0.1).teststat[2] == 0.0
        @test dof_asgn([1,2,2,2,3]) == [1, 3, 1]
    end
    fterm = FunctionTerm(log, x->log(x), (:t, ), :(log(t)), [])
    caterm = CategoricalTerm(:x, StatsModels.ContrastsMatrix(StatsModels.FullDummyCoding(), [1,2,3]))
    global f = FormulaTerm(ContinuousTerm(:y, 0.0, 0.0, 0.0, 0.0), MatrixTerm((InterceptTerm{true}(), caterm, fterm, InteractionTerm((caterm, fterm)))))
    @testset "termIO.jl" begin
        @test factornames(f) == ("y", ["(Intercept)", "x", "log(t)", "x & log(t)"])
        @test factornames(StatsModels.TupleTerm((f.lhs, f.rhs))) == ["y", "(Intercept)", "x", "log(t)", "x & log(t)"]
        @test AnovaBase.select_super_interaction(f.rhs, 1) == Set([1, 2, 3, 4])
        @test AnovaBase.select_super_interaction(f.rhs, 2) == Set([2, 4])
        @test AnovaBase.select_not_super_interaction(f.rhs, 1) == Set(Int[])
        @test AnovaBase.select_not_super_interaction(f.rhs, 2) == Set([1, 3])
        @test AnovaBase.select_sub_interaction(f.rhs, 1) == Set(Int[])
        @test AnovaBase.select_sub_interaction(f.rhs, 2) == Set([1, 2])
        @test AnovaBase.select_not_sub_interaction(f.rhs, 1) == Set([1, 2, 3, 4])
        @test AnovaBase.select_not_sub_interaction(f.rhs, 2) == Set([3, 4])
        @test AnovaBase.subformula(f.lhs, f.rhs, 4, reschema = true).rhs[1:2] == @formula(y ~ 1 + x * log(t)).rhs[1:2]
        @test AnovaBase.subformula(f.lhs, f.rhs, 0, reschema = true).rhs[1] == @formula(y ~ 0).rhs
        @test AnovaBase.subformula(f.lhs, f.rhs, 0, reschema = true).rhs == AnovaBase.subformula(f.lhs, f.rhs, [1, 2, 3, 4], reschema = true).rhs
        @test AnovaBase.subformula(f.lhs, (f.rhs, f.rhs), 0) == FormulaTerm(f.lhs, (MatrixTerm((InterceptTerm{false}(),)), f.rhs))
        @test keys(AnovaBase.extract_contrasts(f)) == keys(Dict(:x => StatsModels.FullDummyCoding()))
        @test !(@test_error test_show(ft))
        @test !(@test_error test_show(lrt))
        @test !(@test_error test_show(lrt2))
    end
end