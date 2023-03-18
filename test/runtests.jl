using AnovaBase
import AnovaBase: dof_residual, nobs, anovatable, prednames, predictors, coeftable, dof_asgn, canonicalgoodnessoffit
using Distributions: Gamma, Binomial
import Base: show
using Test

struct TestTest <: GoodnessOfFit end
test_show(x) = show(IOBuffer(), x)
macro test_error(err, x)
    return quote
        try 
            $x
            false
        catch e
            isa(e, $err) ? true : false
        end
    end
end

macro test_error(x)
    return quote
        try 
            $x
            false
        catch e
            true
        end
    end
end

dof_residual(x::Int) = x
nobs(x) = ntuple(one, length(x))
nobs(x::Int) = one(x)
predictors(::Int) = tuple(Term.(Symbol.(["x$i" for i in 1:7]))...)
predictors(model::StatsModels.TableRegressionModel{Int64, Matrix{Float64}}) = formula(model).rhs
coeftable(model::StatsModels.TableRegressionModel{Int64, Matrix{Float64}}) = []
anovatable(::AnovaResult{<: FullModel{StatsModels.TableRegressionModel{Int64, Matrix{Float64}}}, LikelihoodRatioTest, 7}; rownames = string.(1:7)) = 
    AnovaBase.AnovaTable([
            [1, 1, 1, 1, 1, 1, 1], 
            [NaN, 1, 1e-5, 0.77, 7.7, 77.7, 777.7], 
            [1, 1, 1, 1, 1, 1, 1], 
            [NaN, 1, 1e-5, 0.77, 7.7, 77.7, 777.7], 
            [NaN, 1, 1e-5, 0.77, 7.7, 77.7, 777.7], 
            [NaN, 1, 1e-5, 0.77, 7.7, 77.7, 777.7],  
            [NaN, 1.0, 1e-5, 0.77e-4, 7.71e-4, 7.771e-3, 7.777e-2], 
        ],
    ["DOF", "ΔDOF", "Res.DOF", "Deviance", "ΔDeviance", "Χ²", "Pr(>|Χ²|)"],
    rownames, 7, 6)

@testset "AnovaBase.jl" begin
    global mm = ModelMatrix([1.0 1.0 1.0 1.0 1.0 1.0 1.0;], [1])
    global model1 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{false}(), Term.(Symbol.(["x$i" for i in 1:7]))...)), 
            1, 1, Int
            ), 
        mm)
    global model2 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{true}(), Term.(Symbol.(["x$i" for i in 1:7]))...)), 
            1, 1, Int), 
        mm)
    global model3 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{true}(),)), 
            1, 1, Int), 
        mm)
    global model4 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{false}(),)), 
            1, 1, Int), 
        mm)
    conterm = ContinuousTerm(:y, 0.0, 0.0, 0.0, 0.0)
    fterm = FunctionTerm(log, x->log(x), (:t, ), :(log(t)), [])
    caterm() = CategoricalTerm(:x, StatsModels.ContrastsMatrix(StatsModels.FullDummyCoding(), [1, 2, 3]))
    caterm(i) = CategoricalTerm(Symbol("x$i"), StatsModels.ContrastsMatrix(StatsModels.DummyCoding(), [1, 2, 3]))
    global model5 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{true}(), caterm(), caterm(1), InteractionTerm((caterm(), caterm(1))))), 
            1, 1, Int), 
        mm)
    global model6 = StatsModels.TableRegressionModel(
        1, 
        ModelFrame(
            FormulaTerm(Term(:y), tuple(InterceptTerm{true}(), caterm(), conterm, InteractionTerm((conterm, fterm)))), 
            1, 1, Int), 
        mm)
    @testset "AnovaBase.jl" begin
        @test FullModel(model1, 1, true, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 1, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 1, false, true).pred_id == ntuple(identity, 8)[3:8]
        @test FullModel(model1, 1, false, false).pred_id == ntuple(identity, 8)[3:8]
        @test FullModel(model1, 2, true, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 2, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 2, false, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 2, false, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 3, true, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 3, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 3, false, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model1, 3, false, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 1, true, true).pred_id == ntuple(identity, 8)
        @test FullModel(model2, 1, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 1, false, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 1, false, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 2, true, true).pred_id == ntuple(identity, 8)
        @test FullModel(model2, 2, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 2, false, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 2, false, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 3, true, true).pred_id == ntuple(identity, 8)
        @test FullModel(model2, 3, true, false).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 3, false, true).pred_id == ntuple(identity, 8)[2:8]
        @test FullModel(model2, 3, false, false).pred_id == ntuple(identity, 8)[2:8]
        @test @test_error ArgumentError FullModel(model3, 3, false, true)
        @test @test_error ArgumentError FullModel(model4, 3, true, true)
        @test FullModel(model5, 3, false, true).pred_id == ntuple(identity, 4)[2:4]
        @test FullModel(model6, 3, false, true).pred_id == ntuple(identity, 4)
        @test @test_error ArgumentError NestedModels{Int}(1.5, 2.5, 3)
        @test @test_error ArgumentError NestedModels{Int}((1.5, 2.5, 3))
        @test @test_error !(MixedAovModels{Number}(1, 2.5, 2//3))
        @test @test_error ArgumentError MixedAovModels{Int}((1, 2.5, 2//3))
    end
    global ft = AnovaResult{FTest}(NestedModels{StatsModels.TableRegressionModel}(
                                ntuple(7) do i
                                    StatsModels.TableRegressionModel(
                                        1, 
                                        ModelFrame(
                                            FormulaTerm(Term(:y), tuple(Term.(Symbol.(["x$i" for i in 1:i]))...)), 
                                            1, 1, Int), 
                                        mm)
                                end
                                ), 
                                ntuple(identity, 7), 
                                ntuple(one ∘ float, 7),
                                ntuple(one ∘ float, 7),
                                ntuple(zero ∘ float, 7),
                                NamedTuple())
    global lrt = AnovaResult{LRT}(FullModel(model1, ntuple(identity, 8)[2:8], 3), 
                                ntuple(identity, 7), 
                                ntuple(one ∘ float, 7),
                                ntuple(one ∘ float, 7),
                                ntuple(zero ∘ float, 7),
                                NamedTuple())
    global lrt2 = AnovaResult{LRT}(NestedModels{StatsModels.TableRegressionModel}(
                                ntuple(3) do i
                                    StatsModels.TableRegressionModel(
                                        1, 
                                        ModelFrame(
                                            FormulaTerm(Term(:y), tuple(Term.(Symbol.(["x$i" for i in 0:2i]))...)), 
                                            1, 1, Int), 
                                        mm)
                                end
                                ), 
                                ntuple(identity, 3), 
                                ntuple(one ∘ float, 3),
                                ntuple(one ∘ float, 3),
                                ntuple(zero ∘ float, 3),
                                NamedTuple())
    @testset "attr.jl" begin
        @test nobs(ft) == 1
        @test nobs(lrt) == 1
        @test dof_residual(lrt) == (1, 1, 1, 1, 1, 1, 1)
        @test dof_residual(lrt2) == (1, 1, 1)
        @test deviance(ft) == ft.deviance
        @test teststat(lrt) == lrt.teststat
        @test pval(ft) == ft.pval
        @test anova_test(lrt) == LRT
        @test anova_type(ft) == 1
        @test anova_type(lrt) == lrt.anovamodel.type
    end
    @testset "fit.jl" begin
        @test @test_error ErrorException nestedmodels(model1)
        @test @test_error ErrorException nestedmodels(typeof(model1), @formula(y~x1+x2+x3+x4+x5+x6+x7), [])
        @test @test_error ErrorException anova(model1; test = FTest)
        @test @test_error ErrorException anova(FTest, model1)
        @test @test_error ErrorException anova(LRT, model1, model1)
        @test AnovaBase._diff((1, 2, 3)) == (1, 1)
        @test AnovaBase._diffn((1, 2, 3)) == (-1, -1)
        @test canonicalgoodnessoffit(Gamma()) == FTest
        @test canonicalgoodnessoffit(Binomial()) == LRT
        @test AnovaBase.lrt_nested(NestedModels{StatsModels.TableRegressionModel}(model1, model1), (1,2), (1.5, 1.5), 0.1).teststat[2] == 0.0
        @test AnovaBase.ftest_nested(NestedModels{StatsModels.TableRegressionModel}((model1, model1)), (1,2), (10, 10), (1.5, 1.5), 0.1).teststat[2] == 0.0
        @test dof_asgn([1, 2, 2, 2, 3]) == (1, 3, 1)
    end
    global f = FormulaTerm(conterm, MatrixTerm((InterceptTerm{true}(), caterm(), fterm, InteractionTerm((caterm(), fterm)))))
    @testset "term.jl" begin
        @test AnovaBase.has_intercept(f.rhs)
        @test AnovaBase.any_not_aliased_with_1(f.rhs)
        @test !AnovaBase.any_not_aliased_with_1(MatrixTerm((InterceptTerm{true}(), caterm(), caterm(1), caterm(2), InteractionTerm((caterm(), caterm(1), caterm(2))))))
        @test AnovaBase.any_not_aliased_with_1(MatrixTerm((InterceptTerm{true}(), caterm(), conterm, fterm, InteractionTerm((caterm(), conterm, fterm)))))
        @test AnovaBase.select_super_interaction(f.rhs, 1) == Set([1, 2, 3, 4])
        @test AnovaBase.select_super_interaction(f.rhs, 2) == Set([2, 4])
        @test AnovaBase.select_not_super_interaction(f.rhs, 1) == Set(Int[])
        @test AnovaBase.select_not_super_interaction(f.rhs, 2) == Set([1, 3])
        @test AnovaBase.select_sub_interaction(f.rhs, 1) == Set(Int[])
        @test AnovaBase.select_sub_interaction(f.rhs, 2) == Set([1, 2])
        @test AnovaBase.select_not_sub_interaction(f.rhs, 1) == Set([1, 2, 3, 4])
        @test AnovaBase.select_not_sub_interaction(f.rhs, 2) == Set([3, 4])
        @test AnovaBase.subformula(f, 4; reschema = true).rhs[1:2] == @formula(y ~ 1 + x * log(t)).rhs[1:2]
        @test AnovaBase.subformula(f.lhs, f.rhs, 0; reschema = true).rhs[1] == @formula(y ~ 0).rhs
        @test AnovaBase.subformula(f.lhs, f.rhs, 0; reschema = true).rhs == AnovaBase.subformula(f.lhs, f.rhs, [1, 2, 3, 4], reschema = true).rhs
        @test AnovaBase.subformula(f.lhs, (f.rhs, f.rhs), 0) == FormulaTerm(f.lhs, (MatrixTerm((InterceptTerm{false}(),)), f.rhs))
        @test AnovaBase.subformula(f.lhs, (f.rhs, f.rhs), 0; rhs_id = 2, reschema = true).rhs[2] == FormulaTerm(f.lhs, (@formula(y ~ 1 + x * log(t)).rhs, @formula(y ~ 0).rhs)).rhs[2]
        @test keys(AnovaBase.extract_contrasts(f)) == keys(Dict(:x => StatsModels.FullDummyCoding()))
        @test prednames(f) == ["(Intercept)", "x", "log(t)", "x & log(t)"]
        @test prednames((f.lhs, f.rhs)) == ["y", "(Intercept)", "x", "log(t)", "x & log(t)"]
    end
    @testset "io.jl" begin
        @test prednames(lrt2)[3] == "x5+x6"
        @test @test_error ErrorException anovatable(AnovaResult{FullModel{Int64, 7}, FTest, 7}(
            FullModel(1, ntuple(identity, 7), 3),
            ntuple(identity, 7), 
            ntuple(one ∘ float, 7),
            ntuple(one ∘ float, 7),
            ntuple(zero ∘ float, 7),
            NamedTuple())
        )
        @test @test_error ErrorException anovatable(AnovaResult{NestedModels{Int, 7}, TestTest, 7}(
            NestedModels{Int}(ntuple(identity, 7)),
            ntuple(identity, 7), 
            ntuple(one ∘ float, 7),
            ntuple(one ∘ float, 7),
            ntuple(zero ∘ float, 7),
            NamedTuple())
        )
        @test @test_error ErrorException anovatable(AnovaResult{MixedAovModels{Number, 7}, TestTest, 7}(
            MixedAovModels{Number}(ntuple(identity, 7)),
            ntuple(identity, 7), 
            ntuple(one ∘ float, 7),
            ntuple(one ∘ float, 7),
            ntuple(zero ∘ float, 7),
            NamedTuple())
        )
        @test !(@test_error test_show(ft.anovamodel))
        @test !(@test_error test_show(lrt.anovamodel))
        @test !(@test_error test_show(ft))
        @test !(@test_error test_show(lrt))
        @test !(@test_error test_show(lrt2))
    end
end