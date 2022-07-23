using AnovaBase
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

@testset "AnovaBase.jl" begin
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