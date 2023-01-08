# Function related to I/O
"""
    prednames(aov::AnovaResult)
    prednames(anovamodel::FullModel) 
    prednames(anovamodel::NestedModels)
    prednames(<model>)

Return the name of predictors as a vector of strings.
When there are multiple models, return value is `nothing`.
"""
prednames(aov::AnovaResult) = prednames(aov.anovamodel)
prednames(anovamodel::FullModel) = collect(prednames.(getindex.(Ref(predictors(anovamodel.model)), anovamodel.pred_id)))
function prednames(anovamodel::NestedModels)
    names = collect(prednames.(anovamodel.model))
    names[2:end] = [setdiff(b, a) for (a, b) in @views zip(names[1:end - 1], names[2:end])]
    join.(names, "+")
end
prednames(model::RegressionModel) = vectorize(prednames(predictors(model)))

@deprecate coefnames(aov::AnovaResult) prednames(aov::AnovaResult)
@deprecate coefnames(x, ::Val{:anova}) prednames(x)

"""
    prednames(<term>)

Return the name(s) of predictor(s). Return value is either a `String`, an iterable of `String`s or `nothing`.

# Examples
```julia
julia> iris = dataset("datasets", "iris");

julia> f = formula(lm(@formula(log(SepalLength) ~ SepalWidth + PetalLength * PetalWidth), iris))
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  SepalWidth(continuous)
  PetalLength(continuous)
  PetalWidth(continuous)
  PetalLength(continuous) & PetalWidth(continuous)

julia> prednames(f)
["(Intercept)", "SepalWidth", "PetalLength", "PetalWidth", "PetalLength & PetalWidth"]

julia> prednames(InterceptTerm{false}())

```
"""
prednames(t::FormulaTerm) = filter!(!isnothing, vectorize(prednames(t.rhs)))
prednames(t::MatrixTerm) = prednames(t.terms)
prednames(ts::StatsModels.TupleTerm) = filter!(!isnothing, vectorize(mapreduce(prednames, vcat, ts)))
prednames(::InterceptTerm{H}) where H = H ? "(Intercept)" : nothing
prednames(t::ContinuousTerm) = string(t.sym)
prednames(t::CategoricalTerm) = string(t.sym)
prednames(t::FunctionTerm) = string(t.exorig)
prednames(t::InteractionTerm) = join(prednames.(t.terms), " & ")
prednames(t::Term) = string(t)
prednames(t::ConstantTerm{H}) where H = string(t)
prednames(t) = coefnames(t)

# test name
"""
    testname(::Type{FTest}) = "F test"
    testname(::Type{LRT}) = "Likelihood-ratio test"

Name of tests.
"""
testname(::Type{FTest}) = "F test"
testname(::Type{LRT}) = "Likelihood-ratio test"
#testname(M::AnovaStatsRao) = "Rao score test"
#testname(M::AnovaStatsCp) = "Mallow's Cp"
@deprecate tname testname

# ============================================================================================================================
# AnovaTable, mostly from CoefTable
"""
    AnovaTable

A table with coefficients and related statistics of ANOVA. It is mostly modified from `StatsModels.CoefTable`.

# Fields
* `cols`: values of each statiscics.
* `colnms`: names of statiscics.
* `rownms`: names of each row.
* `pvalcol`: the index of column repressenting p-value.
* `teststatcol`: the index of column representing test statiscics.

# Constructor
    AnovaTable(cols::Vector, colnms::Vector, rownms::Vector, pvalcol::Int = 0, teststatcol::Int = 0)
    AnovaTable(mat::Matrix, colnms::Vector, rownms::Vector, pvalcol::Int = 0, teststatcol::Int = 0)
"""
mutable struct AnovaTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    pvalcol::Int
    teststatcol::Int
    function AnovaTable(cols::Vector, colnms::Vector, rownms::Vector,
                       pvalcol::Int = 0, teststatcol::Int = 0)
        nc = length(cols)
        nrs = map(length, cols)
        nr = nrs[1]
        length(colnms) in [0, nc] || throw(ArgumentError("colnms should have length 0 or $nc"))
        length(rownms) in [0, nr] || throw(ArgumentError("rownms should have length 0 or $nr"))
        all(nrs .== nr) || throw(ArgumentError("Elements of cols should have equal lengths, but got $nrs"))
        pvalcol in 0:nc || throw(ArgumentError("pvalcol should be between 0 and $nc"))
        teststatcol in 0:nc || throw(ArgumentError("teststatcol should be between 0 and $nc"))
        new(cols, colnms, rownms, pvalcol, teststatcol)
    end
end
AnovaTable(mat::Matrix, colnms::Vector, rownms::Vector, pvalcol::Int = 0, teststatcol::Int = 0) = 
        AnovaTable(Any[mat[:, i] for i in 1:size(mat, 2)], colnms, rownms, pvalcol, teststatcol)

"""
Show a p-value using 6 characters, either using the standard 0.XXXX
representation or as <Xe-YY.
"""
struct PValue
    v::Real
    function PValue(v::Real)
        0 <= v <= 1 || isnan(v) || error("p-values must be in [0; 1]")
        new(v)
    end
end

function show(io::IO, pv::PValue)
    v = pv.v
    if isnan(v)
        print(io,"")
    elseif v >= 1e-4
        @printf(io,"%.4f", v)
    else
        @printf(io,"<1e%2.2d", ceil(Integer, max(nextfloat(log10(v)), -99)))
    end
end

"""Show a test statistic using 2 decimal digits"""
struct TestStat <: Real
    v::Real
end

show(io::IO, x::TestStat) = isnan(x.v) ? print(io,"") : @printf(io, "%.4f", x.v)

"""Wrap a string so that show omits quotes"""
struct NoQuote
    s::String
end

show(io::IO, n::NoQuote) = print(io, n.s)


"""Filter NaN"""
struct OtherStat <: Real
    v::Real
end

function show(io::IO, x::OtherStat)
    v = x.v
    if isnan(v) 
        print(io, "")
    elseif floor(v) == v
        print(io, Int(v))
    elseif v < 1e-4 && v > 0
        @printf(io,"<1e%2.2d", ceil(Integer, max(nextfloat(log10(v)), -99)))
    elseif abs(v) < 10
        @printf(io,"%.4f", v)
    else 
        @printf(io,"%.2f", v)
    end
end

function show(io::IO, at::AnovaTable)
    cols = at.cols; rownms = at.rownms; colnms = at.colnms;
    nc = length(cols)
    nr = length(cols[1])
    if length(rownms) == 0
        rownms = [lpad("[$i]",floor(Integer, log10(nr)) + 3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j - 1 == at.pvalcol ? PValue(cols[j - 1][i]) :
           j - 1 in at.teststatcol ? TestStat(cols[j - 1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j - 1][i]) : OtherStat(cols[j - 1][i])
           for i in 1:nr, j in 1:nc + 1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1] + nmswidths[i] - sum(A[i]), A[i][2]) : A[i]
         for i in eachindex(A)]
    totwidth = sum(sum.(A)) + 2 * (length(A) - 1)
    println(io, repeat('─', totwidth))
    print(io, repeat(' ', sum(A[1])))
    for j in 1:length(colnms)
        print(io, "  ", lpad(colnms[j], sum(A[j + 1])))
    end
    println(io, '\n', repeat('─', totwidth))
    for i in 1:size(mat, 1)
        Base.print_matrix_row(io, mat, A, i, 1:size(mat, 2), "  ")
        i != size(mat, 1) && println(io)
    end
    print(io, '\n', repeat('─', totwidth))
    nothing
end

# ====================================================================================================================================
# AnovaTable api
"""
    anovatable(aov::AnovaResult{<: FullModel, Test}; rownames = prednames(aov))
    anovatable(aov::AnovaResult{<: NestedModels, Test}; rownames = string.(1:N))
    anovatable(aov::AnovaResult{<: NestedModels, FTest, N}; rownames = string.(1:N)) where N
    anovatable(aov::AnovaResult{<: NestedModels, LRT, N}; rownames = string.(1:N)) where N

Return a table with coefficients and related statistics of ANOVA.
When displaying `aov` in repl, `rownames` will be `prednames(aov)` for `FullModel` and `string.(1:N)` for `NestedModels`. 

For nested models, there are two default methods for `FTest` and `LRT`; one can also define new methods dispatching on `::NestedModels{M}` where `M` is a model type. 

For a single model, no default api is implemented.

The returned `AnovaTable` object implements the [`Tables.jl`](https://github.com/JuliaData/Tables.jl/) interface, and can be  
converted e.g. to a DataFrame via `using DataFrames; DataFrame(anovatable(aov))`.
"""
function anovatable(aov::AnovaResult{T}; rownames = prednames(aov)) where {T <: FullModel}
    throw(function_arg_error(anovatable, AnovaResult{T}))
end

function anovatable(::AnovaResult{T, S, N}; rownames = string.(1:N)) where {T <: NestedModels, S, N}
    throw(function_arg_error(anovatable, AnovaResult{T}))
end

# default anovatable api for comparing multiple models
function anovatable(aov::AnovaResult{<: NestedModels, FTest, N}; rownames = string.(1:N)) where N
    AnovaTable([
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    [NaN, _diffn(deviance(aov))...], 
                    teststat(aov), 
                    pval(aov)
                ],
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "ΔDeviance", "F value", "Pr(>|F|)"],
              rownames, 7, 6)
end 

function anovatable(aov::AnovaResult{<: NestedModels, LRT, N}; rownames = string.(1:N)) where N
    AnovaTable([
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    teststat(aov), 
                    pval(aov)
                ],
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "χ²", "Pr(>|χ²|)"],
              rownames, 6, 5)
end 

# Show function that delegates to anovatable
function show(io::IO, aov::AnovaResult{<: FullModel, T}) where {T <: GoodnessOfFit}
    at = anovatable(aov)
    println(io, "Analysis of Variance")
    println(io)
    println(io, "Type $(anova_type(aov)) test / $(testname(T))")
    println(io)
    println(io, formula(aov.anovamodel.model))
    println(io)
    println(io, "Table:")
    show(io, at)
end

function show(io::IO, aov::AnovaResult{<: NestedModels, T}) where {T <: GoodnessOfFit}
    at = anovatable(aov)
    println(io,"Analysis of Variance")
    println(io)
    println(io, "Type $(anova_type(aov)) test / $(testname(T))")
    println(io)
    for(id, m) in enumerate(aov.anovamodel.model)
        println(io,"Model $id: ", formula(m))
    end
    println(io)
    println(io,"Table:")
    show(io, at)
end