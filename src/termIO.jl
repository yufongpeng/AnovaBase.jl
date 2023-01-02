# =========================================================================
# Function related to terms, variable names and I/O

"""
    factornames(<term>)

Return string repressentation of a factor.
"""
factornames(t::MatrixTerm) = mapreduce(factornames, vcat, t.terms)
factornames(t::FormulaTerm) = (factornames(t.lhs), vectorize(factornames(t.rhs)))
factornames(::InterceptTerm{H}) where H = H ? ["(Intercept)"] : String[]
factornames(t::ContinuousTerm) = string(t.sym)
factornames(t::CategoricalTerm) = string(t.sym)
factornames(t::FunctionTerm) = string(t.exorig)
factornames(ts::StatsModels.TupleTerm) = mapreduce(factornames, vcat, ts)
factornames(t::InteractionTerm) = join(factornames.(t.terms), " & ")
factornames(t::Term) = string(t)
factornames(t::ConstantTerm{H}) where H = string(t)
factornames(t) = coefnames(t)

# Base.show(io::IO, t::FunctionTerm) = print(io, "$(t.exorig)")

"""
    getterms(<term>)

Get an array of `Expr` or `Symbol` of terms.
"""
getterms(::AbstractTerm) = Symbol[]
getterms(term::Union{Term, CategoricalTerm, ContinuousTerm}) = [term.sym]
#getterms(::InterceptTerm{true}) = [Symbol(1)]
getterms(term::InteractionTerm) = mapreduce(getterms, union, term.terms)
getterms(term::FunctionTerm) = [term.exorig]
getterms(term::StatsModels.TupleTerm) = mapreduce(getterms, union, term)
getterms(term::MatrixTerm) = getterms(term.terms)
#=
"""
    getterm(term::AbstractTerm)
    getterm(term::FunctionTerm)

Get the `Symbol` of a single term.
"""=#
# getterm(term::AbstractTerm) = term.sym
# getterm(term::FunctionTerm) = term.exorig
# getterm(term::InterceptTerm) = Symbol(1)

# Determine selected terms for type 2 ss
"""
    isinteract(f::MatrixTerm, id1::Int, id2::Int)

Determine if `f.terms[id2]` is an interaction term of `f.terms[id1]` and other terms.
"""
isinteract(f::MatrixTerm, id1::Int, id2::Int) = issubset(getterms(f.terms[id1]), getterms(f.terms[id2]))
  
"""
    select_super_interaction(f::MatrixTerm, id::Int)

Return a set of index of `f.terms` which are interaction terms of `f.terms[id]` and other terms.
"""
select_super_interaction(f::MatrixTerm, id::Int) = 
    id == 1 ? Set(eachindex(f.terms)) : Set([idn for idn in eachindex(f.terms) if isinteract(f, id, idn)])

"""
    select_sub_interaction(f::MatrixTerm, id::Int)

Return a set of index of `f.terms` which `f.terms[id]` are interaction terms of those terms and other terms.
"""
select_sub_interaction(f::MatrixTerm, id::Int) = 
    id == 1 ? Set(Int[]) : Set([idn for idn in eachindex(f.terms) if isinteract(f, idn, id)])

"""
    select_not_super_interaction(f::MatrixTerm, id::Int)

Return a set of index of `f.terms` which are not interaction terms of `f.terms[id]` and other terms.
"""
select_not_super_interaction(f::MatrixTerm, id::Int) = 
    id == 1 ? Set(Int[]) : Set([idn for idn in eachindex(f.terms) if !isinteract(f, id, idn)])

"""
    select_not_interaction(f::MatrixTerm, id::Int)

Return a set of index of `f.terms` which `f.terms[id]` are not interaction terms of those terms and other terms.
"""
select_not_sub_interaction(f::MatrixTerm, id::Int) = 
    id == 1 ? Set(eachindex(f.terms)) : Set([idn for idn in eachindex(f.terms) if !isinteract(f, idn, id)])

@deprecate selectcoef(f::MatrixTerm, id::Int) select_super_interaction(f, id)

# Create sub-formula
"""
    subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int; reschema::Bool = false)
    subformula(lhs::AbstractTerm, rhs::MatrixTerm, id; reschema::Bool = false)
    subformula(lhs::AbstractTerm, rhs::NTuple{N, AbstractTerm}, id::Int)

Create formula from existing `lhs` and `rhs` truncated to `1:id` or excluded collection `id`.
"""
subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int; reschema::Bool = false) = 
    id > 0 ? reschema_formula(lhs, rhs.terms[1:id], reschema) : reschema_formula(lhs, (InterceptTerm{false}(),), reschema)
# For nestedmodels & GLM

function subformula(lhs::AbstractTerm, rhs::MatrixTerm, id; reschema::Bool = false)
    terms = rhs.terms[setdiff(eachindex(rhs.terms), id)]
    1 in id && (terms = (InterceptTerm{false}(), terms...))
    reschema_formula(lhs, terms, reschema)
end
# For type 2/3 & GLM

subformula(lhs::AbstractTerm, rhs::NTuple{N, AbstractTerm}, id::Int) where N = 
    id > 0 ? FormulaTerm(lhs, (collect_matrix_terms(first(rhs).terms[1:id]), rhs[2:end]...)) : FormulaTerm(lhs, (collect_matrix_terms((InterceptTerm{false}(),)), rhs[2:end]...))
# For nestedmodels & MixedModels

"""
    extract_contrasts(f::FormulaTerm)

Extract a dictionary of contrasts.
"""
extract_contrasts(f::FormulaTerm) = 
    Dict{Symbol, Any}(t.sym => t.contrasts.contrasts for t in f.rhs.terms if isa(t, CategoricalTerm))

"""
    clearschema(<terms with schema>) = <terms without schema>

Clear any applied schema on terms.
"""
clearschema(::InterceptTerm{true}) = ConstantTerm(1)
clearschema(::InterceptTerm{false}) = ConstantTerm(0)
clearschema(t::FunctionTerm) = t
clearschema(t::Union{CategoricalTerm, ContinuousTerm}) = Term(t.sym)
clearschema(t::InteractionTerm) = InteractionTerm(clearschema.(t.terms))

# reschema only happen when using TupleTerm rather than MatrixTerm
reschema_formula(lhs::AbstractTerm, ts::StatsModels.TupleTerm, reschema::Bool) = 
    reschema ? FormulaTerm(clearschema(lhs), clearschema.(ts)) : FormulaTerm(lhs, collect_matrix_terms(ts))


# test name
tname(::Type{FTest}) = "F test"
tname(::Type{LRT}) = "Likelihood-ratio test"
#tname(M::AnovaStatsRao) = "Rao score test"
#tname(M::AnovaStatsCp) = "Mallow's Cp"

# ============================================================================================================================
# AnovaTable mostly from CoefTable
mutable struct AnovaTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    pvalcol::Int
    teststatcol::Int
    function AnovaTable(cols::Vector, colnms::Vector, rownms::Vector,
                       pvalcol::Int=0, teststatcol::Int = 0)
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

    function AnovaTable(mat::Matrix, colnms::Vector, rownms::Vector,
                       pvalcol::Int=0, teststatcol::Int=0)
        nc = size(mat, 2)
        cols = Any[mat[:, i] for i in 1:nc]
        AnovaTable(cols, colnms, rownms, pvalcol, teststatcol)
    end
end

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
function anova_table(aov::AnovaResult{<: RegressionModel}; kwargs...)
    at = anovatable(aov; kwargs...)
    cfnames = factornames(aov)
    # when first factor is null
    length(at.rownms) == length(cfnames) - 1 && popfirst!(cfnames)
    at.rownms = cfnames
    at
end

# check first and last modeltype
anova_table(aov::AnovaResult{<: Tuple}; kwargs...) = 
    anovatable(aov, typeof(first(aov.model)), typeof(last(aov.model)); kwargs...)

"""
    anovatable(aov::AnovaResult{<: RegressionModel}; kwargs...)
    anovatable(aov::AnovaResult{<: Tuple}, modeltype1, modeltype2; kwargs...)

Return a table with coefficients and related statistics of ANOVA. For nested models, the function will dispatch on the types of the first and the last models. For a single model, no default api was implemented.

The returned AnovaTable object implements the Tables.jl (https://github.com/JuliaData/Tables.jl/) interface, and can be  
converted e.g. to a DataFrame via using DataFrames; DataFrame(anovatable(aov)).
"""
anovatable(aov::AnovaResult{<: Tuple}, modeltype1, modeltype2; kwargs...) = anovatable(aov; kwargs...)

# default anovatable api for comparing multiple models
function anovatable(aov::AnovaResult{<: Tuple, FTest}; kwargs...)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    [NaN, _diffn(deviance(aov))...], 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "ΔDeviance", "F value", "Pr(>|F|)"],
              ["$i" for i in eachindex(pval(aov))], 7, 6)
end 

function anovatable(aov::AnovaResult{<: Tuple, LRT}; kwargs...)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    deviance(aov), 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "ΔDOF", "Res.DOF", "Deviance", "χ²", "Pr(>|χ²|)"],
              ["$i" for i in eachindex(pval(aov))], 6, 5)
end 

# Show function that delegates to anovatable
function show(io::IO, aov::AnovaResult{<: RegressionModel, T}) where {T <: GoodnessOfFit}
    at = anova_table(aov)
    println(io, "Analysis of Variance")
    println(io)
    println(io, "Type $(aov.type) test / $(tname(T))")
    println(io)
    println(io, formula(aov.model))
    println(io)
    println(io, "Table:")
    show(io, at)
end

function show(io::IO, aov::AnovaResult{<: Tuple, T}) where {T <: GoodnessOfFit}
    at = anova_table(aov)
    println(io,"Analysis of Variance")
    println(io)
    println(io, "Type $(aov.type) test / $(tname(T))")
    println(io)
    for(id, m) in enumerate(aov.model)
        println(io,"Model $id: ", formula(m))
    end
    println(io)
    println(io,"Table:")
    show(io, at)
end

