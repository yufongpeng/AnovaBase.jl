# =========================================================================
# Function related to terms, variable names and I/O

# Customize coefnames
StatsBase.coefnames(model::TableRegressionModel{<: GeneralizedLinearModel, T},anova::Val{:anova}) where T =  begin
    v = coefnames(model.mf, anova)
    # push!(v,"(Dispersion)")
    v
end

StatsBase.coefnames(model::TableRegressionModel{<: LinearModel, T},anova::Val{:anova}) where T =  begin
    v = coefnames(model.mf, anova)
    push!(v,"(Residual)")
end

StatsBase.coefnames(mf::ModelFrame, anova::Val{:anova}) = begin
    vectorize(coefnames(mf.f.rhs, anova))
end

StatsBase.coefnames(model::MixedModel, anova::Val{:anova}) = begin 
    v = vectorize(coefnames(model.formula.rhs[1], anova))
    # push!(v, "(Residual)", "(Residual)")
    v
end

StatsBase.coefnames(t::MatrixTerm, anova::Val{:anova}) = mapreduce(coefnames, vcat, t.terms, repeat([anova], length(t.terms)))
StatsBase.coefnames(t::FormulaTerm, ::Val{:anova}) = (coefnames(t.lhs), coefnames(t.rhs))
StatsBase.coefnames(::InterceptTerm{H}, ::Val{:anova}) where {H} = H ? "(Intercept)" : []
StatsBase.coefnames(t::ContinuousTerm, ::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::CategoricalTerm, ::Val{:anova}) = string(t.sym)
StatsBase.coefnames(t::FunctionTerm, ::Val{:anova}) = string(t.exorig)
StatsBase.coefnames(ts::StatsModels.TupleTerm, ::Val{:anova}) = reduce(vcat, coefnames.(ts))

StatsBase.coefnames(t::InteractionTerm, anova::Val{:anova}) = begin
    join(coefnames.(t.terms, anova), " & ")
end
    
# Base.show(io::IO, t::FunctionTerm) = print(io, "$(t.exorig)")

# =============================================================================================================================
# Subsetting coefnames for type 2 anova
getterms(term::AbstractTerm) = Union{Symbol,Expr}[term.sym]
getterms(::InterceptTerm) = Union{Symbol,Expr}[Symbol(1)]
getterms(term::InteractionTerm) = map(i->getterm(i), term.terms)
getterms(term::FunctionTerm) = Union{Symbol,Expr}[term.exorig]
getterms(term::MatrixTerm) = union(getterms.(term.terms)...)
getterm(term::AbstractTerm) = term.sym
getterm(term::FunctionTerm) = term.exorig
getterm(term::InterceptTerm) = Symbol(1)

isinteract(f::MatrixTerm, id1::Int, id2::Int) = issubset(getterms(f.terms[id1]), getterms(f.terms[id2]))
selectcoef(f::MatrixTerm, id::Int) = Set([comp for comp in 1:length(f.terms) if isinteract(f, id, comp)])

# Unify formula api
formula(model::TableRegressionModel) = model.mf.f
formula(model::MixedModel) = model.formula

# Calculate dof from assign
function dof(v::Vector{Int})
    dofv = zeros(Int, v[end])
    prev = 1
    ind = 1
    n = length(v)
    while ind <= n
        v[ind] == prev || (prev = v[ind])
        dofv[prev] += 1
        ind += 1
    end
    dofv
end

# test name
tname(M::AnovaStatsF) = "F test"
tname(M::AnovaStatsLRT) = "Likelihood-ratio test"
tname(M::AnovaStatsRao) = "Rao score test"
tname(M::AnovaStatsCp) = "Mallow's Cp"

# ============================================================================================================================
# AnovaTable mostly from CoefTable
mutable struct AnovaTable
    cols::Vector
    colnms::Vector
    rownms::Vector
    pvalcol::Int
    teststatcol::Int
    function AnovaTable(cols::Vector,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = length(cols)
        nrs = map(length,cols)
        nr = nrs[1]
        length(colnms) in [0,nc] || throw(ArgumentError("colnms should have length 0 or $nc"))
        length(rownms) in [0,nr] || throw(ArgumentError("rownms should have length 0 or $nr"))
        all(nrs .== nr) || throw(ArgumentError("Elements of cols should have equal lengths, but got $nrs"))
        pvalcol in 0:nc || throw(ArgumentError("pvalcol should be between 0 and $nc"))
        teststatcol in 0:nc || throw(ArgumentError("teststatcol should be between 0 and $nc"))
        new(cols,colnms,rownms,pvalcol,teststatcol)
    end

    function AnovaTable(mat::Matrix,colnms::Vector,rownms::Vector,
                       pvalcol::Int=0,teststatcol::Int=0)
        nc = size(mat,2)
        cols = Any[mat[:, i] for i in 1:nc]
        AnovaTable(cols,colnms,rownms,pvalcol,teststatcol)
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
        rownms = [lpad("[$i]",floor(Integer, log10(nr))+3) for i in 1:nr]
    end
    mat = [j == 1 ? NoQuote(rownms[i]) :
           j-1 == at.pvalcol ? PValue(cols[j-1][i]) :
           j-1 in at.teststatcol ? TestStat(cols[j-1][i]) :
           cols[j-1][i] isa AbstractString ? NoQuote(cols[j-1][i]) : OtherStat(cols[j-1][i])
           for i in 1:nr, j in 1:nc+1]
    # Code inspired by print_matrix in Base
    io = IOContext(io, :compact=>true, :limit=>false)
    A = Base.alignment(io, mat, 1:size(mat, 1), 1:size(mat, 2),
                       typemax(Int), typemax(Int), 3)
    nmswidths = pushfirst!(length.(colnms), 0)
    A = [nmswidths[i] > sum(A[i]) ? (A[i][1]+nmswidths[i]-sum(A[i]), A[i][2]) : A[i]
         for i in 1:length(A)]
    totwidth = sum(sum.(A)) + 2 * (length(A) - 1)
    println(io, repeat('─', totwidth))
    print(io, repeat(' ', sum(A[1])))
    for j in 1:length(colnms)
        print(io, "  ", lpad(colnms[j], sum(A[j+1])))
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
# Anovatable api
function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: RegressionModel, S <: AbstractAnovaStats}
    at = anovatable(model.stats, kwargs...)
    cfnames = coefnames(model.model, Val(:anova))
    if length(at.rownms) == length(cfnames)
        at.rownms = cfnames
    end
    at
end

function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: TableRegressionModel{<: LinearModel, U}, S <: AbstractAnovaStats} where U
    at = anovatable(model.stats, kwargs...)
    cfnames = coefnames(model.model, Val(:anova))
    if length(at.rownms) == length(cfnames)
        at.rownms = cfnames
    end
    at
end


function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: Tuple, S <: AbstractAnovaStats} where U
    at = anovatable(model.stats, kwargs...)
    cfnames = coefnames(last(model.model), Val(:anova))
    if length(at.rownms) == length(cfnames)
        at.rownms = cfnames
    end
    at
end

function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: Tuple, S <: NestedAnovaStats}
    anovatable(model.stats, kwargs...)
end

function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: Tuple, S <: NestedAnovaStatsF}
    try
        rs = r2.(model.model)
        Δrs = _diff(rs)
        anovatable(model.stats, rs, Δrs, kwargs...)
    catch 
        anovatable(model.stats, kwargs...)
    end
end

# ----------------------------------------------------------------------------------------------------------------------------------------
# anovatable for AnovaStats
function anovatable(stats::FixedAnovaStatsF{LinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...], [(stats.deviance) ./ stats.dof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "Sum of Squares", "Mean of Squares", "F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    at
end 

function anovatable(stats::FixedAnovaStatsF{GeneralizedLinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...], [(stats.deviance) ./ stats.dof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "Deviance", "Mean of deviance", "F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    at
end 

function anovatable(stats::FixedAnovaStatsLRT{GeneralizedLinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...],  [stats.lrstat...], [stats.pval...]),
              ["DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end 

function anovatable(stats::MixedAnovaStatsF; kwargs...)
    at = AnovaTable(hcat([stats.dof...], [stats.resdof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "Res.DOF", "F value", "Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end

function anovatable(stats::NestedAnovaStatsF; kwargs...)
    at = AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs .+ 1 .- [stats.dof...], [stats.deviance...], [NaN, _diffn(stats.deviance)...], [stats.fstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "Deviance", "ΔDeviance", "F value", "Pr(>|F|)"],
              ["$i" for i = 1:length(stats.dof)], 7, 6)
    at
end 

function anovatable(stats::NestedAnovaStatsF, rs::NTuple{N, Float64}, Δrs::NTuple{M, Float64}; kwargs...) where {N, M}
    at = AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs + 1 .- [stats.dof...], [rs...], [NaN, Δrs...], [stats.deviance...], [NaN, _diffn(stats.deviance)...], [stats.fstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "R²", "ΔR²", "Deviance", "Sum of Squares", "F value", "Pr(>|F|)"],
              ["$i" for i = 1:length(stats.dof)], 9, 8)
    at
end 

function anovatable(stats::NestedAnovaStatsLRT; kwargs...)
    # Δdeviance
    at = AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs + 1 .- [stats.dof...], [stats.deviance...], [stats.lrstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["$i" for i = 1:length(stats.dof)], 6, 5)
    at
end 

# Show function that delegates to anovatable
function show(io::IO, model::AnovaResult{T, S}) where {T <: RegressionModel, S <: AbstractAnovaStats}
    at = anovatable(model)
    println(io, "Analysis of Variance")
    println(io)
    println(io, "Type $(model.stats.type) test / F test")
    println(io)
    println(io, formula(model.model))
    println(io)
    println(io, "Table:")
    show(io, at)
end

function show(io::IO, model::AnovaResult{T, S}) where {T <: Tuple, S <: AbstractAnovaStats}
    at = anovatable(model)
    println(io, "Analysis of Variance")
    println(io)
    println(io, "Type $(model.stats.type) test / F test")
    println(io)
    println(io, formula(last(model.model)))
    println(io)
    println(io, "Table:")
    show(io, at)
end

function show(io::IO, model::AnovaResult{T, S}) where {T <: Tuple, S <: NestedAnovaStats}
    at = anovatable(model)
    println(io,"Analysis of Variance")
    println(io)
    println(io, "$(tname(model.stats))")
    println(io)
    for(id, m) in enumerate(model.model)
        println(io,"Model $id: ", formula(m))
    end
    println(io)
    println(io,"Table:")
    show(io, at)
end

