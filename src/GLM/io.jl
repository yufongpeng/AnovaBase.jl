# ======================================================================================================
# IO
<<<<<<< Updated upstream

"""
    coefnames(<model>, anova::Val{:anova})

Customize coefnames for anova
"""
coefnames(model::TableRegressionModel{<: GeneralizedLinearModel, T}, anova::Val{:anova}) where T =  begin
    v = coefnames(model.mf, anova)
    # push!(v,"(Dispersion)")
    v
end

coefnames(model::TableRegressionModel{<: LinearModel, T}, anova::Val{:anova}) where T =  begin
    coefnames(model.mf, anova)
end

coefnames(mf::ModelFrame, anova::Val{:anova}) = begin
    vectorize(coefnames(mf.f.rhs, anova))
end
=======
function coefnames(aov::AnovaResult{T, FTest}) where {T <: TableRegressionModel{<: Union{LinearModel, GeneralizedLinearModel}}}
    v = coefnames(aov.model, Val(:anova))
    push!(v, "(Residuals)")
    v
end

coefnames(trm::TableRegressionModel{<: Union{LinearModel, GeneralizedLinearModel}}, anova::Val{:anova}) =
    vectorize(coefnames(trm.mf.f.rhs, anova))
>>>>>>> Stashed changes

# anovatable api
function anovatable(model::AnovaResult{T, S}; kwargs...) where {T <: TableRegressionModel{<: LinearModel, <: AbstractArray}, S <: AbstractAnovaStats}
    at = anovatable(model.stats, kwargs...)
    cfnames = coefnames(model.model, Val(:anova))
    S <: AnovaStatsF && push!(cfnames,"(Residual)")
    if length(at.rownms) == length(cfnames)
        at.rownms = cfnames
    end
    at
end

# anovatable api for AnovaStats
function anovatable(stats::FixedAnovaStatsF{LinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...], [(stats.deviance) ./ stats.dof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "Sum of Squares", "Mean of Squares", "F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    at
end 

function anovatable(stats::FixedAnovaStatsF{GeneralizedLinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...], [(stats.deviance) ./ stats.dof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "ΔDeviance", "Mean of deviance", "F value","Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 5, 4)
    at
end 

function anovatable(stats::FixedAnovaStatsLRT{LinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...],  [stats.lrstat...], [stats.pval...]),
              ["DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end 

function anovatable(stats::FixedAnovaStatsLRT{GeneralizedLinearModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...],  [stats.lrstat...], [stats.pval...]),
              ["DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end

function anovatable(stats::NestedAnovaStatsF, rs::NTuple{N, Float64}, Δrs::NTuple{M, Float64}; kwargs...) where {N, M}
    at = AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs + 1 .- [stats.dof...], [rs...], [NaN, Δrs...], [stats.deviance...], [NaN, _diffn(stats.deviance)...], [stats.fstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "R²", "ΔR²", "Deviance", "Sum of Squares", "F value", "Pr(>|F|)"],
              ["$i" for i = 1:length(stats.dof)], 9, 8)
    at
end 

function anovatable(stats::NestedAnovaStatsLRT, rs::NTuple{N, Float64}, Δrs::NTuple{M, Float64}; kwargs...) where {N, M}
    # Δdeviance
    at = AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs + 1 .- [stats.dof...], [rs...], [NaN, Δrs...], [stats.deviance...], [stats.lrstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "R²", "ΔR²","Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["$i" for i = 1:length(stats.dof)], 8, 7)
    at
end 