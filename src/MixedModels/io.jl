# ================================================================================
# IO

"""
    coefnames(<model>, anova::Val{:anova})

Customize coefnames for anova
"""
coefnames(model::MixedModel, anova::Val{:anova}) = begin 
    v = vectorize(coefnames(model.formula.rhs[1], anova))
    # push!(v, "(Residual)", "(Residual)")
    v
end

# anovatable api for AnovaStats
function anovatable(stats::MixedAnovaStatsF; kwargs...)
    at = AnovaTable(hcat([stats.dof...], [stats.resdof...], [stats.fstat...], [stats.pval...]),
              ["DOF", "Res.DOF", "F value", "Pr(>|F|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end

function anovatable(stats::MixedAnovaStatsLRT{LinearMixedModel, N}; kwargs...) where N
    at = AnovaTable(hcat([stats.dof...], [stats.deviance...],  [stats.lrstat...], [stats.pval...]),
              ["DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
    at
end 