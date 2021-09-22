# ================================================================================
# IO
coefnames(model::MixedModel, anova::Val{:anova}) = vectorize(coefnames(model.formula.rhs[1], anova))

# anovatable api
function _anovatable(aov::AnovaResult{<: LinearMixedModel, FTest}; kwargs...)
    AnovaTable(hcat(collect.((
                    dof(aov), 
                    dof_residual(aov), 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "Res.DOF", "F value", "Pr(>|F|)"],
              ["x$i" for i in eachindex(pval(aov))], 4, 3)
end
#=
function anovatable(stats::MixedAnovaStatsLRT{LinearMixedModel, N}; kwargs...) where N
    AnovaTable(hcat([stats.dof...], [stats.deviance...], [stats.lrstat...], [stats.pval...]),
              ["DOF", "Deviance", "Likelihood Ratio", "Pr(>|χ²|)"],
              ["x$i" for i = 1:length(stats.dof)], 4, 3)
end 

function anovatable(stats::NestedAnovaStatsF; modeltype::Type{<: LinearMixedModel}, model::Tuple, kwargs...) where N
    AnovaTable(hcat([stats.dof...], [NaN, _diff(stats.dof)...], stats.nobs .+ 1 .- [stats.dof...], [aic.(model)...], [bic.(model)...], [stats.deviance...], [NaN, _diffn(stats.deviance)...], [stats.fstat...], [stats.pval...]),
              ["DOF", "ΔDOF", "Res. DOF", "AIC", "BIC", "Deviance", "ΔDeviance", "F value", "Pr(>|F|)"],
              ["$i" for i = 1:length(stats.dof)], 9, 8)
end
=#
function _anovatable(aov::AnovaResult{<: Tuple, LRT}, 
                    modeltype1::Union{Type{<: LinearMixedModel}, Type{<: TableRegressionModel{<: GLM.LinearModel}}},
                    modeltype2::Type{<: LinearMixedModel}; 
                    kwargs...)
    if last(aov.model).optsum.REML 
        AnovaTable(hcat(collect.((
                        dof(aov), 
                        [NaN, _diff(dof(aov))...], 
                        dof_residual(aov), 
                        deviance(aov), 
                        teststat(aov), 
                        pval(aov)
                        ))...),
                ["DOF", "ΔDOF", "Res.DOF", "-2 logLik", "χ²", "Pr(>|χ²|)"],
                ["$i" for i in eachindex(pval(aov))], 6, 5)
    else
        AnovaTable(hcat(collect.((
                        dof(aov), 
                        [NaN, _diff(dof(aov))...], 
                        dof_residual(aov), 
                        aic.(aov.model),
                        bic.(aov.model),
                        deviance(aov), 
                        teststat(aov), 
                        pval(aov)
                        ))...),
                ["DOF", "ΔDOF", "Res.DOF", "AIC", "BIC", "-2 logLik", "χ²", "Pr(>|χ²|)"],
                ["$i" for i in eachindex(pval(aov))], 8, 7)
    end
end