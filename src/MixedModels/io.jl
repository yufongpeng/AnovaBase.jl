# ================================================================================
# IO
coefnames(model::MixedModel, anova::Val{:anova}) = 
    vectorize(coefnames(first(model.formula.rhs), anova))
  
# anovatable api
function _anovatable(aov::AnovaResult{<: LinearMixedModel, FTest}; kwargs...)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    dof_residual(aov), 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "Res.DOF", "F value", "Pr(>|F|)"],
              ["x$i" for i in eachindex(pval(aov))], 4, 3)
end

function _anovatable(aov::AnovaResult{<: Tuple, LRT}, 
                    modeltype1::Union{Type{<: LinearMixedModel}, Type{<: TableRegressionModel{<: GLM.LinearModel}}},
                    modeltype2::Type{<: LinearMixedModel}; 
                    kwargs...)
    if last(aov.model).optsum.REML 
        AnovaTable(hcat(vectorize.((
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
        AnovaTable(hcat(vectorize.((
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