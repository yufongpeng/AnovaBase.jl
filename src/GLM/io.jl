# ======================================================================================================
# IO
function coefnames(aov::AnovaResult{T, FTest}; kwargs...) where {T <: TableRegressionModel{<: LinPredModel}}
    v = coefnames(aov.model, Val(:anova))
    push!(v, "(Residuals)")
    v
end

coefnames(trm::TableRegressionModel{<: LinPredModel}, anova::Val{:anova}) =
    coefnames(trm.mf, anova)

coefnames(mf::ModelFrame, anova::Val{:anova}) = vectorize(coefnames(mf.f.rhs, anova))

# anovatable api
function _anovatable(aov::AnovaResult{<: TableRegressionModel{<: LinearModel}, FTest}; kwargs...)
    AnovaTable(hcat(vectorize.((dof(aov), deviance(aov), deviance(aov) ./ dof(aov), teststat(aov), pval(aov)))...),
              ["DOF", "Exp.SS", "Mean Square", "F value","Pr(>|F|)"],
              ["x$i" for i in eachindex(pval(aov))], 5, 4)
end 

function _anovatable(aov::AnovaResult{<: TableRegressionModel{<: GeneralizedLinearModel}, FTest}; kwargs...)
    AnovaTable(hcat(vectorize.((dof(aov), deviance(aov), deviance(aov) ./ dof(aov), teststat(aov), pval(aov)))...),
              ["DOF", "ΔDeviance", "Mean ΔDev", "F value","Pr(>|F|)"],
              ["x$i" for i in eachindex(pval(aov))], 5, 4)
end 

function _anovatable(aov::AnovaResult{<: TableRegressionModel{<: LinearModel}, LRT}; kwargs...)
    AnovaTable(hcat(vectorize.((dof(aov), deviance(aov), teststat(aov), pval(aov)))...),
              ["DOF", "Res.SS", "χ²", "Pr(>|χ²|)"],
              ["x$i" for i in eachindex(pval(aov))], 4, 3)
end 

function _anovatable(aov::AnovaResult{<: TableRegressionModel{<: GeneralizedLinearModel}, LRT}; kwargs...)
    AnovaTable(hcat(vectorize.((dof(aov), deviance(aov), teststat(aov), pval(aov)))...),
              ["DOF", "Deviance", "χ²", "Pr(>|χ²|)"],
              ["x$i" for i in eachindex(pval(aov))], 4, 3)
end

function _anovatable(aov::AnovaResult{<: Tuple, FTest}, 
                    modeltype1::Type{<: TableRegressionModel{<: LinearModel}}, 
                    modeltype2::Type{<: TableRegressionModel{<: LinearModel}};
                    kwargs...)

    rs = r2.(aov.model)
    Δrs = _diff(rs)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    rs,
                    [NaN, Δrs...],
                    deviance(aov), 
                    [NaN, _diffn(deviance(aov))...], 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "ΔDOF", "Res.DOF", "R²", "ΔR²", "Res.SS", "Exp.SS", "F value", "Pr(>|F|)"],
              ["$i" for i in eachindex(pval(aov))], 9, 8)
end 

function _anovatable(aov::AnovaResult{<: Tuple, LRT}, 
                    modeltype1::Type{<: TableRegressionModel{<: LinearModel}},
                    modeltype2::Type{<: TableRegressionModel{<: LinearModel}};
                    kwargs...)
    rs = r2.(aov.model)
    Δrs = _diff(rs)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    rs,
                    [NaN, Δrs...],
                    deviance(aov), 
                    teststat(aov), 
                    pval(aov)
                    ))...),
              ["DOF", "ΔDOF", "Res.DOF", "R²", "ΔR²", "Res.SS", "χ²", "Pr(>|χ²|)"],
              ["$i" for i in eachindex(pval(aov))], 8, 7)
end 