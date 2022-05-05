# ======================================================================================================
# IO
function coefnames(aov::AnovaResult{T, FTest}; kwargs...) where {T <: TableRegressionModel{<: FixedEffectModel}}
    v = coefnames(aov.model, Val(:anova))
    push!(v, "(Residuals)")
    v
end

coefnames(trm::TableRegressionModel{<: FixedEffectModel}, anova::Val{:anova}) =
    vectorize(coefnames(formula(trm).rhs.terms[unique(trm.mm.assign)], anova))

# anovatable api
function _anovatable(aov::AnovaResult{<: TableRegressionModel{<: FixedEffectModel}, FTest}; kwargs...)
    AnovaTable(hcat(vectorize.((dof(aov), deviance(aov), deviance(aov) ./ dof(aov), teststat(aov), pval(aov)))...),
              ["DOF", "Exp.SS", "Mean Square", "F value","Pr(>|F|)"],
              ["x$i" for i in eachindex(pval(aov))], 5, 4)
end 

function _anovatable(aov::AnovaResult{<: Tuple, FTest}, 
                    modeltype1::Type{<: TableRegressionModel{<: FixedEffectModel}}, 
                    modeltype2::Type{<: TableRegressionModel{<: FixedEffectModel}};
                    kwargs...)

    rs = r2.(aov.model)
    rws = ntuple(length(aov.model)) do i 
        aov.model[i].model.r2_within
    end
    Δrs = _diff(rs)
    Δrws = _diff(rws)
    AnovaTable(hcat(vectorize.((
                    dof(aov), 
                    [NaN, _diff(dof(aov))...], 
                    dof_residual(aov), 
                    rs,
                    [NaN, Δrs...],
                    rws,
                    [NaN, Δrws...],
                    deviance(aov), 
                    [NaN, _diffn(deviance(aov))...], 
                    teststat(aov), 
                    pval(aov)
                    ))...),
                ["DOF", "ΔDOF", "Res.DOF", "R²", "ΔR²", "R²_within", "ΔR²_within", "Res.SS", "Exp.SS", "F value", "Pr(>|F|)"],
                ["$i" for i in eachindex(pval(aov))], 11, 10)
end 