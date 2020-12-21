# Old version functions
@deprecate function isbtw(fet::MatrixTerm, assign::Array{Int64,1}, remat::ReMat, X::Matrix)
    n = length(fet.terms)
    between = ones(Bool, n)
    select = 1:length(assign)
    for id in 1:size(remat, 2)
        loc = findall(==(1), view(remat, :, id))
        x = view(X, loc, :)
        for level in select
            if length(unique(x[:, level])) > 1
                factor = assign[level]
                select = select[assign[select] .!= factor]
                between[factor] = false
            end
        end
    end
    between[1] = false
    between
end

# Calculate number of groups
@deprecate nlevels(term::CategoricalTerm) = length(term.contrasts.levels)
@deprecate nlevels(term::ContinuousTerm) = 1 
@deprecate nlevels(term::InterceptTerm) = 1 
@deprecate nlevels(term::InteractionTerm) = prod(nlevels.(term.terms))

@deprecate _diff(v::Vector{T}) where T = cat(v[1], -diff(v), dims = 1)