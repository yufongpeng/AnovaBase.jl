# Function related to terms
"""
    any_not_aliased_with_1(terms)

Return `true` if there are any terms not aliased with the intercept, e.g. `ContinuousTerm` or `FunctionTerm`.

Terms without schema are considered aliased with the intercept.
"""
any_not_aliased_with_1(f::FormulaTerm) = any_not_aliased_with_1(f.rhs)
any_not_aliased_with_1(f::MatrixTerm) = any_not_aliased_with_1(f.terms)
any_not_aliased_with_1(f::TupleTerm) = any(any_not_aliased_with_1, f)
any_not_aliased_with_1(::InterceptTerm) = false
any_not_aliased_with_1(t::ContinuousTerm) = true
any_not_aliased_with_1(t::CategoricalTerm) = false
any_not_aliased_with_1(t::FunctionTerm) = true
any_not_aliased_with_1(t::InteractionTerm) = any(any_not_aliased_with_1, t.terms)
any_not_aliased_with_1(t::Term) = false
any_not_aliased_with_1(t::ConstantTerm) = false

"""
    prednames(term)

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
prednames(ts::TupleTerm) = filter!(!isnothing, vectorize(mapreduce(prednames, vcat, ts)))
prednames(::InterceptTerm{H}) where H = H ? "(Intercept)" : nothing
prednames(t::ContinuousTerm) = string(t.sym)
prednames(t::CategoricalTerm) = string(t.sym)
prednames(t::FunctionTerm) = string(t.exorig)
prednames(t::InteractionTerm) = join(prednames.(t.terms), " & ")
prednames(t::Term) = string(t)
prednames(t::ConstantTerm) = string(t)
prednames(t) = coefnames(t)

"""
    getterms(term)

Return the symbol of term(s) as a vector of `Expr` or `Symbol`.

# Examples
```julia
julia> iris = dataset("datasets", "iris");

julia> f = formula(lm(@formula(log(SepalLength) ~ Species + PetalLength * PetalWidth), iris))
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)
  PetalLength(continuous)
  PetalWidth(continuous)
  PetalLength(continuous) & PetalWidth(continuous)

julia> getterms(f)
(Expr[:(log(SepalLength))], [:Species, :PetalLength, :PetalWidth])

julia> getterms(InterceptTerm{true}())
Symbol[]

```
"""
getterms(::AbstractTerm) = Symbol[]
getterms(t::FormulaTerm) = (getterms(t.lhs), vectorize(getterms(t.rhs)))
getterms(term::MatrixTerm) = getterms(term.terms)
getterms(term::TupleTerm) = mapreduce(getterms, union, term)
getterms(term::Union{Term, CategoricalTerm, ContinuousTerm}) = [term.sym]
#getterms(::InterceptTerm{true}) = [Symbol(1)]
getterms(term::InteractionTerm) = mapreduce(getterms, union, term.terms)
getterms(term::FunctionTerm) = [term.exorig]

# Determine selected terms for type 2 ss
"""
    isinteract(f::Union{MatrixTerm, TupleTerm}, id1::Int, id2::Int)

Determine if `f[id2]` is an interaction term of `f[id1]` and other terms.

# Examples
```julia
julia> iris = dataset("datasets", "iris");

julia> f = formula(lm(@formula(log(SepalLength) ~ Species + PetalLength * PetalWidth), iris))
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)
  PetalLength(continuous)
  PetalWidth(continuous)
  PetalLength(continuous) & PetalWidth(continuous)

julia> isinteract(f.rhs, 1, 2)
true

julia> isinteract(f.rhs, 3, 4)
false

julia> isinteract(f.rhs, 4, 5)
true

```
"""
isinteract(f::MatrixTerm, id1::Int, id2::Int) = isinteract(f.terms, id1, id2)
isinteract(f::TupleTerm, id1::Int, id2::Int) = issubset(getterms(f[id1]), getterms(f[id2]))

const doc_select_interaction = """
    select_super_interaction(f::Union{MatrixTerm, TupleTerm}, id::Int)
    select_sub_interaction(f::Union{MatrixTerm, TupleTerm}, id::Int)
    select_not_super_interaction(f::Union{MatrixTerm, TupleTerm}, id::Int)
    select_not_sub_interaction(f::Union{MatrixTerm, TupleTerm}, id::Int)

Return a set of index of `f`, which

1. returned terms are interaction terms of `f[id]` and other terms.\n
2. `f[id]` is an interaction term of returned terms and other terms.\n
3. returned terms not interaction terms of `f[id]` and other terms.\n
4. `f[id]` is not interaction term of returned terms and other terms.

# Examples
```julia
julia> iris = dataset("datasets", "iris");

julia> f = formula(lm(@formula(log(SepalLength) ~ Species + PetalLength * PetalWidth), iris))
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)
  PetalLength(continuous)
  PetalWidth(continuous)
  PetalLength(continuous) & PetalWidth(continuous)

julia> select_super_interaction(f.rhs, 3)
Set{Int64} with 2 elements:
  5
  3

julia> select_sub_interaction(f.rhs, 3)
Set{Int64} with 2 elements:
  3
  1

julia> select_not_super_interaction(f.rhs, 3)
Set{Int64} with 3 elements:
  4
  2
  1

julia> select_not_sub_interaction(f.rhs, 3)
Set{Int64} with 3 elements:
  5
  4
  2

```
"""
@doc doc_select_interaction
select_super_interaction(f::MatrixTerm, id::Int) = select_super_interaction(f.terms, id) 
function select_super_interaction(f::TupleTerm, id::Int)
    s = id ≡ 1 ? Set(eachindex(f)) : Set([idn for idn in eachindex(f) if isinteract(f, id, idn)])
    hasintercept(f) || filter!(!=(1), s)
    s
end

@doc doc_select_interaction
select_sub_interaction(f::MatrixTerm, id::Int) = select_sub_interaction(f.terms, id)
function select_sub_interaction(f::TupleTerm, id::Int)
    s = id ≡ 1 ? Set(Int[]) : Set([idn for idn in eachindex(f) if isinteract(f, idn, id)])
    hasintercept(f) || filter!(!=(1), s)
    s
end

@doc doc_select_interaction
select_not_super_interaction(f::MatrixTerm, id::Int) = select_not_super_interaction(f.terms, id)
function select_not_super_interaction(f::TupleTerm, id::Int)
    s = id ≡ 1 ? Set(Int[]) : Set([idn for idn in eachindex(f) if !isinteract(f, id, idn)])
    hasintercept(f) || filter!(!=(1), s)
    s
end

@doc doc_select_interaction
select_not_sub_interaction(f::MatrixTerm, id::Int) = select_not_sub_interaction(f.terms, id)
function select_not_sub_interaction(f::TupleTerm, id::Int)
    s = id ≡ 1 ? Set(eachindex(f)) : Set([idn for idn in eachindex(f) if !isinteract(f, idn, id)])
    hasintercept(f) || filter!(!=(1), s)
    s
end

# Create sub-formula
"""
    subformula(f::FormulaTerm, id; kwargs...)
    subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int; reschema::Bool = false)
    subformula(lhs::AbstractTerm, rhs::MatrixTerm, id; reschema::Bool = false)
    subformula(lhs::AbstractTerm, rhs::NTuple{N, AbstractTerm}, id::Int; rhs_id::Int = 1, reschema::Bool = false)

Create formula from existing `lhs` and `rhs` (or `rhs[tuple_id]`) truncated to `1:id` or excluded collection `id`. 
When `id` is 0, all terms in `rhs` (or `rhs[tuple_id]`) will be removed.

If `reschema` is true, all terms' schema will be removed.

# Examples
```julia
julia> iris = dataset("datasets", "iris");

julia> f = formula(lm(@formula(log(SepalLength) ~ Species + PetalLength * PetalWidth), iris))
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)
  PetalLength(continuous)
  PetalWidth(continuous)
  PetalLength(continuous) & PetalWidth(continuous)

julia> subformula(f, 2)
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)

julia> subformula(f, [3, 5]; reschema = true)
FormulaTerm
Response:
  (SepalLength)->log(SepalLength)
Predictors:
  1
  Species(DummyCoding:3→2)
  PetalWidth(unknown)

julia> f = formula(fit(LinearMixedModel, @formula(SepalLength ~ SepalWidth + (SepalWidth|Species)), iris))
FormulaTerm
Response:
  SepalLength(continuous)
Predictors:
  1
  SepalWidth(continuous)
  (1 + SepalWidth | Species)

julia> subformula(f, 0)
FormulaTerm
Response:
  SepalLength(continuous)
Predictors:
  0
  (1 + SepalWidth | Species)

```
"""
subformula(f::FormulaTerm, id; kwargs...) = subformula(f.lhs, f.rhs, id; kwargs...)
subformula(lhs::AbstractTerm, rhs::MatrixTerm, id::Int; reschema::Bool = false) = 
    id > 0 ? reschema_formula(lhs, rhs.terms[1:id], reschema) : reschema_formula(lhs, (InterceptTerm{false}(),), reschema)
# For nestedmodels & GLM

function subformula(lhs::AbstractTerm, rhs::MatrixTerm, id; reschema::Bool = false)
    terms = rhs.terms[setdiff(eachindex(rhs.terms), id)]
    1 in id && (terms = (InterceptTerm{false}(), terms...))
    reschema_formula(lhs, terms, reschema)
end
# For type 2/3 & GLM

function subformula(lhs::AbstractTerm, rhs::NTuple{N, AbstractTerm}, id::Int; rhs_id::Int = 1, reschema::Bool = false) where N
    ts = collect(Any, rhs)
    ts[rhs_id] = id > 0 ? collect_matrix_terms(ts[rhs_id].terms[1:id]) : collect_matrix_terms((InterceptTerm{false}(),))
    if reschema
        lhs = clear_schema(lhs)
        ts = clear_schema.(ts)
    end
    FormulaTerm(lhs, tuple(ts...)) 
end
# For nestedmodels & MixedModels

"""
    extract_contrasts(f::FormulaTerm)

Extract a dictionary of contrasts. The keys are symbols of term; the values are contrasts (`AbstractContrasts`).
"""
extract_contrasts(f::FormulaTerm) = 
    Dict{Symbol, Any}(t.sym => t.contrasts.contrasts for t in f.rhs.terms if isa(t, CategoricalTerm))

"""
    clear_schema(terms_with_schema) -> terms_without_schema

Clear any applied schema on terms.
"""
clear_schema(::InterceptTerm{true}) = ConstantTerm(1)
clear_schema(::InterceptTerm{false}) = ConstantTerm(0)
clear_schema(t::FunctionTerm) = t
clear_schema(t::Union{CategoricalTerm, ContinuousTerm}) = Term(t.sym)
clear_schema(t::InteractionTerm) = InteractionTerm(clear_schema.(t.terms))
function clear_schema(t::MatrixTerm) 
    ts = ntuple(i -> clear_schema(t.terms[i]), length(t.terms))
    length(ts) ≡ 1 ? ts[1] : ts
end

# reschema only happen when using TupleTerm rather than MatrixTerm
reschema_formula(lhs::AbstractTerm, ts::TupleTerm, reschema::Bool) = 
    reschema ? FormulaTerm(clear_schema(lhs), clear_schema.(ts)) : FormulaTerm(lhs, collect_matrix_terms(ts))