"""
    qseries(factors, range; prefactor=(), sqrt_prefactor=false, alternating=false, sign=1)

Build a compact `FactorialSum` from affine factorial triples `(a,b,c)` representing
`[a*z+b]!^c`. Construction retains the rule without evaluating the sum.
Example: `qeval(qseries([(2,0,1)], 1:3); q=1)` is `2!+4!+6! = 746`.
"""
qseries(factors::Union{Tuple,AbstractVector}, r::UnitRange{<:Integer}; kwargs...) =
    FactorialSum(r; factors=factors, kwargs...)

"""
    qeval(rule::FactorialSum; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate a finite factorial rule, classically when no target is supplied (`k` and `q` are exclusive).
Classical and level values use ratio kernels, exact classical values integer arithmetic, and real or
complex q scaled ratios with adaptive precision; `workspace` reuses fixed-q tables. Factorial arguments
must be nonnegative, and individually polar terms are rejected. `qeval(Symbolic(), rule)` returns the
formula as a [`QForm`](@ref); with Nemo, `xform` expands it.
"""
function qeval end

Base.@constprop :aggressive function qeval(s::FactorialSum; k=nothing,q=nothing,exact::Bool=false,
               T::Type{TT}=Float64,workspace=nothing) where {TT}
    _validate_rule(s)
    q = _evaluation_q(k,q,exact)
    if _is_classical(q)
        return exact ? classical_exact(s) : classical_value(s,T;workspace=workspace)
    elseif !isnothing(k)
        k isa Integer || throw(ArgumentError("k must be an integer level"))
        kk=Int(k)
        kk>=0 || throw(DomainError(k,"level must be nonnegative"))
        exact && return qform(s; k = kk)       # see `_level_value`
        return convert(T, value_at_level(s,kk,T;workspace=workspace))
    elseif !isnothing(q)
        exact && !_is_classical(q) && throw(ArgumentError("exact rule evaluation requires k or q=1"))
        return analytic_value(s,q,T;workspace=workspace)
    end
    throw(ArgumentError("specify k or q to evaluate a factorial rule"))
end
