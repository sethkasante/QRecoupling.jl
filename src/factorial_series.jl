"""
    qseries(factors, range; prefactor=(), sqrt_prefactor=false, alternating=false, sign=1)

Build a compact `FactorialSum` from affine factorial triples `(a,b,c)` representing 
`[a*z+b]!^c`. Unlike the callback form of `qseries`, this does not expand a DCR.
Example: `qeval(qseries([(2,0,1)], 1:3); q=1)` is `2!+4!+6! = 746`.
"""
qseries(factors::Union{Tuple,AbstractVector}, r::UnitRange{<:Integer}; kwargs...) =
    FactorialSum(r; factors=factors, kwargs...)

"Explicit DCR view used for exact/analytic projections and root-of-unity cancellations."
function _factorial_dcr(s::FactorialSum)
    is_empty_sum(s) && return ZERO_DCR
    buf = CycloBuffer(max_argument(s))
    d = build_dcr!(buf,
        b -> foreach(p -> add_qfact!(b,Int(p.first),Int(p.second)),s.pre),
        (b,z) -> foreach(f -> add_qfact!(b,_arg(f,z),Int(f.c)),s.fac),
        (b,z) -> begin
            for f in s.fac
                lo,hi,c = _factor_step(f,z)
                c == 0 && continue
                for n in lo:hi
                    add_qint!(b,n,c)
                end
            end
        end,
        s.zlo,s.zhi; extract_radical=s.sqrt_pre, alternating_sign=s.alternating)
    return s.sign0 == 1 ? d : DCR(d.root,d.radical,-d.base,d.ratios,d.z_range,d.max_d)
end

"""
    qeval(rule::FactorialSum; k=nothing, q=nothing, exact=false, T=Float64, workspace=nothing)

Evaluate a finite factorial rule. Classical and in-table level evaluations use ratio
kernels; exact, analytic and out-of-table evaluations use its explicit DCR view.
Factorial arguments must be nonnegative; individually polar terms are rejected even
if the complete sum might have a removable singularity. `k` and `q` are exclusive.
For general functions, use the existing callback `qseries`/DCR interface instead.
"""
function qeval(s::FactorialSum; k=nothing,q=nothing,exact::Bool=false,
               T::Type{TT}=Float64,workspace=nothing) where {TT}
    _validate_rule(s)
    isnothing(k) || isnothing(q) || throw(ArgumentError("specify k or q, not both"))
    if _is_classical(q) && !exact
        return classical_value(s,T;workspace=workspace)
    elseif !isnothing(k)
        k isa Integer || throw(ArgumentError("k must be an integer level"))
        kk=Int(k)
        kk>=0 || throw(DomainError(k,"level must be nonnegative"))
        exact && return project_exact(_factorial_dcr(s),kk)
        return value_at_level(s,kk,T;workspace=workspace,
                             fallback=()->project_discrete(_factorial_dcr(s),kk,T))
    elseif !isnothing(q)
        exact && !_is_classical(q) && throw(ArgumentError("exact rule evaluation requires k or q=1"))
        return qeval(_factorial_dcr(s);q=q,exact=exact,T=T)
    end
    throw(ArgumentError("specify k or q to evaluate a factorial rule"))
end
