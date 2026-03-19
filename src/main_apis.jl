# ==============================================================================
# Global Cache Management
# ==============================================================================

"""
    clear_caches!()

Empties all internal LRU caches (numeric, exact Nemo phases, and cyclo sieves).
Highly recommended to call this when looping over completely different manifolds 
or changing the level `k` drastically to free up RAM.
"""
function clear_caches!()
    empty!(LOGQFACT_CACHE)
    empty!(Q6J_NUMERIC_CACHE)
    empty!(GLOBAL_SIEVE_CACHE)
    empty!(UNIT_CIRCLE_CACHE)
    empty!(ANALYTIC_CACHE)
    empty!(EXACT_PHI_CACHE) # The new Nemo exact sieve cache
    return nothing
end

# Define a single, zero-allocation empty result for topological zeros
const EMPTY_CYCLO_RESULT = CycloResult(
    EMPTY_MONOMIAL, EMPTY_MONOMIAL, EMPTY_MONOMIAL, CycloMonomial[], 0:-1, 0
)

# ==============================================================================
# Quantum 6j Symbol API
# ==============================================================================
# export q6j, q3j

"""
    q6j(j1, j2, j3, j4, j5, j6, [k]; mode=:cyclo, T=Float64, prec=256)

Evaluates the quantum 6j-symbol for the SU(2)_k fusion category.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    # 1. k-Independent Topologies
    if mode == :cyclo || mode == :classical || mode == :classical_exact
        if !δtet(j1, j2, j3, j4, j5, j6)
            mode == :cyclo && return EMPTY_CYCLO_RESULT
            mode == :classical_exact && return ClassicalResult(0, 0//1)
            return 0.0
        end
        
        mode == :cyclo && return q6j_cyclo(j1, j2, j3, j4, j5, j6)
        mode == :classical && return qracah6j_classical(j1, j2, j3, j4, j5, j6)
        mode == :classical_exact && return qracah6j_classical_exact(j1, j2, j3, j4, j5, j6)
    end
    
    # 2. Level k Required
    isnothing(k) && throw(ArgumentError("Mode :$mode requires a level k."))

    # 3. Quantum Admissibility
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        if mode == :exact
            _, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            return ExactResult(k, EMPTY_MONOMIAL, zero(z))
        else
            return zero(T)
        end
    end

    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    # 4. Engine Dispatch
    if mode == :exact
        # The ultimate unification: cyclo -> exact!
        res_cyclo = q6j_cyclo(c_spins[1], c_spins[2], c_spins[3], c_spins[4], c_spins[5], c_spins[6])
        return evaluate_level_exact(res_cyclo, k)
        
    elseif mode == :numeric
        return get!(Q6J_NUMERIC_CACHE, (c_spins, k, prec)) do
            if T === BigFloat
                return setprecision(BigFloat, prec) do
                    model = NumericSU2kModel(k; T=T, prec=prec)
                    _qracah6j_stable(model, c_spins[1], c_spins[2], c_spins[3], c_spins[4], c_spins[5], c_spins[6])
                end
            else
                model = NumericSU2kModel(k; T=T, prec=prec)
                return _qracah6j_stable(model, c_spins[1], c_spins[2], c_spins[3], c_spins[4], c_spins[5], c_spins[6])
            end
        end
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end

# ==============================================================================
# Quantum 3j Symbol API
# ==============================================================================

"""
    q3j(j1, j2, j3, m1, m2, [m3], [k]; mode=:cyclo)
"""
q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin; kwargs...) = 
    q3j(j1, j2, j3, m1, m2, -m1-m2, nothing; kwargs...)

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::OptInt=nothing; 
             mode=nothing, T::Type{<:AbstractFloat}=Float64, prec=256)
    
    mode = isnothing(mode) ? (isnothing(k) ? :cyclo : :numeric) : mode

    if mode == :cyclo || mode == :classical || mode == :classical_exact
        if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
            mode == :cyclo && return EMPTY_CYCLO_RESULT
            mode == :classical_exact && return ClassicalResult(0, 0//1)
            return 0.0
        end
        
        mode == :cyclo && return q3j_cyclo(j1, j2, j3, m1, m2, m3)
        mode == :classical && return qracah3j_classical(j1, j2, j3, m1, m2, m3)
        mode == :classical_exact && return qracah3j_classical_exact(j1, j2, j3, m1, m2, m3)
    end
    
    isnothing(k) && throw(ArgumentError("Mode :$mode requires a level k."))

    if !qδ(j1, j2, j3, k) || !iszero(m1 + m2 + m3) 
        if mode == :exact
            _, z = Nemo.cyclotomic_field(2 * (k + 2), "ζ")
            return ExactResult(k, EMPTY_MONOMIAL, zero(z))
        else
            return zero(T)
        end
    end

    if mode == :exact
        res_cyclo = q3j_cyclo(j1, j2, j3, m1, m2, m3)
        return evaluate_level_exact(res_cyclo, k)
        
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    else
        throw(ArgumentError("Unknown mode: $mode"))
    end
end


export cyclo_to_numeric


"""
    cyclo_to_numeric(res::CycloResult, ::Type{T}=Float64; k=nothing, q=nothing, theta=nothing, prec=512)

The master numeric dispatcher for `QRacahSymbols`.
Takes a compiled `CycloResult` and safely projects it into the requested numeric regime.

# Targeting Options
- `k=val`: Fast discrete evaluation at SU(2)_k root of unity (Topological regime).
- `theta=val`: Fast continuous evaluation on the unit circle q = exp(iθ).
- `q=val`: Full analytic continuation to arbitrary complex parameter q ∈ C (SL(2,C) regime).
"""
function cyclo_to_numeric(res::CycloResult, ::Type{T}=Float64; k=nothing, q=nothing, theta=nothing, prec=512) where {T}
    targets_defined = (!isnothing(k)) + (!isnothing(q)) + (!isnothing(theta))
    if targets_defined != 1
        throw(ArgumentError("Specify exactly one evaluation target: k (integer), theta (real), or q (complex)."))
    end

    if !isnothing(k)
        return evaluate_level(res, Int(k), T; prec=prec)
    elseif !isnothing(theta)
        return evaluate_unit_circle(res, theta, T; prec=prec)
    else
        return evaluate_analytic(res, q, T; prec=prec)
    end
end