# src/QRacahSymbols.jl
module QRacahSymbols

using Nemo
using LRUCache

# Data Structures and Types
export CycloMonomial, ExactSU2kModel, NumericSU2kModel
export GenericResult, ExactResult
export Spin

# Internal Handlers for advanced usage
export qracah6j_generic, qracah6j_exact, qracah6j_numeric, qracah6j_classical
export qracah3j_generic, qracah3j_exact, qracah3j_numeric, qracah3j_classical


# Primary API Functions 
export q6j, q3j, fsymbol, rmatrix, gsymbol, qint, qdim
export evaluate_exact, evaluate_classical, evaluate_symbolic

# Include architecture files
include("types.jl")
include("admissibility.jl")
include("symmetry.jl")
include("generic_engine.jl")
include("exact_engine.jl")
include("numerics_engine.jl")
include("tqft.jl")

# Initialize the global Evaluation & Model Caches
const Q6J_NUMERIC_CACHE = LRU{Tuple{NTuple{6, Float64}, Int, DataType}, Any}(maxsize=50000)
const Q6J_EXACT_CACHE   = LRU{Tuple{NTuple{6, Float64}, Int}, ExactResult}(maxsize=10000)
const EXACT_MODEL_CACHE = LRU{Int, ExactSU2kModel}(maxsize=20)


# ============================================================
# 1. Quantum 6j Symbol
# ============================================================

"""
    q6j(j1, j2, j3, j4, j5, j6; mode=:generic)

Evaluates the k-independent 6j-symbol. 
Valid modes: `:generic` (Symbolic CycloMonomials) or `:classical` (Ponzano-Regge limit).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; mode=:generic)
    # Classical/Generic admissible conditions (independent of k)
    if !δtet(j1, j2, j3, j4, j5, j6)
        if mode == :generic return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[]) end
        if mode == :classical return 0.0 end
    end

    if mode == :classical
        return qracah6j_classical(j1, j2, j3, j4, j5, j6)
    elseif mode == :generic
        return qracah6j_generic(j1, j2, j3, j4, j5, j6)
    else
        error("Mode :$mode requires a level `k`. Call as `q6j(..., k; mode=:$mode)`.")
    end
end

"""
    q6j(j1, j2, j3, j4, j5, j6; k::Int, mode=:numeric)

Evaluates the k-independent 6j-symbol. 
Valid modes: `:generic` (Symbolic CycloMonomials) or `:classical` (Ponzano-Regge limit).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic || mode == :classical
        return q6j(j1, j2, j3, j4, j5, j6; mode=mode) 
    end

    # Quantum admissible conditions
    if !qδtet(j1, j2, j3, j4, j5, j6, k) 
        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return ExactResult(k, model.K(0), model.K(0))
        else
            return zero(T)
        end
    end

    c_spins = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    if mode == :exact
        return get!(Q6J_EXACT_CACHE, (c_spins, k)) do
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            qracah6j_exact(model, c_spins...)
        end
    elseif mode == :numeric
        return get!(Q6J_NUMERIC_CACHE, (c_spins, k, T)) do
            model = NumericSU2kModel(k; T=T, prec=prec)
            _qracah6j_stable(model, c_spins...)
        end
    else
        error("Unknown mode: $mode")
    end
end

# ============================================================
# 2. Quantum 3j Symbol
# ============================================================

# Auto-completes m3 if omitted
q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin; mode=:generic) = 
    q3j(j1, j2, j3, m1, m2, -m1-m2; mode=mode)

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin; mode=:generic)
    # Classical/Generic Gatekeeper
    if !δ(j1, j2, j3) || !iszero(m1 + m2 + m3)
        if mode == :generic return GenericResult(CycloMonomial(0, 0, Int[]), CycloMonomial[]) end
        if mode == :classical return 0.0 end
    end

    if mode == :classical
        return qracah3j_classical(j1, j2, j3, m1, m2, m3)
    elseif mode == :generic
        return qracah3j_generic(j1, j2, j3, m1, m2, m3)
    else
        error("Mode :$mode requires a level `k`. Call as `q3j(..., m3, k; mode=:$mode)`.")
    end
end

function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin, k::Int; mode=:numeric, T::Type{<:AbstractFloat}=Float64, prec=256)
    if mode == :generic || mode == :classical
        return q3j(j1, j2, j3, m1, m2, m3; mode=mode)
    end
    
    # Quantum Gatekeeper
    if !qδ(j1, j2, j3, k) || !iszero(m1 + m2 + m3)
        if mode == :exact
            model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
            return ExactResult(k, model.K(0), model.K(0))
        else
            return zero(T)
        end
    end

    if mode == :exact
        model = get!(() -> ExactSU2kModel(k), EXACT_MODEL_CACHE, k)
        return qracah3j_exact(model, j1, j2, j3, m1, m2)
    elseif mode == :numeric
        model = NumericSU2kModel(k; T=T, prec=prec)
        return _qracah3j_stable(model, j1, j2, j3, m1, m2)
    else
        error("Unknown mode: $mode")
    end
end


# ============================================================
# Quantum Dimensions
# ============================================================

"""
    qint(n::Int; mode=:generic)
    qint(n::Int, k::Int; mode=:exact, T=Float64)

Returns a representation of the quantum integer [n]_q. 
The `mode` determines the choice of representation:
- `:generic` : Returns a symbolic `CycloMonomial` factorization (requires only `n`).
- `:exact`   : Returns the exact algebraic representation in the cyclotomic field Q(ζ).
- `:numeric` : Returns the numerical floating-point value using the stable sine formula.
"""
function qint end



# ============================================================
# Master API: TQFT Symbols
# ============================================================

"""
    qdim(j, [k]; mode=:generic)
    fsymbol(j1, j2, j3, j4, j5, j6, [k]; mode=:numeric)
    rmatrix(j1, j2, j3, [k]; mode=:generic)

Access standard SU(2)_k category data.
"""
function qdim end
function fsymbol end
function rmatrix end
function gsymbol end



end # module