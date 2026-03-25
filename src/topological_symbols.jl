
# ----------------------------------------------------------------------
#                   ---- Topological Symbols -----
# Internal engines for higher-order composite tensors:
# Quantum Integers, Dimensions, R-Matrices, F-Symbols, and G-Symbols.
# ----------------------------------------------------------------------
"""
    _multiply_prefactors(res::CycloResult, extra_sq::CycloMonomial)

Safely fuses a set of extra quantum dimensions into a 6j-symbol's prefactor remainder.
Automatically extracts any newly formed algebraic perfect squares into `root`.
"""
function _multiply_prefactors(res::CycloResult, extra_sq::CycloMonomial)
    res.base_term.sign == 0 && return res # Fast path for topological zeros
    
    max_d = 0
    !isempty(res.radical.exps) && (max_d = max(max_d, res.radical.exps[end].first))
    !isempty(extra_sq.exps) && (max_d = max(max_d, extra_sq.exps[end].first))
    
    buf = SymbolicBuffer(max_d)
    buf.sign = res.radical.sign * extra_sq.sign
    buf.z_pow = res.radical.z_pow + extra_sq.z_pow
    
    @inbounds for (d, e) in res.radical.exps; buf.exps[d] += e; end
    @inbounds for (d, e) in extra_sq.exps; buf.exps[d] += e; end
    
    bonus_root, final_rad = snapshot_square_root(buf)
    final_root = res.root * bonus_root
    
    return CycloResult(final_root, final_rad, res.base_term, res.ratios, res.z_range, max(res.max_d, max_d))
end

# -----------------------------------------------------
#  --- Quantum Integers [n]_q and Dimensions [2j+1]_q
# ------------------------------------------------------

function _qint_cyclo(n::Int)
    n == 0 && return CycloMonomial(0, 0, Pair{Int,Int}[]) 
    n == 1 && return EMPTY_MONOMIAL
    buf = SymbolicBuffer(n)
    add_qint!(buf, n, 1)
    return snapshot(buf)
end

"""
    _qint_exact(n::Int, k::Int)

Internal fast-path for exact cyclotomic integers. Returns a raw `nf_elem`.
Uses EXACT_MODEL_CACHE if available, otherwise computes in isolated O(1) memory 
using Nemo's highly optimized `divexact`.
"""
function _qint_exact(n::Int, k::Int)
    if haskey(EXACT_MODEL_CACHE, k)
        return EXACT_MODEL_CACHE[k].q_ints[n+1]
    end
    
    h = k + 2
    K, z = Nemo.cyclotomic_field(2 * h, "ζ")
    
    # Must use K(0) and K(1) to guarantee type stability!
    n == 0 && return K(0)
    n == 1 && return K(1)
    
    num = z^(2n) - K(1)
    den = z^2 - K(1)
    return Nemo.divexact(num, den) * z^(1-n)
end

"""
    _qint_numeric(n::Int, k::Int, T::Type, prec::Int)

Evaluates [n]_q directly using arbitrary precision.
"""
function _qint_numeric(n::Int, k::Int, ::Type{T}, prec::Int) where {T <: AbstractFloat}
    n == 0 && return zero(T)
    n == 1 && return one(T)
    
    return setprecision(BigFloat, prec) do
        h_bf = BigFloat(k + 2)
        val = sinpi(BigFloat(n) / h_bf) / sinpi(one(BigFloat) / h_bf)
        return T(val)
    end
end

# ------------------------------------------
#  -- R-Matrix (Braiding / Framing Phase)
# ------------------------------------------

function rmatrix_cyclo(j1::Spin, j2::Spin, j3::Spin)
    z_pow_val = round(Int, 2 * (j3*(j3+1) - j1*(j1+1) - j2*(j2+1)))
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    return CycloMonomial(s, z_pow_val, Pair{Int,Int}[])
end

function rmatrix_exact(j1::Spin, j2::Spin, j3::Spin, k::Int)
    h = k + 2
    K, z = Nemo.cyclotomic_field(2 * h, "ζ")
    z_pow_val = round(Int, 2 * (j3*(j3+1) - j1*(j1+1) - j2*(j2+1)))
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    res = s == 1 ? z^z_pow_val : -(z^z_pow_val)
    return CycloExactResult(k, EMPTY_MONOMIAL, res)
end

function rmatrix_numeric(j1::Spin, j2::Spin, j3::Spin, k::Int; T::Type{<:AbstractFloat}=Float64)
    phase_val = 2 * (j3*(j3+1) - j1*(j1+1) - j2*(j2+1))
    s = iseven(round(Int, j1 + j2 - j3)) ? 1 : -1
    return s * cispi(T(phase_val) / T(2 * (k + 2))) 
end


# ------ F-Symbol (Unitary Fusion Crossing) ------

function fsymbol_cyclo(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    res_6j.base_term.sign == 0 && return res_6j 
    
    buf = SymbolicBuffer(round(Int, max(2j3+1, 2j6+1)))
    add_qint!(buf, round(Int, 2j3+1), 1)
    add_qint!(buf, round(Int, 2j6+1), 1)
    dims_mono = snapshot(buf)
    
    res_f = _multiply_prefactors(res_6j, dims_mono)
    
    phase_s = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
    new_base_term = CycloMonomial(res_f.base_term.sign * phase_s, res_f.base_term.z_pow, copy(res_f.base_term.exps))
    
    return CycloResult(res_f.root, res_f.radical, new_base_term, res_f.ratios, res_f.z_range, res_f.max_d)
end

function fsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _q6j_stable(model, j1, j2, j3, j4, j5, j6)
    phase = iseven(round(Int, j1 + j2 + j4 + j5)) ? 1 : -1
    return phase * sqrt(qdim_numeric(model, j3) * qdim_numeric(model, j6)) * val_6j
end


#  ----- G-Symbol (Symmetric Tetrahedral Weight) -----

function gsymbol_cyclo(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)
    res_6j = q6j_cyclo(j1, j2, j3, j4, j5, j6)
    res_6j.base_term.sign == 0 && return res_6j
    
    max_spin = max(j1, j2, j3, j4, j5, j6)
    buf = SymbolicBuffer(round(Int, 2 * max_spin + 1))
    
    add_qint!(buf, round(Int, 2j1 + 1), 1); add_qint!(buf, round(Int, 2j2 + 1), 1); add_qint!(buf, round(Int, 2j3 + 1), 1)
    add_qint!(buf, round(Int, 2j4 + 1), 1); add_qint!(buf, round(Int, 2j5 + 1), 1); add_qint!(buf, round(Int, 2j6 + 1), 1)
    dims_mono = snapshot(buf)
    
    return _multiply_prefactors(res_6j, dims_mono)
end

function gsymbol_numeric(model::NumericSU2kModel{T}, j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin)::T where {T}
    val_6j = _q6j_stable(model, j1, j2, j3, j4, j5, j6)
    dims_prod = qdim_numeric(model, j1) * qdim_numeric(model, j2) * qdim_numeric(model, j3) * qdim_numeric(model, j4) * qdim_numeric(model, j5) * qdim_numeric(model, j6)
    return val_6j * sqrt(dims_prod)
end