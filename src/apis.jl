# ==============================================================================
# File: api.jl
# Unified Public API for SU(2) TQFT Kernels
# ==============================================================================

"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, mode=:discrete, eager=false, T=Float64)
    
Returns the Quantum 6j-symbol. 
If `k` is nothing, returns a `DCR` object. 
If `k` is provided, projects according to `mode`.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
             k=nothing, mode=:discrete, eager=false, T::Type=Float64)
    
    js = canonical_spins(j1, j2, j3, j4, j5, j6)
    dcr = q6j_dcr(js...)
    
    isnothing(k) && return dcr
    
    if mode == :discrete
        return eager ? q6j_direct(js..., k; T=T) : project_discrete(dcr, k, T)
    elseif mode == :exact
        return eager ? q6j_exact(js..., k) : project_exact(dcr, k)
    elseif mode == :classical
        return eager ? project_classical_exact(dcr) : project_classical(dcr, T)
    end
end

"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, mode=:discrete, T=Float64)
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2; 
             k=nothing, mode=:discrete, T::Type=Float64)
    
    dcr = q3j_dcr(j1, j2, j3, m1, m2, m3)
    isnothing(k) && return dcr
    
    if mode == :discrete
        return project_discrete(dcr, k, T)
    elseif mode == :exact
        return project_exact(dcr, k)
    elseif mode == :classical
        return project_classical(dcr, T)
    end
end

"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, T=Float64)
Unitary crossing matrix element: √([d3][d6]) * {6j}.
Returns `DCR` if `k` is nothing.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, T::Type=Float64, mode=:discrete)
    
    dcr = fsymbol_dcr(j1, j2, j3, j4, j5, j6)
    isnothing(k) && return dcr
    
    if mode == :exact
        return project_exact(dcr, k)
    else
        return project_discrete(dcr, k, T)
    end
end

"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, T=Float64)
Tetrahedrally symmetric invariant: √(Π[di]) * {6j}.
Returns `DCR` if `k` is nothing.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, T::Type=Float64, mode=:discrete)
    
    dcr = gsymbol_dcr(j1, j2, j3, j4, j5, j6)
    isnothing(k) && return dcr
    
    if mode == :exact
        return project_exact(dcr, k)
    else
        return project_discrete(dcr, k, T)
    end
end

"""
    theta_value(j1, j2, j3; k=nothing, T=Float64)
Value of the Theta-graph. Returns `CyclotomicMonomial` if `k` is nothing.
"""
function theta_value(j1::Spin, j2::Spin, j3::Spin; k=nothing, T::Type=Float64)
    mono = theta_dcr(j1, j2, j3) # Returns a Monomial
    isnothing(k) && return mono
    return project_discrete(mono, k, T)
end

"""
    qdim(j; k=nothing, T=Float64)
Quantum dimension [2j+1]_q. Returns `CyclotomicMonomial` if `k` is nothing.
"""
function qdim(j::Spin; k=nothing, T::Type=Float64)
    mono = qdim_mono(j)
    isnothing(k) && return mono
    return project_discrete(mono, k, T)
end









# # ==============================================================================
# # File: api.jl
# # Unified Public API for SU(2) TQFT Symbols
# # ==============================================================================

# """
#     q6j(j1, j2, j3, j4, j5, j6; k, mode=:discrete, eager=false, T=Float64)

# Unified API for the Quantum 6j-symbol.
# - `mode=:discrete`  -> Evaluation at root of unity q = exp(iπ/(k+2)).
# - `mode=:classical` -> Evaluation at q=1 (Ponzano-Regge limit).
# - `mode=:exact`     -> Algebraic evaluation in Nemo cyclotomic field.
# - `eager=true`      -> Use direct summation (faster for single calls).
# - `eager=false`     -> Build DCR first (faster for repeated symbolic terms).
# """
# function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
#              k=nothing, mode=:discrete, eager=false, T::Type=Float64)
    
#     # 1. Canonicalize to maximize cache hits
#     js = canonical_spins(j1, j2, j3, j4, j5, j6)
    
#     # 2. Admissibility Guard
#     if mode == :discrete || mode == :exact
#         isnothing(k) && error("Level k required for mode $mode")
#         !qδtet(js..., k) && return zero(T)
#     else
#         !δtet(js...) && return zero(T)
#     end

#     # 3. Routing
#     if mode == :discrete
#         return eager ? q6j_direct(js..., k; T=T) : project_discrete(q6j_dcr(js...), k, T)
#     elseif mode == :exact
#         return eager ? q6j_exact(js..., k) : project_exact(q6j_dcr(js...), k)
#     elseif mode == :classical
#         # Exact rational classical is usually preferred over numerical
#         return eager ? project_classical_exact(q6j_dcr(js...)) : project_classical(q6j_dcr(js...), T)
#     end
# end

# """
#     q3j(j1, j2, j3, m1, m2, m3; k, mode=:discrete)
# Unified API for the Quantum 3j-symbol (Wigner coefficient).
# """
# function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2; 
#              k=nothing, mode=:discrete, T::Type=Float64)
#     # 3j symbols aren't tetrahedral, so we skip canonical_spins here 
#     # unless you implement Regge symmetries for 3j.
#     if mode == :discrete
#         isnothing(k) && error("Level k required")
#         return project_discrete(q3j_dcr(j1, j2, j3, m1, m2, m3), k, T)
#     elseif mode == :exact
#         return project_exact(q3j_dcr(j1, j2, j3, m1, m2, m3), k)
#     end
# end

# """
#     fsymbol(j1, j2, j3, j4, j5, j6; k)
# Unitary crossing matrix element (Normalized 6j).
# """
# function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; k::Int, T::Type=Float64)
#     return project_discrete(fsymbol_dcr(j1, j2, j3, j4, j5, j6), k, T)
# end

# """
#     gsymbol(j1, j2, j3, j4, j5, j6; k)
# Tetrahedrally symmetric invariant (Fully normalized 6j).
# """
# function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; k::Int, T::Type=Float64)
#     return project_discrete(gsymbol_dcr(j1, j2, j3, j4, j5, j6), k, T)
# end

# """
#     theta_value(j1, j2, j3; k)
# Value of the Theta-graph (Triangle dimension).
# """
# function theta_value(j1::Spin, j2::Spin, j3::Spin; k::Int, T::Type=Float64)
#     return project_discrete(theta_dcr(j1, j2, j3), k, T)
# end

# """
#     qdim(j; k)
# Quantum dimension [2j+1]_q.
# """
# function qdim(j::Spin; k::Int, T::Type=Float64)
#     return project_discrete(qdim_mono(j), k, T)
# end