# --------------------------------------------
# Unified Public API for SU(2) TQFT Kernels
# ---------------------------------------------


"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T=Float64)
    
Returns the Quantum 6j-symbol. 
- If `k` and `q` are omitted, returns the abstract `DCR` object. 
- If `k` is provided, projects to the root of unity (Float64 by default, Cyclotomic if `exact=true`).
- If `q=1`, computes the exact classical Ponzano-Regge limit.
- If `eager=true`, bypasses DCR construction for raw speed (only valid for root of unity `k`).
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
             k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T::Type=Float64)
    
    js = canonical_spins(j1, j2, j3, j4, j5, j6)

    # --- eager evaluation circuit ---
    if eager && !isnothing(k) && isnothing(q)
        return exact ? q6j_exact(js..., k) : q6j_direct(js..., k, T)
    end

    # ---  DCR construction ---
    # Check level-k triangle admissibility
    if !isnothing(k) && !_qδtet(js..., k)
        dcr = ZERO_DCR
    else
        dcr = q6j_dcr(js...)
    end
    
    # --- return raw graph if no target is specified ---
    if isnothing(k) && isnothing(q) 
        return dcr
    end

    # --- DCR Projections ---
    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T=Float64)
Returns the quantum Wigner 3j-symbol.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2; 
             k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T::Type=Float64)
    
    Js = doubled(j1, j2, j3, m1, m2, m3)

    # --- eager evaluation circuit ---
    if eager && !isnothing(k) && isnothing(q)
        return exact ? q3j_exact(Js..., k) : q3j_direct(Js..., k, T)
    end

    # --- DCR Construction ---
    if !isnothing(k) && !_qδ(Js[1], Js[2], Js[3], k)
        dcr = ZERO_DCR
    else
        dcr = q3j_dcr(Js...)
    end
    
    # --- Return raw graph if no target is specified ---
    if isnothing(k) && isnothing(q)
        return dcr
    end

    # --- DCR Projections ---
    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Unitary crossing matrix element: √([d3][d6]) * {6j}.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    # F-symbol not fully symmetric! Just double
    Js = doubled(j1, j2, j3, j4, j5, j6)
    
    if !isnothing(k) && !_qδtet(Js..., k)
        dcr = ZERO_DCR
    else
        dcr = fsymbol_dcr(Js...)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Tetrahedrally symmetric invariant: √(Π[di]) * {6j}.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    # G-symbol is fully symmetric.
    Js = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    if !isnothing(k) && !_qδtet(Js..., k)
        dcr = ZERO_DCR
    else
        dcr = gsymbol_dcr(Js...)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Evaluates the standard closed tetrahedron network.
"""
function tetrahedron(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    # Tetrahedron is symmetric.
    Js = canonical_spins(j1, j2, j3, j4, j5, j6)
    
    if !isnothing(k) && !_qδtet(Js..., k)
        dcr = ZERO_DCR
    else
        dcr = tetrahedron_dcr(Js...)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return qeval(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    theta_value(j1, j2, j3; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Value of the Theta-graph.
"""
function theta_value(j1::Spin, j2::Spin, j3::Spin; 
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    Js = doubled(j1, j2, j3)
    
    if !isnothing(k) && !_qδ(Js..., k)
        mono = ZERO_MONOMIAL
    else
        mono = theta_mono(Js...) 
    end
    
    # classical limit 
    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(mono) : project_classical(mono, T)
    end
    
    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono
    
    return exact ? project_exact(mono, k) : project_discrete(mono, k, T)
end


"""
    qdim(j; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Quantum dimension [2j+1]_q. 
"""
function qdim(j::Spin; k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    J = doubled(j)
    mono = qdim_mono(J)
    
    if !isnothing(q) && (q == 1 || q == 1.0)
        return exact ? project_classical_exact(mono) : project_classical(mono, T)
    end
    
    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono
    
    return exact ? project_exact(mono, k) : project_discrete(mono, k, T)
end

#---- clear caches --- 

clear_numeric_caches!() = (empty!(LOGQFACT_CACHE); nothing)

clear_exact_caches!() = begin
    empty!(EXACT_PHI_CACHE)
    empty!(EXACT_MODEL_CACHE)
    nothing
end

clear_sieve_caches!() = (empty!(MAG_SIEVE_CACHE); nothing)

"""
    empty_caches!()
Useful for freeing RAM during long interactive sessions or resetting state for benchmarking.
"""
function empty_caches!()
    clear_numeric_caches!() 
    clear_exact_caches!()
    clear_sieve_caches!()
    
    @info "QRecoupling internal caches have been successfully cleared." maxlog=1
    return nothing
end