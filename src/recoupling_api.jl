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
    if !isnothing(k) && !qδtet(js..., k)
        dcr = ZERO_DCR
    else
        dcr = q6j_dcr(js...)
    end
    
    # --- return raw graph if no target is specified ---
    if isnothing(k) && isnothing(q) 
        return dcr
    end

    # --- DCR Projections ---
    return project_dcr(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T=Float64)
Returns the Wigner 3j-symbol.
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2; 
             k=nothing, q=nothing, exact::Bool=false, eager::Bool=false, T::Type=Float64)
             
    # --- eager evaluation circuit ---
    if eager && !isnothing(k) && isnothing(q)
        return exact ? q3j_exact(j1, j2, j3, m1, m2, m3, k) : q3j_direct(j1, j2, j3, m1, m2, m3, k, T)
    end

    # --- DCR Construction ---
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        dcr = ZERO_DCR
    else
        dcr = q3j_dcr(j1, j2, j3, m1, m2, m3)
    end
    
    # --- Return raw graph if no target is specified ---
    if isnothing(k) && isnothing(q)
        return dcr
    end

    # --- DCR Projections ---
    return project_dcr(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Unitary crossing matrix element: √([d3][d6]) * {6j}.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = fsymbol_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return project_dcr(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Tetrahedrally symmetric invariant: √(Π[di]) * {6j}.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = gsymbol_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return project_dcr(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Evaluates the standard closed tetrahedron network.
"""
function tetrahedron(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = tetrahedron_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    if isnothing(k) && isnothing(q)
        return dcr
    end

    return project_dcr(dcr; k=k, q=q, exact=exact, T=T)
end


"""
    theta_value(j1, j2, j3; k=nothing, q=nothing, exact::Bool=false, T=Float64)
Value of the Theta-graph.
"""
function theta_value(j1::Spin, j2::Spin, j3::Spin; 
                     k=nothing, q=nothing, exact::Bool=false, T::Type=Float64)
    
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        mono = ZERO_MONOMIAL
    else
        mono = theta_mono(j1, j2, j3) 
    end
    
    # Classical limit catch
    if (!isnothing(q) && (q == 1 || q == 1.0)) || (!isnothing(k) && k == Inf)
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
    mono = qdim_mono(j)
    
    if (!isnothing(q) && (q == 1 || q == 1.0)) || (!isnothing(k) && k == Inf)
        return exact ? project_classical_exact(mono) : project_classical(mono, T)
    end
    
    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono
    
    return exact ? project_exact(mono, k) : project_discrete(mono, k, T)
end


"""
    rmatrix(j1, j2, j3; k=nothing, q=nothing, exact::Bool=false, T=ComplexF64)
Returns the R-matrix phase (braiding eigenvalue) for j1, j2 crossing into j3.
"""
function rmatrix(j1::Spin, j2::Spin, j3::Spin; 
                 k=nothing, q=nothing, exact::Bool=false, T::Type=ComplexF64)
    
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        mono = ZERO_MONOMIAL
    else
        mono = rmatrix_mono(j1, j2, j3) 
    end
    
    # Classical limit catch
    if (!isnothing(q) && (q == 1 || q == 1.0)) || (!isnothing(k) && k == Inf)
        return T(mono.sign)
    end
    
    # Analytic
    if !isnothing(q)
        q_C = Complex{typeof(q * 1.0)}(q)
        q_quarter = sqrt(sqrt(q_C)) 
        return T(mono.sign * (q_quarter^mono.q_pow))
    end
    
    isnothing(k) && return mono
    
    # Root of Unity
    if exact
        @info "Note: The exact R-matrix phase involves q^{1/4} and is returned in the expanded cyclotomic field ζ_{$(4k+8)}." maxlog=1
        return project_exact(mono, k)
    else
        h = k + 2
        phase_angle = (pi * mono.q_pow) / (4 * h)
        return T(mono.sign * cis(phase_angle))
    end
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