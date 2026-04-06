# --------------------------------------------
# Unified Public API for SU(2) TQFT Kernels
# ---------------------------------------------


"""
    q6j(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, mode=:discrete, eager=false, T=Float64)
    
Returns the Quantum 6j-symbol. 
- If `k` and `q` are both nothing, returns a `DCR` object. 
- If `k` is provided, projects according to `mode` (:discrete, :exact, :classical).
- If `q` is provided, evaluates analytically for that specific complex/real q.
"""
function q6j(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
             k=nothing, q=nothing, mode=:discrete, eager=false, T::Type=Float64)
    
    js = canonical_spins(j1, j2, j3, j4, j5, j6)

    # Check level-k triangle admissibility
    if !isnothing(k) && !qδtet(js..., k)
        dcr = ZERO_DCR
    else
        dcr = q6j_dcr(js...)
    end
    
    if mode == :classical
        return eager ? project_classical_exact(dcr) : project_classical(dcr, T)
    end

    !isnothing(q) && return project_analytic(dcr, q)
    isnothing(k) && return dcr
    
    if mode == :discrete
        return eager ? q6j_direct(js..., k; T=T) : project_discrete(dcr, k, T)
    elseif mode == :exact
        return eager ? q6j_exact(js..., k) : project_exact(dcr, k)
    end
end

"""
    q3j(j1, j2, j3, m1, m2, m3; k=nothing, q=nothing, mode=:discrete, T=Float64)
"""
function q3j(j1::Spin, j2::Spin, j3::Spin, m1::Spin, m2::Spin, m3::Spin=-m1-m2; 
             k=nothing, q=nothing, mode=:discrete, T::Type=Float64)
    
    # Check level-k triangle admissibility
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        dcr = ZERO_DCR
    else
        dcr = q3j_dcr(j1, j2, j3, m1, m2, m3)
    end
    
    if mode == :classical
        return project_classical(dcr, T)
    end

    !isnothing(q) && return project_analytic(dcr, q)
    isnothing(k) && return dcr
    
    if mode == :discrete
        return project_discrete(dcr, k, T)
    elseif mode == :exact
        return project_exact(dcr, k)
    end
end



"""
    fsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, T=Float64, mode=:discrete)
Unitary crossing matrix element: √([d3][d6]) * {6j}.
"""
function fsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, T::Type=Float64, mode=:discrete)
    
    # Check level-k triangle admissibility
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = fsymbol_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    !isnothing(q) && return project_analytic(dcr, q)
    isnothing(k) && return dcr
    
    if mode == :exact
        return project_exact(dcr, k)
    else
        return project_discrete(dcr, k, T)
    end
end

"""
    gsymbol(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, T=Float64, mode=:discrete)
Tetrahedrally symmetric invariant: √(Π[di]) * {6j}.
"""
function gsymbol(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                 k=nothing, q=nothing, T::Type=Float64, mode=:discrete)
    
    # Check level-k triangle admissibility
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = gsymbol_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    !isnothing(q) && return project_analytic(dcr, q)
    isnothing(k) && return dcr
    
    if mode == :exact
        return project_exact(dcr, k)
    else
        return project_discrete(dcr, k, T)
    end
end


"""
    theta_value(j1, j2, j3; k=nothing, q=nothing, T=Float64)
Value of the Theta-graph. 
"""
function theta_value(j1::Spin, j2::Spin, j3::Spin; k=nothing, q=nothing, T::Type=Float64)
    
    # Check level-k triangle admissibility
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        mono = ZERO_MONOMIAL
    else
        mono = theta_mono(j1, j2, j3) 
    end

    mono = theta_mono(j1, j2, j3) 
    
    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono
    
    return project_discrete(mono, k, T)
end


"""
    qdim(j; k=nothing, q=nothing, T=Float64)
Quantum dimension [2j+1]_q. 
"""
function qdim(j::Spin; k=nothing, q=nothing, T::Type=Float64)
    mono = qdim_mono(j)
    
    !isnothing(q) && return project_analytic(mono, q)
    isnothing(k) && return mono
    
    return project_discrete(mono, k, T)
end


"""
    tetrahedron(j1, j2, j3, j4, j5, j6; k=nothing, q=nothing, mode=:discrete, T=Float64)

Evaluates the standard tetrahedron network. 
(Note: Depending on conventions, this is often identical to the unnormalized `q6j`).
"""
function tetrahedron(j1::Spin, j2::Spin, j3::Spin, j4::Spin, j5::Spin, j6::Spin; 
                     k=nothing, q=nothing, mode=:discrete, T::Type=Float64)
    
    # Check level-k triangle admissibility
    if !isnothing(k) && !qδtet(j1, j2, j3, j4, j5, j6, k)
        dcr = ZERO_DCR
    else
        dcr = tetrahedron_dcr(j1, j2, j3, j4, j5, j6)
    end
    
    !isnothing(q) && return project_analytic(dcr, q)
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
    rmatrix(j1, j2, j3; k=nothing, q=nothing, T=ComplexF64, mode=:discrete)
Returns the R-matrix phase (braiding eigenvalue) for j1, j2 crossing into j3.
Bypasses the projection engines for direct.
"""
function rmatrix(j1::Spin, j2::Spin, j3::Spin; 
                 k=nothing, q=nothing, T::Type=ComplexF64, mode=:discrete)
    
     # Check level-k triangle admissibility
    if !isnothing(k) && !qδ(j1, j2, j3, k)
        mono = ZERO_MONOMIAL
    else
        mono = rmatrix_mono(j1, j2, j3) 
    end
    
    # analytic 
    if !isnothing(q)
        # force complex type to handle phase
        q_C = Complex{typeof(q * 1.0)}(q)
        # Principal 4th root via double sqrt is branch-cut safe and fast
        q_quarter = sqrt(sqrt(q_C)) 
        return T(mono.sign * (q_quarter^mono.q_pow))
    end
    
    # algebraic
    isnothing(k) && return mono
    
    # 3. Mode Projections
    if mode == :discrete
        # Discrete exact bypass!
        # q = exp(i * pi / h)  =>  q^(1/4) = exp(i * pi / 4h)
        # Therefore: (q^(1/4))^P = cis(pi * P / 4h)
        h = k + 2
        phase_angle = (pi * mono.q_pow) / (4 * h)
        return T(mono.sign * cis(phase_angle))
        
    elseif mode == :exact
       @info "Note: The exact R-matrix phase involves q^{1/4} and is returned in the expanded cyclotomic field ζ_{$(4k+8)}." maxlog=1
       return project_exact(mono, k)
        
    elseif mode == :classical
        # At q=1, the phase is just the parity sign
        return T(mono.sign)
    end
end



#---- clear caches --- 

clear_numeric_caches!() = (empty!(LOGQFACT_CACHE); nothing)

clear_exact_caches!() = (empty!(EXACT_PHI_CACHE); nothing)

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