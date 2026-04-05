# ---------------------------------------------------------------------
# File: builder.jl
# Geometric compilation and in-place arithmetic for DCR synthesis.
# ---------------------------------------------------------------------

# --- Buffer Management ---

@inline function ensure_capacity!(buf::CycloBuffer, d_req::Int)
    curr = length(buf.exps)
    if d_req > curr
        new_cap = max(d_req + 10, curr * 2)
        resize!(buf.exps, new_cap)
        @inbounds fill!(view(buf.exps, curr+1:new_cap), 0)
    end
end

@inline function reset!(buf::CycloBuffer, sign::Int=1)
    buf.sign = sign
    buf.q_pow = 0
    @inbounds fill!(view(buf.exps, 1:buf.max_d), 0)
    buf.max_d = 0
end

@inline function update_exps!(buf::CycloBuffer, d::Int, p::Int)
    @inbounds buf.exps[d] += p
    d > buf.max_d && (buf.max_d = d)
end

# --- In-Place Arithmetic ---

function mul!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    (a.sign == 0 || b.sign == 0) && (buf.sign = 0; return buf)
    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow + b.q_pow
    @inbounds for (d, e) in a.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in b.phi_exps; update_exps!(buf, d, e); end
    return buf
end

function div!(buf::CycloBuffer, a::CyclotomicMonomial, b::CyclotomicMonomial)
    b.sign == 0 && throw(DivideError())
    a.sign == 0 && (buf.sign = 0; return buf)
    ensure_capacity!(buf, max(a.max_d, b.max_d))
    reset!(buf, a.sign * b.sign)
    buf.q_pow = a.q_pow - b.q_pow
    @inbounds for (d, e) in a.phi_exps; update_exps!(buf, d, e); end
    @inbounds for (d, e) in b.phi_exps; update_exps!(buf, d, -e); end
    return buf
end

Base.:*(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(mul!(CycloBuffer(max(a.max_d, b.max_d)), a, b))
Base.:/(a::CyclotomicMonomial, b::CyclotomicMonomial) = snapshot(div!(CycloBuffer(max(a.max_d, b.max_d)), a, b))


# --- Factorization Engine ---

@inline function add_qint!(buf::CycloBuffer, n::Int, p::Int=1)
    n <= 1 && return
    buf.q_pow += p * (1 - n) 
    ensure_capacity!(buf, n)
    exps = buf.exps # Local ref for SIMD
    @inbounds for d in 2:n
        (n % d == 0) && (exps[d] += p)
    end
    n > buf.max_d && (buf.max_d = n)
end

@inline function add_qfact!(buf::CycloBuffer, n::Int, p::Int=1)
    n <= 1 && return
    buf.q_pow += p * (n * (1 - n) ÷ 2)
    ensure_capacity!(buf, n)
    exps = buf.exps
    @inbounds @simd for d in 2:n
        exps[d] += p * (n ÷ d)
    end
    n > buf.max_d && (buf.max_d = n)
end

# --- Universal Builder ---

function build_dcr!(buf, pre_func, base_func, ratio_func, z_min::Int, z_max::Int)
    reset!(buf)
    pre_func(buf)
    m_root, m_rad = snapshot_square_root(buf)
    g_max_d = buf.max_d

    z_min > z_max && return DCR(m_root, m_rad, CyclotomicMonomial(0,0,[],0), [], z_min:z_max, g_max_d)

    reset!(buf, iseven(z_min) ? 1 : -1)
    base_func(buf, z_min)
    m_base = snapshot(buf)
    g_max_d = max(g_max_d, buf.max_d)

    ratios = Vector{CyclotomicMonomial}(undef, z_max - z_min)
    for (i, z) in enumerate(z_min:z_max-1)
        reset!(buf, -1) 
        ratio_func(buf, z)
        ratios[i] = snapshot(buf)
        buf.max_d > g_max_d && (g_max_d = buf.max_d)
    end

    return DCR(m_root, m_rad, m_base, ratios, z_min:z_max, g_max_d)
end



# --- Symbols ---

function qtriangle!(buf, j1, j2, j3)
    add_qfact!(buf, Int(j1 + j2 - j3))
    add_qfact!(buf, Int(j1 - j2 + j3))
    add_qfact!(buf, Int(-j1 + j2 + j3))
    add_qfact!(buf, Int(j1 + j2 + j3 + 1), -1)
end

function qtetrahedron!(buf, j1, j2, j3, j4, j5, j6)
    qtriangle!(buf, j1, j2, j3); qtriangle!(buf, j1, j5, j6)
    qtriangle!(buf, j2, j4, j6); qtriangle!(buf, j3, j4, j5)
end

function q3j_dcr(j1, j2, j3, m1, m2, m3 = -m1 - m2)
    α = (Int(j3 - j2 + m1), Int(j3 - j1 - m2))
    β = (Int(j1 + j2 - j3), Int(j1 - m1), Int(j2 + m2))
    z_min, z_max = max(0, -α[1], -α[2]), min(β[1], β[2], β[3])
    buf = CycloBuffer(max(z_max + 2, Int(j1 + j2 + j3 + 1)))

    return let α=α, β=β, j1=j1, j2=j2, j3=j3, m1=m1, m2=m2, m3=m3
        build_dcr!(buf,
            b -> (qtriangle!(b, j1,j2,j3); [add_qfact!(b, Int(j1±m1)) for m1 in (m1,-m1)]; [add_qfact!(b, Int(j2±m2)) for m2 in (m2,-m2)]; [add_qfact!(b, Int(j3±m3)) for m3 in (m3,-m3)]),
            (b, z) -> (add_qfact!(b, z, -1); for a in α add_qfact!(b, a+z, -1) end; for b_val in β add_qfact!(b, b_val-z, -1) end),
            (b, z) -> (for b_val in β add_qint!(b, b_val-z) end; add_qint!(b, z+1, -1); for a in α add_qint!(b, a+z+1, -1) end),
            z_min, z_max
        )
    end
end

function q6j_dcr(j1, j2, j3, j4, j5, j6)
    α = (Int(j1+j2+j3), Int(j1+j5+j6), Int(j2+j4+j6), Int(j3+j4+j5))
    β = (Int(j1+j2+j4+j5), Int(j1+j3+j4+j6), Int(j2+j3+j5+j6))
    z_min, z_max = max(α...), min(β...)
    buf = CycloBuffer(max(z_max + 2, β...))

    return let α=α, β=β, j1=j1, j2=j2, j3=j3, j4=j4, j5=j5, j6=j6
        build_dcr!(buf,
            b -> qtetrahedron!(b, j1, j2, j3, j4, j5, j6),
            (b, z) -> (add_qfact!(b, z+1); for a in α add_qfact!(b, z-a, -1) end; for b_v in β add_qfact!(b, b_v-z, -1) end),
            (b, z) -> (add_qint!(b, z+2); for b_v in β add_qint!(b, b_v-z) end; for a in α add_qint!(b, z+1-a, -1) end),
            z_min, z_max
        )
    end
end