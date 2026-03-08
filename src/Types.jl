# src/Types.jl

# Define a robust type alias for Spins to eliminate method ambiguities
const Spin = Real
const OptInt = Union{Nothing, Int}

"""
    CycloMonomial

Represents a product of cyclotomic polynomials evaluated at q.
"""
struct CycloMonomial
    sign::Int
    z_pow::Int 
    exps::Vector{Int} 
end

function Base.:*(a::CycloMonomial, b::CycloMonomial)
    len_a = length(a.exps)
    len_b = length(b.exps)
    max_len = max(len_a, len_b)
    
    new_exps = zeros(Int, max_len)
    
    @inbounds @simd for i in 1:len_a
        new_exps[i] += a.exps[i]
    end
    @inbounds @simd for i in 1:len_b
        new_exps[i] += b.exps[i]
    end
    
    while !isempty(new_exps) && new_exps[end] == 0
        pop!(new_exps)
    end
    
    return CycloMonomial(a.sign * b.sign, a.z_pow + b.z_pow, new_exps)
end

function Base.:/(a::CycloMonomial, b::CycloMonomial)
    len_a = length(a.exps)
    len_b = length(b.exps)
    max_len = max(len_a, len_b)
    
    new_exps = zeros(Int, max_len)
    
    @inbounds @simd for i in 1:len_a
        new_exps[i] += a.exps[i]
    end
    @inbounds @simd for i in 1:len_b
        new_exps[i] -= b.exps[i]
    end
    
    while !isempty(new_exps) && new_exps[end] == 0
        pop!(new_exps)
    end
    
    return CycloMonomial(a.sign * b.sign, a.z_pow - b.z_pow, new_exps)
end

Base.://(a::CycloMonomial, b::CycloMonomial) = a / b

function Base.inv(M::CycloMonomial)
    M.sign == 0 && throw(DivideError())
    return CycloMonomial(M.sign, -M.z_pow, -M.exps)
end


mutable struct SymbolicBuffer
    sign::Int
    z_pow::Int
    exps::Vector{Int}
end

# Constructor for pre-allocating the vector capacity
SymbolicBuffer(capacity::Int) = SymbolicBuffer(1, 0, zeros(Int, capacity))
# Snapshot function to freeze the buffer into a CycloMonomial
snapshot(buf::SymbolicBuffer) = CycloMonomial(buf.sign, buf.z_pow, copy(buf.exps))

# ----------------------------------------
# Human-Readable Output
# ----------------------------------------

const SUBSCRIPTS = Dict('0'=>'₀', '1'=>'₁', '2'=>'₂', '3'=>'₃', '4'=>'₄', 
                        '5'=>'₅', '6'=>'₆', '7'=>'₇', '8'=>'₈', '9'=>'₉')
const SUPERSCRIPTS = Dict('0'=>'⁰', '1'=>'¹', '2'=>'²', '3'=>'³', '4'=>'⁴', 
                          '5'=>'⁵', '6'=>'⁶', '7'=>'⁷', '8'=>'⁸', '9'=>'⁹', '-'=>'⁻')

to_subscript(n::Int) = map(c -> SUBSCRIPTS[c], string(n))
to_superscript(n::Int) = map(c -> SUPERSCRIPTS[c], string(n))

function Base.show(io::IO, M::CycloMonomial)
    if M.sign == 0
        print(io, "0")
        return
    end

    parts = String[]
    if M.z_pow != 0
        push!(parts, M.z_pow == 1 ? "z" : "z" * to_superscript(M.z_pow))
    end
    
    for d in 1:length(M.exps)
        e = M.exps[d]
        e == 0 && continue
        base_str = "Φ" * to_subscript(d)
        push!(parts, e == 1 ? base_str : base_str * to_superscript(e))
    end
    
    if isempty(parts)
        print(io, M.sign == -1 ? "-1" : "1")
    else
        prefix = M.sign == -1 ? "-" : ""
        print(io, prefix * join(parts, " "))
    end
end

# ----------------------------------------
# Result Structs
# ----------------------------------------

struct GenericResult
    pref_sq::CycloMonomial
    series::Vector{CycloMonomial}
end

function Base.show(io::IO, res::GenericResult)
    print(io, "√(", res.pref_sq, ") × (")
    n_terms = length(res.series)
    if n_terms <= 5
        join(io, res.series, "  +  ")
    else
        print(io, res.series[1], "  +  ", res.series[2], "  +  ... (", n_terms - 4, " more terms) ...  +  ", res.series[end-1], "  +  ", res.series[end])
    end
    print(io, ")")
end

"""
    ExactValue{T}
A self-contained exact algebraic value (like a quantum dimension or integer) 
in the SU(2)_k cyclotomic field.
"""
struct ExactValue{T}
    k::Int
    val::T
end

# Pretty Printing for the REPL
function Base.show(io::IO, ev::ExactValue)
    k_sub = to_subscript(ev.k)
    print(io, "Exact SU(2)$k_sub Value:\n  ")
    print(io, ev.val)
end

struct ExactResult
    k::Int
    pref_sq::nf_elem
    sum_cf::nf_elem
end

function Base.show(io::IO, res::ExactResult)
    k_sub = to_subscript(res.k)
    print(io, "Exact SU(2)$k_sub Symbol:\n")
    print(io, "  Prefactor(Δ²): ", res.pref_sq, "\n")
    print(io, "  Racah Sum(Σ):  ", res.sum_cf)
end

function Base.:*(a::ExactResult, b::ExactResult)
    @assert a.k == b.k "Cannot multiply results from different levels k"
    return ExactResult(a.k, a.pref_sq * b.pref_sq, a.sum_cf * b.sum_cf)
end

function Base.:/(a::ExactResult, b::ExactResult)
    @assert a.k == b.k "Cannot divide results from different levels k"
    return ExactResult(a.k, a.pref_sq * inv(b.pref_sq), a.sum_cf * inv(b.sum_cf))
end


Base.://(a::ExactResult, b::ExactResult) = a / b

# Equality check (Algebraic)
function Base.:(==)(a::ExactResult, b::ExactResult)
    return a.k == b.k && a.pref_sq == b.pref_sq && a.sum_cf == b.sum_cf
end
