module QRecouplingTensorKitExt

# SU(2)_k and real-q SU(2) as TensorKit sectors, with F- and R-symbols from the QRecoupling kernels. It needs only
# TensorKitSectors, so `using TensorKit` loads it, and none of it touches the exact layer or Nemo.
#
# The F-symbols are real and orthogonal at every level and real q, so their scalar type is Float64 and planar tensors can
# stay real. Only the R-symbols are complex. Labels are doubled spins, as everywhere in the kernels.

using PrecompileTools: @compile_workload
using QRecoupling
using QRecoupling: fsymbol_sum, _symbol_value, _qδ, _qδtet
import QRecoupling: SU2k, SU2q

using TensorKitSectors: Sector, SectorValues, SectorProductIterator, SimpleFusion, Anyonic, NoBraiding
import TensorKitSectors: Nsymbol, Fsymbol, Rsymbol, dim, dual, unit, twist, findindex, FusionStyle, BraidingStyle,
    fusionscalartype, braidingscalartype, sectorscalartype, type_repr

abstract type DeformedSU2 <: Sector end

"""
A spin of SU(2) at level `k`: `SU2k(k)(j)`, with `2j` integer and `0 ≤ j ≤ k/2`, at q = e^{iπ/(k+2)}. The field `J` is `2j`.
"""
struct SU2kIrrep{k} <: DeformedSU2
    J::Int
    function SU2kIrrep{k}(j::Real) where {k}
        J = _doubled(j)
        J <= k || throw(DomainError(j, "a spin at level $k is at most $(k // 2)"))
        return new{k}(J)
    end
end

"A spin of SU(2) at real q: `SU2q(q)(j)`, with `2j` a nonnegative integer. The field `J` is `2j`."
struct SU2qIrrep{q} <: DeformedSU2
    J::Int
    SU2qIrrep{q}(j::Real) where {q} = new{q}(_doubled(j))
end

function _doubled(j::Real)
    J = 2j
    isinteger(J) && J >= 0 || throw(DomainError(j, "a spin is a nonnegative multiple of 1/2"))
    return Int(J)
end

function SU2k(k::Integer)
    k >= 0 || throw(DomainError(k, "level must be nonnegative"))
    return SU2kIrrep{Int(k)}
end
function SU2q(q::Real)
    isfinite(q) && q > 0 || throw(DomainError(q, "q must be positive and finite"))
    return SU2qIrrep{Float64(q)}
end

# ---- labels ----

Base.convert(::Type{I}, j::Real) where {I<:DeformedSU2} = I(j)
Base.hash(a::DeformedSU2, h::UInt) = hash(a.J, h)
Base.isless(a::I, b::I) where {I<:DeformedSU2} = isless(a.J, b.J)

type_repr(::Type{SU2kIrrep{k}}) where {k} = "SU2k($k)"
type_repr(::Type{SU2qIrrep{q}}) where {q} = "SU2q($q)"
Base.show(io::IO, ::Type{I}) where {I<:Union{SU2kIrrep{k} where {k},SU2qIrrep{q} where {q}}} =
    isconcretetype(I) ? print(io, type_repr(I)) : Base.show_datatype(io, I)
function Base.show(io::IO, a::I) where {I<:DeformedSU2}
    j = isodd(a.J) ? string(a.J, "//2") : string(a.J ÷ 2)
    get(io, :typeinfo, nothing) === I ? print(io, j) : print(io, type_repr(I), "(", j, ")")
end

unit(::Type{I}) where {I<:DeformedSU2} = I(0)
dual(a::DeformedSU2) = a

# A graded space stores one integer per value of a finite sector in its type, which only suits a small sector.
const SMALL_LEVEL = 32
Base.IteratorSize(::Type{SectorValues{SU2kIrrep{k}}}) where {k} = k <= SMALL_LEVEL ? Base.HasLength() : Base.SizeUnknown()
Base.IteratorSize(::Type{<:SectorValues{<:SU2qIrrep}}) = Base.IsInfinite()
Base.length(::SectorValues{SU2kIrrep{k}}) where {k} = k + 1
Base.iterate(::SectorValues{SU2kIrrep{k}}, J::Int = 0) where {k} = J > k ? nothing : (SU2kIrrep{k}(J // 2), J + 1)
Base.iterate(::SectorValues{I}, J::Int = 0) where {I<:SU2qIrrep} = (I(J // 2), J + 1)
Base.@propagate_inbounds function Base.getindex(::SectorValues{SU2kIrrep{k}}, i::Int) where {k}
    @boundscheck 1 <= i <= k + 1 || throw(BoundsError(values(SU2kIrrep{k}), i))
    return SU2kIrrep{k}((i - 1) // 2)
end
Base.@propagate_inbounds function Base.getindex(::SectorValues{I}, i::Int) where {I<:SU2qIrrep}
    @boundscheck 1 <= i || throw(BoundsError(values(I), i))
    return I((i - 1) // 2)
end
findindex(::SectorValues{I}, a::I) where {I<:DeformedSU2} = a.J + 1

# ---- fusion ----

FusionStyle(::Type{<:DeformedSU2}) = SimpleFusion()
BraidingStyle(::Type{<:SU2kIrrep}) = Anyonic()
BraidingStyle(::Type{<:SU2qIrrep}) = NoBraiding()      # the R-symbols of real q are not unitary, which TensorKit assumes
fusionscalartype(::Type{<:DeformedSU2}) = Float64
braidingscalartype(::Type{<:SU2kIrrep}) = ComplexF64
sectorscalartype(::Type{<:SU2kIrrep}) = ComplexF64
sectorscalartype(::Type{<:SU2qIrrep}) = Float64

_top(a::SU2kIrrep{k}, b::SU2kIrrep{k}) where {k} = min(a.J + b.J, 2k - a.J - b.J)
_top(a::SU2qIrrep{q}, b::SU2qIrrep{q}) where {q} = a.J + b.J

Base.IteratorSize(::Type{<:SectorProductIterator{<:DeformedSU2}}) = Base.HasLength()
Base.length(it::SectorProductIterator{<:DeformedSU2}) = (_top(it.a, it.b) - abs(it.a.J - it.b.J)) ÷ 2 + 1
function Base.iterate(it::SectorProductIterator{I}, J::Int = abs(it.a.J - it.b.J)) where {I<:DeformedSU2}
    return J > _top(it.a, it.b) ? nothing : (I(J // 2), J + 2)
end

Nsymbol(a::SU2kIrrep{k}, b::SU2kIrrep{k}, c::SU2kIrrep{k}) where {k} = _qδ(a.J, b.J, c.J, k)
Nsymbol(a::I, b::I, c::I) where {I<:SU2qIrrep} = abs(a.J - b.J) <= c.J <= a.J + b.J && iseven(a.J + b.J + c.J)

_dim(J::Int, k::Int) = qdim(J // 2; k = k)::Float64
_dim(J::Int, q::Float64) = qdim(J // 2; q = q)::Float64
dim(a::SU2kIrrep{k}) where {k} = _dim(a.J, k)
dim(a::SU2qIrrep{q}) where {q} = _dim(a.J, q)

# ---- topological data ----
#
# TensorKit's F^{abc}_{d,ef} is (−1)^{a+b+c+d} √([2e+1][2f+1]) {a b e; c d f}, which is `fsymbol(a, b, e, c, d, f)`.

# Everything here runs through one function per kind of sector, so a new level compiles only thin wrappers.
function _fsymbol(Js::NTuple{6,Int}, k::Int)
    _qδtet(Js..., k) || return 0.0
    return _symbol_value(fsymbol_sum(Js...), () -> true, k, nothing, false, Float64; family = Val(:f), labels = Js)::Float64
end
function _fsymbol(Js::NTuple{6,Int}, q::Float64)
    return _symbol_value(fsymbol_sum(Js...), () -> true, nothing, q, false, Float64; family = Val(:f), labels = Js)::Float64
end

Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {k,I<:SU2kIrrep{k}} = _fsymbol((a.J, b.J, e.J, c.J, d.J, f.J), k)
function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {q,I<:SU2qIrrep{q}}
    (Nsymbol(a, b, e) && Nsymbol(e, c, d) && Nsymbol(b, c, f) && Nsymbol(a, f, d)) || return 0.0
    return _fsymbol((a.J, b.J, e.J, c.J, d.J, f.J), q)
end

_rsymbol(A::Int, B::Int, C::Int, k::Int) = rmatrix(A // 2, B // 2, C // 2; k = k)::ComplexF64
_twist(J::Int, k::Int) = QRecoupling.twist(J // 2; k = k)::ComplexF64
Rsymbol(a::SU2kIrrep{k}, b::SU2kIrrep{k}, c::SU2kIrrep{k}) where {k} = _rsymbol(a.J, b.J, c.J, k)
twist(a::SU2kIrrep{k}) where {k} = _twist(a.J, k)

@compile_workload begin
    Js = (2, 2, 2, 2, 2, 2)
    _fsymbol(Js, 10); _fsymbol(Js, 0.8); _dim(2, 10); _dim(2, 0.8); _rsymbol(2, 2, 2, 10); _twist(2, 10)
    QRecoupling.empty_caches!()
end

end
