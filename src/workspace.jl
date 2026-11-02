"""
    EvaluationWorkspace()

Reusable Float64 scratch and fixed-q analytic tables. Pass `workspace=work` to `q6j`, `q3j`, `qcg`,
`fsymbol`, `gsymbol`, or `qeval(::FactorialSum)`. Reuse it sequentially, never
across concurrent tasks (batches give each worker its own). Tables are keyed by q, type and precision, with
at most four tiers; CG scratch keeps one sector.
"""
struct EvaluationWorkspace
    ratios::Vector{Float64}
    column::ColumnWork
    analytic::Dict{Any,Any}
    cg::Base.RefValue{Any}
end
EvaluationWorkspace(ratios::Vector{Float64}, column::ColumnWork, analytic::Dict{Any,Any}) =
    EvaluationWorkspace(ratios,column,analytic,Ref{Any}(nothing))
EvaluationWorkspace(ratios::Vector{Float64}, column::ColumnWork) =
    EvaluationWorkspace(ratios,column,Dict{Any,Any}())
EvaluationWorkspace() = EvaluationWorkspace(Float64[],ColumnWork())

_ratio_buffer(::Nothing, ::Type{T}, n::Int) where {T} = Vector{T}(undef,n)
function _ratio_buffer(w::EvaluationWorkspace, ::Type{Float64}, n::Int)
    length(w.ratios) < n && resize!(w.ratios,max(n,2length(w.ratios)))
    return w.ratios
end

# ---- scratch for calls that bring no workspace ----
#
# The `:lazy` policy stores each ratio's low part for the compensated fallback, and a column run needs three arrays as long
# as the column. Calls without a workspace borrow these per thread, each guarded by an atomic flag: a task finding one
# busy allocates, and one never returned (an exception) degrades that thread to allocating, never to sharing.

mutable struct _RatioSlot
    @atomic busy::Bool
    buf::Vector{Float64}
    @atomic column_busy::Bool
    column::ColumnWork{DWord}
end

const _RATIO_SLOTS = Ref{Vector{_RatioSlot}}(_RatioSlot[])
const _RATIO_SLOTS_LOCK = ReentrantLock()

function _ratio_slots()
    slots = _RATIO_SLOTS[]
    Threads.threadid() <= length(slots) && return slots
    @lock _RATIO_SLOTS_LOCK begin               # first use, or a thread added since: grow once
        slots = _RATIO_SLOTS[]
        n = max(Threads.maxthreadid(), Threads.threadid())
        if n > length(slots)
            fresh = vcat(slots, [_RatioSlot(false, Float64[], false, ColumnWork()) for _ in length(slots)+1:n])
            _RATIO_SLOTS[] = fresh
            slots = fresh
        end
    end
    return slots
end

"Borrow this thread's ratio buffer, or allocate one if it is busy. Return it with `_return_ratio_buffer`."
@inline function _borrow_ratio_buffer(n::Int)
    slots = _ratio_slots()
    tid = Threads.threadid()
    if tid <= length(slots)
        slot = @inbounds slots[tid]
        if (@atomicreplace slot.busy false => true).success
            length(slot.buf) < n && resize!(slot.buf, max(n, 2length(slot.buf)))
            return slot, slot.buf
        end
    end
    return nothing, Vector{Float64}(undef, n)
end
@inline _return_ratio_buffer(::Nothing) = nothing
@inline _return_ratio_buffer(slot::_RatioSlot) = (@atomic slot.busy = false; nothing)
"Column scratch for a call: its workspace's, or this thread's, or a new one if that is busy. Return it with `_return_column`."
@inline _borrow_column(w::EvaluationWorkspace) = (nothing, w.column)
@inline function _borrow_column(::Nothing)
    slots = _ratio_slots()
    tid = Threads.threadid()
    if tid <= length(slots)
        slot = @inbounds slots[tid]
        (@atomicreplace slot.column_busy false => true).success && return slot, slot.column
    end
    return nothing, ColumnWork()
end
@inline _return_column(::Nothing) = nothing
function _return_column(slot::_RatioSlot)
    # a very long column gives its memory back: allocating again is small next to the run that needs it
    length(slot.column.w) > COLUMN_KEEP && (slot.column = ColumnWork())
    @atomic slot.column_busy = false
    return nothing
end
const COLUMN_KEEP = 2^16
