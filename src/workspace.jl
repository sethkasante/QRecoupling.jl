"""
    EvaluationWorkspace()

Reusable Float64 ratio/recoupling scratch and bounded fixed-q analytic tables. Pass `workspace=work` to
`q6j`, `q3j_factorial`, `fsymbol`, `gsymbol`, or `qeval(::FactorialSum)`. Storage grows on demand.
One workspace may be reused sequentially, but must not be shared by concurrent tasks.
Batched symbol evaluation creates a separate workspace for each worker automatically.
Analytic tables are keyed by q, numeric type, precision, and capacity; changing q
invalidates them. A workspace retains at most four precision/type tiers.
"""
struct EvaluationWorkspace
    ratios::Vector{Float64}
    column::ColumnWork
    analytic::Dict{Any,Any}
end
EvaluationWorkspace(ratios::Vector{Float64}, column::ColumnWork) =
    EvaluationWorkspace(ratios,column,Dict{Any,Any}())
EvaluationWorkspace() = EvaluationWorkspace(Float64[],ColumnWork())

_ratio_buffer(::Nothing, ::Type{T}, n::Int) where {T} = Vector{T}(undef,n)
function _ratio_buffer(w::EvaluationWorkspace, ::Type{Float64}, n::Int)
    length(w.ratios) < n && resize!(w.ratios,max(n,2length(w.ratios)))
    return w.ratios
end

# ---- a ratio buffer for calls that bring no workspace ----
#
# The `:lazy` policy stores each ratio's low part during the plain pass so that the compensated fallback
# reuses it, which makes the fallback 13–24% cheaper (`q6j(40×6; k = 200)`: 501 ns against 572 ns). Without a
# caller workspace that buffer was a fresh `Vector` per call: two heap allocations on every symbol with two
# or more ratio steps, although almost none of them ever fall back. So a call without a workspace borrows
# one buffer per thread instead.
#
# Borrowing is guarded by an atomic flag, not by the thread id alone. The region that holds the buffer is
# pure arithmetic and never yields, but a guard that relied on that would break silently the day someone
# added a yield (a log message is enough); with the flag a second task that finds the slot busy — because
# the first yielded, migrated, or is re-entering — simply allocates, which is the old behaviour. A slot that
# is never returned (an exception inside the region) degrades that thread to allocating, never to sharing.

mutable struct _RatioSlot
    @atomic busy::Bool
    buf::Vector{Float64}
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
            fresh = vcat(slots, [_RatioSlot(false, Float64[]) for _ in length(slots)+1:n])
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
_column_workspace(::Nothing) = ColumnWork()
_column_workspace(w::EvaluationWorkspace) = w.column
