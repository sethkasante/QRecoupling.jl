"""
    EvaluationWorkspace()

Reusable Float64 ratio and recoupling-recurrence scratch. Pass `workspace=work` to
`q6j`, `q3j`, `fsymbol`, `gsymbol`, or `qeval(::FactorialSum)`. Storage grows on demand.
One workspace may be reused sequentially, but must not be shared by concurrent tasks.
Batched symbol evaluation creates a separate workspace for each worker automatically.
"""
struct EvaluationWorkspace
    ratios::Vector{Float64}
    column::ColumnWork
end
EvaluationWorkspace() = EvaluationWorkspace(Float64[],ColumnWork())

_ratio_buffer(::Nothing, ::Type{T}, n::Int) where {T} = Vector{T}(undef,n)
function _ratio_buffer(w::EvaluationWorkspace, ::Type{Float64}, n::Int)
    length(w.ratios) < n && resize!(w.ratios,max(n,2length(w.ratios)))
    return w.ratios
end
_column_workspace(::Nothing) = ColumnWork()
_column_workspace(w::EvaluationWorkspace) = w.column
