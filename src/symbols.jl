# ---------------------------------------------------------------------------------
#  The symbol interface
#
#  Everything above the factorial rule — kernels, policy, precision tiers, exact projection, structural
#  queries, batching — already works for any rule and knows nothing about SU(2). What tied the package to
#  six specific symbols was not the mathematics but four parallel lookup tables, each keyed on something
#  different: `_rule_for` and `_admissible_at` on the public function, `_family` on the public function
#  again, and `_rule_family` on an internal rule-builder. Adding a symbol meant finding all four.
#
#  This file is the single table. A symbol is a singleton type with three required methods,
#
#      symbol_rule(sym, labels)        the factorial rule, from *undoubled* labels
#      nlabels(sym)                    how many labels it takes
#      level_admissible(sym, k, labels)
#
#  and two optional ones (`symbol_family` for recurrence routing, `symbol_of` to attach a public
#  function). The old helpers are kept, and now dispatch through here, so no call site changed and
#  behaviour is identical — `test/symbols.jl` asserts that against the tables this replaced.
#
#  To add a symbol: define the type, the three methods, and, if a public function should route to it,
#  one `symbol_of`. Nothing else in the package needs to know.
# ---------------------------------------------------------------------------------

"""
    QSymbol

A recoupling symbol represented by a type. Its factorial rule, admissibility checks and recurrence family
are selected by method dispatch.
"""
abstract type QSymbol end

struct SixJ <: QSymbol end
struct ThreeJ <: QSymbol end
struct FSymbol <: QSymbol end
struct GSymbol <: QSymbol end
struct Tetrahedron <: QSymbol end
struct ThetaValue <: QSymbol end

"Number of labels the symbol takes. The 3j accepts 5 or 6: `m₃` is determined by the other two."
nlabels(::SixJ) = 6
nlabels(::ThreeJ) = 6
nlabels(::FSymbol) = 6
nlabels(::GSymbol) = 6
nlabels(::Tetrahedron) = 6
nlabels(::ThetaValue) = 3

"""
    symbol_rule(sym, labels...) -> FactorialSum

The factorial rule of a symbol, from undoubled labels. This is the one method every other layer of the
package consumes, and the only one a new symbol must get right.
"""
symbol_rule(::SixJ, js...) = sixj_sum(doubled(js...)...)
symbol_rule(::FSymbol, js...) = fsymbol_sum(doubled(js...)...)
symbol_rule(::GSymbol, js...) = gsymbol_sum(doubled(js...)...)
symbol_rule(::Tetrahedron, js...) = tetrahedron_sum(doubled(js...)...)
symbol_rule(::ThreeJ, j1, j2, j3, m1, m2, m3 = -m1 - m2) =
    threej_sum(doubled(j1, j2, j3, m1, m2, m3)...)

"""
    level_admissible(sym, k, labels...) -> Bool

Whether the labels carry a representation at level `k`. Labels outside the admissible set are zero by
convention rather than singular, and every layer follows that convention.
"""
level_admissible(::SixJ, k::Int, js...) = _qδtet(doubled(js...)..., k)
level_admissible(::FSymbol, k::Int, js...) = _qδtet(doubled(js...)..., k)
level_admissible(::GSymbol, k::Int, js...) = _qδtet(doubled(js...)..., k)
level_admissible(::Tetrahedron, k::Int, js...) = _qδtet(doubled(js...)..., k)
level_admissible(::ThreeJ, k::Int, j1, j2, j3, m1, m2, m3 = -m1 - m2) =
    _qδ(doubled(j1, j2, j3)..., k)

"""
    symbol_family(sym) -> Val or nothing

Which recurrence family serves the symbol, for the column and near-edge routing. `nothing` means the
symbol has no recurrence and is always evaluated from its rule — a perfectly good answer, and the right
default for a new symbol.
"""
symbol_family(::SixJ) = Val(:sixj)
symbol_family(::FSymbol) = Val(:f)
symbol_family(::GSymbol) = Val(:g)
symbol_family(::ThreeJ) = Val(:threej)
symbol_family(::QSymbol) = nothing

"""
    symbol_of(f) -> QSymbol or nothing

The symbol a public evaluation function computes. Attaching a new function to a symbol is one method.
"""
symbol_of(::typeof(q6j)) = SixJ()
symbol_of(::typeof(q3j_factorial)) = ThreeJ()
symbol_of(::typeof(fsymbol)) = FSymbol()
symbol_of(::typeof(gsymbol)) = GSymbol()
symbol_of(::Any) = nothing
