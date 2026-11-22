import DFL.GraphFlux

/-!
# Two dimensions: Wilson loops and the counting of flux sectors
(paper: Theorem `thm:2dmap`(a), Lemma `lem:2dgauge`(i)–(iii))

* `wilsonLoop z = ∏_{e ∈ z} σᶻ_e` for a cycle `z ∈ ker ∂` (in the `σˣ` basis: it flips the labels of the
  edges of `z`). Wilson loops are involutions, commute with each other, with every Gauss operator and
  with the Hamiltonian; `∏_i G_i = (-1)^{N_f}` (paper: Theorem `thm:2dmap`(a), Eq. `eq:2dparity`).
* On a connected graph, a gauge transformation fixes a sign configuration iff it is constant
  (uniqueness up to a global sign, Lemma `lem:2dgauge`(i)); every gauge orbit has `2^{|V|-1}` elements
  (Lemma `lem:2dgauge`(iii)), and there are `2^β = 2^{|E| - |V| + 1}` orbits, i.e. flux sectors
  (Lemma `lem:2dgauge`(ii) with `flux_eq_iff_gauge`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The Wilson loop of `z ∈ Cfg M`: flips the labels of the edges in `z`. -/
def wilsonLoop (_G : FinGraph N M) (z : Cfg M) : Matrix (Cfg N × Cfg M) (Cfg N × Cfg M) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ (fun x' x => if x' = x + z then 1 else 0)

/-- Wilson loops are involutions. -/
theorem wilsonLoop_mul_self (z : Cfg M) : G.wilsonLoop z * G.wilsonLoop z = 1 := by
  sorry

/-- Wilson loops commute with each other. -/
theorem wilsonLoop_comm (z z' : Cfg M) :
    G.wilsonLoop z * G.wilsonLoop z' = G.wilsonLoop z' * G.wilsonLoop z := by
  sorry

/-- Wilson loops of cycles commute with the Gauss operators. -/
theorem wilsonLoop_comm_gaussOp (z : Cfg M) (hz : G.bdry z = 0) (i : Fin N) :
    G.gaussOp i * G.wilsonLoop z = G.wilsonLoop z * G.gaussOp i := by
  sorry

/-- Wilson loops of cycles commute with the Hamiltonian (paper: Theorem `thm:2dmap`(a)). -/
theorem wilsonLoop_comm_ham (J h : ℝ) (z : Cfg M) (hz : G.bdry z = 0) :
    G.ham J h * G.wilsonLoop z = G.wilsonLoop z * G.ham J h := by
  sorry

/-- `∏_i G_i = (-1)^{N_f}` (paper: Eq. `eq:2dparity`); the Gauss operators are diagonal and this is the
product of their entries. -/
theorem gaussOp_prod (p : Cfg N × Cfg M) : ∏ i, G.gaussOp i p p = sgn (∑ i, p.1 i) := by
  sorry

/-- Adjacency of vertices. -/
def Adj (i j : Fin N) : Prop := ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)

/-- The graph is connected. -/
def Connected : Prop := ∀ i j : Fin N, Relation.ReflTransGen G.Adj i j

/-- **Uniqueness of the gauge transformation up to a global sign** (paper: Lemma `lem:2dgauge`(i)): on
a connected graph, `θ` fixes `s` iff `θ` is constant. -/
theorem gaugeAct_eq_self_iff (hG : G.Connected) (θ : Cfg N) (s : Cfg M) :
    G.gaugeAct θ s = s ↔ ∀ i j, θ i = θ j := by
  sorry

/-- **Gauge orbits** (paper: Lemma `lem:2dgauge`(iii)): on a connected graph with `N ≥ 1` vertices,
every gauge orbit has `2^{N-1}` elements. -/
theorem card_gaugeOrbit (hG : G.Connected) (hN : 0 < N) (s : Cfg M) :
    (univ.image (fun θ : Cfg N => G.gaugeAct θ s)).card = 2 ^ (N - 1) := by
  sorry

/-- **Number of flux sectors** (paper: Lemma `lem:2dgauge`(ii)): on a connected graph with `N ≥ 1`
vertices there are `2^β`, `β = M - N + 1`, gauge orbits of sign configurations. -/
theorem card_fluxClasses (hG : G.Connected) (hN : 0 < N) :
    (univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s))).card =
      2 ^ (M + 1 - N) := by
  sorry

end FinGraph

end DFL
