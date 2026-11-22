import DFL.ECProof

/-!
# Uniformity of (EC) over the volume (paper: Theorem `thm:EC`(ii))

For a window `Λ = [a, a + L)` of the chain `[0, N)` let `Q^Λ(x, y)` be the eigenfunction correlator of the
Bernoulli chain on `Λ` (charges `q|_Λ`), `0` unless `x, y ∈ Λ` (`corrQwin`). Under the large-deviation input:
`E sup_{Λ ∋ x, y} Q^Λ(x, y) ≤ C e^{-η|x-y|}`, with the supremum over all windows of `[0, N)`, uniformly in `N`
(`ec_sup_of_ldt`). Since the charges are i.i.d., this is the statement of the paper for all finite intervals
`Λ ⊂ ℤ` containing `x` and `y`; with windows `[0, L)` it gives `E sup_L Q^{[0,L)}(x, y) ≤ C e^{-η|x-y|}` for one
sequence of charges.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The charges of `q` on the window `[a, a + L)`, relabelled as `Fin L`. -/
def winCfg (q : Cfg N) (a L : ℕ) (h : a + L ≤ N) : Cfg L := fun i => q ⟨a + i.val, by omega⟩

/-- `Q^Λ(x, y)` for the window `Λ = [a, a + L)` (`0` unless `x, y ∈ Λ`). -/
def corrQwin (J h : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (x y : Fin N) : ℝ :=
  if hxy : a ≤ x.val ∧ x.val < a + L ∧ a ≤ y.val ∧ y.val < a + L then
    corrQ (hq J h (winCfg q a L hL)) (hq_isHermitian J h _) ⟨x.val - a, by omega⟩ ⟨y.val - a, by omega⟩
  else 0

/-- `sup_{Λ ∋ x, y} Q^Λ(x, y)`, the supremum over all windows `Λ = [a, a + L)` of `[0, N)`. -/
def corrQsup (J h : ℝ) (q : Cfg N) (x y : Fin N) : ℝ :=
  (univ : Finset (Fin (N + 1) × Fin (N + 1))).sup' univ_nonempty
    (fun p => if hp : p.1.val + p.2.val ≤ N then corrQwin J h q p.1.val p.2.val hp x y else 0)

/-- **(EC) uniformly over the volume** (paper: Theorem `thm:EC`(ii), eq. `eq:ECsup`):
`E sup_{Λ ∋ x, y} Q^Λ(x, y) ≤ C e^{-η|x-y|}`, uniformly in `N`. -/
theorem ec_sup_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q : Cfg N => corrQsup J h q x y) ≤ C * Real.exp (-η * Nat.dist x.val y.val) := by
  sorry

end DFL
