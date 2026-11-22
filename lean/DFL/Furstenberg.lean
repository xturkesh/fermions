import DFL.Transfer

/-!
# The hypotheses of Furstenberg's theorem for the Bernoulli–Anderson model
(paper: Proposition `prop:nocritical`(i),(ii),(iv))

For `λ ≠ 0` and every energy `E`, the two transfer matrices `T_E(±λ)` (the support of the law of
`T_E(λ q₀)`) generate a non-compact, strongly irreducible subgroup of `SL(2, ℝ)`, and they do not
commute (no critical energies in the sense of Jitomirskaya–Schulz-Baldes–Stolz). With Furstenberg's
theorem (not formalized) this gives the positivity of the Lyapunov exponent, which is part of the
literature input `DFL.LDT`.

* `transferT_mul_inv`: `P := T_E(λ) T_E(-λ)⁻¹ = [[1, -2λ], [0, 1]]`.
* `transferT_commutator`: `T_E(v) T_E(v') - T_E(v') T_E(v) = (v - v') [[0, 1], [1, 0]]`.
* `unipotent_pow_unbounded`: the powers of `P` are unbounded (non-compactness, (i)).
* `no_finite_invariant_lines`: no finite non-empty set of lines in `ℝ²` is invariant under both
  `T_E(λ)` and `T_E(-λ)` (strong irreducibility, (ii); invariance under the generators is what has to
  be excluded).
-/

open Matrix

noncomputable section

namespace DFL

/-- `T_E(λ) T_E(-λ)⁻¹ = [[1, -2λ], [0, 1]]` (with the inverse written as the adjugate). -/
theorem transferT_mul_inv (E lam : ℝ) :
    transferT E lam * (transferT E (-lam)).adjugate = !![1, -2 * lam; 0, 1] := by
  sorry

/-- The adjugate is the inverse: `T_E(v) adj(T_E(v)) = 1`. -/
theorem transferT_mul_adjugate (E v : ℝ) : transferT E v * (transferT E v).adjugate = 1 := by
  sorry

/-- **No critical energies** (paper: Proposition `prop:nocritical`(iv)):
`T_E(v) T_E(v') - T_E(v') T_E(v) = (v - v') [[0, 1], [1, 0]]`; for `v ≠ v'` the transfer matrices do
not commute. -/
theorem transferT_commutator (E v v' : ℝ) :
    transferT E v * transferT E v' - transferT E v' * transferT E v = (v - v') • !![0, 1; 1, 0] := by
  sorry

/-- The transfer matrices `T_E(λ)`, `T_E(-λ)` do not commute for `λ ≠ 0`. -/
theorem transferT_not_commute (E lam : ℝ) (hlam : lam ≠ 0) :
    transferT E lam * transferT E (-lam) ≠ transferT E (-lam) * transferT E lam := by
  sorry

/-- **Non-compactness** (paper: Proposition `prop:nocritical`(i)): the powers of
`P = [[1, -2λ], [0, 1]]` have unbounded entries. -/
theorem unipotent_pow_unbounded (lam : ℝ) (hlam : lam ≠ 0) (R : ℝ) :
    ∃ k : ℕ, R < |((!![1, -2 * lam; 0, 1] : Matrix (Fin 2) (Fin 2) ℝ) ^ k) 0 1| := by
  sorry

/-- **Strong irreducibility** (paper: Proposition `prop:nocritical`(ii)): a finite set of non-zero
vectors whose spanned lines are permuted by `T_E(λ)` and by `T_E(-λ)` is empty. Lines are encoded by
non-zero spanning vectors; "the image of the line of `u` is a line of the set" means that `T u` is a
non-zero multiple of some element. -/
theorem no_finite_invariant_lines (E lam : ℝ) (hlam : lam ≠ 0) (L : Finset (Fin 2 → ℝ))
    (hL0 : ∀ u ∈ L, u ≠ 0)
    (hinv : ∀ u ∈ L, ∀ v : ℝ, (v = lam ∨ v = -lam) →
      ∃ u' ∈ L, ∃ c : ℝ, c ≠ 0 ∧ transferT E v *ᵥ u = c • u') :
    L = ∅ := by
  sorry

end DFL
