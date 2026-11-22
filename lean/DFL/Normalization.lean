import DFL.Imbalance

/-!
# Normalization of the one-body operators (paper: Lemma `lem:normalization`)

With `Γ = diag((-1)^x)` (`stag`):

* `Γ h_q Γ⁻¹ = J (Δ + λ diag(q))`, `λ = 2h/J` (`stag_conj_hq`, `stag_conj_hq_eq_smul`): the sign of
  the hopping is irrelevant;
* consequently the eigenfunction correlator is the same for `J` and `-J` (`corrQ_hq_neg`);
* the joint law of `(Q(x, y))_{x,y}` is the same for `λ` and `-λ` (`avgQ_law_neg_h`), because
  `h_q` at `-h` is `h_{-q}` at `h` (`hq_neg_h`) and `-q` has the same law as `q`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Γ² = 1`. -/
theorem stag_mul_stag : stag N * stag N = 1 := by
  sorry

/-- `Γ h_q Γ = h_q` with the opposite hopping sign. -/
theorem stag_conj_hq (J h : ℝ) (q : Cfg N) : stag N * hq J h q * stag N = hq (-J) h q := by
  sorry

/-- **Normalization** (paper: Lemma `lem:normalization`): `Γ h_q Γ⁻¹ = J (Δ + λ diag(q))`,
`λ = 2h/J`. -/
theorem stag_conj_hq_eq_smul (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    stag N * hq J h q * stag N =
      J • (pathAdj N + (2 * h / J) • diagonal (fun x => sgn (q x))) := by
  sorry

/-- The eigenfunction correlator does not depend on the sign of the hopping. -/
theorem corrQ_hq_neg (J h : ℝ) (q : Cfg N) (x y : Fin N) :
    corrQ (hq (-J) h q) (hq_isHermitian (-J) h q) x y = corrQ (hq J h q) (hq_isHermitian J h q) x y := by
  sorry

/-- Flipping the sign of `h` is the same as flipping all charges. -/
theorem hq_neg_h (J h : ℝ) (q : Cfg N) : hq J (-h) q = hq J h (fun x => q x + 1) := by
  sorry

/-- **Symmetry `λ ↦ -λ`** (paper: Lemma `lem:normalization`): the joint law of the correlators
`(Q(x, y))_{x,y}` under i.i.d. uniform charges is the same for `h` and `-h`. -/
theorem avgQ_law_neg_h (J h : ℝ) (f : (Fin N → Fin N → ℝ) → ℝ) :
    avgQ (fun q : Cfg N => f (fun x y => corrQ (hq J (-h) q) (hq_isHermitian J (-h) q) x y)) =
      avgQ (fun q : Cfg N => f (fun x y => corrQ (hq J h q) (hq_isHermitian J h q) x y)) := by
  sorry

end DFL
