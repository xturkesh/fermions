import DFL.AnnealedInterval
import DFL.TransportSharp

/-!
# Number entropy of bulk intervals in finite volume (paper: proof of Theorem `thm:infvol`(e))

For an interval `B = [a, b]` of sites and a Fock state with occupied set `S`, `s_B = |S ∩ B|`:

* `⟨Ψ(t), (N_B - s_B)² Ψ(t)⟩ ≤ 2Θ(F) + 4Σ'₁(F)²` (`lgt_bulk_transport_sq`): the law of total variance,
  `Var_q(N_B) = ∑_{x ∈ B, y ∉ B} |K_{xy}|² ≤ ∑_{x ∈ B, y ∉ B} |K_{xy}|` and `|⟨N_B⟩_q - s_B| ≤ ∑_{x ∈ B, z ∉ B} Q(x, z)`,
  with Lemma `lem:annealed`(iii);
* `S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` (`lgt_numberEntropy_bulk`), uniformly in `N`, `t`, `B`, `S`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

namespace NumberEntropyBulkAux

open TransportSharpAux

/-- `|z|² ≤ |z|` whenever `|z| ≤ 1`. -/
lemma normSq_le_norm_of_norm_le_one {z : ℂ} (hz : ‖z‖ ≤ 1) : Complex.normSq z ≤ ‖z‖ := by
  rw [Complex.normSq_eq_norm_sq]
  have h0 := norm_nonneg z
  nlinarith

/-- `Var_q(N_B) = ∑_{x ∈ B, y ∉ B} |K_q(t)_{xy}|² ≤ sup_s ∑_{x ∈ B, y ∉ B} |K_q(s)_{xy}|`
(since `|K| ≤ 1`). -/
lemma sum_normSq_corrMat_le_iSup (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (a b : ℕ) (t : ℝ) :
    ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, Complex.normSq (corrMat (hq J h q) S t x y) ≤
      ⨆ s : ℝ, ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S s x y‖ := by
  have hbdd : BddAbove (Set.range fun s : ℝ =>
      ∑ x ∈ ivl (N := N) a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S s x y‖) := by
    refine ⟨∑ x ∈ ivl (N := N) a b, ∑ y ∈ univ \ ivl (N := N) a b, (1 : ℝ), ?_⟩
    rintro _ ⟨s, rfl⟩
    exact Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      norm_corrMat_le_one _ (hq_isSymm J h q) S s x y
  calc ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, Complex.normSq (corrMat (hq J h q) S t x y)
      ≤ ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S t x y‖ :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
          normSq_le_norm_of_norm_le_one (norm_corrMat_le_one _ (hq_isSymm J h q) S t x y)
    _ ≤ _ := le_ciSup hbdd t

/-- `(⟨N_B⟩_q - s_B)² ≤ (∑_{x ∈ B, z ∉ B} Q(x, z))²`, from `|⟨N_B⟩_q - s_B| ≤ ∑_{x ∈ B, z ∉ B} |U_{xz}|²`
and `|U|² ≤ Q`. -/
lemma mean_sub_sq_le (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (a b : ℕ) (t : ℝ) :
    (∑ x ∈ ivl a b, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
        ((S ∩ ivl a b).card : ℝ)) ^ 2 ≤
      (∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b, corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2 := by
  have hab := abs_mean_sub_le_symm (evo (hq J h q) t) (evo_unitary _ (hq_isSymm J h q) t)
    (evo_transpose _ (hq_isSymm J h q) t) S (ivl a b)
  have hle : ∑ x ∈ ivl (N := N) a b, ∑ z ∈ univ \ ivl a b,
      Complex.normSq (evo (hq J h q) t x z) ≤
      ∑ x ∈ ivl (N := N) a b, ∑ z ∈ univ \ ivl a b,
        corrQ (hq J h q) (hq_isHermitian J h q) x z :=
    Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun z _ => by
      rw [Complex.normSq_eq_norm_sq]
      exact norm_evo_sq_le_corrQ J h q t x z
  rw [← sq_abs]
  exact pow_le_pow_left₀ (abs_nonneg _) (hab.trans hle) 2

/-- **Averaged second moment of the charge of an interval**:
`E_q ∑_n (N_B(n) - s_B)² |ψ_q(t)(n)|² ≤ 2Θ(F) + 4Σ'₁(F)²`. -/
lemma avgQ_bulk_second_moment_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N))
    (a b : ℕ) (t : ℝ) :
    avgQ (fun q : Cfg N => ∑ n, (count (ivl a b) n - ((S ∩ ivl a b).card : ℝ)) ^ 2 *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) ≤
      2 * thetaF F + 4 * sigma1' F ^ 2 := by
  have e1 : (fun q : Cfg N => ∑ n, (count (ivl a b) n - ((S ∩ ivl a b).card : ℝ)) ^ 2 *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) =
      fun q : Cfg N => (∑ x ∈ ivl (N := N) a b, ∑ y ∈ univ \ ivl (N := N) a b,
        Complex.normSq (corrMat (hq J h q) S t x y)) +
        (∑ x ∈ ivl (N := N) a b, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
          ((S ∩ ivl a b).card : ℝ)) ^ 2 := by
    funext q
    exact sector_second_moment _ (hq_isSymm J h q) S _ t _
  rw [e1, avgQ_add]
  exact add_le_add
    ((avgQ_mono fun q => sum_normSq_corrMat_le_iSup J h q S a b t).trans
      (annealed_ivl_K J h hEC S a b))
    ((avgQ_mono fun q => mean_sub_sq_le J h q S a b t).trans (annealed_ivl_Q_sq J h hEC a b))

end NumberEntropyBulkAux

open TransportSharpAux NumberEntropyBulkAux

/-- **Second moment of the charge in a bulk interval**: `⟨Ψ(t), (N_B - s_B)² Ψ(t)⟩ ≤ 2Θ(F) + 4Σ'₁(F)²`. -/
theorem lgt_bulk_transport_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ)
    (t : ℝ) :
    (expect (matter ((numOp (ivl a b) - ((S ∩ ivl a b).card : ℝ) • 1) *
        (numOp (ivl a b) - ((S ∩ ivl a b).card : ℝ) • 1))) (lgtState J h (fockVec S) t)).re ≤
      2 * thetaF F + 4 * sigma1' F ^ 2 := by
  rw [lgt_second_moment_eq]
  exact avgQ_bulk_second_moment_le hEC S a b t

/-- **Bounded number entropy of a bulk interval**:
`S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))`. -/
theorem lgt_numberEntropy_bulk (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ)
    (t : ℝ) :
    shannon N (lgtNumDist J h (ivl a b) (fockVec S) t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (2 * thetaF F + 4 * sigma1' F ^ 2 + 1 / 12)) := by
  have hp := lgtNumDist_isProb J h (ivl a b) S t
  refine shannon_le_of_pvar_le hp ?_
  calc pvar N (lgtNumDist J h (ivl a b) (fockVec S) t)
      ≤ ∑ k ∈ range (N + 1), lgtNumDist J h (ivl a b) (fockVec S) t k *
          ((k : ℝ) - ((S ∩ ivl a b).card : ℝ)) ^ 2 := pvar_le_sum_sq hp _
    _ = avgQ (fun q => ∑ k ∈ range (N + 1),
          numDist (ivl a b) (evo (dGamma (hq J h q)) t *ᵥ fockVec S) k *
            ((k : ℝ) - ((S ∩ ivl a b).card : ℝ)) ^ 2) := by
        rw [avgQ_sum]
        refine Finset.sum_congr rfl (fun k _ => ?_)
        rw [lgtNumDist_eq_avgQ, mul_comm, ← avgQ_const_mul]
        congr 1
        funext q
        ring
    _ = avgQ (fun q : Cfg N => ∑ n, (count (ivl a b) n - ((S ∩ ivl a b).card : ℝ)) ^ 2 *
          Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) := by
        congr 1
        funext q
        exact sum_numDist_mul_sq _ _ _
    _ ≤ 2 * thetaF F + 4 * sigma1' F ^ 2 := avgQ_bulk_second_moment_le hEC S a b t

end DFL
