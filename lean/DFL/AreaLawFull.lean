import DFL.AreaLawSharp
import DFL.CutLocal

/-!
# Area law for definite-charge initial states with the constants of the paper
(paper: Theorem `thm:arealaw`(b),(c))

For `Ψ₀ = |ξ⟩_σ ⊗ |n⟩` (links in `σˣ` eigenstates) the full matter-plus-gauge entanglement is the free-fermion
entanglement in the sector `q = q(n, ξ)` (`vonNeumann_rdmLGT_xprod`). For uniformly distributed `ξ`:

* `E_ξ S_full(t) = S_G(t) = E_q S(ρ^q_A(t))` for every `t` (`avg_full_eq_SG`);
* `sup_t S_G(t) ≤ 8 ln 2 M_2(F)` (`sup_SG_le`);
* `E_ξ sup_t S_full(t) ≤ 8 ln 2 M_2(F)` (`avg_sup_full_le`);
* `P_ξ(sup_t S_full(t) ≥ s) ≤ 8 ln 2 M_2(F)/s` for `s > 0` (`prob_sup_full_ge_le`).

Here `M_2(F) = ∑_{d ≥ 0} (d+1)² F(d)`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

namespace AreaLawFullAux

open AreaLawAux AreaLawSharpAux

/-- For fixed matter configuration `n`, the map `ξ ↦ chargeOf n ξ` is a bijection, with inverse
`linkOf n`. -/
def chargeEquiv (n : Cfg N) : Cfg N ≃ Cfg N :=
  ⟨chargeOf n, linkOf n, linkOf_chargeOf n, chargeOf_linkOf n⟩

/-- The Fock configuration of the occupied set of `n` is `n`. -/
lemma fockCfg_occupied (n : Cfg N) : fockCfg (univ.filter (fun x => n x = 1)) = n := by
  funext x
  unfold fockCfg
  rcases CutLocalAux.zmod2_cases (n x) with h0 | h0 <;> simp [h0]

/-- The basis vector `|n⟩` is the Fock vector of the occupied set of `n`. -/
lemma single_eq_fockVec (n : Cfg N) :
    (Pi.single n (1 : ℂ) : Cfg N → ℂ) = fockVec (univ.filter (fun x => n x = 1)) := by
  rw [fockVec, fockCfg_occupied]

/-- The entanglement entropy of a Slater determinant is nonnegative. -/
lemma vonNeumann_slater_nonneg (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    0 ≤ vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) := by
  rw [vonNeumann_rdm_slater A hA S t ℓ hℓ]
  exact Finset.sum_nonneg fun i _ =>
    Real.binEntropy_nonneg (eig_corrA_mem A hA S t ℓ hℓ i).1 (eig_corrA_mem A hA S t ℓ hℓ i).2

/-- Time-independent bound: `S(ρ^q_A(t)) ≤ 2 ln 2 ∑_{x<ℓ≤y} ∑_z Q(x,z) Q(y,z)`. -/
lemma vonNeumann_le_corrQ (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (t : ℝ) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) ≤
      2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z :=
  (vonNeumann_slater_sharp (hq J h q) (hq_isSymm J h q) S t ℓ hℓ).trans
    (mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      norm_corrMat_le_corrQ J h q S t x y) two_log_two_nonneg)

/-- The entanglement entropy of a sector is bounded in time. -/
lemma bddAbove_vonNeumann (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    BddAbove (Set.range fun t : ℝ =>
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) :=
  ⟨_, by
    rintro _ ⟨t, rfl⟩
    exact vonNeumann_le_corrQ J h q S ℓ hℓ t⟩

/-- `2 ln 2 Θ(F) ≤ 8 ln 2 M_2(F)`. -/
lemma two_log_two_thetaF_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    2 * Real.log 2 * thetaF F ≤ 8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  calc 2 * Real.log 2 * thetaF F ≤ 2 * Real.log 2 * (4 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d) :=
        mul_le_mul_of_nonneg_left (thetaF_le J h hEC) two_log_two_nonneg
    _ = 8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by ring

/-- **Markov's inequality** for the uniform average over a finite set: if `X ≥ 0` and
`c⁻¹ ∑ X ≤ M`, then `#{X ≥ s} / c ≤ M / s`. -/
lemma markov_card {ι : Type*} [Fintype ι] (X : ι → ℝ) (hX : ∀ i, 0 ≤ X i) (c M s : ℝ)
    (hc : 0 < c) (hs : 0 < s) (hM : (1 / c) * ∑ i, X i ≤ M) :
    ((univ.filter (fun i => s ≤ X i)).card : ℝ) / c ≤ M / s := by
  have hsum : ∑ i, X i ≤ M * c := by
    rw [one_div, inv_mul_le_iff₀ hc] at hM
    linarith
  have hmarkov : ((univ.filter (fun i => s ≤ X i)).card : ℝ) * s ≤ ∑ i, X i := by
    calc ((univ.filter (fun i => s ≤ X i)).card : ℝ) * s =
          ∑ i ∈ univ.filter (fun i => s ≤ X i), s := by
          rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ ∑ i ∈ univ.filter (fun i => s ≤ X i), X i :=
          Finset.sum_le_sum fun i hi => (Finset.mem_filter.1 hi).2
      _ ≤ ∑ i, X i :=
          Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun i _ _ => hX i
  rw [div_le_div_iff₀ hc hs]
  linarith

end AreaLawFullAux

open AreaLawAux AreaLawSharpAux AreaLawFullAux

/-- **`E_ξ S_full(t) = S_G(t)`** (paper: Theorem `thm:arealaw`(c)). -/
theorem avg_full_eq_SG (J h : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ) :
    (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      avgQ (fun q : Cfg N => vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))
        (rdm_isHermitian ℓ hℓ _)) := by
  unfold avgQ
  congr 1
  simp only [vonNeumann_rdmLGT_xprod]
  exact Equiv.sum_comp (chargeEquiv n) (fun q => vonNeumann (rdm ℓ hℓ
    (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1)) (rdm_isHermitian ℓ hℓ _))

/-- **`sup_t S_G(t) ≤ 8 ln 2 M_2(F)`** (paper: Theorem `thm:arealaw`(b)). -/
theorem sup_SG_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ) :
    avgQ (fun q : Cfg N => vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))
      (rdm_isHermitian ℓ hℓ _)) ≤ 8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  rw [single_eq_fockVec n]
  refine le_trans ?_ ((lgt_areaLaw_sharp J h hEC ℓ hℓ (univ.filter (fun x => n x = 1))).trans
    (two_log_two_thetaF_le hEC))
  exact avgQ_mono' fun q => le_ciSup (bddAbove_vonNeumann J h q _ ℓ hℓ) t

/-- **`E_ξ sup_t S_full(t) ≤ 8 ln 2 M_2(F)`** (paper: Theorem `thm:arealaw`(c)). -/
theorem avg_sup_full_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) :
    (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N, ⨆ t : ℝ,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) ≤
      8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  have h1 : ∑ ξ : Cfg N, ⨆ t : ℝ,
      vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      ∑ q : Cfg N, ⨆ t : ℝ, vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ
        fockVec (univ.filter (fun x => n x = 1)))) (rdm_isHermitian ℓ hℓ _) := by
    rw [← single_eq_fockVec n]
    simp only [vonNeumann_rdmLGT_xprod]
    exact Equiv.sum_comp (chargeEquiv n) (fun q => ⨆ t : ℝ, vonNeumann (rdm ℓ hℓ
      (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1)) (rdm_isHermitian ℓ hℓ _))
  rw [h1]
  exact (lgt_areaLaw_sharp J h hEC ℓ hℓ _).trans (two_log_two_thetaF_le hEC)

/-- **Markov bound** (paper: Theorem `thm:arealaw`(c)):
`P_ξ(sup_t S_full(t) ≥ s) ≤ 8 ln 2 M_2(F)/s`. -/
theorem prob_sup_full_ge_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N)
    (s : ℝ) (hs : 0 < s) :
    ((univ.filter (fun ξ : Cfg N => s ≤ ⨆ t : ℝ,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _))).card : ℝ) /
        2 ^ N ≤
      8 * Real.log 2 * (∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d) / s := by
  have hX0 : ∀ ξ : Cfg N, 0 ≤ ⨆ t : ℝ,
      vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) := by
    intro ξ
    refine Real.iSup_nonneg fun t => ?_
    rw [vonNeumann_rdmLGT_xprod, single_eq_fockVec n]
    exact vonNeumann_slater_nonneg _ (hq_isSymm J h _) _ t ℓ hℓ
  exact markov_card (fun ξ : Cfg N => ⨆ t : ℝ,
      vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _))
    hX0 (2 ^ N) _ s (by positivity) hs (avg_sup_full_le J h hEC ℓ hℓ n)

end DFL
