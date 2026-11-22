import DFL.CutLocal

/-!
# Von Neumann entropy: Schur bound, concavity, mixing bound, two-block pinching
(paper: Lemma `lem:pinching`, and the entropy inequalities used in Proposition `prop:entropies`)

`vonNeumann ρ hρ = ∑_i η(λ_i)`, with `η(x) = -x log x` (`Real.negMulLog`) and `λ_i` the eigenvalues
of the Hermitian matrix `ρ` (defined in `DFL.AreaLaw`).

* `vonNeumann_le_diag`: `S(ρ) ≤ ∑_a η(ρ_aa)` for `ρ ≥ 0`. The diagonal of `ρ` is a doubly
  stochastic average of its eigenvalues, `ρ_aa = ∑_i |U_ai|² λ_i`, and `η` is concave.
* `vonNeumann_diagonal`: the entropy of a real diagonal matrix is the Shannon entropy of its diagonal.
* `vonNeumann_le_log_card`: `S(ρ) ≤ log dim` for density matrices.
* `vonNeumann_mul_conjTranspose_comm`: `S(A Aᴴ) = S(Aᴴ A)` (equal nonzero spectra).
* `vonNeumann_concave`: concavity, `∑ w_i S(σ_i) ≤ S(∑ w_i σ_i)`.
* `vonNeumann_mix_le`: the mixing bound `S(∑ w_i σ_i) ≤ ∑ w_i S(σ_i) + H(w)` for density matrices.
* `vonNeumann_le_pinching`, `vonNeumann_pinching_le`: the two-block pinching lemma
  `S(σ) ≤ S(PσP + P'σP') ≤ S(σ) + ln 2`, `P' = 1 - P`.
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

namespace VNEntropyAux

/-- Jensen's inequality for `η = negMulLog` and a doubly stochastic matrix `M`:
`∑_j η(x_j) ≤ ∑_a η(∑_j M_{aj} x_j)`. -/
lemma sum_negMulLog_le_stoch {α β : Type*} [Fintype α] [Fintype β] (M : α → β → ℝ)
    (x : β → ℝ) (hM0 : ∀ a j, 0 ≤ M a j) (hrow : ∀ a, ∑ j, M a j = 1)
    (hcol : ∀ j, ∑ a, M a j = 1) (hx : ∀ j, 0 ≤ x j) :
    ∑ j, Real.negMulLog (x j) ≤ ∑ a, Real.negMulLog (∑ j, M a j * x j) := by
  have hJ : ∀ a, ∑ j, M a j * Real.negMulLog (x j) ≤ Real.negMulLog (∑ j, M a j * x j) := by
    intro a
    simpa only [smul_eq_mul] using Real.concaveOn_negMulLog.le_map_sum (t := univ) (w := M a)
      (p := x) (fun j _ => hM0 a j) (hrow a) (fun j _ => Set.mem_Ici.2 (hx j))
  calc ∑ j, Real.negMulLog (x j) = ∑ a, ∑ j, M a j * Real.negMulLog (x j) := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [← Finset.sum_mul, hcol j, one_mul]
    _ ≤ _ := Finset.sum_le_sum fun a _ => hJ a

/-- The Shannon entropy of a probability vector is at most `log` of the number of outcomes. -/
lemma sum_negMulLog_le_log_card {α : Type*} [Fintype α] (p : α → ℝ) (hp0 : ∀ a, 0 ≤ p a)
    (hp1 : ∑ a, p a = 1) : ∑ a, Real.negMulLog (p a) ≤ Real.log (Fintype.card α) := by
  have hn : (0 : ℝ) < Fintype.card α := by
    rcases isEmpty_or_nonempty α with h | h
    · rw [Finset.univ_eq_empty, Finset.sum_empty] at hp1
      exact absurd hp1 zero_ne_one
    · exact_mod_cast Fintype.card_pos
  have hJ := Real.concaveOn_negMulLog.le_map_sum (t := univ)
    (w := fun _ => (Fintype.card α : ℝ)⁻¹) (p := p) (fun _ _ => inv_nonneg.2 hn.le)
    (by rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_inv_cancel₀ hn.ne'])
    (fun a _ => Set.mem_Ici.2 (hp0 a))
  simp only [smul_eq_mul, ← Finset.mul_sum, hp1, mul_one] at hJ
  have hη : Real.negMulLog (Fintype.card α : ℝ)⁻¹ =
      (Fintype.card α : ℝ)⁻¹ * Real.log (Fintype.card α) := by
    rw [Real.negMulLog, Real.log_inv]; ring
  rw [hη] at hJ
  exact le_of_mul_le_mul_left hJ (inv_pos.2 hn)

/-- Spectral decomposition of a Hermitian matrix, entrywise:
`ρ_{ab} = ∑_j U_{aj} λ_j conj U_{bj}`. -/
lemma spectral_apply {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (a b : ι) :
    ρ a b = ∑ j, (hH.eigenvectorUnitary : Matrix ι ι ℂ) a j * (hH.eigenvalues j : ℂ) *
      star ((hH.eigenvectorUnitary : Matrix ι ι ℂ) b j) := by
  conv_lhs => rw [hH.spectral_theorem]
  rw [Unitary.conjStarAlgAut_apply, Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, Matrix.star_apply, Function.comp_apply]
  rfl

/-- The rows of the eigenvector unitary are normalized. -/
lemma row_normSq {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (a : ι) :
    ∑ j, Complex.normSq ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a j) = 1 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hH.eigenvectorUnitary) a) a
  rw [Matrix.mul_apply, Matrix.one_apply_eq] at h
  have h' : ((∑ j, Complex.normSq ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a j) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Unitary.coe_star, Matrix.star_apply, Complex.star_def, Complex.mul_conj]
  exact_mod_cast h'

/-- The columns of the eigenvector unitary are normalized. -/
lemma col_normSq {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (j : ι) :
    ∑ a, Complex.normSq ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a j) = 1 := by
  have h := congrFun (congrFun (Unitary.coe_star_mul_self hH.eigenvectorUnitary) j) j
  rw [Matrix.mul_apply, Matrix.one_apply_eq] at h
  have h' : ((∑ a, Complex.normSq ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a j) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Matrix.star_apply, Complex.star_def, mul_comm, Complex.mul_conj]
  exact_mod_cast h'

/-- The diagonal of a Hermitian matrix is a stochastic average of its eigenvalues:
`ρ_aa = ∑_j |U_aj|² λ_j`. -/
lemma diag_re_eq {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (a : ι) :
    (ρ a a).re = ∑ j, Complex.normSq ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a j) *
      hH.eigenvalues j := by
  rw [spectral_apply hH, Complex.re_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [mul_right_comm, Complex.star_def, Complex.mul_conj, ← Complex.ofReal_mul,
    Complex.ofReal_re]

end VNEntropyAux

/-- **Schur bound**: the von Neumann entropy of a positive semidefinite matrix is at most the
Shannon entropy of its diagonal. -/
theorem vonNeumann_le_diag (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) :
    vonNeumann ρ hρ.isHermitian ≤ ∑ a, Real.negMulLog (ρ a a).re := by
  unfold vonNeumann
  simp only [VNEntropyAux.diag_re_eq hρ.isHermitian]
  exact VNEntropyAux.sum_negMulLog_le_stoch _ _ (fun a j => Complex.normSq_nonneg _)
    (VNEntropyAux.row_normSq hρ.isHermitian) (VNEntropyAux.col_normSq hρ.isHermitian)
    (fun j => hρ.eigenvalues_nonneg j)

/-- The entropy of a real diagonal matrix is the Shannon entropy of its diagonal. -/
theorem vonNeumann_diagonal (d : ι → ℝ) (hd : (diagonal (fun a => (d a : ℂ))).IsHermitian) :
    vonNeumann (diagonal (fun a => (d a : ℂ))) hd = ∑ a, Real.negMulLog (d a) := by
  have h1 : (diagonal (fun a => (d a : ℂ))).charpoly.roots =
      Multiset.map (fun a => (d a : ℂ)) univ.val := by
    rw [Matrix.charpoly_diagonal, Polynomial.roots_prod]
    · simp
    · simp [Finset.prod_ne_zero_iff, Polynomial.X_sub_C_ne_zero]
  rw [hd.roots_charpoly_eq_eigenvalues] at h1
  have h2 := congrArg (fun s : Multiset ℂ => (s.map (fun z => Real.negMulLog z.re)).sum) h1
  simpa [vonNeumann, Multiset.map_map, Function.comp_def] using h2

/-- The entropy of a density matrix is at most `log` of the dimension. -/
theorem vonNeumann_le_log_card (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (htr : ρ.trace = 1) :
    vonNeumann ρ hρ.isHermitian ≤ Real.log (Fintype.card ι) := by
  refine (vonNeumann_le_diag ρ hρ).trans
    (VNEntropyAux.sum_negMulLog_le_log_card _ (fun a => ?_) ?_)
  · exact (Complex.nonneg_iff.1 (hρ.diag_nonneg (i := a))).1
  · have := congrArg Complex.re htr
    simpa [Matrix.trace, Complex.re_sum] using this

/-- `A Aᴴ` and `Aᴴ A` have the same nonzero eigenvalues, hence the same entropy. -/
theorem vonNeumann_mul_conjTranspose_comm {κ : Type*} [Fintype κ] [DecidableEq κ]
    (A : Matrix ι κ ℂ) :
    vonNeumann (A * Aᴴ) (isHermitian_mul_conjTranspose_self A) =
      vonNeumann (Aᴴ * A) (isHermitian_conjTranspose_mul_self A) :=
  CutLocalAux.vonNeumann_eq_of_charpoly _ _ _ _ _ _ (Matrix.charpoly_mul_comm' A Aᴴ)

namespace VNEntropyAux

/-- Schur bound in a rotated basis: `S(σ) ≤ ∑_j η((Wᴴ σ W)_jj)` for a unitary `W`. -/
lemma vonNeumann_le_diag_conj (σ : Matrix ι ι ℂ) (hσ : σ.PosSemidef) (W : Matrix ι ι ℂ)
    (hW : W * Wᴴ = 1) :
    vonNeumann σ hσ.isHermitian ≤ ∑ j, Real.negMulLog ((Wᴴ * σ * W) j j).re := by
  have hpsd : (Wᴴ * σ * W).PosSemidef := hσ.conjTranspose_mul_mul_same W
  have h1 := CutLocalAux.vonNeumann_conj_isometry Wᴴ (by rwa [conjTranspose_conjTranspose]) σ
    hσ.isHermitian (by rw [conjTranspose_conjTranspose]; exact hpsd.isHermitian)
  calc vonNeumann σ hσ.isHermitian = vonNeumann (Wᴴ * σ * W) hpsd.isHermitian := by
        rw [← h1]; exact CutLocalAux.vonNeumann_congr _ _ (by rw [conjTranspose_conjTranspose])
    _ ≤ _ := vonNeumann_le_diag _ hpsd

/-- The eigenvalues of a Hermitian `τ` are the diagonal entries of `Vᴴ τ V`, `V` the eigenvector
unitary. -/
lemma eigenvalues_eq_diag {τ : Matrix ι ι ℂ} (hτ : τ.IsHermitian) (j : ι) :
    hτ.eigenvalues j = ((((hτ.eigenvectorUnitary : Matrix ι ι ℂ))ᴴ * τ *
      (hτ.eigenvectorUnitary : Matrix ι ι ℂ)) j j).re := by
  have h := hτ.conjStarAlgAut_star_eigenvectorUnitary
  rw [Unitary.conjStarAlgAut_star_apply] at h
  rw [← Matrix.star_eq_conjTranspose, h, diagonal_apply_eq, Function.comp_apply]
  rfl

end VNEntropyAux

/-- **Concavity** of the von Neumann entropy (finite convex combinations of positive semidefinite
matrices). -/
theorem vonNeumann_concave {κ : Type*} [Fintype κ] (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (hw1 : ∑ i, w i = 1) (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).PosSemidef)
    (hmix : (∑ i, (w i : ℂ) • σ i).IsHermitian) :
    ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian ≤ vonNeumann (∑ i, (w i : ℂ) • σ i) hmix := by
  obtain ⟨V, hV⟩ : ∃ V : Matrix ι ι ℂ, V = (hmix.eigenvectorUnitary : Matrix ι ι ℂ) := ⟨_, rfl⟩
  have hVV : V * Vᴴ = 1 := by
    rw [hV, ← Matrix.star_eq_conjTranspose]; exact Unitary.coe_mul_star_self _
  -- `c i j = ⟨v_j, σ_i v_j⟩`
  obtain ⟨c, hc⟩ : ∃ c : κ → ι → ℝ, c = fun i j => ((Vᴴ * σ i * V) j j).re := ⟨_, rfl⟩
  have hc0 : ∀ i j, 0 ≤ c i j := fun i j => by
    rw [hc]; exact (Complex.nonneg_iff.1 ((hσ i).conjTranspose_mul_mul_same V).diag_nonneg).1
  have hq : ∀ j, hmix.eigenvalues j = ∑ i, w i * c i j := by
    intro j
    rw [VNEntropyAux.eigenvalues_eq_diag hmix j, ← hV, Finset.mul_sum, Finset.sum_mul,
      Matrix.sum_apply, Complex.re_sum, hc]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Matrix.mul_smul, Matrix.smul_mul, Matrix.smul_apply, smul_eq_mul, Complex.re_ofReal_mul]
  -- Jensen in each eigen-direction `v_j` of the mixture
  have hJ : ∀ j, ∑ i, w i * Real.negMulLog (c i j) ≤ Real.negMulLog (hmix.eigenvalues j) := by
    intro j
    rw [hq j]
    simpa only [smul_eq_mul] using Real.concaveOn_negMulLog.le_map_sum (t := univ) (w := w)
      (p := fun i => c i j) (fun i _ => hw0 i) hw1 (fun i _ => Set.mem_Ici.2 (hc0 i j))
  calc ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian
      ≤ ∑ i, w i * ∑ j, Real.negMulLog (c i j) := Finset.sum_le_sum fun i _ =>
        mul_le_mul_of_nonneg_left
          (by rw [hc]; exact VNEntropyAux.vonNeumann_le_diag_conj _ (hσ i) V hVV) (hw0 i)
    _ = ∑ j, ∑ i, w i * Real.negMulLog (c i j) := by
        simp only [Finset.mul_sum]; exact Finset.sum_comm
    _ ≤ ∑ j, Real.negMulLog (hmix.eigenvalues j) := Finset.sum_le_sum fun j _ => hJ j
    _ = vonNeumann _ hmix := rfl

namespace VNEntropyAux

variable {κ : Type*} [Fintype κ]

/-- The factor `A` of a mixture, `A Aᴴ = ∑_i w_i σ_i`: its columns are `√(w_i λ_{ik}) u_{ik}`,
with `u_{ik}` and `λ_{ik}` the eigenvectors and eigenvalues of `σ_i`. -/
def mixFactor (w : κ → ℝ) (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).IsHermitian) :
    Matrix ι (κ × ι) ℂ :=
  fun a p => (Real.sqrt (w p.1 * (hσ p.1).eigenvalues p.2) : ℂ) *
    ((hσ p.1).eigenvectorUnitary : Matrix ι ι ℂ) a p.2

/-- `A Aᴴ = ∑_i w_i σ_i` for the factor `A = mixFactor w σ`. -/
lemma mixFactor_mul_conjTranspose (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i) (σ : κ → Matrix ι ι ℂ)
    (hσ : ∀ i, (σ i).PosSemidef) :
    mixFactor w σ (fun i => (hσ i).isHermitian) *
      (mixFactor w σ (fun i => (hσ i).isHermitian))ᴴ = ∑ i, (w i : ℂ) • σ i := by
  ext a b
  rw [Matrix.mul_apply, Fintype.sum_prod_type, Matrix.sum_apply]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Matrix.smul_apply, smul_eq_mul, spectral_apply (hσ i).isHermitian a b, Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  have hs : (Real.sqrt (w i * (hσ i).isHermitian.eigenvalues k) : ℂ) *
      (Real.sqrt (w i * (hσ i).isHermitian.eigenvalues k) : ℂ) =
      (w i : ℂ) * ((hσ i).isHermitian.eigenvalues k : ℂ) := by
    rw [← Complex.ofReal_mul,
      Real.mul_self_sqrt (mul_nonneg (hw0 i) ((hσ i).eigenvalues_nonneg k)), Complex.ofReal_mul]
  simp only [mixFactor, conjTranspose_apply, star_mul', Complex.star_def, Complex.conj_ofReal]
  linear_combination (((hσ i).isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) a k *
    (starRingEnd ℂ) (((hσ i).isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) b k)) * hs

omit [Fintype κ] in
/-- The diagonal of `Aᴴ A` for the factor `A = mixFactor w σ`:
`(AᴴA)_{(i,k),(i,k)} = w_i λ_{ik}`. -/
lemma mixFactor_diag (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i) (σ : κ → Matrix ι ι ℂ)
    (hσ : ∀ i, (σ i).PosSemidef) (i : κ) (k : ι) :
    (((mixFactor w σ (fun i => (hσ i).isHermitian))ᴴ *
      mixFactor w σ (fun i => (hσ i).isHermitian)) (i, k) (i, k)).re =
      w i * (hσ i).isHermitian.eigenvalues k := by
  have hx : 0 ≤ w i * (hσ i).isHermitian.eigenvalues k :=
    mul_nonneg (hw0 i) ((hσ i).eigenvalues_nonneg k)
  have key : ∀ a, ((mixFactor w σ (fun i => (hσ i).isHermitian))ᴴ (i, k) a *
      mixFactor w σ (fun i => (hσ i).isHermitian) a (i, k)).re =
      w i * (hσ i).isHermitian.eigenvalues k *
        Complex.normSq (((hσ i).isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) a k) := by
    intro a
    simp only [mixFactor, conjTranspose_apply, star_mul', Complex.star_def, Complex.conj_ofReal]
    rw [show ∀ (s : ℝ) (u : ℂ), (s : ℂ) * (starRingEnd ℂ) u * ((s : ℂ) * u) =
        ((s * s * Complex.normSq u : ℝ) : ℂ) from fun s u => by
      rw [Complex.ofReal_mul, Complex.ofReal_mul, Complex.normSq_eq_conj_mul_self]; ring]
    rw [Complex.ofReal_re, Real.mul_self_sqrt hx]
  rw [Matrix.mul_apply, Complex.re_sum, Finset.sum_congr rfl fun a _ => key a, ← Finset.mul_sum,
    col_normSq, mul_one]

/-- The mixing bound, for a mixture indexed by a type with decidable equality. -/
lemma vonNeumann_mix_le_aux [DecidableEq κ] (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).PosSemidef) (htr : ∀ i, (σ i).trace = 1)
    (hmix : (∑ i, (w i : ℂ) • σ i).IsHermitian) :
    vonNeumann (∑ i, (w i : ℂ) • σ i) hmix ≤
      ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian + ∑ i, Real.negMulLog (w i) := by
  have h1 : vonNeumann (∑ i, (w i : ℂ) • σ i) hmix =
      vonNeumann ((mixFactor w σ (fun i => (hσ i).isHermitian))ᴴ *
        mixFactor w σ (fun i => (hσ i).isHermitian))
        (isHermitian_conjTranspose_mul_self _) := by
    rw [← vonNeumann_mul_conjTranspose_comm]
    exact CutLocalAux.vonNeumann_congr _ _ (mixFactor_mul_conjTranspose w hw0 σ hσ).symm
  have hlam1 : ∀ i, ∑ k, (hσ i).isHermitian.eigenvalues k = 1 := by
    intro i
    have h := (hσ i).isHermitian.trace_eq_sum_eigenvalues
    rw [htr i] at h
    have h' : ((∑ k, (hσ i).isHermitian.eigenvalues k : ℝ) : ℂ) = 1 := by
      rw [Complex.ofReal_sum]; exact h.symm
    exact_mod_cast h'
  rw [h1]
  refine (vonNeumann_le_diag _ (posSemidef_conjTranspose_mul_self _)).trans (le_of_eq ?_)
  rw [Fintype.sum_prod_type]
  simp only [mixFactor_diag w hw0 σ hσ, Real.negMulLog_mul]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, hlam1 i, one_mul, ← Finset.mul_sum]
  unfold vonNeumann
  ring

end VNEntropyAux

/-- **Mixing bound**: `S(∑ w_i σ_i) ≤ ∑ w_i S(σ_i) + H(w)` for density matrices `σ_i` and weights
`w_i ≥ 0`. -/
theorem vonNeumann_mix_le {κ : Type*} [Fintype κ] (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).PosSemidef) (htr : ∀ i, (σ i).trace = 1)
    (hmix : (∑ i, (w i : ℂ) • σ i).IsHermitian) :
    vonNeumann (∑ i, (w i : ℂ) • σ i) hmix ≤
      ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian + ∑ i, Real.negMulLog (w i) := by
  classical
  exact VNEntropyAux.vonNeumann_mix_le_aux w hw0 σ hσ htr hmix

/-- The two-block pinching `P σ P + (1 - P) σ (1 - P)`. -/
def pinching (P σ : Matrix ι ι ℂ) : Matrix ι ι ℂ := P * σ * P + (1 - P) * σ * (1 - P)

/-- The pinching of a positive semidefinite matrix by a Hermitian `P` is positive semidefinite. -/
theorem pinching_posSemidef (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hσ : σ.PosSemidef) :
    (pinching P σ).PosSemidef := by
  have h1 := hσ.mul_mul_conjTranspose_same P
  have h2 := hσ.mul_mul_conjTranspose_same (1 - P)
  rw [hP.eq] at h1
  rw [conjTranspose_sub, conjTranspose_one, hP.eq] at h2
  exact h1.add h2

namespace VNEntropyAux

omit [Fintype ι] in
/-- The reflection `R = 2P - 1` of a Hermitian `P` is Hermitian. -/
lemma refl_conjTranspose (P : Matrix ι ι ℂ) (hP : P.IsHermitian) : (P + P - 1)ᴴ = P + P - 1 := by
  rw [conjTranspose_sub, conjTranspose_add, conjTranspose_one, hP.eq]

/-- The reflection `R = 2P - 1` of a projection `P` squares to one. -/
lemma refl_mul_self (P : Matrix ι ι ℂ) (hPP : P * P = P) : (P + P - 1) * (P + P - 1) = 1 := by
  have h : (P + P - 1) * (P + P - 1) =
      (P * P - P) + (P * P - P) + (P * P - P) + (P * P - P) + 1 := by noncomm_ring
  rw [h, hPP, sub_self]
  simp

/-- Pinching as an equal-weight mixture: `pinching P σ = ½ σ + ½ R σ R`, `R = 2P - 1`. -/
lemma pinching_eq_mix (P σ : Matrix ι ι ℂ) :
    pinching P σ = (((1 : ℝ) / 2 : ℝ) : ℂ) • σ +
      (((1 : ℝ) / 2 : ℝ) : ℂ) • ((P + P - 1) * σ * (P + P - 1)) := by
  have h : σ + (P + P - 1) * σ * (P + P - 1) = pinching P σ + pinching P σ := by
    unfold pinching; noncomm_ring
  rw [← smul_add, h, ← two_smul ℂ (pinching P σ), smul_smul]
  norm_num

/-- The equal-weight mixture of `σ` and `R σ R`, as a sum over `Fin 2`, is the pinching. -/
lemma sum_fin_two_eq_pinching (P σ : Matrix ι ι ℂ) :
    ∑ i : Fin 2, (((1 : ℝ) / 2 : ℝ) : ℂ) • ![σ, (P + P - 1) * σ * (P + P - 1)] i =
      pinching P σ := by
  rw [Fin.sum_univ_two, pinching_eq_mix]; rfl

/-- Conjugation by the reflection `R = 2P - 1` of an orthogonal projection preserves positivity
and the entropy. -/
lemma refl_conj (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) :
    ∃ h : ((P + P - 1) * σ * (P + P - 1)).PosSemidef,
      vonNeumann ((P + P - 1) * σ * (P + P - 1)) h.isHermitian = vonNeumann σ hσ.isHermitian := by
  obtain ⟨R, hR⟩ : ∃ R : Matrix ι ι ℂ, R = P + P - 1 := ⟨_, rfl⟩
  rw [← hR]
  have hRh : Rᴴ = R := by rw [hR]; exact refl_conjTranspose P hP
  have hRR : R * R = 1 := by rw [hR]; exact refl_mul_self P hPP
  have hRσ : (R * σ * R).PosSemidef := by
    have h := hσ.mul_mul_conjTranspose_same R; rwa [hRh] at h
  refine ⟨hRσ, ?_⟩
  rw [← CutLocalAux.vonNeumann_conj_isometry R (by rw [hRh, hRR]) σ hσ.isHermitian
    (by rw [hRh]; exact hRσ.isHermitian)]
  exact CutLocalAux.vonNeumann_congr _ _ (by rw [hRh])

end VNEntropyAux

/-- **Two-block pinching, lower bound** (paper: Lemma `lem:pinching`): for an orthogonal projection
`P`, pinching does not decrease the entropy. -/
theorem vonNeumann_le_pinching (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) (h' : (pinching P σ).IsHermitian) :
    vonNeumann σ hσ.isHermitian ≤ vonNeumann (pinching P σ) h' := by
  obtain ⟨hRσ, hS⟩ := VNEntropyAux.refl_conj P σ hP hPP hσ
  have hσ2 : ∀ i, (![σ, (P + P - 1) * σ * (P + P - 1)] i).PosSemidef :=
    Fin.forall_fin_two.2 ⟨hσ, hRσ⟩
  have hsum := VNEntropyAux.sum_fin_two_eq_pinching P σ
  have hm : (∑ i : Fin 2, (((1 : ℝ) / 2 : ℝ) : ℂ) •
      ![σ, (P + P - 1) * σ * (P + P - 1)] i).IsHermitian := by
    rw [hsum]; exact h'
  calc vonNeumann σ hσ.isHermitian
      = ∑ i : Fin 2, (1 : ℝ) / 2 *
          vonNeumann (![σ, (P + P - 1) * σ * (P + P - 1)] i) (hσ2 i).isHermitian := by
        rw [Fin.sum_univ_two]
        change _ = 1 / 2 * vonNeumann σ hσ.isHermitian +
          1 / 2 * vonNeumann ((P + P - 1) * σ * (P + P - 1)) hRσ.isHermitian
        rw [hS]; ring
    _ ≤ vonNeumann (∑ i : Fin 2, (((1 : ℝ) / 2 : ℝ) : ℂ) •
          ![σ, (P + P - 1) * σ * (P + P - 1)] i) hm :=
        vonNeumann_concave (fun _ => (1 : ℝ) / 2) (fun _ => by norm_num)
          (by rw [Fin.sum_univ_two]; norm_num) _ hσ2 hm
    _ = vonNeumann (pinching P σ) h' := CutLocalAux.vonNeumann_congr _ _ hsum

/-- **Two-block pinching, upper bound** (paper: Lemma `lem:pinching`): for a density matrix `σ` and
an orthogonal projection `P`, pinching increases the entropy by at most `ln 2`. -/
theorem vonNeumann_pinching_le (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) (htr : σ.trace = 1) (h' : (pinching P σ).IsHermitian) :
    vonNeumann (pinching P σ) h' ≤ vonNeumann σ hσ.isHermitian + Real.log 2 := by
  obtain ⟨hRσ, hS⟩ := VNEntropyAux.refl_conj P σ hP hPP hσ
  have htrR : ((P + P - 1) * σ * (P + P - 1)).trace = 1 := by
    rw [Matrix.trace_mul_comm, ← Matrix.mul_assoc, VNEntropyAux.refl_mul_self P hPP,
      Matrix.one_mul, htr]
  have hσ2 : ∀ i, (![σ, (P + P - 1) * σ * (P + P - 1)] i).PosSemidef :=
    Fin.forall_fin_two.2 ⟨hσ, hRσ⟩
  have htr2 : ∀ i, (![σ, (P + P - 1) * σ * (P + P - 1)] i).trace = 1 :=
    Fin.forall_fin_two.2 ⟨htr, htrR⟩
  have hsum := VNEntropyAux.sum_fin_two_eq_pinching P σ
  have hm : (∑ i : Fin 2, (((1 : ℝ) / 2 : ℝ) : ℂ) •
      ![σ, (P + P - 1) * σ * (P + P - 1)] i).IsHermitian := by
    rw [hsum]; exact h'
  have hη : Real.negMulLog ((1 : ℝ) / 2) + Real.negMulLog ((1 : ℝ) / 2) = Real.log 2 := by
    rw [Real.negMulLog, one_div, Real.log_inv]; ring
  calc vonNeumann (pinching P σ) h'
      = vonNeumann (∑ i : Fin 2, (((1 : ℝ) / 2 : ℝ) : ℂ) •
          ![σ, (P + P - 1) * σ * (P + P - 1)] i) hm :=
        CutLocalAux.vonNeumann_congr _ _ hsum.symm
    _ ≤ ∑ i : Fin 2, (1 : ℝ) / 2 *
          vonNeumann (![σ, (P + P - 1) * σ * (P + P - 1)] i) (hσ2 i).isHermitian +
          ∑ _i : Fin 2, Real.negMulLog ((1 : ℝ) / 2) :=
        vonNeumann_mix_le (fun _ => (1 : ℝ) / 2) (fun _ => by norm_num) _ hσ2 htr2 hm
    _ = vonNeumann σ hσ.isHermitian + Real.log 2 := by
        rw [Fin.sum_univ_two, Fin.sum_univ_two, hη]
        change 1 / 2 * vonNeumann σ hσ.isHermitian +
          1 / 2 * vonNeumann ((P + P - 1) * σ * (P + P - 1)) hRσ.isHermitian + _ = _
        rw [hS]; ring

end DFL
