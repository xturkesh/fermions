import DFL.Entropies

/-!
# Block decomposition of the entropy (paper: Proposition `prop:entropies`(c), the identity)

For a Hermitian matrix that is block diagonal with respect to a labelling `lab : ι → κ`, the von
Neumann entropy is the sum of the entropies of the diagonal blocks (`vonNeumann_block`), and for a
positive semidefinite block `B` of trace `p`, `S(B) = p S(B/p) + η(p)` (`vonNeumann_smul_trace`).
Hence `S(ρ) = H(p) + ∑_k p_k S(ρ^{(k)})`, `ρ^{(k)} = B_k / p_k` (`vonNeumann_block_decomp`). Applied to the
reduced states of the quench, block diagonal in the particle number `N_A` of `A` with weights the
number distribution `p_t`, this is the identity `S_X(t) = S_N(t) + ∑_k p_t(k) S(ρ_X^{(k)}(t))` for
`X ∈ {full, f}` (paper: Proposition `prop:entropies`(c)).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-- The diagonal block of `ρ` on the labels `{a | lab a = k}`. -/
def blockOf (ρ : Matrix ι ι ℂ) (lab : ι → κ) (k : κ) :
    Matrix {a // lab a = k} {a // lab a = k} ℂ :=
  ρ.submatrix (fun a => a.1) (fun a => a.1)

/-- The normalized block `B_k / Tr B_k` (the zero matrix if `Tr B_k = 0`). -/
def normBlock (ρ : Matrix ι ι ℂ) (lab : ι → κ) (k : κ) :
    Matrix {a // lab a = k} {a // lab a = k} ℂ :=
  ((blockOf ρ lab k).trace)⁻¹ • blockOf ρ lab k

set_option linter.unusedSectionVars false in
/-- Blocks of a Hermitian matrix are Hermitian. -/
theorem blockOf_isHermitian (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ) (k : κ) :
    (blockOf ρ lab k).IsHermitian :=
  hρ.submatrix _

/-- Normalized blocks of a Hermitian matrix are Hermitian. -/
theorem normBlock_isHermitian (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ) (k : κ) :
    (normBlock ρ lab k).IsHermitian := by
  have hB := blockOf_isHermitian ρ hρ lab k
  refine hB.smul ?_
  show star ((blockOf ρ lab k).trace)⁻¹ = ((blockOf ρ lab k).trace)⁻¹
  rw [star_inv₀, ← Matrix.trace_conjTranspose, hB.eq]

namespace EntropyBlocksAux

/-- A block-diagonal matrix has characteristic polynomial the product of those of its blocks. -/
lemma charpoly_eq_prod_blockOf (ρ : Matrix ι ι ℂ) (lab : ι → κ)
    (hblock : ∀ a b, lab a ≠ lab b → ρ a b = 0) :
    ρ.charpoly = ∏ k, (blockOf ρ lab k).charpoly := by
  set e := Fintype.equivFin κ
  set b : ι → Fin (Fintype.card κ) := fun a => e (lab a) with hb
  have hBT : ρ.BlockTriangular b := fun i j hij =>
    hblock i j fun h => absurd hij (by simp [hb, h])
  have h1 : ρ.charpoly = ∏ a, (ρ.toSquareBlock b a).charpoly := by
    rw [Matrix.charpoly, hBT.charmatrix.det_fintype]
    simp only [Matrix.charpoly, charmatrix_toSquareBlock]
  rw [h1, ← Equiv.prod_comp e (fun a => (ρ.toSquareBlock b a).charpoly)]
  refine Finset.prod_congr rfl fun k _ => ?_
  rw [← Matrix.charpoly_reindex (Equiv.subtypeEquivRight (fun x => e.apply_eq_iff_eq) :
    {x // b x = e k} ≃ {x // lab x = k})]
  rfl

/-- The entropy of a Hermitian matrix from the roots of its characteristic polynomial. -/
lemma vonNeumann_eq_roots {μ : Type*} [Fintype μ] [DecidableEq μ] (A : Matrix μ μ ℂ)
    (hA : A.IsHermitian) :
    vonNeumann A hA = (A.charpoly.roots.map (fun z => Real.negMulLog (RCLike.re z))).sum := by
  rw [hA.roots_charpoly_eq_eigenvalues, Multiset.map_map, vonNeumann, Finset.sum_eq_multiset_sum]
  rfl

/-- The entropy of a real multiple of a Hermitian matrix: `S(c A) = ∑_i η(c λ_i)`. -/
lemma vonNeumann_ofReal_smul {μ : Type*} [Fintype μ] [DecidableEq μ] (A : Matrix μ μ ℂ)
    (hA : A.IsHermitian) (c : ℝ) (h : ((c : ℂ) • A).IsHermitian) :
    vonNeumann ((c : ℂ) • A) h = ∑ i, Real.negMulLog (c * hA.eigenvalues i) := by
  set U : Matrix μ μ ℂ := (hA.eigenvectorUnitary : Matrix μ μ ℂ)
  have hUU : Uᴴ * U = 1 := by
    rw [← Matrix.star_eq_conjTranspose]
    exact Unitary.coe_star_mul_self hA.eigenvectorUnitary
  have hD : (diagonal (fun i => ((c * hA.eigenvalues i : ℝ) : ℂ))).IsHermitian := by
    refine isHermitian_diagonal_of_self_adjoint _ ?_
    funext i
    simp
  have hdiag : diagonal (fun i => ((c * hA.eigenvalues i : ℝ) : ℂ)) =
      (c : ℂ) • diagonal (RCLike.ofReal ∘ hA.eigenvalues) := by
    ext i j
    by_cases hij : i = j
    · subst hij
      simp
    · simp [hij]
  have heq : (c : ℂ) • A = U * diagonal (fun i => ((c * hA.eigenvalues i : ℝ) : ℂ)) * Uᴴ := by
    rw [hdiag, Matrix.mul_smul, Matrix.smul_mul]
    congr 1
    conv_lhs => rw [hA.spectral_theorem]
    rw [Unitary.conjStarAlgAut_apply, Matrix.star_eq_conjTranspose]
  have hH' : (U * diagonal (fun i => ((c * hA.eigenvalues i : ℝ) : ℂ)) * Uᴴ).IsHermitian :=
    heq ▸ h
  rw [CutLocalAux.vonNeumann_congr h hH' heq, CutLocalAux.vonNeumann_conj_isometry U hUU _ hD]
  exact vonNeumann_diagonal (fun i => c * hA.eigenvalues i) hD

end EntropyBlocksAux

/-- **Additivity over blocks**: the entropy of a block-diagonal Hermitian matrix is the sum of the
entropies of its blocks. -/
theorem vonNeumann_block (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ)
    (hblock : ∀ a b, lab a ≠ lab b → ρ a b = 0) :
    vonNeumann ρ hρ = ∑ k, vonNeumann (blockOf ρ lab k) (blockOf_isHermitian ρ hρ lab k) := by
  have hne : (∏ k, (blockOf ρ lab k).charpoly) ≠ 0 :=
    (Polynomial.monic_prod_of_monic _ _ fun k _ => (blockOf ρ lab k).charpoly_monic).ne_zero
  rw [EntropyBlocksAux.vonNeumann_eq_roots ρ hρ,
    EntropyBlocksAux.charpoly_eq_prod_blockOf ρ lab hblock, Polynomial.roots_prod _ _ hne,
    Multiset.map_bind, Multiset.sum_bind, Finset.sum_eq_multiset_sum]
  congr 1
  refine Multiset.map_congr rfl fun k _ => ?_
  rw [EntropyBlocksAux.vonNeumann_eq_roots]

/-- **Scaling**: for a positive semidefinite `B` of trace `p`, `S(B) = p S(B/p) + η(p)`. -/
theorem vonNeumann_smul_trace {μ : Type*} [Fintype μ] [DecidableEq μ] (B : Matrix μ μ ℂ)
    (hB : B.PosSemidef) (hB' : ((B.trace)⁻¹ • B).IsHermitian) :
    vonNeumann B hB.isHermitian =
      (B.trace).re * vonNeumann ((B.trace)⁻¹ • B) hB' + Real.negMulLog (B.trace).re := by
  set lam := hB.isHermitian.eigenvalues
  set p := (B.trace).re with hp
  have hlam0 : ∀ i, 0 ≤ lam i := hB.eigenvalues_nonneg
  have htr : B.trace = ((∑ i, lam i : ℝ) : ℂ) := by
    rw [hB.isHermitian.trace_eq_sum_eigenvalues, Complex.ofReal_sum]
    rfl
  have hpsum : p = ∑ i, lam i := by rw [hp, htr, Complex.ofReal_re]
  have hBtr : B.trace = (p : ℂ) := by rw [htr, hpsum]
  have hinv : (B.trace)⁻¹ = ((p⁻¹ : ℝ) : ℂ) := by rw [hBtr, Complex.ofReal_inv]
  have hH2 : (((p⁻¹ : ℝ) : ℂ) • B).IsHermitian := by rw [← hinv]; exact hB'
  have h2 : vonNeumann ((B.trace)⁻¹ • B) hB' = ∑ i, Real.negMulLog (p⁻¹ * lam i) := by
    rw [CutLocalAux.vonNeumann_congr hB' hH2 (by rw [hinv]),
      EntropyBlocksAux.vonNeumann_ofReal_smul B hB.isHermitian p⁻¹ hH2]
  rw [h2]
  show ∑ i, Real.negMulLog (lam i) = p * ∑ i, Real.negMulLog (p⁻¹ * lam i) + Real.negMulLog p
  have hp0 : 0 ≤ p := hpsum ▸ Finset.sum_nonneg fun i _ => hlam0 i
  rcases hp0.eq_or_lt with h0 | hpos
  · -- `p = 0`: all eigenvalues vanish
    have hz : ∀ i, lam i = 0 := fun i =>
      (Finset.sum_eq_zero_iff_of_nonneg (fun i _ => hlam0 i)).1 (by rw [← hpsum, ← h0]) i
        (Finset.mem_univ _)
    simp [hz, ← h0]
  · have hpne : p ≠ 0 := hpos.ne'
    have key : ∀ i, p * Real.negMulLog (p⁻¹ * lam i) =
        Real.negMulLog (lam i) + lam i * Real.log p := by
      intro i
      rw [Real.negMulLog_mul]
      simp only [Real.negMulLog, Real.log_inv]
      field_simp
      ring
    rw [Finset.mul_sum]
    simp only [key]
    rw [Finset.sum_add_distrib, ← Finset.sum_mul, ← hpsum, Real.negMulLog]
    ring

/-- **Block decomposition**: for a positive semidefinite block-diagonal `ρ`,
`S(ρ) = ∑_k η(p_k) + ∑_k p_k S(ρ^{(k)})` with `p_k = Tr B_k`. -/
theorem vonNeumann_block_decomp (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (lab : ι → κ)
    (hblock : ∀ a b, lab a ≠ lab b → ρ a b = 0) :
    vonNeumann ρ hρ.isHermitian =
      ∑ k, Real.negMulLog ((blockOf ρ lab k).trace).re +
        ∑ k, ((blockOf ρ lab k).trace).re *
          vonNeumann (normBlock ρ lab k) (normBlock_isHermitian ρ hρ.isHermitian lab k) := by
  rw [vonNeumann_block ρ hρ.isHermitian lab hblock, add_comm, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  exact vonNeumann_smul_trace (blockOf ρ lab k) (hρ.submatrix _)
    (normBlock_isHermitian ρ hρ.isHermitian lab k)

variable {N : ℕ}

/-- The matter particle number of `A` as a label in `{0, …, N}` of the configurations of `A`. -/
def cntA (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ × Cfg ℓ) : Fin (N + 1) :=
  ⟨(univ.filter (fun i => a.1 i = 1)).card, by
    have := Finset.card_le_univ (univ.filter (fun i => a.1 i = 1))
    simp only [Fintype.card_fin] at this
    omega⟩

/-- The particle number of the bare modes of `A` as a label in `{0, …, N}`. -/
def cntBare (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) : Fin (N + 1) :=
  ⟨(univ.filter (fun i => a i = 1)).card, by
    have := Finset.card_le_univ (univ.filter (fun i => a i = 1))
    simp only [Fintype.card_fin] at this
    omega⟩

namespace EntropyBlocksAux

/-- `x = 1 ↔ x ≠ 0` in `ZMod 2`. -/
lemma zmod2_eq_one_iff (x : ZMod 2) : x = 1 ↔ x ≠ 0 := by
  revert x; decide

/-- `cntBare` is the particle number `cnt`. -/
lemma cntBare_val (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) :
    (cntBare ℓ hℓ a : ℕ) = TaxonomyAux.cnt a := by
  show (univ.filter (fun i => a i = 1)).card = (univ.filter (fun i => a i ≠ 0)).card
  exact congrArg Finset.card (Finset.filter_congr fun i _ => zmod2_eq_one_iff (a i))

/-- `cntA` is the particle number `cnt` of the matter part. -/
lemma cntA_val (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ × Cfg ℓ) :
    (cntA ℓ hℓ a : ℕ) = TaxonomyAux.cnt a.1 :=
  cntBare_val (N := N) ℓ hℓ a.1

omit [DecidableEq ι] [Fintype κ] in
/-- The trace of a block is the sum of the diagonal entries with the given label. -/
lemma trace_blockOf (ρ : Matrix ι ι ℂ) (lab : ι → κ) (k : κ) :
    (blockOf ρ lab k).trace = ∑ a ∈ univ.filter (fun a => lab a = k), ρ a a :=
  (Finset.sum_subtype (p := fun a => lab a = k) _ (fun a => by simp) (fun a => ρ a a)).symm

end EntropyBlocksAux

/-- **Entropy taxonomy, the identity for `S_full`** (paper: Proposition `prop:entropies`(c)):
`S_full(t) = S_N(t) + ∑_k p_t(k) S(ρ_full^{(k)}(t))` for Fock initial states. -/
theorem fullEntropy_eq_numberEntropy_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmLGT_isHermitian ℓ hℓ _) =
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) +
        ∑ k : Fin (N + 1), lgtNumDist J h (segment ℓ) (Pi.single n 1) t k *
          vonNeumann (normBlock (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (cntA ℓ hℓ) k)
            (normBlock_isHermitian _ (rdmLGT_isHermitian ℓ hℓ _) _ k) := by
  set Ψ := lgtState J h (Pi.single n 1) t with hΨ
  have hpsd : (rdmLGT ℓ hℓ Ψ).PosSemidef := EntropiesAux.rdmLGT_posSemidef ℓ hℓ Ψ
  have hblock : ∀ a b, cntA ℓ hℓ a ≠ cntA ℓ hℓ b → rdmLGT ℓ hℓ Ψ a b = 0 := fun a b hne =>
    EntropiesAux.rdmLGT_lgtState_block J h n t ℓ hℓ a b fun heq =>
      hne (Fin.ext (by rw [EntropyBlocksAux.cntA_val, EntropyBlocksAux.cntA_val, heq]))
  have htr : ∀ k : Fin (N + 1), ((blockOf (rdmLGT ℓ hℓ Ψ) (cntA ℓ hℓ) k).trace).re =
      lgtNumDist J h (segment ℓ) (Pi.single n 1) t k := by
    intro k
    rw [EntropiesAux.lgtNumDist_eq_rdmLGT J h _ t ℓ hℓ k, EntropyBlocksAux.trace_blockOf, ← hΨ]
    congr 1
    refine Finset.sum_congr (Finset.filter_congr fun a _ => ?_) fun _ _ => rfl
    rw [Fin.ext_iff, EntropyBlocksAux.cntA_val]
  rw [vonNeumann_block_decomp (rdmLGT ℓ hℓ Ψ) hpsd (cntA ℓ hℓ) hblock]
  simp only [htr]
  congr 1
  exact Fin.sum_univ_eq_sum_range
    (fun k => Real.negMulLog (lgtNumDist J h (segment ℓ) (Pi.single n 1) t k)) (N + 1)

/-- **Entropy taxonomy, the identity for `S_f`** (paper: Proposition `prop:entropies`(c)):
`S_f(t) = S_N(t) + ∑_k p_t(k) S(ρ_f^{(k)}(t))` for Fock initial states. -/
theorem bareEntropy_eq_numberEntropy_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmBare_isHermitian ℓ hℓ _) =
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) +
        ∑ k : Fin (N + 1), lgtNumDist J h (segment ℓ) (Pi.single n 1) t k *
          vonNeumann (normBlock (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t)) (cntBare ℓ hℓ) k)
            (normBlock_isHermitian _ (rdmBare_isHermitian ℓ hℓ _) _ k) := by
  set Ψ := lgtState J h (Pi.single n 1) t with hΨ
  have hpsd : (rdmBare ℓ hℓ Ψ).PosSemidef := EntropiesAux.rdmBare_posSemidef ℓ hℓ Ψ
  have hblock : ∀ a b, cntBare ℓ hℓ a ≠ cntBare ℓ hℓ b → rdmBare ℓ hℓ Ψ a b = 0 := fun a b hne =>
    Finset.sum_eq_zero fun c _ => EntropiesAux.rdmLGT_lgtState_block J h n t ℓ hℓ (a, c) (b, c)
      fun heq => hne (Fin.ext (by
        rw [EntropyBlocksAux.cntBare_val, EntropyBlocksAux.cntBare_val]
        exact heq))
  have htr : ∀ k : Fin (N + 1), ((blockOf (rdmBare ℓ hℓ Ψ) (cntBare ℓ hℓ) k).trace).re =
      lgtNumDist J h (segment ℓ) (Pi.single n 1) t k := by
    intro k
    rw [EntropiesAux.lgtNumDist_eq_rdmLGT J h _ t ℓ hℓ k, EntropyBlocksAux.trace_blockOf, ← hΨ]
    have hf : univ.filter (fun a => cntBare ℓ hℓ a = k) =
        univ.filter (fun a : Cfg ℓ => TaxonomyAux.cnt a = k) :=
      Finset.filter_congr fun a _ => by rw [Fin.ext_iff, EntropyBlocksAux.cntBare_val]
    rw [hf]
    congr 1
    simp only [rdmBare]
    rw [Finset.sum_filter, Finset.sum_filter, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun a _ => ?_
    split_ifs <;> simp
  rw [vonNeumann_block_decomp (rdmBare ℓ hℓ Ψ) hpsd (cntBare ℓ hℓ) hblock]
  simp only [htr]
  congr 1
  exact Fin.sum_univ_eq_sum_range
    (fun k => Real.negMulLog (lgtNumDist J h (segment ℓ) (Pi.single n 1) t k)) (N + 1)

end DFL
