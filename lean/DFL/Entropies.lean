import DFL.VNEntropy
import DFL.Taxonomy

/-!
# Entropy taxonomy (paper: Lemma `lem:cutlocal`(b), Proposition `prop:entropies`(a)–(c))

For the initial state `|↑⋯↑⟩_z ⊗ ψ` and the region `A` = (sites `< ℓ`, links `< ℓ`):

* `cutlocal_defect`: the duality `U` is local across the cut up to one controlled gate, so for every
  normalized state `Φ` of the dual variables, `|S_A(UΦ) - S_{A'}(Φ)| ≤ ln 2` (pinching argument);
  `lgt_duality_defect` is the case `Φ = Φ(t) = U†Ψ(t)` (paper: `prop:entropies`(a)).
* `ρ^q_A(t)` (`secRdm`) is the reduced state of the sector-`q` free-fermion state
  `ψ_q(t) = Γ(e^{-ith_q}) ψ`, `ρ_mix(t) = E_q ρ^q_A(t)` (`mixRdm`), `S_G(t) = E_q S(ρ^q_A(t))`
  (`condEntropy`), and `H_t` (`sectorLawEntropy`) is the Shannon entropy of the law of the random
  matrix `q ↦ ρ^q_A(t)`. Then `S_G ≤ S_mix ≤ min(S_G + H_t, ℓ ln 2)`, `H_t ≤ N ln 2`, `H_0 = 0`, and
  `ρ_mix` is the marginal of `ρ_{A'}(Φ(t))` on the modes (paper: `prop:entropies`(b)).
* For Fock `ψ`: `p_t = E_q p^q_t`; the number entropy `S_N` is a lower bound for `S_mix` and for the
  entropy `S_f` of the bare fermion modes (and for `S_full`, `DFL.lgt_numberEntropy_le_full`);
  `S_G ≥ E_q H(p^q_t)`, and `E_q H(p^q_t) ≤ S_N ≤ E_q H(p^q_t) + H_t` (paper: `prop:entropies`(c)).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Cut locality of the duality -/

namespace EntropiesAux

/-! #### Generic facts on reduced density matrices -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `a + a = 0` in `ZMod 2`. -/
lemma zmod2_add_self (a : ZMod 2) : a + a = 0 := CharTwo.add_self_eq_zero a

/-- `z conj z = |z|²`. -/
lemma mul_star_self_eq (z : ℂ) : z * star z = (Complex.normSq z : ℂ) := by
  rw [Complex.star_def, Complex.mul_conj]

/-- There are `2^K` configurations of `K` two-level systems. -/
lemma card_cfg (K : ℕ) : Fintype.card (Cfg K) = 2 ^ K := by
  simp [Cfg, ZMod.card]

/-- The reduced density matrix of an LGT pure state is positive semidefinite (`ρ = X X†`). -/
lemma rdmLGT_posSemidef (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmLGT ℓ hℓ Ψ).PosSemidef := by
  have h : rdmLGT ℓ hℓ Ψ = (Matrix.of fun a b => Ψ ((splitLGT ℓ hℓ).symm (a, b))) *
      (Matrix.of fun a b => Ψ ((splitLGT ℓ hℓ).symm (a, b)))ᴴ := by
    ext a a'
    simp [rdmLGT, Matrix.mul_apply, Matrix.conjTranspose_apply]
  rw [h]
  exact Matrix.posSemidef_self_mul_conjTranspose _

/-- The reduced density matrix of a pure state is positive semidefinite (`ρ = X X†`). -/
lemma rdm_posSemidef (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) : (rdm ℓ hℓ v).PosSemidef := by
  have h : rdm ℓ hℓ v = (Matrix.of fun a b => v ((splitCfg ℓ hℓ).symm (a, b))) *
      (Matrix.of fun a b => v ((splitCfg ℓ hℓ).symm (a, b)))ᴴ := by
    ext a a'
    simp [rdm, Matrix.mul_apply, Matrix.conjTranspose_apply]
  rw [h]
  exact Matrix.posSemidef_self_mul_conjTranspose _

/-- The trace of `ρ_A(Ψ)` is the squared norm of `Ψ`. -/
lemma trace_rdmLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmLGT ℓ hℓ Ψ).trace = ((∑ p, Complex.normSq (Ψ p) : ℝ) : ℂ) := by
  simp only [Matrix.trace, Matrix.diag, rdmLGT, mul_star_self_eq, Complex.ofReal_sum]
  rw [← Fintype.sum_prod_type' (fun a b => (Complex.normSq (Ψ ((splitLGT ℓ hℓ).symm (a, b))) : ℂ))]
  exact Equiv.sum_comp (splitLGT ℓ hℓ).symm (fun p => (Complex.normSq (Ψ p) : ℂ))

/-- The trace of `ρ_A(v)` is the squared norm of `v`. -/
lemma trace_rdm (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) :
    (rdm ℓ hℓ v).trace = ((∑ n, Complex.normSq (v n) : ℝ) : ℂ) := by
  simp only [Matrix.trace, Matrix.diag, rdm, mul_star_self_eq, Complex.ofReal_sum]
  rw [← Fintype.sum_prod_type' (fun a b => (Complex.normSq (v ((splitCfg ℓ hℓ).symm (a, b))) : ℂ))]
  exact Equiv.sum_comp (splitCfg ℓ hℓ).symm (fun n => (Complex.normSq (v n) : ℂ))

/-- Relabelling the basis by a bijection does not change the von Neumann entropy. -/
lemma vonNeumann_submatrix_equiv {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] (M : Matrix ι ι ℂ) (hM : M.IsHermitian) (e : κ ≃ ι)
    (h' : (M.submatrix e e).IsHermitian) : vonNeumann (M.submatrix e e) h' = vonNeumann M hM := by
  refine CutLocalAux.vonNeumann_eq_of_charpoly _ h' M hM 0 0 ?_
  have h := Matrix.charpoly_reindex e.symm M
  rw [Matrix.reindex_apply, Equiv.symm_symm] at h
  rw [pow_zero, one_mul, one_mul, h]

/-! #### The cut link -/

/-- The label of the cut link `b_{ℓ-1}` (the label `0` of the fictitious link `b_{-1}` if
`ℓ = 0`). -/
def cutLabel {ℓ : ℕ} (β : Cfg ℓ) : ZMod 2 := if h : 0 < ℓ then β ⟨ℓ - 1, by omega⟩ else 0

/-- `prevLink` of a recombined link configuration, left of the cut. -/
lemma prevLink_splitCfg_symm_lt (ℓ : ℕ) (hℓ : ℓ ≤ N) (β : Cfg ℓ) (δ : Cfg (N - ℓ)) (k : Fin N)
    (hk : k.val < ℓ) : prevLink ((splitCfg ℓ hℓ).symm (β, δ)) k = prevLink β ⟨k.val, hk⟩ := by
  unfold prevLink
  by_cases h0 : k.val = 0
  · simp [h0]
  · simp only [h0, dite_false]
    exact CutLocalAux.splitCfg_symm_apply_lt _ _ _ _ _ (by simp; omega)

/-- A recombined configuration, right of the cut. -/
lemma splitCfg_symm_natAdd (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) (b : Cfg (N - ℓ)) (i : Fin (N - ℓ))
    (h : ℓ + i.val < N) : (splitCfg ℓ hℓ).symm (a, b) ⟨ℓ + i.val, h⟩ = b i := by
  rw [CutLocalAux.splitCfg_symm_apply_ge _ _ _ _ _ (by simp)]
  congr 1
  ext
  simp

/-- `prevLink` of a recombined link configuration, right of the cut: the first link of `B` sees
the cut link. -/
lemma prevLink_splitCfg_symm_natAdd (ℓ : ℕ) (hℓ : ℓ ≤ N) (β : Cfg ℓ) (δ : Cfg (N - ℓ))
    (i : Fin (N - ℓ)) (h : ℓ + i.val < N) :
    prevLink ((splitCfg ℓ hℓ).symm (β, δ)) ⟨ℓ + i.val, h⟩ =
      prevLink δ i + if i.val = 0 then cutLabel β else 0 := by
  by_cases hi : i.val = 0
  · by_cases hl : ℓ = 0
    · subst hl
      simp [prevLink, cutLabel, hi]
    · have hpos : 0 < ℓ := Nat.pos_of_ne_zero hl
      simp only [prevLink, hi, add_zero, hl, dite_false, dite_true, ite_true, zero_add, cutLabel,
        hpos]
      exact CutLocalAux.splitCfg_symm_apply_lt _ _ _ _ _ (by simp; omega)
  · simp only [prevLink, hi, dite_false, ite_false, add_zero]
    rw [dite_eq_right (by omega)]
    rw [CutLocalAux.splitCfg_symm_apply_ge _ _ _ _ _ (by simp; omega)]
    congr 1
    ext
    simp
    omega

/-- Shifting all link labels by a constant only changes `prevLink` away from the first site. -/
lemma prevLink_add_const {K : ℕ} (δ : Cfg K) (c : ZMod 2) (i : Fin K) :
    prevLink (δ + fun _ => c) i = prevLink δ i + if i.val = 0 then 0 else c := by
  by_cases hi : i.val = 0 <;> simp [prevLink, hi]

/-- **Charges of a recombined LGT configuration**: on `A` they are the charges of the chain
`{0, …, ℓ-1}`; on `B` they are the charges of the chain `{ℓ, …, N-1}` with all link labels of `B`
shifted by the label of the cut link. -/
lemma splitCfg_chargeOf (ℓ : ℕ) (hℓ : ℓ ≤ N) (α β : Cfg ℓ) (γ δ : Cfg (N - ℓ)) :
    splitCfg ℓ hℓ (chargeOf ((splitCfg ℓ hℓ).symm (α, γ)) ((splitCfg ℓ hℓ).symm (β, δ))) =
      (chargeOf α β, chargeOf γ (δ + fun _ => cutLabel β)) := by
  ext i
  · rw [CutLocalAux.splitCfg_apply_fst]
    simp only [chargeOf]
    rw [prevLink_splitCfg_symm_lt ℓ hℓ β δ _ (by simp),
      CutLocalAux.splitCfg_symm_apply_lt _ _ _ _ _ (by simp),
      CutLocalAux.splitCfg_symm_apply_lt _ _ _ _ _ (by simp)]
  · rw [CutLocalAux.splitCfg_apply_snd]
    simp only [chargeOf, Pi.add_apply]
    rw [prevLink_splitCfg_symm_natAdd, splitCfg_symm_natAdd, splitCfg_symm_natAdd,
      prevLink_add_const]
    by_cases hi : i.val = 0
    · simp only [hi, ite_true]
      ring
    · simp only [hi, ite_false]
      linear_combination (-1 : ZMod 2) * zmod2_add_self (cutLabel β)

/-- The dual coordinates of a recombined LGT configuration. -/
lemma dualEquiv_symm_splitLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) (α β : Cfg ℓ) (γ δ : Cfg (N - ℓ)) :
    dualEquiv.symm ((splitLGT ℓ hℓ).symm ((α, β), (γ, δ))) =
      (splitLGT ℓ hℓ).symm ((α, chargeOf α β), (γ, chargeOf γ (δ + fun _ => cutLabel β))) := by
  rw [dualEquiv_symm_apply, CutLocalAux.splitLGT_symm_apply, CutLocalAux.splitLGT_symm_apply]
  refine Prod.ext rfl ?_
  exact (Equiv.eq_symm_apply _).2 (splitCfg_chargeOf ℓ hℓ α β γ δ)

/-- For fixed occupations, `x ↦ chargeOf n x` is a bijection. -/
def chargeEquiv {K : ℕ} (n : Cfg K) : Cfg K ≃ Cfg K :=
  ⟨chargeOf n, linkOf n, linkOf_chargeOf n, chargeOf_linkOf n⟩

/-- **Entries of `ρ_A(UΦ)`**: for equal cut labels they are the entries of `ρ_{A'}(Φ)` at the dual
labels `(α, chargeOf α β)` (substitute `δ ↦ chargeOf γ (δ + c)` in the sum over `B`). -/
lemma rdmLGT_dual_apply (ℓ : ℕ) (hℓ : ℓ ≤ N) (Φ : Cfg N × Cfg N → ℂ) (α β α' β' : Cfg ℓ)
    (hc : cutLabel β = cutLabel β') :
    rdmLGT ℓ hℓ (fun p => Φ (dualEquiv.symm p)) (α, β) (α', β') =
      rdmLGT ℓ hℓ Φ (α, chargeOf α β) (α', chargeOf α' β') := by
  simp only [rdmLGT, Fintype.sum_prod_type, dualEquiv_symm_splitLGT, ← hc]
  refine Finset.sum_congr rfl fun γ _ => ?_
  exact Equiv.sum_comp ((Equiv.addRight fun _ => cutLabel β).trans (chargeEquiv γ))
    (fun q => Φ ((splitLGT ℓ hℓ).symm ((α, chargeOf α β), (γ, q))) *
      star (Φ ((splitLGT ℓ hℓ).symm ((α', chargeOf α' β'), (γ, q)))))

/-! #### The pinching by the cut-link projection -/

/-- The projection `P⁺` of `A` onto the cut-link label `0`. -/
def cutProj (ℓ : ℕ) : Matrix (Cfg ℓ × Cfg ℓ) (Cfg ℓ × Cfg ℓ) ℂ :=
  diagonal fun a => if cutLabel a.2 = 0 then 1 else 0

/-- `P⁺` is Hermitian. -/
lemma cutProj_isHermitian (ℓ : ℕ) : (cutProj ℓ).IsHermitian := by
  refine isHermitian_diagonal_of_self_adjoint _ ?_
  funext a
  simp only [Pi.star_apply]
  split_ifs <;> simp

/-- `P⁺` is idempotent. -/
lemma cutProj_mul_self (ℓ : ℕ) : cutProj ℓ * cutProj ℓ = cutProj ℓ := by
  rw [cutProj, diagonal_mul_diagonal]
  congr 1
  funext a
  split_ifs <;> simp

/-- The pinching by `P⁺` keeps exactly the entries between equal cut labels. -/
lemma pinching_cutProj_apply (ℓ : ℕ) (σ : Matrix (Cfg ℓ × Cfg ℓ) (Cfg ℓ × Cfg ℓ) ℂ)
    (a b : Cfg ℓ × Cfg ℓ) :
    pinching (cutProj ℓ) σ a b = if cutLabel a.2 = cutLabel b.2 then σ a b else 0 := by
  have h1 : (1 : Matrix (Cfg ℓ × Cfg ℓ) (Cfg ℓ × Cfg ℓ) ℂ) - cutProj ℓ =
      diagonal fun a => if cutLabel a.2 = 0 then 0 else 1 := by
    rw [cutProj, ← diagonal_one, diagonal_sub]
    congr 1
    funext a
    split_ifs <;> simp
  rw [pinching, h1, cutProj, Matrix.add_apply]
  simp only [Matrix.mul_diagonal, Matrix.diagonal_mul]
  rcases zmod2_cases (cutLabel a.2) with ha | ha <;>
    rcases zmod2_cases (cutLabel b.2) with hb | hb <;> simp [ha, hb]

/-- The pinchings of `ρ_A(UΦ)` and of the relabelled `ρ_{A'}(Φ)` by `P⁺` coincide. -/
lemma pinching_rdmLGT_dual (ℓ : ℕ) (hℓ : ℓ ≤ N) (Φ : Cfg N × Cfg N → ℂ) :
    pinching (cutProj ℓ) (rdmLGT ℓ hℓ (fun p => Φ (dualEquiv.symm p))) =
      pinching (cutProj ℓ) ((rdmLGT ℓ hℓ Φ).submatrix dualEquiv.symm dualEquiv.symm) := by
  ext ⟨α, β⟩ ⟨α', β'⟩
  rw [pinching_cutProj_apply, pinching_cutProj_apply]
  split_ifs with hc
  · rw [rdmLGT_dual_apply ℓ hℓ Φ α β α' β' hc, submatrix_apply, dualEquiv_symm_apply,
      dualEquiv_symm_apply]
  · rfl

end EntropiesAux

open EntropiesAux in
/-- **Duality defect** (paper: Lemma `lem:cutlocal`(b)): for every normalized state `Φ` of the dual
variables `(n, q)`, the entanglement entropy of `A` in `UΦ`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, and
that of `A'` = (modes `< ℓ`, charges `< ℓ`) in `Φ` differ by at most `ln 2`. -/
theorem cutlocal_defect (ℓ : ℕ) (hℓ : ℓ ≤ N) (Φ : Cfg N × Cfg N → ℂ)
    (hΦ : ∑ p, Complex.normSq (Φ p) = 1) :
    |vonNeumann (rdmLGT ℓ hℓ (fun p => Φ (dualEquiv.symm p))) (rdmLGT_isHermitian ℓ hℓ _) -
      vonNeumann (rdmLGT ℓ hℓ Φ) (rdmLGT_isHermitian ℓ hℓ _)| ≤ Real.log 2 := by
  set ρ := rdmLGT ℓ hℓ (fun p => Φ (dualEquiv.symm p)) with hρdef
  set ρ' := (rdmLGT ℓ hℓ Φ).submatrix (dualEquiv (N := ℓ)).symm dualEquiv.symm with hρ'def
  have hρ : ρ.PosSemidef := rdmLGT_posSemidef ℓ hℓ _
  have hρ' : ρ'.PosSemidef := (rdmLGT_posSemidef ℓ hℓ Φ).submatrix _
  have htr : ρ.trace = 1 := by
    rw [hρdef, trace_rdmLGT,
      Equiv.sum_comp dualEquiv.symm (fun p => Complex.normSq (Φ p)), hΦ, Complex.ofReal_one]
  have htr' : ρ'.trace = 1 := by
    rw [hρ'def]
    simp only [Matrix.trace, Matrix.diag, submatrix_apply]
    rw [Equiv.sum_comp dualEquiv.symm (fun a => rdmLGT ℓ hℓ Φ a a)]
    have := trace_rdmLGT ℓ hℓ Φ
    rw [hΦ, Complex.ofReal_one] at this
    exact this
  have hP := cutProj_isHermitian ℓ
  have hPP := cutProj_mul_self ℓ
  have hS' : vonNeumann ρ' hρ'.isHermitian =
      vonNeumann (rdmLGT ℓ hℓ Φ) (rdmLGT_isHermitian ℓ hℓ Φ) :=
    vonNeumann_submatrix_equiv _ _ _ _
  have hpinH := (pinching_posSemidef (cutProj ℓ) ρ hP hρ).isHermitian
  have hpinH' := (pinching_posSemidef (cutProj ℓ) ρ' hP hρ').isHermitian
  have h1 : vonNeumann ρ (rdmLGT_isHermitian ℓ hℓ _) ≤ vonNeumann (pinching (cutProj ℓ) ρ) hpinH :=
    vonNeumann_le_pinching (cutProj ℓ) ρ hP hPP hρ hpinH
  have h2 : vonNeumann (pinching (cutProj ℓ) ρ') hpinH' ≤
      vonNeumann ρ' hρ'.isHermitian + Real.log 2 :=
    vonNeumann_pinching_le (cutProj ℓ) ρ' hP hPP hρ' htr' hpinH'
  have h3 : vonNeumann ρ' hρ'.isHermitian ≤ vonNeumann (pinching (cutProj ℓ) ρ') hpinH' :=
    vonNeumann_le_pinching (cutProj ℓ) ρ' hP hPP hρ' hpinH'
  have h4 : vonNeumann (pinching (cutProj ℓ) ρ) hpinH ≤
      vonNeumann ρ (rdmLGT_isHermitian ℓ hℓ _) + Real.log 2 :=
    vonNeumann_pinching_le (cutProj ℓ) ρ hP hPP hρ htr hpinH
  have h5 : vonNeumann (pinching (cutProj ℓ) ρ) hpinH =
      vonNeumann (pinching (cutProj ℓ) ρ') hpinH' :=
    CutLocalAux.vonNeumann_congr _ _ (pinching_rdmLGT_dual ℓ hℓ Φ)
  rw [abs_sub_le_iff]
  constructor <;> linarith

namespace EntropiesAux

/-! #### The quench state in the dual frame -/

/-- A unitary matrix preserves `∑ |v_i|²`. -/
lemma sum_normSq_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] (U : Matrix ι ι ℂ)
    (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h := AreaLawAux.nsq_mulVec_unitary U hU v
  simpa [AreaLawAux.nsq, Complex.normSq_eq_norm_sq] using h

/-- The sector propagators preserve the norm. -/
lemma sum_normSq_evo (J h : ℝ) (q : Cfg N) (t : ℝ) (ψ : Cfg N → ℂ) :
    ∑ n, Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ ψ) n) = ∑ n, Complex.normSq (ψ n) :=
  sum_normSq_mulVec _ (evo_unitary _ (dGamma_isSymm _ (hq_isSymm J h q)) t) ψ

/-- `|2^{-N/2}|² = 2^{-N}`. -/
lemma normSq_prefactor :
    Complex.normSq ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹) = 1 / 2 ^ N := by
  rw [map_inv₀, Complex.normSq_ofReal, ← mul_pow, Real.mul_self_sqrt (by norm_num)]
  ring

/-- A phase has modulus one. -/
lemma normSq_eq_one_of_mul_star {w : ℂ} (hw : w * star w = 1) : Complex.normSq w = 1 := by
  have := congrArg Complex.re hw
  rwa [mul_star_self_eq, Complex.ofReal_re, Complex.one_re] at this

/-- **Dual-frame form of the quench state**: `Φ(t)(n, q) = 2^{-N/2} w_q (e^{-itdΓ(h_q)} ψ)(n)`
with phases `|w_q| = 1` (duality and exact sector decomposition). -/
lemma lgtState_dual (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) :
    ∃ w : Cfg N → ℂ, (∀ q, w q * star (w q) = 1) ∧ ∀ n q,
      lgtState J h ψ t (dualEquiv (n, q)) =
        (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * w q * (evo (dGamma (hq J h q)) t *ᵥ ψ) n := by
  choose w hw hevo using fun q : Cfg N =>
    CutLocalAux.evo_sub_smul_one (dGamma (hq J h q)) (h * ∑ j, sgn (q j)) t
  refine ⟨w, hw, fun n q => ?_⟩
  have hU : blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) =
      (evo (H (N := N) J h) t).submatrix dualEquiv dualEquiv := by
    rw [evo_submatrix, duality, ← evo_blockDiagonal]
    rfl
  have hstate : (evo (H J h) t *ᵥ initZ ψ) ∘ dualEquiv =
      blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) *ᵥ (initZ ψ ∘ dualEquiv) := by
    rw [hU, submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm, Function.comp_id]
  have h1 := congrFun hstate (n, q)
  simp only [Function.comp_apply] at h1
  rw [lgtState, h1, sec_blockDiagonal_mulVec]
  have hvec : (fun j => (initZ ψ ∘ dualEquiv) (j, q)) = (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ • ψ := by
    ext j
    simp [initZ]
  rw [hvec, Matrix.mulVec_smul, Hsec, hevo q, Matrix.smul_mulVec]
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

/-- The dual quench state is normalized. -/
lemma sum_normSq_lgtState_dual (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1)
    (t : ℝ) : ∑ p, Complex.normSq (lgtState J h ψ t (dualEquiv p)) = 1 := by
  obtain ⟨w, hw, hΦ⟩ := lgtState_dual J h ψ t
  rw [Fintype.sum_prod_type, Finset.sum_comm]
  simp only [hΦ, Complex.normSq_mul, normSq_prefactor, normSq_eq_one_of_mul_star (hw _), mul_one,
    ← Finset.mul_sum, sum_normSq_evo, hψ, Finset.sum_const, Finset.card_univ, card_cfg,
    nsmul_eq_mul]
  push_cast
  field_simp

end EntropiesAux

/-- **Duality defect for the quench** (paper: Proposition `prop:entropies`(a)):
`|S_full(t) - S_{A'}(Φ(t))| ≤ ln 2`, with `Φ(t) = U†Ψ(t)`. -/
theorem lgt_duality_defect (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1)
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    |vonNeumann (rdmLGT ℓ hℓ (lgtState J h ψ t)) (rdmLGT_isHermitian ℓ hℓ _) -
      vonNeumann (rdmLGT ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p)))
        (rdmLGT_isHermitian ℓ hℓ _)| ≤ Real.log 2 := by
  have h := cutlocal_defect ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p))
    (EntropiesAux.sum_normSq_lgtState_dual J h ψ hψ t)
  simpa only [Equiv.apply_symm_apply] using h

/-! ### Conditional versus mixed entropy -/

/-- The sector-`q` reduced state `ρ^q_A(t)` of `Γ(e^{-ith_q}) ψ` on the modes `< ℓ`. -/
def secRdm (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) :
    Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ ψ)

/-- `ρ_mix(t) = E_q ρ^q_A(t)`. -/
def mixRdm (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  ((1 / 2 ^ N : ℝ) : ℂ) • ∑ q, secRdm J h ψ t ℓ hℓ q

namespace EntropiesAux

/-- `ρ^q_A(t)` is positive semidefinite. -/
lemma secRdm_posSemidef (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) :
    (secRdm J h ψ t ℓ hℓ q).PosSemidef :=
  rdm_posSemidef ℓ hℓ _

/-- `ρ^q_A(t)` has trace one for normalized `ψ`. -/
lemma trace_secRdm (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) : (secRdm J h ψ t ℓ hℓ q).trace = 1 := by
  rw [secRdm, trace_rdm, sum_normSq_evo, hψ, Complex.ofReal_one]

/-- `ρ_mix(t)` is positive semidefinite. -/
lemma mixRdm_posSemidef (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    (mixRdm J h ψ t ℓ hℓ).PosSemidef := by
  unfold mixRdm
  refine Matrix.PosSemidef.smul
    (posSemidef_sum _ fun q _ => secRdm_posSemidef J h ψ t ℓ hℓ q) ?_
  exact_mod_cast (by positivity : (0 : ℝ) ≤ 1 / 2 ^ N)

/-- `ρ_mix(t)` has trace one for normalized `ψ`. -/
lemma trace_mixRdm (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) : (mixRdm J h ψ t ℓ hℓ).trace = 1 := by
  rw [mixRdm, Matrix.trace_smul, Matrix.trace_sum]
  simp only [trace_secRdm J h ψ hψ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul,
    smul_eq_mul]
  push_cast
  field_simp

end EntropiesAux

/-- `ρ_mix(t)` is Hermitian. -/
theorem mixRdm_isHermitian (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    (mixRdm J h ψ t ℓ hℓ).IsHermitian :=
  (EntropiesAux.mixRdm_posSemidef J h ψ t ℓ hℓ).isHermitian

/-- `S_G(t) = E_q S(ρ^q_A(t))`, the entanglement conditional on all Gauss charges. -/
def condEntropy (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : ℝ :=
  avgQ (fun q => vonNeumann (secRdm J h ψ t ℓ hℓ q) (rdm_isHermitian ℓ hℓ _))

open Classical in
/-- `H_t`: the Shannon entropy of the law of the random density matrix `q ↦ ρ^q_A(t)` under the
uniform law of `q`. -/
def sectorLawEntropy (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : ℝ :=
  ∑ σ ∈ univ.image (secRdm J h ψ t ℓ hℓ),
    Real.negMulLog (((univ.filter (fun q => secRdm J h ψ t ℓ hℓ q = σ)).card : ℝ) / 2 ^ N)

/-- `ρ_mix(t)` is the marginal on the modes `< ℓ` of the reduced state `ρ_{A'}(Φ(t))` of the dual
state `Φ(t) = U†Ψ(t)` (paper: Proposition `prop:entropies`(b)). -/
theorem mixRdm_eq_marginal (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (a a' : Cfg ℓ) :
    mixRdm J h ψ t ℓ hℓ a a' =
      ∑ c : Cfg ℓ, rdmLGT ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p)) (a, c) (a', c) := by
  obtain ⟨w, hw, hΦ⟩ := EntropiesAux.lgtState_dual J h ψ t
  set c0 : ℂ := (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ with hc0def
  set φ : Cfg N → Cfg N → ℂ := fun q => evo (dGamma (hq J h q)) t *ᵥ ψ with hφ
  have hc0 : c0 * star c0 = ((1 / 2 ^ N : ℝ) : ℂ) := by
    rw [EntropiesAux.mul_star_self_eq, hc0def, EntropiesAux.normSq_prefactor]
  have hterm : ∀ (q : Cfg N) (x y : Cfg N),
      lgtState J h ψ t (dualEquiv (x, q)) * star (lgtState J h ψ t (dualEquiv (y, q))) =
        ((1 / 2 ^ N : ℝ) : ℂ) * (φ q x * star (φ q y)) := by
    intro q x y
    rw [hΦ, hΦ, star_mul', star_mul',
      show c0 * w q * φ q x * (star c0 * star (w q) * star (φ q y)) =
        (c0 * star c0) * (w q * star (w q)) * (φ q x * star (φ q y)) by ring, hc0, hw q, mul_one]
  have hL : mixRdm J h ψ t ℓ hℓ a a' = ∑ q : Cfg N, ∑ γ : Cfg (N - ℓ),
      ((1 / 2 ^ N : ℝ) : ℂ) * (φ q ((splitCfg ℓ hℓ).symm (a, γ)) *
        star (φ q ((splitCfg ℓ hℓ).symm (a', γ)))) := by
    simp only [mixRdm, Matrix.smul_apply, Matrix.sum_apply, secRdm, rdm, smul_eq_mul,
      Finset.mul_sum, hφ]
  have hR : ∑ c : Cfg ℓ, rdmLGT ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p)) (a, c) (a', c) =
      ∑ q : Cfg N, ∑ γ : Cfg (N - ℓ), ((1 / 2 ^ N : ℝ) : ℂ) *
        (φ q ((splitCfg ℓ hℓ).symm (a, γ)) * star (φ q ((splitCfg ℓ hℓ).symm (a', γ)))) := by
    rw [← Equiv.sum_comp (splitCfg ℓ hℓ).symm, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun c _ => ?_
    simp only [rdmLGT, Fintype.sum_prod_type, CutLocalAux.splitLGT_symm_apply]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun δ _ => Finset.sum_congr rfl fun γ _ => ?_
    exact hterm _ _ _
  rw [hL, hR]

/-- `S_G ≤ S_mix` (paper: Proposition `prop:entropies`(b)). -/
theorem condEntropy_le_mix (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    condEntropy J h ψ t ℓ hℓ ≤ vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) := by
  have hw1 : ∑ _q : Cfg N, (1 / 2 ^ N : ℝ) = 1 := by
    rw [Finset.sum_const, Finset.card_univ, EntropiesAux.card_cfg, nsmul_eq_mul]
    push_cast
    field_simp
  have hmix : (∑ q : Cfg N, (((1 / 2 ^ N : ℝ)) : ℂ) • secRdm J h ψ t ℓ hℓ q) =
      mixRdm J h ψ t ℓ hℓ := by
    rw [mixRdm, Finset.smul_sum]
  have hH : (∑ q : Cfg N, (((1 / 2 ^ N : ℝ)) : ℂ) • secRdm J h ψ t ℓ hℓ q).IsHermitian :=
    hmix ▸ mixRdm_isHermitian J h ψ t ℓ hℓ
  have hmain := vonNeumann_concave (fun _ : Cfg N => (1 / 2 ^ N : ℝ)) (fun _ => by positivity) hw1
    (secRdm J h ψ t ℓ hℓ) (fun q => EntropiesAux.secRdm_posSemidef J h ψ t ℓ hℓ q) hH
  rw [CutLocalAux.vonNeumann_congr hH (mixRdm_isHermitian J h ψ t ℓ hℓ) hmix] at hmain
  unfold condEntropy avgQ
  rw [Finset.mul_sum]
  exact hmain

/-- `S_mix ≤ S_G + H_t` (paper: Proposition `prop:entropies`(b)). -/
theorem mix_le_cond_add (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) ≤
      condEntropy J h ψ t ℓ hℓ + sectorLawEntropy J h ψ t ℓ hℓ := by
  set f := secRdm J h ψ t ℓ hℓ with hf
  set T := univ.image f with hT
  set μ : Matrix (Cfg ℓ) (Cfg ℓ) ℂ → ℝ :=
    fun σ => ((univ.filter (fun q => f q = σ)).card : ℝ) / 2 ^ N with hμ
  -- the entropy as a total function on matrices
  set S' : Matrix (Cfg ℓ) (Cfg ℓ) ℂ → ℝ :=
    fun σ => if hσ : σ.IsHermitian then vonNeumann σ hσ else 0 with hS'
  have hS'apply : ∀ σ (hσ : σ.IsHermitian), S' σ = vonNeumann σ hσ := fun σ hσ => dite_eq_left hσ
  have hmemT : ∀ σ : T, ∃ q, f q = σ.1 := fun σ => by
    obtain ⟨q, -, hq⟩ := Finset.mem_image.1 σ.2
    exact ⟨q, hq⟩
  have hσ : ∀ σ : T, (σ.1).PosSemidef := fun σ => by
    obtain ⟨q, hq⟩ := hmemT σ
    rw [← hq]
    exact EntropiesAux.secRdm_posSemidef J h ψ t ℓ hℓ q
  have htr : ∀ σ : T, (σ.1).trace = 1 := fun σ => by
    obtain ⟨q, hq⟩ := hmemT σ
    rw [← hq]
    exact EntropiesAux.trace_secRdm J h ψ hψ t ℓ hℓ q
  have hcomp : ∑ q, f q = ∑ σ ∈ T, (univ.filter (fun q => f q = σ)).card • σ :=
    Finset.sum_comp (fun σ => σ) f
  have hsum : ∑ σ : T, (μ σ.1 : ℂ) • σ.1 = mixRdm J h ψ t ℓ hℓ := by
    rw [Finset.sum_coe_sort T (fun σ => (μ σ : ℂ) • σ), mixRdm, ← hf, hcomp, Finset.smul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    rw [← Nat.cast_smul_eq_nsmul ℂ, smul_smul]
    congr 1
    simp only [hμ]
    push_cast
    ring
  have hH : (∑ σ : T, (μ σ.1 : ℂ) • σ.1).IsHermitian := hsum ▸ mixRdm_isHermitian J h ψ t ℓ hℓ
  have hmain := vonNeumann_mix_le (fun σ : T => μ σ.1) (fun σ => by positivity) (fun σ : T => σ.1)
    hσ htr hH
  rw [CutLocalAux.vonNeumann_congr hH (mixRdm_isHermitian J h ψ t ℓ hℓ) hsum] at hmain
  refine hmain.trans (le_of_eq ?_)
  congr 1
  · have h1 : ∀ σ : T, μ σ.1 * vonNeumann σ.1 (hσ σ).isHermitian = μ σ.1 * S' σ.1 := fun σ => by
      rw [hS'apply σ.1 (hσ σ).isHermitian]
    rw [Finset.sum_congr rfl fun σ _ => h1 σ, Finset.sum_coe_sort T (fun σ => μ σ * S' σ)]
    unfold condEntropy avgQ
    have h2 : ∀ q, vonNeumann (secRdm J h ψ t ℓ hℓ q) (rdm_isHermitian ℓ hℓ _) = S' (f q) :=
      fun q => (hS'apply (f q) (rdm_isHermitian ℓ hℓ _)).symm
    simp only [h2]
    rw [Finset.sum_comp (s := univ) S' f, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp only [hμ, nsmul_eq_mul]
    ring
  · rw [Finset.sum_coe_sort T (fun σ => Real.negMulLog (μ σ))]
    rfl

/-- `S_mix ≤ ℓ ln 2` (paper: Proposition `prop:entropies`(b)). -/
theorem mix_le_log (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) ≤ ℓ * Real.log 2 := by
  have h := vonNeumann_le_log_card (mixRdm J h ψ t ℓ hℓ)
    (EntropiesAux.mixRdm_posSemidef J h ψ t ℓ hℓ) (EntropiesAux.trace_mixRdm J h ψ hψ t ℓ hℓ)
  rw [EntropiesAux.card_cfg, Nat.cast_pow, Real.log_pow] at h
  push_cast at h
  exact h

namespace EntropiesAux

/-- `η(m / 2^N) ≤ (m / 2^N) N ln 2` for natural `m`. -/
lemma negMulLog_nat_div_le (m : ℕ) :
    Real.negMulLog ((m : ℝ) / 2 ^ N) ≤ (m : ℝ) / 2 ^ N * (N * Real.log 2) := by
  rcases Nat.eq_zero_or_pos m with hm | hm
  · subst hm
    simp
  · have hm1 : (1 : ℝ) ≤ m := by exact_mod_cast hm
    have hlogm : 0 ≤ Real.log m := Real.log_nonneg hm1
    have hx : 0 ≤ (m : ℝ) / 2 ^ N := by positivity
    rw [Real.negMulLog, Real.log_div (by positivity) (by positivity), Real.log_pow]
    nlinarith

/-- At `t = 0` the propagators are the identity. -/
lemma evo_zero {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) : evo A 0 = 1 := by
  unfold evo
  simp

end EntropiesAux

/-- `H_t ≤ N ln 2` (paper: Proposition `prop:entropies`(b)). -/
theorem sectorLawEntropy_le (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    sectorLawEntropy J h ψ t ℓ hℓ ≤ N * Real.log 2 := by
  unfold sectorLawEntropy
  set f := secRdm J h ψ t ℓ hℓ with hf
  have hsum : ∑ σ ∈ univ.image f, ((univ.filter (fun q => f q = σ)).card : ℝ) = 2 ^ N := by
    rw [← Nat.cast_sum, ← Finset.card_eq_sum_card_image f univ, Finset.card_univ,
      EntropiesAux.card_cfg]
    push_cast
    ring
  calc ∑ σ ∈ univ.image f, Real.negMulLog (((univ.filter (fun q => f q = σ)).card : ℝ) / 2 ^ N)
      ≤ ∑ σ ∈ univ.image f, ((univ.filter (fun q => f q = σ)).card : ℝ) / 2 ^ N *
          (N * Real.log 2) :=
        Finset.sum_le_sum fun σ _ => EntropiesAux.negMulLog_nat_div_le _
    _ = N * Real.log 2 := by
        rw [← Finset.sum_mul, ← Finset.sum_div, hsum, div_self (by positivity), one_mul]

/-- `H_0 = 0` (paper: Proposition `prop:entropies`(b)). -/
theorem sectorLawEntropy_zero (J h : ℝ) (ψ : Cfg N → ℂ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    sectorLawEntropy J h ψ 0 ℓ hℓ = 0 := by
  have hf : secRdm J h ψ 0 ℓ hℓ = fun _ => rdm ℓ hℓ ψ := by
    funext q
    rw [secRdm, EntropiesAux.evo_zero, Matrix.one_mulVec]
  unfold sectorLawEntropy
  rw [hf, Finset.image_const Finset.univ_nonempty, Finset.sum_singleton]
  simp

/-! ### Number entropy as a lower bound -/

/-- `p_t = E_q p^q_t`: the LGT number distribution is the average of the sector distributions
(paper: Corollary `cor:sector`(c), used in Proposition `prop:entropies`(c)). -/
theorem lgtNumDist_eq_avg (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) :
    lgtNumDist J h R ψ t k = avgQ (fun q => numDist R (evo (dGamma (hq J h q)) t *ᵥ ψ) k) := by
  rw [lgtNumDist, matter, numProj, sector_average_matter]
  have h2 : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [h2, Complex.re_ofReal_mul, Complex.re_sum]
  unfold avgQ
  congr 1
  refine Finset.sum_congr rfl fun q _ => ?_
  dsimp only
  rw [TaxonomyAux.expect_diag_re, numDist, Finset.sum_filter]
  refine Finset.sum_congr rfl fun m _ => ?_
  split_ifs <;> ring

namespace EntropiesAux

open TaxonomyAux

/-- The number distribution of the segment, read off the reduced density matrix. -/
lemma numDist_segment_eq (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) (k : ℕ) :
    numDist (segment ℓ) v k =
      (∑ α ∈ univ.filter (fun α : Cfg ℓ => cnt α = k), rdm ℓ hℓ v α α).re := by
  rw [numDist, Finset.sum_filter, ← (splitCfg ℓ hℓ).symm.sum_comp, Fintype.sum_prod_type,
    Complex.re_sum, Finset.sum_filter]
  refine Finset.sum_congr rfl fun α _ => ?_
  simp only [count_segment_splitCfg_symm, ← cnt_cast, Nat.cast_inj]
  split_ifs with hk
  · simp only [rdm, Complex.re_sum, re_mul_star_self]
  · simp

/-- The reduced density matrix of an evolved Fock state is block diagonal in the particle number
of the segment. -/
lemma rdm_evo_fock_block (A : Matrix (Fin N) (Fin N) ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ)
    (hℓ : ℓ ≤ N) (α α' : Cfg ℓ) (hne : cnt α ≠ cnt α') :
    rdm ℓ hℓ (evo (dGamma A) t *ᵥ Pi.single n 1) α α' = 0 := by
  refine Finset.sum_eq_zero fun β _ => ?_
  by_cases h1 : (evo (dGamma A) t *ᵥ Pi.single n 1) ((splitCfg ℓ hℓ).symm (α, β)) = 0
  · rw [h1, zero_mul]
  by_cases h2 : (evo (dGamma A) t *ᵥ Pi.single n 1) ((splitCfg ℓ hℓ).symm (α', β)) = 0
  · rw [h2, star_zero, mul_zero]
  exfalso
  apply hne
  have e1 : count univ ((splitCfg ℓ hℓ).symm (α, β)) = count univ n := by
    by_contra hc
    exact h1 (evo_fock_eq_zero A t n _ hc)
  have e2 : count univ ((splitCfg ℓ hℓ).symm (α', β)) = count univ n := by
    by_contra hc
    exact h2 (evo_fock_eq_zero A t n _ hc)
  rw [count_univ_splitCfg_symm] at e1 e2
  have : (cnt α : ℝ) = cnt α' := by rw [cnt_cast, cnt_cast]; linarith
  exact_mod_cast this

/-- The reduced LGT density matrix of the quench from a Fock state is block diagonal in the matter
particle number of `A`. -/
lemma rdmLGT_lgtState_block (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (a a' : Cfg ℓ × Cfg ℓ) (hne : cnt a.1 ≠ cnt a'.1) :
    rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t) a a' = 0 := by
  refine Finset.sum_eq_zero fun b _ => ?_
  by_cases h1 : lgtState J h (Pi.single n 1) t ((splitLGT ℓ hℓ).symm (a, b)) = 0
  · rw [h1, zero_mul]
  by_cases h2 : lgtState J h (Pi.single n 1) t ((splitLGT ℓ hℓ).symm (a', b)) = 0
  · rw [h2, star_zero, mul_zero]
  exfalso
  apply hne
  have e1 : count univ ((splitLGT ℓ hℓ).symm (a, b)).1 = count univ n := by
    by_contra hc
    exact h1 (lgtState_eq_zero J h n t _ hc)
  have e2 : count univ ((splitLGT ℓ hℓ).symm (a', b)).1 = count univ n := by
    by_contra hc
    exact h2 (lgtState_eq_zero J h n t _ hc)
  rw [splitLGT_symm_fst, count_univ_splitCfg_symm] at e1 e2
  have : (cnt a.1 : ℝ) = cnt a'.1 := by rw [cnt_cast, cnt_cast]; linarith
  exact_mod_cast this

/-- The LGT number distribution, read off the reduced LGT density matrix. -/
lemma lgtNumDist_eq_rdmLGT (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (k : ℕ) :
    lgtNumDist J h (segment ℓ) ψ t k =
      (∑ a ∈ univ.filter (fun a : Cfg ℓ × Cfg ℓ => cnt a.1 = k),
        rdmLGT ℓ hℓ (lgtState J h ψ t) a a).re := by
  simp only [lgtNumDist, matter, numProj]
  rw [sec_diagonal_kronecker_one, expect_diag_re, ← (splitLGT ℓ hℓ).symm.sum_comp,
    Fintype.sum_prod_type, Complex.re_sum, Finset.sum_filter]
  refine Finset.sum_congr rfl fun a _ => ?_
  simp only [splitLGT_symm_fst, count_segment_splitCfg_symm, ← cnt_cast, Nat.cast_inj]
  split_ifs with hk
  · simp only [one_mul, rdmLGT, Complex.re_sum, re_mul_star_self]
  · simp

/-- A basis vector is normalized. -/
lemma sum_normSq_single (n : Cfg N) :
    ∑ m, Complex.normSq ((Pi.single n (1 : ℂ) : Cfg N → ℂ) m) = 1 := by
  rw [Finset.sum_eq_single n (fun m _ hm => by rw [Pi.single_eq_of_ne hm, map_zero])
    (fun h => absurd (mem_univ _) h), Pi.single_eq_same, map_one]

end EntropiesAux

open EntropiesAux TaxonomyAux in
/-- `S_N ≤ S_mix` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_mix (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (mixRdm J h (Pi.single n 1) t ℓ hℓ)
        (mixRdm_isHermitian J h (Pi.single n 1) t ℓ hℓ) := by
  have hblock : ∀ α α' : Cfg ℓ, cnt α ≠ cnt α' → mixRdm J h (Pi.single n 1) t ℓ hℓ α α' = 0 := by
    intro α α' hne
    simp only [mixRdm, Matrix.smul_apply, Matrix.sum_apply, secRdm,
      rdm_evo_fock_block _ n t ℓ hℓ α α' hne, Finset.sum_const_zero, smul_zero]
  have hmain := block_bound _ (mixRdm_posSemidef J h (Pi.single n 1) t ℓ hℓ) cnt N
    (fun α => (cnt_le α).trans hℓ) hblock
  have hdist : lgtNumDist J h (segment ℓ) (Pi.single n 1) t = fun k =>
      (∑ α ∈ univ.filter (fun α : Cfg ℓ => cnt α = k),
        mixRdm J h (Pi.single n 1) t ℓ hℓ α α).re := by
    funext k
    rw [lgtNumDist_eq_avg]
    unfold avgQ
    simp only [numDist_segment_eq ℓ hℓ, mixRdm, Matrix.smul_apply, Matrix.sum_apply, smul_eq_mul,
      ← Finset.mul_sum]
    rw [Finset.sum_comm, Complex.re_ofReal_mul, Complex.re_sum]
    rfl
  rw [hdist]
  exact hmain

/-- The reduced state `ρ^f_A` of the bare fermion modes `< ℓ`, with all links traced out. -/
def rdmBare (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  fun a a' => ∑ c : Cfg ℓ, rdmLGT ℓ hℓ Ψ (a, c) (a', c)

/-- `ρ^f_A` is Hermitian. -/
theorem rdmBare_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmBare ℓ hℓ Ψ).IsHermitian := by
  ext a a'
  simp only [Matrix.conjTranspose_apply, rdmBare, star_sum]
  refine Finset.sum_congr rfl fun c _ => ?_
  exact congrFun (congrFun (rdmLGT_isHermitian ℓ hℓ Ψ) (a, c)) (a', c)

namespace EntropiesAux

/-- `ρ^f_A = X X†` is positive semidefinite. -/
lemma rdmBare_posSemidef (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmBare ℓ hℓ Ψ).PosSemidef := by
  have h : rdmBare ℓ hℓ Ψ =
      (Matrix.of fun (a : Cfg ℓ) (b : Cfg ℓ × (Cfg (N - ℓ) × Cfg (N - ℓ))) =>
        Ψ ((splitLGT ℓ hℓ).symm ((a, b.1), b.2))) *
      (Matrix.of fun (a : Cfg ℓ) (b : Cfg ℓ × (Cfg (N - ℓ) × Cfg (N - ℓ))) =>
        Ψ ((splitLGT ℓ hℓ).symm ((a, b.1), b.2)))ᴴ := by
    ext a a'
    simp [rdmBare, rdmLGT, Matrix.mul_apply, Matrix.conjTranspose_apply, Fintype.sum_prod_type]
  rw [h]
  exact Matrix.posSemidef_self_mul_conjTranspose _

end EntropiesAux

open EntropiesAux TaxonomyAux in
/-- `S_N ≤ S_f` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_bare (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t))
        (rdmBare_isHermitian ℓ hℓ _) := by
  set Ψ := lgtState J h (Pi.single n 1) t with hΨ
  have hblock : ∀ a a' : Cfg ℓ, cnt a ≠ cnt a' → rdmBare ℓ hℓ Ψ a a' = 0 := fun a a' hne =>
    Finset.sum_eq_zero fun c _ => rdmLGT_lgtState_block J h n t ℓ hℓ (a, c) (a', c) hne
  have hmain := block_bound _ (rdmBare_posSemidef ℓ hℓ Ψ) cnt N (fun a => (cnt_le a).trans hℓ)
    hblock
  have hdist : lgtNumDist J h (segment ℓ) (Pi.single n 1) t = fun k =>
      (∑ a ∈ univ.filter (fun a : Cfg ℓ => cnt a = k), rdmBare ℓ hℓ Ψ a a).re := by
    funext k
    rw [lgtNumDist_eq_rdmLGT J h _ t ℓ hℓ k, ← hΨ]
    congr 1
    simp only [rdmBare]
    rw [Finset.sum_filter, Finset.sum_filter, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun a _ => ?_
    split_ifs <;> simp
  rw [hdist]
  exact hmain

/-- `S_G ≥ E_q H(p^q_t)` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem avg_sector_shannon_le_cond (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) ≤
      condEntropy J h (Pi.single n 1) t ℓ hℓ := by
  unfold condEntropy avgQ
  refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => ?_) (by positivity)
  exact numDist_shannon_le_vonNeumann_slater (hq J h q) (hq_isSymm J h q) n t ℓ hℓ

set_option linter.unusedVariables false in
/-- `E_q H(p^q_t) ≤ S_N(t)` (concavity; paper: Proposition `prop:entropies`(c)). -/
theorem avg_sector_shannon_le_numberEntropy (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) ≤
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) := by
  have hw1 : ∑ _q : Cfg N, (1 / 2 ^ N : ℝ) = 1 := by
    rw [Finset.sum_const, Finset.card_univ, EntropiesAux.card_cfg, nsmul_eq_mul]
    push_cast
    field_simp
  unfold shannon
  simp only [lgtNumDist_eq_avg]
  unfold avgQ
  simp only [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_le_sum fun k _ => ?_
  have hj := Real.concaveOn_negMulLog.le_map_sum (t := univ) (w := fun _ : Cfg N => 1 / 2 ^ N)
    (p := fun q => numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1) k)
    (fun _ _ => by positivity) hw1
    (fun q _ => Set.mem_Ici.2 (Finset.sum_nonneg fun _ _ => Complex.normSq_nonneg _))
  simpa only [smul_eq_mul] using hj

namespace EntropiesAux

/-- Subadditivity of `η(x) = -x log x` on nonnegative reals. -/
lemma negMulLog_sum_le {ι : Type*} (s : Finset ι) (x : ι → ℝ) (hx : ∀ i ∈ s, 0 ≤ x i) :
    Real.negMulLog (∑ i ∈ s, x i) ≤ ∑ i ∈ s, Real.negMulLog (x i) := by
  rw [Real.negMulLog, neg_mul, Finset.sum_mul, ← Finset.sum_neg_distrib]
  refine Finset.sum_le_sum fun i hi => ?_
  rw [Real.negMulLog]
  rcases (hx i hi).eq_or_lt with h0 | hpos
  · rw [← h0]
    simp
  · have hle : x i ≤ ∑ j ∈ s, x j := Finset.single_le_sum hx hi
    have hlog := Real.log_le_log hpos hle
    have := mul_le_mul_of_nonneg_left hlog hpos.le
    linarith

end EntropiesAux

open EntropiesAux TaxonomyAux in
/-- `S_N(t) ≤ E_q H(p^q_t) + H_t` (classical mixing bound, sectors grouped by the value of
`ρ^q_A(t)`, which determines `p^q_t`; paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_avg_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) +
        sectorLawEntropy J h (Pi.single n 1) t ℓ hℓ := by
  set f := secRdm J h (Pi.single n 1) t ℓ hℓ with hf
  set T := univ.image f with hT
  -- the number distribution of the segment as a function of the reduced density matrix
  set g : Matrix (Cfg ℓ) (Cfg ℓ) ℂ → ℕ → ℝ :=
    fun σ k => (∑ α ∈ univ.filter (fun α : Cfg ℓ => cnt α = k), σ α α).re with hg
  -- the law of `q ↦ ρ^q_A(t)`
  set μ : Matrix (Cfg ℓ) (Cfg ℓ) ℂ → ℝ :=
    fun σ => ((univ.filter (fun q => f q = σ)).card : ℝ) / 2 ^ N with hμ
  have hp' : ∀ q k, numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1) k =
      g (f q) k := fun q k => numDist_segment_eq ℓ hℓ _ k
  have hp : ∀ q, numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1) = g (f q) :=
    fun q => funext fun k => hp' q k
  have hμ0 : ∀ σ, 0 ≤ μ σ := fun σ => by positivity
  have hg0 : ∀ σ ∈ T, ∀ k, 0 ≤ g σ k := by
    intro σ hσ k
    obtain ⟨q, -, rfl⟩ := Finset.mem_image.1 hσ
    rw [← hp' q k]
    exact Finset.sum_nonneg fun _ _ => Complex.normSq_nonneg _
  have hg1 : ∀ σ ∈ T, ∑ k ∈ range (N + 1), g σ k = 1 := by
    intro σ hσ
    obtain ⟨q, -, rfl⟩ := Finset.mem_image.1 hσ
    have htr : ∑ α, f q α α = 1 := trace_secRdm J h (Pi.single n 1) (sum_normSq_single n) t ℓ hℓ q
    simp only [hg]
    rw [← Complex.re_sum, Finset.sum_fiberwise_of_maps_to
      (fun α _ => Finset.mem_range.2 (Nat.lt_succ_of_le ((cnt_le α).trans hℓ))), htr,
      Complex.one_re]
  -- `p_t = ∑_σ μ_σ g_σ`
  have hmix : ∀ k, lgtNumDist J h (segment ℓ) (Pi.single n 1) t k = ∑ σ ∈ T, μ σ * g σ k := by
    intro k
    rw [lgtNumDist_eq_avg]
    unfold avgQ
    simp only [hp']
    rw [Finset.sum_comp (s := univ) (fun σ => g σ k) f, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp only [hμ, nsmul_eq_mul]
    ring
  -- `E_q H(p^q_t) = ∑_σ μ_σ H(g_σ)`
  have havg : avgQ (fun q => shannon N (numDist (segment ℓ)
      (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) = ∑ σ ∈ T, μ σ * shannon N (g σ) := by
    unfold avgQ
    simp only [hp]
    rw [Finset.sum_comp (s := univ) (fun σ => shannon N (g σ)) f, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp only [hμ, nsmul_eq_mul]
    ring
  have hH : sectorLawEntropy J h (Pi.single n 1) t ℓ hℓ = ∑ σ ∈ T, Real.negMulLog (μ σ) := rfl
  rw [havg, hH]
  unfold shannon
  simp only [hmix]
  calc ∑ k ∈ range (N + 1), Real.negMulLog (∑ σ ∈ T, μ σ * g σ k)
      ≤ ∑ k ∈ range (N + 1), ∑ σ ∈ T, Real.negMulLog (μ σ * g σ k) :=
        Finset.sum_le_sum fun k _ =>
          negMulLog_sum_le _ _ fun σ hσ => mul_nonneg (hμ0 σ) (hg0 σ hσ k)
    _ = ∑ σ ∈ T, (μ σ * ∑ k ∈ range (N + 1), Real.negMulLog (g σ k) +
          Real.negMulLog (μ σ) * ∑ k ∈ range (N + 1), g σ k) := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun σ _ => ?_
        rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [Real.negMulLog_mul]
        ring
    _ = ∑ σ ∈ T, μ σ * ∑ k ∈ range (N + 1), Real.negMulLog (g σ k) +
          ∑ σ ∈ T, Real.negMulLog (μ σ) := by
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun σ hσ => ?_
        rw [hg1 σ hσ, mul_one]

end DFL
