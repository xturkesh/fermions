import DFL.Entropies

/-!
# Interior intervals (paper: Remark `rem:interior`)

For an interval `I = [ℓ₁, ℓ₁ + m)` (sites and links; Lean indexing, `m ≥ 1`) of the open chain, the
duality `U`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, changes the entanglement entropy of `I` by at most
`2 ln 2`: `|S_I(UΦ) - S_{I'}(Φ)| ≤ 2 ln 2` (`interior_defect`), with the partial trace in the
occupation basis (the Jordan–Wigner partial trace) on both sides. In configuration language:
* the right gate (control: the last link of `I`, target: the links right of `I`) is removed by a
  substitution in the sum over the right part, up to the pinching by the label of the last link of
  `I`; this costs at most `ln 2` (as in `cutlocal_defect`);
* the left gate (control: the last link left of `I`, target: all links right of it) splits the
  reduced state as `ω₀ + F ω₁ F†`, where `ω_c` is the part with control label `c` and `F` flips all
  link labels of `I`; concavity and the mixing bound give a change of at most `ln 2`.
The value `2 ln 2` is attained for `N = 3`, `I` = site `1` and link `1` (paper: site 2 and link `b_2`):
`exInt` is a product state across `I'` whose image under `U` has maximally mixed reduced state on
`I` (`exInt_entropy`, `exInt_U_entropy`). The value `2 ln 2` is attained already by the matter
vacuum for `N = 4` and a two-site interval (`exVac`, `interior_defect_attained_vacuum`). For states
of definite total fermion parity the reduced state of `I` is block diagonal with respect to the
parity of `I` (`rdmMid_parity_block`).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Split an LGT configuration `(n, x)` into the interval `I = [ℓ₁, ℓ₁ + m)` (sites and links), the
left part `[0, ℓ₁)` and the right part `[ℓ₁ + m, N)`. -/
def splitMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) :
    Cfg N × Cfg N ≃
      (Cfg m × Cfg m) × ((Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) :=
  (splitLGT ℓ₁ (by omega)).trans <|
    (Equiv.prodCongr (Equiv.refl _) (splitLGT (N := N - ℓ₁) m (by omega))).trans
      { toFun := fun p => (p.2.1, (p.1, p.2.2))
        invFun := fun p => (p.2.1, (p.1, p.2.2))
        left_inv := fun _ => rfl
        right_inv := fun _ => rfl }

/-- The reduced density matrix of the interval `[ℓ₁, ℓ₁ + m)`: partial trace over the left and the
right part in the occupation basis (the Jordan–Wigner partial trace). -/
def rdmMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    Matrix (Cfg m × Cfg m) (Cfg m × Cfg m) ℂ :=
  fun a a' => ∑ b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) * star (Ψ ((splitMid ℓ₁ m h).symm (a', b)))

namespace InteriorAux

open EntropiesAux

/-! ### Splitting into three parts -/

/-- Recombining a configuration split by `splitMid`: first the interval and the right part, then
the left part. -/
lemma splitMid_symm_apply (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (a : Cfg m × Cfg m)
    (l : Cfg ℓ₁ × Cfg ℓ₁) (r : Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m)) :
    (splitMid ℓ₁ m h).symm (a, (l, r)) =
      (splitLGT ℓ₁ (by omega)).symm (l, (splitLGT (N := N - ℓ₁) m (by omega)).symm (a, r)) :=
  rfl

/-- The coefficient matrix `X_{a b} = Ψ(a, b)` of a state split by `splitMid`. -/
def coefMat (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    Matrix (Cfg m × Cfg m) ((Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) ℂ :=
  Matrix.of fun a b => Ψ ((splitMid ℓ₁ m h).symm (a, b))

/-- The reduced density matrix of the interval is `X X†`. -/
lemma rdmMid_eq_mul (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    rdmMid ℓ₁ m h Ψ = coefMat ℓ₁ m h Ψ * (coefMat ℓ₁ m h Ψ)ᴴ := by
  ext a a'
  simp [rdmMid, coefMat, Matrix.mul_apply, Matrix.conjTranspose_apply]

/-- The sum of the labels of a recombined configuration is the sum over the two parts. -/
lemma sum_splitCfg_symm (ℓ : ℕ) (hℓ : ℓ ≤ N) (x : Cfg ℓ) (y : Cfg (N - ℓ)) :
    ∑ k, (splitCfg ℓ hℓ).symm (x, y) k = ∑ k, x k + ∑ k, y k := by
  rw [AreaLawAux.sum_fin_split ℓ hℓ]
  congr 1
  · refine Finset.sum_congr rfl fun i _ => ?_
    rw [← AreaLawAux.splitCfg_fst ℓ hℓ ((splitCfg ℓ hℓ).symm (x, y)) i, Equiv.apply_symm_apply]
  · refine Finset.sum_congr rfl fun i _ => ?_
    rw [← AreaLawAux.splitCfg_snd ℓ hℓ ((splitCfg ℓ hℓ).symm (x, y)) i, Equiv.apply_symm_apply]

/-- The fermion parity of a configuration recombined by `splitMid` is the sum of the parities of
the three parts. -/
lemma sum_splitMid_symm (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (a : Cfg m × Cfg m)
    (b : (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) :
    ∑ k, ((splitMid ℓ₁ m h).symm (a, b)).1 k = ∑ k, a.1 k + (∑ k, b.1.1 k + ∑ k, b.2.1 k) := by
  obtain ⟨l, r⟩ := b
  rw [splitMid_symm_apply]
  change ∑ k, (splitCfg ℓ₁ (by omega)).symm
    (l.1, (splitCfg (N := N - ℓ₁) m (by omega)).symm (a.1, r.1)) k = _
  rw [sum_splitCfg_symm, sum_splitCfg_symm]
  ring

end InteriorAux

/-- The reduced density matrix of an interval is Hermitian. -/
theorem rdmMid_isHermitian (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).IsHermitian := by
  ext a a'
  simp only [Matrix.conjTranspose_apply, rdmMid, star_sum, star_mul', star_star]
  exact Finset.sum_congr rfl fun b _ => mul_comm _ _

/-- The reduced density matrix of an interval is positive semidefinite. -/
theorem rdmMid_posSemidef (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).PosSemidef := by
  rw [InteriorAux.rdmMid_eq_mul]
  exact Matrix.posSemidef_self_mul_conjTranspose _

/-- The trace of the reduced density matrix is the squared norm of the state. -/
theorem trace_rdmMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).trace = ((∑ p, Complex.normSq (Ψ p) : ℝ) : ℂ) := by
  simp only [Matrix.trace, Matrix.diag, rdmMid, EntropiesAux.mul_star_self_eq,
    Complex.ofReal_sum]
  rw [← Fintype.sum_prod_type'
    (fun a b => (Complex.normSq (Ψ ((splitMid ℓ₁ m h).symm (a, b))) : ℂ))]
  exact Equiv.sum_comp (splitMid ℓ₁ m h).symm (fun p => (Complex.normSq (Ψ p) : ℂ))

namespace InteriorAux

open EntropiesAux

/-! ### The duality on a split configuration -/

/-- Adding a constant to all labels of a recombined configuration. -/
lemma splitCfg_symm_add_const (ℓ : ℕ) (hℓ : ℓ ≤ N) (β : Cfg ℓ) (δ : Cfg (N - ℓ)) (c : ZMod 2) :
    (splitCfg ℓ hℓ).symm (β, δ) + (fun _ => c) =
      (splitCfg ℓ hℓ).symm (β + fun _ => c, δ + fun _ => c) := by
  ext k
  by_cases hk : k.val < ℓ
  · simp only [Pi.add_apply, CutLocalAux.splitCfg_symm_apply_lt _ _ _ _ _ hk]
  · simp only [Pi.add_apply, CutLocalAux.splitCfg_symm_apply_ge _ _ _ _ _ (not_lt.1 hk)]

/-- Shifting all labels of a nonempty chain shifts the label of its last link. -/
lemma cutLabel_add_const {K : ℕ} (hK : 0 < K) (β : Cfg K) (c : ZMod 2) :
    cutLabel (β + fun _ => c) = cutLabel β + c := by
  simp [cutLabel, hK]

/-- Adding the constant `0` does nothing. -/
lemma add_const_zero {K : ℕ} (β : Cfg K) : (β + fun _ => (0 : ZMod 2)) = β := by
  funext k
  simp

/-- **Charges of a configuration recombined by `splitMid`** (`dualEquiv_symm_splitLGT` at `ℓ₁`
and then at `m` on the remainder): with interval part `(αI, βI)`, left part `(αL, βL)` and right
part `(γR, δR)`, the dual configuration has interval part `(αI, chargeOf αI (βI + c₁))`,
`c₁ = cutLabel βL` (left gate), left part `(αL, chargeOf αL βL)` and right part
`(γR, chargeOf γR (δR + cutLabel βI))` (right gate). -/
lemma dualEquiv_symm_splitMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (αI βI : Cfg m)
    (αL βL : Cfg ℓ₁) (γR δR : Cfg (N - ℓ₁ - m)) :
    dualEquiv.symm ((splitMid ℓ₁ m h).symm ((αI, βI), ((αL, βL), (γR, δR)))) =
      (splitMid ℓ₁ m h).symm ((αI, chargeOf αI (βI + fun _ => cutLabel βL)),
        ((αL, chargeOf αL βL), (γR, chargeOf γR (δR + fun _ => cutLabel βI)))) := by
  rw [splitMid_symm_apply, splitMid_symm_apply, CutLocalAux.splitLGT_symm_apply (N := N - ℓ₁),
    dualEquiv_symm_splitLGT, CutLocalAux.splitLGT_symm_apply (N := N - ℓ₁)]
  have key : chargeOf ((splitCfg (N := N - ℓ₁) m (by omega)).symm (αI, γR))
      ((splitCfg (N := N - ℓ₁) m (by omega)).symm (βI, δR) + fun _ => cutLabel βL) =
      (splitCfg (N := N - ℓ₁) m (by omega)).symm
        (chargeOf αI (βI + fun _ => cutLabel βL), chargeOf γR (δR + fun _ => cutLabel βI)) := by
    rw [splitCfg_symm_add_const]
    refine (Equiv.eq_symm_apply _).2 ?_
    rw [splitCfg_chargeOf, cutLabel_add_const hm, Prod.mk.injEq]
    refine ⟨rfl, ?_⟩
    congr 1
    funext k
    simp only [Pi.add_apply]
    linear_combination zmod2_add_self (cutLabel βL)
  rw [key]

/-- The amplitudes of `UΦ` without the left gate: the interval and the left part in dual
coordinates, and the right part in dual coordinates after shifting its link labels by the label of
the last link of the interval (the right gate). -/
def amp (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) (a : Cfg m × Cfg m)
    (b : (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) : ℂ :=
  Φ ((splitMid ℓ₁ m h).symm (dualEquiv.symm a,
    (dualEquiv.symm b.1, dualEquiv.symm (b.2.1, b.2.2 + fun _ => cutLabel a.2))))

/-- Flipping all link labels: `(α, β) ↦ (α, β + 1)`. -/
def flipLinks (K : ℕ) : Cfg K × Cfg K ≃ Cfg K × Cfg K :=
  Equiv.prodCongr (Equiv.refl _) (Equiv.addRight fun _ => 1)

/-- Flipping all link labels. -/
@[simp] lemma flipLinks_apply (K : ℕ) (a : Cfg K × Cfg K) :
    flipLinks K a = (a.1, a.2 + fun _ => 1) := rfl

/-- Flipping all link labels of the right part. -/
def flipRight (ℓ₁ K : ℕ) :
    (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg K × Cfg K) ≃ (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg K × Cfg K) :=
  Equiv.prodCongr (Equiv.refl _) (flipLinks K)

/-- Flipping all link labels of the right part. -/
@[simp] lemma flipRight_apply (ℓ₁ K : ℕ) (b : (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg K × Cfg K)) :
    flipRight ℓ₁ K b = (b.1, flipLinks K b.2) := rfl

/-- `UΦ` on configurations whose left control label is `0`: no left gate. -/
lemma U_splitMid_of_zero (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (Φ : Cfg N × Cfg N → ℂ)
    (a : Cfg m × Cfg m) (l : Cfg ℓ₁ × Cfg ℓ₁) (r : Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))
    (hl : cutLabel l.2 = 0) :
    Φ (dualEquiv.symm ((splitMid ℓ₁ m h).symm (a, (l, r)))) = amp ℓ₁ m h Φ a (l, r) := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨αL, βL⟩ := l
  obtain ⟨γR, δR⟩ := r
  simp only at hl
  rw [dualEquiv_symm_splitMid ℓ₁ m h hm, hl, add_const_zero]
  rfl

/-- `UΦ` on configurations whose left control label is `1`: the left gate flips all link labels
of the interval, and (through the right gate) all link labels of the right part. -/
lemma U_splitMid_of_one (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (Φ : Cfg N × Cfg N → ℂ)
    (a : Cfg m × Cfg m) (l : Cfg ℓ₁ × Cfg ℓ₁) (r : Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))
    (hl : cutLabel l.2 = 1) :
    Φ (dualEquiv.symm ((splitMid ℓ₁ m h).symm (a, (l, r)))) =
      amp ℓ₁ m h Φ (flipLinks m a) (l, flipLinks _ r) := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨αL, βL⟩ := l
  obtain ⟨γR, δR⟩ := r
  simp only at hl
  rw [dualEquiv_symm_splitMid ℓ₁ m h hm, hl]
  simp only [amp, flipLinks_apply, dualEquiv_symm_apply, cutLabel_add_const hm]
  have e : (δR + fun _ => cutLabel βI) =
      (δR + fun _ => (1 : ZMod 2)) + fun _ => cutLabel βI + 1 := by
    funext k
    simp only [Pi.add_apply]
    linear_combination (-1 : ZMod 2) * zmod2_add_self (1 : ZMod 2)
  rw [e]

/-! ### The control parts `ω_c` -/

/-- The factor `X_c` of the part `ω_c = X_c X_c†` with left control label `c`. -/
def omegaFac (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) (c : ZMod 2) :
    Matrix (Cfg m × Cfg m) ((Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) ℂ :=
  Matrix.of fun a b => if cutLabel b.1.2 = c then amp ℓ₁ m h Φ a b else 0

/-- The part `ω_c` (left control label `c`) of the reduced state of `UΦ` without the left gate. -/
def omegaC (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) (c : ZMod 2) :
    Matrix (Cfg m × Cfg m) (Cfg m × Cfg m) ℂ :=
  omegaFac ℓ₁ m h Φ c * (omegaFac ℓ₁ m h Φ c)ᴴ

/-- `ω_c` is positive semidefinite. -/
lemma omegaC_posSemidef (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) (c : ZMod 2) :
    (omegaC ℓ₁ m h Φ c).PosSemidef :=
  Matrix.posSemidef_self_mul_conjTranspose _

/-- Entries of `ω_c`. -/
lemma omegaC_apply (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) (c : ZMod 2)
    (a a' : Cfg m × Cfg m) :
    omegaC ℓ₁ m h Φ c a a' = ∑ b, if cutLabel b.1.2 = c then
      amp ℓ₁ m h Φ a b * star (amp ℓ₁ m h Φ a' b) else 0 := by
  simp only [omegaC, omegaFac, Matrix.mul_apply, Matrix.conjTranspose_apply, Matrix.of_apply]
  refine Finset.sum_congr rfl fun b _ => ?_
  split_ifs <;> simp

/-- Entries of `F ω₁ F†`, `F` the flip of all link labels of the interval, after flipping all link
labels of the right part in the sum. -/
lemma omegaC_one_flip (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ)
    (a a' : Cfg m × Cfg m) :
    omegaC ℓ₁ m h Φ 1 (flipLinks m a) (flipLinks m a') =
      ∑ b, if cutLabel b.1.2 = 1 then amp ℓ₁ m h Φ (flipLinks m a) (flipRight ℓ₁ _ b) *
        star (amp ℓ₁ m h Φ (flipLinks m a') (flipRight ℓ₁ _ b)) else 0 := by
  rw [omegaC_apply]
  exact (Equiv.sum_comp (flipRight ℓ₁ (N - ℓ₁ - m)) (fun b => if cutLabel b.1.2 = 1 then
    amp ℓ₁ m h Φ (flipLinks m a) b * star (amp ℓ₁ m h Φ (flipLinks m a') b) else 0)).symm

/-- **The left gate**: `ρ_I(UΦ) = ω₀ + F ω₁ F†`, `F` the flip of all link labels of the interval. -/
lemma rdmMid_U_eq (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (Φ : Cfg N × Cfg N → ℂ) :
    rdmMid ℓ₁ m h (fun p => Φ (dualEquiv.symm p)) =
      omegaC ℓ₁ m h Φ 0 + (omegaC ℓ₁ m h Φ 1).submatrix (flipLinks m) (flipLinks m) := by
  ext a a'
  rw [Matrix.add_apply, Matrix.submatrix_apply, omegaC_apply, omegaC_one_flip,
    ← Finset.sum_add_distrib]
  simp only [rdmMid]
  refine Finset.sum_congr rfl fun ⟨l, r⟩ _ => ?_
  rcases zmod2_cases (cutLabel l.2) with hl | hl
  · rw [U_splitMid_of_zero ℓ₁ m h hm Φ a l r hl, U_splitMid_of_zero ℓ₁ m h hm Φ a' l r hl]
    simp [hl]
  · rw [U_splitMid_of_one ℓ₁ m h hm Φ a l r hl, U_splitMid_of_one ℓ₁ m h hm Φ a' l r hl]
    simp [hl]

/-- Entries of `ω₀ + ω₁`: the reduced state of `UΦ` without the left gate. -/
lemma omegaC_add_apply (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ)
    (a a' : Cfg m × Cfg m) :
    (omegaC ℓ₁ m h Φ 0 + omegaC ℓ₁ m h Φ 1) a a' =
      ∑ b, amp ℓ₁ m h Φ a b * star (amp ℓ₁ m h Φ a' b) := by
  rw [Matrix.add_apply, omegaC_apply, omegaC_apply, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun b _ => ?_
  rcases zmod2_cases (cutLabel b.1.2) with hb | hb <;> simp [hb]

/-- **The right gate**: the pinchings of `ω₀ + ω₁` and of the relabelled `ρ_{I'}(Φ)` by the
projection onto the label `0` of the last link of the interval coincide (substitute
`βL ↦ chargeOf αL βL` and `δR ↦ chargeOf γR (δR + c)` in the sum over the complement). -/
lemma pinching_omegaC (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Φ : Cfg N × Cfg N → ℂ) :
    pinching (cutProj m) (omegaC ℓ₁ m h Φ 0 + omegaC ℓ₁ m h Φ 1) =
      pinching (cutProj m) ((rdmMid ℓ₁ m h Φ).submatrix dualEquiv.symm dualEquiv.symm) := by
  ext a a'
  rw [pinching_cutProj_apply, pinching_cutProj_apply]
  split_ifs with hc
  · rw [omegaC_add_apply, submatrix_apply]
    simp only [amp, ← hc]
    exact Equiv.sum_comp (Equiv.prodCongr dualEquiv.symm
      ((Equiv.prodCongr (Equiv.refl _) (Equiv.addRight fun _ => cutLabel a.2)).trans
        dualEquiv.symm))
      (fun b => Φ ((splitMid ℓ₁ m h).symm (dualEquiv.symm a, b)) *
        star (Φ ((splitMid ℓ₁ m h).symm (dualEquiv.symm a', b))))
  · rfl

/-! ### Entropy estimates -/

/-- Relabelling the basis by a bijection does not change the trace. -/
lemma trace_submatrix_equiv {ι κ : Type*} [Fintype ι] [Fintype κ] (A : Matrix ι ι ℂ)
    (e : κ ≃ ι) : (A.submatrix e e).trace = A.trace := by
  simp only [Matrix.trace, Matrix.diag, Matrix.submatrix_apply]
  exact Equiv.sum_comp e (fun i => A i i)

/-- Two density matrices with the same two-block pinching have entropies differing by at most
`ln 2` (as in `cutlocal_defect`). -/
lemma abs_vonNeumann_sub_le_of_pinching {ι : Type*} [Fintype ι] [DecidableEq ι]
    (P σ σ' : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) (hσ' : σ'.PosSemidef) (htr : σ.trace = 1) (htr' : σ'.trace = 1)
    (heq : pinching P σ = pinching P σ') :
    |vonNeumann σ hσ.isHermitian - vonNeumann σ' hσ'.isHermitian| ≤ Real.log 2 := by
  have hpinH := (pinching_posSemidef P σ hP hσ).isHermitian
  have hpinH' := (pinching_posSemidef P σ' hP hσ').isHermitian
  have h1 := vonNeumann_le_pinching P σ hP hPP hσ hpinH
  have h2 := vonNeumann_pinching_le P σ' hP hPP hσ' htr' hpinH'
  have h3 := vonNeumann_le_pinching P σ' hP hPP hσ' hpinH'
  have h4 := vonNeumann_pinching_le P σ hP hPP hσ htr hpinH
  have h5 : vonNeumann (pinching P σ) hpinH = vonNeumann (pinching P σ') hpinH' :=
    CutLocalAux.vonNeumann_congr _ _ heq
  rw [abs_sub_le_iff]
  constructor <;> linarith

/-- Concavity and the mixing bound for a mixture of two density matrices:
`p₀ S(σ₀) + p₁ S(σ₁) ≤ S(p₀ σ₀ + p₁ σ₁) ≤ p₀ S(σ₀) + p₁ S(σ₁) + ln 2`. -/
lemma mix_two_bounds {ι : Type*} [Fintype ι] [DecidableEq ι] (σ₀ σ₁ : Matrix ι ι ℂ)
    (h₀ : σ₀.PosSemidef) (h₁ : σ₁.PosSemidef) (t₀ : σ₀.trace = 1) (t₁ : σ₁.trace = 1)
    (p₀ p₁ : ℝ) (hp₀ : 0 ≤ p₀) (hp₁ : 0 ≤ p₁) (hp : p₀ + p₁ = 1)
    (hM : ((p₀ : ℂ) • σ₀ + (p₁ : ℂ) • σ₁).IsHermitian) :
    p₀ * vonNeumann σ₀ h₀.isHermitian + p₁ * vonNeumann σ₁ h₁.isHermitian ≤
        vonNeumann ((p₀ : ℂ) • σ₀ + (p₁ : ℂ) • σ₁) hM ∧
      vonNeumann ((p₀ : ℂ) • σ₀ + (p₁ : ℂ) • σ₁) hM ≤
        p₀ * vonNeumann σ₀ h₀.isHermitian + p₁ * vonNeumann σ₁ h₁.isHermitian + Real.log 2 := by
  have hσ : ∀ i, (![σ₀, σ₁] i).PosSemidef := Fin.forall_fin_two.2 ⟨h₀, h₁⟩
  have ht : ∀ i, (![σ₀, σ₁] i).trace = 1 := Fin.forall_fin_two.2 ⟨t₀, t₁⟩
  have hw0 : ∀ i, 0 ≤ ![p₀, p₁] i := Fin.forall_fin_two.2 ⟨hp₀, hp₁⟩
  have hw1 : ∑ i, ![p₀, p₁] i = 1 := by rw [Fin.sum_univ_two]; exact hp
  have hsum : ∑ i, ((![p₀, p₁] i : ℝ) : ℂ) • ![σ₀, σ₁] i = (p₀ : ℂ) • σ₀ + (p₁ : ℂ) • σ₁ := by
    rw [Fin.sum_univ_two]; rfl
  have hm : (∑ i, ((![p₀, p₁] i : ℝ) : ℂ) • ![σ₀, σ₁] i).IsHermitian := hsum ▸ hM
  have hS : vonNeumann (∑ i, ((![p₀, p₁] i : ℝ) : ℂ) • ![σ₀, σ₁] i) hm =
      vonNeumann ((p₀ : ℂ) • σ₀ + (p₁ : ℂ) • σ₁) hM := CutLocalAux.vonNeumann_congr _ _ hsum
  have hA : ∑ i, ![p₀, p₁] i * vonNeumann (![σ₀, σ₁] i) (hσ i).isHermitian =
      p₀ * vonNeumann σ₀ h₀.isHermitian + p₁ * vonNeumann σ₁ h₁.isHermitian := by
    rw [Fin.sum_univ_two]; rfl
  have hH : ∑ i, Real.negMulLog (![p₀, p₁] i) ≤ Real.log 2 := by
    have := VNEntropyAux.sum_negMulLog_le_log_card (![p₀, p₁]) hw0 hw1
    simpa using this
  have hc := vonNeumann_concave (![p₀, p₁]) hw0 hw1 (![σ₀, σ₁]) hσ hm
  have hx := vonNeumann_mix_le (![p₀, p₁]) hw0 (![σ₀, σ₁]) hσ ht hm
  rw [hA, hS] at hc hx
  constructor <;> linarith

/-- **Mixing defect**: for positive semidefinite `ω₀, ω₁` with `tr(ω₀ + ω₁) = 1` and a relabelling
`e` of the basis, `|S(ω₀ + ω₁ ∘ e) - S(ω₀ + ω₁)| ≤ ln 2` (both lie in
`[p₀ S(σ₀) + p₁ S(σ₁), p₀ S(σ₀) + p₁ S(σ₁) + ln 2]`, `ω_c = p_c σ_c`). -/
lemma abs_vonNeumann_mix_sub_le {ι : Type*} [Fintype ι] [DecidableEq ι] (ω₀ ω₁ : Matrix ι ι ℂ)
    (h₀ : ω₀.PosSemidef) (h₁ : ω₁.PosSemidef) (htr : (ω₀ + ω₁).trace = 1) (e : ι ≃ ι)
    (hA : (ω₀ + ω₁.submatrix e e).IsHermitian) (hB : (ω₀ + ω₁).IsHermitian) :
    |vonNeumann (ω₀ + ω₁.submatrix e e) hA - vonNeumann (ω₀ + ω₁) hB| ≤ Real.log 2 := by
  have hlog : 0 ≤ Real.log 2 := Real.log_nonneg (by norm_num)
  have hre : ∀ ω : Matrix ι ι ℂ, ω.PosSemidef →
      0 ≤ ω.trace.re ∧ ω.trace = ((ω.trace.re : ℝ) : ℂ) := by
    intro ω hω
    have h := Complex.nonneg_iff.1 hω.trace_nonneg
    exact ⟨h.1, Complex.ext (by simp) (by simp [← h.2])⟩
  obtain ⟨hp₀, ht₀⟩ := hre ω₀ h₀
  obtain ⟨hp₁, ht₁⟩ := hre ω₁ h₁
  generalize ω₀.trace.re = p₀ at hp₀ ht₀
  generalize ω₁.trace.re = p₁ at hp₁ ht₁
  have hp : p₀ + p₁ = 1 := by
    have h := htr
    rw [Matrix.trace_add, ht₀, ht₁, ← Complex.ofReal_add] at h
    exact_mod_cast h
  rcases hp₀.eq_or_lt with h0 | h0
  · -- `ω₀ = 0`
    have hz : ω₀ = 0 := h₀.trace_eq_zero_iff.1 (by rw [ht₀, ← h0, Complex.ofReal_zero])
    subst hz
    have e1 : vonNeumann (0 + ω₁.submatrix e e) hA = vonNeumann ω₁ h₁.isHermitian := by
      rw [CutLocalAux.vonNeumann_congr hA (h₁.submatrix e).isHermitian (zero_add _)]
      exact vonNeumann_submatrix_equiv _ _ _ _
    have e2 : vonNeumann (0 + ω₁) hB = vonNeumann ω₁ h₁.isHermitian :=
      CutLocalAux.vonNeumann_congr _ _ (zero_add _)
    rw [e1, e2, sub_self, abs_zero]
    exact hlog
  rcases hp₁.eq_or_lt with h1 | h1
  · -- `ω₁ = 0`
    have hz : ω₁ = 0 := h₁.trace_eq_zero_iff.1 (by rw [ht₁, ← h1, Complex.ofReal_zero])
    subst hz
    have e1 : vonNeumann (ω₀ + (0 : Matrix ι ι ℂ).submatrix e e) hA = vonNeumann (ω₀ + 0) hB :=
      CutLocalAux.vonNeumann_congr _ _ rfl
    rw [e1, sub_self, abs_zero]
    exact hlog
  -- both parts are nonzero: normalize them
  have hσ₀ : (((p₀⁻¹ : ℝ) : ℂ) • ω₀).PosSemidef :=
    h₀.smul (by exact_mod_cast (inv_nonneg.2 hp₀))
  have hσ₁ : (((p₁⁻¹ : ℝ) : ℂ) • ω₁).PosSemidef :=
    h₁.smul (by exact_mod_cast (inv_nonneg.2 hp₁))
  have hσ₁' : ((((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e).PosSemidef := hσ₁.submatrix e
  have tσ₀ : (((p₀⁻¹ : ℝ) : ℂ) • ω₀).trace = 1 := by
    rw [Matrix.trace_smul, ht₀, smul_eq_mul, ← Complex.ofReal_mul, inv_mul_cancel₀ h0.ne',
      Complex.ofReal_one]
  have tσ₁ : (((p₁⁻¹ : ℝ) : ℂ) • ω₁).trace = 1 := by
    rw [Matrix.trace_smul, ht₁, smul_eq_mul, ← Complex.ofReal_mul, inv_mul_cancel₀ h1.ne',
      Complex.ofReal_one]
  have tσ₁' : ((((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e).trace = 1 := by
    rw [trace_submatrix_equiv, tσ₁]
  have hsc : ∀ (p : ℝ) (ω : Matrix ι ι ℂ), p ≠ 0 → (p : ℂ) • (((p⁻¹ : ℝ) : ℂ) • ω) = ω := by
    intro p ω hp
    rw [smul_smul, ← Complex.ofReal_mul, mul_inv_cancel₀ hp, Complex.ofReal_one, one_smul]
  have eq1 : (p₀ : ℂ) • (((p₀⁻¹ : ℝ) : ℂ) • ω₀) + (p₁ : ℂ) • (((p₁⁻¹ : ℝ) : ℂ) • ω₁) =
      ω₀ + ω₁ := by
    rw [hsc p₀ ω₀ h0.ne', hsc p₁ ω₁ h1.ne']
  have eq2 : (p₀ : ℂ) • (((p₀⁻¹ : ℝ) : ℂ) • ω₀) +
      (p₁ : ℂ) • ((((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e) = ω₀ + ω₁.submatrix e e := by
    rw [show (((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e = ((p₁⁻¹ : ℝ) : ℂ) • ω₁.submatrix e e from rfl,
      hsc p₀ ω₀ h0.ne', hsc p₁ _ h1.ne']
  have hM1 : ((p₀ : ℂ) • (((p₀⁻¹ : ℝ) : ℂ) • ω₀) +
      (p₁ : ℂ) • (((p₁⁻¹ : ℝ) : ℂ) • ω₁)).IsHermitian := eq1 ▸ hB
  have hM2 : ((p₀ : ℂ) • (((p₀⁻¹ : ℝ) : ℂ) • ω₀) +
      (p₁ : ℂ) • ((((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e)).IsHermitian := eq2 ▸ hA
  obtain ⟨a1, a2⟩ := mix_two_bounds _ _ hσ₀ hσ₁ tσ₀ tσ₁ p₀ p₁ hp₀ hp₁ hp hM1
  obtain ⟨b1, b2⟩ := mix_two_bounds _ _ hσ₀ hσ₁' tσ₀ tσ₁' p₀ p₁ hp₀ hp₁ hp hM2
  have hS1 : vonNeumann _ hM1 = vonNeumann (ω₀ + ω₁) hB := CutLocalAux.vonNeumann_congr _ _ eq1
  have hS2 : vonNeumann _ hM2 = vonNeumann (ω₀ + ω₁.submatrix e e) hA :=
    CutLocalAux.vonNeumann_congr _ _ eq2
  have hS3 : vonNeumann ((((p₁⁻¹ : ℝ) : ℂ) • ω₁).submatrix e e) hσ₁'.isHermitian =
      vonNeumann (((p₁⁻¹ : ℝ) : ℂ) • ω₁) hσ₁.isHermitian :=
    vonNeumann_submatrix_equiv _ _ _ _
  rw [hS1] at a1 a2
  rw [hS2, hS3] at b1 b2
  rw [abs_sub_le_iff]
  constructor <;> linarith

end InteriorAux

open InteriorAux EntropiesAux in
/-- **Interior intervals** (paper: Remark `rem:interior`): for every normalized state `Φ` of the dual
variables `(n, q)`, the entanglement entropy of the interval `[ℓ₁, ℓ₁ + m)` (sites and links) in
`UΦ`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, and that of the interval (modes and charges) in `Φ` differ
by at most `2 ln 2`. -/
theorem interior_defect (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (Φ : Cfg N × Cfg N → ℂ)
    (hΦ : ∑ p, Complex.normSq (Φ p) = 1) :
    |vonNeumann (rdmMid ℓ₁ m h (fun p => Φ (dualEquiv.symm p))) (rdmMid_isHermitian ℓ₁ m h _) -
      vonNeumann (rdmMid ℓ₁ m h Φ) (rdmMid_isHermitian ℓ₁ m h _)| ≤ 2 * Real.log 2 := by
  have h₀ := omegaC_posSemidef ℓ₁ m h Φ 0
  have h₁ := omegaC_posSemidef ℓ₁ m h Φ 1
  -- `ρ'' = ρ_{I'}(Φ)` relabelled by `(α, β) ↦ (α, chargeOf α β)`
  have hρ'' : ((rdmMid ℓ₁ m h Φ).submatrix (dualEquiv (N := m)).symm
      dualEquiv.symm).PosSemidef := (rdmMid_posSemidef ℓ₁ m h Φ).submatrix _
  have htr'' :
      ((rdmMid ℓ₁ m h Φ).submatrix (dualEquiv (N := m)).symm dualEquiv.symm).trace = 1 := by
    rw [InteriorAux.trace_submatrix_equiv, trace_rdmMid, hΦ, Complex.ofReal_one]
  have hpin := pinching_omegaC ℓ₁ m h Φ
  -- `ρ̃ = ω₀ + ω₁` has trace one: its diagonal is that of `ρ''`
  have htr : (omegaC ℓ₁ m h Φ 0 + omegaC ℓ₁ m h Φ 1).trace = 1 := by
    rw [← htr'']
    simp only [Matrix.trace, Matrix.diag]
    refine Finset.sum_congr rfl fun a _ => ?_
    have := congrFun (congrFun hpin a) a
    rwa [pinching_cutProj_apply, pinching_cutProj_apply, ite_eq_left rfl, ite_eq_left rfl] at this
  have hρt := h₀.add h₁
  have hρ := h₀.add (h₁.submatrix (flipLinks m))
  -- the left gate costs at most `ln 2`, the right gate at most `ln 2`
  have hc := abs_vonNeumann_mix_sub_le _ _ h₀ h₁ htr (flipLinks m) hρ.isHermitian
    hρt.isHermitian
  have hd := abs_vonNeumann_sub_le_of_pinching (cutProj m) _ _ (cutProj_isHermitian m)
    (cutProj_mul_self m) hρt hρ'' htr htr'' hpin
  have e1 : vonNeumann (rdmMid ℓ₁ m h (fun p => Φ (dualEquiv.symm p)))
      (rdmMid_isHermitian ℓ₁ m h _) = vonNeumann (omegaC ℓ₁ m h Φ 0 +
        (omegaC ℓ₁ m h Φ 1).submatrix (flipLinks m) (flipLinks m)) hρ.isHermitian :=
    CutLocalAux.vonNeumann_congr _ _ (rdmMid_U_eq ℓ₁ m h hm Φ)
  have e2 : vonNeumann ((rdmMid ℓ₁ m h Φ).submatrix (dualEquiv (N := m)).symm dualEquiv.symm)
      hρ''.isHermitian = vonNeumann (rdmMid ℓ₁ m h Φ) (rdmMid_isHermitian ℓ₁ m h _) :=
    vonNeumann_submatrix_equiv _ _ _ _
  rw [e1, ← e2]
  refine (abs_sub_le _
    (vonNeumann (omegaC ℓ₁ m h Φ 0 + omegaC ℓ₁ m h Φ 1) hρt.isHermitian) _).trans ?_
  linarith

/-- **Parity blocks** (paper: Remark `rem:interior`): if `Ψ` has definite total fermion parity, its
reduced state on an interval vanishes between configurations of the interval of different fermion
parity. -/
theorem rdmMid_parity_block (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) (π : ZMod 2)
    (hΨ : ∀ p, Ψ p ≠ 0 → ∑ k, p.1 k = π) (a a' : Cfg m × Cfg m)
    (ha : ∑ k, a.1 k ≠ ∑ k, a'.1 k) : rdmMid ℓ₁ m h Ψ a a' = 0 := by
  unfold rdmMid
  refine Finset.sum_eq_zero fun b _ => ?_
  by_cases h1 : Ψ ((splitMid ℓ₁ m h).symm (a, b)) = 0
  · rw [h1, zero_mul]
  by_cases h2 : Ψ ((splitMid ℓ₁ m h).symm (a', b)) = 0
  · rw [h2, star_zero, mul_zero]
  exfalso
  have e1 := hΨ _ h1
  have e2 := hΨ _ h2
  rw [InteriorAux.sum_splitMid_symm] at e1 e2
  exact ha (by linear_combination e1 - e2)

/-! ### Sharpness -/

namespace InteriorAux

open EntropiesAux

/-- The reduced density matrix of a product state `Ψ(a, b) = [C b] u a`: `|C| u u†`. -/
lemma rdmMid_prod (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ)
    (C : (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m)) → Prop) [DecidablePred C]
    (u : Cfg m × Cfg m → ℂ)
    (hΨ : ∀ a b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) = if C b then u a else 0)
    (a a' : Cfg m × Cfg m) :
    rdmMid ℓ₁ m h Ψ a a' = ((univ.filter C).card : ℂ) * (u a * star (u a')) := by
  have e : ∀ b, (if C b then u a else 0) * star (if C b then u a' else 0) =
      if C b then u a * star (u a') else 0 := by
    intro b
    split_ifs <;> simp
  simp only [rdmMid, hΨ, e]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, nsmul_eq_mul]

/-- The reduced density matrix of a state supported on the graph `b = g a` of a map from the
configurations of the interval to those of the complement. -/
lemma rdmMid_graph (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ)
    (g : Cfg m × Cfg m → (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m)))
    (f : Cfg m × Cfg m → ℂ)
    (hΨ : ∀ a b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) = if b = g a then f a else 0)
    (a a' : Cfg m × Cfg m) :
    rdmMid ℓ₁ m h Ψ a a' = if g a = g a' then f a * star (f a') else 0 := by
  simp only [rdmMid, hΨ, ite_mul, zero_mul]
  rw [Finset.sum_ite_eq' Finset.univ (g a)]
  simp only [Finset.mem_univ, ite_true]
  split_ifs <;> simp

/-- A rank-one projector `v v†` with `‖v‖ = 1` has entropy `0`. -/
lemma vonNeumann_rankOne {ι : Type*} [Fintype ι] [DecidableEq ι] (v : ι → ℂ)
    (hv : ∑ i, Complex.normSq (v i) = 1) (M : Matrix ι ι ℂ) (hM : M.IsHermitian)
    (hMv : ∀ i j, M i j = v i * star (v j)) : vonNeumann M hM = 0 := by
  obtain ⟨A, hA⟩ : ∃ A : Matrix ι Unit ℂ, A = Matrix.of fun i _ => v i := ⟨_, rfl⟩
  have hMA : M = A * Aᴴ := by
    ext i j
    simp [hA, Matrix.mul_apply, hMv]
  have hAA : Aᴴ * A = diagonal (fun _ : Unit => ((1 : ℝ) : ℂ)) := by
    ext u u'
    rw [Subsingleton.elim u' u, diagonal_apply_eq, Matrix.mul_apply]
    simp only [conjTranspose_apply, hA, Matrix.of_apply]
    rw [← hv, Complex.ofReal_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [mul_comm, mul_star_self_eq]
  have hD : (diagonal (fun _ : Unit => ((1 : ℝ) : ℂ))).IsHermitian :=
    hAA ▸ isHermitian_conjTranspose_mul_self A
  rw [CutLocalAux.vonNeumann_congr hM (isHermitian_mul_conjTranspose_self A) hMA,
    vonNeumann_mul_conjTranspose_comm, CutLocalAux.vonNeumann_congr _ hD hAA, vonNeumann_diagonal]
  simp

/-- `|2^{-1/2}|² = 1/2`. -/
lemma invSqrtTwo_mul_star :
    ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * star (((Real.sqrt 2 : ℝ) : ℂ)⁻¹) = ((1 / 2 : ℝ) : ℂ) := by
  rw [← Complex.ofReal_inv, Complex.star_def, Complex.conj_ofReal, ← Complex.ofReal_mul]
  congr 1
  rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  norm_num

/-- The vector `2^{-1/2} [P a]` with two nonzero entries is normalized. -/
lemma sum_normSq_invSqrtTwo {ι : Type*} [Fintype ι] (P : ι → Prop) [DecidablePred P]
    (hP : (univ.filter P).card = 2) :
    ∑ a, Complex.normSq (if P a then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0) = 1 := by
  have e : ∀ a, Complex.normSq (if P a then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0) =
      if P a then 1 / 2 else 0 := by
    intro a
    split_ifs
    · rw [map_inv₀, Complex.normSq_ofReal, Real.mul_self_sqrt (by norm_num)]
      norm_num
    · simp
  simp only [e]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, hP, nsmul_eq_mul]
  norm_num

/-- A state with amplitudes `½` on four configurations (and `0` elsewhere) is normalized. -/
lemma sum_normSq_half {ι : Type*} [Fintype ι] (P : ι → Prop) [DecidablePred P]
    (hP : (univ.filter P).card = 4) :
    ∑ a, Complex.normSq (if P a then (1 / 2 : ℂ) else 0) = 1 := by
  have e : ∀ a, Complex.normSq (if P a then (1 / 2 : ℂ) else 0) = if P a then 1 / 4 else 0 := by
    intro a
    split_ifs
    · rw [Complex.normSq_div, Complex.normSq_one, Complex.normSq_ofNat]
      norm_num
    · simp
  simp only [e]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, hP, nsmul_eq_mul]
  norm_num

/-- The Shannon entropy of four entries `1/4` (and zeros) is `2 ln 2`. -/
lemma sum_negMulLog_quarter {ι : Type*} [Fintype ι] (P : ι → Prop) [DecidablePred P]
    (hP : (univ.filter P).card = 4) :
    ∑ a, Real.negMulLog (if P a then 1 / 4 else 0) = 2 * Real.log 2 := by
  have e : ∀ a, Real.negMulLog (if P a then 1 / 4 else 0) =
      if P a then Real.negMulLog (1 / 4) else 0 := by
    intro a
    split_ifs <;> simp
  simp only [e]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, hP, nsmul_eq_mul]
  have h4 : Real.log 4 = 2 * Real.log 2 := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.log_pow]
    norm_num
  simp only [Real.negMulLog, one_div, Real.log_inv, h4]
  push_cast
  ring

/-- **Product states**: if `Ψ(a, b) = [C b] [P a] / 2` with `|C| = |P| = 2`, the reduced state of
the interval is a rank-one projector and its entropy vanishes. -/
lemma vonNeumann_rdmMid_prod_half (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ)
    (C : (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m)) → Prop) [DecidablePred C]
    (P : Cfg m × Cfg m → Prop) [DecidablePred P]
    (hC : (univ.filter C).card = 2) (hP : (univ.filter P).card = 2)
    (hΨ : ∀ a b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) =
      if C b then (if P a then 1 / 2 else 0) else 0) :
    vonNeumann (rdmMid ℓ₁ m h Ψ) (rdmMid_isHermitian ℓ₁ m h Ψ) = 0 := by
  refine vonNeumann_rankOne (fun a => if P a then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0)
    (sum_normSq_invSqrtTwo P hP) _ _ ?_
  intro a a'
  rw [rdmMid_prod ℓ₁ m h Ψ C (fun a => if P a then 1 / 2 else 0) hΨ a a', hC]
  by_cases ha : P a <;> by_cases ha' : P a' <;>
    simp only [ha, ha', ite_true, ite_false, star_zero, mul_zero, zero_mul]
  rw [invSqrtTwo_mul_star]
  norm_num

/-- **Graph states**: if `Ψ(a, b) = [b = g a] [P a] / 2` with `g` injective on `P` and `|P| = 4`,
the reduced state of the interval is `diag([P a] / 4)`, with entropy `2 ln 2`. -/
lemma vonNeumann_rdmMid_graph_quarter (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ)
    (g : Cfg m × Cfg m → (Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m)))
    (P : Cfg m × Cfg m → Prop) [DecidablePred P] (hP : (univ.filter P).card = 4)
    (hg : ∀ a a', g a = g a' → P a → P a' → a = a')
    (hΨ : ∀ a b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) =
      if b = g a then (if P a then 1 / 2 else 0) else 0) :
    vonNeumann (rdmMid ℓ₁ m h Ψ) (rdmMid_isHermitian ℓ₁ m h Ψ) = 2 * Real.log 2 := by
  have hD : rdmMid ℓ₁ m h Ψ = diagonal (fun a => (((if P a then 1 / 4 else 0 : ℝ)) : ℂ)) := by
    ext a a'
    rw [rdmMid_graph ℓ₁ m h Ψ g (fun a => if P a then 1 / 2 else 0) hΨ, diagonal_apply]
    by_cases ha : a = a'
    · subst ha
      by_cases h0 : P a
      · simp only [h0, ite_true]
        norm_num
      · simp only [h0, ite_true, ite_false, star_zero, mul_zero, Complex.ofReal_zero]
    · by_cases hga : g a = g a'
      · by_cases h0 : P a
        · by_cases h0' : P a'
          · exact absurd (hg a a' hga h0 h0') ha
          · simp only [ha, hga, h0', ite_true, ite_false, star_zero, mul_zero]
        · simp only [ha, hga, h0, ite_true, ite_false, zero_mul]
      · simp only [ha, hga, ite_false]
  have hH : (diagonal (fun a => (((if P a then 1 / 4 else 0 : ℝ)) : ℂ))).IsHermitian :=
    hD ▸ rdmMid_isHermitian ℓ₁ m h Ψ
  rw [CutLocalAux.vonNeumann_congr _ hH hD, vonNeumann_diagonal]
  exact sum_negMulLog_quarter P hP

/-- `|S - 0| = 2 ln 2` if `S = 2 ln 2`. -/
lemma abs_two_log_two_sub_zero : |2 * Real.log 2 - 0| = 2 * Real.log 2 := by
  rw [sub_zero, abs_of_nonneg (mul_nonneg zero_le_two (Real.log_nonneg one_le_two))]

end InteriorAux

/-- `Φ = ½ |n₀ = 0, n₂ = 0⟩ ⊗ (|n₁ = 0⟩ + |n₁ = 1⟩) ⊗ (|q₀ = 0⟩ + |q₀ = 1⟩) ⊗ |q₁ = 0, q₂ = 0⟩`
(`N = 3`), a product state across the interval `I'` = (mode `1`, charge `1`). -/
def exInt : Cfg 3 × Cfg 3 → ℂ :=
  fun p => if p.1 0 = 0 ∧ p.1 2 = 0 ∧ p.2 1 = 0 ∧ p.2 2 = 0 then 1 / 2 else 0

namespace InteriorAux

set_option synthInstance.maxSize 2048

/-- Support of `exInt` in split coordinates: interval part `(n₁, q₁)`, complement
`((n₀, q₀), (n₂, q₂))`. -/
lemma exInt_support : ∀ (αI βI αL βL : Cfg 1) (γR δR : Cfg (3 - 1 - 1)),
    (((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).1 0 = 0 ∧
      ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).1 2 = 0 ∧
      ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).2 1 = 0 ∧
      ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).2 2 = 0) ↔
      ((αL = 0 ∧ γR = 0 ∧ δR = 0) ∧ βI = 0) := by
  decide

/-- `exInt` in split coordinates: `[n₀ = 0, n₂ = 0, q₂ = 0] [q₁ = 0] / 2`. -/
lemma exInt_split (a : Cfg 1 × Cfg 1)
    (b : (Cfg 1 × Cfg 1) × (Cfg (3 - 1 - 1) × Cfg (3 - 1 - 1))) :
    exInt ((splitMid 1 1 (by norm_num)).symm (a, b)) =
      if (b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0) then (if a.2 = 0 then 1 / 2 else 0) else 0 := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨⟨αL, βL⟩, ⟨γR, δR⟩⟩ := b
  have key := exInt_support αI βI αL βL γR δR
  simp only [exInt]
  by_cases hC : (αL = 0 ∧ γR = 0 ∧ δR = 0)
  · by_cases hu : βI = 0
    · rw [ite_eq_left (key.2 ⟨hC, hu⟩), ite_eq_left hC, ite_eq_left hu]
    · rw [ite_eq_right (fun h => hu (key.1 h).2), ite_eq_left hC, ite_eq_right hu]
  · rw [ite_eq_right (fun h => hC (key.1 h).1), ite_eq_right hC]

/-- The complement configuration `((0, x₁ + n₁), (0, x₁))` carried by `U exInt`. -/
def gU (a : Cfg 1 × Cfg 1) : (Cfg 1 × Cfg 1) × (Cfg (3 - 1 - 1) × Cfg (3 - 1 - 1)) :=
  ((0, fun _ => a.2 0 + a.1 0), (0, fun _ => a.2 0))

/-- Support of `U exInt` in split coordinates. -/
lemma exIntU_support : ∀ (αI βI αL βL : Cfg 1) (γR δR : Cfg (3 - 1 - 1)),
    ((dualEquiv.symm ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).1 0
        = 0 ∧
      (dualEquiv.symm ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).1 2
        = 0 ∧
      (dualEquiv.symm ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).2 1
        = 0 ∧
      (dualEquiv.symm ((splitMid 1 1 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).2 2
        = 0) ↔
      (αL = 0 ∧ βL = (fun _ => βI 0 + αI 0) ∧ γR = 0 ∧ δR = (fun _ => βI 0)) := by
  decide

/-- `U exInt` in split coordinates: supported on the graph of `gU`, with amplitude `½`. -/
lemma exIntU_split (a : Cfg 1 × Cfg 1)
    (b : (Cfg 1 × Cfg 1) × (Cfg (3 - 1 - 1) × Cfg (3 - 1 - 1))) :
    exInt (dualEquiv.symm ((splitMid 1 1 (by norm_num)).symm (a, b))) =
      if b = gU a then 1 / 2 else 0 := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨⟨αL, βL⟩, ⟨γR, δR⟩⟩ := b
  have key := exIntU_support αI βI αL βL γR δR
  have hb : ((αL, βL), (γR, δR)) = gU (αI, βI) ↔
      (αL = 0 ∧ βL = (fun _ => βI 0 + αI 0) ∧ γR = 0 ∧ δR = (fun _ => βI 0)) := by
    simp only [gU, Prod.mk.injEq, and_assoc]
  simp only [exInt]
  by_cases hg : ((αL, βL), (γR, δR)) = gU (αI, βI)
  · rw [ite_eq_left (key.2 (hb.1 hg)), ite_eq_left hg]
  · rw [ite_eq_right (fun h => hg (hb.2 (key.1 h))), ite_eq_right hg]

/-- `gU` is injective. -/
lemma gU_injective : ∀ a a' : Cfg 1 × Cfg 1, gU a = gU a' → a = a' := by
  decide

/-- Two complement configurations satisfy `n₀ = 0, n₂ = 0, q₂ = 0`. -/
lemma card_exInt_compl : (univ.filter (fun b : (Cfg 1 × Cfg 1) × (Cfg (3 - 1 - 1) ×
    Cfg (3 - 1 - 1)) => b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0)).card = 2 := by
  decide

/-- Two interval configurations satisfy `q₁ = 0`. -/
lemma card_exInt_int : (univ.filter (fun a : Cfg 1 × Cfg 1 => a.2 = 0)).card = 2 := by
  decide

/-- There are four interval configurations. -/
lemma card_exInt_all : (univ.filter (fun _ : Cfg 1 × Cfg 1 => True)).card = 4 := by
  decide

/-- Four configurations carry the amplitude of `exInt`. -/
lemma card_exInt : (univ.filter (fun p : Cfg 3 × Cfg 3 =>
    p.1 0 = 0 ∧ p.1 2 = 0 ∧ p.2 1 = 0 ∧ p.2 2 = 0)).card = 4 := by
  decide

end InteriorAux

open InteriorAux in
/-- `exInt` is normalized. -/
theorem exInt_norm : ∑ p, Complex.normSq (exInt p) = 1 :=
  sum_normSq_half _ card_exInt

open InteriorAux in
/-- `S_{I'}(Φ) = 0` for `Φ = exInt`, `I'` = (mode `1`, charge `1`). -/
theorem exInt_entropy :
    vonNeumann (rdmMid 1 1 (by norm_num) exInt) (rdmMid_isHermitian 1 1 (by norm_num) _) = 0 :=
  vonNeumann_rdmMid_prod_half 1 1 _ exInt (fun b => b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0)
    (fun a => a.2 = 0) card_exInt_compl card_exInt_int exInt_split

open InteriorAux in
/-- `S_I(UΦ) = 2 ln 2` for `Φ = exInt`, `I` = (site `1`, link `1`): the reduced state is maximally
mixed. -/
theorem exInt_U_entropy :
    vonNeumann (rdmMid 1 1 (by norm_num) (fun p => exInt (dualEquiv.symm p)))
      (rdmMid_isHermitian 1 1 (by norm_num) _) = 2 * Real.log 2 := by
  have hD : rdmMid 1 1 (by norm_num) (fun p => exInt (dualEquiv.symm p)) =
      diagonal (fun a : Cfg 1 × Cfg 1 => (((if True then 1 / 4 else 0 : ℝ)) : ℂ)) := by
    ext a a'
    rw [rdmMid_graph 1 1 _ _ gU (fun _ => 1 / 2) exIntU_split, diagonal_apply]
    by_cases ha : a = a'
    · rw [ite_eq_left (congrArg gU ha), ite_eq_left ha]
      norm_num
    · rw [ite_eq_right (fun hg => ha (gU_injective a a' hg)), ite_eq_right ha]
  have hH : (diagonal (fun a : Cfg 1 × Cfg 1 =>
      (((if True then 1 / 4 else 0 : ℝ)) : ℂ))).IsHermitian :=
    hD ▸ rdmMid_isHermitian 1 1 (by norm_num) _
  rw [CutLocalAux.vonNeumann_congr _ hH hD, vonNeumann_diagonal]
  exact sum_negMulLog_quarter _ card_exInt_all

/-- **The bound `2 ln 2` of `interior_defect` is attained.** -/
theorem interior_defect_attained :
    |vonNeumann (rdmMid 1 1 (by norm_num) (fun p => exInt (dualEquiv.symm p)))
        (rdmMid_isHermitian 1 1 (by norm_num) _) -
      vonNeumann (rdmMid 1 1 (by norm_num) exInt) (rdmMid_isHermitian 1 1 (by norm_num) _)| =
      2 * Real.log 2 := by
  rw [exInt_U_entropy, exInt_entropy]
  exact InteriorAux.abs_two_log_two_sub_zero

/-- `Φ = ½ |0000⟩ ⊗ (|q₀ = 0⟩ + |q₀ = 1⟩) ⊗ |q₁ = 0⟩ ⊗ (|q₂ = 0⟩ + |q₂ = 1⟩) ⊗ |q₃ = 0⟩` (`N = 4`):
the matter vacuum, a product state across the interval `I'` = (modes `1, 2`, charges `1, 2`). Its
image `UΦ = ½ |0000⟩ ⊗ ∑_{x₀, x₂} |x₀, x₀, x₂, x₂⟩` contains a Bell pair across each end of `I`. -/
def exVac : Cfg 4 × Cfg 4 → ℂ :=
  fun p => if p.1 = 0 ∧ p.2 1 = 0 ∧ p.2 3 = 0 then 1 / 2 else 0

namespace InteriorAux

set_option synthInstance.maxSize 2048

/-- Support of `exVac` in split coordinates: interval part `((n₁, n₂), (q₁, q₂))`, complement
`((n₀, q₀), (n₃, q₃))`. -/
lemma exVac_support : ∀ (αI βI : Cfg 2) (αL βL : Cfg 1) (γR δR : Cfg (4 - 1 - 2)),
    (((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).1 = 0 ∧
      ((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).2 1 = 0 ∧
      ((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR)))).2 3 = 0) ↔
      ((αL = 0 ∧ γR = 0 ∧ δR = 0) ∧ (αI = 0 ∧ βI 0 = 0)) := by
  decide

/-- `exVac` in split coordinates: `[n₀ = 0, n₃ = 0, q₃ = 0] [n₁ = n₂ = 0, q₁ = 0] / 2`. -/
lemma exVac_split (a : Cfg 2 × Cfg 2)
    (b : (Cfg 1 × Cfg 1) × (Cfg (4 - 1 - 2) × Cfg (4 - 1 - 2))) :
    exVac ((splitMid 1 2 (by norm_num)).symm (a, b)) =
      if (b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0) then
        (if a.1 = 0 ∧ a.2 0 = 0 then 1 / 2 else 0) else 0 := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨⟨αL, βL⟩, ⟨γR, δR⟩⟩ := b
  have key := exVac_support αI βI αL βL γR δR
  simp only [exVac]
  by_cases hC : (αL = 0 ∧ γR = 0 ∧ δR = 0)
  · by_cases hu : αI = 0 ∧ βI 0 = 0
    · rw [ite_eq_left (key.2 ⟨hC, hu⟩), ite_eq_left hC, ite_eq_left hu]
    · rw [ite_eq_right (fun h => hu (key.1 h).2), ite_eq_left hC, ite_eq_right hu]
  · rw [ite_eq_right (fun h => hC (key.1 h).1), ite_eq_right hC]

/-- The complement configuration `((0, x₁), (0, x₂))` carried by `U exVac`. -/
def gV (a : Cfg 2 × Cfg 2) : (Cfg 1 × Cfg 1) × (Cfg (4 - 1 - 2) × Cfg (4 - 1 - 2)) :=
  ((0, fun _ => a.2 0), (0, fun _ => a.2 1))

/-- Support of `U exVac` in split coordinates. -/
lemma exVacU_support : ∀ (αI βI : Cfg 2) (αL βL : Cfg 1) (γR δR : Cfg (4 - 1 - 2)),
    ((dualEquiv.symm ((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).1
        = 0 ∧
      (dualEquiv.symm ((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).2 1
        = 0 ∧
      (dualEquiv.symm ((splitMid 1 2 (by norm_num)).symm ((αI, βI), ((αL, βL), (γR, δR))))).2 3
        = 0) ↔
      ((αL = 0 ∧ βL = (fun _ => βI 0) ∧ γR = 0 ∧ δR = (fun _ => βI 1)) ∧ αI = 0) := by
  decide

/-- `U exVac` in split coordinates: supported on the graph of `gV`, with amplitude
`½ [n₁ = n₂ = 0]`. -/
lemma exVacU_split (a : Cfg 2 × Cfg 2)
    (b : (Cfg 1 × Cfg 1) × (Cfg (4 - 1 - 2) × Cfg (4 - 1 - 2))) :
    exVac (dualEquiv.symm ((splitMid 1 2 (by norm_num)).symm (a, b))) =
      if b = gV a then (if a.1 = 0 then 1 / 2 else 0) else 0 := by
  obtain ⟨αI, βI⟩ := a
  obtain ⟨⟨αL, βL⟩, ⟨γR, δR⟩⟩ := b
  have key := exVacU_support αI βI αL βL γR δR
  have hb : ((αL, βL), (γR, δR)) = gV (αI, βI) ↔
      (αL = 0 ∧ βL = (fun _ => βI 0) ∧ γR = 0 ∧ δR = (fun _ => βI 1)) := by
    simp only [gV, Prod.mk.injEq, and_assoc]
  simp only [exVac]
  by_cases hg : ((αL, βL), (γR, δR)) = gV (αI, βI)
  · by_cases hα : αI = 0
    · rw [ite_eq_left (key.2 ⟨hb.1 hg, hα⟩), ite_eq_left hg, ite_eq_left hα]
    · rw [ite_eq_right (fun h => hα (key.1 h).2), ite_eq_left hg, ite_eq_right hα]
  · rw [ite_eq_right (fun h => hg (hb.2 (key.1 h).1)), ite_eq_right hg]

/-- On the configurations with empty interval, `gV` is injective. -/
lemma gV_injective : ∀ a a' : Cfg 2 × Cfg 2, gV a = gV a' → a.1 = 0 → a'.1 = 0 → a = a' := by
  decide

/-- Two complement configurations satisfy `n₀ = 0, n₃ = 0, q₃ = 0`. -/
lemma card_exVac_compl : (univ.filter (fun b : (Cfg 1 × Cfg 1) × (Cfg (4 - 1 - 2) ×
    Cfg (4 - 1 - 2)) => b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0)).card = 2 := by
  decide

/-- Two interval configurations satisfy `n₁ = n₂ = 0, q₁ = 0`. -/
lemma card_exVac_int : (univ.filter (fun a : Cfg 2 × Cfg 2 => a.1 = 0 ∧ a.2 0 = 0)).card = 2 := by
  decide

/-- Four interval configurations have empty sites. -/
lemma card_exVac_empty : (univ.filter (fun a : Cfg 2 × Cfg 2 => a.1 = 0)).card = 4 := by
  decide

/-- Four configurations carry the amplitude of `exVac`. -/
lemma card_exVac : (univ.filter (fun p : Cfg 4 × Cfg 4 =>
    p.1 = 0 ∧ p.2 1 = 0 ∧ p.2 3 = 0)).card = 4 := by
  decide

end InteriorAux

open InteriorAux in
/-- `exVac` is normalized. -/
theorem exVac_norm : ∑ p, Complex.normSq (exVac p) = 1 :=
  sum_normSq_half _ card_exVac

open InteriorAux in
/-- `S_{I'}(Φ) = 0` for `Φ = exVac`, `I'` = (modes `1, 2`, charges `1, 2`). -/
theorem exVac_entropy :
    vonNeumann (rdmMid 1 2 (by norm_num) exVac) (rdmMid_isHermitian 1 2 (by norm_num) _) = 0 :=
  vonNeumann_rdmMid_prod_half 1 2 _ exVac (fun b => b.1.1 = 0 ∧ b.2.1 = 0 ∧ b.2.2 = 0)
    (fun a => a.1 = 0 ∧ a.2 0 = 0) card_exVac_compl card_exVac_int exVac_split

open InteriorAux in
/-- `S_I(UΦ) = 2 ln 2` for `Φ = exVac`, `I` = (sites `1, 2`, links `1, 2`). -/
theorem exVac_U_entropy :
    vonNeumann (rdmMid 1 2 (by norm_num) (fun p => exVac (dualEquiv.symm p)))
      (rdmMid_isHermitian 1 2 (by norm_num) _) = 2 * Real.log 2 :=
  vonNeumann_rdmMid_graph_quarter 1 2 _ _ gV (fun a => a.1 = 0) card_exVac_empty gV_injective
    exVacU_split

/-- **The bound `2 ln 2` of `interior_defect` is attained by a state with empty matter.** -/
theorem interior_defect_attained_vacuum :
    |vonNeumann (rdmMid 1 2 (by norm_num) (fun p => exVac (dualEquiv.symm p)))
        (rdmMid_isHermitian 1 2 (by norm_num) _) -
      vonNeumann (rdmMid 1 2 (by norm_num) exVac) (rdmMid_isHermitian 1 2 (by norm_num) _)| =
      2 * Real.log 2 := by
  rw [exVac_U_entropy, exVac_entropy]
  exact InteriorAux.abs_two_log_two_sub_zero

end DFL
