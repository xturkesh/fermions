import DFL.AreaLaw
import DFL.ProductLaw

/-!
# Full entanglement of definite-charge states (paper: Lemma `lem:cutlocal`(c),
Theorem `thm:arealaw`(c))

Let the links start in a `σˣ`-product state `|ξ⟩_σ` and the matter in a Fock state `|n⟩`; then the
LGT state stays in the single charge sector `q = chargeOf n ξ` (`DFL.expect_xprod_fock`). For the
region `A` = (sites `< ℓ`, links `< ℓ`), the full matter + gauge entanglement entropy of the LGT
state equals the entanglement entropy of the free-fermion state `e^{-it dΓ(h_q)} |n⟩` for the modes
`< ℓ`. Reason: in the LGT state the link labels in `A` are a function of the occupations in `A`
(Gauss law), and the link labels outside `A` only depend on the particle-number parity of `A`, which
is the same for all configurations in the support (particle-number conservation). Combined with
the area law of `DFL.AreaLaw`, the average over uniformly random `ξ` of the full entanglement is
bounded uniformly in time and volume.

Proof outline (auxiliary lemmas in `DFL.CutLocalAux`):
* in the dual frame the evolved state is `Ψ(t)(n', x') = [x' = x(n', q)] e^{iθ} φ(n')` with
  `φ = e^{-itdΓ(h_q)} |n⟩` (`evo_xprodInit_apply`);
* `φ` conserves the particle-number parity (`evo_dGamma_eq_zero`);
* hence `ρ_LGT = V ρ_A V†` with the isometry `V : α ↦ (α, L_A(α))` (`rdmLGT_evo_xprod`);
* `X^{|A|} χ(V ρ V†) = X^{|A × A|} χ(ρ)` (Weinstein–Aronszajn), so the nonzero eigenvalues agree
  and the von Neumann entropies coincide (`vonNeumann_conj_isometry`).
-/

open Matrix Finset NormedSpace

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Splitting an LGT configuration (occupations, link labels) into the region
`A` = (sites `< ℓ`, links `< ℓ`) and its complement. -/
def splitLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    Cfg N × Cfg N ≃ (Cfg ℓ × Cfg ℓ) × (Cfg (N - ℓ) × Cfg (N - ℓ)) :=
  (Equiv.prodCongr (splitCfg ℓ hℓ) (splitCfg ℓ hℓ)).trans (Equiv.prodProdProdComm _ _ _ _)

/-- Reduced density matrix of an LGT pure state `Ψ` on the region `A`. -/
def rdmLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    Matrix (Cfg ℓ × Cfg ℓ) (Cfg ℓ × Cfg ℓ) ℂ :=
  fun a a' => ∑ b : Cfg (N - ℓ) × Cfg (N - ℓ),
    Ψ ((splitLGT ℓ hℓ).symm (a, b)) * star (Ψ ((splitLGT ℓ hℓ).symm (a', b)))

theorem rdmLGT_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmLGT ℓ hℓ Ψ).IsHermitian := by
  ext a a'
  simp only [Matrix.conjTranspose_apply, rdmLGT, star_sum, star_mul', star_star]
  exact Finset.sum_congr rfl fun b _ => mul_comm _ _

/-- The `σˣ`-product gauge state `|ξ⟩_σ` with Fock matter `|n⟩`. -/
def xprodInit (ξ n : Cfg N) : Cfg N × Cfg N → ℂ :=
  prodInit (fun j x => if x = ξ j then 1 else 0) (Pi.single n 1)

namespace CutLocalAux

/-! ### Splitting configurations -/

/-- Entries of a recombined configuration in the left segment. -/
lemma splitCfg_symm_apply_lt (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) (b : Cfg (N - ℓ)) (k : Fin N)
    (hk : k.val < ℓ) : (splitCfg ℓ hℓ).symm (a, b) k = a ⟨k.val, hk⟩ := by
  simp only [splitCfg, Equiv.symm_trans_apply, Equiv.arrowCongr_symm, Equiv.arrowCongr_apply,
    Equiv.symm_symm, Equiv.refl_symm, Equiv.coe_refl, Function.comp_apply, id]
  rw [show (finCongr (by omega : N = ℓ + (N - ℓ))) k = Fin.castAdd (N - ℓ) ⟨k.val, hk⟩ from rfl,
    finSumFinEquiv_symm_apply_castAdd]
  rfl

/-- Entries of a recombined configuration in the right part. -/
lemma splitCfg_symm_apply_ge (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) (b : Cfg (N - ℓ)) (k : Fin N)
    (hk : ℓ ≤ k.val) : (splitCfg ℓ hℓ).symm (a, b) k = b ⟨k.val - ℓ, by omega⟩ := by
  simp only [splitCfg, Equiv.symm_trans_apply, Equiv.arrowCongr_symm, Equiv.arrowCongr_apply,
    Equiv.symm_symm, Equiv.refl_symm, Equiv.coe_refl, Function.comp_apply, id]
  have e : (finCongr (by omega : N = ℓ + (N - ℓ))) k = Fin.natAdd ℓ ⟨k.val - ℓ, by omega⟩ := by
    ext; simp; omega
  rw [e, finSumFinEquiv_symm_apply_natAdd]
  rfl

/-- The left part of a split configuration. -/
lemma splitCfg_apply_fst (ℓ : ℕ) (hℓ : ℓ ≤ N) (m : Cfg N) (i : Fin ℓ) :
    ((splitCfg ℓ hℓ) m).1 i = m ⟨i.val, by omega⟩ := by
  simp [splitCfg, Equiv.arrowCongr]
  rfl

/-- The right part of a split configuration. -/
lemma splitCfg_apply_snd (ℓ : ℕ) (hℓ : ℓ ≤ N) (m : Cfg N) (i : Fin (N - ℓ)) :
    ((splitCfg ℓ hℓ) m).2 i = m ⟨ℓ + i.val, by omega⟩ := by
  simp [splitCfg, Equiv.arrowCongr]
  rfl

/-- Recombining an LGT configuration. -/
lemma splitLGT_symm_apply (ℓ : ℕ) (hℓ : ℓ ≤ N) (α β : Cfg ℓ) (γ δ : Cfg (N - ℓ)) :
    (splitLGT ℓ hℓ).symm ((α, β), (γ, δ)) =
      ((splitCfg ℓ hℓ).symm (α, γ), (splitCfg ℓ hℓ).symm (β, δ)) := rfl

/-! ### Parity conservation -/

/-- Particle-number parity of a configuration. -/
def par (m : Cfg N) : ZMod 2 := ∑ k, m k

/-- Emptying an occupied mode changes the parity. -/
lemma par_eq_update_add_one (m : Cfg N) (y : Fin N) (hy : m y = 1) :
    par m = par (Function.update m y 0) + 1 := by
  unfold par
  rw [← Finset.add_sum_erase _ _ (mem_univ y), ← Finset.add_sum_erase _ _ (mem_univ y),
    Function.update_self, hy]
  have : ∑ k ∈ univ.erase y, Function.update m y 0 k = ∑ k ∈ univ.erase y, m k :=
    Finset.sum_congr rfl fun k hk => Function.update_of_ne (Finset.ne_of_mem_erase hk) _ _
  rw [this]
  ring

/-- `c_x† c_y` conserves the parity. -/
lemma cre_mul_ann_eq_zero (x y : Fin N) (m' m : Cfg N) (h : par m' ≠ par m) :
    (cre x * ann y) m' m = 0 := by
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann]
  refine Finset.sum_eq_zero fun p _ => ?_
  by_cases h1 : m' x = 1 ∧ p = Function.update m' x 0
  · by_cases h2 : m y = 1 ∧ p = Function.update m y 0
    · exfalso
      apply h
      rw [par_eq_update_add_one m' x h1.1, par_eq_update_add_one m y h2.1, ← h1.2, ← h2.2]
    · simp [h2]
  · simp [h1]

/-- Quadratic Hamiltonians conserve the parity. -/
lemma dGamma_eq_zero (A : Matrix (Fin N) (Fin N) ℝ) (m' m : Cfg N) (h : par m' ≠ par m) :
    dGamma A m' m = 0 := by
  simp only [dGamma, Matrix.sum_apply, Matrix.smul_apply, cre_mul_ann_eq_zero _ _ m' m h,
    smul_zero, Finset.sum_const_zero]

/-- The free-fermion propagator conserves the parity. -/
lemma evo_dGamma_eq_zero (A : Matrix (Fin N) (Fin N) ℝ) (t : ℝ) (m' m : Cfg N)
    (h : par m' ≠ par m) : evo (dGamma A) t m' m = 0 := by
  set X : Matrix (Cfg N) (Cfg N) ℂ := (-(t : ℂ) * Complex.I) • cplx (dGamma A) with hX
  set d : Cfg N → ℂ := fun m => ((par m).val : ℂ) with hd
  have hdne : ∀ m₁ m₂ : Cfg N, par m₁ ≠ par m₂ → d m₁ ≠ d m₂ := by
    intro m₁ m₂ h12 h'
    exact h12 (ZMod.val_injective _ (Nat.cast_injective h'))
  have hcomm : Commute X (diagonal d) := by
    show X * diagonal d = diagonal d * X
    ext a b
    rw [mul_diagonal, diagonal_mul]
    by_cases hab : par a = par b
    · simp only [hd, hab, mul_comm]
    · simp [hX, cplx, dGamma_eq_zero A a b hab]
  have h2 := congrFun (congrFun hcomm.exp_left.eq m') m
  rw [mul_diagonal, diagonal_mul] at h2
  have h3 : exp X m' m * (d m - d m') = 0 := by rw [mul_sub, h2]; ring
  rcases mul_eq_zero.mp h3 with h4 | h4
  · exact h4
  · exact absurd (sub_eq_zero.mp h4).symm (hdne m' m h)

/-- Applying a matrix to a basis vector. -/
lemma mulVec_single_one_apply {ι : Type*} [Fintype ι] [DecidableEq ι] (M : Matrix ι ι ℂ)
    (i j : ι) : (M *ᵥ Pi.single j 1) i = M i j := by
  simp [Matrix.mulVec, dotProduct, Pi.single_apply]

/-- The support of an evolved Fock state has the parity of the initial configuration. -/
lemma par_eq_of_ne_zero (A : Matrix (Fin N) (Fin N) ℝ) (t : ℝ) (n m : Cfg N)
    (hm : (evo (dGamma A) t *ᵥ Pi.single n 1) m ≠ 0) : par m = par n := by
  by_contra hc
  rw [mulVec_single_one_apply, evo_dGamma_eq_zero A t m n hc] at hm
  exact hm rfl

/-! ### Locality of the dual link labels -/

/-- The link labels left of `ℓ` only depend on the occupations left of `ℓ`. -/
lemma linkOf_head_eq (ℓ : ℕ) (m₁ m₂ q : Cfg N) (hA : ∀ k : Fin N, k.val < ℓ → m₁ k = m₂ k)
    (j : Fin N) (hj : j.val < ℓ) : linkOf m₁ q j = linkOf m₂ q j := by
  unfold linkOf
  refine Finset.sum_congr rfl fun k hk => ?_
  have hkj : k ≤ j := (Finset.mem_filter.1 hk).2
  rw [hA k (lt_of_le_of_lt (Fin.le_def.1 hkj) hj)]

/-- The link labels right of `ℓ` only depend on the occupations right of `ℓ` and the total
parity. -/
lemma linkOf_tail_eq (ℓ : ℕ) (m₁ m₂ q : Cfg N) (hB : ∀ k : Fin N, ℓ ≤ k.val → m₁ k = m₂ k)
    (hp : par m₁ = par m₂) (j : Fin N) (hj : ℓ ≤ j.val) : linkOf m₁ q j = linkOf m₂ q j := by
  unfold linkOf
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib]
  congr 1
  have e1 := Finset.sum_filter_add_sum_filter_not univ (fun k => k ≤ j) m₁
  have e2 := Finset.sum_filter_add_sum_filter_not univ (fun k => k ≤ j) m₂
  have e3 : ∑ k ∈ univ.filter (fun k => ¬ k ≤ j), m₁ k =
      ∑ k ∈ univ.filter (fun k => ¬ k ≤ j), m₂ k := by
    refine Finset.sum_congr rfl fun k hk => hB k ?_
    have := (Finset.mem_filter.1 hk).2
    rw [not_le, Fin.lt_def] at this
    omega
  unfold par at hp
  rw [← e1, ← e2, e3] at hp
  exact add_right_cancel hp

/-! ### The evolved definite-charge state -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- The initial state is a basis vector. -/
lemma xprodInit_eq (ξ n : Cfg N) : xprodInit ξ n = Pi.single (n, ξ) 1 := by
  funext p
  obtain ⟨n', x'⟩ := p
  simp only [xprodInit, prodInit, prodLink]
  rw [Finset.prod_boole]
  by_cases h1 : n' = n <;> by_cases h2 : x' = ξ
  · subst h1; subst h2; simp
  · have : ¬ ∀ i ∈ univ, x' i = ξ i := fun h => h2 (funext fun i => h i (mem_univ i))
    simp [h1, this, Pi.single_apply, h2]
  · simp [h1, Pi.single_apply]
  · simp [h1, Pi.single_apply]

/-- In the dual frame, the evolved state lives in the single charge sector `q = chargeOf n ξ`. -/
lemma evo_xprodInit_apply' (J h : ℝ) (ξ n : Cfg N) (t : ℝ) (n' x' : Cfg N) :
    (evo (H J h) t *ᵥ xprodInit ξ n) (n', x') =
      if x' = linkOf n' (chargeOf n ξ) then evo (Hsec J h (chargeOf n ξ)) t n' n else 0 := by
  rw [xprodInit_eq, mulVec_single_one_apply]
  have hU : (evo (H (N := N) J h) t).submatrix dualEquiv dualEquiv =
      blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) := by
    rw [evo_submatrix, duality, ← evo_blockDiagonal]
    rfl
  have e1 : evo (H J h) t (n', x') (n, ξ) = (evo (H (N := N) J h) t).submatrix dualEquiv
      dualEquiv (n', chargeOf n' x') (n, chargeOf n ξ) := by
    simp [linkOf_chargeOf]
  rw [e1, hU, blockDiagonal_apply']
  by_cases hc : chargeOf n' x' = chargeOf n ξ
  · have hx : x' = linkOf n' (chargeOf n ξ) := by rw [← hc, linkOf_chargeOf]
    simp only [hc, hx, ite_true]
  · have hx : x' ≠ linkOf n' (chargeOf n ξ) := fun hx => hc (by rw [hx, chargeOf_linkOf])
    simp only [hc, hx, ite_false]

/-- The exponential of a scalar matrix is a scalar matrix. -/
lemma exp_smul_one {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ℂ) :
    exp (d • (1 : Matrix ι ι ℂ)) = exp d • (1 : Matrix ι ι ℂ) := by
  rw [smul_one_eq_diagonal, smul_one_eq_diagonal, Matrix.exp_diagonal]
  congr 1
  ext i
  exact Pi.coe_exp _ i

/-- A real constant shift of the generator multiplies the propagator by a phase. -/
lemma evo_sub_smul_one {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) (c t : ℝ) :
    ∃ w : ℂ, w * star w = 1 ∧ evo (A - c • (1 : Matrix ι ι ℝ)) t = w • evo A t := by
  have hc : cplx (A - c • (1 : Matrix ι ι ℝ)) = cplx A + (-(c : ℂ)) • (1 : Matrix ι ι ℂ) := by
    ext i j
    by_cases h : i = j
    · subst h; simp [cplx, sub_eq_add_neg]
    · simp [cplx, one_apply_ne h]
  refine ⟨Complex.exp (((t * c : ℝ) : ℂ) * Complex.I), ?_, ?_⟩
  · rw [Complex.star_def, Complex.mul_conj, Complex.normSq_eq_norm_sq,
      Complex.norm_exp_ofReal_mul_I]
    norm_num
  · unfold evo
    rw [hc, smul_add, smul_smul, Matrix.exp_add_of_commute _ _ ?_]
    · rw [exp_smul_one, mul_smul_one, Complex.exp_eq_exp_ℂ]
      congr 2
      push_cast
      ring
    · exact (Commute.one_right _).smul_right _

/-- The evolved state in the `(occupations, link labels)` basis:
`Ψ(t)(n', x') = [x' = x(n', q)] w φ(n')` with `|w| = 1` and `φ = e^{-itdΓ(h_q)} |n⟩`. -/
lemma evo_xprodInit_apply (J h : ℝ) (ξ n : Cfg N) (t : ℝ) :
    ∃ w : ℂ, w * star w = 1 ∧ ∀ n' x' : Cfg N,
      (evo (H J h) t *ᵥ xprodInit ξ n) (n', x') =
        if x' = linkOf n' (chargeOf n ξ) then
          w * (evo (dGamma (hq J h (chargeOf n ξ))) t *ᵥ Pi.single n 1) n' else 0 := by
  obtain ⟨w, hw, hevo⟩ := evo_sub_smul_one (dGamma (hq J h (chargeOf n ξ)))
    (h * ∑ j, sgn (chargeOf n ξ j)) t
  refine ⟨w, hw, fun n' x' => ?_⟩
  rw [evo_xprodInit_apply', mulVec_single_one_apply, Hsec, hevo, Matrix.smul_apply, smul_eq_mul]

/-! ### Isometries and von Neumann entropy -/

/-- The isometric embedding `α ↦ (α, g α)` as a matrix. -/
def isoMat {ι ι' : Type*} [DecidableEq ι] [DecidableEq ι'] (g : ι → ι') : Matrix (ι × ι') ι ℂ :=
  fun p i => if p.1 = i ∧ p.2 = g i then 1 else 0

variable {ι ι' : Type*} [Fintype ι] [DecidableEq ι] [Fintype ι'] [DecidableEq ι']

/-- `isoMat g` is an isometry. -/
lemma isoMat_conjTranspose_mul (g : ι → ι') : (isoMat g)ᴴ * isoMat g = 1 := by
  ext i j
  simp only [Matrix.mul_apply, conjTranspose_apply, isoMat, Fintype.sum_prod_type]
  by_cases hij : i = j
  · subst hij
    simp [ite_and, Finset.sum_ite_eq', Finset.sum_ite_eq]
  · rw [Matrix.one_apply_ne hij]
    refine Finset.sum_eq_zero fun a _ => Finset.sum_eq_zero fun b _ => ?_
    by_cases h1 : a = i ∧ b = g i
    · have h2 : ¬ (a = j ∧ b = g j) := fun h2 => hij (h1.1.symm.trans h2.1)
      simp [h2]
    · simp [h1]

/-- Entries of `V ρ V†` for `V = isoMat g`. -/
lemma isoMat_conj_apply (g : ι → ι') (ρ : Matrix ι ι ℂ) (a : ι) (b : ι') (a' : ι) (b' : ι') :
    (isoMat g * ρ * (isoMat g)ᴴ) (a, b) (a', b') =
      if b = g a ∧ b' = g a' then ρ a a' else 0 := by
  have h1 : ∀ j, (isoMat g * ρ) (a, b) j = if b = g a then ρ a j else 0 := by
    intro j
    simp only [Matrix.mul_apply, isoMat]
    rw [Finset.sum_eq_single a]
    · by_cases hb : b = g a <;> simp [hb]
    · intro i _ hi
      simp [Ne.symm hi]
    · simp
  rw [Matrix.mul_apply]
  simp only [h1, conjTranspose_apply, isoMat]
  rw [Finset.sum_eq_single a']
  · by_cases hb : b = g a <;> by_cases hb' : b' = g a' <;> simp [hb, hb']
  · intro j _ hj
    simp [Ne.symm hj]
  · simp

/-- Equal matrices have equal von Neumann entropies. -/
lemma vonNeumann_congr {κ : Type*} [Fintype κ] [DecidableEq κ] {ρ σ : Matrix κ κ ℂ}
    (hρ : ρ.IsHermitian) (hσ : σ.IsHermitian) (h : ρ = σ) : vonNeumann ρ hρ = vonNeumann σ hσ := by
  subst h
  rfl

/-- If the characteristic polynomials agree up to powers of `X` (i.e. the nonzero eigenvalues
agree with multiplicity), then the von Neumann entropies agree. -/
lemma vonNeumann_eq_of_charpoly {κ : Type*} [Fintype κ] [DecidableEq κ] (σ : Matrix κ κ ℂ)
    (hσ : σ.IsHermitian) (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (a b : ℕ)
    (h : Polynomial.X ^ a * σ.charpoly = Polynomial.X ^ b * ρ.charpoly) :
    vonNeumann σ hσ = vonNeumann ρ hρ := by
  have h1 : (Polynomial.X ^ a * σ.charpoly : Polynomial ℂ) ≠ 0 :=
    mul_ne_zero (pow_ne_zero _ Polynomial.X_ne_zero) σ.charpoly_monic.ne_zero
  have h2 : (Polynomial.X ^ b * ρ.charpoly : Polynomial ℂ) ≠ 0 :=
    mul_ne_zero (pow_ne_zero _ Polynomial.X_ne_zero) ρ.charpoly_monic.ne_zero
  have hr := congrArg Polynomial.roots h
  rw [Polynomial.roots_mul h1, Polynomial.roots_mul h2, Polynomial.roots_X_pow,
    Polynomial.roots_X_pow, hσ.roots_charpoly_eq_eigenvalues,
    hρ.roots_charpoly_eq_eigenvalues] at hr
  have := congrArg (fun s : Multiset ℂ => (s.map (fun z => Real.negMulLog (RCLike.re z))).sum) hr
  simpa [vonNeumann, Multiset.map_nsmul, Real.negMulLog_zero, Finset.sum_eq_multiset_sum,
    Function.comp_def] using this

/-- Conjugating by an isometry does not change the von Neumann entropy. -/
lemma vonNeumann_conj_isometry {κ : Type*} [Fintype κ] [DecidableEq κ] (V : Matrix κ ι ℂ)
    (hV : Vᴴ * V = 1) (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (hσ : (V * ρ * Vᴴ).IsHermitian) :
    vonNeumann (V * ρ * Vᴴ) hσ = vonNeumann ρ hρ := by
  refine vonNeumann_eq_of_charpoly _ hσ ρ hρ (Fintype.card ι) (Fintype.card κ) ?_
  have := Matrix.charpoly_mul_comm' V (ρ * Vᴴ)
  rwa [← Matrix.mul_assoc, Matrix.mul_assoc ρ, hV, Matrix.mul_one] at this

/-! ### The reduced density matrix of the definite-charge state -/

/-- Sum over the gauge labels of the complement: only one label contributes, and on the support
the two conditions coincide. -/
lemma sum_delta {κ : Type*} [Fintype κ] [DecidableEq κ] (P Q : Prop) [Decidable P]
    [Decidable Q] (L₁ L₂ : κ) (u v w : ℂ) (hw : w * star w = 1)
    (hL : u * star v ≠ 0 → L₁ = L₂) :
    ∑ δ : κ, (if P ∧ δ = L₁ then w * u else 0) * star (if Q ∧ δ = L₂ then w * v else 0) =
      if P ∧ Q then u * star v else 0 := by
  have hwuv : w * u * star (w * v) = u * star v := by
    rw [star_mul', show w * u * (star w * star v) = (w * star w) * (u * star v) by ring, hw,
      one_mul]
  by_cases huv : u * star v = 0
  · have h0 : ∀ δ : κ,
        (if P ∧ δ = L₁ then w * u else 0) * star (if Q ∧ δ = L₂ then w * v else 0) = 0 := by
      intro δ
      split_ifs <;> simp [hwuv, huv]
    simp only [h0, Finset.sum_const_zero, huv, ite_self]
  · obtain rfl := hL huv
    rw [Finset.sum_eq_single L₁]
    · by_cases hP : P <;> by_cases hQ : Q <;> simp [hP, hQ, hwuv]
    · intro δ _ hδ
      simp [hδ]
    · simp

/-- The occupations of `A` determine the link labels of `A`: `L_A(α)`. -/
def linkA (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) (α : Cfg ℓ) : Cfg ℓ :=
  ((splitCfg ℓ hℓ) (linkOf ((splitCfg ℓ hℓ).symm (α, 0)) q)).1

/-- The link labels of `A` only depend on the occupations of `A`. -/
lemma fst_splitCfg_linkOf (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) (α : Cfg ℓ) (γ : Cfg (N - ℓ)) :
    ((splitCfg ℓ hℓ) (linkOf ((splitCfg ℓ hℓ).symm (α, γ)) q)).1 = linkA ℓ hℓ q α := by
  funext i
  rw [linkA, splitCfg_apply_fst, splitCfg_apply_fst]
  refine linkOf_head_eq ℓ _ _ q (fun k hk => ?_) _ (by simp)
  rw [splitCfg_symm_apply_lt _ _ _ _ k hk, splitCfg_symm_apply_lt _ _ _ _ k hk]

/-- On the support, the link labels of the complement agree. -/
lemma snd_splitCfg_linkOf (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) (α α' : Cfg ℓ) (γ : Cfg (N - ℓ))
    (hp : par ((splitCfg ℓ hℓ).symm (α, γ)) = par ((splitCfg ℓ hℓ).symm (α', γ))) :
    ((splitCfg ℓ hℓ) (linkOf ((splitCfg ℓ hℓ).symm (α, γ)) q)).2 =
      ((splitCfg ℓ hℓ) (linkOf ((splitCfg ℓ hℓ).symm (α', γ)) q)).2 := by
  funext i
  rw [splitCfg_apply_snd, splitCfg_apply_snd]
  refine linkOf_tail_eq ℓ _ _ q (fun k hk => ?_) hp _ (by simp)
  rw [splitCfg_symm_apply_ge _ _ _ _ k hk, splitCfg_symm_apply_ge _ _ _ _ k hk]

/-- **Reduced density matrix of the definite-charge state**: `ρ_LGT = V ρ_A V†` with the isometry
`V : α ↦ (α, L_A(α))` and `ρ_A` the reduced density matrix of `e^{-itdΓ(h_q)} |n⟩`. -/
lemma rdmLGT_evo_xprod (J h : ℝ) (ξ n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n) =
      isoMat (linkA ℓ hℓ (chargeOf n ξ)) *
        rdm ℓ hℓ (evo (dGamma (hq J h (chargeOf n ξ))) t *ᵥ Pi.single n 1) *
          (isoMat (linkA ℓ hℓ (chargeOf n ξ)))ᴴ := by
  set q := chargeOf n ξ with hq_def
  set φ := evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1 with hφ
  obtain ⟨w, hw, hΨ⟩ := evo_xprodInit_apply J h ξ n t
  rw [← hq_def, ← hφ] at hΨ
  set sym := (splitCfg ℓ hℓ).symm with hsym
  -- the evolved state, in split coordinates
  have hΨ' : ∀ (α β : Cfg ℓ) (γ δ : Cfg (N - ℓ)),
      (evo (H J h) t *ᵥ xprodInit ξ n) ((splitLGT ℓ hℓ).symm ((α, β), (γ, δ))) =
        if β = linkA ℓ hℓ q α ∧ δ = ((splitCfg ℓ hℓ) (linkOf (sym (α, γ)) q)).2 then
          w * φ (sym (α, γ)) else 0 := by
    intro α β γ δ
    rw [splitLGT_symm_apply, hΨ]
    refine if_congr ?_ rfl rfl
    rw [hsym, Equiv.symm_apply_eq, Prod.ext_iff, fst_splitCfg_linkOf]
  ext ⟨α, β⟩ ⟨α', β'⟩
  rw [isoMat_conj_apply]
  simp only [rdmLGT, Fintype.sum_prod_type, hΨ']
  have hγ : ∀ γ : Cfg (N - ℓ),
      ∑ δ : Cfg (N - ℓ),
        (if β = linkA ℓ hℓ q α ∧ δ = ((splitCfg ℓ hℓ) (linkOf (sym (α, γ)) q)).2 then
          w * φ (sym (α, γ)) else 0) *
        star (if β' = linkA ℓ hℓ q α' ∧ δ = ((splitCfg ℓ hℓ) (linkOf (sym (α', γ)) q)).2 then
          w * φ (sym (α', γ)) else 0) =
      if β = linkA ℓ hℓ q α ∧ β' = linkA ℓ hℓ q α' then
        φ (sym (α, γ)) * star (φ (sym (α', γ))) else 0 := by
    intro γ
    refine sum_delta _ _ _ _ _ _ w hw fun huv => ?_
    have hu : φ (sym (α, γ)) ≠ 0 := left_ne_zero_of_mul huv
    have hv : φ (sym (α', γ)) ≠ 0 := fun h0 => huv (by rw [h0, star_zero, mul_zero])
    exact snd_splitCfg_linkOf ℓ hℓ q α α' γ
      ((par_eq_of_ne_zero _ t n _ hu).trans (par_eq_of_ne_zero _ t n _ hv).symm)
  simp only [hγ]
  by_cases hc : β = linkA ℓ hℓ q α ∧ β' = linkA ℓ hℓ q α'
  · simp only [hc, and_self, ite_true, rdm]
  · simp only [hc, ite_false, Finset.sum_const_zero]

end CutLocalAux

open CutLocalAux in
/-- **Definite-charge states** (paper: Lemma `lem:cutlocal`(c), Theorem `thm:arealaw`(c)): the full
matter + gauge entanglement of the evolved LGT state equals the entanglement of the free-fermion
state of the sector `q = chargeOf n ξ`. -/
theorem vonNeumann_rdmLGT_xprod (J h : ℝ) (ξ n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h (chargeOf n ξ))) t *ᵥ Pi.single n 1))
        (rdm_isHermitian ℓ hℓ _) := by
  have hEq := rdmLGT_evo_xprod J h ξ n t ℓ hℓ
  have hσ := hEq ▸ rdmLGT_isHermitian ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)
  rw [vonNeumann_congr _ hσ hEq]
  exact vonNeumann_conj_isometry _ (isoMat_conjTranspose_mul _) _ _ _

/-- The Fock configuration of the occupied set of `n` is `n`. -/
private lemma fockCfg_occupied (n : Cfg N) : fockCfg (univ.filter (fun x => n x = 1)) = n := by
  funext x
  unfold fockCfg
  rcases CutLocalAux.zmod2_cases (n x) with h0 | h0 <;> simp [h0]

/-- **Area law for the full entanglement of definite-charge states** (paper: Theorem
`thm:arealaw`(c)): under Hypothesis (EC), the average over uniformly random link labels `ξ` of the
full entanglement entropy is bounded uniformly in time, volume, `ℓ` and the Fock state. -/
theorem lgt_areaLaw_full (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ),
      (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) ≤
          K := by
  obtain ⟨K, hK⟩ := lgt_areaLaw J h hEC
  refine ⟨K, fun N ℓ hℓ n t => ?_⟩
  -- the bijection `ξ ↦ chargeOf n ξ`
  let e : Cfg N ≃ Cfg N := ⟨chargeOf n, linkOf n, linkOf_chargeOf n, chargeOf_linkOf n⟩
  have h1 : ∑ ξ : Cfg N,
      vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      ∑ q : Cfg N, vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ
        fockVec (univ.filter (fun x => n x = 1)))) (rdm_isHermitian ℓ hℓ _) := by
    rw [fockVec, fockCfg_occupied]
    simp only [vonNeumann_rdmLGT_xprod]
    exact Equiv.sum_comp e (fun q => vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ
      Pi.single n 1)) (rdm_isHermitian ℓ hℓ _))
  have h2 := hK N ℓ hℓ (univ.filter (fun x => n x = 1)) t
  unfold avgQ at h2
  rw [h1]
  exact h2

end DFL
