import DFL.SlaterGeneral
import DFL.AnnealedInterval

/-!
# Entanglement of Slater determinants for arbitrary regions (paper: Lemma `lem:slater`(a),(b))

For a region `R` (e.g. a bulk interval) the paper identifies `𝓕_Λ ≅ 𝓕_R ⊗ 𝓕_B`, `B = Λ \ R`, by
`c(a) ↦ c(a) ⊗ 1`, `c(b) ↦ (-1)^{N_R} ⊗ c(b)`. In the Jordan–Wigner occupation basis this maps `|n⟩` to
`ε(α, β) |α⟩ ⊗ |β⟩`, where `α = n|_R`, `β = n|_B` and `ε(α, β) = (-1)^{#{(a, b) : a ∈ R, b ∉ R, b < a, α_a = β_b = 1}}`
(`reorderSign`). The reduced density matrix of the modes in `R` (`rdmF`) is therefore
`ρ_R(α, α') = ∑_β ε(α, β) ε(α', β) ψ(α ⊔ β) conj ψ(α' ⊔ β)`.

* `vonNeumann_rdmF_slater`: for `ψ = Γ(e^{-itA})|S⟩`, `S(ρ_R) = ∑_i h(ν_i)`, `ν_i` the eigenvalues of `C_R`;
* `vonNeumann_rdmF_le`: `S(ρ_R) ≤ 2 ln 2 ∑_{x ∈ R, y ∉ R} |K_{xy}|`;
* `vonNeumann_rdmF_segment`: for a left segment this is the reduced state `rdm` of `DFL.AreaLaw`;
* `lgt_areaLaw_bulk`: for an interval `B = [a, b]`, `E_q sup_t S(ρ^q_B(t)) ≤ 4 ln 2 Θ(F)` under (EC)
  (finite-volume form of Theorem `thm:infvol`(d)).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Splitting a configuration into its restrictions to `R` and to the complement of `R`. -/
def splitR (R : Finset (Fin N)) : Cfg N ≃ (R → ZMod 2) × ({x : Fin N // x ∉ R} → ZMod 2) :=
  Equiv.piEquivPiSubtypeProd (fun x : Fin N => x ∈ R) (fun _ => ZMod 2)

/-- The reordering sign `ε(α, β) = (-1)^{#{(a, b) : a ∈ R, b ∉ R, b < a, α_a = β_b = 1}}`. -/
def reorderSign (R : Finset (Fin N)) (α : R → ZMod 2) (β : {x : Fin N // x ∉ R} → ZMod 2) : ℝ :=
  (-1 : ℝ) ^ ((univ : Finset (R × {x : Fin N // x ∉ R})).filter
    (fun p => (p.2 : Fin N) < (p.1 : Fin N) ∧ α p.1 = 1 ∧ β p.2 = 1)).card

/-- The reduced density matrix of the modes in `R`:
`ρ_R(α, α') = ∑_β ε(α, β) ε(α', β) ψ(α ⊔ β) conj ψ(α' ⊔ β)`. -/
def rdmF (R : Finset (Fin N)) (v : Cfg N → ℂ) : Matrix (R → ZMod 2) (R → ZMod 2) ℂ :=
  fun α α' => ∑ β : {x : Fin N // x ∉ R} → ZMod 2,
    ((reorderSign R α β * reorderSign R α' β : ℝ) : ℂ) *
      v ((splitR R).symm (α, β)) * star (v ((splitR R).symm (α', β)))

/-- `ρ_R` is Hermitian. -/
theorem rdmF_isHermitian (R : Finset (Fin N)) (v : Cfg N → ℂ) : (rdmF R v).IsHermitian := by
  ext α α'
  simp only [Matrix.conjTranspose_apply, rdmF, star_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [star_mul', star_mul', star_star, Complex.star_def, Complex.conj_ofReal]
  push_cast
  ring

namespace SlaterBulkAux

open AreaLawAux SlaterExactAux SlaterGeneralAux

/-! ### Restrictions of configurations and the reordering parity -/

/-- The configuration `n` with the modes outside `R` emptied. -/
def cutIn (R : Finset (Fin N)) (n : Cfg N) : Cfg N := fun x => if x ∈ R then n x else 0

/-- The configuration `n` with the modes in `R` emptied. -/
def cutOut (R : Finset (Fin N)) (n : Cfg N) : Cfg N := fun x => if x ∈ R then 0 else n x

/-- The reordering parity `∑_{a ∈ R, b ∉ R, b < a} n_a n_b ∈ ZMod 2`. -/
def rpar (R : Finset (Fin N)) (n : Cfg N) : ZMod 2 :=
  ∑ x ∈ R, ∑ y ∈ univ \ R, if y < x then n x * n y else 0

/-- The fermionic product vector `(X ⋆ Y)(n) = ε(n) X(n|_R) Y(n|_{R^c})` (with `ε` the reordering
sign), the image of `X ⊗ Y` under `𝓕_R ⊗ 𝓕_B ≅ 𝓕_Λ`. -/
def ftens (R : Finset (Fin N)) (X Y : Cfg N → ℂ) : Cfg N → ℂ :=
  fun n => ((sgn (rpar R n) : ℝ) : ℂ) * X (cutIn R n) * Y (cutOut R n)

/-- Emptying a mode of `R` commutes with `cutIn`. -/
lemma cutIn_update_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∈ R) (v : ZMod 2) :
    cutIn R (Function.update n x v) = Function.update (cutIn R n) x v := by
  funext k
  by_cases hk : k = x
  · subst hk; simp [cutIn, hx]
  · simp only [cutIn, Function.update_of_ne hk]

/-- Changing a mode of `R` does not change `cutOut`. -/
lemma cutOut_update_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∈ R) (v : ZMod 2) :
    cutOut R (Function.update n x v) = cutOut R n := by
  funext k
  by_cases hk : k = x
  · subst hk; simp [cutOut, hx]
  · simp only [cutOut, Function.update_of_ne hk]

/-- Changing a mode outside `R` does not change `cutIn`. -/
lemma cutIn_update_not_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∉ R)
    (v : ZMod 2) : cutIn R (Function.update n x v) = cutIn R n := by
  funext k
  by_cases hk : k = x
  · subst hk; simp [cutIn, hx]
  · simp only [cutIn, Function.update_of_ne hk]

/-- Emptying a mode outside `R` commutes with `cutOut`. -/
lemma cutOut_update_not_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∉ R)
    (v : ZMod 2) : cutOut R (Function.update n x v) = Function.update (cutOut R n) x v := by
  funext k
  by_cases hk : k = x
  · subst hk; simp [cutOut, hx]
  · simp only [cutOut, Function.update_of_ne hk]

/-- The Jordan–Wigner string splits into the contributions of `R` and of its complement. -/
lemma jwPar_cut (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) :
    jwPar n x = jwPar (cutIn R n) x + jwPar (cutOut R n) x := by
  unfold jwPar
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  by_cases hk : k ∈ R <;> simp [cutIn, cutOut, hk]

/-- Removing a particle at `x ∈ R` changes the reordering parity by the number of particles of the
complement to the left of `x`. -/
lemma rpar_update_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∈ R) (h1 : n x = 1) :
    rpar R n = rpar R (Function.update n x 0) + jwPar (cutOut R n) x := by
  unfold rpar
  rw [← Finset.add_sum_erase _ _ hx, ← Finset.add_sum_erase _ _ hx]
  have e1 : ∑ y ∈ univ \ R, (if y < x then Function.update n x 0 x * Function.update n x 0 y
      else 0) = 0 := by simp
  have e2 : ∑ a ∈ R.erase x, ∑ y ∈ univ \ R,
      (if y < a then Function.update n x 0 a * Function.update n x 0 y else 0) =
      ∑ a ∈ R.erase x, ∑ y ∈ univ \ R, (if y < a then n a * n y else 0) := by
    refine Finset.sum_congr rfl fun a ha => Finset.sum_congr rfl fun y hy => ?_
    have hax : a ≠ x := Finset.ne_of_mem_erase ha
    have hyx : y ≠ x := fun h => (Finset.mem_sdiff.mp hy).2 (h ▸ hx)
    rw [Function.update_of_ne hax, Function.update_of_ne hyx]
  have e3 : ∑ y ∈ univ \ R, (if y < x then n x * n y else 0) = jwPar (cutOut R n) x := by
    have hs : (univ \ R : Finset (Fin N)) = univ.filter (fun y => y ∉ R) := by ext; simp
    unfold jwPar
    rw [hs, Finset.sum_filter, Finset.sum_filter, h1]
    refine Finset.sum_congr rfl fun k _ => ?_
    by_cases hk : k ∈ R <;> by_cases hkx : k < x <;> simp [cutOut, hk, hkx]
  rw [e1, e2, e3]
  ring

/-- Removing a particle at `x ∉ R` changes the reordering parity by the number of particles of `R`
to the right of `x`. -/
lemma rpar_update_not_mem (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∉ R)
    (h1 : n x = 1) :
    rpar R n = rpar R (Function.update n x 0) + ∑ a ∈ R, (if x < a then n a else 0) := by
  unfold rpar
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun a ha => ?_
  have hax : a ≠ x := fun h => hx (h ▸ ha)
  have hxR : x ∈ univ \ R := by simp [hx]
  rw [← Finset.add_sum_erase _ _ hxR, ← Finset.add_sum_erase _ _ hxR]
  have e2 : ∑ y ∈ (univ \ R).erase x,
      (if y < a then Function.update n x 0 a * Function.update n x 0 y else 0) =
      ∑ y ∈ (univ \ R).erase x, (if y < a then n a * n y else 0) := by
    refine Finset.sum_congr rfl fun y hy => ?_
    have hyx : y ≠ x := Finset.ne_of_mem_erase hy
    rw [Function.update_of_ne hax, Function.update_of_ne hyx]
  rw [e2, Function.update_self, Function.update_of_ne hax, h1]
  by_cases hxa : x < a <;> simp [hxa]
  ring

/-- For `x ∉ R`, the number of particles of `R` splits into those left and right of `x`. -/
lemma sum_cutIn (R : Finset (Fin N)) (n : Cfg N) (x : Fin N) (hx : x ∉ R) :
    ∑ i, cutIn R n i = jwPar (cutIn R n) x + ∑ a ∈ R, (if x < a then n a else 0) := by
  unfold jwPar
  rw [Finset.sum_filter, ← Fintype.sum_extend_by_zero R, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases hiR : i ∈ R
  · have hix : i ≠ x := fun h => hx (h ▸ hiR)
    rcases lt_or_gt_of_ne hix with h | h
    · simp [cutIn, hiR, h, not_lt_of_gt h]
    · simp [cutIn, hiR, h, not_lt_of_gt h]
  · simp [cutIn, hiR]

/-! ### Creation operators acting on fermionic product vectors -/

/-- A creation operator of a mode in `R` acts on the first factor. -/
lemma cre_ftens_mem (R : Finset (Fin N)) (x : Fin N) (hx : x ∈ R) (X Y : Cfg N → ℂ) :
    cplx (cre x) *ᵥ ftens R X Y = ftens R (cplx (cre x) *ᵥ X) Y := by
  funext n
  rw [cre_mulVec_apply]
  simp only [ftens, cre_mulVec_apply]
  have hcx : cutIn R n x = n x := by simp [cutIn, hx]
  rw [hcx]
  split_ifs with h
  · rw [cutIn_update_mem R n x hx, cutOut_update_mem R n x hx]
    have hs : sgn (jwPar n x) * sgn (rpar R (Function.update n x 0)) =
        sgn (rpar R n) * sgn (jwPar (cutIn R n) x) := by
      rw [← FreeDynAux.sgn_add, ← FreeDynAux.sgn_add, rpar_update_mem R n x hx h, jwPar_cut R n x]
      congr 1
      ring
    have hs' : ((sgn (jwPar n x) : ℝ) : ℂ) * ((sgn (rpar R (Function.update n x 0)) : ℝ) : ℂ) =
        ((sgn (rpar R n) : ℝ) : ℂ) * ((sgn (jwPar (cutIn R n) x) : ℝ) : ℂ) := by
      exact_mod_cast hs
    linear_combination (X (Function.update (cutIn R n) x 0) * Y (cutOut R n)) * hs'
  · simp

/-- A creation operator of a mode outside `R` acts on the second factor, with the parity of the
first one. -/
lemma cre_ftens_not_mem (R : Finset (Fin N)) (x : Fin N) (hx : x ∉ R) (X Y : Cfg N → ℂ) :
    cplx (cre x) *ᵥ ftens R X Y = ftens R (parMul X) (cplx (cre x) *ᵥ Y) := by
  funext n
  rw [cre_mulVec_apply]
  simp only [ftens, cre_mulVec_apply, parMul, parC]
  have hcx : cutOut R n x = n x := by simp [cutOut, hx]
  rw [hcx]
  split_ifs with h
  · rw [cutIn_update_not_mem R n x hx, cutOut_update_not_mem R n x hx]
    have hs : sgn (jwPar n x) * sgn (rpar R (Function.update n x 0)) =
        sgn (rpar R n) * (sgn (∑ i, cutIn R n i) * sgn (jwPar (cutOut R n) x)) := by
      rw [← FreeDynAux.sgn_add, ← FreeDynAux.sgn_add, ← FreeDynAux.sgn_add,
        rpar_update_not_mem R n x hx h, jwPar_cut R n x, sum_cutIn R n x hx]
      congr 1
      have h2 := CharTwo.add_self_eq_zero (∑ a ∈ R, if x < a then n a else 0)
      linear_combination -h2
    have hs' : ((sgn (jwPar n x) : ℝ) : ℂ) * ((sgn (rpar R (Function.update n x 0)) : ℝ) : ℂ) =
        ((sgn (rpar R n) : ℝ) : ℂ) * (((sgn (∑ i, cutIn R n i) : ℝ) : ℂ) *
          ((sgn (jwPar (cutOut R n) x) : ℝ) : ℂ)) := by
      exact_mod_cast hs
    linear_combination (X (cutIn R n) * Y (Function.update (cutOut R n) x 0)) * hs'
  · simp

/-- `⋆` is linear in the first factor. -/
lemma ftens_sum_smul_left (R : Finset (Fin N)) {ι : Type*} (s : Finset ι) (c : ι → ℂ)
    (X : ι → Cfg N → ℂ) (Y : Cfg N → ℂ) :
    ftens R (∑ i ∈ s, c i • X i) Y = ∑ i ∈ s, c i • ftens R (X i) Y := by
  funext n
  simp only [ftens, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- `⋆` is linear in the second factor. -/
lemma ftens_sum_smul_right (R : Finset (Fin N)) {ι : Type*} (s : Finset ι) (c : ι → ℂ)
    (X : Cfg N → ℂ) (Y : ι → Cfg N → ℂ) :
    ftens R X (∑ i ∈ s, c i • Y i) = ∑ i ∈ s, c i • ftens R X (Y i) := by
  funext n
  simp only [ftens, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- `(c X) ⋆ Y = c (X ⋆ Y)`. -/
lemma ftens_smul_left (R : Finset (Fin N)) (c : ℂ) (X Y : Cfg N → ℂ) :
    ftens R (c • X) Y = c • ftens R X Y := by
  funext n
  simp only [ftens, Pi.smul_apply, smul_eq_mul]
  ring

/-- `c†(w)` acting on a fermionic product vector. -/
lemma crOp_ftens (R : Finset (Fin N)) (w : Fin N → ℂ) (X Y : Cfg N → ℂ) :
    crOp w *ᵥ ftens R X Y =
      ftens R (crOp (fun x => if x ∈ R then w x else 0) *ᵥ X) Y +
        ftens R (parMul X) (crOp (fun x => if x ∈ R then 0 else w x) *ᵥ Y) := by
  rw [crOp_mulVec_eq, crOp_mulVec_eq, crOp_mulVec_eq, ftens_sum_smul_left, ftens_sum_smul_right,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun x _ => ?_
  by_cases hx : x ∈ R
  · simp [hx, cre_ftens_mem R x hx]
  · simp [hx, cre_ftens_not_mem R x hx]

/-- The vacuum is the product of the vacua. -/
lemma ftens_vac (R : Finset (Fin N)) :
    ftens R (fockVec ∅) (fockVec ∅) = fockVec (∅ : Finset (Fin N)) := by
  funext n
  have h0 : fockCfg (∅ : Finset (Fin N)) = 0 := by funext x; simp [fockCfg]
  simp only [ftens, fockVec, Pi.single_apply, h0]
  by_cases hn : n = 0
  · subst hn
    have e1 : cutIn R (0 : Cfg N) = 0 := by funext x; simp [cutIn]
    have e2 : cutOut R (0 : Cfg N) = 0 := by funext x; simp [cutOut]
    have e3 : rpar R (0 : Cfg N) = 0 := by simp [rpar]
    simp [e1, e2, e3]
  · rw [ite_eq_right hn]
    obtain ⟨x, hx⟩ : ∃ x, n x ≠ 0 := by
      by_contra hc
      push Not at hc
      exact hn (funext hc)
    by_cases hxR : x ∈ R
    · have : cutIn R n ≠ 0 := fun h => hx (by simpa [cutIn, hxR] using congrFun h x)
      simp [this]
    · have : cutOut R n ≠ 0 := fun h => hx (by simpa [cutOut, hxR] using congrFun h x)
      simp [this]

/-- **Expansion of a Slater determinant in the fermionic product basis**: writing
`w_k = P_R w_k + P_{R^c} w_k`,
`c†(w_1) ⋯ c†(w_M) |0⟩ = ∑_σ s_σ (∏_{σ k} c†(P_R w_k) |0⟩) ⋆ (∏_{¬σ k} c†(P_{R^c} w_k) |0⟩)`
with `|s_σ| = 1`. -/
lemma slater_expand_gen (R : Finset (Fin N)) {M : ℕ} (w : Fin M → Fin N → ℂ) :
    ∃ s : (Fin M → Bool) → ℂ, (∀ σ, ‖s σ‖ = 1) ∧
      slaterAlt N M w = ∑ σ, s σ • ftens R (prodSel (fun k x => if x ∈ R then w k x else 0) σ)
        (prodSel (fun k x => if x ∈ R then 0 else w k x) (fun k => !σ k)) := by
  induction M with
  | zero =>
    refine ⟨fun _ => 1, fun _ => by simp, ?_⟩
    rw [Fintype.sum_unique]
    simp [slaterAlt_apply, prodSel_zero, ftens_vac]
  | succ M ih =>
    obtain ⟨s', hs', h'⟩ := ih (fun k => w k.succ)
    refine ⟨fun σ => if σ 0 then s' (fun k : Fin M => σ k.succ) else
      s' (fun k : Fin M => σ k.succ) * ∏ k : Fin M, (if σ k.succ then (-1 : ℂ) else 1), ?_, ?_⟩
    · intro σ
      dsimp only
      split_ifs
      · exact hs' _
      · rw [norm_mul, hs', one_mul, norm_prod]
        exact Finset.prod_eq_one fun k _ => by split_ifs <;> simp
    · have h1 : slaterAlt N (M + 1) w = crOp (w 0) *ᵥ slaterAlt N M (fun k => w k.succ) := by
        rw [slaterAlt_apply, slaterAlt_apply, List.ofFn_succ, List.prod_cons,
          ← Matrix.mulVec_mulVec]
      rw [h1, h', Matrix.mulVec_sum, sum_fun_succ]
      refine Finset.sum_congr rfl fun τ _ => ?_
      have hXt : prodSel (fun k x => if x ∈ R then w k x else 0) (Fin.cons true τ) =
          crOp (fun x => if x ∈ R then w 0 x else 0) *ᵥ
            prodSel (fun k x => if x ∈ R then w k.succ x else 0) τ := by
        rw [prodSel_succ]; simp
      have hYt : prodSel (fun k x => if x ∈ R then 0 else w k x)
          (fun k => !(Fin.cons true τ : Fin (M + 1) → Bool) k) =
          prodSel (fun k x => if x ∈ R then 0 else w k.succ x) (fun k => !τ k) := by
        rw [prodSel_succ]; simp
      have hXf : prodSel (fun k x => if x ∈ R then w k x else 0) (Fin.cons false τ) =
          prodSel (fun k x => if x ∈ R then w k.succ x else 0) τ := by
        rw [prodSel_succ]; simp
      have hYf : prodSel (fun k x => if x ∈ R then 0 else w k x)
          (fun k => !(Fin.cons false τ : Fin (M + 1) → Bool) k) =
          crOp (fun x => if x ∈ R then 0 else w 0 x) *ᵥ
            prodSel (fun k x => if x ∈ R then 0 else w k.succ x) (fun k => !τ k) := by
        rw [prodSel_succ]; simp
      rw [Matrix.mulVec_smul, crOp_ftens, hXt, hYt, hXf, hYf, parMul_prodSel, ftens_smul_left,
        smul_add, smul_smul]
      simp only [Fin.cons_zero, Fin.cons_succ, ite_true, Bool.false_eq_true, ite_false]

/-! ### Supports, splitting and the reordering sign -/

/-- A product of creation operators of orbitals vanishing on `T` has no particle in `T`. -/
lemma isNumR_prodSel_zero (T : Finset (Fin N)) {M : ℕ} (u : Fin M → Fin N → ℂ)
    (hu : ∀ k, ∀ x ∈ T, u k x = 0) (σ : Fin M → Bool) : IsNumR T 0 (prodSel u σ) := by
  induction M with
  | zero =>
    rw [prodSel_zero]
    intro n hn
    simp only [fockVec, Pi.single_apply]
    rw [ite_eq_right]
    rintro rfl
    apply hn
    simp [count, fockCfg, occ]
  | succ M ih =>
    rw [prodSel_succ]
    have ih' := ih (fun k => u k.succ) (fun k => hu k.succ) (fun k => σ k.succ)
    split_ifs
    · exact isNumR_crOp_out T _ (hu 0) _ ih'
    · rw [Matrix.one_mulVec]; exact ih'

/-- `N_T(n) = 0` iff `n` vanishes on `T`. -/
lemma count_eq_zero_iff (T : Finset (Fin N)) (n : Cfg N) : count T n = 0 ↔ ∀ x ∈ T, n x = 0 := by
  unfold count
  rw [Finset.sum_eq_zero_iff_of_nonneg (fun x _ => by unfold occ; split_ifs <;> norm_num)]
  refine forall₂_congr fun x _ => ?_
  unfold occ
  split_ifs with h <;> simp [h]

/-- Entries of `α ⊔ β` in `R`. -/
lemma splitR_symm_mem (R : Finset (Fin N)) (α : R → ZMod 2) (β : {x : Fin N // x ∉ R} → ZMod 2)
    (x : Fin N) (hx : x ∈ R) : (splitR R).symm (α, β) x = α ⟨x, hx⟩ := by
  simp [splitR, hx]

/-- Entries of `α ⊔ β` outside `R`. -/
lemma splitR_symm_not_mem (R : Finset (Fin N)) (α : R → ZMod 2)
    (β : {x : Fin N // x ∉ R} → ZMod 2) (x : Fin N) (hx : x ∉ R) :
    (splitR R).symm (α, β) x = β ⟨x, hx⟩ := by
  simp [splitR, hx]

/-- `(α ⊔ β)|_R = α ⊔ 0`. -/
lemma cutIn_splitR_symm (R : Finset (Fin N)) (α : R → ZMod 2)
    (β : {x : Fin N // x ∉ R} → ZMod 2) :
    cutIn R ((splitR R).symm (α, β)) = (splitR R).symm (α, 0) := by
  funext x
  by_cases hx : x ∈ R
  · simp only [cutIn, hx, ite_true, splitR_symm_mem R _ _ x hx]
  · simp only [cutIn, hx, ite_false, splitR_symm_not_mem R _ _ x hx, Pi.zero_apply]

/-- `(α ⊔ β)|_{R^c} = 0 ⊔ β`. -/
lemma cutOut_splitR_symm (R : Finset (Fin N)) (α : R → ZMod 2)
    (β : {x : Fin N // x ∉ R} → ZMod 2) :
    cutOut R ((splitR R).symm (α, β)) = (splitR R).symm (0, β) := by
  funext x
  by_cases hx : x ∈ R
  · simp only [cutOut, hx, ite_true, splitR_symm_mem R _ _ x hx, Pi.zero_apply]
  · simp only [cutOut, hx, ite_false, splitR_symm_not_mem R _ _ x hx]

/-- `(-1)^k = sgn (k mod 2)`. -/
lemma neg_one_pow_eq_sgn (k : ℕ) : (-1 : ℝ) ^ k = sgn (k : ZMod 2) := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ, ih, Nat.cast_succ, FreeDynAux.sgn_add, FreeDynAux.sgn_one]

/-- The reordering sign is `sgn` of the reordering parity. -/
lemma reorderSign_eq (R : Finset (Fin N)) (α : R → ZMod 2) (β : {x : Fin N // x ∉ R} → ZMod 2) :
    reorderSign R α β = sgn (rpar R ((splitR R).symm (α, β))) := by
  rw [reorderSign, neg_one_pow_eq_sgn, Finset.card_filter]
  congr 1
  push_cast
  rw [Fintype.sum_prod_type, rpar, ← Finset.sum_coe_sort R]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.sum_subtype (univ \ R) (p := fun y => y ∉ R) (fun y => by simp)]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [splitR_symm_mem R α β a a.2, splitR_symm_not_mem R α β b b.2]
  have ha : (⟨(a : Fin N), a.2⟩ : R) = a := rfl
  have hb : (⟨(b : Fin N), b.2⟩ : {x : Fin N // x ∉ R}) = b := rfl
  rw [ha, hb]
  by_cases hlt : (b : Fin N) < (a : Fin N)
  · rcases FreeDynAux.zmod2_cases (α a) with h1 | h1 <;>
      rcases FreeDynAux.zmod2_cases (β b) with h2 | h2 <;> simp [hlt, h1, h2]
  · simp [hlt]

/-- `ε(α, β)² = 1`. -/
lemma reorderSign_mul_self (R : Finset (Fin N)) (α : R → ZMod 2)
    (β : {x : Fin N // x ∉ R} → ZMod 2) : reorderSign R α β * reorderSign R α β = 1 := by
  rw [reorderSign_eq, FreeDynAux.sgn_mul_self]

/-- `(X ⋆ Y)(α ⊔ β) = ε(α, β) X(α ⊔ 0) Y(0 ⊔ β)`. -/
lemma ftens_splitR_symm (R : Finset (Fin N)) (X Y : Cfg N → ℂ) (α : R → ZMod 2)
    (β : {x : Fin N // x ∉ R} → ZMod 2) :
    ftens R X Y ((splitR R).symm (α, β)) =
      (reorderSign R α β : ℂ) * X ((splitR R).symm (α, 0)) * Y ((splitR R).symm (0, β)) := by
  rw [ftens, cutIn_splitR_symm, cutOut_splitR_symm, reorderSign_eq]

/-- Inner products of vectors without particles in `R` are sums over the configurations of the
complement. -/
lemma dot_eq_sum_out (R : Finset (Fin N)) (Y Y' : Cfg N → ℂ) (hY : IsNumR R 0 Y) :
    star Y ⬝ᵥ Y' = ∑ β, star (Y ((splitR R).symm (0, β))) * Y' ((splitR R).symm (0, β)) := by
  rw [dotProduct, ← (splitR R).symm.sum_comp, Fintype.sum_prod_type,
    Finset.sum_eq_single (0 : R → ZMod 2)]
  · rfl
  · intro α _ hα
    refine Finset.sum_eq_zero fun β _ => ?_
    rw [Pi.star_apply, hY _ ?_, star_zero, zero_mul]
    intro hc
    apply hα
    rw [Nat.cast_zero, count_eq_zero_iff] at hc
    funext a
    have := hc a a.2
    rwa [splitR_symm_mem R α β a a.2] at this
  · simp

/-- Inner products of vectors without particles outside `R` are sums over the configurations
of `R`. -/
lemma dot_eq_sum_in (R : Finset (Fin N)) (X X' : Cfg N → ℂ) (hX : IsNumR (univ \ R) 0 X) :
    star X ⬝ᵥ X' = ∑ α, star (X ((splitR R).symm (α, 0))) * X' ((splitR R).symm (α, 0)) := by
  rw [dotProduct, ← (splitR R).symm.sum_comp, Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.sum_eq_single (0 : {x : Fin N // x ∉ R} → ZMod 2)]
  · rfl
  · intro β _ hβ
    rw [Pi.star_apply, hX _ ?_, star_zero, zero_mul]
    intro hc
    apply hβ
    rw [Nat.cast_zero, count_eq_zero_iff] at hc
    funext b
    have := hc b (by simp [b.2])
    rwa [splitR_symm_not_mem R α β b b.2] at this
  · simp

/-! ### The reduced density matrix of a Slater determinant -/

/-- The reduced density matrix of `ψ(α ⊔ β) = ε(α, β) ∑_i c_i X_i(α) Y_i(β)` with orthogonal `Y_i`
and `|c_i| = 1` is `∑_i ‖Y_i‖² |X_i⟩⟨X_i|` (the reordering signs cancel). -/
lemma rdmF_eq_sum (R : Finset (Fin N)) {ι : Type*} [Fintype ι] [DecidableEq ι] (ψ : Cfg N → ℂ)
    (c : ι → ℂ) (X : ι → (R → ZMod 2) → ℂ) (Y : ι → ({x : Fin N // x ∉ R} → ZMod 2) → ℂ)
    (μ : ι → ℂ)
    (hψ : ∀ α β, ψ ((splitR R).symm (α, β)) =
      (reorderSign R α β : ℂ) * ∑ i, c i * (X i α * Y i β))
    (hc : ∀ i, ‖c i‖ = 1) (hY : ∀ i j, ∑ β, star (Y j β) * Y i β = if i = j then μ i else 0)
    (α α' : R → ZMod 2) :
    rdmF R ψ α α' = ∑ i, μ i * (X i α * star (X i α')) := by
  have hterm : ∀ β, ((reorderSign R α β * reorderSign R α' β : ℝ) : ℂ) *
      ψ ((splitR R).symm (α, β)) * star (ψ ((splitR R).symm (α', β))) =
      (∑ i, c i * (X i α * Y i β)) * star (∑ j, c j * (X j α' * Y j β)) := by
    intro β
    have e1 : (reorderSign R α β : ℂ) * (reorderSign R α β : ℂ) = 1 := by
      exact_mod_cast reorderSign_mul_self R α β
    have e2 : (reorderSign R α' β : ℂ) * (reorderSign R α' β : ℂ) = 1 := by
      exact_mod_cast reorderSign_mul_self R α' β
    have e3 : star ((reorderSign R α' β : ℝ) : ℂ) = (reorderSign R α' β : ℂ) := by
      rw [Complex.star_def, Complex.conj_ofReal]
    rw [hψ, hψ, star_mul', e3, Complex.ofReal_mul]
    linear_combination ((reorderSign R α' β : ℂ) * (reorderSign R α' β : ℂ) *
      (∑ i, c i * (X i α * Y i β)) * star (∑ j, c j * (X j α' * Y j β))) * e1 +
      ((∑ i, c i * (X i α * Y i β)) * star (∑ j, c j * (X j α' * Y j β))) * e2
  simp only [rdmF]
  rw [Finset.sum_congr rfl fun β _ => hterm β]
  calc ∑ β, (∑ i, c i * (X i α * Y i β)) * star (∑ j, c j * (X j α' * Y j β))
      = ∑ β, ∑ i, ∑ j, c i * (X i α * Y i β) * star (c j * (X j α' * Y j β)) := by
        simp only [star_sum, Finset.sum_mul_sum]
    _ = ∑ i, ∑ j, ∑ β, c i * (X i α * Y i β) * star (c j * (X j α' * Y j β)) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ i, ∑ j, (c i * star (c j) * X i α * star (X j α')) * (∑ β, star (Y j β) * Y i β) := by
        refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
        simp only [Finset.mul_sum, star_mul']
        exact Finset.sum_congr rfl fun β _ => by ring
    _ = ∑ i, μ i * (X i α * star (X i α')) := by
        simp only [hY, mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
        refine Finset.sum_congr rfl fun i _ => ?_
        have h1 : c i * star (c i) = 1 := by
          rw [Complex.star_def, Complex.mul_conj', hc i]; simp
        linear_combination (μ i * X i α * star (X i α')) * h1

/-- **Schmidt decomposition of a Slater determinant for an arbitrary region** (paper: proof of
Lemma `lem:slater`(a)), summarized: there are canonical occupations `ν_k ∈ [0, 1]` with
`C_R = B Bᴴ`, `Bᴴ B = diag(ν)`, `∑_k √(ν_k(1-ν_k)) ≤ ∑_{x ∈ R, y ∉ R} |K_{xy}|` and
`S(ρ_R) = ∑_k h(ν_k)`. -/
lemma slater_core_bulk (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    ∃ (M : ℕ) (B : Matrix R (Fin M) ℂ) (ν : Fin M → ℝ),
      corrR A S t R = B * Bᴴ ∧ Bᴴ * B = diagonal (fun k => (ν k : ℂ)) ∧
      (∀ k, 0 ≤ ν k ∧ ν k ≤ 1) ∧
      ∑ k, √(ν k * (1 - ν k)) ≤ ∑ x ∈ R, ∑ y ∈ univ \ R, ‖corrMat A S t x y‖ ∧
      vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) =
        ∑ k, Real.binEntropy (ν k) := by
  obtain ⟨M, Φ, c, ν, hc, hψ, hC, h1, h2⟩ := slater_setup_gen A hA S t R
  -- the orbitals restricted to `R` and to its complement
  set a : Fin M → Fin N → ℂ := fun k x => if x ∈ R then Φ x k else 0 with ha_def
  set b : Fin M → Fin N → ℂ := fun k x => if x ∈ R then 0 else Φ x k with hb_def
  have daa : ∀ k l, star (a k) ⬝ᵥ a l = if k = l then (ν k : ℂ) else 0 := by
    intro k l
    rw [← h1 k l]
    have e : ∀ x, star (a k) x * a l x = if x ∈ R then star (Φ x k) * Φ x l else 0 := by
      intro x; by_cases hx : x ∈ R <;> simp [ha_def, hx]
    simp only [dotProduct, e, Finset.sum_ite_mem, Finset.univ_inter]
  have dbb : ∀ k l, star (b k) ⬝ᵥ b l = if k = l then ((1 - ν k : ℝ) : ℂ) else 0 := by
    intro k l
    rw [← h2 k l]
    have e : ∀ x, star (b k) x * b l x =
        if x ∈ univ \ R then star (Φ x k) * Φ x l else 0 := by
      intro x; by_cases hx : x ∈ R <;> simp [hb_def, hx]
    simp only [dotProduct, e, Finset.sum_ite_mem, Finset.univ_inter]
  have hν0 : ∀ k, 0 ≤ ν k := by
    intro k
    have h := daa k k
    rw [ite_eq_left rfl, star_dotProduct_self] at h
    have h' : nsq (a k) = ν k := by exact_mod_cast h
    rw [← h']; exact nsq_nonneg _
  have hν1 : ∀ k, ν k ≤ 1 := by
    intro k
    have h := dbb k k
    rw [ite_eq_left rfl, star_dotProduct_self] at h
    have h' : nsq (b k) = 1 - ν k := by exact_mod_cast h
    linarith [nsq_nonneg (b k)]
  -- `C_R = B Bᴴ` and `Bᴴ B = diag(ν)`
  set B : Matrix R (Fin M) ℂ := Matrix.of fun x k => Φ (x : Fin N) k with hB_def
  have hCB : corrR A S t R = B * Bᴴ := by
    ext x y
    simp [corrR, hC, Matrix.mul_apply, hB_def, Matrix.conjTranspose_apply]
  have hBB : Bᴴ * B = diagonal (fun k => (ν k : ℂ)) := by
    ext k l
    rw [Matrix.mul_apply, Matrix.diagonal_apply, ← h1 k l]
    simp only [hB_def, Matrix.conjTranspose_apply, Matrix.of_apply]
    exact Finset.sum_coe_sort R (fun x => star (Φ x k) * Φ x l)
  have hsqrt := sum_sqrt_nu_le Φ R ν h1 h2
  rw [← hC] at hsqrt
  -- the expansion in the fermionic product basis
  obtain ⟨s, hs, hexp⟩ := slater_expand_gen R (fun k x => Φ x k)
  set X : (Fin M → Bool) → Cfg N → ℂ := fun σ => prodSel a σ with hX_def
  set Y : (Fin M → Bool) → Cfg N → ℂ := fun σ => prodSel b (fun k => !σ k) with hY_def
  set ψ := evo (dGamma A) t *ᵥ fockVec S with hψ_def
  have hψ' : ψ = ∑ σ, (c * s σ) • ftens R (X σ) (Y σ) := by
    rw [hψ, hexp, Finset.smul_sum]
    simp only [smul_smul]
    rfl
  have hcs : ∀ σ, ‖c * s σ‖ = 1 := fun σ => by rw [norm_mul, hc, hs, one_mul]
  have hX : ∀ σ τ, star (X σ) ⬝ᵥ X τ =
      if σ = τ then ∏ k, (if σ k then (ν k : ℂ) else 1) else 0 :=
    fun σ τ => gram_prodSel a _ daa σ τ
  set μ : (Fin M → Bool) → ℝ := fun σ => ∏ k, (if σ k then 1 else 1 - ν k) with hμ_def
  have hYg : ∀ σ τ, star (Y τ) ⬝ᵥ Y σ = if σ = τ then (μ σ : ℂ) else 0 := by
    intro σ τ
    have hg := gram_prodSel b _ dbb (fun k => !τ k) (fun k => !σ k)
    by_cases h : σ = τ
    · subst h
      rw [hY_def]
      simp only [hg, ite_true, hμ_def]
      push_cast
      refine Finset.prod_congr rfl fun k _ => ?_
      cases σ k <;> simp
    · have h' : (fun k => !τ k) ≠ (fun k => !σ k) := by
        intro h'; apply h; funext k; have := congrFun h' k; simpa using this.symm
      rw [hY_def]
      simp only [hg, h, h', ite_false]
  -- supports of the two factors
  have hYn : ∀ σ, IsNumR R 0 (Y σ) := fun σ =>
    isNumR_prodSel_zero R b (fun k x hx => by simp [hb_def, hx]) _
  have hXn : ∀ σ, IsNumR (univ \ R) 0 (X σ) := fun σ =>
    isNumR_prodSel_zero (univ \ R) a
      (fun k x hx => by simp [ha_def, (Finset.mem_sdiff.mp hx).2]) _
  -- the reduced density matrix
  have hv : ∀ α β, ψ ((splitR R).symm (α, β)) = (reorderSign R α β : ℂ) *
      ∑ σ, (c * s σ) * (X σ ((splitR R).symm (α, 0)) * Y σ ((splitR R).symm (0, β))) := by
    intro α β
    rw [hψ', Finset.sum_apply, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    rw [Pi.smul_apply, smul_eq_mul, ftens_splitR_symm]
    ring
  have hYh : ∀ σ τ, ∑ β, star (Y τ ((splitR R).symm (0, β))) * Y σ ((splitR R).symm (0, β)) =
      if σ = τ then (μ σ : ℂ) else 0 := by
    intro σ τ; rw [← dot_eq_sum_out R _ _ (hYn τ), hYg]
  have hρ : ∀ α α', rdmF R ψ α α' = ∑ σ, (μ σ : ℂ) * (X σ ((splitR R).symm (α, 0)) *
      star (X σ ((splitR R).symm (α', 0)))) :=
    rdmF_eq_sum R ψ (fun σ => c * s σ) (fun σ α => X σ ((splitR R).symm (α, 0)))
      (fun σ β => Y σ ((splitR R).symm (0, β))) (fun σ => (μ σ : ℂ)) hv hcs hYh
  have hμ0 : ∀ σ, 0 ≤ μ σ := fun σ =>
    Finset.prod_nonneg fun k _ => by split_ifs <;> linarith [hν1 k]
  have hXnsq : ∀ σ, nsq (X σ) = ∏ k, (if σ k then ν k else 1) := by
    intro σ
    have h := hX σ σ
    rw [ite_eq_left rfl, star_dotProduct_self] at h
    have h' : ((∏ k, (if σ k then ν k else 1) : ℝ) : ℂ) =
        ∏ k, (if σ k then (ν k : ℂ) else 1) := by
      rw [Complex.ofReal_prod]
      exact Finset.prod_congr rfl fun k _ => by split_ifs <;> simp
    rw [← h'] at h
    exact_mod_cast h
  set lam : (Fin M → Bool) → ℝ := fun σ => ∏ k, (if σ k then ν k else 1 - ν k) with hlam_def
  have hlam : ∀ σ, μ σ * nsq (X σ) = lam σ := by
    intro σ
    rw [hXnsq, hμ_def, hlam_def]
    simp only
    rw [← Finset.prod_mul_distrib]
    refine Finset.prod_congr rfl fun k _ => ?_
    cases σ k <;> simp
  -- the entropy
  set Z : Matrix (R → ZMod 2) (Fin M → Bool) ℂ :=
    Matrix.of fun α σ => ((√(μ σ) : ℝ) : ℂ) * X σ ((splitR R).symm (α, 0)) with hZ_def
  have hsq : ∀ σ, ((√(μ σ) : ℝ) : ℂ) * ((√(μ σ) : ℝ) : ℂ) = (μ σ : ℂ) := fun σ => by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hμ0 σ)]
  have hZZ : rdmF R ψ = Z * Zᴴ := by
    ext α α'
    rw [hρ, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp only [hZ_def, Matrix.of_apply, Matrix.conjTranspose_apply, star_mul', Complex.star_def,
      Complex.conj_ofReal]
    rw [← hsq σ]
    ring
  have hZZ' : Zᴴ * Z = diagonal (fun σ => (lam σ : ℂ)) := by
    ext σ τ
    rw [Matrix.mul_apply, Matrix.diagonal_apply]
    have e : ∑ α, Zᴴ σ α * Z α τ =
        ((√(μ σ) : ℝ) : ℂ) * ((√(μ τ) : ℝ) : ℂ) * (star (X σ) ⬝ᵥ X τ) := by
      rw [dot_eq_sum_in R _ _ (hXn σ), Finset.mul_sum]
      refine Finset.sum_congr rfl fun α _ => ?_
      simp only [hZ_def, Matrix.of_apply, Matrix.conjTranspose_apply, star_mul',
        Complex.star_def, Complex.conj_ofReal]
      ring
    rw [e, hX]
    split_ifs with h
    · subst h
      rw [hsq, ← hlam σ, hXnsq, Complex.ofReal_mul,
        Complex.ofReal_prod (s := univ) (f := fun k => if σ k then ν k else 1)]
      congr 1
      exact Finset.prod_congr rfl fun k _ => by split_ifs <;> simp
    · simp
  have hHd : (diagonal (fun σ => (lam σ : ℂ))).IsHermitian :=
    hZZ' ▸ isHermitian_conjTranspose_mul_self Z
  have hent : vonNeumann (rdmF R ψ) (rdmF_isHermitian R ψ) = ∑ k, Real.binEntropy (ν k) := by
    rw [CutLocalAux.vonNeumann_congr _ (isHermitian_mul_conjTranspose_self Z) hZZ,
      vonNeumann_mul_conjTranspose_comm, CutLocalAux.vonNeumann_congr _ hHd hZZ',
      vonNeumann_diagonal, ← sum_negMulLog_prod_sel ν]
  exact ⟨M, B, ν, hCB, hBB, fun k => ⟨hν0 k, hν1 k⟩, hsqrt, hent⟩

/-! ### Left segments -/

/-- The von Neumann entropy is invariant under reindexing. -/
lemma vonNeumann_submatrix_equiv' {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] (M : Matrix ι ι ℂ) (hM : M.IsHermitian) (e : κ ≃ ι)
    (h' : (M.submatrix e e).IsHermitian) : vonNeumann (M.submatrix e e) h' = vonNeumann M hM := by
  refine CutLocalAux.vonNeumann_eq_of_charpoly _ h' M hM 0 0 ?_
  have h := Matrix.charpoly_reindex e.symm M
  rw [Matrix.reindex_apply, Equiv.symm_symm] at h
  rw [pow_zero, one_mul, one_mul, h]

/-- The sites of the left segment, indexed by `Fin ℓ`. -/
def segEquiv (ℓ : ℕ) (hℓ : ℓ ≤ N) : Fin ℓ ≃ (segment (N := N) ℓ) where
  toFun a := ⟨inL ℓ hℓ a, by simp [segment, inL]⟩
  invFun x := ⟨x.1.val, by
    have := x.2
    simp only [segment, Finset.mem_filter] at this
    exact this.2⟩
  left_inv a := rfl
  right_inv x := rfl

/-- The sites of the complement of the left segment, indexed by `Fin (N - ℓ)`. -/
def segCompEquiv (ℓ : ℕ) (hℓ : ℓ ≤ N) : Fin (N - ℓ) ≃ {x : Fin N // x ∉ segment ℓ} where
  toFun b := ⟨inR ℓ hℓ b, by simp [segment, inR]⟩
  invFun x := ⟨x.1.val - ℓ, by
    have h1 := x.2
    have h2 := x.1.isLt
    simp only [segment, Finset.mem_filter, Finset.mem_univ, true_and, not_lt] at h1
    omega⟩
  left_inv b := by ext; simp [inR]
  right_inv x := by
    have h1 := x.2
    simp only [segment, Finset.mem_filter, Finset.mem_univ, true_and, not_lt] at h1
    ext; simp [inR]; omega

/-- Recombining along `splitR (segment ℓ)` is recombining along `splitCfg ℓ`. -/
lemma splitR_segment_symm (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) (b : Cfg (N - ℓ)) :
    (splitR (segment ℓ)).symm (Equiv.arrowCongr (segEquiv ℓ hℓ) (Equiv.refl _) a,
      Equiv.arrowCongr (segCompEquiv ℓ hℓ) (Equiv.refl _) b) = (splitCfg ℓ hℓ).symm (a, b) := by
  funext x
  by_cases hx : x.val < ℓ
  · have hxs : x ∈ segment (N := N) ℓ := by simp [segment, hx]
    rw [splitR_symm_mem _ _ _ x hxs, CutLocalAux.splitCfg_symm_apply_lt ℓ hℓ a b x hx]
    rfl
  · have hxs : x ∉ segment (N := N) ℓ := by simp [segment, hx]
    rw [splitR_symm_not_mem _ _ _ x hxs,
      CutLocalAux.splitCfg_symm_apply_ge ℓ hℓ a b x (not_lt.mp hx)]
    rfl

/-- For a left segment there are no reordering signs. -/
lemma reorderSign_segment (ℓ : ℕ) (α : segment (N := N) ℓ → ZMod 2)
    (β : {x : Fin N // x ∉ segment ℓ} → ZMod 2) : reorderSign (segment ℓ) α β = 1 := by
  rw [reorderSign, Finset.filter_false_of_mem, Finset.card_empty, pow_zero]
  rintro ⟨x, y⟩ _ ⟨hyx, -⟩
  have h1 := x.2
  have h2 := y.2
  simp only [segment, Finset.mem_filter, Finset.mem_univ, true_and, not_lt] at h1 h2
  dsimp only at hyx
  have := Fin.lt_def.mp hyx
  omega

/-- For a left segment, `ρ_R` is `rdm` up to reindexing. -/
lemma rdmF_segment_submatrix (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) :
    (rdmF (segment ℓ) v).submatrix (Equiv.arrowCongr (segEquiv ℓ hℓ) (Equiv.refl _))
      (Equiv.arrowCongr (segEquiv ℓ hℓ) (Equiv.refl _)) = rdm ℓ hℓ v := by
  ext a a'
  simp only [Matrix.submatrix_apply, rdmF, rdm, reorderSign_segment, mul_one, Complex.ofReal_one,
    one_mul]
  rw [← (Equiv.arrowCongr (segCompEquiv ℓ hℓ) (Equiv.refl (ZMod 2))).sum_comp]
  simp only [splitR_segment_symm]

end SlaterBulkAux

/-- **Entanglement entropy of a Slater determinant for any region** (paper: Lemma `lem:slater`(a)):
`S(ρ_R) = ∑_i h(ν_i)`. -/
theorem vonNeumann_rdmF_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) =
      ∑ i, Real.binEntropy ((corrR_isHermitian A S t R).eigenvalues i) := by
  obtain ⟨M, B, ν, hCB, hBB, -, -, hent⟩ := SlaterBulkAux.slater_core_bulk A hA S t R
  have hm := SlaterExactAux.mset_eigenvalues _ (corrR_isHermitian A S t R) B hCB ν hBB
  rw [hent]
  exact (SlaterExactAux.sum_eq_of_mset hm Real.binEntropy Real.binEntropy_zero).symm

/-- **Area law for Slater determinants and any region** (paper: Lemma `lem:slater`(b), Lemma `lem:binent`):
`S(ρ_R) ≤ 2 ln 2 ∑_{x ∈ R, y ∉ R} |K_{xy}|`. -/
theorem vonNeumann_rdmF_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) ≤
      2 * Real.log 2 * ∑ x ∈ R, ∑ y ∈ univ \ R, ‖corrMat A S t x y‖ := by
  obtain ⟨M, B, ν, hCB, hBB, -, hsq, -⟩ := SlaterBulkAux.slater_core_bulk A hA S t R
  have hm := SlaterExactAux.mset_eigenvalues _ (corrR_isHermitian A S t R) B hCB ν hBB
  rw [vonNeumann_rdmF_slater A hA S t R]
  calc ∑ i, Real.binEntropy ((corrR_isHermitian A S t R).eigenvalues i)
      ≤ ∑ i, 2 * Real.log 2 * √((corrR_isHermitian A S t R).eigenvalues i *
          (1 - (corrR_isHermitian A S t R).eigenvalues i)) :=
        Finset.sum_le_sum fun i _ => binEntropy_le_two_log_two_sqrt _
          (eig_corrR_mem A hA S t R i).1 (eig_corrR_mem A hA S t R i).2
    _ = 2 * Real.log 2 * ∑ i, √((corrR_isHermitian A S t R).eigenvalues i *
          (1 - (corrR_isHermitian A S t R).eigenvalues i)) := by
        rw [Finset.mul_sum]
    _ = 2 * Real.log 2 * ∑ k, √(ν k * (1 - ν k)) := by
        rw [SlaterExactAux.sum_eq_of_mset hm (fun x => √(x * (1 - x))) (by simp)]
    _ ≤ 2 * Real.log 2 * ∑ x ∈ R, ∑ y ∈ univ \ R, ‖corrMat A S t x y‖ :=
        mul_le_mul_of_nonneg_left hsq AreaLawSharpAux.two_log_two_nonneg

/-- For a left segment, `ρ_R` is the reduced state `rdm` of `DFL.AreaLaw` (no reordering signs). -/
theorem vonNeumann_rdmF_segment (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) :
    vonNeumann (rdmF (segment ℓ) v) (rdmF_isHermitian _ v) =
      vonNeumann (rdm ℓ hℓ v) (rdm_isHermitian ℓ hℓ v) := by
  have hsub := SlaterBulkAux.rdmF_segment_submatrix ℓ hℓ v
  have hH : ((rdmF (segment ℓ) v).submatrix
      (Equiv.arrowCongr (SlaterBulkAux.segEquiv ℓ hℓ) (Equiv.refl _))
      (Equiv.arrowCongr (SlaterBulkAux.segEquiv ℓ hℓ) (Equiv.refl _))).IsHermitian :=
    hsub ▸ rdm_isHermitian ℓ hℓ v
  rw [← CutLocalAux.vonNeumann_congr hH (rdm_isHermitian ℓ hℓ v) hsub]
  exact (SlaterBulkAux.vonNeumann_submatrix_equiv' _ (rdmF_isHermitian _ v) _ hH).symm

/-- **Area law for bulk intervals** (paper: finite-volume form of Theorem `thm:infvol`(d)):
`E_q sup_t S(ρ^q_B(t)) ≤ 2 ln 2 · 2Θ(F) = 4 ln 2 Θ(F)` for every interval `B = [a, b]`. -/
theorem lgt_areaLaw_bulk (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      vonNeumann (rdmF (ivl a b) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdmF_isHermitian _ _)) ≤
      4 * Real.log 2 * thetaF F := by
  -- pointwise in `q`, a time-independent majorant: the cut sums at `a` and at `b + 1`
  have hpt : ∀ q : Cfg N, (⨆ t : ℝ,
      vonNeumann (rdmF (ivl a b) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
        (rdmF_isHermitian _ _)) ≤
      2 * Real.log 2 * ((∑ x ∈ segment (N := N) a, ∑ y ∈ univ \ segment a,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z) +
      ∑ x ∈ segment (N := N) (b + 1), ∑ y ∈ univ \ segment (b + 1),
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z) := by
    intro q
    refine ciSup_le fun t =>
      (vonNeumann_rdmF_le (hq J h q) (hq_isSymm J h q) S t (ivl a b)).trans ?_
    refine mul_le_mul_of_nonneg_left ?_ AreaLawSharpAux.two_log_two_nonneg
    refine (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      AreaLawSharpAux.norm_corrMat_le_corrQ J h q S t x y).trans ?_
    exact AnnealedIntervalAux.sum_ivl_le_two_cuts
      (fun x y => ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z)
      (fun x y => Finset.sum_nonneg fun z _ =>
        mul_nonneg (corrQ_nonneg _ _ x z) (corrQ_nonneg _ _ y z))
      (fun x y => Finset.sum_congr rfl fun z _ => mul_comm _ _) a b
  refine (AreaLawAux.avgQ_mono' hpt).trans ?_
  rw [AreaLawAux.avgQ_const_mul', AnnealedIntervalAux.avgQ_add_eq]
  have h1 := AreaLawSharpAux.avgQ_cut_corrQ_le (N := N) hEC a
  have h2 := AreaLawSharpAux.avgQ_cut_corrQ_le (N := N) hEC (b + 1)
  calc 2 * Real.log 2 * (_ + _) ≤ 2 * Real.log 2 * (thetaF F + thetaF F) :=
        mul_le_mul_of_nonneg_left (add_le_add h1 h2) AreaLawSharpAux.two_log_two_nonneg
    _ = 4 * Real.log 2 * thetaF F := by ring

end DFL
