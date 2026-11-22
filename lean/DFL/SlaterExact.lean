import DFL.AreaLaw
import DFL.VNEntropy

/-!
# Reduced states of Slater determinants (paper: Lemmas `lem:quasifree`, `lem:slater`)

For the evolved Fock state `ψ = Γ(e^{-itA}) |S⟩` (a Slater determinant with one-body projection
`K = corrMat A S t`) and the left segment `{0, …, ℓ-1}`, let `C_A = P_A K P_A` (`corrA`), with
eigenvalues `ν_1, …, ν_ℓ ∈ [0, 1]`. The canonical orbitals give a Schmidt decomposition of `ψ` across
the cut, with Schmidt coefficients `∏_{j ∈ T} ν_j ∏_{j ∉ T} (1 - ν_j)`, whence

* `vonNeumann_rdm_slater`: `S(ρ_A) = ∑_i h(ν_i)`, `h` the binary entropy (paper: `lem:slater`(a));
* `numDist_slater`: the law of `N_A` is the Poisson-binomial law of the `ν_i`;
* `genfun_slater`: its generating function is `det(1 + (z-1) C_A)` (paper: `lem:quasifree`(b));
* `sum_sqrt_eig_le`: `∑_i √(ν_i (1 - ν_i)) ≤ ∑_{x ∈ A, y ∉ A} |K_{xy}|` (paper: `lem:slater`(b)).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `C_A = P_A K P_A`: the correlation matrix restricted to the left segment `{0, …, ℓ-1}`. -/
def corrA (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    Matrix (Fin ℓ) (Fin ℓ) ℂ :=
  (corrMat A S t).submatrix (AreaLawAux.inL ℓ hℓ) (AreaLawAux.inL ℓ hℓ)

/-- The Poisson-binomial law of the parameters `ν`. -/
def poissonBinomial {ι : Type*} [Fintype ι] [DecidableEq ι] (ν : ι → ℝ) (k : ℕ) : ℝ :=
  ∑ T ∈ (univ : Finset (Finset ι)).filter (fun T => T.card = k),
    (∏ i ∈ T, ν i) * ∏ i ∈ univ \ T, (1 - ν i)

namespace SlaterExactAux

open AreaLawAux

variable {K : ℕ}

/-! ### Creation and annihilation operators of vectors -/

/-- General CAR: `c(v) c†(u) + c†(u) c(v) = ⟨v, u⟩`. -/
lemma crOp_CAR_gen (v u : Fin K → ℂ) :
    (crOp v)ᴴ * crOp u + crOp u * (crOp v)ᴴ =
      (star v ⬝ᵥ u) • (1 : Matrix (Cfg K) (Cfg K) ℂ) := by
  have hcar : ∀ x y : Fin K, cplx (ann x) * cplx (cre y) + cplx (cre y) * cplx (ann x) =
      if x = y then 1 else 0 := by
    intro x y
    rw [← cplx_mul', ← cplx_mul']
    have h := congrArg cplx (ann_mul_cre_add x y)
    have hadd : cplx (ann x * cre y + cre y * ann x) =
        cplx (ann x * cre y) + cplx (cre y * ann x) := by simp [cplx, Matrix.map_add]
    rw [← hadd, h]
    split_ifs
    · exact cplx_one'
    · simp [cplx]
  rw [crOp_conjTranspose]
  calc (∑ x, star (v x) • cplx (ann x)) * crOp u + crOp u * (∑ x, star (v x) • cplx (ann x))
      = ∑ x, ∑ y, ((star (v x) * u y) • (cplx (ann x) * cplx (cre y)) +
          (u y * star (v x)) • (cplx (cre y) * cplx (ann x))) := by
        simp only [crOp, Finset.sum_mul_sum, smul_mul_smul_comm, Finset.sum_add_distrib]
        rw [Finset.sum_comm (f := fun x y => (u x * star (v y)) • (cplx (cre x) * cplx (ann y)))]
    _ = ∑ x, ∑ y, (star (v x) * u y) • (if x = y then (1 : Matrix (Cfg K) (Cfg K) ℂ) else 0) := by
        refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
        rw [mul_comm (u y), ← smul_add, hcar]
    _ = (star v ⬝ᵥ u) • 1 := by
        simp_rw [smul_ite, smul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true,
          ← Finset.sum_smul]
        rfl

/-- `c(v) |0⟩ = 0`. -/
lemma crOp_adj_vac (v : Fin K → ℂ) : (crOp v)ᴴ *ᵥ fockVec (∅ : Finset (Fin K)) = 0 := by
  rw [crOp_conjTranspose, Matrix.sum_mulVec]
  simp [Matrix.smul_mulVec, ann_mulVec_vac]

/-- `⟨c†(u) X, Y⟩ = ⟨X, c(u) Y⟩`. -/
lemma dot_crOp_left (u : Fin K → ℂ) (X Y : Cfg K → ℂ) :
    star (crOp u *ᵥ X) ⬝ᵥ Y = star X ⬝ᵥ ((crOp u)ᴴ *ᵥ Y) := by
  rw [Matrix.star_mulVec, Matrix.dotProduct_mulVec]

/-- `⟨X, c†(u) Y⟩ = ⟨c(u) X, Y⟩`. -/
lemma dot_crOp_right (u : Fin K → ℂ) (X Y : Cfg K → ℂ) :
    star X ⬝ᵥ (crOp u *ᵥ Y) = star ((crOp u)ᴴ *ᵥ X) ⬝ᵥ Y := by
  rw [Matrix.star_mulVec, Matrix.conjTranspose_conjTranspose, Matrix.dotProduct_mulVec]

/-! ### Products of selected creation operators -/

/-- The ordered product `∏_{k : σ k} c†(u_k) |0⟩` of the selected creation operators applied to
the vacuum. -/
def prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) (σ : Fin M → Bool) : Cfg K → ℂ :=
  (List.ofFn fun k => if σ k then crOp (u k) else 1).prod *ᵥ fockVec ∅

/-- `prodSel` with no orbitals is the vacuum. -/
lemma prodSel_zero (u : Fin 0 → Fin K → ℂ) (σ : Fin 0 → Bool) : prodSel u σ = fockVec ∅ := by
  simp [prodSel]

/-- Peeling off the first orbital of `prodSel`. -/
lemma prodSel_succ {M : ℕ} (u : Fin (M + 1) → Fin K → ℂ) (σ : Fin (M + 1) → Bool) :
    prodSel u σ = (if σ 0 then crOp (u 0) else 1) *ᵥ
      prodSel (fun k => u k.succ) (fun k => σ k.succ) := by
  rw [prodSel, prodSel, List.ofFn_succ, List.prod_cons, ← Matrix.mulVec_mulVec]

/-- `c(v)` kills `prodSel u σ` if `v` is orthogonal to the selected orbitals. -/
lemma adj_prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) (σ : Fin M → Bool) (v : Fin K → ℂ)
    (hv : ∀ k, σ k = true → star v ⬝ᵥ u k = 0) : (crOp v)ᴴ *ᵥ prodSel u σ = 0 := by
  induction M with
  | zero => rw [prodSel_zero, crOp_adj_vac]
  | succ M ih =>
    have ih' := ih (fun k => u k.succ) (fun k => σ k.succ) (fun k hk => hv k.succ hk)
    rw [prodSel_succ]
    by_cases h0 : σ 0 = true
    · rw [ite_eq_left h0, Matrix.mulVec_mulVec]
      have hcar := crOp_CAR_gen v (u 0)
      rw [hv 0 h0, zero_smul, add_eq_zero_iff_eq_neg] at hcar
      rw [hcar, Matrix.neg_mulVec, ← Matrix.mulVec_mulVec, ih', Matrix.mulVec_zero, neg_zero]
    · rw [ite_eq_right h0, Matrix.one_mulVec, ih']

/-- **CAR Gram identity**: for an orthogonal family `u`, the vectors `prodSel u σ` are orthogonal,
with `‖prodSel u σ‖² = ∏_{k : σ k} ‖u_k‖²`. -/
lemma gram_prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) (d : Fin M → ℂ)
    (hu : ∀ k l, star (u k) ⬝ᵥ u l = if k = l then d k else 0) (σ τ : Fin M → Bool) :
    star (prodSel u σ) ⬝ᵥ prodSel u τ = if σ = τ then ∏ k, (if σ k then d k else 1) else 0 := by
  induction M with
  | zero =>
    have h : σ = τ := Subsingleton.elim _ _
    simp [prodSel_zero, h, fockVec]
  | succ M ih =>
    have hu' : ∀ k l : Fin M, star (u k.succ) ⬝ᵥ u l.succ = if k = l then d k.succ else 0 := by
      intro k l; rw [hu]; simp [Fin.succ_inj]
    have ih' := ih (fun k => u k.succ) (fun k => d k.succ) hu'
    have hkill : ∀ ρ : Fin M → Bool, (crOp (u 0))ᴴ *ᵥ prodSel (fun k => u k.succ) ρ = 0 :=
      fun ρ => adj_prodSel _ _ _ (fun k _ => by rw [hu]; simp [(Fin.succ_ne_zero k).symm])
    have hσ : σ = τ ↔ σ 0 = τ 0 ∧ (fun k : Fin M => σ k.succ) = (fun k : Fin M => τ k.succ) := by
      constructor
      · rintro rfl; exact ⟨rfl, rfl⟩
      · rintro ⟨h0, hs⟩
        funext k
        exact Fin.cases h0 (fun k => congrFun hs k) k
    rw [prodSel_succ, prodSel_succ, Fin.prod_univ_succ]
    simp only [hσ]
    cases h0 : σ 0 <;> cases h1 : τ 0
    · simp only [Bool.false_eq_true, ite_false, Matrix.one_mulVec, ih', true_and, one_mul]
    · simp only [ite_false, ite_true, Matrix.one_mulVec, dot_crOp_right, hkill, star_zero,
        zero_dotProduct, false_and, Bool.false_eq, Bool.true_eq_false]
    · simp only [Bool.false_eq_true, ite_false, ite_true, Matrix.one_mulVec, dot_crOp_left,
        hkill, dotProduct_zero, false_and, Bool.true_eq_false]
    · simp only [ite_true, true_and]
      rw [dot_crOp_left, Matrix.mulVec_mulVec]
      have hcar := crOp_CAR_gen (u 0) (u 0)
      rw [hu, ite_eq_left rfl, ← eq_sub_iff_add_eq] at hcar
      rw [hcar, Matrix.sub_mulVec, ← Matrix.mulVec_mulVec, hkill, Matrix.mulVec_zero, sub_zero,
        Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, ih', smul_eq_mul]
      split_ifs <;> ring

/-! ### Parity -/

/-- Multiplication by the parity `(-1)^{N}`. -/
def parMul (Z : Cfg K → ℂ) : Cfg K → ℂ := fun n => parC n * Z n

/-- Emptying an occupied mode flips the parity sign. -/
lemma parC_update (n : Cfg K) (x : Fin K) (hx : n x = 1) :
    parC n = -parC (Function.update n x 0) := by
  have h := CutLocalAux.par_eq_update_add_one n x hx
  simp only [CutLocalAux.par] at h
  simp only [parC, h, sgn_add, FreeDynAux.sgn_one]
  push_cast; ring

/-- The parity anticommutes with `c_x†`. -/
lemma parMul_cre (x : Fin K) (Z : Cfg K → ℂ) :
    parMul (cplx (cre x) *ᵥ Z) = -(cplx (cre x) *ᵥ parMul Z) := by
  funext n
  simp only [parMul, Pi.neg_apply, cre_mulVec_apply]
  split_ifs with h
  · rw [parC_update n x h]; ring
  · simp

/-- The parity anticommutes with `c†(u)`. -/
lemma parMul_crOp (u : Fin K → ℂ) (Z : Cfg K → ℂ) :
    parMul (crOp u *ᵥ Z) = -(crOp u *ᵥ parMul Z) := by
  rw [crOp_mulVec_eq, crOp_mulVec_eq]
  funext n
  simp only [parMul, Finset.sum_apply, Pi.smul_apply, Pi.neg_apply, smul_eq_mul, Finset.mul_sum,
    ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun x _ => ?_
  have h := congrFun (parMul_cre x Z) n
  simp only [parMul, Pi.neg_apply] at h
  rw [mul_left_comm, h]; ring

/-- The vacuum is even. -/
lemma parMul_vac : parMul (fockVec (∅ : Finset (Fin K))) = fockVec ∅ := by
  funext n
  by_cases hn : n = fockCfg ∅
  · subst hn
    have h0 : ∑ i, fockCfg (∅ : Finset (Fin K)) i = 0 := by simp [fockCfg]
    simp [parMul, parC, h0, sgn]
  · simp [parMul, fockVec, hn]

/-- `parMul` is linear. -/
lemma parMul_smul (c : ℂ) (Z : Cfg K → ℂ) : parMul (c • Z) = c • parMul Z := by
  funext n; simp only [parMul, Pi.smul_apply, smul_eq_mul]; ring

/-- `prodSel u σ` has the parity of the number of selected orbitals. -/
lemma parMul_prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) (σ : Fin M → Bool) :
    parMul (prodSel u σ) = (∏ k, if σ k then (-1 : ℂ) else 1) • prodSel u σ := by
  induction M with
  | zero => simp [prodSel_zero, parMul_vac]
  | succ M ih =>
    rw [prodSel_succ, Fin.prod_univ_succ]
    cases h : σ 0
    · simp only [Bool.false_eq_true, ite_false, Matrix.one_mulVec, one_mul]
      exact ih _ _
    · simp only [ite_true, parMul_crOp, ih, Matrix.mulVec_smul, neg_one_mul, neg_smul]

/-! ### Particle number -/

/-- `Z` is supported on configurations with exactly `m` particles. -/
def IsNum (m : ℕ) (Z : Cfg K → ℂ) : Prop := ∀ n, count univ n ≠ (m : ℝ) → Z n = 0

/-- Emptying an occupied mode removes one particle. -/
lemma count_update (n : Cfg K) (x : Fin K) (hx : n x = 1) :
    count univ n = count univ (Function.update n x 0) + 1 := by
  unfold count
  rw [← Finset.add_sum_erase _ _ (mem_univ x), ← Finset.add_sum_erase _ _ (mem_univ x),
    Function.update_self, hx]
  have : ∑ k ∈ univ.erase x, occ (Function.update n x 0 k) = ∑ k ∈ univ.erase x, occ (n k) :=
    Finset.sum_congr rfl fun k hk => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hk)]
  rw [this]
  simp [occ]
  ring

/-- `c†(u)` raises the particle number by one. -/
lemma isNum_crOp {m : ℕ} (u : Fin K → ℂ) (Z : Cfg K → ℂ) (hZ : IsNum m Z) :
    IsNum (m + 1) (crOp u *ᵥ Z) := by
  intro n hn
  rw [crOp_mulVec_eq, Finset.sum_apply]
  refine Finset.sum_eq_zero fun x _ => ?_
  rw [Pi.smul_apply, cre_mulVec_apply]
  split_ifs with h
  · rw [hZ _ ?_, mul_zero, smul_zero]
    intro hc; apply hn; rw [count_update n x h, hc]; push_cast; ring
  · simp

/-- `prodSel u σ` has exactly as many particles as selected orbitals. -/
lemma isNum_prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) (σ : Fin M → Bool) :
    IsNum (∑ k, if σ k then 1 else 0) (prodSel u σ) := by
  induction M with
  | zero =>
    intro n hn
    rw [prodSel_zero]
    simp only [fockVec, Pi.single_apply]
    rw [ite_eq_right]
    rintro rfl
    apply hn
    simp [count, fockCfg, occ]
  | succ M ih =>
    rw [prodSel_succ, Fin.sum_univ_succ]
    cases h : σ 0
    · simp only [Bool.false_eq_true, ite_false, Matrix.one_mulVec, zero_add]
      exact ih _ _
    · simp only [ite_true]
      rw [add_comm]
      exact isNum_crOp _ _ (ih _ _)

/-! ### Expansion of a Slater determinant across the cut -/

/-- Splitting a sum over `Fin (M+1) → Bool` according to the first value. -/
lemma sum_fun_succ {β : Type*} [AddCommMonoid β] {M : ℕ} (f : (Fin (M + 1) → Bool) → β) :
    ∑ σ, f σ = ∑ τ : Fin M → Bool, (f (Fin.cons true τ) + f (Fin.cons false τ)) := by
  rw [← (Fin.consEquiv fun _ => Bool).sum_comp, Fintype.sum_prod_type, Fintype.sum_bool,
    Finset.sum_add_distrib]
  rfl

/-- **Expansion of a Slater determinant across the cut**: writing `w_k = a_k ⊕ b_k`,
`c†(w_1) ⋯ c†(w_M) |0⟩ = ∑_σ s_σ (∏_{σ k} c†(a_k) |0⟩) ⊗ (∏_{¬σ k} c†(b_k) |0⟩)` with `|s_σ| = 1`. -/
lemma slater_expand (ℓ : ℕ) (hℓ : ℓ ≤ N) {M : ℕ} (w : Fin M → Fin N → ℂ) :
    ∃ s : (Fin M → Bool) → ℂ, (∀ σ, ‖s σ‖ = 1) ∧
      slaterAlt N M w = ∑ σ, s σ • tens ℓ hℓ (prodSel (fun k a => w k (inL ℓ hℓ a)) σ)
        (prodSel (fun k b => w k (inR ℓ hℓ b)) (fun k => !σ k)) := by
  induction M with
  | zero =>
    refine ⟨fun _ => 1, fun _ => by simp, ?_⟩
    rw [Fintype.sum_unique]
    simp [slaterAlt_apply, prodSel_zero, tens_vac]
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
      have hXt : prodSel (fun k a => w k (inL ℓ hℓ a)) (Fin.cons true τ) =
          crOp (fun a => w 0 (inL ℓ hℓ a)) *ᵥ prodSel (fun k a => w k.succ (inL ℓ hℓ a)) τ := by
        rw [prodSel_succ]; simp
      have hYt : prodSel (fun k b => w k (inR ℓ hℓ b))
          (fun k => !(Fin.cons true τ : Fin (M + 1) → Bool) k) =
          prodSel (fun k b => w k.succ (inR ℓ hℓ b)) (fun k => !τ k) := by
        rw [prodSel_succ]; simp
      have hXf : prodSel (fun k a => w k (inL ℓ hℓ a)) (Fin.cons false τ) =
          prodSel (fun k a => w k.succ (inL ℓ hℓ a)) τ := by
        rw [prodSel_succ]; simp
      have hYf : prodSel (fun k b => w k (inR ℓ hℓ b))
          (fun k => !(Fin.cons false τ : Fin (M + 1) → Bool) k) =
          crOp (fun b => w 0 (inR ℓ hℓ b)) *ᵥ
            prodSel (fun k b => w k.succ (inR ℓ hℓ b)) (fun k => !τ k) := by
        rw [prodSel_succ]; simp
      have hpar : (fun a => parC a * prodSel (fun k a => w k.succ (inL ℓ hℓ a)) τ a) =
          (∏ k, if τ k then (-1 : ℂ) else 1) • prodSel (fun k a => w k.succ (inL ℓ hℓ a)) τ :=
        parMul_prodSel _ _
      rw [Matrix.mulVec_smul, crOp_tens, hXt, hYt, hXf, hYf, hpar, tens_smul_left, smul_add,
        smul_smul]
      simp only [Fin.cons_zero, Fin.cons_succ, ite_true, Bool.false_eq_true, ite_false]

/-! ### The reduced density matrix of a sum of product vectors -/

/-- The reduced density matrix of `∑_i c_i X_i ⊗ Y_i` with orthogonal `Y_i` and `|c_i| = 1` is
`∑_i ‖Y_i‖² |X_i⟩⟨X_i|`. -/
lemma rdm_eq_sum (ℓ : ℕ) (hℓ : ℓ ≤ N) {ι : Type*} [Fintype ι] [DecidableEq ι] (ψ : Cfg N → ℂ)
    (c : ι → ℂ) (X : ι → Cfg ℓ → ℂ) (Y : ι → Cfg (N - ℓ) → ℂ) (μ : ι → ℂ)
    (hψ : ψ = ∑ i, c i • tens ℓ hℓ (X i) (Y i)) (hc : ∀ i, ‖c i‖ = 1)
    (hY : ∀ i j, star (Y j) ⬝ᵥ Y i = if i = j then μ i else 0) (a a' : Cfg ℓ) :
    rdm ℓ hℓ ψ a a' = ∑ i, μ i * (X i a * star (X i a')) := by
  simp only [rdm, hψ, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, tens, Equiv.apply_symm_apply]
  calc ∑ b, (∑ i, c i * (X i a * Y i b)) * star (∑ j, c j * (X j a' * Y j b))
      = ∑ b, ∑ i, ∑ j, c i * (X i a * Y i b) * star (c j * (X j a' * Y j b)) := by
        simp only [star_sum, Finset.sum_mul_sum]
    _ = ∑ i, ∑ j, ∑ b, c i * (X i a * Y i b) * star (c j * (X j a' * Y j b)) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ i, ∑ j, (c i * star (c j) * X i a * star (X j a')) * (star (Y j) ⬝ᵥ Y i) := by
        refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
        simp only [dotProduct, Pi.star_apply, Finset.mul_sum, star_mul']
        exact Finset.sum_congr rfl fun b _ => by ring
    _ = ∑ i, μ i * (X i a * star (X i a')) := by
        simp only [hY, mul_ite, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]
        refine Finset.sum_congr rfl fun i _ => ?_
        have h1 : c i * star (c i) = 1 := by
          rw [Complex.star_def, Complex.mul_conj', hc i]; simp
        linear_combination (μ i * X i a * star (X i a')) * h1

/-- The number of particles in the left segment is the number of particles of the left factor. -/
lemma count_segment (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) :
    count (segment ℓ) n = count univ (splitCfg ℓ hℓ n).1 := by
  unfold count
  rw [sum_segment ℓ hℓ]
  rfl

/-- `P(N_A = k)` as a partial trace of the reduced density matrix. -/
lemma numDist_eq_sum_rdm (ℓ : ℕ) (hℓ : ℓ ≤ N) (ψ : Cfg N → ℂ) (k : ℕ) :
    numDist (segment ℓ) ψ k =
      ∑ a ∈ univ.filter (fun a : Cfg ℓ => count univ a = k), (rdm ℓ hℓ ψ a a).re := by
  unfold numDist
  rw [Finset.sum_filter, Finset.sum_filter, ← (splitCfg ℓ hℓ).symm.sum_comp,
    Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun a _ => ?_
  simp only [count_segment ℓ hℓ, Equiv.apply_symm_apply]
  split_ifs with h
  · simp only [rdm, Complex.re_sum]
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [Complex.star_def, Complex.mul_conj, Complex.ofReal_re]
  · simp

/-! ### Product laws -/

/-- The product weights `∏_k (σ_k ? p_k : 1 - p_k)` sum to one. -/
lemma sum_prod_sel {M : ℕ} (p : Fin M → ℝ) :
    ∑ σ : Fin M → Bool, ∏ k, (if σ k then p k else 1 - p k) = 1 := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [sum_fun_succ]
    simp only [Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, ite_true, Bool.false_eq_true,
      ite_false]
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, ih (fun k => p k.succ)]
    ring

/-- **Entropy of a product law**: `∑_σ η(∏_k (σ_k ? p_k : 1 - p_k)) = ∑_k h(p_k)`. -/
lemma sum_negMulLog_prod_sel {M : ℕ} (p : Fin M → ℝ) :
    ∑ σ : Fin M → Bool, Real.negMulLog (∏ k, (if σ k then p k else 1 - p k)) =
      ∑ k, Real.binEntropy (p k) := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [sum_fun_succ, Fin.sum_univ_succ, ← ih (fun k => p k.succ)]
    simp only [Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, ite_true, Bool.false_eq_true,
      ite_false, Real.negMulLog_mul]
    have h1 := sum_prod_sel (fun k => p k.succ)
    rw [Real.binEntropy_eq_negMulLog_add_negMulLog_one_sub]
    simp only [Finset.sum_add_distrib, ← Finset.sum_mul, ← Finset.mul_sum, h1]
    ring

/-- Subsets of `Fin M` as Boolean functions. -/
def boolEquiv (M : ℕ) : (Fin M → Bool) ≃ Finset (Fin M) where
  toFun σ := univ.filter fun k => σ k = true
  invFun T k := decide (k ∈ T)
  left_inv σ := by funext k; simp
  right_inv T := by ext k; simp

/-- The law of the number of selected orbitals is the Poisson-binomial law. -/
lemma sum_sel_eq_poissonBinomial {M : ℕ} (ν : Fin M → ℝ) (n : ℕ) :
    ∑ σ : Fin M → Bool, (if (∑ j, if σ j then 1 else 0) = n then
      ∏ j, (if σ j then ν j else 1 - ν j) else 0) = poissonBinomial ν n := by
  rw [poissonBinomial, Finset.sum_filter, ← (boolEquiv M).sum_comp]
  refine Finset.sum_congr rfl fun σ _ => ?_
  have hcard : ((boolEquiv M) σ).card = ∑ j, if σ j then 1 else 0 := by
    show (univ.filter fun k => σ k = true).card = _
    rw [Finset.card_filter]
  have hprod : (∏ i ∈ (boolEquiv M) σ, ν i) * ∏ i ∈ univ \ (boolEquiv M) σ, (1 - ν i) =
      ∏ j, (if σ j then ν j else 1 - ν j) := by
    rw [Finset.prod_ite, Finset.filter_not]
    rfl
  rw [hcard, hprod]

open Polynomial in
/-- The generating polynomial `∏_i (ν_i X + 1 - ν_i)` of the Poisson-binomial law. -/
def pbPoly {ι : Type*} [Fintype ι] (ν : ι → ℝ) : ℝ[X] := ∏ i, (C (ν i) * X + C (1 - ν i))

open Polynomial in
/-- The coefficients of `pbPoly ν` are the Poisson-binomial probabilities. -/
lemma coeff_pbPoly {ι : Type*} [Fintype ι] [DecidableEq ι] (ν : ι → ℝ) (n : ℕ) :
    (pbPoly ν).coeff n = poissonBinomial ν n := by
  rw [pbPoly, Finset.prod_add, Finset.powerset_univ, poissonBinomial, Finset.sum_filter,
    Polynomial.finsetSum_coeff]
  refine Finset.sum_congr rfl fun T _ => ?_
  have h : (∏ i ∈ T, C (ν i) * X) * ∏ i ∈ univ \ T, C (1 - ν i) =
      C ((∏ i ∈ T, ν i) * ∏ i ∈ univ \ T, (1 - ν i)) * X ^ T.card := by
    simp only [Finset.prod_mul_distrib, Finset.prod_const, map_mul, map_prod]
    ring
  rw [h, Polynomial.coeff_C_mul_X_pow]
  rcases eq_or_ne T.card n with h' | h'
  · simp [h']
  · simp [h', Ne.symm h']

open Polynomial in
/-- `pbPoly ν` has degree at most the number of parameters. -/
lemma natDegree_pbPoly {ι : Type*} [Fintype ι] (ν : ι → ℝ) :
    (pbPoly ν).natDegree ≤ Fintype.card ι := by
  unfold pbPoly
  refine (Polynomial.natDegree_prod_le univ _).trans ?_
  refine (Finset.sum_le_sum fun i _ => Polynomial.natDegree_linear_le).trans ?_
  simp

open Polynomial in
/-- **Generating function of the Poisson-binomial law**:
`∑_{k ≤ n} P(k) z^k = ∏_i (ν_i z + 1 - ν_i)` for `n ≥ #parameters`. -/
lemma genfun_poissonBinomial {ι : Type*} [Fintype ι] [DecidableEq ι] (ν : ι → ℝ) (n : ℕ)
    (hn : Fintype.card ι ≤ n) (z : ℂ) :
    ∑ k ∈ range (n + 1), (poissonBinomial ν k : ℂ) * z ^ k =
      ∏ i, ((ν i : ℂ) * z + (1 - ν i)) := by
  set q : ℂ[X] := (pbPoly ν).map Complex.ofRealHom with hq
  have hdeg : q.natDegree < n + 1 :=
    lt_of_le_of_lt ((Polynomial.natDegree_map_le).trans (natDegree_pbPoly ν)) (by omega)
  have heval : q.eval z = ∏ i, ((ν i : ℂ) * z + (1 - ν i)) := by
    simp [hq, pbPoly, Polynomial.map_prod, Polynomial.eval_prod]
  rw [← heval, Polynomial.eval_eq_sum_range' hdeg]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [hq, Polynomial.coeff_map, coeff_pbPoly]
  rfl

/-! ### Equal nonzero spectra -/

/-- Equal matrices have equal eigenvalues. -/
lemma eigenvalues_congr {ι : Type*} [Fintype ι] [DecidableEq ι] {ρ σ : Matrix ι ι ℂ}
    (hρ : ρ.IsHermitian) (hσ : σ.IsHermitian) (h : ρ = σ) : hρ.eigenvalues = hσ.eigenvalues := by
  subst h; rfl

/-- `B Bᴴ` and `Bᴴ B = diag(ν)` have the same nonzero spectrum: an identity of multisets padded
with zeros. -/
lemma mset_eigenvalues {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
    (C : Matrix ι ι ℂ) (hC : C.IsHermitian) (B : Matrix ι κ ℂ) (hCB : C = B * Bᴴ) (ν : κ → ℝ)
    (hG : Bᴴ * B = diagonal (fun k => (ν k : ℂ))) :
    Fintype.card κ • ({0} : Multiset ℝ) + univ.val.map hC.eigenvalues =
      Fintype.card ι • {0} + univ.val.map ν := by
  subst hCB
  have h := Matrix.charpoly_mul_comm' B Bᴴ
  rw [hG] at h
  have h1 : (Polynomial.X ^ Fintype.card κ * (B * Bᴴ).charpoly : Polynomial ℂ) ≠ 0 :=
    mul_ne_zero (pow_ne_zero _ Polynomial.X_ne_zero) (Matrix.charpoly_monic _).ne_zero
  have h2 : (Polynomial.X ^ Fintype.card ι *
      (diagonal (fun k => (ν k : ℂ))).charpoly : Polynomial ℂ) ≠ 0 :=
    mul_ne_zero (pow_ne_zero _ Polynomial.X_ne_zero) (Matrix.charpoly_monic _).ne_zero
  have hr := congrArg Polynomial.roots h
  rw [Polynomial.roots_mul h1, Polynomial.roots_mul h2, Polynomial.roots_X_pow,
    Polynomial.roots_X_pow, hC.roots_charpoly_eq_eigenvalues] at hr
  have hd : (diagonal (fun k => (ν k : ℂ))).charpoly.roots =
      univ.val.map (fun k => (ν k : ℂ)) := by
    rw [Matrix.charpoly_diagonal, Polynomial.roots_prod]
    · simp
    · simp [Finset.prod_ne_zero_iff, Polynomial.X_sub_C_ne_zero]
  rw [hd] at hr
  have := congrArg (Multiset.map Complex.re) hr
  simpa [Multiset.map_add, Multiset.map_nsmul, Multiset.map_map, Function.comp_def] using this

/-- Transfer of sums `∑ f` with `f 0 = 0` along a padded multiset identity. -/
lemma sum_eq_of_mset {ι κ : Type*} [Fintype ι] [Fintype κ] {e : ι → ℝ} {ν : κ → ℝ} {a b : ℕ}
    (h : a • ({0} : Multiset ℝ) + univ.val.map e = b • {0} + univ.val.map ν) {R : Type*}
    [AddCommMonoid R] (f : ℝ → R) (hf : f 0 = 0) : ∑ i, f (e i) = ∑ k, f (ν k) := by
  have := congrArg (fun s : Multiset ℝ => (s.map f).sum) h
  simpa [Multiset.map_add, Multiset.sum_add, Multiset.map_nsmul, Multiset.sum_nsmul, hf,
    Multiset.map_map, Function.comp_def] using this

/-- Transfer of products `∏ F` with `F 0 = 1` along a padded multiset identity. -/
lemma prod_eq_of_mset {ι κ : Type*} [Fintype ι] [Fintype κ] {e : ι → ℝ} {ν : κ → ℝ} {a b : ℕ}
    (h : a • ({0} : Multiset ℝ) + univ.val.map e = b • {0} + univ.val.map ν) {R : Type*}
    [CommMonoid R] (F : ℝ → R) (hF : F 0 = 1) : ∏ i, F (e i) = ∏ k, F (ν k) := by
  have := congrArg (fun s : Multiset ℝ => (s.map F).prod) h
  simpa [Multiset.map_add, Multiset.prod_add, Multiset.map_nsmul, Multiset.prod_nsmul, hF,
    Multiset.map_map, Function.comp_def] using this

/-- **Weinstein–Aronszajn**: `det(1 + (z-1) B Bᴴ) = ∏_k (ν_k z + 1 - ν_k)` if `Bᴴ B = diag(ν)`. -/
lemma det_one_add_smul {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
    (B : Matrix ι κ ℂ) (ν : κ → ℝ) (hG : Bᴴ * B = diagonal (fun k => (ν k : ℂ))) (z : ℂ) :
    (1 + (z - 1) • (B * Bᴴ)).det = ∏ k, ((ν k : ℂ) * z + (1 - ν k)) := by
  rw [← Matrix.smul_mul, Matrix.det_one_add_mul_comm, Matrix.mul_smul, hG]
  have hD : (1 : Matrix κ κ ℂ) + (z - 1) • diagonal (fun k => (ν k : ℂ)) =
      diagonal (fun k => (ν k : ℂ) * z + (1 - ν k)) := by
    ext k l
    by_cases h : k = l
    · subst h; simp; ring
    · simp [h]
  rw [hD, Matrix.det_diagonal]

/-! ### The Schmidt decomposition of the evolved Fock state -/

/-- Setup: the evolved Fock state is `c ·` the Slater determinant of orthonormal canonical
orbitals `Φ`, with `K = ΦΦ†` and diagonal overlaps on the left segment and on its complement. -/
lemma slater_setup (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (ℓ : ℕ) :
    ∃ (M : ℕ) (Φ : Matrix (Fin N) (Fin M) ℂ) (c : ℂ) (ν : Fin M → ℝ), ‖c‖ = 1 ∧
      evo (dGamma A) t *ᵥ fockVec S = c • slaterAlt N M (fun k x => Φ x k) ∧
      corrMat A S t = Φ * Φᴴ ∧
      (∀ k l, ∑ x ∈ segment ℓ, star (Φ x k) * Φ x l = if k = l then (ν k : ℂ) else 0) ∧
      (∀ k l, ∑ x ∈ univ \ segment ℓ, star (Φ x k) * Φ x l =
        if k = l then ((1 - ν k : ℝ) : ℂ) else 0) := by
  classical
  obtain ⟨M, z, c, hz, hzS, hc, hv⟩ := evo_fockVec_eq A hA t S
  set Φ : Matrix (Fin N) (Fin M) ℂ := Matrix.of fun x k => evo A t x (z k) with hΦdef
  have hΦ : Φᴴ * Φ = 1 := by
    ext k l
    have h := congrFun (congrFun (evo_unitary A hA t) (z k)) (z l)
    simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, Matrix.one_apply] at h ⊢
    simp only [hΦdef, Matrix.of_apply]
    rw [h]
    exact if_congr hz.eq_iff rfl rfl
  have hC : corrMat A S t = Φ * Φᴴ := by
    ext x y
    rw [corrMat, Matrix.mul_apply, Matrix.mul_apply]
    simp only [Matrix.mul_diagonal, Matrix.conjTranspose_apply, mul_ite, mul_one, mul_zero,
      ite_mul, zero_mul, hΦdef, Matrix.of_apply]
    rw [Finset.sum_ite_mem, Finset.univ_inter]
    have hS : S = Finset.univ.image z := by
      ext x; simp [hzS x]
    rw [hS, Finset.sum_image (fun k _ l _ h => hz h)]
  obtain ⟨V, hVu, ν, h1, h2⟩ := canonical_orbitals Φ hΦ (segment ℓ)
  have hV1 : star V * V = 1 := Matrix.mem_unitaryGroup_iff'.mp hVu
  have hV2 : V * star V = 1 := Matrix.mem_unitaryGroup_iff.mp hVu
  have hCC : (Φ * V) * (Φ * V)ᴴ = Φ * Φᴴ := by
    rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc V,
      ← Matrix.star_eq_conjTranspose, hV2, Matrix.one_mul]
  have hdet : ‖V.det‖ = 1 := by
    have h := congrArg Matrix.det hV1
    rw [Matrix.det_mul, Matrix.star_eq_conjTranspose, Matrix.det_conjTranspose, Matrix.det_one,
      Complex.star_def, Complex.conj_mul'] at h
    have h' : ‖V.det‖ ^ 2 = 1 := by exact_mod_cast h
    exact (pow_eq_one_iff_of_nonneg (norm_nonneg _) two_ne_zero).mp h'
  have hdet0 : V.det ≠ 0 := by
    intro h0; rw [h0, norm_zero] at hdet; exact zero_ne_one hdet
  refine ⟨M, Φ * V, c * (V.det)⁻¹, ν, ?_, ?_, ?_, h1, h2⟩
  · rw [norm_mul, norm_inv, hc, hdet]; norm_num
  · rw [hv, slaterAlt_mul Φ V, smul_smul, mul_assoc, inv_mul_cancel₀ hdet0, mul_one]
    rfl
  · rw [hC, hCC]

/-- **Schmidt decomposition of a Slater determinant across the cut** (paper: proof of
Lemma `lem:slater`(a)), summarized: there are canonical occupations `ν_k ∈ [0, 1]` with
`C_A = B Bᴴ`, `Bᴴ B = diag(ν)` (so `C_A` and `diag(ν)` have the same nonzero spectrum),
`S(ρ_A) = ∑_k h(ν_k)`, the law of `N_A` is the Poisson-binomial law of `ν`, and
`∑_k √(ν_k(1-ν_k)) ≤ ∑_{x ∈ A, y ∉ A} |K_{xy}|`. -/
lemma slater_core (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    ∃ (M : ℕ) (B : Matrix (Fin ℓ) (Fin M) ℂ) (ν : Fin M → ℝ),
      corrA A S t ℓ hℓ = B * Bᴴ ∧ Bᴴ * B = diagonal (fun k => (ν k : ℂ)) ∧
      (∀ k, 0 ≤ ν k ∧ ν k ≤ 1) ∧
      ∑ k, √(ν k * (1 - ν k)) ≤ ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ ∧
      vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) =
        ∑ k, Real.binEntropy (ν k) ∧
      ∀ n : ℕ, numDist (segment ℓ) (evo (dGamma A) t *ᵥ fockVec S) n = poissonBinomial ν n := by
  classical
  obtain ⟨M, Φ, c, ν, hc, hψ, hC, h1, h2⟩ := slater_setup A hA S t ℓ
  -- the orbitals on the two sides of the cut
  set a : Fin M → Fin ℓ → ℂ := fun k x => Φ (inL ℓ hℓ x) k with ha_def
  set b : Fin M → Fin (N - ℓ) → ℂ := fun k y => Φ (inR ℓ hℓ y) k with hb_def
  have ha : ∀ k l, star (a k) ⬝ᵥ a l = if k = l then (ν k : ℂ) else 0 := by
    intro k l; rw [← h1 k l, sum_segment ℓ hℓ]; rfl
  have hb : ∀ k l, star (b k) ⬝ᵥ b l = if k = l then ((1 - ν k : ℝ) : ℂ) else 0 := by
    intro k l; rw [← h2 k l, sum_compl_segment ℓ hℓ]; rfl
  have hν0 : ∀ k, 0 ≤ ν k := by
    intro k
    have h := ha k k
    rw [ite_eq_left rfl, star_dotProduct_self] at h
    have h' : nsq (a k) = ν k := by exact_mod_cast h
    rw [← h']; exact nsq_nonneg _
  have hν1 : ∀ k, ν k ≤ 1 := by
    intro k
    have h := hb k k
    rw [ite_eq_left rfl, star_dotProduct_self] at h
    have h' : nsq (b k) = 1 - ν k := by exact_mod_cast h
    linarith [nsq_nonneg (b k)]
  set B : Matrix (Fin ℓ) (Fin M) ℂ := Matrix.of fun x k => Φ (inL ℓ hℓ x) k with hB_def
  have hCB : corrA A S t ℓ hℓ = B * Bᴴ := by
    ext x y
    simp [corrA, hC, Matrix.mul_apply, hB_def, Matrix.conjTranspose_apply]
  have hBB : Bᴴ * B = diagonal (fun k => (ν k : ℂ)) := by
    ext k l
    rw [Matrix.mul_apply, Matrix.diagonal_apply, ← h1 k l, sum_segment ℓ hℓ]
    rfl
  have hsqrt := sum_sqrt_nu_le Φ (segment ℓ) ν h1 h2
  rw [← hC] at hsqrt
  -- the Schmidt decomposition
  obtain ⟨s, hs, hexp⟩ := slater_expand ℓ hℓ (fun k x => Φ x k)
  set X : (Fin M → Bool) → Cfg ℓ → ℂ := fun σ => prodSel a σ with hX_def
  set Y : (Fin M → Bool) → Cfg (N - ℓ) → ℂ := fun σ => prodSel b (fun k => !σ k) with hY_def
  set ψ := evo (dGamma A) t *ᵥ fockVec S with hψ_def
  have hψ' : ψ = ∑ σ, (c * s σ) • tens ℓ hℓ (X σ) (Y σ) := by
    rw [hψ, hexp, Finset.smul_sum]
    simp only [smul_smul]
    rfl
  have hcs : ∀ σ, ‖c * s σ‖ = 1 := fun σ => by rw [norm_mul, hc, hs, one_mul]
  have hX : ∀ σ τ, star (X σ) ⬝ᵥ X τ =
      if σ = τ then ∏ k, (if σ k then (ν k : ℂ) else 1) else 0 :=
    fun σ τ => gram_prodSel a _ ha σ τ
  set μ : (Fin M → Bool) → ℝ := fun σ => ∏ k, (if σ k then 1 else 1 - ν k) with hμ_def
  have hY : ∀ σ τ, star (Y τ) ⬝ᵥ Y σ = if σ = τ then (μ σ : ℂ) else 0 := by
    intro σ τ
    have hg := gram_prodSel b _ hb (fun k => !τ k) (fun k => !σ k)
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
  have hρ : ∀ x x', rdm ℓ hℓ ψ x x' = ∑ σ, (μ σ : ℂ) * (X σ x * star (X σ x')) :=
    rdm_eq_sum ℓ hℓ ψ (fun σ => c * s σ) X Y (fun σ => (μ σ : ℂ)) hψ' hcs hY
  have hμ0 : ∀ σ, 0 ≤ μ σ := fun σ =>
    Finset.prod_nonneg fun k _ => by split_ifs <;> linarith [hν1 k]
  have hXn : ∀ σ, nsq (X σ) = ∏ k, (if σ k then ν k else 1) := by
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
    rw [hXn, hμ_def, hlam_def]
    simp only
    rw [← Finset.prod_mul_distrib]
    refine Finset.prod_congr rfl fun k _ => ?_
    cases σ k <;> simp
  -- the entropy
  set Z : Matrix (Cfg ℓ) (Fin M → Bool) ℂ :=
    Matrix.of fun x σ => ((√(μ σ) : ℝ) : ℂ) * X σ x with hZ_def
  have hsq : ∀ σ, ((√(μ σ) : ℝ) : ℂ) * ((√(μ σ) : ℝ) : ℂ) = (μ σ : ℂ) := fun σ => by
    rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hμ0 σ)]
  have hZZ : rdm ℓ hℓ ψ = Z * Zᴴ := by
    ext x x'
    rw [hρ, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp only [hZ_def, Matrix.of_apply, Matrix.conjTranspose_apply, star_mul', Complex.star_def,
      Complex.conj_ofReal]
    rw [← hsq σ]
    ring
  have hZZ' : Zᴴ * Z = diagonal (fun σ => (lam σ : ℂ)) := by
    ext σ τ
    rw [Matrix.mul_apply, Matrix.diagonal_apply]
    have e : ∑ x, Zᴴ σ x * Z x τ =
        ((√(μ σ) : ℝ) : ℂ) * ((√(μ τ) : ℝ) : ℂ) * (star (X σ) ⬝ᵥ X τ) := by
      rw [dotProduct, Finset.mul_sum]
      refine Finset.sum_congr rfl fun x _ => ?_
      simp only [hZ_def, Matrix.of_apply, Matrix.conjTranspose_apply, star_mul', Pi.star_apply,
        Complex.star_def, Complex.conj_ofReal]
      ring
    rw [e, hX]
    split_ifs with h
    · subst h
      rw [hsq, ← hlam σ, hXn, Complex.ofReal_mul, Complex.ofReal_prod]
      congr 1
      exact Finset.prod_congr rfl fun k _ => by split_ifs <;> simp
    · simp
  have hHd : (diagonal (fun σ => (lam σ : ℂ))).IsHermitian :=
    hZZ' ▸ isHermitian_conjTranspose_mul_self Z
  have hent : vonNeumann (rdm ℓ hℓ ψ) (rdm_isHermitian ℓ hℓ ψ) = ∑ k, Real.binEntropy (ν k) := by
    rw [CutLocalAux.vonNeumann_congr _ (isHermitian_mul_conjTranspose_self Z) hZZ,
      vonNeumann_mul_conjTranspose_comm, CutLocalAux.vonNeumann_congr _ hHd hZZ',
      vonNeumann_diagonal, ← sum_negMulLog_prod_sel ν]
  -- the law of `N_A`
  have hnum : ∀ n : ℕ, numDist (segment ℓ) ψ n = poissonBinomial ν n := by
    intro n
    rw [numDist_eq_sum_rdm ℓ hℓ ψ n, ← sum_sel_eq_poissonBinomial]
    simp only [hρ, Complex.re_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun σ _ => ?_
    have hre : ∀ x, ((μ σ : ℂ) * (X σ x * star (X σ x))).re = μ σ * ‖X σ x‖ ^ 2 := by
      intro x
      rw [Complex.star_def, Complex.mul_conj', ← Complex.ofReal_pow, ← Complex.ofReal_mul,
        Complex.ofReal_re]
    simp only [hre, ← Finset.mul_sum]
    have hsupp : IsNum (∑ k, if σ k then 1 else 0) (X σ) := isNum_prodSel a σ
    by_cases hn : (∑ j, if σ j then 1 else 0) = n
    · rw [ite_eq_left hn, ← hlam σ, nsq]
      congr 1
      rw [Finset.sum_filter_of_ne]
      intro x _ hx
      by_contra hcount
      apply hx
      simp [hsupp x (by rw [hn]; exact hcount)]
    · rw [ite_eq_right hn]
      apply mul_eq_zero_of_right
      apply Finset.sum_eq_zero
      intro x hx
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
      simp [hsupp x (by rw [hx]; exact_mod_cast (Ne.symm hn))]
  exact ⟨M, B, ν, hCB, hBB, fun k => ⟨hν0 k, hν1 k⟩, hsqrt, hent, hnum⟩

end SlaterExactAux

/-- `C_A` is Hermitian. -/
theorem corrA_isHermitian (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : (corrA A S t ℓ hℓ).IsHermitian := by
  sorry

/-- The eigenvalues of `C_A` lie in `[0, 1]`. -/
theorem eig_corrA_mem (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (i : Fin ℓ) :
    0 ≤ (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i ∧
      (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i ≤ 1 := by
  sorry

/-- **Entanglement entropy of a Slater determinant** (paper: Lemma `lem:slater`(a)):
`S(ρ_A) = ∑_i h(ν_i)`. -/
theorem vonNeumann_rdm_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) =
      ∑ i, Real.binEntropy ((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i) := by
  sorry

/-- **The law of `N_A` is Poisson-binomial** (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b)). -/
theorem numDist_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (k : ℕ) :
    numDist (segment ℓ) (evo (dGamma A) t *ᵥ fockVec S) k =
      poissonBinomial (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues k := by
  sorry

/-- **Generating function** (paper: Lemma `lem:quasifree`(b)):
`∑_k P(N_A = k) z^k = det(1 + (z - 1) C_A)`. -/
theorem genfun_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (z : ℂ) :
    ∑ k ∈ range (N + 1), (numDist (segment ℓ) (evo (dGamma A) t *ᵥ fockVec S) k : ℂ) * z ^ k =
      (1 + (z - 1) • corrA A S t ℓ hℓ).det := by
  sorry

/-- **Trace-norm bound** (paper: Lemma `lem:slater`(b)):
`∑_i √(ν_i (1 - ν_i)) = ‖P_B K P_A‖₁ ≤ ∑_{x ∈ A, y ∉ A} |K_{xy}|`. -/
theorem sum_sqrt_eig_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    ∑ i, √((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i *
        (1 - (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i)) ≤
      ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  sorry

end DFL
