import DFL.SlaterExact

/-!
# The law of `N_R` for Slater determinants and arbitrary regions (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b))

For the evolved Fock state `ψ = Γ(e^{-itA}) |S⟩` (a Slater determinant with one-body projection
`K = corrMat A S t`) and an **arbitrary** set of sites `R` (e.g. a bulk interval), let `C_R = P_R K P_R`
(`corrR`), with eigenvalues `ν_i ∈ [0, 1]`. Then the law of `N_R` is the Poisson-binomial law of the `ν_i`
(`numDist_slater_gen`), with generating function `det(1 + (z - 1) C_R)` (`genfun_slater_gen`). For left
segments this is `numDist_slater`, `genfun_slater` of `DFL.SlaterExact`.

## Proof outline (auxiliary lemmas in `DFL.SlaterGeneralAux`)

1. `slater_setup_gen`: `ψ = c · c†(w_1) ⋯ c†(w_M) |0⟩` (`|c| = 1`) for canonical orbitals `w_k`
   (`AreaLawAux.canonical_orbitals` for the region `R`) with `K = ∑_k |w_k⟩⟨w_k|`,
   `⟨P_R w_k, P_R w_l⟩ = ν_k δ_kl` and `⟨P_{R^c} w_k, P_{R^c} w_l⟩ = (1 - ν_k) δ_kl`.
2. Multilinearity (`AlternatingMap.map_add_univ`) with `w_k = a_k + b_k`, `a_k = P_R w_k`,
   `b_k = P_{R^c} w_k`, expands `ψ = c ∑_{s ⊆ {1..M}} T_s`, where `T_s` is the Slater determinant of the
   orbitals `a_k` (`k ∈ s`) and `b_k` (`k ∉ s`).
3. `gram_slaterAlt_diag`: if `⟨u_k, v_l⟩ = 0` for `k ≠ l`, the Slater determinants of `u` and `v` have
   overlap `∏_k ⟨u_k, v_k⟩`; hence the `T_s` are orthogonal with
   `‖T_s‖² = ∏_{k ∈ s} ν_k ∏_{k ∉ s} (1 - ν_k)`.
4. `isNumR_slaterAlt`: `T_s` has exactly `|s|` particles in `R`. Hence `P(N_R = n) = ∑_{|s| = n} ‖T_s‖²`
   is the Poisson-binomial law of `ν` (`slater_core_gen`).
5. `C_R = B B†` with `B† B = diag(ν)`, so (as in `DFL.SlaterExact`) `C_R` and `diag(ν)` have the same
   nonzero spectrum, which transfers the law and gives `det(1 + (z - 1) C_R) = ∏_k (1 + (z - 1) ν_k)`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `C_R = P_R K P_R`: the correlation matrix restricted to the sites of `R`. -/
def corrR (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (R : Finset (Fin N)) :
    Matrix R R ℂ :=
  (corrMat A S t).submatrix (fun x : R => (x : Fin N)) (fun x : R => (x : Fin N))

namespace SlaterGeneralAux

open AreaLawAux SlaterExactAux

variable {K : ℕ}

/-! ### Slater determinants: a Gram identity -/

/-- `slaterAlt` is `prodSel` with every orbital selected. -/
lemma slaterAlt_eq_prodSel {M : ℕ} (u : Fin M → Fin K → ℂ) :
    slaterAlt K M u = prodSel u (fun _ => true) := by
  rw [slaterAlt_apply, prodSel]
  simp

/-- Peeling off the first orbital of a Slater determinant. -/
lemma slaterAlt_succ {M : ℕ} (u : Fin (M + 1) → Fin K → ℂ) :
    slaterAlt K (M + 1) u = crOp (u 0) *ᵥ slaterAlt K M (fun k => u k.succ) := by
  rw [slaterAlt_apply, slaterAlt_apply, List.ofFn_succ, List.prod_cons, ← Matrix.mulVec_mulVec]

/-- `c(v)` kills a Slater determinant of orbitals orthogonal to `v`. -/
lemma adj_slaterAlt {M : ℕ} (u : Fin M → Fin K → ℂ) (v : Fin K → ℂ)
    (hv : ∀ k, star v ⬝ᵥ u k = 0) : (crOp v)ᴴ *ᵥ slaterAlt K M u = 0 := by
  rw [slaterAlt_eq_prodSel]
  exact adj_prodSel u _ v (fun k _ => hv k)

/-- **Gram identity for diagonal overlaps**: if `⟨u_k, v_l⟩ = 0` for `k ≠ l`, then
`⟨c†(u_1) ⋯ c†(u_M) |0⟩, c†(v_1) ⋯ c†(v_M) |0⟩⟩ = ∏_k ⟨u_k, v_k⟩`. -/
lemma gram_slaterAlt_diag {M : ℕ} (u v : Fin M → Fin K → ℂ)
    (huv : ∀ k l, k ≠ l → star (u k) ⬝ᵥ v l = 0) :
    star (slaterAlt K M u) ⬝ᵥ slaterAlt K M v = ∏ k, star (u k) ⬝ᵥ v k := by
  induction M with
  | zero =>
    simp only [slaterAlt_apply, List.ofFn_zero, List.prod_nil, Matrix.one_mulVec,
      Finset.univ_eq_empty, Finset.prod_empty]
    rw [star_dotProduct_self, fockVec, nsq_pi_single, Complex.ofReal_one]
  | succ M ih =>
    have ih' := ih (fun k => u k.succ) (fun k => v k.succ)
      (fun k l hkl => huv _ _ fun h => hkl (Fin.succ_injective _ h))
    have hkill : (crOp (u 0))ᴴ *ᵥ slaterAlt K M (fun k => v k.succ) = 0 :=
      adj_slaterAlt _ _ fun k => huv 0 k.succ (Fin.succ_ne_zero k).symm
    have hcar := crOp_CAR_gen (u 0) (v 0)
    rw [← eq_sub_iff_add_eq] at hcar
    rw [slaterAlt_succ, slaterAlt_succ, Fin.prod_univ_succ, dot_crOp_left, Matrix.mulVec_mulVec,
      hcar, Matrix.sub_mulVec, ← Matrix.mulVec_mulVec, hkill, Matrix.mulVec_zero, sub_zero,
      Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, ih', smul_eq_mul]

/-! ### Particle number in a region -/

/-- `Z` is supported on configurations with exactly `m` particles in `R`. -/
def IsNumR (R : Finset (Fin K)) (m : ℕ) (Z : Cfg K → ℂ) : Prop :=
  ∀ n, count R n ≠ (m : ℝ) → Z n = 0

/-- Emptying an occupied mode of `R` removes one particle from `R`. -/
lemma count_update_mem (R : Finset (Fin K)) (n : Cfg K) (x : Fin K) (hx : n x = 1)
    (hxR : x ∈ R) : count R n = count R (Function.update n x 0) + 1 := by
  unfold count
  rw [← Finset.add_sum_erase _ _ hxR, ← Finset.add_sum_erase _ _ hxR, Function.update_self, hx]
  have : ∑ k ∈ R.erase x, occ (Function.update n x 0 k) = ∑ k ∈ R.erase x, occ (n k) :=
    Finset.sum_congr rfl fun k hk => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hk)]
  rw [this]
  simp [occ]
  ring

/-- Emptying a mode outside `R` does not change the number of particles in `R`. -/
lemma count_update_not_mem (R : Finset (Fin K)) (n : Cfg K) (x : Fin K) (hxR : x ∉ R) :
    count R (Function.update n x 0) = count R n := by
  unfold count
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [Function.update_of_ne (fun h : k = x => hxR (h ▸ hk))]

/-- `c†(u)` with `u` supported in `R` raises the particle number in `R` by one. -/
lemma isNumR_crOp_in (R : Finset (Fin K)) {m : ℕ} (u : Fin K → ℂ) (hu : ∀ x ∉ R, u x = 0)
    (Z : Cfg K → ℂ) (hZ : IsNumR R m Z) : IsNumR R (m + 1) (crOp u *ᵥ Z) := by
  intro n hn
  rw [crOp_mulVec_eq, Finset.sum_apply]
  refine Finset.sum_eq_zero fun x _ => ?_
  rw [Pi.smul_apply, cre_mulVec_apply]
  by_cases hxR : x ∈ R
  · split_ifs with h
    · rw [hZ _ ?_, mul_zero, smul_zero]
      intro hc; apply hn; rw [count_update_mem R n x h hxR, hc]; push_cast; ring
    · simp
  · rw [hu x hxR, zero_smul]

/-- `c†(u)` with `u` supported outside `R` preserves the particle number in `R`. -/
lemma isNumR_crOp_out (R : Finset (Fin K)) {m : ℕ} (u : Fin K → ℂ) (hu : ∀ x ∈ R, u x = 0)
    (Z : Cfg K → ℂ) (hZ : IsNumR R m Z) : IsNumR R m (crOp u *ᵥ Z) := by
  intro n hn
  rw [crOp_mulVec_eq, Finset.sum_apply]
  refine Finset.sum_eq_zero fun x _ => ?_
  rw [Pi.smul_apply, cre_mulVec_apply]
  by_cases hxR : x ∈ R
  · rw [hu x hxR, zero_smul]
  · split_ifs with h
    · rw [hZ _ ?_, mul_zero, smul_zero]
      rwa [count_update_not_mem R n x hxR]
    · simp

/-- A Slater determinant of orbitals each supported in `R` (for `p k`) or outside `R` (for `¬ p k`)
has exactly `#{k | p k}` particles in `R`. -/
lemma isNumR_slaterAlt (R : Finset (Fin K)) {M : ℕ} (u : Fin M → Fin K → ℂ) (p : Fin M → Prop)
    [DecidablePred p] (hin : ∀ k, p k → ∀ x ∉ R, u k x = 0)
    (hout : ∀ k, ¬ p k → ∀ x ∈ R, u k x = 0) :
    IsNumR R (∑ k, if p k then 1 else 0) (slaterAlt K M u) := by
  induction M with
  | zero =>
    intro n hn
    rw [slaterAlt_apply]
    simp only [List.ofFn_zero, List.prod_nil, Matrix.one_mulVec, fockVec, Pi.single_apply]
    rw [ite_eq_right]
    rintro rfl
    apply hn
    simp [count, fockCfg, occ]
  | succ M ih =>
    have ih' := ih (fun k => u k.succ) (fun k => p k.succ) (fun k hk => hin _ hk)
      (fun k hk => hout _ hk)
    rw [slaterAlt_succ, Fin.sum_univ_succ]
    by_cases h0 : p 0
    · rw [ite_eq_left h0, add_comm]
      exact isNumR_crOp_in R _ (hin 0 h0) _ ih'
    · rw [ite_eq_right h0, zero_add]
      exact isNumR_crOp_out R _ (hout 0 h0) _ ih'

/-- `∏_k (k ∈ s ? ν_k : 1 - ν_k) = ∏_{k ∈ s} ν_k ∏_{k ∉ s} (1 - ν_k)`. -/
lemma prod_ite_mem_eq {M : ℕ} (ν : Fin M → ℝ) (s : Finset (Fin M)) :
    ∏ k, (if k ∈ s then ν k else 1 - ν k) = (∏ k ∈ s, ν k) * ∏ k ∈ univ \ s, (1 - ν k) := by
  rw [← Finset.prod_mul_prod_compl s, Finset.compl_eq_univ_sdiff]
  congr 1
  · exact Finset.prod_congr rfl fun k hk => ite_eq_left hk
  · exact Finset.prod_congr rfl fun k hk => ite_eq_right (Finset.mem_sdiff.mp hk).2

/-! ### The law of `N_R` -/

/-- Setup (as `SlaterExactAux.slater_setup`, for an arbitrary region `R`): the evolved Fock state is
`c ·` the Slater determinant of orthonormal canonical orbitals `Φ`, with `K = ΦΦ†` and diagonal
overlaps on `R` and on its complement. -/
lemma slater_setup_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    ∃ (M : ℕ) (Φ : Matrix (Fin N) (Fin M) ℂ) (c : ℂ) (ν : Fin M → ℝ), ‖c‖ = 1 ∧
      evo (dGamma A) t *ᵥ fockVec S = c • slaterAlt N M (fun k x => Φ x k) ∧
      corrMat A S t = Φ * Φᴴ ∧
      (∀ k l, ∑ x ∈ R, star (Φ x k) * Φ x l = if k = l then (ν k : ℂ) else 0) ∧
      (∀ k l, ∑ x ∈ univ \ R, star (Φ x k) * Φ x l =
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
  obtain ⟨V, hVu, ν, h1, h2⟩ := canonical_orbitals Φ hΦ R
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

/-- **The law of `N_R` for a Slater determinant**, `R` arbitrary: there are canonical occupations
`ν_k ∈ [0, 1]` with `C_R = B Bᴴ`, `Bᴴ B = diag(ν)`, and the law of `N_R` is the Poisson-binomial law
of `ν`. -/
lemma slater_core_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    ∃ (M : ℕ) (B : Matrix R (Fin M) ℂ) (ν : Fin M → ℝ),
      corrR A S t R = B * Bᴴ ∧ Bᴴ * B = diagonal (fun k => (ν k : ℂ)) ∧
      (∀ k, 0 ≤ ν k ∧ ν k ≤ 1) ∧
      ∀ n : ℕ, numDist R (evo (dGamma A) t *ᵥ fockVec S) n = poissonBinomial ν n := by
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
  have dab : ∀ k l, star (a k) ⬝ᵥ b l = 0 := by
    intro k l
    simp only [dotProduct]
    refine Finset.sum_eq_zero fun x _ => ?_
    by_cases hx : x ∈ R <;> simp [ha_def, hb_def, hx]
  have dba : ∀ k l, star (b k) ⬝ᵥ a l = 0 := by
    intro k l
    simp only [dotProduct]
    refine Finset.sum_eq_zero fun x _ => ?_
    by_cases hx : x ∈ R <;> simp [ha_def, hb_def, hx]
  -- the canonical occupations lie in `[0, 1]`
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
  -- expansion of the Slater determinant
  set T : Finset (Fin M) → Cfg N → ℂ := fun s => slaterAlt N M (s.piecewise a b) with hT_def
  set w : Finset (Fin M) → ℝ := fun s => (∏ k ∈ s, ν k) * ∏ k ∈ univ \ s, (1 - ν k) with hw_def
  have hab : (fun k x => Φ x k) = a + b := by
    funext k x
    by_cases hx : x ∈ R <;> simp [ha_def, hb_def, hx]
  set ψ := evo (dGamma A) t *ᵥ fockVec S with hψ_def
  have hψ' : ψ = c • ∑ s, T s := by
    rw [hψ, hab, AlternatingMap.map_add_univ]
  -- the Gram matrix of the `T s`
  have huv : ∀ (s s' : Finset (Fin M)) (k l : Fin M), k ≠ l →
      star (s.piecewise a b k) ⬝ᵥ s'.piecewise a b l = 0 := by
    intro s s' k l hkl
    by_cases hk : k ∈ s <;> by_cases hl : l ∈ s'
    · rw [Finset.piecewise_eq_of_mem _ _ _ hk, Finset.piecewise_eq_of_mem _ _ _ hl, daa,
        ite_eq_right hkl]
    · rw [Finset.piecewise_eq_of_mem _ _ _ hk, Finset.piecewise_eq_of_notMem _ _ _ hl, dab]
    · rw [Finset.piecewise_eq_of_notMem _ _ _ hk, Finset.piecewise_eq_of_mem _ _ _ hl, dba]
    · rw [Finset.piecewise_eq_of_notMem _ _ _ hk, Finset.piecewise_eq_of_notMem _ _ _ hl, dbb,
        ite_eq_right hkl]
  have hkk_same : ∀ (s : Finset (Fin M)) (k : Fin M),
      star (s.piecewise a b k) ⬝ᵥ s.piecewise a b k =
        ((if k ∈ s then ν k else 1 - ν k : ℝ) : ℂ) := by
    intro s k
    by_cases hk : k ∈ s
    · rw [Finset.piecewise_eq_of_mem _ _ _ hk, daa, ite_eq_left rfl, ite_eq_left hk]
    · rw [Finset.piecewise_eq_of_notMem _ _ _ hk, dbb, ite_eq_left rfl, ite_eq_right hk]
  have hkk_diff : ∀ (s s' : Finset (Fin M)) (k : Fin M), ¬(k ∈ s ↔ k ∈ s') →
      star (s.piecewise a b k) ⬝ᵥ s'.piecewise a b k = 0 := by
    intro s s' k hk
    by_cases hk1 : k ∈ s <;> by_cases hk2 : k ∈ s'
    · exact absurd (iff_of_true hk1 hk2) hk
    · rw [Finset.piecewise_eq_of_mem _ _ _ hk1, Finset.piecewise_eq_of_notMem _ _ _ hk2, dab]
    · rw [Finset.piecewise_eq_of_notMem _ _ _ hk1, Finset.piecewise_eq_of_mem _ _ _ hk2, dba]
    · exact absurd (iff_of_false hk1 hk2) hk
  have hT : ∀ s s' : Finset (Fin M), star (T s) ⬝ᵥ T s' = if s = s' then (w s : ℂ) else 0 := by
    intro s s'
    simp only [hT_def]
    rw [gram_slaterAlt_diag _ _ (huv s s')]
    split_ifs with hss
    · subst hss
      simp only [hw_def]
      rw [← prod_ite_mem_eq, Complex.ofReal_prod]
      exact Finset.prod_congr rfl fun k _ => hkk_same s k
    · obtain ⟨k, hk⟩ : ∃ k, ¬(k ∈ s ↔ k ∈ s') := by
        by_contra h
        exact hss (Finset.ext fun k => not_not.mp fun hk => h ⟨k, hk⟩)
      exact Finset.prod_eq_zero (Finset.mem_univ k) (hkk_diff s s' k hk)
  -- `T s` has `|s|` particles in `R`
  have hTn : ∀ s : Finset (Fin M), IsNumR R s.card (T s) := by
    intro s
    have hin : ∀ k, k ∈ s → ∀ x ∉ R, s.piecewise a b k x = 0 := by
      intro k hk x hx
      rw [Finset.piecewise_eq_of_mem _ _ _ hk]
      simp [ha_def, hx]
    have hout : ∀ k, k ∉ s → ∀ x ∈ R, s.piecewise a b k x = 0 := by
      intro k hk x hx
      rw [Finset.piecewise_eq_of_notMem _ _ _ hk]
      simp [hb_def, hx]
    have h := isNumR_slaterAlt R (s.piecewise a b) (fun k => k ∈ s) hin hout
    have hcard : (∑ k, if k ∈ s then 1 else 0 : ℕ) = s.card := by simp
    intro n hn
    exact h n (by rw [hcard]; exact hn)
  -- the law of `N_R`
  have hnum : ∀ n : ℕ, numDist R ψ n = poissonBinomial ν n := by
    intro n
    set F := (univ : Finset (Finset (Fin M))).filter (fun s => s.card = n) with hF
    set X : Cfg N → ℂ := ∑ s ∈ F, T s with hX
    have hXsupp : ∀ m, count R m ≠ (n : ℝ) → X m = 0 := by
      intro m hm
      rw [hX, Finset.sum_apply]
      refine Finset.sum_eq_zero fun s hs => ?_
      have hs' : s.card = n := (Finset.mem_filter.mp hs).2
      exact hTn s m (by rw [hs']; exact hm)
    have hψX : ∀ m, count R m = (n : ℝ) → ψ m = c * X m := by
      intro m hm
      rw [hψ', Pi.smul_apply, smul_eq_mul, hX, Finset.sum_apply, Finset.sum_apply]
      congr 1
      symm
      refine Finset.sum_subset (Finset.subset_univ F) fun s _ hs => ?_
      have hs' : s.card ≠ n := fun h => hs (Finset.mem_filter.mpr ⟨Finset.mem_univ s, h⟩)
      exact hTn s m (by rw [hm]; exact_mod_cast hs'.symm)
    have h3 : numDist R ψ n = nsq X := by
      unfold numDist nsq
      calc ∑ m ∈ univ.filter (fun m => count R m = n), Complex.normSq (ψ m)
          = ∑ m ∈ univ.filter (fun m => count R m = n), ‖X m‖ ^ 2 := by
            refine Finset.sum_congr rfl fun m hm => ?_
            have hm' : count R m = n := (Finset.mem_filter.mp hm).2
            rw [hψX m hm', Complex.normSq_mul, Complex.normSq_eq_norm_sq c, hc, one_pow, one_mul,
              Complex.normSq_eq_norm_sq]
        _ = ∑ m, ‖X m‖ ^ 2 := by
            refine Finset.sum_subset (Finset.subset_univ _) fun m _ hm => ?_
            have hm' : count R m ≠ n := fun h =>
              hm (Finset.mem_filter.mpr ⟨Finset.mem_univ m, h⟩)
            rw [hXsupp m hm', norm_zero]
            norm_num
    have h4 : ((nsq X : ℝ) : ℂ) = ((∑ s ∈ F, w s : ℝ) : ℂ) := by
      rw [← star_dotProduct_self, hX, star_sum, sum_dotProduct, Complex.ofReal_sum]
      refine Finset.sum_congr rfl fun s hs => ?_
      rw [dotProduct_sum]
      simp only [hT]
      rw [Finset.sum_ite_eq F s (fun _ => (w s : ℂ)), ite_eq_left hs]
    have h5 : nsq X = ∑ s ∈ F, w s := by exact_mod_cast h4
    rw [h3, h5]
    rfl
  exact ⟨M, B, ν, hCB, hBB, fun k => ⟨hν0 k, hν1 k⟩, hnum⟩

end SlaterGeneralAux

/-- `C_R` is Hermitian. -/
theorem corrR_isHermitian (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) : (corrR A S t R).IsHermitian := by
  unfold Matrix.IsHermitian corrR
  rw [Matrix.conjTranspose_submatrix, FreeDynAux.corrMat_conjTranspose]

namespace SlaterGeneralAux

open SlaterExactAux

/-- The law of `N_R` is the Poisson-binomial law of the eigenvalues of `C_R`. -/
lemma numDist_eq_pb_eig (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) (k : ℕ) :
    numDist R (evo (dGamma A) t *ᵥ fockVec S) k =
      poissonBinomial (corrR_isHermitian A S t R).eigenvalues k := by
  obtain ⟨M, B, ν, hCB, hBB, -, hnum⟩ := slater_core_gen A hA S t R
  have hm := mset_eigenvalues _ (corrR_isHermitian A S t R) B hCB ν hBB
  have hp := prod_eq_of_mset hm
    (fun x => Polynomial.C x * Polynomial.X + Polynomial.C (1 - x)) (by simp)
  rw [hnum k, ← coeff_pbPoly, ← coeff_pbPoly, pbPoly, pbPoly, hp]

end SlaterGeneralAux

/-- The eigenvalues of `C_R` lie in `[0, 1]`. -/
theorem eig_corrR_mem (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (i : R) :
    0 ≤ (corrR_isHermitian A S t R).eigenvalues i ∧ (corrR_isHermitian A S t R).eigenvalues i ≤ 1 := by
  obtain ⟨M, B, ν, hCB, hBB, hν, -⟩ := SlaterGeneralAux.slater_core_gen A hA S t R
  have hm := SlaterExactAux.mset_eigenvalues _ (corrR_isHermitian A S t R) B hCB ν hBB
  simp only [Fintype.card_fin] at hm
  have hi : (corrR_isHermitian A S t R).eigenvalues i ∈
      M • ({0} : Multiset ℝ) + univ.val.map (corrR_isHermitian A S t R).eigenvalues :=
    Multiset.mem_add.mpr (Or.inr (Multiset.mem_map_of_mem _ (Finset.mem_univ_val i)))
  rw [hm, Multiset.mem_add] at hi
  rcases hi with h | h
  · rw [Multiset.mem_singleton.mp (Multiset.mem_of_mem_nsmul h)]
    norm_num
  · obtain ⟨k, -, hk⟩ := Multiset.mem_map.mp h
    rw [← hk]
    exact hν k

/-- **Generating function of `N_R`** (paper: Lemma `lem:quasifree`(b)), any region `R`:
`∑_k P(N_R = k) z^k = det(1 + (z - 1) C_R)`. -/
theorem genfun_slater_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (z : ℂ) :
    ∑ k ∈ range (N + 1), (numDist R (evo (dGamma A) t *ᵥ fockVec S) k : ℂ) * z ^ k =
      (1 + (z - 1) • corrR A S t R).det := by
  obtain ⟨M, B, ν, hCB, hBB, -, -⟩ := SlaterGeneralAux.slater_core_gen A hA S t R
  have hm := SlaterExactAux.mset_eigenvalues _ (corrR_isHermitian A S t R) B hCB ν hBB
  simp only [SlaterGeneralAux.numDist_eq_pb_eig A hA S t R]
  rw [SlaterExactAux.genfun_poissonBinomial _ N (by simpa using Finset.card_le_univ R) z,
    SlaterExactAux.prod_eq_of_mset hm (fun x => (x : ℂ) * z + (1 - x)) (by simp), hCB,
    SlaterExactAux.det_one_add_smul B ν hBB z]

/-- **The law of `N_R` is Poisson-binomial** (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b)), any
region `R`. -/
theorem numDist_slater_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (k : ℕ) :
    numDist R (evo (dGamma A) t *ᵥ fockVec S) k =
      poissonBinomial (corrR_isHermitian A S t R).eigenvalues k :=
  SlaterGeneralAux.numDist_eq_pb_eig A hA S t R k

end DFL
