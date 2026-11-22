import DFL.GaussLaw
import DFL.FreeDyn

/-!
# The duality dictionary, dressed fermions and their dynamics
(paper: Theorem `thm:duality`(c),(d),(e))

`U` is the permutation `|n, q⟩ ↦ |n⟩ ⊗ |x(n, q)⟩` (`dualEquiv`); `U† X U` is
`X.submatrix dualEquiv dualEquiv`. In the dual frame `c_j = ann j ⊗ 1`, `τˣ_j = diag(sgn q_j)` and `τᶻ_j`
flips the charge `q_j` (`τᶻ_N := 1`). `flipSet T` flips the labels in `T`.

* (c) **Dictionary**: `U† f_j U = c_j τᶻ_j` (`ann_dual`), `U† σᶻ_{b_j} U = τᶻ_j τᶻ_{j+1}` (`sz_dual`),
  `U† ∏_{k ≥ j} σᶻ_{b_k} U = τᶻ_j` (`szTail_dual`), `U† σˣ_{b_j} U = ∏_{k ≤ j} τˣ_k (-1)^{n_k}`
  (`sx_dual`), `U† σˣ_{b_{j-1}} σˣ_{b_j} U = τˣ_j (-1)^{n_j}` (`elec_dual`),
  `U† σᶻ_{b_j} f_j† f_{j+1} U = c_j† c_{j+1}` (`hopTerm_dual`), `U† (-1)^{N_f} U = (-1)^{N_c}`
  (`parity_dual`).
* (d) **Dressed fermions** `c̃_j = U c_j U† = f_j ∏_{k ≥ j} σᶻ_{b_k}` (`dressed`, `dressed_eq`): they
  satisfy the CAR (`dressed_anticomm`, `dressed_CAR`), commute with every `G_k` (`dressed_comm_gauss`),
  are odd (`dressed_odd`), `c̃_j† c̃_j = n_j` (`dressed_number`), `c̃_x† c̃_y = f_x† (∏_{x ≤ k < y} σᶻ_{b_k}) f_y`
  for `x < y` (`dressed_hop`), and `H = -J ∑_j (c̃_j† c̃_{j+1} + h.c.) - h ∑_j G_j (-1)^{n_j}`
  (`H_dressed`).
* (e) **Dynamics**: `U† e^{-itH} U = ⊕_q e^{iht ∑_j q_j} Γ(e^{-ith_q})` (`evo_dual`) and the Heisenberg
  evolution of the dressed fermions `e^{itH} c̃_x e^{-itH} = ∑_q P_q ∑_y (e^{-ith_q})_{xy} c̃_y`
  (`heisenberg_dressed`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The flip of the labels in `T`: `|y⟩ ↦ |y + 1_T⟩`. -/
def flipSet (T : Finset (Fin N)) : Matrix (Cfg N) (Cfg N) ℝ := fun y' y =>
  if y' = y + fockCfg T then 1 else 0

/-- The fermion parity `(-1)^{N_f}` on the matter Fock space. -/
def parityF : Matrix (Cfg N) (Cfg N) ℝ := diagonal (fun n => sgn (∑ k, n k))

/-- The sites `k ≥ j`. -/
def tailOf (j : Fin N) : Finset (Fin N) := univ.filter (fun k => j ≤ k)

/-- The sites `x ≤ k < y`. -/
def between (x y : Fin N) : Finset (Fin N) := univ.filter (fun k => x ≤ k ∧ k < y)

/-! ### Auxiliary lemmas -/

namespace DualityFullAux

/-- `a + a = 0` in `ZMod 2`. -/
lemma zmod2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `v + v = 0` for configurations. -/
lemma cfg_add_self (v : Cfg N) : v + v = 0 := by
  funext k; exact zmod2_add_self (v k)

/-- `y' = y + v ↔ y = y' + v` for configurations (characteristic `2`). -/
lemma cfg_eq_add_comm (y' y v : Cfg N) : y' = y + v ↔ y = y' + v := by
  constructor <;> rintro rfl <;> rw [add_assoc, cfg_add_self, add_zero]

/-- `1_T + 1_{T'} = 1_{T ∆ T'}`. -/
lemma fockCfg_add (T T' : Finset (Fin N)) :
    fockCfg T + fockCfg T' = fockCfg (symmDiff T T') := by
  funext k
  simp only [Pi.add_apply, fockCfg, Finset.mem_symmDiff]
  by_cases h : k ∈ T <;> by_cases h' : k ∈ T' <;> simp [h, h', GaussLawAux.zmod2_one_add_one]

/-- `1_∅ = 0`. -/
lemma fockCfg_empty : fockCfg (∅ : Finset (Fin N)) = 0 := by
  funext k; simp [fockCfg]

/-- Flipping one label: `Function.update x j (x j + 1) = x + 1_{{j}}`. -/
lemma update_add_one_eq (x : Cfg N) (j : Fin N) :
    Function.update x j (x j + 1) = x + fockCfg {j} := by
  funext k
  by_cases hk : k = j
  · subst hk; simp [fockCfg]
  · simp [fockCfg, hk]

/-- Emptying an occupied site: `Function.update n j 0 = n + 1_{{j}}` if `n j = 1`. -/
lemma update_zero_eq {n : Cfg N} {j : Fin N} (hn : n j = 1) :
    Function.update n j 0 = n + fockCfg {j} := by
  funext k
  by_cases hk : k = j
  · subst hk; simp [fockCfg, hn, GaussLawAux.zmod2_one_add_one]
  · simp [fockCfg, hk]

/-- The product of two flips is the flip of the symmetric difference. -/
lemma flipSet_mul (T T' : Finset (Fin N)) : flipSet T * flipSet T' = flipSet (symmDiff T T') := by
  ext y' y
  simp only [flipSet, Matrix.mul_apply, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq',
    Finset.mem_univ, ite_true]
  rw [add_assoc, fockCfg_add, symmDiff_comm]

/-- Flips are symmetric. -/
lemma flipSet_transpose (T : Finset (Fin N)) : (flipSet T)ᵀ = flipSet T := by
  ext y' y
  simp only [transpose_apply, flipSet]
  by_cases h : y' = y + fockCfg T
  · rw [ite_eq_left h, ite_eq_left ((cfg_eq_add_comm _ _ _).1 h)]
  · rw [ite_eq_right h, ite_eq_right (fun h' => h ((cfg_eq_add_comm _ _ _).2 h'))]

/-- The empty flip is the identity. -/
lemma flipSet_empty : flipSet (∅ : Finset (Fin N)) = 1 := by
  ext y' y
  simp [flipSet, fockCfg_empty, Matrix.one_apply]

/-- Flips are involutions. -/
lemma flipSet_mul_self (T : Finset (Fin N)) : flipSet T * flipSet T = 1 := by
  rw [flipSet_mul, symmDiff_self, Finset.bot_eq_empty, flipSet_empty]

/-- `{k ≥ x} ∆ {k ≥ y} = {x ≤ k < y}` for `x < y`. -/
lemma tailOf_symmDiff_of_lt {x y : Fin N} (hxy : x < y) :
    symmDiff (tailOf x) (tailOf y) = between x y := by
  ext k
  simp only [tailOf, between, Finset.mem_symmDiff, Finset.mem_filter, Finset.mem_univ, true_and,
    Fin.le_def, Fin.lt_def]
  rw [Fin.lt_def] at hxy
  omega

/-- `{k ≥ j} ∆ {k ≥ j+1} = {j}`. -/
lemma tailOf_symmDiff_succ (j j' : Fin N) (h : j'.val = j.val + 1) :
    symmDiff (tailOf j) (tailOf j') = {j} := by
  ext k
  simp only [tailOf, Finset.mem_symmDiff, Finset.mem_filter, Finset.mem_univ, true_and,
    Finset.mem_singleton, Fin.le_def, Fin.ext_iff]
  omega

/-- The link configuration determines the charges: `x(n, q) = x ↔ q = q(n, x)`. -/
lemma linkOf_eq_iff (n q x : Cfg N) : linkOf n q = x ↔ q = chargeOf n x := by
  constructor
  · rintro rfl; rw [chargeOf_linkOf]
  · rintro rfl; rw [linkOf_chargeOf]

/-- `q ↦ x(n, q)` is injective. -/
lemma linkOf_inj' (n q' q : Cfg N) : linkOf n q' = linkOf n q ↔ q' = q := by
  rw [linkOf_eq_iff, chargeOf_linkOf]

/-- `prevLink` is additive. -/
lemma prevLink_add (x z : Cfg N) (j : Fin N) :
    prevLink (x + z) j = prevLink x j + prevLink z j := by
  unfold prevLink
  split_ifs <;> simp

/-- Changing the links by `z` changes the charges by the discrete derivative of `z`. -/
lemma chargeOf_add_right (n x z : Cfg N) :
    chargeOf n (x + z) = chargeOf n x + fun k => prevLink z k + z k := by
  funext k
  simp only [chargeOf, Pi.add_apply, prevLink_add]
  ring

/-- Changing the occupations by `m` changes the charges by `m`. -/
lemma chargeOf_add_left (n m x : Cfg N) : chargeOf (n + m) x = chargeOf n x + m := by
  funext k
  simp only [chargeOf, Pi.add_apply]
  ring

/-- The discrete derivative of `1_{{j}}` is `1_{{j, j+1}}`. -/
lemma delta_single (j : Fin N) :
    (fun k => prevLink (fockCfg {j}) k + fockCfg {j} k) =
      fockCfg (univ.filter (fun k : Fin N => k = j ∨ k.val = j.val + 1)) := by
  funext k
  simp only [prevLink, fockCfg, Finset.mem_singleton, Finset.mem_filter, Finset.mem_univ, true_and]
  split_ifs <;> simp_all [Fin.ext_iff, GaussLawAux.zmod2_one_add_one] <;> omega

/-- The discrete derivative of `1_{k ≥ j}` is `1_{{j}}`. -/
lemma delta_tail (j : Fin N) :
    (fun k => prevLink (fockCfg (tailOf j)) k + fockCfg (tailOf j) k) = fockCfg {j} := by
  funext k
  simp only [prevLink, fockCfg, tailOf, Finset.mem_singleton, Finset.mem_filter, Finset.mem_univ,
    true_and, Fin.le_def, Fin.ext_iff]
  split_ifs <;> simp_all [GaussLawAux.zmod2_one_add_one] <;> omega

/-- Support of `c_j† c_{j'}` for `j ≠ j'`: `n_j ≠ n_{j'}` and `n' = n ∘ swap j j'`. -/
lemma cre_mul_ann_ne_zero {j j' : Fin N} (hjj' : j ≠ j') {n' n : Cfg N}
    (h : (cre j * ann j') n' n ≠ 0) : n j ≠ n j' ∧ n' = n ∘ Equiv.swap j j' := by
  rw [FreeDynAux.mul_ann_apply, FreeDynAux.cre_apply] at h
  by_cases h1 : n j' = 1
  · rw [ite_eq_left h1] at h
    by_cases h2 : Function.update n j' 0 j = 0 ∧
        n' = Function.update (Function.update n j' 0) j 1
    · obtain ⟨h2, rfl⟩ := h2
      rw [Function.update_of_ne hjj'] at h2
      refine ⟨by rw [h1, h2]; decide, ?_⟩
      funext k
      by_cases hkj : k = j
      · subst hkj; simp [h1]
      · by_cases hkj' : k = j'
        · subst hkj'; simp [Function.update_of_ne hkj, h2]
        · simp [Function.update_of_ne hkj, Function.update_of_ne hkj',
            Equiv.swap_apply_of_ne_of_ne hkj hkj']
    · rw [ite_eq_right h2, mul_zero] at h; exact absurd rfl h
  · rw [ite_eq_right h1] at h; exact absurd rfl h

end DualityFullAux

open DualityFullAux

/-! ### (c) The dictionary -/

/-- `σᶻ_{b_j}` is the flip of the link `b_j`. -/
theorem sz_eq_flipSet (j : Fin N) : sz j = flipSet {j} := by
  ext x' x
  simp only [sz, flipSet]
  rw [update_add_one_eq]

/-- `U† f_j U = c_j τᶻ_j`. -/
theorem ann_dual (j : Fin N) :
    (ann j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      ann j ⊗ₖ flipSet {j} := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply, kronecker_apply, kronecker_apply]
  by_cases ha : n j = 1 ∧ n' = Function.update n j 0
  · obtain ⟨hn, rfl⟩ := ha
    congr 1
    simp only [Matrix.one_apply, flipSet, linkOf_eq_iff, update_zero_eq hn, chargeOf_add_left,
      chargeOf_linkOf]
  · have : ann j n' n = 0 := by rw [FreeDynAux.ann_apply, ite_eq_right ha]
    rw [this, zero_mul, zero_mul]

/-- `U† σᶻ_{b_j} U = τᶻ_j τᶻ_{j+1}` (`τᶻ_N := 1`). -/
theorem sz_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sz j).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (univ.filter (fun k => k = j ∨ k.val = j.val + 1)) := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply, kronecker_apply, kronecker_apply]
  by_cases hn : n' = n
  · subst hn
    congr 1
    rw [sz_eq_flipSet]
    simp only [flipSet, linkOf_eq_iff, chargeOf_add_right, chargeOf_linkOf, delta_single]
  · simp [hn]

/-- `U† (∏_{k ≥ j} σᶻ_{b_k}) U = τᶻ_j`. -/
theorem szTail_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (tailOf j)).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet {j} := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply, kronecker_apply, kronecker_apply]
  by_cases hn : n' = n
  · subst hn
    congr 1
    simp only [flipSet, linkOf_eq_iff, chargeOf_add_right, chargeOf_linkOf, delta_tail]
  · simp [hn]

/-- `U† σˣ_{b_j} U = ∏_{k ≤ j} τˣ_k (-1)^{n_k}`. -/
theorem sx_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sx j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun a : Cfg N × Cfg N => ∏ k ∈ univ.filter (· ≤ j), sgn (a.2 k) * sgn (a.1 k)) := by
  rw [sx, GaussLawAux.one_kronecker_diagonal, submatrix_diagonal_equiv]
  congr 1
  funext a
  simp only [Function.comp_apply, dualEquiv_apply, linkOf]
  rw [GaussLawAux.sgn_sum]
  exact Finset.prod_congr rfl fun k _ => GaussLawAux.sgn_add _ _

/-- `U† σˣ_{b_{j-1}} σˣ_{b_j} U = τˣ_j (-1)^{n_j}`. -/
theorem elec_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elec j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun a : Cfg N × Cfg N => sgn (a.2 j) * sgn (a.1 j)) := by
  rw [elec, GaussLawAux.one_kronecker_diagonal, submatrix_diagonal_equiv]
  congr 1
  funext a
  simp only [Function.comp_apply, dualEquiv_apply]
  rw [prevLink_linkOf_add, GaussLawAux.sgn_add]

/-- `U† (σᶻ_{b_j} f_j† f_{j+1}) U = c_j† c_{j+1}` (`j < N - 1`). -/
theorem hopTerm_dual (j : Fin N) (hj : j.val + 1 < N) :
    ((cre j * ann ⟨j.val + 1, hj⟩) ⊗ₖ sz j).submatrix (dualEquiv (N := N)) dualEquiv =
      (cre j * ann ⟨j.val + 1, hj⟩) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  have hjj' : j ≠ ⟨j.val + 1, hj⟩ := by
    intro h; have := congrArg Fin.val h; simp at this
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply, kronecker_apply, kronecker_apply]
  by_cases h0 : (cre j * ann ⟨j.val + 1, hj⟩) n' n = 0
  · rw [h0, zero_mul, zero_mul]
  · obtain ⟨hne, rfl⟩ := cre_mul_ann_ne_zero hjj' h0
    congr 1
    simp only [sz, Matrix.one_apply]
    rw [← linkOf_swap n q j hj hne]
    simp only [linkOf_inj']

/-- `U† (-1)^{N_f} U = (-1)^{N_c}`: `U` is even. -/
theorem parity_dual :
    (parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  rw [parityF, ← diagonal_one, diagonal_kronecker_diagonal, submatrix_diagonal_equiv]
  rfl

/-! ### (d) Dressed fermions -/

/-- The dressed fermion `c̃_j = f_j ∏_{k ≥ j} σᶻ_{b_k}` (Wilson line to the boundary). -/
def dressed (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ := ann j ⊗ₖ flipSet (tailOf j)

/-- **`c̃_j = U c_j U†`** (paper: eq. `eq:dressed`). -/
theorem dressed_eq (j : Fin N) :
    dressed j =
      (ann j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)).symm dualEquiv.symm := by
  sorry

/-- CAR: `{c̃_x, c̃_y} = 0`. -/
theorem dressed_anticomm (x y : Fin N) : dressed x * dressed y + dressed y * dressed x = 0 := by
  sorry

/-- CAR: `{c̃_x, c̃_y†} = δ_{xy}`. -/
theorem dressed_CAR (x y : Fin N) :
    dressed x * (dressed y)ᵀ + (dressed y)ᵀ * dressed x = if x = y then 1 else 0 := by
  sorry

/-- The dressed fermions are gauge invariant. -/
theorem dressed_comm_gauss (j k : Fin N) : dressed j * gauss k = gauss k * dressed j := by
  sorry

/-- The dressed fermions are odd. -/
theorem dressed_odd (j : Fin N) :
    dressed j * (parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) =
      -((parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) * dressed j) := by
  sorry

/-- `c̃_j† c̃_j = n_j`. -/
theorem dressed_number (j : Fin N) :
    (dressed j)ᵀ * dressed j = numOp {j} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

/-- **Gauge-invariant hopping** (paper: Theorem `thm:duality`(d)): for `x < y`,
`c̃_x† c̃_y = f_x† (∏_{x ≤ k < y} σᶻ_{b_k}) f_y`. -/
theorem dressed_hop (x y : Fin N) (hxy : x < y) :
    (dressed x)ᵀ * dressed y =
      (cre x ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
        ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (between x y)) *
          (ann y ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  sorry

/-- **`H` in terms of dressed fermions** (paper: eq. `eq:Hdressed`):
`H = -J ∑_j (c̃_j† c̃_{j+1} + c̃_{j+1}† c̃_j) - h ∑_j G_j (-1)^{n_j}`. -/
theorem H_dressed (J h : ℝ) :
    H J h =
      (-J) • ∑ j : Fin N, (if hj : j.val + 1 < N then
          (dressed j)ᵀ * dressed ⟨j.val + 1, hj⟩ + (dressed ⟨j.val + 1, hj⟩)ᵀ * dressed j else 0) -
        h • ∑ j : Fin N, gauss j * (diagonal (fun n : Cfg N => sgn (n j)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  sorry

/-! ### (e) Dynamics -/

/-- **The dynamics in the dual frame** (paper: Theorem `thm:duality`(e)):
`U† e^{-itH} U = ⊕_q e^{iht ∑_j q_j} Γ(e^{-ith_q})`. -/
theorem evo_dual (J h : ℝ) (t : ℝ) :
    (evo (H J h) t).submatrix (dualEquiv (N := N)) dualEquiv =
      blockDiagonal (fun q : Cfg N =>
        Complex.exp (Complex.I * (h * t * ∑ j, sgn (q j) : ℝ)) • evo (dGamma (hq J h q)) t) := by
  sorry

/-- **Heisenberg evolution of the dressed fermions** (paper: eq. `eq:heisenberg`):
`e^{itH} c̃_x e^{-itH} = ∑_q P_q ∑_y (e^{-ith_q})_{xy} c̃_y`. -/
theorem heisenberg_dressed (J h : ℝ) (x : Fin N) (t : ℝ) :
    (evo (H J h) t)ᴴ * cplx (dressed x) * evo (H J h) t =
      ∑ q : Cfg N, cplx (gaussProj q) * ∑ y : Fin N, evo (hq J h q) t x y • cplx (dressed y) := by
  sorry

end DFL
