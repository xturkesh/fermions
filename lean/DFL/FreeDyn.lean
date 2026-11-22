import DFL.Basic
import DFL.Evolution

/-!
# Free-fermion dynamics

Canonical anticommutation relations, the Heisenberg evolution of the annihilation operators under
a quadratic Hamiltonian `dΓ(A)`, and the density and number-variance formulas for Fock initial
states (Wick's theorem for Fock states).

Main results:
* `DFL.ann_mul_ann_add`, `DFL.ann_mul_cre_add`: CAR;
* `DFL.dGamma_comm_ann`: `[dΓ(A), c_x] = -∑_y A_{xy} c_y`;
* `DFL.heisenberg_ann`: `W(t)† c_x W(t) = ∑_y (e^{-itA})_{xy} c_y`, with `W(t) = e^{-itdΓ(A)}`;
* `DFL.two_point`: `⟨c_x† c_y⟩_t = C(t)_{yx}` with `C(t) = U(t) P_S U(t)†`;
* `DFL.density`: `⟨n_x⟩_t = ∑_{z ∈ S} |U(t)_{xz}|²`;
* `DFL.number_variance`: `Var N_R = ∑_{x ∈ R, y ∉ R} |C(t)_{xy}|²`.

All auxiliary lemmas live in the namespace `DFL.FreeDynAux`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

namespace FreeDynAux

/-! ### Signs, Jordan–Wigner parities and matrix elements -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
theorem zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `sgn` is multiplicative. -/
theorem sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  fin_cases a <;> fin_cases b <;> simp [sgn] <;> rfl

/-- `sgn a ^ 2 = 1`. -/
theorem sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases zmod2_cases a with h | h <;> simp [sgn, h]

/-- `sgn 0 = 1`. -/
@[simp] theorem sgn_zero : sgn 0 = 1 := by simp [sgn]

/-- `sgn 1 = -1`. -/
@[simp] theorem sgn_one : sgn 1 = -1 := by simp [sgn]

/-- `sgn (-1) = -1`. -/
@[simp] theorem sgn_neg_one : sgn (-1) = -1 := by
  have : (-1 : ZMod 2) = 1 := by decide
  rw [this, sgn_one]

/-- Jordan–Wigner parity of an updated configuration. -/
theorem jwPar_update (n : Cfg N) (y : Fin N) (a : ZMod 2) (x : Fin N) :
    jwPar (Function.update n y a) x = jwPar n x + if y < x then a - n y else 0 := by
  have h : Function.update n y a = n + Pi.single y (a - n y) := by
    ext k
    by_cases hk : k = y
    · subst hk; simp
    · simp [hk]
  unfold jwPar
  rw [h]
  simp only [Pi.add_apply, Finset.sum_add_distrib, Finset.sum_pi_single']
  simp [Finset.mem_filter]

/-- Updating the mode `x` itself does not change the parity to the left of `x`. -/
theorem jwPar_update_self (n : Cfg N) (x : Fin N) (a : ZMod 2) :
    jwPar (Function.update n x a) x = jwPar n x := by
  simp [jwPar_update]

/-- Matrix elements of the annihilation operator. -/
theorem ann_apply (x : Fin N) (n' n : Cfg N) :
    ann x n' n = if n x = 1 ∧ n' = Function.update n x 0 then sgn (jwPar n x) else 0 := rfl

/-- Matrix elements of the creation operator. -/
theorem cre_apply (x : Fin N) (n' n : Cfg N) :
    cre x n' n = if n x = 0 ∧ n' = Function.update n x 1 then sgn (jwPar n x) else 0 := by
  simp only [cre, transpose_apply, ann_apply]
  by_cases h1 : n' x = 1 ∧ n = Function.update n' x 0
  · obtain ⟨h1, rfl⟩ := h1
    simp [jwPar_update_self, ← h1]
  · have h2 : ¬ (n x = 0 ∧ n' = Function.update n x 1) := by
      rintro ⟨h2, rfl⟩
      exact h1 ⟨by simp, by simp [← h2]⟩
    simp only [h1, h2, ite_false]

/-- Right multiplication by an annihilation operator. -/
theorem mul_ann_apply (M : Matrix (Cfg N) (Cfg N) ℝ) (y : Fin N) (n' n : Cfg N) :
    (M * ann y) n' n = if n y = 1 then sgn (jwPar n y) * M n' (Function.update n y 0) else 0 := by
  simp only [Matrix.mul_apply, ann_apply]
  split_ifs with h
  · simp [h, mul_comm]
  · simp [h]

/-- Right multiplication by a creation operator. -/
theorem mul_cre_apply (M : Matrix (Cfg N) (Cfg N) ℝ) (y : Fin N) (n' n : Cfg N) :
    (M * cre y) n' n = if n y = 0 then sgn (jwPar n y) * M n' (Function.update n y 1) else 0 := by
  simp only [Matrix.mul_apply, cre_apply]
  split_ifs with h
  · simp [h, mul_comm]
  · simp [h]

end FreeDynAux

/-! ### Canonical anticommutation relations -/

theorem ann_mul_ann_add (x y : Fin N) : ann x * ann y + ann y * ann x = 0 := by
  ext n' n
  simp only [Matrix.add_apply, FreeDynAux.mul_ann_apply, FreeDynAux.ann_apply, Matrix.zero_apply]
  by_cases hxy : x = y
  · subst hxy; simp
  · have hyx : y ≠ x := Ne.symm hxy
    rw [Function.update_of_ne hxy, Function.update_of_ne hyx, Function.update_comm hyx]
    by_cases hx : n x = 1
    · by_cases hy : n y = 1
      · by_cases hn : n' = Function.update (Function.update n x 0) y 0
        · simp only [hx, hy, hn, and_self, ite_true, FreeDynAux.jwPar_update]
          rcases lt_or_gt_of_ne hxy with h | h
          · simp [h, not_lt.mpr h.le, FreeDynAux.sgn_add]
            ring
          · simp [h, not_lt.mpr h.le, FreeDynAux.sgn_add]
            ring
        · simp [hn]
      · simp [hy]
    · simp [hx]

theorem ann_mul_cre_add (x y : Fin N) :
    ann x * cre y + cre y * ann x = if x = y then 1 else 0 := by
  ext n' n
  simp only [Matrix.add_apply, FreeDynAux.mul_cre_apply, FreeDynAux.mul_ann_apply,
    FreeDynAux.ann_apply, FreeDynAux.cre_apply]
  by_cases hxy : x = y
  · subst hxy
    simp only [ite_true, Function.update_self, Function.update_idem,
      FreeDynAux.jwPar_update_self, true_and, Matrix.one_apply]
    rcases FreeDynAux.zmod2_cases (n x) with h | h
    · have hn : Function.update n x 0 = n := by rw [← h]; exact Function.update_eq_self x n
      simp [h, hn, FreeDynAux.sgn_mul_self]
    · have hn : Function.update n x 1 = n := by rw [← h]; exact Function.update_eq_self x n
      simp [h, hn, FreeDynAux.sgn_mul_self]
  · have hyx : y ≠ x := Ne.symm hxy
    simp only [hxy, ite_false, Matrix.zero_apply]
    rw [Function.update_of_ne hxy, Function.update_of_ne hyx, Function.update_comm hyx]
    by_cases hx : n x = 1
    · by_cases hy : n y = 0
      · by_cases hn : n' = Function.update (Function.update n x 0) y 1
        · simp only [hx, hy, hn, and_self, ite_true, FreeDynAux.jwPar_update]
          rcases lt_or_gt_of_ne hxy with h | h
          · simp [h, not_lt.mpr h.le, FreeDynAux.sgn_add]
            ring
          · simp [h, not_lt.mpr h.le, FreeDynAux.sgn_add]
            ring
        · simp [hn]
      · simp [hy]
    · simp [hx]

/-- Commutator of `c_z† c_w` with `c_x`: `[c_z† c_w, c_x] = -δ_{xz} c_w`. -/
theorem FreeDynAux.comm_cre_mul_ann (z w x : Fin N) :
    cre z * ann w * ann x - ann x * (cre z * ann w) = -(if x = z then ann w else 0) := by
  have h1 : ann w * ann x = -(ann x * ann w) := eq_neg_of_add_eq_zero_left (ann_mul_ann_add w x)
  have h2 : ann x * cre z = (if x = z then 1 else 0) - cre z * ann x :=
    eq_sub_of_add_eq (ann_mul_cre_add x z)
  calc cre z * ann w * ann x - ann x * (cre z * ann w)
      = cre z * (ann w * ann x) - (ann x * cre z) * ann w := by simp only [Matrix.mul_assoc]
    _ = -(if x = z then ann w else 0) := by
        rw [h1, h2]
        split_ifs <;> noncomm_ring

/-- `[dΓ(A), c_x] = -∑_y A_{xy} c_y`. -/
theorem dGamma_comm_ann (A : Matrix (Fin N) (Fin N) ℝ) (x : Fin N) :
    dGamma A * ann x - ann x * dGamma A = -∑ y, A x y • ann y := by
  have key : ∀ z w, (A z w • (cre z * ann w)) * ann x - ann x * (A z w • (cre z * ann w)) =
      -(if x = z then A z w • ann w else 0) := by
    intro z w
    rw [Matrix.smul_mul, Matrix.mul_smul, ← smul_sub, FreeDynAux.comm_cre_mul_ann]
    split_ifs <;> simp
  simp only [dGamma, Finset.sum_mul, Finset.mul_sum, ← Finset.sum_sub_distrib, key,
    Finset.sum_neg_distrib]
  rw [Finset.sum_comm]
  simp [Finset.sum_ite_eq]

/-- `dΓ` of a symmetric matrix is symmetric. -/
theorem dGamma_isSymm (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) : (dGamma A).IsSymm := by
  unfold Matrix.IsSymm dGamma
  simp only [Matrix.transpose_sum, Matrix.transpose_smul, Matrix.transpose_mul, cre,
    Matrix.transpose_transpose]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun b _ => ?_))
  rw [hA.apply a b]

/-! ### Heisenberg evolution -/

namespace FreeDynAux

open NormedSpace

section Analysis

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

omit [DecidableEq ι] in
/-- Entries of a differentiable matrix-valued function are differentiable. -/
theorem hasDerivAt_entry (f : ℂ → Matrix ι ι ℂ) (f' : Matrix ι ι ℂ) (u : ℂ)
    (hf : HasDerivAt f f' u) (i j : ι) : HasDerivAt (fun v => f v i j) (f' i j) u := by
  have h1 := (hasDerivAt_pi (φ := fun v => (f v : ι → ι → ℂ)) (φ' := (f' : ι → ι → ℂ))).1 hf i
  exact (hasDerivAt_pi.1 h1) j

open scoped Matrix.Norms.Operator in
/-- Abstract Heisenberg evolution: if `[X, c_x] = -∑_y Y_{xy} c_y` for all `x`, then
`e^{uX} (∑_y (e^{uY})_{xy} c_y) e^{-uX} = c_x`. -/
theorem heis_aux (X : Matrix κ κ ℂ) (Y : Matrix ι ι ℂ) (c : ι → Matrix κ κ ℂ)
    (hc : ∀ x, X * c x - c x * X = -∑ y, Y x y • c y) (x : ι) (u : ℂ) :
    exp (u • X) * (∑ y, exp (u • Y) x y • c y) * exp (u • (-X)) = c x := by
  let M : ℂ → Matrix κ κ ℂ := fun v => ∑ y, exp (v • Y) x y • c y
  let M' : ℂ → Matrix κ κ ℂ := fun v => ∑ y, (exp (v • Y) * Y) x y • c y
  let F : ℂ → Matrix κ κ ℂ := fun v => exp (v • X) * M v * exp (v • (-X))
  have hM : ∀ v, HasDerivAt M (M' v) v := by
    intro v
    apply HasDerivAt.fun_sum
    intro y _
    exact (hasDerivAt_entry _ _ v (hasDerivAt_exp_smul_const Y v) x y).smul_const (c y)
  have hkey : ∀ v, X * M v - M v * X + M' v = 0 := by
    intro v
    have e1 : X * M v - M v * X = ∑ y, exp (v • Y) x y • (X * c y - c y * X) := by
      simp only [M, Finset.mul_sum, Finset.sum_mul, Matrix.mul_smul, Matrix.smul_mul, smul_sub,
        Finset.sum_sub_distrib]
    rw [e1]
    simp only [M', hc, smul_neg, Finset.smul_sum, smul_smul, Finset.sum_neg_distrib,
      Matrix.mul_apply, Finset.sum_smul]
    rw [Finset.sum_comm]
    exact neg_add_cancel _
  have hF : ∀ v, HasDerivAt F 0 v := by
    intro v
    have h := ((hasDerivAt_exp_smul_const X v).mul (hM v)).mul
      (hasDerivAt_exp_smul_const' (-X) v)
    have h' : HasDerivAt F ((exp (v • X) * X * M v + exp (v • X) * M' v) * exp (v • -X) +
        exp (v • X) * M v * (-X * exp (v • -X))) v := h
    have e2 : (exp (v • X) * X * M v + exp (v • X) * M' v) * exp (v • -X) +
        exp (v • X) * M v * (-X * exp (v • -X)) =
        exp (v • X) * (X * M v - M v * X + M' v) * exp (v • -X) := by noncomm_ring
    rw [e2, hkey v, Matrix.mul_zero, Matrix.zero_mul] at h'
    exact h'
  have h0 := is_const_of_deriv_eq_zero (fun v => (hF v).differentiableAt)
    (fun v => (hF v).deriv) u 0
  simpa [F, M, Matrix.one_apply] using h0

end Analysis

section Cplx

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- Entries of `cplx A`. -/
@[simp] theorem cplx_apply (A : Matrix ι ι ℝ) (i j : ι) : cplx A i j = (A i j : ℂ) := rfl

omit [DecidableEq ι] in
/-- `cplx` is multiplicative. -/
theorem cplx_mul (A B : Matrix ι ι ℝ) : cplx (A * B) = cplx A * cplx B := by
  ext i j; simp [Matrix.mul_apply]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` commutes with subtraction. -/
theorem cplx_sub (A B : Matrix ι ι ℝ) : cplx (A - B) = cplx A - cplx B := by
  ext i j; simp

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` commutes with negation. -/
theorem cplx_neg (A : Matrix ι ι ℝ) : cplx (-A) = -cplx A := by
  ext i j; simp

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` commutes with real scalar multiplication. -/
theorem cplx_smul (r : ℝ) (A : Matrix ι ι ℝ) : cplx (r • A) = (r : ℂ) • cplx A := by
  ext i j; simp

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` commutes with finite sums. -/
theorem cplx_sum {α : Type*} (s : Finset α) (f : α → Matrix ι ι ℝ) :
    cplx (∑ a ∈ s, f a) = ∑ a ∈ s, cplx (f a) := by
  ext i j; simp [Matrix.sum_apply]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` turns transposes into conjugate transposes. -/
theorem cplx_conjTranspose (A : Matrix ι ι ℝ) : (cplx A)ᴴ = cplx Aᵀ := by
  ext i j; simp [conjTranspose_apply]

omit [Fintype ι] in
/-- `cplx` of a diagonal matrix. -/
theorem cplx_diagonal (d : ι → ℝ) : cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
  ext i j; by_cases h : i = j <;> simp [h]

omit [Fintype ι] [DecidableEq ι] in
/-- The generator of `evo` is skew-Hermitian for symmetric `A`. -/
theorem gen_conjTranspose (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) :
    ((-(t : ℂ) * Complex.I) • cplx A)ᴴ = -((-(t : ℂ) * Complex.I) • cplx A) := by
  rw [conjTranspose_smul, cplx_conjTranspose, hA.eq, ← neg_smul]
  congr 1
  simp

/-- `e^{-itA}` is unitary for symmetric `A` (`U U† = 1`). -/
theorem evo_mul_conjTranspose (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) :
    evo A t * (evo A t)ᴴ = 1 := by
  unfold evo
  rw [← Matrix.exp_conjTranspose, gen_conjTranspose A hA t,
    ← Matrix.exp_add_of_commute _ _ (Commute.refl ((-(t : ℂ) * Complex.I) • cplx A)).neg_right,
    add_neg_cancel, exp_zero]

/-- `e^{-itA}` is unitary for symmetric `A` (`U† U = 1`). -/
theorem conjTranspose_mul_evo (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) :
    (evo A t)ᴴ * evo A t = 1 := by
  unfold evo
  rw [← Matrix.exp_conjTranspose, gen_conjTranspose A hA t,
    ← Matrix.exp_add_of_commute _ _ (Commute.refl ((-(t : ℂ) * Complex.I) • cplx A)).neg_left,
    neg_add_cancel, exp_zero]

/-- Conjugating a product by a unitary. -/
theorem conj_mul (W M M' : Matrix ι ι ℂ) (hW : W * Wᴴ = 1) :
    Wᴴ * (M * M') * W = (Wᴴ * M * W) * (Wᴴ * M' * W) := by
  calc Wᴴ * (M * M') * W = Wᴴ * M * (W * Wᴴ) * M' * W := by
        rw [hW, Matrix.mul_one]; simp only [Matrix.mul_assoc]
    _ = (Wᴴ * M * W) * (Wᴴ * M' * W) := by simp only [Matrix.mul_assoc]

end Cplx

end FreeDynAux

open NormedSpace in
/-- **Heisenberg evolution of the annihilation operators** under a quadratic Hamiltonian:
`e^{itdΓ(A)} c_x e^{-itdΓ(A)} = ∑_y (e^{-itA})_{xy} c_y`. -/
theorem heisenberg_ann (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (x : Fin N) (t : ℝ) :
    (evo (dGamma A) t)ᴴ * cplx (ann x) * evo (dGamma A) t = ∑ y, evo A t x y • cplx (ann y) := by
  set X : Matrix (Cfg N) (Cfg N) ℂ := (-Complex.I) • cplx (dGamma A) with hXdef
  set Y : Matrix (Fin N) (Fin N) ℂ := (-Complex.I) • cplx A with hYdef
  have hW : evo (dGamma A) t = exp ((t : ℂ) • X) := by
    simp only [evo, X, smul_smul, mul_neg, neg_mul]
  have hWh : (evo (dGamma A) t)ᴴ = exp ((t : ℂ) • (-X)) := by
    rw [hW, ← Matrix.exp_conjTranspose]
    congr 1
    rw [conjTranspose_smul, hXdef, conjTranspose_smul, FreeDynAux.cplx_conjTranspose,
      (dGamma_isSymm A hA).eq]
    simp
  have hU : ∀ y, evo A t x y = exp ((t : ℂ) • Y) x y := by
    intro y
    simp only [evo, Y, smul_smul, mul_neg, neg_mul]
  have hc : ∀ z, X * cplx (ann z) - cplx (ann z) * X = -∑ y, Y z y • cplx (ann y) := by
    intro z
    have h := congrArg cplx (dGamma_comm_ann A z)
    rw [FreeDynAux.cplx_sub, FreeDynAux.cplx_mul, FreeDynAux.cplx_mul, FreeDynAux.cplx_neg,
      FreeDynAux.cplx_sum] at h
    simp only [FreeDynAux.cplx_smul] at h
    simp only [X, Y, Matrix.smul_mul, Matrix.mul_smul, ← smul_sub, h, smul_neg, Finset.smul_sum,
      smul_smul, Matrix.smul_apply, FreeDynAux.cplx_apply, smul_eq_mul]
  have key := FreeDynAux.heis_aux X Y (fun y => cplx (ann y)) hc x (t : ℂ)
  have inv1 : exp ((t : ℂ) • -X) * exp ((t : ℂ) • X) = 1 := by
    rw [← Matrix.exp_add_of_commute, smul_neg, neg_add_cancel, exp_zero]
    rw [smul_neg]
    exact (Commute.refl _).neg_left
  rw [hWh, hW, ← key]
  calc exp ((t : ℂ) • -X) * (exp ((t : ℂ) • X) * (∑ y, exp ((t : ℂ) • Y) x y • cplx (ann y)) *
        exp ((t : ℂ) • -X)) * exp ((t : ℂ) • X)
      = (exp ((t : ℂ) • -X) * exp ((t : ℂ) • X)) * (∑ y, exp ((t : ℂ) • Y) x y • cplx (ann y)) *
        (exp ((t : ℂ) • -X) * exp ((t : ℂ) • X)) := by simp only [Matrix.mul_assoc]
    _ = ∑ y, evo A t x y • cplx (ann y) := by
        rw [inv1, Matrix.one_mul, Matrix.mul_one]
        simp only [hU]

/-! ### Fock initial states -/

/-- The Fock basis vector with occupied set `S`. -/
def fockVec (S : Finset (Fin N)) : Cfg N → ℂ := Pi.single (fockCfg S) 1

/-- The one-body correlation matrix `C(t) = U(t) P_S U(t)†` of the evolved Fock state, where
`U(t) = e^{-itA}` and `P_S` is the projection onto `S`. -/
def corrMat (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) :
    Matrix (Fin N) (Fin N) ℂ :=
  evo A t * diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) * (evo A t)ᴴ

namespace FreeDynAux

/-- The diagonal matrix elements of `c_z† c_w`. -/
theorem cre_mul_ann_diag (z w : Fin N) (s : Cfg N) :
    (cre z * ann w) s s = if z = w ∧ s w = 1 then 1 else 0 := by
  rw [mul_ann_apply, cre_apply]
  by_cases hw : s w = 1
  · by_cases hzw : z = w
    · subst hzw
      have h1 : Function.update s z 1 = s := by rw [← hw]; exact Function.update_eq_self z s
      simp [hw, h1, jwPar_update_self, sgn_mul_self]
    · have h1 : s ≠ Function.update (Function.update s w 0) z 1 := by
        intro h
        have := congrFun h w
        rw [Function.update_of_ne (Ne.symm hzw), Function.update_self, hw] at this
        exact one_ne_zero this
      simp [hw, hzw, h1]
  · simp [hw]

/-- `c_x† c_x` is the occupation operator. -/
theorem cre_mul_ann_self (x : Fin N) : cre x * ann x = diagonal (fun n => occ (n x)) := by
  ext n' n
  rw [mul_ann_apply, cre_apply, diagonal_apply]
  by_cases hn : n' = n
  · subst hn
    rcases zmod2_cases (n' x) with h | h
    · simp [h, occ]
    · have h1 : Function.update n' x 1 = n' := by rw [← h]; exact Function.update_eq_self x n'
      simp [h, h1, occ, jwPar_update_self, sgn_mul_self]
  · rcases zmod2_cases (n x) with h | h
    · simp [h, hn]
    · have h1 : Function.update n x 1 = n := by rw [← h]; exact Function.update_eq_self x n
      simp [h, h1, hn]

/-- The Fock configuration takes the value `1` exactly on `S`. -/
theorem fockCfg_eq_one (S : Finset (Fin N)) (z : Fin N) : fockCfg S z = 1 ↔ z ∈ S := by
  by_cases h : z ∈ S <;> simp [fockCfg, h]

/-- The Fock configuration takes the value `0` exactly off `S`. -/
theorem fockCfg_eq_zero (S : Finset (Fin N)) (z : Fin N) : fockCfg S z = 0 ↔ z ∉ S := by
  by_cases h : z ∈ S <;> simp [fockCfg, h]

section Expect

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
/-- Expectations in a transformed state. -/
theorem expect_mulVec (M W : Matrix ι ι ℂ) (v : ι → ℂ) :
    expect M (W *ᵥ v) = expect (Wᴴ * M * W) v := by
  unfold expect
  rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec]
  simp only [Matrix.mulVec_mulVec, Matrix.mul_assoc]

/-- Expectation in a basis vector. -/
theorem expect_single (M : Matrix ι ι ℂ) (i : ι) : expect M (Pi.single i 1) = M i i := by
  unfold expect
  simp [single_dotProduct]

end Expect

/-- Heisenberg evolution of the creation operators. -/
theorem heisenberg_cre (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (x : Fin N) (t : ℝ) :
    (evo (dGamma A) t)ᴴ * cplx (cre x) * evo (dGamma A) t =
      ∑ y, star (evo A t x y) • cplx (cre y) := by
  have h := congrArg conjTranspose (heisenberg_ann A hA x t)
  simp only [Matrix.conjTranspose_mul, Matrix.conjTranspose_conjTranspose,
    Matrix.conjTranspose_sum, Matrix.conjTranspose_smul, cplx_conjTranspose, Matrix.mul_assoc]
    at h
  simpa only [cre, Matrix.mul_assoc] using h

/-- Expansion of the correlation matrix. -/
theorem corrMat_apply (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (i j : Fin N) :
    corrMat A S t i j = ∑ z ∈ S, evo A t i z * star (evo A t j z) := by
  rw [corrMat, Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, conjTranspose_apply, mul_ite, mul_one, mul_zero, ite_mul,
    zero_mul]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

end FreeDynAux

/-- The two-point function of an evolved Fock state: `⟨c_x† c_y⟩_t = C(t)_{yx}`. -/
theorem two_point (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x y : Fin N) (t : ℝ) :
    expect (cplx (cre x * ann y)) (evo (dGamma A) t *ᵥ fockVec S) = corrMat A S t y x := by
  have hW := FreeDynAux.evo_mul_conjTranspose (dGamma A) (dGamma_isSymm A hA) t
  rw [FreeDynAux.expect_mulVec, fockVec, FreeDynAux.expect_single, FreeDynAux.cplx_mul,
    FreeDynAux.conj_mul _ _ _ hW, FreeDynAux.heisenberg_cre A hA, heisenberg_ann A hA,
    FreeDynAux.corrMat_apply]
  simp only [Finset.sum_mul, Finset.mul_sum, Matrix.smul_mul, Matrix.mul_smul, smul_smul,
    ← FreeDynAux.cplx_mul, Matrix.sum_apply, Matrix.smul_apply, FreeDynAux.cplx_apply,
    FreeDynAux.cre_mul_ann_diag, smul_eq_mul, FreeDynAux.fockCfg_eq_one]
  simp only [apply_ite (fun r : ℝ => (r : ℂ)), Complex.ofReal_one, Complex.ofReal_zero, ite_and,
    mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

/-- The density of an evolved Fock state: `⟨n_x⟩_t = ∑_{z ∈ S} |U(t)_{xz}|²`. -/
theorem density (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (x : Fin N)
    (t : ℝ) :
    expect (cplx (diagonal fun n => occ (n x))) (evo (dGamma A) t *ᵥ fockVec S) =
      ((∑ z ∈ S, Complex.normSq (evo A t x z) : ℝ) : ℂ) := by
  rw [← FreeDynAux.cre_mul_ann_self, two_point A hA S x x t, FreeDynAux.corrMat_apply]
  push_cast
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [Complex.star_def, Complex.mul_conj]

/-- The variance of the particle number of a region `R` in an evolved Fock state:
`⟨N_R²⟩ - ⟨N_R⟩² = ∑_{x ∈ R} ∑_{y ∉ R} |C(t)_{xy}|²`. -/
theorem number_variance (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    expect (cplx (numOp R * numOp R)) (evo (dGamma A) t *ᵥ fockVec S) -
        expect (cplx (numOp R)) (evo (dGamma A) t *ᵥ fockVec S) ^ 2 =
      ((∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y) : ℝ) : ℂ) := by
  sorry

/-! ### The distribution of the particle number -/

/-- Distribution of the particle number `N_R` in the state `v`: `p(k) = ∑_{n : N_R(n) = k} |v(n)|²`. -/
def numDist (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) : ℝ :=
  ∑ n ∈ univ.filter (fun n => count R n = k), Complex.normSq (v n)

/-- The distribution `numDist` is the spectral measure of `N_R`: its moments are expectations. -/
theorem numDist_eq_expect (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) :
    ((numDist R v k : ℝ) : ℂ) = expect (cplx (numProj R k)) v := by
  sorry

end DFL
