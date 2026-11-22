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

namespace FreeDynAux

/-- `N_R = ∑_{x ∈ R} c_x† c_x`. -/
theorem numOp_eq_sum (R : Finset (Fin N)) : numOp R = ∑ x ∈ R, cre x * ann x := by
  simp only [cre_mul_ann_self, numOp]
  ext n' n
  rw [Matrix.sum_apply]
  by_cases h : n' = n
  · subst h; simp [count]
  · simp [h]

/-- Column `s` of `c_a† c_b`: it is supported on the single configuration obtained from `s` by
emptying `b` and filling `a`. -/
theorem cre_mul_ann_col (a b : Fin N) (m s : Cfg N) :
    (cre a * ann b) m s = if m = Function.update (Function.update s b 0) a 1 then
      (if s b = 1 ∧ Function.update s b 0 a = 0 then
        sgn (jwPar s b) * sgn (jwPar (Function.update s b 0) a) else 0) else 0 := by
  rw [mul_ann_apply, cre_apply]
  by_cases h1 : s b = 1 <;> by_cases h2 : Function.update s b 0 a = 0 <;>
    by_cases h3 : m = Function.update (Function.update s b 0) a 1 <;> simp [h1, h2, h3]

/-- Two "hopping" configurations `s - e_z + e_w` and `s - e_v + e_u` coincide only if
`(z, w) = (v, u)`. -/
theorem upd_upd_eq_iff {s : Cfg N} {z w u v : Fin N} (hzw : z ≠ w) (huv : u ≠ v)
    (hz : s z = 1) (hw : s w = 0) (hv : s v = 1) (hu : s u = 0) :
    Function.update (Function.update s z 0) w 1 = Function.update (Function.update s v 0) u 1 ↔
      z = v ∧ w = u := by
  constructor
  · intro h
    have hzv : z = v := by
      by_contra hzv
      have h' := congrFun h z
      rw [Function.update_of_ne hzw, Function.update_self] at h'
      by_cases hzu : z = u
      · subst hzu
        rw [Function.update_self] at h'
        exact zero_ne_one h'
      · rw [Function.update_of_ne hzu, Function.update_of_ne hzv, hz] at h'
        exact zero_ne_one h'
    subst hzv
    refine ⟨rfl, ?_⟩
    by_contra hwu
    have h' := congrFun h w
    rw [Function.update_self, Function.update_of_ne hwu, Function.update_of_ne hzw.symm, hw] at h'
    exact one_ne_zero h'
  · rintro ⟨rfl, rfl⟩
    rfl

/-- Diagonal matrix elements of `c_z† c_w c_u† c_v` (Wick's theorem for a basis state). -/
theorem wick4_entry (z w u v : Fin N) (s : Cfg N) :
    (cre z * ann w * (cre u * ann v)) s s =
      (if z = w ∧ s z = 1 then 1 else 0) * (if u = v ∧ s u = 1 then 1 else 0) +
      (if z = v ∧ w = u ∧ s z = 1 ∧ s w = 0 then 1 else 0) := by
  have hT : ∀ m, (cre z * ann w) s m = (cre w * ann z) m s := by
    intro m
    have h : (cre z * ann w)ᵀ = cre w * ann z := by
      simp [Matrix.transpose_mul, cre]
    rw [← h, transpose_apply]
  rw [Matrix.mul_apply]
  simp only [hT, cre_mul_ann_col, ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ,
    ite_true]
  rcases eq_or_ne z w with rfl | hzw
  · simp only [Function.update_self, Function.update_idem, jwPar_update_self, sgn_mul_self,
      and_true]
    rcases eq_or_ne u v with rfl | huv
    · simp only [Function.update_self, Function.update_idem, jwPar_update_self, sgn_mul_self,
        and_true]
      by_cases hz : s z = 1
      · have h1 : Function.update s z 1 = s := by rw [← hz]; exact Function.update_eq_self z s
        by_cases hu : s u = 1
        · have h2 : Function.update s u 1 = s := by rw [← hu]; exact Function.update_eq_self u s
          simp [hz, hu, h1, h2]
        · simp [hz, hu]
      · simp [hz]
    · rw [Function.update_of_ne huv]
      by_cases hz : s z = 1
      · have h1 : Function.update s z 1 = s := by rw [← hz]; exact Function.update_eq_self z s
        have hne' : ¬ (z = v ∧ z = u) := fun h => huv (h.2.symm.trans h.1)
        by_cases hc : s v = 1 ∧ s u = 0
        · have hne : s ≠ Function.update (Function.update s v 0) u 1 := by
            intro h
            have h' := congrFun h u
            rw [Function.update_self, hc.2] at h'
            exact zero_ne_one h'
          simp [hz, h1, hne, huv]
        · simp [hz, h1, hc, huv]
      · simp [hz]
  · rw [Function.update_of_ne (Ne.symm hzw)]
    rcases eq_or_ne u v with rfl | huv
    · simp only [Function.update_self, Function.update_idem, jwPar_update_self, sgn_mul_self,
        and_true]
      have hne' : ¬ (z = u ∧ w = u) := fun h => hzw (h.1.trans h.2.symm)
      by_cases hc : s z = 1 ∧ s w = 0
      · by_cases hu : s u = 1
        · have h2 : Function.update s u 1 = s := by rw [← hu]; exact Function.update_eq_self u s
          have hne : Function.update (Function.update s z 0) w 1 ≠ s := by
            intro h
            have h' := congrFun h w
            rw [Function.update_self, hc.2] at h'
            exact one_ne_zero h'
          simp [hc, h2, hne, hzw, hne']
        · simp [hc, hu, hzw, hne']
      · simp [hc, hzw]
    · rw [Function.update_of_ne huv]
      by_cases hc : s z = 1 ∧ s w = 0 ∧ s v = 1 ∧ s u = 0
      · obtain ⟨hz, hw, hv, hu⟩ := hc
        simp only [upd_upd_eq_iff hzw huv hz hw hv hu]
        by_cases he : z = v ∧ w = u
        · obtain ⟨rfl, rfl⟩ := he
          simp only [hz, hw, and_self, ite_true, hzw, false_and, zero_add, ite_false]
          rw [mul_mul_mul_comm, sgn_mul_self, sgn_mul_self, mul_one]
        · simp [hz, hw, hzw, he]
      · have hR : ¬ (z = v ∧ w = u ∧ s z = 1 ∧ s w = 0) := by
          rintro ⟨rfl, rfl, hz, hw⟩
          exact hc ⟨hz, hw, hz, hw⟩
        by_cases h1 : s z = 1 ∧ s w = 0
        · have h2 : ¬ (s v = 1 ∧ s u = 0) := fun h2 => hc ⟨h1.1, h1.2, h2.1, h2.2⟩
          have hR' : ¬ (z = v ∧ w = u) := fun h => hR ⟨h.1, h.2, h1.1, h1.2⟩
          simp [h1, h2, hzw]
          exact fun a b => hR' ⟨a, b⟩
        · simp [h1, hzw]

/-- Expectation of a quadratic form `∑ K_{zw} c_z† c_w` in a basis state. -/
theorem wick2_sum (K : Matrix (Fin N) (Fin N) ℂ) (s : Cfg N) :
    (∑ z, ∑ w, K z w • cplx (cre z * ann w)) s s = ∑ z, if s z = 1 then K z z else 0 := by
  simp only [Matrix.sum_apply, Matrix.smul_apply, cplx_apply, cre_mul_ann_diag, smul_eq_mul,
    apply_ite (fun r : ℝ => (r : ℂ)), Complex.ofReal_one, Complex.ofReal_zero, ite_and, mul_ite,
    mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]

/-- Expectation of a product of two quadratic forms in a basis state (Wick's theorem). -/
theorem wick4_sum (K L : Matrix (Fin N) (Fin N) ℂ) (s : Cfg N) :
    ((∑ z, ∑ w, K z w • cplx (cre z * ann w)) * (∑ u, ∑ v, L u v • cplx (cre u * ann v))) s s =
      (∑ z, if s z = 1 then K z z else 0) * (∑ u, if s u = 1 then L u u else 0) +
      ∑ z, ∑ w, if s z = 1 ∧ s w = 0 then K z w * L w z else 0 := by
  have expand : (∑ z, ∑ w, K z w • cplx (cre z * ann w)) *
      (∑ u, ∑ v, L u v • cplx (cre u * ann v)) =
      ∑ z, ∑ w, ∑ u, ∑ v, (K z w * L u v) • cplx (cre z * ann w * (cre u * ann v)) := by
    simp only [Finset.sum_mul]
    simp only [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun z _ => Finset.sum_congr rfl (fun w _ =>
      Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))))
    rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, ← cplx_mul]
  rw [expand]
  simp only [Matrix.sum_apply, Matrix.smul_apply, cplx_apply, wick4_entry, smul_eq_mul,
    Complex.ofReal_add, Complex.ofReal_mul, apply_ite (fun r : ℝ => (r : ℂ)), Complex.ofReal_one,
    Complex.ofReal_zero, mul_add, Finset.sum_add_distrib]
  congr 1
  · have hK : ∀ M : Matrix (Fin N) (Fin N) ℂ,
        (∑ z, ∑ w, M z w * (if z = w ∧ s z = 1 then (1 : ℂ) else 0)) =
          ∑ z, if s z = 1 then M z z else 0 := by
      intro M
      simp only [ite_and, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ,
        ite_true]
    rw [← hK K, ← hK L]
    simp only [Finset.sum_mul]
    simp only [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun z _ => Finset.sum_congr rfl (fun w _ =>
      Finset.sum_congr rfl (fun u _ => Finset.sum_congr rfl (fun v _ => ?_))))
    ring
  · simp only [ite_and, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]

/-- Matrix entries of `U† P_R U`. -/
theorem conj_proj_apply (U : Matrix (Fin N) (Fin N) ℂ) (R : Finset (Fin N)) (z w : Fin N) :
    (Uᴴ * diagonal (fun x => if x ∈ R then (1 : ℂ) else 0) * U) z w =
      ∑ x ∈ R, star (U x z) * U x w := by
  rw [Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, conjTranspose_apply, mul_ite, mul_one, mul_zero, ite_mul,
    zero_mul]
  rw [Finset.sum_ite_mem, Finset.univ_inter]

/-- The Heisenberg-evolved number operator of a region is the quadratic form of `U† P_R U`. -/
theorem conj_numOp (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (R : Finset (Fin N)) (t : ℝ) :
    (evo (dGamma A) t)ᴴ * cplx (numOp R) * evo (dGamma A) t =
      ∑ z, ∑ w, ((evo A t)ᴴ * diagonal (fun x => if x ∈ R then (1 : ℂ) else 0) * evo A t) z w •
        cplx (cre z * ann w) := by
  have hW := evo_mul_conjTranspose (dGamma A) (dGamma_isSymm A hA) t
  have h1 : ∀ x, (evo (dGamma A) t)ᴴ * cplx (cre x * ann x) * evo (dGamma A) t =
      ∑ z, ∑ w, (star (evo A t x z) * evo A t x w) • cplx (cre z * ann w) := by
    intro x
    rw [cplx_mul, conj_mul _ _ _ hW, heisenberg_cre A hA, heisenberg_ann A hA]
    simp only [Finset.sum_mul]
    simp only [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun z _ => Finset.sum_congr rfl (fun w _ => ?_))
    rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul, ← cplx_mul]
  rw [numOp_eq_sum, cplx_sum, Finset.mul_sum, Finset.sum_mul]
  simp only [h1, conj_proj_apply, Finset.sum_smul]
  rw [Finset.sum_comm (s := R)]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [Finset.sum_comm (s := R)]

/-- The projection `P_S = diag(1_S)` is idempotent. -/
theorem proj_mul_self (S : Finset (Fin N)) :
    diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) *
      diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) =
      diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) := by
  rw [diagonal_mul_diagonal]
  congr 1
  ext z
  by_cases h : z ∈ S <;> simp [h]

/-- The correlation matrix is Hermitian. -/
theorem corrMat_conjTranspose (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) :
    (corrMat A S t)ᴴ = corrMat A S t := by
  have hd : (diagonal (fun z => if z ∈ S then (1 : ℂ) else 0))ᴴ =
      diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) := by
    rw [diagonal_conjTranspose]
    congr 1
    ext z
    by_cases h : z ∈ S <;> simp [h]
  simp only [corrMat, conjTranspose_mul, conjTranspose_conjTranspose, hd, Matrix.mul_assoc]

/-- A "boundary" double sum as a trace: `∑_{z ∈ S, w ∉ S} M_{zw} M_{wz} = tr(P M (1 - P) M)`. -/
theorem sum_eq_trace (M : Matrix (Fin N) (Fin N) ℂ) (S : Finset (Fin N)) :
    (∑ z, ∑ w, if z ∈ S ∧ w ∉ S then M z w * M w z else 0) =
      trace (diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) * M *
        (1 - diagonal (fun z => if z ∈ S then (1 : ℂ) else 0)) * M) := by
  have h1 : (1 : Matrix (Fin N) (Fin N) ℂ) - diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) =
      diagonal (fun z => if z ∈ S then 0 else 1) := by
    rw [← diagonal_one, diagonal_sub]
    congr 1
    ext z
    by_cases h : z ∈ S <;> simp [h]
  rw [h1, Matrix.trace]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [Matrix.diag_apply, Matrix.mul_apply]
  refine Finset.sum_congr rfl (fun w _ => ?_)
  simp only [Matrix.mul_diagonal, Matrix.diagonal_mul]
  by_cases hz : z ∈ S <;> by_cases hw : w ∈ S <;> simp [hz, hw]

/-- Rewriting a double sum over `R × Rᶜ` as a sum with an indicator. -/
theorem sum_sdiff_eq_ite (R : Finset (Fin N)) (g : Fin N → Fin N → ℂ) :
    (∑ x ∈ R, ∑ y ∈ univ \ R, g x y) = ∑ x, ∑ y, if x ∈ R ∧ y ∉ R then g x y else 0 := by
  simp only [ite_and, Finset.sum_ite_irrel, Finset.sum_const_zero]
  rw [Finset.sum_ite_mem, Finset.univ_inter]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Finset.sdiff_eq_filter, Finset.sum_filter]

/-- The trace identity behind the number variance:
`tr(P K (1-P) K) = tr(Q C (1-Q) C)` with `K = U† Q U`, `C = U P U†`, `U` unitary and `P`, `Q`
idempotent. -/
theorem trace_identity (U P Q : Matrix (Fin N) (Fin N) ℂ) (hU1 : U * Uᴴ = 1) (hU2 : Uᴴ * U = 1)
    (hP : P * P = P) (hQ : Q * Q = Q) :
    trace (P * (Uᴴ * Q * U) * (1 - P) * (Uᴴ * Q * U)) =
      trace (Q * (U * P * Uᴴ) * (1 - Q) * (U * P * Uᴴ)) := by
  have hU1' : ∀ M : Matrix (Fin N) (Fin N) ℂ, U * (Uᴴ * M) = M := fun M => by
    rw [← Matrix.mul_assoc, hU1, Matrix.one_mul]
  have hU2' : ∀ M : Matrix (Fin N) (Fin N) ℂ, Uᴴ * (U * M) = M := fun M => by
    rw [← Matrix.mul_assoc, hU2, Matrix.one_mul]
  have hP' : ∀ M : Matrix (Fin N) (Fin N) ℂ, P * (P * M) = P * M := fun M => by
    rw [← Matrix.mul_assoc, hP]
  have hQ' : ∀ M : Matrix (Fin N) (Fin N) ℂ, Q * (Q * M) = Q * M := fun M => by
    rw [← Matrix.mul_assoc, hQ]
  have e1 : trace (P * (Uᴴ * (Q * U))) = trace (Q * (U * (P * Uᴴ))) := by
    rw [← Matrix.mul_assoc P, Matrix.trace_mul_comm (P * Uᴴ), Matrix.mul_assoc]
  have e2 : trace (P * (Uᴴ * (Q * (U * (P * (Uᴴ * (Q * U))))))) =
      trace (Q * (U * (P * (Uᴴ * (Q * (U * (P * Uᴴ))))))) := by
    rw [← Matrix.mul_assoc P, Matrix.trace_mul_comm (P * Uᴴ)]
    simp only [Matrix.mul_assoc]
  simp only [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, Matrix.mul_assoc, hU1', hU2', hP',
    hQ', Matrix.trace_sub]
  rw [e1, e2]

end FreeDynAux

/-- The variance of the particle number of a region `R` in an evolved Fock state:
`⟨N_R²⟩ - ⟨N_R⟩² = ∑_{x ∈ R} ∑_{y ∉ R} |C(t)_{xy}|²`. -/
theorem number_variance (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    expect (cplx (numOp R * numOp R)) (evo (dGamma A) t *ᵥ fockVec S) -
        expect (cplx (numOp R)) (evo (dGamma A) t *ᵥ fockVec S) ^ 2 =
      ((∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y) : ℝ) : ℂ) := by
  have hW := FreeDynAux.evo_mul_conjTranspose (dGamma A) (dGamma_isSymm A hA) t
  have hU1 := FreeDynAux.evo_mul_conjTranspose A hA t
  have hU2 := FreeDynAux.conjTranspose_mul_evo A hA t
  rw [FreeDynAux.expect_mulVec, FreeDynAux.expect_mulVec, fockVec, FreeDynAux.expect_single,
    FreeDynAux.expect_single, FreeDynAux.cplx_mul, FreeDynAux.conj_mul _ _ _ hW,
    FreeDynAux.conj_numOp A hA R t, FreeDynAux.wick4_sum, FreeDynAux.wick2_sum, sq,
    add_sub_cancel_left]
  simp only [FreeDynAux.fockCfg_eq_one, FreeDynAux.fockCfg_eq_zero]
  have hnorm : ∀ x y, (Complex.normSq (corrMat A S t x y) : ℂ) =
      corrMat A S t x y * corrMat A S t y x := by
    intro x y
    rw [← Complex.mul_conj, ← Complex.star_def]
    congr 1
    rw [← conjTranspose_apply, FreeDynAux.corrMat_conjTranspose]
  push_cast
  simp only [hnorm]
  rw [FreeDynAux.sum_sdiff_eq_ite, FreeDynAux.sum_eq_trace, FreeDynAux.sum_eq_trace]
  exact FreeDynAux.trace_identity _ _ _ hU1 hU2 (FreeDynAux.proj_mul_self S)
    (FreeDynAux.proj_mul_self R)

/-! ### The distribution of the particle number -/

/-- Distribution of the particle number `N_R` in the state `v`: `p(k) = ∑_{n : N_R(n) = k} |v(n)|²`. -/
def numDist (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) : ℝ :=
  ∑ n ∈ univ.filter (fun n => count R n = k), Complex.normSq (v n)

/-- The distribution `numDist` is the spectral measure of `N_R`: its moments are expectations. -/
theorem numDist_eq_expect (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) :
    ((numDist R v k : ℝ) : ℂ) = expect (cplx (numProj R k)) v := by
  unfold numDist expect numProj
  rw [FreeDynAux.cplx_diagonal]
  simp only [Matrix.mulVec_diagonal, dotProduct, Pi.star_apply]
  push_cast
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  split_ifs with h
  · rw [Complex.normSq_eq_conj_mul_self, Complex.star_def, Complex.ofReal_one, one_mul]
  · simp

end DFL
