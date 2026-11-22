import DFL.Evolution

/-!
# Finite propagation speed for one-dimensional Schrödinger operators (paper: Lemma `lem:dimerTL`)

For an interval `Λ = [a, b] ∋ 0` and a real potential `w`, let `H_Λ = -Δ_Λ + diag(w)` (Dirichlet ends).

* (i) `|⟨δ_y, e^{-isH_Λ} δ_0⟩| ≤ ∑_{k ≥ |y|} (2|s|)^k / k! ≤ (2|s|)^{|y|}/|y|! · e^{2|s|}`, uniformly in the
  potential (`dimerTL_i`, `dimerTL_i'`);
* (ii) for `Λ ⊆ Λ' = [a', b']` and `|w| ≤ W` on `Λ'`,
  `|⟨δ_y, (e^{-isH_Λ} - e^{-isH_{Λ'}}) δ_0⟩| ≤ 2 ∑_{k ≥ R} ((2 + W)|s|)^k / k!` with
  `R = dist(0, ℤ \ Λ)` (`dimerTL_ii`; here `Λ'` is finite);
* (iii) for `0 < p ≤ 2`, `∑_{y ∈ Λ} |y|^p |⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ (2|s|)^p` (`dimerTL_iii`).
-/

open Matrix

noncomputable section

namespace DFL

/-- `H_Λ = -Δ_Λ + diag(w)` on `ℓ²([a, b])`. -/
def schrH (a b : ℤ) (w : ℤ → ℝ) : Matrix (Set.Icc a b) (Set.Icc a b) ℝ := fun x y =>
  if x.1 = y.1 then w x.1 else if |x.1 - y.1| = 1 then -1 else 0

/-- The origin as a site of `[a, b]`. -/
def origin (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) : Set.Icc a b := ⟨0, ha, hb⟩

namespace DimerTLAux

open NormedSpace

section General

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The generator of `evo` as a real multiple of a fixed matrix. -/
lemma evo_eq_exp_smul (A : Matrix ι ι ℝ) (t : ℝ) :
    evo A t = exp (t • ((-Complex.I) • cplx A)) := by
  unfold evo
  congr 1
  ext i j
  simp only [Matrix.smul_apply, smul_eq_mul, Complex.real_smul]
  ring

/-- `evo A 0 = 1`. -/
lemma evo_zero (A : Matrix ι ι ℝ) : evo A 0 = 1 := by
  simp [evo]

section Deriv

open scoped Matrix.Norms.Operator

set_option backward.isDefEq.respectTransparency false

/-- Derivative of the propagator (generator on the right). -/
lemma hasDerivAt_evo (A : Matrix ι ι ℝ) (t : ℝ) :
    HasDerivAt (fun u : ℝ => evo A u) (evo A t * ((-Complex.I) • cplx A)) t := by
  simp only [evo_eq_exp_smul]
  exact hasDerivAt_exp_smul_const _ t

/-- Derivative of the propagator (generator on the left). -/
lemma hasDerivAt_evo' (A : Matrix ι ι ℝ) (t : ℝ) :
    HasDerivAt (fun u : ℝ => evo A u) (((-Complex.I) • cplx A) * evo A t) t := by
  simp only [evo_eq_exp_smul]
  exact hasDerivAt_exp_smul_const' _ t

/-- Derivative of the Heisenberg-evolved observable `r ↦ U(-r) X U(r)`. -/
lemma hasDerivAt_heis (A X : Matrix ι ι ℝ) (s : ℝ) :
    HasDerivAt (fun r : ℝ => evo A (-r) * cplx X * evo A r)
      (evo A (-s) * (Complex.I • (cplx A * cplx X - cplx X * cplx A)) * evo A s) s := by
  have h1 := (hasDerivAt_evo A (-s)).scomp s (hasDerivAt_neg s)
  have h2 := (h1.mul_const (cplx X)).mul (hasDerivAt_evo' A s)
  convert h2 using 1
  · rfl
  · rfl
  · simp only [Function.comp_apply, neg_smul, one_smul, smul_neg, Matrix.mul_neg, Matrix.neg_mul,
      Matrix.mul_smul, Matrix.smul_mul, smul_sub, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc,
      neg_neg]
    abel

end Deriv

/-- Column extraction as a linear map. -/
def colLin (z : ι) : Matrix ι ι ℂ →ₗ[ℝ] EuclideanSpace ℂ ι where
  toFun M := WithLp.toLp 2 (fun x => M x z)
  map_add' M M' := by
    rw [← WithLp.toLp_add]
    rfl
  map_smul' c M := by
    rw [RingHom.id_apply, ← WithLp.toLp_smul]
    rfl

/-- Column extraction as a continuous linear map. -/
def colL (z : ι) : Matrix ι ι ℂ →L[ℝ] EuclideanSpace ℂ ι :=
  LinearMap.toContinuousLinearMap (colLin z)

omit [DecidableEq ι] in
/-- Column extraction. -/
lemma colL_apply (z : ι) (M : Matrix ι ι ℂ) : colL z M = WithLp.toLp 2 (fun x => M x z) := rfl

section Deriv2

open scoped Matrix.Norms.Operator

set_option backward.isDefEq.respectTransparency false

/-- Derivative of the column `z` of the Heisenberg-evolved observable. -/
lemma hasDerivAt_heis_col (A X : Matrix ι ι ℝ) (s : ℝ) (z : ι) :
    HasDerivAt (fun r : ℝ => colL z (evo A (-r) * cplx X * evo A r))
      (colL z (evo A (-s) * (Complex.I • (cplx A * cplx X - cplx X * cplx A)) * evo A s)) s :=
  (colL z).hasFDerivAt.comp_hasDerivAt s (hasDerivAt_heis A X s)

end Deriv2

omit [DecidableEq ι] in
/-- The column of a product. -/
lemma col_mul (U W : Matrix ι ι ℂ) (z : ι) :
    (fun x => (U * W) x z) = U *ᵥ (fun y => W y z) := by
  funext x
  simp [Matrix.mul_apply, Matrix.mulVec, dotProduct]

/-- A unitary matrix preserves the sum of squared moduli. -/
lemma sum_norm_sq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ x, ‖(U *ᵥ v) x‖ ^ 2 = ∑ x, ‖v x‖ ^ 2 := by
  have h : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec, Matrix.mulVec_mulVec, hU,
      Matrix.one_mulVec]
  have e : ∀ w : ι → ℂ, star w ⬝ᵥ w = ((∑ x, ‖w x‖ ^ 2 : ℝ) : ℂ) := by
    intro w
    simp only [dotProduct, Pi.star_apply, Complex.ofReal_sum, Complex.ofReal_pow]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [Complex.star_def, mul_comm, Complex.mul_conj, Complex.normSq_eq_norm_sq]
    push_cast
    ring
  rw [e, e] at h
  exact_mod_cast h

/-- A unitary matrix preserves the Euclidean norm. -/
lemma norm_toLp_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ‖(WithLp.toLp 2 (U *ᵥ v) : EuclideanSpace ℂ ι)‖ =
      ‖(WithLp.toLp 2 v : EuclideanSpace ℂ ι)‖ := by
  rw [EuclideanSpace.norm_eq, EuclideanSpace.norm_eq]
  rw [sum_norm_sq_mulVec U hU v]

/-- The columns of a unitary propagator are unit vectors. -/
lemma sum_norm_sq_evo (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) (z : ι) :
    ∑ x, ‖evo A t x z‖ ^ 2 = 1 := by
  have h := sum_norm_sq_mulVec (evo A t) (evo_unitary A hA t) (fun y => if y = z then 1 else 0)
  have e : evo A t *ᵥ (fun y => if y = z then (1 : ℂ) else 0) = fun x => evo A t x z := by
    funext x
    simp [Matrix.mulVec, dotProduct]
  rw [e] at h
  rw [h]
  simp [apply_ite]

omit [DecidableEq ι] in
/-- **Schur test**: row and column `ℓ¹`-norms at most `r` give `‖C v‖ ≤ r ‖v‖`. -/
lemma schur (C : Matrix ι ι ℂ) (r : ℝ) (hr : 0 ≤ r) (hrow : ∀ i, ∑ j, ‖C i j‖ ≤ r)
    (hcol : ∀ j, ∑ i, ‖C i j‖ ≤ r) (v : ι → ℂ) :
    ∑ i, ‖(C *ᵥ v) i‖ ^ 2 ≤ r ^ 2 * ∑ i, ‖v i‖ ^ 2 := by
  have h1 : ∀ i, ‖(C *ᵥ v) i‖ ^ 2 ≤ r * ∑ j, ‖C i j‖ * ‖v j‖ ^ 2 := by
    intro i
    have e1 : ‖(C *ᵥ v) i‖ ≤ ∑ j, ‖C i j‖ * ‖v j‖ := by
      simp only [Matrix.mulVec, dotProduct]
      refine (norm_sum_le _ _).trans (le_of_eq ?_)
      simp
    have e2 : (∑ j, ‖C i j‖ * ‖v j‖) ^ 2 ≤ (∑ j, ‖C i j‖) * ∑ j, ‖C i j‖ * ‖v j‖ ^ 2 :=
      Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul _ (fun j _ => norm_nonneg _)
        (fun j _ => by positivity) (fun j _ => le_of_eq (by ring))
    calc ‖(C *ᵥ v) i‖ ^ 2 ≤ (∑ j, ‖C i j‖ * ‖v j‖) ^ 2 := by gcongr
      _ ≤ (∑ j, ‖C i j‖) * ∑ j, ‖C i j‖ * ‖v j‖ ^ 2 := e2
      _ ≤ r * ∑ j, ‖C i j‖ * ‖v j‖ ^ 2 := by
          gcongr
          exact hrow i
  calc ∑ i, ‖(C *ᵥ v) i‖ ^ 2 ≤ ∑ i, r * ∑ j, ‖C i j‖ * ‖v j‖ ^ 2 := Finset.sum_le_sum fun i _ => h1 i
    _ = r * ∑ j, (∑ i, ‖C i j‖) * ‖v j‖ ^ 2 := by
        rw [← Finset.mul_sum, Finset.sum_comm]
        simp [Finset.sum_mul]
    _ ≤ r * ∑ j, r * ‖v j‖ ^ 2 := by
        gcongr with j
        exact hcol j
    _ = r ^ 2 * ∑ i, ‖v i‖ ^ 2 := by rw [← Finset.mul_sum]; ring

end General

/-- `cplx` of a diagonal matrix. -/
lemma cplx_diagonal {ι : Type*} [DecidableEq ι] (d : ι → ℝ) :
    cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
  ext i j; by_cases h : i = j <;> simp [cplx, h]

/-- `schrH` is symmetric. -/
lemma schrH_isSymm (a b : ℤ) (w : ℤ → ℝ) : (schrH a b w).IsSymm := by
  refine Matrix.IsSymm.ext fun x y => ?_
  simp only [schrH]
  by_cases h : x.1 = y.1
  · rw [ite_eq_left h.symm, ite_eq_left h, h]
  · rw [ite_eq_right (Ne.symm h), ite_eq_right h, abs_sub_comm]

/-- Two integers at distance one. -/
lemma abs_sub_eq_one_iff (c v : ℤ) : |c - v| = 1 ↔ v = c + 1 ∨ v = c - 1 := by
  rw [abs_eq (by norm_num)]
  omega

/-- At most one site of `[a, b]` has a given position. -/
lemma sum_ite_val_le_one (a b d : ℤ) :
    ∑ v : Set.Icc a b, (if v.1 = d then (1 : ℝ) else 0) ≤ 1 := by
  rw [Finset.sum_boole]
  exact_mod_cast Finset.card_le_one.mpr fun x hx y hy => by
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx hy
    exact Subtype.ext (hx.trans hy.symm)

/-- Every site of `[a, b]` has at most two neighbours. -/
lemma sum_adj_le_two (a b c : ℤ) :
    ∑ v : Set.Icc a b, (if |c - v.1| = 1 then (1 : ℝ) else 0) ≤ 2 := by
  calc ∑ v : Set.Icc a b, (if |c - v.1| = 1 then (1 : ℝ) else 0)
      ≤ ∑ v : Set.Icc a b, ((if v.1 = c + 1 then (1 : ℝ) else 0) +
          (if v.1 = c - 1 then (1 : ℝ) else 0)) := by
        refine Finset.sum_le_sum fun v _ => ?_
        simp only [abs_sub_eq_one_iff]
        split_ifs <;> norm_num
        omega
    _ = ∑ v : Set.Icc a b, (if v.1 = c + 1 then (1 : ℝ) else 0) +
          ∑ v : Set.Icc a b, (if v.1 = c - 1 then (1 : ℝ) else 0) := Finset.sum_add_distrib
    _ ≤ 1 + 1 := add_le_add (sum_ite_val_le_one a b _) (sum_ite_val_le_one a b _)
    _ = 2 := by norm_num

/-- The position operator on `[a, b]`. -/
def posX (a b : ℤ) : Matrix (Set.Icc a b) (Set.Icc a b) ℝ := diagonal (fun y => ((y.1 : ℤ) : ℝ))

/-- The entries of `i [H, X]` have modulus at most one, and vanish unless `|u - v| = 1`. -/
lemma norm_comm_apply_le (a b : ℤ) (w : ℤ → ℝ) (u v : Set.Icc a b) :
    ‖(Complex.I • (cplx (schrH a b w) * cplx (posX a b) - cplx (posX a b) * cplx (schrH a b w)))
      u v‖ ≤ if |u.1 - v.1| = 1 then 1 else 0 := by
  have e : (Complex.I • (cplx (schrH a b w) * cplx (posX a b) -
      cplx (posX a b) * cplx (schrH a b w))) u v =
      Complex.I * ((schrH a b w u v * (((v.1 : ℤ) : ℝ) - ((u.1 : ℤ) : ℝ)) : ℝ) : ℂ) := by
    rw [posX, cplx_diagonal]
    simp only [Matrix.smul_apply, Matrix.sub_apply, Matrix.mul_diagonal,
      Matrix.diagonal_mul, smul_eq_mul]
    simp only [cplx, Matrix.map_apply]
    push_cast
    ring
  rw [e, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul]
  simp only [schrH]
  by_cases h1 : u.1 = v.1
  · rw [h1]
    simp
  · rw [ite_eq_right h1]
    by_cases h2 : |u.1 - v.1| = 1
    · rw [ite_eq_left h2, ite_eq_left h2, abs_neg, abs_one, one_mul, ← Int.cast_sub, ← Int.cast_abs,
        abs_sub_comm, h2, Int.cast_one]
    · rw [ite_eq_right h2, ite_eq_right h2, abs_zero, zero_mul]

section Ballistic

set_option backward.isDefEq.respectTransparency false

/-- **Ballistic bound** for the second moment: `‖X e^{-isH} δ_0‖² ≤ (2|s|)²`. -/
lemma second_moment_le (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) :
    ∑ y : Set.Icc a b, ((y.1 : ℤ) : ℝ) ^ 2 * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      (2 * |s|) ^ 2 := by
  set H := schrH a b w with hH
  set o := origin a b ha hb with ho
  set X := posX a b with hX
  have hsym : H.IsSymm := schrH_isSymm a b w
  set C := Complex.I • (cplx H * cplx X - cplx X * cplx H) with hC
  -- the derivative is bounded by `2`
  have hder : ∀ r : ℝ, ‖colL o (evo H (-r) * C * evo H r)‖ ≤ 2 := by
    intro r
    rw [colL_apply, Matrix.mul_assoc, col_mul, norm_toLp_mulVec _ (evo_unitary H hsym (-r)),
      col_mul, EuclideanSpace.norm_eq]
    have hrow : ∀ u, ∑ v, ‖C u v‖ ≤ 2 := fun u =>
      (Finset.sum_le_sum fun v _ => norm_comm_apply_le a b w u v).trans (sum_adj_le_two a b u.1)
    have hcol : ∀ v, ∑ u, ‖C u v‖ ≤ 2 := by
      intro v
      refine (Finset.sum_le_sum fun u _ => norm_comm_apply_le a b w u v).trans ?_
      simp_rw [abs_sub_comm _ v.1]
      exact sum_adj_le_two a b v.1
    have h1 := schur C 2 (by norm_num) hrow hcol (fun y => evo H r y o)
    rw [sum_norm_sq_evo H hsym r o, mul_one] at h1
    calc Real.sqrt (∑ x, ‖(WithLp.toLp 2 (C *ᵥ fun y => evo H r y o) : EuclideanSpace ℂ _) x‖ ^ 2)
        ≤ Real.sqrt (2 ^ 2) := Real.sqrt_le_sqrt h1
      _ = 2 := Real.sqrt_sq (by norm_num)
  -- mean value inequality
  have key := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := fun r : ℝ => colL o (evo H (-r) * cplx X * evo H r)) (s := Set.univ) (x := 0)
    (y := s) (fun r _ => (hasDerivAt_heis_col H X r o).hasDerivWithinAt)
    (fun r _ => hder r) convex_univ trivial trivial
  have hX0 : (fun x => cplx X x o) = 0 := by
    funext x
    rw [hX, posX, cplx_diagonal, Matrix.diagonal_apply]
    split_ifs with h
    · rw [h, ho, origin]; simp
    · rfl
  have h0 : colL o (evo H (-0) * cplx X * evo H 0) = 0 := by
    rw [neg_zero, evo_zero, Matrix.one_mul, Matrix.mul_one, colL_apply, hX0, WithLp.toLp_zero]
  simp only [h0, sub_zero, Real.norm_eq_abs] at key
  have hs : ‖colL o (evo H (-s) * cplx X * evo H s)‖ ^ 2 =
      ∑ y, ((y.1 : ℤ) : ℝ) ^ 2 * ‖evo H s y o‖ ^ 2 := by
    rw [colL_apply, Matrix.mul_assoc, col_mul, norm_toLp_mulVec _ (evo_unitary H hsym (-s)),
      EuclideanSpace.norm_eq, Real.sq_sqrt (Finset.sum_nonneg fun _ _ => by positivity)]
    refine Finset.sum_congr rfl fun y _ => ?_
    rw [hX, posX, cplx_diagonal]
    simp only [Matrix.diagonal_mul]
    rw [norm_mul, mul_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs]
  rw [← hs]
  exact pow_le_pow_left₀ (norm_nonneg _) key 2

end Ballistic

/-- `(x²)^{p/2} = |x|^p`. -/
lemma sq_rpow_half (x p : ℝ) : (x ^ 2) ^ (p / 2) = |x| ^ p := by
  rw [← sq_abs, ← Real.rpow_natCast, ← Real.rpow_mul (abs_nonneg x)]
  congr 1
  push_cast
  ring


section Series

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

open scoped Matrix.Norms.Operator in
/-- Entrywise exponential series of the propagator. -/
lemma hasSum_evo_apply (A : Matrix ι ι ℝ) (t : ℝ) (i j : ι) :
    HasSum (fun k : ℕ => ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) *
      (((A ^ k) i j : ℝ) : ℂ)) (evo A t i j) := by
  have h := NormedSpace.exp_series_hasSum_exp' (𝕂 := ℂ) ((-(t : ℂ) * Complex.I) • cplx A)
  have h2 := Pi.hasSum.mp (Pi.hasSum.mp h i) j
  convert h2 using 1
  swap
  · rfl
  funext k
  have hc : (cplx A) ^ k = cplx (A ^ k) := by
    unfold cplx
    exact (Matrix.map_pow A Complex.ofRealHom k).symm
  rw [smul_pow, hc]
  simp only [Matrix.smul_apply, cplx, Matrix.map_apply, smul_eq_mul]
  ring

/-- The coefficients of the exponential series: `‖(-it)^k/k!‖ = |t|^k/k!`. -/
lemma norm_coef (t : ℝ) (k : ℕ) :
    ‖(-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)‖ = |t| ^ k / k.factorial := by
  rw [norm_div, norm_pow, norm_mul, norm_neg, Complex.norm_real, Complex.norm_I, mul_one,
    Complex.norm_natCast, Real.norm_eq_abs]

/-- A term of the exponential series with an entry bounded by `c^k`. -/
lemma norm_term_le (t : ℝ) (k : ℕ) (x c : ℝ) (hx : |x| ≤ c ^ k) :
    ‖((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((x : ℝ) : ℂ)‖ ≤
      (c * |t|) ^ k / k.factorial := by
  rw [norm_mul, norm_coef, Complex.norm_real, Real.norm_eq_abs, mul_pow]
  calc |t| ^ k / k.factorial * |x| ≤ |t| ^ k / k.factorial * c ^ k :=
        mul_le_mul_of_nonneg_left hx (by positivity)
    _ = c ^ k * |t| ^ k / k.factorial := by ring

/-- Entries of powers of a matrix whose rows have `ℓ¹`-norm at most `c`. -/
lemma abs_pow_apply_le (A : Matrix ι ι ℝ) (c : ℝ) (hc : 0 ≤ c)
    (hA : ∀ i, ∑ l, |A i l| ≤ c) (k : ℕ) (i j : ι) : |(A ^ k) i j| ≤ c ^ k := by
  induction k generalizing i with
  | zero =>
    rw [pow_zero, pow_zero, Matrix.one_apply]
    split_ifs <;> simp
  | succ k ih =>
    rw [pow_succ', Matrix.mul_apply, pow_succ']
    calc |∑ l, A i l * (A ^ k) l j| ≤ ∑ l, |A i l * (A ^ k) l j| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ l, |A i l| * c ^ k := by
          refine Finset.sum_le_sum fun l _ => ?_
          rw [abs_mul]
          exact mul_le_mul_of_nonneg_left (ih l) (abs_nonneg _)
      _ = (∑ l, |A i l|) * c ^ k := by rw [Finset.sum_mul]
      _ ≤ c * c ^ k := mul_le_mul_of_nonneg_right (hA i) (by positivity)

end Series

/-- The rows of `H_Λ` have `ℓ¹`-norm at most `2 + |w(y)|`. -/
lemma sum_abs_schrH_le (a b : ℤ) (w : ℤ → ℝ) (y : Set.Icc a b) :
    ∑ u, |schrH a b w y u| ≤ |w y.1| + 2 := by
  have h1 : ∀ u : Set.Icc a b, |schrH a b w y u| ≤
      |w y.1| * (if u.1 = y.1 then (1 : ℝ) else 0) + (if |y.1 - u.1| = 1 then (1 : ℝ) else 0) := by
    intro u
    simp only [schrH]
    by_cases h : y.1 = u.1
    · rw [ite_eq_left h, ite_eq_left h.symm, h, sub_self]
      simp
    · rw [ite_eq_right h, ite_eq_right (Ne.symm h), mul_zero, zero_add]
      split_ifs <;> simp
  calc ∑ u, |schrH a b w y u|
      ≤ ∑ u : Set.Icc a b, (|w y.1| * (if u.1 = y.1 then (1 : ℝ) else 0) +
          (if |y.1 - u.1| = 1 then (1 : ℝ) else 0)) := Finset.sum_le_sum fun u _ => h1 u
    _ = |w y.1| * ∑ u : Set.Icc a b, (if u.1 = y.1 then (1 : ℝ) else 0) +
          ∑ u : Set.Icc a b, (if |y.1 - u.1| = 1 then (1 : ℝ) else 0) := by
        rw [Finset.sum_add_distrib, Finset.mul_sum]
    _ ≤ |w y.1| * 1 + 2 := by
        gcongr
        · exact sum_ite_val_le_one a b y.1
        · exact sum_adj_le_two a b y.1
    _ = |w y.1| + 2 := by ring

/-- The inclusion `[a, b] → [a', b']`. -/
def incl {a b a' b' : ℤ} (ha' : a' ≤ a) (hb' : b ≤ b') (y : Set.Icc a b) : Set.Icc a' b' :=
  ⟨y.1, ha'.trans y.2.1, y.2.2.trans hb'⟩

/-- The inclusion is injective. -/
lemma incl_injective {a b a' b' : ℤ} (ha' : a' ≤ a) (hb' : b ≤ b') :
    Function.Injective (incl ha' hb') := fun x y h => by
  have h' := congrArg Subtype.val h
  exact Subtype.ext h'

/-- The entries of `H_{Λ'}` on `Λ` are those of `H_Λ`. -/
lemma schrH_incl {a b a' b' : ℤ} (ha' : a' ≤ a) (hb' : b ≤ b') (w : ℤ → ℝ) (u x : Set.Icc a b) :
    schrH a' b' w (incl ha' hb' u) (incl ha' hb' x) = schrH a b w u x := rfl

/-- Nonzero entries of `H` join sites at distance at most one. -/
lemma schrH_eq_zero {a b : ℤ} (w : ℤ → ℝ) (u x : Set.Icc a b) (h : 1 < |u.1 - x.1|) :
    schrH a b w u x = 0 := by
  simp only [schrH]
  rw [ite_eq_right (fun e => by rw [e, sub_self, abs_zero] at h; norm_num at h),
    ite_eq_right (fun e => by rw [e] at h; norm_num at h)]

/-- **Finite propagation**: the powers of `H_Λ` and `H_{Λ'}` agree near the origin. -/
lemma pow_schrH_incl (a b a' b' : ℤ) (ha' : a' ≤ a) (hb' : b ≤ b') (w : ℤ → ℝ) (k : ℕ) :
    ∀ x : Set.Icc a b, |x.1| + k ≤ min b (-a) → ∀ y : Set.Icc a b,
      (schrH a b w ^ k) y x = (schrH a' b' w ^ k) (incl ha' hb' y) (incl ha' hb' x) := by
  induction k with
  | zero =>
    intro x _ y
    simp only [pow_zero, Matrix.one_apply]
    by_cases h : y = x
    · rw [ite_eq_left h, ite_eq_left (congrArg _ h)]
    · rw [ite_eq_right h, ite_eq_right (fun e => h (incl_injective ha' hb' e))]
  | succ k ih =>
    intro x hx y
    have hx1 : |x.1| + 1 ≤ min b (-a) := by push_cast at hx; linarith
    have hxa : a + 1 ≤ x.1 := by
      have := neg_abs_le x.1
      have := min_le_right b (-a)
      linarith
    have hxb : x.1 + 1 ≤ b := by
      have := le_abs_self x.1
      have := min_le_left b (-a)
      linarith
    rw [pow_succ, pow_succ, Matrix.mul_apply, Matrix.mul_apply]
    refine Fintype.sum_of_injective (incl ha' hb') (incl_injective ha' hb') _ _ ?_ ?_
    · intro u hu
      have hux : 1 < |u.1 - (incl ha' hb' x).1| := by
        by_contra hc
        have h1 := abs_le.mp (not_lt.mp hc)
        exact hu ⟨⟨u.1, by simp only [incl] at h1; linarith, by simp only [incl] at h1; linarith⟩,
          rfl⟩
      rw [schrH_eq_zero w u (incl ha' hb' x) hux, mul_zero]
    · intro u
      rw [schrH_incl]
      by_cases hu : 1 < |u.1 - x.1|
      · rw [schrH_eq_zero w u x hu, mul_zero, mul_zero]
      · have h1 := abs_le.mp (not_lt.mp hu)
        rw [ih u ?_ y]
        have := abs_sub_abs_le_abs_sub u.1 x.1
        push_cast at hx
        linarith

end DimerTLAux

/-- **Finite speed, uniformly in the potential** (paper: Lemma `lem:dimerTL`(i)). -/
theorem dimerTL_i (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ≤
      ∑' k : ℕ, if (y.1).natAbs ≤ k then (2 * |s|) ^ k / (k.factorial : ℝ) else 0 := by
  sorry

/-- The tail bound in (i): `∑_{k ≥ m} c^k/k! ≤ (c^m / m!) e^c`. -/
theorem dimerTL_i' (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ≤
      (2 * |s|) ^ (y.1).natAbs / ((y.1).natAbs.factorial : ℝ) * Real.exp (2 * |s|) := by
  sorry

/-- **Thermodynamic limit** (paper: Lemma `lem:dimerTL`(ii), finite `Λ'`). -/
theorem dimerTL_ii (a b a' b' : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (ha' : a' ≤ a) (hb' : b ≤ b')
    (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, a' ≤ x → x ≤ b' → |w x| ≤ W) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb) -
        evo (schrH a' b' w) s ⟨y.1, ha'.trans y.2.1, y.2.2.trans hb'⟩
          (origin a' b' (ha'.trans ha) (hb.trans hb'))‖ ≤
      2 * ∑' k : ℕ, if (min b (-a)).toNat + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ) else 0 := by
  set o := origin a b ha hb with ho
  set m := (min b (-a)).toNat with hm
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0 (ha'.trans ha) (hb.trans hb'))
  have hy' : (⟨y.1, ha'.trans y.2.1, y.2.2.trans hb'⟩ : Set.Icc a' b') =
      DimerTLAux.incl ha' hb' y := rfl
  have ho' : origin a' b' (ha'.trans ha) (hb.trans hb') = DimerTLAux.incl ha' hb' o := rfl
  rw [hy', ho']
  have h3 := (DimerTLAux.hasSum_evo_apply (schrH a b w) s y o).sub
    (DimerTLAux.hasSum_evo_apply (schrH a' b' w) s (DimerTLAux.incl ha' hb' y)
      (DimerTLAux.incl ha' hb' o))
  set f : ℕ → ℝ := fun k => if m + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ) else 0
    with hf
  have hf0 : ∀ k, 0 ≤ f k := fun k => by
    simp only [hf]
    split_ifs <;> positivity
  have hfs : Summable f := by
    refine Summable.of_nonneg_of_le hf0 (fun k => ?_) (Real.summable_pow_div_factorial ((2 + W) * |s|))
    simp only [hf]
    split_ifs
    · exact le_rfl
    · positivity
  -- row bounds
  have hrow : ∀ u : Set.Icc a b, ∑ l, |schrH a b w u l| ≤ 2 + W := fun u =>
    (DimerTLAux.sum_abs_schrH_le a b w u).trans
      (by linarith [hw u.1 (ha'.trans u.2.1) (u.2.2.trans hb')])
  have hrow' : ∀ u : Set.Icc a' b', ∑ l, |schrH a' b' w u l| ≤ 2 + W := fun u =>
    (DimerTLAux.sum_abs_schrH_le a' b' w u).trans (by linarith [hw u.1 u.2.1 u.2.2])
  refine (h3.norm_le_of_bounded (hfs.mul_left 2).hasSum fun k => ?_).trans
    (le_of_eq tsum_mul_left)
  by_cases hk : m + 1 ≤ k
  · have e1 := DimerTLAux.norm_term_le s k _ (2 + W)
      (DimerTLAux.abs_pow_apply_le _ _ (by linarith) hrow k y o)
    have e2 := DimerTLAux.norm_term_le s k _ (2 + W)
      (DimerTLAux.abs_pow_apply_le _ _ (by linarith) hrow' k (DimerTLAux.incl ha' hb' y)
        (DimerTLAux.incl ha' hb' o))
    simp only [hf]
    rw [ite_eq_left hk]
    refine (norm_sub_le _ _).trans ?_
    linarith
  · have hk' : |o.1| + (k : ℤ) ≤ min b (-a) := by
      have h0 : (0 : ℤ) ≤ min b (-a) := le_min hb (by linarith)
      have : (k : ℤ) ≤ m := by exact_mod_cast (by omega : k ≤ m)
      rw [hm, Int.toNat_of_nonneg h0] at this
      simp only [ho, origin, abs_zero, zero_add]
      exact this
    rw [DimerTLAux.pow_schrH_incl a b a' b' ha' hb' w k o hk' y, sub_self, norm_zero]
    positivity

/-- **Ballistic bound** (paper: Lemma `lem:dimerTL`(iii)). -/
theorem dimerTL_iii (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (p : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) :
    ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      (2 * |s|) ^ p := by
  set φ := fun y : Set.Icc a b => ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 with hφ
  have hsum : ∑ y, φ y = 1 :=
    DimerTLAux.sum_norm_sq_evo _ (DimerTLAux.schrH_isSymm a b w) s _
  have h2 := DimerTLAux.second_moment_le a b ha hb w s
  have hconc := Real.concaveOn_rpow (p := p / 2) (by positivity) (by linarith)
  have hJ := hconc.le_map_sum (t := Finset.univ) (w := φ)
    (p := fun y : Set.Icc a b => ((y.1 : ℤ) : ℝ) ^ 2)
    (fun y _ => by positivity) hsum (fun y _ => Set.mem_Ici.mpr (by positivity))
  calc ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2
      = ∑ y, φ y • (((y.1 : ℤ) : ℝ) ^ 2) ^ (p / 2) := by
        refine Finset.sum_congr rfl fun y _ => ?_
        rw [DimerTLAux.sq_rpow_half, smul_eq_mul, mul_comm]
    _ ≤ (∑ y, φ y • ((y.1 : ℤ) : ℝ) ^ 2) ^ (p / 2) := hJ
    _ ≤ ((2 * |s|) ^ 2) ^ (p / 2) := by
        refine Real.rpow_le_rpow (Finset.sum_nonneg fun y _ => by positivity) ?_ (by positivity)
        refine le_of_eq_of_le (Finset.sum_congr rfl fun y _ => ?_) h2
        rw [smul_eq_mul, mul_comm]
    _ = (2 * |s|) ^ p := by
        rw [DimerTLAux.sq_rpow_half, abs_of_nonneg (by positivity)]

end DFL
