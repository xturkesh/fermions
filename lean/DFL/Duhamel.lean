import DFL.OneBody

/-!
# Duhamel sensitivity to distant charges (paper: Lemma `lem:duhamel`)

Flipping the charges on a set `D` changes the one-body propagator, acting on a site `z`, by an
amount controlled by the propagation from `z` to `D`:
`‖(U_q(s) - U_{q'}(s)) δ_z‖ ≤ 4|h| |s| ∑_{x ∈ D} Q_{q'}(x, z)`, where `U_q(s) = e^{-ish_q}`.
Under Hypothesis (EC) the Bernoulli average is at most `4|h| |s| ∑_{x ∈ D} F(|x - z|)`: charges at
distance `R` from `z` influence the dynamics near `z` only after times of order `e^{ηR}`.

Vector norms are Euclidean (`EuclideanSpace ℂ (Fin N)`).

Proof strategy (no integrals): for `f(r) = U_q(s - r) U_{q'}(r) δ_z` we have
`f'(r) = i U_q(s - r) (h_q - h_{q'}) U_{q'}(r) δ_z`, whose norm is at most
`4|h| ∑_{x ∈ D} |U_{q'}(r)_{xz}| ≤ 4|h| ∑_{x ∈ D} Q_{q'}(x, z)` (unitarity of `U_q`), and we apply
the mean value inequality between `r = 0` and `r = s`.
-/

open Matrix Finset NormedSpace

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Flipping a configuration on a set `D`. -/
def flipOn (D : Finset (Fin N)) (q : Cfg N) : Cfg N := fun x => if x ∈ D then q x + 1 else q x

/-- The Euclidean norm of the column `U δ_z` of a matrix. -/
def colNorm (U : Matrix (Fin N) (Fin N) ℂ) (z : Fin N) : ℝ :=
  ‖(WithLp.toLp 2 (fun x => U x z) : EuclideanSpace ℂ (Fin N))‖

namespace DuhamelAux

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

/-- Derivative of the Duhamel interpolation `r ↦ U_A(s - r) U_B(r)`. -/
lemma hasDerivAt_prod (A B : Matrix ι ι ℝ) (s r : ℝ) :
    HasDerivAt (fun r : ℝ => evo A (s - r) * evo B r)
      (evo A (s - r) * (Complex.I • (cplx A - cplx B)) * evo B r) r := by
  have h1 := (hasDerivAt_evo A (s - r)).scomp r ((hasDerivAt_id r).const_sub s)
  have h2 := hasDerivAt_evo' B r
  have h3 := h1.mul h2
  convert h3 using 1
  · rfl
  · rfl
  · simp only [Function.comp_apply, neg_smul, one_smul, smul_neg, Matrix.mul_neg, Matrix.neg_mul,
      Matrix.mul_smul, Matrix.smul_mul, smul_sub, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc]
    abel

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

/-- Derivative of the column `z` of the Duhamel interpolation. -/
lemma hasDerivAt_col (A B : Matrix ι ι ℝ) (s r : ℝ) (z : ι) :
    HasDerivAt (fun r : ℝ => colL z (evo A (s - r) * evo B r))
      (colL z (evo A (s - r) * (Complex.I • (cplx A - cplx B)) * evo B r)) r :=
  (colL z).hasFDerivAt.comp_hasDerivAt r (hasDerivAt_prod A B s r)

end Deriv

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

omit [DecidableEq ι] in
/-- The Euclidean norm is bounded by the `ℓ¹` norm. -/
lemma norm_toLp_le_sum (v : ι → ℂ) :
    ‖(WithLp.toLp 2 v : EuclideanSpace ℂ ι)‖ ≤ ∑ x, ‖v x‖ := by
  rw [EuclideanSpace.norm_eq, Real.sqrt_le_iff]
  exact ⟨Finset.sum_nonneg fun x _ => norm_nonneg _,
    Finset.sum_sq_le_sq_sum_of_nonneg fun x _ => norm_nonneg _⟩

/-- Entries of `i (h_q - h_{q'}) W`: the difference of the generators is diagonal. -/
lemma diff_mul_apply (J h : ℝ) (q q' : Cfg N) (W : Matrix (Fin N) (Fin N) ℂ) (x z : Fin N) :
    ((Complex.I • (cplx (hq J h q) - cplx (hq J h q'))) * W) x z =
      Complex.I * ((2 * h * (sgn (q x) - sgn (q' x)) : ℝ) : ℂ) * W x z := by
  rw [Matrix.mul_apply, Finset.sum_eq_single x]
  · simp only [Matrix.smul_apply, Matrix.sub_apply, cplx, Matrix.map_apply, hq, Matrix.add_apply,
      diagonal_apply_eq, smul_eq_mul]
    push_cast
    ring
  · intro y _ hy
    simp [cplx, hq, diagonal_apply_ne _ (Ne.symm hy)]
  · simp

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- The generators differ by `±4h` on `D` and agree off `D`. -/
lemma abs_diff_le (h : ℝ) (q : Cfg N) (D : Finset (Fin N)) (x : Fin N) :
    |2 * h * (sgn (q x) - sgn (flipOn D q x))| ≤ if x ∈ D then 4 * |h| else 0 := by
  unfold flipOn
  by_cases hx : x ∈ D
  · have h11 : (1 : ZMod 2) + 1 = 0 := by decide
    rcases zmod2_cases (q x) with h0 | h0
    · simp only [hx, ite_true, h0, zero_add, sgn]
      norm_num [abs_mul]
      linarith
    · simp only [hx, ite_true, h0, h11, sgn]
      norm_num [abs_mul]
      linarith
  · simp [hx]

/-- Pointwise bound on the derivative of the Duhamel interpolation. -/
lemma norm_deriv_le (J h : ℝ) (q : Cfg N) (D : Finset (Fin N)) (z : Fin N) (s r : ℝ) :
    ‖colL z (evo (hq J h q) (s - r) *
        (Complex.I • (cplx (hq J h q) - cplx (hq J h (flipOn D q)))) *
          evo (hq J h (flipOn D q)) r)‖ ≤
      4 * |h| * ∑ x ∈ D, corrQ (hq J h (flipOn D q)) (hq_isHermitian J h (flipOn D q)) x z := by
  set B := hq J h (flipOn D q)
  rw [colL_apply, Matrix.mul_assoc, col_mul,
    norm_toLp_mulVec _ (evo_unitary _ (hq_isSymm J h q) _)]
  refine (norm_toLp_le_sum _).trans ?_
  have hpt : ∀ x, ‖((Complex.I • (cplx (hq J h q) - cplx B)) * evo B r) x z‖ ≤
      if x ∈ D then 4 * |h| * corrQ B (hq_isHermitian J h (flipOn D q)) x z else 0 := by
    intro x
    rw [diff_mul_apply, norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real,
      Real.norm_eq_abs]
    have h1 := abs_diff_le h q D x
    have h2 := norm_evo_le_corrQ B (hq_isHermitian J h (flipOn D q)) r x z
    have h3 := norm_nonneg (evo B r x z)
    split_ifs with hx
    · simp only [hx, ↓reduceIte] at h1
      exact mul_le_mul h1 h2 h3 (by positivity)
    · simp only [hx, ↓reduceIte] at h1
      have : |2 * h * (sgn (q x) - sgn (flipOn D q x))| = 0 :=
        le_antisymm h1 (abs_nonneg _)
      rw [this, zero_mul]
  refine (Finset.sum_le_sum fun x _ => hpt x).trans (le_of_eq ?_)
  rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.mul_sum]

end DuhamelAux

open DuhamelAux in
/-- **Duhamel bound** (deterministic): if `q'` is `q` flipped on `D`, then
`‖(U_q(s) - U_{q'}(s)) δ_z‖ ≤ 4|h| |s| ∑_{x ∈ D} Q_{q'}(x, z)`. -/
theorem colNorm_evo_sub_le (J h : ℝ) (q : Cfg N) (D : Finset (Fin N)) (z : Fin N) (s : ℝ) :
    colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z ≤
      4 * |h| * |s| *
        ∑ x ∈ D, corrQ (hq J h (flipOn D q)) (hq_isHermitian J h (flipOn D q)) x z := by
  set A := hq J h q with hA
  set B := hq J h (flipOn D q) with hB
  set C := 4 * |h| * ∑ x ∈ D, corrQ B (hq_isHermitian J h (flipOn D q)) x z with hC
  have key := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le
    (f := fun r : ℝ => colL z (evo A (s - r) * evo B r)) (s := Set.univ) (x := 0) (y := s)
    (fun r _ => (hasDerivAt_col A B s r z).hasDerivWithinAt)
    (fun r _ => norm_deriv_le J h q D z s r) convex_univ trivial trivial
  simp only [sub_zero, evo_zero, Matrix.mul_one, Matrix.one_mul, sub_self] at key
  have hcol : colNorm (evo A s - evo B s) z = ‖colL z (evo B s) - colL z (evo A s)‖ := by
    rw [colNorm, ← map_sub, colL_apply, ← norm_neg, ← WithLp.toLp_neg]
    congr 2
    funext x
    simp
  rw [hcol]
  calc ‖colL z (evo B s) - colL z (evo A s)‖ ≤ C * ‖s‖ := key
    _ = 4 * |h| * |s| * ∑ x ∈ D, corrQ B (hq_isHermitian J h (flipOn D q)) x z := by
      rw [hC, Real.norm_eq_abs]
      ring

/-- Flipping on `D` is an involution. -/
private lemma flipOn_flipOn (D : Finset (Fin N)) (q : Cfg N) : flipOn D (flipOn D q) = q := by
  funext x
  unfold flipOn
  by_cases hx : x ∈ D
  · simp only [hx, ite_true, add_assoc]
    have h11 : (1 : ZMod 2) + 1 = 0 := by decide
    rw [h11, add_zero]
  · simp [hx]

/-- **Averaged Duhamel bound** under Hypothesis (EC), uniformly for `|s| ≤ t`. -/
theorem avg_colNorm_evo_sub_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (D : Finset (Fin N))
    (z : Fin N) (t s : ℝ) (hs : |s| ≤ t) :
    avgQ (fun q : Cfg N => colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z) ≤
      4 * |h| * t * ∑ x ∈ D, F (Nat.dist x.val z.val) := by
  have hpos : (0 : ℝ) ≤ 1 / 2 ^ N := by positivity
  -- pointwise Duhamel bound
  have h1 : avgQ (fun q : Cfg N => colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z) ≤
      avgQ (fun q : Cfg N => 4 * |h| * |s| *
        ∑ x ∈ D, corrQ (hq J h (flipOn D q)) (hq_isHermitian J h (flipOn D q)) x z) := by
    unfold avgQ
    exact mul_le_mul_of_nonneg_left
      (Finset.sum_le_sum fun q _ => colNorm_evo_sub_le J h q D z s) hpos
  -- reindex by the involution `flipOn D`
  let e : Cfg N ≃ Cfg N := ⟨flipOn D, flipOn D, flipOn_flipOn D, flipOn_flipOn D⟩
  have h2 : avgQ (fun q : Cfg N => 4 * |h| * |s| *
        ∑ x ∈ D, corrQ (hq J h (flipOn D q)) (hq_isHermitian J h (flipOn D q)) x z) =
      4 * |h| * |s| * ∑ x ∈ D, avgQ (fun q : Cfg N =>
        corrQ (hq J h q) (hq_isHermitian J h q) x z) := by
    unfold avgQ
    simp only [Finset.mul_sum]
    conv_lhs => rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun x _ => ?_
    conv_rhs => rw [← Equiv.sum_comp e]
    refine Finset.sum_congr rfl fun q _ => ?_
    simp only [e, Equiv.coe_fn_mk]
    ring
  -- Hypothesis (EC)
  have h3 : ∑ x ∈ D, avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤
      ∑ x ∈ D, F (Nat.dist x.val z.val) :=
    Finset.sum_le_sum fun x _ => hEC.bound N x z
  have hF : 0 ≤ ∑ x ∈ D, F (Nat.dist x.val z.val) :=
    Finset.sum_nonneg fun x _ => hEC.nonneg _
  have hh : 0 ≤ 4 * |h| := by positivity
  calc avgQ (fun q : Cfg N => colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z)
      ≤ 4 * |h| * |s| * ∑ x ∈ D, avgQ (fun q : Cfg N =>
          corrQ (hq J h q) (hq_isHermitian J h q) x z) := h1.trans (le_of_eq h2)
    _ ≤ 4 * |h| * |s| * ∑ x ∈ D, F (Nat.dist x.val z.val) :=
        mul_le_mul_of_nonneg_left h3 (by positivity)
    _ ≤ 4 * |h| * t * ∑ x ∈ D, F (Nat.dist x.val z.val) :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hs hh) hF

end DFL
