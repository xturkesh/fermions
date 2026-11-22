import Mathlib.Data.Matrix.Basis
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Algebra.Order.Star.Basic
import Mathlib.Tactic.Module
import Mathlib.Tactic.FieldSimp

/-!
# Projective measurement of a site occupation on a Slater determinant

A pure, particle-number-conserving fermionic Gaussian state (Slater determinant) on a finite
set of sites `n` is encoded in its correlation matrix `C j k = ⟨c_j† c_k⟩`, an orthogonal
projector (`Cᴴ = C`, `C * C = C`) whose trace is the particle number.

Measuring the occupation `n_i` gives `1` with Born probability `C i i`, after which the
correlation matrix is

  `measOcc C i = C - (C i i)⁻¹ • (C * P i * C) + P i`,   `P i = single i i 1`,

and gives `0` with probability `1 - C i i`, after which it is the particle–hole conjugate
update `measEmp C i = 1 - measOcc (1 - C) i`.

Main results (purely algebraic):
* `idem_update`: the abstract identity behind the update rule;
* `measOcc_mul_self`, `measEmp_mul_self`: the updated matrices are again projectors
  (Gaussianity and purity are preserved by the measurement);
* `measOcc_isHermitian`, `measEmp_isHermitian`;
* `measOcc_apply_self`, `measOcc_apply_row`, `measOcc_apply_col` (and `measEmp_…`): the
  measured site is left in a definite occupation and decoupled from the rest;
* `trace_measOcc`, `trace_measEmp`: the particle number is conserved;
* `born_prob_nonneg`, `born_prob_le_one`: Born probabilities lie in `[0, 1]`;
* `born_average`: the Born average of the post-measurement correlation matrix is the dephased
  matrix `(1 - P i) * C * (1 - P i) + P i * C * P i`, i.e. averaged two-point functions follow
  the *linear* dephasing (Lindblad) evolution: a measurement-induced transition can only show
  up in quantities that are nonlinear in the state (entanglement, cumulants, …).
-/

open Matrix

namespace MonitoredFermions

/-! ### An abstract identity -/

section Abstract

variable {𝕜 R : Type*} [Field 𝕜] [Ring R] [Module 𝕜 R] [IsScalarTower 𝕜 R R]
  [SMulCommClass 𝕜 R R]

/-- If `C`, `P` are idempotents with `P * C * P = a • P`, `a ≠ 0`, then
`C - a⁻¹ • (C * P * C) + P` is idempotent. -/
theorem idem_update {C P : R} {a : 𝕜} (ha : a ≠ 0) (hC : C * C = C) (hP : P * P = P)
    (hPCP : P * C * P = a • P) :
    (C - a⁻¹ • (C * P * C) + P) * (C - a⁻¹ • (C * P * C) + P) =
      C - a⁻¹ • (C * P * C) + P := by
  have hC' : ∀ Y : R, C * (C * Y) = C * Y := fun Y => by rw [← mul_assoc, hC]
  have hPCP₀ : P * (C * P) = a • P := by rw [← mul_assoc, hPCP]
  have hPCP' : ∀ Y : R, P * (C * (P * Y)) = a • (P * Y) := fun Y => by
    rw [← mul_assoc, ← mul_assoc, hPCP, smul_mul_assoc]
  simp only [add_mul, mul_add, sub_mul, mul_sub, smul_mul_assoc, mul_smul_comm, mul_assoc,
    hC, hP, hC', hPCP₀, hPCP', smul_smul]
  match_scalars <;> field_simp <;> ring

lemma one_sub_mul_self {Q : R} (hQ : Q * Q = Q) : (1 - Q) * (1 - Q) = 1 - Q := by
  simp only [sub_mul, mul_sub, one_mul, mul_one, hQ]
  abel

end Abstract

/-! ### Matrices -/

section Matrices

variable {n 𝕜 : Type*} [Fintype n] [DecidableEq n] [Field 𝕜]

/-- Projector onto the single-particle orbital localized at site `i`. -/
def siteProj (i : n) : Matrix n n 𝕜 := single i i 1

omit [Fintype n] in
lemma siteProj_apply (i j k : n) :
    (siteProj i : Matrix n n 𝕜) j k = if i = j ∧ i = k then 1 else 0 := rfl

lemma siteProj_mul_self (i : n) : (siteProj i : Matrix n n 𝕜) * siteProj i = siteProj i := by
  simp [siteProj]

lemma siteProj_sandwich (i : n) (C : Matrix n n 𝕜) :
    siteProj i * C * siteProj i = C i i • siteProj i := by
  simp [siteProj]

/-- Entries of `M * single i i 1 * N` (rectangular version). -/
lemma mul_single_mul_apply {l m o : Type*} [Fintype m] [DecidableEq m]
    (M : Matrix l m 𝕜) (N : Matrix m o 𝕜) (i : m) (j : l) (k : o) :
    (M * single i i (1 : 𝕜) * N : Matrix l o 𝕜) j k = M j i * N i k := by
  rw [Matrix.mul_apply, Finset.sum_eq_single i]
  · simp
  · intro x _ hx
    simp [hx]
  · simp

lemma mul_siteProj_mul_apply (C D : Matrix n n 𝕜) (i j k : n) :
    (C * siteProj i * D : Matrix n n 𝕜) j k = C j i * D i k :=
  mul_single_mul_apply C D i j k

/-- Correlation matrix after the measurement outcome `n_i = 1`. -/
noncomputable def measOcc (C : Matrix n n 𝕜) (i : n) : Matrix n n 𝕜 :=
  C - (C i i)⁻¹ • (C * siteProj i * C) + siteProj i

/-- Correlation matrix after the measurement outcome `n_i = 0` (particle–hole conjugate). -/
noncomputable def measEmp (C : Matrix n n 𝕜) (i : n) : Matrix n n 𝕜 :=
  1 - measOcc (1 - C) i

omit [Fintype n] in
lemma one_sub_apply_self (C : Matrix n n 𝕜) (i : n) : (1 - C) i i = 1 - C i i := by
  simp

/-- Outcome `1`: the post-measurement correlation matrix is again a projector. -/
theorem measOcc_mul_self {C : Matrix n n 𝕜} (hC : C * C = C) {i : n} (h : C i i ≠ 0) :
    measOcc C i * measOcc C i = measOcc C i :=
  idem_update h hC (siteProj_mul_self i) (siteProj_sandwich i C)

/-- Outcome `0`: the post-measurement correlation matrix is again a projector. -/
theorem measEmp_mul_self {C : Matrix n n 𝕜} (hC : C * C = C) {i : n} (h : C i i ≠ 1) :
    measEmp C i * measEmp C i = measEmp C i := by
  have h' : (1 - C) i i ≠ 0 := by
    rw [one_sub_apply_self]
    exact sub_ne_zero.mpr (Ne.symm h)
  exact one_sub_mul_self (measOcc_mul_self (one_sub_mul_self hC) h')

/-- After the outcome `1`, site `i` is occupied with certainty. -/
theorem measOcc_apply_self (C : Matrix n n 𝕜) {i : n} (h : C i i ≠ 0) :
    measOcc C i i i = 1 := by
  have e : (C i i)⁻¹ * (C i i * C i i) = C i i := by field_simp
  simp only [measOcc, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    mul_siteProj_mul_apply, smul_eq_mul, siteProj_apply, and_self, ite_true, e]
  ring

/-- After the outcome `1`, site `i` is decoupled from every other site (row). -/
theorem measOcc_apply_row (C : Matrix n n 𝕜) {i j : n} (h : C i i ≠ 0) (hij : i ≠ j) :
    measOcc C i i j = 0 := by
  have e : (C i i)⁻¹ * (C i i * C i j) = C i j := by field_simp
  simp only [measOcc, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    mul_siteProj_mul_apply, smul_eq_mul, siteProj_apply, hij, and_false, ite_false, e]
  ring

/-- After the outcome `1`, site `i` is decoupled from every other site (column). -/
theorem measOcc_apply_col (C : Matrix n n 𝕜) {i j : n} (h : C i i ≠ 0) (hij : i ≠ j) :
    measOcc C i j i = 0 := by
  have e : (C i i)⁻¹ * (C j i * C i i) = C j i := by field_simp
  simp only [measOcc, Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply,
    mul_siteProj_mul_apply, smul_eq_mul, siteProj_apply, hij, false_and, ite_false, e]
  ring

/-- After the outcome `0`, site `i` is empty with certainty. -/
theorem measEmp_apply_self (C : Matrix n n 𝕜) {i : n} (h : C i i ≠ 1) :
    measEmp C i i i = 0 := by
  have h' : (1 - C) i i ≠ 0 := by
    rw [one_sub_apply_self]
    exact sub_ne_zero.mpr (Ne.symm h)
  simp [measEmp, Matrix.sub_apply, measOcc_apply_self (1 - C) h']

theorem measEmp_apply_row (C : Matrix n n 𝕜) {i j : n} (h : C i i ≠ 1) (hij : i ≠ j) :
    measEmp C i i j = 0 := by
  have h' : (1 - C) i i ≠ 0 := by
    rw [one_sub_apply_self]
    exact sub_ne_zero.mpr (Ne.symm h)
  simp [measEmp, Matrix.sub_apply, measOcc_apply_row (1 - C) h' hij, Matrix.one_apply_ne hij]

theorem measEmp_apply_col (C : Matrix n n 𝕜) {i j : n} (h : C i i ≠ 1) (hij : i ≠ j) :
    measEmp C i j i = 0 := by
  have h' : (1 - C) i i ≠ 0 := by
    rw [one_sub_apply_self]
    exact sub_ne_zero.mpr (Ne.symm h)
  simp [measEmp, Matrix.sub_apply, measOcc_apply_col (1 - C) h' hij, Matrix.one_apply_ne (Ne.symm hij)]

/-- The particle number `trace C` is conserved by the outcome `1`. -/
theorem trace_measOcc {C : Matrix n n 𝕜} (hC : C * C = C) {i : n} (h : C i i ≠ 0) :
    trace (measOcc C i) = trace C := by
  have hX : trace (C * siteProj i * C) = C i i := by
    rw [trace_mul_cycle, hC]
    simp [siteProj, trace_mul_single]
  have hP : trace (siteProj i : Matrix n n 𝕜) = 1 := by
    simp [siteProj, trace_single_eq_same]
  rw [measOcc, trace_add, trace_sub, trace_smul, hX, hP, smul_eq_mul, inv_mul_cancel₀ h]
  ring

/-- The particle number `trace C` is conserved by the outcome `0`. -/
theorem trace_measEmp {C : Matrix n n 𝕜} (hC : C * C = C) {i : n} (h : C i i ≠ 1) :
    trace (measEmp C i) = trace C := by
  have h' : (1 - C) i i ≠ 0 := by
    rw [one_sub_apply_self]
    exact sub_ne_zero.mpr (Ne.symm h)
  rw [measEmp, trace_sub, trace_measOcc (one_sub_mul_self hC) h', trace_sub]
  ring

/-- **Born average = dephasing.** Averaging the two post-measurement correlation matrices with
the Born weights `C i i` and `1 - C i i` gives the dephased correlation matrix. -/
theorem born_average {C : Matrix n n 𝕜} {i : n} (h0 : C i i ≠ 0) (h1 : C i i ≠ 1) :
    C i i • measOcc C i + (1 - C i i) • measEmp C i =
      (1 - siteProj i) * C * (1 - siteProj i) + siteProj i * C * siteProj i := by
  have h1' : 1 - C i i ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have hPCP : siteProj i * (C * siteProj i) = C i i • siteProj i := by
    rw [← mul_assoc, siteProj_sandwich]
  simp only [measEmp, measOcc, one_sub_apply_self]
  simp only [sub_mul, mul_sub, one_mul, mul_one, mul_assoc, hPCP, smul_smul, smul_add,
    smul_sub]
  match_scalars <;> field_simp <;> ring

end Matrices

/-! ### Hermiticity and Born probabilities -/

section Star

variable {n 𝕜 : Type*} [Fintype n] [DecidableEq n] [Field 𝕜] [StarRing 𝕜]

omit [Fintype n] in
lemma siteProj_conjTranspose (i : n) : (siteProj i : Matrix n n 𝕜)ᴴ = siteProj i := by
  simp [siteProj, conjTranspose_single]

theorem measOcc_isHermitian {C : Matrix n n 𝕜} (hC : C.IsHermitian) (i : n) :
    (measOcc C i).IsHermitian := by
  have hd : star (C i i) = C i i := hC.apply i i
  unfold measOcc IsHermitian
  rw [conjTranspose_add, conjTranspose_sub, conjTranspose_smul, conjTranspose_mul,
    conjTranspose_mul, hC.eq, siteProj_conjTranspose, star_inv₀, hd, ← mul_assoc]

theorem measEmp_isHermitian {C : Matrix n n 𝕜} (hC : C.IsHermitian) (i : n) :
    (measEmp C i).IsHermitian := by
  have h1C : (1 - C).IsHermitian := isHermitian_one.sub hC
  exact isHermitian_one.sub (measOcc_isHermitian h1C i)

omit [DecidableEq n] in
/-- The diagonal entries of a Hermitian projector are sums of squared moduli. -/
lemma diag_eq_sum {C : Matrix n n 𝕜} (hH : C.IsHermitian) (hC : C * C = C) (i : n) :
    C i i = ∑ k, C i k * star (C i k) := by
  conv_lhs => rw [← hC]
  rw [Matrix.mul_apply]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [hH.apply k i]

variable [PartialOrder 𝕜] [StarOrderedRing 𝕜]

omit [DecidableEq n] in
/-- Born probability of the outcome `n_i = 1` is nonnegative. -/
theorem born_prob_nonneg {C : Matrix n n 𝕜} (hH : C.IsHermitian) (hC : C * C = C) (i : n) :
    0 ≤ C i i := by
  rw [diag_eq_sum hH hC i]
  exact Finset.sum_nonneg fun k _ => mul_star_self_nonneg (C i k)

/-- Born probability of the outcome `n_i = 1` is at most one. -/
theorem born_prob_le_one {C : Matrix n n 𝕜} (hH : C.IsHermitian) (hC : C * C = C) (i : n) :
    C i i ≤ 1 := by
  have h := born_prob_nonneg (isHermitian_one.sub hH) (one_sub_mul_self hC) i
  rw [one_sub_apply_self] at h
  exact sub_nonneg.mp h

end Star

end MonitoredFermions
