import DFL.OneBody
import DFL.FreeDyn

/-!
# Sector-averaged zero-velocity Lieb–Robinson bounds (paper: Proposition `prop:LR`)

In the dual frame the dressed fermions of the LGT are the free fermions `c_x` of each charge
sector, evolving by `c_x(t) = ∑_y (e^{-ith_q})_{xy} c_y` (`DFL.heisenberg_ann`). Under Hypothesis
(EC) their anticommutators and their commutators with local observables decay exponentially in the
distance *on average over the sectors*, uniformly in time. (No state-independent operator-norm
bound holds: the clean sector `q ≡ 0` is ballistic.)

Norms of matrices are `ℓ²` operator norms (`Matrix.Norms.L2Operator`).
-/

open Matrix Finset
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

variable {N : ℕ}

namespace LiebRobinsonAux

/-- The sector average is monotone. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (hfg : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  unfold avgQ
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => hfg q) (by positivity)

/-- The sector average commutes with a constant factor and a finite sum. -/
lemma avgQ_const_mul_sum (c : ℝ) (Y : Finset (Fin N)) (f : Cfg N → Fin N → ℝ) :
    avgQ (fun q => c * ∑ y ∈ Y, f q y) = c * ∑ y ∈ Y, avgQ (fun q => f q y) := by
  unfold avgQ
  simp_rw [Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun q _ => ?_
  ring

/-- Complexification is multiplicative. -/
lemma cplx_mul (A B : Matrix (Cfg N) (Cfg N) ℝ) : cplx (A * B) = cplx A * cplx B :=
  Matrix.map_mul (f := Complex.ofRealHom)

/-- The conjugate transpose of a complexified real matrix is its complexified transpose. -/
lemma cplx_conjTranspose_eq (A : Matrix (Cfg N) (Cfg N) ℝ) : (cplx A)ᴴ = cplx Aᵀ := by
  ext i j
  simp [cplx, conjTranspose_apply]

/-- Complexification of a diagonal matrix. -/
lemma cplx_diagonal_eq (d : Cfg N → ℝ) : cplx (diagonal d) = diagonal (fun n => (d n : ℂ)) :=
  Matrix.diagonal_map (by simp)

/-- `sgn a * sgn a = 1`. -/
lemma sgn_mul_self_eq (a : ZMod 2) : sgn a * sgn a = 1 := by
  unfold sgn; split_ifs <;> norm_num

/-- Every element of `ZMod 2` other than `1` is `0`. -/
lemma zmod2_eq_zero_of_ne_one {a : ZMod 2} (ha : a ≠ 1) : a = 0 := by
  revert a; decide

/-- Two configurations with the same value at `x` that agree after emptying `x` are equal. -/
lemma eq_of_update_zero_eq {x : Fin N} {n n' : Cfg N} (hx : n' x = n x)
    (h : Function.update n x 0 = Function.update n' x 0) : n' = n := by
  funext k
  by_cases hk : k = x
  · subst hk; exact hx
  · have := congrFun h k
    simpa [Function.update_of_ne hk] using this.symm

/-- `c_y† c_y` is the diagonal occupation operator `n_y`. -/
lemma cre_mul_ann_diag (y : Fin N) : cre y * ann y = diagonal (fun n => occ (n y)) := by
  ext n' n
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann, diagonal_apply]
  rw [Finset.sum_eq_single (Function.update n y 0)]
  · by_cases hn : n y = 1
    · by_cases hn' : n' = n
      · subst hn'
        simp [hn, sgn_mul_self_eq, occ]
      · have h' : ¬ (n' y = 1 ∧ Function.update n y 0 = Function.update n' y 0) := by
          rintro ⟨h1, h2⟩
          exact hn' (eq_of_update_zero_eq (h1.trans hn.symm) h2)
        simp only [h', hn', ite_false, zero_mul]
    · simp only [hn, false_and, ite_false, mul_zero]
      split_ifs with hn'
      · subst hn'
        simp [occ, zmod2_eq_zero_of_ne_one hn]
      · rfl
  · intro m _ hm
    have : ¬ (n y = 1 ∧ m = Function.update n y 0) := fun h => hm h.2
    simp only [this, ite_false, mul_zero]
  · intro h
    exact absurd (Finset.mem_univ _) h

end LiebRobinsonAux

/-- **Sector-averaged anticommutator bound**: `avg_q |{c_x(t), c_y†}_q| ≤ F(|x - y|)` for every
time, where `{c_x(t), c_y†}_q = (e^{-ith_q})_{xy}`. -/
theorem avg_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (LiebRobinsonAux.avgQ_mono fun q =>
    norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans (hEC.bound N x y)

/-- Uniform-in-time version, with the supremum over `t` inside the average. -/
theorem avg_iSup_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (LiebRobinsonAux.avgQ_mono fun q => ciSup_le fun t =>
    norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans (hEC.bound N x y)

/-- The Heisenberg-evolved annihilation operator in the sector `q`. -/
def annT (J h : ℝ) (q : Cfg N) (x : Fin N) (t : ℝ) : Matrix (Cfg N) (Cfg N) ℂ :=
  (evo (dGamma (hq J h q)) t)ᴴ * cplx (ann x) * evo (dGamma (hq J h q)) t

/-- Annihilation operators have operator norm at most one. -/
theorem norm_cplx_ann_le (y : Fin N) : ‖cplx (ann (N := N) y)‖ ≤ 1 := by
  -- `c_y† c_y = n_y` is a diagonal `0/1` matrix
  have h1 : (cplx (ann (N := N) y))ᴴ * cplx (ann y) =
      diagonal (fun n : Cfg N => ((occ (n y) : ℝ) : ℂ)) := by
    rw [LiebRobinsonAux.cplx_conjTranspose_eq, ← LiebRobinsonAux.cplx_mul,
      show (ann (N := N) y)ᵀ * ann y = cre y * ann y from rfl, LiebRobinsonAux.cre_mul_ann_diag,
      LiebRobinsonAux.cplx_diagonal_eq]
  have h2 : ‖(cplx (ann (N := N) y))ᴴ * cplx (ann y)‖ ≤ 1 := by
    rw [h1, l2_opNorm_diagonal]
    refine (pi_norm_le_iff_of_nonneg zero_le_one).mpr fun n => ?_
    unfold occ
    split_ifs <;> simp
  -- C*-identity `‖A‖² = ‖A† A‖`
  rw [l2_opNorm_conjTranspose_mul_self] at h2
  nlinarith [norm_nonneg (cplx (ann (N := N) y))]

/-- **Sector-averaged zero-velocity Lieb–Robinson bound**: if `B` commutes with every `c_y`,
`y ∉ Y` (e.g. an even observable supported in `Y`), then for every time
`avg_q ‖[c_x(t), B]‖ ≤ 2 ‖B‖ ∑_{y ∈ Y} F(|x - y|)`. -/
theorem avg_norm_comm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (Y : Finset (Fin N))
    (B : Matrix (Cfg N) (Cfg N) ℂ) (hB : ∀ y ∉ Y, cplx (ann y) * B = B * cplx (ann y))
    (x : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖annT J h q x t * B - B * annT J h q x t‖) ≤
      2 * ‖B‖ * ∑ y ∈ Y, F (Nat.dist x.val y.val) := by
  -- `‖[c_y, B]‖ ≤ 2 ‖B‖`
  have hcy : ∀ y : Fin N, ‖cplx (ann y) * B - B * cplx (ann y)‖ ≤ 2 * ‖B‖ := by
    intro y
    have h1 : ‖cplx (ann y) * B‖ ≤ ‖B‖ :=
      (norm_mul_le _ _).trans (mul_le_of_le_one_left (norm_nonneg B) (norm_cplx_ann_le y))
    have h2 : ‖B * cplx (ann y)‖ ≤ ‖B‖ :=
      (norm_mul_le _ _).trans (mul_le_of_le_one_right (norm_nonneg B) (norm_cplx_ann_le y))
    linarith [norm_sub_le (cplx (ann y) * B) (B * cplx (ann y))]
  -- the bound in each sector: `‖[c_x(t), B]‖ ≤ 2 ‖B‖ ∑_{y ∈ Y} |(e^{-ith_q})_{xy}|`
  have hpt : ∀ q : Cfg N, ‖annT J h q x t * B - B * annT J h q x t‖ ≤
      2 * ‖B‖ * ∑ y ∈ Y, ‖evo (hq J h q) t x y‖ := by
    intro q
    have hU : annT J h q x t = ∑ y, evo (hq J h q) t x y • cplx (ann y) :=
      heisenberg_ann (hq J h q) (hq_isSymm J h q) x t
    have hcomm : annT J h q x t * B - B * annT J h q x t =
        ∑ y ∈ Y, evo (hq J h q) t x y • (cplx (ann y) * B - B * cplx (ann y)) := by
      rw [hU, Matrix.sum_mul, Matrix.mul_sum, ← Finset.sum_sub_distrib]
      symm
      refine (Finset.sum_subset (Finset.subset_univ Y) (fun y _ hy => ?_)).trans ?_
      · rw [hB y hy, sub_self, smul_zero]
      · refine Finset.sum_congr rfl fun y _ => ?_
        rw [smul_mul_assoc, mul_smul_comm, smul_sub]
    rw [hcomm, Finset.mul_sum]
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun y _ => ?_)
    refine (norm_smul_le _ _).trans ?_
    rw [mul_comm (2 * ‖B‖)]
    exact mul_le_mul_of_nonneg_left (hcy y) (norm_nonneg _)
  -- average over the sectors
  calc avgQ (fun q : Cfg N => ‖annT J h q x t * B - B * annT J h q x t‖)
      ≤ avgQ (fun q : Cfg N => 2 * ‖B‖ * ∑ y ∈ Y, ‖evo (hq J h q) t x y‖) :=
        LiebRobinsonAux.avgQ_mono hpt
    _ = 2 * ‖B‖ * ∑ y ∈ Y, avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖) :=
        LiebRobinsonAux.avgQ_const_mul_sum _ _ _
    _ ≤ 2 * ‖B‖ * ∑ y ∈ Y, F (Nat.dist x.val y.val) := by
        refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun y _ => ?_) (by positivity)
        exact avg_norm_evo_le J h hEC x y t

end DFL
