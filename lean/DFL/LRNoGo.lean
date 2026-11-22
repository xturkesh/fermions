import DFL.LRNorm
import DFL.DimerTL
import DFL.InfVolNumber

/-!
# No operator-norm zero-velocity bound (paper: Proposition `prop:LR`(d))

There is no function `G ≥ 0` with `∑_d d² G(d)² < ∞` such that `sup_t ‖[n_x(t), n_y]‖ ≤ G(|x - y|)` for all `N` and
all `x, y`; the same holds with `[n_x(t), n_y]` replaced by `{c̃_x(t), c̃_y†}` or by `[c̃_x(t), n_y]`
(`no_uniform_LR_nn`, `no_uniform_LR_anticomm`, `no_uniform_LR_cn`).

The operator norms are maxima over the sectors (`norm_comm_nn_eq`, `norm_dressedAnticomm_eq`, `norm_comm_number_eq`) and
include the clean sector, where the dressed fermion is a free particle on the chain. It moves ballistically:
`∑_r r² |U(t)_{r, c}|² ≥ (√2 |J t| - ε)²` in the middle `c` of a long chain (Heisenberg equations `Ẋ = V`, `V̇ = i[H, V]`,
with `[H, V]` supported at the two ends, which the particle reaches only with exponentially small amplitude), which is
incompatible with `|U(t)_{rc}| ≤ √2 G(|r - c|)` for large `|r - c|`.

## Proof (finite-volume version of the paper's argument; auxiliary lemmas in `DFL.LRNoGoAux`)

* `no_clean_bound`: in the clean sector `q = 0` no profile `F` with `∑_d d² F(d)² < ∞` bounds
  `|U_0(t)_{xy}|` for `|x - y| ≥ d₀`. On `[0, 2L]` the second moment about the centre `L` is then bounded
  uniformly in `L` and `t` (`sum_dist_le`), while `ballistic_lower` gives
  `|t| |J| ≤ (∑_r (r - L)² |U_0(t)_{rL}|²)^{1/2} + 3 J² w_{K|t|}(L) t²` with the propagation weight
  `w_c(L) = c^L/L! e^c → 0` (`InfVolNumberAux.norm_evo_hq_le`): the double commutator `[h_0, [h_0, X]]` is
  diagonal and vanishes in the bulk (`comm_comm_eq_diagonal`), and `taylor_two` is the second-order mean
  value inequality.
* For `[n_x(t), n_y]` the bound `g(|U|) ≤ G` with `g(s) = s (1 - s²)^{1/2}` gives `|U| ≤ √2 G` once
  `G < 1/2`, by the intermediate value theorem (`|U_0(0)_{xy}| = 0` and `g(2^{-1/2}) = 1/2`).
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator

noncomputable section

namespace DFL

namespace LRNoGoAux

open DimerTLAux

variable {N : ℕ}

/-- The superdiagonal shift `S_{ij} = [j = i + 1]`. -/
def shiftS (N : ℕ) : Matrix (Fin N) (Fin N) ℝ := fun i j => if j.val = i.val + 1 then 1 else 0

/-- The position operator `X = diag(r - c)`. -/
def posC (N c : ℕ) : Matrix (Fin N) (Fin N) ℝ := diagonal (fun r => (r.val : ℝ) - c)

/-- `Δ = S + Sᵀ`. -/
lemma pathAdj_eq_shift : pathAdj N = shiftS N + (shiftS N)ᵀ := by
  ext i j
  simp only [pathAdj, shiftS, Matrix.add_apply, Matrix.transpose_apply]
  by_cases h1 : j.val = i.val + 1
  · rw [ite_eq_left (Or.inl h1.symm), ite_eq_left h1, ite_eq_right (by omega), add_zero]
  · by_cases h2 : i.val = j.val + 1
    · rw [ite_eq_left (Or.inr h2.symm), ite_eq_right h1, ite_eq_left h2, zero_add]
    · rw [ite_eq_right (by omega), ite_eq_right h1, ite_eq_right h2, add_zero]

/-- `[S, X] = S`. -/
lemma shiftS_comm_posC (c : ℕ) : shiftS N * posC N c - posC N c * shiftS N = shiftS N := by
  ext i j
  simp only [posC, Matrix.sub_apply, Matrix.mul_diagonal, Matrix.diagonal_mul, shiftS]
  split_ifs with h
  · have : (j.val : ℝ) = i.val + 1 := by exact_mod_cast h
    rw [this]; ring
  · ring

/-- `[Sᵀ, X] = -Sᵀ`. -/
lemma shiftST_comm_posC (c : ℕ) :
    (shiftS N)ᵀ * posC N c - posC N c * (shiftS N)ᵀ = -(shiftS N)ᵀ := by
  ext i j
  simp only [posC, Matrix.sub_apply, Matrix.neg_apply, Matrix.mul_diagonal, Matrix.diagonal_mul,
    Matrix.transpose_apply, shiftS]
  split_ifs with h
  · have : (i.val : ℝ) = j.val + 1 := by exact_mod_cast h
    rw [this]; ring
  · ring

/-- `SᵀS = diag([1 ≤ i])`. -/
lemma shiftST_mul_shiftS :
    (shiftS N)ᵀ * shiftS N = diagonal (fun i : Fin N => if 1 ≤ i.val then (1 : ℝ) else 0) := by
  ext i k
  rw [Matrix.mul_apply, Matrix.diagonal_apply]
  simp only [Matrix.transpose_apply, shiftS]
  by_cases hik : i = k
  · subst hik
    rw [ite_eq_left rfl]
    by_cases hi : 1 ≤ i.val
    · rw [ite_eq_left hi, Finset.sum_eq_single ⟨i.val - 1, by omega⟩]
      · rw [ite_eq_left (by simp only; omega), one_mul]
      · intro j _ hj
        rw [ite_eq_right, zero_mul]
        intro h
        exact hj (Fin.ext (by simp only; omega))
      · intro h
        exact absurd (Finset.mem_univ _) h
    · rw [ite_eq_right hi]
      refine Finset.sum_eq_zero fun j _ => ?_
      rw [ite_eq_right (by omega), zero_mul]
  · rw [ite_eq_right hik]
    refine Finset.sum_eq_zero fun j _ => ?_
    by_cases h1 : i.val = j.val + 1
    · rw [ite_eq_left h1, ite_eq_right (fun h2 => hik (Fin.ext (by omega))), mul_zero]
    · rw [ite_eq_right h1, zero_mul]

/-- `SSᵀ = diag([i + 1 < N])`. -/
lemma shiftS_mul_shiftST :
    shiftS N * (shiftS N)ᵀ = diagonal (fun i : Fin N => if i.val + 1 < N then (1 : ℝ) else 0) := by
  ext i k
  rw [Matrix.mul_apply, Matrix.diagonal_apply]
  simp only [Matrix.transpose_apply, shiftS]
  by_cases hik : i = k
  · subst hik
    rw [ite_eq_left rfl]
    by_cases hi : i.val + 1 < N
    · rw [ite_eq_left hi, Finset.sum_eq_single ⟨i.val + 1, hi⟩]
      · simp
      · intro j _ hj
        rw [ite_eq_right (fun h => hj (Fin.ext h)), zero_mul]
      · intro h
        exact absurd (Finset.mem_univ _) h
    · rw [ite_eq_right hi]
      refine Finset.sum_eq_zero fun j _ => ?_
      rw [ite_eq_right (by have := j.isLt; omega), zero_mul]
  · rw [ite_eq_right hik]
    refine Finset.sum_eq_zero fun j _ => ?_
    by_cases h1 : j.val = i.val + 1
    · rw [ite_eq_left h1, ite_eq_right (fun h2 => hik (Fin.ext (by omega))), mul_zero]
    · rw [ite_eq_right h1, zero_mul]

/-- The clean-sector Hamiltonian `h_0 = -J Δ + 2h`. -/
lemma hq_zero_eq (J h : ℝ) :
    hq J h (0 : Cfg N) = (-J) • pathAdj N + (2 * h) • (1 : Matrix (Fin N) (Fin N) ℝ) := by
  have e : diagonal (fun x : Fin N => sgn ((0 : Cfg N) x)) = 1 := by
    simp [sgn]
  rw [hq, e]

/-- `[h_0, X] = -J (S - Sᵀ)`. -/
lemma comm_hq_posC (J h : ℝ) (c : ℕ) :
    hq J h (0 : Cfg N) * posC N c - posC N c * hq J h 0 = (-J) • (shiftS N - (shiftS N)ᵀ) := by
  have h1 := shiftS_comm_posC (N := N) c
  have h2 := shiftST_comm_posC (N := N) c
  rw [hq_zero_eq, pathAdj_eq_shift]
  calc ((-J) • (shiftS N + (shiftS N)ᵀ) + (2 * h) • (1 : Matrix (Fin N) (Fin N) ℝ)) * posC N c -
        posC N c * ((-J) • (shiftS N + (shiftS N)ᵀ) + (2 * h) • 1)
      = (-J) • ((shiftS N * posC N c - posC N c * shiftS N) +
          ((shiftS N)ᵀ * posC N c - posC N c * (shiftS N)ᵀ)) := by
        simp only [Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.mul_smul, Matrix.one_mul,
          Matrix.mul_one, smul_add, smul_sub]
        abel
    _ = (-J) • (shiftS N - (shiftS N)ᵀ) := by rw [h1, h2, sub_eq_add_neg]

/-- `[h_0, [h_0, X]] = 2J² (SᵀS - SSᵀ)`. -/
lemma comm_hq_comm_hq_posC (J h : ℝ) (c : ℕ) :
    hq J h (0 : Cfg N) * (hq J h 0 * posC N c - posC N c * hq J h 0) -
      (hq J h 0 * posC N c - posC N c * hq J h 0) * hq J h 0 =
      (2 * J ^ 2) • ((shiftS N)ᵀ * shiftS N - shiftS N * (shiftS N)ᵀ) := by
  rw [comm_hq_posC, hq_zero_eq, pathAdj_eq_shift]
  simp only [Matrix.add_mul, Matrix.mul_add, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul,
    Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, smul_add, smul_sub, smul_smul]
  module

/-- The double commutator is diagonal and vanishes in the bulk. -/
lemma comm_comm_eq_diagonal (J h : ℝ) (c : ℕ) :
    hq J h (0 : Cfg N) * (hq J h 0 * posC N c - posC N c * hq J h 0) -
      (hq J h 0 * posC N c - posC N c * hq J h 0) * hq J h 0 =
      diagonal (fun i : Fin N => 2 * J ^ 2 *
        ((if 1 ≤ i.val then (1 : ℝ) else 0) - (if i.val + 1 < N then 1 else 0))) := by
  rw [comm_hq_comm_hq_posC, shiftST_mul_shiftS, shiftS_mul_shiftST, Matrix.diagonal_sub,
    ← Matrix.diagonal_smul]
  rfl

/-! ### Analytic helpers -/

/-- `cplx` respects subtraction. -/
lemma cplx_sub_eq (A B : Matrix (Fin N) (Fin N) ℝ) : cplx (A - B) = cplx A - cplx B := by
  ext i j
  simp [cplx]

/-- `cplx` respects products. -/
lemma cplx_mul_eq (A B : Matrix (Fin N) (Fin N) ℝ) : cplx (A * B) = cplx A * cplx B := by
  ext i j
  simp [cplx, Matrix.mul_apply]

/-- Column extraction commutes with complex scalars. -/
lemma colL_smul (z : Fin N) (a : ℂ) (M : Matrix (Fin N) (Fin N) ℂ) :
    colL z (a • M) = a • colL z M := by
  rw [colL_apply, colL_apply, ← WithLp.toLp_smul]
  rfl

/-- The squared norm of a column of `U(-r) M U(r)`. -/
lemma norm_colL_heis_sq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (M : Matrix (Fin N) (Fin N) ℂ) (r : ℝ) (z : Fin N) :
    ‖colL z (evo A (-r) * M * evo A r)‖ ^ 2 = ∑ x, ‖(M *ᵥ fun y => evo A r y z) x‖ ^ 2 := by
  rw [colL_apply, Matrix.mul_assoc, col_mul, norm_toLp_mulVec _ (evo_unitary A hA (-r)), col_mul,
    EuclideanSpace.norm_sq_eq]

/-- **Second-order Taylor estimate** from a bound on the second derivative. -/
lemma taylor_two {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] (f f' f'' : ℝ → E)
    (hf : ∀ r, HasDerivAt f (f' r) r) (hf' : ∀ r, HasDerivAt f' (f'' r) r) (M t : ℝ)
    (hM : ∀ r, |r| ≤ |t| → ‖f'' r‖ ≤ M) :
    ‖f t - f 0 - t • f' 0‖ ≤ M * |t| * |t| := by
  have hT0 : 0 ≤ |t| := abs_nonneg t
  have hM0 : 0 ≤ M := (norm_nonneg _).trans (hM 0 (by simp))
  have hs : Convex ℝ (Set.Icc (-|t|) |t|) := convex_Icc _ _
  have h0 : (0 : ℝ) ∈ Set.Icc (-|t|) |t| := ⟨by linarith, hT0⟩
  have ht : t ∈ Set.Icc (-|t|) |t| := ⟨neg_abs_le t, le_abs_self t⟩
  have ha : ∀ r ∈ Set.Icc (-|t|) |t|, ‖f' r - f' 0‖ ≤ M * |t| := by
    intro r hr
    have key := hs.norm_image_sub_le_of_norm_hasDerivWithin_le (f := f') (f' := f'')
      (fun x _ => (hf' x).hasDerivWithinAt) (fun x hx => hM x (abs_le.mpr hx)) h0 hr
    rw [sub_zero, Real.norm_eq_abs] at key
    exact key.trans (mul_le_mul_of_nonneg_left (abs_le.mpr hr) hM0)
  have hg : ∀ r, HasDerivAt (fun r => f r - r • f' 0) (f' r - f' 0) r := by
    intro r
    have := (hf r).fun_sub ((hasDerivAt_id r).smul_const (f' 0))
    simpa using this
  have key := hs.norm_image_sub_le_of_norm_hasDerivWithin_le (f := fun r => f r - r • f' 0)
    (f' := fun r => f' r - f' 0) (fun x _ => (hg x).hasDerivWithinAt) ha h0 ht
  simp only [zero_smul, sub_zero, Real.norm_eq_abs] at key
  calc ‖f t - f 0 - t • f' 0‖ = ‖f t - t • f' 0 - f 0‖ := by congr 1; abel
    _ ≤ M * |t| * |t| := key

/-- `|a - b|²` as the square of `Nat.dist`. -/
lemma cast_dist_sq (a b : ℕ) : ((Nat.dist a b : ℕ) : ℝ) ^ 2 = ((a : ℝ) - b) ^ 2 := by
  rcases le_total a b with hab | hab
  · rw [Nat.dist_eq_sub_of_le hab, Nat.cast_sub hab]
    ring
  · rw [Nat.dist_eq_sub_of_le_right hab, Nat.cast_sub hab]

/-- Each distance is attained by at most two sites: `∑_r g(|r - c|) ≤ 2 ∑_d g(d)`. -/
lemma sum_dist_le (c : ℕ) (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) (hgs : Summable g) :
    ∑ r : Fin N, g (Nat.dist r.val c) ≤ 2 * ∑' d, g d := by
  have key : ∀ (s : Finset (Fin N)) (e : Fin N → ℕ), Set.InjOn e s →
      ∑ r ∈ s, g (e r) ≤ ∑' d, g d := by
    intro s e he
    rw [← Finset.sum_image he]
    exact hgs.sum_le_tsum _ (fun d _ => hg d)
  rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun r : Fin N => c ≤ r.val), two_mul]
  refine add_le_add ?_ ?_
  · calc ∑ r ∈ Finset.univ.filter (fun r : Fin N => c ≤ r.val), g (Nat.dist r.val c)
        = ∑ r ∈ Finset.univ.filter (fun r : Fin N => c ≤ r.val), g (r.val - c) := by
          refine Finset.sum_congr rfl fun r hr => ?_
          rw [Finset.mem_filter] at hr
          rw [Nat.dist_eq_sub_of_le_right hr.2]
      _ ≤ ∑' d, g d := by
          refine key _ _ ?_
          intro a ha b hb hab
          simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb hab
          exact Fin.ext (by omega)
  · calc ∑ r ∈ Finset.univ.filter (fun r : Fin N => ¬ c ≤ r.val), g (Nat.dist r.val c)
        = ∑ r ∈ Finset.univ.filter (fun r : Fin N => ¬ c ≤ r.val), g (c - r.val) := by
          refine Finset.sum_congr rfl fun r hr => ?_
          rw [Finset.mem_filter] at hr
          rw [Nat.dist_eq_sub_of_le (by omega)]
      _ ≤ ∑' d, g d := by
          refine key _ _ ?_
          intro a ha b hb hab
          simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb hab
          exact Fin.ext (by omega)

/-! ### The ballistic lower bound in the clean sector -/

/-- **Ballistic lower bound** in the middle `c = L` of the clean chain `[0, 2L]`: with `K = 2|J| + 2|h|`,
`|t| |J| ≤ (∑_r (r - L)² |U_0(t)_{rL}|²)^{1/2} + 3 J² w_{K|t|}(L) t²`. The Heisenberg position
`X(t) δ_L = U(t)^* X U(t) δ_L` has `X(0) δ_L = 0`, velocity `i[h_0, X] δ_L` of norm `≥ |J|` at `t = 0`, and
acceleration `-U(t)^* [h_0, [h_0, X]] U(t) δ_L`, which only sees the two ends of the chain. -/
lemma ballistic_lower (J h : ℝ) (L : ℕ) (hL : 1 ≤ L) (t : ℝ) :
    |t| * |J| ≤ √(∑ r : Fin (2 * L + 1), ((r.val : ℝ) - L) ^ 2 *
        ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ^ 2) +
      3 * J ^ 2 * InfVolNumberAux.wt (|t| * (2 * |J| + 2 * |h|)) L * |t| * |t| := by
  set c : Fin (2 * L + 1) := ⟨L, by omega⟩ with hc
  set c₁ : Fin (2 * L + 1) := ⟨L + 1, by omega⟩ with hc₁
  set A := hq J h (0 : Cfg (2 * L + 1)) with hA
  have hAs : A.IsSymm := hq_isSymm J h 0
  set X := posC (2 * L + 1) L with hX
  set D := A * X - X * A with hD
  set β : Fin (2 * L + 1) → ℝ := fun i => 2 * J ^ 2 *
    ((if 1 ≤ i.val then (1 : ℝ) else 0) - (if i.val + 1 < 2 * L + 1 then 1 else 0)) with hβ
  have hQ : A * D - D * A = diagonal β := comm_comm_eq_diagonal J h L
  set w := InfVolNumberAux.wt (|t| * (2 * |J| + 2 * |h|)) L with hw
  have hw0 : 0 ≤ w := InfVolNumberAux.wt_nonneg (by positivity) L
  set f : ℝ → EuclideanSpace ℂ (Fin (2 * L + 1)) :=
    fun r => colL c (evo A (-r) * cplx X * evo A r) with hf
  set f' : ℝ → EuclideanSpace ℂ (Fin (2 * L + 1)) :=
    fun r => colL c (evo A (-r) * (Complex.I • (cplx A * cplx X - cplx X * cplx A)) * evo A r)
    with hf'
  set f'' : ℝ → EuclideanSpace ℂ (Fin (2 * L + 1)) :=
    fun r => Complex.I • colL c (evo A (-r) * (Complex.I • (cplx A * cplx D - cplx D * cplx A)) *
      evo A r) with hf''
  have hAX : cplx A * cplx X - cplx X * cplx A = cplx D := by
    rw [hD, cplx_sub_eq, cplx_mul_eq, cplx_mul_eq]
  have hAD : cplx A * cplx D - cplx D * cplx A = cplx (diagonal β) := by
    rw [← hQ, cplx_sub_eq (A * D) (D * A), cplx_mul_eq A D, cplx_mul_eq D A]
  have hf'eq : f' = fun r => Complex.I • colL c (evo A (-r) * cplx D * evo A r) := by
    funext r
    simp only [hf']
    rw [hAX, Matrix.mul_smul, Matrix.smul_mul, colL_smul]
  have hd1 : ∀ r, HasDerivAt f (f' r) r := fun r => hasDerivAt_heis_col A X r c
  have hd2 : ∀ r, HasDerivAt f' (f'' r) r := by
    intro r
    rw [hf'eq]
    exact (hasDerivAt_heis_col A D r c).const_smul Complex.I
  -- the acceleration only sees the two ends of the chain
  have hM : ∀ r, |r| ≤ |t| → ‖f'' r‖ ≤ 3 * J ^ 2 * w := by
    intro r hr
    have hwr : ∀ x : Fin (2 * L + 1), Nat.dist x.val L = L → ‖evo A r x c‖ ≤ w := by
      intro x hx
      have h1 := InfVolNumberAux.norm_evo_hq_le J h r 0 x c
      rw [show (c : Fin (2 * L + 1)).val = L from rfl, hx] at h1
      refine h1.trans ?_
      rw [hw]
      unfold InfVolNumberAux.wt
      gcongr
    have hterm : ∀ x : Fin (2 * L + 1),
        ‖((Complex.I • cplx (diagonal β)) *ᵥ fun y => evo A r y c) x‖ ^ 2 ≤
          4 * J ^ 4 * w ^ 2 * ((if x.val = 0 then 1 else 0) + (if x.val = 2 * L then 1 else 0)) := by
      intro x
      have e : ‖((Complex.I • cplx (diagonal β)) *ᵥ fun y => evo A r y c) x‖ =
          |β x| * ‖evo A r x c‖ := by
        rw [Matrix.smul_mulVec, Pi.smul_apply, cplx_diagonal, Matrix.mulVec_diagonal, norm_smul,
          Complex.norm_I, one_mul, norm_mul, Complex.norm_real, Real.norm_eq_abs]
      rw [e, mul_pow, sq_abs]
      have hJ4 : 0 ≤ 4 * J ^ 4 * w ^ 2 := by positivity
      by_cases h1 : 1 ≤ x.val
      · by_cases h2 : x.val + 1 < 2 * L + 1
        · have hb : β x = 0 := by
            simp only [hβ]
            rw [ite_eq_left h1, ite_eq_left h2]
            ring
          rw [hb]
          have : (0 : ℝ) ≤ (if x.val = 0 then 1 else 0) + (if x.val = 2 * L then 1 else 0) := by
            positivity
          nlinarith
        · have hx : x.val = 2 * L := by omega
          have hb : β x ^ 2 = 4 * J ^ 4 := by
            simp only [hβ]
            rw [ite_eq_left h1, ite_eq_right h2]
            ring
          have hu := hwr x (by rw [hx, Nat.dist_eq_sub_of_le_right (by omega)]; omega)
          have hu2 : ‖evo A r x c‖ ^ 2 ≤ w ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hu 2
          rw [hb, ite_eq_right (by omega), ite_eq_left hx]
          nlinarith
      · have hx : x.val = 0 := by omega
        have hb : β x ^ 2 = 4 * J ^ 4 := by
          simp only [hβ]
          rw [ite_eq_right h1, ite_eq_left (by omega)]
          ring
        have hu := hwr x (by rw [hx, Nat.dist_zero_left])
        have hu2 : ‖evo A r x c‖ ^ 2 ≤ w ^ 2 := pow_le_pow_left₀ (norm_nonneg _) hu 2
        rw [hb, ite_eq_left hx, ite_eq_right (by omega)]
        nlinarith
    have hsq : ‖f'' r‖ ^ 2 ≤ (3 * J ^ 2 * w) ^ 2 := by
      simp only [hf'']
      rw [norm_smul, Complex.norm_I, one_mul, hAD, norm_colL_heis_sq A hAs]
      calc ∑ x, ‖((Complex.I • cplx (diagonal β)) *ᵥ fun y => evo A r y c) x‖ ^ 2
          ≤ ∑ x : Fin (2 * L + 1), 4 * J ^ 4 * w ^ 2 *
              ((if x.val = 0 then 1 else 0) + (if x.val = 2 * L then 1 else 0)) :=
            Finset.sum_le_sum fun x _ => hterm x
        _ ≤ 4 * J ^ 4 * w ^ 2 * 2 := by
            rw [← Finset.mul_sum, Finset.sum_add_distrib]
            have e1 := InfVolAux.sum_ite_val_le_one (N := 2 * L + 1) 0
            have e2 := InfVolAux.sum_ite_val_le_one (N := 2 * L + 1) (2 * L)
            have : 0 ≤ 4 * J ^ 4 * w ^ 2 := by positivity
            nlinarith
        _ ≤ (3 * J ^ 2 * w) ^ 2 := by
            have : 0 ≤ J ^ 4 * w ^ 2 := by positivity
            nlinarith
    exact (pow_le_pow_iff_left₀ (norm_nonneg _) (by positivity) two_ne_zero).mp hsq
  have htaylor := taylor_two f f' f'' hd1 hd2 (3 * J ^ 2 * w) t hM
  -- the initial position and velocity
  have hf0 : f 0 = 0 := by
    simp only [hf, neg_zero, evo_zero, Matrix.one_mul, Matrix.mul_one]
    rw [colL_apply]
    have : (fun x => cplx X x c) = 0 := by
      funext x
      simp only [hX, posC, cplx, Matrix.map_apply, Matrix.diagonal_apply, Pi.zero_apply]
      split_ifs with hxc
      · rw [hxc]
        simp [hc]
      · simp
    rw [this, WithLp.toLp_zero]
  have hDc : D c₁ c = J := by
    rw [hD, hA, hX, comm_hq_posC, hc, hc₁]
    simp only [Matrix.smul_apply, Matrix.sub_apply, Matrix.transpose_apply, shiftS, smul_eq_mul]
    split_ifs <;> first | omega | ring
  have hf'0 : |J| ≤ ‖f' 0‖ := by
    rw [hf'eq]
    simp only [neg_zero, evo_zero, Matrix.one_mul, Matrix.mul_one]
    rw [norm_smul, Complex.norm_I, one_mul]
    refine le_trans ?_ (PiLp.norm_apply_le _ c₁)
    rw [colL_apply]
    simp [cplx, hDc]
  have hft : ‖f t‖ = √(∑ r : Fin (2 * L + 1), ((r.val : ℝ) - L) ^ 2 * ‖evo A t r c‖ ^ 2) := by
    rw [← Real.sqrt_sq (norm_nonneg (f t))]
    congr 1
    simp only [hf]
    rw [norm_colL_heis_sq A hAs]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [hX, posC, cplx_diagonal, Matrix.mulVec_diagonal, norm_mul, mul_pow, Complex.norm_real,
      Real.norm_eq_abs, sq_abs]
  rw [hf0, sub_zero] at htaylor
  have h1 : ‖t • f' 0‖ ≤ ‖f t‖ + ‖f t - t • f' 0‖ := by
    calc ‖t • f' 0‖ = ‖f t - (f t - t • f' 0)‖ := by rw [sub_sub_cancel]
      _ ≤ ‖f t‖ + ‖f t - t • f' 0‖ := norm_sub_le _ _
  rw [norm_smul, Real.norm_eq_abs] at h1
  rw [← hft]
  calc |t| * |J| ≤ |t| * ‖f' 0‖ := mul_le_mul_of_nonneg_left hf'0 (abs_nonneg t)
    _ ≤ ‖f t‖ + 3 * J ^ 2 * w * |t| * |t| := by linarith

/-- **The clean sector is ballistic**: no profile `F` with `∑_d d² F(d)² < ∞` bounds `|U_0(t)_{xy}|` for all
`N`, all `t` and all `x, y` at distance at least `d₀`. -/
lemma no_clean_bound (J h : ℝ) (hJ : J ≠ 0) (F : ℕ → ℝ)
    (hF : Summable (fun d : ℕ => (d : ℝ) ^ 2 * F d ^ 2)) (d₀ : ℕ)
    (hb : ∀ (N : ℕ) (x y : Fin N) (t : ℝ), d₀ ≤ Nat.dist x.val y.val →
      ‖evo (hq J h (0 : Cfg N)) t x y‖ ≤ F (Nat.dist x.val y.val)) : False := by
  set g : ℕ → ℝ := fun d => (d : ℝ) ^ 2 * F d ^ 2 + if d < d₀ then (d : ℝ) ^ 2 else 0 with hg
  have hg0 : ∀ d, 0 ≤ g d := fun d => by
    simp only [hg]
    have : 0 ≤ (d : ℝ) ^ 2 * F d ^ 2 := by positivity
    split_ifs <;> positivity
  have hgs : Summable g := by
    refine hF.add (summable_of_ne_finset_zero (s := Finset.range d₀) fun d hd => ?_)
    rw [Finset.mem_range] at hd
    rw [ite_eq_right hd]
  set S := 2 * ∑' d, g d with hS
  have hS0 : 0 ≤ S := by
    have : 0 ≤ ∑' d, g d := tsum_nonneg hg0
    linarith
  -- the second moment is bounded uniformly in `L` and `t`
  have hup : ∀ (L : ℕ) (t : ℝ), ∑ r : Fin (2 * L + 1), ((r.val : ℝ) - L) ^ 2 *
      ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ^ 2 ≤ S := by
    intro L t
    refine le_trans (Finset.sum_le_sum fun r _ => ?_) (sum_dist_le L g hg0 hgs)
    rw [← cast_dist_sq]
    have hU1 : ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ≤ 1 :=
      norm_evo_le_one _ (hq_isSymm J h 0) t _ _
    have hd2 : (0 : ℝ) ≤ ((Nat.dist r.val L : ℕ) : ℝ) ^ 2 := by positivity
    simp only [hg]
    by_cases hdd : d₀ ≤ Nat.dist r.val L
    · have hU : ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ≤ F (Nat.dist r.val L) :=
        hb (2 * L + 1) r ⟨L, by omega⟩ t hdd
      rw [ite_eq_right (by omega)]
      have : ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ^ 2 ≤ F (Nat.dist r.val L) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) hU 2
      nlinarith
    · rw [ite_eq_left (by omega)]
      have : ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ^ 2 ≤ 1 := by
        have := pow_le_pow_left₀ (norm_nonneg _) hU1 2
        rwa [one_pow] at this
      have : 0 ≤ ((Nat.dist r.val L : ℕ) : ℝ) ^ 2 * F (Nat.dist r.val L) ^ 2 := by positivity
      nlinarith
  -- choose a large time, then a long chain
  have hJ0 : 0 < |J| := abs_pos.mpr hJ
  set K0 := S + 1 with hK0
  have hK0pos : 0 < K0 := by linarith
  set t : ℝ := 2 * K0 / |J| with ht
  have ht0 : 0 < t := by positivity
  have htJ : |t| * |J| = 2 * K0 := by
    rw [abs_of_pos ht0, ht]
    field_simp
  have hlim : Filter.Tendsto
      (fun L : ℕ => 3 * J ^ 2 * InfVolNumberAux.wt (|t| * (2 * |J| + 2 * |h|)) L * |t| * |t|)
      Filter.atTop (nhds 0) := by
    have h1 := ((((FloorSemiring.tendsto_pow_div_factorial_atTop (|t| * (2 * |J| + 2 * |h|))).mul_const
      (Real.exp (|t| * (2 * |J| + 2 * |h|)))).const_mul (3 * J ^ 2)).mul_const |t|).mul_const |t|
    rw [zero_mul, mul_zero, zero_mul, zero_mul] at h1
    exact h1
  obtain ⟨L, hLsmall, hL1⟩ :=
    ((hlim.eventually (gt_mem_nhds hK0pos)).and (Filter.eventually_ge_atTop 1)).exists
  have hlow := ballistic_lower J h L hL1 t
  have hup' := hup L t
  set m := ∑ r : Fin (2 * L + 1), ((r.val : ℝ) - L) ^ 2 *
    ‖evo (hq J h (0 : Cfg (2 * L + 1))) t r ⟨L, by omega⟩‖ ^ 2 with hm
  have hm0 : 0 ≤ m := Finset.sum_nonneg fun r _ => by positivity
  have hK0le : K0 ≤ √m := by
    rw [htJ] at hlow
    linarith
  have hK0sq : K0 ^ 2 ≤ m := by
    have := pow_le_pow_left₀ hK0pos.le hK0le 2
    rwa [Real.sq_sqrt hm0] at this
  nlinarith

end LRNoGoAux

open LRNoGoAux

/-- **No operator-norm zero-velocity bound for `[n_x(t), n_y]`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_nn (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ), ‖numComm (N := N) J h x y t‖ ≤ G (Nat.dist x.val y.val) := by
  have _ := hG0  -- the nonnegativity of `G` is not needed (it follows from the bound)
  intro hbound
  -- `G(d) < 1/2` for large `d`
  obtain ⟨d₁, hd₁⟩ := Filter.eventually_atTop.mp
    (hG.tendsto_atTop_zero.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 4)))
  have hGsmall : ∀ d, max d₁ 1 ≤ d → G d < 1 / 2 := by
    intro d hd
    have h1 : (d : ℝ) ^ 2 * G d ^ 2 < 1 / 4 := hd₁ d (le_trans (le_max_left _ _) hd)
    have hd1 : (1 : ℝ) ≤ d := by exact_mod_cast le_trans (le_max_right _ _) hd
    have h2 : G d ^ 2 < 1 / 4 := by nlinarith [sq_nonneg (G d)]
    nlinarith
  refine no_clean_bound J h hJ (fun d => √2 * G d) ?_ (max d₁ 1) fun N x y t hd => ?_
  · have e : (fun d : ℕ => (d : ℝ) ^ 2 * (√2 * G d) ^ 2) =
        fun d : ℕ => 2 * ((d : ℝ) ^ 2 * G d ^ 2) := by
      funext d
      rw [mul_pow, Real.sq_sqrt (by norm_num)]
      ring
    rw [e]
    exact hG.mul_left 2
  · -- continuity argument: `|U_0(t)_{xy}|` never reaches `2^{-1/2}`
    set d := Nat.dist x.val y.val with hd'
    set u : ℝ → ℝ := fun s => ‖evo (hq J h (0 : Cfg N)) s x y‖ with hu
    have hucont : Continuous u :=
      (LRNormAux.continuous_evo_apply _ (hq_isHermitian J h 0) x y).norm
    have hxy : x ≠ y := by
      intro e
      rw [e, Nat.dist_self] at hd'
      have := le_max_right d₁ 1
      omega
    have hu0 : u 0 = 0 := by
      simp only [hu]
      rw [LRNormAux.evo_zero_eq, Matrix.one_apply_ne hxy, norm_zero]
    have hg : ∀ s, u s * √(1 - u s ^ 2) ≤ G d := by
      intro s
      refine le_trans ?_ (hbound N x y s)
      rw [norm_comm_nn_eq]
      exact Finset.le_sup' (fun q : Cfg N =>
        ‖evo (hq J h q) s x y‖ * √(1 - ‖evo (hq J h q) s x y‖ ^ 2)) (Finset.mem_univ 0)
    have hGd : G d < 1 / 2 := hGsmall d hd
    have hsmall : ∀ s, u s ^ 2 < 1 / 2 := by
      intro s
      by_contra hcon
      have hs0 : √(1 / 2) ≤ u s := by
        calc √(1 / 2) ≤ √(u s ^ 2) := Real.sqrt_le_sqrt (not_lt.mp hcon)
          _ = u s := Real.sqrt_sq (norm_nonneg _)
      have hmem : √(1 / 2) ∈ Set.uIcc (u 0) (u s) := by
        rw [hu0]
        exact Set.mem_uIcc.mpr (Or.inl ⟨Real.sqrt_nonneg _, hs0⟩)
      obtain ⟨s', _, hs'⟩ := intermediate_value_uIcc hucont.continuousOn hmem
      have h1 := hg s'
      rw [hs', Real.sq_sqrt (by norm_num), show (1 : ℝ) - 1 / 2 = 1 / 2 by norm_num,
        Real.mul_self_sqrt (by norm_num)] at h1
      linarith
    have ht := hsmall t
    have hut0 : 0 ≤ u t := norm_nonneg _
    have h1 : √(1 / 2) ≤ √(1 - u t ^ 2) := Real.sqrt_le_sqrt (by linarith)
    have h2 : u t * √(1 / 2) ≤ G d := le_trans (mul_le_mul_of_nonneg_left h1 hut0) (hg t)
    have h3 : √(1 / 2) * √2 = 1 := by
      rw [← Real.sqrt_mul (by norm_num)]
      norm_num
    show u t ≤ √2 * G d
    calc u t = u t * √(1 / 2) * √2 := by rw [mul_assoc, h3, mul_one]
      _ ≤ G d * √2 := mul_le_mul_of_nonneg_right h2 (Real.sqrt_nonneg 2)
      _ = √2 * G d := mul_comm _ _

/-- **No operator-norm zero-velocity bound for `{c̃_x(t), c̃_y†}`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_anticomm (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ),
      ‖dressedAnticomm (N := N) J h x y t‖ ≤ G (Nat.dist x.val y.val) := by
  have _ := hG0  -- the nonnegativity of `G` is not needed (it follows from the bound)
  intro hbound
  refine no_clean_bound J h hJ G hG 0 fun N x y t _ => ?_
  refine le_trans ?_ (hbound N x y t)
  rw [norm_dressedAnticomm_eq]
  exact Finset.le_sup' (fun q : Cfg N => ‖evo (hq J h q) t x y‖) (Finset.mem_univ 0)

/-- **No operator-norm zero-velocity bound for `[c̃_x(t), n_y]`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_cn (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ),
      ‖heisT J h t (cplx (dressed x)) * cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) -
          cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) * heisT J h t (cplx (dressed x))‖ ≤
        G (Nat.dist x.val y.val) := by
  have _ := hG0  -- the nonnegativity of `G` is not needed (it follows from the bound)
  intro hbound
  refine no_clean_bound J h hJ G hG 0 fun N x y t _ => ?_
  refine le_trans ?_ (hbound N x y t)
  rw [norm_comm_number_eq]
  exact Finset.le_sup' (fun q : Cfg N => ‖evo (hq J h q) t x y‖) (Finset.mem_univ 0)

end DFL
