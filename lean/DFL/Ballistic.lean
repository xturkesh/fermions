import DFL.InfVol

/-!
# The clean sector is ballistic on `ℤ` (paper: Remark `rem:gaussproj`(iii))

For the clean sector on `ℤ`, `h_ℤ = -JΔ + 2h` (`UZ J h (fun _ => 0)`, the exponential series of the
tridiagonal `hZ` of `InfVol.lean`), the mean-square displacement of a particle started at `x` is exactly
`∑_y (y - x)² |⟨δ_y, e^{-ith_ℤ} δ_x⟩|² = 2J²t²` for all `t` (`ballistic_clean`): on finitely supported
vectors `e^{ith_ℤ} X e^{-ith_ℤ} = X + tV`, with `V = i[h_ℤ, X]`, because `V` commutes with `h_ℤ`, and
`‖V δ_x‖² = 2J²`. Equivalently, in the power series of `t ↦ ⟨δ_x, e^{ith} (X - x)² e^{-ith} δ_x⟩` the
coefficients are `⟨δ_x, ad_h^n((X - x)²) δ_x⟩ (i^n / n!)`, and `ad_h³((X - x)²) = 0`.

Consequently there is no bound `sup_t |⟨δ_x, e^{-ith_ℤ} δ_y⟩| ≤ F(|x - y|)` with `∑_d d² F(d)² < ∞`
(`no_uniform_decay_clean`).

## Proof

* `BallisticAux.hZpow_comm`: the commutator identity `[X, h^{k+1}] = (k+1) [X, h] h^k` in matrix
  elements, `(y - z) (h^{k+1})_{yz} = -J (k+1) ((h^k)_{y,z+1} - (h^k)_{y,z-1})` (by induction, using
  the translation invariance of the clean `h_ℤ`). Summing the exponential series,
  `(y - x) U_{yx}(t) = iJt (U_{y,x+1}(t) - U_{y,x-1}(t))` (`BallisticAux.UZ_comm`).
* `BallisticAux.hasSum_conj_UZ_mul`: the columns of `U^ℤ_q(t)` are orthonormal (limit of the
  finite-volume unitarity, by dominated convergence with the uniform bound
  `|U^{(M)}_{ya}(t)| ≤ e^{2r} 2^{-|y-a|}` coming from finite propagation).
* Hence `∑_y (y - x)² |U_{yx}|² = J²t² ∑_y |U_{y,x+1} - U_{y,x-1}|² = 2J²t²`.
-/

open Filter Topology

noncomputable section

namespace DFL

namespace BallisticAux

/-! ### Finite volume: finite propagation and unitarity -/

/-- Finite propagation in finite volume: `(h_q^k)_{ij} = 0` if `|i - j| > k`. -/
lemma pow_hq_eq_zero {N : ℕ} (J h : ℝ) (q : Cfg N) (k : ℕ) :
    ∀ i j : Fin N, k < Nat.dist i.val j.val → ((hq J h q) ^ k) i j = 0 := by
  sorry

/-- Entry bound for the finite-volume propagator, with exponential decay in the distance. -/
lemma norm_evo_hq_le {N : ℕ} (J h : ℝ) (q : Cfg N) (t : ℝ) (i j : Fin N) :
    ‖evo (hq J h q) t i j‖ ≤
      (∑' k : ℕ, (2 * (|t| * (2 * |J| + 2 * |h|))) ^ k / k.factorial) *
        (1 / 2 : ℝ) ^ Nat.dist i.val j.val := by
  sorry

/-- Entry bound for the finite-volume propagators on `ℤ`, uniform in the volume. -/
lemma norm_evoZ_le (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (y a : ℤ) :
    ‖evoZ J h q t M y a‖ ≤
      (∑' k : ℕ, (2 * (|t| * (2 * |J| + 2 * |h|))) ^ k / k.factorial) *
        (1 / 2 : ℝ) ^ (y - a).natAbs := by
  sorry

/-- `y ↦ 2^{-|y - a|}` is summable over `ℤ`. -/
lemma summable_half_pow (a : ℤ) : Summable (fun y : ℤ => (1 / 2 : ℝ) ^ (y - a).natAbs) := by
  sorry

/-- Unitarity of the finite-volume propagator, written as a sum over `ℤ`. -/
lemma tsum_star_evoZ_mul (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ) (ha : |a| ≤ M)
    (hb : |b| ≤ M) :
    ∑' y : ℤ, star (evoZ J h q t M y a) * evoZ J h q t M y b = if a = b then 1 else 0 := by
  sorry

/-- **Unitarity in infinite volume**: the columns of `U^ℤ_q(t)` are orthonormal. -/
lemma hasSum_star_UZ_mul (J h : ℝ) (q : CfgZ) (t : ℝ) (a b : ℤ) :
    HasSum (fun y : ℤ => star (UZ J h q t y a) * UZ J h q t y b) (if a = b then 1 else 0) := by
  sorry

/-- Two distinct columns of `U^ℤ_q(t)`: `∑_y |U_{ya} - U_{yb}|² = 2`. -/
lemma hasSum_norm_sub_sq (J h : ℝ) (q : CfgZ) (t : ℝ) (a b : ℤ) (hab : a ≠ b) :
    HasSum (fun y : ℤ => ‖UZ J h q t y a - UZ J h q t y b‖ ^ 2) 2 := by
  sorry

/-! ### The commutator identity in the clean sector -/

/-- The clean charge configuration `q ≡ 0` (all signs `+1`). -/
abbrev q0 : CfgZ := fun _ => 0

/-- Three-term recursion for the powers of `h^ℤ` in the clean sector. -/
lemma hZpow_succ (J h : ℝ) (k : ℕ) (y z : ℤ) :
    hZpow J h q0 (k + 1) y z =
      -J * (hZpow J h q0 k (y - 1) z + hZpow J h q0 k (y + 1) z) +
        2 * h * hZpow J h q0 k y z := by
  have hI : Finset.Icc (y - 1) (y + 1) = {y - 1, y, y + 1} := by
    ext w; simp only [Finset.mem_Icc, Finset.mem_insert, Finset.mem_singleton]; omega
  have e1 : hZ J h q0 y (y - 1) = -J := by
    unfold hZ
    rw [ite_eq_left (show y = y - 1 + 1 ∨ y - 1 = y + 1 by omega),
      ite_eq_right (show ¬ y = y - 1 by omega)]
    ring
  have e2 : hZ J h q0 y y = 2 * h := by
    unfold hZ sgn
    rw [ite_eq_right (show ¬ (y = y + 1 ∨ y = y + 1) by omega), ite_eq_left rfl,
      ite_eq_left rfl]
    ring
  have e3 : hZ J h q0 y (y + 1) = -J := by
    unfold hZ
    rw [ite_eq_left (show y = y + 1 + 1 ∨ y + 1 = y + 1 by omega),
      ite_eq_right (show ¬ y = y + 1 by omega)]
    ring
  simp only [hZpow]
  rw [hI, Finset.sum_insert (by simp only [Finset.mem_insert, Finset.mem_singleton]; omega),
    Finset.sum_pair (by omega), e1, e2, e3]
  ring

/-- Translation invariance of the powers of the clean `h^ℤ`. -/
lemma hZpow_shift (J h : ℝ) (k : ℕ) :
    ∀ y z a : ℤ, hZpow J h q0 k (y + a) (z + a) = hZpow J h q0 k y z := by
  induction k with
  | zero => intro y z a; simp [hZpow]
  | succ k ih =>
    intro y z a
    rw [hZpow_succ, hZpow_succ, show y + a - 1 = (y - 1) + a by ring,
      show y + a + 1 = (y + 1) + a by ring, ih, ih, ih]

/-- The commutator identity `[X, h^{k+1}] = (k+1) [X, h] h^k` in matrix elements. -/
lemma hZpow_comm (J h : ℝ) (k : ℕ) :
    ∀ y z : ℤ, ((y - z : ℤ) : ℝ) * hZpow J h q0 (k + 1) y z =
      -J * ((k : ℝ) + 1) * (hZpow J h q0 k y (z + 1) - hZpow J h q0 k y (z - 1)) := by
  induction k with
  | zero =>
    intro y z
    rw [hZpow_succ]
    simp only [hZpow]
    rcases (show y = z + 1 ∨ y = z - 1 ∨ (y ≠ z + 1 ∧ y ≠ z - 1) by omega) with
      rfl | rfl | ⟨h1, h2⟩
    all_goals split_ifs <;> first | omega | (subst_vars; push_cast; ring)
  | succ k ih =>
    intro y z
    have T1 : hZpow J h q0 (k + 1) (y - 1) z =
        -J * (hZpow J h q0 k (y - 1) (z + 1) + hZpow J h q0 k (y + 1) (z + 1)) +
          2 * h * hZpow J h q0 k y (z + 1) := by
      have := hZpow_shift J h (k + 1) (y - 1) z 1
      rw [show y - 1 + 1 = y by ring] at this
      rw [← this, hZpow_succ]
    have T2 : hZpow J h q0 (k + 1) (y + 1) z =
        -J * (hZpow J h q0 k (y - 1) (z - 1) + hZpow J h q0 k (y + 1) (z - 1)) +
          2 * h * hZpow J h q0 k y (z - 1) := by
      have := hZpow_shift J h (k + 1) (y + 1) z (-1)
      rw [show y + 1 + -1 = y by ring, show z + -1 = z - 1 by ring] at this
      rw [← this, hZpow_succ]
    have I1 := ih (y - 1) z
    have I2 := ih (y + 1) z
    have I3 := ih y z
    rw [hZpow_succ J h (k + 1) y z, hZpow_succ J h k y (z + 1), hZpow_succ J h k y (z - 1)]
    push_cast at I1 I2 I3 ⊢
    linear_combination (-J) * I1 + (-J) * I2 + 2 * h * I3 + (-J) * T1 + J * T2

/-- The commutator identity for the propagator: `(y - x) U_{yx} = iJt (U_{y,x+1} - U_{y,x-1})`. -/
lemma UZ_comm (J h t : ℝ) (x y : ℤ) :
    ((y - x : ℤ) : ℂ) * UZ J h q0 t y x =
      Complex.I * J * t * (UZ J h q0 t y (x + 1) - UZ J h q0 t y (x - 1)) := by
  sorry

/-- The mean-square displacement series of the clean sector has sum `2J²t²`. -/
lemma hasSum_ballistic (J h t : ℝ) (x : ℤ) :
    HasSum (fun y : ℤ => ((y - x : ℤ) : ℝ) ^ 2 * ‖UZ J h q0 t y x‖ ^ 2) (2 * J ^ 2 * t ^ 2) := by
  sorry

end BallisticAux

/-- **The clean sector is ballistic** (paper: `rem:gaussproj`(iii)):
`∑_y (y - x)² |⟨δ_y, e^{-ith_ℤ} δ_x⟩|² = 2J²t²`. -/
theorem ballistic_clean (J h t : ℝ) (x : ℤ) :
    ∑' y : ℤ, ((y - x : ℤ) : ℝ) ^ 2 * ‖UZ J h (fun _ => 0) t y x‖ ^ 2 = 2 * J ^ 2 * t ^ 2 := by
  exact (BallisticAux.hasSum_ballistic J h t x).tsum_eq

/-- The series in `ballistic_clean` converges. -/
theorem summable_ballistic_clean (J h t : ℝ) (x : ℤ) :
    Summable (fun y : ℤ => ((y - x : ℤ) : ℝ) ^ 2 * ‖UZ J h (fun _ => 0) t y x‖ ^ 2) := by
  exact (BallisticAux.hasSum_ballistic J h t x).summable

/-- **No uniform decay in the clean sector** (paper: `rem:gaussproj`(iii)): for `J ≠ 0` there is no `F` with
`∑_d d² F(d)² < ∞` and `|⟨δ_y, e^{-ith_ℤ} δ_x⟩| ≤ F(|x - y|)` for all `t`, `x`, `y`. -/
theorem no_uniform_decay_clean (J h : ℝ) (hJ : J ≠ 0) (F : ℕ → ℝ)
    (hF : Summable (fun d : ℕ => (d : ℝ) ^ 2 * F d ^ 2)) :
    ¬ ∀ (t : ℝ) (x y : ℤ), ‖UZ J h (fun _ => 0) t y x‖ ≤ F (y - x).natAbs := by
  sorry

end DFL
