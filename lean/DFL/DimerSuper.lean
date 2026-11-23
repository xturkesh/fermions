import DFL.DimerSecond
import DFL.DimerState

/-!
# Superdiffusion of the random dimer model (paper: Proposition `prop:dimersuper`, Theorem `thm:dimer`(c),(d))

For the random dimer potential `v = v^{(r)}` on `ℤ` (values `±λ`, `λ = 2|h|/J ∈ (0, 1)`):

* a deterministic bound (paper: eq. `eq:superdet`): with `z = E + i/(2T)`, `L = ⌊T⌋` and
  `X_j(E) = ‖T^z(j+1, 1)‖₂² + ‖T^z(0, -j)‖₂²`,
  `M^ℤ_{v,p}(T; 0) ≥ c₀²/(4πT) ∫_{λ-δ}^{λ+δ} ∑_{j ≤ L} j^p / X_j(E) dE` (Plancherel `dimerMom_plancherel` and the
  Green-function propagation of `DimerLower.lean` on `[-L', L']`, with `‖A⁻¹‖ ≤ ‖A‖` in `SL(2, ℂ)`, and the limit
  `L' → ∞`, `tendsto_dimerMom_Z`);
* `dimer_super_annealed` (paper: `prop:dimersuper`(a)): `E M^ℤ_{v,p}(T; 0) ≥ c T^{p - 1/2}` for `T ≥ T₀`, by
  `δ = T^{-1/2}`, `lem:dimersecond` (`integral_hs2_tmProd_le`) and Jensen's inequality `E[1/X] ≥ 1/E[X]`;
* `dimer_super_ae` (paper: `prop:dimersuper`(b)): almost surely `M^ℤ_{v,p}(T; 0) ≥ T^{p - 1/2 - α}` eventually, for
  every `α > 0` (Cauchy–Schwarz, Markov, Borel–Cantelli along `T = 2^k`);
* `dimerMInf_super`, `dimerMInf_exponent` (paper: `thm:dimer`(c)): `M^{∞,r}_p(T) ≥ c (JT)^{p - 1/2}`;
* `dimerD_unbounded_half`, `dimer_LGT_unbounded_half`, `not_EC_dimer_half` (paper: `thm:dimer`(d) for every
  `p > 1/2`).
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

namespace DimerSuperAux

open DimerLowerAux

/-! ### Hilbert–Schmidt norms of `2 × 2` matrices -/

/-- `‖A‖₂² ≥ 0`. -/
lemma hs2_nonneg (A : Matrix (Fin 2) (Fin 2) ℂ) : 0 ≤ hs2 A := by
  unfold hs2
  positivity

/-- Cauchy–Schwarz in `ℂ²`: `|a x + b y|² ≤ (|a|² + |b|²)(|x|² + |y|²)`. -/
lemma norm_sq_add_mul_le (a b x y : ℂ) :
    ‖a * x + b * y‖ ^ 2 ≤ (‖a‖ ^ 2 + ‖b‖ ^ 2) * (‖x‖ ^ 2 + ‖y‖ ^ 2) := by
  have h1 : ‖a * x + b * y‖ ≤ ‖a‖ * ‖x‖ + ‖b‖ * ‖y‖ := by
    calc ‖a * x + b * y‖ ≤ ‖a * x‖ + ‖b * y‖ := norm_add_le _ _
      _ = ‖a‖ * ‖x‖ + ‖b‖ * ‖y‖ := by rw [norm_mul, norm_mul]
  calc ‖a * x + b * y‖ ^ 2 ≤ (‖a‖ * ‖x‖ + ‖b‖ * ‖y‖) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ (‖a‖ ^ 2 + ‖b‖ ^ 2) * (‖x‖ ^ 2 + ‖y‖ ^ 2) := by
        nlinarith [sq_nonneg (‖a‖ * ‖y‖ - ‖b‖ * ‖x‖)]

/-- `‖A‖₂²` in coordinates. -/
lemma hs2_eq (A : Matrix (Fin 2) (Fin 2) ℂ) :
    hs2 A = ‖A 0 0‖ ^ 2 + ‖A 0 1‖ ^ 2 + (‖A 1 0‖ ^ 2 + ‖A 1 1‖ ^ 2) := by
  simp only [hs2, Fin.sum_univ_two]

open scoped Matrix.Norms.L2Operator in
/-- The `ℓ²` operator norm is bounded by the Hilbert–Schmidt norm: `‖A‖² ≤ ‖A‖₂²`. -/
lemma norm_sq_le_hs2 (A : Matrix (Fin 2) (Fin 2) ℂ) : ‖A‖ ^ 2 ≤ hs2 A := by
  refine DimerTransferAux.norm_sq_le_of_forall (hs2_nonneg A) fun u => ?_
  have e0 : (A *ᵥ u) 0 = A 0 0 * u 0 + A 0 1 * u 1 := by
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  have e1 : (A *ᵥ u) 1 = A 1 0 * u 0 + A 1 1 * u 1 := by
    simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  rw [e0, e1, hs2_eq]
  have h0 := norm_sq_add_mul_le (A 0 0) (A 0 1) (u 0) (u 1)
  have h1 := norm_sq_add_mul_le (A 1 0) (A 1 1) (u 0) (u 1)
  nlinarith [h0, h1]

/-- For `det A = 1`: `‖A‖₂² ≥ 2`. -/
lemma two_le_hs2 (A : Matrix (Fin 2) (Fin 2) ℂ) (hA : A.det = 1) : 2 ≤ hs2 A := by
  rw [Matrix.det_fin_two] at hA
  have h1 : (1 : ℝ) ≤ ‖A 0 0‖ * ‖A 1 1‖ + ‖A 0 1‖ * ‖A 1 0‖ := by
    calc (1 : ℝ) = ‖A 0 0 * A 1 1 - A 0 1 * A 1 0‖ := by rw [hA, norm_one]
      _ ≤ ‖A 0 0 * A 1 1‖ + ‖A 0 1 * A 1 0‖ := norm_sub_le _ _
      _ = ‖A 0 0‖ * ‖A 1 1‖ + ‖A 0 1‖ * ‖A 1 0‖ := by rw [norm_mul, norm_mul]
  rw [hs2_eq]
  nlinarith [sq_nonneg (‖A 0 0‖ - ‖A 1 1‖), sq_nonneg (‖A 0 1‖ - ‖A 1 0‖)]

/-- `X_j(z) = ‖T^z(j+1, 1)‖₂² + ‖T^z(0, -j)‖₂²` (paper: proof of `prop:dimersuper`). -/
def Xj (v : ℤ → ℝ) (z : ℂ) (j : ℕ) : ℝ := hs2 (tmProd v z 1 j) + hs2 (tmProd v z (-(j : ℤ)) j)

/-- `X_j ≥ 4`, since the transfer matrices have determinant one. -/
lemma four_le_Xj (v : ℤ → ℝ) (z : ℂ) (j : ℕ) : 4 ≤ Xj v z j := by
  have h1 := two_le_hs2 _ (det_tmProd v z 1 j)
  have h2 := two_le_hs2 _ (det_tmProd v z (-(j : ℤ)) j)
  unfold Xj
  linarith

/-- The tangent line of `1/X` at `a`: `2/a - X/a² ≤ 1/X`. -/
lemma two_div_sub_le_inv (X a : ℝ) (hX : 0 < X) (ha : 0 < a) : 2 / a - X / a ^ 2 ≤ 1 / X := by
  have e : 1 / X - (2 / a - X / a ^ 2) = (a - X) ^ 2 / (a ^ 2 * X) := by
    field_simp
    ring
  have : 0 ≤ (a - X) ^ 2 / (a ^ 2 * X) := by positivity
  linarith

/-! ### Steps 2–3 of `prop:dimerlower` with `X_j` (paper: proof of `prop:dimersuper`) -/

section Green

open scoped Matrix.Norms.L2Operator

/-- **Steps 2–3 with `X_j` in place of the uniform transfer-matrix bound**: at the origin of
`[a, b] ⊇ [-L, L]`, `c₀² ∑_{j ≤ L} j^p / X_j(z) ≤ 2 ∑_y |y|^p |G_z(y)|²`, `c₀ = 1/(3 + 2λ)`. -/
lemma green_moment_ge_X (lam p : ℝ) (hl : 0 ≤ lam) (hp : 0 < p) (a b : ℤ) (ha : a ≤ 0)
    (hb : 0 ≤ b) (v : ℤ → ℝ) (L : ℕ) (haL : a ≤ -(L : ℤ)) (hbL : (L : ℤ) ≤ b) (z : ℂ)
    (hR : (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) *
      (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ))⁻¹ = 1)
    (hvz : ‖(v 0 : ℂ) - z‖ ≤ 2 * lam + 1) :
    (1 / (3 + 2 * lam)) ^ 2 * ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p / Xj v z (i + 1) ≤
      2 * ∑ y : Set.Icc a b, |((y.1 - (origin a b ha hb).1 : ℤ) : ℝ)| ^ p *
        ‖dimerGreen a b v z (origin a b ha hb) y‖ ^ 2 := by
  set x := origin a b ha hb with hxdef
  have hx0 : x.1 = 0 := rfl
  set G := DimerLowerAux.ext a b (dimerGreen a b v z x) with hG
  set S := ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖dimerGreen a b v z x y‖ ^ 2
    with hS
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  have hc₀0 : 0 < c₀ := by positivity
  have hc₀eq : (3 + 2 * lam) * c₀ = 1 := by rw [hc₀]; field_simp
  -- the equation `(H - z) G = δ_0`
  have hrec : ∀ k, a ≤ k → k ≤ b → k ≠ 0 →
      G (k + 1) = ((v k : ℂ) - z) * G k - G (k - 1) := by
    intro k hka hkb hk0
    have h := green_eq a b v z hR x ⟨k, hka, hkb⟩
    rw [ite_eq_right (fun e => hk0 ((congrArg Subtype.val e).trans hx0))] at h
    linear_combination -h
  have hX : ((v 0 : ℂ) - z) * G 0 - G 1 - G (-1) = 1 := by
    have h := green_eq a b v z hR x x
    rw [ite_eq_left rfl, hx0, zero_add, zero_sub] at h
    exact h
  have hthree : 1 ≤ (2 * lam + 1) * ‖G 0‖ + ‖G 1‖ + ‖G (-1)‖ := by
    have h := congrArg norm hX
    rw [norm_one] at h
    calc (1 : ℝ) = ‖((v 0 : ℂ) - z) * G 0 - G 1 - G (-1)‖ := h.symm
      _ ≤ ‖((v 0 : ℂ) - z) * G 0 - G 1‖ + ‖G (-1)‖ := norm_sub_le _ _
      _ ≤ ‖((v 0 : ℂ) - z) * G 0‖ + ‖G 1‖ + ‖G (-1)‖ := by
          gcongr
          exact norm_sub_le _ _
      _ = ‖(v 0 : ℂ) - z‖ * ‖G 0‖ + ‖G 1‖ + ‖G (-1)‖ := by rw [norm_mul]
      _ ≤ (2 * lam + 1) * ‖G 0‖ + ‖G 1‖ + ‖G (-1)‖ := by gcongr
  have hXpos : ∀ j, 0 < Xj v z j := fun j => by linarith [four_le_Xj v z j]
  -- `c₀²/X ≤ ‖Φ'‖²` from `c₀² ≤ ‖Φ‖² ≤ X ‖Φ'‖²`
  have hdiv : ∀ (X n m : ℝ), 0 < X → c₀ ^ 2 ≤ n → n ≤ X * m → c₀ ^ 2 / X ≤ m := by
    intro X n m hX h1 h2
    rw [div_le_iff₀ hX]
    calc c₀ ^ 2 ≤ n := h1
      _ ≤ X * m := h2
      _ = m * X := mul_comm _ _
  have hmain : ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1)) ≤ 2 * S := by
    by_cases hcase : c₀ ≤ ‖G 1‖ ∨ c₀ ≤ ‖G 0‖
    · -- `‖Φ(1)‖ ≥ c₀`: propagate to the right
      have hΦ : c₀ ^ 2 ≤ nsq ![G 1, G (1 - 1)] := by
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one, sub_self]
        rcases hcase with h | h
        · nlinarith [sq_nonneg ‖G 0‖]
        · nlinarith [sq_nonneg ‖G 1‖]
      have hj : ∀ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1)) ≤
          |((1 + ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (1 + ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 +
            |((1 + ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (1 + ((i + 1 : ℕ) : ℤ))‖ ^ 2 := by
        intro i hi
        rw [Finset.mem_range] at hi
        have htr := transfer v z G 1 (i + 1) (fun i' hi' =>
          hrec (1 + i') (by omega) (by omega) (by omega))
        have h1' := nsq_le_of_det _ (det_tmProd v z 1 (i + 1)) ![G 1, G (1 - 1)]
        rw [← htr] at h1'
        have hA : ‖tmProd v z 1 (i + 1)‖ ^ 2 ≤ Xj v z (i + 1) := by
          have := norm_sq_le_hs2 (tmProd v z 1 (i + 1))
          have := hs2_nonneg (tmProd v z (-((i + 1 : ℕ) : ℤ)) (i + 1))
          unfold Xj
          linarith
        have h2' : nsq ![G 1, G (1 - 1)] ≤ Xj v z (i + 1) *
            nsq ![G (1 + ((i + 1 : ℕ) : ℤ)), G (1 + ((i + 1 : ℕ) : ℤ) - 1)] :=
          h1'.trans (mul_le_mul_of_nonneg_right hA (nsq_nonneg _))
        have h3' := hdiv _ _ _ (hXpos (i + 1)) hΦ h2'
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one] at h3'
        simp only [hx0, sub_zero]
        rw [abs_cast_rpow _ ((i : ℝ) + 1) (by push_cast; ring) (by positivity),
          abs_cast_rpow _ ((i : ℝ) + 2) (by push_cast; ring) (by positivity)]
        exact weighted_le _ _ _ p i hp (by positivity) (by linarith)
      calc ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1))
          ≤ ∑ i ∈ Finset.range L,
            (|((1 + ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (1 + ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 +
              |((1 + ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (1 + ((i + 1 : ℕ) : ℤ))‖ ^ 2) := Finset.sum_le_sum hj
        _ = ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => 1 + ((i + 1 : ℕ) : ℤ) - 1) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => 1 + ((i + 1 : ℕ) : ℤ) - 1) i)‖ ^ 2 +
            ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => 1 + ((i + 1 : ℕ) : ℤ)) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => 1 + ((i + 1 : ℕ) : ℤ)) i)‖ ^ 2 :=
            Finset.sum_add_distrib
        _ ≤ S + S := add_le_add
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
        _ = 2 * S := by ring
    · -- otherwise `|G(-1)| ≥ c₀`, so `‖Φ(0)‖ ≥ c₀`: propagate to the left
      rw [not_or, not_le, not_le] at hcase
      have hxm : c₀ ≤ ‖G (-1)‖ := by
        have := mul_le_mul_of_nonneg_left hcase.2.le (by linarith : (0 : ℝ) ≤ 2 * lam + 1)
        nlinarith [norm_nonneg (G 0)]
      have hΦ : c₀ ^ 2 ≤ nsq ![G 0, G (0 - 1)] := by
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one, zero_sub]
        nlinarith [sq_nonneg ‖G 0‖]
      have hj : ∀ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1)) ≤
          |((-((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (-((i + 1 : ℕ) : ℤ))‖ ^ 2 +
            |((-((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (-((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 := by
        intro i hi
        rw [Finset.mem_range] at hi
        have htr := transfer v z G (-((i + 1 : ℕ) : ℤ)) (i + 1) (fun i' hi' =>
          hrec (-((i + 1 : ℕ) : ℤ) + i') (by omega) (by omega) (by omega))
        have e1 : -((i + 1 : ℕ) : ℤ) + ((i + 1 : ℕ) : ℤ) = 0 := by ring
        rw [e1] at htr
        have h1' := nsq_mulVec_le (tmProd v z (-((i + 1 : ℕ) : ℤ)) (i + 1))
          ![G (-((i + 1 : ℕ) : ℤ)), G (-((i + 1 : ℕ) : ℤ) - 1)]
        rw [← htr] at h1'
        have hA : ‖tmProd v z (-((i + 1 : ℕ) : ℤ)) (i + 1)‖ ^ 2 ≤ Xj v z (i + 1) := by
          have := norm_sq_le_hs2 (tmProd v z (-((i + 1 : ℕ) : ℤ)) (i + 1))
          have := hs2_nonneg (tmProd v z 1 (i + 1))
          unfold Xj
          linarith
        have h2' : nsq ![G 0, G (0 - 1)] ≤ Xj v z (i + 1) *
            nsq ![G (-((i + 1 : ℕ) : ℤ)), G (-((i + 1 : ℕ) : ℤ) - 1)] :=
          h1'.trans (mul_le_mul_of_nonneg_right hA (nsq_nonneg _))
        have h3' := hdiv _ _ _ (hXpos (i + 1)) hΦ h2'
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one] at h3'
        simp only [hx0, sub_zero]
        rw [abs_cast_rpow' _ ((i : ℝ) + 1) (by push_cast; ring) (by positivity),
          abs_cast_rpow' _ ((i : ℝ) + 2) (by push_cast; ring) (by positivity)]
        exact weighted_le _ _ _ p i hp (by positivity) h3'
      calc ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1))
          ≤ ∑ i ∈ Finset.range L,
            (|((-((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (-((i + 1 : ℕ) : ℤ))‖ ^ 2 +
              |((-((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (-((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2) := Finset.sum_le_sum hj
        _ = ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => -((i + 1 : ℕ) : ℤ)) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => -((i + 1 : ℕ) : ℤ)) i)‖ ^ 2 +
            ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => -((i + 1 : ℕ) : ℤ) - 1) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => -((i + 1 : ℕ) : ℤ) - 1) i)‖ ^ 2 :=
            Finset.sum_add_distrib
        _ ≤ S + S := add_le_add
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
        _ = 2 * S := by ring
  calc c₀ ^ 2 * ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p / Xj v z (i + 1)
      = ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * (c₀ ^ 2 / Xj v z (i + 1)) := by
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl fun i _ => ?_
        ring
    _ ≤ 2 * S := hmain

end Green

/-! ### The deterministic bound (paper: eq. `eq:superdet`) -/

/-- `a ↦ T_a` is continuous. -/
lemma continuous_tmC : Continuous (fun a : ℂ => tmC a) := by
  refine continuous_matrix fun i j => ?_
  fin_cases i <;> fin_cases j <;> simp [tmC] <;> fun_prop

/-- `z ↦ T^z(m + L, m)` is continuous. -/
lemma continuous_tmProd (v : ℤ → ℝ) (m : ℤ) (L : ℕ) :
    Continuous (fun z : ℂ => tmProd v z m L) := by
  induction L with
  | zero => exact continuous_const
  | succ L ih =>
    exact (continuous_tmC.comp (continuous_const.sub continuous_id)).matrix_mul ih

/-- `‖·‖₂²` is continuous. -/
lemma continuous_hs2 : Continuous hs2 := by
  unfold hs2
  exact continuous_finsetSum _ fun i _ => continuous_finsetSum _ fun j _ =>
    ((continuous_id.matrix_elem i j).norm).pow 2

/-- The energy integrand `Y(E) = ∑_{i<L} (i+1)^p X_{i+1}(E + i/(2T))`. -/
def Yfun (v : ℤ → ℝ) (p T : ℝ) (L : ℕ) (E : ℝ) : ℝ :=
  ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1)

/-- `Y` is continuous in the energy. -/
lemma continuous_Yfun (v : ℤ → ℝ) (p T : ℝ) (L : ℕ) : Continuous (Yfun v p T L) := by
  have hz : Continuous (fun E : ℝ => (E : ℂ) + Complex.I / (2 * T)) :=
    Complex.continuous_ofReal.add continuous_const
  unfold Yfun Xj
  exact continuous_finsetSum _ fun i _ => continuous_const.mul
    ((continuous_hs2.comp ((continuous_tmProd v 1 (i + 1)).comp hz)).add
      (continuous_hs2.comp ((continuous_tmProd v _ (i + 1)).comp hz)))

/-- `B = ∫_{λ-δ}^{λ+δ} Y(E) dE`. -/
def Bfun (v : ℤ → ℝ) (p T δ lam : ℝ) (L : ℕ) : ℝ :=
  ∫ E in Set.Icc (lam - δ) (lam + δ), Yfun v p T L E

/-- `|w - z| ≤ 2λ + 1` for `w = ±λ` and `z = E + i/(2T)` with `|E - λ| ≤ 1/2`, `T ≥ 1`. -/
lemma norm_v_sub_z_le (lam : ℝ) (hl : 0 < lam) (w : ℝ) (hw : w = lam ∨ w = -lam) (E T : ℝ)
    (hT : 1 ≤ T) (hE : |E - lam| ≤ 1 / 2) :
    ‖(w : ℂ) - ((E : ℂ) + Complex.I / (2 * T))‖ ≤ 2 * lam + 1 := by
  have hT0 : 0 < T := by linarith
  have e : (w : ℂ) - ((E : ℂ) + Complex.I / (2 * T)) =
      ((w - E : ℝ) : ℂ) + ((-(1 / (2 * T)) : ℝ) : ℂ) * Complex.I := by
    push_cast
    ring
  have h1 : |w - E| ≤ |w - lam| + |E - lam| := by
    rw [abs_sub_comm E lam]
    exact abs_sub_le w lam E
  have h2 : |w - lam| ≤ 2 * lam := by
    rcases hw with hw | hw <;> subst hw <;> rw [abs_le] <;> constructor <;> linarith
  have h3 : 1 / (2 * T) ≤ 1 / 2 := one_div_le_one_div_of_le (by norm_num) (by linarith)
  rw [e]
  calc ‖((w - E : ℝ) : ℂ) + ((-(1 / (2 * T)) : ℝ) : ℂ) * Complex.I‖
      ≤ ‖((w - E : ℝ) : ℂ)‖ + ‖((-(1 / (2 * T)) : ℝ) : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |w - E| + 1 / (2 * T) := by
        rw [Complex.norm_real, norm_mul, Complex.norm_real, Complex.norm_I, mul_one,
          Real.norm_eq_abs, Real.norm_eq_abs, abs_neg,
          abs_of_pos (by positivity : (0 : ℝ) < 1 / (2 * T))]
    _ ≤ 2 * lam + 1 := by linarith

/-- **The deterministic bound on `[-n, n]`** (Plancherel and Steps 2–3, with the tangent-line
bound `1/X ≥ 2/a - X/a²`). -/
lemma det_bound_fin (lam : ℝ) (hl : 0 < lam) (v : ℤ → ℝ) (hv0 : v 0 = lam ∨ v 0 = -lam)
    (p T δ : ℝ) (hp : 0 < p) (hT : 1 ≤ T) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1 / 2) (L n : ℕ)
    (hLn : L ≤ n) (a : ℝ) (ha : 0 < a) :
    1 / (2 * Real.pi * T) * ((1 / (3 + 2 * lam)) ^ 2 / 2 *
      (2 * δ * (2 / a) * ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p - Bfun v p T δ lam L / a ^ 2)) ≤
    dimerMom (-(n : ℤ)) n v p T (origin (-(n : ℤ)) n (by omega) (by omega)) := by
  have hT0 : 0 < T := by linarith
  have ha0 : -(n : ℤ) ≤ 0 := by omega
  have hb0 : (0 : ℤ) ≤ n := by omega
  rw [dimerMom_plancherel _ _ v p T hp hT0]
  set o := origin (-(n : ℤ)) n ha0 hb0 with ho
  set F : ℝ → ℝ := fun E => ∑ y : Set.Icc (-(n : ℤ)) n, |((y.1 - o.1 : ℤ) : ℝ)| ^ p *
    ‖dimerGreen (-(n : ℤ)) n v (E + Complex.I / (2 * T)) o y‖ ^ 2 with hF
  have hF0 : ∀ E, 0 ≤ F E := fun E => Finset.sum_nonneg fun y _ => by positivity
  have hFint : Integrable F :=
    DimerLowerAux.integrable_green_sum (schrH _ _ v) (DimerTLAux.schrH_isSymm _ _ v)
      (fun y => |((y.1 - o.1 : ℤ) : ℝ)| ^ p) T hT0 o
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  set SL := ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p with hSL
  set g : ℝ → ℝ := fun E => c₀ ^ 2 / 2 * (2 / a * SL - Yfun v p T L E / a ^ 2) with hg
  have hgF : ∀ E ∈ Set.Icc (lam - δ) (lam + δ), g E ≤ F E := by
    intro E hE
    have hzim : ((E : ℂ) + Complex.I / (2 * T)).im ≠ 0 := by
      rw [DimerLowerAux.im_add_I_div]
      positivity
    have hR := DimerLowerAux.resolvent_mul_inv (schrH (-(n : ℤ)) n v)
      (DimerTLAux.schrH_isSymm _ _ v) _ hzim
    have hE' : |E - lam| ≤ 1 / 2 := by
      rw [abs_le]
      constructor <;> linarith [hE.1, hE.2]
    have hvz := norm_v_sub_z_le lam hl (v 0) hv0 E T hT hE'
    have h : c₀ ^ 2 * ∑ i ∈ Finset.range L,
        ((i : ℝ) + 1) ^ p / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) ≤ 2 * F E :=
      green_moment_ge_X lam p hl.le hp _ _ ha0 hb0 v L (by omega) (by omega) _ hR hvz
    have h2 : 2 / a * SL - Yfun v p T L E / a ^ 2 ≤ ∑ i ∈ Finset.range L,
        ((i : ℝ) + 1) ^ p / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) := by
      have e : 2 / a * SL - Yfun v p T L E / a ^ 2 = ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p *
          (2 / a - Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) / a ^ 2) := by
        simp only [hSL, Yfun, Finset.mul_sum, Finset.sum_div, ← Finset.sum_sub_distrib]
        refine Finset.sum_congr rfl fun i _ => ?_
        ring
      rw [e]
      refine Finset.sum_le_sum fun i _ => ?_
      have hX := four_le_Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1)
      calc ((i : ℝ) + 1) ^ p * (2 / a - Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) / a ^ 2)
          ≤ ((i : ℝ) + 1) ^ p * (1 / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1)) :=
            mul_le_mul_of_nonneg_left (two_div_sub_le_inv _ _ (by linarith) ha) (by positivity)
        _ = ((i : ℝ) + 1) ^ p / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) := mul_one_div _ _
    calc g E = c₀ ^ 2 / 2 * (2 / a * SL - Yfun v p T L E / a ^ 2) := rfl
      _ ≤ c₀ ^ 2 / 2 * ∑ i ∈ Finset.range L,
          ((i : ℝ) + 1) ^ p / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1) := by gcongr
      _ = (c₀ ^ 2 * ∑ i ∈ Finset.range L,
          ((i : ℝ) + 1) ^ p / Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1)) / 2 := by ring
      _ ≤ (2 * F E) / 2 := by gcongr
      _ = F E := by ring
  have hYint : IntegrableOn (Yfun v p T L) (Set.Icc (lam - δ) (lam + δ)) :=
    (continuous_Yfun v p T L).integrableOn_Icc
  have hgint : IntegrableOn g (Set.Icc (lam - δ) (lam + δ)) :=
    (continuous_const.mul (continuous_const.sub
      ((continuous_Yfun v p T L).div_const _))).integrableOn_Icc
  have h1 : ∫ E in Set.Icc (lam - δ) (lam + δ), g E ≤ ∫ E in Set.Icc (lam - δ) (lam + δ), F E :=
    setIntegral_mono_on hgint hFint.integrableOn measurableSet_Icc hgF
  have h2 : ∫ E in Set.Icc (lam - δ) (lam + δ), F E ≤ ∫ E, F E :=
    setIntegral_le_integral hFint (ae_of_all _ hF0)
  have hvol : volume.real (Set.Icc (lam - δ) (lam + δ)) = 2 * δ := by
    rw [Measure.real, Real.volume_Icc, ENNReal.toReal_ofReal (by linarith)]
    ring
  have h3 : ∫ E in Set.Icc (lam - δ) (lam + δ), g E =
      c₀ ^ 2 / 2 * (2 * δ * (2 / a) * SL - Bfun v p T δ lam L / a ^ 2) := by
    simp only [hg]
    rw [integral_const_mul, integral_sub (continuous_const.integrableOn_Icc)
      (hYint.div_const _), setIntegral_const, integral_div, hvol]
    simp only [Bfun, smul_eq_mul]
    ring
  calc 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / a) * SL - Bfun v p T δ lam L / a ^ 2))
      = 1 / (2 * Real.pi * T) * ∫ E in Set.Icc (lam - δ) (lam + δ), g E := by rw [h3]
    _ ≤ 1 / (2 * Real.pi * T) * ∫ E, F E := by
        gcongr
        exact h1.trans h2

/-- **The deterministic bound** (paper: eq. `eq:superdet`, with the tangent-line bound): for a
bounded potential with `v(0) = ±λ`, all `a > 0`, `T ≥ 1`, `0 < δ ≤ 1/2` and `L`,
`M^ℤ_{v,p}(T; 0) ≥ (2πT)⁻¹ (c₀²/2) (2δ (2/a) ∑_{j ≤ L} j^p - B/a²)`. -/
lemma det_bound (lam : ℝ) (hl : 0 < lam) (v : ℤ → ℝ) (W : ℝ) (hW : ∀ y, |v y| ≤ W)
    (hv0 : v 0 = lam ∨ v 0 = -lam) (p T δ : ℝ) (hp : 0 < p) (hT : 1 ≤ T) (hδ0 : 0 < δ)
    (hδ1 : δ ≤ 1 / 2) (L : ℕ) (a : ℝ) (ha : 0 < a) :
    1 / (2 * Real.pi * T) * ((1 / (3 + 2 * lam)) ^ 2 / 2 *
      (2 * δ * (2 / a) * ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p - Bfun v p T δ lam L / a ^ 2)) ≤
    dimerMomZ v p T 0 := by
  have hlim := tendsto_dimerMom_Z v W hW (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ))
    (fun n => by simp) (fun n => by positivity)
    (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop p T hp
    (by linarith)
  refine ge_of_tendsto hlim ?_
  filter_upwards [eventually_ge_atTop L] with n hn
  exact det_bound_fin lam hl v hv0 p T δ hp hT hδ0 hδ1 L n hn a ha

/-! ### Functions of finitely many signs -/

/-- The restriction of a sign configuration to a finite set of sites. -/
def restr (S : Finset ℤ) (η : CfgZ) : S → ZMod 2 := fun k => η k.1

/-- The extension by `0` of a sign configuration on a finite set of sites. -/
def liftC (S : Finset ℤ) (c : S → ZMod 2) : CfgZ := fun k => if hk : k ∈ S then c ⟨k, hk⟩ else 0

/-- The restriction is measurable. -/
lemma measurable_restr (S : Finset ℤ) : Measurable (restr S) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- `liftC ∘ restr` does not change the signs in `S`. -/
lemma liftC_restr (S : Finset ℤ) (η : CfgZ) (k : ℤ) (hk : k ∈ S) :
    liftC S (restr S η) k = η k := by
  simp [liftC, restr, hk]

/-- A function of the signs in `S` factors through the restriction to `S`. -/
lemma dep_factor {β : Type*} (S : Finset ℤ) (F : CfgZ → β)
    (hF : ∀ η η' : CfgZ, (∀ k ∈ S, η k = η' k) → F η = F η') :
    F = (F ∘ liftC S) ∘ restr S :=
  funext fun η => hF _ _ fun k hk => (liftC_restr S η k hk).symm

/-- A function of finitely many signs is measurable. -/
lemma measurable_dep {β : Type*} [MeasurableSpace β] (S : Finset ℤ) (F : CfgZ → β)
    (hF : ∀ η η' : CfgZ, (∀ k ∈ S, η k = η' k) → F η = F η') : Measurable F := by
  rw [dep_factor S F hF]
  exact (measurable_of_finite _).comp (measurable_restr S)

/-- A real function of finitely many signs is integrable. -/
lemma integrable_dep (S : Finset ℤ) (F : CfgZ → ℝ)
    (hF : ∀ η η' : CfgZ, (∀ k ∈ S, η k = η' k) → F η = F η') : Integrable F muZ := by
  rw [dep_factor S F hF]
  exact (integrable_map_measure (measurable_of_finite _).aestronglyMeasurable
    (measurable_restr S).aemeasurable).mp Integrable.of_finite

/-- **Fubini for functions of finitely many signs**: `∫ dμ_ℤ` and `∫_s dE` commute. -/
lemma integral_swap_dep (S : Finset ℤ) (f : CfgZ → ℝ → ℝ)
    (hf : ∀ η η' : CfgZ, (∀ k ∈ S, η k = η' k) → f η = f η') (s : Set ℝ)
    (hint : ∀ η, IntegrableOn (f η) s) :
    ∫ η, (∫ E in s, f η E) ∂muZ = ∫ E in s, (∫ η, f η E ∂muZ) := by
  classical
  set ν := muZ.map (restr S) with hν
  have hm := measurable_restr S
  have e1 : ∀ η, f η = f (liftC S (restr S η)) := fun η =>
    hf _ _ fun k hk => (liftC_restr S η k hk).symm
  have hL : ∫ η, (∫ E in s, f η E) ∂muZ = ∑ c, ν.real {c} • ∫ E in s, f (liftC S c) E := by
    calc ∫ η, (∫ E in s, f η E) ∂muZ
        = ∫ η, (fun c => ∫ E in s, f (liftC S c) E) (restr S η) ∂muZ := by
          congr 1
          funext η
          rw [e1 η]
      _ = ∫ c, (fun c => ∫ E in s, f (liftC S c) E) c ∂ν :=
          (integral_map (f := fun c => ∫ E in s, f (liftC S c) E) hm.aemeasurable
            (measurable_of_finite _).aestronglyMeasurable).symm
      _ = ∑ c, ν.real {c} • ∫ E in s, f (liftC S c) E := integral_fintype Integrable.of_finite
  have hR : ∀ E, ∫ η, f η E ∂muZ = ∑ c, ν.real {c} • f (liftC S c) E := by
    intro E
    calc ∫ η, f η E ∂muZ = ∫ η, (fun c => f (liftC S c) E) (restr S η) ∂muZ := by
          congr 1
          funext η
          rw [e1 η]
      _ = ∫ c, (fun c => f (liftC S c) E) c ∂ν :=
          (integral_map (f := fun c => f (liftC S c) E) hm.aemeasurable
            (measurable_of_finite _).aestronglyMeasurable).symm
      _ = ∑ c, ν.real {c} • f (liftC S c) E := integral_fintype Integrable.of_finite
  rw [hL]
  simp_rw [hR]
  rw [integral_finsetSum]
  · refine Finset.sum_congr rfl fun c _ => ?_
    rw [integral_smul]
  · intro c _
    exact Integrable.smul (ν.real {c}) (hint (liftC S c))

/-- The sites `{⌊(y + r)/2⌋ : |y| ≤ n}` carrying the signs of `v^{(r)}` on `[-n, n]`. -/
def sitesS (r : ℤ) (n : ℕ) : Finset ℤ := (Finset.Icc (-(n : ℤ)) n).image (fun y => (y + r) / 2)

/-- `v^{(r)}` on `[-n, n]` depends only on the signs in `sitesS r n`. -/
lemma dimerPotZ_congr (J h : ℝ) (r : ℤ) (n : ℕ) (η η' : CfgZ)
    (hη : ∀ k ∈ sitesS r n, η k = η' k) (y : ℤ) (hy1 : -(n : ℤ) ≤ y) (hy2 : y ≤ n) :
    dimerPotZ J h r η y = dimerPotZ J h r η' y := by
  have hk : (y + r) / 2 ∈ sitesS r n :=
    Finset.mem_image.mpr ⟨y, Finset.mem_Icc.mpr ⟨hy1, hy2⟩, rfl⟩
  simp only [dimerPotZ, hη _ hk]

/-- `T^z(m + L, m)` depends only on `v(m), …, v(m + L - 1)`. -/
lemma tmProd_congr (v w : ℤ → ℝ) (z : ℂ) (m : ℤ) (L : ℕ)
    (h : ∀ y, m ≤ y → y < m + L → v y = w y) : tmProd v z m L = tmProd w z m L := by
  induction L with
  | zero => rfl
  | succ L ih =>
    simp only [tmProd]
    rw [h (m + L) (by omega) (by push_cast; omega),
      ih fun y h1 h2 => h y h1 (by push_cast; omega)]

/-- `X_j` depends only on `v` on `[-j, j]`. -/
lemma Xj_congr (v w : ℤ → ℝ) (z : ℂ) (j : ℕ) (h : ∀ y, -(j : ℤ) ≤ y → y ≤ j → v y = w y) :
    Xj v z j = Xj w z j := by
  unfold Xj
  rw [tmProd_congr v w z 1 j fun y h1 h2 => h y (by omega) (by omega),
    tmProd_congr v w z (-(j : ℤ)) j fun y h1 h2 => h y (by omega) (by omega)]

/-- `Y` depends only on `v` on `[-L, L]`. -/
lemma Yfun_congr (v w : ℤ → ℝ) (p T : ℝ) (L : ℕ) (h : ∀ y, -(L : ℤ) ≤ y → y ≤ L → v y = w y) :
    Yfun v p T L = Yfun w p T L := by
  funext E
  unfold Yfun
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.mem_range] at hi
  rw [Xj_congr v w _ (i + 1) fun y h1 h2 => h y (by omega) (by omega)]

/-- `Y ≥ 0`. -/
lemma Yfun_nonneg (v : ℤ → ℝ) (p T : ℝ) (L : ℕ) (E : ℝ) : 0 ≤ Yfun v p T L E :=
  Finset.sum_nonneg fun i _ => mul_nonneg (by positivity)
    (by linarith [four_le_Xj v ((E : ℂ) + Complex.I / (2 * T)) (i + 1)])

/-- `H_v` on `[a, b]` depends only on `v` on `[a, b]`. -/
lemma schrH_congr (a b : ℤ) (v w : ℤ → ℝ) (h : ∀ y, a ≤ y → y ≤ b → v y = w y) :
    schrH a b v = schrH a b w := by
  ext x y
  simp only [schrH, h x.1 x.2.1 x.2.2]

/-- `M_p(T; x)` on `[a, b]` depends only on `v` on `[a, b]`. -/
lemma dimerMom_congr (a b : ℤ) (v w : ℤ → ℝ) (h : ∀ y, a ≤ y → y ≤ b → v y = w y) (p T : ℝ)
    (x : Set.Icc a b) : dimerMom a b v p T x = dimerMom a b w p T x := by
  unfold dimerMom
  rw [schrH_congr a b v w h]

/-- The finite-volume approximations of `M^ℤ_{v^{(r)},p}(T; 0)` on `[-n, n]`. -/
lemma tendsto_dimerMom_pot (J h : ℝ) (r : ℤ) (η : CfgZ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    Tendsto (fun n : ℕ => dimerMom (-(n : ℤ)) n (dimerPotZ J h r η) p T
      (origin (-(n : ℤ)) n (by omega) (by omega))) atTop
      (𝓝 (dimerMomZ (dimerPotZ J h r η) p T 0)) :=
  tendsto_dimerMom_Z (dimerPotZ J h r η) |2 * h / J| (DimerInfVolAux.abs_dimerPotZ_le J h r η)
    (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ)) (fun n => by simp) (fun n => by positivity)
    (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop p T hp
    hT

/-- `η ↦ M^ℤ_{v^{(r)},p}(T; 0)` is measurable (a limit of functions of finitely many signs). -/
lemma measurable_M (J h : ℝ) (r : ℤ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    Measurable (fun η => dimerMomZ (dimerPotZ J h r η) p T 0) := by
  refine measurable_of_tendsto_metrizable (f := fun (n : ℕ) η => dimerMom (-(n : ℤ)) n
    (dimerPotZ J h r η) p T (origin (-(n : ℤ)) n (by omega) (by omega))) (fun n => ?_) ?_
  · exact measurable_dep (sitesS r n) _ fun η η' hηη' =>
      dimerMom_congr _ _ _ _ (fun y hy1 hy2 => dimerPotZ_congr J h r n η η' hηη' y hy1 hy2) p T _
  · rw [tendsto_pi_nhds]
    intro η
    exact tendsto_dimerMom_pot J h r η p T hp hT

/-- `|M^ℤ_{v,p}(T; 0)| ≤ ∫ momBound`, uniformly in the bounded potential. -/
lemma norm_dimerMomZ_le (v : ℤ → ℝ) (W : ℝ) (hW : ∀ y, |v y| ≤ W) (p T : ℝ) (hp : 0 < p)
    (hT : 0 < T) :
    ‖dimerMomZ v p T 0‖ ≤ ∫ t in Set.Ioi (0 : ℝ), DimerInfVolAux.momBound p T t := by
  have hlim := tendsto_dimerMom_Z v W hW (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ))
    (fun n => by simp) (fun n => by positivity)
    (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop p T hp
    hT
  exact le_of_tendsto hlim.norm (Eventually.of_forall fun n =>
    DimerInfVolAux.norm_dimerMom_le _ _ _ _ v p T hp hT)

/-- `η ↦ M^ℤ_{v^{(r)},p}(T; 0)` is integrable. -/
lemma integrable_M (J h : ℝ) (r : ℤ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    Integrable (fun η => dimerMomZ (dimerPotZ J h r η) p T 0) muZ :=
  Integrable.mono' (integrable_const (∫ t in Set.Ioi (0 : ℝ), DimerInfVolAux.momBound p T t))
    (measurable_M J h r p T hp hT).aestronglyMeasurable
    (ae_of_all _ fun η => norm_dimerMomZ_le _ _ (DimerInfVolAux.abs_dimerPotZ_le J h r η) p T hp hT)

/-- `B` for `v^{(r)}` depends only on the signs in `sitesS r L`. -/
lemma Bfun_dep (J h : ℝ) (r : ℤ) (p T δ lam : ℝ) (L : ℕ) (η η' : CfgZ)
    (hη : ∀ k ∈ sitesS r L, η k = η' k) :
    Bfun (dimerPotZ J h r η) p T δ lam L = Bfun (dimerPotZ J h r η') p T δ lam L := by
  unfold Bfun
  rw [Yfun_congr _ _ p T L fun y hy1 hy2 => dimerPotZ_congr J h r L η η' hη y hy1 hy2]

/-- `η ↦ B` is integrable. -/
lemma integrable_Bfun (J h : ℝ) (r : ℤ) (p T δ lam : ℝ) (L : ℕ) :
    Integrable (fun η => Bfun (dimerPotZ J h r η) p T δ lam L) muZ :=
  integrable_dep (sitesS r L) _ fun η η' hη => Bfun_dep J h r p T δ lam L η η' hη

/-- **`E B ≤ 2δ (∑_{j ≤ L} j^p) C'`** if `E X_j(E + i/(2T)) ≤ C'` for `1 ≤ j ≤ L` and
`E ∈ [λ - δ, λ + δ]` (Fubini). -/
lemma integral_Bfun_le (J h : ℝ) (r : ℤ) (p T δ lam C' : ℝ) (hδ : 0 ≤ δ) (L : ℕ)
    (hX : ∀ E ∈ Set.Icc (lam - δ) (lam + δ), ∀ j : ℕ, 1 ≤ j → j ≤ L →
      Integrable (fun η => Xj (dimerPotZ J h r η) ((E : ℂ) + Complex.I / (2 * T)) j) muZ ∧
      ∫ η, Xj (dimerPotZ J h r η) ((E : ℂ) + Complex.I / (2 * T)) j ∂muZ ≤ C') :
    ∫ η, Bfun (dimerPotZ J h r η) p T δ lam L ∂muZ ≤
      2 * δ * (∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p) * C' := by
  unfold Bfun
  rw [integral_swap_dep (sitesS r L) (fun η E => Yfun (dimerPotZ J h r η) p T L E)
    (fun η η' hηη' => Yfun_congr _ _ p T L fun y hy1 hy2 =>
      dimerPotZ_congr J h r L η η' hηη' y hy1 hy2)
    _ (fun η => (continuous_Yfun _ p T L).integrableOn_Icc)]
  have hpt : ∀ E ∈ Set.Icc (lam - δ) (lam + δ), ∫ η, Yfun (dimerPotZ J h r η) p T L E ∂muZ ≤
      (∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p) * C' := by
    intro E hE
    unfold Yfun
    rw [integral_finsetSum _ fun i hi =>
      ((hX E hE (i + 1) (by omega) (by rw [Finset.mem_range] at hi; omega)).1.const_mul _),
      Finset.sum_mul]
    refine Finset.sum_le_sum fun i hi => ?_
    rw [integral_const_mul]
    exact mul_le_mul_of_nonneg_left
      (hX E hE (i + 1) (by omega) (by rw [Finset.mem_range] at hi; omega)).2 (by positivity)
  have hvol : volume.real (Set.Icc (lam - δ) (lam + δ)) = 2 * δ := by
    rw [Measure.real, Real.volume_Icc, ENNReal.toReal_ofReal (by linarith)]
    ring
  calc ∫ E in Set.Icc (lam - δ) (lam + δ), ∫ η, Yfun (dimerPotZ J h r η) p T L E ∂muZ
      ≤ ∫ E in Set.Icc (lam - δ) (lam + δ), (∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p) * C' :=
        integral_mono_of_nonneg (ae_of_all _ fun E => integral_nonneg fun η => Yfun_nonneg _ _ _ _ _)
          continuous_const.integrableOn_Icc (ae_restrict_of_forall_mem measurableSet_Icc hpt)
    _ = 2 * δ * (∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p) * C' := by
        rw [setIntegral_const, hvol, smul_eq_mul]
        ring

/-! ### The second moment of `X_j` (paper: `lem:dimersecond`) -/

/-- `Re(E + i/(2T)) = E`. -/
lemma re_add_I_div (E T : ℝ) : ((E : ℂ) + Complex.I / (2 * T)).re = E := by
  rw [show (Complex.I / (2 * T) : ℂ) = ((1 / (2 * T) : ℝ) : ℂ) * Complex.I by push_cast; ring]
  simp

/-- `|E + i/(2T) - λ| ≤ |E - λ| + 1/(2T)`. -/
lemma norm_z_sub_le' (E lam T : ℝ) (hT : 0 < T) :
    ‖((E : ℂ) + Complex.I / (2 * T)) - lam‖ ≤ |E - lam| + 1 / (2 * T) := by
  have e : ((E : ℂ) + Complex.I / (2 * T)) - lam =
      ((E - lam : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I := by
    push_cast
    ring
  rw [e]
  calc ‖((E - lam : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I‖
      ≤ ‖((E - lam : ℝ) : ℂ)‖ + ‖((1 / (2 * T) : ℝ) : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |E - lam| + 1 / (2 * T) := by
        rw [Complex.norm_real, norm_mul, Complex.norm_real, Complex.norm_I, mul_one,
          Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos (by positivity : (0 : ℝ) < 1 / (2 * T))]

/-- **`E X_j ≤ C' = 2C e^{3K/2}`** (paper: proof of `prop:dimersuper`(a), from `lem:dimersecond`): for
`T ≥ max(1, 4/δ₀²)`, `|E - λ| ≤ T^{-1/2}` and `j ≤ T`. -/
lemma integral_Xj_le (J h : ℝ) (r : ℤ) (C K δ₀ : ℝ) (hC : 0 < C) (hK : 0 ≤ K) (hδ₀ : 0 < δ₀)
    (hbb : ∀ z : ℂ, ‖z - ((2 * |h| / J : ℝ) : ℂ)‖ ≤ δ₀ → ∀ (m : ℤ) (L : ℕ),
      Integrable (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) muZ ∧
      ∫ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ∂muZ ≤
        C * Real.exp (K * L * ((z.re - 2 * |h| / J) ^ 2 + |z.im|)))
    (T : ℝ) (hT1 : 1 ≤ T) (hTδ : 4 / δ₀ ^ 2 ≤ T) (E : ℝ)
    (hE : |E - 2 * |h| / J| ≤ 1 / Real.sqrt T) (j : ℕ) (hj : (j : ℝ) ≤ T) :
    Integrable (fun η => Xj (dimerPotZ J h r η) ((E : ℂ) + Complex.I / (2 * T)) j) muZ ∧
    ∫ η, Xj (dimerPotZ J h r η) ((E : ℂ) + Complex.I / (2 * T)) j ∂muZ ≤
      2 * C * Real.exp (3 * K / 2) := by
  have hT0 : 0 < T := by linarith
  set s := Real.sqrt T with hs
  have hs0 : 0 < s := Real.sqrt_pos.mpr hT0
  have hs2 : s ^ 2 = T := Real.sq_sqrt hT0.le
  have hs1 : 1 ≤ s := by
    have h := Real.sqrt_le_sqrt hT1
    rwa [Real.sqrt_one] at h
  set z := (E : ℂ) + Complex.I / (2 * T) with hz
  set lam := 2 * |h| / J with hlam
  have hzδ : ‖z - ((lam : ℝ) : ℂ)‖ ≤ δ₀ := by
    have h1 := norm_z_sub_le' E lam T hT0
    have h2 : 1 / (2 * T) ≤ 1 / s := by
      rw [div_le_div_iff₀ (by positivity) hs0]
      nlinarith
    have h3 : 2 / s ≤ δ₀ := by
      rw [div_le_iff₀ hs0]
      have h4 : 4 ≤ (δ₀ * s) ^ 2 := by
        rw [mul_pow, hs2]
        rw [div_le_iff₀ (by positivity)] at hTδ
        linarith
      nlinarith [mul_pos hδ₀ hs0]
    calc ‖z - ((lam : ℝ) : ℂ)‖ ≤ |E - lam| + 1 / (2 * T) := h1
      _ ≤ 1 / s + 1 / s := add_le_add hE h2
      _ = 2 / s := by ring
      _ ≤ δ₀ := h3
  obtain ⟨hi1, hb1⟩ := hbb z hzδ 1 j
  obtain ⟨hi2, hb2⟩ := hbb z hzδ (-(j : ℤ)) j
  have hexp : K * (j : ℝ) * ((z.re - lam) ^ 2 + |z.im|) ≤ 3 * K / 2 := by
    rw [hz, re_add_I_div, DimerLowerAux.im_add_I_div,
      abs_of_pos (by positivity : (0 : ℝ) < 1 / (2 * T))]
    have h5 : (E - lam) ^ 2 ≤ 1 / T := by
      have := pow_le_pow_left₀ (abs_nonneg _) hE 2
      rw [sq_abs, div_pow, one_pow, hs2] at this
      exact this
    have h6 : (j : ℝ) * ((E - lam) ^ 2 + 1 / (2 * T)) ≤ T * (1 / T + 1 / (2 * T)) :=
      mul_le_mul hj (by linarith) (by positivity) hT0.le
    have h7 : T * (1 / T + 1 / (2 * T)) = 3 / 2 := by
      field_simp
      ring
    calc K * (j : ℝ) * ((E - lam) ^ 2 + 1 / (2 * T))
        = K * ((j : ℝ) * ((E - lam) ^ 2 + 1 / (2 * T))) := by ring
      _ ≤ K * (3 / 2) := by rw [← h7]; exact mul_le_mul_of_nonneg_left h6 hK
      _ = 3 * K / 2 := by ring
  have hb : C * Real.exp (K * (j : ℝ) * ((z.re - lam) ^ 2 + |z.im|)) ≤
      C * Real.exp (3 * K / 2) :=
    mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr hexp) hC.le
  refine ⟨hi1.add hi2, ?_⟩
  unfold Xj
  rw [integral_add hi1 hi2]
  linarith

/-- **`E B ≤ 2δ (∑_{j ≤ L} j^p) C'`** for `δ = T^{-1/2}`, `L = ⌊T⌋`, `T ≥ max(4, 4/δ₀²)`. -/
lemma integral_Bfun_le' (J h : ℝ) (r : ℤ) (p : ℝ) (C K δ₀ : ℝ) (hC : 0 < C) (hK : 0 ≤ K)
    (hδ₀ : 0 < δ₀)
    (hbb : ∀ z : ℂ, ‖z - ((2 * |h| / J : ℝ) : ℂ)‖ ≤ δ₀ → ∀ (m : ℤ) (L : ℕ),
      Integrable (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) muZ ∧
      ∫ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ∂muZ ≤
        C * Real.exp (K * L * ((z.re - 2 * |h| / J) ^ 2 + |z.im|)))
    (T : ℝ) (hT4 : 4 ≤ T) (hTδ : 4 / δ₀ ^ 2 ≤ T) :
    ∫ η, Bfun (dimerPotZ J h r η) p T (1 / Real.sqrt T) (2 * |h| / J) ⌊T⌋₊ ∂muZ ≤
      2 * (1 / Real.sqrt T) * (∑ i ∈ Finset.range ⌊T⌋₊, ((i : ℝ) + 1) ^ p) *
        (2 * C * Real.exp (3 * K / 2)) := by
  have hT0 : 0 < T := by linarith
  refine integral_Bfun_le J h r p T _ _ _ (by positivity) _ fun E hE j hj1 hjL => ?_
  exact integral_Xj_le J h r C K δ₀ hC hK hδ₀ hbb T (by linarith) hTδ E
    (abs_le.mpr ⟨by linarith [hE.1], by linarith [hE.2]⟩) j
    (le_trans (Nat.cast_le.mpr hjL) (Nat.floor_le hT0.le))

/-- `√T ≥ 2` for `T ≥ 4`. -/
lemma two_le_sqrt (T : ℝ) (hT : 4 ≤ T) : 2 ≤ Real.sqrt T := by
  have h := Real.sqrt_le_sqrt hT
  rwa [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)] at h

/-- **The constant** (paper: proof of `prop:dimersuper`(a)): for `T ≥ 2` and `δ = T^{-1/2}`,
`(2πT)⁻¹ (c₀²/2) 2δ ∑_{j ≤ ⌊T⌋} j^p / C' ≥ c₀²/(2^{p+2} π (p+1) C') T^{p - 1/2}`. -/
lemma lower_const (lam p T C' : ℝ) (hp : 0 < p) (hT : 2 ≤ T) (hC' : 0 < C') (hl : 0 ≤ lam) :
    (1 / (3 + 2 * lam)) ^ 2 / (2 ^ (p + 2) * Real.pi * (p + 1) * C') * T ^ (p - 1 / 2) ≤
      1 / (2 * Real.pi * T) * ((1 / (3 + 2 * lam)) ^ 2 / 2 *
        (2 * (1 / Real.sqrt T) * ∑ i ∈ Finset.range ⌊T⌋₊, ((i : ℝ) + 1) ^ p)) / C' := by
  have hT0 : 0 < T := by linarith
  set s := Real.sqrt T with hs
  have hs0 : 0 < s := Real.sqrt_pos.mpr hT0
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  have hc₀0 : 0 < c₀ := by positivity
  set L := ⌊T⌋₊ with hL
  set SL := ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p with hSL
  have hLT2 : T / 2 ≤ L := by
    have := Nat.lt_floor_add_one T
    linarith
  have hSLge : (T / 2) ^ (p + 1) / (p + 1) ≤ SL :=
    (div_le_div_of_nonneg_right (Real.rpow_le_rpow (by positivity) hLT2 (by linarith))
      (by linarith)).trans (DimerLowerAux.sum_rpow_ge p hp L)
  have e0 : T ^ (p + 1) = T ^ (p - 1 / 2) * T * s := by
    have h1 : T ^ (p + 1) = T ^ ((p - 1 / 2) + 1 + 1 / 2) := by
      congr 1
      ring
    rw [h1, Real.rpow_add hT0, Real.rpow_add hT0, Real.rpow_one, hs, Real.sqrt_eq_rpow]
  have e1 : (T / 2) ^ (p + 1) = T ^ (p - 1 / 2) * T * s / 2 ^ (p + 1) := by
    rw [Real.div_rpow hT0.le (by norm_num), e0]
  have e2 : (2 : ℝ) ^ (p + 2) = 2 ^ (p + 1) * 2 := by
    rw [show p + 2 = (p + 1) + 1 by ring, Real.rpow_add two_pos, Real.rpow_one]
  have h2p : (0 : ℝ) < 2 ^ (p + 1) := by positivity
  have hpi : 0 < Real.pi := Real.pi_pos
  calc c₀ ^ 2 / (2 ^ (p + 2) * Real.pi * (p + 1) * C') * T ^ (p - 1 / 2)
      = c₀ ^ 2 * ((T / 2) ^ (p + 1) / (p + 1)) / (2 * Real.pi * T * s * C') := by
        rw [e1, e2]
        field_simp
    _ ≤ c₀ ^ 2 * SL / (2 * Real.pi * T * s * C') := by gcongr
    _ = 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * (1 / s) * SL)) / C' := by
        field_simp

/-! ### Almost-sure superdiffusion (paper: proof of `prop:dimersuper`(b)) -/

/-- `B ≥ 0`. -/
lemma Bfun_nonneg (v : ℤ → ℝ) (p T δ lam : ℝ) (L : ℕ) : 0 ≤ Bfun v p T δ lam L :=
  setIntegral_nonneg measurableSet_Icc fun E _ => Yfun_nonneg v p T L E

/-- **Dyadic interpolation**: `M^ℤ(T₁) ≤ 2 M^ℤ(T)` for `T₁ ≤ T ≤ 2T₁`, since
`T₁⁻¹ e^{-t/T₁} ≤ 2 T⁻¹ e^{-t/T}` for `t > 0`. -/
lemma dimerMomZ_le_two_mul (v : ℤ → ℝ) (W : ℝ) (hW : ∀ y, |v y| ≤ W) (p T₁ T : ℝ) (hp : 0 < p)
    (hT₁ : 0 < T₁) (h1 : T₁ ≤ T) (h2 : T ≤ 2 * T₁) :
    dimerMomZ v p T₁ 0 ≤ 2 * dimerMomZ v p T 0 := by
  have hT : 0 < T := by linarith
  have hlim : ∀ T' : ℝ, 0 < T' → Tendsto (fun n : ℕ => dimerMom (-(n : ℤ)) n v p T'
      (origin (-(n : ℤ)) n (by omega) (by omega))) atTop (𝓝 (dimerMomZ v p T' 0)) :=
    fun T' hT' => tendsto_dimerMom_Z v W hW (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ))
      (fun n => by simp) (fun n => by positivity)
      (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop p T'
      hp hT'
  refine le_of_tendsto_of_tendsto' (hlim T₁ hT₁) ((hlim T hT).const_mul 2) fun n => ?_
  have ha0 : -(n : ℤ) ≤ 0 := by omega
  have hb0 : (0 : ℤ) ≤ n := by omega
  have hint : Integrable (fun t => (1 / T) * Real.exp (-t / T) *
      ∑ y : Set.Icc (-(n : ℤ)) n, |((y.1 - (origin (-(n : ℤ)) n ha0 hb0).1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH (-(n : ℤ)) n v) t y (origin (-(n : ℤ)) n ha0 hb0)‖ ^ 2)
      (volume.restrict (Set.Ioi 0)) :=
    Integrable.mono' (DimerInfVolAux.integrableOn_momBound p T hp hT)
      (DimerInfVolAux.continuous_integrand _ _ _ T).aestronglyMeasurable
      (ae_of_all _ fun t => DimerInfVolAux.integrand_le _ _ ha0 hb0 v p T hp hT t)
  unfold dimerMom
  rw [← integral_const_mul]
  refine integral_mono_of_nonneg (ae_of_all _ fun t => ?_) (hint.const_mul 2)
    (ae_restrict_of_forall_mem measurableSet_Ioi fun t ht => ?_)
  · have : 0 ≤ ∑ y : Set.Icc (-(n : ℤ)) n, |((y.1 - (origin (-(n : ℤ)) n ha0 hb0).1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH (-(n : ℤ)) n v) t y (origin (-(n : ℤ)) n ha0 hb0)‖ ^ 2 :=
      Finset.sum_nonneg fun y _ => by positivity
    positivity
  · have ht0 : 0 < t := ht
    set S := ∑ y : Set.Icc (-(n : ℤ)) n, |((y.1 - (origin (-(n : ℤ)) n ha0 hb0).1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH (-(n : ℤ)) n v) t y (origin (-(n : ℤ)) n ha0 hb0)‖ ^ 2 with hS
    have hS0 : 0 ≤ S := Finset.sum_nonneg fun y _ => by positivity
    have e1 : 1 / T₁ ≤ 2 * (1 / T) := by
      rw [div_le_iff₀ hT₁]
      field_simp
      linarith
    have e2 : Real.exp (-t / T₁) ≤ Real.exp (-t / T) := by
      apply Real.exp_le_exp.mpr
      rw [neg_div, neg_div, neg_le_neg_iff]
      exact div_le_div_of_nonneg_left ht0.le hT₁ h1
    calc 1 / T₁ * Real.exp (-t / T₁) * S ≤ 2 * (1 / T) * Real.exp (-t / T) * S := by gcongr
      _ = 2 * (1 / T * Real.exp (-t / T) * S) := by ring

/-- **From the dyadic times to all times**: if `M(2^N) ≥ c (2^N)^q/(N+1)²` for `N ≥ N₀` and
`M(T₁) ≤ 2M(T)` for `T₁ ≤ T ≤ 2T₁`, then `M(T) ≥ T^{q - α}` for all large `T`. -/
lemma eventually_ge_of_dyadic (M : ℝ → ℝ) (q c α : ℝ) (hc : 0 < c) (hα : 0 < α) (N₀ : ℕ)
    (hdy : ∀ N : ℕ, N₀ ≤ N → c * ((2 : ℝ) ^ N) ^ q / ((N : ℝ) + 1) ^ 2 ≤ M (2 ^ N))
    (hint : ∀ T₁ T : ℝ, 0 < T₁ → T₁ ≤ T → T ≤ 2 * T₁ → M T₁ ≤ 2 * M T) :
    ∀ᶠ T in atTop, T ^ (q - α) ≤ M T := by
  set D := max 1 ((2 : ℝ) ^ q) with hD
  have hD1 : 1 ≤ D := le_max_left _ _
  set r := (2 : ℝ) ^ α with hr
  have hr1 : 1 < r := Real.one_lt_rpow (by norm_num) hα
  have hlim := tendsto_pow_const_div_const_pow_of_one_lt 2 hr1
  have hev : ∀ᶠ N : ℕ in atTop, (N : ℝ) ^ 2 / r ^ N < c / (8 * D) :=
    (tendsto_order.1 hlim).2 _ (by positivity)
  obtain ⟨N₃, hN₃⟩ := eventually_atTop.mp hev
  set N₂ := max (max N₀ N₃) 1 with hN₂
  filter_upwards [eventually_ge_atTop ((2 : ℝ) ^ N₂)] with T hT
  have hT1 : 1 ≤ T := le_trans (one_le_pow₀ (by norm_num)) hT
  have hT0 : 0 < T := by linarith
  set N := Nat.log 2 ⌊T⌋₊ with hN
  have hfl : 1 ≤ ⌊T⌋₊ := Nat.le_floor (by exact_mod_cast hT1)
  have hlow : (2 : ℝ) ^ N ≤ T := by
    have h1 : 2 ^ N ≤ ⌊T⌋₊ := Nat.pow_log_le_self 2 (by omega)
    calc (2 : ℝ) ^ N = ((2 ^ N : ℕ) : ℝ) := by push_cast; ring
      _ ≤ (⌊T⌋₊ : ℝ) := by exact_mod_cast h1
      _ ≤ T := Nat.floor_le hT0.le
  have hup : T ≤ 2 * 2 ^ N := by
    have h1 : ⌊T⌋₊ < 2 ^ (N + 1) := Nat.lt_pow_succ_log_self (by norm_num) _
    have h2 : T < ⌊T⌋₊ + 1 := Nat.lt_floor_add_one T
    have h3 : ((⌊T⌋₊ : ℕ) : ℝ) + 1 ≤ ((2 ^ (N + 1) : ℕ) : ℝ) := by exact_mod_cast h1
    push_cast at h3
    rw [pow_succ] at h3
    linarith
  have hNge : N₂ ≤ N := by
    have h1 : 2 ^ N₂ ≤ ⌊T⌋₊ := Nat.le_floor (by push_cast; exact hT)
    exact Nat.le_log_of_pow_le (by norm_num) h1
  have hN0 : N₀ ≤ N := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hNge
  have hN3 : N₃ ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hNge
  have hN1 : 1 ≤ N := le_trans (le_max_right _ _) hNge
  set TN := (2 : ℝ) ^ N with hTN
  have hTN0 : 0 < TN := by positivity
  have hM1 := hdy N hN0
  have hM2 := hint TN T hTN0 hlow hup
  have hTq : T ^ q ≤ D * TN ^ q := by
    rcases le_or_gt 0 q with hq | hq
    · calc T ^ q ≤ (2 * TN) ^ q := Real.rpow_le_rpow hT0.le hup hq
        _ = 2 ^ q * TN ^ q := Real.mul_rpow (by norm_num) hTN0.le
        _ ≤ D * TN ^ q := by gcongr; exact le_max_right _ _
    · calc T ^ q ≤ TN ^ q := Real.rpow_le_rpow_of_nonpos hTN0 hlow hq.le
        _ = 1 * TN ^ q := (one_mul _).symm
        _ ≤ D * TN ^ q := by gcongr
  have hTα : TN ^ α ≤ T ^ α := Real.rpow_le_rpow hTN0.le hlow hα.le
  have hTNα : TN ^ α = r ^ N := by
    rw [hTN, ← Real.rpow_natCast, ← Real.rpow_mul (by norm_num), mul_comm,
      Real.rpow_mul (by norm_num), Real.rpow_natCast]
  have h3 := hN₃ N hN3
  have hrN : 0 < r ^ N := by positivity
  have h4 : 2 * D * ((N : ℝ) + 1) ^ 2 ≤ c * r ^ N := by
    rw [div_lt_div_iff₀ hrN (by positivity)] at h3
    have h5 : ((N : ℝ) + 1) ^ 2 ≤ 4 * (N : ℝ) ^ 2 := by
      have : (1 : ℝ) ≤ N := by exact_mod_cast hN1
      nlinarith
    nlinarith
  have hTq0 : 0 < TN ^ q := by positivity
  calc T ^ (q - α) = T ^ q / T ^ α := Real.rpow_sub hT0 q α
    _ ≤ D * TN ^ q / TN ^ α := by gcongr
    _ ≤ c * TN ^ q / ((N : ℝ) + 1) ^ 2 / 2 := by
        rw [hTNα, div_div, div_le_div_iff₀ hrN (by positivity)]
        calc D * TN ^ q * (((N : ℝ) + 1) ^ 2 * 2) = TN ^ q * (2 * D * ((N : ℝ) + 1) ^ 2) := by
              ring
          _ ≤ TN ^ q * (c * r ^ N) := by gcongr
          _ = c * TN ^ q * r ^ N := by ring
    _ ≤ M TN / 2 := by gcongr
    _ ≤ M T := by linarith

end DimerSuperAux

/-- **Superdiffusion on average** (paper: Proposition `prop:dimersuper`(a)):
`E M^ℤ_{v^{(r)},p}(T; 0) ≥ c T^{p - 1/2}` for `T ≥ T₀`. -/
theorem dimer_super_annealed (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ)
    (p : ℝ) (hp : 0 < p) :
    ∃ c T₀ : ℝ, 0 < c ∧ 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ T →
      c * T ^ (p - 1 / 2) ≤ ∫ η, dimerMomZ (dimerPotZ J h r η) p T 0 ∂muZ := by
  open DimerSuperAux in
  obtain ⟨C, K, δ₀, hC, hK, hδ₀, hbb⟩ := integral_hs2_tmProd_le J h hJ hh0 hh1 r
  have hl0 : 0 < 2 * |h| / J := by positivity
  set lam := 2 * |h| / J with hlam
  set C' := 2 * C * Real.exp (3 * K / 2) with hC'
  have hC'0 : 0 < C' := by positivity
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  refine ⟨c₀ ^ 2 / (2 ^ (p + 2) * Real.pi * (p + 1) * C'), max 4 (4 / δ₀ ^ 2), by positivity,
    by positivity, fun T hT => ?_⟩
  have hT4 : 4 ≤ T := le_trans (le_max_left _ _) hT
  have hTδ : 4 / δ₀ ^ 2 ≤ T := le_trans (le_max_right _ _) hT
  have hT0 : 0 < T := by linarith
  set δ := 1 / Real.sqrt T with hδ
  have hδ0 : 0 < δ := by positivity
  have hδ1 : δ ≤ 1 / 2 := one_div_le_one_div_of_le (by norm_num) (two_le_sqrt T hT4)
  set L := ⌊T⌋₊ with hL
  set SL := ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p with hSL
  have hEB : ∫ η, Bfun (dimerPotZ J h r η) p T δ lam L ∂muZ ≤ 2 * δ * SL * C' :=
    integral_Bfun_le' J h r p C K δ₀ hC hK hδ₀ hbb T hT4 hTδ
  have hpt : ∀ η, 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / C') * SL -
      Bfun (dimerPotZ J h r η) p T δ lam L / C' ^ 2)) ≤ dimerMomZ (dimerPotZ J h r η) p T 0 :=
    fun η => det_bound lam hl0 (dimerPotZ J h r η) |2 * h / J|
      (DimerInfVolAux.abs_dimerPotZ_le J h r η) ((isDimerCfgZ_dimerPotZ J h hJ r η).1 0) p T δ hp
      (by linarith) hδ0 hδ1 L C' hC'0
  have hBint := integrable_Bfun J h r p T δ lam L
  have hlow_int : Integrable (fun η => 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 *
      (2 * δ * (2 / C') * SL - Bfun (dimerPotZ J h r η) p T δ lam L / C' ^ 2))) muZ :=
    (((integrable_const _).sub (hBint.div_const _)).const_mul _).const_mul _
  have h1 := integral_mono hlow_int (integrable_M J h r p T hp hT0) hpt
  have h2 : ∫ η, 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 *
      (2 * δ * (2 / C') * SL - Bfun (dimerPotZ J h r η) p T δ lam L / C' ^ 2)) ∂muZ =
      1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / C') * SL -
        (∫ η, Bfun (dimerPotZ J h r η) p T δ lam L ∂muZ) / C' ^ 2)) := by
    rw [integral_const_mul, integral_const_mul,
      integral_sub (integrable_const _) (hBint.div_const _), integral_const, integral_div]
    simp
  have h3 := lower_const lam p T C' hp (by linarith) hC'0 hl0.le
  calc c₀ ^ 2 / (2 ^ (p + 2) * Real.pi * (p + 1) * C') * T ^ (p - 1 / 2)
      ≤ 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * SL)) / C' := h3
    _ = 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 *
          (2 * δ * (2 / C') * SL - 2 * δ * SL * C' / C' ^ 2)) := by
        field_simp
        ring
    _ ≤ 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / C') * SL -
          (∫ η, Bfun (dimerPotZ J h r η) p T δ lam L ∂muZ) / C' ^ 2)) := by
        gcongr
    _ = _ := h2.symm
    _ ≤ _ := h1

/-- **Superdiffusion almost surely** (paper: Proposition `prop:dimersuper`(b)): almost surely, for every
`α > 0`, `M^ℤ_{v^{(r)},p}(T; 0) ≥ T^{p - 1/2 - α}` for all large `T`. -/
theorem dimer_super_ae (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) :
    ∀ᵐ η ∂muZ, ∀ α : ℝ, 0 < α → ∀ᶠ T in atTop,
      T ^ (p - 1 / 2 - α) ≤ dimerMomZ (dimerPotZ J h r η) p T 0 := by
  open DimerSuperAux in
  obtain ⟨C, K, δ₀, hC, hK, hδ₀, hbb⟩ := integral_hs2_tmProd_le J h hJ hh0 hh1 r
  have hl0 : 0 < 2 * |h| / J := by positivity
  set lam := 2 * |h| / J with hlam
  set C' := 2 * C * Real.exp (3 * K / 2) with hC'
  have hC'0 : 0 < C' := by positivity
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  set c := c₀ ^ 2 / (2 ^ (p + 2) * Real.pi * (p + 1) * C') with hc
  have hc0 : 0 < c := by positivity
  -- the dyadic times `T_N = 2^N ≥ max(4, 4/δ₀²)` for `N ≥ N₀`
  obtain ⟨N₀, hN₀⟩ := pow_unbounded_of_one_lt (max 4 (4 / δ₀ ^ 2)) (by norm_num : (1 : ℝ) < 2)
  have hTN : ∀ N : ℕ, N₀ ≤ N → 4 ≤ (2 : ℝ) ^ N ∧ 4 / δ₀ ^ 2 ≤ (2 : ℝ) ^ N := fun N hN => by
    have : (2 : ℝ) ^ N₀ ≤ 2 ^ N := pow_le_pow_right₀ (by norm_num) hN
    exact ⟨by linarith [le_max_left 4 (4 / δ₀ ^ 2)], by linarith [le_max_right 4 (4 / δ₀ ^ 2)]⟩
  -- `A_N = 2δ_N ∑_{j ≤ L_N} j^p`, `B_N`, and the bad events `{B_N ≥ (N+1)² C' A_N}`
  set A : ℕ → ℝ := fun N => 2 * (1 / Real.sqrt (2 ^ N)) *
    ∑ i ∈ Finset.range ⌊(2 : ℝ) ^ N⌋₊, ((i : ℝ) + 1) ^ p with hA
  set B : ℕ → CfgZ → ℝ := fun N η =>
    Bfun (dimerPotZ J h r η) p (2 ^ N) (1 / Real.sqrt (2 ^ N)) lam ⌊(2 : ℝ) ^ N⌋₊ with hB
  set Bad : ℕ → Set CfgZ := fun N => {η | ((N : ℝ) + 1) ^ 2 * C' * A N ≤ B N η} with hBad
  have hApos : ∀ N, N₀ ≤ N → 0 < A N := by
    intro N hN
    have h4 := (hTN N hN).1
    have hL1 : 1 ≤ ⌊(2 : ℝ) ^ N⌋₊ := Nat.le_floor (by push_cast; linarith)
    have hS : 0 < ∑ i ∈ Finset.range ⌊(2 : ℝ) ^ N⌋₊, ((i : ℝ) + 1) ^ p :=
      Finset.sum_pos (fun i _ => by positivity) ⟨0, Finset.mem_range.mpr (by omega)⟩
    simp only [hA]
    positivity
  -- Markov's inequality
  have hmeas : ∀ N, N₀ ≤ N → muZ.real (Bad N) ≤ 1 / ((N : ℝ) + 1) ^ 2 := by
    intro N hN
    obtain ⟨h4, hδ⟩ := hTN N hN
    have hEB : ∫ η, B N η ∂muZ ≤ A N * C' :=
      integral_Bfun_le' J h r p C K δ₀ hC hK hδ₀ hbb (2 ^ N) h4 hδ
    have hM := mul_meas_ge_le_integral_of_nonneg
      (ae_of_all _ fun η => Bfun_nonneg (dimerPotZ J h r η) p (2 ^ N) (1 / Real.sqrt (2 ^ N)) lam
        ⌊(2 : ℝ) ^ N⌋₊)
      (integrable_Bfun J h r p (2 ^ N) (1 / Real.sqrt (2 ^ N)) lam ⌊(2 : ℝ) ^ N⌋₊)
      (((N : ℝ) + 1) ^ 2 * C' * A N)
    have hpos : 0 < C' * A N := mul_pos hC'0 (hApos N hN)
    have key : (C' * A N) * (muZ.real (Bad N) * ((N : ℝ) + 1) ^ 2) ≤ (C' * A N) * 1 := by
      calc (C' * A N) * (muZ.real (Bad N) * ((N : ℝ) + 1) ^ 2)
          = ((N : ℝ) + 1) ^ 2 * C' * A N * muZ.real (Bad N) := by ring
        _ ≤ ∫ η, B N η ∂muZ := hM
        _ ≤ A N * C' := hEB
        _ = (C' * A N) * 1 := by ring
    rw [le_div_iff₀ (by positivity)]
    exact le_of_mul_le_mul_left key hpos
  have hmeas' : ∀ N, muZ.real (Bad N) ≤ ((N₀ : ℝ) + 1) ^ 2 / ((N : ℝ) + 1) ^ 2 := by
    intro N
    rcases le_or_gt N₀ N with hN | hN
    · refine (hmeas N hN).trans (div_le_div_of_nonneg_right ?_ (by positivity))
      have : (0 : ℝ) ≤ N₀ := Nat.cast_nonneg _
      nlinarith
    · refine measureReal_le_one.trans ?_
      rw [one_le_div (by positivity)]
      have : (N : ℝ) + 1 ≤ (N₀ : ℝ) + 1 := by exact_mod_cast (by omega : N + 1 ≤ N₀ + 1)
      gcongr
  -- Borel–Cantelli
  have hsum : ∑' N, muZ (Bad N) ≠ ⊤ := by
    have h1 : Summable (fun N : ℕ => 1 / ((N : ℝ) + 1) ^ 2) := by
      have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
      simpa using this
    have hg : Summable (fun N : ℕ => ((N₀ : ℝ) + 1) ^ 2 / ((N : ℝ) + 1) ^ 2) :=
      (h1.mul_left (((N₀ : ℝ) + 1) ^ 2)).congr fun N => by ring
    refine ne_top_of_le_ne_top (ENNReal.ofReal_ne_top
      (r := ∑' N : ℕ, ((N₀ : ℝ) + 1) ^ 2 / ((N : ℝ) + 1) ^ 2)) ?_
    rw [ENNReal.ofReal_tsum_of_nonneg (fun N => by positivity) hg]
    refine tsum_le_tsum fun N => ?_
    rw [← ofReal_measureReal (measure_ne_top muZ _)]
    exact ENNReal.ofReal_le_ofReal (hmeas' N)
  filter_upwards [ae_eventually_notMem hsum] with η hη
  intro α hα
  obtain ⟨N₁, hN₁⟩ := eventually_atTop.mp hη
  refine eventually_ge_of_dyadic (fun T => dimerMomZ (dimerPotZ J h r η) p T 0) (p - 1 / 2) c α
    hc0 hα (max N₀ N₁) (fun N hN => ?_) (fun T₁ T hT₁ h1 h2 => dimerMomZ_le_two_mul _ |2 * h / J|
      (DimerInfVolAux.abs_dimerPotZ_le J h r η) p T₁ T hp hT₁ h1 h2)
  -- the dyadic bound `M(2^N) ≥ c 2^{N(p - 1/2)}/(N+1)²` off the bad events
  show c * ((2 : ℝ) ^ N) ^ (p - 1 / 2) / ((N : ℝ) + 1) ^ 2 ≤
    dimerMomZ (dimerPotZ J h r η) p (2 ^ N) 0
  have hN0 : N₀ ≤ N := le_trans (le_max_left _ _) hN
  have hN1 : N₁ ≤ N := le_trans (le_max_right _ _) hN
  obtain ⟨h4, hδ⟩ := hTN N hN0
  have hnot : η ∉ Bad N := hN₁ N hN1
  have hBlt : B N η < ((N : ℝ) + 1) ^ 2 * C' * A N := not_le.mp hnot
  have hT0 : 0 < (2 : ℝ) ^ N := by positivity
  have hδ1 : 1 / Real.sqrt (2 ^ N) ≤ 1 / 2 :=
    one_div_le_one_div_of_le (by norm_num) (two_le_sqrt _ h4)
  have hdet := det_bound lam hl0 (dimerPotZ J h r η) |2 * h / J|
    (DimerInfVolAux.abs_dimerPotZ_le J h r η) ((isDimerCfgZ_dimerPotZ J h hJ r η).1 0) p (2 ^ N)
    (1 / Real.sqrt (2 ^ N)) hp (by linarith) (by positivity) hδ1 ⌊(2 : ℝ) ^ N⌋₊
    (((N : ℝ) + 1) ^ 2 * C') (by positivity)
  have hlc := lower_const lam p (2 ^ N) C' hp (by linarith) hC'0 hl0.le
  set T := (2 : ℝ) ^ N with hT
  set SL := ∑ i ∈ Finset.range ⌊T⌋₊, ((i : ℝ) + 1) ^ p with hSL
  set δ := 1 / Real.sqrt T with hδdef
  have hA' : A N = 2 * δ * SL := rfl
  have hB' : B N η = Bfun (dimerPotZ J h r η) p T δ lam ⌊T⌋₊ := rfl
  rw [hA', hB'] at hBlt
  have hN2 : 0 < ((N : ℝ) + 1) ^ 2 := by positivity
  calc c * T ^ (p - 1 / 2) / ((N : ℝ) + 1) ^ 2
      ≤ (1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * SL)) / C') / ((N : ℝ) + 1) ^ 2 :=
        div_le_div_of_nonneg_right hlc hN2.le
    _ = 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / (((N : ℝ) + 1) ^ 2 * C')) * SL -
          ((N : ℝ) + 1) ^ 2 * C' * (2 * δ * SL) / (((N : ℝ) + 1) ^ 2 * C') ^ 2)) := by
        field_simp
        ring
    _ ≤ 1 / (2 * Real.pi * T) * (c₀ ^ 2 / 2 * (2 * δ * (2 / (((N : ℝ) + 1) ^ 2 * C')) * SL -
          Bfun (dimerPotZ J h r η) p T δ lam ⌊T⌋₊ / (((N : ℝ) + 1) ^ 2 * C') ^ 2)) := by
        gcongr
    _ ≤ dimerMomZ (dimerPotZ J h r η) p T 0 := hdet

/-- **Theorem `thm:dimer`(c)**: `M^{∞,r}_p(T) ≥ c (JT)^{p - 1/2}` for `JT ≥ T₀`. -/
theorem dimerMInf_super (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) :
    ∃ c T₀ : ℝ, 0 < c ∧ 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ J * T →
      c * (J * T) ^ (p - 1 / 2) ≤ dimerMInf J h r p T := by
  sorry

/-- **Theorem `thm:dimer`(c)**, the annealed exponent: for every `α > 0`, `M^{∞,r}_p(T) ≥ T^{p - 1/2 - α}` for all
large `T`. -/
theorem dimerMInf_exponent (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) (α : ℝ) (hα : 0 < α) :
    ∀ᶠ T in atTop, T ^ (p - 1 / 2 - α) ≤ dimerMInf J h r p T := by
  sorry

/-- **Theorem `thm:dimer`(d)** for every `p > 1/2`: `sup_{N even, x, t} D^{(N)}_p(t; x) = ∞`. -/
theorem dimerD_unbounded_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerD J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

/-- **Theorem `thm:dimer`(d)** for every `p > 1/2`, for the single-fermion dimer state. -/
theorem dimer_LGT_unbounded_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerDLGT J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

/-- **Hypothesis (EC) fails for the dimer law** (paper: `thm:dimer`(d)), for every `p > 1/2`: there is no `F ≥ 0`
with `∑_d (d+1)^p F(d) < ∞` and `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N`. -/
theorem not_EC_dimer_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d)
    (hFs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d)) :
    ¬ ∀ (N : ℕ), Even N → ∀ x y : Fin N,
      dimerAvg (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        F (Nat.dist x.val y.val) := by
  sorry

end DFL
