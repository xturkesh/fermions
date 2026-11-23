import DFL.DimerInfVol

/-!
# Second moment of the dimer transfer matrices near the critical energy (paper: Lemma `lem:dimersecond`)

Let `0 < λ < 1` and `D_±(z) = T_{±λ - z}²` (`dimerDp`, `dimerDm`). In the basis `S` of `lem:dimercritical`(iii)
(`dimerS`), `M_±(z) = S⁻¹ D_±(z) S` are, at `z = λ`, the rotations `-1` and `R_{2θ}`. The averaged map
`𝓛_z(X) = ½ ∑_± M_±(z)* X M_±(z)` on Hermitian `2 × 2` matrices, in the coordinates
`X = x 1 + y σ_y + w₁ σ_z + w₂ σ_x`, is the identity on `(x, y)` and the contraction `½(1 + R_{4θ})` (norm
`|cos 2θ| = |2λ² - 1| < 1`) on `w` at `z = λ`; for real `z` the `y`-coordinate decouples (the `M_±` are real
symplectic), and the first-order change of the `(x, x)`-entry vanishes (`det M_± ≡ 1`). In the weighted norm
`|x| + |y| + c|w|`, `c ≍ |z - λ|`, this gives `‖𝓛_z‖ ≤ 1 + K((Re z - λ)² + |Im z|)`, hence
`E ‖D_{ε_n}(z) ⋯ D_{ε_1}(z)‖₂² ≤ C e^{K n ((Re z - λ)² + |Im z|)}` (`dimer_second_moment`).

`integral_hs2_tmProd_le` is the consequence for the transfer matrices `T^z(m + L, m)` (`tmProd`) of the random dimer
potential `v^{(r)}` over any window: complete dimers are i.i.d., and at most two further sites contribute a bounded
factor.
-/

open Matrix Finset MeasureTheory

noncomputable section

namespace DFL

/-- The squared Hilbert–Schmidt norm of a `2 × 2` complex matrix. -/
def hs2 (A : Matrix (Fin 2) (Fin 2) ℂ) : ℝ := ∑ i, ∑ j, ‖A i j‖ ^ 2

/-! ### Proof: a Lyapunov quadratic form

Instead of the weighted norm of the paper we use the explicit Hermitian form `G_u(x) = x* Y_u x`
with `Y_u = !![1, λ; λ, 1 - 2λu]`, `u = Re z - λ` (`Gq`). At `z = λ` it is invariant under
`D_±(λ)` (it is `(S S*)⁻¹` up to a factor), and the first-order correction `-2λu E₂₂` is chosen
so that `G_u(D₊(z)w) + G_u(D₋(z)w) = 2 G_u(w) + 2 w* E w` with every entry of `E` of size
`O((Re z - λ)² + |Im z|)` (`key_identity`, `Eform_le`). Since `G_u ≍ |·|²`, this gives
`½ ∑_± G_u(D_± w) ≤ (1 + K ε) G_u(w)` (`one_step`); iterating over the words (`sum_iter`) and
summing over the columns of the products gives the bound.
-/

namespace DimerSecondAux

/-- The quadratic form `x* Y x`, `Y = !![1, λ; λ, 1 - 2λu]`, in the real coordinates of
`x = (a + ib, c + id)`. -/
def Gf (lam u a b c d : ℝ) : ℝ :=
  a ^ 2 + b ^ 2 + c ^ 2 + d ^ 2 + 2 * lam * (a * c + b * d) - 2 * lam * u * (c ^ 2 + d ^ 2)

/-- The quadratic form `G_u(x) = x* Y x` on `ℂ²`. -/
def Gq (lam u : ℝ) (x : Fin 2 → ℂ) : ℝ := Gf lam u (x 0).re (x 0).im (x 1).re (x 1).im

/-- Real part of the first entry of `!![b² - 1, b; -b, -1] (p + iq, r + is)`, `b = U + iv`. -/
def c0 (U v p q r s : ℝ) : ℝ := (U ^ 2 - v ^ 2 - 1) * p - 2 * U * v * q + U * r - v * s

/-- Imaginary part of the first entry. -/
def c1 (U v p q r s : ℝ) : ℝ := (U ^ 2 - v ^ 2 - 1) * q + 2 * U * v * p + U * s + v * r

/-- Real part of the second entry. -/
def c2 (U v p q r _s : ℝ) : ℝ := -(U * p - v * q) - r

/-- Imaginary part of the second entry. -/
def c3 (U v p q _r s : ℝ) : ℝ := -(U * q + v * p) - s

/-- The error form `w* (𝓛(Y) - Y) w` in real coordinates (`u = Re z - λ`, `v = Im z`). -/
def Eform (lam u v p q r s : ℝ) : ℝ :=
  ((u ^ 2 + v ^ 2) ^ 2 + (2 * lam ^ 2 - 1) * u ^ 2 + (2 * lam ^ 2 + 3) * v ^ 2) * (p ^ 2 + q ^ 2) +
    (u ^ 2 + v ^ 2) * (r ^ 2 + s ^ 2) +
    2 * ((lam * (v ^ 2 - u ^ 2) + u * (u ^ 2 + v ^ 2)) * (p * r + q * s) -
      (2 * lam * u * v - v * (u ^ 2 + v ^ 2) - 2 * v) * (p * s - q * r))

/-- The exact one-step identity `G(D₊w) + G(D₋w) = 2 G(w) + 2 w* E w`. -/
lemma key_identity (lam u v p q r s : ℝ) :
    Gf lam u (c0 u v p q r s) (c1 u v p q r s) (c2 u v p q r s) (c3 u v p q r s) +
      Gf lam u (c0 (u + 2 * lam) v p q r s) (c1 (u + 2 * lam) v p q r s)
        (c2 (u + 2 * lam) v p q r s) (c3 (u + 2 * lam) v p q r s) =
    2 * Gf lam u p q r s + 2 * Eform lam u v p q r s := by
  unfold Gf c0 c1 c2 c3 Eform
  ring

/-- Elementary bounds for `u² + v² ≤ 1`. -/
lemma uv_bounds (u v : ℝ) (huv : u ^ 2 + v ^ 2 ≤ 1) :
    |u| ≤ 1 ∧ |v| ≤ 1 ∧ v ^ 2 ≤ |v| ∧ (u ^ 2 + v ^ 2) ^ 2 ≤ u ^ 2 + v ^ 2 := by
  have hu2 : u ^ 2 ≤ 1 := by nlinarith [sq_nonneg v]
  have hv2 : v ^ 2 ≤ 1 := by nlinarith [sq_nonneg u]
  have hv1 : |v| ≤ 1 := (sq_le_one_iff_abs_le_one v).mp hv2
  refine ⟨(sq_le_one_iff_abs_le_one u).mp hu2, hv1, ?_, ?_⟩
  · rw [← sq_abs]; nlinarith [abs_nonneg v]
  · have h0 : 0 ≤ u ^ 2 + v ^ 2 := by positivity
    nlinarith

/-- `2|pr + qs| ≤ p² + q² + r² + s²` and `2|ps - qr| ≤ p² + q² + r² + s²`. -/
lemma Z_bounds (p q r s : ℝ) :
    |2 * (p * r + q * s)| ≤ p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2) ∧
      |2 * (p * s - q * r)| ≤ p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2) :=
  ⟨abs_le.mpr ⟨by nlinarith [sq_nonneg (p + r), sq_nonneg (q + s)],
      by nlinarith [sq_nonneg (p - r), sq_nonneg (q - s)]⟩,
    abs_le.mpr ⟨by nlinarith [sq_nonneg (p + s), sq_nonneg (q - r)],
      by nlinarith [sq_nonneg (p - s), sq_nonneg (q + r)]⟩⟩

/-- The `(x, x)`-type coefficient of the error form is `O(u² + |v|)`. -/
lemma E00_abs_le (lam u v : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (huv : u ^ 2 + v ^ 2 ≤ 1) :
    |(u ^ 2 + v ^ 2) ^ 2 + (2 * lam ^ 2 - 1) * u ^ 2 + (2 * lam ^ 2 + 3) * v ^ 2| ≤
      6 * (u ^ 2 + |v|) := by
  obtain ⟨-, -, hvv, hsq⟩ := uv_bounds u v huv
  have hl : 0 ≤ 1 - lam ^ 2 := by nlinarith
  have k1 := mul_nonneg hl (sq_nonneg u)
  have k2 := mul_nonneg hl (sq_nonneg v)
  have k3 := mul_nonneg (sq_nonneg lam) (sq_nonneg u)
  have k4 := mul_nonneg (sq_nonneg lam) (sq_nonneg v)
  have k5 := sq_nonneg (u ^ 2 + v ^ 2)
  have k6 := sq_nonneg u
  have k7 := sq_nonneg v
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- The real part of the off-diagonal coefficient is `O(u² + |v|)`. -/
lemma E01re_abs_le (lam u v : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (huv : u ^ 2 + v ^ 2 ≤ 1) :
    |lam * (v ^ 2 - u ^ 2) + u * (u ^ 2 + v ^ 2)| ≤ 2 * (u ^ 2 + |v|) := by
  obtain ⟨hu1, -, hvv, -⟩ := uv_bounds u v huv
  have hu := abs_le.mp hu1
  have hS : 0 ≤ u ^ 2 + v ^ 2 := by positivity
  have k1 := mul_nonneg h0.le (sq_nonneg u)
  have k2 := mul_nonneg h0.le (sq_nonneg v)
  have k3 := mul_nonneg (show 0 ≤ 1 - lam by linarith) (sq_nonneg u)
  have k4 := mul_nonneg (show 0 ≤ 1 - lam by linarith) (sq_nonneg v)
  have k5 := mul_nonneg (show 0 ≤ 1 - u by linarith) hS
  have k6 := mul_nonneg (show 0 ≤ 1 + u by linarith) hS
  exact abs_le.mpr ⟨by linarith, by linarith⟩

/-- The imaginary part of the off-diagonal coefficient is `O(|v|)`. -/
lemma E01im_abs_le (lam u v : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (huv : u ^ 2 + v ^ 2 ≤ 1) :
    |2 * lam * u * v - v * (u ^ 2 + v ^ 2) - 2 * v| ≤ 5 * (u ^ 2 + |v|) := by
  obtain ⟨hu1, -, -, -⟩ := uv_bounds u v huv
  have hu := abs_le.mp hu1
  have e : 2 * lam * u * v - v * (u ^ 2 + v ^ 2) - 2 * v =
      v * (2 * lam * u - (u ^ 2 + v ^ 2) - 2) := by ring
  have k1 := mul_nonneg h0.le (show 0 ≤ 1 - u by linarith)
  have k2 := mul_nonneg h0.le (show 0 ≤ 1 + u by linarith)
  have e1 : |2 * lam * u - (u ^ 2 + v ^ 2) - 2| ≤ 5 :=
    abs_le.mpr ⟨by nlinarith [sq_nonneg u, sq_nonneg v], by nlinarith [sq_nonneg u, sq_nonneg v]⟩
  rw [e, abs_mul]
  have := mul_le_mul_of_nonneg_left e1 (abs_nonneg v)
  nlinarith [sq_nonneg u]

/-- The error form is `O((Re z - λ)² + |Im z|)`. -/
lemma Eform_le (lam u v p q r s : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (huv : u ^ 2 + v ^ 2 ≤ 1) :
    Eform lam u v p q r s ≤ 13 * (u ^ 2 + |v|) * (p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2)) := by
  obtain ⟨hZ1, hZ2⟩ := Z_bounds p q r s
  have hE00 := E00_abs_le lam u v h0 h1 huv
  have ha := E01re_abs_le lam u v h0 h1 huv
  have hb := E01im_abs_le lam u v h0 h1 huv
  obtain ⟨-, -, hvv, -⟩ := uv_bounds u v huv
  have hX0 : 0 ≤ p ^ 2 + q ^ 2 := by positivity
  have hY0 : 0 ≤ r ^ 2 + s ^ 2 := by positivity
  have hε : 0 ≤ u ^ 2 + |v| := by positivity
  have k1 := mul_le_mul_of_nonneg_right
    ((le_abs_self _).trans hE00) hX0
  have k2 : (u ^ 2 + v ^ 2) * (r ^ 2 + s ^ 2) ≤ (u ^ 2 + |v|) * (r ^ 2 + s ^ 2) :=
    mul_le_mul_of_nonneg_right (by linarith) hY0
  have k3 : (lam * (v ^ 2 - u ^ 2) + u * (u ^ 2 + v ^ 2)) * (2 * (p * r + q * s)) ≤
      2 * (u ^ 2 + |v|) * (p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2)) := by
    refine (le_abs_self _).trans ?_
    rw [abs_mul]
    exact mul_le_mul ha hZ1 (abs_nonneg _) (by positivity)
  have k4 : -((2 * lam * u * v - v * (u ^ 2 + v ^ 2) - 2 * v) * (2 * (p * s - q * r))) ≤
      5 * (u ^ 2 + |v|) * (p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2)) := by
    refine (neg_le_abs _).trans ?_
    rw [abs_mul]
    exact mul_le_mul hb hZ2 (abs_nonneg _) (by positivity)
  have k5 : 0 ≤ (u ^ 2 + |v|) * (r ^ 2 + s ^ 2) := by positivity
  have e : Eform lam u v p q r s =
      ((u ^ 2 + v ^ 2) ^ 2 + (2 * lam ^ 2 - 1) * u ^ 2 + (2 * lam ^ 2 + 3) * v ^ 2) *
        (p ^ 2 + q ^ 2) + (u ^ 2 + v ^ 2) * (r ^ 2 + s ^ 2) +
        (lam * (v ^ 2 - u ^ 2) + u * (u ^ 2 + v ^ 2)) * (2 * (p * r + q * s)) +
        -((2 * lam * u * v - v * (u ^ 2 + v ^ 2) - 2 * v) * (2 * (p * s - q * r))) := by
    unfold Eform; ring
  rw [e]
  nlinarith

/-- Lower bound: `G_u(x) ≥ (1 - λ)/2 · |x|²` for `|u| ≤ (1 - λ)/4`. -/
lemma Gf_ge (lam u p q r s : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (hu : |u| ≤ (1 - lam) / 4) :
    (1 - lam) / 2 * (p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2) ≤ Gf lam u p q r s := by
  unfold Gf
  have hZ : -(p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2) ≤ 2 * (p * r + q * s) := by
    nlinarith [sq_nonneg (p + r), sq_nonneg (q + s)]
  have hu' := abs_le.mp hu
  have hY : 0 ≤ r ^ 2 + s ^ 2 := by positivity
  have h2 : 2 * lam * u * (r ^ 2 + s ^ 2) ≤ (1 - lam) / 2 * (r ^ 2 + s ^ 2) := by
    have : 2 * lam * u ≤ (1 - lam) / 2 := by nlinarith
    exact mul_le_mul_of_nonneg_right this hY
  nlinarith [mul_le_mul_of_nonneg_left hZ h0.le]

/-- `T_{-b}² = !![b² - 1, b; -b, -1]`. -/
lemma tmC_neg_sq (b : ℂ) : tmC (-b) * tmC (-b) = !![b ^ 2 - 1, b; -b, -1] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [tmC, Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- `G_u(T_{-b}² w)` in real coordinates. -/
lemma Gq_sq (lam u : ℝ) (b : ℂ) (w : Fin 2 → ℂ) :
    Gq lam u ((tmC (-b) * tmC (-b)) *ᵥ w) =
      Gf lam u (c0 b.re b.im (w 0).re (w 0).im (w 1).re (w 1).im)
        (c1 b.re b.im (w 0).re (w 0).im (w 1).re (w 1).im)
        (c2 b.re b.im (w 0).re (w 0).im (w 1).re (w 1).im)
        (c3 b.re b.im (w 0).re (w 0).im (w 1).re (w 1).im) := by
  rw [tmC_neg_sq]
  simp only [Gq, Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.of_apply, Matrix.cons_val',
    Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.empty_val', Matrix.cons_val_fin_one,
    Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im, Complex.sub_re,
    Complex.sub_im, Complex.neg_re, Complex.neg_im, Complex.one_re, Complex.one_im, sq]
  unfold c0 c1 c2 c3
  ring_nf

/-- **One step**: `G(D₊(z) w) + G(D₋(z) w) ≤ 2 (1 + K ε) G(w)`, `ε = (Re z - λ)² + |Im z|`,
`K = 26/(1 - λ)`, for `|z - λ|` small. -/
lemma one_step (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (z : ℂ)
    (huv : (z.re - lam) ^ 2 + z.im ^ 2 ≤ 1) (hu : |z.re - lam| ≤ (1 - lam) / 4)
    (w : Fin 2 → ℂ) :
    Gq lam (z.re - lam) (dimerDp lam z *ᵥ w) + Gq lam (z.re - lam) (dimerDm lam z *ᵥ w) ≤
      2 * (1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|)) * Gq lam (z.re - lam) w := by
  have hp : dimerDp lam z = tmC (-(z - lam)) * tmC (-(z - lam)) := by rw [dimerDp, neg_sub]
  have hm : dimerDm lam z = tmC (-(z + lam)) * tmC (-(z + lam)) := by
    rw [dimerDm, neg_add', sub_eq_add_neg, add_comm (-(lam : ℂ)) (-z), ← sub_eq_add_neg]
  rw [hp, hm, Gq_sq, Gq_sq]
  have e1 : (z - (lam : ℂ)).re = z.re - lam := by simp
  have e2 : (z - (lam : ℂ)).im = z.im := by simp
  have e3 : (z + (lam : ℂ)).re = (z.re - lam) + 2 * lam := by simp; ring
  have e4 : (z + (lam : ℂ)).im = z.im := by simp
  rw [e1, e2, e3, e4, key_identity]
  set u := z.re - lam with hudef
  set v := z.im with hvdef
  set p := (w 0).re
  set q := (w 0).im
  set r := (w 1).re
  set s := (w 1).im
  have hE := Eform_le lam u v p q r s h0 h1 huv
  have hG := Gf_ge lam u p q r s h0 h1 hu
  have hGq : Gq lam u w = Gf lam u p q r s := rfl
  rw [hGq]
  have hl : 0 < 1 - lam := by linarith
  have hε : 0 ≤ u ^ 2 + |v| := by positivity
  have hN : p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2) ≤ 2 / (1 - lam) * Gf lam u p q r s := by
    have : p ^ 2 + q ^ 2 + (r ^ 2 + s ^ 2) =
        2 / (1 - lam) * ((1 - lam) / 2 * (p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2)) := by
      field_simp; ring
    rw [this]
    exact mul_le_mul_of_nonneg_left hG (by positivity)
  have hN' := mul_le_mul_of_nonneg_left hN (show 0 ≤ 13 * (u ^ 2 + |v|) by positivity)
  have e5 : 13 * (u ^ 2 + |v|) * (2 / (1 - lam) * Gf lam u p q r s) =
      26 / (1 - lam) * (u ^ 2 + |v|) * Gf lam u p q r s := by ring
  nlinarith

/-- **Iteration**: if `∑_a G(D_a w) ≤ c G(w)` for all `w`, then
`∑_{s ∈ αⁿ} G(D_{s₀} ⋯ D_{s_{n-1}} w) ≤ cⁿ G(w)`. -/
lemma sum_iter {α : Type*} [Fintype α] (D : α → Matrix (Fin 2) (Fin 2) ℂ)
    (G : (Fin 2 → ℂ) → ℝ) (c : ℝ) (hc : 0 ≤ c) (hstep : ∀ w, ∑ a, G (D a *ᵥ w) ≤ c * G w) :
    ∀ (n : ℕ) (w : Fin 2 → ℂ),
      ∑ s : Fin n → α, G ((List.ofFn fun i => D (s i)).prod *ᵥ w) ≤ c ^ n * G w := by
  intro n
  induction n with
  | zero => intro w; simp
  | succ n ih =>
    intro w
    rw [← (Fin.consEquiv fun _ => α).sum_comp, Fintype.sum_prod_type]
    simp only [Fin.consEquiv_apply, List.ofFn_succ, Fin.cons_zero, Fin.cons_succ, List.prod_cons,
      ← Matrix.mulVec_mulVec]
    rw [Finset.sum_comm]
    calc ∑ t : Fin n → α, ∑ a, G (D a *ᵥ ((List.ofFn fun i => D (t i)).prod *ᵥ w))
        ≤ ∑ t : Fin n → α, c * G ((List.ofFn fun i => D (t i)).prod *ᵥ w) :=
          Finset.sum_le_sum fun t _ => hstep _
      _ = c * ∑ t : Fin n → α, G ((List.ofFn fun i => D (t i)).prod *ᵥ w) := by
          rw [Finset.mul_sum]
      _ ≤ c * (c ^ n * G w) := mul_le_mul_of_nonneg_left (ih w) hc
      _ = c ^ (n + 1) * G w := by ring

/-- The squared norm of a vector of `ℂ²` in real coordinates. -/
def Nrm (x : Fin 2 → ℂ) : ℝ := (x 0).re ^ 2 + (x 0).im ^ 2 + ((x 1).re ^ 2 + (x 1).im ^ 2)

/-- `‖A‖₂² = |A e₀|² + |A e₁|²`. -/
lemma hs2_eq (A : Matrix (Fin 2) (Fin 2) ℂ) :
    hs2 A = Nrm (A *ᵥ ![1, 0]) + Nrm (A *ᵥ ![0, 1]) := by
  simp only [hs2, Nrm, Fin.sum_univ_two, Matrix.mulVec, dotProduct, Complex.sq_norm,
    Complex.normSq_apply, Matrix.cons_val_zero, Matrix.cons_val_one, mul_one, mul_zero, add_zero,
    zero_add]
  ring

/-- `|x|² ≤ 2/(1-λ) · G_u(x)`. -/
lemma Nrm_le (lam u : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (hu : |u| ≤ (1 - lam) / 4)
    (x : Fin 2 → ℂ) : Nrm x ≤ 2 / (1 - lam) * Gq lam u x := by
  have hG := Gf_ge lam u (x 0).re (x 0).im (x 1).re (x 1).im h0 h1 hu
  have hl : 0 < 1 - lam := by linarith
  have : Nrm x = 2 / (1 - lam) * ((1 - lam) / 2 *
      ((x 0).re ^ 2 + (x 0).im ^ 2 + (x 1).re ^ 2 + (x 1).im ^ 2)) := by
    unfold Nrm; field_simp; ring
  rw [this]
  exact mul_le_mul_of_nonneg_left hG (by positivity)

/-- `G_u(e₀) + G_u(e₁) = 2 - 2λu ≤ 4`. -/
lemma Gq_e_le (lam u : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (hu : |u| ≤ 1) :
    Gq lam u ![1, 0] + Gq lam u ![0, 1] ≤ 4 := by
  simp [Gq, Gf]
  have := abs_le.mp hu
  nlinarith

/-- **The finite second-moment bound**: for letters `D_a` (`a ∈ α`) with
`∑_a G_u(D_a w) ≤ 2ρ G_u(w)`, `2^{-n} ∑_s ‖D_{s₀} ⋯ D_{s_{n-1}}‖₂² ≤ 8/(1-λ) · ρⁿ`. -/
lemma avg_hs2_le {α : Type*} [Fintype α] (D : α → Matrix (Fin 2) (Fin 2) ℂ) (lam u ρ : ℝ)
    (h0 : 0 < lam) (h1 : lam < 1) (hu : |u| ≤ (1 - lam) / 4) (hρ : 0 ≤ ρ)
    (hstep : ∀ w, ∑ a, Gq lam u (D a *ᵥ w) ≤ 2 * ρ * Gq lam u w) (n : ℕ) :
    (∑ s : Fin n → α, hs2 (List.ofFn fun i => D (s i)).prod) / 2 ^ n ≤
      8 / (1 - lam) * ρ ^ n := by
  have hl : 0 < 1 - lam := by linarith
  have hu1 : |u| ≤ 1 := hu.trans (by linarith)
  have h2 := sum_iter D (Gq lam u) (2 * ρ) (by positivity) hstep n
  have hsum : ∑ s : Fin n → α, hs2 (List.ofFn fun i => D (s i)).prod ≤
      2 / (1 - lam) * (∑ s : Fin n → α, Gq lam u ((List.ofFn fun i => D (s i)).prod *ᵥ ![1, 0]) +
        ∑ s : Fin n → α, Gq lam u ((List.ofFn fun i => D (s i)).prod *ᵥ ![0, 1])) := by
    rw [← Finset.sum_add_distrib, Finset.mul_sum]
    refine Finset.sum_le_sum fun s _ => ?_
    rw [hs2_eq, mul_add]
    exact add_le_add (Nrm_le lam u h0 h1 hu _) (Nrm_le lam u h0 h1 hu _)
  have h3 := Gq_e_le lam u h0 h1 hu1
  have h4 := add_le_add (h2 ![1, 0]) (h2 ![0, 1])
  have h5 : ∑ s : Fin n → α, hs2 (List.ofFn fun i => D (s i)).prod ≤
      2 / (1 - lam) * ((2 * ρ) ^ n * 4) := by
    refine hsum.trans (mul_le_mul_of_nonneg_left ?_ (by positivity))
    calc _ ≤ (2 * ρ) ^ n * Gq lam u ![1, 0] + (2 * ρ) ^ n * Gq lam u ![0, 1] := h4
      _ = (2 * ρ) ^ n * (Gq lam u ![1, 0] + Gq lam u ![0, 1]) := by ring
      _ ≤ (2 * ρ) ^ n * 4 := mul_le_mul_of_nonneg_left h3 (by positivity)
  rw [div_le_iff₀ (by positivity)]
  calc _ ≤ 2 / (1 - lam) * ((2 * ρ) ^ n * 4) := h5
    _ = 8 / (1 - lam) * ρ ^ n * 2 ^ n := by rw [mul_pow]; ring

/-- `|z - λ| ≤ δ ≤ 1` gives `(Re z - λ)² + (Im z)² ≤ 1` and `|Re z - λ| ≤ δ`. -/
lemma near (lam δ : ℝ) (hδ : δ ≤ 1) (z : ℂ) (hz : ‖z - (lam : ℂ)‖ ≤ δ) :
    (z.re - lam) ^ 2 + z.im ^ 2 ≤ 1 ∧ |z.re - lam| ≤ δ := by
  have e : ‖z - (lam : ℂ)‖ ^ 2 = (z.re - lam) ^ 2 + z.im ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; simp; ring
  have hn := norm_nonneg (z - (lam : ℂ))
  refine ⟨?_, ?_⟩
  · rw [← e]; nlinarith
  · have := Complex.abs_re_le_norm (z - (lam : ℂ))
    simp only [Complex.sub_re, Complex.ofReal_re] at this
    linarith

/-- A sum over `ZMod 2`. -/
lemma sum_zmod2 (F : ZMod 2 → ℝ) : ∑ a : ZMod 2, F a = F 0 + F 1 := by
  rw [show (Finset.univ : Finset (ZMod 2)) = {0, 1} by decide, Finset.sum_pair (by decide)]

/-- The dimer transfer matrix `T_{c sgn(a) - z}²` for the dimer sign `a`. -/
def Dz (c : ℝ) (z : ℂ) (a : ZMod 2) : Matrix (Fin 2) (Fin 2) ℂ :=
  tmC (((c * sgn a : ℝ) : ℂ) - z) * tmC (((c * sgn a : ℝ) : ℂ) - z)

/-- For `c = ±λ`, the two letters `Dz c z a` are `D₊(z)` and `D₋(z)`. -/
lemma sum_Dz (lam c : ℝ) (hc : c = lam ∨ c = -lam) (z : ℂ) (F : Matrix (Fin 2) (Fin 2) ℂ → ℝ) :
    ∑ a : ZMod 2, F (Dz c z a) = F (dimerDp lam z) + F (dimerDm lam z) := by
  rw [sum_zmod2]
  have s0 : sgn 0 = 1 := by simp [sgn]
  have s1 : sgn 1 = -1 := by simp [sgn]
  rcases hc with rfl | rfl
  · simp only [Dz, s0, s1, dimerDp, dimerDm, mul_one, mul_neg, Complex.ofReal_neg]
  · simp only [Dz, s0, s1, dimerDp, dimerDm, mul_one, mul_neg, neg_neg, Complex.ofReal_neg]
    ring

/-- The window map `η ↦ (η(k + n - 1 - i))_{i < n}` (the dimer signs from right to left). -/
def winR (k : ℤ) (n : ℕ) (η : CfgZ) : Fin n → ZMod 2 := fun i => η (k + n - 1 - ((i : ℕ) : ℤ))

/-- The extension of a window configuration to `ℤ` (by `0` outside `[k, k + n)`). -/
def extR (k : ℤ) (n : ℕ) (t : Fin n → ZMod 2) (j : ℤ) : ZMod 2 :=
  if h : k ≤ j ∧ j < k + n then t ⟨(k + n - 1 - j).toNat, by omega⟩ else 0

/-- `winR` is measurable. -/
lemma measurable_winR (k : ℤ) (n : ℕ) : Measurable (winR k n) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- `extR` inverts `winR` on the window. -/
lemma extR_winR (k : ℤ) (n : ℕ) (η : CfgZ) (j : ℤ) (hj : k ≤ j ∧ j < k + n) :
    extR k n (winR k n η) j = η j := by
  simp only [extR, dite_eq_left hj, winR]
  congr 1
  omega

/-- The preimage of a window configuration is a cylinder set. -/
lemma winR_preimage (k : ℤ) (n : ℕ) (t : Fin n → ZMod 2) :
    winR k n ⁻¹' {t} = Set.pi ↑(Finset.Ico k (k + n)) (fun j => {extR k n t j}) := by
  ext η
  simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_pi, Finset.coe_Ico, Set.mem_Ico]
  constructor
  · rintro rfl j hj
    exact (extR_winR k n η j hj).symm
  · intro H
    funext i
    have hi := i.isLt
    have h1 := H (k + n - 1 - ((i : ℕ) : ℤ)) ⟨by omega, by omega⟩
    simp only [winR]
    rw [h1, extR, dite_eq_left ⟨by omega, by omega⟩]
    congr 1
    ext
    simp only
    omega

/-- Under `muZ`, each window configuration has probability `2^{-n}`. -/
lemma muZ_winR_preimage (k : ℤ) (n : ℕ) (t : Fin n → ZMod 2) :
    muZ (winR k n ⁻¹' {t}) = ((2 : ENNReal) ^ n)⁻¹ := by
  rw [winR_preimage, muZ, Measure.infinitePi_pi _ (fun _ _ => MeasurableSet.of_discrete)]
  simp only [InfVolAux.unifZ2_singleton]
  rw [Finset.prod_const, Int.card_Ico, ENNReal.inv_pow]
  congr 2
  omega

/-- Integrals of functions of a window are uniform averages. -/
lemma integral_comp_winR (k : ℤ) (n : ℕ) (f : (Fin n → ZMod 2) → ℝ) :
    ∫ η, f (winR k n η) ∂muZ = (∑ t, f t) / 2 ^ n := by
  rw [← integral_map (measurable_winR k n).aemeasurable
    (measurable_of_countable f).aestronglyMeasurable, integral_fintype Integrable.of_finite,
    Finset.sum_div]
  refine Finset.sum_congr rfl fun t _ => ?_
  rw [Measure.real, Measure.map_apply (measurable_winR k n) (measurableSet_singleton t),
    muZ_winR_preimage, smul_eq_mul, ENNReal.toReal_inv, ENNReal.toReal_pow]
  norm_num
  ring

/-- Functions of a window are integrable. -/
lemma integrable_comp_winR (k : ℤ) (n : ℕ) (f : (Fin n → ZMod 2) → ℝ) :
    Integrable (fun η => f (winR k n η)) muZ :=
  (Integrable.of_finite (μ := Measure.map (winR k n) muZ)).comp_measurable (measurable_winR k n)

/-- `T^z(m + L, m)` only depends on the potential on the sites `m, …, m + L - 1`. -/
lemma tmProd_congr (v v' : ℤ → ℝ) (z : ℂ) (m : ℤ) :
    ∀ L : ℕ, (∀ i : ℕ, i < L → v (m + i) = v' (m + i)) → tmProd v z m L = tmProd v' z m L
  | 0, _ => rfl
  | L + 1, h => by
    rw [tmProd, tmProd, h L (by omega), tmProd_congr v v' z m L (fun i hi => h i (by omega))]

/-- `‖T^z(m + L, m)‖₂²` is a function of the dimer signs in a window. -/
lemma hs2_tmProd_eq_winR (J h : ℝ) (r m : ℤ) (L : ℕ) (z : ℂ) (η : CfgZ) :
    hs2 (tmProd (dimerPotZ J h r η) z m L) =
      hs2 (tmProd (dimerPotZ J h r (extR ((m + r) / 2) (L + 1)
        (winR ((m + r) / 2) (L + 1) η))) z m L) := by
  congr 1
  refine tmProd_congr _ _ z m L fun i hi => ?_
  simp only [dimerPotZ]
  rw [extR_winR _ _ _ _ ⟨by omega, by omega⟩]

/-- `‖T^z(m + L, m)‖₂²` is integrable. -/
lemma integrable_hs2_tmProd (J h : ℝ) (r m : ℤ) (L : ℕ) (z : ℂ) :
    Integrable (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) muZ := by
  have e : (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) = fun η =>
      (fun t => hs2 (tmProd (dimerPotZ J h r (extR ((m + r) / 2) (L + 1) t)) z m L))
        (winR ((m + r) / 2) (L + 1) η) :=
    funext fun η => hs2_tmProd_eq_winR J h r m L z η
  rw [e]
  exact integrable_comp_winR ((m + r) / 2) (L + 1)
    (fun t => hs2 (tmProd (dimerPotZ J h r (extR ((m + r) / 2) (L + 1) t)) z m L))

/-- Over `n` complete dimers starting at `2k - r`: `T^z = D_{η_{k+n-1}} ⋯ D_{η_k}`. -/
lemma tmProd_dimers (J h : ℝ) (r k : ℤ) (z : ℂ) (η : CfgZ) : ∀ n : ℕ,
    tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n) =
      (List.ofFn fun i => Dz (2 * h / J) z (winR k n η i)).prod
  | 0 => by simp [tmProd]
  | n + 1 => by
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, tmProd, tmProd, tmProd_dimers J h r k z η n,
      List.ofFn_succ, List.prod_cons, ← Matrix.mul_assoc]
    congr 1
    · have e1 : (2 * k - r + ((2 * n + 1 : ℕ) : ℤ) + r) / 2 =
          k + (n + 1 : ℕ) - 1 - ((0 : ℕ) : ℤ) := by
        push_cast; omega
      have e2 : (2 * k - r + ((2 * n : ℕ) : ℤ) + r) / 2 = k + (n + 1 : ℕ) - 1 - ((0 : ℕ) : ℤ) := by
        push_cast; omega
      simp only [dimerPotZ, Dz, winR, Fin.val_zero]
      rw [e1, e2]
    · congr 2
      funext i
      simp only [winR, Fin.val_succ]
      congr 2
      push_cast
      ring

/-- `(1 + Kε)ⁿ ≤ e^{K L ε}` for `K, ε ≥ 0` and `n ≤ L`. -/
lemma rho_pow_le (K ε : ℝ) (hK : 0 ≤ K) (hε : 0 ≤ ε) (n L : ℕ) (hn : n ≤ L) :
    (1 + K * ε) ^ n ≤ Real.exp (K * L * ε) := by
  have h1 : 1 + K * ε ≤ Real.exp (K * ε) := by linarith [Real.add_one_le_exp (K * ε)]
  have h2 : 0 ≤ K * ε := mul_nonneg hK hε
  calc (1 + K * ε) ^ n ≤ Real.exp (K * ε) ^ n := pow_le_pow_left₀ (by positivity) h1 n
    _ = Real.exp (n * (K * ε)) := (Real.exp_nat_mul _ n).symm
    _ ≤ Real.exp (K * L * ε) := by
      rw [Real.exp_le_exp]
      have : (n : ℝ) ≤ L := by exact_mod_cast hn
      nlinarith

/-- Cauchy–Schwarz in `ℂ²`. -/
lemma norm_sq_add_mul_le (a b c d : ℂ) :
    ‖a * c + b * d‖ ^ 2 ≤ (‖a‖ ^ 2 + ‖b‖ ^ 2) * (‖c‖ ^ 2 + ‖d‖ ^ 2) := by
  have h1 : ‖a * c + b * d‖ ≤ ‖a‖ * ‖c‖ + ‖b‖ * ‖d‖ := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, norm_mul]
  have h2 := pow_le_pow_left₀ (norm_nonneg _) h1 2
  nlinarith [sq_nonneg (‖a‖ * ‖d‖ - ‖b‖ * ‖c‖)]

/-- The Hilbert–Schmidt norm is submultiplicative. -/
lemma hs2_mul_le (A B : Matrix (Fin 2) (Fin 2) ℂ) : hs2 (A * B) ≤ hs2 A * hs2 B := by
  simp only [hs2, Matrix.mul_apply, Fin.sum_univ_two]
  have k00 := norm_sq_add_mul_le (A 0 0) (A 0 1) (B 0 0) (B 1 0)
  have k01 := norm_sq_add_mul_le (A 0 0) (A 0 1) (B 0 1) (B 1 1)
  have k10 := norm_sq_add_mul_le (A 1 0) (A 1 1) (B 0 0) (B 1 0)
  have k11 := norm_sq_add_mul_le (A 1 0) (A 1 1) (B 0 1) (B 1 1)
  linarith

/-- `‖A‖₂² ≥ 0`. -/
lemma hs2_nonneg (A : Matrix (Fin 2) (Fin 2) ℂ) : 0 ≤ hs2 A := by
  unfold hs2; positivity

/-- `‖T_a‖₂² = |a|² + 2`. -/
lemma hs2_tmC (a : ℂ) : hs2 (tmC a) = ‖a‖ ^ 2 + 2 := by
  simp [hs2, tmC, Fin.sum_univ_two]
  ring

/-- At most one site: `‖T^z(x + L, x)‖₂² ≤ 11` for `L ≤ 1`, `v = ±λ` and `|z - λ| ≤ 1`. -/
lemma hs2_tmProd_short (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (v : ℤ → ℝ)
    (hv : ∀ y, v y = lam ∨ v y = -lam) (z : ℂ) (hz : ‖z - lam‖ ≤ 1) (x : ℤ) (L : ℕ)
    (hL : L ≤ 1) : hs2 (tmProd v z x L) ≤ 11 := by
  interval_cases L
  · simp [tmProd, hs2, Fin.sum_univ_two, Matrix.one_apply]
    norm_num
  · simp only [tmProd, Matrix.mul_one]
    rw [hs2_tmC]
    have hw := DimerTransferAux.norm_w_sub_le lam h0 (hv (x + ((0 : ℕ) : ℤ)))
    have htri : ‖(v (x + ((0 : ℕ) : ℤ)) : ℂ) - z‖ ≤
        ‖(v (x + ((0 : ℕ) : ℤ)) : ℂ) - lam‖ + ‖(lam : ℂ) - z‖ :=
      norm_sub_le_norm_sub_add_norm_sub _ _ _
    rw [norm_sub_rev (lam : ℂ) z] at htri
    have h3 : ‖(v (x + ((0 : ℕ) : ℤ)) : ℂ) - z‖ ≤ 3 := by linarith
    have h4 := norm_nonneg ((v (x + ((0 : ℕ) : ℤ)) : ℂ) - z)
    nlinarith

/-- Splitting a window into at most one site, complete dimers, and at most one site. -/
lemma window_split (r m : ℤ) (L : ℕ) : ∃ (k : ℤ) (n L₁ L₂ : ℕ), L₁ ≤ 1 ∧ L₂ ≤ 1 ∧ n ≤ L ∧
    ∀ (v : ℤ → ℝ) (z : ℂ), tmProd v z m L =
      tmProd v z (2 * k - r + 2 * n) L₂ * tmProd v z (2 * k - r) (2 * n) * tmProd v z m L₁ := by
  rcases Nat.eq_zero_or_pos L with rfl | hL
  · exact ⟨0, 0, 0, 0, by omega, by omega, le_rfl, fun v z => by simp [tmProd]⟩
  obtain ⟨k, hk⟩ := Int.even_or_odd' (m + r)
  rcases hk with hk | hk
  · obtain ⟨n, L₂, hnL, hL₂⟩ : ∃ n L₂ : ℕ, L = 2 * n + L₂ ∧ L₂ ≤ 1 :=
      ⟨L / 2, L % 2, by omega, by omega⟩
    subst hnL
    refine ⟨k, n, 0, L₂, by omega, hL₂, by omega, fun v z => ?_⟩
    have hm : 2 * k - r = m := by omega
    rw [hm, DimerTransferAux.tmProd_add, show tmProd v z m 0 = 1 from rfl, Matrix.mul_one]
    push_cast
    rfl
  · obtain ⟨n, L₂, hnL, hL₂⟩ : ∃ n L₂ : ℕ, L = 1 + (2 * n + L₂) ∧ L₂ ≤ 1 :=
      ⟨(L - 1) / 2, (L - 1) % 2, by omega, by omega⟩
    subst hnL
    refine ⟨k + 1, n, 1, L₂, le_rfl, hL₂, by omega, fun v z => ?_⟩
    rw [DimerTransferAux.tmProd_add, DimerTransferAux.tmProd_add]
    have e1 : m + ((1 : ℕ) : ℤ) + ((2 * n : ℕ) : ℤ) = 2 * (k + 1) - r + 2 * (n : ℤ) := by
      push_cast; omega
    have e2 : m + ((1 : ℕ) : ℤ) = 2 * (k + 1) - r := by push_cast; omega
    rw [e1, e2]

end DimerSecondAux

open DimerSecondAux

/-- **Second moment near the critical energy** (paper: Lemma `lem:dimersecond`): the average over the `2^n` dimer
words of `‖D_{ε_1}(z) ⋯ D_{ε_n}(z)‖₂²` is at most `C e^{K n ((Re z - λ)² + |Im z|)}` for `|z - λ| ≤ δ₀`. -/
theorem dimer_second_moment (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ∃ C K δ₀ : ℝ, 0 < C ∧ 0 ≤ K ∧ 0 < δ₀ ∧ ∀ z : ℂ, ‖z - (lam : ℂ)‖ ≤ δ₀ → ∀ n : ℕ,
      (∑ s : Fin n → Bool,
          hs2 (List.ofFn fun i => if s i then dimerDp lam z else dimerDm lam z).prod) / 2 ^ n ≤
        C * Real.exp (K * n * ((z.re - lam) ^ 2 + |z.im|)) := by
  have hl : 0 < 1 - lam := by linarith
  refine ⟨8 / (1 - lam), 26 / (1 - lam), (1 - lam) / 4, by positivity, by positivity,
    by positivity, fun z hz n => ?_⟩
  obtain ⟨huv, hu⟩ := near lam ((1 - lam) / 4) (by linarith) z hz
  have hε : 0 ≤ (z.re - lam) ^ 2 + |z.im| := by positivity
  have hρ : 0 ≤ 1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|) := by positivity
  have hstep : ∀ w, ∑ b : Bool, Gq lam (z.re - lam)
      ((if b then dimerDp lam z else dimerDm lam z) *ᵥ w) ≤
        2 * (1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|)) * Gq lam (z.re - lam) w := by
    intro w
    simpa only [Fintype.sum_bool, ite_true, Bool.false_eq_true, ite_false] using
      one_step lam h0 h1 z huv hu w
  refine (avg_hs2_le (fun b : Bool => if b then dimerDp lam z else dimerDm lam z) lam
    (z.re - lam) _ h0 h1 hu hρ hstep n).trans ?_
  exact mul_le_mul_of_nonneg_left (rho_pow_le _ _ (by positivity) hε n n le_rfl) (by positivity)

/-- **Second moment over windows of the random dimer potential**: for `λ = 2|h|/J ∈ (0, 1)`, every alignment `r`,
every window `[m, m + L - 1]` and `z` near `λ`, `E ‖T^z(m + L, m)‖₂² ≤ C e^{K L ((Re z - λ)² + |Im z|)}`. -/
theorem integral_hs2_tmProd_le (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) :
    ∃ C K δ₀ : ℝ, 0 < C ∧ 0 ≤ K ∧ 0 < δ₀ ∧ ∀ z : ℂ, ‖z - ((2 * |h| / J : ℝ) : ℂ)‖ ≤ δ₀ →
      ∀ (m : ℤ) (L : ℕ),
        Integrable (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) muZ ∧
        ∫ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ∂muZ ≤
          C * Real.exp (K * L * ((z.re - 2 * |h| / J) ^ 2 + |z.im|)) := by
  have h0 : 0 < 2 * |h| / J := by positivity
  have h1 : 2 * |h| / J < 1 := by rw [div_lt_one hJ]; exact hh1
  set lam := 2 * |h| / J with hlam
  have hl : 0 < 1 - lam := by linarith
  refine ⟨121 * (8 / (1 - lam)), 26 / (1 - lam), (1 - lam) / 4, by positivity, by positivity,
    by positivity, fun z hz m L => ⟨integrable_hs2_tmProd J h r m L z, ?_⟩⟩
  obtain ⟨huv, hu⟩ := near lam ((1 - lam) / 4) (by linarith) z hz
  have hε : 0 ≤ (z.re - lam) ^ 2 + |z.im| := by positivity
  have hρ : 0 ≤ 1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|) := by positivity
  have hc : 2 * h / J = lam ∨ 2 * h / J = -lam := by
    rcases abs_choice h with e | e
    · left; rw [hlam, e]
    · right; rw [hlam, e]; ring
  have hv : ∀ η y, dimerPotZ J h r η y = lam ∨ dimerPotZ J h r η y = -lam :=
    fun η y => (isDimerCfgZ_dimerPotZ J h hJ r η).1 y
  have hz1 : ‖z - (lam : ℂ)‖ ≤ 1 := hz.trans (by linarith)
  obtain ⟨k, n, L₁, L₂, hL₁, hL₂, hn, hsplit⟩ := window_split r m L
  have hpt : ∀ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ≤
      121 * hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) := by
    intro η
    rw [hsplit]
    have hA := hs2_tmProd_short lam h0 h1 _ (hv η) z hz1 (2 * k - r + 2 * n) L₂ hL₂
    have hC := hs2_tmProd_short lam h0 h1 _ (hv η) z hz1 m L₁ hL₁
    have hB := hs2_nonneg (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n))
    have hA0 := hs2_nonneg (tmProd (dimerPotZ J h r η) z (2 * k - r + 2 * n) L₂)
    have hC0 := hs2_nonneg (tmProd (dimerPotZ J h r η) z m L₁)
    have k1 := hs2_mul_le (tmProd (dimerPotZ J h r η) z (2 * k - r + 2 * n) L₂ *
      tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) (tmProd (dimerPotZ J h r η) z m L₁)
    have k2 := hs2_mul_le (tmProd (dimerPotZ J h r η) z (2 * k - r + 2 * n) L₂)
      (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n))
    have k3 := hs2_nonneg (tmProd (dimerPotZ J h r η) z (2 * k - r + 2 * n) L₂ *
      tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n))
    calc _ ≤ _ := k1
      _ ≤ (hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r + 2 * n) L₂) *
            hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n))) *
            hs2 (tmProd (dimerPotZ J h r η) z m L₁) := mul_le_mul_of_nonneg_right k2 hC0
      _ ≤ (11 * hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n))) * 11 := by
          gcongr
      _ = 121 * hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) := by ring
  have hmid : ∫ η, hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) ∂muZ ≤
      8 / (1 - lam) * (1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|)) ^ n := by
    simp_rw [tmProd_dimers J h r k z _ n]
    calc ∫ η, hs2 (List.ofFn fun i => Dz (2 * h / J) z (winR k n η i)).prod ∂muZ
        = (∑ t : Fin n → ZMod 2, hs2 (List.ofFn fun i => Dz (2 * h / J) z (t i)).prod) / 2 ^ n :=
          integral_comp_winR k n (fun t => hs2 (List.ofFn fun i => Dz (2 * h / J) z (t i)).prod)
      _ ≤ _ := by
          refine avg_hs2_le (Dz (2 * h / J) z) lam (z.re - lam) _ h0 h1 hu hρ (fun w => ?_) n
          rw [sum_Dz lam (2 * h / J) hc z (fun D => Gq lam (z.re - lam) (D *ᵥ w))]
          exact one_step lam h0 h1 z huv hu w
  calc ∫ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ∂muZ
      ≤ ∫ η, 121 * hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) ∂muZ :=
        integral_mono (integrable_hs2_tmProd J h r m L z)
          ((integrable_hs2_tmProd J h r _ _ z).const_mul 121) hpt
    _ = 121 * ∫ η, hs2 (tmProd (dimerPotZ J h r η) z (2 * k - r) (2 * n)) ∂muZ :=
        integral_const_mul _ _
    _ ≤ 121 * (8 / (1 - lam) * (1 + 26 / (1 - lam) * ((z.re - lam) ^ 2 + |z.im|)) ^ n) := by
        gcongr
    _ ≤ 121 * (8 / (1 - lam)) * Real.exp (26 / (1 - lam) * L * ((z.re - lam) ^ 2 + |z.im|)) := by
        rw [← mul_assoc]
        exact mul_le_mul_of_nonneg_left (rho_pow_le _ _ (by positivity) hε n L hn) (by positivity)

end DFL
