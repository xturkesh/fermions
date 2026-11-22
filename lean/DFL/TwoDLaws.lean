import DFL.TwoD
import DFL.Evolution

/-!
# Checkerboard propagator and the diagonal of `K_s²` (paper: Proposition `prop:2dlaws`(c)(iii),(iv))

For a bipartite graph with sublattice signs `ε` and the checkerboard charges `q⁰ = -ε`, the sector
operator is `h = h_{q⁰,s} = K_s - 2h Γ_ε` with `h² = K_s² + 4h²` (`FinGraph.checkerboard_sq`). With
`R = (K_s² + 4h²)^{1/2}`:

* (iii) `e^{-ith} = cos(tR) - i h R⁻¹ sin(tR)`: on every eigenvector `v` of `K_s²` with eigenvalue `μ`,
  `e^{-ith} v = cos(tω) v - i ω⁻¹ sin(tω) h v`, `ω = √(μ + 4h²)` (`checkerboard_evo_eigvec`); as an
  operator identity, with `cos(tR) = ∑_k (-1)^k t^{2k} (R²)^k / (2k)!` and
  `R⁻¹ sin(tR) = ∑_k (-1)^k t^{2k+1} (R²)^k / (2k+1)!` (`checkerboard_evo`);
* (iv) on a simple graph and for signs `s ∈ {±1}^E`, `(K_s²)_{ii} = J² deg(i)` (`Ks_sq_diag`).
-/

open Matrix Finset

noncomputable section

namespace DFL

namespace TwoDLawsAux

/-! ### Power series for `exp (-itH)` when `H² = X` -/

section Series

open scoped Matrix.Norms.Operator

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The even terms of the exponential series of `-itH`, for `H² = X`. -/
lemma expSeries_even (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) (k : ℕ) :
    (((2 * k).factorial : ℂ)⁻¹) • ((-(t : ℂ) * Complex.I) • H) ^ (2 * k) =
      (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k := by
  have hc2 : (-(t : ℂ) * Complex.I) ^ (2 * k) = (-1) ^ k * (t : ℂ) ^ (2 * k) := by
    rw [pow_mul, pow_mul, mul_pow, neg_sq, Complex.I_sq]
    ring
  rw [smul_pow, smul_smul, hc2, pow_mul H 2 k, sq H, hHX]
  congr 1
  ring

/-- The odd terms of the exponential series of `-itH`, for `H² = X`. -/
lemma expSeries_odd (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) (k : ℕ) :
    (((2 * k + 1).factorial : ℂ)⁻¹) • ((-(t : ℂ) * Complex.I) • H) ^ (2 * k + 1) =
      (-Complex.I) • (H * ((((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) /
        ((2 * k + 1).factorial : ℂ))) • X ^ k)) := by
  have hc2 : (-(t : ℂ) * Complex.I) ^ (2 * k) = (-1) ^ k * (t : ℂ) ^ (2 * k) := by
    rw [pow_mul, pow_mul, mul_pow, neg_sq, Complex.I_sq]
    ring
  have hH2 : H ^ (2 * k) = X ^ k := by rw [pow_mul, sq, hHX]
  rw [smul_pow, smul_smul, pow_succ', pow_succ', hc2, hH2, Matrix.mul_smul, smul_smul]
  congr 1
  ring

/-- The norms of the terms of the series `∑ (-1)^k t^{2k}/(2k)! X^k` are summable (`X = H²`). -/
lemma summable_norm_cosSeries (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) :
    Summable (fun k : ℕ =>
      ‖(((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k‖) := by
  have hnorm := NormedSpace.norm_expSeries_summable' (𝕂 := ℂ) ((-(t : ℂ) * Complex.I) • H)
  have h := hnorm.comp_injective (i := fun k => 2 * k) (mul_right_injective₀ (two_ne_zero' ℕ))
  simp only [Function.comp_def, expSeries_even H X hHX t] at h
  exact h

/-- The series `cos(tR) = ∑ (-1)^k t^{2k}/(2k)! X^k` converges (`X = R² = H²`). -/
lemma summable_cosSeries (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) :
    Summable (fun k : ℕ =>
      (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k) :=
  (summable_norm_cosSeries H X hHX t).of_norm

/-- The series `R⁻¹ sin(tR) = ∑ (-1)^k t^{2k+1}/(2k+1)! X^k` converges (`X = R² = H²`). -/
lemma summable_sincSeries (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) :
    Summable (fun k : ℕ =>
      (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) / ((2 * k + 1).factorial : ℂ))) • X ^ k) := by
  refine Summable.of_norm_bounded ((summable_norm_cosSeries H X hHX t).mul_left |t|)
    (fun k => ?_)
  have hb : (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) / ((2 * k + 1).factorial : ℂ))) • X ^ k =
      ((t : ℂ) * (((2 * k + 1 : ℕ) : ℂ))⁻¹) •
        ((((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k) := by
    rw [smul_smul]
    congr 1
    rw [Nat.factorial_succ]
    push_cast
    field_simp
    ring
  rw [hb, norm_smul, norm_mul, norm_inv, Complex.norm_real, Complex.norm_natCast,
    Real.norm_eq_abs]
  have h1 : ((2 * k + 1 : ℕ) : ℝ)⁻¹ ≤ 1 := inv_le_one_of_one_le₀ (by norm_cast; omega)
  have h2 := norm_nonneg ((((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k)
  calc |t| * ((2 * k + 1 : ℕ) : ℝ)⁻¹ *
        ‖(((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k‖
      ≤ |t| * 1 * ‖(((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k‖ := by
        gcongr
    _ = |t| * ‖(((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k‖ := by
        ring

/-- **`e^{-itH} = cos(tR) - i H R⁻¹ sin(tR)`** for `H² = R² = X`, with `cos(tR)` and
`R⁻¹ sin(tR)` given by their power series in `X`. -/
lemma exp_eq_cos_sub_sinc (H X : Matrix ι ι ℂ) (hHX : H * H = X) (t : ℝ) :
    NormedSpace.exp ((-(t : ℂ) * Complex.I) • H) =
      (∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k) -
        Complex.I • (H * ∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) /
          ((2 * k + 1).factorial : ℂ))) • X ^ k) := by
  have hexp : HasSum (fun n : ℕ => ((n.factorial : ℂ)⁻¹) • ((-(t : ℂ) * Complex.I) • H) ^ n)
      (NormedSpace.exp ((-(t : ℂ) * Complex.I) • H)) :=
    NormedSpace.exp_series_hasSum_exp' (𝕂 := ℂ) ((-(t : ℂ) * Complex.I) • H)
  have hE : HasSum
      (fun k : ℕ => (((2 * k).factorial : ℂ)⁻¹) • ((-(t : ℂ) * Complex.I) • H) ^ (2 * k))
      (∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) • X ^ k) := by
    simp only [expSeries_even H X hHX t]
    exact (summable_cosSeries H X hHX t).hasSum
  have hO : HasSum
      (fun k : ℕ => (((2 * k + 1).factorial : ℂ)⁻¹) •
        ((-(t : ℂ) * Complex.I) • H) ^ (2 * k + 1))
      ((-Complex.I) • (H * ∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) /
          ((2 * k + 1).factorial : ℂ))) • X ^ k)) := by
    simp only [expSeries_odd H X hHX t]
    exact ((summable_sincSeries H X hHX t).hasSum.mul_left H).const_smul (-Complex.I)
  rw [hexp.unique (hE.even_add_odd hO), neg_smul, sub_eq_add_neg]

end Series

/-! ### Power series on eigenvectors -/

/-- A convergent power series `∑ c_k X^k` acts on an eigenvector `X v = m v` as the scalar
`∑ c_k m^k`. -/
lemma tsum_smul_pow_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] (X : Matrix ι ι ℂ)
    (c : ℕ → ℂ) (hc : Summable (fun k => c k • X ^ k)) (m S : ℂ)
    (hS : HasSum (fun k => c k * m ^ k) S) (v : ι → ℂ) (hv : X *ᵥ v = m • v) :
    (∑' k, c k • X ^ k) *ᵥ v = S • v := by
  have hpow : ∀ k : ℕ, X ^ k *ᵥ v = m ^ k • v := by
    intro k
    induction k with
    | zero => simp
    | succ k ih =>
      rw [pow_succ', ← Matrix.mulVec_mulVec, ih, Matrix.mulVec_smul, hv, smul_smul, pow_succ]
  let g : Matrix ι ι ℂ →+ (ι → ℂ) :=
    { toFun := fun A => A *ᵥ v
      map_zero' := Matrix.zero_mulVec v
      map_add' := fun A B => Matrix.add_mulVec A B v }
  have hg : Continuous g := continuous_id.matrix_mulVec continuous_const
  have h1 := hc.hasSum.map g hg
  have h3 : (g ∘ fun k => c k • X ^ k) = fun k => (c k * m ^ k) • v := by
    funext k
    show (c k • X ^ k) *ᵥ v = (c k * m ^ k) • v
    rw [Matrix.smul_mulVec, hpow, smul_smul]
  rw [h3] at h1
  exact h1.unique (hS.smul_const v)

/-! ### The checkerboard operator -/

variable {N M : ℕ}

/-- `K_s` is symmetric (entrywise form). -/
lemma Ks_symm (G : FinGraph N M) (J : ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.Ks J s j i = G.Ks J s i j := by
  simp only [FinGraph.Ks]
  congr 2
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  tauto

/-- `K_s` is symmetric. -/
lemma Ks_transpose (G : FinGraph N M) (J : ℝ) (s : Fin M → ℝ) : (G.Ks J s)ᵀ = G.Ks J s := by
  ext i j
  exact Ks_symm G J s i j

/-- Complexification of a sum of real matrices. -/
lemma cplx_add (A B : Matrix (Fin N) (Fin N) ℝ) : cplx (A + B) = cplx A + cplx B :=
  Matrix.map_add _ (fun a b => Complex.ofReal_add a b) A B

/-- Complexification of a product of real matrices. -/
lemma cplx_mul (A B : Matrix (Fin N) (Fin N) ℝ) : cplx (A * B) = cplx A * cplx B :=
  TwoDAux.map_ofReal_mul A B

/-- Complexification of a real multiple of the identity. -/
lemma cplx_smul_one (c : ℝ) :
    cplx (c • (1 : Matrix (Fin N) (Fin N) ℝ)) = (c : ℂ) • (1 : Matrix (Fin N) (Fin N) ℂ) :=
  TwoDAux.map_ofReal_smul_one c

/-- The checkerboard identity `h² = K_s² + 4h²` for the complexified operators. -/
lemma cplx_hqs_sq (G : FinGraph N M) (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    cplx (G.hqs J h (fun i => -ε i) s) * cplx (G.hqs J h (fun i => -ε i) s) =
      cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) := by
  rw [← G.checkerboard_sq J h s ε hε1 hε, cplx_mul]

/-- Eigenvalues of `K_s²` (with complex eigenvectors) are nonnegative. -/
lemma eigval_nonneg (G : FinGraph N M) (J : ℝ) (s : Fin M → ℝ) (μ : ℝ) (v : Fin N → ℂ)
    (hv0 : v ≠ 0) (hv : cplx (G.Ks J s * G.Ks J s) *ᵥ v = (μ : ℂ) • v) : 0 ≤ μ := by
  have hKH := TwoDAux.conjTranspose_map_ofReal _ (Ks_transpose G J s)
  have h2 := TwoDAux.star_mulVec_dotProduct _ hKH v
  have hv' : ((G.Ks J s).map (fun r : ℝ => (r : ℂ)) * (G.Ks J s).map (fun r : ℝ => (r : ℂ))) *ᵥ v =
      (μ : ℂ) • v := by
    rw [← TwoDAux.map_ofReal_mul]
    exact hv
  rw [hv', dotProduct_smul, TwoDAux.star_dotProduct_self, TwoDAux.star_dotProduct_self,
    smul_eq_mul] at h2
  have h3 : (∑ i, ‖((G.Ks J s).map (fun r : ℝ => (r : ℂ)) *ᵥ v) i‖ ^ 2) =
      μ * ∑ i, ‖v i‖ ^ 2 := by
    exact_mod_cast h2
  have hSpos : 0 < ∑ i, ‖v i‖ ^ 2 := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hv0
    exact Finset.sum_pos' (fun j _ => by positivity)
      ⟨i, Finset.mem_univ _, pow_pos (norm_pos_iff.mpr hi) 2⟩
  have hT : 0 ≤ ∑ i, ‖((G.Ks J s).map (fun r : ℝ => (r : ℂ)) *ᵥ v) i‖ ^ 2 :=
    Finset.sum_nonneg (fun j _ => by positivity)
  nlinarith

end TwoDLawsAux

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- **Checkerboard propagator on eigenvectors** (paper: Proposition `prop:2dlaws`(c)(iii)): if
`K_s² v = μ v`, then `e^{-ith} v = cos(tω) v - i ω⁻¹ sin(tω) h v` with `ω = √(μ + 4h²)`. -/
theorem checkerboard_evo_eigvec (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (t μ : ℝ) (v : Fin N → ℂ) (hv : cplx (G.Ks J s * G.Ks J s) *ᵥ v = (μ : ℂ) • v) :
    evo (G.hqs J h (fun i => -ε i) s) t *ᵥ v =
      ((Real.cos (t * √(μ + 4 * h ^ 2)) : ℝ) : ℂ) • v -
        (Complex.I * ((Real.sin (t * √(μ + 4 * h ^ 2)) / √(μ + 4 * h ^ 2) : ℝ) : ℂ)) •
          (cplx (G.hqs J h (fun i => -ε i) s) *ᵥ v) := by
  by_cases hv0 : v = 0
  · subst hv0
    simp only [Matrix.mulVec_zero, smul_zero, sub_zero]
  have hμ : 0 ≤ μ := TwoDLawsAux.eigval_nonneg G J s μ v hv0 hv
  have hpos : 0 < μ + 4 * h ^ 2 := by positivity
  set ω := √(μ + 4 * h ^ 2) with hω
  have hω0 : ω ≠ 0 := (Real.sqrt_pos.mpr hpos).ne'
  have hω2 : ω ^ 2 = μ + 4 * h ^ 2 := Real.sq_sqrt hpos.le
  have hXv : cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) *ᵥ v =
      ((μ + 4 * h ^ 2 : ℝ) : ℂ) • v := by
    rw [TwoDLawsAux.cplx_add, TwoDLawsAux.cplx_smul_one, Matrix.add_mulVec, hv,
      Matrix.smul_mulVec, Matrix.one_mulVec, ← add_smul]
    push_cast
    rfl
  have hcos : HasSum (fun k : ℕ => ((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) /
      ((2 * k).factorial : ℂ)) * ((μ + 4 * h ^ 2 : ℝ) : ℂ) ^ k) ((Real.cos (t * ω) : ℝ) : ℂ) := by
    rw [Complex.ofReal_cos]
    convert Complex.hasSum_cos ((t * ω : ℝ) : ℂ) using 1
    funext k
    rw [← hω2]
    push_cast
    ring
  have hsin : HasSum (fun k : ℕ => ((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) /
      ((2 * k + 1).factorial : ℂ)) * ((μ + 4 * h ^ 2 : ℝ) : ℂ) ^ k)
      ((Real.sin (t * ω) / ω : ℝ) : ℂ) := by
    rw [Complex.ofReal_div, Complex.ofReal_sin]
    convert (Complex.hasSum_sin ((t * ω : ℝ) : ℂ)).div_const (ω : ℂ) using 1
    funext k
    rw [← hω2]
    have hω0' : (ω : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hω0
    push_cast
    field_simp
    ring
  have hsq := TwoDLawsAux.cplx_hqs_sq G J h s ε hε1 hε
  unfold evo
  rw [TwoDLawsAux.exp_eq_cos_sub_sinc _ _ hsq t, Matrix.sub_mulVec, Matrix.smul_mulVec,
    ← Matrix.mulVec_mulVec,
    TwoDLawsAux.tsum_smul_pow_mulVec _ _ (TwoDLawsAux.summable_cosSeries _ _ hsq t) _ _ hcos v hXv,
    TwoDLawsAux.tsum_smul_pow_mulVec _ _ (TwoDLawsAux.summable_sincSeries _ _ hsq t) _ _ hsin v
      hXv,
    Matrix.mulVec_smul, smul_smul]

/-- **Checkerboard propagator** (paper: Proposition `prop:2dlaws`(c)(iii)):
`e^{-ith} = cos(tR) - i h R⁻¹ sin(tR)`, `R² = K_s² + 4h²`, with `cos(tR)` and `R⁻¹ sin(tR)` given by
their (even) power series in `R²`. -/
theorem checkerboard_evo (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (t : ℝ) :
    evo (G.hqs J h (fun i => -ε i) s) t =
      (∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) •
          cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) ^ k) -
        Complex.I • (cplx (G.hqs J h (fun i => -ε i) s) *
          ∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) / ((2 * k + 1).factorial : ℂ))) •
            cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) ^ k) :=
  TwoDLawsAux.exp_eq_cos_sub_sinc _ _ (TwoDLawsAux.cplx_hqs_sq G J h s ε hε1 hε) t

/-- **Diagonal of `K_s²`** (paper: Proposition `prop:2dlaws`(c)(iv)): on a simple graph (no parallel
edges) and for signs `s_e ∈ {±1}`, `(K_s²)_{ii} = J² deg(i)`. -/
theorem Ks_sq_diag (J : ℝ) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i : Fin N) :
    (G.Ks J s * G.Ks J s) i i =
      J ^ 2 * ((univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card : ℝ) := by
  -- `F j`: the edges joining `i` and `j` (at most one, by simplicity)
  set F : Fin N → Finset (Fin M) := fun j => univ.filter (fun e =>
    (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) with hF
  have hK : ∀ j, G.Ks J s i j = -J * ∑ e ∈ F j, s e := fun j => rfl
  have hcard : ∀ j, (F j).card ≤ 1 := by
    intro j
    refine Finset.card_le_one.mpr fun e he e' he' => ?_
    simp only [hF, Finset.mem_filter, Finset.mem_univ, true_and] at he he'
    rcases he with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> rcases he' with ⟨h3, h4⟩ | ⟨h3, h4⟩
    · exact hsimple e e' (Or.inl ⟨h1.trans h3.symm, h2.trans h4.symm⟩)
    · exact hsimple e e' (Or.inr ⟨h1.trans h3.symm, h2.trans h4.symm⟩)
    · exact hsimple e e' (Or.inr ⟨h2.trans h4.symm, h1.trans h3.symm⟩)
    · exact hsimple e e' (Or.inl ⟨h2.trans h4.symm, h1.trans h3.symm⟩)
  have hsq : ∀ j, (∑ e ∈ F j, s e) ^ 2 = ((F j).card : ℝ) := by
    intro j
    rcases Nat.le_one_iff_eq_zero_or_eq_one.mp (hcard j) with h0 | h1
    · rw [Finset.card_eq_zero.mp h0]
      simp
    · obtain ⟨a, ha⟩ := Finset.card_eq_one.mp h1
      rw [ha, Finset.sum_singleton, Finset.card_singleton]
      rcases hs a with h | h <;> rw [h] <;> norm_num
  -- the other endpoint of an edge at `i`; its fibres are the `F j`
  set o : Fin M → Fin N := fun e => if G.src e = i then G.tgt e else G.src e with ho
  have hfib : ∀ j, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).filter (fun e => o e = j)
      = F j := by
    intro j
    ext e
    simp only [hF, ho, Finset.mem_filter, Finset.mem_univ, true_and]
    have hne := G.src_ne_tgt e
    by_cases h1 : G.src e = i
    · rw [ite_eq_left h1]
      constructor
      · rintro ⟨-, h⟩
        exact Or.inl ⟨h1, h⟩
      · rintro (⟨-, h⟩ | ⟨h2, -⟩)
        · exact ⟨Or.inl h1, h⟩
        · exact absurd (h1.trans h2.symm) hne
    · rw [ite_eq_right h1]
      simp only [h1, false_or, false_and]
  rw [Matrix.mul_apply]
  calc ∑ j, G.Ks J s i j * G.Ks J s j i = ∑ j, J ^ 2 * ((F j).card : ℝ) := by
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [TwoDLawsAux.Ks_symm G J s i j, hK, ← hsq]
        ring
    _ = J ^ 2 * ((univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card : ℝ) := by
        rw [← Finset.mul_sum, Finset.card_eq_sum_card_fiberwise (f := o) (t := univ)
          (fun e _ => Finset.mem_univ _)]
        push_cast
        simp_rw [hfib]

end FinGraph

end DFL
