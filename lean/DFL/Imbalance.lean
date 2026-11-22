import DFL.OneBody

/-!
# The density-wave imbalance (paper: Theorem `thm:imbalance`)

One-body part of the imbalance theorem: the chiral anticommutation identity, simplicity of the
spectrum of the Jacobi matrix `h_q`, the infinite-time (time-averaged) imbalance
`Ī_q = N⁻¹ ∑_α a_α²` with `a_α = ⟨φ_α, Γ φ_α⟩`, and lower bounds on its Bernoulli average:
a finite-`N` bound for every `h ≠ 0`, and a bound uniform in `N` under Hypothesis (EC).
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Staggered weight `a_α = ⟨φ_α, Γ φ_α⟩ = ∑_x (-1)^x φ_α(x)²` of the eigenvector `φ_α` of `h_q`. -/
def stagWeight (J h : ℝ) (q : Cfg N) (α : Fin N) : ℝ :=
  ∑ x, (-1 : ℝ) ^ x.val * ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2

/-- The one-body imbalance in sector `q`: `N⁻¹ Re Tr(Γ U(t) Γ U(t)†)`, `U(t) = e^{-ith_q}`. -/
def imbOne (J h : ℝ) (q : Cfg N) (t : ℝ) : ℝ :=
  (1 / N : ℝ) *
    (trace (cplx (stag N) * evo (hq J h q) t * cplx (stag N) * (evo (hq J h q) t)ᴴ)).re

/-- The infinite-time imbalance in sector `q`: `N⁻¹ ∑_α a_α²`. -/
def imbInf (J h : ℝ) (q : Cfg N) : ℝ := (1 / N : ℝ) * ∑ α, stagWeight J h q α ^ 2

/-- `Γ` anticommutes with the path adjacency matrix (the path is bipartite). -/
theorem stag_pathAdj_anticomm (N : ℕ) : stag N * pathAdj N + pathAdj N * stag N = 0 := by
  ext x y
  simp only [stag, Matrix.add_apply, Matrix.diagonal_mul, Matrix.mul_diagonal, pathAdj,
    Matrix.zero_apply]
  split_ifs with hxy
  · rcases hxy with h | h
    · rw [← h, pow_succ]; ring
    · rw [← h, pow_succ]; ring
  · ring

/-- Chiral identity `Γ h_q + h_q Γ = 4h V_q Γ`, with `V_q = diag(q)`. -/
theorem stag_hq_anticomm (J h : ℝ) (q : Cfg N) :
    stag N * hq J h q + hq J h q * stag N =
      (4 * h) • (diagonal (fun x => sgn (q x)) * stag N) := by
  ext x y
  have kxy := congrFun (congrFun (stag_pathAdj_anticomm N) x) y
  simp only [stag, Matrix.add_apply, Matrix.diagonal_mul, Matrix.mul_diagonal,
    Matrix.zero_apply] at kxy
  simp only [stag, hq, Matrix.add_apply, Matrix.diagonal_mul, Matrix.mul_diagonal,
    Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply]
  by_cases hxy : x = y
  · subst hxy
    have : pathAdj N x x = 0 := by simp [pathAdj]
    simp only [this, ite_true]
    ring
  · simp only [hxy, ite_false]
    linear_combination (-J) * kxy

/-- The eigenvalue identity `E_α a_α = 2h ⟨φ_α, V_q Γ φ_α⟩`. -/
theorem eigenvalue_mul_stagWeight (J h : ℝ) (q : Cfg N) (α : Fin N) :
    (hq_isHermitian J h q).eigenvalues α * stagWeight J h q α =
      2 * h * ∑ x, sgn (q x) * (-1 : ℝ) ^ x.val *
        ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2 := by
  set hH := hq_isHermitian J h q with hHdef
  set v : Fin N → ℝ := ⇑(hH.eigenvectorBasis α) with hvdef
  set E := hH.eigenvalues α with hE
  have hv : hq J h q *ᵥ v = E • v := hH.mulVec_eigenvectorBasis α
  have hsymm : v ᵥ* hq J h q = hq J h q *ᵥ v := by
    rw [← Matrix.mulVec_transpose, (hq_isSymm J h q).eq]
  have hG : v ⬝ᵥ (stag N *ᵥ v) = stagWeight J h q α := by
    simp only [stagWeight, dotProduct, stag, Matrix.mulVec_diagonal]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [← hHdef]
    ring
  have h1 : v ⬝ᵥ ((stag N * hq J h q + hq J h q * stag N) *ᵥ v) =
      2 * E * stagWeight J h q α := by
    rw [Matrix.add_mulVec, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, hv,
      dotProduct_add, Matrix.dotProduct_mulVec v (hq J h q), hsymm, hv, Matrix.mulVec_smul,
      dotProduct_smul, smul_dotProduct, hG]
    simp only [smul_eq_mul]
    ring
  have h2 : v ⬝ᵥ (((4 * h) • (diagonal (fun x => sgn (q x)) * stag N)) *ᵥ v) =
      4 * h * ∑ x, sgn (q x) * (-1 : ℝ) ^ x.val * (v x) ^ 2 := by
    simp only [stag]
    rw [Matrix.diagonal_mul_diagonal, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]
    simp only [dotProduct, Matrix.mulVec_diagonal, Finset.mul_sum]
    refine Finset.sum_congr rfl fun x _ => ?_
    ring
  rw [stag_hq_anticomm] at h1
  rw [h1] at h2
  linarith

namespace ImbalanceAux

/-- Orthonormality of the eigenvectors: `∑_x φ_α(x) φ_β(x) = δ_{αβ}`. -/
lemma eig_orth (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α β : Fin N) :
    ∑ x, hA.eigenvectorBasis α x * hA.eigenvectorBasis β x = if α = β then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_star_mul_self hA.eigenvectorUnitary) α) β
  simp only [Matrix.mul_apply, Matrix.one_apply, Matrix.star_apply,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial] at h
  exact h

/-- Completeness of the eigenvectors: `∑_α φ_α(x) φ_α(y) = δ_{xy}`. -/
lemma eig_compl (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    ∑ α, hA.eigenvectorBasis α x * hA.eigenvectorBasis α y = if x = y then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hA.eigenvectorUnitary) x) y
  simp only [Matrix.mul_apply, Matrix.one_apply, Unitary.coe_star, Matrix.star_apply,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial] at h
  exact h

/-- The eigenvalue equation, entrywise: `∑_y A_{xy} φ_α(y) = E_α φ_α(x)`. -/
lemma eig_eq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α x : Fin N) :
    ∑ y, A x y * hA.eigenvectorBasis α y = hA.eigenvalues α * hA.eigenvectorBasis α x := by
  have h := congrFun (hA.mulVec_eigenvectorBasis α) x
  simpa [Matrix.mulVec, dotProduct] using h

/-- Superdiagonal entries of `h_q`. -/
lemma hq_apply_succ (J h : ℝ) (q : Cfg N) (x y : Fin N) (hxy : x.val + 1 = y.val) :
    hq J h q x y = -J := by
  have hne : x ≠ y := by intro e; rw [e] at hxy; omega
  simp [hq, pathAdj, hxy, diagonal_apply_ne _ hne]

/-- `h_q` vanishes above the superdiagonal. -/
lemma hq_apply_far (J h : ℝ) (q : Cfg N) (x y : Fin N) (h1 : x.val + 1 < y.val) :
    hq J h q x y = 0 := by
  have hne : x ≠ y := by intro e; rw [e] at h1; omega
  have : ¬ (x.val + 1 = y.val ∨ y.val + 1 = x.val) := by omega
  simp [hq, pathAdj, this, diagonal_apply_ne _ hne]

/-- An eigenvector of the Jacobi matrix `h_q` (`J ≠ 0`) vanishing at the first site vanishes. -/
lemma hq_eigvec_zero (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (E : ℝ) (v : Fin N → ℝ)
    (hv : ∀ x, ∑ y, hq J h q x y * v y = E * v x)
    (h0 : ∀ x : Fin N, x.val = 0 → v x = 0) : ∀ x, v x = 0 := by
  have key : ∀ k : ℕ, ∀ x : Fin N, x.val ≤ k → v x = 0 := by
    intro k
    induction k with
    | zero => intro x hx; exact h0 x (by omega)
    | succ k ih =>
      intro x hx
      rcases Nat.lt_or_ge x.val (k + 1) with hlt | hge
      · exact ih x (by omega)
      · have hxk : x.val = k + 1 := by omega
        have hk : k < N := by omega
        have hz := hv ⟨k, hk⟩
        have hvz : v ⟨k, hk⟩ = 0 := ih ⟨k, hk⟩ le_rfl
        rw [hvz, mul_zero, Finset.sum_eq_single x] at hz
        · rw [hq_apply_succ J h q ⟨k, hk⟩ x (by simp; omega)] at hz
          rcases mul_eq_zero.mp hz with h' | h'
          · exact absurd (neg_eq_zero.mp h') hJ
          · exact h'
        · intro y _ hyx
          rcases Nat.lt_or_ge y.val (k + 1) with hy | hy
          · rw [ih y (by omega), mul_zero]
          · have : k + 1 < y.val := by
              have : y.val ≠ x.val := fun e => hyx (Fin.ext e)
              omega
            rw [hq_apply_far J h q ⟨k, hk⟩ y (by simp; omega), zero_mul]
        · intro hx'; exact absurd (Finset.mem_univ x) hx'
  intro x
  exact key x.val x le_rfl

end ImbalanceAux

/-- Jacobi matrices with non-zero off-diagonal entries have simple spectrum. -/
theorem hq_eigenvalues_injective (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    Function.Injective (hq_isHermitian J h q).eigenvalues := by
  intro α β hαβ
  by_contra hne
  set hH := hq_isHermitian J h q with hHdef
  have hN : 0 < N := Fin.pos α
  set z : Fin N := ⟨0, hN⟩ with hz
  set a := hH.eigenvectorBasis α z
  set b := hH.eigenvectorBasis β z
  set w : Fin N → ℝ := fun x => b * hH.eigenvectorBasis α x - a * hH.eigenvectorBasis β x
    with hwdef
  have hw : ∀ x, ∑ y, hq J h q x y * w y = hH.eigenvalues α * w x := by
    intro x
    have e1 := ImbalanceAux.eig_eq (hq J h q) hH α x
    have e2 := ImbalanceAux.eig_eq (hq J h q) hH β x
    have : ∀ y, hq J h q x y * w y = b * (hq J h q x y * hH.eigenvectorBasis α y) -
        a * (hq J h q x y * hH.eigenvectorBasis β y) := fun y => by simp only [w]; ring
    rw [Finset.sum_congr rfl (fun y _ => this y), Finset.sum_sub_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum, e1, e2, ← hαβ]
    simp only [w]; ring
  have hw0 : ∀ x, w x = 0 := ImbalanceAux.hq_eigvec_zero J h hJ q _ w hw (fun x hx => by
    have : x = z := Fin.ext hx
    subst this; simp only [w]; ring)
  have hb : b = 0 := by
    have h1 : ∑ x, hH.eigenvectorBasis α x * w x = 0 := by simp [hw0]
    have h2 : ∑ x, hH.eigenvectorBasis α x * w x =
        b * (∑ x, hH.eigenvectorBasis α x * hH.eigenvectorBasis α x) -
          a * (∑ x, hH.eigenvectorBasis α x * hH.eigenvectorBasis β x) := by
      simp only [w, Finset.mul_sum, ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl fun x _ => by ring
    rw [ImbalanceAux.eig_orth, ImbalanceAux.eig_orth, ite_eq_left rfl, ite_eq_right hne] at h2
    linarith
  have hβ0 : ∀ x, hH.eigenvectorBasis β x = 0 :=
    ImbalanceAux.hq_eigvec_zero J h hJ q (hH.eigenvalues β) _ (ImbalanceAux.eig_eq _ hH β) (fun x hx => by
      have : x = z := Fin.ext hx
      rw [this]; exact hb)
  have hnorm := ImbalanceAux.eig_orth (hq J h q) hH β β
  simp [hβ0] at hnorm

namespace ImbalanceAux

/-- Four nested finite sums can be reordered. -/
lemma sum_comm4 {ι κ M : Type*} [Fintype ι] [Fintype κ] [AddCommMonoid M]
    (f : ι → ι → κ → κ → M) :
    ∑ x, ∑ y, ∑ a, ∑ b, f x y a b = ∑ a, ∑ b, ∑ x, ∑ y, f x y a b := by
  calc ∑ x, ∑ y, ∑ a, ∑ b, f x y a b = ∑ x, ∑ a, ∑ y, ∑ b, f x y a b :=
        Finset.sum_congr rfl fun x _ => Finset.sum_comm
    _ = ∑ a, ∑ x, ∑ y, ∑ b, f x y a b := Finset.sum_comm
    _ = ∑ a, ∑ x, ∑ b, ∑ y, f x y a b :=
        Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun x _ => Finset.sum_comm
    _ = ∑ a, ∑ b, ∑ x, ∑ y, f x y a b := Finset.sum_congr rfl fun a _ => Finset.sum_comm

/-- The phase product `e^{-itE_α} conj(e^{-itE_β})` has real part `cos(t(E_β - E_α))`. -/
lemma re_phase (t a b : ℝ) :
    (Complex.exp (-(t : ℂ) * Complex.I * (a : ℂ)) *
      star (Complex.exp (-(t : ℂ) * Complex.I * (b : ℂ)))).re = Real.cos (t * (b - a)) := by
  rw [Complex.star_def, ← Complex.exp_conj, ← Complex.exp_add]
  have : -(t : ℂ) * Complex.I * (a : ℂ) + (starRingEnd ℂ) (-(t : ℂ) * Complex.I * (b : ℂ)) =
      ((t * (b - a) : ℝ) : ℂ) * Complex.I := by
    simp only [map_mul, map_neg, Complex.conj_ofReal, Complex.conj_I]
    push_cast; ring
  rw [this, Complex.exp_ofReal_mul_I_re]

/-- Trace formula `Re Tr(Γ U Γ U†) = ∑_{αβ} G_{αβ}² cos(t(E_β - E_α))` for a propagator in
spectral form. -/
lemma trace_formula_gen (U : Matrix (Fin N) (Fin N) ℂ) (φ : Fin N → Fin N → ℝ)
    (E : Fin N → ℝ) (t : ℝ)
    (hU : ∀ x y, U x y = ∑ α, Complex.exp (-(t : ℂ) * Complex.I * (E α : ℂ)) *
      ((φ α x : ℂ) * (φ α y : ℂ))) :
    (trace (cplx (stag N) * U * cplx (stag N) * Uᴴ)).re =
      ∑ α, ∑ β, (∑ x, (-1 : ℝ) ^ x.val * φ α x * φ β x) ^ 2 * Real.cos (t * (E β - E α)) := by
  have hstag : cplx (stag N) = diagonal (fun x => (((-1 : ℝ) ^ x.val : ℝ) : ℂ)) := by
    unfold cplx stag; rw [Matrix.diagonal_map (by simp)]
  have h1 : trace (cplx (stag N) * U * cplx (stag N) * Uᴴ) =
      ∑ x, ∑ y, ((((-1 : ℝ) ^ x.val) * ((-1 : ℝ) ^ y.val) : ℝ) : ℂ) * (U x y * star (U x y)) := by
    rw [Matrix.trace]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [Matrix.diag_apply, Matrix.mul_apply]
    refine Finset.sum_congr rfl fun y _ => ?_
    rw [hstag, Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.conjTranspose_apply]
    push_cast; ring
  have hstar : ∀ x y, star (U x y) = ∑ β, star (Complex.exp (-(t : ℂ) * Complex.I * (E β : ℂ))) *
      ((φ β x : ℂ) * (φ β y : ℂ)) := by
    intro x y
    rw [hU, star_sum]
    refine Finset.sum_congr rfl fun β _ => ?_
    rw [star_mul, star_mul, Complex.star_def, Complex.conj_ofReal, Complex.conj_ofReal]
    ring
  have h2 : ∀ x y, ((((-1 : ℝ) ^ x.val) * ((-1 : ℝ) ^ y.val) : ℝ) : ℂ) * (U x y * star (U x y)) =
      ∑ α, ∑ β, Complex.exp (-(t : ℂ) * Complex.I * (E α : ℂ)) *
        star (Complex.exp (-(t : ℂ) * Complex.I * (E β : ℂ))) *
        ((((-1 : ℝ) ^ x.val * φ α x * φ β x) * ((-1 : ℝ) ^ y.val * φ α y * φ β y) : ℝ) : ℂ) := by
    intro x y
    rw [hstar x y, hU x y, Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun α _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun β _ => ?_
    push_cast; ring
  have h3 : ∀ α β, (∑ x, ∑ y, Complex.exp (-(t : ℂ) * Complex.I * (E α : ℂ)) *
        star (Complex.exp (-(t : ℂ) * Complex.I * (E β : ℂ))) *
        ((((-1 : ℝ) ^ x.val * φ α x * φ β x) * ((-1 : ℝ) ^ y.val * φ α y * φ β y) : ℝ) : ℂ)).re =
      (∑ x, (-1 : ℝ) ^ x.val * φ α x * φ β x) ^ 2 * Real.cos (t * (E β - E α)) := by
    intro α β
    have : ∑ x, ∑ y, Complex.exp (-(t : ℂ) * Complex.I * (E α : ℂ)) *
        star (Complex.exp (-(t : ℂ) * Complex.I * (E β : ℂ))) *
        ((((-1 : ℝ) ^ x.val * φ α x * φ β x) * ((-1 : ℝ) ^ y.val * φ α y * φ β y) : ℝ) : ℂ) =
        Complex.exp (-(t : ℂ) * Complex.I * (E α : ℂ)) *
          star (Complex.exp (-(t : ℂ) * Complex.I * (E β : ℂ))) *
          (((∑ x, (-1 : ℝ) ^ x.val * φ α x * φ β x) ^ 2 : ℝ) : ℂ) := by
      rw [sq, Finset.sum_mul_sum, Complex.ofReal_sum, Finset.mul_sum]
      refine Finset.sum_congr rfl fun x _ => ?_
      rw [Complex.ofReal_sum, Finset.mul_sum]
    rw [this, Complex.re_mul_ofReal, re_phase]
    ring
  rw [h1, Finset.sum_congr rfl (fun x _ => Finset.sum_congr rfl (fun y _ => h2 x y)), sum_comm4,
    Complex.re_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  exact h3 α β

/-- Time average of `cos(ωt)`: tends to `1` if `ω = 0` and to `0` otherwise. -/
lemma tendsto_avg_cos (ω : ℝ) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, Real.cos (ω * t)) atTop
      (𝓝 (if ω = 0 then 1 else 0)) := by
  by_cases hω : ω = 0
  · subst hω
    simp only [zero_mul, Real.cos_zero, intervalIntegral.integral_const, sub_zero, smul_eq_mul,
      mul_one, ite_true]
    refine tendsto_const_nhds.congr' ?_
    filter_upwards [eventually_gt_atTop 0] with T hT
    field_simp
  · rw [ite_eq_right hω]
    have hI : ∀ T, ∫ t in (0 : ℝ)..T, Real.cos (ω * t) = Real.sin (ω * T) / ω := by
      intro T
      rw [intervalIntegral.integral_comp_mul_left (fun t => Real.cos t) hω, integral_cos,
        mul_zero, Real.sin_zero, sub_zero, smul_eq_mul]
      field_simp
    simp_rw [hI]
    have hω' : 0 < |ω| := abs_pos.mpr hω
    have hbound : ∀ᶠ T in atTop, ‖(1 / T) * (Real.sin (ω * T) / ω)‖ ≤ (1 / |ω|) * T⁻¹ := by
      filter_upwards [eventually_gt_atTop 0] with T hT
      rw [Real.norm_eq_abs, abs_mul, abs_of_pos (by positivity : 0 < 1 / T), abs_div]
      have hs := Real.abs_sin_le_one (ω * T)
      calc 1 / T * (|Real.sin (ω * T)| / |ω|) ≤ 1 / T * (1 / |ω|) := by gcongr
        _ = (1 / |ω|) * T⁻¹ := by ring
    have hlim0 : Tendsto (fun T : ℝ => (1 / |ω|) * T⁻¹) atTop (𝓝 0) := by
      simpa using tendsto_inv_atTop_zero.const_mul (1 / |ω|)
    exact squeeze_zero_norm' hbound hlim0

end ImbalanceAux

/-- The time average of the one-body imbalance converges to `Ī_q = N⁻¹ ∑_α a_α²`. -/
theorem imbOne_timeAvg (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h q t) atTop
      (𝓝 (imbInf J h q)) := by
  set hH := hq_isHermitian J h q with hHdef
  set E := hH.eigenvalues with hE
  set G : Fin N → Fin N → ℝ := fun α β =>
    ∑ x, (-1 : ℝ) ^ x.val * hH.eigenvectorBasis α x * hH.eigenvectorBasis β x with hG
  have hform : ∀ t, imbOne J h q t =
      ∑ α, ∑ β, (1 / N : ℝ) * G α β ^ 2 * Real.cos ((E β - E α) * t) := by
    intro t
    rw [imbOne, ImbalanceAux.trace_formula_gen (evo (hq J h q) t)
      (fun α x => hH.eigenvectorBasis α x) E t (OneBodyAux.evo_apply_eq_sum _ hH t), Finset.mul_sum]
    refine Finset.sum_congr rfl fun α _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun β _ => ?_
    simp only [G]
    rw [mul_comm t (E β - E α)]
    ring
  have hcont : ∀ (c ω : ℝ), Continuous (fun t : ℝ => c * Real.cos (ω * t)) := by
    intro c ω; fun_prop
  have hint : ∀ T, ∫ t in (0 : ℝ)..T, imbOne J h q t =
      ∑ α, ∑ β, (1 / N : ℝ) * G α β ^ 2 * ∫ t in (0 : ℝ)..T, Real.cos ((E β - E α) * t) := by
    intro T
    simp_rw [hform]
    rw [intervalIntegral.integral_finsetSum (fun α _ => ?_)]
    · refine Finset.sum_congr rfl fun α _ => ?_
      rw [intervalIntegral.integral_finsetSum (fun β _ => ?_)]
      · refine Finset.sum_congr rfl fun β _ => ?_
        exact intervalIntegral.integral_const_mul _ _
      · exact (hcont _ _).intervalIntegrable _ _
    · have : Continuous (fun t : ℝ => ∑ β, (1 / N : ℝ) * G α β ^ 2 * Real.cos ((E β - E α) * t)) := by
        fun_prop
      exact this.intervalIntegrable _ _
  have hlim : imbInf J h q =
      ∑ α, ∑ β, (1 / N : ℝ) * G α β ^ 2 * (if E β - E α = 0 then 1 else 0) := by
    rw [imbInf, Finset.mul_sum]
    refine Finset.sum_congr rfl fun α _ => ?_
    have hinj := hq_eigenvalues_injective J h hJ q
    have : ∀ β, (if E β - E α = 0 then (1 : ℝ) else 0) = if β = α then 1 else 0 := by
      intro β
      by_cases hβ : β = α
      · simp [hβ]
      · have : E β - E α ≠ 0 := fun e => hβ (hinj (sub_eq_zero.mp e))
        simp [hβ, this]
    simp only [this, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    simp only [G, stagWeight]
    congr 2
    refine Finset.sum_congr rfl fun x _ => ?_
    ring
  rw [hlim]
  simp_rw [hint, Finset.mul_sum]
  refine tendsto_finsetSum _ fun α _ => tendsto_finsetSum _ fun β _ => ?_
  have := (ImbalanceAux.tendsto_avg_cos (E β - E α)).const_mul ((1 / N : ℝ) * G α β ^ 2)
  refine this.congr fun T => ?_
  ring

namespace ImbalanceAux

/-! ### Averages over charge configurations -/

/-- There are `2^N` charge configurations. -/
lemma card_cfg : (Fintype.card (Cfg N) : ℝ) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- `∑_q 1 = 2^N`. -/
lemma sum_one_cfg : ∑ _q : Cfg N, (1 : ℝ) = 2 ^ N := by
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one, card_cfg]

/-- The average of a constant. -/
lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  rw [avgQ, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, card_cfg]
  field_simp

/-- Monotonicity of the average. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (hfg : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  unfold avgQ
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => hfg q) (by positivity)

/-- Additivity of the average. -/
lemma avgQ_add (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  unfold avgQ; rw [Finset.sum_add_distrib]; ring

/-- The average commutes with subtraction. -/
lemma avgQ_sub (f g : Cfg N → ℝ) : avgQ (fun q => f q - g q) = avgQ f - avgQ g := by
  unfold avgQ; rw [Finset.sum_sub_distrib]; ring

/-- Homogeneity of the average. -/
lemma avgQ_mul_left (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  unfold avgQ; rw [← Finset.mul_sum]; ring

/-- The average commutes with finite sums. -/
lemma avgQ_sum {ι : Type*} (s : Finset ι) (f : ι → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  unfold avgQ; rw [Finset.sum_comm, Finset.mul_sum]

/-- Jensen (Cauchy–Schwarz) for the average: `(E f)² ≤ E f²`. -/
lemma sq_avgQ_le (f : Cfg N → ℝ) : avgQ f ^ 2 ≤ avgQ (fun q => f q ^ 2) := by
  have h0 : 0 ≤ avgQ (fun q => (f q - avgQ f) ^ 2) := by
    have := avgQ_mono (f := fun _ => (0 : ℝ)) (g := fun q => (f q - avgQ f) ^ 2)
      (fun q => sq_nonneg _)
    rwa [avgQ_const] at this
  have : avgQ (fun q => (f q - avgQ f) ^ 2) = avgQ (fun q => f q ^ 2) - avgQ f ^ 2 := by
    have e : ∀ q, (f q - avgQ f) ^ 2 = f q ^ 2 + ((-2 * avgQ f) * f q + avgQ f ^ 2) :=
      fun q => by ring
    simp_rw [e]
    rw [avgQ_add, avgQ_add, avgQ_mul_left, avgQ_const]
    ring
  linarith

/-! ### Rademacher variables `r_x(q) = sgn(q x) (-1)^x` -/

/-- Flipping a label flips its sign. -/
lemma sgn_add_one (a : ZMod 2) : sgn (a + 1) = -sgn a := by
  have key : ∀ b : ZMod 2, (b + 1 = 0 ↔ ¬ b = 0) := by decide
  unfold sgn
  by_cases ha : a = 0
  · rw [ite_eq_right ((key a).not.mpr (not_not.mpr ha)), ite_eq_left ha]
  · rw [ite_eq_left ((key a).mpr ha), ite_eq_right ha, neg_neg]

/-- `sgn a ^ 2 = 1`. -/
lemma sgn_sq (a : ZMod 2) : sgn a ^ 2 = 1 := by
  unfold sgn; split_ifs <;> norm_num

/-- `|sgn a| = 1`. -/
lemma abs_sgn (a : ZMod 2) : |sgn a| = 1 := by
  unfold sgn; split_ifs <;> norm_num

/-- The Rademacher variable `r_x(q) = sgn(q x) (-1)^x`. -/
def rad (q : Cfg N) (x : Fin N) : ℝ := sgn (q x) * (-1 : ℝ) ^ x.val

/-- `r_x(q)² = 1`. -/
lemma rad_sq (q : Cfg N) (x : Fin N) : rad q x ^ 2 = 1 := by
  rw [rad, mul_pow, sgn_sq, ← pow_mul, mul_comm x.val 2, pow_mul]; norm_num

/-- `|r_x(q)| = 1`. -/
lemma abs_rad (q : Cfg N) (x : Fin N) : |rad q x| = 1 := by
  rw [rad, abs_mul, abs_sgn, abs_pow, abs_neg, abs_one, one_pow, one_mul]

/-- Flipping `q` at `x` changes the sign of `r_x`. -/
lemma rad_flip_self (q : Cfg N) (x : Fin N) : rad (q + Pi.single x 1) x = -rad q x := by
  simp only [rad, Pi.add_apply, Pi.single_eq_same, sgn_add_one]; ring

/-- Flipping `q` at `x` does not change `r_y` for `y ≠ x`. -/
lemma rad_flip_ne (q : Cfg N) {x y : Fin N} (hyx : y ≠ x) : rad (q + Pi.single x 1) y = rad q y := by
  simp [rad, Pi.single_eq_of_ne hyx]

/-- If `g` does not depend on the label at `x`, then `∑_q r_x(q) g(q) = 0`. -/
lemma sum_rad_mul_eq_zero (x : Fin N) (g : Cfg N → ℝ) (hg : ∀ q, g (q + Pi.single x 1) = g q) :
    ∑ q, rad q x * g q = 0 := by
  have h := Equiv.sum_comp (Equiv.addRight (Pi.single x 1 : Cfg N)) (fun q => rad q x * g q)
  simp only [Equiv.coe_addRight] at h
  have h2 : ∀ q : Cfg N, rad (q + Pi.single x 1) x * g (q + Pi.single x 1) = -(rad q x * g q) := by
    intro q
    rw [hg, rad_flip_self]; ring
  simp only [h2, Finset.sum_neg_distrib] at h
  linarith

/-- The block sum `S_B(q) = ∑_{x ∈ B} r_x(q)` does not depend on labels outside `B`. -/
lemma sum_rad_flip (B : Finset (Fin N)) {x : Fin N} (hx : x ∉ B) (q : Cfg N) :
    ∑ y ∈ B, rad (q + Pi.single x 1) y = ∑ y ∈ B, rad q y :=
  Finset.sum_congr rfl fun _ hy => rad_flip_ne q (fun e => hx (e ▸ hy))

/-- Second moment: `∑_q S_B(q)² = 2^N |B|`. -/
lemma moment2 (B : Finset (Fin N)) :
    ∑ q : Cfg N, (∑ x ∈ B, rad q x) ^ 2 = 2 ^ N * B.card := by
  induction B using Finset.induction_on with
  | empty => simp
  | insert x B hx ih =>
    rw [Finset.card_insert_of_notMem hx]
    have e : ∀ q : Cfg N, (∑ y ∈ insert x B, rad q y) ^ 2 =
        1 + 2 * (rad q x * (∑ y ∈ B, rad q y) ^ 1) + (∑ y ∈ B, rad q y) ^ 2 := by
      intro q
      rw [Finset.sum_insert hx]
      linear_combination rad_sq q x
    simp_rw [e]
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib, sum_one_cfg, ← Finset.mul_sum,
      sum_rad_mul_eq_zero x _ (fun q => by rw [sum_rad_flip B hx q]), ih]
    push_cast; ring

/-- Fourth moment: `∑_q S_B(q)⁴ = 2^N (3|B|² - 2|B|)`. -/
lemma moment4 (B : Finset (Fin N)) :
    ∑ q : Cfg N, (∑ x ∈ B, rad q x) ^ 4 = 2 ^ N * (3 * (B.card : ℝ) ^ 2 - 2 * B.card) := by
  induction B using Finset.induction_on with
  | empty => simp
  | insert x B hx ih =>
    rw [Finset.card_insert_of_notMem hx]
    have e : ∀ q : Cfg N, (∑ y ∈ insert x B, rad q y) ^ 4 =
        1 + 4 * (rad q x * (∑ y ∈ B, rad q y) ^ 1) + 6 * (∑ y ∈ B, rad q y) ^ 2 +
          4 * (rad q x * (∑ y ∈ B, rad q y) ^ 3) + (∑ y ∈ B, rad q y) ^ 4 := by
      intro q
      rw [Finset.sum_insert hx]
      linear_combination (rad q x ^ 2 + 1 + 4 * rad q x * (∑ y ∈ B, rad q y) +
        6 * (∑ y ∈ B, rad q y) ^ 2) * rad_sq q x
    simp_rw [e]
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
      Finset.sum_add_distrib, sum_one_cfg, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum,
      sum_rad_mul_eq_zero x _ (fun q => by rw [sum_rad_flip B hx q]),
      sum_rad_mul_eq_zero x _ (fun q => by rw [sum_rad_flip B hx q]), moment2, ih]
    push_cast; ring

/-! ### Row sums of `h_q²` -/

/-- At most one site lies to the right of `x`. -/
lemma sum_succ_le_one (x : Fin N) :
    ∑ y : Fin N, (if x.val + 1 = y.val then (1 : ℝ) else 0) ≤ 1 := by
  rw [Finset.sum_boole]
  have : (Finset.univ.filter (fun y : Fin N => x.val + 1 = y.val)).card ≤ 1 :=
    Finset.card_le_one.mpr fun a ha b hb => by
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
      exact Fin.ext (by omega)
  exact_mod_cast this

/-- At most one site lies to the left of `x`. -/
lemma sum_pred_le_one (x : Fin N) :
    ∑ y : Fin N, (if y.val + 1 = x.val then (1 : ℝ) else 0) ≤ 1 := by
  rw [Finset.sum_boole]
  have : (Finset.univ.filter (fun y : Fin N => y.val + 1 = x.val)).card ≤ 1 :=
    Finset.card_le_one.mpr fun a ha b hb => by
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
      exact Fin.ext (by omega)
  exact_mod_cast this

/-- Row sums of squares of `h_q`: `(h_q²)_{xx} = ∑_y (h_q)_{xy}² ≤ 2J² + 4h²`. -/
lemma hq_row_sq_le (J h : ℝ) (q : Cfg N) (x : Fin N) :
    ∑ y, (hq J h q x y) ^ 2 ≤ 2 * J ^ 2 + 4 * h ^ 2 := by
  have hpt : ∀ y, (hq J h q x y) ^ 2 ≤
      J ^ 2 * ((if x.val + 1 = y.val then 1 else 0) + (if y.val + 1 = x.val then 1 else 0)) +
        4 * h ^ 2 * (if x = y then 1 else 0) := by
    intro y
    by_cases hxy : x = y
    · subst hxy
      have h1 : pathAdj N x x = 0 := by simp [pathAdj]
      have h2 : ¬ (x.val + 1 = x.val) := by omega
      simp only [hq, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, h1,
        Matrix.diagonal_apply_eq, h2, ite_false, ite_true]
      nlinarith [sgn_sq (q x), sq_nonneg J]
    · have hd : diagonal (fun x => sgn (q x)) x y = 0 := diagonal_apply_ne _ hxy
      simp only [hq, pathAdj, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, hd, mul_zero,
        add_zero, hxy, ite_false]
      split_ifs <;> first | nlinarith [sq_nonneg J] | (exfalso; omega)
  calc ∑ y, (hq J h q x y) ^ 2
      ≤ ∑ y, (J ^ 2 * ((if x.val + 1 = y.val then 1 else 0) +
          (if y.val + 1 = x.val then 1 else 0)) + 4 * h ^ 2 * (if x = y then 1 else 0)) :=
        Finset.sum_le_sum fun y _ => hpt y
    _ = J ^ 2 * (∑ y : Fin N, (if x.val + 1 = y.val then (1 : ℝ) else 0) +
          ∑ y : Fin N, (if y.val + 1 = x.val then (1 : ℝ) else 0)) + 4 * h ^ 2 := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_add_distrib]
        simp
    _ ≤ J ^ 2 * (1 + 1) + 4 * h ^ 2 := by
        gcongr
        · exact sum_succ_le_one x
        · exact sum_pred_le_one x
    _ = 2 * J ^ 2 + 4 * h ^ 2 := by ring

/-- `∑_α (E_α φ_α(x))² = ∑_y A_{xy}²` (the diagonal of `A²` in the eigenbasis). -/
lemma sum_eig_sq_mul (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x : Fin N) :
    ∑ α, (hA.eigenvalues α * hA.eigenvectorBasis α x) ^ 2 = ∑ y, (A x y) ^ 2 := by
  have e : ∀ α, (hA.eigenvalues α * hA.eigenvectorBasis α x) ^ 2 =
      ∑ y, ∑ z, A x y * A x z * (hA.eigenvectorBasis α y * hA.eigenvectorBasis α z) := by
    intro α
    rw [← eig_eq A hA α x, sq, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun z _ => by ring
  simp_rw [e]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [Finset.sum_comm]
  simp_rw [← Finset.mul_sum, eig_compl A hA]
  simp [sq]

/-- `∑_α φ_α(x)² = 1`. -/
lemma sum_eig_sq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x : Fin N) :
    ∑ α, (hA.eigenvectorBasis α x) ^ 2 = 1 := by
  have := eig_compl A hA x x
  simp only [ite_true] at this
  rw [← this]
  exact Finset.sum_congr rfl fun α _ => sq _

/-- `∑_x φ_α(x)² = 1`. -/
lemma sum_eig_sq' (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α : Fin N) :
    ∑ x, (hA.eigenvectorBasis α x) ^ 2 = 1 := by
  have := eig_orth A hA α α
  simp only [ite_true] at this
  rw [← this]
  exact Finset.sum_congr rfl fun x _ => sq _

/-- `∑_α E_α a_α = 2h ∑_x r_x`. -/
lemma sum_eig_mul_stag (J h : ℝ) (q : Cfg N) :
    ∑ α, (hq_isHermitian J h q).eigenvalues α * stagWeight J h q α = 2 * h * ∑ x, rad q x := by
  simp_rw [eigenvalue_mul_stagWeight, ← Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [← Finset.mul_sum, sum_eig_sq _ (hq_isHermitian J h q) x, mul_one, rad]

/-- `∑_α E_α² ≤ N (2J² + 4h²)`. -/
lemma sum_eig_sq_le (J h : ℝ) (q : Cfg N) :
    ∑ α, (hq_isHermitian J h q).eigenvalues α ^ 2 ≤ N * (2 * J ^ 2 + 4 * h ^ 2) := by
  set hH := hq_isHermitian J h q
  have e : ∀ α, hH.eigenvalues α ^ 2 = ∑ x, (hH.eigenvalues α * hH.eigenvectorBasis α x) ^ 2 := by
    intro α
    calc hH.eigenvalues α ^ 2 = hH.eigenvalues α ^ 2 * ∑ x, (hH.eigenvectorBasis α x) ^ 2 := by
          rw [sum_eig_sq' _ hH α, mul_one]
      _ = ∑ x, (hH.eigenvalues α * hH.eigenvectorBasis α x) ^ 2 := by
          rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun x _ => by ring
  simp_rw [e]
  rw [Finset.sum_comm]
  calc ∑ x, ∑ α, (hH.eigenvalues α * hH.eigenvectorBasis α x) ^ 2
      = ∑ x, ∑ y, (hq J h q x y) ^ 2 := Finset.sum_congr rfl fun x _ => sum_eig_sq_mul _ hH x
    _ ≤ ∑ _x : Fin N, (2 * J ^ 2 + 4 * h ^ 2) :=
        Finset.sum_le_sum fun x _ => hq_row_sq_le J h q x
    _ = N * (2 * J ^ 2 + 4 * h ^ 2) := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-- `E (∑_x r_x)² = N`. -/
lemma avgQ_sq_rad_univ : avgQ (fun q : Cfg N => (∑ x, rad q x) ^ 2) = N := by
  rw [avgQ, moment2, Finset.card_univ, Fintype.card_fin]
  field_simp

end ImbalanceAux

/-- **Finite-volume lower bound** for every `h ≠ 0`:
`2^{-N} ∑_q Ī_q ≥ h² / (N (|J| + |h|)²)`. -/
theorem avg_imbInf_ge (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) (hN : 1 ≤ N) :
    h ^ 2 / (N * (|J| + |h|) ^ 2) ≤ avgQ (fun q : Cfg N => imbInf J h q) := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  set M : ℝ := N * (2 * J ^ 2 + 4 * h ^ 2) with hMdef
  have hq : ∀ q : Cfg N, 4 * h ^ 2 * (∑ x, ImbalanceAux.rad q x) ^ 2 ≤
      M * ∑ α, stagWeight J h q α ^ 2 := by
    intro q
    have cs := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ
      (fun α => (hq_isHermitian J h q).eigenvalues α) (fun α => stagWeight J h q α)
    rw [ImbalanceAux.sum_eig_mul_stag] at cs
    have hE := ImbalanceAux.sum_eig_sq_le J h q
    have ha : 0 ≤ ∑ α, stagWeight J h q α ^ 2 := Finset.sum_nonneg fun _ _ => sq_nonneg _
    calc 4 * h ^ 2 * (∑ x, ImbalanceAux.rad q x) ^ 2 = (2 * h * ∑ x, ImbalanceAux.rad q x) ^ 2 := by ring
      _ ≤ (∑ α, (hq_isHermitian J h q).eigenvalues α ^ 2) * ∑ α, stagWeight J h q α ^ 2 := cs
      _ ≤ M * ∑ α, stagWeight J h q α ^ 2 := mul_le_mul_of_nonneg_right hE ha
  have havg : 4 * h ^ 2 * N ≤ M * avgQ (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2) := by
    have := ImbalanceAux.avgQ_mono hq
    rwa [ImbalanceAux.avgQ_mul_left, ImbalanceAux.avgQ_mul_left, ImbalanceAux.avgQ_sq_rad_univ] at this
  have hA0 : 0 ≤ avgQ (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2) := by
    have := ImbalanceAux.avgQ_mono (f := fun _ : Cfg N => (0 : ℝ))
      (g := fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2)
      (fun q => Finset.sum_nonneg fun _ _ => sq_nonneg _)
    rwa [ImbalanceAux.avgQ_const] at this
  have heq : avgQ (fun q : Cfg N => imbInf J h q) =
      (1 / N) * avgQ (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2) := by
    simp only [imbInf]; exact ImbalanceAux.avgQ_mul_left _ _
  rw [heq]
  set A := avgQ (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2)
  have hJh : 2 * J ^ 2 + 4 * h ^ 2 ≤ 4 * (|J| + |h|) ^ 2 := by
    nlinarith [abs_nonneg J, abs_nonneg h, sq_abs J, sq_abs h, mul_nonneg (abs_nonneg J) (abs_nonneg h)]
  have h4 : 4 * h ^ 2 ≤ (2 * J ^ 2 + 4 * h ^ 2) * A := by
    have : N * (4 * h ^ 2) ≤ N * ((2 * J ^ 2 + 4 * h ^ 2) * A) := by
      rw [hMdef] at havg; linarith
    exact le_of_mul_le_mul_left this hN'
  have hpos : 0 < (N : ℝ) * (|J| + |h|) ^ 2 := by
    have : 0 < |J| := abs_pos.mpr hJ
    positivity
  rw [div_le_iff₀ hpos]
  have : 1 / (N : ℝ) * A * (N * (|J| + |h|) ^ 2) = A * (|J| + |h|) ^ 2 := by
    field_simp
  rw [this]
  nlinarith

/-- **Uniform lower bound** under Hypothesis (EC): the averaged infinite-time imbalance is bounded
below by a positive constant independent of the volume. -/
theorem avg_imbInf_uniform (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), 1 ≤ N → c ≤ avgQ (fun q : Cfg N => imbInf J h q) := by
  sorry

end DFL
