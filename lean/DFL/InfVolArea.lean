import DFL.InfVolNumber
import DFL.SlaterBulk

/-!
# Area law for bulk intervals in the thermodynamic limit (paper: Theorem `thm:infvol`(d))

For a finite interval `A = [a, b] ⊂ ℤ`, a Fock configuration `n` on `ℤ` and the chain `Λ_M = [-M, M]`, let
`S^{(M)}_q(t)` be the entropy of the reduced state of the modes of `A` (`rdmF`, the fermionic reduced state of
`DFL.SlaterBulk`) in the sector state `ψ^{(M)}_q(t) = Γ(e^{-ith^{(M)}_q}) |n|_{Λ_M}⟩` (`entM`), and
`S^ℤ_q(t) = ∑_i h(ν_i(K^ℤ_q(t)|_A))` (`entZ`). Then

* `S^{(M)}_q(t) → S^ℤ_q(t)` uniformly in `q` (`tendstoUniformly_entM`);
* the bulk conditional entanglement converges: `S_G^{(M)}(t) = E S^{(M)}_q(t) → E S^ℤ_q(t)` (`tendsto_SGM`);
* `E sup_t S^ℤ_q(t) ≤ 4 ln 2 Θ(F)` under (EC) (`lintegral_iSup_entZ_le`);
* for the half-infinite regions `A_M = [-M, m]`, `E S(ρ^q_{A_M}(t)) ≤ 2 ln 2 Θ(F)` for every `M` and `t`
  (`areaLaw_halfInfinite`), in particular for the `limsup`.

## Proof outline (auxiliary lemmas in `DFL.InfVolAreaAux`)

1. Quantitative finite propagation: for `x, s ∈ Λ_M`,
   `|U^{(M)}_q(t)_{xs} - U^ℤ_q(t)_{xs}| ≤ 2 (c^m/m!) e^c`, `m = M - |x| + 1`, uniformly in `q`
   (`norm_evoZ_sub_UZ_le`), so the propagator entries converge uniformly in `q`
   (`tendstoUniformly_evoZ`); dominated convergence for series, uniformly in a parameter
   (`tendstoUniformly_tsum`), gives `K^{(M)}_q(t)|_A → K^ℤ_q(t)|_A` uniformly in `q`
   (`tendstoUniformly_KMA`, paper: Step 1).
2. For a polynomial `p`, `∑_i p(ν_i(C)) = G_p(C) := ∑_k p_k Re Tr C^k`
   (`sum_eval_eigenvalues`), a continuous function of the entries, invariant under relabelling
   (`Gp_submatrix`). It is uniformly continuous on the compact set of matrices with entries in the
   unit disc, which contains `K^{(M)}_q(t)|_A` and `K^ℤ_q(t)|_A`, so
   `G_p(K^{(M)}_q(t)|_A) → G_p(K^ℤ_q(t)|_A)` uniformly in `q` (`tendstoUniformly_Gp`).
3. The eigenvalues of `K^ℤ_q(t)|_A` lie in `[0, 1]` (`eig_KZA_mem`): `0 ≤ Re⟨v, C v⟩ ≤ ‖v‖²`
   passes to the limit. By Weierstrass, `h` is uniformly approximated by polynomials on `[0, 1]`,
   and with `S^{(M)}_q(t) = ∑_i h(ν_i(C_{bulk}))` (`vonNeumann_rdmF_slater`) this gives the uniform
   convergence.
4. The averages are integrals over `muZ` of functions of `q|_{Λ_M}`; dominated convergence gives
   `tendsto_SGM`. Fatou's lemma and the finite-volume bulk area law `lgt_areaLaw_bulk` give
   `lintegral_iSup_entZ_le`; `areaLaw_halfInfinite` is `lgt_areaLaw_sharp` for the left segment
   `[-M, m]`.
-/

open Matrix Finset MeasureTheory Filter Topology

noncomputable section

namespace DFL

/-- `S^{(M)}_q(t)`: entanglement entropy of the modes of `A = [a, b]` in the sector state on `Λ_M`. -/
def entM (J h : ℝ) (n q : CfgZ) (a b : ℤ) (t : ℝ) (M : ℕ) : ℝ :=
  vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M q))) t *ᵥ fockVec (occM M n)))
    (rdmF_isHermitian _ _)

/-- `S^ℤ_q(t) = ∑_i h(ν_i(K^ℤ_q(t)|_A))`. -/
def entZ (J h : ℝ) (n q : CfgZ) (a b : ℤ) (t : ℝ) : ℝ :=
  ∑ i, Real.binEntropy ((KZA_isHermitian J h q n t a b).eigenvalues i)

namespace InfVolAreaAux

open InfVolNumberAux

/-- `|U^ℤ_q(t)_{xs}| ≤ w_c(|x - s|)`. -/
lemma norm_UZ_le_wt (J h : ℝ) (q : CfgZ) (t : ℝ) (x s : ℤ) :
    ‖UZ J h q t x s‖ ≤ wt (|t| * (2 * |J| + 2 * |h|)) (x - s).natAbs :=
  le_of_tendsto' (tendsto_norm.comp (tendsto_evoZ J h q t x s)) fun M =>
    norm_evoZ_le_wt J h q t M x s

/-- **Quantitative finite propagation**: for `x, s ∈ Λ_M`,
`|U^{(M)}_q(t)_{xs} - U^ℤ_q(t)_{xs}| ≤ 2 (c^m/m!) e^c`, `m = M - |x| + 1`, uniformly in `q`. -/
lemma norm_evoZ_sub_UZ_le (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x s : ℤ) (hx : |x| ≤ M)
    (hs : |s| ≤ M) :
    ‖evoZ J h q t M x s - UZ J h q t x s‖ ≤
      2 * ((|t| * (2 * |J| + 2 * |h|)) ^ (M - x.natAbs + 1) /
        ((M - x.natAbs + 1).factorial : ℝ) * Real.exp (|t| * (2 * |J| + 2 * |h|))) := by
  set c := |t| * (2 * |J| + 2 * |h|) with hc
  have hc0 : 0 ≤ c := by positivity
  set m := M - x.natAbs + 1 with hm
  have h1 : HasSum (fun k : ℕ => ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) *
      ((((hq J h (restrictZ M q)) ^ k) (InfVolAux.toFin M x) (InfVolAux.toFin M s) : ℝ) : ℂ))
      (evoZ J h q t M x s) := by
    have := InfVolAux.hasSum_evo_apply (hq J h (restrictZ M q)) t (InfVolAux.toFin M x)
      (InfVolAux.toFin M s)
    rw [evoZ, dite_eq_left ⟨hx, hs⟩, InfVolAux.idxZ_eq_toFin, InfVolAux.idxZ_eq_toFin]
    exact this
  have h2 := (summable_UZ J h q t x s).hasSum
  have h3 := h1.sub h2
  refine (h3.norm_le_of_bounded ((summable_tail c hc0 m).hasSum.mul_left 2) fun k => ?_).trans
    (mul_le_mul_of_nonneg_left (tail_exp_le c hc0 m) (by norm_num))
  by_cases hk : m ≤ k
  · rw [ite_eq_left hk]
    refine ((norm_sub_le _ _).trans (add_le_add
      (InfVolAux.norm_term_le t k _ _ (InfVolAux.abs_pow_hq_le J h _ k _ _))
      (InfVolAux.norm_term_le t k _ _ (InfVolAux.abs_hZpow_le J h q k x s)))).trans
      (le_of_eq ?_)
    ring
  · have hxk : |x| + k ≤ M := by
      rw [Int.abs_eq_natAbs]
      have : x.natAbs ≤ M := by rw [Int.abs_eq_natAbs] at hx; exact_mod_cast hx
      omega
    rw [ite_eq_right hk, mul_zero, InfVolAux.pow_hq_apply J h q M s hs k x hxk, sub_self, norm_zero]

/-- **Uniform convergence of the propagator entries** (paper: eq. `eq:main-conv`):
`U^{(M)}_q(t)_{xs} → U^ℤ_q(t)_{xs}` uniformly in `q`. -/
lemma tendstoUniformly_evoZ (J h : ℝ) (t : ℝ) (x s : ℤ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => evoZ J h q t M x s)
      (fun q => UZ J h q t x s) atTop := by
  set c := |t| * (2 * |J| + 2 * |h|) with hc
  have hlim : Tendsto (fun M : ℕ => 2 * (c ^ (M - x.natAbs + 1) /
      ((M - x.natAbs + 1).factorial : ℝ) * Real.exp c)) atTop (𝓝 0) := by
    have h0 := FloorSemiring.tendsto_pow_div_factorial_atTop c
    have hm : Tendsto (fun M : ℕ => M - x.natAbs + 1) atTop atTop :=
      (tendsto_add_atTop_nat 1).comp (tendsto_sub_atTop_nat x.natAbs)
    simpa using ((h0.comp hm).mul_const (Real.exp c)).const_mul 2
  rw [Metric.tendstoUniformly_iff]
  intro ε hε
  filter_upwards [hlim.eventually (gt_mem_nhds hε), eventually_ge_atTop (x.natAbs + s.natAbs)]
    with M hM1 hM2 q
  have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  have hs : |s| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  rw [dist_comm, dist_eq_norm]
  exact (norm_evoZ_sub_UZ_le J h q t M x s hx hs).trans_lt hM1

/-- **Dominated convergence for series, uniformly in a parameter**: if the terms converge
uniformly and are dominated by a summable `b`, the sums converge uniformly. -/
lemma tendstoUniformly_tsum {α β : Type*} {f : ℕ → α → β → ℂ} {g : α → β → ℂ} {b : β → ℝ}
    (hb : Summable b) (hf : ∀ M x s, ‖f M x s‖ ≤ b s) (hg : ∀ x s, ‖g x s‖ ≤ b s)
    (hlim : ∀ s, TendstoUniformly (fun M x => f M x s) (fun x => g x s) atTop) :
    TendstoUniformly (fun M x => ∑' s, f M x s) (fun x => ∑' s, g x s) atTop := by
  rw [tendstoUniformly_iff_tendsto, tendsto_uniformity_iff_dist_tendsto_zero]
  have hsf : ∀ M x, Summable (f M x) := fun M x => Summable.of_norm_bounded hb (hf M x)
  have hsg : ∀ x, Summable (g x) := fun x => Summable.of_norm_bounded hb (hg x)
  have key : Tendsto (fun p : ℕ × α => ∑' s, (g p.2 s - f p.1 p.2 s)) (atTop ×ˢ ⊤)
      (𝓝 (∑' _ : β, (0 : ℂ))) := by
    refine tendsto_tsum_of_dominated_convergence (bound := fun s => 2 * b s) (hb.mul_left 2)
      (fun s => ?_) (Eventually.of_forall fun p s => ?_)
    · have := tendstoUniformly_iff_tendsto.1 (hlim s)
      rw [tendsto_uniformity_iff_dist_tendsto_zero] at this
      rw [tendsto_zero_iff_norm_tendsto_zero]
      simpa [dist_eq_norm] using this
    · refine (norm_sub_le _ _).trans ?_
      linarith [hf p.1 p.2 s, hg p.2 s]
  rw [tsum_zero] at key
  have := (continuous_norm.tendsto 0).comp key
  simp only [norm_zero] at this
  refine this.congr fun p => ?_
  simp only [Function.comp_apply, dist_eq_norm]
  rw [(hsg p.2).tsum_sub (hsf p.1 p.2)]

/-- Uniform convergence of `F_n conj(G_n)` for uniformly convergent sequences bounded by one. -/
lemma tendstoUniformly_mul_star {α : Type*} {F G : ℕ → α → ℂ} {f g : α → ℂ}
    (hF : TendstoUniformly F f atTop) (hG : TendstoUniformly G g atTop)
    (hF1 : ∀ n x, ‖F n x‖ ≤ 1) (hg1 : ∀ x, ‖g x‖ ≤ 1) :
    TendstoUniformly (fun n x => F n x * star (G n x)) (fun x => f x * star (g x)) atTop := by
  rw [Metric.tendstoUniformly_iff] at hF hG ⊢
  intro ε hε
  filter_upwards [hF (ε / 2) (by positivity), hG (ε / 2) (by positivity)] with n h1 h2 x
  rw [dist_eq_norm]
  have e : f x * star (g x) - F n x * star (G n x) =
      (f x - F n x) * star (g x) + F n x * star (g x - G n x) := by
    rw [star_sub]; ring
  have h1' := h1 x
  have h2' := h2 x
  rw [dist_eq_norm] at h1' h2'
  rw [e]
  calc ‖(f x - F n x) * star (g x) + F n x * star (g x - G n x)‖
      ≤ ‖f x - F n x‖ * ‖g x‖ + ‖F n x‖ * ‖g x - G n x‖ := by
        refine (norm_add_le _ _).trans ?_
        rw [norm_mul, norm_mul, norm_star, norm_star]
    _ ≤ ‖f x - F n x‖ * 1 + 1 * ‖g x - G n x‖ :=
        add_le_add (mul_le_mul_of_nonneg_left (hg1 x) (norm_nonneg _))
          (mul_le_mul_of_nonneg_right (hF1 n x) (norm_nonneg _))
    _ < ε := by linarith

/-- **Convergence of the correlation matrices, uniformly in `q`** (paper: Step 1 of the proof of
`thm:infvol`): `K^{(M)}_q(t)_{xy} → K^ℤ_q(t)_{xy}` uniformly in `q`. -/
lemma tendstoUniformly_KM (J h : ℝ) (n : CfgZ) (t : ℝ) (x y : ℤ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => KM J h q n t M x y)
      (fun q => KZ J h q n t x y) atTop := by
  set c := |t| * (2 * |J| + 2 * |h|) with hc'
  have hc : 0 ≤ c := by positivity
  unfold KM KZ
  refine tendstoUniformly_tsum
    (f := fun M q s => if n s = 1 then evoZ J h q t M x s * star (evoZ J h q t M y s) else 0)
    (g := fun q s => if n s = 1 then UZ J h q t x s * star (UZ J h q t y s) else 0)
    (b := fun s => wt c (x - s).natAbs) (summable_wt_int c x) (fun M q s => ?_) (fun q s => ?_)
    (fun s => ?_)
  · split_ifs
    · rw [norm_mul, norm_star]
      calc ‖evoZ J h q t M x s‖ * ‖evoZ J h q t M y s‖ ≤ wt c (x - s).natAbs * 1 :=
            mul_le_mul (norm_evoZ_le_wt J h q t M x s) (InfVolAux.norm_evoZ_le_one J h q t M y s)
              (norm_nonneg _) (wt_nonneg hc _)
        _ = _ := mul_one _
    · rw [norm_zero]; exact wt_nonneg hc _
  · split_ifs
    · rw [norm_mul, norm_star]
      calc ‖UZ J h q t x s‖ * ‖UZ J h q t y s‖ ≤ wt c (x - s).natAbs * 1 :=
            mul_le_mul (norm_UZ_le_wt J h q t x s) (norm_UZ_le_one J h q t y s)
              (norm_nonneg _) (wt_nonneg hc _)
        _ = _ := mul_one _
    · rw [norm_zero]; exact wt_nonneg hc _
  · by_cases hs : n s = 1
    · simp only [hs, ite_true]
      exact tendstoUniformly_mul_star (tendstoUniformly_evoZ J h t x s)
        (tendstoUniformly_evoZ J h t y s) (fun M q => InfVolAux.norm_evoZ_le_one J h q t M x s)
        (fun q => norm_UZ_le_one J h q t y s)
    · simp only [hs, ite_false]
      exact tendsto_const_nhds.tendstoUniformly_const

/-- Uniform convergence of matrices follows from uniform convergence of the entries. -/
lemma tendstoUniformly_matrix {ι α m n' : Type*} {p : Filter ι} {F : ι → α → Matrix m n' ℂ}
    {f : α → Matrix m n' ℂ}
    (h : ∀ i j, TendstoUniformly (fun k x => F k x i j) (fun x => f x i j) p) :
    TendstoUniformly F f p := by
  rw [tendstoUniformly_iff_tendsto, Matrix.uniformity]
  simp only [tendsto_iInf, tendsto_comap_iff]
  intro i j
  exact tendstoUniformly_iff_tendsto.1 (h i j)

/-- `K^{(M)}_q(t)|_A → K^ℤ_q(t)|_A` uniformly in `q`. -/
lemma tendstoUniformly_KMA (J h : ℝ) (n : CfgZ) (t : ℝ) (a b : ℤ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => KMA J h q n t M a b)
      (fun q => KZA J h q n t a b) atTop :=
  tendstoUniformly_matrix fun x y => tendstoUniformly_KM J h n t x y

/-! ### Spectral calculus for Hermitian matrices -/

section Spectral

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-- `Tr C^k = ∑_i ν_i^k` for a Hermitian `C`. -/
lemma trace_pow_eq (C : Matrix ι ι ℂ) (hC : C.IsHermitian) (k : ℕ) :
    (C ^ k).trace = ∑ i, ((hC.eigenvalues i : ℂ)) ^ k := by
  have hspec := hC.spectral_theorem
  have hpow : C ^ k = Unitary.conjStarAlgAut ℂ _ hC.eigenvectorUnitary
      (diagonal (RCLike.ofReal ∘ hC.eigenvalues) ^ k) := by
    conv_lhs => rw [hspec]
    rw [map_pow]
  rw [hpow, Unitary.conjStarAlgAut_apply, trace_mul_cycle, Unitary.coe_star_mul_self, one_mul,
    diagonal_pow, trace_diagonal]
  rfl

/-- The polynomial functional `G_p(C) = ∑_k p_k Re Tr C^k`, which equals `∑_i p(ν_i)` for
Hermitian `C`. -/
def Gp (p : Polynomial ℝ) (C : Matrix ι ι ℂ) : ℝ :=
  ∑ k ∈ range (p.natDegree + 1), p.coeff k * ((C ^ k).trace).re

/-- `∑_i p(ν_i(C)) = G_p(C)` for a Hermitian `C`. -/
lemma sum_eval_eigenvalues (p : Polynomial ℝ) (C : Matrix ι ι ℂ) (hC : C.IsHermitian) :
    ∑ i, p.eval (hC.eigenvalues i) = Gp p C := by
  unfold Gp
  simp_rw [Polynomial.eval_eq_sum_range, trace_pow_eq C hC]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [← Finset.mul_sum, Complex.re_sum]
  congr 1
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [← Complex.ofReal_pow, Complex.ofReal_re]

/-- `G_p` is continuous in the entries. -/
lemma continuous_Gp (p : Polynomial ℝ) : Continuous (fun C : Matrix ι ι ℂ => Gp p C) := by
  unfold Gp
  exact continuous_finsetSum _ fun k _ =>
    continuous_const.mul (Complex.continuous_re.comp ((continuous_pow k).matrix_trace))

/-- `G_p` is invariant under relabelling rows and columns. -/
lemma Gp_submatrix (p : Polynomial ℝ) (C : Matrix ι ι ℂ) (e : κ ≃ ι) :
    Gp p (C.submatrix e e) = Gp p C := by
  unfold Gp
  refine Finset.sum_congr rfl fun k _ => ?_
  have hpow : (C.submatrix e e) ^ k = (C ^ k).submatrix e e := by
    have := map_pow (Matrix.reindexAlgEquiv ℂ ℂ e.symm) C k
    simpa [Matrix.coe_reindexAlgEquiv, Matrix.reindex_apply] using this.symm
  rw [hpow]
  congr 2
  simp only [Matrix.trace, Matrix.diag, Matrix.submatrix_apply]
  exact e.sum_comp (fun j => (C ^ k) j j)

/-- For a unitary `U`, `⟨U w, U y⟩ = ⟨w, y⟩`. -/
lemma dotProduct_unitary (U : Matrix ι ι ℂ) (hU : star U * U = 1) (w y : ι → ℂ) :
    star (U *ᵥ w) ⬝ᵥ (U *ᵥ y) = star w ⬝ᵥ y := by
  rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, ← star_eq_conjTranspose, hU, vecMul_one]

/-- Spectral form of the quadratic form of a Hermitian matrix: with `w = U^* v`,
`Re ⟨v, C v⟩ = ∑_i ν_i |w_i|²` and `Re ⟨v, v⟩ = ∑_i |w_i|²`. -/
lemma re_dotProduct_spectral (C : Matrix ι ι ℂ) (hC : C.IsHermitian) (v : ι → ℂ) :
    (star v ⬝ᵥ (C *ᵥ v)).re = ∑ i, hC.eigenvalues i *
        Complex.normSq ((star (hC.eigenvectorUnitary : Matrix ι ι ℂ) *ᵥ v) i) ∧
      (star v ⬝ᵥ v).re =
        ∑ i, Complex.normSq ((star (hC.eigenvectorUnitary : Matrix ι ι ℂ) *ᵥ v) i) := by
  set U : Matrix ι ι ℂ := (hC.eigenvectorUnitary : Matrix ι ι ℂ) with hUdef
  have hU1 : U * star U = 1 := Unitary.mul_star_self_of_mem hC.eigenvectorUnitary.2
  have hU2 : star U * U = 1 := Unitary.star_mul_self_of_mem hC.eigenvectorUnitary.2
  set w := star U *ᵥ v with hw
  have hv : v = U *ᵥ w := by rw [hw, mulVec_mulVec, hU1, one_mulVec]
  have hC' : C = U * diagonal (fun i => (hC.eigenvalues i : ℂ)) * star U :=
    hC.spectral_theorem.trans rfl
  have hCv : C *ᵥ v = U *ᵥ (diagonal (fun i => (hC.eigenvalues i : ℂ)) *ᵥ w) := by
    conv_lhs => rw [hC']
    rw [← mulVec_mulVec, ← mulVec_mulVec]
  constructor
  · rw [hCv]
    conv_lhs => rw [hv]
    rw [dotProduct_unitary U hU2]
    simp only [dotProduct, mulVec_diagonal, Pi.star_apply, Complex.re_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Complex.normSq_apply]
    simp only [Complex.mul_re, Complex.mul_im, Complex.star_def, Complex.conj_re,
      Complex.conj_im, Complex.ofReal_re, Complex.ofReal_im]
    ring
  · conv_lhs => rw [hv]
    rw [dotProduct_unitary U hU2]
    simp only [dotProduct, Pi.star_apply, Complex.re_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Complex.normSq_apply]
    simp only [Complex.mul_re, Complex.star_def, Complex.conj_re, Complex.conj_im]
    ring

/-- If the eigenvalues of a Hermitian `C` lie in `[0, 1]`, then `0 ≤ Re⟨v, C v⟩ ≤ Re⟨v, v⟩`. -/
lemma re_dotProduct_mem (C : Matrix ι ι ℂ) (hC : C.IsHermitian)
    (hν : ∀ i, 0 ≤ hC.eigenvalues i ∧ hC.eigenvalues i ≤ 1) (v : ι → ℂ) :
    0 ≤ (star v ⬝ᵥ (C *ᵥ v)).re ∧ (star v ⬝ᵥ (C *ᵥ v)).re ≤ (star v ⬝ᵥ v).re := by
  obtain ⟨h1, h2⟩ := re_dotProduct_spectral C hC v
  rw [h1, h2]
  exact ⟨Finset.sum_nonneg fun i _ => mul_nonneg (hν i).1 (Complex.normSq_nonneg _),
    Finset.sum_le_sum fun i _ => mul_le_of_le_one_left (Complex.normSq_nonneg _) (hν i).2⟩

/-- Conversely, if `0 ≤ Re⟨v, C v⟩ ≤ Re⟨v, v⟩` for all `v`, the eigenvalues lie in `[0, 1]`. -/
lemma eigenvalues_mem_of_re_dotProduct (C : Matrix ι ι ℂ) (hC : C.IsHermitian)
    (hv : ∀ v : ι → ℂ, 0 ≤ (star v ⬝ᵥ (C *ᵥ v)).re ∧
      (star v ⬝ᵥ (C *ᵥ v)).re ≤ (star v ⬝ᵥ v).re) (i : ι) :
    0 ≤ hC.eigenvalues i ∧ hC.eigenvalues i ≤ 1 := by
  set U : Matrix ι ι ℂ := (hC.eigenvectorUnitary : Matrix ι ι ℂ) with hUdef
  have hU2 : star U * U = 1 := Unitary.star_mul_self_of_mem hC.eigenvectorUnitary.2
  have hw : star U *ᵥ (U *ᵥ Pi.single i 1) = Pi.single i 1 := by
    rw [mulVec_mulVec, hU2, one_mulVec]
  obtain ⟨h1, h2⟩ := re_dotProduct_spectral C hC (U *ᵥ Pi.single i 1)
  rw [← hUdef, hw] at h1 h2
  have e1 : ∑ j, hC.eigenvalues j * Complex.normSq ((Pi.single i (1 : ℂ) : ι → ℂ) j) =
      hC.eigenvalues i := by
    rw [Finset.sum_eq_single i (fun j _ hj => by simp [hj]) (by simp)]
    simp
  have e2 : ∑ j, Complex.normSq ((Pi.single i (1 : ℂ) : ι → ℂ) j) = 1 := by
    rw [Finset.sum_eq_single i (fun j _ hj => by simp [hj]) (by simp)]
    simp
  have := hv (U *ᵥ Pi.single i 1)
  rw [h1, h2, e1, e2] at this
  exact this

omit [DecidableEq ι] [DecidableEq κ] in
/-- Relabelling the quadratic form. -/
lemma dotProduct_submatrix (C : Matrix ι ι ℂ) (e : κ ≃ ι) (v : κ → ℂ) :
    star v ⬝ᵥ ((C.submatrix e e) *ᵥ v) = star (v ∘ e.symm) ⬝ᵥ (C *ᵥ (v ∘ e.symm)) := by
  rw [submatrix_mulVec_equiv]
  simp only [dotProduct, Pi.star_apply, Function.comp_apply]
  exact (e.sum_comp (fun j => star (v (e.symm j)) * (C *ᵥ (v ∘ e.symm)) j)).symm.trans
    (by simp only [Equiv.symm_apply_apply]) |>.symm

omit [DecidableEq ι] [DecidableEq κ] in
/-- Relabelling the norm. -/
lemma dotProduct_self_comp (e : κ ≃ ι) (v : κ → ℂ) :
    star v ⬝ᵥ v = star (v ∘ e.symm) ⬝ᵥ (v ∘ e.symm) := by
  simp only [dotProduct, Pi.star_apply, Function.comp_apply]
  exact (e.sum_comp (fun j => star (v (e.symm j)) * v (e.symm j))).symm.trans
    (by simp only [Equiv.symm_apply_apply]) |>.symm

end Spectral

/-! ### The spectra of `K^{(M)}_q(t)|_A` and `K^ℤ_q(t)|_A` -/

/-- For `[a, b] ⊆ Λ_M`, `0 ≤ Re⟨v, K^{(M)}_q(t)|_A v⟩ ≤ Re⟨v, v⟩`. -/
lemma re_dotProduct_KMA_mem (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) (v : Set.Icc a b → ℂ) :
    0 ≤ (star v ⬝ᵥ (KMA J h q n t M a b *ᵥ v)).re ∧
      (star v ⬝ᵥ (KMA J h q n t M a b *ᵥ v)).re ≤ (star v ⬝ᵥ v).re := by
  rw [← corrR_submatrix_eqvM J h q n t M a b hM, dotProduct_submatrix,
    dotProduct_self_comp (eqvM M a b hM) v]
  exact re_dotProduct_mem _ (corrR_isHermitian _ _ _ _)
    (eig_corrR_mem _ (hq_isSymm J h _) _ t _) _

/-- The eigenvalues of `K^ℤ_q(t)|_A` lie in `[0, 1]`. -/
lemma eig_KZA_mem (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) (i : Set.Icc a b) :
    0 ≤ (KZA_isHermitian J h q n t a b).eigenvalues i ∧
      (KZA_isHermitian J h q n t a b).eigenvalues i ≤ 1 := by
  refine eigenvalues_mem_of_re_dotProduct _ _ (fun v => ?_) i
  have hcont : Continuous (fun C : Matrix (Set.Icc a b) (Set.Icc a b) ℂ =>
      (star v ⬝ᵥ (C *ᵥ v)).re) :=
    Complex.continuous_re.comp (continuous_const.dotProduct
      (continuous_id.matrix_mulVec continuous_const))
  have hlim := (hcont.tendsto _).comp (tendsto_KMA J h q n t a b)
  have hev : ∀ᶠ M in atTop, 0 ≤ (star v ⬝ᵥ (KMA J h q n t M a b *ᵥ v)).re ∧
      (star v ⬝ᵥ (KMA J h q n t M a b *ᵥ v)).re ≤ (star v ⬝ᵥ v).re := by
    filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
    exact re_dotProduct_KMA_mem J h q n t M a b
      ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩ v
  exact ⟨ge_of_tendsto hlim (hev.mono fun M hM => hM.1),
    le_of_tendsto hlim (hev.mono fun M hM => hM.2)⟩

/-- For `[a, b] ⊆ Λ_M`, the entries of `K^{(M)}_q(t)|_A` have modulus at most one. -/
lemma norm_KMA_le_one (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) (x y : Set.Icc a b) :
    ‖KMA J h q n t M a b x y‖ ≤ 1 := by
  rw [← corrR_submatrix_eqvM J h q n t M a b hM]
  exact TransportSharpAux.norm_corrMat_le_one _ (hq_isSymm J h _) _ t _ _

/-- The entries of `K^ℤ_q(t)|_A` have modulus at most one. -/
lemma norm_KZA_le_one (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) (x y : Set.Icc a b) :
    ‖KZA J h q n t a b x y‖ ≤ 1 := by
  refine le_of_tendsto (tendsto_norm.comp (tendsto_KM J h q n t x y)) ?_
  filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
  exact norm_KMA_le_one J h q n t M a b ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩ x y

/-- **Uniform convergence of `G_p`**: `G_p(K^{(M)}_q(t)|_A) → G_p(K^ℤ_q(t)|_A)` uniformly in `q`,
since `G_p` is uniformly continuous on the compact set of matrices with entries in the unit disc. -/
lemma tendstoUniformly_Gp (p : Polynomial ℝ) (J h : ℝ) (n : CfgZ) (t : ℝ) (a b : ℤ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => Gp p (KMA J h q n t M a b))
      (fun q => Gp p (KZA J h q n t a b)) atTop := by
  have hB : IsCompact ((Metric.closedBall (0 : ℂ) 1).matrix :
      Set (Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) := (isCompact_closedBall 0 1).matrix
  have hG : UniformContinuousOn (Gp p) ((Metric.closedBall (0 : ℂ) 1).matrix :
      Set (Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) :=
    hB.uniformContinuousOn_of_continuous (continuous_Gp p).continuousOn
  refine hG.comp_tendstoUniformly_eventually ?_ (fun q => ?_) (tendstoUniformly_KMA J h n t a b)
  · filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM q
    exact Set.mem_matrix.2 fun x y => mem_closedBall_zero_iff.2
      (norm_KMA_le_one J h q n t M a b ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩ x y)
  · exact Set.mem_matrix.2 fun x y => mem_closedBall_zero_iff.2 (norm_KZA_le_one J h q n t a b x y)

/-- Uniform polynomial approximation of `∑_i h(ν_i)`. -/
lemma abs_sum_binEntropy_sub_le {ι : Type*} [Fintype ι] (ν : ι → ℝ)
    (hν : ∀ i, 0 ≤ ν i ∧ ν i ≤ 1) (p : Polynomial ℝ) (δ : ℝ)
    (hp : ∀ x ∈ Set.Icc (0 : ℝ) 1, |p.eval x - Real.binEntropy x| < δ) :
    |∑ i, Real.binEntropy (ν i) - ∑ i, p.eval (ν i)| ≤ Fintype.card ι * δ := by
  rw [← Finset.sum_sub_distrib]
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  calc ∑ i, |Real.binEntropy (ν i) - p.eval (ν i)| ≤ ∑ _i : ι, δ :=
        Finset.sum_le_sum fun i _ => by
          rw [abs_sub_comm]; exact (hp _ ⟨(hν i).1, (hν i).2⟩).le
    _ = Fintype.card ι * δ := by rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

/-- `0 ≤ S(ρ_R) ≤ |R| ln 2` for the evolved Fock states. -/
lemma vonNeumann_rdmF_mem {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) (R : Finset (Fin N)) :
    0 ≤ vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) ∧
      vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) ≤
        R.card * Real.log 2 := by
  rw [vonNeumann_rdmF_slater A hA S t R]
  have hν := eig_corrR_mem A hA S t R
  refine ⟨Finset.sum_nonneg fun i _ => Real.binEntropy_nonneg (hν i).1 (hν i).2, ?_⟩
  calc ∑ i, Real.binEntropy ((corrR_isHermitian A S t R).eigenvalues i)
      ≤ ∑ _i : R, Real.log 2 := Finset.sum_le_sum fun i _ => Real.binEntropy_le_log_two
    _ = R.card * Real.log 2 := by
        rw [Finset.sum_const, Finset.card_univ, Fintype.card_coe, nsmul_eq_mul]

/-- For `[a, b] ⊆ Λ_M`, the bulk interval of `Λ_M` is an interval of `Fin (2M+1)`. -/
lemma bulkM_eq_ivl (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    bulkM M a b = ivl (a + M).toNat (b + M).toNat := by
  obtain ⟨h3, h4⟩ := hM
  ext i
  simp only [bulkM, ivl, Finset.mem_filter, Finset.mem_univ, true_and]
  omega

/-- For `[a, b] ⊆ Λ_M`, `|bulkM M a b| = |[a, b]|`. -/
lemma card_bulkM (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    (bulkM M a b).card = Fintype.card (Set.Icc a b) := by
  rw [← Fintype.card_coe, Fintype.card_congr (eqvM M a b hM)]

/-- Finite-volume form of the bound: `E sup_t S^{(M)}_q(t) ≤ 4 ln 2 Θ(F)` as a lower integral. -/
lemma lintegral_iSup_entM_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (a b : ℤ) (M : ℕ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    ∫⁻ q, ⨆ t : ℝ, ENNReal.ofReal (entM J h n q a b t M) ∂muZ ≤
      ENNReal.ofReal (4 * Real.log 2 * thetaF F) := by
  set f : Cfg (2 * M + 1) → ℝ → ℝ := fun q' t =>
    vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
      (rdmF_isHermitian _ _) with hf
  have hbdd : ∀ q', BddAbove (Set.range (f q')) := fun q' =>
    ⟨(bulkM M a b).card * Real.log 2, by
      rintro _ ⟨t, rfl⟩
      exact (vonNeumann_rdmF_mem _ (hq_isSymm J h q') _ t _).2⟩
  have hnn : ∀ q', 0 ≤ ⨆ t : ℝ, f q' t := fun q' =>
    Real.iSup_nonneg fun t => (vonNeumann_rdmF_mem _ (hq_isSymm J h q') _ t _).1
  calc ∫⁻ q, ⨆ t : ℝ, ENNReal.ofReal (entM J h n q a b t M) ∂muZ
      = ∫⁻ q, (fun q' => ⨆ t : ℝ, ENNReal.ofReal (f q' t)) (restrictZ M q) ∂muZ := rfl
    _ = ∑ q', (⨆ t : ℝ, ENNReal.ofReal (f q' t)) * ((2 : ENNReal) ^ (2 * M + 1))⁻¹ :=
        InfVolAux.lintegral_comp_restrictZ M (fun q' => ⨆ t : ℝ, ENNReal.ofReal (f q' t))
    _ ≤ ∑ q' : Cfg (2 * M + 1), ENNReal.ofReal (⨆ t : ℝ, f q' t) *
          ((2 : ENNReal) ^ (2 * M + 1))⁻¹ := by
        refine Finset.sum_le_sum fun q' _ => mul_le_mul_left ?_ _
        exact iSup_le fun t => ENNReal.ofReal_le_ofReal (le_ciSup (hbdd q') t)
    _ = ENNReal.ofReal (avgQ fun q' => ⨆ t : ℝ, f q' t) := by
        unfold avgQ
        rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_sum_of_nonneg fun q' _ => hnn q',
          Finset.mul_sum]
        refine Finset.sum_congr rfl fun q' _ => ?_
        rw [mul_comm]
        congr 1
        rw [one_div, ENNReal.ofReal_inv_of_pos (by positivity), ENNReal.ofReal_pow (by norm_num),
          ENNReal.ofReal_ofNat]
    _ ≤ ENNReal.ofReal (4 * Real.log 2 * thetaF F) := by
        refine ENNReal.ofReal_le_ofReal ?_
        have := lgt_areaLaw_bulk J h hEC (occM M n) (a + M).toNat (b + M).toNat
        rw [← bulkM_eq_ivl M a b hM] at this
        exact this

/-- The entanglement entropy of a left segment is bounded in time (by a time-independent
majorant). -/
lemma bddAbove_vonNeumann_rdm {N : ℕ} (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (ℓ : ℕ)
    (hℓ : ℓ ≤ N) :
    BddAbove (Set.range fun t : ℝ =>
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) :=
  ⟨2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z, by
    rintro _ ⟨t, rfl⟩
    exact (vonNeumann_slater_sharp (hq J h q) (hq_isSymm J h q) S t ℓ hℓ).trans
      (mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
        AreaLawSharpAux.norm_corrMat_le_corrQ J h q S t x y) AreaLawSharpAux.two_log_two_nonneg)⟩

end InfVolAreaAux

open InfVolNumberAux InfVolAreaAux in
/-- **`S^{(M)}_q(t) → S^ℤ_q(t)` uniformly in `q`** (paper: Theorem `thm:infvol`(d)). -/
theorem tendstoUniformly_entM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => entM J h n q a b t M)
      (fun q : CfgZ => entZ J h n q a b t) atTop := by
  rw [Metric.tendstoUniformly_iff]
  intro ε hε
  set d : ℝ := (Fintype.card (Set.Icc a b) : ℝ) with hd
  have hd0 : 0 ≤ d := Nat.cast_nonneg _
  set δ := ε / (3 * (d + 1)) with hδ
  have hδ0 : 0 < δ := by positivity
  have hdδ : d * δ < ε / 3 := by
    have e : d * δ = ε / 3 * (d / (d + 1)) := by rw [hδ]; field_simp
    rw [e]
    exact mul_lt_of_lt_one_right (by positivity) ((div_lt_one (by positivity)).2 (by linarith))
  obtain ⟨p, hp⟩ := exists_polynomial_near_of_continuousOn 0 1 Real.binEntropy
    Real.binEntropy_continuous.continuousOn δ hδ0
  filter_upwards [Metric.tendstoUniformly_iff.1 (tendstoUniformly_Gp p J h n t a b) (ε / 3)
    (by positivity), eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM1 hM2 q
  have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM2, le_of_max_le_right hM2⟩
  have hcard : (Fintype.card ↥(bulkM M a b) : ℝ) = d := by
    rw [hd, Fintype.card_congr (eqvM M a b hM')]
  have e1 : entM J h n q a b t M = ∑ i, Real.binEntropy
      ((corrR_isHermitian (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)).eigenvalues i) :=
    vonNeumann_rdmF_slater _ (hq_isSymm J h _) _ t _
  have e2 : Gp p (KMA J h q n t M a b) = ∑ i, p.eval
      ((corrR_isHermitian (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)).eigenvalues i) := by
    rw [← corrR_submatrix_eqvM J h q n t M a b hM', Gp_submatrix, sum_eval_eigenvalues]
  have e3 : Gp p (KZA J h q n t a b) =
      ∑ i, p.eval ((KZA_isHermitian J h q n t a b).eigenvalues i) :=
    (sum_eval_eigenvalues _ _ _).symm
  have b1 := abs_sum_binEntropy_sub_le _ (eig_corrR_mem (hq J h (restrictZ M q))
    (hq_isSymm J h _) (occM M n) t (bulkM M a b)) p δ hp
  have b2 := abs_sum_binEntropy_sub_le _ (eig_KZA_mem J h q n t a b) p δ hp
  rw [hcard] at b1
  rw [← hd] at b2
  have b3 := hM1 q
  rw [Real.dist_eq, e2, e3] at b3
  rw [Real.dist_eq, e1]
  unfold entZ
  rw [abs_lt] at b3 ⊢
  rw [abs_le] at b1 b2
  constructor <;> linarith [b1.1, b1.2, b2.1, b2.2, b3.1, b3.2]

open InfVolNumberAux InfVolAreaAux in
/-- **Convergence of the bulk conditional entanglement** (paper: Theorem `thm:infvol`(d)):
`S_G^{(M)}(t) = E_q S^{(M)}_q(t) → E S^ℤ_q(t)`. -/
theorem tendsto_SGM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) :
    Tendsto (fun M : ℕ => avgQ (fun q' : Cfg (2 * M + 1) =>
        vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
          (rdmF_isHermitian _ _)))
      atTop (𝓝 (∫ q, entZ J h n q a b t ∂muZ)) := by
  refine Tendsto.congr (fun M => InfVolAux.integral_comp_restrictZ M (fun q' =>
    vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
      (rdmF_isHermitian _ _))) ?_
  refine tendsto_integral_filter_of_dominated_convergence
    (fun _ => (Fintype.card (Set.Icc a b) : ℝ) * Real.log 2)
    (Eventually.of_forall fun M => (InfVolAux.measurable_comp_restrictZ M (fun q' =>
      vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
        (rdmF_isHermitian _ _))).aestronglyMeasurable)
    ?_ (integrable_const _) (ae_of_all _ fun q => (tendstoUniformly_entM J h n a b t).tendsto_at q)
  filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
  have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
  refine ae_of_all _ fun q => ?_
  have h1 := vonNeumann_rdmF_mem (hq J h (restrictZ M q)) (hq_isSymm J h _) (occM M n) t
    (bulkM M a b)
  rw [card_bulkM M a b hM'] at h1
  rw [Real.norm_eq_abs, abs_le]
  exact ⟨by linarith [h1.1, Real.log_nonneg (by norm_num : (1 : ℝ) ≤ 2),
    (Nat.cast_nonneg _ : (0 : ℝ) ≤ Fintype.card (Set.Icc a b))], h1.2⟩

open InfVolNumberAux InfVolAreaAux in
/-- **Area law for bulk intervals on `ℤ`** (paper: Theorem `thm:infvol`(d)): `E sup_t S^ℤ_q(t) ≤ 4 ln 2 Θ(F)`. -/
theorem lintegral_iSup_entZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (a b : ℤ) :
    ∫⁻ q, ⨆ t : ℝ, ENNReal.ofReal (entZ J h n q a b t) ∂muZ ≤
      ENNReal.ofReal (4 * Real.log 2 * thetaF F) := by
  set G : ℕ → CfgZ → ENNReal := fun M q => ⨆ t : ℝ, ENNReal.ofReal (entM J h n q a b t M) with hG
  have hpt : ∀ q, ⨆ t : ℝ, ENNReal.ofReal (entZ J h n q a b t) ≤
      liminf (fun M => G M q) atTop := by
    intro q
    refine iSup_le fun t => ?_
    have hlim : Tendsto (fun M => ENNReal.ofReal (entM J h n q a b t M)) atTop
        (𝓝 (ENNReal.ofReal (entZ J h n q a b t))) :=
      (ENNReal.continuous_ofReal.tendsto _).comp ((tendstoUniformly_entM J h n a b t).tendsto_at q)
    rw [← hlim.liminf_eq]
    exact liminf_le_liminf (Eventually.of_forall fun M =>
      le_iSup (fun t => ENNReal.ofReal (entM J h n q a b t M)) t)
  have hmeas : ∀ M, Measurable (G M) := fun M =>
    InfVolAux.measurable_comp_restrictZ M (fun q' => ⨆ t : ℝ, ENNReal.ofReal
      (vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
        (rdmF_isHermitian _ _)))
  calc ∫⁻ q, ⨆ t : ℝ, ENNReal.ofReal (entZ J h n q a b t) ∂muZ
      ≤ ∫⁻ q, liminf (fun M => G M q) atTop ∂muZ := lintegral_mono hpt
    _ ≤ liminf (fun M => ∫⁻ q, G M q ∂muZ) atTop := lintegral_liminf_le hmeas
    _ ≤ ENNReal.ofReal (4 * Real.log 2 * thetaF F) := by
        refine liminf_le_of_frequently_le' (Eventually.frequently ?_)
        filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
        exact lintegral_iSup_entM_le J h hEC n a b M ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩

/-- **Half-infinite regions** (paper: Theorem `thm:infvol`(d)): for `A_M = [-M, m]`,
`E_q S(ρ^q_{A_M}(t)) ≤ 2 ln 2 Θ(F)` for every `M ≥ |m|` and every `t`. -/
theorem areaLaw_halfInfinite (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ)
    (M : ℕ) (hm : |m| ≤ M) :
    avgQ (fun q' : Cfg (2 * M + 1) =>
      vonNeumann (rdm (m + M + 1).toNat (by have h1 := (abs_le.mp hm).1; have h2 := (abs_le.mp hm).2; omega)
        (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
        (rdm_isHermitian _ _ _)) ≤ 2 * Real.log 2 * thetaF F := by
  have hℓ : (m + M + 1).toNat ≤ 2 * M + 1 := by
    have h1 := (abs_le.mp hm).1
    have h2 := (abs_le.mp hm).2
    omega
  refine le_trans ?_ (lgt_areaLaw_sharp J h hEC (m + M + 1).toNat hℓ (occM M n))
  exact AreaLawAux.avgQ_mono' fun q' =>
    le_ciSup (InfVolAreaAux.bddAbove_vonNeumann_rdm J h q' (occM M n) _ hℓ) t

end DFL
