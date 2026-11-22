import DFL.DualityFull
import DFL.SectorNorm
import DFL.AreaLawSharp

/-!
# Bounded transport and number entropy with the constants of the paper
(paper: Theorems `thm:transport`(a)–(d), `thm:numberentropy`)

Constants: `Σ₁(F) = ∑_{d ≥ 1} d F(d)` (`sigma1`), `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}` (`sigma1'`),
`Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)` (`thetaF`), `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)` (`momentF`).

* (a) **Propagator**: `{c̃_x(t), c̃_y†} = ∑_q U_q(t)_{xy} P_q` (`dressed_anticomm_t`), so
  `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|` (`norm_sector_anticomm`),
  `E_q sup_t ‖P_q {c̃_x(t), c̃_y†} P_q‖ ≤ F(|x-y|)` (`avg_sup_anticomm_le`) and
  `⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)` (`expect_abs_anticomm_le`).
* (b) **One particle**: `∑_x (|x-y|+1)^p ⟨n_x⟩_t ≤ E_q sup_t ∑_x (|x-y|+1)^p |U_q(t)_{xy}|² ≤ 2 M_p(F)`
  for `0 ≤ p ≤ 2` (`one_particle_moments`).
* (c) **Transport across a cut**: `|⟨N_{≤m}⟩_{q,t} - s_m| ≤ ∑_{x ≤ m < z} |U_q(t)_{xz}|²`
  (`sector_transport_le`), `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (`avg_sup_transport_le`),
  `|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)` (`lgt_transport_sharp`) and
  `⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²` (`lgt_transport_sq`).
* (d) **Dressed correlations**: `|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)` (`dressed_corr_le`).
* **Number entropy**: `S_N(t) ≤ ½ ln(2πe(Θ(F) + Σ'₁(F)² + 1/12))` (`lgt_numberEntropy_sharp`) and, in a
  single sector, `E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))` (`avg_sup_sector_numberEntropy_le`).

Sites are 0-indexed: `leftOf m = {0, …, m}` is the paper's `[1, m+1]`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Σ₁(F) = ∑_{d ≥ 1} d F(d)`. -/
def sigma1 (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * F d

/-- `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}`. -/
def sigma1' (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * √(F d)

/-- `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)`. -/
def momentF (p : ℝ) (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d

/-- The dressed fermion at time `t`, `c̃_x(t) = e^{itH} c̃_x e^{-itH}`. -/
def dressedT (J h : ℝ) (x : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (H J h) t)ᴴ * cplx (dressed x) * evo (H J h) t

/-- The anticommutator `{c̃_x(t), c̃_y†}`. -/
def dressedAnticomm (J h : ℝ) (x y : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  dressedT J h x t * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * dressedT J h x t


namespace TransportSharpAux

/-! ### Auxiliary lemmas -/

section Cplx

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
/-- `cplx` is multiplicative. -/
lemma cplx_mul' (A B : Matrix ι ι ℝ) : cplx (A * B) = cplx A * cplx B := by
  ext i j; simp [cplx, Matrix.mul_apply]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` is additive. -/
lemma cplx_add' (A B : Matrix ι ι ℝ) : cplx (A + B) = cplx A + cplx B := by
  ext i j; simp [cplx]

omit [Fintype ι] in
/-- `cplx 1 = 1`. -/
lemma cplx_one' : cplx (1 : Matrix ι ι ℝ) = 1 := by
  ext i j; by_cases h : i = j <;> simp [cplx, Matrix.one_apply, h]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx 0 = 0`. -/
lemma cplx_zero' : cplx (0 : Matrix ι ι ℝ) = 0 := by
  ext i j; simp [cplx]

/-- `e^{-itA}` is a symmetric matrix when `A` is symmetric. -/
lemma evo_transpose (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) : (evo A t)ᵀ = evo A t := by
  unfold evo
  rw [← Matrix.exp_transpose, transpose_smul]
  congr 2
  unfold cplx
  rw [← Matrix.transpose_map, hA.eq]

end Cplx

/-- A real matrix commuting with every Gauss operator only connects configurations with the same
charges. -/
lemma charge_eq_of_comm_gauss (M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (hM : ∀ k, M * gauss k = gauss k * M) (a b : Cfg N × Cfg N) (hab : M a b ≠ 0) :
    chargeOf a.1 a.2 = chargeOf b.1 b.2 := by
  funext k
  have h := congrFun (congrFun (hM k) a) b
  rw [GaussLawAux.gauss_eq, mul_diagonal, diagonal_mul] at h
  by_contra hne
  have hs := sec_sgn_ne hne
  apply hab
  have h0 : M a b * (sgn (chargeOf b.1 b.2 k) - sgn (chargeOf a.1 a.2 k)) = 0 := by
    rw [mul_sub, h]; ring
  rcases mul_eq_zero.mp h0 with h1 | h1
  · exact h1
  · exact absurd (sub_eq_zero.mp h1).symm hs

/-- A real matrix commuting with every Gauss operator commutes with every sector projection. -/
lemma gaussProj_comm_of_comm_gauss (M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (hM : ∀ k, M * gauss k = gauss k * M) (q : Cfg N) :
    gaussProj q * M = M * gaussProj q := by
  unfold gaussProj
  exact GaussLawAux.diagonal_comm_of_support _ _ fun a b hab => by
    rw [charge_eq_of_comm_gauss M hM a b hab]

/-- `c̃_j†` commutes with every Gauss operator. -/
lemma dressed_transpose_comm_gauss (j k : Fin N) :
    (dressed j)ᵀ * gauss k = gauss k * (dressed j)ᵀ := by
  have h := congrArg Matrix.transpose (dressed_comm_gauss j k)
  rw [transpose_mul, transpose_mul, (gauss_isSymm k).eq] at h
  exact h.symm

/-- `P_q` commutes with `c̃_y†`. -/
lemma cplx_gaussProj_comm_dressed_transpose (q : Cfg N) (y : Fin N) :
    cplx (gaussProj q) * cplx (dressed y)ᵀ = cplx (dressed y)ᵀ * cplx (gaussProj q) := by
  rw [← cplx_mul', ← cplx_mul',
    gaussProj_comm_of_comm_gauss _ (dressed_transpose_comm_gauss y) q]

/-- The CAR `{c̃_z, c̃_y†} = δ_{zy}`, complexified. -/
lemma cplx_dressed_CAR (z y : Fin N) :
    cplx (dressed z) * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * cplx (dressed z) =
      if z = y then 1 else 0 := by
  rw [← cplx_mul', ← cplx_mul', ← cplx_add', dressed_CAR]
  split_ifs
  · exact cplx_one'
  · exact cplx_zero'

/-- The Bernoulli average is monotone. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  unfold avgQ
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Averaged propagator bound from Hypothesis (EC). -/
lemma avgQ_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (avgQ_mono fun q => norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans
    (hEC.bound N x y)

/-- Uniform-in-time averaged propagator bound from Hypothesis (EC). -/
lemma avgQ_iSup_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (avgQ_mono fun q => ciSup_le fun t =>
    norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans (hEC.bound N x y)


/-- There are `2^N` charge configurations. -/
lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Average of a constant. -/
lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  simp only [avgQ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
  push_cast
  field_simp

/-- Additivity of the average. -/
lemma avgQ_add (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  simp only [avgQ, Finset.sum_add_distrib]
  ring

/-- The average commutes with finite sums. -/
lemma avgQ_sum {κ : Type*} (s : Finset κ) (f : κ → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  simp only [avgQ]
  rw [Finset.sum_comm, Finset.mul_sum]

/-- Homogeneity of the average. -/
lemma avgQ_const_mul (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  simp only [avgQ, ← Finset.mul_sum]
  ring

/-- Average of a function minus a constant. -/
lemma avgQ_sub_const (f : Cfg N → ℝ) (c : ℝ) : avgQ (fun q => f q - c) = avgQ f - c := by
  have : avgQ (fun q => f q - c) = avgQ (fun q => f q + (-c)) := by
    simp only [sub_eq_add_neg]
  rw [this, avgQ_add, avgQ_const]
  ring

/-- `|avg f| ≤ avg |f|`. -/
lemma abs_avgQ_le (f : Cfg N → ℝ) : |avgQ f| ≤ avgQ (fun q => |f q|) := by
  simp only [avgQ]
  rw [abs_mul, abs_of_pos (by positivity)]
  exact mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) (by positivity)

/-- Real part of a complex Bernoulli average of real numbers. -/
lemma re_avg (f : Cfg N → ℝ) :
    ((1 / 2 ^ N : ℂ) * ∑ q : Cfg N, ((f q : ℝ) : ℂ)).re = avgQ f := by
  simp only [avgQ]
  rw [← Complex.ofReal_sum]
  have : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_mul, Complex.ofReal_re]

/-- `‖U_{xz}‖² ≤ ‖U_{xz}‖ ≤ Q(x, z)` for the propagator of `h_q`. -/
lemma norm_evo_sq_le_corrQ (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    ‖evo (hq J h q) t x z‖ ^ 2 ≤ corrQ (hq J h q) (hq_isHermitian J h q) x z := by
  have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x z
  have h0 := norm_nonneg (evo (hq J h q) t x z)
  have h2 := norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x z
  nlinarith


section Unitary

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Expectation of a complexified real diagonal matrix: `⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
lemma expect_cplx_diagonal (d : ι → ℝ) (v : ι → ℂ) :
    expect (cplx (diagonal d)) v = ((∑ i, d i * Complex.normSq (v i) : ℝ) : ℂ) := by
  have hc : cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
    unfold cplx; exact Matrix.diagonal_map (by simp)
  rw [expect, hc]
  simp only [dotProduct, mulVec_diagonal, Pi.star_apply, Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  simp only [Complex.star_def]
  ring

omit [DecidableEq ι] in
/-- `∑ |v_i|² = star v ⬝ᵥ v`. -/
lemma sum_normSq_eq_dotProduct (w : ι → ℂ) :
    ((∑ i, Complex.normSq (w i) : ℝ) : ℂ) = star w ⬝ᵥ w := by
  simp only [Complex.ofReal_sum, dotProduct, Pi.star_apply]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  rfl

/-- An isometry preserves `∑ |v_i|²`. -/
lemma sum_normSq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h1 : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, hU, vecMul_one]
  exact_mod_cast (sum_normSq_eq_dotProduct (U *ᵥ v)).trans
    (h1.trans (sum_normSq_eq_dotProduct v).symm)

/-- Rows of a unitary matrix are normalized. -/
lemma sum_normSq_row (U : Matrix ι ι ℂ) (hU : U * Uᴴ = 1) (x : ι) :
    ∑ z, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU x) x
  simp only [one_apply_eq] at h
  have : ((∑ z, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun z _ => ?_)
    rw [← Complex.mul_conj]
    rfl
  exact_mod_cast this

/-- Columns of a unitary matrix are normalized. -/
lemma sum_normSq_col (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (z : ι) :
    ∑ x, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU z) z
  simp only [one_apply_eq] at h
  have : ((∑ x, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [Complex.normSq_eq_conj_mul_self]
    rfl
  exact_mod_cast this

end Unitary

/-- `⟨n_x⟩` of an evolved Fock state, as a real sum over configurations. -/
lemma sum_occ_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x : Fin N) (t : ℝ) :
    ∑ n, occ (n x) * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  have := density A hA S x t
  rw [expect_cplx_diagonal] at this
  exact_mod_cast this

/-- `⟨N_R⟩` of an evolved Fock state: `∑_{x ∈ R} ∑_{z ∈ S} |U_{xz}|²`. -/
lemma sum_count_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  simp only [count, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun x _ => sum_occ_normSq A hA S x t)

/-- The evolved Fock state is normalized. -/
lemma sum_normSq_evo_fock (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) :
    ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) = 1 := by
  rw [sum_normSq_mulVec _ (evo_unitary _ (dGamma_isSymm A hA) t), fockVec,
    Finset.sum_eq_single (fockCfg S) (fun n _ hn => by rw [Pi.single_eq_of_ne hn, map_zero])
      (fun h => absurd (mem_univ _) h), Pi.single_eq_same, map_one]

/-- LGT expectation of a diagonal matter observable, as a Bernoulli average of sector sums. -/
lemma lgt_expect_diag_re (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    (expect (matter (diagonal d)) (lgtState J h ψ t)).re =
      avgQ (fun q => ∑ n, d n * Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ ψ) n)) := by
  rw [matter, sector_average_matter]
  simp_rw [expect_cplx_diagonal]
  exact re_avg _

/-- Sector expectation of `N_R` in an evolved Fock state. -/
lemma expect_numOp_re (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    (expect (cplx (numOp R)) (evo (dGamma A) t *ᵥ fockVec S)).re =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  rw [numOp, expect_cplx_diagonal, Complex.ofReal_re, sum_count_normSq A hA S R t]

/-- **Deviation of the mean particle number** for a symmetric unitary `U`:
`|∑_{x ∈ R, z ∈ S} |U_{xz}|² - |S ∩ R|| ≤ ∑_{x ∈ R, z ∉ R} |U_{xz}|²`. -/
lemma abs_mean_sub_le_symm (U : Matrix (Fin N) (Fin N) ℂ) (hU : Uᴴ * U = 1) (hT : Uᵀ = U)
    (S R : Finset (Fin N)) :
    |∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ)| ≤
      ∑ x ∈ R, ∑ z ∈ univ \ R, Complex.normSq (U x z) := by
  have ha : ∀ x z, 0 ≤ Complex.normSq (U x z) := fun x z => Complex.normSq_nonneg _
  have hsym : ∀ x z, U x z = U z x := fun x z => by
    have := congrFun (congrFun hT z) x
    rwa [transpose_apply] at this
  have hcol : ∀ z,
      ∑ x ∈ R, Complex.normSq (U x z) + ∑ x ∈ univ \ R, Complex.normSq (U x z) = 1 := by
    intro z
    rw [add_comm, Finset.sum_sdiff (subset_univ R)]
    exact sum_normSq_col U hU z
  set D := ∑ x ∈ R, ∑ z ∈ univ \ R, Complex.normSq (U x z) with hD
  set a := ∑ z ∈ S \ R, ∑ x ∈ R, Complex.normSq (U x z) with ha_def
  set b := ∑ z ∈ S ∩ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) with hb_def
  have hab : ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ) = a - b := by
    rw [Finset.sum_comm, ← Finset.sum_inter_add_sum_sdiff S R]
    have h1 : ∑ z ∈ S ∩ R, ∑ x ∈ R, Complex.normSq (U x z) = ((S ∩ R).card : ℝ) - b := by
      rw [show ((S ∩ R).card : ℝ) = ∑ z ∈ S ∩ R, (1 : ℝ) by simp, hb_def,
        ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl fun z _ => by linarith [hcol z]
    rw [h1]
    ring
  have ha0 : 0 ≤ a := Finset.sum_nonneg fun z _ => Finset.sum_nonneg fun x _ => ha x z
  have hb0 : 0 ≤ b := Finset.sum_nonneg fun z _ => Finset.sum_nonneg fun x _ => ha x z
  have haD : a ≤ D := by
    calc a ≤ ∑ z ∈ univ \ R, ∑ x ∈ R, Complex.normSq (U x z) :=
          Finset.sum_le_sum_of_subset_of_nonneg (Finset.sdiff_subset_sdiff (subset_univ S) le_rfl)
            (fun z _ _ => Finset.sum_nonneg fun x _ => ha x z)
      _ = D := Finset.sum_comm
  have hbD : b ≤ D := by
    calc b ≤ ∑ z ∈ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) :=
          Finset.sum_le_sum_of_subset_of_nonneg Finset.inter_subset_right
            (fun z _ _ => Finset.sum_nonneg fun x _ => ha x z)
      _ = D := Finset.sum_congr rfl fun z _ => Finset.sum_congr rfl fun x _ => by rw [hsym]
  rw [hab, abs_sub_le_iff]
  constructor <;> linarith

/-! #### Counting lemmas -/

/-- Bounding `∑_{a ∈ s} g(φ a)` by the cardinalities of the fibres of `φ`. -/
lemma sum_comp_le_of_card_fiber {α : Type*} (s : Finset α) (φ : α → ℕ) (M : ℕ)
    (hφ : ∀ a ∈ s, φ a < M) (c : ℕ → ℝ)
    (hc : ∀ d, ((s.filter (fun a => φ a = d)).card : ℝ) ≤ c d)
    (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) :
    ∑ a ∈ s, g (φ a) ≤ ∑ d ∈ range M, c d * g d := by
  rw [← Finset.sum_fiberwise_of_maps_to' (t := range M) (fun a ha => mem_range.2 (hφ a ha))]
  refine Finset.sum_le_sum (fun d _ => ?_)
  rw [Finset.sum_const, nsmul_eq_mul]
  exact mul_le_mul_of_nonneg_right (hc d) (hg d)

/-- At most `d` pairs `(x, y)` with `x ≤ m < y` are at distance `d`. -/
lemma card_cross_fiber_le (m d : ℕ) :
    (((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card : ℝ) ≤ d := by
  have h : ((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card ≤ (range d).card := by
    apply Finset.card_le_card_of_injOn (fun p : Fin N × Fin N => m - p.1.val)
    · intro p hp
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp
      simp only [coe_range, Set.mem_Iio]
      unfold Nat.dist at hp
      omega
    · intro p hp p' hp' hpp'
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp hp'
      simp only at hpp'
      unfold Nat.dist at hp hp'
      ext <;> omega
  rw [card_range] at h
  exact_mod_cast h

/-- **Crossing pairs**: `∑_{x ≤ m < y} g(|x - y|) ≤ ∑_d d g(d)`. -/
lemma sum_cross_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d)
    (hs : Summable (fun d : ℕ => (d : ℝ) * g d)) (m : ℕ) :
    ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, g (Nat.dist x.val y.val) ≤
      ∑' d : ℕ, (d : ℝ) * g d := by
  rw [← Finset.sum_product' (f := fun x y : Fin N => g (Nat.dist x.val y.val))]
  calc _ ≤ ∑ d ∈ range N, (d : ℝ) * g d :=
        sum_comp_le_of_card_fiber _ (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val) N
          (fun p _ => by have := p.1.isLt; have := p.2.isLt; unfold Nat.dist; omega)
          (fun d => d) (card_cross_fiber_le m) g hg
    _ ≤ _ := hs.sum_le_tsum _ (fun d _ => mul_nonneg (Nat.cast_nonneg d) (hg d))

/-! #### Consequences of Hypothesis (EC) -/

/-- Summability of `d ↦ d F(d)`. -/
lemma ec_summable_mul {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * F d) :=
  hEC.summable.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (hEC.nonneg d))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (hEC.nonneg d))

/-- Summability of `d ↦ d √F(d)`. -/
lemma ec_summable_mul_sqrt {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * √(F d)) :=
  hEC.summable_sqrt.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (Real.sqrt_nonneg _))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (Real.sqrt_nonneg _))

/-- `E ∑_{x ≤ m < z} Q(x, z) ≤ Σ₁(F)`. -/
lemma avgQ_cross_corrQ_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ sigma1 F := by
  calc _ = ∑ x ∈ leftOf (N := N) m, ∑ z ∈ univ \ leftOf (N := N) m,
        avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x z) := by
        rw [avgQ_sum]
        exact Finset.sum_congr rfl fun x _ => avgQ_sum _ _
    _ ≤ ∑ x ∈ leftOf (N := N) m, ∑ z ∈ univ \ leftOf (N := N) m, F (Nat.dist x.val z.val) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun z _ => hEC.bound N x z
    _ ≤ sigma1 F := sum_cross_le F hEC.nonneg (ec_summable_mul hEC) m

end TransportSharpAux

open TransportSharpAux

/-! ### (a) The dressed propagator -/

/-- **The dressed anticommutator is the one-body propagator in each sector** (paper: eq. `eq:main-anticomm`). -/
theorem dressed_anticomm_t (J h : ℝ) (x y : Fin N) (t : ℝ) :
    dressedAnticomm J h x y t = ∑ q : Cfg N, evo (hq J h q) t x y • cplx (gaussProj q) := by
  unfold dressedAnticomm dressedT
  rw [heisenberg_dressed]
  calc (∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z • cplx (dressed z)) *
        cplx (dressed y)ᵀ + cplx (dressed y)ᵀ *
        (∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z • cplx (dressed z))
      = ∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z •
          (cplx (dressed z) * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * cplx (dressed z)) := by
        rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun q _ => ?_
        rw [Matrix.mul_assoc, ← Matrix.mul_assoc (cplx (dressed y)ᵀ) (cplx (gaussProj q)),
          ← cplx_gaussProj_comm_dressed_transpose q y, Matrix.mul_assoc, ← Matrix.mul_add]
        congr 1
        rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun z _ => ?_
        rw [Matrix.smul_mul, Matrix.mul_smul, ← smul_add]
    _ = ∑ q : Cfg N, evo (hq J h q) t x y • cplx (gaussProj q) := by
        refine Finset.sum_congr rfl fun q _ => ?_
        simp only [cplx_dressed_CAR, smul_ite, smul_zero, Finset.sum_ite_eq', Finset.mem_univ,
          ite_true]
        rw [Matrix.mul_smul, Matrix.mul_one]

/-- `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|`. -/
theorem norm_sector_anticomm (J h : ℝ) (x y : Fin N) (t : ℝ) (q : Cfg N) :
    ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖ = ‖evo (hq J h q) t x y‖ := by
  have hdiag : cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q) =
      diagonal (fun p : Cfg N × Cfg N =>
        if chargeOf p.1 p.2 = q then evo (hq J h q) t x y else 0) := by
    rw [dressed_anticomm_t, SectorNormAux.sum_smul_gaussProj, SectorNormAux.cplx_gaussProj,
      diagonal_mul_diagonal, diagonal_mul_diagonal]
    congr 1
    funext p
    split_ifs with hp
    · rw [hp]; ring
    · ring
  rw [hdiag, l2_opNorm_diagonal]
  apply le_antisymm
  · refine (pi_norm_le_iff_of_nonneg (norm_nonneg _)).mpr fun p => ?_
    split_ifs
    · exact le_rfl
    · simp
  · set v : Cfg N × Cfg N → ℂ := fun p =>
      if chargeOf p.1 p.2 = q then evo (hq J h q) t x y else 0 with hv
    have h1 : v ((0 : Cfg N), linkOf 0 q) = evo (hq J h q) t x y := by
      simp only [hv]
      exact ite_eq_left (chargeOf_linkOf 0 q)
    calc ‖evo (hq J h q) t x y‖ = ‖v ((0 : Cfg N), linkOf 0 q)‖ := by rw [h1]
      _ ≤ ‖v‖ := norm_le_pi_norm v _

/-- **Localization of the dressed propagator** (paper: eq. `eq:main-prop`). -/
theorem avg_sup_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N =>
      ⨆ t : ℝ, ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖) ≤
      F (Nat.dist x.val y.val) := by
  simp_rw [norm_sector_anticomm]
  exact avgQ_iSup_norm_evo_le J h hEC x y

/-- **`⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)`** (paper: Theorem `thm:transport`(a)): the operator
`∑_q |U_q(t)_{xy}| P_q` is positive semidefinite with square `{c̃_x(t), c̃_y†}^* {c̃_x(t), c̃_y†}` (so it is the
absolute value), and its expectation in `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` is at most `F(|x-y|)`. -/
theorem expect_abs_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ψ : Cfg N → ℂ)
    (hψ : star ψ ⬝ᵥ ψ = 1) (x y : Fin N) (t : ℝ) :
    (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)).PosSemidef ∧
      (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) *
          (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) =
        (dressedAnticomm J h x y t)ᴴ * dressedAnticomm J h x y t ∧
      (expect (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) (initZ ψ)).re ≤
        F (Nat.dist x.val y.val) := by
  obtain ⟨h1, h2⟩ := sectorDiag_abs (fun q : Cfg N => evo (hq J h q) t x y)
  refine ⟨h1, ?_, ?_⟩
  · rw [h2, dressed_anticomm_t]
  · rw [expect_sectorDiag_abs ψ hψ (fun q : Cfg N => evo (hq J h q) t x y), Complex.ofReal_re]
    exact avgQ_norm_evo_le J h hEC x y t

/-! ### (b) One particle -/

/-- **Moments of a single particle** (paper: Theorem `thm:transport`(b)), `0 ≤ p ≤ 2`. -/
theorem one_particle_moments (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (y : Fin N) (p : ℝ)
    (hp0 : 0 ≤ p) (hp2 : p ≤ 2) (t : ℝ) :
    ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec {y}) t)).re ≤
      avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ∧
    avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ≤
      2 * momentF p F := by
  have hw0 : ∀ x : Fin N, 0 ≤ ((Nat.dist x.val y.val : ℝ) + 1) ^ p := fun x =>
    Real.rpow_nonneg (by positivity) p
  constructor
  · have e : ∀ x : Fin N,
        (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec {y}) t)).re =
          avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖ ^ 2) := by
      intro x
      rw [lgt_density]
      simp only [Finset.sum_singleton, Complex.normSq_eq_norm_sq]
    simp_rw [e, ← avgQ_const_mul]
    rw [← avgQ_sum]
    refine avgQ_mono fun q => ?_
    refine le_ciSup (f := fun s : ℝ => ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
      ‖evo (hq J h q) s x y‖ ^ 2) ⟨∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p, ?_⟩ t
    rintro _ ⟨s, rfl⟩
    refine Finset.sum_le_sum fun x _ => mul_le_of_le_one_right (hw0 x) ?_
    exact pow_le_one₀ (norm_nonneg _) (norm_evo_le_one _ (hq_isSymm J h q) s x y)
  · have hpt : ∀ q : Cfg N, (⨆ s : ℝ, ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        ‖evo (hq J h q) s x y‖ ^ 2) ≤ ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        corrQ (hq J h q) (hq_isHermitian J h q) x y := fun q =>
      ciSup_le fun s => Finset.sum_le_sum fun x _ =>
        mul_le_mul_of_nonneg_left (norm_evo_sq_le_corrQ J h q s x y) (hw0 x)
    have hg0 : ∀ d : ℕ, 0 ≤ ((d : ℝ) + 1) ^ p * F d := fun d =>
      mul_nonneg (Real.rpow_nonneg (by positivity) p) (hEC.nonneg d)
    have hgs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d) := by
      refine hEC.summable.of_nonneg_of_le hg0 fun d => ?_
      refine mul_le_mul_of_nonneg_right ?_ (hEC.nonneg d)
      have h1 : (1 : ℝ) ≤ (d : ℝ) + 1 := by linarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]
      calc ((d : ℝ) + 1) ^ p ≤ ((d : ℝ) + 1) ^ (2 : ℝ) := Real.rpow_le_rpow_of_exponent_le h1 hp2
        _ = ((d : ℝ) + 1) ^ 2 := Real.rpow_two _
    calc avgQ (fun q : Cfg N => ⨆ s : ℝ, ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          ‖evo (hq J h q) s x y‖ ^ 2)
        ≤ avgQ (fun q : Cfg N => ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          corrQ (hq J h q) (hq_isHermitian J h q) x y) := avgQ_mono hpt
      _ = ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) := by
          rw [avgQ_sum]
          exact Finset.sum_congr rfl fun x _ => avgQ_const_mul _ _
      _ ≤ ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * F (Nat.dist x.val y.val) :=
          Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (hEC.bound N x y) (hw0 x)
      _ ≤ 2 * momentF p F :=
          AreaLawAux.sum_dist_le (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d) hg0 hgs y

/-! ### (c) Transport across a cut -/

/-- **Transport in one sector** (paper: Theorem `thm:transport`(c)). -/
theorem sector_transport_le (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (m : ℕ) (t : ℝ) :
    |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m, ‖evo (hq J h q) t x z‖ ^ 2 := by
  rw [expect_numOp_re _ (hq_isSymm J h q)]
  refine (abs_mean_sub_le_symm _ (evo_unitary _ (hq_isSymm J h q) t)
    (evo_transpose _ (hq_isSymm J h q) t) S _).trans (le_of_eq ?_)
  simp only [Complex.normSq_eq_norm_sq]

/-- `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (paper: eq. `eq:main-transport`). -/
theorem avg_sup_transport_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)|) ≤ sigma1 F := by
  refine (avgQ_mono fun q => ciSup_le fun t => ?_).trans (avgQ_cross_corrQ_le hEC m)
  exact (sector_transport_le J h q S m t).trans (Finset.sum_le_sum fun x _ =>
    Finset.sum_le_sum fun z _ => norm_evo_sq_le_corrQ J h q t x z)

/-- **Bounded transport with the constant of the paper** (paper: eq. `eq:main-transport`):
`|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)`. -/
theorem lgt_transport_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ sigma1 F := by
  rw [numOp, lgt_expect_diag_re, ← avgQ_sub_const]
  refine (abs_avgQ_le _).trans ((avgQ_mono fun q => ?_).trans (avgQ_cross_corrQ_le hEC m))
  have h1 := sector_transport_le J h q S m t
  rw [numOp, expect_cplx_diagonal, Complex.ofReal_re] at h1
  exact h1.trans (Finset.sum_le_sum fun x _ =>
    Finset.sum_le_sum fun z _ => norm_evo_sq_le_corrQ J h q t x z)

/-- **Second moment of the transported charge** (paper: eq. `eq:main-transport2`):
`⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²`. -/
theorem lgt_transport_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    (expect (matter ((numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1) *
        (numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1))) (lgtState J h (fockVec S) t)).re ≤
      thetaF F + sigma1' F ^ 2 := by
  sorry

/-! ### (d) Dressed correlations -/

/-- **Dressed correlations** (paper: Theorem `thm:transport`(d)): for a Fock state and `x ≠ y`,
`|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)`. -/
theorem dressed_corr_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (x y : Fin N)
    (hxy : x ≠ y) (t : ℝ) :
    ‖expect (cplx ((dressed x)ᵀ * dressed y)) (lgtState J h (fockVec S) t)‖ ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  sorry

/-! ### Number entropy -/

/-- **Bounded number entropy with the constants of the paper** (paper: Theorem `thm:numberentropy`):
`S_N(t) ≤ ½ ln(2πe(K₁ + K₂ + 1/12))`, `K₁ = Θ(F)`, `K₂ = Σ'₁(F)²`. -/
theorem lgt_numberEntropy_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + sigma1' F ^ 2 + 1 / 12)) := by
  sorry

/-- **Number entropy in a single sector** (paper: Theorem `thm:numberentropy`):
`E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))`. -/
theorem avg_sup_sector_numberEntropy_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N))
    (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      shannon N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + 1 / 12)) := by
  sorry

end DFL
