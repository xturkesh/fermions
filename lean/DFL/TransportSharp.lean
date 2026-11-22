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
  sorry

/-- `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (paper: eq. `eq:main-transport`). -/
theorem avg_sup_transport_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)|) ≤ sigma1 F := by
  sorry

/-- **Bounded transport with the constant of the paper** (paper: eq. `eq:main-transport`):
`|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)`. -/
theorem lgt_transport_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ sigma1 F := by
  sorry

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
