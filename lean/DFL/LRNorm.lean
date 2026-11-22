import DFL.LRFull
import DFL.Imbalance

/-!
# The operator norm sees the clean sector (paper: Proposition `prop:LR`(c))

* `‖[n_x(t), n_y]‖ = max_q |U_q(t)_{xy}| (1 - |U_q(t)_{xy}|²)^{1/2}` (`norm_comm_nn_eq`): in the sector `q` the
  commutator is `dΓ(Z)` with `Z = [|w⟩⟨w|, |δ_y⟩⟨δ_y|]`, `w = U_q(t)^* δ_x`, whose eigenvalues are `±iμ, 0`,
  `μ = |U_{xy}| (1 - |U_{xy}|²)^{1/2}`;
* in the clean sector `q = (1, …, 1)` the eigenvectors of `h_q = -JΔ + 2h` are the sine modes, and
  `sup_t |U(t)_{1N}|² ≥ ∑_k φ_k(1)⁴ = 3/(2(N+1))` for `N ≥ 2`, whence
  `sup_t ‖{c̃_1(t), c̃_N†}‖ ≥ (3/(2(N+1)))^{1/2}` and `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}`
  (`sup_norm_dressedAnticomm_ge`, `sup_norm_comm_nn_ge`), while the sector averages are at most `F(N-1)`.

Sites are 0-indexed: the paper's sites `1` and `N` are `0` and `N - 1`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The commutator `[n_x(t), n_y]`. -/
def numComm (J h : ℝ) (x y : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  heisT J h t (cplx (numOp {x} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) *
      cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) -
    cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
      heisT J h t (cplx (numOp {x} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))

namespace LRNormAux

open LRFullAux

/-! ### Two orthogonal fermionic modes -/

section Alg

variable {𝒜 : Type*} [Ring 𝒜] [Algebra ℂ 𝒜]

omit [Algebra ℂ 𝒜] in
/-- Right-associated form of a product relation. -/
lemma mulX {x y z : 𝒜} (h : x * y = z) (X : 𝒜) : x * (y * X) = z * X := by
  rw [← mul_assoc, h]

/-- The canonical anticommutation relations of two orthogonal fermionic modes `e` (normalized) and
`r` (`{r, r†} = s`), as normal-ordering rules for the order `e† < r† < e < r`. -/
structure TwoModes (e ed r rd : 𝒜) (s : ℂ) : Prop where
  ee : e * e = 0
  dd : ed * ed = 0
  e_ed : e * ed = 1 - ed * e
  rr : r * r = 0
  rdrd : rd * rd = 0
  r_rd : r * rd = s • 1 - rd * r
  r_e : r * e = -(e * r)
  e_rd : e * rd = -(rd * e)
  r_ed : r * ed = -(ed * r)
  rd_ed : rd * ed = -(ed * rd)

/-- Normal ordering of polynomials in two orthogonal fermionic modes. -/
macro "car_norm" H:term : tactic => `(tactic| (simp only [mul_add, add_mul, mul_sub, sub_mul,
    mul_neg, neg_mul, smul_mul_assoc, mul_smul_comm, mul_assoc, one_mul, mul_one, zero_mul,
    mul_zero, smul_add, smul_sub, smul_neg, neg_neg, smul_zero, smul_smul, neg_zero, add_zero,
    zero_add, sub_zero, zero_sub, ($H).ee, ($H).dd, ($H).e_ed, ($H).rr, ($H).rdrd, ($H).r_rd,
    ($H).r_e, ($H).e_rd, ($H).r_ed, ($H).rd_ed, mulX ($H).ee, mulX ($H).dd, mulX ($H).e_ed,
    mulX ($H).rr, mulX ($H).rdrd, mulX ($H).r_rd, mulX ($H).r_e, mulX ($H).e_rd, mulX ($H).r_ed,
    mulX ($H).rd_ed]))

variable {e ed r rd : 𝒜} {s : ℂ}

set_option linter.unusedSimpArgs false

/-- `[A†A, e†e] = a r†e - ā e†r` for `A = a e + r`. -/
lemma TwoModes.comm_number (H : TwoModes e ed r rd s) (a : ℂ) :
    (star a • ed + rd) * (a • e + r) * (ed * e) - (ed * e) * ((star a • ed + rd) * (a • e + r)) =
      a • (rd * e) - star a • (ed * r) := by
  car_norm H
  module

/-- `K†K = |a|² D` for `K = a r†e - ā e†r`, with `D = r†r + s e†e + 2 e†r†er`. -/
lemma TwoModes.kdk (H : TwoModes e ed r rd s) (a : ℂ) :
    (star a • (ed * r) - a • (rd * e)) * (a • (rd * e) - star a • (ed * r)) =
      (star a * a) • (rd * r + s • (ed * e) + (2 : ℂ) • (ed * (rd * (e * r)))) := by
  car_norm H
  module

/-- `D² = s D`. -/
lemma TwoModes.dsq (H : TwoModes e ed r rd s) :
    (rd * r + s • (ed * e) + (2 : ℂ) • (ed * (rd * (e * r)))) *
      (rd * r + s • (ed * e) + (2 : ℂ) • (ed * (rd * (e * r)))) =
      s • (rd * r + s • (ed * e) + (2 : ℂ) • (ed * (rd * (e * r)))) := by
  car_norm H
  module

/-- `D e† Ω = s e† Ω` when `r Ω = 0`. -/
lemma TwoModes.dvec (H : TwoModes e ed r rd s) (Q : 𝒜) (hrQ : r * Q = 0) :
    (rd * r + s • (ed * e) + (2 : ℂ) • (ed * (rd * (e * r)))) * (ed * Q) = s • (ed * Q) := by
  car_norm H
  simp only [hrQ, mul_zero, smul_zero, add_zero, zero_add, sub_zero, neg_zero]

end Alg

/-! ### Norm of a sector block -/

/-- If the orthogonal projection `P` commutes with `K`, `(K†K)² = λ K†K`, and `K†K T = λ T` for some
`T ≠ 0` with `P T = T`, then `‖P K P‖ = √λ`. -/
lemma norm_block_eq {ι : Type*} [Fintype ι] [DecidableEq ι] (P K T : Matrix ι ι ℂ) (l : ℝ)
    (hl : 0 ≤ l) (hPP : P * P = P) (hPh : Pᴴ = P) (hPK : P * K = K * P)
    (hX : (Kᴴ * K) * (Kᴴ * K) = (l : ℂ) • (Kᴴ * K)) (hT : T ≠ 0) (hPT : P * T = T)
    (hXT : (Kᴴ * K) * T = (l : ℂ) • T) : ‖P * K * P‖ = √l := by
  have hPKh : P * Kᴴ = Kᴴ * P := by
    have h1 := congrArg conjTranspose hPK
    rw [conjTranspose_mul, conjTranspose_mul, hPh] at h1
    exact h1.symm
  have hPX : P * (Kᴴ * K) = (Kᴴ * K) * P := by
    rw [← Matrix.mul_assoc, hPKh, Matrix.mul_assoc, hPK, Matrix.mul_assoc]
  have hY : P * K * P = P * K := by
    rw [Matrix.mul_assoc, ← hPK, ← Matrix.mul_assoc, hPP]
  have hYC : (P * K)ᴴ * (P * K) = P * (Kᴴ * K) := by
    rw [conjTranspose_mul, hPh, Matrix.mul_assoc, ← Matrix.mul_assoc P P K, hPP, hPK,
      ← Matrix.mul_assoc, ← hPX]
  have hCh : (P * (Kᴴ * K))ᴴ = P * (Kᴴ * K) := by
    rw [conjTranspose_mul, conjTranspose_mul, conjTranspose_conjTranspose, hPh, ← hPX]
  have hCC : P * (Kᴴ * K) * (P * (Kᴴ * K)) = (l : ℂ) • (P * (Kᴴ * K)) := by
    calc P * (Kᴴ * K) * (P * (Kᴴ * K)) = P * ((Kᴴ * K) * P * (Kᴴ * K)) := by
          simp only [Matrix.mul_assoc]
      _ = P * (P * ((Kᴴ * K) * (Kᴴ * K))) := by rw [← hPX]; simp only [Matrix.mul_assoc]
      _ = (l : ℂ) • (P * (Kᴴ * K)) := by
          rw [← Matrix.mul_assoc, hPP, hX, Matrix.mul_smul]
  have hCT : P * (Kᴴ * K) * T = (l : ℂ) • T := by
    rw [Matrix.mul_assoc, hXT, Matrix.mul_smul, hPT]
  set C := P * (Kᴴ * K) with hC
  have hnormY : ‖P * K‖ * ‖P * K‖ = ‖C‖ := by
    rw [← hYC, l2_opNorm_conjTranspose_mul_self]
  have hnormC : ‖C‖ * ‖C‖ = l * ‖C‖ := by
    rw [← l2_opNorm_conjTranspose_mul_self, hCh, hCC, norm_smul, Complex.norm_real,
      Real.norm_of_nonneg hl]
  have hlC : l ≤ ‖C‖ := by
    have h1 : l * ‖T‖ ≤ ‖C‖ * ‖T‖ := by
      calc l * ‖T‖ = ‖C * T‖ := by
            rw [hCT, norm_smul, Complex.norm_real, Real.norm_of_nonneg hl]
        _ ≤ ‖C‖ * ‖T‖ := norm_mul_le _ _
    exact le_of_mul_le_mul_right h1 (norm_pos_iff.mpr hT)
  have hCl : ‖C‖ = l := by
    have h2 : ‖C‖ * (‖C‖ - l) = 0 := by rw [mul_sub, hnormC]; ring
    rcases mul_eq_zero.mp h2 with h | h
    · rw [h] at hlC ⊢
      linarith
    · linarith
  rw [hY, ← Real.sqrt_mul_self (norm_nonneg (P * K)), hnormY, hCl]

/-! ### Mode operators -/

/-- The mode operator `∑_z f_z c̃_z`. -/
def modeOp (f : Fin N → ℂ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  ∑ z, f z • cplx (dressed z)

/-- `(∑_z f_z c̃_z)† = ∑_z f̄_z c̃_z†`. -/
lemma modeOp_conjTranspose (f : Fin N → ℂ) :
    (modeOp f)ᴴ = ∑ z, star (f z) • cplx (dressed z)ᵀ := by
  rw [modeOp, conjTranspose_sum]
  exact Finset.sum_congr rfl fun z _ => by rw [conjTranspose_smul, cplx_conjTranspose']

/-- `modeOp` is additive. -/
lemma modeOp_add (f g : Fin N → ℂ) : modeOp (f + g) = modeOp f + modeOp g := by
  simp only [modeOp, Pi.add_apply, add_smul, Finset.sum_add_distrib]

/-- `modeOp` is homogeneous. -/
lemma modeOp_smul (c : ℂ) (f : Fin N → ℂ) : modeOp (c • f) = c • modeOp f := by
  simp only [modeOp, Pi.smul_apply, smul_eq_mul, mul_smul, Finset.smul_sum]

/-- `modeOp δ_y = c̃_y`. -/
lemma modeOp_single (y : Fin N) : modeOp (Pi.single y 1) = cplx (dressed y) := by
  rw [modeOp, Finset.sum_eq_single y]
  · rw [Pi.single_eq_same, one_smul]
  · intro z _ hz
    rw [Pi.single_eq_of_ne hz, zero_smul]
  · intro h
    exact absurd (Finset.mem_univ y) h

/-- `{c̃_z, c̃_{z'}} = 0` (complexified). -/
lemma cplx_dressed_anticomm (z z' : Fin N) :
    cplx (dressed z) * cplx (dressed z') + cplx (dressed z') * cplx (dressed z) = 0 := by
  rw [← TransportSharpAux.cplx_mul', ← TransportSharpAux.cplx_mul',
    ← TransportSharpAux.cplx_add', dressed_anticomm, TransportSharpAux.cplx_zero']

/-- `{b(f), b(g)} = 0`. -/
lemma modeOp_anticomm (f g : Fin N → ℂ) : modeOp f * modeOp g + modeOp g * modeOp f = 0 := by
  unfold modeOp
  rw [Finset.sum_mul_sum, Finset.sum_mul_sum, Finset.sum_comm (s := univ) (t := univ)
    (f := fun z z' => g z • cplx (dressed z) * (f z' • cplx (dressed z'))),
    ← Finset.sum_add_distrib]
  refine Finset.sum_eq_zero fun z _ => ?_
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_eq_zero fun z' _ => ?_
  rw [smul_mul_smul_comm, smul_mul_smul_comm, mul_comm (g z') (f z), ← smul_add,
    cplx_dressed_anticomm, smul_zero]

/-- `{b(f), b(g)†} = ⟨g, f⟩`. -/
lemma modeOp_CAR (f g : Fin N → ℂ) :
    modeOp f * (modeOp g)ᴴ + (modeOp g)ᴴ * modeOp f = (∑ z, f z * star (g z)) • 1 := by
  rw [modeOp_conjTranspose]
  unfold modeOp
  rw [Finset.sum_mul_sum, Finset.sum_mul_sum, Finset.sum_comm (s := univ) (t := univ)
    (f := fun z z' => star (g z) • cplx (dressed z)ᵀ * (f z' • cplx (dressed z'))),
    ← Finset.sum_add_distrib, Finset.sum_smul]
  refine Finset.sum_congr rfl fun z _ => ?_
  rw [← Finset.sum_add_distrib]
  have e1 : ∀ z', f z • cplx (dressed z) * (star (g z') • cplx (dressed z')ᵀ) +
      star (g z') • cplx (dressed z')ᵀ * (f z • cplx (dressed z)) =
      (f z * star (g z')) • (if z = z' then (1 : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) else 0) := by
    intro z'
    rw [smul_mul_smul_comm, smul_mul_smul_comm, mul_comm (star (g z')) (f z), ← smul_add,
      TransportSharpAux.cplx_dressed_CAR]
  rw [Finset.sum_congr rfl fun z' _ => e1 z']
  simp only [smul_ite, smul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true]

/-- Two orthogonal mode operators, the first one normalized, satisfy `TwoModes`. -/
lemma twoModes_modeOp (f g : Fin N → ℂ) (hf : ∑ z, f z * star (f z) = 1)
    (hfg : ∑ z, f z * star (g z) = 0) :
    TwoModes (modeOp f) (modeOp f)ᴴ (modeOp g) (modeOp g)ᴴ (∑ z, g z * star (g z)) := by
  have hsq : ∀ M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ, M + M = 0 → M = 0 := fun M hM => by
    have h2 : (2 : ℂ) • M = 0 := by rw [two_smul]; exact hM
    exact (smul_eq_zero.mp h2).resolve_left two_ne_zero
  have hgf : ∑ z, g z * star (f z) = 0 := by
    have h1 := congrArg star hfg
    rw [star_sum, star_zero] at h1
    rw [← h1]
    exact Finset.sum_congr rfl fun z _ => by rw [star_mul, star_star, mul_comm]
  have hff := hsq _ (modeOp_anticomm f f)
  have hgg := hsq _ (modeOp_anticomm g g)
  have hfg' : modeOp g * modeOp f = -(modeOp f * modeOp g) :=
    eq_neg_of_add_eq_zero_right (modeOp_anticomm f g)
  refine ⟨hff, ?_, ?_, hgg, ?_, ?_, hfg', ?_, ?_, ?_⟩
  · rw [← conjTranspose_mul, hff, conjTranspose_zero]
  · have h1 := modeOp_CAR f f
    rw [hf, one_smul] at h1
    exact eq_sub_of_add_eq h1
  · rw [← conjTranspose_mul, hgg, conjTranspose_zero]
  · exact eq_sub_of_add_eq (modeOp_CAR g g)
  · have h1 := modeOp_CAR f g
    rw [hfg, zero_smul] at h1
    exact eq_neg_of_add_eq_zero_left h1
  · have h1 := modeOp_CAR g f
    rw [hgf, zero_smul] at h1
    exact eq_neg_of_add_eq_zero_left h1
  · rw [← conjTranspose_mul, ← conjTranspose_mul, hfg', conjTranspose_neg, neg_neg]

/-- `P_q` commutes with every mode operator. -/
lemma comm_gaussProj_modeOp (q : Cfg N) (f : Fin N → ℂ) :
    cplx (gaussProj q) * modeOp f = modeOp f * cplx (gaussProj q) :=
  comm_gaussProj_sum_dressedSharp q f false

/-- `P_q` commutes with the adjoint of every mode operator. -/
lemma comm_gaussProj_modeOp_conjTranspose (q : Cfg N) (f : Fin N → ℂ) :
    cplx (gaussProj q) * (modeOp f)ᴴ = (modeOp f)ᴴ * cplx (gaussProj q) := by
  rw [modeOp_conjTranspose]
  exact comm_gaussProj_sum_dressedSharp q (fun z => star (f z)) true

/-! ### The vacuum of a sector -/

/-- The configuration without fermions in the sector `q`. -/
def vac (q : Cfg N) : Cfg N × Cfg N := (0, linkOf 0 q)

/-- The projection onto the vacuum of the sector `q`. -/
def vacProj (q : Cfg N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  diagonal (fun p => if p = vac q then 1 else 0)

/-- `c̃_z Ω_q = 0`. -/
lemma cplx_dressed_mul_vacProj (z : Fin N) (q : Cfg N) : cplx (dressed z) * vacProj q = 0 := by
  ext a b
  rw [vacProj, mul_diagonal, Matrix.zero_apply]
  split_ifs with hb
  · subst hb
    obtain ⟨a1, a2⟩ := a
    have h01 : (0 : ZMod 2) ≠ 1 := by decide
    simp [cplx, dressed, vac, ann, h01]
  · rw [mul_zero]

/-- `b(f) Ω_q = 0`. -/
lemma modeOp_mul_vacProj (f : Fin N → ℂ) (q : Cfg N) : modeOp f * vacProj q = 0 := by
  rw [modeOp, Finset.sum_mul]
  exact Finset.sum_eq_zero fun z _ => by rw [smul_mul_assoc, cplx_dressed_mul_vacProj, smul_zero]

/-- The vacuum projection is idempotent. -/
lemma vacProj_mul_self (q : Cfg N) : vacProj q * vacProj q = vacProj q := by
  rw [vacProj, diagonal_mul_diagonal]
  congr 1
  funext p
  split_ifs <;> simp

/-- The vacuum projection is self-adjoint. -/
lemma vacProj_conjTranspose (q : Cfg N) : (vacProj q)ᴴ = vacProj q := by
  rw [vacProj, diagonal_conjTranspose]
  congr 1
  funext p
  simp only [Pi.star_apply]
  split_ifs <;> simp

/-- The vacuum of the sector `q` lies in the sector `q`. -/
lemma gaussProj_mul_vacProj (q : Cfg N) : cplx (gaussProj q) * vacProj q = vacProj q := by
  rw [SectorNormAux.cplx_gaussProj, vacProj, diagonal_mul_diagonal]
  congr 1
  funext p
  split_ifs with h1 h2 h2
  · rw [mul_one]
  · rw [mul_zero]
  · exact absurd (by rw [h2]; exact chargeOf_linkOf 0 q) h1
  · rw [mul_zero]

/-- The vacuum projection is nonzero. -/
lemma vacProj_ne_zero (q : Cfg N) : vacProj q ≠ 0 := by
  intro h
  have h1 := congrFun (congrFun h (vac q)) (vac q)
  simp [vacProj] at h1

/-- `c̃_y† Ω_q ≠ 0`. -/
lemma dressed_transpose_mul_vacProj_ne_zero (y : Fin N) (q : Cfg N) :
    (cplx (dressed y))ᴴ * vacProj q ≠ 0 := by
  intro h
  have hCAR : cplx (dressed y) * (cplx (dressed y))ᴴ = 1 - (cplx (dressed y))ᴴ * cplx (dressed y) := by
    have h2 := TransportSharpAux.cplx_dressed_CAR y y
    rw [ite_eq_left rfl, ← cplx_conjTranspose'] at h2
    exact eq_sub_of_add_eq h2
  have h1 : ((cplx (dressed y))ᴴ * vacProj q)ᴴ * ((cplx (dressed y))ᴴ * vacProj q) = vacProj q := by
    rw [conjTranspose_mul, conjTranspose_conjTranspose, vacProj_conjTranspose, Matrix.mul_assoc,
      ← Matrix.mul_assoc (cplx (dressed y)), hCAR, Matrix.sub_mul, Matrix.one_mul,
      Matrix.mul_assoc, cplx_dressed_mul_vacProj, Matrix.mul_zero, sub_zero, vacProj_mul_self]
  rw [h, Matrix.mul_zero] at h1
  exact vacProj_ne_zero q h1.symm

/-! ### The number commutator in a sector -/

/-- The row `z ↦ U_q(t)_{xz}` with the entry `z = y` removed. -/
def rowRest (J h : ℝ) (q : Cfg N) (t : ℝ) (x y : Fin N) : Fin N → ℂ :=
  fun z => if z = y then 0 else evo (hq J h q) t x z

/-- `∑_{z ≠ y} |U_q(t)_{xz}|² = 1 - |U_q(t)_{xy}|²`. -/
lemma sum_rowRest (J h : ℝ) (q : Cfg N) (t : ℝ) (x y : Fin N) :
    ∑ z, rowRest J h q t x y z * star (rowRest J h q t x y z) =
      1 - ((‖evo (hq J h q) t x y‖ ^ 2 : ℝ) : ℂ) := by
  have hrow : ∑ z, evo (hq J h q) t x z * star (evo (hq J h q) t x z) = 1 := by
    have h1 := congrFun (congrFun (evo_unitary' (hq J h q) (hq_isSymm J h q) t) x) x
    simpa [Matrix.mul_apply, conjTranspose_apply] using h1
  have hsplit : ∀ z, evo (hq J h q) t x z * star (evo (hq J h q) t x z) =
      rowRest J h q t x y z * star (rowRest J h q t x y z) +
        (if z = y then evo (hq J h q) t x y * star (evo (hq J h q) t x y) else 0) := by
    intro z
    by_cases hz : z = y
    · subst hz
      simp [rowRest]
    · simp [rowRest, hz]
  have ha : evo (hq J h q) t x y * star (evo (hq J h q) t x y) =
      ((‖evo (hq J h q) t x y‖ ^ 2 : ℝ) : ℂ) := by
    rw [Complex.star_def, Complex.mul_conj']
    push_cast
    ring
  rw [Finset.sum_congr rfl (fun z _ => hsplit z), Finset.sum_add_distrib, Finset.sum_ite_eq'] at hrow
  simp only [Finset.mem_univ, ite_true] at hrow
  rw [ha] at hrow
  linear_combination hrow

/-- `∑_z (δ_y)_z g_z = g_y`. -/
lemma sum_single_mul (y : Fin N) (g : Fin N → ℂ) :
    ∑ z, (Pi.single y 1 : Fin N → ℂ) z * g z = g y := by
  rw [Finset.sum_eq_single y (fun z _ hz => by rw [Pi.single_eq_of_ne hz, zero_mul])
    (fun h => absurd (Finset.mem_univ y) h), Pi.single_eq_same, one_mul]

/-- `c̃_y` and `r = ∑_{z ≠ y} U_q(t)_{xz} c̃_z` are two orthogonal modes, `{r, r†} = ∑_{z ≠ y} |U_{xz}|²`. -/
lemma twoModes_row (J h : ℝ) (q : Cfg N) (t : ℝ) (x y : Fin N) :
    TwoModes (cplx (dressed y)) (cplx (dressed y))ᴴ (modeOp (rowRest J h q t x y))
      (modeOp (rowRest J h q t x y))ᴴ
      (∑ z, rowRest J h q t x y z * star (rowRest J h q t x y z)) := by
  have H0 := twoModes_modeOp (Pi.single y 1) (rowRest J h q t x y)
    (by rw [sum_single_mul, Pi.single_eq_same, star_one])
    (by rw [sum_single_mul]; simp [rowRest])
  rwa [modeOp_single] at H0

/-- `P [X, n] = P [M, n]` if `P X = P M` and `P n = n P`. -/
lemma mul_comm_sub_eq {ι : Type*} [Fintype ι] [DecidableEq ι] (P X M n : Matrix ι ι ℂ)
    (key : P * X = P * M) (hPn : P * n = n * P) : P * (X * n - n * X) = P * (M * n - n * M) := by
  calc P * (X * n - n * X) = P * X * n - (P * n) * X := by
        rw [Matrix.mul_sub, ← Matrix.mul_assoc, ← Matrix.mul_assoc]
    _ = P * M * n - n * (P * X) := by rw [hPn, Matrix.mul_assoc n P X, key]
    _ = P * M * n - (P * n) * M := by rw [key, ← Matrix.mul_assoc, ← hPn]
    _ = P * (M * n - n * M) := by rw [Matrix.mul_sub, Matrix.mul_assoc, Matrix.mul_assoc]

/-- `[n_x(t), n_y]` commutes with every sector projection. -/
lemma gaugeInv_numComm (J h : ℝ) (x y : Fin N) (t : ℝ) : GaugeInv (numComm J h x y t) := by
  refine gaugeInv_of_comm_gaussProj _ fun q => ?_
  have c1 : Commute (cplx (gaussProj q)) (heisT J h t (cplx (dressed x)ᵀ)) :=
    comm_gaussProj_heis_dressedSharp J h q t x true
  have c2 : Commute (cplx (gaussProj q)) (heisT J h t (cplx (dressed x))) :=
    comm_gaussProj_heis_dressedSharp J h q t x false
  have c3 : Commute (cplx (gaussProj q)) (cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) :=
    comm_gaussProj_number q y
  have c12 := c1.mul_right c2
  unfold numComm
  rw [← cplx_dressed_number x, heis_mul]
  exact ((c12.mul_right c3).sub_right (c3.mul_right c12)).eq

/-- **The number commutator in a sector**: `P_q [n_x(t), n_y] P_q = P_q K P_q` with
`K = a r†c̃_y - ā c̃_y† r`, `a = U_q(t)_{xy}` and `r = ∑_{z ≠ y} U_q(t)_{xz} c̃_z`. -/
lemma block_numComm (J h : ℝ) (x y : Fin N) (t : ℝ) (q : Cfg N) :
    cplx (gaussProj q) * numComm J h x y t * cplx (gaussProj q) =
      cplx (gaussProj q) * (evo (hq J h q) t x y •
          ((modeOp (rowRest J h q t x y))ᴴ * cplx (dressed y)) -
        star (evo (hq J h q) t x y) • ((cplx (dressed y))ᴴ * modeOp (rowRest J h q t x y))) *
      cplx (gaussProj q) := by
  set P := cplx (gaussProj q) with hP
  set u : Fin N → ℂ := fun z => evo (hq J h q) t x z with hu_def
  set a := evo (hq J h q) t x y with ha_def
  set v := rowRest J h q t x y with hv_def
  set ny := cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) with hny_def
  have hu : u = a • Pi.single y 1 + v := by
    funext z
    by_cases hz : z = y
    · subst hz
      simp [hu_def, ha_def, hv_def, rowRest]
    · simp [hu_def, hv_def, rowRest, hz]
  have hA : modeOp u = a • cplx (dressed y) + modeOp v := by
    rw [hu, modeOp_add, modeOp_smul, modeOp_single]
  have hAd : (modeOp u)ᴴ = star a • (cplx (dressed y))ᴴ + (modeOp v)ᴴ := by
    rw [hA, conjTranspose_add, conjTranspose_smul]
  have hny : ny = (cplx (dressed y))ᴴ * cplx (dressed y) := by
    rw [hny_def, ← cplx_dressed_number, cplx_conjTranspose']
  have hPa : P * heisT J h t (cplx (dressed x)) = P * modeOp u :=
    (gaussProj_heis_dressedSharp J h q t x false).1
  have hPd : P * heisT J h t (cplx (dressed x)ᵀ) = P * (modeOp u)ᴴ := by
    rw [modeOp_conjTranspose]
    exact (gaussProj_heis_dressedSharp J h q t x true).1
  have hPu : P * (modeOp u)ᴴ = (modeOp u)ᴴ * P := comm_gaussProj_modeOp_conjTranspose q u
  have hPn : P * ny = ny * P := comm_gaussProj_number q y
  have key : P * (heisT J h t (cplx (dressed x)ᵀ) * heisT J h t (cplx (dressed x))) =
      P * ((modeOp u)ᴴ * modeOp u) := by
    rw [← Matrix.mul_assoc, hPd, hPu, Matrix.mul_assoc, hPa, ← Matrix.mul_assoc, ← hPu,
      Matrix.mul_assoc]
  have hnum : numComm J h x y t =
      heisT J h t (cplx (dressed x)ᵀ) * heisT J h t (cplx (dressed x)) * ny -
        ny * (heisT J h t (cplx (dressed x)ᵀ) * heisT J h t (cplx (dressed x))) := by
    rw [hny_def]
    unfold numComm
    rw [← cplx_dressed_number x, heis_mul]
  have hK : (modeOp u)ᴴ * modeOp u * ny - ny * ((modeOp u)ᴴ * modeOp u) =
      a • ((modeOp v)ᴴ * cplx (dressed y)) - star a • ((cplx (dressed y))ᴴ * modeOp v) := by
    rw [hAd, hA, hny]
    exact (twoModes_row J h q t x y).comm_number a
  rw [hnum, mul_comm_sub_eq P _ _ ny key hPn, hK]

/-- `‖P_q K P_q‖ = |a| (1 - |a|²)^{1/2}` for the sector operator `K` of `block_numComm`. -/
lemma norm_block_K (J h : ℝ) (x y : Fin N) (t : ℝ) (q : Cfg N) :
    ‖cplx (gaussProj q) * (evo (hq J h q) t x y •
          ((modeOp (rowRest J h q t x y))ᴴ * cplx (dressed y)) -
        star (evo (hq J h q) t x y) • ((cplx (dressed y))ᴴ * modeOp (rowRest J h q t x y))) *
      cplx (gaussProj q)‖ = ‖evo (hq J h q) t x y‖ * √(1 - ‖evo (hq J h q) t x y‖ ^ 2) := by
  have H := twoModes_row J h q t x y
  set P := cplx (gaussProj q) with hP
  set a := evo (hq J h q) t x y with ha_def
  set v := rowRest J h q t x y with hv_def
  set e := cplx (dressed y) with he_def
  set r := modeOp v with hr_def
  set s := ∑ z, v z * star (v z) with hs_def
  set K := a • (rᴴ * e) - star a • (eᴴ * r) with hK
  have hKh : Kᴴ = star a • (eᴴ * r) - a • (rᴴ * e) := by
    rw [hK, conjTranspose_sub, conjTranspose_smul, conjTranspose_smul, conjTranspose_mul,
      conjTranspose_mul, conjTranspose_conjTranspose, conjTranspose_conjTranspose, star_star]
  set D := rᴴ * r + s • (eᴴ * e) + (2 : ℂ) • (eᴴ * (rᴴ * (e * r))) with hD
  have hKK : Kᴴ * K = (star a * a) • D := by
    rw [hKh, hK]
    exact H.kdk a
  have hs : s = 1 - ((‖a‖ ^ 2 : ℝ) : ℂ) := sum_rowRest J h q t x y
  have ha2 : star a * a = ((‖a‖ ^ 2 : ℝ) : ℂ) := by
    rw [Complex.star_def, Complex.conj_mul']
    push_cast
    ring
  set l : ℝ := ‖a‖ ^ 2 * (1 - ‖a‖ ^ 2) with hl_def
  have hl : (l : ℂ) = star a * a * s := by
    rw [hs, ha2, hl_def]
    push_cast
    ring
  have ha1 : ‖a‖ ≤ 1 := norm_evo_le_one _ (hq_isSymm J h q) t x y
  have hl0 : 0 ≤ l := mul_nonneg (sq_nonneg _) (by nlinarith [norm_nonneg a])
  have c1 : Commute P e := comm_gaussProj_dressedSharp q y false
  have c2 : Commute P eᴴ := by
    show P * eᴴ = eᴴ * P
    rw [he_def, cplx_conjTranspose']
    exact comm_gaussProj_dressedSharp q y true
  have c3 : Commute P r := comm_gaussProj_modeOp q v
  have c4 : Commute P rᴴ := comm_gaussProj_modeOp_conjTranspose q v
  have hPK : P * K = K * P :=
    (((c4.mul_right c1).smul_right a).sub_right ((c2.mul_right c3).smul_right (star a))).eq
  have hT := dressed_transpose_mul_vacProj_ne_zero y q
  have hPT : P * (eᴴ * vacProj q) = eᴴ * vacProj q := by
    rw [← Matrix.mul_assoc, c2.eq, Matrix.mul_assoc, gaussProj_mul_vacProj]
  have hX : (Kᴴ * K) * (Kᴴ * K) = (l : ℂ) • (Kᴴ * K) := by
    rw [hKK, smul_mul_smul_comm, H.dsq, smul_smul, hl, smul_smul]
    congr 1
    ring
  have hXT : (Kᴴ * K) * (eᴴ * vacProj q) = (l : ℂ) • (eᴴ * vacProj q) := by
    rw [hKK, smul_mul_assoc, H.dvec _ (modeOp_mul_vacProj v q), smul_smul, hl]
  rw [norm_block_eq P K _ l hl0 (SectorNormAux.cplx_gaussProj_mul_self q)
    (SectorNormAux.cplx_gaussProj_conjTranspose q) hPK hX hT hPT hXT, hl_def,
    Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq (norm_nonneg _)]

end LRNormAux

open LRNormAux

/-- **`‖[n_x(t), n_y]‖ = max_q |U_q(t)_{xy}| (1 - |U_q(t)_{xy}|²)^{1/2}`** (paper: eq. `eq:main-LRnorm`). -/
theorem norm_comm_nn_eq (J h : ℝ) (x y : Fin N) (t : ℝ) :
    ‖numComm J h x y t‖ =
      univ.sup' univ_nonempty (fun q : Cfg N =>
        ‖evo (hq J h q) t x y‖ * √(1 - ‖evo (hq J h q) t x y‖ ^ 2)) := by
  rw [norm_gaugeInv _ (gaugeInv_numComm J h x y t)]
  exact Finset.sup'_congr univ_nonempty rfl fun q _ => by
    rw [block_numComm, norm_block_K]

/-- **The clean sector forbids a uniform operator-norm bound** (paper: Proposition `prop:LR`(c)):
`sup_t ‖{c̃_1(t), c̃_N†}‖ ≥ (3/(2(N+1)))^{1/2}` for `N ≥ 2`. -/
theorem sup_norm_dressedAnticomm_ge (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    √(3 / (2 * ((N : ℝ) + 1))) ≤
      ⨆ t : ℝ, ‖dressedAnticomm (N := N) J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t‖ := by
  sorry

/-- `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}` for `N ≥ 2` (paper: Proposition `prop:LR`(c)). -/
theorem sup_norm_comm_nn_ge (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    √(3 / (4 * ((N : ℝ) + 1))) ≤ ⨆ t : ℝ, ‖numComm (N := N) J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t‖ := by
  sorry

/-- The sector average is small: `E_q sup_t ‖P_q {c̃_1(t), c̃_N†} P_q‖ ≤ F(N-1)` (paper: Proposition
`prop:LR`(c)). -/
theorem avg_sup_anticomm_corner_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (hN : 1 ≤ N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖cplx (gaussProj q) *
      dressedAnticomm J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t * cplx (gaussProj q)‖) ≤ F (N - 1) := by
  have hb := avg_sup_anticomm_le J h hEC (⟨0, by omega⟩ : Fin N) ⟨N - 1, by omega⟩
  rwa [Nat.dist_zero_left] at hb

end DFL
