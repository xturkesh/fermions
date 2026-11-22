import DFL.LRFull
import DFL.Imbalance
import DFL.Normalization

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
local macro "car_norm" H:term : tactic => `(tactic| (simp only [mul_add, add_mul, mul_sub, sub_mul,
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
      (f z * star (g z')) •
        (if z = z' then (1 : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) else 0) := by
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
  have hCAR : cplx (dressed y) * (cplx (dressed y))ᴴ =
      1 - (cplx (dressed y))ᴴ * cplx (dressed y) := by
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
  rw [Finset.sum_congr rfl (fun z _ => hsplit z), Finset.sum_add_distrib,
    Finset.sum_ite_eq'] at hrow
  simp only [Finset.mem_univ, ite_true] at hrow
  rw [ha] at hrow
  linear_combination hrow

/-- `∑_z (δ_y)_z g_z = g_y`. -/
lemma sum_single_mul (y : Fin N) (g : Fin N → ℂ) :
    ∑ z, (Pi.single y 1 : Fin N → ℂ) z * g z = g y := by
  rw [Finset.sum_eq_single y (fun z _ hz => by rw [Pi.single_eq_of_ne hz, zero_mul])
    (fun h => absurd (Finset.mem_univ y) h), Pi.single_eq_same, one_mul]

/-- `c̃_y` and `r = ∑_{z ≠ y} U_q(t)_{xz} c̃_z` are two orthogonal modes, with
`{r, r†} = ∑_{z ≠ y} |U_{xz}|²`. -/
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

/-! ### Long-time averages -/

section TimeAvg

open Filter Topology

/-- `|U(t)_{xy}|² = ∑_{α,β} p_α p_β cos((E_β - E_α) t)` with `p_α = φ_α(x) φ_α(y)`. -/
lemma normSq_evo_eq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (t : ℝ) (x y : Fin N) :
    ‖evo A t x y‖ ^ 2 = ∑ α, ∑ β, (hA.eigenvectorBasis α x * hA.eigenvectorBasis α y) *
      (hA.eigenvectorBasis β x * hA.eigenvectorBasis β y) *
        Real.cos ((hA.eigenvalues β - hA.eigenvalues α) * t) := by
  have h1 : ‖evo A t x y‖ ^ 2 = (evo A t x y * star (evo A t x y)).re := by
    rw [Complex.star_def, Complex.mul_conj', ← Complex.ofReal_pow, Complex.ofReal_re]
  have hsr : ∀ r : ℝ, star (r : ℂ) = r := fun r => Complex.conj_ofReal r
  rw [h1, OneBodyAux.evo_apply_eq_sum A hA t x y, star_sum, Finset.sum_mul_sum, Complex.re_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  have e1 : Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) *
      ((hA.eigenvectorBasis α x : ℂ) * (hA.eigenvectorBasis α y : ℂ)) *
      star (Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues β : ℂ)) *
        ((hA.eigenvectorBasis β x : ℂ) * (hA.eigenvectorBasis β y : ℂ))) =
      (Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) *
        star (Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues β : ℂ)))) *
      (((hA.eigenvectorBasis α x * hA.eigenvectorBasis α y) *
        (hA.eigenvectorBasis β x * hA.eigenvectorBasis β y) : ℝ) : ℂ) := by
    rw [star_mul, star_mul, hsr, hsr]
    push_cast
    ring
  rw [e1, Complex.re_mul_ofReal, ImbalanceAux.re_phase, mul_comm t]
  ring

/-- `t ↦ U(t)_{xy}` is continuous. -/
lemma continuous_evo_apply (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    Continuous (fun t : ℝ => evo A t x y) := by
  have e : (fun t : ℝ => evo A t x y) = fun t : ℝ =>
      ∑ α, Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) *
        ((hA.eigenvectorBasis α x : ℂ) * (hA.eigenvectorBasis α y : ℂ)) :=
    funext fun t => OneBodyAux.evo_apply_eq_sum A hA t x y
  rw [e]
  refine continuous_finsetSum _ fun α _ => Continuous.mul ?_ continuous_const
  exact Complex.continuous_exp.comp
    (((Complex.continuous_ofReal.neg).mul continuous_const).mul continuous_const)

/-- **Long-time average**: for simple spectrum, `T⁻¹ ∫_0^T |U(t)_{xy}|² dt → R(x, y)`. -/
lemma tendsto_avg_normSq_evo (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues) (x y : Fin N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, ‖evo A t x y‖ ^ 2) atTop
      (𝓝 (transProb A hA x y)) := by
  set E := hA.eigenvalues with hE
  set p : Fin N → ℝ := fun α => hA.eigenvectorBasis α x * hA.eigenvectorBasis α y with hp
  have hcont : ∀ (c ω : ℝ), Continuous (fun t : ℝ => c * Real.cos (ω * t)) := by
    intro c ω
    fun_prop
  have hint : ∀ T, ∫ t in (0 : ℝ)..T, ‖evo A t x y‖ ^ 2 =
      ∑ α, ∑ β, p α * p β * ∫ t in (0 : ℝ)..T, Real.cos ((E β - E α) * t) := by
    intro T
    simp_rw [normSq_evo_eq A hA]
    rw [intervalIntegral.integral_finsetSum (fun α _ => ?_)]
    · refine Finset.sum_congr rfl fun α _ => ?_
      rw [intervalIntegral.integral_finsetSum (fun β _ => ?_)]
      · refine Finset.sum_congr rfl fun β _ => ?_
        exact intervalIntegral.integral_const_mul _ _
      · exact (hcont _ _).intervalIntegrable _ _
    · have : Continuous (fun t : ℝ => ∑ β, p α * p β * Real.cos ((E β - E α) * t)) := by
        fun_prop
      exact this.intervalIntegrable _ _
  have hlim : transProb A hA x y =
      ∑ α, ∑ β, p α * p β * (if E β - E α = 0 then 1 else 0) := by
    rw [transProb]
    refine Finset.sum_congr rfl fun α _ => ?_
    have : ∀ β, (if E β - E α = 0 then (1 : ℝ) else 0) = if β = α then 1 else 0 := by
      intro β
      by_cases hβ : β = α
      · simp [hβ]
      · have : E β - E α ≠ 0 := fun e => hβ (hinj (sub_eq_zero.mp e))
        simp [hβ, this]
    simp only [this, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true, hp]
    ring
  rw [hlim]
  simp_rw [hint, Finset.mul_sum]
  refine tendsto_finsetSum _ fun α _ => tendsto_finsetSum _ fun β _ => ?_
  have := (ImbalanceAux.tendsto_avg_cos (E β - E α)).const_mul (p α * p β)
  refine this.congr fun T => ?_
  ring

/-- `sup_t |U(t)_{xy}|² ≥ R(x, y)` for simple spectrum. -/
lemma transProb_le_of_forall (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues) (x y : Fin N) (M : ℝ)
    (hM : ∀ t, ‖evo A t x y‖ ^ 2 ≤ M) : transProb A hA x y ≤ M := by
  refine le_of_tendsto (tendsto_avg_normSq_evo A hA hinj x y) ?_
  filter_upwards [eventually_gt_atTop 0] with T hT
  have hc : Continuous (fun t : ℝ => ‖evo A t x y‖ ^ 2) :=
    ((continuous_evo_apply A hA x y).norm).pow 2
  have h1 : ∫ t in (0 : ℝ)..T, ‖evo A t x y‖ ^ 2 ≤ ∫ t in (0 : ℝ)..T, M :=
    intervalIntegral.integral_mono_on hT.le (hc.intervalIntegrable _ _)
      intervalIntegrable_const (fun t _ => hM t)
  rw [intervalIntegral.integral_const, sub_zero, smul_eq_mul] at h1
  rw [one_div, inv_mul_le_iff₀ hT]
  linarith

end TimeAvg

/-! ### The clean sector `q = 0` (all `G_j = +1`) -/

section Clean

open Real

/-- `∑_{j<n} cos(2π m j / n) = 0` if `n ∤ m`. -/
lemma sum_range_cos_eq_zero (n m : ℕ) (hn : 0 < n) (hm : ¬ n ∣ m) :
    ∑ j ∈ range n, Real.cos (2 * π * m * j / n) = 0 := by
  have hprim := Complex.isPrimitiveRoot_exp n hn.ne'
  set ζ := Complex.exp (2 * π * Complex.I / n) with hζ
  have h1 : ζ ^ m ≠ 1 := fun h => hm ((hprim.pow_eq_one_iff_dvd m).mp h)
  have h2 : (ζ ^ m) ^ n = 1 := by rw [← pow_mul, mul_comm, pow_mul, hprim.pow_eq_one, one_pow]
  have h3 : ∑ j ∈ range n, (ζ ^ m) ^ j = 0 := by rw [geom_sum_eq h1, h2, sub_self, zero_div]
  have h4 : ∀ j : ℕ, ((ζ ^ m) ^ j).re = Real.cos (2 * π * m * j / n) := by
    intro j
    rw [← pow_mul, hζ, ← Complex.exp_nat_mul]
    have e : ((m * j : ℕ) : ℂ) * (2 * π * Complex.I / n) =
        ((2 * π * m * j / n : ℝ) : ℂ) * Complex.I := by
      push_cast
      ring
    rw [e, Complex.exp_ofReal_mul_I_re]
  calc ∑ j ∈ range n, Real.cos (2 * π * m * j / n) = (∑ j ∈ range n, (ζ ^ m) ^ j).re := by
        rw [Complex.re_sum]
        exact Finset.sum_congr rfl fun j _ => (h4 j).symm
    _ = 0 := by rw [h3, Complex.zero_re]

/-- `∑_{j=1}^{n} cos(2π m j/(n+1)) = -1` if `n + 1 ∤ m`. -/
lemma sum_cos_succ (n m : ℕ) (hm : ¬ (n + 1) ∣ m) :
    ∑ j : Fin n, Real.cos (2 * π * m * ((j : ℝ) + 1) / ((n : ℝ) + 1)) = -1 := by
  have h := sum_range_cos_eq_zero (n + 1) m (Nat.succ_pos n) hm
  rw [Finset.sum_range_succ'] at h
  rw [Fin.sum_univ_eq_sum_range (fun j => Real.cos (2 * π * m * ((j : ℝ) + 1) / ((n : ℝ) + 1))) n]
  push_cast at h
  simp only [mul_zero, zero_div, Real.cos_zero] at h
  linarith

/-- `sin⁴θ = 3/8 - cos(2θ)/2 + cos(4θ)/8`. -/
lemma sin_pow_four (θ : ℝ) :
    Real.sin θ ^ 4 = 3 / 8 - Real.cos (2 * θ) / 2 + Real.cos (4 * θ) / 8 := by
  have h2 : Real.cos (2 * θ) = 1 - 2 * Real.sin θ ^ 2 := by
    rw [Real.cos_two_mul, Real.cos_sq']
    ring
  have h4 : Real.cos (4 * θ) = 2 * Real.cos (2 * θ) ^ 2 - 1 := by
    rw [show 4 * θ = 2 * (2 * θ) by ring, Real.cos_two_mul]
  rw [h4, h2]
  ring

/-- The wave number angle `θ_k = π k/(N+1)`, `k = κ + 1 ∈ {1, …, N}`. -/
def sineAngle (κ : Fin N) : ℝ := π * ((κ : ℝ) + 1) / ((N : ℝ) + 1)

/-- The sine mode `φ_k(x) = (2/(N+1))^{1/2} sin(θ_k (x+1))`. -/
def sineMode (κ x : Fin N) : ℝ := √(2 / ((N : ℝ) + 1)) * Real.sin (sineAngle κ * ((x : ℝ) + 1))

/-- The eigenvalue `-2J cos θ_k + 2h` of the sine mode `φ_k` in the clean sector. -/
def sineEig (J h : ℝ) (κ : Fin N) : ℝ := -2 * J * Real.cos (sineAngle κ) + 2 * h

/-- `∑_y A_{xy} f(y+1) = f(x) + f(x+2)` for the path adjacency matrix `A`, if
`f(0) = f(N+1) = 0`. -/
lemma sum_pathAdj_mul (f : ℕ → ℝ) (h0 : f 0 = 0) (hN : f (N + 1) = 0) (x : Fin N) :
    ∑ y : Fin N, pathAdj N x y * f (y.val + 1) = f x.val + f (x.val + 2) := by
  have hsplit : ∀ y : Fin N, pathAdj N x y * f (y.val + 1) =
      (if y.val = x.val + 1 then f (y.val + 1) else 0) +
        (if y.val + 1 = x.val then f (y.val + 1) else 0) := by
    intro y
    simp only [pathAdj]
    by_cases h1 : y.val = x.val + 1
    · have h2 : ¬ y.val + 1 = x.val := by omega
      rw [ite_eq_left (Or.inl h1.symm), ite_eq_left h1, ite_eq_right h2]
      ring
    · by_cases h2 : y.val + 1 = x.val
      · rw [ite_eq_left (Or.inr h2), ite_eq_right h1, ite_eq_left h2]
        ring
      · have h3 : ¬ (x.val + 1 = y.val ∨ y.val + 1 = x.val) := by omega
        rw [ite_eq_right h3, ite_eq_right h1, ite_eq_right h2]
        ring
  rw [Finset.sum_congr rfl fun y _ => hsplit y, Finset.sum_add_distrib,
    Fin.sum_univ_eq_sum_range (fun j => if j = x.val + 1 then f (j + 1) else 0) N,
    Fin.sum_univ_eq_sum_range (fun j => if j + 1 = x.val then f (j + 1) else 0) N,
    Finset.sum_ite_eq', add_comm]
  congr 1
  · rcases Nat.eq_zero_or_pos x.val with hx | hx
    · rw [hx, h0]
      exact Finset.sum_eq_zero fun j _ => ite_eq_right (by omega)
    · obtain ⟨k, hk⟩ : ∃ k, x.val = k + 1 := ⟨x.val - 1, by omega⟩
      rw [hk]
      have e : ∀ i : ℕ, (if i + 1 = k + 1 then f (i + 1) else 0) =
          if i = k then f (k + 1) else 0 := by
        intro i
        by_cases hi : i = k
        · rw [ite_eq_left (by omega), ite_eq_left hi, hi]
        · rw [ite_eq_right (by omega), ite_eq_right hi]
      rw [Finset.sum_congr rfl fun i _ => e i, Finset.sum_ite_eq', ite_eq_left]
      have := x.isLt
      simp only [Finset.mem_range]
      omega
  · split_ifs with hx
    · rfl
    · have : x.val + 2 = N + 1 := by
        have := x.isLt
        simp only [Finset.mem_range] at hx
        omega
      rw [this, hN]

/-- In the clean sector `h_0 = -J A + 2h`. -/
lemma sum_hq_zero_mul (J h : ℝ) (v : Fin N → ℝ) (x : Fin N) :
    ∑ y, hq J h (0 : Cfg N) x y * v y = -J * ∑ y, pathAdj N x y * v y + 2 * h * v x := by
  simp [hq, Matrix.add_apply, Matrix.one_apply, add_mul, Finset.sum_add_distrib, Finset.mul_sum,
    mul_assoc, sgn]

/-- The sine modes are eigenvectors of `h_0`. -/
lemma sineMode_eig (J h : ℝ) (κ x : Fin N) :
    ∑ y, hq J h (0 : Cfg N) x y * sineMode κ y = sineEig J h κ * sineMode κ x := by
  set c := √(2 / ((N : ℝ) + 1)) with hc
  set θ := sineAngle κ with hθ
  set f : ℕ → ℝ := fun j => c * Real.sin (θ * j) with hf
  have hmode : ∀ y : Fin N, sineMode κ y = f (y.val + 1) := by
    intro y
    simp only [sineMode, hf, hc, hθ]
    push_cast
    ring_nf
  have h0 : f 0 = 0 := by simp [hf]
  have hN : f (N + 1) = 0 := by
    simp only [hf]
    have e : θ * ((N + 1 : ℕ) : ℝ) = ((κ.val + 1 : ℕ) : ℝ) * π := by
      rw [hθ, sineAngle]
      push_cast
      field_simp
    rw [e, Real.sin_nat_mul_pi, mul_zero]
  rw [sum_hq_zero_mul]
  simp only [hmode]
  rw [sum_pathAdj_mul f h0 hN x]
  simp only [hf, sineEig, ← hθ]
  have e1 : θ * ((x.val : ℕ) : ℝ) = θ * ((x.val + 1 : ℕ) : ℝ) - θ := by push_cast; ring
  have e2 : θ * ((x.val + 2 : ℕ) : ℝ) = θ * ((x.val + 1 : ℕ) : ℝ) + θ := by push_cast; ring
  rw [e1, e2, Real.sin_sub, Real.sin_add]
  ring

/-- The sine modes are normalized. -/
lemma sineMode_norm (κ : Fin N) : ∑ x, sineMode κ x ^ 2 = 1 := by
  have hN1 : (0 : ℝ) < (N : ℝ) + 1 := by positivity
  have hc2 : √(2 / ((N : ℝ) + 1)) ^ 2 = 2 / ((N : ℝ) + 1) := Real.sq_sqrt (by positivity)
  have hdvd : ¬ (N + 1) ∣ (κ.val + 1) := fun hd => by
    have := Nat.le_of_dvd (Nat.succ_pos _) hd
    have := κ.isLt
    omega
  have hcos := sum_cos_succ N (κ.val + 1) hdvd
  have e : ∀ x : Fin N, sineMode κ x ^ 2 = 2 / ((N : ℝ) + 1) *
      (1 / 2 - Real.cos (2 * π * ((κ.val + 1 : ℕ) : ℝ) * ((x : ℝ) + 1) / ((N : ℝ) + 1)) / 2) := by
    intro x
    have a : 2 * (π * ((κ : ℝ) + 1) / ((N : ℝ) + 1) * ((x : ℝ) + 1)) =
        2 * π * ((κ.val + 1 : ℕ) : ℝ) * ((x : ℝ) + 1) / ((N : ℝ) + 1) := by
      push_cast
      ring
    rw [sineMode, mul_pow, hc2, Real.sin_sq_eq_half_sub, sineAngle, a]
  rw [Finset.sum_congr rfl fun x _ => e x, ← Finset.mul_sum, Finset.sum_sub_distrib,
    ← Finset.sum_div, ← Finset.sum_div, hcos]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
  field_simp
  ring

/-- Distinct sine modes have distinct eigenvalues (`J ≠ 0`). -/
lemma sineEig_injective (J h : ℝ) (hJ : J ≠ 0) : Function.Injective (sineEig (N := N) J h) := by
  intro κ κ' e
  have hN1 : (0 : ℝ) < (N : ℝ) + 1 := by positivity
  have hmem : ∀ κ : Fin N, sineAngle κ ∈ Set.Icc 0 π := by
    intro κ
    have hk : (κ : ℝ) + 1 ≤ (N : ℝ) + 1 := by
      have := κ.isLt
      have : (κ : ℝ) ≤ N := by exact_mod_cast this.le
      linarith
    refine ⟨by unfold sineAngle; positivity, ?_⟩
    unfold sineAngle
    rw [div_le_iff₀ hN1]
    nlinarith [Real.pi_pos]
  have hcos : Real.cos (sineAngle κ) = Real.cos (sineAngle κ') := by
    unfold sineEig at e
    have : -2 * J * (Real.cos (sineAngle κ) - Real.cos (sineAngle κ')) = 0 := by linarith
    rcases mul_eq_zero.mp this with h1 | h1
    · exact absurd (by linarith : J = 0) hJ
    · linarith
  have hθ := Real.injOn_cos (hmem κ) (hmem κ') hcos
  unfold sineAngle at hθ
  have h2 : ((κ : ℕ) : ℝ) = ((κ' : ℕ) : ℝ) := by
    field_simp at hθ
    nlinarith [Real.pi_pos]
  exact Fin.ext (by exact_mod_cast h2)

/-- `∑_{k=1}^{N} sin⁴(π k/(N+1)) = 3(N+1)/8` for `N ≥ 2`. -/
lemma sum_sin_pow_four (hN : 2 ≤ N) :
    ∑ κ : Fin N, Real.sin (sineAngle κ) ^ 4 = 3 * ((N : ℝ) + 1) / 8 := by
  have h1 : ¬ (N + 1) ∣ 1 := fun hd => by have := Nat.le_of_dvd one_pos hd; omega
  have h2 : ¬ (N + 1) ∣ 2 := fun hd => by have := Nat.le_of_dvd two_pos hd; omega
  have hc1 := sum_cos_succ N 1 h1
  have hc2 := sum_cos_succ N 2 h2
  have e : ∀ κ : Fin N, Real.sin (sineAngle κ) ^ 4 =
      3 / 8 - Real.cos (2 * π * ((1 : ℕ) : ℝ) * ((κ : ℝ) + 1) / ((N : ℝ) + 1)) / 2 +
        Real.cos (2 * π * ((2 : ℕ) : ℝ) * ((κ : ℝ) + 1) / ((N : ℝ) + 1)) / 8 := by
    intro κ
    have a2 : 2 * (π * ((κ : ℝ) + 1) / ((N : ℝ) + 1)) =
        2 * π * ((1 : ℕ) : ℝ) * ((κ : ℝ) + 1) / ((N : ℝ) + 1) := by
      push_cast
      ring
    have a4 : 4 * (π * ((κ : ℝ) + 1) / ((N : ℝ) + 1)) =
        2 * π * ((2 : ℕ) : ℝ) * ((κ : ℝ) + 1) / ((N : ℝ) + 1) := by
      push_cast
      ring
    rw [sin_pow_four, sineAngle, a2, a4]
  rw [Finset.sum_congr rfl fun κ _ => e κ, Finset.sum_add_distrib, Finset.sum_sub_distrib,
    ← Finset.sum_div, ← Finset.sum_div, ← Finset.sum_div, hc1, hc2]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  ring

/-- **The clean-sector transition probability** `R(1, N) = ∑_k φ_k(1)⁴ = 3/(2(N+1))` for `N ≥ 2`. -/
lemma transProb_clean (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    transProb (hq J h (0 : Cfg N)) (hq_isHermitian J h 0) ⟨0, by omega⟩ ⟨N - 1, by omega⟩ =
      3 / (2 * ((N : ℝ) + 1)) := by
  set hA := hq_isHermitian J h (0 : Cfg N) with hAdef
  have hinj := hq_eigenvalues_injective J h hJ (0 : Cfg N)
  have key : ∀ κ : Fin N, ∃ β, hA.eigenvalues β = sineEig J h κ ∧
      ∀ x, |sineMode κ x| = |hA.eigenvectorBasis β x| := fun κ =>
    NormalizationAux.exists_abs_eq_eigvec _ hA hinj _ _ (sineMode_eig J h κ) (sineMode_norm κ)
  choose σ hσ1 hσ2 using key
  have hσinj : Function.Injective σ := by
    intro κ κ' e
    exact sineEig_injective J h hJ (by rw [← hσ1 κ, ← hσ1 κ', e])
  have hN1 : (0 : ℝ) < (N : ℝ) + 1 := by positivity
  have hc2 : √(2 / ((N : ℝ) + 1)) ^ 2 = 2 / ((N : ℝ) + 1) := Real.sq_sqrt (by positivity)
  set x0 : Fin N := ⟨0, by omega⟩ with hx0
  set x1 : Fin N := ⟨N - 1, by omega⟩ with hx1
  have hterm : ∀ κ : Fin N,
      (hA.eigenvectorBasis (σ κ) x0) ^ 2 * (hA.eigenvectorBasis (σ κ) x1) ^ 2 =
      (2 / ((N : ℝ) + 1)) ^ 2 * Real.sin (sineAngle κ) ^ 4 := by
    intro κ
    rw [← sq_abs (hA.eigenvectorBasis (σ κ) x0), ← hσ2 κ x0, sq_abs,
      ← sq_abs (hA.eigenvectorBasis (σ κ) x1), ← hσ2 κ x1, sq_abs, sineMode, sineMode]
    have e0 : sineAngle κ * ((x0 : ℝ) + 1) = sineAngle κ := by simp [hx0]
    have e1 : sineAngle κ * ((x1 : ℝ) + 1) = ((κ.val + 1 : ℕ) : ℝ) * π - sineAngle κ := by
      have hc : ((x1 : ℕ) : ℝ) + 1 = N := by
        rw [hx1]
        simp only
        rw [Nat.cast_sub (by omega : 1 ≤ N)]
        push_cast
        ring
      rw [hc, sineAngle]
      push_cast
      field_simp
      ring
    have hsq : ((-1 : ℝ) ^ (κ.val + 1)) ^ 2 = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, neg_one_sq, one_pow]
    rw [e0, e1, Real.sin_nat_mul_pi_sub, mul_pow, mul_pow, hc2, neg_sq, mul_pow, hsq, one_mul]
    ring
  rw [transProb, ← hσinj.bijective_of_finite.sum_comp
    (fun β => (hA.eigenvectorBasis β x0) ^ 2 * (hA.eigenvectorBasis β x1) ^ 2)]
  simp only
  rw [Finset.sum_congr rfl fun κ _ => hterm κ, ← Finset.mul_sum, sum_sin_pow_four hN]
  field_simp
  ring

/-- `U(0) = 1`. -/
lemma evo_zero_eq {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) : evo A 0 = 1 := by
  simp [evo]

end Clean

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
  set x0 : Fin N := ⟨0, by omega⟩ with hx0
  set x1 : Fin N := ⟨N - 1, by omega⟩ with hx1
  have hbdd : BddAbove (Set.range fun t : ℝ => ‖dressedAnticomm (N := N) J h x0 x1 t‖) := by
    refine ⟨1, ?_⟩
    rintro _ ⟨t, rfl⟩
    dsimp only
    rw [norm_dressedAnticomm_eq]
    exact Finset.sup'_le _ _ fun q _ => norm_evo_le_one _ (hq_isSymm J h q) t _ _
  have hle : ∀ t : ℝ, ‖evo (hq J h (0 : Cfg N)) t x0 x1‖ ≤
      ⨆ t : ℝ, ‖dressedAnticomm (N := N) J h x0 x1 t‖ := by
    intro t
    refine le_trans ?_ (le_ciSup hbdd t)
    rw [norm_dressedAnticomm_eq]
    exact Finset.le_sup' (fun q : Cfg N => ‖evo (hq J h q) t x0 x1‖) (Finset.mem_univ 0)
  set S := ⨆ t : ℝ, ‖dressedAnticomm (N := N) J h x0 x1 t‖ with hS
  have hS0 : 0 ≤ S := (norm_nonneg _).trans (hle 0)
  have htp := transProb_le_of_forall _ (hq_isHermitian J h 0) (hq_eigenvalues_injective J h hJ 0)
    x0 x1 (S ^ 2) fun t => pow_le_pow_left₀ (norm_nonneg _) (hle t) 2
  rw [transProb_clean J h hJ hN] at htp
  calc √(3 / (2 * ((N : ℝ) + 1))) ≤ √(S ^ 2) := Real.sqrt_le_sqrt htp
    _ = S := Real.sqrt_sq hS0

/-- `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}` for `N ≥ 2` (paper: Proposition `prop:LR`(c)). -/
theorem sup_norm_comm_nn_ge (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    √(3 / (4 * ((N : ℝ) + 1))) ≤ ⨆ t : ℝ, ‖numComm (N := N) J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t‖ := by
  set x0 : Fin N := ⟨0, by omega⟩ with hx0
  set x1 : Fin N := ⟨N - 1, by omega⟩ with hx1
  set u : ℝ → ℝ := fun t => ‖evo (hq J h (0 : Cfg N)) t x0 x1‖ with hu
  have hucont : Continuous u := (continuous_evo_apply _ (hq_isHermitian J h 0) x0 x1).norm
  have hx01 : x0 ≠ x1 := by
    rw [hx0, hx1, Ne, Fin.ext_iff]
    simp only
    omega
  have hu0 : u 0 = 0 := by
    simp only [hu]
    rw [evo_zero_eq, Matrix.one_apply_ne hx01, norm_zero]
  have hbdd : BddAbove (Set.range fun t : ℝ => ‖numComm (N := N) J h x0 x1 t‖) := by
    refine ⟨1, ?_⟩
    rintro _ ⟨t, rfl⟩
    dsimp only
    rw [norm_comm_nn_eq]
    refine Finset.sup'_le _ _ fun q _ => ?_
    have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x0 x1
    have h2 : √(1 - ‖evo (hq J h q) t x0 x1‖ ^ 2) ≤ 1 :=
      Real.sqrt_le_one.mpr (by nlinarith [norm_nonneg (evo (hq J h q) t x0 x1)])
    calc ‖evo (hq J h q) t x0 x1‖ * √(1 - ‖evo (hq J h q) t x0 x1‖ ^ 2) ≤ 1 * 1 :=
          mul_le_mul h1 h2 (Real.sqrt_nonneg _) zero_le_one
      _ = 1 := one_mul 1
  set S := ⨆ t : ℝ, ‖numComm (N := N) J h x0 x1 t‖ with hS
  have hg : ∀ t, u t * √(1 - u t ^ 2) ≤ S := by
    intro t
    refine le_trans ?_ (le_ciSup hbdd t)
    rw [norm_comm_nn_eq]
    exact Finset.le_sup' (fun q : Cfg N =>
      ‖evo (hq J h q) t x0 x1‖ * √(1 - ‖evo (hq J h q) t x0 x1‖ ^ 2)) (Finset.mem_univ 0)
  have hS0 : 0 ≤ S := by
    have := hg 0
    rwa [hu0, zero_mul] at this
  set s₀ := √(3 / (2 * ((N : ℝ) + 1))) with hs₀
  have hN1 : (0 : ℝ) < N + 1 := by positivity
  have hs₀0 : 0 ≤ s₀ := Real.sqrt_nonneg _
  have hs₀sq : s₀ ^ 2 = 3 / (2 * ((N : ℝ) + 1)) := Real.sq_sqrt (by positivity)
  have hs₀half : s₀ ^ 2 ≤ 1 / 2 := by
    rw [hs₀sq, div_le_div_iff₀ (by positivity) (by norm_num)]
    have : (2 : ℝ) ≤ N := by exact_mod_cast hN
    nlinarith
  have hreach : ∀ s, 0 ≤ s → s < s₀ → ∃ t, u t = s := by
    intro s hs0 hs
    obtain ⟨t₁, ht₁⟩ : ∃ t₁, s ≤ u t₁ := by
      by_contra hcon
      simp only [not_exists, not_le] at hcon
      have htp := transProb_le_of_forall _ (hq_isHermitian J h 0)
        (hq_eigenvalues_injective J h hJ 0) x0 x1 (s ^ 2) fun t =>
          pow_le_pow_left₀ (norm_nonneg _) (hcon t).le 2
      rw [transProb_clean J h hJ hN, ← hs₀sq] at htp
      nlinarith
    have hmem : s ∈ Set.uIcc (u 0) (u t₁) := by
      rw [hu0]
      exact Set.mem_uIcc.mpr (Or.inl ⟨hs0, ht₁⟩)
    obtain ⟨t, _, ht⟩ := intermediate_value_uIcc hucont.continuousOn hmem
    exact ⟨t, ht⟩
  have hsqrt2 : 0 < √2 := Real.sqrt_pos.mpr two_pos
  have key : ∀ c, c < s₀ / √2 → c ≤ S := by
    intro c hc
    rcases le_or_gt c 0 with hc0 | hc0
    · exact hc0.trans hS0
    · set s := c * √2 with hs_def
      have hs : s < s₀ := by
        rw [lt_div_iff₀ hsqrt2] at hc
        exact hc
      have hs0 : 0 ≤ s := by positivity
      obtain ⟨t, ht⟩ := hreach s hs0 hs
      have h1 := hg t
      rw [ht] at h1
      have hs2 : s ^ 2 ≤ 1 / 2 := by nlinarith
      have h2 : √(1 / 2) ≤ √(1 - s ^ 2) := Real.sqrt_le_sqrt (by linarith)
      have h3 : c = s * √(1 / 2) := by
        rw [hs_def, mul_assoc, ← Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 2),
          show (2 : ℝ) * (1 / 2) = 1 by norm_num, Real.sqrt_one, mul_one]
      calc c = s * √(1 / 2) := h3
        _ ≤ s * √(1 - s ^ 2) := mul_le_mul_of_nonneg_left h2 hs0
        _ ≤ S := h1
  have hfin : s₀ / √2 ≤ S := le_of_forall_lt_imp_le_of_dense key
  have he : √(3 / (4 * ((N : ℝ) + 1))) = s₀ / √2 := by
    rw [hs₀, ← Real.sqrt_div (by positivity)]
    congr 1
    field_simp
    ring
  rw [he]
  exact hfin

/-- The sector average is small: `E_q sup_t ‖P_q {c̃_1(t), c̃_N†} P_q‖ ≤ F(N-1)` (paper: Proposition
`prop:LR`(c)). -/
theorem avg_sup_anticomm_corner_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (hN : 1 ≤ N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖cplx (gaussProj q) *
      dressedAnticomm J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t * cplx (gaussProj q)‖) ≤ F (N - 1) := by
  have hb := avg_sup_anticomm_le J h hEC (⟨0, by omega⟩ : Fin N) ⟨N - 1, by omega⟩
  rwa [Nat.dist_zero_left] at hb

end DFL
