import DFL.Main
import DFL.Variants

/-!
# Persistent imbalance on the open chain (paper: Corollary `cor:imbalanceopen`)

The model of Smith–Knolle–Kovrizhin–Moessner without the boundary link is the open chain of
`DFL/Variants.lean` (`openInit`: links `0, …, N-2` in `|↑⟩_z`, boundary link in `|+⟩`). Its sector
law is uniform on the parity hyperplane `{q : ∑_j q_j = π}`, `π` the parity of the matter state, so
the imbalance is the parity-hyperplane average `E'_q` of the sector imbalances.

Three properties of `E'` replace the independence of the charges in `avg_imbInf_uniform`: the
charges of any `N - 1` sites are i.i.d. uniform under `E'`; `E' f ≤ 2 E f` for `f ≥ 0`; and the
charges are pairwise independent (`N ≥ 3`). Hence, under (EC), the infinite-time imbalance of the
open chain is bounded below uniformly in the volume.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The average over the parity hyperplane, `E'_q f = 2^{-(N-1)} ∑_{∑ q = π} f(q)`. -/
def avgPar (π : ZMod 2) (f : Cfg N → ℝ) : ℝ :=
  (1 / 2 ^ (N - 1) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π), f q

/-- The density-wave imbalance of the open chain. -/
def lgtImbalanceOpen (J h : ℝ) (N : ℕ) (hN : 0 < N) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (matter (diagonal fun n => occ (n x)))
      (evo (H J h) t *ᵥ openInit hN (fockVec (cdw N)))).re

namespace OpenImbalanceAux

/-! ### The imbalance of the open chain as a parity-hyperplane average -/

/-- Rows of a unitary matrix are normalized: `∑_z |U_{xz}|² = 1`. -/
lemma sum_normSq_row {ι : Type*} [Fintype ι] [DecidableEq ι] (U : Matrix ι ι ℂ)
    (hU : U * Uᴴ = 1) (x : ι) : ∑ z, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU x) x
  simp only [one_apply_eq] at h
  have : ((∑ z, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum, Matrix.mul_apply]
    refine Finset.sum_congr rfl (fun z _ => ?_)
    rw [conjTranspose_apply, ← Complex.mul_conj]
    rfl
  exact_mod_cast this

/-- `∑_z (-1)^z a_z = 2 ∑_{z even} a_z - ∑_z a_z`. -/
lemma sum_neg_one_pow_mul (a : Fin N → ℝ) :
    ∑ z : Fin N, (-1 : ℝ) ^ z.val * a z = 2 * ∑ z ∈ cdw N, a z - ∑ z, a z := by
  rw [cdw, Finset.sum_filter, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rcases Nat.even_or_odd z.val with he | ho
  · rw [he.neg_one_pow, ite_eq_left he]; ring
  · rw [ho.neg_one_pow, ite_eq_right (Nat.not_even_iff_odd.2 ho)]; ring

/-- `∑_{x < N} (-1)^x = 0` for even `N`. -/
lemma sum_neg_one_pow_fin (hN : Even N) : ∑ x : Fin N, (-1 : ℝ) ^ x.val = 0 := by
  rw [Fin.sum_univ_eq_sum_range (fun i => (-1 : ℝ) ^ i) N]
  obtain ⟨k, rfl⟩ := hN
  induction k with
  | zero => simp
  | succ k ih =>
    rw [show k + 1 + (k + 1) = (k + k) + 1 + 1 by ring, Finset.sum_range_succ,
      Finset.sum_range_succ, ih, pow_succ, Even.neg_one_pow ⟨k, rfl⟩]
    ring

/-- `Re Tr(Γ U Γ U†) = ∑_{x,z} (-1)^x (-1)^z |U_{xz}|²`. -/
lemma trace_stag_re (U : Matrix (Fin N) (Fin N) ℂ) :
    (trace (cplx (stag N) * U * cplx (stag N) * Uᴴ)).re =
      ∑ x : Fin N, ∑ z : Fin N,
        (-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) := by
  have hs : cplx (stag N) = diagonal (fun x : Fin N => (((-1 : ℝ) ^ x.val : ℝ) : ℂ)) := by
    unfold cplx stag; exact diagonal_map (by simp)
  rw [hs, trace, Complex.re_sum]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Matrix.diag_apply, Matrix.mul_apply, Complex.re_sum]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [mul_diagonal, diagonal_mul, conjTranspose_apply]
  have : (((-1 : ℝ) ^ x.val : ℝ) : ℂ) * U x z * (((-1 : ℝ) ^ z.val : ℝ) : ℂ) * star (U x z) =
      (((-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) : ℝ) : ℂ) := by
    rw [Complex.ofReal_mul, Complex.ofReal_mul, ← Complex.mul_conj]
    simp only [Complex.star_def]
    ring
  rw [this, Complex.ofReal_re]

/-- The one-body imbalance through the propagator (`N` even):
`imbOne = (2/N) ∑_x (-1)^x ∑_{z even} |U_{xz}|²`. -/
lemma imbOne_eq (J h : ℝ) (hN : Even N) (q : Cfg N) (t : ℝ) :
    imbOne J h q t = (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
      ∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z) := by
  have hrow := sum_normSq_row (evo (hq J h q) t) (evo_unitary' _ (hq_isSymm J h q) t)
  have h1 : ∀ x : Fin N, ∑ z : Fin N, (-1 : ℝ) ^ x.val *
      ((-1 : ℝ) ^ z.val * Complex.normSq (evo (hq J h q) t x z)) =
      2 * ((-1 : ℝ) ^ x.val * ∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z)) -
        (-1 : ℝ) ^ x.val := by
    intro x
    rw [← Finset.mul_sum, sum_neg_one_pow_mul, hrow x]
    ring
  rw [imbOne, trace_stag_re, Finset.sum_congr rfl (fun x _ => h1 x), Finset.sum_sub_distrib,
    ← Finset.mul_sum, sum_neg_one_pow_fin hN]
  ring

/-- The density of the open-chain state at site `x`: the parity-hyperplane sum of the densities of
the evolved density-wave Fock state, `2^{-(N-1)} ∑_{∑ q = π} ∑_{z even} |(e^{-ith_q})_{xz}|²`. -/
lemma open_density_re (J h : ℝ) (hN0 : 0 < N) (x : Fin N) (t : ℝ) :
    (expect (matter (diagonal fun n => occ (n x)))
      (evo (H J h) t *ᵥ openInit hN0 (fockVec (cdw N)))).re =
      (1 / 2 ^ (N - 1) : ℝ) *
        ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = ∑ j, fockCfg (cdw N) j),
          ∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z) := by
  have hψ : ∀ n, fockVec (cdw N) n ≠ 0 → ∑ j, n j = ∑ j, fockCfg (cdw N) j := by
    intro n hn
    by_cases hne : n = fockCfg (cdw N)
    · rw [hne]
    · exact (hn (by rw [fockVec, Pi.single_eq_of_ne hne])).elim
  rw [sector_average_open J h hN0 _ _ hψ _ (gaugeInv_matter _) t]
  have hterm : ∀ q : Cfg N, expect (sectorBlock (matter (diagonal fun n => occ (n x))) q)
      (evo (Hsec J h q) t *ᵥ fockVec (cdw N)) =
      ((∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z) : ℝ) : ℂ) := by
    intro q
    rw [sectorBlock_matter, Hsec, expect_evo_sub_smul_one, density _ (hq_isSymm J h q)]
  rw [Finset.sum_congr rfl (fun q _ => hterm q), ← Complex.ofReal_sum,
    show (1 / 2 ^ (N - 1) : ℂ) = ((1 / 2 ^ (N - 1) : ℝ) : ℂ) by push_cast; ring,
    ← Complex.ofReal_mul, Complex.ofReal_re]

set_option backward.isDefEq.respectTransparency false in
open scoped Matrix.Norms.Operator in
/-- The propagator `t ↦ e^{-itA}` is continuous. -/
lemma continuous_evo {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    Continuous (fun t : ℝ => evo A t) := by
  unfold evo
  exact NormedSpace.exp_continuous.comp (by fun_prop)

/-- The one-body imbalance is continuous in time. -/
lemma continuous_imbOne (J h : ℝ) (q : Cfg N) : Continuous (fun t => imbOne J h q t) := by
  have hU := continuous_evo (hq J h q)
  unfold imbOne
  exact continuous_const.mul (Complex.continuous_re.comp
    ((((continuous_const.matrix_mul hU).matrix_mul continuous_const).matrix_mul
      hU.matrix_conjTranspose).matrix_trace))

end OpenImbalanceAux

/-- The imbalance of the open chain is the parity-hyperplane average of the sector imbalances. -/
theorem lgtImbalanceOpen_eq (J h : ℝ) (hN0 : 0 < N) (hN : Even N) (t : ℝ) :
    lgtImbalanceOpen J h N hN0 t =
      avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbOne J h q t) := by
  simp only [lgtImbalanceOpen, OpenImbalanceAux.open_density_re J h hN0, avgPar,
    OpenImbalanceAux.imbOne_eq J h hN, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun q _ => Finset.sum_congr rfl fun x _ =>
    Finset.sum_congr rfl fun z _ => ?_
  ring

/-- The time average of the imbalance of the open chain converges to the parity-hyperplane average of
the infinite-time sector imbalances. -/
theorem lgtImbalanceOpen_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN0 : 0 < N) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalanceOpen J h N hN0 t) atTop
      (𝓝 (avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbInf J h q))) := by
  have heq : (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalanceOpen J h N hN0 t) =
      fun T => (1 / 2 ^ (N - 1) : ℝ) *
        ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = ∑ j, fockCfg (cdw N) j),
          ((1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h q t) := by
    funext T
    simp_rw [lgtImbalanceOpen_eq J h hN0 hN, avgPar]
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_finsetSum
      (fun q _ => (OpenImbalanceAux.continuous_imbOne J h q).intervalIntegrable 0 T),
      Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    ring
  rw [heq]
  unfold avgPar
  exact (tendsto_finsetSum _ (fun q _ => imbOne_timeAvg J h hJ q)).const_mul _

namespace OpenImbalanceAux

/-! ### Elementary properties of the parity-hyperplane average -/

/-- `a + a = 0` in `ZMod 2`. -/
lemma zmod2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `a + 1 ≠ a` in `ZMod 2`. -/
lemma zmod2_add_one_ne (a : ZMod 2) : a + 1 ≠ a := by
  revert a; decide

/-- In `ZMod 2`, `a + 1 = b ↔ a ≠ b`. -/
lemma zmod2_add_one_eq_iff (a b : ZMod 2) : a + 1 = b ↔ ¬ a = b := by
  revert a b; decide

/-- `2^N = 2 · 2^{N-1}` for `N ≥ 1`. -/
lemma two_pow_eq_two_mul (hN : 0 < N) : (2 : ℝ) ^ N = 2 * 2 ^ (N - 1) := by
  rw [← pow_succ']
  congr 1
  omega

/-- The parity-hyperplane average of a constant (`N ≥ 1`). -/
lemma avgPar_const (hN : 0 < N) (π : ZMod 2) (c : ℝ) : avgPar π (fun _ : Cfg N => c) = c := by
  unfold avgPar
  rw [Finset.sum_const, card_parity_sector hN π, nsmul_eq_mul]
  push_cast
  field_simp

/-- Monotonicity of the parity-hyperplane average. -/
lemma avgPar_mono (π : ZMod 2) {f g : Cfg N → ℝ} (hfg : ∀ q, f q ≤ g q) :
    avgPar π f ≤ avgPar π g := by
  unfold avgPar
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => hfg q) (by positivity)

/-- Additivity of the parity-hyperplane average. -/
lemma avgPar_add (π : ZMod 2) (f g : Cfg N → ℝ) :
    avgPar π (fun q => f q + g q) = avgPar π f + avgPar π g := by
  unfold avgPar; rw [Finset.sum_add_distrib]; ring

/-- The parity-hyperplane average commutes with subtraction. -/
lemma avgPar_sub (π : ZMod 2) (f g : Cfg N → ℝ) :
    avgPar π (fun q => f q - g q) = avgPar π f - avgPar π g := by
  unfold avgPar; rw [Finset.sum_sub_distrib]; ring

/-- Homogeneity of the parity-hyperplane average. -/
lemma avgPar_mul_left (π : ZMod 2) (c : ℝ) (f : Cfg N → ℝ) :
    avgPar π (fun q => c * f q) = c * avgPar π f := by
  unfold avgPar; rw [← Finset.mul_sum]; ring

/-- The parity-hyperplane average commutes with finite sums. -/
lemma avgPar_sum {ι : Type*} (π : ZMod 2) (s : Finset ι) (f : ι → Cfg N → ℝ) :
    avgPar π (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgPar π (f i) := by
  unfold avgPar; rw [Finset.sum_comm, Finset.mul_sum]

/-- Jensen (Cauchy–Schwarz) for the parity-hyperplane average: `(E' f)² ≤ E' f²`. -/
lemma sq_avgPar_le (hN : 0 < N) (π : ZMod 2) (f : Cfg N → ℝ) :
    avgPar π f ^ 2 ≤ avgPar π (fun q => f q ^ 2) := by
  have h0 : 0 ≤ avgPar π (fun q => (f q - avgPar π f) ^ 2) := by
    have := avgPar_mono π (f := fun _ => (0 : ℝ)) (g := fun q => (f q - avgPar π f) ^ 2)
      (fun q => sq_nonneg _)
    rwa [avgPar_const hN] at this
  have : avgPar π (fun q => (f q - avgPar π f) ^ 2) =
      avgPar π (fun q => f q ^ 2) - avgPar π f ^ 2 := by
    have e : ∀ q, (f q - avgPar π f) ^ 2 = f q ^ 2 + ((-2 * avgPar π f) * f q + avgPar π f ^ 2) :=
      fun q => by ring
    simp_rw [e]
    rw [avgPar_add, avgPar_add, avgPar_mul_left, avgPar_const hN]
    ring
  linarith

/-- **Marginals.** If `g` is invariant under flipping the charge at some site `x₀` (e.g. `g` depends
only on the charges outside `x₀`), its parity-hyperplane average is its Bernoulli average. -/
lemma avgPar_eq_avgQ_of_flip (π : ZMod 2) (x₀ : Fin N) (g : Cfg N → ℝ)
    (hg : ∀ q, g (q + Pi.single x₀ 1) = g q) : avgPar π g = avgQ g := by
  have hsum : ∀ q : Cfg N, ∑ j, (q + Pi.single x₀ 1 : Cfg N) j = ∑ j, q j + 1 := by
    intro q
    simp only [Pi.add_apply, Finset.sum_add_distrib, Fintype.sum_pi_single']
  have hinv : ∀ q : Cfg N, q + Pi.single x₀ 1 + Pi.single x₀ 1 = q := by
    intro q
    funext j
    simp only [Pi.add_apply]
    rw [add_assoc, zmod2_add_self, add_zero]
  have hswap : ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π), g q =
      ∑ q ∈ univ.filter (fun q : Cfg N => ¬ ∑ j, q j = π), g q := by
    refine Finset.sum_nbij' (· + Pi.single x₀ 1) (· + Pi.single x₀ 1) ?_ ?_ ?_ ?_ ?_
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hq ⊢
      rw [hsum, hq]
      exact zmod2_add_one_ne π
    · intro q hq
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hq ⊢
      rw [hsum]
      exact (zmod2_add_one_eq_iff _ _).2 hq
    · intro q _
      exact hinv q
    · intro q _
      exact hinv q
    · intro q _
      exact (hg q).symm
  have htot := Finset.sum_filter_add_sum_filter_not univ (fun q : Cfg N => ∑ j, q j = π) g
  unfold avgPar avgQ
  rw [← htot, ← hswap, two_pow_eq_two_mul x₀.pos]
  field_simp
  ring

/-- **Density.** `E' f ≤ 2 E f` for `f ≥ 0`. -/
lemma avgPar_le_two_avgQ (hN : 0 < N) (π : ZMod 2) {f : Cfg N → ℝ} (hf : ∀ q, 0 ≤ f q) :
    avgPar π f ≤ 2 * avgQ f := by
  unfold avgPar avgQ
  have h1 : ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π), f q ≤ ∑ q, f q :=
    Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) (fun q _ _ => hf q)
  rw [two_pow_eq_two_mul hN]
  calc (1 / 2 ^ (N - 1) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π), f q
      ≤ (1 / 2 ^ (N - 1) : ℝ) * ∑ q, f q := mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = 2 * (1 / (2 * 2 ^ (N - 1)) * ∑ q, f q) := by field_simp

/-- **Pairs.** For `N ≥ 3`, `E' (∑_x r_x)² = N` (the charges are pairwise independent). -/
lemma avgPar_sq_rad_univ (hN : 3 ≤ N) (π : ZMod 2) :
    avgPar π (fun q : Cfg N => (∑ x, ImbalanceAux.rad q x) ^ 2) = N := by
  have hN0 : 0 < N := by omega
  have hpair : ∀ x y : Fin N,
      avgPar π (fun q : Cfg N => ImbalanceAux.rad q x * ImbalanceAux.rad q y) =
        if x = y then 1 else 0 := by
    intro x y
    by_cases hxy : x = y
    · subst hxy
      rw [ite_eq_left rfl]
      have e : ∀ q : Cfg N, ImbalanceAux.rad q x * ImbalanceAux.rad q x = 1 := fun q => by
        rw [← sq, ImbalanceAux.rad_sq]
      simp_rw [e]
      exact avgPar_const hN0 π 1
    · rw [ite_eq_right hxy]
      obtain ⟨z, hzx, hzy⟩ : ∃ z : Fin N, z ≠ x ∧ z ≠ y := by
        have hcard : ({x, y} : Finset (Fin N)).card < (univ : Finset (Fin N)).card := by
          rw [Finset.card_univ, Fintype.card_fin]
          exact Finset.card_le_two.trans_lt (by omega)
        obtain ⟨z, -, hz⟩ := Finset.exists_mem_notMem_of_card_lt_card hcard
        simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hz
        exact ⟨z, hz⟩
      rw [avgPar_eq_avgQ_of_flip π z (fun q : Cfg N => ImbalanceAux.rad q x * ImbalanceAux.rad q y)
        (fun q => by
          rw [ImbalanceAux.rad_flip_ne q (Ne.symm hzx), ImbalanceAux.rad_flip_ne q (Ne.symm hzy)])]
      unfold avgQ
      rw [ImbalanceAux.sum_rad_mul_eq_zero x (fun q => ImbalanceAux.rad q y)
        (fun q => ImbalanceAux.rad_flip_ne q (Ne.symm hxy)), mul_zero]
  have e : ∀ q : Cfg N, (∑ x, ImbalanceAux.rad q x) ^ 2 =
      ∑ x, ∑ y, ImbalanceAux.rad q x * ImbalanceAux.rad q y := fun q => by
    rw [sq, Finset.sum_mul_sum]
  simp_rw [e]
  rw [avgPar_sum]
  simp_rw [avgPar_sum, hpair]
  simp

/-! ### The block argument under the parity-hyperplane average -/

/-- The boundary term `C_B` is nonnegative. -/
lemma Cb_nonneg (J h : ℝ) (q : Cfg N) (B : Finset (Fin N)) : 0 ≤ ImbalanceAux.Cb J h q B :=
  Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => ImbalanceAux.transProb_nonneg _ _ _ _

/-- Block estimate for the parity-hyperplane average: `E' X_I² ≥ h² m²` for an interval `I` of
length `3m²` missing at least one site, with `m ≥ 8κ` (the factor `2` of the density property is
absorbed in the larger block). -/
lemma avgPar_Xb_sq_ge (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (π : ZMod 2) (a m : ℕ)
    (hm : 1 ≤ m) (hκ : 8 * (2 * (∑' d, Real.sqrt (F d)) ^ 2) ≤ m) (hN : a + 3 * m ^ 2 ≤ N)
    (x₀ : Fin N) (hx₀ : x₀ ∉ ImbalanceAux.Ival (N := N) a (3 * m ^ 2)) :
    h ^ 2 * m ^ 2 ≤
      avgPar π (fun q : Cfg N => ImbalanceAux.Xb J h q (ImbalanceAux.Ival a (3 * m ^ 2)) ^ 2) := by
  have hN0 : 0 < N := x₀.pos
  set I : Finset (Fin N) := ImbalanceAux.Ival a (3 * m ^ 2) with hI
  have hS : (m : ℝ) ≤ avgPar π (fun q : Cfg N => |∑ x ∈ I, ImbalanceAux.rad q x|) := by
    rw [avgPar_eq_avgQ_of_flip π x₀ (fun q : Cfg N => |∑ x ∈ I, ImbalanceAux.rad q x|)
      (fun q => by rw [ImbalanceAux.sum_rad_flip I hx₀ q])]
    exact ImbalanceAux.avgQ_abs_S_ge I m hm (ImbalanceAux.card_Ival hN)
  have hC : avgPar π (fun q : Cfg N => ImbalanceAux.Cb J h q I) ≤
      2 * (2 * (∑' d, Real.sqrt (F d)) ^ 2) :=
    (avgPar_le_two_avgQ hN0 π (fun q => Cb_nonneg J h q I)).trans
      (mul_le_mul_of_nonneg_left (ImbalanceAux.avgQ_Cb_le J h hEC a (3 * m ^ 2)) (by norm_num))
  have hX : 2 * |h| * (avgPar π (fun q : Cfg N => |∑ x ∈ I, ImbalanceAux.rad q x|) -
      2 * avgPar π (fun q : Cfg N => ImbalanceAux.Cb J h q I)) ≤
      avgPar π (fun q : Cfg N => |ImbalanceAux.Xb J h q I|) := by
    have := avgPar_mono π (fun q : Cfg N => ImbalanceAux.abs_Xb_ge J h q I)
    rwa [avgPar_mul_left, avgPar_sub, avgPar_mul_left] at this
  have hh0 : 0 ≤ |h| := abs_nonneg h
  have hpos : |h| * m ≤ avgPar π (fun q : Cfg N => |ImbalanceAux.Xb J h q I|) := by
    nlinarith [mul_nonneg hh0 (sub_nonneg.mpr hS), mul_nonneg hh0 (sub_nonneg.mpr hC),
      mul_nonneg hh0 (sub_nonneg.mpr hκ)]
  have hsq := sq_avgPar_le hN0 π (fun q : Cfg N => |ImbalanceAux.Xb J h q I|)
  simp only [sq_abs] at hsq
  calc h ^ 2 * (m : ℝ) ^ 2 = (|h| * m) ^ 2 := by rw [mul_pow, sq_abs]
    _ ≤ (avgPar π fun q : Cfg N => |ImbalanceAux.Xb J h q I|) ^ 2 :=
        pow_le_pow_left₀ (by positivity) hpos 2
    _ ≤ _ := hsq

/-- **Global bound** for the parity-hyperplane average (`N ≥ 3`):
`E' Ī_q ≥ 4h² / (N (2J² + 4h²))`. -/
lemma avgPar_imbInf_ge_global (J h : ℝ) (hJ : J ≠ 0) (hN : 3 ≤ N) (π : ZMod 2) :
    4 * h ^ 2 / (N * (2 * J ^ 2 + 4 * h ^ 2)) ≤ avgPar π (fun q : Cfg N => imbInf J h q) := by
  have hN0 : 0 < N := by omega
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN0
  have hM : 0 < 2 * J ^ 2 + 4 * h ^ 2 := by
    have : 0 < J ^ 2 := by positivity
    positivity
  have hq : ∀ q : Cfg N, 4 * h ^ 2 * (∑ x, ImbalanceAux.rad q x) ^ 2 ≤
      (N * (2 * J ^ 2 + 4 * h ^ 2)) * ∑ α, stagWeight J h q α ^ 2 := by
    intro q
    have cs := Finset.sum_mul_sq_le_sq_mul_sq Finset.univ
      (fun α => (hq_isHermitian J h q).eigenvalues α) (fun α => stagWeight J h q α)
    rw [ImbalanceAux.sum_eig_mul_stag] at cs
    have hE := ImbalanceAux.sum_eig_sq_le J h q
    have ha : 0 ≤ ∑ α, stagWeight J h q α ^ 2 := Finset.sum_nonneg fun _ _ => sq_nonneg _
    calc 4 * h ^ 2 * (∑ x, ImbalanceAux.rad q x) ^ 2
        = (2 * h * ∑ x, ImbalanceAux.rad q x) ^ 2 := by ring
      _ ≤ (∑ α, (hq_isHermitian J h q).eigenvalues α ^ 2) * ∑ α, stagWeight J h q α ^ 2 := cs
      _ ≤ (N * (2 * J ^ 2 + 4 * h ^ 2)) * ∑ α, stagWeight J h q α ^ 2 :=
          mul_le_mul_of_nonneg_right hE ha
  have havg := avgPar_mono π hq
  rw [avgPar_mul_left, avgPar_mul_left, avgPar_sq_rad_univ hN π] at havg
  have heq : avgPar π (fun q : Cfg N => imbInf J h q) =
      (1 / N) * avgPar π (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2) := by
    simp only [imbInf]; exact avgPar_mul_left _ _ _
  rw [heq, div_le_iff₀ (by positivity)]
  set A := avgPar π (fun q : Cfg N => ∑ α, stagWeight J h q α ^ 2)
  have : 1 / (N : ℝ) * A * (N * (2 * J ^ 2 + 4 * h ^ 2)) = A * (2 * J ^ 2 + 4 * h ^ 2) := by
    field_simp
  rw [this]
  have h4 : (N : ℝ) * (4 * h ^ 2) ≤ N * (A * (2 * J ^ 2 + 4 * h ^ 2)) := by linarith
  exact le_of_mul_le_mul_left h4 hN'

end OpenImbalanceAux

/-- **Uniform lower bound for the parity-hyperplane average** under Hypothesis (EC). -/
theorem avgPar_imbInf_uniform (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), 4 ≤ N → ∀ π : ZMod 2,
      c ≤ avgPar π (fun q : Cfg N => imbInf J h q) := by
  set κ : ℝ := 2 * (∑' d, Real.sqrt (F d)) ^ 2 with hκdef
  set m : ℕ := ⌊8 * κ⌋₊ + 1 with hmdef
  have hm1 : 1 ≤ m := Nat.le_add_left 1 _
  have hmκ : 8 * κ ≤ (m : ℝ) := by
    have := Nat.lt_floor_add_one (8 * κ)
    rw [hmdef]; push_cast; linarith
  set ℓ : ℕ := 3 * m ^ 2 with hℓdef
  have hℓ : 0 < ℓ := by positivity
  have hℓR : (ℓ : ℝ) = 3 * (m : ℝ) ^ 2 := by rw [hℓdef]; push_cast; ring
  set M : ℝ := 2 * J ^ 2 + 4 * h ^ 2 with hMdef
  have hJ2 : 0 < J ^ 2 := by positivity
  have hM : 0 < M := by positivity
  have hm' : (0 : ℝ) < m := by exact_mod_cast hm1
  have hh2 : 0 < h ^ 2 := by positivity
  refine ⟨h ^ 2 / (6 * ℓ * M), by positivity, ?_⟩
  intro N hN4 π
  have hN0 : 0 < N := by omega
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN0
  by_cases hNℓ : N < 2 * ℓ
  · -- small volumes: the global bound
    refine le_trans ?_ (OpenImbalanceAux.avgPar_imbInf_ge_global J h hJ (by omega) π)
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    have : (N : ℝ) ≤ 2 * ℓ := by exact_mod_cast hNℓ.le
    nlinarith [mul_pos hh2 hM]
  · -- large volumes: blocks of length `ℓ ≤ N / 2`
    replace hNℓ : 2 * ℓ ≤ N := not_lt.mp hNℓ
    set K := N / ℓ with hKdef
    have hK1 : 1 ≤ K := Nat.div_pos (by omega) hℓ
    have hKN : K * ℓ ≤ N := Nat.div_mul_le_self N ℓ
    have hNK : N ≤ 2 * K * ℓ := by
      have h1 := Nat.div_add_mod N ℓ
      have h2 := Nat.mod_lt N hℓ
      have h3 : ℓ ≤ K * ℓ := Nat.le_mul_of_pos_left ℓ hK1
      rw [← hKdef] at h1
      have : ℓ * K = K * ℓ := Nat.mul_comm _ _
      rw [mul_assoc]
      omega
    have havg := OpenImbalanceAux.avgPar_mono π
      (fun q : Cfg N => ImbalanceAux.imbInf_ge_blocks J h q ℓ K hKN hM hℓ)
    rw [OpenImbalanceAux.avgPar_mul_left, OpenImbalanceAux.avgPar_sum] at havg
    have hblock : ∀ k ∈ Finset.range K, h ^ 2 * m ^ 2 / (ℓ * M) ≤
        avgPar π (fun q : Cfg N => ImbalanceAux.Xb J h q (ImbalanceAux.Ival (k * ℓ) ℓ) ^ 2 /
          (ℓ * M)) := by
      intro k hk
      have hk' : k * ℓ + 3 * m ^ 2 ≤ N := by
        have : (k + 1) * ℓ ≤ K * ℓ :=
          Nat.mul_le_mul_right ℓ (by rw [Finset.mem_range] at hk; omega)
        rw [add_mul, one_mul] at this
        omega
      obtain ⟨x₀, hx₀⟩ : ∃ x₀ : Fin N, x₀ ∉ ImbalanceAux.Ival (N := N) (k * ℓ) ℓ := by
        have hcard : (ImbalanceAux.Ival (N := N) (k * ℓ) ℓ).card <
            (univ : Finset (Fin N)).card := by
          rw [ImbalanceAux.card_Ival hk', Finset.card_univ, Fintype.card_fin]
          omega
        obtain ⟨x, -, hx⟩ := Finset.exists_mem_notMem_of_card_lt_card hcard
        exact ⟨x, hx⟩
      have hb := OpenImbalanceAux.avgPar_Xb_sq_ge J h hEC π (k * ℓ) m hm1 hmκ hk' x₀ hx₀
      have e : avgPar π (fun q : Cfg N =>
          ImbalanceAux.Xb J h q (ImbalanceAux.Ival (k * ℓ) ℓ) ^ 2 / (ℓ * M)) =
          (1 / (ℓ * M)) * avgPar π
            (fun q : Cfg N => ImbalanceAux.Xb J h q (ImbalanceAux.Ival (k * ℓ) ℓ) ^ 2) := by
        rw [← OpenImbalanceAux.avgPar_mul_left]
        congr 1; ext q; ring
      rw [e, div_eq_mul_one_div, mul_comm]
      exact mul_le_mul_of_nonneg_left hb (by positivity)
    have hsum := Finset.sum_le_sum hblock
    rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul] at hsum
    refine le_trans ?_ havg
    refine le_trans ?_ (mul_le_mul_of_nonneg_left hsum (by positivity))
    have hNK' : (N : ℝ) ≤ 2 * K * ℓ := by exact_mod_cast hNK
    have hK' : (0 : ℝ) < K := by exact_mod_cast hK1
    rw [hℓR] at hNK' ⊢
    rw [show 1 / (N : ℝ) * (K * (h ^ 2 * m ^ 2 / (3 * m ^ 2 * M))) = K * h ^ 2 / (3 * N * M) by
      field_simp]
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_pos hh2 hM]

/-- **Persistent imbalance on the open chain** (paper: Corollary `cor:imbalanceopen`): under
Hypothesis (EC), the time-averaged density-wave imbalance of the model of Smith–Knolle–Kovrizhin–
Moessner without boundary link is eventually above a positive constant independent of the volume. -/
theorem lgtImbalanceOpen_lower (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ) (hN0 : 0 < N), Even N → 4 ≤ N →
      ∀ ε > 0, ∀ᶠ T in atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalanceOpen J h N hN0 t := by
  obtain ⟨c, hc, hcN⟩ := avgPar_imbInf_uniform J h hJ hh hEC
  refine ⟨c, hc, fun N hN0 hN hN4 ε hε => ?_⟩
  have hge : c ≤ avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbInf J h q) :=
    hcN N hN4 _
  filter_upwards [(lgtImbalanceOpen_timeAvg J h hJ hN0 hN).eventually
    (lt_mem_nhds (show avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbInf J h q) - ε <
      avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbInf J h q) by linarith))] with T hT
  linarith

end DFL
