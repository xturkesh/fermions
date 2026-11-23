import DFL.DimerInfVol
import DFL.DualityFull
import DFL.GaussLaw

/-!
# Dimer product states (paper: Definition `def:dimerstate`, Proposition `prop:dimerlaw`,
Theorem `thm:dimer` for the single-fermion dimer state)

Sites and links are 0-indexed: the paper's even links `b_{2k}` are the odd links `b = 2k - 1` here, and the
dimers are the pairs `{2k, 2k+1}`. Gauge states are written in the `σˣ` basis (labels in `ZMod 2`); a
`z`-polarized link is the uniform superposition of its two labels.

* `prodStateXi n ξ` is `|n⟩ ⊗ χ_ξ`: matter `|n⟩`, odd links in the `σˣ` eigenstates `ξ_b`, even links
  `z`-polarized (unnormalized: amplitude `1` on every allowed label configuration).
* `dimerState n` (paper: `Ψ^dim(n)`, Definition `def:dimerstate`) is `prodStateXi n (parityLe n)`: the odd link
  `b` carries the matter parity `N_{≤ b}(n)` to its left. `dimerState 0` is the dimer vacuum `Ω^dim` and
  `dimerState (fockCfg {x})` is the single-fermion state `Ψ^x`.
* `prodStateXi_dual`, `pairDefect_eq_zero_iff`, `gauss_pair_prodStateXi` (paper: `prop:dimerlaw`(e)).
* `dimerState_dual`, `normSq_dimerState`, `normSq_gaussProj_dimerState`, `dimerState_cdw`
  (paper: `prop:dimerlaw`(a)).
* `gauss_pair_dimerState`, `gauss_prefix_dimerState` (paper: `prop:dimerlaw`(b)).
* `dressed_cre_dimerState`, `dimerState_eq_dressed`, `dimerState_single_eq` (paper: `prop:dimerlaw`(c)).
* `sector_average_dimerState`, `dimer_density_single` (paper: `prop:dimerlaw`(d), eqs. `eq:dimeravg` and
  `eq:dimerdensity`).
* `dimerDLGT_eq`, `dimer_LGT_moment_lower`, `dimer_LGT_TL_D`, `dimer_LGT_TL_M`, `dimer_LGT_unbounded`
  (paper: Theorem `thm:dimer`(a),(b),(d) for the moments `D^{(N)}_p`, `M^{(N)}_p` of eq. `eq:LGTmoments`).
-/

open Matrix Finset MeasureTheory Filter Topology
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The matter parity `N_{≤ b}(n)` of the sites `0, …, b` (in `ZMod 2`). -/
def parityLe (n : Cfg N) (b : Fin N) : ZMod 2 := ∑ j ∈ univ.filter (fun j : Fin N => j ≤ b), n j

/-- The product state `|n⟩ ⊗ χ_ξ`: matter `|n⟩`, odd links in the `σˣ` eigenstates `ξ_b`, even links
`z`-polarized (unnormalized). -/
def prodStateXi (n : Cfg N) (ξ : Fin N → ZMod 2) : Cfg N × Cfg N → ℂ :=
  fun p => if p.1 = n ∧ ∀ b : Fin N, b.val % 2 = 1 → p.2 b = ξ b then 1 else 0

/-- **Dimer product states** (paper: Definition `def:dimerstate`): `Ψ^dim(n) = |n⟩ ⊗ χ^dim(n)`, where the odd
link `b` carries the matter parity to its left. -/
def dimerState (n : Cfg N) : Cfg N × Cfg N → ℂ := prodStateXi n (parityLe n)

/-- `ξ_{2k-1}` (`0` for `k = 0`): the label of the odd link to the left of the pair `{2k, 2k+1}`. -/
def prevOdd (ξ : Fin N → ZMod 2) (k : ℕ) : ZMod 2 :=
  if h : 0 < k ∧ 2 * k - 1 < N then ξ ⟨2 * k - 1, h.2⟩ else 0

/-- `d_k(n, ξ) = ξ_{2k+1} + ξ_{2k-1} + n_{2k} + n_{2k+1}` for the pair `{2k, 2k+1}` (paper:
eq. `eq:dimersigns`, additively). -/
def pairDefect (n : Cfg N) (ξ : Fin N → ZMod 2) (k : ℕ) (hk : 2 * k + 1 < N) : ZMod 2 :=
  ξ ⟨2 * k + 1, hk⟩ + prevOdd ξ k + n ⟨2 * k, by omega⟩ + n ⟨2 * k + 1, hk⟩

/-- The moment `D^{(N)}_p(t; x) = ⟨Ψ^x(t), ∑_y |y - x|^p n_y Ψ^x(t)⟩` of the normalized single-fermion dimer
state (paper: eq. `eq:LGTmoments`; `‖Ψ^x‖² = 2^{N/2}` for the unnormalized `dimerState`). -/
def dimerDLGT (J h : ℝ) (N : ℕ) (p t : ℝ) (x : Fin N) : ℝ :=
  (1 / 2 ^ (N / 2) : ℝ) *
    (expect (cplx (diagonal (fun n : Cfg N => ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * occ (n y)) ⊗ₖ
      (1 : Matrix (Cfg N) (Cfg N) ℝ))) (evo (H J h) t *ᵥ dimerState (fockCfg {x}))).re

/-! ### Auxiliary lemmas -/

namespace DimerStateAux

/-- In `ZMod 2`, `a + b = c ↔ a = c + b`. -/
lemma zmod2_add_eq_iff (a b c : ZMod 2) : a + b = c ↔ a = c + b := by
  revert a b c; decide

/-- Two indicators agree when their conditions are equivalent (for any decidability instances). -/
lemma ite_congr_iff {α : Type*} {P Q : Prop} {hP : Decidable P} {hQ : Decidable Q} (h : P ↔ Q)
    (a b : α) : @ite α P hP a b = @ite α Q hQ a b := by
  by_cases hp : P
  · exact (ite_eq_left hp).trans (ite_eq_left (h.mp hp)).symm
  · exact (ite_eq_right hp).trans (ite_eq_right (mt h.mpr hp)).symm

/-- The product of two `0/1` indicators is the indicator of the conjunction. -/
lemma ite_mul_ite {P Q : Prop} {hP : Decidable P} {hQ : Decidable Q} :
    @ite ℂ P hP 1 0 * @ite ℂ Q hQ 1 0 = @ite ℂ (P ∧ Q) instDecidableAnd 1 0 := by
  by_cases hp : P <;> by_cases hq : Q <;> simp [hp, hq]

/-- The product of a real and a complex `0/1` indicator is the indicator of the conjunction. -/
lemma ofReal_ite_mul_ite {P Q : Prop} {hP : Decidable P} {hQ : Decidable Q} :
    ((@ite ℝ P hP 1 0 : ℝ) : ℂ) * @ite ℂ Q hQ 1 0 = @ite ℂ (P ∧ Q) instDecidableAnd 1 0 := by
  by_cases hp : P <;> by_cases hq : Q <;> simp [hp, hq]

/-- `|1_c|² = 1_c`. -/
lemma normSq_ite (c : Prop) [Decidable c] :
    Complex.normSq (if c then (1 : ℂ) else 0) = if c then 1 else 0 := by
  split_ifs <;> simp

/-- `prevOdd` vanishes for the first pair. -/
lemma prevOdd_zero (ξ : Fin N → ZMod 2) : prevOdd ξ 0 = 0 := by
  unfold prevOdd
  exact dite_eq_right (by omega)

/-- `prevOdd ξ (k + 1) = ξ_{2k+1}`. -/
lemma prevOdd_succ (ξ : Fin N → ZMod 2) (k : ℕ) (hk : 2 * k + 1 < N) :
    prevOdd ξ (k + 1) = ξ ⟨2 * k + 1, hk⟩ := by
  have hc : 0 < k + 1 ∧ 2 * (k + 1) - 1 < N := ⟨Nat.succ_pos k, by omega⟩
  unfold prevOdd
  rw [dite_eq_left hc]
  exact congrArg ξ (Fin.ext (by show 2 * (k + 1) - 1 = 2 * k + 1; omega))

/-- `prevOdd` only sees the labels of the odd links. -/
lemma prevOdd_congr {x ξ : Fin N → ZMod 2} (h : ∀ b : Fin N, b.val % 2 = 1 → x b = ξ b)
    (k : ℕ) : prevOdd x k = prevOdd ξ k := by
  unfold prevOdd
  by_cases hk : 0 < k ∧ 2 * k - 1 < N
  · rw [dite_eq_left hk, dite_eq_left hk]
    exact h _ (by show (2 * k - 1) % 2 = 1; omega)
  · rw [dite_eq_right hk, dite_eq_right hk]

/-- `pairDefect` only sees the labels of the odd links. -/
lemma pairDefect_congr (n : Cfg N) {x ξ : Fin N → ZMod 2}
    (h : ∀ b : Fin N, b.val % 2 = 1 → x b = ξ b) (k : ℕ) (hk : 2 * k + 1 < N) :
    pairDefect n x k hk = pairDefect n ξ k hk := by
  unfold pairDefect
  rw [h ⟨2 * k + 1, hk⟩ (by show (2 * k + 1) % 2 = 1; omega), prevOdd_congr h k]

/-- `prevLink` at a site given by its value. -/
lemma prevLink_mk (x : Cfg N) (i : ℕ) (hi : i < N) :
    prevLink x ⟨i, hi⟩ = if h : i = 0 then 0 else x ⟨i - 1, by omega⟩ := rfl

/-- The pair defect of a link configuration is the total charge `q_{2k} + q_{2k+1}` of the pair. -/
lemma pairDefect_eq_chargeOf (n x : Cfg N) (k : ℕ) (hk : 2 * k + 1 < N) :
    pairDefect n x k hk =
      chargeOf n x ⟨2 * k, by omega⟩ + chargeOf n x ⟨2 * k + 1, hk⟩ := by
  have h1 : prevLink x ⟨2 * k, by omega⟩ = prevOdd x k := by
    rw [prevLink_mk]
    unfold prevOdd
    split_ifs <;> first | rfl | omega
  have h2 : prevLink x ⟨2 * k + 1, hk⟩ = x ⟨2 * k, by omega⟩ := by
    rw [prevLink_mk]
    split_ifs <;> first | omega | rfl
  have key : ∀ a b c d e : ZMod 2, a + b + c + d = b + e + c + (e + a + d) := by decide
  simp only [pairDefect, chargeOf, h1, h2]
  exact key _ _ _ _ _

/-- The pair defect of the dual link labels `x(n, q)` is the charge `q_{2k} + q_{2k+1}` of the pair. -/
lemma pairDefect_linkOf (n q : Cfg N) (k : ℕ) (hk : 2 * k + 1 < N) :
    pairDefect n (linkOf n q) k hk = q ⟨2 * k, by omega⟩ + q ⟨2 * k + 1, hk⟩ := by
  rw [pairDefect_eq_chargeOf, chargeOf_linkOf]

/-- Equal pair defects force equal labels on the odd links (induction along the chain). -/
lemma odd_eq_of_pairDefect_eq (n : Cfg N) {x ξ : Fin N → ZMod 2}
    (h : ∀ (k : ℕ) (hk : 2 * k + 1 < N), pairDefect n x k hk = pairDefect n ξ k hk) :
    ∀ (k : ℕ) (hk : 2 * k + 1 < N), x ⟨2 * k + 1, hk⟩ = ξ ⟨2 * k + 1, hk⟩ := by
  intro k
  induction k with
  | zero =>
    intro hk
    have h0 := h 0 hk
    simp only [pairDefect, prevOdd_zero, add_zero, add_left_inj] at h0
    exact h0
  | succ k ih =>
    intro hk
    have h1 := h (k + 1) hk
    simp only [pairDefect, prevOdd_succ x k (by omega), prevOdd_succ ξ k (by omega),
      ih (by omega), add_left_inj] at h1
    exact h1

/-- A statement about all odd links follows from the statements about the links `2k+1`. -/
lemma forall_odd {P : Fin N → Prop} (h : ∀ (k : ℕ) (hk : 2 * k + 1 < N), P ⟨2 * k + 1, hk⟩) :
    ∀ b : Fin N, b.val % 2 = 1 → P b := by
  intro b hb
  have hb' := b.isLt
  have e : b = ⟨2 * (b.val / 2) + 1, by omega⟩ :=
    Fin.ext (by show b.val = 2 * (b.val / 2) + 1; omega)
  rw [e]
  exact h _ _

/-- Two link configurations agree on the odd links iff they have the same pair defects. -/
lemma odd_eq_iff_pairDefect_eq (n : Cfg N) (x ξ : Fin N → ZMod 2) :
    (∀ b : Fin N, b.val % 2 = 1 → x b = ξ b) ↔
      ∀ (k : ℕ) (hk : 2 * k + 1 < N), pairDefect n x k hk = pairDefect n ξ k hk :=
  ⟨fun h k hk => pairDefect_congr n h k hk,
    fun h => forall_odd (P := fun b => x b = ξ b) (odd_eq_of_pairDefect_eq n h)⟩

/-- `N_{≤ j}(n)` are the dual link labels of `n` in the neutral sector `q = 0`. -/
lemma parityLe_eq_linkOf (n : Cfg N) : parityLe n = linkOf n 0 := by
  funext j
  simp [parityLe, linkOf]

/-- The dimer orientation has no pair defects. -/
lemma pairDefect_parityLe (n : Cfg N) (k : ℕ) (hk : 2 * k + 1 < N) :
    pairDefect n (parityLe n) k hk = 0 := by
  rw [parityLe_eq_linkOf, pairDefect_linkOf]
  simp

/-- Dimer charge configurations: every pair `{2k, 2k+1}` is neutral. -/
lemma isDimer_iff (q : Cfg N) :
    IsDimer q ↔ ∀ (k : ℕ) (hk : 2 * k + 1 < N), q ⟨2 * k, by omega⟩ + q ⟨2 * k + 1, hk⟩ = 0 := by
  rw [isDimer_iff_pairs]
  refine forall_congr' fun k => ⟨fun h hk => ?_, fun h hk => ?_⟩
  · have := h hk
    rwa [dimerExt_of_lt _ _ (by omega), dimerExt_of_lt _ _ hk] at this
  · rw [dimerExt_of_lt _ _ (by omega), dimerExt_of_lt _ _ hk]
    exact h hk

/-- `|n⟩ ⊗ χ_ξ` vanishes off its support. -/
lemma prodStateXi_apply_of_not {n : Cfg N} {ξ : Fin N → ZMod 2} {p : Cfg N × Cfg N}
    (hp : ¬ (p.1 = n ∧ ∀ b : Fin N, b.val % 2 = 1 → p.2 b = ξ b)) : prodStateXi n ξ p = 0 :=
  ite_eq_right hp

/-- Values of the dimer state. -/
lemma dimerState_apply (n m x : Cfg N) :
    dimerState n (m, x) =
      if m = n ∧ ∀ b : Fin N, b.val % 2 = 1 → x b = parityLe n b then 1 else 0 := rfl

/-- The dimer state `Ψ^dim(n)` vanishes off the matter configuration `n`. -/
lemma dimerState_apply_of_ne {n m : Cfg N} (hm : m ≠ n) (x : Cfg N) : dimerState n (m, x) = 0 := by
  rw [dimerState_apply]
  exact ite_eq_right (not_and_of_not_left _ hm)

/-- The one-dimensional Gauss law: `∑_{k ≤ b} q_k(n, x) = x_b - N_{≤ b}(n)`. -/
lemma sum_chargeOf_le (m x : Cfg N) (b : Fin N) :
    ∑ k ∈ univ.filter (· ≤ b), chargeOf m x k = x b - parityLe m b := by
  have h1 := congrFun (linkOf_chargeOf m x) b
  simp only [linkOf, Finset.sum_add_distrib] at h1
  rw [← h1, parityLe, add_sub_cancel_right]

/-- Adding a particle on an empty site: `Function.update n s 1 = n + 1_{{s}}` if `n s = 0`. -/
lemma update_one_eq {n : Cfg N} {s : Fin N} (hn : n s = 0) :
    Function.update n s 1 = n + fockCfg {s} := by
  funext k
  by_cases hk : k = s
  · subst hk; simp [fockCfg, hn]
  · simp [fockCfg, hk]

/-- `N_{≤ b}(n + e_s) = N_{≤ b}(n) + 1_{s ≤ b}`. -/
lemma parityLe_add_single (n : Cfg N) (s b : Fin N) :
    parityLe (n + fockCfg {s}) b = parityLe n b + fockCfg (tailOf s) b := by
  simp only [parityLe, Pi.add_apply, Finset.sum_add_distrib]
  congr 1
  simp [fockCfg, tailOf]

/-- Matrix elements of `c_s†` on a configuration that is empty up to `s` (no Jordan–Wigner sign). -/
lemma cre_apply_of_empty {n : Cfg N} {s : Fin N} (hs : ∀ k, k ≤ s → n k = 0) (n' : Cfg N) :
    cre s n' n = if n' = n + fockCfg {s} then 1 else 0 := by
  have hj : jwPar n s = 0 :=
    Finset.sum_eq_zero fun k hk => hs k (Finset.mem_filter.mp hk).2.le
  rw [FreeDynAux.cre_apply, update_one_eq (hs s le_rfl), hj, FreeDynAux.sgn_zero]
  exact ite_congr_iff (and_iff_right (hs s le_rfl)) 1 0

/-- `(flipSet T)_{x', x' + 1_T} = 1`. -/
lemma flipSet_apply_add (T : Finset (Fin N)) (x' : Cfg N) : flipSet T x' (x' + fockCfg T) = 1 :=
  ite_eq_left (by rw [add_assoc, DualityFullAux.cfg_add_self, add_zero])

/-- `(flipSet T)_{x', x} = 0` unless `x = x' + 1_T`. -/
lemma flipSet_apply_of_ne (T : Finset (Fin N)) {x' x : Cfg N} (hx : x ≠ x' + fockCfg T) :
    flipSet T x' x = 0 :=
  ite_eq_right (fun h => hx ((DualityFullAux.cfg_eq_add_comm x' x _).mp h))

/-- Matrix elements of `c̃_s† = c_s† ⊗ ∏_{k ≥ s} σᶻ_{b_k}`. -/
lemma dressedT_apply (s : Fin N) (n' x' m x : Cfg N) :
    cplx (dressed s)ᵀ (n', x') (m, x) = ((cre s n' m * flipSet (tailOf s) x' x : ℝ) : ℂ) := by
  rw [DualityFullAux.dressed_transpose]
  rfl

/-- `1_S + 1_{{s}} = 1_{insert s S}` for `s ∉ S`. -/
lemma fockCfg_insert {S : Finset (Fin N)} {s : Fin N} (hs : s ∉ S) :
    fockCfg S + fockCfg {s} = fockCfg (insert s S) := by
  funext k
  by_cases hk : k = s
  · subst hk; simp [fockCfg, hs]
  · simp [fockCfg, hk]

/-- Expectation of a diagonal matter observable: `Re ⟨v, (D ⊗ 1) v⟩ = ∑_p d(p.1) |v_p|²`. -/
lemma expect_matter_re (d : Cfg N → ℝ) (v : Cfg N × Cfg N → ℂ) :
    (expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) v).re =
      ∑ p, d p.1 * Complex.normSq (v p) := by
  rw [sec_diagonal_kronecker_one, sec_cplx_diagonal]
  unfold expect
  rw [dotProduct, Complex.re_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [mulVec_diagonal, Pi.star_apply]
  have e : star (v p) * ((d p.1 : ℂ) * v p) = ((d p.1 * Complex.normSq (v p) : ℝ) : ℂ) := by
    rw [Complex.ofReal_mul, Complex.normSq_eq_conj_mul_self, Complex.star_def]
    ring
  rw [e, Complex.ofReal_re]

end DimerStateAux

/-! ### Other orientations of the odd links (paper: `prop:dimerlaw`(e)) -/

open Classical in
/-- In the dual frame, `|n⟩ ⊗ χ_ξ` is the uniform superposition of `|n⟩ ⊗ |q⟩` over the charges with
`q_{2k} + q_{2k+1} = d_k(n, ξ)` for every pair (paper: `prop:dimerlaw`(e)). -/
theorem prodStateXi_dual (n : Cfg N) (ξ : Fin N → ZMod 2) (m q : Cfg N) :
    prodStateXi n ξ (dualEquiv (m, q)) =
      if m = n ∧ ∀ (k : ℕ) (hk : 2 * k + 1 < N),
          q ⟨2 * k, by omega⟩ + q ⟨2 * k + 1, hk⟩ = pairDefect n ξ k hk then 1 else 0 := by
  simp only [prodStateXi, dualEquiv_apply]
  refine DimerStateAux.ite_congr_iff (and_congr_right fun hm => ?_) 1 0
  rw [DimerStateAux.odd_eq_iff_pairDefect_eq n (linkOf m q) ξ, hm]
  simp only [DimerStateAux.pairDefect_linkOf]

/-- The sector law of `|n⟩ ⊗ χ_ξ` is the dimer law iff `ξ` is the dimer orientation `ξ_b = N_{≤ b}(n)`
(paper: `prop:dimerlaw`(e)). -/
theorem pairDefect_eq_zero_iff (n : Cfg N) (ξ : Fin N → ZMod 2) :
    (∀ (k : ℕ) (hk : 2 * k + 1 < N), pairDefect n ξ k hk = 0) ↔
      ∀ b : Fin N, b.val % 2 = 1 → ξ b = parityLe n b := by
  rw [DimerStateAux.odd_eq_iff_pairDefect_eq n]
  simp only [DimerStateAux.pairDefect_parityLe]

/-- `G_{2k} G_{2k+1}` acts on `|n⟩ ⊗ χ_ξ` as the number `(-1)^{d_k}` (paper: `prop:dimerlaw`(e)). -/
theorem gauss_pair_prodStateXi (n : Cfg N) (ξ : Fin N → ZMod 2) (k : ℕ) (hk : 2 * k + 1 < N) :
    cplx (gauss (⟨2 * k, by omega⟩ : Fin N) * gauss ⟨2 * k + 1, hk⟩) *ᵥ prodStateXi n ξ =
      ((sgn (pairDefect n ξ k hk) : ℝ) : ℂ) • prodStateXi n ξ := by
  rw [GaussLawAux.gauss_eq, GaussLawAux.gauss_eq, diagonal_mul_diagonal, sec_cplx_diagonal]
  funext p
  rw [mulVec_diagonal, Pi.smul_apply, smul_eq_mul, ← GaussLawAux.sgn_add,
    ← DimerStateAux.pairDefect_eq_chargeOf p.1 p.2 k hk]
  by_cases hp : p.1 = n ∧ ∀ b : Fin N, b.val % 2 = 1 → p.2 b = ξ b
  · rw [hp.1, DimerStateAux.pairDefect_congr n hp.2 k hk]
  · rw [DimerStateAux.prodStateXi_apply_of_not hp, mul_zero, mul_zero]

/-! ### The dimer sector law (paper: `prop:dimerlaw`(a),(b)) -/

/-- **The dimer sector law** (paper: `prop:dimerlaw`(a)): in the dual frame `Ψ^dim(n)` is the uniform
superposition of `|n⟩ ⊗ |q⟩` over the dimer charge configurations. -/
theorem dimerState_dual (n m q : Cfg N) :
    dimerState n (dualEquiv (m, q)) = if m = n ∧ IsDimer q then 1 else 0 := by
  rw [dimerState, prodStateXi_dual]
  refine DimerStateAux.ite_congr_iff (and_congr_right fun _ => ?_) 1 0
  simp only [DimerStateAux.pairDefect_parityLe]
  exact (DimerStateAux.isDimer_iff q).symm

/-- `‖Ψ^dim(n)‖² = 2^{N/2}` (`N` even). -/
theorem normSq_dimerState (hN : Even N) (n : Cfg N) :
    ∑ p, Complex.normSq (dimerState n p) = 2 ^ (N / 2) := by
  rw [← (dualEquiv (N := N)).sum_comp, Fintype.sum_prod_type]
  simp only [dimerState_dual]
  rw [Finset.sum_eq_single n (fun m _ hm => by simp [hm]) (by simp)]
  have e1 : ∀ y : Cfg N, Complex.normSq (if n = n ∧ IsDimer y then (1 : ℂ) else 0) =
      if IsDimer y then 1 else 0 := by
    intro y
    by_cases hy : IsDimer y <;> simp [hy]
  simp only [e1]
  rw [Finset.sum_boole, card_isDimer hN]
  push_cast
  rfl

/-- The sector weights (paper: `prop:dimerlaw`(a)): `‖P_q Ψ^dim(n)‖² = 1` for a dimer sector `q` and `0`
otherwise (unnormalized; divide by `‖Ψ^dim(n)‖² = 2^{N/2}`). -/
theorem normSq_gaussProj_dimerState (n q : Cfg N) :
    ∑ p, Complex.normSq ((cplx (gaussProj q) *ᵥ dimerState n) p) = if IsDimer q then 1 else 0 := by
  rw [gaussProj, sec_cplx_diagonal, ← (dualEquiv (N := N)).sum_comp, Fintype.sum_prod_type]
  simp only [mulVec_diagonal, dimerState_dual]
  simp only [dualEquiv_apply, chargeOf_linkOf]
  rw [Finset.sum_eq_single n (fun m _ hm => by simp [hm]) (by simp)]
  rw [Finset.sum_eq_single q (fun q' _ hq' => by simp [hq']) (by simp)]
  by_cases hq : IsDimer q <;> simp [hq]

/-- The CDW dimer state of `Dimer.lean` is `Ψ^dim` of the density wave on the even sites. -/
theorem dimerState_cdw : dimerState (fockCfg (univ.filter (fun x : Fin N => Even x.val))) = dimerInit N := by
  funext p
  obtain ⟨m, x⟩ := p
  have hlab : ∀ b : Fin N, b.val % 2 = 1 →
      parityLe (fockCfg (univ.filter (fun x : Fin N => Even x.val))) b = dimerLabel b := by
    intro b hb
    rw [DimerStateAux.parityLe_eq_linkOf]
    exact (dimerLinks_linkOf_cdw 0).mpr (fun k _ hk => rfl) b hb
  rw [DimerStateAux.dimerState_apply]
  simp only [dimerInit, dimerGauge]
  rw [DimerStateAux.ite_mul_ite]
  refine DimerStateAux.ite_congr_iff (and_congr_right fun _ => ?_) 1 0
  exact forall_congr' fun b => imp_congr_right fun hb => by rw [hlab b hb]

/-- **Neutral dimers** (paper: `prop:dimerlaw`(b)): `G_{2k} G_{2k+1} Ψ^dim(n) = Ψ^dim(n)`. -/
theorem gauss_pair_dimerState (n : Cfg N) (k : ℕ) (hk : 2 * k + 1 < N) :
    cplx (gauss (⟨2 * k, by omega⟩ : Fin N) * gauss ⟨2 * k + 1, hk⟩) *ᵥ dimerState n =
      dimerState n := by
  rw [dimerState, gauss_pair_prodStateXi, DimerStateAux.pairDefect_parityLe, FreeDynAux.sgn_zero,
    Complex.ofReal_one, one_smul]

/-- (paper: `prop:dimerlaw`(b)): `∏_{j ≤ b} G_j Ψ^dim(n) = Ψ^dim(n)` for every odd link `b`. -/
theorem gauss_prefix_dimerState (n : Cfg N) (b : Fin N) (hb : b.val % 2 = 1) :
    cplx ((((List.finRange N).filter (fun k => decide (k ≤ b))).map gauss).prod) *ᵥ dimerState n =
      dimerState n := by
  rw [GaussLawAux.list_prod_gauss _ ((List.nodup_finRange N).filter _), sec_cplx_diagonal]
  have hfs : ((List.finRange N).filter (fun k => decide (k ≤ b))).toFinset =
      univ.filter (· ≤ b) := by
    ext k
    simp
  rw [hfs]
  funext p
  obtain ⟨m, x⟩ := p
  rw [mulVec_diagonal, DimerStateAux.sum_chargeOf_le]
  by_cases hp : m = n ∧ ∀ c : Fin N, c.val % 2 = 1 → x c = parityLe n c
  · rw [hp.1, hp.2 b hb, sub_self, FreeDynAux.sgn_zero, Complex.ofReal_one, one_mul]
  · rw [DimerStateAux.dimerState_apply, ite_eq_right hp, mul_zero]

/-! ### Dressed fermions (paper: `prop:dimerlaw`(c)) -/

/-- Adding a dressed fermion to the left of all particles: `c̃_s† Ψ^dim(n) = Ψ^dim(n + e_s)` if `n_k = 0` for
all `k ≤ s`. -/
theorem dressed_cre_dimerState (n : Cfg N) (s : Fin N) (hs : ∀ k, k ≤ s → n k = 0) :
    cplx (dressed s)ᵀ *ᵥ dimerState n = dimerState (n + fockCfg {s}) := by
  funext p
  obtain ⟨n', x'⟩ := p
  simp only [mulVec, dotProduct, Fintype.sum_prod_type, DimerStateAux.dressedT_apply]
  rw [Finset.sum_eq_single n (fun m _ hm => Finset.sum_eq_zero fun x _ => by
      rw [DimerStateAux.dimerState_apply_of_ne hm, mul_zero]) (fun h => absurd (Finset.mem_univ n) h)]
  rw [Finset.sum_eq_single (x' + fockCfg (tailOf s)) (fun x _ hx => by
      rw [DimerStateAux.flipSet_apply_of_ne _ hx, mul_zero, Complex.ofReal_zero, zero_mul])
      (fun h => absurd (Finset.mem_univ _) h)]
  rw [DimerStateAux.flipSet_apply_add, mul_one, DimerStateAux.cre_apply_of_empty hs,
    DimerStateAux.dimerState_apply, DimerStateAux.dimerState_apply, DimerStateAux.ofReal_ite_mul_ite]
  refine DimerStateAux.ite_congr_iff
    (and_congr_right fun _ => (and_iff_right rfl).trans ?_) 1 0
  exact forall_congr' fun b => imp_congr_right fun _ => by
    rw [Pi.add_apply, DimerStateAux.parityLe_add_single, DimerStateAux.zmod2_add_eq_iff]

/-- **Dressed fermions** (paper: `prop:dimerlaw`(c)): `Ψ^dim(n) = c̃_{s_1}† ⋯ c̃_{s_m}† Ω^dim` for the
occupied sites `s_1 < ⋯ < s_m`. -/
theorem dimerState_eq_dressed (l : List (Fin N)) (hl : l.Pairwise (· < ·)) :
    (l.map (fun s => cplx (dressed s)ᵀ)).prod *ᵥ dimerState 0 = dimerState (fockCfg l.toFinset) := by
  induction l with
  | nil =>
    rw [List.map_nil, List.prod_nil, one_mulVec, List.toFinset_nil, DualityFullAux.fockCfg_empty]
  | cons s l ih =>
    rw [List.pairwise_cons] at hl
    have hs : ∀ k, k ≤ s → fockCfg l.toFinset k = 0 := by
      intro k hk
      rw [FreeDynAux.fockCfg_eq_zero, List.mem_toFinset]
      intro hkl
      exact absurd hk (not_le.mpr (hl.1 k hkl))
    have hsl : s ∉ l.toFinset := by
      rw [List.mem_toFinset]
      intro hsl
      exact lt_irrefl s (hl.1 s hsl)
    rw [List.map_cons, List.prod_cons, ← mulVec_mulVec, ih hl.2, dressed_cre_dimerState _ s hs,
      DimerStateAux.fockCfg_insert hsl, List.toFinset_cons]

/-- In particular `Ψ^x = c̃_x† Ω^dim` (paper: `prop:dimerlaw`(c)). -/
theorem dimerState_single_eq (x : Fin N) :
    cplx (dressed x)ᵀ *ᵥ dimerState 0 = dimerState (fockCfg {x}) := by
  rw [dressed_cre_dimerState 0 x (fun _ _ => rfl), zero_add]

/-! ### Dynamics (paper: `prop:dimerlaw`(d)) -/

/-- **Sector averaging for dimer states** (paper: `prop:dimerlaw`(d), eq. `eq:dimeravg`). -/
theorem sector_average_dimerState (J h : ℝ) (n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ dimerState n) =
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ Pi.single n 1) := by
  rw [dimer_expect_evo_dual J h _ O hO t, Finset.sum_filter]
  refine Finset.sum_congr rfl fun q _ => ?_
  by_cases hq : IsDimer q
  · rw [ite_eq_left hq]
    congr 2
    ext m
    rw [dimerState_dual, Pi.single_apply]
    simp [hq]
  · rw [ite_eq_right hq]
    have hvec : (fun m => dimerState n (dualEquiv (m, q))) = 0 := by
      ext m
      rw [dimerState_dual]
      simp [hq]
    rw [hvec, mulVec_zero]
    simp [expect]

/-- **The density of the single-fermion dimer state** (paper: eq. `eq:dimerdensity`), normalized:
`⟨Ψ^x(t), n_y Ψ^x(t)⟩ = E_dim |⟨δ_y, e^{-ith_q} δ_x⟩|²`. -/
theorem dimer_density_single (J h : ℝ) (hN : Even N) (x y : Fin N) (t : ℝ) :
    (1 / 2 ^ (N / 2) : ℝ) *
      (expect (cplx (diagonal (fun n : Cfg N => occ (n y)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))
        (evo (H J h) t *ᵥ dimerState (fockCfg {x}))).re =
      dimerAvg (fun q : Cfg N => ‖evo (hq J h q) t y x‖ ^ 2) := by
  unfold dimerAvg
  rw [card_isDimer hN, sector_average_dimerState J h (fockCfg {x}) _ (gaugeInv_matter _) t,
    Complex.re_sum, div_eq_inv_mul, one_div]
  push_cast
  congr 1
  refine Finset.sum_congr rfl fun q _ => ?_
  rw [sectorBlock_matter, Hsec, expect_evo_sub_smul_one]
  have hv : (Pi.single (fockCfg {x}) 1 : Cfg N → ℂ) = fockVec {x} := rfl
  rw [hv, density _ (hq_isSymm J h q), Complex.ofReal_re, Finset.sum_singleton, Complex.sq_norm]

/-! ### Transport of the single-fermion dimer state (paper: Theorem `thm:dimer`) -/

/-- `D^{(N)}_p(t; x) = E_dim ∑_y |y - x|^p |⟨δ_y, e^{-ith_q} δ_x⟩|²` (paper: eq. `eq:dimerdensity`). -/
theorem dimerDLGT_eq (J h : ℝ) (hN : Even N) (p t : ℝ) (x : Fin N) :
    dimerDLGT J h N p t x = dimerD J h N p t x := by
  unfold dimerDLGT dimerD
  have e := DimerLGTAux.dimerAvg_sum_mul (N := N) (fun y : Fin N => |((y.val : ℝ) - x.val)| ^ p)
    (fun y q => ‖evo (hq J h q) t y x‖ ^ 2)
  simp only at e
  rw [e]
  simp only [← dimer_density_single J h hN x, DimerStateAux.expect_matter_re]
  simp only [Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun r _ => ?_
  ring

/-- **Theorem `thm:dimer`(a)** for the single-fermion dimer state: `M^{(N)}_p(T; x) ≥ c_{λ,p} (JT)^{p-1}` for
`1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites). -/
theorem dimer_LGT_moment_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (hN : Even N) (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T) (hT2 : J * T ≤ x.val)
    (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤
      ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * dimerDLGT J h N p t x := by
  have e : (fun t => (1 / T) * Real.exp (-t / T) * dimerDLGT J h N p t x) =
      fun t => (1 / T) * Real.exp (-t / T) * dimerD J h N p t x :=
    funext fun t => by rw [dimerDLGT_eq J h hN p t x]
  rw [e]
  exact dimerMN_lower J h hJ hh0 hh1 x p T hp hT1 hT2 hT3

/-- **Theorem `thm:dimer`(b)**, first limit, for the single-fermion dimer state. -/
theorem dimer_LGT_TL_D (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p t : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) :
    Tendsto (fun n => dimerDLGT J h (Ns n) p t (xs n)) atTop
      (𝓝 (∫ η, ∑' y : ℤ, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 ∂muZ)) := by
  have e : (fun n => dimerDLGT J h (Ns n) p t (xs n)) = fun n => dimerD J h (Ns n) p t (xs n) :=
    funext fun n => dimerDLGT_eq J h (hN n) p t (xs n)
  rw [e]
  exact dimer_TL_D J h hJ r Ns xs hN hpar hx1 hx2 p t hp0 hp2

/-- **Theorem `thm:dimer`(b)**, second limit, for the single-fermion dimer state. -/
theorem dimer_LGT_TL_M (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p T : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) (hT : 0 < T) :
    Tendsto (fun n => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        dimerDLGT J h (Ns n) p t (xs n)) atTop (𝓝 (dimerMInf J h r p T)) := by
  have e : (fun n => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
      dimerDLGT J h (Ns n) p t (xs n)) = fun n => dimerMN J h (Ns n) p T (xs n) := by
    funext n
    simp only [dimerDLGT_eq J h (hN n)]
    rfl
  rw [e]
  exact dimer_TL_M J h hJ r Ns xs hN hpar hx1 hx2 p T hp0 hp2 hT

/-- **Theorem `thm:dimer`(d)**, eq. `eq:dimerLGTd` for `p > 1`, for the single-fermion dimer state:
`sup_{N even, x, t} D^{(N)}_p(t; x) = ∞`. -/
theorem dimer_LGT_unbounded (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerDLGT J h (2 * r.1) p r.2.2 r.2.1)) := by
  have e : (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) => dimerDLGT J h (2 * r.1) p r.2.2 r.2.1) =
      fun r => dimerD J h (2 * r.1) p r.2.2 r.2.1 :=
    funext fun r => dimerDLGT_eq J h (even_two_mul r.1) p r.2.2 r.2.1
  rw [e]
  exact dimerD_unbounded J h hJ hh0 hh1 p hp

end DFL
