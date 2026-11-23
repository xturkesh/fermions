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

/-! ### Other orientations of the odd links (paper: `prop:dimerlaw`(e)) -/

open Classical in
/-- In the dual frame, `|n⟩ ⊗ χ_ξ` is the uniform superposition of `|n⟩ ⊗ |q⟩` over the charges with
`q_{2k} + q_{2k+1} = d_k(n, ξ)` for every pair (paper: `prop:dimerlaw`(e)). -/
theorem prodStateXi_dual (n : Cfg N) (ξ : Fin N → ZMod 2) (m q : Cfg N) :
    prodStateXi n ξ (dualEquiv (m, q)) =
      if m = n ∧ ∀ (k : ℕ) (hk : 2 * k + 1 < N),
          q ⟨2 * k, by omega⟩ + q ⟨2 * k + 1, hk⟩ = pairDefect n ξ k hk then 1 else 0 := by
  sorry

/-- The sector law of `|n⟩ ⊗ χ_ξ` is the dimer law iff `ξ` is the dimer orientation `ξ_b = N_{≤ b}(n)`
(paper: `prop:dimerlaw`(e)). -/
theorem pairDefect_eq_zero_iff (n : Cfg N) (ξ : Fin N → ZMod 2) :
    (∀ (k : ℕ) (hk : 2 * k + 1 < N), pairDefect n ξ k hk = 0) ↔
      ∀ b : Fin N, b.val % 2 = 1 → ξ b = parityLe n b := by
  sorry

/-- `G_{2k} G_{2k+1}` acts on `|n⟩ ⊗ χ_ξ` as the number `(-1)^{d_k}` (paper: `prop:dimerlaw`(e)). -/
theorem gauss_pair_prodStateXi (n : Cfg N) (ξ : Fin N → ZMod 2) (k : ℕ) (hk : 2 * k + 1 < N) :
    cplx (gauss (⟨2 * k, by omega⟩ : Fin N) * gauss ⟨2 * k + 1, hk⟩) *ᵥ prodStateXi n ξ =
      ((sgn (pairDefect n ξ k hk) : ℝ) : ℂ) • prodStateXi n ξ := by
  sorry

/-! ### The dimer sector law (paper: `prop:dimerlaw`(a),(b)) -/

/-- **The dimer sector law** (paper: `prop:dimerlaw`(a)): in the dual frame `Ψ^dim(n)` is the uniform
superposition of `|n⟩ ⊗ |q⟩` over the dimer charge configurations. -/
theorem dimerState_dual (n m q : Cfg N) :
    dimerState n (dualEquiv (m, q)) = if m = n ∧ IsDimer q then 1 else 0 := by
  sorry

/-- `‖Ψ^dim(n)‖² = 2^{N/2}` (`N` even). -/
theorem normSq_dimerState (hN : Even N) (n : Cfg N) :
    ∑ p, Complex.normSq (dimerState n p) = 2 ^ (N / 2) := by
  sorry

/-- The sector weights (paper: `prop:dimerlaw`(a)): `‖P_q Ψ^dim(n)‖² = 1` for a dimer sector `q` and `0`
otherwise (unnormalized; divide by `‖Ψ^dim(n)‖² = 2^{N/2}`). -/
theorem normSq_gaussProj_dimerState (n q : Cfg N) :
    ∑ p, Complex.normSq ((cplx (gaussProj q) *ᵥ dimerState n) p) = if IsDimer q then 1 else 0 := by
  sorry

/-- The CDW dimer state of `Dimer.lean` is `Ψ^dim` of the density wave on the even sites. -/
theorem dimerState_cdw : dimerState (fockCfg (univ.filter (fun x : Fin N => Even x.val))) = dimerInit N := by
  sorry

/-- **Neutral dimers** (paper: `prop:dimerlaw`(b)): `G_{2k} G_{2k+1} Ψ^dim(n) = Ψ^dim(n)`. -/
theorem gauss_pair_dimerState (n : Cfg N) (k : ℕ) (hk : 2 * k + 1 < N) :
    cplx (gauss (⟨2 * k, by omega⟩ : Fin N) * gauss ⟨2 * k + 1, hk⟩) *ᵥ dimerState n =
      dimerState n := by
  sorry

/-- (paper: `prop:dimerlaw`(b)): `∏_{j ≤ b} G_j Ψ^dim(n) = Ψ^dim(n)` for every odd link `b`. -/
theorem gauss_prefix_dimerState (n : Cfg N) (b : Fin N) (hb : b.val % 2 = 1) :
    cplx ((((List.finRange N).filter (fun k => decide (k ≤ b))).map gauss).prod) *ᵥ dimerState n =
      dimerState n := by
  sorry

/-! ### Dressed fermions (paper: `prop:dimerlaw`(c)) -/

/-- Adding a dressed fermion to the left of all particles: `c̃_s† Ψ^dim(n) = Ψ^dim(n + e_s)` if `n_k = 0` for
all `k ≤ s`. -/
theorem dressed_cre_dimerState (n : Cfg N) (s : Fin N) (hs : ∀ k, k ≤ s → n k = 0) :
    cplx (dressed s)ᵀ *ᵥ dimerState n = dimerState (n + fockCfg {s}) := by
  sorry

/-- **Dressed fermions** (paper: `prop:dimerlaw`(c)): `Ψ^dim(n) = c̃_{s_1}† ⋯ c̃_{s_m}† Ω^dim` for the
occupied sites `s_1 < ⋯ < s_m`. -/
theorem dimerState_eq_dressed (l : List (Fin N)) (hl : l.Pairwise (· < ·)) :
    (l.map (fun s => cplx (dressed s)ᵀ)).prod *ᵥ dimerState 0 = dimerState (fockCfg l.toFinset) := by
  sorry

/-- In particular `Ψ^x = c̃_x† Ω^dim` (paper: `prop:dimerlaw`(c)). -/
theorem dimerState_single_eq (x : Fin N) :
    cplx (dressed x)ᵀ *ᵥ dimerState 0 = dimerState (fockCfg {x}) := by
  sorry

/-! ### Dynamics (paper: `prop:dimerlaw`(d)) -/

/-- **Sector averaging for dimer states** (paper: `prop:dimerlaw`(d), eq. `eq:dimeravg`). -/
theorem sector_average_dimerState (J h : ℝ) (n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ dimerState n) =
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ Pi.single n 1) := by
  sorry

/-- **The density of the single-fermion dimer state** (paper: eq. `eq:dimerdensity`), normalized:
`⟨Ψ^x(t), n_y Ψ^x(t)⟩ = E_dim |⟨δ_y, e^{-ith_q} δ_x⟩|²`. -/
theorem dimer_density_single (J h : ℝ) (hN : Even N) (x y : Fin N) (t : ℝ) :
    (1 / 2 ^ (N / 2) : ℝ) *
      (expect (cplx (diagonal (fun n : Cfg N => occ (n y)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))
        (evo (H J h) t *ᵥ dimerState (fockCfg {x}))).re =
      dimerAvg (fun q : Cfg N => ‖evo (hq J h q) t y x‖ ^ 2) := by
  sorry

/-! ### Transport of the single-fermion dimer state (paper: Theorem `thm:dimer`) -/

/-- `D^{(N)}_p(t; x) = E_dim ∑_y |y - x|^p |⟨δ_y, e^{-ith_q} δ_x⟩|²` (paper: eq. `eq:dimerdensity`). -/
theorem dimerDLGT_eq (J h : ℝ) (hN : Even N) (p t : ℝ) (x : Fin N) :
    dimerDLGT J h N p t x = dimerD J h N p t x := by
  sorry

/-- **Theorem `thm:dimer`(a)** for the single-fermion dimer state: `M^{(N)}_p(T; x) ≥ c_{λ,p} (JT)^{p-1}` for
`1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites). -/
theorem dimer_LGT_moment_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (hN : Even N) (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T) (hT2 : J * T ≤ x.val)
    (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤
      ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * dimerDLGT J h N p t x := by
  sorry

/-- **Theorem `thm:dimer`(b)**, first limit, for the single-fermion dimer state. -/
theorem dimer_LGT_TL_D (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p t : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) :
    Tendsto (fun n => dimerDLGT J h (Ns n) p t (xs n)) atTop
      (𝓝 (∫ η, ∑' y : ℤ, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 ∂muZ)) := by
  sorry

/-- **Theorem `thm:dimer`(b)**, second limit, for the single-fermion dimer state. -/
theorem dimer_LGT_TL_M (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p T : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) (hT : 0 < T) :
    Tendsto (fun n => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        dimerDLGT J h (Ns n) p t (xs n)) atTop (𝓝 (dimerMInf J h r p T)) := by
  sorry

/-- **Theorem `thm:dimer`(d)**, eq. `eq:dimerLGTd` for `p > 1`, for the single-fermion dimer state:
`sup_{N even, x, t} D^{(N)}_p(t; x) = ∞`. -/
theorem dimer_LGT_unbounded (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerDLGT J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

end DFL
