import DFL.Resonance

/-!
# Large deviations and resonance elimination (paper: Proposition `prop:LDT`, Lemma `lem:entry`,
Lemma `lem:double`)

Probability is the uniform measure on the sign configurations `q ∈ Cfg N` (i.i.d. Bernoulli
charges); probabilities are written as cardinalities, `P(S) = #S / 2^N`.

* `LDT λ` is the literature input: a function `γ` (the Lyapunov exponent) bounded below on `I_λ` by
  some `γ_* > 0`, and uniform large deviations for `log ‖Φ_m(E) u‖ / m` around `γ(E)`, for all unit
  vectors `u`, `E ∈ I_λ` and `m ≥ 1`. For the Bernoulli–Anderson model this is
  Bucaj–Damanik–Fillman–Gerbuz–VandenBoom–Wang–Zhang, Prop. 3.6 (vectorwise uniform LDT), together with
  Furstenberg's theorem and the continuity of `γ` (their Thms. 2.1 and 2.6); see the paper,
  Propositions `prop:nocritical` and `prop:LDT`. Only the vector form of the LDT is used.
* `entry_ldt` (paper: Lemma `lem:entry`): matrix entries `⟨w, M e⟩` of products `M` of `m` transfer
  matrices over distinct sites are smaller than `e^{κ m}` with probability `≤ C₁ e^{-c₁ m}`.
* `single_box`, `two_box` (paper: Lemma `lem:double`): a fixed energy is resonant for a box with
  probability `≤ 2 C₁ e^{-c₁ n}`, and two disjoint boxes share a resonant energy with probability
  `≤ 4 (8n+6) C₁ e^{-c₁ n}`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The product of the transfer matrices `T_E(λ q_z)` along the list of sites `zs` (first site
leftmost). -/
def tlist (lam E : ℝ) (q : Cfg N) (zs : List (Fin N)) : Matrix (Fin 2) (Fin 2) ℝ :=
  (zs.map (fun z => transferT E (lam * sgn (q z)))).prod

/-- **Uniform large deviations** (literature input; paper: Propositions `prop:nocritical`, `prop:LDT`,
Eq. `eq:LDTvec`). -/
def LDT (lam : ℝ) : Prop :=
  ∃ γ : ℝ → ℝ, ∃ γs : ℝ, 0 < γs ∧ (∀ E ∈ Ilam lam, γs ≤ γ E) ∧
    ∀ ε : ℝ, 0 < ε → ∃ C c : ℝ, 0 < c ∧ ∀ m : ℕ, 1 ≤ m → ∀ E ∈ Ilam lam, ∀ u : Fin 2 → ℝ,
      u 0 ^ 2 + u 1 ^ 2 = 1 →
        ((univ.filter (fun q : Cfg m =>
          ε ≤ |Real.log (enorm2 (tlist lam E q (List.finRange m) *ᵥ u)) / m - γ E|)).card : ℝ) ≤
          C * Real.exp (-c * m) * 2 ^ m

/-- The conclusion of Lemma `lem:entry` for given constants `κ, C₁, c₁, m₀`. -/
def EntryBound (lam κ C1 c1 : ℝ) (m0 : ℕ) : Prop :=
  ∀ (N : ℕ) (zs : List (Fin N)), zs.Nodup → m0 ≤ zs.length → ∀ E ∈ Ilam lam,
    ∀ e w : Fin 2 → ℝ, e 0 ^ 2 + e 1 ^ 2 = 1 → w 0 ^ 2 + w 1 ^ 2 = 1 →
      ((univ.filter (fun q : Cfg N =>
        |w ⬝ᵥ (tlist lam E q zs *ᵥ e)| ≤ Real.exp (κ * zs.length))).card : ℝ) ≤
        C1 * Real.exp (-c1 * zs.length) * 2 ^ N

/-- **Large deviations for entries in arbitrary directions** (paper: Lemma `lem:entry`). -/
theorem entry_ldt (lam : ℝ) (hL : LDT lam) :
    ∃ κ C1 c1 : ℝ, ∃ m0 : ℕ, 0 < κ ∧ 0 < c1 ∧ 0 ≤ C1 ∧ EntryBound lam κ C1 c1 m0 := by
  sorry

open Classical in
/-- **Single-box estimate** (paper: Lemma `lem:double`, Eq. `eq:single`): for a box
`[a, a+L-1] ⊆ [0, N-1]` and `x = a + i` whose selected ends are at distance `n` from `x`, a fixed
`E ∈ I_λ` is resonant with probability at most `2 C₁ e^{-c₁ n}` (threshold `e^{-κn}`). -/
theorem single_box (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (a : ℤ)
    (L i n : ℕ) (endL endR : Bool) (E : ℝ) (ha : 0 ≤ a) (haL : a + L ≤ N) (hi : i < L)
    (hn : m0 ≤ n) (hR : endR = true → L - 1 - i = n) (hLe : endL = true → i = n)
    (hE : E ∈ Ilam lam) :
    ((univ.filter (fun q : Cfg N =>
      E ∈ resSet (potZ lam q) lam a L i endL endR (Real.exp (-(κ * n))))).card : ℝ) ≤
      2 * C1 * Real.exp (-c1 * n) * 2 ^ N := by
  sorry

open Classical in
/-- **Two-box estimate** (paper: Lemma `lem:double`, Eq. `eq:double`): for two disjoint boxes of
length at most `2n + 1` as in `single_box`, the probability that they have a common resonant energy
is at most `4 (8n+6) C₁ e^{-c₁ n}`. -/
theorem two_box (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (hC1 : 0 ≤ C1)
    (a1 : ℤ) (L1 i1 : ℕ) (endL1 endR1 : Bool) (a2 : ℤ) (L2 i2 : ℕ) (endL2 endR2 : Bool) (n : ℕ)
    (ha1 : 0 ≤ a1) (h12 : a1 + L1 ≤ a2) (ha2 : a2 + L2 ≤ N) (hi1 : i1 < L1) (hi2 : i2 < L2)
    (hL1 : L1 ≤ 2 * n + 1) (hL2 : L2 ≤ 2 * n + 1) (hn : m0 ≤ n)
    (hR1 : endR1 = true → L1 - 1 - i1 = n) (hLe1 : endL1 = true → i1 = n)
    (hR2 : endR2 = true → L2 - 1 - i2 = n) (hLe2 : endL2 = true → i2 = n) :
    ((univ.filter (fun q : Cfg N =>
      (resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 (Real.exp (-(κ * n))) ∩
        resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 (Real.exp (-(κ * n)))).Nonempty)).card : ℝ) ≤
      4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N := by
  sorry

end DFL
