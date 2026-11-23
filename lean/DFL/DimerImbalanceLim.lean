import DFL.DimerImbalance

/-!
# The dimer imbalance vanishes in the thermodynamic limit (paper: Lemma `lem:fewapprox`,
Proposition `prop:dimerimbalance`(iii))

* `card_eig_near_le` (paper: Lemma `lem:fewapprox`): for every charge configuration `q` (no average), `h_q` has
  at most `2 (N/L + 1)` eigenvalues in `[E - ε, E + ε]` whenever `ε L (K + 1)^L < |J|`, with
  `K = (2|h| + |E|)/|J|`. A vector in the span of the corresponding eigenvectors that vanishes on the first
  two sites of every block of length `L` is `0`: on each block the eigenvalue recursion propagates the two
  zeros with errors controlled by `(K + 1)^L ‖(h_q - E) v‖ / |J|`.
* `dimerImbInf_le`: `Ī_N ≤ 8 L² (4|h|/|J| + 1)^{2L} / N + 4/L + 4/N` for every `L ≥ 1`, from the dimer identity
  of `prop:dimerimbalance`(i) (`dimer_stagWeight_identity`): eigenvectors away from the energies `±2h` have
  small staggered weight, and there are few eigenvalues near `±2h`.
* `dimerImbInf_tendsto_zero` (paper: `prop:dimerimbalance`(iii)): `Ī_N → 0` as `N → ∞` through even values.
  (Choosing `L ≍ log N` gives `Ī_N = O(1 / log N)`.)

This replaces the paper's earlier route through the density of states of the infinite-volume dimer model.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **Few eigenvalues near any energy** (paper: Lemma `lem:fewapprox`): if `ε L (K + 1)^L < |J|` with
`K = (2|h| + |E|)/|J|`, then `h_q` has at most `2 (N/L + 1)` eigenvalues in `[E - ε, E + ε]`. -/
theorem card_eig_near_le (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (E ε : ℝ) (L : ℕ) (hL : 0 < L)
    (hε : ε * L * ((2 * |h| + |E|) / |J| + 1) ^ L < |J|) :
    ((univ.filter (fun α : Fin N => |(hq_isHermitian J h q).eigenvalues α - E| ≤ ε)).card : ℝ) ≤
      2 * ((N : ℝ) / L + 1) := by
  sorry

/-- **Quantitative vanishing of the dimer imbalance**: for every `L ≥ 1`,
`Ī_N ≤ 8 L² (4|h|/|J| + 1)^{2L} / N + 4/L + 4/N`. -/
theorem dimerImbInf_le (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (hN0 : 0 < N) (L : ℕ) (hL : 0 < L) :
    (1 / 2 ^ (N / 2) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbInf J h q ≤
      8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / N + 4 / L + 4 / N := by
  sorry

/-- **The imbalance of the dimer CDW state vanishes** (paper: Proposition `prop:dimerimbalance`(iii)):
`Ī_N → 0` as `N → ∞` through even values. -/
theorem dimerImbInf_tendsto_zero (J h : ℝ) (hJ : J ≠ 0) :
    Tendsto (fun M : ℕ => (1 / 2 ^ (2 * M / 2) : ℝ) *
      ∑ q ∈ univ.filter (fun q : Cfg (2 * M) => IsDimer q), imbInf J h q) atTop (𝓝 0) := by
  sorry

end DFL
