import DFL.TwoD

/-!
# Strong coupling: cluster confinement up to times of order `|h|/J²`
(paper: Proposition `prop:cluster`)

On a finite graph with at most four edges at every vertex (e.g. a box in `ℤ²`), with link signs
`s = ±1` and charges `q = ±1`, let `h = h_{q,s} = K_s + 2h diag(q)` and let `h⁰` be `h` with the
hopping across all bonds joining opposite charges removed (`hqs0`, the pinching of `h` by the
projections `P_± = 1_{q = ±1}`). The equal-charge clusters are the connected components of the
subgraphs induced on `{q = +1}` and `{q = -1}`.

* `hqs0_cluster` (paper: (a)): `h⁰` leaves `ℓ²(C)` invariant for every equal-charge cluster `C`.
* `cluster_bound` (paper: (b)): for `|h| > 2J` and `g = 4|h| - 8J`,
  `‖e^{-ith} - e^{-ith⁰}‖ ≤ 8J/g + 48|t|J²/g` (operator norm), uniformly in the graph, `q` and `s`.
* `cluster_leak` (paper: (b)): the weight of `e^{-ith} δ_x` outside the cluster of `x` is at most
  `(8J/g + 48|t|J²/g)²`.
* `cluster_avg` (paper: (c), up to the percolation input): for any law of `q` with density at most
  `K` with respect to the uniform law, the averaged weight outside a region `B ∋ x` is at most
  `(8J/g + 48|t|J²/g)² + K · P_unif(C_q(x) ⊄ B)`. The exponential bound on `P_unif(C_q(x) ⊄ B_R(x))`
  for i.i.d. charges on `ℤ²` is the percolation estimate of Li–Liang (paper: proof of (c)); it is
  not formalized here.
-/

open Matrix Finset
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- `h⁰_{q,s}`: the hopping across bonds joining opposite charges removed. -/
def hqs0 (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  fun i j => if q i = q j then G.hqs J h q s i j else 0

/-- Adjacent vertices of equal charge. -/
def eqChargeAdj (q : Fin N → ℝ) (i j : Fin N) : Prop :=
  (∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) ∧ q i = q j

/-- `y` lies in the equal-charge cluster `C_q(x)` of `x`. -/
def InCluster (q : Fin N → ℝ) (x y : Fin N) : Prop :=
  Relation.ReflTransGen (G.eqChargeAdj q) x y

/-- **Clusters are invariant** (paper: Proposition `prop:cluster`(a)): if `z ∈ C_q(x)` and
`h⁰_{yz} ≠ 0`, then `y ∈ C_q(x)`. -/
theorem hqs0_cluster (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (x y z : Fin N)
    (hz : G.InCluster q x z) (hyz : G.hqs0 J h q s y z ≠ 0) : G.InCluster q x y := by
  sorry

/-- **Cluster confinement, operator-norm bound** (paper: Proposition `prop:cluster`(b)). -/
theorem cluster_bound (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ) :
    ‖evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t‖ ≤
      8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J) := by
  sorry

open Classical in
/-- **Leakage out of the cluster** (paper: Proposition `prop:cluster`(b)): the weight of
`e^{-ith} δ_x` outside `C_q(x)` is at most `(8J/g + 48|t|J²/g)²`. -/
theorem cluster_leak (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) :
    ∑ y ∈ univ.filter (fun y => ¬ G.InCluster q x y), ‖evo (G.hqs J h q s) t y x‖ ^ 2 ≤
      (8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 := by
  sorry

open Classical in
/-- **Averaged confinement** (paper: Proposition `prop:cluster`(c), with the percolation estimate
left as the input `P_unif(C_q(x) ⊄ B)`): charges `q ∈ {±1}^V` encoded as `sgn ∘ q`, `q : Cfg N`,
distributed with weights `w ≤ K · 2^{-N}`. -/
theorem cluster_avg (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) (B : Finset (Fin N)) (w : Cfg N → ℝ) (K : ℝ) (hw0 : ∀ q, 0 ≤ w q)
    (hw1 : ∑ q, w q = 1) (hwK : ∀ q, w q ≤ K / 2 ^ N) :
    ∑ q, w q * ∑ y ∈ univ.filter (fun y => y ∉ B),
        ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2 ≤
      (8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 +
        K * ((univ.filter (fun q : Cfg N =>
          ∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y)).card : ℝ) / 2 ^ N := by
  sorry

end FinGraph

end DFL
