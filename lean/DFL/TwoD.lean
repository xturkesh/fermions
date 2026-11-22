import DFL.Graph

/-!
# Two dimensions: no compactly supported eigenvectors; the checkerboard chiral gap
(paper: Proposition `prop:nocompact`; Proposition `prop:2dlaws`(c)(i),(ii), and the gauge covariance
`D_θ h_{q,s} D_θ = h_{q,θ·s}` below Eq. `eq:hqs`)

* `no_compact_eigenvector`: a nearest-neighbour operator on `ℤ²` with nonvanishing hoppings has no
  finitely supported eigenvector (look at the lexicographically extremal point of the support on its
  top antidiagonal). In particular no `h_{q,s}` on `ℓ²(ℤ²)` has one (`no_compact_eigenvector_hqs`).
* On a finite graph with link signs `s` (binary fluxes), `K_s = -J ∑_e s_e (|i⟩⟨j| + |j⟩⟨i|)` and
  `h_{q,s} = K_s + 2h diag(q)`. For the checkerboard state on a bipartite graph the charges are
  `q⁰ = -ε` (sublattice signs), `h_{q⁰,s} = K_s - 2h Γ_ε` with `K_s Γ_ε = -Γ_ε K_s`, hence
  `h_{q⁰,s}² = K_s² + 4h²`: every eigenvalue has modulus at least `2|h|`, and the eigenspace of
  `h²` at `4h²` (the band edges `±2|h|`) is `ker K_s`.
-/

open Matrix Finset

noncomputable section

namespace DFL

/-! ### No compactly supported eigenvectors on `ℤ²` -/

/-- Nearest neighbours in `ℤ²`. -/
def IsNN (x y : ℤ × ℤ) : Prop := (x.1 - y.1).natAbs + (x.2 - y.2).natAbs = 1

/-- **No compactly supported eigenvectors on `ℤ²`** (paper: Proposition `prop:nocompact`): if
`(Hu)(x) = ∑_{y ∼ x} a_{xy} u(y) + V(x) u(x)` with all nearest-neighbour `a_{xy} ≠ 0`, a finitely
supported solution of `Hu = Eu` vanishes. -/
theorem no_compact_eigenvector (a : ℤ × ℤ → ℤ × ℤ → ℂ) (V : ℤ × ℤ → ℂ) (E : ℂ)
    (ha : ∀ x y, IsNN x y → a x y ≠ 0) (u : ℤ × ℤ → ℂ) (hu : (Function.support u).Finite)
    (heig : ∀ x, a x (x + (1, 0)) * u (x + (1, 0)) + a x (x - (1, 0)) * u (x - (1, 0)) +
        a x (x + (0, 1)) * u (x + (0, 1)) + a x (x - (0, 1)) * u (x - (0, 1)) + V x * u x =
      E * u x) :
    u = 0 := by
  sorry

/-- No `h_{q,s} = K_s + 2h diag(q)` on `ℓ²(ℤ²)` (signs `s_{xy} = ±1`, `J ≠ 0`, any potential `q`) has
a finitely supported eigenvector. -/
theorem no_compact_eigenvector_hqs (J h : ℝ) (hJ : J ≠ 0) (s : ℤ × ℤ → ℤ × ℤ → ℝ)
    (hs : ∀ x y, s x y = 1 ∨ s x y = -1) (q : ℤ × ℤ → ℝ) (E : ℂ) (u : ℤ × ℤ → ℂ)
    (hu : (Function.support u).Finite)
    (heig : ∀ x, ((-J : ℝ) : ℂ) * ((s x (x + (1, 0)) : ℂ) * u (x + (1, 0)) +
        (s x (x - (1, 0)) : ℂ) * u (x - (1, 0)) + (s x (x + (0, 1)) : ℂ) * u (x + (0, 1)) +
        (s x (x - (0, 1)) : ℂ) * u (x - (0, 1))) + ((2 * h * q x : ℝ) : ℂ) * u x = E * u x) :
    u = 0 := by
  sorry

/-! ### Binary fluxes on a finite graph; the checkerboard -/

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The hopping operator with link signs `s`: `K_s = -J ∑_e s_e (|s(e)⟩⟨t(e)| + |t(e)⟩⟨s(e)|)`. -/
def Ks (J : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  fun i j => -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
    (G.tgt e = i ∧ G.src e = j)), s e

/-- `h_{q,s} = K_s + 2h diag(q)`. -/
def hqs (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  G.Ks J s + (2 * h) • diagonal q

/-- With all signs `+1`, `h_{q,+}` is the Anderson–Bernoulli operator `h_q` of the trivial flux
sector. -/
theorem hqs_one (J h : ℝ) (q : Cfg N) :
    G.hqs J h (fun i => sgn (q i)) (fun _ => 1) = G.hq J h q := by
  sorry

/-- **Gauge covariance**: `D_θ h_{q,s} D_θ = h_{q,θ·s}`, `(θ·s)_e = θ_{s(e)} s_e θ_{t(e)}`, for
`θ ∈ {±1}^V`. -/
theorem gauge_hqs (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (θ : Fin N → ℝ)
    (hθ : ∀ i, θ i = 1 ∨ θ i = -1) :
    diagonal θ * G.hqs J h q s * diagonal θ =
      G.hqs J h q (fun e => θ (G.src e) * s e * θ (G.tgt e)) := by
  sorry

/-- On a bipartite graph (sublattice signs `ε`, every edge joins opposite signs), `K_s`
anticommutes with `Γ_ε = diag(ε)`. -/
theorem Ks_anticomm (J : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.Ks J s * diagonal ε + diagonal ε * G.Ks J s = 0 := by
  sorry

/-- **Checkerboard identity** (paper: Proposition `prop:2dlaws`(c)(i)): for the checkerboard charges
`q⁰ = -ε`, `h_{q⁰,s}² = K_s² + 4h²`. -/
theorem checkerboard_sq (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.hqs J h (fun i => -ε i) s * G.hqs J h (fun i => -ε i) s =
      G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ) := by
  sorry

/-- **Chiral gap** (paper: Proposition `prop:2dlaws`(c)(ii)): every eigenvalue `E` of the checkerboard
operator `h_{q⁰,s}` (real signs `s`) has `|E| ≥ 2|h|`. -/
theorem checkerboard_gap (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (E : ℂ) (v : Fin N → ℂ) (hv : v ≠ 0)
    (heig : (G.hqs J h (fun i => -ε i) s).map (fun r => (r : ℂ)) *ᵥ v = E • v) :
    2 * |h| ≤ ‖E‖ := by
  sorry

/-- **Band edges = chiral zero modes** (paper: Proposition `prop:2dlaws`(c)(ii)): the eigenspace of
`h_{q⁰,s}²` at `4h²` (the spectral subspace of `|h_{q⁰,s}|` at `2|h|`) is `ker K_s`. -/
theorem checkerboard_edge_iff (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (v : Fin N → ℝ) :
    G.Ks J s *ᵥ v = 0 ↔
      (G.hqs J h (fun i => -ε i) s * G.hqs J h (fun i => -ε i) s) *ᵥ v = (4 * h ^ 2) • v := by
  sorry

end FinGraph

end DFL
