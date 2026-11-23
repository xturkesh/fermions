import DFL.Cluster
import DFL.TwoDLaws
import DFL.TwoDLattices

/-!
# The checkerboard state at strong coupling (paper: Remark `rem:checkerstrong`)

For the checkerboard charges `q⁰ = -ε` of a bipartite graph (sublattice signs `ε`, every edge joins
opposite signs):

* every equal-charge cluster is a single site (`checker_inCluster_iff`), `h⁰_{q⁰,s} = 2h diag(q⁰)`
  (`hqs0_checker`), and Proposition `prop:cluster`(b) gives the deterministic return amplitude
  `|⟨δ_x, e^{-ith_{q⁰,s}} δ_x⟩| ≥ 1 - 8J/g - 16|t|J²/g` in every flux sector (`checker_return`);
* `R_s = (K_s² + 4h²)^{1/2}` (`checkerR`, continuous functional calculus of the Hermitian matrix
  `K_s`) satisfies `R_s² = K_s² + 4h²` (`checkerR_sq`) and, by the scalar inequality
  `0 ≤ 1 + y/2 - √(1+y) ≤ y²/8` (`one_add_half_sub_sqrt_bounds`) and `‖K_s‖ ≤ 4J`
  (`norm_cplx_Ks_le`),
  `‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ ‖K_s‖⁴/(64|h|³) ≤ 4J⁴/|h|³` (`norm_checkerR_sub_le`,
  `norm_checkerR_sub_le_J`); hence `‖e^{-itR_s} - e^{-it(2|h| + K_s²/(4|h|))}‖ ≤ |t| ‖K_s‖⁴/(64|h|³)`
  (`norm_exp_checkerR_sub_le`, `norm_exp_checkerR_sub_le_J`);
* `K_s²` connects only sites of the same sublattice (`Ks_sq_sublattice`); between the two ends of a
  diagonal of a plaquette it equals `J² s s (1 + W_p)` (`Ks_sq_plaquette`), between the ends of a
  straight path of length two `J² s s` (`Ks_sq_line`), and it vanishes between distinct sites without
  a common neighbour (`Ks_sq_far`). On the planar box (`box_simple`, `box_deg_le_four`) this gives
  the effective sublattice hopping of the remark: `±2J²` (resp. `0`) on the plaquette diagonals if
  `W_p = +1` (resp. `-1`), and `±J²` on the straight next-nearest-neighbour bonds
  (`box_Ks_sq_diag`, `box_Ks_sq_antidiag`, `box_Ks_sq_horiz`, `box_Ks_sq_vert`).
-/

open Matrix Finset NormedSpace
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-! ### Clusters and the return amplitude -/

/-- For the checkerboard charges `q⁰ = -ε` every equal-charge cluster is a single site. -/
theorem checker_inCluster_iff (ε : Fin N → ℝ) (hε1 : ∀ i, ε i = 1 ∨ ε i = -1)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (x y : Fin N) :
    G.InCluster (fun i => -ε i) x y ↔ x = y := by
  sorry

/-- `h⁰_{q⁰,s} = 2h diag(q⁰)` for the checkerboard charges. -/
theorem hqs0_checker (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.hqs0 J h (fun i => -ε i) s = (2 * h) • diagonal (fun i => -ε i) := by
  sorry

/-- **Return amplitude of the checkerboard state** (paper: Remark `rem:checkerstrong`): in every
flux sector, `|⟨δ_x, e^{-ith_{q⁰,s}} δ_x⟩| ≥ 1 - 8J/g - 16|t|J²/g`, `g = 4|h| - 8J`. -/
theorem checker_return (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (ε : Fin N → ℝ) (hε1 : ∀ i, ε i = 1 ∨ ε i = -1)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) :
    1 - 8 * J / (4 * |h| - 8 * J) - 16 * |t| * J ^ 2 / (4 * |h| - 8 * J) ≤
      ‖evo (G.hqs J h (fun i => -ε i) s) t x x‖ := by
  sorry

/-! ### The slow generator -/

/-- The scalar inequality `0 ≤ 1 + y/2 - √(1+y) ≤ y²/8` for `y ≥ 0`. -/
theorem one_add_half_sub_sqrt_bounds (y : ℝ) (hy : 0 ≤ y) :
    0 ≤ 1 + y / 2 - √(1 + y) ∧ 1 + y / 2 - √(1 + y) ≤ y ^ 2 / 8 := by
  sorry

/-- `‖K_s‖ ≤ 4J` on graphs of degree at most four (Schur test). -/
theorem norm_cplx_Ks_le (J : ℝ) (hJ : 0 ≤ J) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) :
    ‖cplx (G.Ks J s)‖ ≤ 4 * J := by
  sorry

/-- `R_s = (K_s² + 4h²)^{1/2}`, by the continuous functional calculus of the Hermitian matrix
`K_s` (complexified). -/
def checkerR (J h : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℂ :=
  cfc (fun x : ℝ => √(x ^ 2 + 4 * h ^ 2)) (cplx (G.Ks J s))

/-- The effective slow generator `2|h| + K_s²/(4|h|)`. -/
def checkerEff (J h : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  (2 * |h|) • (1 : Matrix (Fin N) (Fin N) ℝ) + (1 / (4 * |h|)) • (G.Ks J s * G.Ks J s)

/-- `R_s` is Hermitian. -/
theorem checkerR_isHermitian (J h : ℝ) (s : Fin M → ℝ) : (G.checkerR J h s).IsHermitian := by
  sorry

/-- `R_s² = K_s² + 4h²`. -/
theorem checkerR_sq (J h : ℝ) (s : Fin M → ℝ) :
    G.checkerR J h s * G.checkerR J h s =
      cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) := by
  sorry

/-- **Slow generator** (paper: Remark `rem:checkerstrong`):
`‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ ‖K_s‖⁴/(64|h|³)`. -/
theorem norm_checkerR_sub_le (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) :
    ‖G.checkerR J h s - cplx (G.checkerEff J h s)‖ ≤ ‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3) := by
  sorry

/-- `‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ 4J⁴/|h|³` on graphs of degree at most four. -/
theorem norm_checkerR_sub_le_J (J h : ℝ) (hJ : 0 ≤ J) (hh : h ≠ 0) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) :
    ‖G.checkerR J h s - cplx (G.checkerEff J h s)‖ ≤ 4 * J ^ 4 / |h| ^ 3 := by
  sorry

/-- **The slow dynamics** (paper: Remark `rem:checkerstrong`): up to `|t| ‖K_s‖⁴/(64|h|³)`, the
dynamics generated by `R_s` is that generated by `2|h| + K_s²/(4|h|)`. -/
theorem norm_exp_checkerR_sub_le (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) (t : ℝ) :
    ‖exp ((-(t : ℂ) * Complex.I) • G.checkerR J h s) - evo (G.checkerEff J h s) t‖ ≤
      |t| * (‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3)) := by
  sorry

/-- The slow dynamics up to `4|t|J⁴/|h|³` on graphs of degree at most four. -/
theorem norm_exp_checkerR_sub_le_J (J h : ℝ) (hJ : 0 ≤ J) (hh : h ≠ 0) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ) :
    ‖exp ((-(t : ℂ) * Complex.I) • G.checkerR J h s) - evo (G.checkerEff J h s) t‖ ≤
      |t| * (4 * J ^ 4 / |h| ^ 3) := by
  sorry

/-! ### Structure of `K_s²` -/

/-- `K_s²` connects only sites of the same sublattice. -/
theorem Ks_sq_sublattice (J : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (i k : Fin N) (hik : ε i ≠ ε k) :
    (G.Ks J s * G.Ks J s) i k = 0 := by
  sorry

/-- `K_s²` vanishes between sites without a common neighbour. -/
theorem Ks_sq_far (J : ℝ) (s : Fin M → ℝ) (i k : Fin N)
    (hfar : ∀ j, ¬ ((∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) ∧
      (∃ e, (G.src e = j ∧ G.tgt e = k) ∨ (G.tgt e = j ∧ G.src e = k)))) :
    (G.Ks J s * G.Ks J s) i k = 0 := by
  sorry

/-- **Plaquette diagonal**: on a simple graph, if `j₁ ≠ j₂` are the only common neighbours of `i` and
`k`, joined by the edges `e₁ = {i, j₁}`, `e₂ = {j₁, k}`, `e₃ = {k, j₂}`, `e₄ = {j₂, i}`, then
`(K_s²)_{ik} = J² s_{e₁} s_{e₂} (1 + W)` with the plaquette flux `W = s_{e₁} s_{e₂} s_{e₃} s_{e₄}`. -/
theorem Ks_sq_plaquette (J : ℝ) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i k j₁ j₂ : Fin N) (hj : j₁ ≠ j₂) (e₁ e₂ e₃ e₄ : Fin M)
    (h₁ : (G.src e₁ = i ∧ G.tgt e₁ = j₁) ∨ (G.tgt e₁ = i ∧ G.src e₁ = j₁))
    (h₂ : (G.src e₂ = j₁ ∧ G.tgt e₂ = k) ∨ (G.tgt e₂ = j₁ ∧ G.src e₂ = k))
    (h₃ : (G.src e₃ = k ∧ G.tgt e₃ = j₂) ∨ (G.tgt e₃ = k ∧ G.src e₃ = j₂))
    (h₄ : (G.src e₄ = j₂ ∧ G.tgt e₄ = i) ∨ (G.tgt e₄ = j₂ ∧ G.src e₄ = i))
    (hcommon : ∀ j, (∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) →
      (∃ e, (G.src e = j ∧ G.tgt e = k) ∨ (G.tgt e = j ∧ G.src e = k)) → j = j₁ ∨ j = j₂) :
    (G.Ks J s * G.Ks J s) i k = J ^ 2 * s e₁ * s e₂ * (1 + s e₁ * s e₂ * s e₃ * s e₄) := by
  sorry

/-- **Straight path of length two**: on a simple graph, if `j` is the only common neighbour of `i`
and `k`, joined by the edges `e₁ = {i, j}` and `e₂ = {j, k}`, then `(K_s²)_{ik} = J² s_{e₁} s_{e₂}`. -/
theorem Ks_sq_line (J : ℝ) (s : Fin M → ℝ)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i k j : Fin N) (e₁ e₂ : Fin M)
    (h₁ : (G.src e₁ = i ∧ G.tgt e₁ = j) ∨ (G.tgt e₁ = i ∧ G.src e₁ = j))
    (h₂ : (G.src e₂ = j ∧ G.tgt e₂ = k) ∨ (G.tgt e₂ = j ∧ G.src e₂ = k))
    (hcommon : ∀ j', (∃ e, (G.src e = i ∧ G.tgt e = j') ∨ (G.tgt e = i ∧ G.src e = j')) →
      (∃ e, (G.src e = j' ∧ G.tgt e = k) ∨ (G.tgt e = j' ∧ G.src e = k)) → j' = j) :
    (G.Ks J s * G.Ks J s) i k = J ^ 2 * s e₁ * s e₂ := by
  sorry

/-! ### The planar box -/

open TwoDLatticesAux

/-- The corner `(x + a, y + b)` of the box plaquette `p = (x, y)` (`a, b ∈ {0, 1}`). -/
def plaqCorner (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) (a b : Fin 2) : Fin (L₁ * L₂) :=
  boxVert ⟨p.1.val + a.val, by have := p.1.isLt; have := a.isLt; omega⟩
    ⟨p.2.val + b.val, by have := p.2.isLt; have := b.isLt; omega⟩

/-- The plaquette flux `W_p = ∏_{e ∈ ∂p} s_e` of a box plaquette. -/
def boxFlux (L₁ L₂ : ℕ) (s : Fin (boxM L₁ L₂) → ℝ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : ℝ :=
  s (boxBot L₁ L₂ p) * s (boxRight L₁ L₂ p) * s (boxTop L₁ L₂ p) * s (boxLeft L₁ L₂ p)

/-- The box is a simple graph. -/
theorem box_simple (L₁ L₂ : ℕ) (e e' : Fin (boxM L₁ L₂))
    (h : ((box L₁ L₂).src e = (box L₁ L₂).src e' ∧ (box L₁ L₂).tgt e = (box L₁ L₂).tgt e') ∨
      ((box L₁ L₂).src e = (box L₁ L₂).tgt e' ∧ (box L₁ L₂).tgt e = (box L₁ L₂).src e')) :
    e = e' := by
  sorry

/-- Every vertex of the box has degree at most four. -/
theorem box_deg_le_four (L₁ L₂ : ℕ) (i : Fin (L₁ * L₂)) :
    (univ.filter (fun e => (box L₁ L₂).src e = i ∨ (box L₁ L₂).tgt e = i)).card ≤ 4 := by
  sorry

/-- **Plaquette diagonals of the box** (paper: Remark `rem:checkerstrong`):
`(K_s²)_{(x,y),(x+1,y+1)} = J² s_bot s_right (1 + W_p)`. -/
theorem box_Ks_sq_diag (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s) (plaqCorner L₁ L₂ p 0 0) (plaqCorner L₁ L₂ p 1 1) =
      J ^ 2 * s (boxBot L₁ L₂ p) * s (boxRight L₁ L₂ p) * (1 + boxFlux L₁ L₂ s p) := by
  sorry

/-- The other diagonal: `(K_s²)_{(x+1,y),(x,y+1)} = J² s_bot s_left (1 + W_p)`. -/
theorem box_Ks_sq_antidiag (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s) (plaqCorner L₁ L₂ p 1 0) (plaqCorner L₁ L₂ p 0 1) =
      J ^ 2 * s (boxBot L₁ L₂ p) * s (boxLeft L₁ L₂ p) * (1 + boxFlux L₁ L₂ s p) := by
  sorry

/-- **Straight next-nearest-neighbour bonds of the box**, horizontal:
`(K_s²)_{(x,y),(x+2,y)} = J² s_{(x,y)—(x+1,y)} s_{(x+1,y)—(x+2,y)}`. -/
theorem box_Ks_sq_horiz (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (x : Fin (L₁ - 2)) (y : Fin L₂) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s)
        (boxVert ⟨x.val, by have := x.isLt; omega⟩ y)
        (boxVert ⟨x.val + 2, by have := x.isLt; omega⟩ y) =
      J ^ 2 * s (boxEdge L₁ L₂ (Sum.inl (⟨x.val, by have := x.isLt; omega⟩, y))) *
        s (boxEdge L₁ L₂ (Sum.inl (⟨x.val + 1, by have := x.isLt; omega⟩, y))) := by
  sorry

/-- Vertical: `(K_s²)_{(x,y),(x,y+2)} = J² s_{(x,y)—(x,y+1)} s_{(x,y+1)—(x,y+2)}`. -/
theorem box_Ks_sq_vert (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (x : Fin L₁) (y : Fin (L₂ - 2)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s)
        (boxVert x ⟨y.val, by have := y.isLt; omega⟩)
        (boxVert x ⟨y.val + 2, by have := y.isLt; omega⟩) =
      J ^ 2 * s (boxEdge L₁ L₂ (Sum.inr (x, ⟨y.val, by have := y.isLt; omega⟩))) *
        s (boxEdge L₁ L₂ (Sum.inr (x, ⟨y.val + 1, by have := y.isLt; omega⟩))) := by
  sorry

end FinGraph

end DFL
