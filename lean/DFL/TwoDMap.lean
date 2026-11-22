import DFL.GraphWilson

/-!
# The exact mapping on graphs: joint eigenvectors of the Gauss operators and Wilson loops
(paper: Theorem `thm:2dmap`(b))

For a connected finite graph, charges `q ∈ Cfg N`, link signs `s ∈ Cfg M` and occupations `n` with
`∑ n = ∑ q`, the vector
`Ψ_{n,q,s} = |ker ∂|^{-1/2} ∑_{x : n + ∂x = q} (-1)^{⟨s,x⟩} |n, x⟩`   (`fluxVec q s n`)
is a normalized joint eigenvector of the occupations, of every Gauss operator (eigenvalue `q_i`) and
of every Wilson loop `W_z`, `z ∈ ker ∂` (eigenvalue `(-1)^{⟨s,z⟩}`, i.e. the flux of `s`). Vectors with
different labels `(q, n)` or different fluxes are orthogonal, gauge-equivalent signs give the same
vector up to the sign `q(θ)(-1)^{⟨θ,n⟩}` (paper: `V_{q,θ·s} = q(θ) V_{q,s} Γ(D_θ)`), and there are
`2^N · 2^β · 2^{N-1} = dim` labels, so one representative `s` per flux class gives an orthonormal
basis. In this basis `H` acts as `dΓ(h_{q,s}) - h ∑ q_i` (paper: Eqs. `eq:2ddict`, `eq:2dsector`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The number of cycles, `|ker ∂|`. -/
def cycleCard : ℕ := (univ.filter (fun z : Cfg M => G.bdry z = 0)).card

/-- The flux classes (gauge orbits of sign configurations). -/
def fluxClasses : Finset (Finset (Cfg M)) :=
  univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s))

/-- The joint eigenvector `Ψ_{n,q,s}`. -/
def fluxVec (q : Cfg N) (s : Cfg M) (n : Cfg N) : Cfg N × Cfg M → ℂ :=
  fun p => if p.1 = n ∧ n + G.bdry p.2 = q then
    ((Real.sqrt G.cycleCard : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) else 0

/-- `Ψ_{n,q,s}` is an eigenvector of the occupations with eigenvalues `n`. -/
theorem fluxVec_occ (q : Cfg N) (s : Cfg M) (n : Cfg N) (i : Fin N) :
    cplx (diagonal (fun p : Cfg N × Cfg M => occ (p.1 i))) *ᵥ G.fluxVec q s n =
      ((occ (n i) : ℝ) : ℂ) • G.fluxVec q s n := by
  sorry

/-- `G_i Ψ_{n,q,s} = q_i Ψ_{n,q,s}`. -/
theorem fluxVec_gauss (q : Cfg N) (s : Cfg M) (n : Cfg N) (i : Fin N) :
    cplx (G.gaussOp i) *ᵥ G.fluxVec q s n = ((sgn (q i) : ℝ) : ℂ) • G.fluxVec q s n := by
  sorry

/-- `W_z Ψ_{n,q,s} = (-1)^{⟨s,z⟩} Ψ_{n,q,s}` for every cycle `z`. -/
theorem fluxVec_wilson (q : Cfg N) (s : Cfg M) (n : Cfg N) (z : Cfg M) (hz : G.bdry z = 0) :
    cplx (G.wilsonLoop z) *ᵥ G.fluxVec q s n =
      ((sgn (∑ e, s e * z e) : ℝ) : ℂ) • G.fluxVec q s n := by
  sorry

/-- **Gauge covariance** (paper: `V_{q,θ·s} = q(θ) V_{q,s} Γ(D_θ)`). -/
theorem fluxVec_gaugeAct (q : Cfg N) (s : Cfg M) (n : Cfg N) (θ : Cfg N) :
    G.fluxVec q (G.gaugeAct θ s) n = ((sgn (∑ i, θ i * (q i + n i)) : ℝ) : ℂ) • G.fluxVec q s n := by
  sorry

/-- `Ψ_{n,q,s}` is normalized when `∑ n = ∑ q` (connected graph). -/
theorem fluxVec_norm (hG : G.Connected) (q : Cfg N) (s : Cfg M) (n : Cfg N)
    (hn : ∑ i, n i = ∑ i, q i) :
    Star.star (G.fluxVec q s n) ⬝ᵥ G.fluxVec q s n = 1 := by
  sorry

/-- Vectors with different charges, occupations or fluxes are orthogonal. -/
theorem fluxVec_orthogonal (q q' : Cfg N) (s s' : Cfg M) (n n' : Cfg N)
    (h : q' ≠ q ∨ n' ≠ n ∨ ¬ ∃ θ : Cfg N, s' = G.gaugeAct θ s) :
    Star.star (G.fluxVec q' s' n') ⬝ᵥ G.fluxVec q s n = 0 := by
  sorry

/-- **Completeness**: `2^N` charge configurations, `2^β` flux classes and `2^{N-1}` occupations of the
right parity give `2^{N+M}` labels, the dimension of the Hilbert space (connected graph, `N ≥ 1`). -/
theorem card_fluxLabels (hG : G.Connected) (hN : 0 < N) :
    2 ^ N * G.fluxClasses.card * 2 ^ (N - 1) = Fintype.card (Cfg N × Cfg M) := by
  sorry

/-- **Sector Hamiltonians** (paper: Eqs. `eq:2ddict`, `eq:2dsector`): `H Ψ_{n,q,s}` is the combination
of the `Ψ_{n',q,s}` with the matrix elements of `dΓ(h_{q,s}) - h ∑_i q_i`. -/
theorem ham_fluxVec (J h : ℝ) (q : Cfg N) (s : Cfg M) (n : Cfg N) :
    cplx (G.ham J h) *ᵥ G.fluxVec q s n =
      ∑ n', (((dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) -
        (h * ∑ i, sgn (q i)) • 1) n' n : ℝ) : ℂ) • G.fluxVec q s n' := by
  sorry

end FinGraph

end DFL
