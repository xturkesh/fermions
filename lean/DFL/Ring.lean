import DFL.Variants

/-!
# The ring: Wilson loop and twisted boundary condition (paper: Proposition `prop:variants`(ii))

On the ring (0-indexed sites `0, …, N-1`, links `b_j = (j, j+1 mod N)`), the hopping sum includes
the bond `(N-1, 0)` carried by the link `b_{N-1}`, and the Gauss law at site `0` involves
`σˣ_{b_{N-1}}`. In the dual variables `(n, q)` of the open chain (`dualEquiv`), with the same map,
* `G_0^ring ↦ (-1)^{N_f} ∏_{k ≠ 0} τˣ_k`, `G_j^ring ↦ τˣ_j` (`j ≠ 0`), `W = ∏_j σᶻ_{b_j} ↦ τᶻ_0`;
* the ring bond becomes `-J τᶻ_0 (c†_{N-1} c_0 + h.c.)`.
Hence on the joint eigenspace `{G_j^ring = q_j ∀ j, W = w}` the Hamiltonian acts as
`dΓ(h_{q,w}) - h ∑_j q_j` on the matter states of parity `∑_j q_j`, where `h_{q,w}` is `h_q` with the
bond `(N-1, 0)` added with amplitude `-J w`: the twist is exactly the Wilson loop eigenvalue.

Conventions: `ω : ZMod 2` labels `w = sgn ω`; `q : Cfg N` labels `q_j = sgn (q j)`.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ} [NeZero N]

/-- The last site `N - 1`; the ring bond `(N - 1, 0)` is carried by the link `b_{N-1}`. -/
def lastSite : Fin N := ⟨N - 1, Nat.sub_lt (NeZero.pos N) Nat.one_pos⟩

/-- The link to the left of site `j` on the ring: `b_{j-1}`, with `b_{-1} := b_{N-1}`. -/
def prevLinkRing (x : Cfg N) (j : Fin N) : ZMod 2 :=
  if j.val = 0 then x lastSite else prevLink x j

/-- Ring electric term `σˣ_{b_{j-1}} σˣ_{b_j}`. -/
def elecRing (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun x => sgn (prevLinkRing x j + x j))

/-- Ring Gauss operator `G_j^ring = (-1)^{n_j} σˣ_{b_{j-1}} σˣ_{b_j}`. -/
def gaussRing (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun p => sgn (p.1 j) * sgn (prevLinkRing p.2 j + p.2 j))

/-- Hopping across the ring bond `(N - 1, 0)`. -/
def hopWrap : Matrix (Cfg N) (Cfg N) ℝ :=
  cre lastSite * ann 0 + cre 0 * ann lastSite

/-- The ring Hamiltonian
`H^ring = -J ∑_{j} σᶻ_{b_j} (c†_j c_{j+1} + h.c.) - h ∑_j σˣ_{b_{j-1}} σˣ_{b_j}` (indices mod `N`). -/
def Hring (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • (∑ j, hop j ⊗ₖ sz j + hopWrap ⊗ₖ sz lastSite) -
    h • ∑ j, (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elecRing j

/-- `∏_j σᶻ_{b_j}` in the `σˣ` basis: flips every link label. -/
def flipAll : Matrix (Cfg N) (Cfg N) ℝ := fun x' x => if x' = x + 1 then 1 else 0

/-- The Wilson loop `W = ∏_j σᶻ_{b_j}`. -/
def wilson : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipAll

/-- `τᶻ_0` in the dual variables: flips the charge label `q_0`. -/
def flip0 : Matrix (Cfg N) (Cfg N) ℝ :=
  fun q' q => if q' = Function.update q 0 (q 0 + 1) then 1 else 0

/-- The ring Hamiltonian in the dual variables `(n, q)`. -/
def HringDual (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • (∑ j, hop j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) + hopWrap ⊗ₖ flip0) -
    h • (diagonal (fun p : Cfg N × Cfg N => sgn (p.1 0 + ∑ k, p.1 k + ∑ k ∈ univ.erase 0, p.2 k)) +
      ∑ j ∈ univ.erase (0 : Fin N), diagonal (fun p : Cfg N × Cfg N => sgn (p.1 j + p.2 j)))

/-- The adjacency matrix of the ring bond `(N - 1, 0)`. -/
def wrapAdj : Matrix (Fin N) (Fin N) ℝ :=
  fun x y => if (x = lastSite ∧ y = 0) ∨ (x = 0 ∧ y = lastSite) then 1 else 0

/-- The twisted one-body Hamiltonian `h_{q,w} = h_q - J w (|N-1⟩⟨0| + |0⟩⟨N-1|)`, `w = sgn ω`. -/
def hqw (J h : ℝ) (q : Cfg N) (ω : ZMod 2) : Matrix (Fin N) (Fin N) ℝ :=
  hq J h q + (-J * sgn ω) • wrapAdj

/-! ### The dual frame -/

/-- **Ring duality** (paper: Proposition `prop:variants`(ii)): the same map `U` that dualizes the
open chain turns the ring Hamiltonian into `HringDual`. -/
theorem ring_duality (J h : ℝ) :
    (Hring J h).submatrix (dualEquiv (N := N)) dualEquiv = HringDual J h := by
  sorry

/-- `U† G_0^ring U = (-1)^{N_f} ∏_{k ≠ 0} τˣ_k`. -/
theorem gaussRing_dual_zero :
    (gaussRing (0 : Fin N)).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun p : Cfg N × Cfg N => sgn (∑ k, p.1 k + ∑ k ∈ univ.erase 0, p.2 k)) := by
  sorry

/-- `U† G_j^ring U = τˣ_j` for `j ≠ 0`. -/
theorem gaussRing_dual (j : Fin N) (hj : j ≠ 0) :
    (gaussRing j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun p : Cfg N × Cfg N => sgn (p.2 j)) := by
  sorry

/-- `U† W U = τᶻ_0`. -/
theorem wilson_dual :
    (wilson (N := N)).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flip0 := by
  sorry

/-- The Gauss operators commute with the ring Hamiltonian. -/
theorem Hring_comm_gaussRing (J h : ℝ) (j : Fin N) :
    Hring J h * gaussRing j = gaussRing j * Hring J h := by
  sorry

/-- The Wilson loop commutes with the ring Hamiltonian. -/
theorem Hring_comm_wilson (J h : ℝ) :
    Hring J h * wilson (N := N) = wilson (N := N) * Hring J h := by
  sorry

/-- `∏_j G_j^ring = (-1)^{N_f}` (the Gauss operators are diagonal; this is the product of their
entries). -/
theorem gaussRing_prod (p : Cfg N × Cfg N) :
    ∏ j, gaussRing j p p = sgn (∑ k, p.1 k) := by
  sorry

/-! ### Joint eigenvectors and the twisted sector Hamiltonians -/

/-- The joint eigenvector of all `G_j^ring` (eigenvalues `sgn (q j)`) and of `W` (eigenvalue
`sgn ω`) with matter `|n⟩`, for `∑ n = ∑ q`. In the dual frame it is
`|n⟩ ⊗ (|q⁰⟩ + sgn ω |q¹⟩)/√2`, where `q⁰, q¹` agree with `q` off site `0` and have `q_0 = 0, 1`. -/
def ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) : Cfg N × Cfg N → ℂ :=
  fun p => if p.1 = n ∧ ∀ k, k ≠ 0 → chargeOf p.1 p.2 k = q k then
    ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * ((sgn (ω * chargeOf p.1 p.2 0) : ℝ) : ℂ) else 0

/-- The vectors `ringVec q ω n`, `∑ n = ∑ q`, are orthonormal. -/
theorem ringVec_orthonormal (q q' : Cfg N) (ω ω' : ZMod 2) (n n' : Cfg N)
    (hn : ∑ k, n k = ∑ k, q k) (hn' : ∑ k, n' k = ∑ k, q' k) :
    star (ringVec q' ω' n') ⬝ᵥ ringVec q ω n = if q' = q ∧ ω' = ω ∧ n' = n then 1 else 0 := by
  sorry

/-- There are `4^N = dim` labels `(q, ω, n)` with `∑ n = ∑ q`: with orthonormality, the vectors
`ringVec` form an orthonormal basis, so every joint eigenspace `{G^ring = q, W = w}` is non-zero and
spanned by the `ringVec q ω n`. -/
theorem card_ringLabels :
    (univ.filter (fun a : Cfg N × ZMod 2 × Cfg N => ∑ k, a.2.2 k = ∑ k, a.1 k)).card =
      Fintype.card (Cfg N × Cfg N) := by
  sorry

/-- `G_j^ring ringVec = sgn (q j) ringVec`. -/
theorem gaussRing_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) (hn : ∑ k, n k = ∑ k, q k)
    (j : Fin N) :
    cplx (gaussRing j) *ᵥ ringVec q ω n = ((sgn (q j) : ℝ) : ℂ) • ringVec q ω n := by
  sorry

/-- `W ringVec = sgn ω ringVec`. -/
theorem wilson_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) :
    cplx wilson *ᵥ ringVec q ω n = ((sgn ω : ℝ) : ℂ) • ringVec q ω n := by
  sorry

/-- `ringVec q ω n` is a joint eigenvector of the occupations, with eigenvalues `n`. -/
theorem occ_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) (i : Fin N) :
    cplx (diagonal (fun p : Cfg N × Cfg N => occ (p.1 i))) *ᵥ ringVec q ω n =
      ((occ (n i) : ℝ) : ℂ) • ringVec q ω n := by
  sorry

/-- **Twisted sector Hamiltonians** (paper: Proposition `prop:variants`(ii), Eq. `eq:ringtwist`): the
span of the `ringVec q ω n`, `∑ n = ∑ q`, is invariant under `H^ring`, which acts there as
`dΓ(h_{q,w}) - h ∑_j q_j` in the basis `|n⟩`. -/
theorem Hring_ringVec (hN : 2 ≤ N) (J h : ℝ) (q : Cfg N) (ω : ZMod 2) (n : Cfg N)
    (hn : ∑ k, n k = ∑ k, q k) :
    cplx (Hring J h) *ᵥ ringVec q ω n =
      ∑ n' ∈ univ.filter (fun n' : Cfg N => ∑ k, n' k = ∑ k, q k),
        (((dGamma (hqw J h q ω) - (h * ∑ j, sgn (q j)) • 1) n' n : ℝ) : ℂ) • ringVec q ω n' := by
  sorry

/-! ### The `z`-polarized initial state -/

/-- The `z`-polarized gauge state has Wilson loop `W = +1`. -/
theorem wilson_initZ (ψ : Cfg N → ℂ) : cplx wilson *ᵥ initZ ψ = initZ ψ := by
  sorry

/-- **Ring, sector law**: for a matter state of definite parity `π` and `z`-polarized links, the
dynamics of every matter observable diagonal in the occupations is the uniform average over the
`2^{N-1}` charge configurations with `∑_j q_j = π` of free-fermion dynamics with the periodic
one-body Hamiltonians `h_{q,+1}`. -/
theorem sector_average_ring (hN : 2 ≤ N) (J h : ℝ) (ψ : Cfg N → ℂ) (π : ZMod 2)
    (hψ : ∀ n, ψ n ≠ 0 → ∑ j, n j = π) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) (evo (Hring J h) t *ᵥ initZ ψ) =
      (1 / 2 ^ (N - 1) : ℂ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π),
        expect (cplx (diagonal d)) (evo (dGamma (hqw J h q 0)) t *ᵥ ψ) := by
  sorry

end DFL
