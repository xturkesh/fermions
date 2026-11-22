import Mathlib

/-!
# Disorder-free localization in the 1D ℤ₂ lattice gauge theory — basic definitions

This file fixes the objects shared by the whole formalization. It contains **definitions only**;
statements and proofs live in the other files of the `DFL` library.

## Conventions (0-indexed version of `paper/BRIEF.md`, Section 2)

* Sites `j : Fin N` carry spinless fermions; links `b_j : Fin N` join `j` and `j+1`; the last link
  `b_{N-1}` is a *boundary link* (no hopping crosses it); the fictitious link `b_{-1}` has
  `σˣ_{b_{-1}} := 1`.
* All two-level labels live in `ZMod 2`:
  * fermion occupations `n_j ∈ {0,1}`;
  * link labels `x_b`, where `x_b = 0 ↔ σˣ_b = +1` and `x_b = 1 ↔ σˣ_b = -1` (we work in the
    `σˣ` eigenbasis of the links);
  * Gauss charges `q_j`, where `q_j = 0 ↔ G_j = +1` and `q_j = 1 ↔ G_j = -1`.
* The LGT Hilbert space is `Cfg N × Cfg N → ℂ` (occupations × link labels). The dual space is
  `Cfg N × Cfg N → ℂ` (occupations × charges). Operators are matrices indexed by configurations;
  all Hamiltonians are real.
* Fermions use the Jordan–Wigner convention `c_x |n⟩ = (-1)^{∑_{k<x} n_k} |n - e_x⟩`.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

/-! ### Signs and configurations -/

/-- The sign `(-1)^a` of `a : ZMod 2`. -/
def sgn (a : ZMod 2) : ℝ := if a = 0 then 1 else -1

/-- The occupation number of `a : ZMod 2` as a real number (`0 ↦ 0`, `1 ↦ 1`). -/
def occ (a : ZMod 2) : ℝ := if a = 0 then 0 else 1

/-- Configurations of `N` two-level systems (occupations, link labels or charges). -/
abbrev Cfg (N : ℕ) := Fin N → ZMod 2

variable {N : ℕ}

/-! ### Fermions -/

/-- Jordan–Wigner parity `∑_{k<x} n_k` of the modes strictly to the left of `x`. -/
def jwPar (n : Cfg N) (x : Fin N) : ZMod 2 := ∑ k ∈ univ.filter (fun k => k < x), n k

/-- Matrix of the annihilation operator `c_x` in the occupation basis: entry `(n', n)` is
`⟨n'| c_x |n⟩`, i.e. `c_x |n⟩ = (-1)^{∑_{k<x} n_k} |n - e_x⟩` if `n_x = 1` and `0` otherwise. -/
def ann (x : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := fun n' n =>
  if n x = 1 ∧ n' = Function.update n x 0 then sgn (jwPar n x) else 0

/-- Creation operator `c_x† = (c_x)ᵀ` (all entries are real). -/
def cre (x : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := (ann x)ᵀ

/-- Second quantization `dΓ(A) = ∑_{x,y} A_{xy} c_x† c_y` of a one-body matrix `A`. -/
def dGamma (A : Matrix (Fin N) (Fin N) ℝ) : Matrix (Cfg N) (Cfg N) ℝ :=
  ∑ x, ∑ y, A x y • (cre x * ann y)

/-- Number of particles of the configuration `n` in the region `A`. -/
def count (A : Finset (Fin N)) (n : Cfg N) : ℝ := ∑ x ∈ A, occ (n x)

/-- Number operator `N_A = ∑_{x ∈ A} c_x† c_x` of a region `A` (diagonal). -/
def numOp (A : Finset (Fin N)) : Matrix (Cfg N) (Cfg N) ℝ := diagonal (count A)

/-- Spectral projection of `N_A` onto the eigenvalue `k`. -/
def numProj (A : Finset (Fin N)) (k : ℕ) : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun n => if count A n = k then 1 else 0)

/-- The Fock configuration whose occupied set is `S`. -/
def fockCfg (S : Finset (Fin N)) : Cfg N := fun x => if x ∈ S then 1 else 0

/-! ### One-body operators -/

/-- Adjacency matrix of the path `0 — 1 — ⋯ — (N-1)` (Dirichlet ends). -/
def pathAdj (N : ℕ) : Matrix (Fin N) (Fin N) ℝ := fun x y =>
  if x.val + 1 = y.val ∨ y.val + 1 = x.val then 1 else 0

/-- The dual one-body Hamiltonian in the charge sector `q`:
`h_q = -J Δ_N + 2h diag(q)` with `q_x = sgn (q x) ∈ {±1}`. -/
def hq (J h : ℝ) (q : Cfg N) : Matrix (Fin N) (Fin N) ℝ :=
  (-J) • pathAdj N + (2 * h) • diagonal (fun x => sgn (q x))

/-- Staggering matrix `Γ = diag((-1)^x)`. -/
def stag (N : ℕ) : Matrix (Fin N) (Fin N) ℝ := diagonal (fun x => (-1 : ℝ) ^ x.val)

/-! ### Links -/

/-- `σˣ` on link `b` in the `σˣ` eigenbasis: diagonal with entries `sgn (x b)`. -/
def sx (b : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := diagonal (fun x => sgn (x b))

/-- `σᶻ` on link `b` in the `σˣ` eigenbasis: flips the label of link `b`. -/
def sz (b : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := fun x' x =>
  if x' = Function.update x b (x b + 1) then 1 else 0

/-- Label of the link `b_{j-1}` to the left of site `j`; the fictitious link `b_{-1}` has label `0`
(i.e. `σˣ_{b_{-1}} := 1`). -/
def prevLink (x : Cfg N) (j : Fin N) : ZMod 2 :=
  if hj : j.val = 0 then 0 else x ⟨j.val - 1, by have := j.isLt; omega⟩

/-- Electric term `σˣ_{b_{j-1}} σˣ_{b_j}` at site `j` (diagonal in the `σˣ` basis). -/
def elec (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun x => sgn (prevLink x j + x j))

/-- Nearest-neighbour hopping `c_j† c_{j+1} + c_{j+1}† c_j` across the bond `(j, j+1)`;
it is `0` for the last site `j = N-1` (the boundary link carries no hopping). -/
def hop (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ :=
  if hj : j.val + 1 < N then
    cre j * ann ⟨j.val + 1, hj⟩ + cre ⟨j.val + 1, hj⟩ * ann j
  else 0

/-! ### The lattice gauge theory and its dual -/

/-- The Hamiltonian of the 1D ℤ₂ lattice gauge theory of Smith–Knolle–Kovrizhin–Moessner,
with boundary link and `σˣ_{b_{-1}} := 1`, on `Cfg N × Cfg N` (occupations × link labels):
`H = -J ∑_j σᶻ_{b_j} (c_j† c_{j+1} + h.c.) - h ∑_j σˣ_{b_{j-1}} σˣ_{b_j}`. -/
def H (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • ∑ j, hop j ⊗ₖ sz j - h • ∑ j, (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elec j

/-- Gauss operator `G_j = σˣ_{b_{j-1}} (-1)^{n_j} σˣ_{b_j}` (diagonal). -/
def gauss (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun p => sgn (p.1 j) * sgn (prevLink p.2 j + p.2 j))

/-- The sector Hamiltonian `dΓ(h_q) - h ∑_j q_j` acting on the matter Fock space. -/
def Hsec (J h : ℝ) (q : Cfg N) : Matrix (Cfg N) (Cfg N) ℝ :=
  dGamma (hq J h q) - (h * ∑ j, sgn (q j)) • (1 : Matrix (Cfg N) (Cfg N) ℝ)

/-- The dual Hamiltonian on `Cfg N × Cfg N` (occupations × charges): block diagonal in the
charges, `⊕_q (dΓ(h_q) - h ∑_j q_j)`. -/
def Hdual (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  blockDiagonal (fun q => Hsec J h q)

/-- Link labels dual to `(n, q)`: `x_{b_j} = ∑_{k ≤ j} (q_k + n_k)` in `ZMod 2`. -/
def linkOf (n q : Cfg N) : Cfg N := fun j => ∑ k ∈ univ.filter (fun k => k ≤ j), (q k + n k)

/-- Gauss charges of `(n, x)`: `q_j = x_{b_{j-1}} + x_{b_j} + n_j` in `ZMod 2`. -/
def chargeOf (n x : Cfg N) : Cfg N := fun j => prevLink x j + x j + n j

end DFL
