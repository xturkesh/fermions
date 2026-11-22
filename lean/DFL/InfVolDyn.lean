import DFL.QuasiFree
import DFL.InfVolArea

/-!
# Limit dynamics of local gauge-invariant observables (paper: Theorem `thm:infvol`(a))

Let `X = [a, b] ⊂ ℤ` with `L = b - a + 1` sites, labelled `i ↦ a + i` by `Fin L`. A local gauge-invariant observable
`O ∈ 𝔊(X)` acts in the dual frame as `U†OU = ∑_q O_{q|_X} ⊗ |q⟩⟨q|`, where `O_{q|_X}` is an operator on the Fock space
`𝓕_X` of the modes of `X` depending only on the charges in `X`. On `Λ_M = [-M, M]` we build it from a family
`O : Cfg L → Matrix (Cfg L) (Cfg L) ℂ` (`locObsM`): the operator `O_{q|_X}` on `𝓕_X` is embedded into the Fock space of
`Λ_M` through the identification `𝓕 ≅ 𝓕_X ⊗ 𝓕_B` of Lemma `lem:slater` (`embF`, with the reordering sign; then
`⟨ψ, embF(O) ψ⟩ = Tr(ρ_X(ψ) O)`, `expect_embF`).

With `K^ℤ_q(t)|_X` the restricted correlation matrix (`KZL`) and `ρ[·]` the quasi-free state (`qfState`):

* `tendsto_omega`: for every `t`, `⟨Ψ^{(M)}(t), O Ψ^{(M)}(t)⟩ → ω_t(O) := E Tr(ρ[K^ℤ_q(t)|_X] O_{q|_X})`;
* `tendstoUniformlyOn_sector`: sector by sector, `Tr(ρ_X(ψ^{(M)}_q(t)) O_{q|_X}) → Tr(ρ[K^ℤ_q(t)|_X] O_{q|_X})` uniformly
  in `q` and in `t ∈ [-T, T]`.
-/

open Matrix Finset MeasureTheory Filter Topology

noncomputable section

namespace DFL

/-- The embedding of an operator on the modes of `R` into the Fock space, through `𝓕 ≅ 𝓕_R ⊗ 𝓕_B` (`O ⊗ 1`). -/
def embF {N : ℕ} (R : Finset (Fin N)) (O : Matrix (R → ZMod 2) (R → ZMod 2) ℂ) : Matrix (Cfg N) (Cfg N) ℂ :=
  fun n' n => if (splitR R n').2 = (splitR R n).2 then
    ((reorderSign R (splitR R n').1 (splitR R n).2 * reorderSign R (splitR R n).1 (splitR R n).2 : ℝ) : ℂ) *
      O (splitR R n').1 (splitR R n).1 else 0

/-- **Expectations of embedded operators are traces against the reduced state**: `⟨ψ, (O ⊗ 1) ψ⟩ = Tr(ρ_R(ψ) O)`. -/
theorem expect_embF {N : ℕ} (R : Finset (Fin N)) (O : Matrix (R → ZMod 2) (R → ZMod 2) ℂ) (v : Cfg N → ℂ) :
    expect (embF R O) v = (rdmF R v * O).trace := by
  sorry

/-- The number of sites of `[a, b]`. -/
def ivlLen (a b : ℤ) : ℕ := (b - a + 1).toNat

/-- The site `a + i` of `X = [a, b]` as a site of `Λ_M` (relabelled in `Fin (2M+1)`). -/
def siteM (M : ℕ) (a b : ℤ) (i : Fin (ivlLen a b)) : Fin (2 * M + 1) := InfVolAux.toFin M (a + i.val)

/-- The charges of `q ∈ Cfg (2M+1)` on `X`. -/
def chargesX (M : ℕ) (a b : ℤ) (q : Cfg (2 * M + 1)) : Cfg (ivlLen a b) := fun i => q (siteM M a b i)

/-- An operator on `𝓕_X` (modes labelled by `Fin L`), as an operator on the modes of `bulkM M a b`. -/
def relabelX (M : ℕ) (a b : ℤ) (O : Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    Matrix (bulkM M a b → ZMod 2) (bulkM M a b → ZMod 2) ℂ :=
  O.submatrix (fun α i => if h : siteM M a b i ∈ bulkM M a b then α ⟨siteM M a b i, h⟩ else 0)
    (fun α i => if h : siteM M a b i ∈ bulkM M a b then α ⟨siteM M a b i, h⟩ else 0)

/-- **The local gauge-invariant observable on `Λ_M`**: `U (∑_q (O_{q|_X} ⊗ 1) ⊗ |q⟩⟨q|) U†`. -/
def locObsM (M : ℕ) (a b : ℤ) (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    Matrix (Cfg (2 * M + 1) × Cfg (2 * M + 1)) (Cfg (2 * M + 1) × Cfg (2 * M + 1)) ℂ :=
  (blockDiagonal (fun q : Cfg (2 * M + 1) => embF (bulkM M a b) (relabelX M a b (O (chargesX M a b q))))).submatrix
    (dualEquiv (N := 2 * M + 1)).symm (dualEquiv (N := 2 * M + 1)).symm

/-- `K^ℤ_q(t)|_X` with the modes of `X` labelled by `Fin L`. -/
def KZL (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) : Matrix (Fin (ivlLen a b)) (Fin (ivlLen a b)) ℂ :=
  fun i j => KZ J h q n t (a + i.val) (a + j.val)

/-- **Limit dynamics** (paper: Theorem `thm:infvol`(a)): `⟨Ψ^{(M)}(t), O Ψ^{(M)}(t)⟩ → E Tr(ρ[K^ℤ_q(t)|_X] O_{q|_X})`. -/
theorem tendsto_omega (J h : ℝ) (n : CfgZ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) (t : ℝ) :
    Tendsto (fun M : ℕ => expect (locObsM M a b O) (lgtState J h (fockVec (occM M n)) t)) atTop
      (𝓝 (∫ q, (qfState (KZL J h q n t a b) * O (fun i => q (a + i.val))).trace ∂muZ)) := by
  sorry

/-- **Sector by sector, uniformly in `q` and in `t ∈ [-T, T]`** (paper: Theorem `thm:infvol`(a)). -/
theorem tendstoUniformlyOn_sector (J h : ℝ) (n : CfgZ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) (T : ℝ) :
    TendstoUniformlyOn
      (fun (M : ℕ) (p : CfgZ × ℝ) =>
        (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M p.1))) p.2 *ᵥ fockVec (occM M n)) *
          relabelX M a b (O (fun i => p.1 (a + i.val)))).trace)
      (fun p : CfgZ × ℝ => (qfState (KZL J h p.1 n p.2 a b) * O (fun i => p.1 (a + i.val))).trace)
      atTop {p | |p.2| ≤ T} := by
  sorry

end DFL
