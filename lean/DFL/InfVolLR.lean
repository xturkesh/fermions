import DFL.LRFull
import DFL.InfVol

/-!
# The Lieb–Robinson bound in the thermodynamic limit (paper: Theorem `thm:infvol`(c))

Let `X ⊂ ℤ` be finite, `Y ⊂ ℤ` finite, and `A = ∑_S a_S γ_S` an expansion in ordered Majorana monomials of the dressed
fermions with sites in `X` (`|A|₁ ≤ ∑_S |a_S|`); the coefficients do not depend on the volume. On `Λ_M = [-M, M]`
(relabelled as `Fin (2M+1)`, `InfVolAux.toFin`), let `A_M` be the corresponding operator and `B` any observable localized in
`Y` (`LocB`, e.g. `B ∈ 𝔊(Y)`). Then, under (EC), for every `t` and every `M` with `X ∪ Y ⊂ Λ_M`,
`E_q ‖P_q [A_M(t), B]‖ ≤ 8 |A|₁ ‖B‖ ∑_{x ∈ X, y ∈ Y} F(|x - y|)` (`infvol_LR`): the bound of the paper holds uniformly in
`M`, in particular for the `limsup`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator

noncomputable section

namespace DFL

/-- A Majorana index list with sites in `ℤ`, relabelled into `Λ_M`. -/
def majListM (M : ℕ) (l : List (ℤ × Bool)) : List (Fin (2 * M + 1) × Bool) :=
  l.map (fun p => (InfVolAux.toFin M p.1, p.2))

/-- **Infinite-volume Lieb–Robinson bound** (paper: Theorem `thm:infvol`(c)), uniformly in the volume. -/
theorem infvol_LR (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (X Y : Finset ℤ)
    (s : Finset (List (ℤ × Bool))) (hs : ∀ l ∈ s, l.Nodup ∧ ∀ p ∈ l, p.1 ∈ X)
    (a : List (ℤ × Bool) → ℂ) (M : ℕ) (hM : ∀ z ∈ X ∪ Y, |z| ≤ M)
    (B : Matrix (Cfg (2 * M + 1) × Cfg (2 * M + 1)) (Cfg (2 * M + 1) × Cfg (2 * M + 1)) ℂ)
    (hB : LocB (Y.image (InfVolAux.toFin M)) B) (t : ℝ) :
    avgQ (fun q : Cfg (2 * M + 1) => ‖cplx (gaussProj q) *
        (heisT J h t (∑ l ∈ s, a l • majMono (majListM M l)) * B -
          B * heisT J h t (∑ l ∈ s, a l • majMono (majListM M l)))‖) ≤
      8 * (∑ l ∈ s, ‖a l‖) * ‖B‖ * ∑ x ∈ X, ∑ y ∈ Y, F (x - y).natAbs := by
  sorry

end DFL
