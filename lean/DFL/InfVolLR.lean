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

namespace InfVolLRAux

/-- `toFin` is injective on the box `[-M, M]`. -/
lemma toFin_inj (M : ℕ) {z w : ℤ} (hz : |z| ≤ M) (hw : |w| ≤ M)
    (h : InfVolAux.toFin M z = InfVolAux.toFin M w) : z = w := by
  have h1 := InfVolAux.toFin_val M z hz
  have h2 := InfVolAux.toFin_val M w hw
  rw [h] at h1
  omega

/-- `toFin` preserves distances in the box `[-M, M]`. -/
lemma dist_toFin (M : ℕ) {x y : ℤ} (hx : |x| ≤ M) (hy : |y| ≤ M) :
    Nat.dist (InfVolAux.toFin M x).val (InfVolAux.toFin M y).val = (x - y).natAbs := by
  rw [← InfVolAux.idxZ_eq_toFin M x hx, ← InfVolAux.idxZ_eq_toFin M y hy]
  exact InfVolAux.dist_idxZ M x y hx hy

/-- A relabelled list without repetitions and with sites in `X ⊂ [-M, M]` has no repetitions and
sites in the relabelled region. -/
lemma majListM_nodup (M : ℕ) (X : Finset ℤ) (hX : ∀ z ∈ X, |z| ≤ M) (l : List (ℤ × Bool))
    (hl : l.Nodup ∧ ∀ p ∈ l, p.1 ∈ X) :
    (majListM M l).Nodup ∧ ∀ p ∈ majListM M l, p.1 ∈ X.image (InfVolAux.toFin M) := by
  unfold majListM
  refine ⟨hl.1.map_on ?_, fun p hp => ?_⟩
  · intro p hp p' hp' hpp
    simp only [Prod.mk.injEq] at hpp
    exact Prod.ext (toFin_inj M (hX _ (hl.2 p hp)) (hX _ (hl.2 p' hp')) hpp.1) hpp.2
  · obtain ⟨p', hp', rfl⟩ := List.mem_map.mp hp
    exact Finset.mem_image_of_mem _ (hl.2 p' hp')

/-- Regrouping a Majorana expansion along a relabelling `g` of the index lists:
`∑_{l ∈ s} a_l γ_{g(l)} = ∑_{l' ∈ g(s)} (∑_{g(l) = l'} a_l) γ_{l'}`. -/
lemma sum_smul_majMono_comp {N : ℕ} {ι : Type*} (s : Finset ι) (g : ι → List (Fin N × Bool))
    (a : ι → ℂ) :
    ∑ l ∈ s, a l • majMono (g l) =
      ∑ l' ∈ s.image g, (∑ l ∈ s with g l = l', a l) • majMono l' := by
  rw [← Finset.sum_fiberwise_of_maps_to (fun l hl => Finset.mem_image_of_mem g hl)
    (fun l => a l • majMono (g l))]
  refine Finset.sum_congr rfl fun l' _ => ?_
  rw [Finset.sum_smul]
  refine Finset.sum_congr rfl fun l hl => ?_
  rw [(Finset.mem_filter.mp hl).2]

/-- The regrouped coefficients have total weight at most `∑_l |a_l|`. -/
lemma sum_norm_fiber_le {ι κ : Type*} [DecidableEq κ] (s : Finset ι) (g : ι → κ) (a : ι → ℂ) :
    ∑ l' ∈ s.image g, ‖∑ l ∈ s with g l = l', a l‖ ≤ ∑ l ∈ s, ‖a l‖ := by
  rw [← Finset.sum_fiberwise_of_maps_to (fun l hl => Finset.mem_image_of_mem g hl)
    (fun l => ‖a l‖)]
  exact Finset.sum_le_sum fun l' _ => norm_sum_le _ _

end InfVolLRAux

open InfVolLRAux in
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
  have hX : ∀ z ∈ X, |z| ≤ M := fun z hz => hM z (Finset.mem_union_left _ hz)
  have hY : ∀ z ∈ Y, |z| ≤ M := fun z hz => hM z (Finset.mem_union_right _ hz)
  have hA := sum_smul_majMono_comp s (majListM M) a
  have hs' : ∀ l' ∈ s.image (majListM M),
      l'.Nodup ∧ ∀ p ∈ l', p.1 ∈ X.image (InfVolAux.toFin M) := by
    intro l' hl'
    obtain ⟨l, hl, rfl⟩ := Finset.mem_image.mp hl'
    exact majListM_nodup M X hX l (hs l hl)
  -- the sector-wise bound at time `t`
  have hpt : ∀ q : Cfg (2 * M + 1), ‖cplx (gaussProj q) *
      (heisT J h t (∑ l ∈ s, a l • majMono (majListM M l)) * B -
        B * heisT J h t (∑ l ∈ s, a l • majMono (majListM M l)))‖ ≤
      8 * (∑ l ∈ s, ‖a l‖) * ‖B‖ * ∑ x ∈ X.image (InfVolAux.toFin M),
        ∑ y ∈ Y.image (InfVolAux.toFin M), ‖evo (hq J h q) t x y‖ := by
    intro q
    rw [hA]
    refine (LR_sector_maj J h (X.image (InfVolAux.toFin M)) (Y.image (InfVolAux.toFin M))
      (s.image (majListM M)) hs' (fun l' => ∑ l ∈ s with majListM M l = l', a l) B hB q t).trans ?_
    refine mul_le_mul_of_nonneg_right ?_
      (Finset.sum_nonneg fun x _ => Finset.sum_nonneg fun y _ => norm_nonneg _)
    refine mul_le_mul_of_nonneg_right ?_ (norm_nonneg B)
    exact mul_le_mul_of_nonneg_left (sum_norm_fiber_le s (majListM M) a) (by norm_num)
  -- the sector average under (EC)
  have hinjX : Set.InjOn (InfVolAux.toFin M) (X : Set ℤ) := fun x hx x' hx' hxx =>
    toFin_inj M (hX x hx) (hX x' hx') hxx
  have hinjY : Set.InjOn (InfVolAux.toFin M) (Y : Set ℤ) := fun y hy y' hy' hyy =>
    toFin_inj M (hY y hy) (hY y' hy') hyy
  refine (TransportSharpAux.avgQ_mono hpt).trans ?_
  rw [TransportSharpAux.avgQ_const_mul]
  refine mul_le_mul_of_nonneg_left ?_ (by positivity)
  calc avgQ (fun q : Cfg (2 * M + 1) => ∑ x ∈ X.image (InfVolAux.toFin M),
        ∑ y ∈ Y.image (InfVolAux.toFin M), ‖evo (hq J h q) t x y‖)
      = ∑ x ∈ X.image (InfVolAux.toFin M), ∑ y ∈ Y.image (InfVolAux.toFin M),
          avgQ (fun q : Cfg (2 * M + 1) => ‖evo (hq J h q) t x y‖) := by
        rw [TransportSharpAux.avgQ_sum]
        refine Finset.sum_congr rfl fun x _ => ?_
        rw [TransportSharpAux.avgQ_sum]
    _ ≤ ∑ x ∈ X.image (InfVolAux.toFin M), ∑ y ∈ Y.image (InfVolAux.toFin M),
          F (Nat.dist x.val y.val) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
          TransportSharpAux.avgQ_norm_evo_le J h hEC x y t
    _ = ∑ x ∈ X, ∑ y ∈ Y, F (x - y).natAbs := by
        rw [Finset.sum_image hinjX]
        refine Finset.sum_congr rfl fun x hx => ?_
        rw [Finset.sum_image hinjY]
        refine Finset.sum_congr rfl fun y hy => ?_
        rw [dist_toFin M (hX x hx) (hY y hy)]

end DFL
