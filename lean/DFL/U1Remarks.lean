import DFL.QLM
import DFL.Wilson

/-!
# Remarks on the U(1) models (paper: Remarks `rem:adhh`(a), `rem:wilson`(c), `rem:gausscap`(b))

* `energy_conservation_bound` (`rem:adhh`(a)): if `H = ν 𝒩 + H'` with `ν > 0` and `|⟨H'⟩_v| ≤ C ‖v‖²`, then
  `|⟨𝒩⟩_t - ⟨𝒩⟩_0| ≤ 2C/ν` for every normalized initial state and all `t`, since `ν⟨𝒩⟩_t + ⟨H'⟩_t` is
  conserved.
* `spinPlus`, `spinPlus_opNorm`, `max_spinPlus_weight_even`, `max_spinPlus_weight_odd`, `cS_bounds`
  (`rem:wilson`(c)): for spin `S = n/2` the raising operator `S⁺` has `‖S⁺‖ = max_m ((S-m)(S+m+1))^{1/2}`; the
  maximum is `S(S+1)` for integer `S` and `(S+1/2)²` for half-integer `S`, so
  `1 ≤ c_S = ‖S⁺‖/√(S(S+1)) ≤ 2/√3`.
* `shannon_uniform`, `qlm_numberEntropy_attained` (`rem:gausscap`(b)): the equal-weight superposition of states
  of one sector with `r` distinct values of `N_A` has number entropy `ln r`, so the caps of `prop:gausscap` are
  attained whenever a sector contains such states.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

/-! ### Energy conservation (paper: `rem:adhh`(a)) -/

/-- **An elementary bound from energy conservation** (paper: `rem:adhh`(a)). -/
theorem energy_conservation_bound {ι : Type*} [Fintype ι] [DecidableEq ι]
    (Nop Hp : Matrix ι ι ℝ) (hNop : Nop.IsSymm) (hHp : Hp.IsSymm) (ν C : ℝ) (hν : 0 < ν)
    (hC : ∀ v : ι → ℂ, |(expect (cplx Hp) v).re| ≤ C * ∑ i, Complex.normSq (v i))
    (ψ : ι → ℂ) (hψ : ∑ i, Complex.normSq (ψ i) = 1) (t : ℝ) :
    |(expect (cplx Nop) (evo (ν • Nop + Hp) t *ᵥ ψ)).re - (expect (cplx Nop) ψ).re| ≤ 2 * C / ν := by
  sorry

/-! ### The norm of the spin raising operator (paper: `rem:wilson`(c)) -/

/-- The spin raising operator for spin `S = n/2` on the `n + 1` states `m = -S + k`, `k = 0, …, n`:
`S⁺ |k⟩ = ((n - k)(k + 1))^{1/2} |k + 1⟩`. -/
def spinPlus (n : ℕ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ := fun i j =>
  if i.val = j.val + 1 then Real.sqrt (((n - j.val : ℕ) : ℝ) * ((j.val + 1 : ℕ) : ℝ)) else 0

open scoped Matrix.Norms.L2Operator in
/-- `‖S⁺‖ = max_k ((n - k)(k + 1))^{1/2}` (operator norm). -/
theorem spinPlus_opNorm (n : ℕ) (hn : 0 < n) :
    ‖spinPlus n‖ = Real.sqrt ((Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ)))) := by
  sorry

/-- Integer spin `S = m`: `max_k (2m - k)(k + 1) = m (m + 1) = S(S+1)`. -/
theorem max_spinPlus_weight_even (m : ℕ) (hm : 0 < m) :
    (Finset.range (2 * m)).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((2 * m - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) = (m : ℝ) * (m + 1) := by
  sorry

/-- Half-integer spin `S = m + 1/2`: `max_k (2m + 1 - k)(k + 1) = (m + 1)² = (S + 1/2)²`. -/
theorem max_spinPlus_weight_odd (m : ℕ) :
    (Finset.range (2 * m + 1)).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((2 * m + 1 - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) = ((m : ℝ) + 1) ^ 2 := by
  sorry

/-- `1 ≤ c_S ≤ 2/√3` for every spin `S = n/2 > 0`, where `c_S² = max_k (n - k)(k + 1) / (S(S+1))`. -/
theorem cS_bounds (n : ℕ) (hn : 0 < n) :
    1 ≤ (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) / ((n : ℝ) / 2 * ((n : ℝ) / 2 + 1)) ∧
      (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) / ((n : ℝ) / 2 * ((n : ℝ) / 2 + 1)) ≤
        4 / 3 := by
  sorry

/-! ### The caps of `prop:gausscap` are attained (paper: `rem:gausscap`(b)) -/

/-- The uniform distribution on `r` values has Shannon entropy `ln r`. -/
theorem shannon_uniform (M : ℕ) (K : Finset ℕ) (hK : ∀ k ∈ K, k ≤ M) (hne : K.Nonempty) :
    shannon M (fun k => if k ∈ K then 1 / (K.card : ℝ) else 0) = Real.log K.card := by
  sorry

/-- **The caps are attained** (paper: `rem:gausscap`(b)): if `φ_k`, `k ∈ K`, are normalized states with
`N_{[a,b]} φ_k = k φ_k`, their equal-weight superposition is normalized and has number entropy
`ln |K|`. -/
theorem qlm_numberEntropy_attained {N : ℕ} (S : ℕ) (a b : ℕ) (K : Finset ℕ)
    (hK : ∀ k ∈ K, k ≤ N) (hne : K.Nonempty)
    (φ : ℕ → Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hφn : ∀ k ∈ K, ∑ c, Complex.normSq (φ k c) = 1)
    (hφs : ∀ k ∈ K, ∀ c, φ k c ≠ 0 → countInt a b c.1 = k) :
    (∑ c, Complex.normSq ((Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c) = 1) ∧
      shannon N (qlmNumDist S a b (fun c => (Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c)) =
        Real.log K.card := by
  sorry

end DFL
