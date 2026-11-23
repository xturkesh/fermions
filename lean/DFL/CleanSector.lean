import DFL.Imbalance
import DFL.Normalization
import DFL.DimerTransfer

/-!
# Single sectors (paper: Remarks `rem:gaussproj`(i),(ii), `rem:main-imbalance`(ii), `rem:dimersigns`)

The clean sector is `q ≡ 1` (all Gauss charges `+1`; the label `0 : ZMod 2`), with
`h_{1_N} = -J Δ_N + 2h`.

* `corrQ_clean_ends` (`rem:gaussproj`(i)): no decay of the eigenfunction correlator in the clean sector,
  `Q(1, N) = ∑_α |φ_α(1) φ_α(N)| = 1` for every `N`, the largest possible value. (The eigenvectors of the
  reflection-symmetric `h_{1_N}` are even or odd, since the spectrum is simple.)
* `stagWeight_clean`, `imbInf_clean`, `imbOne_clean_timeAvg` (`rem:gaussproj`(ii)): for even `N`,
  `Γ h + h Γ = 4h Γ` and `2h` is not an eigenvalue, so every eigenvector has `a_α = ⟨φ_α, Γ φ_α⟩ = 0` and the
  infinite-time imbalance of the clean sector vanishes.
* `stagWeight_h_zero` (`rem:main-imbalance`(ii)): for `h = 0` and even `N`, the chiral symmetry
  `Γ h_q Γ = -h_q` forces `a_α = 0` in every sector.
* `hq_flip_all`, `norm_evo_flip_all` (`rem:dimersigns`(i)): the global reversal `h_{-q} = -Γ h_q Γ` is harmless,
  `|e^{-ith_{-q}}_{xy}| = |e^{-ith_q}_{xy}|`; `trace_hq_dimer_example`, `not_similar_dimer_example`: reversing the
  charges of one dimer out of two is not a unitary equivalence (`Tr h_{(1,1,1,1)} = 8h ≠ 0 = Tr h_{(1,1,-1,-1)}`).
* `antidimer_commutator`, `antidimer_comm_ne` (`rem:dimersigns`(ii)): the anti-dimer transfer matrices
  `T_a T_b`, `T_b T_a` (`a = -λ - E`, `b = λ - E`) satisfy
  `[T_a T_b, T_b T_a] = [[b² - a², (a-b)ab], [-(a-b)ab, a² - b²]]` and never commute for `λ ≠ 0`.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### The clean sector (paper: `rem:gaussproj`) -/

/-- **No decay of the eigenfunction correlator in the clean sector** (paper: `rem:gaussproj`(i)):
`Q_{1_N}(1, N) = 1`. -/
theorem corrQ_clean_ends (J h : ℝ) (hJ : J ≠ 0) (hN : 0 < N) :
    corrQ (hq J h (0 : Cfg N)) (hq_isHermitian J h 0) ⟨0, hN⟩ ⟨N - 1, by omega⟩ = 1 := by
  sorry

/-- **No memory of a density wave in the clean sector** (paper: `rem:gaussproj`(ii)): for even `N`, every
eigenvector of `h_{1_N}` has staggered weight `a_α = 0`. -/
theorem stagWeight_clean (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (α : Fin N) :
    stagWeight J h (0 : Cfg N) α = 0 := by
  sorry

/-- The infinite-time imbalance of the clean sector is `N⁻¹ ∑_α a_α² = 0` (paper: `rem:gaussproj`(ii)). -/
theorem imbInf_clean (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) : imbInf J h (0 : Cfg N) = 0 := by
  sorry

/-- The time-averaged imbalance of the clean sector tends to `0` (paper: `rem:gaussproj`(ii)). -/
theorem imbOne_clean_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h (0 : Cfg N) t) atTop (𝓝 0) := by
  sorry

/-- **Chiral symmetry at `h = 0`** (paper: `rem:main-imbalance`(ii)): for even `N`, `a_α = 0` in every
sector. -/
theorem stagWeight_h_zero (J : ℝ) (hJ : J ≠ 0) (hN : Even N) (q : Cfg N) (α : Fin N) :
    stagWeight J 0 q α = 0 := by
  sorry

/-! ### Deterministic signs (paper: `rem:dimersigns`) -/

/-- The global reversal of the charges: `h_{-q} = -Γ h_q Γ` (paper: `rem:dimersigns`(i)). -/
theorem hq_flip_all (J h : ℝ) (q : Cfg N) :
    hq J h (fun x => q x + 1) = -(stag N * hq J h q * stag N) := by
  sorry

/-- `|e^{-ith_{-q}}_{xy}| = |e^{-ith_q}_{xy}|` (paper: `rem:dimersigns`(i)). -/
theorem norm_evo_flip_all (J h : ℝ) (q : Cfg N) (t : ℝ) (x y : Fin N) :
    ‖evo (hq J h (fun x => q x + 1)) t x y‖ = ‖evo (hq J h q) t x y‖ := by
  sorry

/-- `Tr h_{(1,1,1,1)} = 8h` and `Tr h_{(1,1,-1,-1)} = 0` (paper: `rem:dimersigns`(i)). -/
theorem trace_hq_dimer_example (J h : ℝ) :
    trace (hq J h (0 : Cfg 4)) = 8 * h ∧ trace (hq J h (![0, 0, 1, 1] : Cfg 4)) = 0 := by
  sorry

/-- Reversing the charges of one dimer out of two is not a similarity transformation for `h ≠ 0`
(paper: `rem:dimersigns`(i)). -/
theorem not_similar_dimer_example (J h : ℝ) (hh : h ≠ 0) (U : Matrix (Fin 4) (Fin 4) ℝ)
    (hU : IsUnit U) : U * hq J h (0 : Cfg 4) * U⁻¹ ≠ hq J h (![0, 0, 1, 1] : Cfg 4) := by
  sorry

/-- The commutator of the anti-dimer transfer matrices (paper: `rem:dimersigns`(ii)). -/
theorem antidimer_commutator (a b : ℂ) :
    tmC a * tmC b * (tmC b * tmC a) - tmC b * tmC a * (tmC a * tmC b) =
      !![b ^ 2 - a ^ 2, (a - b) * a * b; -((a - b) * a * b), a ^ 2 - b ^ 2] := by
  sorry

/-- **The anti-dimer law has no critical energy** (paper: `rem:dimersigns`(ii)): for `λ ≠ 0` the two
anti-dimer transfer matrices never commute. -/
theorem antidimer_comm_ne (lam : ℝ) (hl : lam ≠ 0) (E : ℂ) :
    tmC (-lam - E) * tmC (lam - E) * (tmC (lam - E) * tmC (-lam - E)) ≠
      tmC (lam - E) * tmC (-lam - E) * (tmC (-lam - E) * tmC (lam - E)) := by
  sorry

end DFL
