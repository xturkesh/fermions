import DFL.LRNorm

/-!
# The constants of (EC) degenerate at weak disorder (paper: Remark `rem:rate`(iii))

* `corrQ_free_mid`: for the free chain on `N = 2m + 1` sites (0-indexed), the eigenfunction correlator between
  the first site and the middle site is `Q_{Δ_N}(1, (N+1)/2) = 2 / ((N+1) sin(π/(N+1)))`: with the sine
  eigenvectors `φ_k(x) = (2/(N+1))^{1/2} sin(πkx/(N+1))`, only odd `k` contribute, and
  `∑_{k odd} sin(πk/(N+1)) = 1/sin(π/(N+1))`.
* `tendsto_corrQ_free_mid`: hence `Q_{Δ_N}(1, (N+1)/2) → 2/π`.
* `no_uniform_EC_weak`: consequently no constants `C, η > 0` give `E_q Q^{[1,N]}_q(x, y) ≤ C e^{-η|x-y|}` for all
  `0 < λ ≤ λ₁`, all `N` and all `x, y`. For fixed `N` and `q`, `Q_q` is lower semicontinuous at `λ = 0`
  (`|p(h_q)_{xy}| ≤ sup|p| · Q_q(x, y)` for polynomials `p`, and `p(h_q)_{xy}` is a polynomial in `λ`), so the
  bound would pass to the free chain, contradicting the previous item.
-/

open Filter Topology Real

noncomputable section

namespace DFL

/-- **The free-chain correlator between the end and the middle** (paper: `rem:rate`(iii)): for `N = 2m+1`,
`Q_{Δ_N}(1, (N+1)/2) = 2 / ((N+1) sin(π/(N+1)))`. -/
theorem corrQ_free_mid (J : ℝ) (hJ : J ≠ 0) (m : ℕ) :
    corrQ (hq J 0 (0 : Cfg (2 * m + 1))) (hq_isHermitian J 0 0) ⟨0, by omega⟩ ⟨m, by omega⟩ =
      2 / ((2 * (m : ℝ) + 2) * Real.sin (π / (2 * (m : ℝ) + 2))) := by
  sorry

/-- `Q_{Δ_N}(1, (N+1)/2) → 2/π` along odd `N`. -/
theorem tendsto_corrQ_free_mid (J : ℝ) (hJ : J ≠ 0) :
    Tendsto (fun m : ℕ => corrQ (hq J 0 (0 : Cfg (2 * m + 1))) (hq_isHermitian J 0 0) ⟨0, by omega⟩
      ⟨m, by omega⟩) atTop (𝓝 (2 / π)) := by
  sorry

/-- **The constants of (EC) must degenerate at weak disorder** (paper: `rem:rate`(iii)): for `J ≠ 0` there are no
`C, η > 0` with `E_q Q^{[1,N]}_q(x, y) ≤ C e^{-η|x-y|}` for all `0 < λ = 2h/J ≤ λ₁`, all `N` and all `x, y`. -/
theorem no_uniform_EC_weak (J : ℝ) (hJ : J ≠ 0) (lam1 : ℝ) (hlam1 : 0 < lam1) (C η : ℝ)
    (hη : 0 < η) :
    ¬ ∀ h : ℝ, 0 < 2 * h / J → 2 * h / J ≤ lam1 → ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C * Real.exp (-η * (Nat.dist x.val y.val : ℝ)) := by
  sorry

end DFL
