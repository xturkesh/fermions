import DFL.Cluster
import DFL.Graph

/-!
# Strong-coupling confinement in the gauge theory (paper: Remark `rem:clusterLGT`)

For `σᶻ`-polarized links and a Fock state `|n⟩` with occupied set `S_n`, the density of the evolved LGT state is
`⟨Ψ(t), n_y Ψ(t)⟩ = |sectors|⁻¹ ∑_{q ∈ n + im ∂} ∑_{a ∈ S_n} |(e^{-ith_q})_{ya}|²` (`sector_average_matter`, `density`),
and `h_q = h_{q,+}` is the zero-flux operator of `Cluster.lean`. Hence, by `cluster_avg`, the expected number of
particles outside given regions `B a ∋ a` is at most `|S_n| [(8J/g + 16|t|J²/g)² + K P]` (`cluster_LGT`), where
`K = 2^N / |sectors|` bounds the density of the sector law with respect to the uniform law and `P` bounds the
probability that the equal-charge cluster of an occupied site leaves its region (the percolation input of
`prop:cluster`(c)).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

open Classical in
/-- **Meaning for the gauge theory** (paper: Remark `rem:clusterLGT`). -/
theorem cluster_LGT (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (n : Cfg N) (t : ℝ)
    (B : Fin N → Finset (Fin N)) (P : ℝ)
    (hP : ∀ a, n a = 1 → ((univ.filter (fun q : Cfg N =>
      ∃ y, y ∉ B a ∧ G.InCluster (fun i => sgn (q i)) a y)).card : ℝ) / 2 ^ N ≤ P) :
    ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
        (expect (cplx (diagonal (fun m : Cfg N => occ (m y)) ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
          (evo (G.ham J h) t *ᵥ G.initZ (Pi.single n 1))).re ≤
      ((univ.filter (fun a => n a = 1)).card : ℝ) *
        ((8 * J / (4 * |h| - 8 * J) + 16 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 +
          (2 ^ N / ((G.sectors n).card : ℝ)) * P) := by
  sorry

end FinGraph

end DFL
