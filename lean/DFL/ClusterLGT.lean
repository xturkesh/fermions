import DFL.Cluster
import DFL.Graph
import DFL.FreeDyn

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

namespace ClusterLGTAux

/-- The adjacency matrix is symmetric (entrywise). -/
lemma adj_symm (i j : Fin N) : G.adj i j = G.adj j i := by
  simp only [adj]
  congr 2
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  tauto

/-- The one-body Hamiltonian `h_q` of a charge sector is symmetric. -/
lemma hq_isSymm (J h : ℝ) (q : Cfg N) : (G.hq J h q).IsSymm := by
  refine Matrix.IsSymm.ext fun i j => ?_
  simp only [hq, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply,
    adj_symm G j i]
  by_cases hij : i = j
  · subst hij; rfl
  · rw [ite_eq_right (Ne.symm hij), ite_eq_right hij]

/-- The Fock configuration of the occupied set of `n` is `n`. -/
lemma fockCfg_filter_eq (n : Cfg N) [DecidablePred (fun a => n a = 1)] :
    fockCfg (univ.filter (fun a => n a = 1)) = n := by
  funext x
  simp only [fockCfg, Finset.mem_filter, Finset.mem_univ, true_and]
  rcases FreeDynAux.zmod2_cases (n x) with hx | hx <;> simp [hx]

end ClusterLGTAux

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
  set S : Finset (Fin N) := univ.filter (fun a => n a = 1) with hS
  set E : ℝ := (8 * J / (4 * |h| - 8 * J) + 16 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 with hE
  set c : ℝ := ((G.sectors n).card : ℝ) with hc
  have hne : (G.sectors n).Nonempty :=
    ⟨_, Finset.mem_image_of_mem (fun x => n + G.bdry x) (Finset.mem_univ (0 : Cfg M))⟩
  have hc0 : 0 < c := by
    rw [hc]
    exact_mod_cast Finset.card_pos.mpr hne
  -- the density at each site: an average of free-fermion densities
  have hdens : ∀ y, (expect (cplx (diagonal (fun m : Cfg N => occ (m y)) ⊗ₖ
      (1 : Matrix (Cfg M) (Cfg M) ℝ))) (evo (G.ham J h) t *ᵥ G.initZ (Pi.single n 1))).re =
      (1 / c) * ∑ q ∈ G.sectors n, ∑ a ∈ S, ‖evo (G.hq J h q) t y a‖ ^ 2 := by
    intro y
    have hfock : (Pi.single n 1 : Cfg N → ℂ) = fockVec S := by
      rw [fockVec, hS, ClusterLGTAux.fockCfg_filter_eq]
    rw [G.sector_average_matter J h n (fun m => occ (m y)) t, hfock]
    simp only [density _ (ClusterLGTAux.hq_isSymm G J h _) S y t, Complex.normSq_eq_norm_sq]
    have hre : (1 / ((G.sectors n).card : ℂ)) * ∑ q ∈ G.sectors n,
        (((∑ a ∈ S, ‖evo (G.hq J h q) t y a‖ ^ 2 : ℝ)) : ℂ) =
        (((1 / c) * ∑ q ∈ G.sectors n, ∑ a ∈ S, ‖evo (G.hq J h q) t y a‖ ^ 2 : ℝ) : ℂ) := by
      rw [hc]
      push_cast
      ring
    rw [hre, Complex.ofReal_re]
  -- the uniform law on the sectors, as a weight on all charge configurations
  set w : Cfg N → ℝ := fun q => if q ∈ G.sectors n then 1 / c else 0 with hw
  have hc1 : 0 ≤ 1 / c := le_of_lt (one_div_pos.mpr hc0)
  have hw0 : ∀ q, 0 ≤ w q := by
    intro q
    simp only [hw]
    split_ifs
    · exact hc1
    · exact le_rfl
  have hw1 : ∑ q, w q = 1 := by
    simp only [hw]
    rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.sum_const, nsmul_eq_mul]
    exact mul_one_div_cancel hc0.ne'
  have hK : (2 ^ N / c) / 2 ^ N = 1 / c := by
    field_simp
  have hwK : ∀ q, w q ≤ (2 ^ N / c) / 2 ^ N := by
    intro q
    rw [hK]
    simp only [hw]
    split_ifs
    · exact le_rfl
    · exact hc1
  have hK0 : 0 ≤ 2 ^ N / c := div_nonneg (by positivity) hc0.le
  -- the bound for each occupied site
  have hbound : ∀ a ∈ S, (1 / c) * ∑ q ∈ G.sectors n,
      ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a), ‖evo (G.hq J h q) t y a‖ ^ 2 ≤
      E + (2 ^ N / c) * P := by
    intro a ha
    have ha1 : n a = 1 := (Finset.mem_filter.mp ha).2
    have hsub : univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a) ⊆
        univ.filter (fun y => y ∉ B a) := by
      intro y hy
      rw [Finset.mem_filter] at hy ⊢
      exact ⟨Finset.mem_univ _, hy.2 a ha1⟩
    have havg := G.cluster_avg J h hJ hh (fun _ => 1) (fun _ => Or.inl rfl) hdeg t a (B a) w
      (2 ^ N / c) hw0 hw1 hwK
    simp only [hqs_one] at havg
    have hleft : ∑ q, w q * ∑ y ∈ univ.filter (fun y => y ∉ B a),
        ‖evo (G.hq J h q) t y a‖ ^ 2 =
        (1 / c) * ∑ q ∈ G.sectors n, ∑ y ∈ univ.filter (fun y => y ∉ B a),
          ‖evo (G.hq J h q) t y a‖ ^ 2 := by
      simp only [hw, ite_mul, zero_mul]
      rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.mul_sum]
    calc (1 / c) * ∑ q ∈ G.sectors n,
          ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a), ‖evo (G.hq J h q) t y a‖ ^ 2
        ≤ (1 / c) * ∑ q ∈ G.sectors n, ∑ y ∈ univ.filter (fun y => y ∉ B a),
          ‖evo (G.hq J h q) t y a‖ ^ 2 := by
          refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => ?_) hc1
          exact Finset.sum_le_sum_of_subset_of_nonneg hsub (fun _ _ _ => by positivity)
      _ = ∑ q, w q * ∑ y ∈ univ.filter (fun y => y ∉ B a), ‖evo (G.hq J h q) t y a‖ ^ 2 :=
          hleft.symm
      _ ≤ E + (2 ^ N / c) * ((univ.filter (fun q : Cfg N =>
          ∃ y, y ∉ B a ∧ G.InCluster (fun i => sgn (q i)) a y)).card : ℝ) / 2 ^ N := havg
      _ ≤ E + (2 ^ N / c) * P := by
          rw [mul_div_assoc]
          exact add_le_add le_rfl (mul_le_mul_of_nonneg_left (hP a ha1) hK0)
  -- swap the sums and add up the bounds
  calc ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
        (expect (cplx (diagonal (fun m : Cfg N => occ (m y)) ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
          (evo (G.ham J h) t *ᵥ G.initZ (Pi.single n 1))).re
      = ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
          (1 / c) * ∑ q ∈ G.sectors n, ∑ a ∈ S, ‖evo (G.hq J h q) t y a‖ ^ 2 :=
        Finset.sum_congr rfl fun y _ => hdens y
    _ = ∑ a ∈ S, (1 / c) * ∑ q ∈ G.sectors n,
          ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a), ‖evo (G.hq J h q) t y a‖ ^ 2 := by
        rw [← Finset.mul_sum, ← Finset.mul_sum]
        congr 1
        calc ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
              ∑ q ∈ G.sectors n, ∑ a ∈ S, ‖evo (G.hq J h q) t y a‖ ^ 2
            = ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
                ∑ a ∈ S, ∑ q ∈ G.sectors n, ‖evo (G.hq J h q) t y a‖ ^ 2 :=
              Finset.sum_congr rfl fun y _ => Finset.sum_comm
          _ = ∑ a ∈ S, ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
                ∑ q ∈ G.sectors n, ‖evo (G.hq J h q) t y a‖ ^ 2 := Finset.sum_comm
          _ = ∑ a ∈ S, ∑ q ∈ G.sectors n,
                ∑ y ∈ univ.filter (fun y => ∀ a, n a = 1 → y ∉ B a),
                  ‖evo (G.hq J h q) t y a‖ ^ 2 :=
              Finset.sum_congr rfl fun a _ => Finset.sum_comm
    _ ≤ ∑ a ∈ S, (E + (2 ^ N / c) * P) := Finset.sum_le_sum hbound
    _ = (S.card : ℝ) * (E + (2 ^ N / c) * P) := by rw [Finset.sum_const, nsmul_eq_mul]

end FinGraph

end DFL
