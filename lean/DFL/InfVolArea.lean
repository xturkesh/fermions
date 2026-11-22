import DFL.InfVolNumber
import DFL.SlaterBulk

/-!
# Area law for bulk intervals in the thermodynamic limit (paper: Theorem `thm:infvol`(d))

For a finite interval `A = [a, b] ⊂ ℤ`, a Fock configuration `n` on `ℤ` and the chain `Λ_M = [-M, M]`, let
`S^{(M)}_q(t)` be the entropy of the reduced state of the modes of `A` (`rdmF`, the fermionic reduced state of
`DFL.SlaterBulk`) in the sector state `ψ^{(M)}_q(t) = Γ(e^{-ith^{(M)}_q}) |n|_{Λ_M}⟩` (`entM`), and
`S^ℤ_q(t) = ∑_i h(ν_i(K^ℤ_q(t)|_A))` (`entZ`). Then

* `S^{(M)}_q(t) → S^ℤ_q(t)` uniformly in `q` (`tendstoUniformly_entM`);
* the bulk conditional entanglement converges: `S_G^{(M)}(t) = E S^{(M)}_q(t) → E S^ℤ_q(t)` (`tendsto_SGM`);
* `E sup_t S^ℤ_q(t) ≤ 4 ln 2 Θ(F)` under (EC) (`lintegral_iSup_entZ_le`);
* for the half-infinite regions `A_M = [-M, m]`, `E S(ρ^q_{A_M}(t)) ≤ 2 ln 2 Θ(F)` for every `M` and `t`
  (`areaLaw_halfInfinite`), in particular for the `limsup`.
-/

open Matrix Finset MeasureTheory Filter Topology

noncomputable section

namespace DFL

/-- `S^{(M)}_q(t)`: entanglement entropy of the modes of `A = [a, b]` in the sector state on `Λ_M`. -/
def entM (J h : ℝ) (n q : CfgZ) (a b : ℤ) (t : ℝ) (M : ℕ) : ℝ :=
  vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M q))) t *ᵥ fockVec (occM M n)))
    (rdmF_isHermitian _ _)

/-- `S^ℤ_q(t) = ∑_i h(ν_i(K^ℤ_q(t)|_A))`. -/
def entZ (J h : ℝ) (n q : CfgZ) (a b : ℤ) (t : ℝ) : ℝ :=
  ∑ i, Real.binEntropy ((KZA_isHermitian J h q n t a b).eigenvalues i)

namespace InfVolAreaAux

/-- The entanglement entropy of a left segment is bounded in time (by a time-independent
majorant). -/
lemma bddAbove_vonNeumann_rdm {N : ℕ} (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (ℓ : ℕ)
    (hℓ : ℓ ≤ N) :
    BddAbove (Set.range fun t : ℝ =>
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) :=
  ⟨2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z, by
    rintro _ ⟨t, rfl⟩
    exact (vonNeumann_slater_sharp (hq J h q) (hq_isSymm J h q) S t ℓ hℓ).trans
      (mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
        AreaLawSharpAux.norm_corrMat_le_corrQ J h q S t x y) AreaLawSharpAux.two_log_two_nonneg)⟩

end InfVolAreaAux

/-- **`S^{(M)}_q(t) → S^ℤ_q(t)` uniformly in `q`** (paper: Theorem `thm:infvol`(d)). -/
theorem tendstoUniformly_entM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) :
    TendstoUniformly (fun (M : ℕ) (q : CfgZ) => entM J h n q a b t M)
      (fun q : CfgZ => entZ J h n q a b t) atTop := by
  sorry

/-- **Convergence of the bulk conditional entanglement** (paper: Theorem `thm:infvol`(d)):
`S_G^{(M)}(t) = E_q S^{(M)}_q(t) → E S^ℤ_q(t)`. -/
theorem tendsto_SGM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) :
    Tendsto (fun M : ℕ => avgQ (fun q' : Cfg (2 * M + 1) =>
        vonNeumann (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
          (rdmF_isHermitian _ _)))
      atTop (𝓝 (∫ q, entZ J h n q a b t ∂muZ)) := by
  sorry

/-- **Area law for bulk intervals on `ℤ`** (paper: Theorem `thm:infvol`(d)): `E sup_t S^ℤ_q(t) ≤ 4 ln 2 Θ(F)`. -/
theorem lintegral_iSup_entZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (a b : ℤ) :
    ∫⁻ q, ⨆ t : ℝ, ENNReal.ofReal (entZ J h n q a b t) ∂muZ ≤
      ENNReal.ofReal (4 * Real.log 2 * thetaF F) := by
  sorry

/-- **Half-infinite regions** (paper: Theorem `thm:infvol`(d)): for `A_M = [-M, m]`,
`E_q S(ρ^q_{A_M}(t)) ≤ 2 ln 2 Θ(F)` for every `M ≥ |m|` and every `t`. -/
theorem areaLaw_halfInfinite (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ)
    (M : ℕ) (hm : |m| ≤ M) :
    avgQ (fun q' : Cfg (2 * M + 1) =>
      vonNeumann (rdm (m + M + 1).toNat (by have h1 := (abs_le.mp hm).1; have h2 := (abs_le.mp hm).2; omega)
        (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)))
        (rdm_isHermitian _ _ _)) ≤ 2 * Real.log 2 * thetaF F := by
  have hℓ : (m + M + 1).toNat ≤ 2 * M + 1 := by
    have h1 := (abs_le.mp hm).1
    have h2 := (abs_le.mp hm).2
    omega
  refine le_trans ?_ (lgt_areaLaw_sharp J h hEC (m + M + 1).toNat hℓ (occM M n))
  exact AreaLawAux.avgQ_mono' fun q' =>
    le_ciSup (InfVolAreaAux.bddAbove_vonNeumann_rdm J h q' (occM M n) _ hℓ) t

end DFL
