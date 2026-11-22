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
  have key : ∀ n' : Cfg N, (embF R O *ᵥ v) n' = ∑ α, ((reorderSign R (splitR R n').1 (splitR R n').2 *
      reorderSign R α (splitR R n').2 : ℝ) : ℂ) * O (splitR R n').1 α *
        v ((splitR R).symm (α, (splitR R n').2)) := by
    intro n'
    simp only [mulVec, dotProduct, embF]
    rw [← (splitR R).symm.sum_comp]
    simp only [Equiv.apply_symm_apply, Fintype.sum_prod_type]
    rw [Finset.sum_comm, Finset.sum_eq_single (splitR R n').2]
    · simp
    · intro β _ hβ
      simp [Ne.symm hβ]
    · simp
  unfold expect
  simp only [dotProduct, key, Pi.star_apply, Matrix.trace, Matrix.diag, Matrix.mul_apply]
  rw [← (splitR R).symm.sum_comp]
  simp only [Equiv.apply_symm_apply, Fintype.sum_prod_type]
  unfold rdmF
  simp only [Finset.mul_sum, Finset.sum_mul]
  conv_rhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun α' _ => ?_
  conv_rhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun β _ => ?_
  refine Finset.sum_congr rfl fun α _ => ?_
  push_cast
  ring

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

namespace InfVolDynAux

open InfVolNumberAux InfVolAreaAux

/-! ### The increasing labelling of the bulk interval -/

/-- The labelling `i ↦ a + i` of `[a, b]` by `Fin (ivlLen a b)`. -/
def finIcc (a b : ℤ) : Fin (ivlLen a b) ≃ Set.Icc a b where
  toFun i := ⟨a + i.val, by
    have := i.isLt
    unfold ivlLen at this
    constructor <;> omega⟩
  invFun x := ⟨(x.1 - a).toNat, by
    obtain ⟨h1, h2⟩ := x.2
    unfold ivlLen
    omega⟩
  left_inv i := by
    apply Fin.ext
    simp only
    omega
  right_inv x := by
    obtain ⟨h1, h2⟩ := x.2
    apply Subtype.ext
    simp only
    omega

/-- For `[a, b] ⊆ Λ_M`, the labelling of the bulk interval `bulkM M a b` by `Fin (ivlLen a b)`. -/
def gM (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) : Fin (ivlLen a b) ≃ bulkM M a b :=
  (finIcc a b).trans (eqvM M a b hM)

/-- `gM` is the site map `siteM`. -/
lemma gM_val (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) (i : Fin (ivlLen a b)) :
    ((gM M a b hM i : bulkM M a b) : Fin (2 * M + 1)) = siteM M a b i := by
  simp only [gM, Equiv.trans_apply, eqvM, finIcc, Equiv.coe_fn_mk, siteM]
  exact InfVolAux.idxZ_eq_toFin _ _ _

/-- `gM` is increasing. -/
lemma gM_strictMono (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    StrictMono (fun i => ((gM M a b hM i : bulkM M a b) : Fin (2 * M + 1))) := by
  intro i j hij
  simp only [gM, Equiv.trans_apply, eqvM, finIcc, Equiv.coe_fn_mk, idxZ]
  rw [Fin.mk_lt_mk]
  have h1 := Fin.lt_def.mp hij
  obtain ⟨h3, h4⟩ := hM
  omega

/-- Relabelling the configurations of `R` along `g : Fin L ≃ R`. -/
def cfgEquiv {N L : ℕ} {R : Finset (Fin N)} (g : Fin L ≃ R) : Cfg L ≃ (R → ZMod 2) where
  toFun n r := n (g.symm r)
  invFun α i := α (g i)
  left_inv n := by funext i; simp
  right_inv α := by funext r; simp

/-- **`lem:slater`(a) for any increasing labelling** of the region `R`. -/
lemma rdmF_eq_qfState_of_strictMono {N L : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) (R : Finset (Fin N)) (g : Fin L ≃ R)
    (hg : StrictMono (fun i => (g i : Fin N))) :
    (rdmF R (evo (dGamma A) t *ᵥ fockVec S)).submatrix (cfgEquiv g) (cfgEquiv g) =
      qfState ((corrR A S t R).submatrix g g) := by
  have hL : L = R.card := by
    have := Fintype.card_congr g
    simpa using this
  subst hL
  have h1 : ∀ i, g i = R.orderIsoOfFin rfl i := by
    have := Finset.orderEmbOfFin_unique (s := R) (h := rfl) (f := fun i => (g i : Fin N))
      (fun i => (g i).2) hg
    intro i
    apply Subtype.ext
    rw [Finset.coe_orderIsoOfFin_apply]
    exact congrFun this i
  have h2 : ∀ r, g.symm r = (R.orderIsoOfFin rfl).symm r := by
    intro r
    apply g.injective
    rw [Equiv.apply_symm_apply, h1, OrderIso.apply_symm_apply]
  have e1 : (⇑(cfgEquiv g)) = fun n : Cfg R.card => fun r : R => n ((R.orderIsoOfFin rfl).symm r) := by
    funext n r
    exact congrArg n (h2 r)
  have e2 : (⇑g) = fun i => R.orderIsoOfFin rfl i := funext h1
  rw [e1, e2]
  exact rdmF_eq_qfState A hA S t R

/-- `relabelX` is the relabelling along `gM`. -/
lemma relabelX_eq (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M)
    (X : Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    relabelX M a b X = X.submatrix (cfgEquiv (gM M a b hM)).symm (cfgEquiv (gM M a b hM)).symm := by
  have e : (fun (α : bulkM M a b → ZMod 2) (i : Fin (ivlLen a b)) =>
      if h : siteM M a b i ∈ bulkM M a b then α ⟨siteM M a b i, h⟩ else 0) =
        ⇑(cfgEquiv (gM M a b hM)).symm := by
    funext α i
    have hmem : siteM M a b i ∈ bulkM M a b := gM_val M a b hM i ▸ (gM M a b hM i).2
    rw [dite_eq_left hmem]
    show α ⟨siteM M a b i, hmem⟩ = α (gM M a b hM i)
    congr 1
    exact Subtype.ext (gM_val M a b hM i).symm
  unfold relabelX
  rw [e]

/-- `K^{(M)}_q(t)|_X` with the modes of `X` labelled by `Fin L`. -/
def KML (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ) :
    Matrix (Fin (ivlLen a b)) (Fin (ivlLen a b)) ℂ :=
  fun i j => KM J h q n t M (a + i.val) (a + j.val)

/-- For `[a, b] ⊆ Λ_M`, `C_R` of the bulk interval, relabelled along `gM`, is `K^{(M)}_q(t)|_X`. -/
lemma corrR_submatrix_gM (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    (corrR (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)).submatrix (gM M a b hM)
      (gM M a b hM) = KML J h q n t M a b := by
  have e : (⇑(gM M a b hM)) = ⇑(eqvM M a b hM) ∘ ⇑(finIcc a b) := rfl
  rw [e, ← submatrix_submatrix, corrR_submatrix_eqvM]
  rfl

/-- The trace is invariant under relabelling. -/
lemma trace_submatrix_equiv {ι κ : Type*} [Fintype ι] [Fintype κ] (A : Matrix ι ι ℂ) (e : κ ≃ ι) :
    (A.submatrix e e).trace = A.trace := by
  simp only [Matrix.trace, Matrix.diag, submatrix_apply]
  exact e.sum_comp (fun i => A i i)

/-- `Tr(ρ X') = Tr(ρ' X)` for the relabelled `ρ' = ρ ∘ e`, `X' = X ∘ e⁻¹`. -/
lemma trace_mul_submatrix_symm {ι κ : Type*} [Fintype ι] [Fintype κ] (ρ : Matrix ι ι ℂ)
    (X : Matrix κ κ ℂ) (e : κ ≃ ι) :
    (ρ * X.submatrix e.symm e.symm).trace = (ρ.submatrix e e * X).trace := by
  rw [← trace_submatrix_equiv (ρ * X.submatrix e.symm e.symm) e,
    ← submatrix_mul_equiv ρ (X.submatrix e.symm e.symm) e e e, submatrix_submatrix,
    Equiv.symm_comp_self, submatrix_id_id]

/-- **Finite volume** (paper: proof of `thm:infvol`(a)): for `[a, b] ⊆ Λ_M`,
`Tr(ρ_X(ψ^{(M)}_q(t)) O) = Tr(ρ[K^{(M)}_q(t)|_X] O)`. -/
lemma trace_rdmF_relabelX (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) (X : Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M q))) t *ᵥ fockVec (occM M n)) *
      relabelX M a b X).trace = (qfState (KML J h q n t M a b) * X).trace := by
  rw [relabelX_eq M a b hM X, trace_mul_submatrix_symm,
    rdmF_eq_qfState_of_strictMono _ (hq_isSymm J h _) _ t _ (gM M a b hM) (gM_strictMono M a b hM),
    corrR_submatrix_gM]

/-! ### `0 ≤ K ≤ 1` -/

section Order

open scoped ComplexOrder

/-- A diagonal `0/1` matrix `P` satisfies `P = P P†`. -/
lemma diag01_eq {ι : Type*} [Fintype ι] [DecidableEq ι] (p : ι → Prop) [DecidablePred p] :
    diagonal (fun z => if p z then (1 : ℂ) else 0) =
      diagonal (fun z => if p z then (1 : ℂ) else 0) *
        (diagonal (fun z => if p z then (1 : ℂ) else 0))ᴴ := by
  rw [diagonal_conjTranspose, diagonal_mul_diagonal]
  congr 1
  funext z
  simp only [Pi.star_apply]
  split_ifs <;> simp

/-- `0 ≤ C(t) ≤ 1` for the correlation matrix of an evolved Fock state. -/
lemma corrMat_mem {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) : (corrMat A S t).PosSemidef ∧ (1 - corrMat A S t).PosSemidef := by
  have hP := diag01_eq (fun z => z ∈ S)
  have hQ := diag01_eq (fun z => z ∉ S)
  have hPQ : (1 : Matrix (Fin N) (Fin N) ℂ) - diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) =
      diagonal (fun z => if z ∉ S then (1 : ℂ) else 0) := by
    rw [← diagonal_one, diagonal_sub]
    congr 1
    funext z
    split_ifs <;> simp_all
  have hUU : evo A t * (evo A t)ᴴ = 1 := evo_unitary' A hA t
  have hC : corrMat A S t = evo A t * diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) *
      (evo A t)ᴴ := rfl
  set U := evo A t
  set P := diagonal (fun z => if z ∈ S then (1 : ℂ) else 0)
  set Q := diagonal (fun z => if z ∉ S then (1 : ℂ) else 0)
  constructor
  · have e : corrMat A S t = (U * P) * (U * P)ᴴ := by
      calc corrMat A S t = U * (P * Pᴴ) * Uᴴ := by rw [← hP, hC]
        _ = (U * P) * (U * P)ᴴ := by rw [conjTranspose_mul]; simp only [Matrix.mul_assoc]
    rw [e]
    exact posSemidef_self_mul_conjTranspose _
  · have e : 1 - corrMat A S t = (U * Q) * (U * Q)ᴴ := by
      calc 1 - corrMat A S t = U * (1 - P) * Uᴴ := by
            rw [hC, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, hUU]
        _ = U * (Q * Qᴴ) * Uᴴ := by rw [hPQ, ← hQ]
        _ = (U * Q) * (U * Q)ᴴ := by rw [conjTranspose_mul]; simp only [Matrix.mul_assoc]
    rw [e]
    exact posSemidef_self_mul_conjTranspose _

/-- `0 ≤ C ≤ 1` passes to principal submatrices. -/
lemma mem_submatrix {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]
    (C : Matrix ι ι ℂ) (e : κ → ι) (he : Function.Injective e)
    (hC : C.PosSemidef ∧ (1 - C).PosSemidef) :
    (C.submatrix e e).PosSemidef ∧ (1 - C.submatrix e e).PosSemidef := by
  refine ⟨hC.1.submatrix _, ?_⟩
  have e1 : 1 - C.submatrix e e = (1 - C).submatrix e e := by
    ext i j
    simp [Matrix.one_apply, he.eq_iff]
  rw [e1]
  exact hC.2.submatrix _

/-- For `[a, b] ⊆ Λ_M`, `0 ≤ K^{(M)}_q(t)|_X ≤ 1`. -/
lemma KML_mem (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    (KML J h q n t M a b).PosSemidef ∧ (1 - KML J h q n t M a b).PosSemidef := by
  rw [← corrR_submatrix_gM J h q n t M a b hM]
  refine mem_submatrix _ _ (gM M a b hM).injective ?_
  exact mem_submatrix _ _ Subtype.val_injective (corrMat_mem _ (hq_isSymm J h _) _ t)

/-- `K^{(M)}_q(t)|_X → K^ℤ_q(t)|_X` for fixed `q` and `t`. -/
lemma tendsto_KML (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) :
    Tendsto (fun M => KML J h q n t M a b) atTop (𝓝 (KZL J h q n t a b)) :=
  tendsto_pi_nhds.2 fun _ => tendsto_pi_nhds.2 fun _ => tendsto_KM J h q n t _ _

/-- The set `{0 ≤ C ≤ 1}` is closed. -/
lemma isClosed_mem (L : ℕ) :
    IsClosed {C : Matrix (Fin L) (Fin L) ℂ | C.PosSemidef ∧ (1 - C).PosSemidef} :=
  isClosed_setOfPred_posSemidef.inter
    (isClosed_setOfPred_posSemidef.preimage (continuous_const.sub continuous_id))

/-- The set `{0 ≤ C ≤ 1}` is compact. -/
lemma isCompact_mem (L : ℕ) :
    IsCompact {C : Matrix (Fin L) (Fin L) ℂ | C.PosSemidef ∧ (1 - C).PosSemidef} := by
  have hK : IsCompact (Set.univ.pi fun _ : Fin L => Set.univ.pi fun _ : Fin L =>
      Metric.closedBall (0 : ℂ) L) :=
    isCompact_univ_pi fun _ => isCompact_univ_pi fun _ => isCompact_closedBall _ _
  refine hK.of_isClosed_subset (isClosed_mem L) ?_
  rintro C ⟨hC0, hC1⟩
  simp only [Set.mem_pi, Set.mem_univ, true_implies, Metric.mem_closedBall, dist_zero_right]
  intro i j
  exact QuasiFreeAux.norm_entry_le hC0 hC1 i j

/-- `0 ≤ K^ℤ_q(t)|_X ≤ 1`. -/
lemma KZL_mem (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) :
    (KZL J h q n t a b).PosSemidef ∧ (1 - KZL J h q n t a b).PosSemidef := by
  refine (isClosed_mem _).mem_of_tendsto (tendsto_KML J h q n t a b) ?_
  filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
  exact KML_mem J h q n t M a b ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩

/-- `Tr(ρ[C] O_c)` is bounded on `{0 ≤ C ≤ 1}`, uniformly in the finitely many `c`. -/
lemma exists_bound_trace (L : ℕ) (O : Cfg L → Matrix (Cfg L) (Cfg L) ℂ) :
    ∃ B : ℝ, ∀ C : Matrix (Fin L) (Fin L) ℂ, C.PosSemidef → (1 - C).PosSemidef →
      ∀ c, ‖(qfState C * O c).trace‖ ≤ B := by
  have hcont : Continuous (fun C : Matrix (Fin L) (Fin L) ℂ => ∑ c, ‖(qfState C * O c).trace‖) :=
    continuous_finsetSum _ fun c _ =>
      ((continuous_qfState.matrix_mul continuous_const).matrix_trace).norm
  obtain ⟨B, hB⟩ := (isCompact_mem L).exists_bound_of_continuousOn hcont.continuousOn
  refine ⟨B, fun C h0 h1 c => ?_⟩
  have h2 := hB C ⟨h0, h1⟩
  rw [Real.norm_eq_abs, abs_of_nonneg (Finset.sum_nonneg fun _ _ => norm_nonneg _)] at h2
  exact (Finset.single_le_sum (f := fun c => ‖(qfState C * O c).trace‖)
    (fun _ _ => norm_nonneg _) (Finset.mem_univ c)).trans h2

/-! ### Convergence uniformly in `q` and in `t ∈ [-T, T]` -/

/-- Charges and times with `|t| ≤ T`. -/
abbrev CT (T : ℝ) := ↥({p : CfgZ × ℝ | |p.2| ≤ T} : Set (CfgZ × ℝ))

/-- The propagation weight is monotone in `c`. -/
lemma wt_mono {c c' : ℝ} (hc : 0 ≤ c) (hcc' : c ≤ c') (d : ℕ) : wt c d ≤ wt c' d := by
  have hc' : 0 ≤ c' := hc.trans hcc'
  unfold wt
  have h1 : 0 ≤ c' ^ d / (d.factorial : ℝ) := by positivity
  gcongr

/-- `c_t ≤ c_T` for `|t| ≤ T`. -/
lemma ct_le (T J h : ℝ) (p : CT T) :
    |p.1.2| * (2 * |J| + 2 * |h|) ≤ |T| * (2 * |J| + 2 * |h|) := by
  have h1 : |p.1.2| ≤ T := p.2
  have h2 : |p.1.2| ≤ |T| := h1.trans (le_abs_self T)
  gcongr

/-- **Uniform convergence of the propagator entries, uniformly in `q` and `|t| ≤ T`.** -/
lemma tendstoUniformly_evoZ_T (J h T : ℝ) (x s : ℤ) :
    TendstoUniformly (fun (M : ℕ) (p : CT T) => evoZ J h p.1.1 p.1.2 M x s)
      (fun p => UZ J h p.1.1 p.1.2 x s) atTop := by
  set c := |T| * (2 * |J| + 2 * |h|) with hc
  have hlim : Tendsto (fun M : ℕ => 2 * (c ^ (M - x.natAbs + 1) /
      ((M - x.natAbs + 1).factorial : ℝ) * Real.exp c)) atTop (𝓝 0) := by
    have h0 := FloorSemiring.tendsto_pow_div_factorial_atTop c
    have hm : Tendsto (fun M : ℕ => M - x.natAbs + 1) atTop atTop :=
      (tendsto_add_atTop_nat 1).comp (tendsto_sub_atTop_nat x.natAbs)
    simpa using ((h0.comp hm).mul_const (Real.exp c)).const_mul 2
  rw [Metric.tendstoUniformly_iff]
  intro ε hε
  filter_upwards [hlim.eventually (gt_mem_nhds hε), eventually_ge_atTop (x.natAbs + s.natAbs)]
    with M hM1 hM2 p
  have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  have hs : |s| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  rw [dist_comm, dist_eq_norm]
  refine ((norm_evoZ_sub_UZ_le J h p.1.1 p.1.2 M x s hx hs).trans ?_).trans_lt hM1
  have hw := wt_mono (by positivity) (ct_le T J h p) (M - x.natAbs + 1)
  unfold wt at hw
  linarith

/-- **`K^{(M)}_q(t)_{xy} → K^ℤ_q(t)_{xy}` uniformly in `q` and `|t| ≤ T`** (paper: Step 1 of the
proof of `thm:infvol`, "the bound is also uniform for `|t|` in compact sets"). -/
lemma tendstoUniformly_KM_T (J h T : ℝ) (n : CfgZ) (x y : ℤ) :
    TendstoUniformly (fun (M : ℕ) (p : CT T) => KM J h p.1.1 n p.1.2 M x y)
      (fun p => KZ J h p.1.1 n p.1.2 x y) atTop := by
  set c := |T| * (2 * |J| + 2 * |h|) with hc'
  have hc : 0 ≤ c := by positivity
  unfold KM KZ
  refine tendstoUniformly_tsum
    (f := fun M (p : CT T) s => if n s = 1 then
      evoZ J h p.1.1 p.1.2 M x s * star (evoZ J h p.1.1 p.1.2 M y s) else 0)
    (g := fun (p : CT T) s => if n s = 1 then
      UZ J h p.1.1 p.1.2 x s * star (UZ J h p.1.1 p.1.2 y s) else 0)
    (b := fun s => wt c (x - s).natAbs) (summable_wt_int c x) (fun M p s => ?_) (fun p s => ?_)
    (fun s => ?_)
  · split_ifs
    · rw [norm_mul, norm_star]
      calc ‖evoZ J h p.1.1 p.1.2 M x s‖ * ‖evoZ J h p.1.1 p.1.2 M y s‖
          ≤ wt c (x - s).natAbs * 1 :=
            mul_le_mul ((norm_evoZ_le_wt J h p.1.1 p.1.2 M x s).trans
                (wt_mono (by positivity) (ct_le T J h p) _))
              (InfVolAux.norm_evoZ_le_one J h p.1.1 p.1.2 M y s) (norm_nonneg _) (wt_nonneg hc _)
        _ = _ := mul_one _
    · rw [norm_zero]; exact wt_nonneg hc _
  · split_ifs
    · rw [norm_mul, norm_star]
      calc ‖UZ J h p.1.1 p.1.2 x s‖ * ‖UZ J h p.1.1 p.1.2 y s‖ ≤ wt c (x - s).natAbs * 1 :=
            mul_le_mul ((norm_UZ_le_wt J h p.1.1 p.1.2 x s).trans
                (wt_mono (by positivity) (ct_le T J h p) _))
              (norm_UZ_le_one J h p.1.1 p.1.2 y s) (norm_nonneg _) (wt_nonneg hc _)
        _ = _ := mul_one _
    · rw [norm_zero]; exact wt_nonneg hc _
  · by_cases hs : n s = 1
    · simp only [hs, ite_true]
      exact tendstoUniformly_mul_star (tendstoUniformly_evoZ_T J h T x s)
        (tendstoUniformly_evoZ_T J h T y s)
        (fun M p => InfVolAux.norm_evoZ_le_one J h p.1.1 p.1.2 M x s)
        (fun p => norm_UZ_le_one J h p.1.1 p.1.2 y s)
    · simp only [hs, ite_false]
      exact tendsto_const_nhds.tendstoUniformly_const

/-- `K^{(M)}_q(t)|_X → K^ℤ_q(t)|_X` uniformly in `q` and `|t| ≤ T`. -/
lemma tendstoUniformly_KML_T (J h T : ℝ) (n : CfgZ) (a b : ℤ) :
    TendstoUniformly (fun (M : ℕ) (p : CT T) => KML J h p.1.1 n p.1.2 M a b)
      (fun p => KZL J h p.1.1 n p.1.2 a b) atTop :=
  tendstoUniformly_matrix fun _ _ => tendstoUniformly_KM_T J h T n _ _

/-- `X ↦ Tr(X B)` is uniformly continuous. -/
lemma uniformContinuous_trace_mul {ι : Type*} [Fintype ι] [DecidableEq ι] (B : Matrix ι ι ℂ) :
    UniformContinuous (fun X : Matrix ι ι ℂ => (X * B).trace) := by
  let φ : Matrix ι ι ℂ →+ ℂ := (Matrix.traceAddMonoidHom ι ℂ).comp (AddMonoidHom.mulRight B)
  exact uniformContinuous_addMonoidHom_of_continuous (f := φ)
    ((continuous_id.matrix_mul continuous_const).matrix_trace)

/-- **`Tr(ρ[K^{(M)}_q(t)|_X] O_{q|_X}) → Tr(ρ[K^ℤ_q(t)|_X] O_{q|_X})` uniformly in `q` and
`|t| ≤ T`** (uniform continuity of `C ↦ ρ[C]` on `{0 ≤ C ≤ 1}`, finitely many `O_c`). -/
lemma tendstoUniformlyOn_trace_qfState (J h T : ℝ) (n : CfgZ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    TendstoUniformlyOn
      (fun (M : ℕ) (p : CfgZ × ℝ) =>
        (qfState (KML J h p.1 n p.2 M a b) * O (fun i => p.1 (a + i.val))).trace)
      (fun p => (qfState (KZL J h p.1 n p.2 a b) * O (fun i => p.1 (a + i.val))).trace)
      atTop {p | |p.2| ≤ T} := by
  rw [tendstoUniformlyOn_iff_tendstoUniformly_comp_coe]
  have hq : TendstoUniformly (fun (M : ℕ) (p : CT T) => qfState (KML J h p.1.1 n p.1.2 M a b))
      (fun p => qfState (KZL J h p.1.1 n p.1.2 a b)) atTop := by
    refine uniformContinuousOn_qfState.comp_tendstoUniformly_eventually ?_
      (fun p => KZL_mem J h p.1.1 n p.1.2 a b) (tendstoUniformly_KML_T J h T n a b)
    filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM p
    exact KML_mem J h p.1.1 n p.1.2 M a b ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
  rw [Metric.tendstoUniformly_iff]
  intro ε hε
  have hall : ∀ c : Cfg (ivlLen a b), ∀ᶠ M in atTop, ∀ p : CT T,
      dist ((qfState (KZL J h p.1.1 n p.1.2 a b) * O c).trace)
        ((qfState (KML J h p.1.1 n p.1.2 M a b) * O c).trace) < ε := fun c =>
    Metric.tendstoUniformly_iff.1 ((uniformContinuous_trace_mul (O c)).comp_tendstoUniformly hq) ε hε
  filter_upwards [Filter.eventually_all.2 hall] with M hM p
  exact hM _ p

end Order

/-! ### The finite-volume expectation as a sector average -/

/-- `locObsM` is gauge invariant. -/
lemma gaugeInv_locObsM (M : ℕ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) :
    GaugeInv (locObsM M a b O) := by
  intro j
  rw [sec_cplx_gauss]
  ext x y
  rw [mul_diagonal, diagonal_mul]
  unfold locObsM
  rw [submatrix_apply, blockDiagonal_apply]
  split_ifs with hq
  · have hx := sec_gauss_dual_apply j (dualEquiv.symm x)
    have hy := sec_gauss_dual_apply j (dualEquiv.symm y)
    rw [Equiv.apply_symm_apply] at hx hy
    rw [hx, hy, hq]
    ring
  · simp

/-- The sector blocks of `locObsM` are the embedded `O_{q|_X}`. -/
lemma sectorBlock_locObsM (M : ℕ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) (q : Cfg (2 * M + 1)) :
    sectorBlock (locObsM M a b O) q = embF (bulkM M a b) (relabelX M a b (O (chargesX M a b q))) := by
  ext n' n
  simp only [sectorBlock, locObsM, submatrix_apply, Equiv.symm_apply_apply, blockDiagonal_apply_eq]

/-- **Sector averaging for `locObsM`** (paper: `cor:sector`(b) and `lem:slater`(a)):
`⟨Ψ^{(M)}(t), O Ψ^{(M)}(t)⟩ = 2^{-N} ∑_q Tr(ρ_X(ψ_q(t)) O_{q|_X})`. -/
lemma expect_locObsM (J h : ℝ) (M : ℕ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ)
    (ψ : Cfg (2 * M + 1) → ℂ) (t : ℝ) :
    expect (locObsM M a b O) (lgtState J h ψ t) = (1 / 2 ^ (2 * M + 1) : ℂ) *
      ∑ q : Cfg (2 * M + 1), (rdmF (bulkM M a b) (evo (dGamma (hq J h q)) t *ᵥ ψ) *
        relabelX M a b (O (chargesX M a b q))).trace := by
  rw [sector_average J h ψ _ (gaugeInv_locObsM M a b O) t]
  refine congrArg (fun s => (1 / 2 ^ (2 * M + 1) : ℂ) * s) (Finset.sum_congr rfl fun q _ => ?_)
  rw [sectorBlock_locObsM, ← expect_embF]
  unfold Hsec
  exact expect_evo_sub_smul_one _ _ _ _ _

/-- Integrals of complex functions of `q|_{Λ_M}` are the uniform averages. -/
lemma integral_comp_restrictZ_cplx (M : ℕ) (g : Cfg (2 * M + 1) → ℂ) :
    ∫ q, g (restrictZ M q) ∂muZ = (1 / 2 ^ (2 * M + 1) : ℂ) * ∑ q', g q' := by
  rw [← integral_map (InfVolAux.measurable_restrictZ M).aemeasurable
    (measurable_of_countable g).aestronglyMeasurable, integral_fintype Integrable.of_finite,
    Finset.mul_sum]
  refine Finset.sum_congr rfl fun q' _ => ?_
  rw [Measure.real, Measure.map_apply (InfVolAux.measurable_restrictZ M) (measurableSet_singleton q'),
    InfVolAux.muZ_restrictZ_preimage, ENNReal.toReal_inv, ENNReal.toReal_pow, Complex.real_smul]
  norm_num

/-- For `[a, b] ⊆ Λ_M`, the charges of `q|_{Λ_M}` on `X` are `i ↦ q (a + i)`. -/
lemma chargesX_restrictZ (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) (q : CfgZ) :
    chargesX M a b (restrictZ M q) = fun i => q (a + i.val) := by
  funext i
  simp only [chargesX, restrictZ, siteM]
  have hi := i.isLt
  unfold ivlLen at hi
  obtain ⟨h3, h4⟩ := hM
  have habs : |a + (i.val : ℤ)| ≤ M := by rw [abs_le]; constructor <;> omega
  rw [InfVolAux.toFin_val M _ habs]
  congr 1
  ring

end InfVolDynAux

/-- **Limit dynamics** (paper: Theorem `thm:infvol`(a)): `⟨Ψ^{(M)}(t), O Ψ^{(M)}(t)⟩ → E Tr(ρ[K^ℤ_q(t)|_X] O_{q|_X})`. -/
theorem tendsto_omega (J h : ℝ) (n : CfgZ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) (t : ℝ) :
    Tendsto (fun M : ℕ => expect (locObsM M a b O) (lgtState J h (fockVec (occM M n)) t)) atTop
      (𝓝 (∫ q, (qfState (KZL J h q n t a b) * O (fun i => q (a + i.val))).trace ∂muZ)) := by
  have h1 : ∀ M : ℕ, expect (locObsM M a b O) (lgtState J h (fockVec (occM M n)) t) =
      ∫ q, (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M q))) t *ᵥ fockVec (occM M n)) *
        relabelX M a b (O (chargesX M a b (restrictZ M q)))).trace ∂muZ := by
    intro M
    rw [InfVolDynAux.expect_locObsM]
    exact (InfVolDynAux.integral_comp_restrictZ_cplx M (fun q' => (rdmF (bulkM M a b)
      (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)) *
        relabelX M a b (O (chargesX M a b q'))).trace)).symm
  refine Tendsto.congr (fun M => (h1 M).symm) ?_
  obtain ⟨B, hB⟩ := InfVolDynAux.exists_bound_trace (ivlLen a b) O
  refine tendsto_integral_filter_of_dominated_convergence (fun _ => B)
    (Eventually.of_forall fun M => (InfVolAux.measurable_comp_restrictZ M (fun q' =>
      (rdmF (bulkM M a b) (evo (dGamma (hq J h q')) t *ᵥ fockVec (occM M n)) *
        relabelX M a b (O (chargesX M a b q'))).trace)).aestronglyMeasurable)
    ?_ (integrable_const _) (ae_of_all _ fun q => ?_)
  · filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
    have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
    refine ae_of_all _ fun q => ?_
    rw [InfVolDynAux.chargesX_restrictZ M a b hM' q,
      InfVolDynAux.trace_rdmF_relabelX J h q n t M a b hM']
    obtain ⟨h0, h1⟩ := InfVolDynAux.KML_mem J h q n t M a b hM'
    exact hB _ h0 h1 _
  · have hc : Continuous (fun C => (qfState C * O (fun i => q (a + i.val))).trace) :=
      (continuous_qfState.matrix_mul continuous_const).matrix_trace
    refine ((hc.tendsto _).comp (InfVolDynAux.tendsto_KML J h q n t a b)).congr' ?_
    filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
    have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
    rw [Function.comp_apply, InfVolDynAux.chargesX_restrictZ M a b hM' q,
      InfVolDynAux.trace_rdmF_relabelX J h q n t M a b hM']

/-- **Sector by sector, uniformly in `q` and in `t ∈ [-T, T]`** (paper: Theorem `thm:infvol`(a)). -/
theorem tendstoUniformlyOn_sector (J h : ℝ) (n : CfgZ) (a b : ℤ)
    (O : Cfg (ivlLen a b) → Matrix (Cfg (ivlLen a b)) (Cfg (ivlLen a b)) ℂ) (T : ℝ) :
    TendstoUniformlyOn
      (fun (M : ℕ) (p : CfgZ × ℝ) =>
        (rdmF (bulkM M a b) (evo (dGamma (hq J h (restrictZ M p.1))) p.2 *ᵥ fockVec (occM M n)) *
          relabelX M a b (O (fun i => p.1 (a + i.val)))).trace)
      (fun p : CfgZ × ℝ => (qfState (KZL J h p.1 n p.2 a b) * O (fun i => p.1 (a + i.val))).trace)
      atTop {p | |p.2| ≤ T} := by
  refine (InfVolDynAux.tendstoUniformlyOn_trace_qfState J h T n a b O).congr ?_
  filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM p _
  exact (InfVolDynAux.trace_rdmF_relabelX J h p.1 n p.2 M a b
    ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩ _).symm

end DFL
