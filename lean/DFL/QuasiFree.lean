import DFL.SlaterBulk

/-!
# Quasi-free density matrices (paper: Lemma `lem:quasifree`, and `lem:slater`(a): `ρ_A = ρ[C_A]`)

For a matrix `C` on `ℓ²(Fin L)` let `Γ(C)` be its second quantization on the Fock space `𝓕_L`,
`Γ(C) c†_{x_1} ⋯ c†_{x_k} |0⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` (`secondQ`), and (eq. `eq:main-rhoformula`)
`ρ[C] = ∑_{j=0}^{L} ((-1)^j / j!) ∑_{x_1, …, x_j} c_{x_1} ⋯ c_{x_j} Γ(C) c†_{x_j} ⋯ c†_{x_1}` (`qfState`).

* (a) For Hermitian `C` with orthonormal eigenvectors `v_i` and eigenvalues `ν_i`,
  `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))` (`qfState_eq_prod`): `ρ[C]` does not depend on the
  eigenbasis; its matrix elements are polynomials in the entries of `C` (`continuous_qfState`), and `C ↦ ρ[C]` is
  uniformly continuous on `{0 ≤ C ≤ 1}` (`uniformContinuousOn_qfState`).
* (b) For `0 ≤ C ≤ 1`, `ρ[C]` is a density matrix (`qfState_posSemidef`, `trace_qfState`) with
  `S(ρ[C]) = ∑_i h(ν_i)` (`vonNeumann_qfState`); its generating function is
  `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)` (`trace_qfState_zN`).
* `lem:slater`(a): the reduced state of a Slater determinant on any region `R` is the quasi-free state with symbol
  `C_R` (`rdmF_eq_qfState`), after relabelling the modes of `R` increasingly by `Fin |R|`.

## Proof outline (auxiliary lemmas in `DFL.QuasiFreeAux`)

1. `secondQ_mulVec_slaterAlt`: `Γ(C) c†(w_1) ⋯ c†(w_m)|0⟩ = c†(Cw_1) ⋯ c†(Cw_m)|0⟩` (both sides are alternating
   maps agreeing on sorted basis tuples), so `Γ(C)` multiplies a Slater determinant of eigenvectors by the product
   of the eigenvalues.
2. `sum_list_eq_xiOp`, `sum_ann_cre_eq`: the inner sums of `ρ[C]` are `Ξ_j = Ξ_1^j` with
   `Ξ_1(Y) = ∑_x c_x Y c†_x = ∑_i c(v_i) Y c†(v_i)` for any orthonormal basis `v`; on a Slater determinant `|s⟩` of
   eigenvectors, `Ξ_j(Γ(C))|s⟩ = ∏_k ν_{s_k} S_j(s)|s⟩` (`xiOp_mulVec_slat`), where `S_j(U)` (`sFun`) is the `j`-th
   derivative at `0` of `∏_{i ∉ U} (1 + ν_i X)`, so Taylor's formula gives the inclusion–exclusion identity
   `∑_j (-1)^j/j! S_j(U) = ∏_{i ∉ U} (1 - ν_i)` (`sum_sFun`).
3. `qfState_eq_conj`: in the orthonormal Slater basis `|T⟩` of any orthonormal eigenbasis, `ρ[C] = W diag(w) W†`
   with `w_T = ∏_{i ∈ T} ν_i ∏_{i ∉ T} (1 - ν_i)`; this gives (a), (b) (product law, `sum_negMulLog_prod_sel`).
4. `rdmF_slater_sum`, `prodSel_embR`: `ρ_R = ∑_σ μ_σ |X_σ⟩⟨X_σ|` (from `SlaterBulkAux`), transported to `Fin |R|`
   (the increasing relabelling preserves the Jordan–Wigner strings); `qfState_eq_sum_prodSel` identifies this
   mixture with `ρ[C_R]` in an eigenbasis of `C_R` adapted to the canonical orbitals (`exists_adapted_basis`, an
   orthonormal extension of the normalized nonzero restricted orbitals).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

open AreaLawAux

variable {L : ℕ}

/-- The occupied sites of `n`. -/
def occSet (n : Cfg L) : Finset (Fin L) := univ.filter (fun x => n x = 1)

/-- **Second quantization** `Γ(C)`: `Γ(C)|n⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` for the occupied sites
`x_1 < ⋯ < x_k` of `n`. -/
def secondQ (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ := fun n' n =>
  slaterAlt L (occSet n).card
    (fun j : Fin (occSet n).card => fun x => C x ((occSet n).orderEmbOfFin rfl j)) n'

/-- **The quasi-free density matrix** `ρ[C]` (paper: eq. `eq:main-rhoformula`). -/
def qfState (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ :=
  ∑ j ∈ range (L + 1), ((-1 : ℂ) ^ j / (j.factorial : ℂ)) •
    ∑ x : Fin j → Fin L,
      (List.ofFn fun i => cplx (ann (x i))).prod * secondQ C *
        (List.ofFn fun i => cplx (cre (x (Fin.rev i)))).prod

/-- The operator `z^{N}` on `𝓕_L`. -/
def zPowN (z : ℂ) : Matrix (Cfg L) (Cfg L) ℂ := diagonal (fun n => z ^ (occSet n).card)

namespace QuasiFreeAux

open SlaterExactAux SlaterGeneralAux

/-! ### Second quantization of Slater determinants -/

/-- `c†(δ_x) = c_x†`. -/
lemma crOp_single (x : Fin L) : crOp (Pi.single x (1 : ℂ)) = cplx (cre x) := by
  simp [crOp, Pi.single_apply]

/-- The occupied set of the Fock configuration of `S` is `S`. -/
lemma occSet_fockCfg (S : Finset (Fin L)) : occSet (fockCfg S) = S := by
  ext x; simp [occSet, FreeDynAux.fockCfg_eq_one]

/-- The image of `Fin.cons a s`. -/
lemma image_cons {α : Type*} [DecidableEq α] {m : ℕ} (a : α) (s : Fin m → α) :
    univ.image (Fin.cons a s : Fin (m + 1) → α) = insert a (univ.image s) := by
  ext x
  simp only [mem_image, mem_univ, true_and, mem_insert, Fin.exists_fin_succ, Fin.cons_zero,
    Fin.cons_succ]
  constructor
  · rintro (h | h)
    · exact Or.inl h.symm
    · exact Or.inr h
  · rintro (h | h)
    · exact Or.inl h.symm
    · exact Or.inr h

/-- A Slater determinant of increasing basis orbitals is a Fock basis vector. -/
lemma slaterAlt_single_of_strictMono : ∀ {m : ℕ} (u : Fin m → Fin L), StrictMono u →
    slaterAlt L m (fun j => Pi.single (u j) (1 : ℂ)) = fockVec (univ.image u)
  | 0, u, _ => by simp [slaterAlt_apply]
  | m + 1, u, hu => by
    have hu' : StrictMono (fun k : Fin m => u k.succ) := hu.comp (Fin.strictMono_succ)
    have himg : univ.image u = insert (u 0) (univ.image (fun k : Fin m => u k.succ)) := by
      have h := image_cons (u 0) (fun k : Fin m => u k.succ)
      have e : (Fin.cons (u 0) (fun k : Fin m => u k.succ) : Fin (m + 1) → Fin L) = u :=
        Fin.cons_self_tail u
      rwa [e] at h
    rw [slaterAlt_succ, slaterAlt_single_of_strictMono _ hu', crOp_single, himg,
      AreaLawAux.fockVec_insert]
    intro x hx
    obtain ⟨k, -, rfl⟩ := mem_image.mp hx
    exact hu (Fin.succ_pos k)

/-- The column of `Γ(C)` at `n`, written via any finset equal to the occupied set of `n`. -/
lemma secondQ_apply_of_eq (C : Matrix (Fin L) (Fin L) ℂ) (n' n : Cfg L) (T : Finset (Fin L))
    (hT : occSet n = T) :
    secondQ C n' n = slaterAlt L T.card (fun j x => C x (T.orderEmbOfFin rfl j)) n' := by
  subst hT; rfl

/-- `Γ(C) |S⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩`. -/
lemma secondQ_mulVec_fockVec (C : Matrix (Fin L) (Fin L) ℂ) (S : Finset (Fin L)) :
    secondQ C *ᵥ fockVec S = slaterAlt L S.card (fun j x => C x (S.orderEmbOfFin rfl j)) := by
  funext n'
  rw [fockVec, Matrix.mulVec_single_one]
  exact secondQ_apply_of_eq C n' _ S (occSet_fockCfg S)

/-- Changing the length argument of `orderEmbOfFin` along an equality. -/
lemma slaterAlt_orderEmb_cast (C : Matrix (Fin L) (Fin L) ℂ) (T : Finset (Fin L)) {m : ℕ}
    (h : T.card = m) :
    slaterAlt L T.card (fun j x => C x (T.orderEmbOfFin rfl j)) =
      slaterAlt L m (fun j x => C x (T.orderEmbOfFin h j)) := by
  subst h; rfl

/-- **Multiplicativity of the second quantization**:
`Γ(C) c†(w_1) ⋯ c†(w_m) |0⟩ = c†(C w_1) ⋯ c†(C w_m) |0⟩`. -/
lemma secondQ_mulVec_slaterAlt (C : Matrix (Fin L) (Fin L) ℂ) {m : ℕ} (w : Fin m → Fin L → ℂ) :
    secondQ C *ᵥ slaterAlt L m w = slaterAlt L m (fun j => C *ᵥ w j) := by
  have hFG : (secondQ C).mulVecLin.compAlternatingMap (slaterAlt L m) =
      (slaterAlt L m).compLinearMap C.mulVecLin := by
    refine Module.Basis.ext_alternating (Pi.basisFun ℂ (Fin L)) fun v hv => ?_
    set σ := Tuple.sort v
    have hmono : StrictMono (v ∘ σ) :=
      (Tuple.monotone_sort v).strictMono_of_injective (hv.comp σ.injective)
    rw [AlternatingMap.map_congr_perm _ _ σ,
      AlternatingMap.map_congr_perm ((slaterAlt L m).compLinearMap C.mulVecLin) _ σ]
    congr 1
    have hc : (univ.image (v ∘ σ)).card = m := by
      rw [card_image_of_injective _ hmono.injective, card_univ, Fintype.card_fin]
    have hs := Finset.orderEmbOfFin_unique hc (fun i => mem_image_of_mem (v ∘ σ) (mem_univ i))
      hmono
    have e1 : ((fun i => (Pi.basisFun ℂ (Fin L)) (v i)) ∘ σ) =
        fun j => Pi.single ((v ∘ σ) j) (1 : ℂ) := by
      funext j; simp
    rw [LinearMap.compAlternatingMap_apply, AlternatingMap.compLinearMap_apply, e1,
      slaterAlt_single_of_strictMono (v ∘ σ) hmono, Matrix.mulVecLin_apply,
      secondQ_mulVec_fockVec, slaterAlt_orderEmb_cast C _ hc, ← hs]
    congr 1
    funext j x
    simp
  exact congrArg (fun f => f w) hFG


/-! ### Slater determinants of an orthonormal family of eigenvectors -/

section Eig

variable (v : Fin L → Fin L → ℂ)

/-- The Slater determinant `c†(v_{s_0}) ⋯ c†(v_{s_{m-1}}) |0⟩` of orbitals of the family `v`. -/
def slat {m : ℕ} (s : Fin m → Fin L) : Cfg L → ℂ := slaterAlt L m (fun j => v (s j))

/-- `c†(v_a)` prepends the orbital `v_a`. -/
lemma crOp_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    crOp (v a) *ᵥ slat v s = slat v (Fin.cons a s : Fin (m + 1) → Fin L) := by
  rw [slat, slat, slaterAlt_succ]
  simp only [Fin.cons_zero, Fin.cons_succ]

/-- A Slater determinant with a repeated orbital vanishes. -/
lemma slat_cons_eq_zero {m : ℕ} (s : Fin m → Fin L) (a : Fin L) (ha : a ∈ univ.image s) :
    slat v (Fin.cons a s : Fin (m + 1) → Fin L) = 0 := by
  obtain ⟨k, -, hk⟩ := mem_image.mp ha
  exact (slaterAlt L (m + 1)).map_eq_zero_of_eq _ (i := 0) (j := k.succ) (by simp [hk])
    (Fin.succ_ne_zero k).symm

/-- `Γ(C)` acts on a Slater determinant of eigenvectors by the product of the eigenvalues. -/
lemma secondQ_mulVec_slat (C : Matrix (Fin L) (Fin L) ℂ) (ν : Fin L → ℂ)
    (hCv : ∀ i, C *ᵥ v i = ν i • v i) {m : ℕ} (s : Fin m → Fin L) :
    secondQ C *ᵥ slat v s = (∏ k, ν (s k)) • slat v s := by
  rw [slat, secondQ_mulVec_slaterAlt]
  simp only [hCv]
  exact AlternatingMap.map_smul_univ _ _ _

variable {v}
variable (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
include hv

/-- `c(v_a)` kills a Slater determinant of the other orbitals. -/
lemma adj_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) (ha : a ∉ univ.image s) :
    (crOp (v a))ᴴ *ᵥ slat v s = 0 :=
  adj_slaterAlt _ _ fun k => by
    rw [hv, ite_eq_right_iff]
    intro h
    exact absurd (mem_image.mpr ⟨k, mem_univ k, h.symm⟩) ha

/-- CAR for a normalized orbital. -/
lemma car_self (a : Fin L) :
    (crOp (v a))ᴴ * crOp (v a) + crOp (v a) * (crOp (v a))ᴴ = 1 := by
  rw [crOp_CAR_gen, hv, ite_eq_left rfl, one_smul]

/-- `c(v_a) c†(v_a) |s⟩ = [a ∉ s] |s⟩`, written with the prepended orbital. -/
lemma adj_mulVec_slat_cons {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    (crOp (v a))ᴴ *ᵥ slat v (Fin.cons a s : Fin (m + 1) → Fin L) =
      if a ∈ univ.image s then 0 else slat v s := by
  split_ifs with ha
  · rw [slat_cons_eq_zero v s a ha, mulVec_zero]
  · rw [← crOp_mulVec_slat, mulVec_mulVec]
    have h := car_self hv a
    rw [← eq_sub_iff_add_eq] at h
    rw [h, sub_mulVec, one_mulVec, ← mulVec_mulVec, adj_mulVec_slat hv s a ha, mulVec_zero,
      sub_zero]

/-- The number operator `n(v_a) = c†(v_a) c(v_a)` on a Slater determinant of the family. -/
lemma num_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    (crOp (v a) * (crOp (v a))ᴴ) *ᵥ slat v s = if a ∈ univ.image s then slat v s else 0 := by
  split_ifs with ha
  · have h := car_self hv a
    rw [← eq_sub_iff_add_eq'] at h
    rw [h, sub_mulVec, one_mulVec, ← mulVec_mulVec, crOp_mulVec_slat,
      slat_cons_eq_zero v s a ha, mulVec_zero, sub_zero]
  · rw [← mulVec_mulVec, adj_mulVec_slat hv s a ha, mulVec_zero]

/-- The hole operator `1 - n(v_a) = c(v_a) c†(v_a)` on a Slater determinant of the family. -/
lemma hole_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    ((crOp (v a))ᴴ * crOp (v a)) *ᵥ slat v s = if a ∈ univ.image s then 0 else slat v s := by
  rw [← mulVec_mulVec, crOp_mulVec_slat, adj_mulVec_slat_cons hv]

end Eig

/-! ### The operators `Ξ_j` -/

/-- `Ξ_j(Y) = ∑_{x_1, …, x_j} c_{x_1} ⋯ c_{x_j} Y c†_{x_j} ⋯ c†_{x_1}`, recursively. -/
def xiOp (Y : Matrix (Cfg L) (Cfg L) ℂ) : ℕ → Matrix (Cfg L) (Cfg L) ℂ
  | 0 => Y
  | j + 1 => ∑ x : Fin L, cplx (ann x) * xiOp Y j * cplx (cre x)

/-- The explicit sum of eq. `eq:main-rhoformula` is `Ξ_j`. -/
lemma sum_list_eq_xiOp (Y : Matrix (Cfg L) (Cfg L) ℂ) : ∀ j : ℕ,
    ∑ x : Fin j → Fin L, (List.ofFn fun i => cplx (ann (x i))).prod * Y *
      (List.ofFn fun i => cplx (cre (x (Fin.rev i)))).prod = xiOp Y j
  | 0 => by simp [xiOp]
  | j + 1 => by
    rw [← (Fin.consEquiv fun _ => Fin L).sum_comp, Fintype.sum_prod_type]
    simp only [Fin.consEquiv_apply, xiOp, ← sum_list_eq_xiOp Y j, Finset.mul_sum,
      Finset.sum_mul]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun x _ => ?_
    rw [List.ofFn_succ, List.prod_cons, List.ofFn_succ', List.prod_concat]
    simp only [Fin.cons_zero, Fin.cons_succ, Fin.rev_castSucc, Fin.rev_last]
    simp only [mul_assoc]

/-- **Completeness**: `∑_x c_x Z c†_x = ∑_i c(v_i) Z c†(v_i)` for an orthonormal basis `v`. -/
lemma sum_ann_cre_eq {v : Fin L → Fin L → ℂ}
    (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
    (Z : Matrix (Cfg L) (Cfg L) ℂ) :
    ∑ x : Fin L, cplx (ann x) * Z * cplx (cre x) = ∑ i, (crOp (v i))ᴴ * Z * crOp (v i) := by
  calc ∑ x : Fin L, cplx (ann x) * Z * cplx (cre x)
      = ∑ x : Fin L, ∑ y : Fin L, (if x = y then (1 : ℂ) else 0) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        refine Finset.sum_congr rfl fun x _ => ?_
        simp [ite_smul, Finset.sum_ite_eq]
    _ = ∑ x : Fin L, ∑ y : Fin L, ∑ i, (star (v i x) * v i y) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        simp only [← hcomp, Finset.sum_smul]
    _ = ∑ i, ∑ x : Fin L, ∑ y : Fin L, (star (v i x) * v i y) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        rw [Finset.sum_congr rfl fun x _ => Finset.sum_comm]
        exact Finset.sum_comm
    _ = ∑ i, (crOp (v i))ᴴ * Z * crOp (v i) := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [crOp_conjTranspose]
        simp only [crOp, Finset.sum_mul, Finset.mul_sum, smul_mul_assoc, mul_smul_comm,
          Finset.smul_sum, smul_smul]
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => ?_
        rw [mul_comm (star (v i x))]

/-! ### The inclusion–exclusion numbers -/

/-- `S_0(U) = 1`, `S_{j+1}(U) = ∑_{i ∉ U} ν_i S_j(U ∪ {i})`; `S_j(U) = j! e_j(ν|_{U^c})`. -/
def sFun (ν : Fin L → ℂ) : ℕ → Finset (Fin L) → ℂ
  | 0, _ => 1
  | j + 1, U => ∑ i ∈ univ \ U, ν i * sFun ν j (insert i U)

open Polynomial in
/-- `S_j(U)` is the `j`-th derivative at `0` of `∏_{i ∉ U} (1 + ν_i X)`. -/
lemma eval_iterate_derivative_prod (ν : Fin L → ℂ) (j : ℕ) : ∀ U : Finset (Fin L),
    (derivative^[j] (∏ i ∈ univ \ U, (C (ν i) * X + C 1))).eval 0 = sFun ν j U := by
  induction j with
  | zero => intro U; simp [sFun, eval_prod]
  | succ j ih =>
    intro U
    have hd : derivative (∏ i ∈ univ \ U, (C (ν i) * X + C 1)) =
        ∑ i ∈ univ \ U, C (ν i) * ∏ k ∈ univ \ insert i U, (C (ν k) * X + C 1) := by
      rw [derivative_prod_finset]
      refine Finset.sum_congr rfl fun i hi => ?_
      have he : (univ \ U).erase i = univ \ insert i U := by
        ext k
        simp only [mem_erase, mem_sdiff, mem_univ, true_and, mem_insert, not_or]
      rw [he, mul_comm]
      congr 1
      simp
    rw [Function.iterate_succ_apply, hd, iterate_derivative_sum, eval_finsetSum, sFun]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [iterate_derivative_C_mul, eval_mul, eval_C, ih]

open Polynomial in
/-- **Inclusion–exclusion**: `∑_{j ≤ L} (-1)^j / j! S_j(U) = ∏_{i ∉ U} (1 - ν_i)`. -/
lemma sum_sFun (ν : Fin L → ℂ) (U : Finset (Fin L)) :
    ∑ j ∈ range (L + 1), ((-1 : ℂ) ^ j / (j.factorial : ℂ)) * sFun ν j U =
      ∏ i ∈ univ \ U, (1 - ν i) := by
  set P : ℂ[X] := ∏ i ∈ univ \ U, (C (ν i) * X + C 1) with hP
  have hdeg : P.natDegree < L + 1 := by
    have h1 : P.natDegree ≤ ∑ i ∈ univ \ U, (C (ν i) * X + C 1).natDegree :=
      natDegree_prod_le _ _
    have h2 : ∑ i ∈ univ \ U, (C (ν i) * X + C 1).natDegree ≤ ∑ i ∈ univ \ U, 1 :=
      Finset.sum_le_sum fun i _ => natDegree_linear_le
    have h3 : ∑ i ∈ univ \ U, (1 : ℕ) ≤ L := by
      rw [Finset.sum_const, smul_eq_mul, mul_one]
      exact (card_le_univ _).trans (by simp)
    omega
  have hev : P.eval (-1) = ∏ i ∈ univ \ U, (1 - ν i) := by
    rw [hP, eval_prod]
    refine Finset.prod_congr rfl fun i _ => ?_
    simp only [eval_add, eval_mul, eval_C, eval_X]
    ring
  rw [← hev, eval_eq_sum_range' hdeg]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hc := coeff_iterate_derivative (k := j) P 0
  rw [zero_add, Nat.descFactorial_self, coeff_zero_eq_eval_zero,
    eval_iterate_derivative_prod, nsmul_eq_mul] at hc
  have hf : (j.factorial : ℂ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero j
  rw [hc]
  field_simp

/-- `∑_i [i ∉ U] f_i = ∑_{i ∈ U^c} f_i`. -/
lemma sum_ite_mem_zero {M : Type*} [AddCommMonoid M] (U : Finset (Fin L)) (f : Fin L → M) :
    ∑ i, (if i ∈ U then 0 else f i) = ∑ i ∈ univ \ U, f i := by
  rw [Finset.sdiff_eq_filter, Finset.sum_filter]
  exact Finset.sum_congr rfl fun i _ => by by_cases h : i ∈ U <;> simp [h]

/-! ### Action of `ρ[C]` on Slater determinants of eigenvectors -/

section Main

variable {v : Fin L → Fin L → ℂ}
  (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
  (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
  {C : Matrix (Fin L) (Fin L) ℂ} {ν : Fin L → ℂ} (hCv : ∀ i, C *ᵥ v i = ν i • v i)
include hv hcomp hCv

/-- `Ξ_j(Γ(C)) |s⟩ = ∏_k ν_{s_k} S_j(s) |s⟩` for a Slater determinant of eigenvectors. -/
lemma xiOp_mulVec_slat (j : ℕ) : ∀ {m : ℕ} (s : Fin m → Fin L),
    xiOp (secondQ C) j *ᵥ slat v s = ((∏ k, ν (s k)) * sFun ν j (univ.image s)) • slat v s := by
  induction j with
  | zero =>
    intro m s
    simp only [xiOp, sFun, mul_one]
    exact secondQ_mulVec_slat v C ν hCv s
  | succ j ih =>
    intro m s
    have h1 : ∀ i, ((crOp (v i))ᴴ * xiOp (secondQ C) j * crOp (v i)) *ᵥ slat v s =
        (if i ∈ univ.image s then 0 else
          ν i * ((∏ k, ν (s k)) * sFun ν j (insert i (univ.image s)))) • slat v s := by
      intro i
      rw [← mulVec_mulVec, ← mulVec_mulVec, crOp_mulVec_slat, ih, mulVec_smul,
        adj_mulVec_slat_cons hv, Fin.prod_univ_succ, image_cons]
      simp only [Fin.cons_zero, Fin.cons_succ]
      split_ifs <;> simp [mul_assoc]
    rw [xiOp, sum_ann_cre_eq hcomp, sum_mulVec]
    simp only [h1, ← Finset.sum_smul]
    congr 1
    rw [sum_ite_mem_zero, sFun, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring

/-- **The quasi-free state on Slater determinants of eigenvectors** (paper: proof of
Lemma `lem:quasifree`(a)): `ρ[C] |s⟩ = ∏_k ν_{s_k} ∏_{i ∉ s} (1 - ν_i) |s⟩`. -/
lemma qfState_mulVec_slat {m : ℕ} (s : Fin m → Fin L) :
    qfState C *ᵥ slat v s =
      ((∏ k, ν (s k)) * ∏ i ∈ univ \ univ.image s, (1 - ν i)) • slat v s := by
  rw [qfState, sum_mulVec]
  simp only [smul_mulVec, sum_list_eq_xiOp, xiOp_mulVec_slat hv hcomp hCv, smul_smul,
    ← Finset.sum_smul]
  rw [← sum_sFun ν (univ.image s), Finset.mul_sum]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by ring

end Main

/-! ### The Slater eigenbasis -/

/-- The Slater basis vector `|T⟩ = ∏_{i ∈ T} c†(v_i) |0⟩` (increasing order). -/
def bvec (v : Fin L → Fin L → ℂ) (T : Finset (Fin L)) : Cfg L → ℂ := slat v (T.orderEmbOfFin rfl)

/-- The increasing enumeration of `T` has image `T`. -/
lemma image_orderEmb (T : Finset (Fin L)) : univ.image (T.orderEmbOfFin rfl) = T :=
  image_orderEmbOfFin_univ T rfl

/-- `∏_k f(t_k) = ∏_{i ∈ T} f(i)` along the increasing enumeration. -/
lemma prod_orderEmb {M : Type*} [CommMonoid M] (T : Finset (Fin L)) (f : Fin L → M) :
    ∏ k, f ((T.orderEmbOfFin rfl) k) = ∏ i ∈ T, f i := by
  conv_rhs => rw [← image_orderEmb T]
  rw [Finset.prod_image fun x _ y _ h => (T.orderEmbOfFin rfl).injective h]

/-- The matrix of the Slater basis vectors (columns). -/
def bmat (v : Fin L → Fin L → ℂ) : Matrix (Cfg L) (Finset (Fin L)) ℂ :=
  Matrix.of fun n T => bvec v T n

/-- `#Cfg L = #Finset (Fin L) = 2^L`. -/
lemma card_cfg_eq : Fintype.card (Cfg L) = Fintype.card (Finset (Fin L)) := by
  simp [Fintype.card_finset, ZMod.card]

section Basis

variable {v : Fin L → Fin L → ℂ} (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
include hv

/-- **Orthonormality of the Slater basis**. -/
lemma bvec_orth (T T' : Finset (Fin L)) :
    star (bvec v T) ⬝ᵥ bvec v T' = if T = T' then 1 else 0 := by
  split_ifs with h
  · subst h
    rw [bvec, slat, gram_slaterAlt_diag]
    · exact Finset.prod_eq_one fun k _ => by rw [hv, ite_eq_left rfl]
    · intro k l hkl
      rw [hv, ite_eq_right_iff]
      intro h
      exact absurd ((T.orderEmbOfFin rfl).injective h) hkl
  · obtain ⟨i, hi⟩ : ∃ i, ¬ (i ∈ T ↔ i ∈ T') := by
      by_contra hc
      exact h (Finset.ext fun k => not_not.mp fun hk => hc ⟨k, hk⟩)
    by_cases hiT : i ∈ T
    · have hiT' : i ∉ T' := fun h' => hi ⟨fun _ => h', fun _ => hiT⟩
      have h1 := num_mulVec_slat hv (T.orderEmbOfFin rfl) i
      rw [image_orderEmb, ite_eq_left hiT] at h1
      rw [bvec, ← h1, ← mulVec_mulVec, dot_crOp_left, bvec,
        adj_mulVec_slat hv (T'.orderEmbOfFin rfl) i (by rwa [image_orderEmb]), dotProduct_zero]
    · have hiT' : i ∈ T' := by tauto
      have h1 := num_mulVec_slat hv (T'.orderEmbOfFin rfl) i
      rw [image_orderEmb, ite_eq_left hiT'] at h1
      rw [bvec, bvec, ← h1, ← mulVec_mulVec, dot_crOp_right,
        adj_mulVec_slat hv (T.orderEmbOfFin rfl) i (by rwa [image_orderEmb]), star_zero,
        zero_dotProduct]

/-- The Slater basis matrix is an isometry. -/
lemma bmat_ct_mul : (bmat v)ᴴ * bmat v = 1 := by
  refine Matrix.ext fun T T' => ?_
  rw [Matrix.mul_apply, Matrix.one_apply, ← bvec_orth hv T T']
  simp [bmat, dotProduct, conjTranspose_apply]

/-- The Slater basis matrix is unitary. -/
lemma bmat_mul_ct : bmat v * (bmat v)ᴴ = 1 :=
  (Matrix.mul_eq_one_comm_of_equiv (Fintype.equivOfCardEq card_cfg_eq.symm)).mp (bmat_ct_mul hv)

/-- An operator acting diagonally on the Slater basis is `W diag(d) W†`. -/
lemma eq_of_mulVec_bvec (X : Matrix (Cfg L) (Cfg L) ℂ) (d : Finset (Fin L) → ℂ)
    (h : ∀ T, X *ᵥ bvec v T = d T • bvec v T) : X = bmat v * diagonal d * (bmat v)ᴴ := by
  have hXW : X * bmat v = bmat v * diagonal d := by
    refine Matrix.ext fun n T => ?_
    rw [mul_diagonal]
    have := congrFun (h T) n
    simp only [mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at this
    rw [Matrix.mul_apply]
    simp only [bmat, of_apply]
    rw [this, mul_comm]
  calc X = X * (bmat v * (bmat v)ᴴ) := by rw [bmat_mul_ct hv, mul_one]
    _ = bmat v * diagonal d * (bmat v)ᴴ := by rw [← Matrix.mul_assoc, hXW]

end Basis

/-- The weights `w_T = ∏_{i ∈ T} ν_i ∏_{i ∉ T} (1 - ν_i)`. -/
def wT (ν : Fin L → ℂ) (T : Finset (Fin L)) : ℂ := ∏ i, if i ∈ T then ν i else 1 - ν i

/-- `w_T` as a product over `T` and its complement. -/
lemma wT_eq (ν : Fin L → ℂ) (T : Finset (Fin L)) :
    wT ν T = (∏ i ∈ T, ν i) * ∏ i ∈ univ \ T, (1 - ν i) := by
  rw [wT, ← Finset.prod_mul_prod_compl T, Finset.compl_eq_univ_sdiff]
  congr 1
  · exact Finset.prod_congr rfl fun k hk => ite_eq_left hk
  · exact Finset.prod_congr rfl fun k hk => ite_eq_right (Finset.mem_sdiff.mp hk).2

/-- **Diagonalization of `ρ[C]`** in the Slater eigenbasis of any orthonormal eigenbasis of `C`:
`ρ[C] = W diag(w) W†`. -/
theorem qfState_eq_conj {v : Fin L → Fin L → ℂ}
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
    (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
    {C : Matrix (Fin L) (Fin L) ℂ} {ν : Fin L → ℂ} (hCv : ∀ i, C *ᵥ v i = ν i • v i) :
    qfState C = bmat v * diagonal (wT ν) * (bmat v)ᴴ :=
  eq_of_mulVec_bvec hv _ _ fun T => by
    rw [bvec, qfState_mulVec_slat hv hcomp hCv, image_orderEmb, wT_eq, prod_orderEmb]

/-! ### Mathlib's eigenbasis -/

section Hermitian

variable {C : Matrix (Fin L) (Fin L) ℂ} (hC : C.IsHermitian)
include hC

/-- Orthonormality of the eigenvectors. -/
lemma eig_orth (i j : Fin L) :
    star (⇑(hC.eigenvectorBasis i) : Fin L → ℂ) ⬝ᵥ ⇑(hC.eigenvectorBasis j) =
      if i = j then 1 else 0 := by
  have h := orthonormal_iff_ite.mp hC.eigenvectorBasis.orthonormal i j
  rw [EuclideanSpace.inner_eq_star_dotProduct] at h
  rw [dotProduct_comm]
  exact h

/-- Completeness of the eigenvectors. -/
lemma eig_comp (x y : Fin L) :
    ∑ i, star ((⇑(hC.eigenvectorBasis i) : Fin L → ℂ) x) * ⇑(hC.eigenvectorBasis i) y =
      if x = y then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hC.eigenvectorUnitary) y) x
  rw [Matrix.mul_apply, Matrix.one_apply] at h
  simp only [IsHermitian.eigenvectorUnitary_apply] at h
  calc ∑ i, star ((⇑(hC.eigenvectorBasis i) : Fin L → ℂ) x) * ⇑(hC.eigenvectorBasis i) y
      = ∑ i, (hC.eigenvectorBasis i).ofLp y * star ((hC.eigenvectorBasis i).ofLp x) :=
        Finset.sum_congr rfl fun i _ => mul_comm _ _
    _ = if y = x then 1 else 0 := h
    _ = if x = y then 1 else 0 := by
        by_cases hxy : x = y
        · subst hxy; simp
        · simp [hxy, Ne.symm hxy]

/-- The eigenvalue equation with complex scalars. -/
lemma eig_eq (i : Fin L) :
    C *ᵥ ⇑(hC.eigenvectorBasis i) = ((hC.eigenvalues i : ℝ) : ℂ) • ⇑(hC.eigenvectorBasis i) := by
  rw [hC.mulVec_eigenvectorBasis]
  funext x
  simp [Complex.real_smul]

/-- `ρ[C]` in Mathlib's eigenbasis. -/
lemma qfState_eq_conj_eig :
    qfState C = bmat (fun i => ⇑(hC.eigenvectorBasis i)) *
      diagonal (wT (fun i => (hC.eigenvalues i : ℂ))) * (bmat (fun i => ⇑(hC.eigenvectorBasis i)))ᴴ :=
  qfState_eq_conj (eig_orth hC) (eig_comp hC) (eig_eq hC)

/-- `det(1 + (z - 1) C) = ∏_i (ν_i z + 1 - ν_i)`. -/
lemma det_one_add_smul_eig (z : ℂ) :
    (1 + (z - 1) • C).det = ∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i)) := by
  set U : Matrix (Fin L) (Fin L) ℂ := (hC.eigenvectorUnitary : Matrix (Fin L) (Fin L) ℂ)
  have hU1 : U * star U = 1 := Unitary.coe_mul_star_self _
  have hU2 : star U * U = 1 := Unitary.coe_star_mul_self _
  have hspec : C = U * diagonal (RCLike.ofReal ∘ hC.eigenvalues) * star U := by
    have := hC.spectral_theorem
    rwa [Unitary.conjStarAlgAut_apply] at this
  have hD : (1 : Matrix (Fin L) (Fin L) ℂ) + (z - 1) • diagonal (RCLike.ofReal ∘ hC.eigenvalues) =
      diagonal (fun i => (hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i)) := by
    ext i j
    by_cases h : i = j
    · subst h; simp; ring
    · simp [h]
  have h1 : 1 + (z - 1) • C =
      U * (1 + (z - 1) • diagonal (RCLike.ofReal ∘ hC.eigenvalues)) * star U := by
    rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hU1, Matrix.mul_smul, Matrix.smul_mul,
      ← hspec]
  have h2 : (star U).det * U.det = 1 := by rw [← det_mul, hU2, det_one]
  rw [h1, det_mul, det_mul, hD, det_diagonal]
  calc U.det * (∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i))) * (star U).det
      = (∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i))) * ((star U).det * U.det) := by
        ring
    _ = _ := by rw [h2, mul_one]

end Hermitian

/-- The eigenvalues of `C` are at most `1` if `1 - C ≥ 0`. -/
lemma eig_le_one {C : Matrix (Fin L) (Fin L) ℂ} (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef)
    (i : Fin L) : hC0.isHermitian.eigenvalues i ≤ 1 := by
  have h := hC1.dotProduct_mulVec_nonneg (⇑(hC0.isHermitian.eigenvectorBasis i))
  rw [sub_mulVec, one_mulVec, eig_eq hC0.isHermitian, dotProduct_sub, dotProduct_smul,
    eig_orth hC0.isHermitian, ite_eq_left rfl, smul_eq_mul, mul_one] at h
  have h' : (0 : ℂ) ≤ ((1 - hC0.isHermitian.eigenvalues i : ℝ) : ℂ) := by
    push_cast; exact h
  exact sub_nonneg.mp (Complex.zero_le_real.mp h')

/-- Real weights. -/
lemma wT_real (ν : Fin L → ℝ) (T : Finset (Fin L)) :
    wT (fun i => (ν i : ℂ)) T = ((∏ i, if i ∈ T then ν i else 1 - ν i : ℝ) : ℂ) := by
  rw [wT, Complex.ofReal_prod]
  exact Finset.prod_congr rfl fun i _ => by split_ifs <;> simp

/-- The weights sum to one. -/
lemma sum_wT (ν : Fin L → ℂ) : ∑ T, wT ν T = 1 := by
  have h := Finset.prod_add ν (fun i => 1 - ν i) univ
  rw [powerset_univ] at h
  rw [Finset.sum_congr rfl fun T _ => wT_eq ν T, ← h]
  simp

/-- Generating function of the weights. -/
lemma sum_wT_pow (ν : Fin L → ℂ) (z : ℂ) :
    ∑ T, wT ν T * z ^ T.card = ∏ i, (ν i * z + (1 - ν i)) := by
  have h := Finset.prod_add (fun i => ν i * z) (fun i => 1 - ν i) univ
  rw [powerset_univ] at h
  rw [h]
  refine Finset.sum_congr rfl fun T _ => ?_
  rw [wT_eq, Finset.prod_mul_distrib, Finset.prod_const]
  ring

/-- The product of commuting-type factors acting diagonally. -/
lemma list_prod_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] {n : ℕ} (F : Fin n → Matrix ι ι ℂ)
    (c : Fin n → ℂ) (φ : ι → ℂ) (h : ∀ i, F i *ᵥ φ = c i • φ) :
    (List.ofFn F).prod *ᵥ φ = (∏ i, c i) • φ := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.ofFn_succ, List.prod_cons, ← mulVec_mulVec,
      ih (fun i => F i.succ) (fun i => c i.succ) (fun i => h i.succ), mulVec_smul, h 0, smul_smul,
      Fin.prod_univ_succ, mul_comm]

/-- The factor `ν_i n_i + (1 - ν_i)(1 - n_i)` on the Slater basis. -/
lemma factor_mulVec_bvec {v : Fin L → Fin L → ℂ}
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) (ν : Fin L → ℝ) (i : Fin L)
    (T : Finset (Fin L)) :
    (((ν i : ℝ) : ℂ) • (crOp (v i) * (crOp (v i))ᴴ) +
        ((1 - ν i : ℝ) : ℂ) • ((crOp (v i))ᴴ * crOp (v i))) *ᵥ bvec v T =
      (if i ∈ T then (ν i : ℂ) else 1 - ν i) • bvec v T := by
  rw [add_mulVec, smul_mulVec, smul_mulVec, bvec, num_mulVec_slat hv, hole_mulVec_slat hv,
    image_orderEmb]
  split_ifs <;> simp

/-- `N = #occupied` as a real number. -/
lemma count_univ_eq (n : Cfg L) : count univ n = ((occSet n).card : ℝ) := by
  rw [count, occSet, card_filter, Nat.cast_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rcases FreeDynAux.zmod2_cases (n x) with h | h <;> simp [occ, h]

/-- `z^N |T⟩ = z^{|T|} |T⟩`. -/
lemma zPowN_mulVec_bvec (v : Fin L → Fin L → ℂ) (z : ℂ) (T : Finset (Fin L)) :
    zPowN z *ᵥ bvec v T = (z ^ T.card) • bvec v T := by
  funext n
  rw [zPowN, mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  by_cases hn : (occSet n).card = T.card
  · rw [hn]
  · have hnum := isNumR_slaterAlt univ (fun j => v (T.orderEmbOfFin rfl j)) (fun _ => True)
      (fun _ _ x hx => absurd (mem_univ x) hx) (fun _ h => absurd trivial h)
    have h0 : bvec v T n = 0 := hnum n (by rw [count_univ_eq]; simpa using hn)
    rw [h0, mul_zero, mul_zero]

/-- `c†(·)` is continuous. -/
lemma continuous_crOp : Continuous (crOp : (Fin L → ℂ) → Matrix (Cfg L) (Cfg L) ℂ) := by
  unfold crOp
  exact continuous_finsetSum _ fun x _ => (continuous_apply x).smul continuous_const

/-- Slater determinants depend continuously on the orbitals. -/
lemma continuous_slaterAlt_comp {X : Type*} [TopologicalSpace X] {m : ℕ}
    (w : Fin m → X → Fin L → ℂ) (hw : ∀ j, Continuous (w j)) :
    Continuous fun a => slaterAlt L m (fun j => w j a) := by
  simp only [slaterAlt_apply]
  refine Continuous.matrix_mulVec ?_ continuous_const
  simp only [List.ofFn_eq_map]
  exact continuous_list_prod _ fun j _ => continuous_crOp.comp (hw j)

/-- The second quantization is continuous. -/
lemma continuous_secondQ :
    Continuous (secondQ : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ) := by
  refine continuous_matrix fun n' n => ?_
  have h := continuous_slaterAlt_comp (m := (occSet n).card)
    (fun j (C : Matrix (Fin L) (Fin L) ℂ) x => C x ((occSet n).orderEmbOfFin rfl j))
    (fun j => continuous_pi fun x => Continuous.matrix_elem continuous_id x _)
  exact (continuous_apply n').comp h

/-- Entries of `U diag(d) U†` with `U` unitary and `|d_k| ≤ 1` are bounded by the dimension. -/
lemma norm_conj_diag_le {n : ℕ} (U : Matrix (Fin n) (Fin n) ℂ)
    (hU : U ∈ Matrix.unitaryGroup (Fin n) ℂ) (d : Fin n → ℂ) (hd : ∀ k, ‖d k‖ ≤ 1)
    (i j : Fin n) : ‖(U * diagonal d * star U) i j‖ ≤ n := by
  rw [Matrix.mul_apply]
  refine (norm_sum_le _ _).trans ?_
  have hb : ∀ k, ‖(U * diagonal d) i k * star U k j‖ ≤ 1 := by
    intro k
    rw [Matrix.mul_diagonal, Matrix.star_apply, norm_mul, norm_mul, norm_star]
    have h1 := entry_norm_bound_of_unitary hU i k
    have h2 := entry_norm_bound_of_unitary hU j k
    have h3 := hd k
    have h4 : 0 ≤ ‖U i k‖ := norm_nonneg _
    have h5 : 0 ≤ ‖d k‖ := norm_nonneg _
    calc ‖U i k‖ * ‖d k‖ * ‖U j k‖ ≤ 1 * 1 * 1 := by gcongr
      _ = 1 := by ring
  calc ∑ k, ‖(U * diagonal d) i k * star U k j‖ ≤ ∑ k : Fin n, (1 : ℝ) :=
        Finset.sum_le_sum fun k _ => hb k
    _ = n := by simp

/-- Entries of `0 ≤ C ≤ 1` are bounded by `L`. -/
lemma norm_entry_le {C : Matrix (Fin L) (Fin L) ℂ} (hC0 : C.PosSemidef)
    (hC1 : (1 - C).PosSemidef) (i j : Fin L) : ‖C i j‖ ≤ L := by
  have hspec := hC0.isHermitian.spectral_theorem
  rw [Unitary.conjStarAlgAut_apply] at hspec
  rw [hspec]
  refine norm_conj_diag_le _ hC0.isHermitian.eigenvectorUnitary.2 _ (fun k => ?_) i j
  simp only [Function.comp_apply]
  rw [RCLike.norm_ofReal, abs_le]
  constructor <;> linarith [hC0.eigenvalues_nonneg k, eig_le_one hC0 hC1 k]

end QuasiFreeAux

open QuasiFreeAux

/-- **Eigenmode form** (paper: Lemma `lem:quasifree`(a)): `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))`
for Hermitian `C` (the factors commute); in particular `ρ[C]` does not depend on the choice of eigenbasis. -/
theorem qfState_eq_prod (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) :
    qfState C = (List.ofFn fun i : Fin L =>
      ((hC.eigenvalues i : ℝ) : ℂ) • (crOp (hC.eigenvectorBasis i) * (crOp (hC.eigenvectorBasis i))ᴴ) +
        ((1 - hC.eigenvalues i : ℝ) : ℂ) •
          ((crOp (hC.eigenvectorBasis i))ᴴ * crOp (hC.eigenvectorBasis i))).prod := by
  rw [qfState_eq_conj_eig hC]
  symm
  refine eq_of_mulVec_bvec (eig_orth hC) _ _ fun T => ?_
  rw [list_prod_mulVec _ (fun i => if i ∈ T then (hC.eigenvalues i : ℂ) else 1 - hC.eigenvalues i) _
    (fun i => factor_mulVec_bvec (eig_orth hC) hC.eigenvalues i T)]
  rfl

/-- The matrix elements of `ρ[C]` are polynomials in the entries of `C`: `C ↦ ρ[C]` is continuous. -/
theorem continuous_qfState : Continuous (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ) := by
  unfold qfState
  refine continuous_finsetSum _ fun j _ =>
    Continuous.fun_const_smul (continuous_finsetSum _ fun x _ => ?_) _
  exact (continuous_const.matrix_mul continuous_secondQ).matrix_mul continuous_const

/-- **Uniform continuity** (paper: Lemma `lem:quasifree`(a)) on `{0 ≤ C ≤ 1}`. -/
theorem uniformContinuousOn_qfState :
    UniformContinuousOn (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ)
      {C | C.PosSemidef ∧ (1 - C).PosSemidef} := by
  refine IsCompact.uniformContinuousOn_of_continuous ?_ continuous_qfState.continuousOn
  have hclosed : IsClosed {C : Matrix (Fin L) (Fin L) ℂ | C.PosSemidef ∧ (1 - C).PosSemidef} :=
    isClosed_setOfPred_posSemidef.inter
      (isClosed_setOfPred_posSemidef.preimage (continuous_const.sub continuous_id))
  have hK : IsCompact (Set.univ.pi fun _ : Fin L => Set.univ.pi fun _ : Fin L =>
      Metric.closedBall (0 : ℂ) L) :=
    isCompact_univ_pi fun _ => isCompact_univ_pi fun _ => isCompact_closedBall _ _
  refine hK.of_isClosed_subset hclosed ?_
  rintro C ⟨hC0, hC1⟩
  simp only [Set.mem_pi, Set.mem_univ, true_implies, Metric.mem_closedBall, dist_zero_right]
  intro i j
  exact norm_entry_le hC0 hC1 i j

/-- `ρ[C] ≥ 0` for `0 ≤ C ≤ 1` (paper: Lemma `lem:quasifree`(b)). -/
theorem qfState_posSemidef (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    (qfState C).PosSemidef := by
  rw [qfState_eq_conj_eig hC0.isHermitian]
  refine PosSemidef.mul_mul_conjTranspose_same ?_ _
  rw [posSemidef_diagonal_iff]
  intro T
  rw [wT_real]
  refine Complex.zero_le_real.mpr (Finset.prod_nonneg fun i _ => ?_)
  split_ifs
  · exact hC0.eigenvalues_nonneg i
  · linarith [eig_le_one hC0 hC1 i]

/-- `Tr ρ[C] = 1` for Hermitian `C` (paper: Lemma `lem:quasifree`(b)). -/
theorem trace_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) : (qfState C).trace = 1 := by
  rw [qfState_eq_conj_eig hC, trace_mul_comm, ← Matrix.mul_assoc, bmat_ct_mul (eig_orth hC),
    Matrix.one_mul, trace_diagonal]
  exact sum_wT _

/-- **Entropy of a quasi-free state** (paper: Lemma `lem:quasifree`(b)): `S(ρ[C]) = ∑_i h(ν_i)`. -/
theorem vonNeumann_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    vonNeumann (qfState C) (qfState_posSemidef C hC0 hC1).isHermitian =
      ∑ i, Real.binEntropy (hC0.isHermitian.eigenvalues i) := by
  set hC := hC0.isHermitian
  set W := bmat (fun i => ⇑(hC.eigenvectorBasis i))
  set d : Finset (Fin L) → ℝ := fun T => ∏ i, if i ∈ T then hC.eigenvalues i else 1 - hC.eigenvalues i
  have hw : wT (fun i => (hC.eigenvalues i : ℂ)) = fun T => (d T : ℂ) :=
    funext fun T => wT_real _ T
  have heq : qfState C = W * diagonal (fun T => (d T : ℂ)) * Wᴴ := by
    rw [qfState_eq_conj_eig hC, hw]
  have hdH : (diagonal (fun T => (d T : ℂ))).IsHermitian := by
    refine isHermitian_diagonal_of_self_adjoint _ (funext fun T => ?_)
    simp [Complex.conj_ofReal]
  have hWH : (W * diagonal (fun T => (d T : ℂ)) * Wᴴ).IsHermitian :=
    heq ▸ (qfState_posSemidef C hC0 hC1).isHermitian
  rw [CutLocalAux.vonNeumann_congr _ hWH heq,
    CutLocalAux.vonNeumann_conj_isometry W (bmat_ct_mul (eig_orth hC)) _ hdH hWH,
    vonNeumann_diagonal, ← SlaterExactAux.sum_negMulLog_prod_sel hC.eigenvalues,
    ← (SlaterExactAux.boolEquiv L).sum_comp]
  refine Finset.sum_congr rfl fun σ _ => ?_
  congr 1
  refine Finset.prod_congr rfl fun k _ => ?_
  simp [SlaterExactAux.boolEquiv]

/-- **Generating function** (paper: Lemma `lem:quasifree`(b)): `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)`. -/
theorem trace_qfState_zN (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) (z : ℂ) :
    (qfState C * zPowN z).trace = (1 + (z - 1) • C).det := by
  have hv := eig_orth hC
  have hWW := bmat_ct_mul hv
  have hZ := eq_of_mulVec_bvec hv (zPowN (L := L) z) (fun T : Finset (Fin L) => z ^ T.card)
    fun T => zPowN_mulVec_bvec _ z T
  rw [qfState_eq_conj_eig hC, hZ, det_one_add_smul_eig hC z, ← sum_wT_pow]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc (bmat _)ᴴ (bmat _), hWW, Matrix.one_mul, trace_mul_comm,
    Matrix.mul_assoc, Matrix.mul_assoc, hWW, Matrix.mul_one, diagonal_mul_diagonal,
    trace_diagonal]

namespace QuasiFreeAux

open SlaterExactAux SlaterGeneralAux SlaterBulkAux

/-! ### Relabelling the modes of a region -/

section Transport

variable {N : ℕ} (R : Finset (Fin N))

/-- The configuration of `Λ` equal to `n` (relabelled by `R.orderIsoOfFin`) on `R` and empty outside. -/
def embR (n : Cfg R.card) : Cfg N :=
  (splitR R).symm (fun a : R => n ((R.orderIsoOfFin rfl).symm a), 0)

/-- `embR n` takes the value `n i` at the `i`-th site of `R`. -/
lemma embR_apply_mem (n : Cfg R.card) (i : Fin R.card) :
    embR R n (((R.orderIsoOfFin rfl) i : R) : Fin N) = n i := by
  rw [embR, splitR_symm_mem R _ _ _ ((R.orderIsoOfFin rfl) i).2]
  show n ((R.orderIsoOfFin rfl).symm ⟨_, _⟩) = n i
  rw [Subtype.coe_eta, OrderIso.symm_apply_apply]

/-- `embR n` is empty outside `R`. -/
lemma embR_apply_not_mem (n : Cfg R.card) (x : Fin N) (hx : x ∉ R) : embR R n x = 0 := by
  rw [embR, splitR_symm_not_mem R _ _ x hx]
  rfl

/-- Every site of `R` is the image of its relabelling. -/
lemma exists_orderIso_eq (x : Fin N) (hx : x ∈ R) :
    ∃ j, (((R.orderIsoOfFin rfl) j : R) : Fin N) = x :=
  ⟨(R.orderIsoOfFin rfl).symm ⟨x, hx⟩, by simp⟩

/-- Updating a site of `R` commutes with `embR`. -/
lemma embR_update (n : Cfg R.card) (i : Fin R.card) (b : ZMod 2) :
    Function.update (embR R n) (((R.orderIsoOfFin rfl) i : R) : Fin N) b =
      embR R (Function.update n i b) := by
  funext x
  by_cases hx : x ∈ R
  · obtain ⟨j, rfl⟩ := exists_orderIso_eq R x hx
    rw [embR_apply_mem]
    by_cases hij : j = i
    · subst hij; simp
    · have hne : (((R.orderIsoOfFin rfl) j : R) : Fin N) ≠ (((R.orderIsoOfFin rfl) i : R) : Fin N) :=
        fun h => hij ((R.orderIsoOfFin rfl).injective (Subtype.ext h))
      rw [Function.update_of_ne hne, embR_apply_mem, Function.update_of_ne hij]
  · have hne : x ≠ (((R.orderIsoOfFin rfl) i : R) : Fin N) :=
      fun h => hx (h ▸ ((R.orderIsoOfFin rfl) i).2)
    rw [Function.update_of_ne hne, embR_apply_not_mem R _ x hx, embR_apply_not_mem R _ x hx]

/-- The relabelling is increasing, so it preserves the Jordan–Wigner strings. -/
lemma jwPar_embR (n : Cfg R.card) (i : Fin R.card) :
    jwPar (embR R n) (((R.orderIsoOfFin rfl) i : R) : Fin N) = jwPar n i := by
  unfold jwPar
  rw [Finset.sum_filter, Finset.sum_filter,
    ← Finset.sum_subset (Finset.subset_univ R) (fun x _ hx => by
      rw [embR_apply_not_mem R n x hx, ite_self]),
    ← Finset.sum_coe_sort R, ← (R.orderIsoOfFin rfl).toEquiv.sum_comp]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hlt : ((((R.orderIsoOfFin rfl) j : R) : Fin N) < (((R.orderIsoOfFin rfl) i : R) : Fin N)) ↔
      j < i := Subtype.coe_lt_coe.trans (R.orderIsoOfFin rfl).lt_iff_lt
  simp only [OrderIso.coe_toEquiv]
  rw [embR_apply_mem]
  exact if_congr hlt rfl rfl

/-- `c_x†` for `x ∈ R` corresponds to `c_i†` on the relabelled modes. -/
lemma cre_embR (i : Fin R.card) (Y : Cfg N → ℂ) (n : Cfg R.card) :
    (cplx (cre (((R.orderIsoOfFin rfl) i : R) : Fin N)) *ᵥ Y) (embR R n) =
      (cplx (cre i) *ᵥ (fun n' => Y (embR R n'))) n := by
  rw [AreaLawAux.cre_mulVec_apply, AreaLawAux.cre_mulVec_apply, embR_apply_mem, jwPar_embR,
    embR_update]

/-- `c†(w)` for `w` supported in `R` corresponds to `c†(w ∘ ι)` on the relabelled modes. -/
lemma crOp_embR (w : Fin N → ℂ) (hw : ∀ x ∉ R, w x = 0) (Y : Cfg N → ℂ) (n : Cfg R.card) :
    (crOp w *ᵥ Y) (embR R n) =
      (crOp (fun i => w (((R.orderIsoOfFin rfl) i : R) : Fin N)) *ᵥ (fun n' => Y (embR R n'))) n := by
  rw [crOp_mulVec_eq, crOp_mulVec_eq]
  simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  rw [← Finset.sum_subset (Finset.subset_univ R) (fun x _ hx => by rw [hw x hx, zero_mul]),
    ← Finset.sum_coe_sort R, ← (R.orderIsoOfFin rfl).toEquiv.sum_comp]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [OrderIso.coe_toEquiv]
  rw [cre_embR]

/-- `embR n` is the vacuum iff `n` is. -/
lemma embR_eq_zero_iff (n : Cfg R.card) : embR R n = 0 ↔ n = 0 := by
  constructor
  · intro h
    funext i
    rw [← embR_apply_mem R n i, h]
    rfl
  · rintro rfl
    funext x
    by_cases hx : x ∈ R
    · obtain ⟨j, rfl⟩ := exists_orderIso_eq R x hx
      rw [embR_apply_mem]; rfl
    · rw [embR_apply_not_mem R _ x hx]; rfl

/-- Products of creation operators of orbitals supported in `R`, relabelled. -/
lemma prodSel_embR : ∀ {M : ℕ} (a : Fin M → Fin N → ℂ), (∀ k, ∀ x ∉ R, a k x = 0) →
    ∀ (σ : Fin M → Bool) (n : Cfg R.card),
    prodSel a σ (embR R n) = prodSel (fun k i => a k (((R.orderIsoOfFin rfl) i : R) : Fin N)) σ n
  | 0, _, _, σ, n => by
    rw [prodSel_zero, prodSel_zero]
    have h0 : fockCfg (∅ : Finset (Fin N)) = 0 := by funext x; simp [fockCfg]
    have h0' : fockCfg (∅ : Finset (Fin R.card)) = 0 := by funext x; simp [fockCfg]
    simp only [fockVec, Pi.single_apply, h0, h0', embR_eq_zero_iff]
  | M + 1, a, ha, σ, n => by
    rw [prodSel_succ, prodSel_succ]
    have ih : ∀ n', prodSel (fun k => a k.succ) (fun k => σ k.succ) (embR R n') =
        prodSel (fun k i => a k.succ (((R.orderIsoOfFin rfl) i : R) : Fin N))
          (fun k => σ k.succ) n' :=
      fun n' => prodSel_embR (fun k => a k.succ) (fun k => ha k.succ) (fun k => σ k.succ) n'
    split_ifs
    · rw [crOp_embR R _ (ha 0)]
      congr 2
      funext n'
      exact ih n'
    · rw [Matrix.one_mulVec, Matrix.one_mulVec]
      exact ih n

end Transport

/-! ### The reduced state of a Slater determinant as a mixture -/

/-- **The reduced density matrix of a Slater determinant** (from the proof of
`SlaterBulkAux.slater_core_bulk`): `ρ_R = ∑_σ μ_σ |X_σ⟩⟨X_σ|`, `X_σ = ∏_{σ k} c†(P_R w_k) |0⟩`. -/
lemma rdmF_slater_sum {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    ∃ (M : ℕ) (Φ : Matrix (Fin N) (Fin M) ℂ) (ν : Fin M → ℝ),
      corrMat A S t = Φ * Φᴴ ∧
      (∀ k l, ∑ x ∈ R, star (Φ x k) * Φ x l = if k = l then (ν k : ℂ) else 0) ∧
      ∀ α α', rdmF R (evo (dGamma A) t *ᵥ fockVec S) α α' =
        ∑ σ : Fin M → Bool, ((∏ k, if σ k then 1 else 1 - ν k : ℝ) : ℂ) *
          (prodSel (fun k x => if x ∈ R then Φ x k else 0) σ ((splitR R).symm (α, 0)) *
            star (prodSel (fun k x => if x ∈ R then Φ x k else 0) σ ((splitR R).symm (α', 0)))) := by
  obtain ⟨M, Φ, c, ν, hc, hψ, hC, h1, h2⟩ := slater_setup_gen A hA S t R
  refine ⟨M, Φ, ν, hC, h1, ?_⟩
  set a : Fin M → Fin N → ℂ := fun k x => if x ∈ R then Φ x k else 0 with ha_def
  set b : Fin M → Fin N → ℂ := fun k x => if x ∈ R then 0 else Φ x k with hb_def
  have dbb : ∀ k l, star (b k) ⬝ᵥ b l = if k = l then ((1 - ν k : ℝ) : ℂ) else 0 := by
    intro k l
    rw [← h2 k l]
    have e : ∀ x, star (b k) x * b l x =
        if x ∈ univ \ R then star (Φ x k) * Φ x l else 0 := by
      intro x; by_cases hx : x ∈ R <;> simp [hb_def, hx]
    simp only [dotProduct, e, Finset.sum_ite_mem, Finset.univ_inter]
  obtain ⟨s, hs, hexp⟩ := slater_expand_gen R (fun k x => Φ x k)
  set X : (Fin M → Bool) → Cfg N → ℂ := fun σ => prodSel a σ with hX_def
  set Y : (Fin M → Bool) → Cfg N → ℂ := fun σ => prodSel b (fun k => !σ k) with hY_def
  set ψ := evo (dGamma A) t *ᵥ fockVec S with hψ_def
  have hψ' : ψ = ∑ σ, (c * s σ) • ftens R (X σ) (Y σ) := by
    rw [hψ, hexp, Finset.smul_sum]
    simp only [smul_smul]
    rfl
  have hcs : ∀ σ, ‖c * s σ‖ = 1 := fun σ => by rw [norm_mul, hc, hs, one_mul]
  set μ : (Fin M → Bool) → ℝ := fun σ => ∏ k, (if σ k then 1 else 1 - ν k) with hμ_def
  have hYg : ∀ σ τ, star (Y τ) ⬝ᵥ Y σ = if σ = τ then (μ σ : ℂ) else 0 := by
    intro σ τ
    have hg := gram_prodSel b _ dbb (fun k => !τ k) (fun k => !σ k)
    by_cases h : σ = τ
    · subst h
      rw [hY_def]
      simp only [hg, ite_true, hμ_def]
      push_cast
      refine Finset.prod_congr rfl fun k _ => ?_
      cases σ k <;> simp
    · have h' : (fun k => !τ k) ≠ (fun k => !σ k) := by
        intro h'; apply h; funext k; have := congrFun h' k; simpa using this.symm
      rw [hY_def]
      simp only [hg, h, h', ite_false]
  have hYn : ∀ σ, IsNumR R 0 (Y σ) := fun σ =>
    isNumR_prodSel_zero R b (fun k x hx => by simp [hb_def, hx]) _
  have hv : ∀ α β, ψ ((splitR R).symm (α, β)) = (reorderSign R α β : ℂ) *
      ∑ σ, (c * s σ) * (X σ ((splitR R).symm (α, 0)) * Y σ ((splitR R).symm (0, β))) := by
    intro α β
    rw [hψ', Finset.sum_apply, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    rw [Pi.smul_apply, smul_eq_mul, ftens_splitR_symm]
    ring
  have hYh : ∀ σ τ, ∑ β, star (Y τ ((splitR R).symm (0, β))) * Y σ ((splitR R).symm (0, β)) =
      if σ = τ then (μ σ : ℂ) else 0 := by
    intro σ τ; rw [← dot_eq_sum_out R _ _ (hYn τ), hYg]
  exact fun α α' => rdmF_eq_sum R ψ (fun σ => c * s σ) (fun σ α => X σ ((splitR R).symm (α, 0)))
    (fun σ β => Y σ ((splitR R).symm (0, β))) (fun σ => (μ σ : ℂ)) hv hcs hYh α α'


/-! ### The abstract identity `∑_σ μ_σ |P_σ⟩⟨P_σ| = ρ[∑_k |u_k⟩⟨u_k|]` -/

/-- An orthonormal family of `L` vectors in `ℂ^L` is complete. -/
lemma comp_of_orth {L : ℕ} (v : Fin L → Fin L → ℂ)
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) (x y : Fin L) :
    ∑ i, star (v i x) * v i y = if x = y then 1 else 0 := by
  set V : Matrix (Fin L) (Fin L) ℂ := Matrix.of fun x i => v i x with hV
  have h1 : Vᴴ * V = 1 := by
    refine Matrix.ext fun i j => ?_
    rw [Matrix.mul_apply, Matrix.one_apply, ← hv i j]
    simp [hV, dotProduct, conjTranspose_apply]
  have h2 : V * Vᴴ = 1 := (Matrix.mul_eq_one_comm_of_equiv (Equiv.refl _)).mp h1
  have h := congrFun (congrFun h2 y) x
  rw [Matrix.mul_apply, Matrix.one_apply] at h
  simp only [hV, conjTranspose_apply, of_apply] at h
  calc ∑ i, star (v i x) * v i y = ∑ i, v i y * star (v i x) :=
        Finset.sum_congr rfl fun i _ => mul_comm _ _
    _ = if y = x then 1 else 0 := h
    _ = if x = y then 1 else 0 := by
        by_cases hxy : x = y
        · subst hxy; simp
        · simp [hxy, Ne.symm hxy]

/-- A mixture of rank-one operators acting on a vector. -/
lemma mulVec_sum_rankOne {ι κ : Type*} [Fintype ι] [Fintype κ] (c : κ → ℂ) (P : κ → ι → ℂ)
    (φ : ι → ℂ) :
    (Matrix.of fun n n' => ∑ σ, c σ * (P σ n * star (P σ n'))) *ᵥ φ =
      ∑ σ, (c σ * (star (P σ) ⬝ᵥ φ)) • P σ := by
  funext n
  simp only [mulVec, dotProduct, of_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul, Finset.mul_sum, Pi.star_apply]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun σ _ => Finset.sum_congr rfl fun n' _ => by ring

/-- `∑_k |u_k⟩⟨u_k|` acting on a vector. -/
lemma mulVec_sum_proj {M ℓ : ℕ} (u : Fin M → Fin ℓ → ℂ) (w : Fin ℓ → ℂ) :
    (Matrix.of fun i j => ∑ k, u k i * star (u k j)) *ᵥ w = ∑ k, (star (u k) ⬝ᵥ w) • u k := by
  funext i
  simp only [mulVec, dotProduct, of_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul,
    Finset.sum_mul, Pi.star_apply]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun j _ => by ring

/-- Products of selected creation operators are Slater determinants of the selected orbitals. -/
lemma prodSel_eq_slaterAlt {K : ℕ} : ∀ {M : ℕ} (w : Fin M → Fin K → ℂ) (σ : Fin M → Bool),
    ∃ (m : ℕ) (s : Fin m → Fin M), Function.Injective s ∧
      (∀ k, σ k = true ↔ k ∈ univ.image s) ∧ prodSel w σ = slaterAlt K m (fun j => w (s j))
  | 0, w, σ => ⟨0, Fin.elim0, fun i => i.elim0, fun k => k.elim0, by
      rw [prodSel_zero, slaterAlt_apply]; simp⟩
  | M + 1, w, σ => by
    obtain ⟨m, s, hs, hσ, hP⟩ := prodSel_eq_slaterAlt (fun k => w k.succ) (fun k => σ k.succ)
    rw [prodSel_succ, hP]
    by_cases h0 : σ 0 = true
    · refine ⟨m + 1, Fin.cons 0 (fun j => (s j).succ), ?_, ?_, ?_⟩
      · rw [Fin.cons_injective_iff]
        exact ⟨fun ⟨j, hj⟩ => Fin.succ_ne_zero _ hj, fun a b h => hs (Fin.succ_injective _ h)⟩
      · intro k
        refine Fin.cases ?_ (fun k => ?_) k
        · simp only [h0, true_iff]
          exact mem_image.mpr ⟨0, mem_univ _, rfl⟩
        · rw [hσ k]
          simp only [mem_image, mem_univ, true_and, Fin.exists_fin_succ, Fin.cons_zero,
            Fin.cons_succ]
          constructor
          · rintro ⟨j, hj⟩; exact Or.inr ⟨j, by rw [hj]⟩
          · rintro (h | ⟨j, hj⟩)
            · exact absurd h.symm (Fin.succ_ne_zero k)
            · exact ⟨j, Fin.succ_injective _ hj⟩
      · rw [ite_eq_left h0, slaterAlt_succ]
        simp only [Fin.cons_zero, Fin.cons_succ]
    · refine ⟨m, fun j => (s j).succ, fun a b h => hs (Fin.succ_injective _ h), ?_, ?_⟩
      · intro k
        refine Fin.cases ?_ (fun k => ?_) k
        · simp only [mem_image, mem_univ, true_and]
          constructor
          · intro h; exact absurd h h0
          · rintro ⟨j, hj⟩; exact absurd hj (Fin.succ_ne_zero _)
        · rw [hσ k]
          simp only [mem_image, mem_univ, true_and]
          exact ⟨fun ⟨j, hj⟩ => ⟨j, by rw [hj]⟩, fun ⟨j, hj⟩ => ⟨j, Fin.succ_injective _ hj⟩⟩
      · rw [ite_eq_right h0, one_mulVec]

/-- Changing the length argument of `orderEmbOfFin` in a Slater basis vector. -/
lemma slat_orderEmb_cast {L : ℕ} (v : Fin L → Fin L → ℂ) (T : Finset (Fin L)) {m : ℕ}
    (h : T.card = m) : slat v (T.orderEmbOfFin h) = bvec v T := by
  subst h; rfl

/-- A Slater determinant of distinct basis orbitals is `±` the Slater basis vector. -/
lemma slat_eq_smul_bvec {L : ℕ} (v : Fin L → Fin L → ℂ) {m : ℕ} (t : Fin m → Fin L)
    (ht : Function.Injective t) :
    ∃ ε : ℂ, star ε * ε = 1 ∧ slat v t = ε • bvec v (univ.image t) := by
  set τ := Tuple.sort t
  have hmono : StrictMono (t ∘ τ) :=
    (Tuple.monotone_sort t).strictMono_of_injective (ht.comp τ.injective)
  have hc : (univ.image t).card = m := by
    rw [card_image_of_injective _ ht, card_univ, Fintype.card_fin]
  have hs := Finset.orderEmbOfFin_unique hc (fun x => mem_image_of_mem t (mem_univ (τ x))) hmono
  have h1 : slat v t = Equiv.Perm.sign τ • slat v (fun x => t (τ x)) :=
    (slaterAlt L m).map_congr_perm _ τ
  rw [hs, slat_orderEmb_cast] at h1
  rcases Int.units_eq_one_or (Equiv.Perm.sign τ) with h | h
  · exact ⟨1, by simp, by rw [h1, h, one_smul, one_smul]⟩
  · exact ⟨-1, by simp, by rw [h1, h]; simp⟩


/-- **Adapted eigenbasis**: for orthogonal `u_k` with `‖u_k‖² = ν_k`, an orthonormal basis `v` of `ℂ^ℓ`
containing the normalized nonzero `u_k` (at the indices `π k`), with the eigenvalues `λ` of
`∑_k |u_k⟩⟨u_k|` (`λ_{π k} = ν_k`, `λ = 0` off the range of `π`). -/
lemma exists_adapted_basis {M ℓ : ℕ} (u : Fin M → Fin ℓ → ℂ) (ν : Fin M → ℝ)
    (hu : ∀ k l, star (u k) ⬝ᵥ u l = if k = l then (ν k : ℂ) else 0) :
    ∃ (v : Fin ℓ → Fin ℓ → ℂ) (lam : Fin ℓ → ℝ) (π : {k : Fin M // ν k ≠ 0} → Fin ℓ),
      (∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) ∧ Function.Injective π ∧
      (∀ k : {k : Fin M // ν k ≠ 0}, u k = ((√(ν k) : ℝ) : ℂ) • v (π k)) ∧
      (∀ k : {k : Fin M // ν k ≠ 0}, lam (π k) = ν k) ∧
      (∀ i, (∀ k, π k ≠ i) → lam i = 0) := by
  classical
  have hν0 : ∀ k, 0 ≤ ν k := fun k => by
    have h := hu k k
    rw [ite_eq_left rfl, AreaLawAux.star_dotProduct_self] at h
    have h' : AreaLawAux.nsq (u k) = ν k := by exact_mod_cast h
    rw [← h']
    exact AreaLawAux.nsq_nonneg _
  have hsq : ∀ k, ν k ≠ 0 → ((√(ν k) : ℝ) : ℂ) ≠ 0 := fun k hk => by
    have : 0 < √(ν k) := Real.sqrt_pos.mpr (lt_of_le_of_ne (hν0 k) (Ne.symm hk))
    exact_mod_cast this.ne'
  set g : Fin M → Fin ℓ → ℂ := fun k => (((√(ν k) : ℝ) : ℂ))⁻¹ • u k with hg
  have hug : ∀ k, ν k ≠ 0 → u k = ((√(ν k) : ℝ) : ℂ) • g k := fun k hk => by
    rw [hg, smul_smul, mul_inv_cancel₀ (hsq k hk), one_smul]
  have hgo : ∀ k l, ν k ≠ 0 → ν l ≠ 0 → star (g k) ⬝ᵥ g l = if k = l then 1 else 0 := by
    intro k l hk hl
    simp only [hg, star_smul, smul_dotProduct, dotProduct_smul, hu, smul_eq_mul]
    split_ifs with hkl
    · subst hkl
      have h1 : ((√(ν k) : ℝ) : ℂ) * ((√(ν k) : ℝ) : ℂ) = (ν k : ℂ) := by
        rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν0 k)]
      rw [Complex.star_def, map_inv₀, Complex.conj_ofReal, ← h1]
      field_simp [hsq k hk]
    · simp
  set m := Fintype.card {k : Fin M // ν k ≠ 0} with hm_def
  set fam0 : {k : Fin M // ν k ≠ 0} → EuclideanSpace ℂ (Fin ℓ) := fun k => WithLp.toLp 2 (g k)
    with hfam0_def
  have hfam0 : Orthonormal ℂ fam0 := by
    rw [orthonormal_iff_ite]
    intro k l
    rw [EuclideanSpace.inner_eq_star_dotProduct, dotProduct_comm]
    simp only [hfam0_def, WithLp.ofLp_toLp]
    rw [hgo k l k.2 l.2]
    exact if_congr Subtype.ext_iff.symm rfl rfl
  have hm : m ≤ ℓ := by
    have := hfam0.linearIndependent.fintype_card_le_finrank
    rwa [finrank_euclideanSpace_fin] at this
  set eJ := Fintype.equivFin {k : Fin M // ν k ≠ 0} with heJ
  set π : {k : Fin M // ν k ≠ 0} → Fin ℓ := fun k => Fin.castLE hm (eJ k) with hπ_def
  have hπ : Function.Injective π := fun a b h => eJ.injective (Fin.castLE_injective hm h)
  have hπk : ∀ k, (⟨(π k).val, (eJ k).isLt⟩ : Fin m) = eJ k := fun k => Fin.ext rfl
  set fam : Fin ℓ → EuclideanSpace ℂ (Fin ℓ) := fun i =>
    if h : i.val < m then fam0 (eJ.symm ⟨i.val, h⟩) else 0 with hfam_def
  have hfam : Orthonormal ℂ ({i : Fin ℓ | i.val < m}.domRestrict fam) := by
    rw [orthonormal_iff_ite]
    rintro ⟨i, hi⟩ ⟨j, hj⟩
    have hi' : i.val < m := hi
    have hj' : j.val < m := hj
    simp only [Set.domRestrict_apply, hfam_def, dite_eq_left hi', dite_eq_left hj']
    rw [orthonormal_iff_ite.mp hfam0]
    refine if_congr ⟨fun h => ?_, fun h => ?_⟩ rfl rfl
    · have h2 := eJ.symm.injective h
      simp only [Fin.mk.injEq] at h2
      exact Subtype.ext (Fin.ext h2)
    · cases h; rfl
  obtain ⟨b, hb⟩ := hfam.exists_orthonormalBasis_extension_of_card_eq (by simp)
  refine ⟨fun i => (b i).ofLp, fun i => if h : i.val < m then ν (eJ.symm ⟨i.val, h⟩).val else 0,
    π, ?_, hπ, ?_, ?_, ?_⟩
  · intro i j
    have h := orthonormal_iff_ite.mp b.orthonormal i j
    rw [EuclideanSpace.inner_eq_star_dotProduct] at h
    rw [dotProduct_comm]
    exact h
  · intro k
    have hk : (π k).val < m := (eJ k).isLt
    have h1 : b (π k) = fam (π k) := hb (π k) hk
    have h2 : fam (π k) = fam0 k := by
      simp only [hfam_def, dite_eq_left hk]
      rw [hπk, Equiv.symm_apply_apply]
    simp only [h1, h2, hfam0_def, WithLp.ofLp_toLp]
    exact hug k k.2
  · intro k
    have hk : (π k).val < m := (eJ k).isLt
    simp only [dite_eq_left hk]
    rw [hπk, Equiv.symm_apply_apply]
  · intro i hi
    dsimp only
    split_ifs with h
    · exfalso
      refine hi (eJ.symm ⟨i.val, h⟩) (Fin.ext ?_)
      simp [hπ_def]
    · rfl

/-- The adapted basis is an eigenbasis of `∑_k |u_k⟩⟨u_k|`. -/
lemma adapted_eig {M ℓ : ℕ} (u : Fin M → Fin ℓ → ℂ) (ν : Fin M → ℝ) (v : Fin ℓ → Fin ℓ → ℂ)
    (lam : Fin ℓ → ℝ) (π : {k : Fin M // ν k ≠ 0} → Fin ℓ)
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) (hπ : Function.Injective π)
    (hup : ∀ k : {k : Fin M // ν k ≠ 0}, u k = ((√(ν k) : ℝ) : ℂ) • v (π k))
    (hu0 : ∀ k, ν k = 0 → u k = 0) (hν0 : ∀ k, 0 ≤ ν k)
    (hlamπ : ∀ k : {k : Fin M // ν k ≠ 0}, lam (π k) = ν k)
    (hlam0 : ∀ i, (∀ k, π k ≠ i) → lam i = 0) (i : Fin ℓ) :
    (Matrix.of fun i j => ∑ k, u k i * star (u k j)) *ᵥ v i = (lam i : ℂ) • v i := by
  classical
  rw [mulVec_sum_proj]
  have hterm : ∀ k, (star (u k) ⬝ᵥ v i) • u k =
      (if h : ν k = 0 then 0 else if π ⟨k, h⟩ = i then (ν k : ℂ) else 0) • v i := by
    intro k
    by_cases hk : ν k = 0
    · rw [dite_eq_left hk, hu0 k hk, smul_zero, zero_smul]
    · rw [dite_eq_right hk, hup ⟨k, hk⟩, star_smul, smul_dotProduct, hv]
      split_ifs with h
      · subst h
        have h1 : ((√(ν k) : ℝ) : ℂ) * ((√(ν k) : ℝ) : ℂ) = (ν k : ℂ) := by
          rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν0 k)]
        rw [smul_smul, Complex.star_def, Complex.conj_ofReal, smul_eq_mul, mul_one, h1]
      · simp
  simp only [hterm, ← Finset.sum_smul]
  congr 1
  by_cases hi : ∃ k0, π k0 = i
  · obtain ⟨k0, rfl⟩ := hi
    rw [hlamπ, Finset.sum_eq_single k0.val]
    · rw [dite_eq_right k0.2, ite_eq_left rfl]
    · intro k _ hk
      split_ifs with h1 h2
      · rfl
      · exact absurd (congrArg Subtype.val (hπ h2)) hk
      · rfl
    · intro h; exact absurd (mem_univ _) h
  · push Not at hi
    rw [hlam0 i hi, Complex.ofReal_zero]
    refine Finset.sum_eq_zero fun k _ => ?_
    split_ifs with h1 h2
    · rfl
    · exact absurd h2 (hi _)
    · rfl


open Classical in
/-- The weight of `|T⟩⟨T|` contributed by the term `σ` of the mixture `∑_σ μ_σ |P_σ⟩⟨P_σ|`. -/
noncomputable def coefQ {M ℓ : ℕ} (ν : Fin M → ℝ) (π : {k : Fin M // ν k ≠ 0} → Fin ℓ)
    (T : Finset (Fin ℓ)) (σ : Fin M → Bool) : ℂ :=
  if (∀ k, σ k = true → ν k ≠ 0) ∧
      (univ.filter (fun k : {k : Fin M // ν k ≠ 0} => σ k.val = true)).image π = T then
    ((∏ k, if σ k then ν k else 1 - ν k : ℝ) : ℂ)
  else 0

/-- Membership in the image set of `coefQ`. -/
lemma mem_image_filter {M ℓ : ℕ} (ν : Fin M → ℝ) (π : {k : Fin M // ν k ≠ 0} → Fin ℓ)
    (σ : Fin M → Bool) (i : Fin ℓ) :
    i ∈ (univ.filter (fun k : {k : Fin M // ν k ≠ 0} => σ k.val = true)).image π ↔
      ∃ k : {k : Fin M // ν k ≠ 0}, σ k.val = true ∧ π k = i := by
  simp only [mem_image, mem_filter, mem_univ, true_and]

/-- The term `σ` of the mixture, applied to the basis vector `|T⟩` of the adapted basis. -/
lemma term_eq_coefQ {M ℓ : ℕ} (u : Fin M → Fin ℓ → ℂ) (ν : Fin M → ℝ) (v : Fin ℓ → Fin ℓ → ℂ)
    (π : {k : Fin M // ν k ≠ 0} → Fin ℓ)
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) (hπ : Function.Injective π)
    (hup : ∀ k : {k : Fin M // ν k ≠ 0}, u k = ((√(ν k) : ℝ) : ℂ) • v (π k))
    (hu0 : ∀ k, ν k = 0 → u k = 0) (hν0 : ∀ k, 0 ≤ ν k) (σ : Fin M → Bool)
    (T : Finset (Fin ℓ)) :
    (((∏ k, if σ k then 1 else 1 - ν k : ℝ) : ℂ) * (star (prodSel u σ) ⬝ᵥ bvec v T)) •
        prodSel u σ = coefQ ν π T σ • bvec v T := by
  classical
  obtain ⟨m', s, hs, hσs, hP⟩ := prodSel_eq_slaterAlt u σ
  unfold coefQ
  by_cases hsel : ∀ k, σ k = true → ν k ≠ 0
  · have hsJ : ∀ j, ν (s j) ≠ 0 := fun j =>
      hsel (s j) ((hσs (s j)).mpr (mem_image_of_mem s (mem_univ j)))
    set tt : Fin m' → Fin ℓ := fun j => π ⟨s j, hsJ j⟩ with htt_def
    have htt : Function.Injective tt := fun a b h => hs (congrArg Subtype.val (hπ h))
    have hP2 : prodSel u σ = ((∏ j, √(ν (s j)) : ℝ) : ℂ) • slat v tt := by
      rw [hP, slat, Complex.ofReal_prod]
      have e : (fun j => u (s j)) = fun j => ((√(ν (s j)) : ℝ) : ℂ) • v (tt j) :=
        funext fun j => hup ⟨s j, hsJ j⟩
      rw [e, AlternatingMap.map_smul_univ]
    obtain ⟨ε, hε, hslat⟩ := slat_eq_smul_bvec v tt htt
    have hmem : ∀ i, i ∈ univ.image tt ↔ ∃ k : {k : Fin M // ν k ≠ 0}, σ k.val = true ∧ π k = i := by
      intro i
      simp only [mem_image, mem_univ, true_and]
      constructor
      · rintro ⟨j, rfl⟩
        exact ⟨⟨s j, hsJ j⟩, (hσs (s j)).mpr (mem_image_of_mem s (mem_univ j)), rfl⟩
      · rintro ⟨k, hk, rfl⟩
        obtain ⟨j, -, hj⟩ := mem_image.mp ((hσs k.val).mp hk)
        exact ⟨j, by simp only [htt_def]; congr 1; exact Subtype.ext hj⟩
    have hiffT : (univ.filter (fun k : {k : Fin M // ν k ≠ 0} => σ k.val = true)).image π = T ↔
        univ.image tt = T := by
      simp only [Finset.ext_iff, mem_image_filter, hmem]
    have hprodR : (∏ k, if σ k then 1 else 1 - ν k) * ∏ j, ν (s j) =
        ∏ k, if σ k then ν k else 1 - ν k := by
      have h1 : ∏ j, ν (s j) = ∏ k, if σ k then ν k else 1 := by
        rw [← Finset.prod_image (fun a _ b _ h => hs h)]
        have h2 : univ.image s = univ.filter (fun k => σ k = true) := by
          ext k; simp only [mem_filter, mem_univ, true_and]; exact (hσs k).symm
        rw [h2, Finset.prod_filter]
      rw [h1, ← Finset.prod_mul_distrib]
      exact Finset.prod_congr rfl fun k _ => by cases σ k <;> simp
    have hcc : ((∏ j, √(ν (s j)) : ℝ) : ℂ) * ((∏ j, √(ν (s j)) : ℝ) : ℂ) = ((∏ j, ν (s j) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, ← Finset.prod_mul_distrib]
      congr 1
      exact Finset.prod_congr rfl fun j _ => Real.mul_self_sqrt (hν0 _)
    have hsc : star ((∏ j, √(ν (s j)) : ℝ) : ℂ) = ((∏ j, √(ν (s j)) : ℝ) : ℂ) :=
      Complex.conj_ofReal _
    rw [hP2, hslat, star_smul, star_smul, smul_dotProduct, smul_dotProduct, bvec_orth hv, smul_smul,
      smul_smul, hsc]
    by_cases hT : univ.image tt = T
    · rw [ite_eq_left hT, ite_eq_left ⟨hsel, hiffT.mpr hT⟩, hT]
      congr 1
      simp only [smul_eq_mul, mul_one]
      calc _ = ((∏ k, if σ k then 1 else 1 - ν k : ℝ) : ℂ) *
            (((∏ j, √(ν (s j)) : ℝ) : ℂ) * ((∏ j, √(ν (s j)) : ℝ) : ℂ)) * (star ε * ε) := by
            ring
        _ = ((∏ k, if σ k then ν k else 1 - ν k : ℝ) : ℂ) := by
            rw [hcc, hε, mul_one, ← Complex.ofReal_mul, hprodR]
    · rw [ite_eq_right hT, ite_eq_right (fun h => hT (hiffT.mp h.2))]
      simp
  · push Not at hsel
    obtain ⟨k, hk1, hk2⟩ := hsel
    obtain ⟨j, -, hj⟩ := mem_image.mp ((hσs k).mp hk1)
    have hP0 : prodSel u σ = 0 := by
      rw [hP]
      exact (slaterAlt ℓ m').map_coord_zero j (by simp only [hj, hu0 k hk2])
    rw [hP0, smul_zero, ite_eq_right (fun h => h.1 k hk1 hk2), zero_smul]

/-- **The weights match**: `∑_σ (weight of |T⟩⟨T| in term σ) = w_T` for the adapted eigenvalues. -/
lemma sum_coefQ {M ℓ : ℕ} (ν : Fin M → ℝ) (lam : Fin ℓ → ℝ) (π : {k : Fin M // ν k ≠ 0} → Fin ℓ)
    (hπ : Function.Injective π) (hlamπ : ∀ k, lam (π k) = ν k)
    (hlam0 : ∀ i, (∀ k, π k ≠ i) → lam i = 0) (T : Finset (Fin ℓ)) :
    ∑ σ, coefQ ν π T σ = wT (fun i => (lam i : ℂ)) T := by
  classical
  by_cases hT : ∀ i ∈ T, ∃ k, π k = i
  · set σT : Fin M → Bool := fun k => decide (∃ h : ν k ≠ 0, π ⟨k, h⟩ ∈ T) with hσT
    have hσTk : ∀ k : {k : Fin M // ν k ≠ 0}, σT k.val = true ↔ π k ∈ T := by
      intro k
      simp only [hσT, decide_eq_true_iff]
      exact ⟨fun ⟨_, h⟩ => h, fun h => ⟨k.2, h⟩⟩
    have hiff : ∀ σ, ((∀ k, σ k = true → ν k ≠ 0) ∧
        (univ.filter (fun k : {k : Fin M // ν k ≠ 0} => σ k.val = true)).image π = T) ↔
          σ = σT := by
      intro σ
      constructor
      · rintro ⟨hsel, himg⟩
        funext k
        by_cases hk : ν k = 0
        · have h1 : σ k = false := by
            cases h : σ k
            · rfl
            · exact absurd hk (hsel k h)
          have h2 : σT k = false := by
            simp only [hσT, decide_eq_false_iff_not]
            rintro ⟨h, -⟩
            exact h hk
          rw [h1, h2]
        · have h3 : σ k = true ↔ π ⟨k, hk⟩ ∈ T := by
            rw [← himg, mem_image_filter]
            constructor
            · intro h; exact ⟨⟨k, hk⟩, h, rfl⟩
            · rintro ⟨k', hk', hπk⟩
              have := congrArg Subtype.val (hπ hπk)
              simp only at this
              rw [← this]
              exact hk'
          exact Bool.eq_iff_iff.mpr (h3.trans (hσTk ⟨k, hk⟩).symm)
      · rintro rfl
        refine ⟨fun k hk => ?_, ?_⟩
        · simp only [hσT, decide_eq_true_iff] at hk
          exact hk.1
        · ext i
          rw [mem_image_filter]
          constructor
          · rintro ⟨k, hk, rfl⟩
            exact (hσTk k).mp hk
          · intro hi
            obtain ⟨k, rfl⟩ := hT i hi
            exact ⟨k, (hσTk k).mpr hi, rfl⟩
    have hc : ∀ σ, coefQ ν π T σ =
        if σ = σT then ((∏ k, if σT k then ν k else 1 - ν k : ℝ) : ℂ) else 0 := by
      intro σ
      unfold coefQ
      by_cases h : σ = σT
      · subst h
        rw [ite_eq_left ((hiff _).mpr rfl), ite_eq_left rfl]
      · rw [ite_eq_right (fun h' => h ((hiff σ).mp h')), ite_eq_right h]
    rw [Finset.sum_congr rfl fun σ _ => hc σ, Finset.sum_ite_eq' univ σT, ite_eq_left (mem_univ _),
      wT_real]
    congr 1
    have hL : ∏ k, (if σT k then ν k else 1 - ν k) =
        ∏ k : {k : Fin M // ν k ≠ 0}, (if π k ∈ T then ν k else 1 - ν k) := by
      rw [← Finset.prod_filter_mul_prod_filter_not univ (fun k => ν k ≠ 0)]
      have h1 : ∏ k ∈ univ.filter (fun k => ¬ ν k ≠ 0), (if σT k then ν k else 1 - ν k) = 1 := by
        refine Finset.prod_eq_one fun k hk => ?_
        have hk' : ν k = 0 := by simpa using (mem_filter.mp hk).2
        have h2 : σT k = false := by
          simp only [hσT, decide_eq_false_iff_not]
          rintro ⟨h, -⟩
          exact h hk'
        rw [h2]
        simp [hk']
      rw [h1, mul_one, Finset.prod_subtype (p := fun k : Fin M => ν k ≠ 0)
        (univ.filter (fun k => ν k ≠ 0)) (fun k => by simp)]
      refine Finset.prod_congr rfl fun k _ => ?_
      by_cases hk : π k ∈ T
      · rw [ite_eq_left ((hσTk k).mpr hk), ite_eq_left hk]
      · rw [ite_eq_right (fun h => hk ((hσTk k).mp h)), ite_eq_right hk]
    have hR : ∏ i, (if i ∈ T then lam i else 1 - lam i) =
        ∏ k : {k : Fin M // ν k ≠ 0}, (if π k ∈ T then ν k else 1 - ν k) := by
      rw [← Finset.prod_mul_prod_compl (univ.image π)]
      have h1 : ∏ i ∈ (univ.image π)ᶜ, (if i ∈ T then lam i else 1 - lam i) = 1 := by
        refine Finset.prod_eq_one fun i hi => ?_
        have hi' : ∀ k, π k ≠ i := fun k hk =>
          (mem_compl.mp hi) (mem_image.mpr ⟨k, mem_univ _, hk⟩)
        have hiT : i ∉ T := fun h => by
          obtain ⟨k, hk⟩ := hT i h
          exact hi' k hk
        rw [ite_eq_right hiT, hlam0 i hi', sub_zero]
      rw [h1, mul_one, Finset.prod_image (fun a _ b _ h => hπ h)]
      exact Finset.prod_congr rfl fun k _ => by rw [hlamπ]
    rw [hL, hR]
  · push Not at hT
    obtain ⟨i0, hi0T, hi0⟩ := hT
    have h0 : ∀ σ, coefQ ν π T σ = 0 := by
      intro σ
      unfold coefQ
      rw [ite_eq_right]
      rintro ⟨-, h⟩
      rw [← h, mem_image_filter] at hi0T
      obtain ⟨k, -, hk⟩ := hi0T
      exact hi0 k hk
    rw [Finset.sum_congr rfl fun σ _ => h0 σ, Finset.sum_const_zero, wT]
    symm
    refine Finset.prod_eq_zero (mem_univ i0) ?_
    rw [ite_eq_left hi0T, hlam0 i0 hi0, Complex.ofReal_zero]

/-- **Reduced states of Slater determinants are quasi-free** (abstract form, paper: proof of
Lemma `lem:slater`(a)): for orthogonal `u_k` with `‖u_k‖² = ν_k`,
`∑_σ μ_σ |P_σ⟩⟨P_σ| = ρ[∑_k |u_k⟩⟨u_k|]`, `P_σ = ∏_{σ k} c†(u_k) |0⟩`. -/
theorem qfState_eq_sum_prodSel {M ℓ : ℕ} (u : Fin M → Fin ℓ → ℂ) (ν : Fin M → ℝ)
    (hu : ∀ k l, star (u k) ⬝ᵥ u l = if k = l then (ν k : ℂ) else 0) :
    qfState (Matrix.of fun i j => ∑ k, u k i * star (u k j)) =
      Matrix.of fun n n' => ∑ σ : Fin M → Bool, ((∏ k, if σ k then 1 else 1 - ν k : ℝ) : ℂ) *
        (prodSel u σ n * star (prodSel u σ n')) := by
  classical
  have hν0 : ∀ k, 0 ≤ ν k := fun k => by
    have h := hu k k
    rw [ite_eq_left rfl, AreaLawAux.star_dotProduct_self] at h
    have h' : AreaLawAux.nsq (u k) = ν k := by exact_mod_cast h
    rw [← h']
    exact AreaLawAux.nsq_nonneg _
  have hu0 : ∀ k, ν k = 0 → u k = 0 := fun k hk => by
    have h := hu k k
    rw [ite_eq_left rfl, AreaLawAux.star_dotProduct_self, hk] at h
    have h' : AreaLawAux.nsq (u k) = 0 := by exact_mod_cast h
    funext x
    have h2 := (Finset.sum_eq_zero_iff_of_nonneg (fun x _ => by positivity)).mp h' x (mem_univ x)
    simpa using h2
  obtain ⟨v, lam, π, hv, hπ, hup, hlamπ, hlam0⟩ := exists_adapted_basis u ν hu
  have hCv := adapted_eig u ν v lam π hv hπ hup hu0 hν0 hlamπ hlam0
  rw [qfState_eq_conj hv (comp_of_orth v hv) hCv]
  symm
  refine eq_of_mulVec_bvec hv _ _ fun T => ?_
  rw [mulVec_sum_rankOne]
  simp only [term_eq_coefQ u ν v π hv hπ hup hu0 hν0, ← Finset.sum_smul]
  rw [sum_coefQ ν lam π hπ hlamπ hlam0 T]

end QuasiFreeAux

/-- **Reduced states of Slater determinants are quasi-free** (paper: Lemma `lem:slater`(a), `ρ_A = ρ[C_A]`): for any
region `R`, relabelled increasingly by `Fin |R|`. -/
theorem rdmF_eq_qfState {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) :
    (rdmF R (evo (dGamma A) t *ᵥ fockVec S)).submatrix
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a))
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a)) =
      qfState ((corrR A S t R).submatrix (fun i => (R.orderIsoOfFin rfl) i) (fun i => (R.orderIsoOfFin rfl) i)) := by
  obtain ⟨M, Φ, ν, hC, h1, hρ⟩ := rdmF_slater_sum A hA S t R
  set u : Fin M → Fin R.card → ℂ := fun k i => Φ (((R.orderIsoOfFin rfl) i : R) : Fin N) k
    with hu_def
  have hu : ∀ k l, star (u k) ⬝ᵥ u l = if k = l then (ν k : ℂ) else 0 := by
    intro k l
    rw [← h1 k l, ← Finset.sum_coe_sort R, ← (R.orderIsoOfFin rfl).toEquiv.sum_comp]
    simp only [dotProduct, Pi.star_apply, hu_def, OrderIso.coe_toEquiv]
  have hCu : (corrR A S t R).submatrix (fun i => (R.orderIsoOfFin rfl) i)
      (fun i => (R.orderIsoOfFin rfl) i) = Matrix.of fun i j => ∑ k, u k i * star (u k j) := by
    refine Matrix.ext fun i j => ?_
    simp [corrR, hC, Matrix.mul_apply, conjTranspose_apply, hu_def]
  have hu' : (fun k i => (fun k x => if x ∈ R then Φ x k else 0) k
      (((R.orderIsoOfFin rfl) i : R) : Fin N)) = u := by
    funext k i
    simp [hu_def]
  rw [hCu, qfState_eq_sum_prodSel u ν hu]
  refine Matrix.ext fun n n' => ?_
  rw [submatrix_apply, hρ, of_apply]
  refine Finset.sum_congr rfl fun σ _ => ?_
  have e1 := prodSel_embR R (fun k x => if x ∈ R then Φ x k else 0)
    (fun k x hx => by simp [hx]) σ n
  have e2 := prodSel_embR R (fun k x => if x ∈ R then Φ x k else 0)
    (fun k x hx => by simp [hx]) σ n'
  rw [hu'] at e1 e2
  simp only [embR] at e1 e2
  rw [e1, e2]

end DFL
