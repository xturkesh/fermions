# Number entropy in U(1)-symmetric disordered chains: rigorous statements, proofs, and checks

Date: 2026-11-22.

Scripts: `/tmp/claude-0/-home-user-fermions/e9e4eeb3-4dfa-5769-bb28-5161ded2b166/scratchpad/code/`.

**Status tags.**
- **[P]**: complete proof given in this note.
- **[P-c]**: complete proof, with routine estimates condensed.
- **[S]**: proof sketch; the gap is named explicitly.
- **[H]**: heuristic.
- **[N]**: numerics run for this note (scripts listed in §9).

Literature claims keep the tags of the colleagues' reports: [snippet] means only abstract-level search excerpts were seen; [memory] means not checked.

**Access limitations.** No paper could be opened in this session: arxiv.org does not resolve, and the shared web-search budget was exhausted. In particular, the gate parametrisation of Aceituno Chávez et al. (arXiv:2312.13420, their Eqs. (8)–(10)) was **not** checked. Everything below that concerns their model uses only the abstract-level facts listed in the literature report. Those facts are: U is a finite-depth circuit of U(1)-preserving two-site gates, H′ is diagonal, and the initial states are Néel states.

---

## 0. Results at a glance

| # | Statement (informal) | Status | Novelty (§8) |
|---|---|---|---|
| T1 (Thm 3.1, Cor 3.2) | Take H = U D U†, with D diagonal in the occupation basis and U U(1)-invariant. Then ‖N_A(t) − U N_A U†‖ = δ_A := ‖U†N_A U − N_A‖ for all t. Consequently Var N_A(t) ≤ (σ₀ + 2δ_A)² and S_N(t) ≤ ½ ln[2πe((σ₀+2δ_A)² + 1/12)], uniformly in t and L. Explicit δ_A: finite circuits δ_A ≤ d_u·\|∂A\| (and then *all* Rényi S_N, Hartley included, are ≤ ln(4d_u\|∂A\|+1)); quasi-local generators δ_A ≤ J ξ² \|∂A\|; random products of local rotations, in probability or expectation. | [P] | Ingredients standard; the combined statement is folklore, apparently never written as a theorem |
| T1′ (Thm 3.3) | **U need not be U(1)-invariant.** It suffices that the LIOMs τ_x = U n_x U† commute with N, which is automatic if [H,N] = 0 and D has simple spectrum. Bound: ‖N_A(t) − N_A‖ ≤ 24\|∂A\| Σ_j (j+1) η(j), where η are the LIOM tails. | [P] | Apparently new in this form |
| Tails (Rem 3.5) | First-moment tails Σ(j+1)η(j) < ∞ suffice: exponential, stretched-exponential, power law j^{-β} with β > 2. Not necessary. In a Gaussian example the variance stays bounded for hopping tails r^{-1.5} (heuristically for anything faster than r^{-1}) and diverges for r^{-1} or slower. | [P] + [N] | New, elementary |
| T2 (Prop 4.1) | For circuit l-bit models (Chávez-type) the *dephased* (infinite-time) number distribution p̄ is **independent of D** (given simple spectrum) and depends only on the gates within 4d_u − 2 sites of the cut. Its law is therefore *exactly L-independent* for L ≥ 8d_u − 4. Moreover the time average Ŝ_N(L) ≤ H(p̄) holds for every realisation. Numerically, Ŝ_N(L) rises with L, E H(p̄) is L-independent, and the gap E[H(p̄) − Ŝ_N] halves with every two added sites. | [P] + [N] | Apparently new |
| T2 verdict | Chávez et al.'s model is in the T1 class. Their L-growing saturation values must stay below an L-independent bound: deterministically ≤ 4.88 for d_u = 16, and in the dephased sense ≤ E H(p̄). There is no contradiction. "Saturation increasing with L" is the generic finite-size signature of a *localised* l-bit model. | [P] + [N] + [H] | |
| T3 (Thm 5.1, Prop 5.2) | Quasi-free U(1) states: ½ ln(8 Var/π) ≤ S_N ≤ ½ ln(2πe(Var+1/12)), Var ≤ S/(4 ln 2), and 2Var ≤ S^{(2)} ≤ 4 ln2·Var. Hence S_N ≤ ½ ln(2πe(S/(4ln2)+1/12)). Also S_N = ½ ln(2πe Var) + o(1) as Var → ∞. | [P] / [P-c] | Upper bounds folklore; lower bound and Rényi-2 sandwich elementary, possibly new in this form |
| T3 (Prop 5.5) | General interacting U(1) states: S_N ≤ S is optimal. For every s ≤ ln(min(\|A\|,\|A^c\|,N)+1) there are states with S_N = S = s. So no sublinear f exists. | [P] | Elementary |
| T4 (Thm 6.2) | Random t-V/XXZ chain, any L. With probability ≥ 1 − \|∂A\| 2^{8n−2}γ/W: ‖N_A(t) − N_A‖ ≤ \|∂A\|(4J/γ + n!(10J/γ)^n J\|t\|) for all t. Hence S_N stays O(1) up to t\* ≥ J⁻¹ exp[(ln(pW/10J\|∂A\|))²/(32 ln 2)·(1 − o(1))] with probability 1 − p. | [P] | Elementary single-cut analogue of known quasi-conservation results; constants useless at physical W |
| T5 (§7) | Logic: unbounded S_N ⇒ no U(1)-invariant LIOMs with E ln(1+δ_A) bounded. The converse is false. Finite-window ln ln t growth is compatible with T1. Sharp diagnostics: max-over-eigenstates Var_E(N_A) and the dephased variance Var_{p̄}(L). | [P] | |

---

## 1. Setting and notation

**Chain and occupation basis.**
- Λ = {1,…,L} is an open chain. On a ring, read distances as ring distances and let ∂A count cut bonds.
- The local space is C², with occupation n_x, N = Σ_x n_x, and N_X = Σ_{x∈X} n_x.
- Statements are made for hard-core bosons (spin-½). For spinless fermions in 1D, every operator used here is U(1)-invariant, hence even. The Jordan–Wigner map sends an even operator supported in X to one supported in the interval hull of X, and all supports used below are intervals, so nothing changes. In Theorem 3.3, U is assumed parity-preserving.
- Occupation (Fock) basis: |c⟩, c ∈ {0,1}^Λ. "Diagonal" means diagonal in this basis.

**The maps E_diag and Π_Y.**
- E_diag(X) := Σ_c |c⟩⟨c|X|c⟩⟨c| = ∫ e^{iθ·n} X e^{−iθ·n} dθ/(2π)^L.
  - It is a product of single-site maps.
  - It is a unital CP contraction, preserves supports, and commutes with partial traces.
- Π_Y is the normalised partial trace onto Y: Π_Y(X) = 2^{−|Y^c|} tr_{Y^c}(X) ⊗ 1_{Y^c}.
  - It is a unital CP contraction.
  - Π_Y X = X if X is supported in Y.
  - Π_Y X = tr(X)·1 if X is supported in Y^c, where tr is the normalised trace.

**Region, boundary, transfer quantities.**
- A ⊂ Λ. ∂A is the set of bonds {x, x+1} with exactly one endpoint in A. For a half chain |∂A| = 1; for a bulk interval |∂A| = 2.
- d_x := dist(x, Λ∖A) for x ∈ A. Then #{x∈A : d_x = d} ≤ |∂A| for every d ≥ 1.
- m_X(A) := min(|X∩A|, |X∖A|). This is the maximal charge an operation on X can move across ∂A.
- B_j(x) := {y : |y − x| ≤ j}.

**Dynamics and number statistics.**
- Heisenberg picture: N_A(t) = e^{iHt} N_A e^{−iHt}.
- For an initial state ρ, p_t(n) := Tr(e^{−iHt}ρe^{iHt} P_n^A), where P_n^A is the spectral projection of N_A. Equivalently, p_t is the law of the integer-valued observable N_A(t) in ρ.
- S_N(t) := −Σ_n p_t(n) ln p_t(n), in nats. S_N^{(α)} is the Rényi version. σ₀² := Var_ρ(N_A).
- For Fock states σ₀ = 0. For Néel states n₀ := N_A is sharp.

---

## 2. Tools

**Lemma 2.1 (integer max-entropy bound) [P].** Let X be integer-valued with Var X = σ² < ∞. Then H(X) ≤ ½ ln[2πe(σ² + 1/12)].

*Proof.*
1. Let V ~ Unif[0,1) be independent of X, and set Y = X + V.
2. Y has density f(y) = p_{⌊y⌋}. Hence h(Y) = −∫ f ln f = H(X), and Var Y = σ² + 1/12.
3. The Gaussian maximises differential entropy at fixed variance: h(Y) ≤ ½ ln(2πe Var Y). ∎

*Reference:* Cover & Thomas, *Elements of Information Theory*, 2nd ed. (2006), §17.3, the theorem bounding H(p₁,p₂,…) by ½ log 2πe(Σp_i i² − (Σ i p_i)² + 1/12). The numbering, Thm 17.3.3 (Thm 16.3.3 in the 1st ed.), is from [memory]. It is also often attributed to Massey (1988) [memory, not verified]. The proof above makes the citation non-load-bearing.

*Remarks.*
- Rényi α ≥ 1 entropies are ≤ Shannon, so they are bounded too.
- Rényi α < 1, including Hartley, is **not** controlled by σ². Example: p = (1−ε)δ₀ + (ε/M)Σ_{k=1}^M δ_k has Var ≤ εM² but S^{(0)} = ln(M+1).
- [N] check_inequalities.py: over 2·10⁴ random integer laws, max(H − bound) = −0.045 < 0.

**Lemma 2.2 (variance under a bounded perturbation) [P].** For self-adjoint X, Y and any state ρ:
- std_ρ(X+Y) ≤ std_ρ(X) + std_ρ(Y);
- std_ρ(Y) ≤ ½(λ_max(Y) − λ_min(Y)) ≤ ‖Y‖.

*Proof.*
1. std_ρ(X) = ‖(X − ⟨X⟩_ρ)ρ^{1/2}‖_HS. The map X ↦ X − ⟨X⟩ is linear, so std is a seminorm.
2. Var_ρ Y = min_c Tr ρ(Y−c)² ≤ min_c ‖Y − c‖². ∎

**Lemma 2.3 (charge moved by U(1)-invariant local operations) [P].** Let X ⊂ Λ.

- (i) Let k = k† be supported in X with [k, N_X] = 0. Then ‖[k, N_A]‖ ≤ m_X(A)‖k‖. In particular [k, N_A] = 0 if X does not straddle ∂A.
- (ii) Let O be unitary, supported in X, with [O, N_X] = 0. Then O†N_A O − N_A is supported in X, and ‖O†N_A O − N_A‖ ≤ m_X(A) · min{1, ‖O − 1‖}.
- (iii) Circuits. Let U = O_K ⋯ O_1 with each O_j as in (ii) on X_j. Then
  δ_A(U) := ‖U†N_A U − N_A‖ ≤ Σ_j m_{X_j}(A) min{1, ‖O_j − 1‖}.
- (iv) Generators. Let U = 𝒯exp(−i∫₀¹K(s)ds) (time-ordered) with K(s) = Σ_X k_X(s), each k_X(s) as in (i). Then
  δ_A(U) ≤ ∫₀¹ Σ_X m_X(A)‖k_X(s)‖ ds.
- (v) Support. In (iii), suppose ψ ∈ Ran 1_I(N_A). Then Uψ ∈ Ran 1_{I+[−κ,κ]}(N_A), where κ := Σ_{j: m_{X_j}(A)>0} m_{X_j}(A).
- (vi) ‖U N_A U† − N_A‖ = δ_A(U).

*Proof.*
- (i)
  - [k, N_A] = [k, N_{X∩A}], and also [k, N_{X∩A}] = −[k, N_{X∖A}].
  - Let M be either N_{X∩A} or N_{X∖A}, with spectrum in [0, m]. Then [k, M] = [k, M − m/2], so ‖[k,M]‖ ≤ m‖k‖. Take the smaller m.
- (ii)
  - O†N_A O − N_A = O†[N_{X∩A}, O] = −(O†N_{X∖A}O − N_{X∖A}).
  - For M as in (i): O†[M, O] = O†[M − m/2, O − 1], which gives the bound m‖O−1‖.
  - Also O†MO and M both have spectrum in [0, m], so their difference has norm ≤ m.
- (iii) Telescoping with W_j := O_{j−1}⋯O₁:
  U†N_AU − N_A = Σ_j W_j†(O_j†N_A O_j − N_A)W_j.
- (iv) ∂_s(U(s)†N_AU(s)) = U(s)† i[K(s),N_A] U(s). Integrate and apply (i).
- (v)
  - A gate with m_{X_j}(A) = 0 commutes with N_A.
  - Otherwise O_j preserves the joint eigenspaces of (N_{A∖X}, N_X). On such an eigenspace, with eigenvalue s for N_X, N_A − N_{A∖X} = N_{X∩A} ranges over [max(0, s−|X∖A|), min(s, |X∩A|)]. This interval has length ≤ m_X(A).
- (vi) U N_A U† − N_A = U(N_A − U†N_AU)U†. ∎

**Lemma 2.4 (dephasing) [P].** If D is diagonal with simple spectrum, then for every X,
lim_{τ→∞} τ⁻¹∫₀^τ e^{iDt} X e^{−iDt} dt = E_diag(X).

**Lemma 2.5 (boundary lemma) [P].** Let {g_x}_{x∈Λ} be operators with tr g_x = 0 and Σ_{x∈Λ} g_x = 0. Put
ε(j) := max_x ‖g_x − Π_{B_j(x)} g_x‖ for j ≥ 0, and ε(−1) := max_x ‖g_x‖.
Then for every A:
‖Σ_{x∈A} g_x‖ ≤ 6 |∂A| Σ_{j≥0} (j+1) ε(j).
If in addition every g_x is supported in B_R(x), then ‖Σ_{x∈A} g_x‖ ≤ 2R|∂A|·max_x‖g_x‖.

*Proof.*
1. Decompose g_x = Σ_{k=0}^{L} g_x^{(k)}, where
   g_x^{(k)} := Π_{B_k(x)}g_x − Π_{B_{k−1}(x)}g_x and Π_{B_{−1}} := tr.
   - Each g_x^{(k)} is traceless and supported in B_k(x).
   - ‖g_x^{(k)}‖ ≤ ε(k−1) + ε(k).
2. For x ∈ A, split the sum at k = d_x:
   - I_A := Σ_{x∈A} Σ_{k<d_x} g_x^{(k)} is supported in A.
   - C_A := Σ_{x∈A} Σ_{k≥d_x} g_x^{(k)}.
   - Define I_{A^c} and C_{A^c} in the same way.
3. Bound C_A by counting. Using #{x∈A : d_x = d} ≤ |∂A|:
   ‖C_A‖ ≤ |∂A| Σ_{k≥1} k(ε(k−1) + ε(k)) ≤ 2|∂A| Σ_{j≥0}(j+1)ε(j).
   The same bound holds for C_{A^c}.
4. Bound I_A. From Σ_x g_x = 0 we get I_A + C_A + I_{A^c} + C_{A^c} = 0. Also Π_A I_A = I_A and Π_A I_{A^c} = tr(I_{A^c})·1 = 0. Hence I_A = −Π_A(C_A + C_{A^c}), and
   ‖I_A‖ ≤ ‖C_A‖ + ‖C_{A^c}‖.
5. Therefore ‖Σ_{x∈A} g_x‖ = ‖I_A + C_A‖ ≤ 2‖C_A‖ + ‖C_{A^c}‖, which gives the first bound.
6. Strict support.
   - Σ_{x∈A} g_x = −Σ_{x∉A} g_x is supported in S := {dist(·,A) ≤ R} ∩ {dist(·,A^c) ≤ R}.
   - Hence Σ_{x∈A} g_x = Σ_{x∈A} Π_S g_x.
   - Π_S g_x = 0 unless B_R(x) meets S, which forces d_x ≤ 2R. There are at most 2R|∂A| such x. ∎

---

## 3. T1: U(1) LIOMs imply a time- and volume-uniform bound on number fluctuations

**Theorem 3.1 [P].** Let D be diagonal, let U be unitary with [U, N] = 0, and let H := U D U†. Then [H, N] = 0. Write δ_A := ‖U†N_AU − N_A‖. For all t ∈ ℝ and every state ρ:

- (a) ‖N_A(t) − U N_A U†‖ = δ_A, and ‖N_A(t) − N_A‖ ≤ 2δ_A.
- (b) Var_ρ N_A(t) ≤ (σ₀ + 2δ_A)².
- (c) S_N^{(α)}(t) ≤ S_N(t) ≤ ½ ln[2πe((σ₀+2δ_A)² + 1/12)] for all α ≥ 1. For Fock states this reads S_N(t) ≤ ½ ln[2πe(4δ_A² + 1/12)].
- (d) Every eigenvector U|c⟩ of H has Var(N_A) ≤ δ_A². For a Fock state, the infinite-time average of ⟨(N_A(t) − n₀)²⟩ is ≤ 4δ_A².

*Proof.*
1. Let Δ := U†N_AU − N_A. Since [D, N_A] = 0,
   N_A(t) = U e^{iDt}(N_A + Δ)e^{−iDt} U† = U N_A U† + U Δ(t) U†, where Δ(t) = e^{iDt}Δe^{−iDt}.
   This gives the first claim of (a).
2. N_A(t) − N_A = (U N_A U† − N_A) + UΔ(t)U†. With Lemma 2.3(vi) this gives the second claim of (a).
3. (b): The operator N_A(t) − N_A has spectral range ≤ 2δ_A + 2δ_A. Apply Lemma 2.2.
4. (c): Apply Lemma 2.1 to p_t.
5. (d): Since (N_A − n_c)|c⟩ = 0, Var_{U|c⟩}(N_A) = Var_{|c⟩}(N_A + Δ) = Var_{|c⟩}(Δ) ≤ ‖Δ‖². Also ⟨(N_A(t) − n₀)²⟩ = ‖(N_A(t) − N_A)ψ‖² ≤ 4δ_A². ∎

**Corollary 3.2 (explicit constants by class) [P].**

**(a) Finite-depth U(1) circuits.**
- For any product of U(1)-invariant gates, δ_A ≤ Σ over gates straddling ∂A of m_X(A)·min{1, ‖O − 1‖}.
- For a brickwork of two-site gates acting d_u times on each bond: δ_A ≤ d_u|∂A|. Hence S_N(t) ≤ ½ ln[2πe(4d_u²|∂A|² + 1/12)] for Fock states, for all t, L and realisations.
- Moreover, for a Fock state p_t is supported in {|n − n₀| ≤ 2d_u|∂A|}. To see this, apply Lemma 2.3(v) to U†, note that e^{−iDt} conserves N_A, then apply 2.3(v) to U. Hence **every** Rényi entropy, including Hartley, satisfies S_N^{(α)}(t) ≤ ln(4d_u|∂A| + 1).
- Numbers for a half chain with d_u = 16 (Chávez et al.'s largest runs [snippet]): δ_A ≤ 16, S_N ≤ ½ ln(2πe·1024.08) = **4.88**, S_N^{(0)} ≤ ln 65 = **4.17**.
- [N] check_T1.py (L = 8, d_u = 3, f = 0.6):
  - exact δ_A = 0.824, against the circuit bound 2.32;
  - max over 40 times of ‖N_A(t) − N_A‖ is 1.245 ≤ 2δ_A = 1.648;
  - the support of p_t stays within n₀ ± 2.

**(b) Quasi-local U(1) generators.** Assume

(K) for all x, y and s: Σ_{X∋x,y} ‖k_X(s)‖ ≤ J F(|x−y|), with F non-increasing.

Then δ_A ≤ J |∂A| Σ_{r≥1} r F(r).

*Proof.*
1. m_X(A) ≤ |X∩A|·|X∖A|.
2. So Σ_X m_X‖k_X‖ ≤ Σ_{x∈A, y∉A} Σ_{X∋x,y} ‖k_X‖ ≤ J Σ_{x∈A, y∉A} F(|x−y|).
3. Each such pair is separated by at least one bond of ∂A, and a given bond separates exactly r pairs at distance r. ∎

Special cases:
- F(r) = e^{−r/ξ}: Σ r e^{−r/ξ} = [4 sinh²(1/2ξ)]⁻¹ ≤ ξ², so **δ_A ≤ J ξ² |∂A|**.
- Terms supported on intervals with ‖k_{[a,b]}‖ ≤ J₀e^{−(b−a)/ξ}: (K) holds with J = J₀(1 − e^{−1/ξ})⁻² ≤ J₀(1+ξ)².
- If the decay ‖k_X‖ ≤ J₀ e^{−diam X/ξ} is imposed on *arbitrary* subsets X, (K) holds only for ξ < 1/ln 2. The reason is that there are about 2^{D−1} subsets of diameter D through two given points.

**(c) Random products of local U(1) rotations, including resonant blocks.**
- Let U = O_K⋯O_1 be random. Each O_j is U(1)-invariant on X_j, and ‖O_j − 1‖ is arbitrary (e.g. exact diagonalisation of a resonant block).
- Deterministically, δ_A ≤ Ξ_A := Σ_j m_{X_j}(A) min{1, ‖O_j − 1‖}.
- If the X_j are intervals, Ξ_A ≤ Σ_{b∈∂A} Ξ_b, where Ξ_b := ½ Σ_{j: X_j ⊃ b} |X_j| min{1, ‖O_j−1‖}. Hence E Ξ_A² ≤ |∂A|² max_b E Ξ_b² by Minkowski.
- Consequences:
  - E sup_t Var N_A(t) ≤ E(σ₀ + 2Ξ_A)².
  - E sup_t S_N(t) ≤ E ½ ln[2πe((σ₀+2Ξ_A)² + 1/12)] ≤ ½ ln[2πe((σ₀ + 2‖Ξ_A‖_{L²})² + 1/12)] by Jensen.
  - Tail bound: P(sup_t S_N > s) ≤ P(Ξ_A > ½(√(e^{2s}/2πe − 1/12) − σ₀)).
- The S_N bound needs only **E ln(1 + Ξ_A) < ∞**. The variance bound needs E Ξ_A² < ∞.
- When is sup_L E Ξ_b² < ∞? It holds whenever the non-resonant rotations covering b at scale D have norms summable against D (e.g. ≤ ε^D, with O(D) of them per scale), and b lies in a resonant block of size ≥ D with probability decaying faster than D^{−3}.
  - This is the structure of Imbrie-type constructions [memory]. However, Imbrie's model has no U(1) symmetry, and I have not checked his estimates against this hypothesis.
  - So (c) is a theorem about a *hypothetical* U(1) analogue of his construction.

**Theorem 3.3 (U(1)-invariance of U is not needed) [P].** Let D be diagonal, U unitary, and H = UDU†. Assume

(I) the LIOMs τ_x := U n_x U† commute with N for all x. Equivalently, U†NU is diagonal.

(I) holds automatically if [H, N] = 0 and D has simple spectrum.

Let m_x := U†n_xU and η(j) := max_x ‖m_x − Π_{B_j(x)}m_x‖. Then for all t and ρ:

‖N_A(t) − N_A‖ ≤ 2‖B_A‖, ‖B_A‖ ≤ 12|∂A| Σ_{j≥0}(j+1)η(j), where B_A := Σ_{x∈A}(m_x − E_diag m_x).

Hence Var_ρ N_A(t) ≤ (σ₀ + 2‖B_A‖)², the S_N bound of Theorem 3.1(c) holds with δ_A → ‖B_A‖, and Var_{U|c⟩}(N_A) ≤ ‖B_A‖². If U is a circuit of light-cone radius R, then ‖B_A‖ ≤ 2R|∂A|.

*Proof.*
1. (I) is automatic. [D, U†NU] = U†[H,N]U = 0, and simple spectrum then forces U†NU to be diagonal. Conversely, U†NU is diagonal iff [U†NU, n_x] = 0 for all x, iff [N, τ_x] = 0 for all x.
2. Exact identity. Put f_x := E_diag m_x and g_x := m_x − f_x. Then U†N_AU = F_A + B_A with F_A := Σ_{x∈A} f_x diagonal. Hence [D, F_A] = 0 and
   N_A(t) − N_A = U(B_A(t) − B_A)U†.
3. Apply Lemma 2.5 to g_x:
   - tr g_x = tr m_x − tr E_diag m_x = 0.
   - Σ_x g_x = U†NU − E_diag(U†NU) = 0 by (I).
   - ‖g_x − Π_{B_j}g_x‖ = ‖(1−E_diag)(m_x − Π_{B_j}m_x)‖ ≤ 2η(j), because E_diag commutes with Π.
   - ‖g_x‖ ≤ 1, because m_x − ½ has norm ½. ∎

[N] check_T1.py tests this with U = V·C, where V is a U(1) brickwork and C a random classical reversible circuit of X, CNOT and Toffoli gates (C does **not** conserve N). Results:
- ‖[H,N]‖ = 0 and ‖[U,N]‖ = 6;
- ‖Σ_x g_x‖ = 6·10⁻¹⁵;
- max_t ‖N_A(t) − N_A‖ = 0.991 ≤ 2‖B_A‖ = 0.996.

**Remark 3.4 (Is a U(1)-invariant U an extra assumption beyond "FMBL + U(1)-symmetric H"?)**

1. *For the conclusions: no, up to degeneracies.* Theorem 3.3 uses only (I). (I) is automatic for simple spectrum of D, which is generic: a.s. for continuous i.i.d. l-bit fields. The price is the constant: 12Σ(j+1)η(j), or 2R, instead of δ_A.
2. *"U symmetric" is a gauge choice.*
   - U ↦ UW, with W a configuration permutation times diagonal phases, preserves the form H = U D′ U† with D′ = W†DW diagonal.
   - Example: W = ∏_{x∈S}σ^x_x turns a symmetric U into a non-symmetric one. The LIOMs stay U(1)-invariant (τ_x ↦ 1 − τ_x).
   - The invariant content is the U(1)-invariance of the LIOMs.
3. *Can one always re-choose U symmetric and quasi-local?*
   - This is equivalent to finding quasi-local projector LIOMs τ′_x with Σ_x τ′_x = N.
   - If U†NU = Σ_x n_x∘π for a local reversible relabelling π (a finite-depth classical reversible circuit or reversible cellular automaton), then yes: take U′ = UW_π^{−1}.
   - Whether every quasi-local integer-valued "isospectral counting function" is of this form is a combinatorial question I cannot settle.
   - In every constructive scheme U comes out symmetric automatically. These schemes are Schrieffer–Wolff/Wegner flows with U(1)-invariant generators, and Imbrie-type schemes in which each rotation either is generated by U(1)-invariant terms or diagonalises a U(1)-invariant block Hamiltonian H_B. In the latter case block eigenvectors can always be chosen inside N_B-sectors, because [H_B, N_B] = 0 and the sector dimensions match.
4. *Degeneracies.* If D is degenerate, [H, N] = 0 only makes U†NU block-diagonal with respect to the eigenspaces of D. Then (I) may fail for a given U, and Theorem 3.3 does not apply as stated. Generic random l-bit Hamiltonians have this with probability 0. Exact non-abelian symmetries are a different matter: they are believed to forbid FMBL (Potter–Vasseur [memory]). This is not the U(1) case.

**Remark 3.5 (Which LIOM tails keep C finite?) [P] + [N].**

The sufficient conditions are first moments: Σ_j (j+1) η(j) < ∞ (Thm 3.3) and Σ_r r F(r) < ∞ (Cor 3.2(b)).
- Exponential tails η(j) ≲ e^{−j/ξ} give C ∝ ξ².
- Stretched exponential e^{−(j/ξ)^β} gives C ∝ ξ²Γ(2/β)/β.
- Power laws η(j) ≲ j^{−β} require β > 2.

These are **not necessary**. Example [N] gaussian_tails.py: take a Gaussian U(1) dressing u = e^{iθ}, θ_xy = 0.3|x−y|^{−α}, with a Néel state and a half cut, L = 64…1024. Heuristically the LIOM tails are η(j) ~ j^{1/2−α}, since τ_x = c†(ue_x)c(ue_x) and ‖ue_x restricted outside B_j(x)‖ ~ (Σ_{r>j} r^{−2α})^{1/2}.

| α | Var_ψ(U†N_AU), L = 64 → 1024 | ‖Δ‖ = ‖U†N_AU − N_A‖ | first moment Σ_{A×A^c}\|θ\| |
|---|---|---|---|
| 0.75 | 0.38 → 0.81 (power-law growth) | 1.13 → 2.13 (grows) | diverges |
| 1.0 | 0.186 → 0.283 (+0.024 per doubling: logarithmic) | 0.78 → 1.17 (log) | diverges |
| 1.5 | 0.0926 → 0.0935 (converged) | 0.461 → 0.489 (converging) | diverges |
| 2.0 | 0.0796 | 0.3643 | diverges (log) |
| 3.0 | 0.0760 | 0.3109 | converges |

Reading the table:
- In this family the variance and ‖Δ‖ stay bounded for α ≥ 1.5 (tested) and diverge for α ≤ 1 (tested). The heuristic Σ_r r θ_r² < ∞ puts the threshold at α = 1, i.e. η(j) ~ j^{−1/2}. That is far weaker than either sufficient condition: α > 2 via the generator (Cor 3.2(b)), or α > 5/2 via the LIOM tails (Thm 3.3, β > 2).
- At and below the threshold the variance *does* diverge. For Gaussian states this forces S_N → ∞ (Theorem 5.1(c)).
- So: **T1 survives exponential and stretched-exponential tails, and power-law LIOM tails faster than j^{−2} in operator norm. In the Gaussian example the conclusion (bounded Var, bounded S_N) survives for much slower decay: tested at hopping tails r^{−1.5}, heuristically for any decay faster than r^{−1}. It genuinely fails for hopping tails ~ r^{−1} or slower.**

**Remark 3.6 (random constants).** For random l-bit data, the S_N statement needs only a logarithmic moment of the constant (Cor 3.2(c)). The variance statement needs a second moment. S_N is therefore an intrinsically *weak* diagnostic of LIOM structure.

---

## 4. T2: reconciling T1 with Aceituno Chávez et al. (arXiv:2312.13420)

### 4.1 Their model is in the T1 class

According to the literature report [snippet]:
- H = U H′ U†, where H′ is diagonal in the τ^z basis, with random fields and interactions decaying exponentially with distance.
- U is a finite-depth circuit of two-site gates, "each built to preserve U(1)". The depth is d_u = 16 for the large runs.
- They quench from Néel states, evolve exactly in the l-bit basis, and use L = 12–28.

Consequences:
- This is class (a) of Cor. 3.2, with hypotheses verified at the level of the abstract. Deterministically, for every realisation, every t and every L:
  - ‖N_A(t) − N_A‖ ≤ 32;
  - S_N(t) ≤ 4.88;
  - all Rényi S_N^{(α)}(t) ≤ ln 65 = 4.17, Hartley included.
- At their sizes the trivial bound ln(L/2 + 1) ≤ ln 15 = 2.71 is stronger. So **T1 places no constraint at L ≤ 28**, and nothing they report can contradict it.
- What T1 does imply is lim sup_{L→∞} sup_t S_N(t) ≤ 4.88 for every realisation. Any growth of the saturation value with L must therefore stop.
- The light cone is too wide to bite at their sizes (4d_u − 2 = 62 sites per side). The effective constant is set by the near-identity gates through min{1, ‖O − 1‖} in Cor 3.2(a). Its value depends on f, which I could not check.

### 4.2 An exact structural fact: the dephased distribution is local and independent of D

**Proposition 4.1 [P].** Setting:
- U = G_T⋯G_1, where each layer G_k is a product of U(1)-invariant gates on disjoint nearest-neighbour bonds;
- D is diagonal with simple spectrum;
- ψ is a Fock state with N_Aψ = n₀ψ;
- A = {1,…,ℓ}, with cut bond b = (ℓ, ℓ+1).

Then p̄(n) := lim_{τ→∞} τ⁻¹∫₀^τ p_t(n) dt exists, and:

- (i) p̄(n) = Σ_c |⟨c|U†ψ⟩|² ⟨c|U†P_n^A U|c⟩. This is **independent of D**.
- (ii) p̄(n₀ + ·) depends only on the gates contained in R := {x : dist(x,b) ≤ 2T−1} and on ψ|_R. For the brickwork of §1 (layers alternate, the first contains b, T = 2d_u), R can be taken to be the 4d_u − 2 sites on each side of b.
- (iii) Suppose the gate ensemble near the cut and the initial-state ensemble are defined relative to the cut. Then the law of p̄ is exactly L-independent once ℓ and L − ℓ are ≥ 4d_u − 2, i.e. L ≥ 8d_u − 4 for the half chain. This applies in particular to E H(p̄) and E Var_{p̄}. The law of p̄ is also independent of the law of D.
- (iv) For every realisation:
  - lim sup_τ τ⁻¹∫₀^τ S_N(t)dt ≤ H(p̄), by concavity of Shannon entropy and Jensen;
  - Var_{p̄} ≤ E_{p̄}(n−n₀)² ≤ 4δ_A².

*Proof.*

(i) p_t(n) = ⟨ψ|U e^{iDt}(U†P_n^AU) e^{−iDt}U†|ψ⟩. Apply Lemma 2.4.

(ii)
1. Let φ := U†ψ. The generating function is χ̄(λ) := Σ_n p̄(n)e^{iλn} = Σ_c |φ_c|² ⟨c|U†e^{iλN_A}U|c⟩.
2. Let U^λ := e^{−iλN_A}Ue^{iλN_A}. Every gate not on b commutes with N_A, so U^λ differs from U only in the gates on b, g ↦ e^{−iλn_ℓ} g e^{iλn_ℓ}.
3. Then U†e^{iλN_A}U = e^{iλN_A}W_λ, with W_λ := (U^λ)†U.
4. Peel W_λ = G₁^{λ†}⋯G_T^{λ†}G_T⋯G₁ from the inside. This gives W_λ = V₁(G₁†V₂G₁)(G₁†G₂†V₃G₂G₁)⋯, where V_k = G_k^{λ†}G_k is supported on b. So W_λ is supported in S := {dist(x,b) ≤ T−1}.
5. Hence χ̄(λ) = E_{c∼|φ|²}[e^{iλN_A(c)} ⟨c|W_λ|c⟩], and ⟨c|W_λ|c⟩ depends only on c_S.
6. The joint law of (c_S, N_A(c)) under |φ|² is fixed by E[X_S e^{iμN_A}] for diagonal X_S on S. Compute:
   ⟨φ|X_S e^{iμN_A}|φ⟩ = ⟨ψ|U e^{iμN_A}X_S U†|ψ⟩ = e^{iμn₀}⟨ψ|U^μX_SU†|ψ⟩ = e^{iμn₀}⟨ψ|U(W_μ†X_S)U†|ψ⟩.
7. The operator W_μ†X_S is supported in S. Conjugating it by the depth-T circuit gives an operator supported in R that is built only from gates inside R. Its expectation in the product state ψ depends only on ψ|_R.
8. For the brickwork, an explicit count gives the sharper R: S has 2T − 2 sites; the first conjugating layer does not enlarge S, because both edge bonds of S belong to the other sublattice of bonds; each later layer enlarges it by one site per side.

(iii) and (iv) follow immediately. ∎

[N] check_T1.py and check_T2_threshold.py verify this exactly:
- p̄ agrees to 10 digits for d_u = 1 at L = 4, 6, 8, 10, and for d_u = 2 at L = 12, 14.
- For d_u = 2 it differs at L = 6, 8, 10, below the threshold 8d_u − 4 = 12.
- Two very different D (W = 1, κ = 0.5 versus W = 2, κ = 1.5) give the same time-averaged p, matching the exact p̄ to 4 digits.

**Interpretation [P for the "size" part; H for the "time" part].**
- The interaction tails of D do **not** enter p̄ at all. They set only (1) the time scales on which p_t approaches p̄, and (2) the size of late-time temporal fluctuations of p_t at finite L. This makes the literature's "two tails" remark precise: the U-tails fix the size of the fluctuations, the D-tails fix the clock.
- **Conjecture 4.2 [H].** Suppose that in infinite volume D couples the cut region to infinitely many l-bits with generic couplings, and that U†ψ has non-degenerate far-l-bit statistics. Then p_t → p̄ as t → ∞: the far l-bits act as a dephasing bath for the number coherences near the cut. Hence S_N(t) → H(p̄).
  - At finite L, Ŝ_N(L) < H(p̄), and the gap closes as L grows.
  - "Saturation increasing with L towards an L-independent limit from below" is therefore exactly what a **localised** l-bit model predicts.
  - Gap to a proof: one needs decay of Σ_{c_out} ρ(c_out) e^{itΣ_j J_{ij}(c_i−c′_i)c_j}-type averages for the *joint* (not product) far-l-bit statistics of U†ψ.

### 4.3 Numerics [N] (lbit_run.py, analyze_lbit.py)

**Model.**
- Open chain with half cut.
- U: U(1) brickwork aligned with the cut, d_u gates per bond. Each gate is exp(−i f M), with M a random U(1)-symmetric Hermitian 4×4 matrix (GUE-like blocks).
- D_c = Σ_X J_X ∏_{i∈X} c_i over **all** subsets X:
  - J_{{i}} = h_i ~ U[−1, 1];
  - J_X = u_X e^{−(diam X−1)/κ} 2^{−(|X|−2)} with u_X ~ U[−1,1].
  - All k-body terms are included, so the D_c are a.s. rationally independent, and infinite-time averages equal averages over independent uniform phases (Kronecker–Weyl).
- Initial state: a Néel state, one of the two chosen at random.
- Observables:
  - S_N(t) on a log grid t = 10^{−1}…10^{13}, with phases computed in extended precision;
  - Ŝ_N := infinite-time average of S_N, estimated as the mean over 96 random phase vectors;
  - E H(p̄), where H(p̄) is the plug-in entropy of the phase-averaged p (slightly biased low, by about 0.006 for set C against the exact value).

**Set A** (d_u = 8, f = 0.25, κ = 0.5; gate depth exceeds the half chain, as in Chávez et al.). Disorder-averaged S_N(t):

| L | samples | t=1 | 10 | 10² | 10⁴ | 10⁸ | 10¹³ | Ŝ_N (∞-time) | E H(p̄) | E[H(p̄)−Ŝ_N] |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 4000 | 0.354 | 0.596 | 0.604 | 0.606 | 0.607 | 0.606 | 0.6055(23) | 0.6240(24) | 0.0185(2) |
| 10 | 4000 | 0.361 | 0.602 | 0.611 | 0.614 | 0.614 | 0.616 | 0.6152(23) | 0.6248(24) | 0.0096(1) |
| 12 | 2000 | 0.358 | 0.602 | 0.612 | 0.619 | 0.617 | 0.620 | 0.6193(33) | 0.6240(34) | 0.0048(1) |
| 14 | 1000 | 0.361 | 0.600 | 0.615 | 0.623 | 0.622 | 0.621 | 0.6217(49) | 0.6243(49) | 0.0026(1) |
| 16 | 400 | 0.357 | 0.603 | 0.616 | 0.628 | 0.629 | 0.627 | 0.6272(75) | 0.6285(75) | 0.0013(1) |

Further facts for set A:
- E Var_{p̄} = 0.207–0.210 for all L.
- The circuit constant is ⟨δ_circ⟩ = 3.6, so the T1 bound is ⟨½ ln(2πe(4δ²+1/12))⟩ = 3.39. The largest sup_t S_N over all samples is 1.12. The minimum over samples of (bound − sup_t S_N) is 2.07.

**Set C** (d_u = 1, f = 0.8, κ = 0.25; exact L-independence of the law of p̄ for L ≥ 4):

| L | samples | t=1 | 10 | 10² | 10⁴ | 10⁸ | 10¹³ | Ŝ_N | E H(p̄) (MC) | E[H(p̄)−Ŝ_N] |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 4000 | 0.225 | 0.446 | 0.458 | 0.464 | 0.462 | 0.462 | 0.4633(29) | 0.4805(30) | 0.0172(2) |
| 10 | 4000 | 0.227 | 0.447 | 0.458 | 0.469 | 0.471 | 0.470 | 0.4694(30) | 0.4799(31) | 0.0105(2) |
| 12 | 2000 | 0.236 | 0.444 | 0.461 | 0.474 | 0.477 | 0.479 | 0.4779(42) | 0.4828(43) | 0.0049(1) |
| 14 | 1000 | 0.234 | 0.454 | 0.456 | 0.475 | 0.483 | 0.482 | 0.4825(58) | 0.4857(59) | 0.0032(1) |
| 16 | 400 | 0.224 | 0.447 | 0.452 | 0.467 | 0.475 | 0.477 | 0.4767(93) | 0.4782(93) | 0.0015(1) |

Further facts for set C:
- Exact E H(p̄), computed from 4000 dense L = 8 samples (exact_pbar_C.py; the law is L-independent here): **0.4867(30)**. E Var_{p̄} = 0.150.
- After the fast rise (t ≲ 10), S_N grows slowly. The slope d⟨S_N⟩/d ln ln t is about 0.014–0.019 on [10, 10⁴] for every L. On [10⁴, 10⁸] it is ≈ 0 (−0.003, +0.002) for L = 8, 10 and 0.005–0.013 for L = 12–16. The growth window therefore lengthens with L, and the saturation value rises with L.

**Set B** (d_u = 8, f = 0.4, κ = 0.25; stronger dressing):

| L | samples | t=1 | 10 | 10² | 10⁴ | 10⁸ | 10¹³ | Ŝ_N | E H(p̄) | E[H(p̄)−Ŝ_N] |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 4000 | 0.536 | 0.816 | 0.821 | 0.822 | 0.825 | 0.823 | 0.8221(15) | 0.8398(16) | 0.0177(1) |
| 10 | 4000 | 0.563 | 0.829 | 0.832 | 0.836 | 0.838 | 0.837 | 0.8380(16) | 0.8458(16) | 0.0078(1) |
| 12 | 2000 | 0.565 | 0.833 | 0.838 | 0.843 | 0.845 | 0.843 | 0.8450(23) | 0.8484(23) | 0.0034(0) |
| 14 | 1000 | 0.566 | 0.827 | 0.832 | 0.840 | 0.842 | 0.841 | 0.8406(33) | 0.8421(33) | 0.0015(0) |
| 16 | 400 | 0.570 | 0.838 | 0.831 | 0.847 | 0.849 | 0.849 | 0.8497(50) | 0.8503(50) | 0.0007(0) |

Further facts for set B:
- ⟨δ_circ⟩ = 5.6, so the T1 bound is 3.83. The largest sup_t S_N over all samples is 1.29.
- E H(p̄) at L = 8 lies 2.6σ below L = 10. This is the expected drift of p̄_L itself below the light-cone threshold (8d_u − 4 = 60), visible here because the dressing is stronger. From L = 10 on, E H(p̄) is flat within errors.
- The gap again shrinks geometrically, by a factor of about 2.2 per two sites.

**What the numerics show.**
1. A model that is *exactly* in the T1 class reproduces the two qualitative features reported by Chávez et al.:
   - slow growth of S_N over many decades after the initial rise;
   - a long-time value Ŝ_N(L) that **increases with L**.
2. In the same model:
   - the per-realisation inequality Ŝ_N ≤ H(p̄) is a theorem (Prop 4.1(iv)). The estimator satisfies it automatically, because both sides use the same phase samples; this is not independent evidence;
   - E H(p̄) does not depend on L: exactly by Prop 4.1 in set C (L ≥ 4), and within errors in set A and in set B from L = 10 on (where exact independence would need L ≥ 60);
   - the gap E[H(p̄) − Ŝ_N] **halves with every two added sites** (0.0185 → 0.0096 → 0.0048 → 0.0026 → 0.0013 in set A; 0.0177 → 0.0078 → 0.0034 → 0.0015 → 0.0007 in set B; similarly in set C).
   So the L-dependence of the saturation value is a finite-size transient converging exponentially fast to an L-independent limit. That limit lies far below the T1 bound.
3. The slow-growth amplitudes here (≈ 0.01–0.03 in S_N after t ≈ 10) are smaller than those reported for microscopic XXZ data. I did not tune parameters to match. The argument is structural, not a fit.
4. **Verdict.** There is no contradiction. Chávez et al.'s L-dependence is a finite-size effect that must stay bounded: deterministically by T1, and in the dephased sense by E H(p̄), which is L-independent for L ≥ 8d_u − 4.
   - For their d_u = 16 (and a cut-aligned brickwork), the law of p̄ becomes *exactly* L-independent only at L ≥ 124.
   - Below that, p̄_L itself can drift. The drift comes from gates near the chain ends inside the light cone, and is exponentially small in L/ξ_τ for near-identity gates.
   - This drift could produce the reported non-monotonic L-dependence (rise up to L ≈ 16, then slow decrease [snippet]). I cannot verify their figure.

---

## 5. T3: Gaussian (quasi-free) states

**Theorem 5.1 [P].** Let ω be a gauge-invariant quasi-free state (pure or mixed) of fermions on Λ. Let C_A be its correlation matrix on A, with C_xy = ω(c_y†c_x), and let ν₁,…,ν_{|A|} ∈ [0,1] be its eigenvalues. Write Var := Var_ω N_A, S := S(ρ_A), and S^{(2)} := −ln Tr ρ_A².

- (a) N_A is Poisson-binomial with parameters ν_k: ω(e^{iλN_A}) = det(1 − C_A + C_A e^{iλ}) = ∏_k (1 − ν_k + ν_k e^{iλ}). Hence Var = Σ_k ν_k(1−ν_k) = Tr C_A(1−C_A).
- (b) S = Σ_k h(ν_k) ≥ 4 ln 2 · Var, and 2 Var ≤ S^{(2)} ≤ 4 ln 2 · Var.
- (c) ½ ln(8 Var/π) ≤ −ln max_n p(n) ≤ S_N ≤ ½ ln[2πe(Var + 1/12)].
- (d) Consequently:
  - S_N ≤ ½ ln[2πe(S/(4 ln 2) + 1/12)];
  - S_N ≤ ½ ln[πe(S^{(2)} + 1/6)];
  - S_N ≥ ½ ln[2S^{(2)}/(π ln 2)].

*Proof.*

(a) Write C_A = Σ_k ν_k|v_k⟩⟨v_k| and b_k := c(v_k). Then ρ_A = ⊗_k[(1−ν_k)|0⟩⟨0| + ν_k|1⟩⟨1|] and N_A = Σ_k b_k†b_k.

(b) Binary entropy.
- With y := (1−2ν)²:
  h(ν) = ln 2 − Σ_{k≥1} y^k/(2k(2k−1)) ≥ ln 2 − y Σ_{k≥1} 1/(2k(2k−1)) = ln 2 (1−y) = 4 ln 2 · ν(1−ν).
  This uses Σ_{k≥1} 1/(2k(2k−1)) = ln 2. The inequality is also Topsøe 2001 [memory].
- For Rényi-2, each mode contributes −ln(1 − x) with x = 2ν(1−ν) ∈ [0, ½]. By convexity, x ≤ −ln(1−x) ≤ 2 ln 2 · x.

(c)
- Upper bound: Lemma 2.1.
- Lower bound:
  1. |ω(e^{iθN_A})|² = ∏_k [1 − 2ν_k(1−ν_k)(1−cos θ)] ≤ e^{−2Var(1−cos θ)}.
  2. On [−π, π], 1 − cos θ ≥ 2θ²/π².
  3. Hence max_n p(n) ≤ (2π)⁻¹∫_{−π}^{π}|ω(e^{iθN_A})| dθ ≤ (2π)⁻¹∫_ℝ e^{−2Var θ²/π²}dθ = √(π/8)/√Var.

(d) Combine (b) and (c). ∎

[N] check_inequalities.py:
- min_ν h/(ν(1−ν)) = 2.7726 = 4 ln 2, attained at ν = ½;
- −ln(1−x)/x ranges over [1, 2 ln 2];
- over 5000 random Poisson-binomial laws: min(H_∞ − ½ln(8Var/π)) = 0.399 ≥ 0 and max(H − upper) = −0.0012 ≤ 0;
- empirically, max p ≤ 0.67·√(π/8)/σ.

**Proposition 5.2 (entropy asymptotics) [P-c].** For Poisson-binomial laws, S_N − ½ ln(2πe Var) → 0 as Var → ∞, uniformly in the parameters ν_k. Numerically the error is O(1/Var): it is −8·10⁻¹⁰ at Var = 2500 with ν ≡ ½, and −1.4·10⁻⁴ at Var = 475 with ν ≡ 0.05.

*Proof.* Let σ² = Var and μ = mean, and let g_k := (2πσ²)^{−1/2} e^{−(k−μ)²/2σ²}.

1. *Local CLT: sup_k |p_k − g_k| ≤ C₀/σ².*
   - Let f_ν(θ) := ln(1−ν+νe^{iθ}) − iνθ. Then f″ = −w(1−w) and f‴ = −iw(1−w)(1−2w), where w = νe^{iθ}/(1−ν+νe^{iθ}).
   - For |θ| ≤ π/2: |1−ν+νe^{iθ}|² ≥ ½, so |f‴| ≤ 8ν(1−ν). Hence |ln φ̃(θ) + σ²θ²/2| ≤ (4/3)σ²|θ|³, where φ̃(θ) = E e^{iθ(N_A−μ)}.
   - On |θ| ≤ 3/16, |φ̃ − e^{−σ²θ²/2}| ≤ (4/3)σ²|θ|³e^{−σ²θ²/4}. Its integral is 64/(3σ²).
   - On 3/16 ≤ |θ| ≤ π, use |φ̃| ≤ e^{−2σ²θ²/π²}. That contribution is exponentially small in σ².
2. *Central region.* Let K := {|k−μ| ≤ σ√(ln σ)}.
   - On K, g_k ≥ (2π)^{−1/2}σ^{−3/2}, so p_k = g_k(1 + O(σ^{−1/2})).
   - Hence Σ_K −p_k ln p_k = Σ_K p_k[½ ln(2πσ²) + (k−μ)²/2σ²] + O(σ^{−1/2}).
3. *Tails.* Bernstein's inequality P(|N_A−μ| ≥ s) ≤ 2e^{−s²/(2(σ²+s/3))} gives:
   - P(K^c) ln σ → 0;
   - E[(N_A−μ)²/σ²; K^c] → 0;
   - Σ_{K^c} −p ln p → 0. For this, split at p = σ^{−3}, and use −p ln p ≤ (2/e)√p below the split.

Summing, S_N = ½ ln(2πσ²) + ½ + o(1). ∎

(Background: Petrov, *Sums of Independent Random Variables*, Ch. VII, for lattice local limit theorems [memory].)

**Corollary 5.3 (Anderson insulators, time-uniform bound) [P].** Setting: a random one-particle h on Λ, and a Fock initial state with occupied set Z. Assume

(EL) E sup_t |⟨x|e^{−iht}|z⟩| ≤ C₁ e^{−|x−z|/ξ}, uniformly in Λ.

Then

E sup_t Var N_A(t) ≤ C₂|∂A|, with C₂ = C₁ Σ_{r≥1} r (r + 1 + 2/(e^{1/ξ}−1)) e^{−r/2ξ},

and E sup_t S_N(t) ≤ ½ ln[2πe(C₂|∂A| + 1/12)].

*Proof.*
1. Var_t = Σ_{x∈A, y∉A} |C_xy(t)|², with |C_xy|² ≤ |C_xy| ≤ Σ_{z∈Z} |U_xz|^{1/2}|U_yz|^{1/2}.
2. By Cauchy–Schwarz in expectation, E sup_t(·) ≤ C₁ Σ_z e^{−(|x−z|+|y−z|)/2ξ}.
3. Sum over pairs as in Cor 3.2(b). ∎

This is ANSS-type (Abdul-Rahman–Nachtergaele–Sims–Stolz bound entanglement, and S_N ≤ S) [memory]. The direct version needs only (EL).

**5.4 The random-hopping chain (the user's code) [P + N].**

Setting: J_i = u_i^δ, zero potential. At E = 0 the chain is critical, so (EL) fails there and Cor 5.3 does not apply. Theorem 5.1 applies at every time, deterministically. Predictions:
- If S(t) ~ a ln ln t, then Var(t) ≤ (a/4 ln 2) ln ln t + O(1). Hence S_N(t) ≤ ½ ln ln ln t + O(1).
- If S ~ a ln t, then S_N ≤ ½ ln ln t + O(1).
- Conversely, S_N → ∞ iff Var → ∞, and then S_N = ½ ln(2πe Var) + o(1).
- At t = 10¹², ln ln t = 3.32, so ½ ln ln ln t ≈ 0.6: number-entropy growth in this model is essentially invisible on any accessible window.

[N] random_hopping_T3.py. Setting as in the code: a ring, A = the first L/2 sites (two cuts), and the code's initial state with occupied sites 4, 8, …, L. L = 256.

| δ | t | ⟨S⟩ | ⟨Var⟩ | ⟨S_N⟩ | ½ln(2πe(⟨Var⟩+1/12)) | ½ln(2πe(⟨S⟩/4ln2+1/12)) |
|---|---|---|---|---|---|---|
| 1 | 10² | 2.10 | 0.638 | 1.150 | 1.256 | 1.332 |
| 1 | 10⁶ | 2.55 | 0.790 | 1.263 | 1.351 | 1.420 |
| 1 | 10¹² | 2.57 | 0.797 | 1.264 | 1.355 | 1.424 |
| 2 | 10² | 0.83 | 0.255 | 0.673 | 0.876 | 0.937 |
| 2 | 10⁶ | 1.15 | 0.360 | 0.831 | 1.012 | 1.069 |
| 2 | 10¹² | 1.22 | 0.384 | 0.882 | 1.039 | 1.096 |

Observations:
- The run used 60 samples for δ = 1 and 40 for δ = 2. Over all samples and times there were no violations of any of the three inequalities.
- At L = 256, finite-size saturation (t_L ~ e^{c δ√L}) sets in within the window for δ = 1.
- Remarks on /home/user/fermions/src. These are observations only; the code was not modified.
  1. `init_state` fills only u(2i, i) for *even* i. The initial state is therefore the period-4, quarter-filling density wave |0001 0001 …⟩ with L/4 particles, not a Néel state. The odd columns are zero, which is harmless for C = UU†. This may be intentional.
  2. `build_ham` writes hmat(ls, ls+1) and hmat(ls+1, ls) for i = ls, outside the array bounds. The wrap-around elements are then set separately when bc = .true.
  3. Var = Σν(1−ν) and S_N (by Poisson-binomial convolution of the ν's) can be obtained at negligible cost from the eigenvalues already computed in `entropy`.

**Proposition 5.5 (no sublinear S_N ≤ f(S) for interacting U(1) states) [P].**

- (a) Let ψ be pure with definite total N. Then S_N ≤ S (Wiseman–Vaccaro, Lukin et al. [memory]). Also S_N ≤ ln(min(|A|, |A^c|, N) + 1).
- (b) Suppose N ≤ |A^c|. For every m ≤ min(|A|, N) and every s ∈ [0, ln(m+1)], there is such a ψ with S_N = S = s. The construction:
  - choose q on {0,…,m} with H(q) = s;
  - let |a_n⟩ := |1^n 0^{|A|−n}⟩_A and |b_n⟩ := a fixed configuration of N − n particles in A^c;
  - set ψ = Σ_n √q_n |a_n⟩|b_n⟩. This is a Schmidt decomposition with N_A ↔ n bijective, so S_conf = 0.
  Hence the optimal universal f is f(s) = s, and **no sublinear f exists** as |A| → ∞.
  - Physical version: a rigid M-particle block (e.g. a magnon bound state of the large-Δ XXZ chain) in a superposition of positions straddling the cut also has S_N = S.
- (c) Without definite N, even S_N ≤ S fails. Example: a product state (|0⟩+|1⟩)^{⊗|A|} ⊗ … has S = 0 and S_N > 0.

The relation S_N ≲ ½ ln S is therefore a free-fermion theorem (Thm 5.1). It must fail asymptotically in a T1-class FMBL phase, where S ∝ ln t from dephasing [standard H] while S_N stays bounded.

---

## 6. T4: rigorous finite-time bounds for the strongly disordered XXZ / t-V chain

**Model.** H = H₀ + J T on any L (open or ring), where
- H₀ = Σ_x h_x n_x + JΔ Σ_x n_x n_{x+1}, with h_x i.i.d. U[−W, W];
- T = Σ_x τ_x, with τ_x = c_x†c_{x+1} + h.c.

(For hard-core bosons, or for fermions with open boundary conditions, T is a real symmetric matrix in the occupation basis.)

**A priori ballistic bound.** dN_A/dt = i[H, N_A] = J × (bond currents across ∂A). Each bond current has norm 1, so ‖N_A(t) − N_A‖ ≤ J|∂A||t|.

**Lemma 6.1 (approximately conserved dressed charge) [P].** Suppose Ñ = Ñ† satisfies ‖Ñ − N_A‖ ≤ δ and ‖[H, Ñ]‖ ≤ ε. Then for all t and ρ:
- ‖N_A(t) − N_A‖ ≤ 2δ + ε|t|;
- Var_ρ N_A(t) ≤ (σ₀ + 2δ + ε|t|)²;
- S_N(t) ≤ ½ ln[2πe((σ₀+2δ+ε|t|)² + 1/12)].

*Proof.* Write N_A(t) − N_A = [Ñ(t) − Ñ] + [Ñ − N_A] − [Ñ − N_A](t). Then use ‖Ñ(t) − Ñ‖ ≤ |t|·‖[H,Ñ]‖ and Lemmas 2.1–2.2. ∎

Theorem 3.1 is the case ε = 0 with Ñ = U N_A U†. Note that only **one** quasi-conserved operator, localised at ∂A, is needed. Resonances away from the cut are irrelevant.

**Theorem 6.2 [P].** Fix n ≥ 1 and γ ≥ 20nJ. There is an event Ω_{n,γ} with the following properties:
- it depends only on the fields within distance 2n of ∂A;
- P(Ω_{n,γ}) ≥ 1 − |∂A| 2^{8n−2} γ/W;
- on Ω_{n,γ} there is a self-adjoint Ñ_n with
  ‖Ñ_n − N_A‖ ≤ 2|∂A| J/γ and ‖[H, Ñ_n]‖ ≤ |∂A| J n! (10J/γ)^n.

Hence on Ω_{n,γ}, for all t and uniformly in L:

‖N_A(t) − N_A‖ ≤ |∂A| (4J/γ + n!(10J/γ)^n J|t|),

with the variance and S_N bounds of Lemma 6.1.

*Proof.*

*Construction.* Set Y₀ = N_A. For k ≥ 1 let R_k := J[T, Y_{k−1}] and ⟨c|Y_k|c′⟩ := −⟨c|R_k|c′⟩/(E_c − E_{c′}), where E_c is the H₀-energy. Then Ñ_n := Σ_{k=0}^n Y_k satisfies [H, Ñ_n] = J[T, Y_n]. (Telescoping, using [H₀, Y_k] = −R_k.)

*Well-definedness (reality).* By induction:
- Y_{k−1} is real symmetric.
- T is real symmetric, so R_k is real antisymmetric and has zero diagonal. The only obstruction to solving [H₀,Y] = −R would be a non-zero diagonal, so none arises.
- Y_k is again real symmetric, since R and the denominators are both antisymmetric.
- Non-degeneracy of the relevant denominators is supplied by Ω.

*Locality.* Consider one boundary bond b, and measure dist to b.
- R_k is supported in {dist ≤ 2k−2}.
- Y_k is supported in I_k := {dist ≤ 2k−1}, with |I_k| = 4k. The diagonal Δ-dependence adds one site per side.
- The number of bonds touching I_k is 4k+1.
- For several boundary bonds the recursion is linear: Y_k = Σ_b Y_k^{(b)}.

*Norms.* Use |||X||| := max_c Σ_{c′}|⟨c|X|c′⟩|. It is submultiplicative, and ‖X‖ ≤ |||X||| for (anti)symmetric X by the Schur test. We have |||τ_x||| = 1 and |||[τ_b, N_A]||| = 1. On Ω, |||Y_k||| ≤ |||R_k|||/γ. Hence
- |||Y_k||| ≤ (J/γ) ∏_{j=1}^{k−1} 2(4j+1)J/γ ≤ (J/γ)(10(k−1)J/γ)^{k−1} ≤ (J/γ) 2^{−(k−1)} for k ≤ n;
- ‖[H,Ñ_n]‖ ≤ J·2(4n+1)·|||Y_n||| ≤ J n!(10J/γ)^n.

*Probability.*
- The denominators met up to order n are E_c − E_{c′} for pairs (c, c′) that differ only on {dist ≤ 2n−2} (4n − 2 sites), with arbitrary values on the two sites at distance 2n−1. There are ≤ 2^{8n−2} such ordered pairs.
- If c ≠ c′, some h_x has coefficient ±1 in E_c − E_{c′}. Conditionally on the other fields, P(|E_c − E_{c′}| < γ) ≤ γ/W.
- A union bound over pairs and over b ∈ ∂A finishes the proof. ∎

**Corollary 6.3 (superpolynomial t\*) [P].** Let p ∈ (0,1) and x := pW/(10|∂A|J). Choose n ≈ ln x/(16 ln 2) and γ = pW/(|∂A| 2^{8n−2}). Then with probability ≥ 1 − p,

sup_{|t| ≤ t\*} ‖N_A(t) − N_A‖ ≤ |∂A|(1 + 4J/γ), where t\* ≥ J⁻¹ exp[(ln x)²/(32 ln 2) − O(ln x · ln ln x)].

On this time window S_N stays O(1), uniformly in L.

*Proof.* With 10J/γ = 2^{8n−2}/x, we get ln(1/(n!(10J/γ)^n)) = n ln x − n(8n−2) ln 2 − ln n!. Optimise over n. The constraint γ ≥ 20nJ amounts to x ≥ 2n 2^{8n−2}, which holds for large x. ∎

[N] check_T4.py (L = 10, W/J = 6, 12, 24):
- the diagonal of every R_k vanishes to machine precision, confirming the reality argument;
- ‖N_A(t) − N_A‖ ≤ 2δ_n + ε_n t holds in all cases;
- at these W the series is spoiled by resonances (smallest denominators ~10⁻³) after one to a few orders, as the theorem's failure probability predicts.

**Honest assessment.**
1. The constants are hopeless at physical disorder. With the stated constants, the n = 1 bound is non-vacuous with reasonable probability and beats the ballistic one only for W ≳ 10³J. (Counting only the three denominators that actually occur at n = 1 gives failure probability ≤ 3γ/W per bond, which improves this only modestly.) Cor 6.3 needs ln x ≫ 1. Better combinatorics would improve constants, not the structure: exp(c ln²(W/J)) is the natural output of fixed-order perturbation theory with a union bound.
2. The literature report says De Roeck–Huveneers–Meeus–Prośniak (2023) [snippet] prove quasi-conservation of a full set of quasi-local quantities up to exp{c(log W)^{2−ε}}. Theorem 6.2 is an elementary single-cut analogue of the same type. I have not checked whether their model class includes U(1) hopping.
3. *What is missing for an all-time (T1-type) unconditional statement:*
   - (i) A multiscale control of resonances at **all** orders near the cut, uniformly in time. This is an Imbrie-type construction for a U(1) chain together with a level-statistics hypothesis (LLA). It does not exist; Imbrie's proof covers a non-U(1) chain and is conditional on LLA.
   - (ii) Alternatively, a logarithmic light cone at positive energy density. Elgart–Klein have this only at the bottom of the spectrum [snippet].
   - Suppose (ii) holds in the form ‖n_x(t) − Π_{B_j(x)}n_x(t)‖ ≤ C e^{−(j−ℓ(t))/ξ′} with ℓ(t) = v ln t. Then Lemma 2.5 applied to g_x = n_x(t) − n_x (traceless, Σ_x g_x = 0) gives ‖N_A(t) − N_A‖ ≤ 6|∂A|(ℓ(t)² + O(ℓ)). Hence **S_N(t) ≤ 2 ln ln t + O(1)** rigorously. The heuristic linear charge-transfer count of the literature report gives ln ln t.
4. Beyond t\*, Lemma 6.1 yields only S_N ≲ ln(εt), i.e. ln t growth. That is the generic envelope of linear charge-transfer growth.

---

## 7. T5: logical consequences and sharp diagnostics

**Logic.**
- **(L1) [P]**
  - Suppose E sup_t S_N(t; L) → ∞ along L → ∞ (or, for L = ∞, as t → ∞). Then there is no family of representations H_L = U_L D_L U_L† (D_L diagonal, U(1)-invariant LIOMs) with sup_L E ln(1 + K_A(U_L)) < ∞. Here K_A = δ_A for symmetric U, or ‖B_A‖ ≤ 12|∂A| Σ(j+1)η(j) in general.
  - In particular, no such family has LIOM tails with uniformly summable first moment.
  - If E sup_t Var N_A → ∞, the same holds with E K_A² in place of the log-moment.
- **(L2) [P] The converse is false.** Bounded S_N does not imply localisation. Example: the U(1) charge is carried by immobile particles that couple only density–density to an ergodic neutral sector. Then S_N ≡ 0 for Fock states while S grows linearly. Fragmentation and strong-Stark (prethermal) systems give less trivial examples. Bounded S_N is necessary for T1-class FMBL, not sufficient.
- **(L3) [P + N] Finite windows decide nothing.**
  - Slow S_N growth over any finite window, ln ln t-like and with L-increasing saturation, is compatible with T1. Section 4 exhibits it in an exactly T1-class model.
  - Finite-(L, t) data can never falsify T1-class FMBL. Only an extrapolation in L can, and only under assumptions on convergence rates.
  - By Lemma 2.1, observed values S_N ≈ 0.3, 0.5 and 1.0 require only Var ≥ 0.023, 0.076 and 0.35, i.e. δ_A ≥ 0.08, 0.14 and 0.30. These are far below any natural bound.
- **(L4) [P] S_N is the weakest of the natural observables.**
  - S_N ≈ ½ ln(2πe Var) compresses: (ln t)^ν growth of Var appears as (ν/2) ln ln t in S_N.
  - For random LIOM constants, S_N needs only a log-moment, whereas Var needs a second moment (Rem 3.6).
  - Rényi α < 1 and Hartley are *not* controlled by Var. They are controlled only for strictly finite circuits (Cor 3.2(a)).
- **(L5) [P]** S_N ~ ½ ln S is a free-fermion theorem (Thm 5.1), with no interacting analogue (Prop 5.5). Its persistence in interacting data at accessible times is evidence neither for nor against MBL.

**Sharp numerical diagnostics** (each is bounded uniformly in L by a T1/T1′ constant):

- **(D1) Static eigenstate test [P].** M(L) := E max_{eigenstates E} Var_E(N_A), taken in a fixed-N sector.
  - T1/T1′ give M(L) ≤ E K_A² (Thm 3.1(d), Thm 3.3).
  - It needs no time evolution and is blind to dephasing time scales. It is accessible by full ED up to L ≈ 20.
  - Ergodic phase: M ~ L. Growth of M(L) beyond any constant (e.g. ∝ ln L) at strong disorder would rule out U(1)-invariant LIOMs with E K² < ∞.
  - Caveat: resonant pairs give Var_E ≈ ¼ at any L. Only the *scaling* is informative, and quantiles of the eigenstate distribution are more robust than the max. Bipartite eigenstate fluctuations have been studied before (Singh–Bardarson–Pollmann 2016 [memory]); the max/uniform-bound angle seems unexploited.
- **(D2) Dephased (diagonal-ensemble) statistics [P].** p̄_L(n) = Σ_E |⟨E|ψ⟩|²⟨E|P_n^A|E⟩, for simple spectrum.
  - Var_{p̄_L} ≤ 4K_A² (Prop 4.1(iv) / Thm 3.1(d)).
  - Use Var_{p̄_L} and H(p̄_L) instead of late-time S_N(t; L). Section 4.3 shows that Ŝ_N(L) < H(p̄_L) and that the gap closes with L in localised models. The rise of saturation values with L is precisely this artefact.
  - In circuit l-bit null models with matched local data, the law of p̄_L is exactly L-independent for L ≥ 8d_u − 4 (Prop 4.1). This gives a sharp calibration.
- **(D3) Worst case over initial states [P].** sup over Fock ψ and over t of Var_ψ N_A(t), or the operator norm ‖N_A(t) − N_A‖, or its infinite-temperature version 2^{−L}Tr(N_A(t) − N_A)². All are ≤ 4K_A².

**What numerics cannot decide.** Boundedness or unboundedness as L → ∞ from L ≤ 28, without convergence-rate assumptions. The T1 constant is not known a priori. Dephasing scales grow like e^{r/κ} and exceed accessible times after a few sites. The light cone of Chávez-type circuits (L ≥ 124 for d_u = 16) is far beyond ED.

---

## 8. Novelty assessment (relative to the two literature reports, which are snippet/memory-level)

- **T1 (Thm 3.1, Cor 3.2).**
  - The key lemma (charge moved by a U(1)-invariant locality-preserving unitary is a bounded boundary operator) is standard in charge-transport mathematical physics (Bachmann–Bols–De Roeck–Fraas; Kapustin–Sopenko [snippet]).
  - The integer max-entropy step is textbook.
  - The statement "LIOMs ⇒ bounded S_N" is folklore heuristic in the MBL literature (KUFS 2020–21, Luitz–Bar Lev 2020 [snippet]).
  - Apparently new as a written theorem with explicit constants. Also apparently new in detail: the support/Hartley bound for circuits, the eigenstate-variance corollary, and the log-moment remark. Assessment: *folklore made rigorous; low novelty, useful*.
- **T1′ (Thm 3.3) and Lemma 2.5.** The partial-trace "boundary lemma" for a locally decomposed zero-sum family is a standard style of argument. Its use to remove U(1)-invariance of the LIOM unitary (only invariant LIOMs needed, automatic for simple spectrum) appears new. *Moderate.*
- **Tails (Rem 3.5).** The first-moment criterion is elementary. The Gaussian numerics (bounded for hopping tails r^{−1.5}, divergent for r^{−1} and slower) are new here, but easy.
- **Prop 4.1 (D-independence and exact locality of the dephased distribution in circuit l-bit models) and the per-realisation bound Ŝ_N ≤ H(p̄).** Not found in the reports. Apparently new, and it directly sharpens the Chávez et al. discussion. The exponential closing of the gap is new numerics here. Conjecture 4.2 is open. *Moderate to high relevance for the debate; low technical difficulty.*
- **T3.**
  - S ≥ 4 ln 2·Var, the ½ ln(2πe(·)) upper bounds, and S_N ~ ln S for free fermions are known in spirit (KUFS SciPost 2020 [snippet]).
  - The explicit two-sided bounds, in particular ½ ln(8Var/π) ≤ S_N and the Rényi-2 sandwich 2Var ≤ S^{(2)} ≤ 4 ln 2·Var, are elementary and possibly new in this form.
  - Prop 5.2 is a routine local-CLT consequence.
  - Prop 5.5 (no sublinear f) is elementary but worth stating.
- **T4.** Lemma 6.1 is trivial. Theorem 6.2 is a self-contained single-cut version of standard perturbative quasi-conservation; the DHMP-type results [snippet] are stronger in scope. *Low novelty.* Its value is clarity about what is missing.
- **T5.** The logic is standard once T1 is stated. Diagnostics D1–D2, in the uniform-bound form, appear unexploited.

---

## 9. Reproducibility

All scripts are in `/tmp/claude-0/-home-user-fermions/e9e4eeb3-4dfa-5769-bb28-5161ded2b166/scratchpad/code/`:
- `check_inequalities.py`: scalar inequalities of Lemmas 2.1 and Thm 5.1, plus Poisson-binomial asymptotics.
- `lbit.py`: fixed-N sector, U(1) gates, cut-aligned brickwork, all-k-body l-bit energies.
- `check_T1.py`, `check_T2_threshold.py`: exact checks of Thm 3.1, Thm 3.3 and Prop 4.1.
- `lbit_run.py` and `analyze_lbit.py`: T2 production runs (sets A, B, C; results in `res_*.npz`). `exact_pbar_C.py`: exact E H(p̄) for set C.
- `gaussian_tails.py`: Rem 3.5.
- `random_hopping_T3.py`: §5.4.
- `check_T4.py`: Thm 6.2 recursion.
