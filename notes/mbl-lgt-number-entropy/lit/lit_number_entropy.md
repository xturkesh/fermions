# Number entropy in MBL: literature map (as of 2026-11-22)

## How this was sourced, and the tags used
- **arxiv.org and ar5iv.labs.arxiv.org could not be reached from this session (DNS blocked by the network allowlist), and the `web-fetch` subagent tool was not available.** So no paper page was opened directly. Everything below that is not from memory comes from search-engine excerpts of the arXiv, APS, SciPost or abstract pages. Those excerpts often quote abstract or body sentences word for word.
- Tags:
  - **[snippet]**: seen this session only as search-engine excerpts or summaries of the actual page. Treat it as "abstract-level read", not a full read.
  - **[memory]**: not checked this session.
  - **[this report]**: my own derivation, with its type given.
  - No item is tagged [read], because no full text was opened.
- Claim types: (P) proven, (N) numerical, (H) heuristic, (E) experiment.
- To do the requested careful read of 2312.13420 (the circuit Eqs. (8)–(10), f values, ξ_τ, and the S_N(L) figure), arxiv.org must be allowed in the session's network settings.

---

## (a) Bottom line, 2026

1. **The numerical fact is not disputed (N).**
   - After quenches from product states in the random-field XXZ / t-V chain, S_N(t) keeps growing slowly over all accessible windows.
   - It is well fit by S_N ≈ (ν/2) ln ln t, i.e. Var N_A ~ (ln t)^ν.
   - This holds for exact diagonalization (ED) up to L≈24, for TDVP at L=50 and W=10 (Sierant–Zakrzewski, per the 2026 review [snippet]), and for quasiperiodic potentials (Sharifan–Kiefer-Emmanouilidis–Sirker 2026 [snippet]).
   - The S_N values involved are small, a few tenths or less at strong disorder [memory].
2. **The interpretation is still split; nothing is settled.**
   - The Sirker group (2020–2026) reads the growth as slow particle delocalization and hence the absence of true MBL. In 2026 they attribute it to an energy-independent "path-entropy obstruction" rather than avalanches [snippet].
   - Luitz–Bar Lev (2020), Ghosh–Žnidarič (2022), Aceituno Chávez et al. (PRL 2024), Chen–Chen–Wang (PRB 2025?) and Faulend–Buljan–Štrkalj (2026) argue the growth is compatible with localization. Their mechanisms are: single-particle fluctuations across the cut, two-level resonances, l-bit tails plus dephasing, and many-body quantum beats.
   - The most decisive item is Chávez et al.: a model localized *by construction* reproduces ultraslow S_N growth and an L-dependent saturation value at L ≤ 28. Their own conclusion is that such observations are *not sufficient* to rule out MBL [snippet].
3. **What is rigorous is very little.**
   - No theorem covers interacting U(1) chains in either direction. Imbrie's LIOM proof is for a non-U(1) Ising-type chain under the limited-level-attraction assumption [memory].
   - Rigorous results that do exist:
     - S = S_N + S_conf.
     - The discrete max-entropy bound.
     - Free-fermion bounds relating S, S_N and Var N_A (Kiefer-Emmanouilidis et al. 2020).
     - Time-uniform bounds for Anderson / disordered XY chains via dynamical localization (Abdul-Rahman–Nachtergaele–Sims–Stolz 2016).
     - The *conditional* statement in (d): U(1) LIOMs with uniform quasi-local tails imply S_N is bounded uniformly in t and L. Its ingredients are standard. In the MBL literature it appears only as a heuristic ("since the conserved charges are local, S_N has to be bounded"), stated as such by Kiefer-Emmanouilidis et al. themselves [snippet]. I did not find it written as a theorem.
4. **The key quantitative point (this report).**
   - Any reasonable LIOM bound C on charge transfer across the cut is far above the observed S_N values.
   - So finite-L, finite-t S_N data cannot falsify the LIOM picture. Only the asymptotic *trend* (unbounded growth with t at L=∞, or of the saturation value with L) could.
   - Chávez et al. show that at L ≤ 28 that trend is not yet diagnostic.
5. **Broader context [memory].** Since 2020 the consensus has shifted: the XXZ chain at W ≈ 5–15 is in a "prethermal MBL" or resonance-dominated regime, and asymptotic MBL (if it exists) only sets in at W ≳ 18–20 (avalanche arguments). Slow S_N growth at W ≈ 5–10 therefore says little about the asymptotic phase either way.

---

## (b) Paper by paper

### Kiefer-Emmanouilidis, Unanyan, Fleischhauer, Sirker (KUFS) and later Sirker-group papers

- **PRL 124, 243601 (2020), arXiv:2003.04849, "Evidence for unbounded growth of the number entropy in MBL phases" (N).**
  - Claim: where S ~ ln t, S_N ~ ln ln t, "indicating continuing particle transport at a very slow rate". Fits S_N ≈ (ν/2) ln ln t for D ∈ [5,30]. "No indications for a saturation" [snippet].
  - Contains the heuristic "since the conserved charges η_j are local, the number entropy S_N has to be bounded" [snippet].
  - Also reports a "close to linear scaling" of S_N(t→∞) with system size [snippet]; whether this is in L or in ln L was not checked.
  - Model: XXZ / t-V chain with random fields. Method: ED. Sizes up to L≈24; the companion comment uses L = 8,…,16, 20, 24 [snippet]. Times: not checked.
  - Status: the data are not disputed; the interpretation is.
- **PRB 103, 024203 (2021), "Slow delocalization of particles in many-body localized phases", arXiv:2010.00565 (arXiv title "Absence of true localization in many-body localized phases") (N + H).**
  - Introduces the Hartley number entropy S_H: the α→0 Rényi limit, counting n with p(n) > p_c.
  - Finds a slow, steady increase in the number of particles crossing the cut in the interacting case. In the Anderson case it saturates for any cutoff.
  - The median S_N grows like the mean, so the growth is "typical, not rare-region" behavior.
  - States (H) that in a localized phase the truncated S_H "for any cutoff p_c>0 has to saturate in the thermodynamic limit".
  - Says a regime where S_N saturates while S keeps growing is never observed [snippet].
  - Reply to Luitz–Bar Lev: S_N is not sensitive enough to exclude slow transport, hence the switch to S_H and to the tails of p(n) [snippet].
- **Ann. Phys. 435, 168481 (2021), arXiv:2012.12436, "Unlimited growth of particle fluctuations in MBL phases" (N + P for free fermions).**
  - Tests alternative readings and argues that none holds. The alternatives are: (i) a transient, with saturation only at later times; (ii) fluctuations "very slow in building up but strictly limited in space", which is exactly the scenario of (d).
  - Derives strict bounds for non-interacting fermions, plus "renormalized" bounds for interacting systems (N/H).
  - States that unlimited growth of the Rényi S_N implies growth of ΔN [snippet].
- **SciPost Phys. 12, 034 (2022), arXiv:2108.03142, "Particle fluctuations and the failure of simple effective models for MBL phases" (N).**
  - Effective models that keep only interactions diagonal in the Anderson basis do *not* reproduce the microscopic S_N growth. Assisted and pair hopping cannot be neglected [snippet].
  - Concludes it is "questionable whether the microscopic model possesses an exponential number of exactly conserved local charges" [snippet].
  - Note: such effective models are a special case of (d), with U a Gaussian U(1) unitary with exponential tails. So their S_N saturation is *rigorously* expected.
- **SciPost Phys. 8, 083 (2020), arXiv:2003.03112, "Bounds on the entanglement entropy by the number entropy in non-interacting fermionic systems" (P, for Gaussian U(1) states).**
  - ln S^(2) is bounded above and below by the Rényi-2 number entropy S_N^(2) [snippet]. Hence S_N ~ ln S for free fermions: growth of S always comes with logarithmically slower growth of S_N, and vice versa.
  - Numerics cover clean, Anderson and bond-disordered (critical) chains [snippet].
  - This is the origin of the "S_N ~ ln S" relation. It is proven only for free fermions; for interacting systems it is an observation (N).
- **Comment arXiv:2203.06689 on Ghosh–Žnidarič, by KUFS (N).**
  - Argues: Ghosh–Žnidarič's initial states have only local dynamics and are unsuitable; the fitted power-law exponent depends on L and is not 1; ln ln t fits hold across sizes; the saturation value grows ~ ln L; the von Neumann entropy behaves similarly, so the dynamics look single-particle-like [snippet].
  - Journal status not found.
- **Zhao, Feng, Hu, Guo, Sirker, PRB 102, 195132 (2020), arXiv:2010.06678, "Entanglement dynamics in the 3D Anderson model" (N).**
  - S_N bounds S for 3D Anderson.
  - At strong disorder there is *no* window with S ~ ln t and S_N ~ ln ln t, in contrast with the Heisenberg chain [snippet].
- **Weisse, Gerstner, Sirker, Phys. Rev. Research (2025), arXiv:2401.08031, "Operator growth in disordered spin chains…" (N).**
  - Finds no convergence of the LIOM-generating unitary U in the interacting case. Suggests U does not exist, so no MBL [snippet].
  - Directly relevant: it attacks the *assumption* of (d).
- **arXiv:2606.07774 (2026), "Length-resolved operator growth and path-entropy obstructions to MBL" (H/N).** Authors not checked; it is cited by Sharifan et al. as the Sirker group's mechanism [snippet].
- **Sharifan, Kiefer-Emmanouilidis, Sirker, arXiv:2610.25557 (Oct 2026), "Signatures of slow particle transport in the interacting Aubry–André chain" (N).**
  - S ~ ln t and S_N ~ ln ln t, with rates following a power of the disorder strength, as for random disorder.
  - The S_H saturation value grows ~ ln L even where MBL was argued to be stable.
  - The quasiperiodic and random cases look alike, which they read as a common origin not based on avalanches.
  - Another excerpt says that in the ℓ-bit models S_N and S_H saturate with only weak L-dependence ("ℓ-bit-like"), and stresses finite-size scaling.
  - Figures at L=16 with 10^4 samples [snippet; excerpts partly inconsistent, check].
- **Not found this session:** dedicated Sirker-group papers on number entropy in Stark MBL or in 2D. There is a DAMOP 2023 talk "The number entropy as an accessible measure to detect transport and topological phases" [snippet, abstract only].

### Luitz & Bar Lev, PRB 102, 100202(R) (2020), arXiv:2007.13767 (arXiv title "Is there slow particle transport in the MBL phase?") (N)
- Method: saturation (long-time) values of S and S_N from *random product states* deep in the MBL phase.
- Results:
  - The *distributions* of saturated S_N converge with L.
  - They have a sharp cutoff at the value corresponding to one particle fluctuating across the cut (≈ the maximum entropy of a two-outcome fluctuation).
  - "Fluctuations of the particle number are strictly limited at strong disorder and preclude an indefinite growth of S_N."
  - Slow entanglement growth is "uniquely configurational" [snippet].
- Sizes not checked. Status: answered by KUFS through Hartley entropy and the tails of p(n). The two analyses probe different parts of p(n): its bulk versus tails of tiny probability.

### Ghosh & Žnidarič, PRB 105, 144203 (2022), arXiv:2112.12987 (v1 title "Theory of growth of number entropy in disordered systems") (N + H)
- Disordered isotropic Heisenberg chain: confirms small S_N growth but argues it is compatible with localization.
- A two-level resonance model (H) predicts:
  - a trivial initial rise;
  - an intermediate power-law approach to saturation;
  - mean (median) saturation S_N ∝ 1/W (1/W²);
  - strong dependence on the initial state.
- Saturation shows no significant increase with L, and at large W "if anything decreases with L" [snippet].
- Reply arXiv:2205.05768 [snippet]:
  - the comment's data do not refute them;
  - W=5 may be ergodic or near-critical;
  - the comment over-focuses on the Néel state;
  - the power-law approach holds over a longer span than ln ln t;
  - the saturation value does not grow with L.

### Aceituno Chávez, Artiaco, Klein Kvorning, Herviou, Bardarson, PRL 133, 126502 (2024), arXiv:2312.13420 (v1 Dec 2023, v2 20 Nov 2024) (N; model localized by construction)
All items [snippet]. The full text could not be opened, so verify the formulas.
- **Model**
  - H = U H′ U†, with H′ diagonal in the occupation (τ^z) basis: random on-site fields h_i plus interactions between ℓ-bits that decay exponentially with distance.
  - U is a **finite-depth random circuit of two-site gates, each built to preserve U(1)**. The paper states that the unitary is quasilocal and U(1)-preserving.
  - Gates are built from random Hermitian matrices M_i with weights w_i and a **"mixing factor" f**. The circuit is in Eqs. (8)–(10), which I could not see.
  - f sets the ℓ-bit localization length ξ_τ indirectly.
  - Example in the main text: depth d_u = 2 (U acts twice on each bond). Large runs use **d_u = 16, "which exceeds the longest half-chain considered"**.
  - ℓ-bit supports are checked numerically: averaged over 10^4 realizations, they fit a stretched exponential |Ō(x)| = C e^{−(x/ξ_τ)^β} with β ≈ 1, i.e. on-average exponential tails.
  - "The circuit structure accounts for the first exponential."
  - So: U(1)-conserving, yes; quasi-local with exponential tails, yes on average, and strictly finite depth.
- **Computed**
  - Disorder-averaged S_N(t) (Shannon) and the Hartley-type S_H (α = 10^−3, cutoff p_c = 10^−8) after a quench from a Néel state (one of the two at random).
  - Bipartition: half-chain (implied by "longest half-chain").
  - L = 12 to 28, with 10^4 realizations (another excerpt says ≥10^5).
  - Time evolution is exact in the ℓ-bit basis, so arbitrarily long times are reachable.
- **Results**
  - Ultraslow growth consistent with S̄_N ~ ln ln t "within a time interval that is exponential in system size".
  - Saturation at a **subextensive value that increases with L**.
  - One excerpt (likely v2): "Following an initial rise when L ≤ 16, the disorder-averaged saturation value then decreases slowly until L=28, suggesting that the number entropy could eventually saturate with system size."
  - Another excerpt: "monotonic and subextensive growth of saturation values with system size up to L=28".
  - The two excerpts may refer to different versions or to S_N versus S_H. **Check this.**
  - The approach to saturation is fitted as y = a − κ ln ln t (supplement).
  - No closed formula for the long-time value was seen.
- **Interpretation**
  - The growth comes from "the two exponential tails in the ℓ-bit model": the support tails of U and the interaction tails of H′, "both of which contain particle fluctuations".
  - Since the model is localized by construction, such growth is "inherent to MBL systems, at least for the system sizes we could simulate". Observing it in microscopic models is not sufficient evidence for delocalization.
  - They note that an earlier study (KUFS SciPost 2022) concluded that ℓ-bit models fail to explain the growth.

### Other theory and numerics
- **Sierant & Zakrzewski, PRB 105, 224203 (2022), arXiv:2109.13608 (N).**
  - Mainly about imbalance decay.
  - Includes S_N in TDVP: an ln ln t fit over t ∈ [20, 5000]. The 2026 review credits it with confirming ln ln t at L=50, W=10.
  - Quasiperiodic case: S_N ~ ln t, a plateau near t≈300, then renewed growth [snippet].
- **Sierant, Lewenstein, Scardicchio, Vidmar, Zakrzewski, Rep. Prog. Phys. 88, 026502 (2025), arXiv:2403.07111 (review).**
  - Expected S_n ∝ (1/z) ln t for dynamical exponent z ≥ 2 (citing KUFS); saturation in MBL; observed ln ln t up to L=28.
  - Restates the free-fermion bounds [snippet].
- **Review chapter "Many-body localization", arXiv:2604.12464 (2026; authors not checked).** Summarizes the debate:
  - KUFS read the data as slow transport;
  - others see a transient from single-particle motion across the cut;
  - Ghosh–Žnidarič tie it to rare resonances, with saturation at long times [snippet].
- **Faulend, Buljan, Štrkalj, arXiv:2603.28721 (2026), "Quantum beats in many-body localized systems" (H + N).**
  - Analytic model: two nearly resonant single-particle hoppings interact and modulate the Rabi amplitude ("many-body quantum beats").
  - This explains S_N growth and quasiparticle spreading in quasiperiodic potentials without thermalization.
  - Proposes engineered weak-disorder pairs for cold-atom tests [snippet].
- **Chen, Chen, Wang, PRB (DOI 10.1103/88sq-cm27), arXiv:2303.14825 (N).**
  - Symmetry- and energy-resolved dynamics in the disordered Bose–Hubbard chain.
  - Channel-resolved analysis suggests the ln ln t growth of S_N "may not signal the breakdown of MBL" in spin/fermion chains [snippet].
- **Chen, Chen, Wang, arXiv:2401.11339 (N).** Uses S = S_N + S_conf to separate two transition points, a thermalization breakdown and a localization transition, with Harris-compatible ν ≈ 2 (random) and ≈1.5 (quasiperiodic) [snippet].
- **Zeng, Hamma et al., arXiv:2305.08334, "Logarithmic light cone, slow entanglement growth, and quantum memory" (P under the log-light-cone assumption).** S(t) ≤ c ln t plus a sub-leading ln ln t term. They note that the ln ln t term "correlates with the number entropy" in XXZ numerics [snippet].
- **Toniolo & Bose, arXiv:2408.02016.** A log t law for α-Rényi entropies (α ≈ 1) from a log light cone [snippet].
- **Evers, Modak, Bera, PRB 108, 134204 (2023), arXiv:2302.11384.** Uses the entropy as an "internal clock". Argues the XXZ chain is not strictly MBL in reasonable W windows. S_N is not central [snippet].
- **Barghathi, Casiano-Diaz, Del Maestro, PRA 100, 022324 (2019), arXiv:1905.03312.** Operationally accessible entanglement (= S_conf) of spinless fermions [snippet id; journal from memory].

### Experiments (E)
- **Lukin et al., Science 364, 256 (2019), arXiv:1805.09819.**
  - Interacting Aubry–André Bose–Hubbard chain with site-resolved full counting statistics.
  - S = S_N + S_conf: number entanglement saturates quickly and is area-law; configurational entanglement grows slowly (logarithmically) toward volume law and is absent without interactions [snippet].
  - Chain length and maximum time (≈10²–10³ tunneling times?) [memory, not checked].
  - Times are far too short to separate ln ln t from saturation.
- **Gong et al., Phys. Rev. Research 3, 033043 (2021), arXiv:2012.11521.** 12 transmons. Half-chain S_N(t), autocorrelation and Hamming distance all change smoothly; the transition is found through sample-to-sample sensitivity [snippet].
- **Rispoli et al., Nature 573, 385 (2019), arXiv:1812.06959.** Subsystem number statistics across the MBL crossover [memory for content].

### Free fermions, full counting statistics, symmetry resolution
- **Klich & Levitov, PRL 102, 100502 (2009), arXiv:0804.1377** (P, Gaussian). Song, Flindt, Rachel, Klich, Le Hur, PRB 83, 161408(R) (2011), arXiv:1008.5191. Song, Rachel, Le Hur, PRB 82, 012405 (2010), id?. Song et al., PRB 85, 035409 (2012), id?.
  - S = Σ_k α_2k κ_2k (even cumulants of N_A), with α_2 = π²/3, α_4 = π⁴/45, α_6 = 2π⁶/945 [memory]. The series has convergence caveats.
- **Calabrese, Mintchev, Vicari, EPL 98, 20003 (2012), arXiv:1111.4836** (P). Exact Rényi–cumulant relations in any dimension. With impurities, all cumulants contribute [snippet; journal from memory].
- **Goldstein & Sela, PRL 120, 200602 (2018), arXiv:1711.09418. Xavier, Alcaraz, Sierra, PRB 98, 041106(R) (2018), arXiv:1804.06357. Bonsignori, Ruggiero, Calabrese, J. Phys. A 52, 475302 (2019), arXiv:1907.02084.** Symmetry-resolved entanglement and equipartition (CFT and free fermions) [snippet].
- **Turkeshi, Ruggiero, Alba, Calabrese, "Entanglement equipartition in critical random spin chains", PRB 102, 014455 (2020), arXiv:2005.03331.** Random-singlet phase: equipartition holds, broken at order s²/ln ℓ [snippet; journal from memory].
- **Abdul-Rahman, Nachtergaele, Sims, Stolz, Lett. Math. Phys. 106, 649 (2016), arXiv:1510.00262** (P). Disordered XY chain under the eigenfunction-correlator localization assumption: bipartite entanglement from a broad class of product states is bounded uniformly in t and L. Transport bounds for XX [snippet]. By S_N ≤ S, S_N is bounded too.
- **Hamza, Sims, Stolz, CMP 315, 215 (2012)** (P). Zero-velocity Lieb–Robinson bound for the disordered XY chain [memory].

---

## (c) Exact logic: S_N versus transport, LIOMs and dephasing

### Rigorous facts
- **F1 (P).** If ψ has definite total N, then ρ_A = ⊕_n p_n ρ_{A,n} and S = S_N + S_conf, with S_conf = Σ p_n S(ρ_{A,n}) ≥ 0. Hence S_N ≤ S. (Wiseman–Vaccaro PRL 91, 097902 (2003) [memory]; Lukin 2019; Goldstein–Sela 2018.)
- **F2 (P).** For an integer-valued variable with variance σ²: H ≤ ½ ln[2πe(σ² + 1/12)] (Massey 1988; Cover & Thomas, discrete-entropy bound [memory]). Rényi α ≥ 1 ≤ Shannon, so these are bounded too.
  - Rényi α < 1 and Hartley are **not** controlled by σ², because tails of tiny probability dominate them.
  - Consequence: S_N unbounded implies Var N_A unbounded. The converse is false, e.g. p = (1−ε)δ_0 + εδ_M.
- **F3 (P, Gaussian fermionic states only).**
  - S ≥ 4 ln2 · Var N_A [this report: from h(ν) ≥ 4 ln2·ν(1−ν) applied to the eigenvalues of the correlation matrix].
  - The cumulant series as above.
  - The KUFS bounds: ln S^(2) is sandwiched by S_N^(2).
  - So S_N ~ ln S is a free-fermion theorem, not a general law.
- **F4 (P, Anderson).**
  - For a Slater product state, Var N_A(t) = Σ_{x∈A, y∈B} |C_xy(t)|².
  - Dynamical localization (E sup_t |⟨x|e^{−iht}|z⟩| ≤ C e^{−|x−z|/ξ}) gives E sup_t Var N_A(t) ≤ C′, uniform in L [this report, elementary; also implied by ANSS 2016 together with F1].
- **F5 (P, conditional).** See (d): U(1) LIOMs with uniform quasi-local tails imply sup_{t,L} S_N ≤ ½ ln(2πe(4C² + 1/12)).
- **F6 (H, this report).** If only a log light cone holds (spreading length ℓ(t) = v ln t), then at most O(ℓ(t)) particles can cross the cut. That gives Var ≲ (v ln t)² and S_N ≲ ln ln t + O(1).
  - So **ln ln t is the maximal S_N growth compatible with a log light cone**. The observed S_N ≈ (ν/2) ln ln t with ν ≤ 2 sits inside that envelope; compare Zeng–Hamma.
  - ln ln t growth is therefore compatible with LIOM-type slow dynamics in general. It is incompatible only with the stronger structure [D, n_x] = 0 plus a U(1) dressing U with uniform tails.

### (i) S_N versus particle transport
- S_N(t) is the Shannon entropy of the full counting statistics of the transferred charge Q(t) = N_A(t) − N_A in the state ψ. It measures *fluctuations* of transferred charge, not net current.
- Scale: diffusion gives Var Q ~ √t, so S_N ≈ ¼ ln t (Gaussian approximation). Subdiffusion with dynamical exponent z gives S_N ≈ (1/2z) ln t. The review quotes (1/z) ln t; prefactor conventions differ.
- So S_N ~ ln ln t means Var ~ (ln t)^ν: zero diffusion constant, effectively z = ∞. Calling it "transport" means logarithmic charge spreading, not a nonzero conductivity.
- Valid:
  - "S_N → ∞ (t→∞ at L=∞, or the t=∞ value as L→∞)" implies "Var Q unbounded" (F2).
- Invalid:
  - S_N growth at finite (L, t) implies anything asymptotic.
  - Var growth implies S_N growth.
  - S_N growth implies σ_dc > 0.
  - Bounded S_N implies localization. Bounded S_N is necessary, not sufficient: it does not even bound Var, and fragmentation-type systems can also bound it.

### (ii) S_N versus the LIOM picture
- Valid (contrapositive of F5): unbounded S_N excludes H = U D U† with U U(1)-symmetric and with tails quasi-local uniformly in L. For random tails, it excludes the case E[C²] < ∞; average with Jensen, S̄_N ≤ ½ ln(2πe(4E[C²]+1/12)).
- What the data actually constrain is only C ≳ √(e^{2S_N}/(2πe) − 1/12)/2.
  - Example: S_N = 0.5 gives 2C ≳ 0.28.
  - This is trivially satisfied, so current data cannot falsify LIOMs.
- KUFS SciPost 2022 ("simple effective models fail") is F5 with U Gaussian: there, saturation is a theorem, so its failure to reproduce growth is expected.
- Chávez et al. show that a *generic* U(1) dressing still produces visible growth at L ≤ 28. That is consistent with F5 because F5 bounds only the asymptotic sup.
- The S_N ~ ln S relation (F3) must fail asymptotically in a true U(1)-LIOM phase: S ∝ ln t from dephasing while S_N stays bounded. Its persistence at accessible times is not evidence either way.
- Caveat on the premise: LIOMs for U(1) chains are unproven, and Weisse–Gerstner–Sirker 2025 argue numerically that U does not converge.

### (iii) S_N versus dephasing
- [D, N_A] = 0, so pure dephasing (U = 1) gives S_N ≡ 0 for product states (P). All S_N comes from the dressing U.
- Dephasing by distant ℓ-bits, at times ~ e^{r/ξ_J}, only reshuffles the *interference* between the two bounded dressing terms Δ′ and UΔ(t)U† defined in (d). It can generate slow time dependence of S_N, but its total amplitude is ≤ 4C² (P for the bound; the mechanism is H).
- This is the precise content of Chávez's "two tails" remark:
  - only the U-tails set the *size* of the fluctuations;
  - the interaction tails set the *times* (an ultraslow approach that looks like ln ln t over limited windows).
- The inference "S_N grows on dephasing time scales, so there is transport beyond the ℓ-bit picture" is invalid.
- S_conf ∝ ln t from dephasing is the standard LIOM prediction (Žnidarič–Prosen–Prelovšek 2008; Bardarson–Pollmann–Moore PRL 109, 017202 (2012), arXiv:1202.5532; Serbyn–Papić–Abanin 2013 [memory]). It is consistent with Lukin 2019 (E).

---

## (d) The proposed statement: is it known, and how does it square with Chávez et al.?

**Statement:** H = U D U†, with D diagonal in the occupation basis and U U(1)-symmetric and quasi-local. Then ‖N_A(t) − U N_A U†‖ ≤ C(ξ)|∂A|; Var N_A(t) is bounded; and S_N(t) ≤ ½ ln(2πe(Var + 1/12)), uniformly in t and L.

**Proof (this report; P, elementary).**
- Let Δ := U† N_A U − N_A and Δ′ := U N_A U† − N_A.
- Since [D, N_A] = 0, N_A(t) = U N_A U† + U Δ(t) U†, with Δ(t) = e^{iDt} Δ e^{−iDt}. Hence ‖N_A(t) − U N_A U†‖ = ‖Δ‖, exactly as stated.
- Then N_A(t) − N_A = Δ′ + UΔ(t)U†.
- For ψ with N_A ψ = n₀ ψ (any occupation product state):
  - Var_ψ N_A(t) ≤ ‖(N_A(t) − N_A)ψ‖² ≤ (‖Δ‖ + ‖Δ′‖)².
  - **Correction to the statement:** the bound is (‖Δ‖+‖Δ′‖)² ≤ 4C²|∂A|², not C².
- Bounding ‖Δ′‖:
  - If U = T exp(−i∫G(s)ds) with G = Σ_X g_X, [g_X, N] = 0 and ‖g_X‖ ≤ J e^{−diam X/ξ}: [g_X, N_A] = 0 unless X straddles ∂A, because [g_X, N_X] = 0. So ‖Δ′‖ ≤ ∫‖[G(s), N_A]‖ ds ≤ Σ_{X straddling} 2‖g_X‖ min(|X∩A|, |X∩B|) ≤ C(J,ξ)|∂A|.
  - For a brickwork of U(1) two-site gates: ‖Δ′‖ ≤ number of gates straddling the cut. Moreover p_t(n) is supported in n₀ ± (straddling gates in U and U†). So **all** Rényi S_N, Hartley included, are bounded by ln(2K+1).
  - With only exponential tails (no strict light cone), F2 bounds Shannon and α ≥ 1. Truncated Hartley needs tail estimates. Heuristically, moving m particles costs ~ e^{−m²/2ξ}, so S_H ≲ ½ ln(ξ ln(1/p_c)) + O(1): bounded, but depending on the cutoff (H).
- Then F2 gives S_N ≤ ½ ln(2πe(4C²|∂A|² + 1/12)).
- Implicit assumption: U commutes with N, i.e. LIOMs are dressed number operators. This is natural, but it is an assumption.

**Is it in the literature?**
- As a **heuristic**, yes:
  - KUFS PRL 2020: "since the conserved charges η_j are local, S_N has to be bounded".
  - KUFS PRB 2021: truncated S_H "has to saturate in the thermodynamic limit" in a localized phase.
  - Luitz–Bar Lev 2020: ℓ-bits predict no particle transport.
  - KUFS Ann. Phys. 2021: the "strictly limited in space" scenario, which they reject numerically.
  - All [snippet].
- The **key lemma** (charge moved by a U(1)-symmetric locality-preserving unitary, U*Q_A U − Q_A, is a bounded operator localized at ∂A) is standard in mathematical physics of charge transport:
  - Bachmann, Bols, De Roeck, Fraas, "A many-body index for quantum charge transport", Commun. Math. Phys. 375, 1249 (2020), arXiv:1810.07351 [snippet for paper; the exact lemma from memory].
  - Kapustin & Sopenko, "Local Noether theorem…", J. Math. Phys. 63, 091903 (2022), arXiv:2201.01327 [snippet].
  - The discrete-entropy step is textbook.
- I did **not** find the combined statement (LIOM, then uniform Var bound, then S_N ≤ ½ ln(2πe(·+1/12))) stated as a theorem in any MBL paper. The search was snippet-limited, so this is "not found", not "absent".

**How it squares with Chávez et al.**
1. Their U is exactly in the theorem's class: a U(1), finite-depth circuit with d_u = 16 fixed.
   - Deterministically, ‖N_A(t) − N_A‖ ≤ 2d_u = 32. That gives S_N ≤ ½ ln(2πe·1024) ≈ 4.9 and S_H ≤ ln(4d_u+1) ≈ 4.2.
   - Both exceed the trivial finite-size maximum ln(L/2+1) ≤ ln 15 ≈ 2.7 at L ≤ 28. The strict light cone bites only for L ≳ 4d_u ≈ 64, so it gives no constraint at their sizes.
   - The binding constraint is the f-controlled exponential decay, with an L-independent C(ξ_τ). Its value relative to the observed S_N is unknown to me, since f and ξ_τ were not extracted.
2. Hence "saturation value grows with L" can only be a **finite-size statement**: the theorem implies lim_{L→∞} sup_t S_N < ∞ for their model.
   - This matches their own (v2) excerpt: rise up to L≈16, slow decrease up to L=28, "could eventually saturate with system size".
   - The excerpt claiming monotonic growth up to L=28 needs checking against the figure.
3. For S_H with p_c = 10^−8, the cutoff heuristic above already exceeds ln 15 for modest ξ_τ. So L-growth of S_H at L ≤ 28 is expected even in this strictly localized model. This also undercuts using S_H(L) at L ≤ 28 (Sharifan et al. 2026) as evidence of delocalization.
4. Net:
   - The theorem and Chávez et al. are consistent.
   - The theorem says the ℓ-bit growth must stop; Chávez says it has not visibly stopped at L ≤ 28.
   - Together they make the point sharper: only **asymptotic unboundedness** of S_N (or of Var N_A) could falsify U(1)-LIOM MBL, and present numerics cannot reach that regime.

---

## (e) References
Tags: [s] = [snippet], [m] = [memory]. "id?" means not checked.

**Core debate**
- Kiefer-Emmanouilidis, Unanyan, Fleischhauer, Sirker. PRL 124, 243601 (2020). arXiv:2003.04849 [s]
- Same authors. PRB 103, 024203 (2021). arXiv:2010.00565 [s]
- Same authors. Ann. Phys. 435, 168481 (2021). arXiv:2012.12436 [s]
- Same authors. SciPost Phys. 12, 034 (2022). arXiv:2108.03142 [s]
- Kiefer-Emmanouilidis, Unanyan, Sirker, Fleischhauer. SciPost Phys. 8, 083 (2020). arXiv:2003.03112 [s]
- KUFS. Comment on Ghosh–Žnidarič. arXiv:2203.06689 [s]
- Ghosh, Žnidarič. Reply. arXiv:2205.05768 [s]
- Luitz, Bar Lev. PRB 102, 100202(R) (2020). arXiv:2007.13767 [s]
- Ghosh, Žnidarič. PRB 105, 144203 (2022). arXiv:2112.12987 [s]
- Aceituno Chávez, Artiaco, Klein Kvorning, Herviou, Bardarson. PRL 133, 126502 (2024). arXiv:2312.13420 [s]
- Sharifan, Kiefer-Emmanouilidis, Sirker. arXiv:2610.25557 (2026) [s]
- Faulend, Buljan, Štrkalj. arXiv:2603.28721 (2026) [s]
- "Length-resolved operator growth and path-entropy obstructions to MBL". arXiv:2606.07774 (authors not checked) [s]
- Weisse, Gerstner, Sirker. Phys. Rev. Research (2025). arXiv:2401.08031 [s]
- Zhao, Feng, Hu, Guo, Sirker. PRB 102, 195132 (2020). arXiv:2010.06678 [s]
- Chen, Chen, Wang. PRB (DOI 10.1103/88sq-cm27). arXiv:2303.14825 [s]
- Chen, Chen, Wang. arXiv:2401.11339 [s]
- Sierant, Zakrzewski. PRB 105, 224203 (2022). arXiv:2109.13608 [s]
- Sierant, Lewenstein, Scardicchio, Vidmar, Zakrzewski. Rep. Prog. Phys. 88, 026502 (2025). arXiv:2403.07111 [s; volume from memory]
- MBL review chapter. arXiv:2604.12464 [s]
- Evers, Modak, Bera. PRB 108, 134204 (2023). arXiv:2302.11384 [s]
- Zeng, Hamma et al. arXiv:2305.08334 [s]
- Toniolo, Bose. arXiv:2408.02016 [s]

**Experiments**
- Lukin et al. Science 364, 256 (2019). arXiv:1805.09819 [s]
- Gong et al. Phys. Rev. Research 3, 033043 (2021). arXiv:2012.11521 [s]
- Rispoli et al. Nature 573, 385 (2019). arXiv:1812.06959 [m]

**Free fermions and symmetry resolution**
- Klich, Levitov. PRL 102, 100502 (2009). arXiv:0804.1377 [s]
- Song, Flindt, Rachel, Klich, Le Hur. PRB 83, 161408(R) (2011). arXiv:1008.5191 [s]
- Song, Rachel, Le Hur. PRB 82, 012405 (2010). id? [s/m]
- Song et al. PRB 85, 035409 (2012). id? [m]
- Calabrese, Mintchev, Vicari. EPL 98, 20003 (2012). arXiv:1111.4836 [s/m]
- Goldstein, Sela. PRL 120, 200602 (2018). arXiv:1711.09418 [s]
- Xavier, Alcaraz, Sierra. PRB 98, 041106(R) (2018). arXiv:1804.06357 [s]
- Bonsignori, Ruggiero, Calabrese. J. Phys. A 52, 475302 (2019). arXiv:1907.02084 [s]
- Turkeshi, Ruggiero, Alba, Calabrese. PRB 102, 014455 (2020). arXiv:2005.03331 [s/m]
- Barghathi, Casiano-Diaz, Del Maestro. PRA 100, 022324 (2019). arXiv:1905.03312 [s/m]
- Wiseman, Vaccaro. PRL 91, 097902 (2003). id? [m]

**Rigorous and LIOM background**
- Abdul-Rahman, Nachtergaele, Sims, Stolz. LMP 106, 649 (2016). arXiv:1510.00262 [s]
- Hamza, Sims, Stolz. CMP 315, 215 (2012). arXiv:1108.3811 [m]
- Imbrie. J. Stat. Phys. 163, 998 (2016). arXiv:1403.7837 [m]
- Kim, Chandran, Abanin. arXiv:1412.3073 [s]
- Friesdorf et al. PRL 114, 170505 (2015). arXiv:1409.1252 [m]
- Serbyn, Papić, Abanin. PRL 111, 127201 (2013). arXiv:1305.5554 [m]
- Huse, Nandkishore, Oganesyan. PRB 90, 174202 (2014). arXiv:1305.4915 [m]
- Bardarson, Pollmann, Moore. PRL 109, 017202 (2012). arXiv:1202.5532 [s]
- Žnidarič, Prosen, Prelovšek. PRB 77, 064426 (2008). id? [m]
- Bachmann, Bols, De Roeck, Fraas. CMP 375, 1249 (2020). arXiv:1810.07351 [s]
- Kapustin, Sopenko. JMP 63, 091903 (2022). arXiv:2201.01327 [s]
- J. L. Massey. "On the entropy of integer-valued random variables", 1988 Beijing Int. Workshop on Information Theory [m]
- Cover & Thomas. *Elements of Information Theory* [m]

**MBL-stability context**
- De Roeck, Huveneers. PRB 95, 155129 (2017). arXiv:1608.01815 [m]
- Morningstar, Colmenarez, Khemani, Luitz, Huse. PRB 105, 174205 (2022). arXiv:2107.05642 [m]
- Long, Crowley, Khemani, Chandran. PRL 131, 106301 (2023). arXiv:2207.05761 [m]
- Sels, Polkovnikov. PRE 104, 054105 (2021). arXiv:2009.04501 [m]
- Šuntajs, Bonča, Prosen, Vidmar. PRE 102, 062144 (2020). arXiv:1905.06345 [m]
