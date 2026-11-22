# Ergodicity breaking, localization and slow dynamics in lattice gauge theories: what is rigorous (literature map, up to Nov 2026)

## 0. Method and how far to trust each tag
- No full paper was read this session. arXiv.org, ar5iv and export.arxiv were blocked by this session's network policy, and there was no `web-fetch` subagent tool. Everything below comes from WebSearch results, i.e. abstract-level excerpts pulled by the search engine.
- Tags:
  - **[read-abs]**: I read the abstract or abstract-level excerpts this session. Bibliographic data was confirmed by a repository or journal listing.
  - **[memory]**: not verified this session.
  - **[derived]**: my own algebra, done this session.
- Evidence types:
  - **(P)** proven theorem
  - **(N)** numerics
  - **(H)** heuristic or phenomenological argument
  - **(E)** experiment
  - **(X)** exact algebraic identity or mapping. This is P-level, but it is not a dynamical theorem by itself.
- **Not done:** the EXACT Smith–Knolle–Kovrizhin–Moessner (SKKM) mapping could not be read from the paper. Section (c) gives a reconstruction that I checked algebraically. The form of the spin term is [memory] and should be verified against arXiv:1701.04748, Eqs. 1–3.
- **arXiv ids to correct in your source list:**
  - SKKM PRL 118 is arXiv:**1701.04748**, not 1703.04583.
  - The "Papaefstathiou–Knolle–Bañuls PRB 2020" paper is really Papaefstathiou–**Smith**–Knolle, PRB 102, 165132 (2020).

## (a) Bottom line
1. **No dynamical localization theorem exists for any LGT model as such.** I found no math-ph paper proving localization, ergodicity breaking or slow thermalization *for an LGT*. The rigorous LGT+fermion math-ph work I found is about equilibrium: the Lieb flux phase, and Goller–Porta 2025 on stability of the π-flux phase.
2. **The exactly solvable DFL models are rigorous by inheritance.** The 1D SKKM model maps exactly (X) onto free fermions with an i.i.d. Bernoulli(1/2) potential ±2h, averaged over all 2^N configurations, for every gauge-invariant (sector-diagonal) observable. Combined with existing theorems on the 1D Bernoulli Anderson model, this gives rigorous averaged dynamical localization of the fermions:
   - spectral localization: Carmona–Klein–Martinelli 1987;
   - exponential dynamical localization in expectation: Bucaj et al. 2019, Jitomirskaya–Zhu 2019, Ge–Zhao 2020.
   - I found no paper that states this corollary explicitly. It is an easy, open "free" theorem.
3. **2D free-fermion DFL is only partly covered by mathematics.**
   - Binary potential (Bernoulli on Z²): only near the band edge (Ding–Smart 2020), or at large disorder away from finitely many exceptional energies (L. Li, arXiv:2002.11580).
   - Binary random flux: no theorem. Random-flux theorems need regular flux distributions and work near band edges (Klopp–Nakamura–Nakano–Nomura; Erdős–Hasler).
   - Expected obstructions:
     - chiral (bipartite) E=0 physics for pure random-sign hopping;
     - a class-D thermal metal when the fermions are Majoranas (Kitaev honeycomb).
4. **Interacting DFL is not more robust than MBL in 1D; it *is* MBL.**
   - With nearest-neighbour density interactions, 1D SKKM is exactly the XXZ chain in a binary random field (X).
   - It inherits every open question about MBL, including avalanches. Binary disorder is a worst case: uniform domains appear with finite density 2^-ℓ, and the exact resonances are not covered by Imbrie's proof, which needs continuous distributions.
   - Tian et al. (arXiv:2607.26156, 2026; N/H) argue that binary disorder does not suppress resonances in the 1D model behind the Google experiment, whereas 4-level disorder still shows localization signatures. This is the most direct challenge to interacting DFL in that setting.
5. **The LGT-specific mechanisms beyond MBL are classical-combinatorial.**
   - Hilbert space fragmentation and real-space percolation of kinetically active regions: Karpov et al. 2021; Chakraborty et al. 2022/2023; Jeyaretnam et al. 2025; Ciavarella et al. 2025; Cataldi et al. 2026.
   - Fragmentation is exact (X/P) wherever the constraint is exact. Localization then follows rigorously from finite clusters (classical percolation input).
   - Recent reinterpretations: Schwinger-model DFL is better explained by approximate fragmentation than by MBL (Jeyaretnam et al.), and quantum-link DFL persists as a classical cellular-automaton effect (Sala–Giudici–Halimeh).
6. **Rigorous slow-dynamics statements exist only through prethermalization theorems.**
   - Abanin–De Roeck–Ho–Huveneers (ADHH): CMP 2017.
   - Yin–Lucas: PRL 2023.
   - Yin–Surace–Lucas: PRX 2025. Its false-vacuum lifetime bound is ≥ exp(c ε^-d), and it proves that short-range-entangled metastable states are scars of a nearby Hamiltonian.
   - Halimeh–Hauke-type gauge protection is rigorous only where ADHH applies (an integer-spectrum, "compliant" penalty). There the authors adapt ADHH in appendices, which are sketches rather than stand-alone theorems.
   - Other protection claims are finite-size quantum-Zeno bounds (which depend on system size N) or numerics.
7. **Experiments.**
   - DFL observed in 1D and 2D on a Google superconducting processor (Science 2026).
   - Statistical localization from fragmentation observed in a Rydberg U(1) LGT (Nat. Phys. 2026).
   - String dynamics and string breaking imaged in (2+1)D (Nature 2025).
   - All of these are finite-time, finite-size (E) results.

## (b) Paper by paper
### B1. Disorder-free localization (DFL): core papers
- **Smith, Knolle, Kovrizhin, Moessner, PRL 118, 266601 (2017).**
  - Model: 1D spinless fermions on sites, spin-1/2 on bonds, local Z2 gauge symmetry, translation invariant.
  - Mechanism: an extensive set of conserved Z2 charges, so that "the time evolution can be described by a non-interacting Hamiltonian with effective binary disorder".
  - Initial states: product states |S⟩⊗|ψ⟩ with an arbitrary fermion state ψ. The randomness is "generated dynamically" by averaging over charge sectors [read-abs].
  - Type: X + N (free-fermion numerics; system sizes not verified) [memory].
  - Exact mapping: see (c).
- **Smith et al., PRL 119, 176601 (2017).**
  - Adding density interactions gives an exact mapping to the XXZ chain in a random binary field, averaged over all charge configurations.
  - Claims: MBL with log entanglement growth in a disorder-free system, and a state with combined area/volume-law entanglement, i.e. a quantum disentangled liquid (QDL).
  - Type: X + N (ED/tDMRG-type numerics; sizes not verified) [read-abs].
- **Smith et al., PRB 97, 245137 (2018).**
  - Model: 2D fermions coupled to dynamical Z2 gauge fields.
  - "Simplest models effectively non-interacting". Binary disorder changes the strong-disorder limit and allows "a delocalization transition via quantum percolation" for d>1.
  - Interactions that make the charges dynamical give only transient localization, or quasi-MBL.
  - Type: X + N + H [read-abs]. My reconstruction of the free 2D structure is in (c) [derived].
- **Brenes, Dalmonte, Heyl, Scardicchio, PRL 120, 030601 (2018).**
  - Model: lattice Schwinger model, U(1).
  - Mechanism: superselection sectors (static background charges) act as disorder; claimed MBL-like dynamics and confinement–entanglement link.
  - Claimed sub-logarithmic (ln ln t) entanglement growth at strong coupling.
  - Type: N (ED) + H [read-abs].
- **Papaefstathiou, Smith, Knolle, PRB 102, 165132 (2020).**
  - Model: (1+1)D QED with modified gauge-field dynamics.
  - Staggered discretization: free-fermion-soluble model with emergent DFL. Wilson discretization: interacting.
  - Type: X + N [read-abs]. The distribution of the emergent disorder was not verified.
- **Karpov, Verdel, Huang, Schmitt, Heyl, PRL 126, 130401 (2021).**
  - Model: 2D U(1) quantum link model (pure gauge, genuinely interacting). Both localized and ergodic phases exist without disorder.
  - Method: bounds the transition via classical correlated percolation, which implies Hilbert space fragmentation. Numerics use a classical-spin variational network.
  - Type: P-conditional (finite clusters imply localization, exactly) + N (percolation threshold, dynamics) + H (claim of greater robustness than MBL in 2D) [read-abs].
- **Chakraborty, Heyl, Karpov, Moessner, PRB 106, L060308 (2022).**
  - Studies the quantum localization transition of the same 2D QLM: continuous transition; scaling collapse with ν=4/3, β=5/36 (2D percolation exponents).
  - Quantum interference matters on top of classical fragmentation.
  - Type: N [read-abs].
  - Follow-up, **PRL 131, 220402 (2023)**: spectral response dominated by discrete peaks from finite active clusters, and halted information spreading. Type: N + H [read-abs].
- **Hart, Gopalakrishnan, Castelnovo, PRL 126, 227202 (2021).**
  - Model: quasi-1D compass and plaquette-Ising ladders, mapped onto free fermions plus strictly localized spins (one conserved Z2 charge per rung).
  - Log entanglement growth obtained from free-fermion Loschmidt echoes.
  - Type: X + N [read-abs].
- **Osborne, McCulloch, Halimeh, arXiv:2301.07720 (journal not verified).**
  - Model: regularized 2+1D QED with dynamical matter (tensor networks).
  - DFL persists for superpositions of many sectors; the imbalance vanishes for a single homogeneous sector. DFL strength depends on matter propagation directionality, and DFL is expected to weaken with dimension.
  - Type: N [read-abs].
- **Cataldi, Calajó, Silvi, Montangero, Halimeh, PRL 136, 170401 (2026).**
  - Model: 1+1D SU(2) hardcore-gluon LGT with matter. Static SU(2) background charges label the superselection sectors.
  - Three regimes: ergodic; fragmented (nonthermal but delocalized); disorder-free MBL.
  - Type: N [read-abs].
- **Not located:** "Giudici et al." on DFL and "Hart–De Tomasi–Castelnovo, fragmentation to DFL". Closest matches:
  - Giudici, Surace, Ebot, Scardicchio, Dalmonte, PRR 2, 032034 (2020): *disordered* U(1) LGTs; confinement reduces resonances, so ergodicity breaks at weaker disorder; weak finite-size effects. N [read-abs].
  - Hart, De Tomasi, Castelnovo, arXiv:2005.03036: random quantum comb; compact localized states turn into scars. N + H [read-abs].

### B2. Protecting or enhancing DFL, and variants
- **Halimeh, Zhao, Hauke, Knolle, "Stabilizing disorder-free localization", arXiv:2111.02427 (journal not verified).**
  - Terms linear in the gauge generator suppress gauge-breaking errors.
  - Protection holds "for times at least polynomial in the protection strength".
  - Type: N + H (perturbative) [read-abs].
- **Halimeh, Homeier, Zhao, Bohrdt, Grusdt, Hauke, Knolle, PRX Quantum 3, 020345 (2022).**
  - Through the quantum Zeno effect, a renormalized gauge theory with an enhanced local symmetry emerges dynamically and enhances DFL (Z2 pseudogenerator; Rydberg proposal).
  - Type: N + H [read-abs].
- **Lang, Hauke, Knolle, Grusdt, Halimeh, PRB 106, 174305 (2022).**
  - "Stark gauge protection": local (pseudo)generators weighted by a Stark potential; works for U(1) and Z2.
  - Type: N + H [read-abs].
- **Halimeh, Hauke, Knolle, Grusdt, "Temperature-induced DFL", arXiv:2206.11273.**
  - Thermal ensembles spanning many sectors also produce DFL.
  - Type: N/H [read-abs title/authors; content memory].
- **Sala, Giudici, Halimeh, PRB 109, L060305 (2024) (journal from citations), arXiv:2302.00681.**
  - Cellular-automaton circuits reproduce quantum-link-model (QLM) DFL quench dynamics and match exact quantum results at small sizes.
  - DFL persists in the thermodynamic limit as a *classical* effect of the finite-S regularization of the QLM. Quantum interference is not necessary.
  - Type: N + H [read-abs].
- **Jeyaretnam, Bhore, Osborne, Halimeh, Papić, Commun. Phys. 8, 172 (2025).**
  - Strong-coupling Schwinger DFL comes from *approximate* Hilbert space fragmentation, which causes jumps in entanglement within each sector.
  - The sector-averaged growth is ln t or a weak power law, not ln ln t. There is a single regime that mimics MBL at accessible sizes but is distinct from it.
  - Type: N (ED, MPS) + H (degenerate perturbation theory) [read-abs].
- **Further variants and probes:**
  - Smith et al., arXiv:1812.07981: logarithmic OTOC spreading without MBL in the free DFL model. Journal: PRL 123, 086602 (2019) [memory].
  - Yarloo, Mohseni-Rajaee, Langari, arXiv:1810.08434: "statistical bubble localization" in a Z2 LGT with a cluster spin chain. N [read-abs].
  - Najafzadeh, Langari, PRB 113, 144308 (2026), arXiv:2602.19586: "reversed QDL" in a disorder-free ladder. N [read-abs].
  - Yang, Jiang, Pan, arXiv:2509.06642: dissipation *enhances* DFL memory. N [read-abs].
  - Bhakuni, Cataldi et al., arXiv:2607.28730: quantum resources (magic) in U(1)/SU(2) DFL. N [read-abs].
  - Monitored Z2 LGT: Chakrabarti, Ara, Nirbhan, Bhattacharyya, Raychowdhury, arXiv:2603.28877 (no-click projective measurements, tensor networks up to 256 sites; framed as entanglement and measurement-induced transitions, not DFL); Ara, arXiv:2603.29900. N [read-abs].

### B3. DFL experiments
- **Gyawali, …, Roushan, Knolle, Halimeh, Moessner et al. (Google Quantum AI), Science 393, 71 (2026), arXiv:2410.06557.**
  - Platform: superconducting qubits; Z2 LGT with matter qubits on vertices and gauge qubits on links; Gauss law built in; second-order Trotter circuit.
  - Initial states are tunable superpositions over all disorder (gauge-sector) configurations. Perturbations fail to diffuse in 1D (search excerpt: a 38-qubit ring) and in 2D.
  - The paper proposes a sampling algorithm with polynomial speedup.
  - Rényi entropies show that superposition-prepared states differ from direct disorder sampling, which is the which-sector entropy point made in (c).
  - Type: E [read-abs].
- **Halimeh, Khodaeva, Kovrizhin, Moessner, Knolle, arXiv:2410.08268.**
  - Gate-based implementation of the SKKM model using an efficient three-qubit-gate encoding. Proposed benchmark: absent correlation spreading plus tunable growth of entanglement toward a volume law.
  - Type: N/proposal [read-abs].
- **Datla, Zhao, Ho, Klco, Loh, Nat. Phys. 22, 355 (2026), arXiv:2505.18143.**
  - First observation of statistical localization from nonlocal conserved quantities (strong fragmentation) in a Rydberg U(1) LGT; persists at infinite temperature.
  - Type: E [read-abs].

### B4. Critiques: is interacting DFL stable, and is it MBL in disguise?
- **Tian et al., arXiv:2607.26156 (2026).**
  - Flavor-extended Z2 LGT, mapped onto a mixed-field Ising chain with n-level bond disorder (the setting of the Google experiment).
  - Binary disorder, despite its larger variance, lacks the local amplitude diversity needed to suppress resonances; n=4 shows consistent localization signatures.
  - Localization is controlled by the local disorder spectrum and the resonant connectivity of the many-body Hilbert space, not by disorder strength alone.
  - Type: N + H [read-abs].
- **My synthesis (H).**
  - (i) Free DFL is Anderson localization, rigorous in 1D. It is fragile to any term that makes the charges dynamical. Then one gets heavy–light physics (QDL, quasi-MBL, asymptotic thermalization), as anticipated in Smith et al. 2018 and in the general literature on translation-invariant quasi-MBL [memory]:
    - De Roeck–Huveneers, CMP 332, 1017 (2014), which is rigorous: "asymptotic" localization from thermal disorder;
    - Yao et al., PRL 117, 240601 (2016);
    - Schiulaz–Silva–Müller, PRB 91, 184202 (2015).
  - (ii) Interacting 1D DFL is MBL with binary disorder, not "MBL in disguise". Avalanche arguments apply with finite-density clean domains (probability 2^-ℓ).
    - I found no paper doing an explicit avalanche analysis of DFL.
    - Note: at large h, binary disorder does *not* become "strong" in the usual sense. Hopping inside equal-q domains is unsuppressed, so the limit is a granular chain of finite clean segments with weak links.
  - (iii) Exact fragmentation is genuinely more robust than MBL, but only while gauge invariance and the constraint are exact. Karpov et al. argue it survives in 2D, where MBL is believed unstable (H).

### B5. Hilbert space fragmentation in LGTs
- **Yang, Liu, Gorshkov, Iadecola, PRL 124, 207602 (2020).**
  - Strict confinement: two U(1) conservation laws (domain-wall number, magnetization) give exponentially many Krylov sectors, some thermal and some Bethe-integrable.
  - Arises perturbatively in the confining limit of a Z2 LGT with fermions.
  - Type: X/P (exact fragmentation in the strict model) + H (perturbative LGT connection) [read-abs].
- **Ciavarella, Bauer, Halimeh, arXiv:2502.03533 (v2 Oct 2025; journal not verified).**
  - Generic fragmentation in Kogut–Susskind LGTs in d+1 dimensions: transitions with large electric-energy change are suppressed relative to plaquette and kinetic matrix elements, leaving an effective block-diagonal Hamiltonian.
  - Type: H (effective Hamiltonian) + N (L=4, electric-field truncation |n|<5; U(1) and SU(2)) [read-abs].
- **Stephen, Hart, Nandkishore, PRL 132, 040401 (2024).**
  - Fragmentation with exponentially many frozen states that is "provably" stable, to all orders in perturbation theory, against arbitrary k-body perturbations with k/L→0.
  - Not an LGT result per se, but the cleanest rigorous template.
  - Type: P (all-orders perturbative statement) [read-abs].

### B6. Quantum many-body scars in LGTs
- **Surace, Mazza, Giudici, Lerose, Gambassi, Dalmonte, PRX 10, 021041 (2020).**
  - Exact mapping of Rydberg-blockade (PXP) dynamics onto the U(1) QLM (spin-1/2 Schwinger); scar oscillations reinterpreted as string inversion.
  - Type: X + N [read-abs].
- **Lin, Motrunich, PRL 122, 173401 (2019).** Exact finite-bond-dimension MPS eigenstates at E=0 in PXP, i.e. the S=1/2 QLM. Type: P (exact eigenstates) [read-abs].
- **Banerjee, Sen, PRL 126, 220601 (2021); Biswas, Banerjee, Sen, SciPost Phys. 12, 148 (2022).**
  - U(1) QLM on ladders (and quantum dimer models): exponentially many mid-spectrum zero modes protected by an index theorem; special combinations survive as scars at finite λ.
  - Type: P (index-theorem zero modes) + N (scars at finite λ) [read-abs].
- **Desaules, Banerjee, Hudomal, Papić, Sen, Halimeh, PRB 107, L201105 (2023).** Scars persist for spin-S QLMs with S>1/2 (resonant and detuned regimes). Type: N [read-abs].
- **Desaules et al., PRB 107, 205112 (2023).** Scars in the truncated Schwinger model; finite-S scaling suggests they survive as S→∞. Type: N [read-abs].
- **Aramthottil et al., PRB 106, L041101 (2022).** Gauged Kitaev chain (Z2): scars exist in the deconfined phase and are absent in the confined phase. Type: N + analytic limit [read-abs].
- **Halimeh, Barbiero, Hauke, Grusdt, Bohrdt, Quantum 7, 1004 (2023).** Scars made robust to gauge-breaking errors by linear gauge protection (quantum Zeno). Type: N + H [read-abs].
- **Experiments:**
  - Bernien et al., Nature 551, 579 (2017): Rydberg PXP [memory].
  - Su et al., PRR 5, 023010 (2023): tilted Bose-Hubbard emulating PXP; low-entropy scar subspace [read-abs].
  - Both are E.
- **James, Konik, Robinson, PRL 122, 130603 (2019).** Confinement in the 1D and 2D Ising model gives nonthermal eigenstates (ordered phase only). Type: N + H [read-abs].

### B7. Confinement-induced slow dynamics and false-vacuum decay
- **Kormos, Collura, Takács, Calabrese, Nat. Phys. 13, 246 (2017).** Ising chain with transverse and longitudinal fields: confinement (mesons) after a quench suppresses the light cone. Type: N + H [read-abs].
- **Lerose, Surace, Mazza, Perfetto, Collura, Gambassi, PRB 102, 041118(R) (2020).**
  - Quasilocalized dynamics from confinement without disorder.
  - Effective Hamiltonians with Stark localization of excitations and weak entanglement growth; slow breakdown related to the Schwinger effect.
  - Type: H + N [read-abs].
- **Birnkammer, Bastianello, Knap, Nat. Commun. 13 (2022), arXiv:2202.12908.**
  - Multistage thermalization: a prethermal Gibbs ensemble with conserved meson number, reached through rare-event scattering cascades; true equilibration via Schwinger-like meson-number violation.
  - Type: H + N [read-abs].
- **Yin, Lucas, PRL 131, 050402 (2023).**
  - Prethermalization is generic for gapped local systems under small perturbations, in any dimension, with time scale ≥ exp[(Δ/ε)^a]. Includes a proof that false-vacuum decay is non-perturbatively slow.
  - Type: P [read-abs].
- **Yin, Surace, Lucas, PRX 15, 011064 (2025).**
  - Definition: a metastable pure state is one that all sufficiently local operators either stabilize or raise in energy.
  - Theorem: such states have non-perturbatively long lifetimes. The false-vacuum lifetime is ≥ exp(c ε^-d) in certain d-dimensional models, matching QFT.
  - Short-range-entangled metastable states are scars of a perturbatively close Hamiltonian.
  - Assumptions: finite local dimension, locality; translation invariance can be relaxed.
  - Type: P [read-abs].

### B8. Gauge-invariance protection (Halimeh–Hauke and others)
- **Halimeh, Hauke, PRL 125, 030503 (2020).**
  - An energy penalty V suppresses gauge-breaking errors of strength λ. Z2 and U(1) numerics "as well as analytic proofs".
  - Type: N + perturbative estimates (P at the level of time-dependent perturbation theory) [read-abs].
- **Halimeh, Lang, Mildenberger, Jiang, Hauke, PRX Quantum 2, 040311 (2021).**
  - Single-body penalties; an emergent global symmetry protects the gauge subspace; exponentially long protection for suitable coefficient sequences.
  - Analytics via a quantum-Zeno "gauge protection theorem". Numerics up to t > 10^10/J.
  - Type: N + P-applied (see note) [read-abs].
- **Halimeh, Homeier, Schweizer, Aidelsburger, Hauke, Grusdt, arXiv:2108.02203 (PRR 2022 [memory]).** Adapts ADHH in an appendix: a "compliant" local-pseudogenerator protection stabilizes gauge invariance up to times exponential in V. Type: P-applied, sketch [read-abs].
- **Van Damme, Lang, Hauke, Halimeh, arXiv:2104.07040.** Thermodynamic limit (iMPS): a renormalized gauge theory up to τ_ren ∝ exp(V/V0)/V0, with an ADHH-based analysis. Type: N + P-applied [read-abs].
- **Halimeh, Hauke, arXiv:2004.07248 ("staircase prethermalization").** Violation timescale delayed exponentially in system size. Type: N [read-abs].
- **Verdict on rigour:**
  - Exponential-in-V bounds are rigorous exactly when the penalty Σ c_j G_j has commensurate (integer) spectrum, so that ADHH applies. ADHH: CMP 354, 809 (2017), time scale "almost exponential", ~exp[c(V/λ)/log³(V/λ)] [read-abs + memory for the log power].
  - The papers apply ADHH in appendices rather than state stand-alone theorems with constants.
  - Non-compliant (incommensurate or random) sequences rely on numerics.
  - Quantum-Zeno bounds are rigorous for finite systems, but their constants grow with system size [memory/H].
  - The N-dependence in "staircase" results is numerical.

### B9. Math-ph background actually needed
- **1D Bernoulli Anderson model:**
  - Carmona–Klein–Martinelli, CMP 108, 41 (1987): localization for any non-degenerate distribution with a finite moment, including Bernoulli. (P) [read-abs]
  - Bucaj et al., Trans. AMS (2019), arXiv:1706.06135: positivity plus large-deviation proof. (P) [title read; journal memory]
  - Jitomirskaya–Zhu, CMP 370, 311 (2019), arXiv:1803.10697: arbitrary, e.g. Bernoulli, disorder without multiscale analysis. (P) [read-abs]
  - Ge–Zhao, J. Spectral Theory (2020), arXiv id?: exponential dynamical localization *in expectation* for any bounded, non-degenerate single-site distribution. (P) [via a secondary source]
- **2D Bernoulli:**
  - Ding–Smart, Invent. Math. 219, 467 (2020), arXiv:1809.09041: near the spectral edge. (P) [read-abs]
  - L. Li, arXiv:2002.11580: large disorder, localization outside small neighbourhoods (measure O(V^-c)) of finitely many energies, namely Dirichlet eigenvalues of -Δ on finite sets ("cluster" resonances). (P) [read-abs]
- **Random magnetic flux on Z²:**
  - Klopp–Nakamura–Nakano–Nomura, Ann. Henri Poincaré 4, 795 (2003): localization near the bottom of the spectrum, with a fixed domino-tiling restriction. (P) [read-abs; authors memory]
  - Erdős–Hasler, arXiv:1012.5185 and arXiv:1101.2139: i.i.d. fluxes, band edges, Wegner estimate without monotonicity. (P) [read-abs]
  - Neither covers binary (0/π) flux.
- **Free fermions to spins:** Hamza–Sims–Stolz, CMP 315, 215 (2012): one-particle dynamical localization implies a zero-velocity Lieb–Robinson bound for the disordered XY chain. (P) [memory]
- **MBL:**
  - Imbrie, J. Stat. Phys. 163, 998 (2016): needs continuous disorder plus a level-attraction assumption. (P, conditional) [memory]
  - De Roeck–Huveneers–Meeus–Prośniak, Physica A 631 (2023), arXiv:2308.08462: quasi-localization up to non-perturbatively long times at strong disorder. (P) [read-abs; disorder assumptions not verified, likely a regular density]
- **LGT + fermions, equilibrium:**
  - Lieb, PRL 73, 2158 (1994): flux phase. (P) [memory]
  - Goller–Porta, arXiv:2501.10065: stability of the π-flux phase for a 2D Z2 LGT with fermions, via reflection positivity. (P) [read-abs]
- **Symmetry-class caveat (N/H):** Kitaev honeycomb at T>0 means free Majoranas in random Z2 flux, which gives a delocalized thermal metal (Self et al., arXiv:1807.07926 [read-abs]). Binary random gauge fields therefore do NOT generically localize in 2D.

## (c) Which models are EXACTLY non-interacting, and which are interacting
**1D SKKM: explicit mapping.**
- Model form [memory; verify against the paper]:
  - H = -J Σ_j σ^z_{j,j+1}(f†_j f_{j+1} + h.c.) - h Σ_j σ^x_{j-1,j} σ^x_{j,j+1}.
- Algebra [derived]:
  - Gauss generators q_j = σ^x_{j-1,j} (-1)^{n_j} σ^x_{j,j+1}, with q_j = ±1 and [H, q_j] = 0. The abstract-level text writes the charges as "τ^z_j (-1)^{n_j}" with τ^z_j = σ^x_{j-1,j}σ^x_{j,j+1}, which is consistent.
  - Then σ^x_{j-1,j}σ^x_{j,j+1} = q_j (1 - 2n_j).
  - Gauge-invariant dressed fermions on an open chain: c_j = (Π_{k<j} σ^z_{k,k+1}) f_j. They obey canonical anticommutation relations and commute with all q_k. Moreover c†_j c_{j+1} = f†_j σ^z_{j,j+1} f_{j+1} and c†_j c_j = n_j.
- Result:
  - H = ⊕_q [ -J Σ_j (c†_j c_{j+1} + h.c.) + 2h Σ_j q_j c†_j c_j - h Σ_j q_j ].
  - So in each sector one has free fermions in a binary potential ±2h.
- Initial state |Ψ0⟩ = |↑…↑⟩_z ⊗ |ψ⟩, with ψ arbitrary:
  - each sector has weight ||P_q Ψ0||² = 2^-N, i.e. the q_j are i.i.d. Bernoulli(1/2);
  - the string acts trivially on Ψ0;
  - so for every gauge-invariant O (all c-bilinears and densities): ⟨O(t)⟩ = 2^-N Σ_q ⟨ψ| e^{iH_q t} O e^{-iH_q t} |ψ⟩.
  - This is an exact quenched binary-disorder average. A different initial spin state gives a different (possibly biased or correlated) disorder law. A single-sector initial state gives one realization, and the homogeneous sector is clean and ballistic.
- What is NOT a disorder average:
  - Spin observables: σ^z_{j,j+1} flips q_j and q_{j+1}, so its correlators are Loschmidt echoes between H_q and H_q' with q' differing on two sites (the origin of slow, log-type spin and OTOC dynamics).
  - The entanglement of the full state, which contains which-sector entropy (the source of the QDL area/volume structure and of the Google superposition-versus-sampling difference).
- Check against the alternative spin term [derived]: with a single-bond field -h Σ σ^x_{j,j+1} instead, the model is the standard interacting 1D Z2 LGT. In a sector it becomes a spin chain whose frustration graph contains a claw, so it is not Jordan–Wigner-free. The paper's "non-interacting" claim and the XXZ mapping in PRL 119 both point to the σ^xσ^x form above.
- Adding V n_j n_{j+1} (PRL 119): the Jordan–Wigner transform of the c's gives an XXZ chain with Δ ∝ V in the field 2h q_j. Genuinely interacting: binary-disorder MBL.

**2D Z2 (Smith 2018).** My reconstruction [derived; the exact Hamiltonian is not verified]:
- With a star term -h Σ_i Π_{j∈+i} σ^x_ij, both q_i and the plaquette fluxes W_p = Π_{∂p} σ^z are conserved, so each sector is free fermions with a binary potential 2h q_i and binary π-flux W_p.
- A σ^z-polarized spin initial state fixes W_p = +1 and randomizes q, giving the 2D Bernoulli Anderson model. A σ^x-polarized state randomizes the fluxes, giving the random π-flux problem.
- Large-h heuristic [derived, H]:
  - At first order in J, equal-q site clusters at p = 1/2 on the square lattice do not percolate, since site p_c ≈ 0.593.
  - Second-order J²/h tunnelling through an opposite site couples next-nearest neighbours, whose site-percolation threshold is ≈ 0.407 < 1/2.
  - So the strong-h limit is a quantum-percolation problem, not trivially localized. This is consistent with the abstract's "quantum percolation" remark.

| Model | Exact reduction | Type |
|---|---|---|
| SKKM 1D (PRL 118) | 1D Anderson–Bernoulli (i.i.d., p = 1/2, ±2h), averaged over 2^N sectors | Free |
| SKKM 2D free variants (PRB 97) | 2D Anderson–Bernoulli and/or binary random π-flux, depending on the initial state | Free |
| Hart–Gopalakrishnan–Castelnovo compass/plaquette ladders | Free fermions + one static Z2 charge per rung | Free (binary disorder) |
| Papaefstathiou et al., staggered U(1) | Free fermions with emergent disorder (law not verified) | Free |
| Kitaev honeycomb (same structure, not called DFL) | Free Majoranas in random Z2 flux | Free; class D: thermal metal, not localized |
| SKKM + density interaction (PRL 119) | XXZ in binary random field | Interacting |
| Brenes et al. Schwinger / U(1) QLM | Background charges give an effective random potential plus long-range (confining) electric interactions | Interacting |
| Papaefstathiou et al., Wilson discretization | — | Interacting |
| Karpov / Chakraborty 2D U(1) QLM | Exact fragmentation into real-space clusters (classical correlated percolation); dynamics inside clusters | Interacting |
| Osborne et al. 2+1D QED with matter | — | Interacting |
| Cataldi et al. SU(2) | — | Interacting |
| Google (2+1)D Z2 with matter | In 1D, a mixed-field Ising chain with binary bond disorder (per Tian et al.) | Interacting |

## (d) Existing rigorous results, and gaps that look within reach
**Existing (P/X):**
1. Exact free-fermion mappings and exact sector decompositions (SKKM; compass ladder; Kitaev-type). Exact PXP ↔ U(1) QLM equivalence.
2. Exact fragmentation:
   - strict confinement;
   - index-theorem zero modes in QLM/QDM ladders;
   - exact scar eigenstates in PXP/QLM;
   - perturbation-robust fragmentation (Stephen–Hart–Nandkishore).
3. Prethermalization theorems:
   - ADHH;
   - Yin–Lucas;
   - Yin–Surace–Lucas (false vacuum; metastable states are scars);
   - ADHH-based gauge-protection arguments (appendix-level).
4. Math-ph theorems for the effective one-body problems:
   - 1D Bernoulli, complete including dynamical localization in expectation;
   - 2D Bernoulli, partial (edge; large disorder away from finitely many energies);
   - continuous random flux (band edges);
   - subcritical percolation, where finite clusters mean trivial localization given exponential cluster-size tails (Menshikov; Aizenman–Barsky) [memory].

**Gaps where a theorem looks within reach:**
- **G1 (easy).** Rigorous DFL for 1D SKKM, as a corollary of the exact mapping plus Ge–Zhao or Jitomirskaya–Zhu. Concretely:
  - E_q sup_t |⟨δ_x, e^{-iH_q t} δ_y⟩| ≤ C e^{-η|x-y|};
  - hence, for any ψ and any gauge-invariant local A, B: sup_t |⟨[A(t), B]⟩_Ψ0| ≤ C e^{-η d(A,B)} (Hamza–Sims–Stolz-type zero-velocity Lieb–Robinson bound);
  - and a positive time-averaged return probability (inverse participation ratio) at every site, uniformly in L.
  - Proving a nonzero density-wave *imbalance* for all h>0 needs a small extra argument; it is immediate only at large h/J.
  - No explicit statement of this exists in the DFL literature as far as I found.
- **G2.** Spin (Loschmidt) observables in 1D SKKM. One needs absence of Anderson orthogonality catastrophe for a local two-site perturbation in the localized regime, i.e. a nonzero long-time spin autocorrelation. Existing AOC results for random Schrödinger operators (Gebert, Küttler, Müller and coauthors) [memory] look adaptable to Bernoulli 1D.
- **G3.** 2D SKKM with a binary potential at large h. L. Li's theorem gives spectral localization outside small energy windows. Turning this into a quench statement (most spectral weight of a CDW initial state is localized) needs dynamical localization in expectation in those windows. This is plausible but harder.
  - Binary π-flux: open. One must exclude E=0 (chiral/Gade physics on bipartite lattices).
  - Majorana (class D) versions are expected NOT to localize.
- **G4.** Prethermal fragmentation in DFL at strong binary field.
  - Q = Σ_j q_j n_j has integer spectrum, so ADHH applies directly to H = 2h Q + J·(hopping) + V·(density–density). V commutes with Q, so it can sit in the commuting part D, at least for V ≲ J.
  - Result: quasi-conservation of Q up to times ~exp(c h/J), in any dimension, with interactions. This is a rigorous slow-dynamics theorem for interacting DFL, short of localization.
  - The same structure applies to strong-coupling electric-energy fragmentation (Ciavarella et al.; Jeyaretnam et al.) with a truncated electric field.
- **G5.** String-breaking lifetime in LGTs with heavy matter, via Yin–Surace–Lucas. A lattice string with 2m > σ·r_min seems locally metastable in their sense, which would give a lower bound exponential in the critical pair separation 2m/σ. This needs care with truncated or unbounded electric fields.
- **G6.** Karpov-type 2D QLM. Make the "classical correlated percolation bound" a theorem:
  - (i) the inactive set is dynamically invariant (combinatorial);
  - (ii) exponential cluster-size tails for the correlated percolation, e.g. by stochastic domination by independent percolation.
  - Together these give rigorous localization on the classically subcritical side. The quantum transition itself is out of reach.
- **G7.** Gauge protection. Write a stand-alone theorem: an ADHH-type bound for the compliant penalty Σ c_j G_j with dynamical matter, explicit constants and locality. Truncated U(1) is fine. Kogut–Susskind with unbounded electric field needs extra work.
- **Out of reach / contested:**
  - interacting 1D DFL (binary-disorder MBL; not even Imbrie's framework covers it);
  - the 2D quantum percolation transition;
  - Schwinger-model DFL beyond prethermal times;
  - non-Abelian DFL;
  - any infinite-time statement once the charges are dynamical.

## (e) References
"Read" means read at abstract level (section 0); otherwise [memory].

**DFL core and variants**
- A. Smith, J. Knolle, D. L. Kovrizhin, R. Moessner, "Disorder-free localization", PRL 118, 266601 (2017), arXiv:1701.04748 [read-abs]
- A. Smith, J. Knolle, R. Moessner, D. L. Kovrizhin, "Absence of ergodicity without quenched disorder: from quantum disentangled liquids to many-body localization", PRL 119, 176601 (2017), arXiv:1705.09143 [read-abs]
- A. Smith, J. Knolle, R. Moessner, D. L. Kovrizhin, "Dynamical localization in Z2 lattice gauge theories", PRB 97, 245137 (2018), arXiv:1803.06574 [read-abs]
- A. Smith et al., "Logarithmic spreading of out-of-time-ordered correlators without many-body localization", PRL 123, 086602 (2019) [journal memory], arXiv:1812.07981 [title read]
- M. Brenes, M. Dalmonte, M. Heyl, A. Scardicchio, "Many-body localization dynamics from gauge invariance", PRL 120, 030601 (2018), arXiv:1706.05878 [read-abs]
- I. Papaefstathiou, A. Smith, J. Knolle, "Disorder-free localization in a simple U(1) lattice gauge theory", PRB 102, 165132 (2020), arXiv:2003.12497 [read-abs]
- P. Karpov, R. Verdel, Y.-P. Huang, M. Schmitt, M. Heyl, "Disorder-free localization in an interacting 2D lattice gauge theory", PRL 126, 130401 (2021), arXiv:2003.04901 [read-abs]
- N. Chakraborty, M. Heyl, P. Karpov, R. Moessner, "Disorder-free localization transition in a two-dimensional lattice gauge theory", PRB 106, L060308 (2022), arXiv:2203.06198 [read-abs]
- N. Chakraborty, M. Heyl, P. Karpov, R. Moessner, "Spectral response of disorder-free localized lattice gauge theories", PRL 131, 220402 (2023), arXiv:2211.14328 [read-abs]
- O. Hart, S. Gopalakrishnan, C. Castelnovo, "Logarithmic entanglement growth from disorder-free localization in the two-leg compass ladder", PRL 126, 227202 (2021), arXiv:2009.00618 [read-abs]
- O. Hart, G. De Tomasi, C. Castelnovo, "From compact localized states to many-body scars in the random quantum comb", arXiv:2005.03036 (PRR 2020 [memory]) [read-abs]
- J. Osborne, I. P. McCulloch, J. C. Halimeh, "Disorder-free localization in 2+1D lattice gauge theories with dynamical matter", arXiv:2301.07720 (journal?) [read-abs]
- J. C. Halimeh, H. Zhao, P. Hauke, J. Knolle, "Stabilizing disorder-free localization", arXiv:2111.02427 (journal?) [read-abs]
- J. C. Halimeh, L. Homeier, H. Zhao, A. Bohrdt, F. Grusdt, P. Hauke, J. Knolle, "Enhancing disorder-free localization through dynamically emergent local symmetries", PRX Quantum 3, 020345 (2022), arXiv:2111.08715 [read-abs]
- H. Lang, P. Hauke, J. Knolle, F. Grusdt, J. C. Halimeh, "Disorder-free localization with Stark gauge protection", PRB 106, 174305 (2022), arXiv:2203.01338 [read-abs]
- J. C. Halimeh, P. Hauke, J. Knolle, F. Grusdt, "Temperature-induced disorder-free localization", arXiv:2206.11273 (journal?) [title/authors read]
- P. Sala, G. Giudici, J. C. Halimeh, "Disorder-free localization as a purely classical effect", PRB 109, L060305 (2024) [journal from citing papers], arXiv:2302.00681 [read-abs]
- J. Jeyaretnam, T. Bhore, J. J. Osborne, J. C. Halimeh, Z. Papić, "Hilbert space fragmentation at the origin of disorder-free localization in the lattice Schwinger model", Commun. Phys. 8, 172 (2025), arXiv:2409.08320 [read-abs]
- G. Cataldi, G. Calajó, P. Silvi, S. Montangero, J. C. Halimeh, "Disorder-free localization and fragmentation in a non-Abelian lattice gauge theory", PRL 136, 170401 (2026), arXiv:2505.04704 [read-abs]
- Tian et al., "Role of flavor degrees of freedom in quantum simulations of disorder-free localization", arXiv:2607.26156 (2026) [read-abs]
- D. S. Bhakuni, G. Cataldi et al., "Quantum resources in disorder-free localization dynamics of gauge theories", arXiv:2607.28730 (2026) [read-abs]
- H. Yarloo, M. Mohseni-Rajaee, A. Langari, "Emergent statistical bubble localization in a Z2 lattice gauge theory", arXiv:1810.08434 (journal?) [read-abs]
- H. Najafzadeh, A. Langari, "Crossover from quantum chaos to a reversed quantum disentangled liquid in a disorder-free spin ladder", PRB 113, 144308 (2026), arXiv:2602.19586 [read-abs]
- X. Yang, X.-P. Jiang, L. Pan, "Dissipation-enhanced localization in a disorder-free Z2 lattice gauge system", arXiv:2509.06642 [read-abs]
- G. Giudici, F. M. Surace, J. E. Ebot, A. Scardicchio, M. Dalmonte, "Breakdown of ergodicity in disordered U(1) lattice gauge theories", PRR 2, 032034 (2020), arXiv:1912.09403 [read-abs]
- I. Sau, D. Banerjee, A. Sen, "Fate of many-body localization in an Abelian lattice gauge theory", PRB 113, 024204 (2026), arXiv:2405.20379 [read-abs]
- T. Grover, M. P. A. Fisher, "Quantum disentangled liquids", J. Stat. Mech. (2014) P10010 [memory], arXiv:1307.2288

**Experiments**
- G. Gyawali et al. (Google Quantum AI), "Observation of disorder-free localization using a (2+1)D lattice gauge theory on a quantum processor", Science 393, 71 (2026), arXiv:2410.06557 [read-abs]
- J. C. Halimeh, U. E. Khodaeva, D. L. Kovrizhin, R. Moessner, J. Knolle, "Disorder-free localization for benchmarking quantum computers", arXiv:2410.08268 [read-abs]
- P. R. Datla, L. Zhao, W. W. Ho, N. Klco, H. Loh, "Statistical localization of U(1) lattice gauge theory in a Rydberg simulator", Nat. Phys. 22, 355 (2026), arXiv:2505.18143 [read-abs]
- T. A. Cochran et al., "Visualizing dynamics of charges and strings in (2+1)D lattice gauge theories", Nature 642, 315 (2025), arXiv:2409.17142 [read-abs]
- D. González-Cuadra et al., "Observation of string breaking on a (2+1)D Rydberg quantum simulator", Nature (2025) [memory], arXiv:2410.16558 [title seen]
- G.-X. Su et al., "Observation of many-body scarring in a Bose-Hubbard quantum simulator", PRR 5, 023010 (2023), arXiv:2201.00821 [read-abs]
- H. Bernien et al., "Probing many-body dynamics on a 51-atom quantum simulator", Nature 551, 579 (2017), arXiv:1707.04344 [memory]

**Scars, fragmentation, confinement**
- F. M. Surace, P. P. Mazza, G. Giudici, A. Lerose, A. Gambassi, M. Dalmonte, "Lattice gauge theories and string dynamics in Rydberg atom quantum simulators", PRX 10, 021041 (2020), arXiv:1902.09551 [read-abs]
- C.-J. Lin, O. I. Motrunich, "Exact quantum many-body scar states in the Rydberg-blockaded atom chain", PRL 122, 173401 (2019), arXiv:1810.00888 [read-abs]
- D. Banerjee, A. Sen, "Quantum scars from zero modes in an Abelian lattice gauge theory on ladders", PRL 126, 220601 (2021), arXiv:2012.08540 [read-abs]
- S. Biswas, D. Banerjee, A. Sen, "Scars from protected zero modes and beyond in U(1) quantum link and quantum dimer models", SciPost Phys. 12, 148 (2022), arXiv:2202.03451 [read-abs]
- J.-Y. Desaules, D. Banerjee, A. Hudomal, Z. Papić, A. Sen, J. C. Halimeh, "Weak ergodicity breaking in the Schwinger model", PRB 107, L201105 (2023), arXiv:2203.08830 [read-abs]
- J.-Y. Desaules, A. Hudomal, D. Banerjee, A. Sen, Z. Papić, J. C. Halimeh, "Prominent quantum many-body scars in a truncated Schwinger model", PRB 107, 205112 (2023), arXiv:2204.01745 [read-abs]
- A. S. Aramthottil et al., "Scar states in deconfined Z2 lattice gauge theories", PRB 106, L041101 (2022), arXiv:2201.10260 [read-abs]
- J. C. Halimeh, L. Barbiero, P. Hauke, F. Grusdt, A. Bohrdt, "Robust quantum many-body scars in lattice gauge theories", Quantum 7, 1004 (2023), arXiv:2203.08828 [read-abs]
- Z.-C. Yang, F. Liu, A. V. Gorshkov, T. Iadecola, "Hilbert-space fragmentation from strict confinement", PRL 124, 207602 (2020), arXiv:1912.04300 [read-abs]
- A. N. Ciavarella, C. W. Bauer, J. C. Halimeh, "Generic Hilbert space fragmentation in Kogut–Susskind lattice gauge theories", arXiv:2502.03533 (journal?) [read-abs]
- D. T. Stephen, O. Hart, R. M. Nandkishore, "Ergodicity breaking provably robust to arbitrary perturbations", PRL 132, 040401 (2024), arXiv:2209.03966 [read-abs]
- M. Kormos, M. Collura, G. Takács, P. Calabrese, "Real-time confinement following a quantum quench to a non-integrable model", Nat. Phys. 13, 246 (2017), arXiv:1604.03571 [read-abs]
- A. Lerose, F. M. Surace, P. P. Mazza, G. Perfetto, M. Collura, A. Gambassi, "Quasilocalized dynamics from confinement of quantum excitations", PRB 102, 041118(R) (2020), arXiv:1911.07877 [read-abs]
- A. J. A. James, R. M. Konik, N. J. Robinson, "Nonthermal states arising from confinement in one and two dimensions", PRL 122, 130603 (2019), arXiv:1804.09990 [read-abs]
- S. Birnkammer, A. Bastianello, M. Knap, "Prethermalization in one-dimensional quantum many-body systems with confinement", Nat. Commun. 13 (2022), arXiv:2202.12908 [read-abs]
- C. Yin, A. Lucas, "Prethermalization and the local robustness of gapped systems", PRL 131, 050402 (2023), arXiv:2209.11242 [read-abs]
- C. Yin, F. M. Surace, A. Lucas, "Theory of metastable states in many-body quantum systems", PRX 15, 011064 (2025), arXiv:2408.05261 [read-abs]

**Gauge protection**
- J. C. Halimeh, P. Hauke, "Reliability of lattice gauge theories", PRL 125, 030503 (2020), arXiv:2001.00024 [read-abs]
- J. C. Halimeh, H. Lang, J. Mildenberger, Z. Jiang, P. Hauke, "Gauge-symmetry protection using single-body terms", PRX Quantum 2, 040311 (2021), arXiv:2007.00668 [read-abs]
- J. C. Halimeh, L. Homeier, C. Schweizer, M. Aidelsburger, P. Hauke, F. Grusdt, "Stabilizing lattice gauge theories through simplified local pseudogenerators", arXiv:2108.02203 (PRR 4, 033120 (2022) [memory]) [read-abs]
- M. Van Damme, H. Lang, P. Hauke, J. C. Halimeh, "Reliability of lattice gauge theories in the thermodynamic limit", arXiv:2104.07040 (journal?) [read-abs]
- J. C. Halimeh, P. Hauke, "Staircase prethermalization and constrained dynamics in lattice gauge theories", arXiv:2004.07248 [read-abs]
- D. Abanin, W. De Roeck, W. W. Ho, F. Huveneers, "A rigorous theory of many-body prethermalization for periodically driven and closed quantum systems", CMP 354, 809 (2017), arXiv:1509.05386 [read-abs]
- J. C. Halimeh, Z. Davoudi et al., "Quantum simulation of out-of-equilibrium dynamics in gauge theories" (review), arXiv:2509.03586 [read-abs]

**Mathematics**
- R. Carmona, A. Klein, F. Martinelli, "Anderson localization for Bernoulli and other singular potentials", CMP 108, 41 (1987) [read-abs]
- V. Bucaj, D. Damanik, J. Fillman, V. Gerbuz, T. VandenBoom, F. Wang, Z. Zhang, "Localization for the one-dimensional Anderson model via positivity and large deviations for the Lyapunov exponent", Trans. AMS (2019) [memory], arXiv:1706.06135 [title read]
- S. Jitomirskaya, X. Zhu, "Large deviations of the Lyapunov exponent and localization for the 1D Anderson model", CMP 370, 311 (2019), arXiv:1803.10697 [read-abs]
- L. Ge, X. Zhao, "Exponential dynamical localization in expectation for the one dimensional Anderson model", J. Spectral Theory (2020) [memory], arXiv id? [result via a secondary source]
- J. Ding, C. K. Smart, "Localization near the edge for the Anderson Bernoulli model on the two dimensional lattice", Invent. Math. 219, 467 (2020), arXiv:1809.09041 [read-abs]
- L. Li, "Anderson-Bernoulli localization at large disorder on the 2D lattice", arXiv:2002.11580 (journal?) [read-abs]
- F. Klopp, S. Nakamura, F. Nakano, Y. Nomura, "Anderson localization for 2D discrete Schrödinger operators with random magnetic fields", Ann. Henri Poincaré 4, 795 (2003) [read-abs; authors memory]
- L. Erdős, D. Hasler, "Wegner estimate and Anderson localization for random magnetic fields", arXiv:1012.5185; "Anderson localization for random magnetic Laplacian on Z²", arXiv:1101.2139 (journals [memory]) [read-abs]
- E. Hamza, R. Sims, G. Stolz, "Dynamical localization in disordered quantum spin systems", CMP 315, 215 (2012), arXiv id? [memory]
- J. Z. Imbrie, "On many-body localization for quantum spin chains", J. Stat. Phys. 163, 998 (2016), arXiv:1403.7837 [memory]
- W. De Roeck, F. Huveneers, B. Meeus, O. Prośniak, "Rigorous and simple results on very slow thermalization, or quasi-localization, of the disordered quantum chain", Physica A 631 (2023), arXiv:2308.08462 [read-abs]
- W. De Roeck, F. Huveneers, "Asymptotic quantum many-body localization from thermal disorder", CMP 332, 1017 (2014), arXiv id? [title seen; content memory]
- E. H. Lieb, "Flux phase of the half-filled band", PRL 73, 2158 (1994) [memory]
- Goller, Porta, "Stability of the π-flux phase for Z2 lattice gauge theory coupled to fermionic matter", arXiv:2501.10065 [read-abs]
- C. N. Self et al., "Thermally induced metallic phase in a gapped quantum spin liquid: a Monte Carlo study of the Kitaev model with parity projection", arXiv:1807.07926 (PRB 2019 [memory]) [read-abs]

**Monitored LGT (adjacent)**
- N. Chakrabarti, N. Ara, N. Nirbhan, A. Bhattacharyya, I. Raychowdhury, "Effects of monitoring on entanglement dynamics for 1+1D Z2 lattice gauge theory", arXiv:2603.28877 [read-abs]
- N. Ara, "Dynamics of entanglement entropy for a locally monitored lattice gauge theory", arXiv:2603.29900 [read-abs]
