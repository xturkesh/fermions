# Monitored free (and weakly interacting) fermions: where a measurement-induced transition is predicted to EXIST (beyond 1D short-range)

Coverage: literature up to Nov 2026 as indexed by web search. Verification method and its limits: direct fetching of arxiv.org and export.arxiv.org was blocked in this session (host not on the allowlist), and the session's shared web-search budget (200 calls) ran out before everything was checked. So arXiv IDs and journal references were checked against arXiv abstract/HTML listings and publisher pages (APS, Quantum) **as they appeared in web-search results**, not by opening the pages directly. Tags used below:
- **[V]**: arXiv ID and journal reference both seen on an arXiv listing and/or a publisher page in the search results.
- **[V-arXiv]**: arXiv listing seen; no journal reference found, so the paper may be unpublished or the reference is unverified.
- **[V-sec]**: the detail was seen only in a secondary source (a citing paper's reference list, an institutional repository or an aggregator).
- **[T]**: only the title or ID was seen; the content comes from a search-engine summary of the arXiv HTML and must be checked.

Dates are arXiv first-submission months (read from the ID) plus the journal year. The evidence-strength gradings in the "Inferences" subsections are my own assessments.

Corrections to the assignment's premises, found during verification:
1. The field-theory paper is by Guo, Foster, Jian and **Ludwig**, not Lucas.
2. The free-fermion adaptive paper is Ravindranath–**Yang**–Chen (Quantum 2025). Ravindranath–**Han**–Yang–Chen (PRB 2023) covers generic adaptive random circuits with feedback and is not described as free-fermion in the sources seen.
3. The "Lévy flights" paper (Poboiko, Szyniszewski, Turner, Gornyi, Mirlin, Pal, PRL 2025) studies a **short-range** chain with misaligned pair measurements, not long-range hopping.

---

## 1. d >= 2: monitored free fermions above one dimension (nature of the metallic phase and of the transition; exponents)

### Takeaway
For number-conserving (U(1)) monitored free fermions, a genuine measurement-induced phase transition (MIPT) exists for d > 1. It separates a "metallic" (gapless) phase, where entanglement and the particle-number second cumulant scale as ℓ^{d−1} ln ℓ, from an area-law phase (∝ ℓ^{d−1}). The transition is an Anderson-type transition of an SU(R) replica nonlinear sigma model (NLSM) in d+1 dimensions with R→1 (class AIII in Chahine–Buchhold's language). The theory is controlled only near d = 1+ε. 2D numerics give ν = 1.4 ± 0.1, strong multifractality, emergent conformal-like scaling, and conductance-like mesoscopic fluctuations. For 2D Majorana systems (measurement-only Kitaev-type circuits), a Majorana-liquid phase with S ~ L log L and continuous transitions with quantum Lifshitz scaling are established.

### Cited Findings
**U(1) (complex) fermions**
- **Poboiko, Gornyi, Mirlin, "Measurement-induced phase transition for free fermions above one dimension"**: arXiv:2309.12405 (Sep 2023); Phys. Rev. Lett. 132, 110403 (2024) [V].
  - Mapping and control: the problem maps onto an SU(R) replica NLSM in d+1 dimensions with R→1. A one-loop RG is controlled for d = 1+ε (ε ≪ 1), so the analytic exponents are only parametrically justified near d = 1.
  - Phases: a gapless phase with ℓ^{d−1} ln ℓ scaling of both the particle-number second cumulant and the entanglement entropy, and an area-law phase with ℓ^{d−1} scaling.
  - 2D numerics: on a square-lattice model the critical point was located numerically, and a scaling collapse gives **ν = 1.4 ± 0.1**.
  - Sources: [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.110403); [arXiv](https://arxiv.org/abs/2309.12405); [arXiv HTML](https://arxiv.org/html/2309.12405)
- **Poboiko & Mirlin, "Quantum dynamics of monitored free fermions: Evolution of quantum correlations and scaling at measurement-induced phase transition"**: arXiv:2512.01772 (1 Dec 2025); Phys. Rev. B 113, 144311 (2026) [V-sec: the journal reference comes from the KIT repository record; the APS page was not seen].
  - It studies, analytically and numerically, the dynamics of many-body free fermions under local density measurements.
  - It states the general NLSM relation ν = [β′(G_c)]^{−1} for d > 1, where β(G_c) = 0.
  - Sources: [arXiv v1](https://arxiv.org/abs/2512.01772v1); [arXiv HTML v2](https://arxiv.org/html/2512.01772v2); [KIT record](https://publikationen.bibliothek.kit.edu/1000193037/180599078)
- **Chahine & Buchhold, "Entanglement phases, localization, and multifractality of monitored free fermions in two dimensions"**: arXiv:2309.12391 (Sep 2023); Phys. Rev. B 110, 054313 (2024) [V].
  - Model and method: continuously monitored U(1) free fermions in 2D, studied with exact numerical simulations of the entanglement entropy and the wave-function inverse participation ratio. The authors derive an exact fermion replica master equation and outline two analytic routes: an NLSM in **symmetry class AIII**, and a third-quantization / non-Hermitian Hubbard-model route for bipartite lattices.
  - Results: the multifractality matches the NLSM of the Anderson transition. The critical point also shows emergent conformal invariance, which is "not typically associated with Anderson localization". The purification time is linked to the multifractal exponent.
  - Sources: [APS](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.054313); [arXiv](https://arxiv.org/abs/2309.12391)
- **Poboiko, Gornyi, Mirlin, "Mesoscopic fluctuations and multifractality at and across measurement-induced phase transition(s)"**: arXiv:2507.11312 (15 Jul 2025); Phys. Rev. B 112, L220202 (2025) [V].
  - Model: 2D free fermions under projective monitoring of local charge; the paper studies fluctuations over the ensemble of quantum trajectories across the MIPT.
  - Result: there is an analogy to Anderson localization. A quantity G_AB corresponds to the two-terminal conductance and C(r) to the two-point conductance, "albeit with different replica limits and unconventional symmetry class, geometry, and boundary conditions".
  - Sources: [APS](https://journals.aps.org/prb/abstract/10.1103/2k3q-h6lz); [arXiv](https://web3.arxiv.org/abs/2507.11312)
- **Chahine & Buchhold, "Spectral transitions of the entanglement Hamiltonian in monitored free fermions"**: arXiv:2509.08584 (Sep 2025) [V-arXiv; the dimension studied is not confirmed].
  - Method: a random-matrix analysis of the entanglement Hamiltonian. The adjacent-gap ratio and a Kullback–Leibler divergence (KL₁) register a transition out of ergodic behaviour, and the spectral form factor and Thouless time support it.
  - A second divergence, KL₂, hints at a possible non-ergodic extended (multifractal) regime at intermediate monitoring. The abstract also mentions Fermi-liquid and Lifshitz-type features.
  - A secondary search summary adds: an ergodic to non-ergodic transition at a quantum Lifshitz fixed point, with estimates of the critical point and of ν (values not retrieved).
  - Sources: [arXiv HTML](https://arxiv.org/html/2509.08584v1); [arXiv](https://arxiv.org/abs/2509.08584)
- **Fava, Piroli, Bernard, Nahum, "Monitored fermions with conserved U(1) charge"** (v1 title: "A tractable model of…"): arXiv:2407.08045 (Jul 2024); Phys. Rev. Research 6, 043246 (Dec 2024) [V].
  - Method: a replica field theory with the flavor number as a large parameter, giving an SU(N) NLSM coupled to fluctuating charge hydrodynamics.
  - In d = 1: area law at all monitoring rates, with nontrivial scaling when the correlation length is large.
  - In higher d: the NLSM has a transition between a nontrivial phase and an area-law phase. The correlation length can be infinite there, and the semiclassical treatment can hold over arbitrarily large scales.
  - Open point flagged by the authors: whether the coupling of charge fluctuations to the entanglement sector is relevant at the higher-d critical point.
  - Sources: [APS](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.043246); [arXiv](https://arxiv.org/abs/2407.08045); [arXiv HTML](https://arxiv.org/html/2407.08045); [arXiv HTML v2](https://arxiv.org/html/2407.08045v2)
- **"Entanglement dynamics of monitored noninteracting fermions on graphics processing units"**: arXiv:2508.18468 (Aug 2025) [T; authors and venue not verified]. It reports a 2D MIPT at a protocol-dependent critical monitoring strength, while 1D Dirac fermions remain area-law. [arXiv HTML](https://arxiv.org/html/2508.18468)
- *One-line 1D context, which other researchers cover:* in 1D short-range U(1) fermions there is no transition; the area law holds at any rate, and the apparent transition is a sharp finite-size crossover. Poboiko, Pöpperl, Gornyi, Mirlin, arXiv:2304.03138; Phys. Rev. X 13, 041046 (2023) [V]. [APS](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.13.041046); [arXiv](https://arxiv.org/pdf/2304.03138)

**2D Majorana / Kitaev-honeycomb measurement-only circuits (free Majoranas in a Z₂ background)**
- **Lavasani, Luo, Vijay, "Monitored quantum dynamics and the Kitaev spin liquid"**: arXiv:2207.02877 (Jul 2022); Phys. Rev. B 108, 115135 (2023) [V].
  - Model: (2+1)D circuits of local projective measurements on the honeycomb lattice.
  - Phases: either a topologically ordered phase or a "critical" phase with a logarithmic violation of the area law and long-range tripartite entanglement.
  - Method: a Majorana-parton description, confirmed by extensive numerics.
  - The explicit "L log L" form is attributed to this phase by later papers (Klocke et al. 2024; arXiv:2509.16758).
  - Sources: [APS accepted](https://link.aps.org/accepted/10.1103/PhysRevB.108.115135); [arXiv](https://arxiv.org/pdf/2207.02877)
- **Sriram, Rakovszky, Khemani, Ippoliti, "Topology, criticality, and dynamically generated qubits in a stochastic measurement-only Kitaev model"**: arXiv:2207.07096 (14 Jul 2022; v2 16 Sep 2023); Phys. Rev. B 108, 094304 (2023) [V-sec: the journal reference comes from a citing PRR reference list].
  - Model: random measurements of two-qubit bond operators.
  - When one bond type dominates: area-law phases that protect two topological qubits on a torus for a time exponential in system size.
  - At comparable rates: a critical phase with a log violation of the area law, which protects the topological qubits only for a time polynomial in system size.
  - Purification has z = 1/2 ("super-ballistic").
  - Sources: [arXiv](https://arxiv.org/abs/2207.07096); [citing PRR](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.L042063)
- **Klocke, Simm, Zhu, Trebst, Buchhold, "Entanglement dynamics in monitored Kitaev circuits: loop models, symmetry classification, and quantum Lifshitz scaling"**: arXiv:2409.02171 (Sep 2024); Phys. Rev. B 111, 224301 (2025) [journal reference from a search summary of the arXiv/group pages; APS page not seen].
  - Model: measurement-only circuits in classes BDI and D (Kitaev and Yao–Kivelson), mapped onto Majorana loop models (loop fugacity n = 1 for monitored circuits).
  - Field theory: a CP^{n−1} or RP^{n−1} NLSM, depending on loop orientability.
  - Phases: a delocalized, highly entangled Majorana-liquid phase with S(A) ~ L log L.
  - Transitions: continuous, with quantum Lifshitz scaling but critical exponents from two distinct universality classes.
  - Sources: [arXiv](https://www.arxiv.org/abs/2409.02171); [Berkeley Kavli](https://kavli.berkeley.edu/node/2120)

### Inferences
- **Existence of a d ≥ 2 transition: strong evidence.** Three independent lines agree:
  - analytic NLSM work: Poboiko et al., controlled at d = 1+ε, and Fava et al., controlled at large flavor number;
  - large-scale 2D numerics from at least two groups (KIT and Cologne), plus a 2025 GPU study;
  - an Anderson-transition-like phenomenology (multifractality, conductance-like fluctuations).
- **Exponents: weak independent support.** The only quantitative ν found is 1.4 ± 0.1, from a single group. I found no independent numerical estimate. It is not established that the d = 1+ε expansion reaches d = 2 reliably.
- Poboiko et al.'s SU(R) principal-chiral-type NLSM and Chahine–Buchhold's "class AIII" NLSM appear to be the same chiral-class structure, since class AIII has a U(n)-type principal-chiral target. This is my inference; I did not check it in the texts.
- In NLSM language, d > 1 has an unstable fixed point G_c. Flows to large G give the ℓ^{d−1} ln ℓ metal; flows to small G give the area law. This is inferred from ν = 1/β′(G_c) and from the stated phase structure.
- The ℓ^{d−1} ln ℓ metallic scaling (U(1)) and the L log L Majorana liquid (2D Kitaev-type circuits) are the same "Fermi-surface-like" log-violation of the area law, in two different symmetry classes. Chahine–Buchhold's mention of "Fermi-liquid features" is consistent with this, but that is an interpretive link.

### Gaps
- Values of ν, η or the multifractal exponents from Chahine–Buchhold (2024, 2025) were not retrieved, nor the quantum-Lifshitz exponents of Klocke et al.
- No 3D numerics were found.
- Whether charge hydrodynamics modifies the U(1) universality class in d ≥ 2 is open, as Fava et al. themselves flag.
- The authors of arXiv:2508.18468 were not verified, and the dimension studied in arXiv:2509.08584 was not confirmed.
- The search budget ran out before I could look for 2026 follow-ups giving independent ν estimates.

---

## 2. Long-range hopping, long-range measurements, and Lévy-flight-like spreading (dependence on the power-law exponent)

### Takeaway
Long-range hopping ∝ r^{−α} is the clearest 1D route to non-area-law steady states for free fermions. Recurring thresholds are α ≈ 1 (no transition below it per Minato et al.; toward volume law per Fan & García-García) and α = 3/2 (above it the system behaves short-ranged). For 1/2 < α < 3/2 a "fractal" phase with S ∝ L^{3/2−α} and dynamical exponent z = (2α−1)/2 < 1 appears at large N. Groups disagree on whether an area-law transition exists in 1 < α < 3/2. The 2022 claim of a log phase for α > 3/2 is superseded, for finite-N complex fermions, by 2025–2026 work. A short-range model with misaligned pair measurements gives Lévy-flight spreading and S ∝ ℓ^{1/3}.

### Cited Findings
- **Minato, Sugimoto, Kuwahara, Saito, "Fate of measurement-induced phase transition in long-range interactions"**: arXiv:2104.09118 (Apr 2021); Phys. Rev. Lett. 128, 010603 (2022) [V].
  - Model: fermion chains with long-range hopping/interactions ∝ r^{−α}.
  - Result: for α > α_c a MIPT is observed, with an α-dependent universality class. For α < α_c no MIPT occurs even for arbitrarily strong measurement. **α_c = 1 for free fermions** (arXiv v1).
  - Method: γ_c(α) from finite-size scaling of the effective central charge.
  - **Conflict between versions:** arXiv v2 states sufficient conditions α > d + 1/2 for general bilinear systems and α > d + 1 for nonintegrable systems. The PubMed (published) abstract gives α > d/2 + 1 for bilinear systems. The two coincide (3/2) at d = 1; I could not resolve this.
  - Sources: [arXiv v2](https://arxiv.org/abs/2104.09118v2); [arXiv HTML v1](https://arxiv.org/html/2104.09118v1); [PubMed](https://pubmed.ncbi.nlm.nih.gov/35061459/)
- **Müller, Diehl, Buchhold, "Measurement-induced dark state phase transitions in long-ranged fermion systems"**: arXiv:2105.08076 (May 2021); Phys. Rev. Lett. 128, 010605 (5 Jan 2022) [V].
  - Model: free fermions with long-range hopping under continuous local density measurement.
  - New phase: an algebraic phase for **1 < p ≲ 3/2, independent of the measurement rate**, with algebraic entanglement growth and algebraically decaying density–density correlations, both with fractional exponents.
  - Phase boundaries: two critical lines separate this phase from a log-entanglement "critical phase" at small rates and an area-law phase at large rates.
  - Method: perturbative RG giving a modified sine-Gordon theory. The MIPT is framed as a transition in the dark state of an effective non-Hermitian Hamiltonian. Data are archived on [Zenodo](https://zenodo.org/records/10455005).
  - Sources: [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.010605); [arXiv](https://arxiv.org/abs/2105.08076v2)
  - **Partly superseded:** see the Müller–Buchhold–Diehl 2025 and Fan & García-García 2026 entries below.
- **Zhang, Liu, Jian, Chen, "Universal entanglement transitions of free fermions with long-range non-unitary dynamics"**: arXiv:2105.08895 (May 2021); Quantum 6, 723 (27 May 2022) [V].
  - Models: solvable Brownian models, namely a large-N SYK₂ chain and a single-flavor chain, which share the same phase diagram; the static SYK₂ chain gives the same structure, which the authors read as universality.
  - Phases for α > 1/2: a log phase for α > 3/2 with z = 1, and a fractal phase for 1/2 < α < 3/2 with **z = (2α−1)/2** and **S ∝ L_A^{1−z} = L_A^{3/2−α}**.
  - Purification cannot distinguish the two phases (the entropy decays as L/T). Higher-d phase diagrams are also discussed.
  - Sources: [Quantum](https://quantum-journal.org/papers/q-2022-05-27-723/); [arXiv](https://arxiv.org/abs/2105.08895); [arXiv HTML](https://arxiv.org/html/2105.08895)
- **Sahu, Jian, Bentsen, Swingle, "Entanglement phases in large-N hybrid Brownian circuits with long-range couplings"**: arXiv:2109.00013 (Sep 2021) [V-arXiv; journal reference not verified]. Long-range couplings are irrelevant for α > 3/2 in 1D, and non-interacting models show fractal-entangled phases. [arXiv](https://arxiv.org/pdf/2109.00013)
- **Russomanno, Piccitto, Rossini, "Entanglement transitions and quantum bifurcations under continuous long-range monitoring"**: arXiv:2307.05685 (Jul 2023) [V-arXiv; **journal reference not verified**; the Pisa repository record gives none].
  - Model: steady-state bipartite entanglement of the quantum trajectories of a free-fermion chain under continuous **nonlocal** monitoring.
  - Claim: qualitative changes ("bifurcations") of the trajectory distribution are connected to the entanglement-transition points.
  - Sources: [arXiv](https://arxiv.org/pdf/2307.05685); [Pisa repository](https://arpi.unipi.it/handle/11568/1215470)
- **Poboiko, Szyniszewski, Turner, Gornyi, Mirlin, Pal, "Measurement-induced Lévy flights of quantum information"**: arXiv:2501.12903 (Jan 2025; v3 Oct 2025); Phys. Rev. Lett. 135, 170403 (21 Oct 2025) [V].
  - Model: a short-range chain monitored by measurements on pairs of adjacent sites. A misalignment angle θ makes neighbouring measurements non-commuting.
  - Result: at maximal misalignment θ = π, **S ∝ ℓ^{1/3}**. The cause is a measurement-induced quasiparticle decay rate that vanishes at one point of the Brillouin zone, which gives superdiffusive, Lévy-flight-like spreading. The behaviour persists in the measurement-only limit.
  - Method: an effective NLSM with long-range couplings, checked with large-scale numerics.
  - Sources: [APS](https://journals.aps.org/prl/abstract/10.1103/tx71-1cd9); [arXiv](https://arxiv.org/abs/2501.12903v3); [KIT](https://www.tkm.kit.edu/english/publications/4301_4367.php)
- **Fan & García-García, "Entanglement transitions and multifractality in monitored free-fermions with random long-range hopping"**: arXiv:2607.23871 (Jul 2026; one aggregator says 8 Sep 2026, a date conflict) [V-arXiv; preprint].
  - Model: a 1D chain of monitored complex fermions with **random** power-law hopping.
  - Phases: for α ≲ 1 the entanglement grows faster than log, tending toward volume law at small α. **For α > 3/2 the area law holds at any monitoring strength.**
  - At the critical point: S ∝ log L, with multifractal density–density correlations.
  - Sources: [arXiv](https://arxiv.org/abs/2607.23871); [arXiv HTML](https://arxiv.org/html/2607.23871)
- **Müller, Buchhold, Diehl (2025)**, the same group as the 2022 PRL: for free (non-interacting) Dirac fermions, "no entanglement transition at nonzero measurement strength in the thermodynamic limit". arXiv:2502.02645 [V-arXiv]; details under Question 3. [arXiv](https://arxiv.org/abs/2502.02645)
- Long-range dissipators: Piccitto, Chiriacò, Rossini, Russomanno, "Entanglement behavior and localization properties in monitored fermion systems", arXiv:2502.09401 (Feb 2025), studies among other things a Kitaev chain with long-range dissipators [T]. [arXiv HTML](https://arxiv.org/html/2502.09401v2)
- Contrast for bosons: continuously monitored free bosons with long-range couplings show **no** MIPT under local measurements, because long-range-coupling entanglement outweighs the measurements. arXiv:2405.19768 [T]. [arXiv HTML](https://arxiv.org/html/2405.19768)

### Inferences
- **Superseded claim:** the log-entanglement phase at small rates for p (α) > 3/2 (Müller et al. 2022; and Zhang et al. 2022 at N = ∞) is, for finite-N complex fermions, superseded by:
  - the 1D short-range area-law result (PRX 2023);
  - the same group's 2025 statement that free Dirac fermions have no entanglement transition;
  - Fan & García-García 2026 (area law at any monitoring for α > 3/2, though with random hopping).
- **Unresolved conflict in the window 1 < α < 3/2:**
  - Müller et al. 2022 report an algebraic phase independent of the measurement rate, i.e. no area-law transition.
  - Minato et al. 2022 report a transition for every α > 1.
  - Fan & García-García 2026 report a transition with log scaling at criticality.
  - The models differ (continuous vs projective monitoring; clean vs random hopping), and none of the sources I saw reconciles them.
- **Unifying picture (my inference):**
  - z < 1 for α < 3/2 (Zhang et al.), Lévy-flight spreading from a vanishing decay rate (Poboiko et al. 2025) and z = 1/2 purification in the 2D Kitaev measurement-only model (Sriram et al.) all suggest the same thing: effectively long-range information spreading is what stabilizes sub-volume "fractal" entanglement against the measurement-induced localization that otherwise gives an area law in 1D.
  - The ℓ^{1/3} law of Poboiko et al. is an example where this arises without long-range hopping.
- **Evidence strength:**
  - Large-N results (Zhang et al.): exact at N = ∞.
  - Finite-N results: numerics plus perturbative RG, with group-to-group disagreements. I rate the evidence moderate.

### Gaps
- The fractional exponents of Müller et al. 2022 and the α-dependent exponents of Minato et al. were not retrieved.
- The behaviour for α < 1/2 in Zhang et al. (presumably volume law) is not confirmed in the retrieved text.
- I don't know whether the ℓ^{1/3} regime of the Lévy-flight paper is a point (θ = π) or a phase.
- The journal references for Russomanno et al. and Sahu et al. are unverified.
- The PRL vs arXiv-v2 discrepancy in Minato et al.'s sufficient condition is unresolved.
- arXiv:2609.04091 ("Effective Hamiltonian description on monitored Majorana chains: correlated power-law hoppings and unconventional entanglement scaling", 2026) was seen as a title only [T]: [arXiv](https://arxiv.org/pdf/2609.04091).

---

## 3. Interactions: do they restore a genuine transition, including in 1D, with volume law at weak monitoring?

### Takeaway
Yes, on the theory side. Independent field theories from 2024–2025 (KIT; Guo–Foster–Jian–Ludwig) predict that interactions allow measurement-induced transitions **in any d, including d = 1**. In theories that include chaos, a volume-law phase appears at weak monitoring: Jian et al. 2021 (large-N SYK₄), Guo et al. 2025 (charge sharpening inside the volume-law phase) and Foster et al. 2025 (volume law via a dangerously irrelevant mass). Direct numerical confirmation of interacting-fermion entanglement transitions is still limited; Müller–Buchhold–Diehl call the interacting entanglement entropy "an open problem".

### Cited Findings
- **Poboiko, Pöpperl, Gornyi, Mirlin, "Measurement-induced transitions for interacting fermions"**: arXiv:2410.07334 (9 Oct 2024; v2 10 Jan 2025) [V-arXiv; **journal reference not verified**].
  - Model: number-conserving interacting fermions with a time-independent Hamiltonian, under measurements. The paper studies how entanglement and charge fluctuations change when interactions are switched on.
  - Method: an NLSM derived for interacting fermions, analysed by RG.
  - Result: "in the presence of interaction, there are measurement transitions in any d, at variance with free fermions, for which a d = 1 system is always in the area-law phase."
  - Sources: [arXiv](https://arxiv.org/abs/2410.07334v1); [arXiv HTML](https://arxiv.org/html/2410.07334v1)
- **Guo, Foster, Jian, Ludwig, "Field theory of monitored, interacting fermion dynamics with charge conservation"**: arXiv:2410.07317 (9 Oct 2024; v2 14 Apr 2025); Phys. Rev. B 112, 064304 (12 Aug 2025) [V].
  - Method: a unified effective field theory for all phases and transitions.
  - Result: interactions reduce the replica symmetry of the non-interacting theory, which allows additional phases and transitions. The theory explains why the charge-sharpening transition occurs only inside the volume-law phase.
  - One search summary described this APS page as concerning "measurement-induced heating"; this seems to be a description of the same paper, not a different one.
  - Sources: [APS](https://journals.aps.org/prb/abstract/10.1103/gdxd-pw8v); [arXiv](https://arxiv.org/abs/2410.07317); [APS March Meeting 2025 abstract](https://meetings-archive.aps.org/smt/2025/mar-n29/10/)
- **Foster, Guo, Jian, Ludwig, "Free-fermion measurement-induced volume- to area-law entanglement transition in the presence of fermion interactions"**: arXiv:2510.23706 (Oct 2025) [V-arXiv; preprint].
  - Model: 1D monitored **interacting Majorana fermions in class DIII**, with no conserved quantities.
  - Conjecture: the volume-to-area-law MIPT is the **non-interacting DIII critical point**. The volume-law phase arises through a dangerously irrelevant "mass" whose vanishing at the transition suppresses chaos.
  - The authors identify a candidate non-interacting critical point analytically and propose numerical tests.
  - Source: [arXiv](https://arxiv.org/abs/2510.23706)
- **Müller, Buchhold, Diehl, "Monitored interacting Dirac fermions"**: arXiv:2502.02645 (4 Feb 2025; v2 21 Feb 2025) [V-arXiv; journal reference not verified].
  - Method: a bosonic replica field theory analysed by RG.
  - Result: for nonzero **attractive** interaction there is a **BKT-type** transition at γ_c. For γ < γ_c, density correlations decay algebraically; for γ > γ_c the system is localized, with exponentially decaying correlations. The BKT character shows in the scaling of the correlation length.
  - Free limit: γ_c → 0, so the algebraic phase shrinks to a point ("implicit double fine-tuning" of weak free-fermion measurements to the BKT endpoint). Along the non-interacting line there is no entanglement transition at γ > 0.
  - The interacting entanglement entropy is left open.
  - Source: [arXiv](https://arxiv.org/abs/2502.02645)
- **Large-N interacting reference point:** Jian, Liu, Chen, Swingle, Zhang (PRL 127, 140601 (2021) [V]; details under Question 5). Interactions lower the continuous O(2) replica symmetry of SYK₂ to a discrete C₄, giving **volume-law** entanglement in the symmetry-broken (weak-monitoring) phase. [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.127.140601)
- A review (arXiv:2203.01513) notes that the absence of interactions is essential for the O(2)×O(2) symmetry, which an on-site Brownian SYK₄ term breaks [T]. [arXiv](https://arxiv.org/pdf/2203.01513)
- **2D interacting Majoranas:** Zhu, Tantivasadakarn, Trebst, "Structured volume-law entanglement in an interacting, monitored Majorana spin liquid", arXiv:2303.17627 (Mar 2023); Phys. Rev. Research 6, L042063 (2024) [V]. Random measurement-only circuits combining the Kitaev bond and plaquette couplings give a volume-law phase with structured entanglement, unlike the area-law or critical phases of the free (bond-only) versions. [arXiv](https://arxiv.org/pdf/2303.17627); [APS](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.L042063)

### Inferences
- **Two independent groups now agree that interactions qualitatively change the 1D answer.** KIT's NLSM with interactions and Ludwig et al.'s reduced replica symmetry both make 1D transitions possible. The mechanism is that interactions break the large continuous replica symmetry responsible for the free-fermion Goldstone/NLSM physics; a concrete large-N example is O(2) → C₄ in SYK.
- **Evidence strength:**
  - Field theory: analytic but perturbative/RG, and mostly uncontrolled at strong coupling.
  - Large-N SYK: exact at N = ∞.
  - Numerics: strong for 2D interacting Majorana circuits (Zhu et al.), limited for 1D U(1) interacting fermions.
  - Foster et al. 2025 is explicitly a conjecture.
- The Müller et al. BKT transition concerns density correlations, not entanglement. Whether it coincides with an entanglement transition is open.

### Gaps
- Not verified this session (search budget exhausted): Xing, Turkeshi, Schirò, Fazio, Rossini, "Interactions and integrability in weakly monitored Hamiltonian systems" (I believe arXiv:2308.09133 / PRB 2024, unconfirmed), and Fuji–Ashida (arXiv:2004.11957, seen only as a title in results).
- The journal references for arXiv:2410.07334 and arXiv:2502.02645 are unverified.
- The exponents or numerics for the interacting volume-to-area transition in the KIT and Ludwig papers were not retrieved.

---

## 4. Disorder and quasiperiodicity: interplay of Anderson localization and monitoring

### Takeaway
In 1D, disorder and quasiperiodic studies (2022–2026) report a "critical/log" phase that survives up to a finite potential strength, with BKT-like criticality. These claims sit inside the contested 1D U(1) picture (asymptotically area law in the clean case), so their thermodynamic-limit status is uncertain. Monitoring suppresses Anderson localization. In 2D the MIPT itself is an Anderson-type transition with multifractality (Question 1). Random long-range hopping gives multifractal density correlations at the critical point.

### Cited Findings
- **Szyniszewski, Lunt, Pal, "Disordered monitored free fermions"**: arXiv:2211.02534 (4 Nov 2022; v3 21 Oct 2023); Phys. Rev. B 108, 165126 (17 Oct 2023; Editors' Suggestion) [V].
  - Model: 1D spinless fermions in a random potential under continuous monitoring.
  - Results: the critical phase is stable up to a finite disorder strength, and the criticality is consistent with **BKT** universality. Monitoring suppresses Anderson localization. The area-law phase at weak dissipation shows power-law decay of single-particle wavefunctions.
  - Sources: [APS](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.165126); [arXiv](https://arxiv.org/abs/2211.02534v1); [UCL PDF](https://discovery.ucl.ac.uk/id/eprint/10182449/1/PhysRevB.108.165126.pdf)
- **Yamamoto & Hamazaki (authors inferred; the search result did not list them directly), "Localization properties in disordered quantum many-body dynamics under continuous measurement"**: arXiv:2301.07290 (Jan 2023); Phys. Rev. B 107, L220201 (2023) [V-sec].
  - Method: the fidelity between random trajectories.
  - Results: disorder and measurement both suppress entanglement spreading but give different dynamics. Localization cannot be read off from the steady-state entanglement alone. The authors propose ultracold-atom tests.
  - Sources: [INSPIRE](https://inspirehep.net/arxiv/2301.07290); [arXiv](https://arxiv.org/pdf/2301.07290); [ICTP talk](https://indico.ictp.it/event/10202/session/22/contribution/99.pdf)
- **Matsubara, Yamamoto, Koga, "Measurement-induced phase transitions for free fermions in a quasiperiodic potential"**: arXiv:2503.23807 (31 Mar 2025) [V-arXiv; journal reference not verified].
  - Model: the Aubry–André–Harper chain, studied with quantum trajectories.
  - Results: robust log scaling up to a critical potential **V_c/J ≈ 2.3**. For weak potential, measurement drives a log to area-law transition; for stronger potential, any finite measurement gives the area law.
  - Sources: [arXiv](https://arxiv.org/abs/2503.23807); [arXiv HTML](https://arxiv.org/html/2503.23807v2)
- **arXiv:2509.09538 (Sep 2025), "Entanglement phases and phase transitions in monitored free fermion system(s) due to localization(s)"** (the title varies between versions) [T; authors not verified]. It reports that the BKT transition of monitored free-fermion chains is robust to localization, with the boundary shifting as the localization potential grows. [arXiv HTML](https://arxiv.org/html/2509.09538)
- **"Controlled Zeno-induced localization of free fermions in a quasiperiodic chain"**: arXiv:2602.12340 (Feb 2026); Phys. Rev. B (DOI 10.1103/27g7-tvqh; volume/page not seen) [V-arXiv + APS DOI; authors not verified].
  - Model: the monitored Aubry–André–Harper chain in the quantum-Zeno regime.
  - Method: a self-consistent effective potential built from single trajectories, without post-selection. It equals the quasiperiodic potential plus a uniform imaginary measurement term (plus fluctuation corrections), so the problem reduces to non-Hermitian localization.
  - Sources: [arXiv HTML](https://arxiv.org/html/2602.12340); [APS](https://journals.aps.org/prb/abstract/10.1103/27g7-tvqh)
- arXiv:2607.03799 (Jul 2026), "Measurement-induced spatially nonuniform fluctuations of the local particle number and their crossover in a quasiperiodic free-fermion chain" [T]. [arXiv](https://arxiv.org/pdf/2607.03799)
- Random long-range hopping (Fan & García-García 2026, Question 2) gives multifractal density–density correlations at the transition. [arXiv](https://arxiv.org/abs/2607.23871)
- 2D: monitoring-induced Anderson-type physics, i.e. multifractality and conductance analogies (Chahine–Buchhold 2024; Poboiko–Gornyi–Mirlin 2025; Question 1). [APS PRB 110, 054313](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.110.054313); [APS PRB 112, L220202](https://journals.aps.org/prb/abstract/10.1103/2k3q-h6lz)

### Inferences
- **Disputed status, not superseded:** the 1D "critical phase stable to finite disorder" (2023) and the "log phase up to V_c ≈ 2.3" (2025) claims rest on finite-size numerics. Asymptotic arguments for the clean 1D case predict an area law at all rates (PRX 2023; Fava et al. 2024), and I found no source extending them to disorder. A search summary paraphrased one paper as saying disorder "may stabilize the critical phase, even if it is absent in the clean case". This remains open.
- **Evidence strength:**
  - 1D disorder: numerics plus BKT-phenomenology fits; moderate to weak for thermodynamic-limit claims.
  - 2D Anderson analogy: strong (analytic plus numerics).

### Gaps
- No d ≥ 2 study of disordered or quasiperiodic monitored free fermions (i.e. disorder added on top of the 2D MIPT) was found.
- The authors of arXiv:2509.09538 and arXiv:2602.12340 were not verified, and the journal references for Matsubara et al. are unverified.
- Critical exponents for the disordered 1D transitions were not retrieved, beyond the statement that they are BKT-consistent.

---

## 5. Large-N, all-to-all, Brownian SYK₂, N-flavor chains and trees: where is the transition controlled or exactly solvable?

### Takeaway
At N = ∞, monitored Brownian SYK₂ chains are exactly solvable via a replica saddle. Spontaneous breaking of an O(2) (O(2)×O(2) → O(2)) replica symmetry gives a Goldstone log phase that turns into an area-law phase at a critical rate; interactions reduce O(2) to C₄ and give volume law. Large flavor number also controls the U(1) NLSM (Fava et al.): in d = 1 it still ends in an area law at all rates, while in d > 1 there is a genuine transition. Quantum-tree MIPTs with real (Born-rule) measurements are analysed for qubits; **I found no free-fermion tree model**.

### Cited Findings
- **Jian, Liu, Chen, Swingle, Zhang, "Measurement-induced phase transition in the monitored Sachdev–Ye–Kitaev model"** (v1 title "SYK meets non-Hermiticity II: measurement-induced phase transition"): arXiv:2104.08270 (Apr 2021); Phys. Rev. Lett. 127, 140601 (2021) [V].
  - Model: Brownian SYK chains under continuous monitoring, with the effective action derived at large N.
  - SYK₂ (non-interacting): a continuous O(2) replica symmetry. The symmetry-broken phase at low rate has a Goldstone mode and **log entanglement**; the symmetric phase at high rate is **area law**.
  - Interacting case: C₄ symmetry and volume law.
  - The large-N critical exponents are checked by solving the Schwinger–Dyson equations numerically.
  - Sources: [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.127.140601); [arXiv](https://arxiv.org/abs/2104.08270); [PubMed](https://pubmed.ncbi.nlm.nih.gov/34652178/)
- **Zhang, Jian, Liu, Chen, "Emergent replica conformal symmetry in non-Hermitian SYK₂ chains"**: arXiv:2104.04088 (Apr 2021; v2 28 Oct 2021); Quantum 5, 579 (16 Nov 2021) [V].
  - Models: static and Brownian non-Hermitian SYK₂ chains.
  - Mechanism: O(2)×O(2) replica symmetry broken to O(2) by the saddle, giving an emergent CFT of Goldstone modes.
  - Result: S(L_A) equals the energy of a half-vortex pair, so it scales logarithmically with a coefficient set by the Goldstone stiffness.
  - Sources: [Quantum](https://quantum-journal.org/papers/q-2021-11-16-579/); [arXiv](https://arxiv.org/pdf/2104.04088)
- **Zhang, Liu, Jian, Chen, Quantum 6, 723 (2022)**, arXiv:2105.08895 [V]: the large-N long-range version (Question 2), with the same phase diagram for the single-flavor chain. [Quantum](https://quantum-journal.org/papers/q-2022-05-27-723/)
- **Fava, Piroli, Bernard, Nahum, Phys. Rev. Research 6, 043246 (2024)**, arXiv:2407.08045 [V]: the large-flavor SU(N) NLSM plus hydrodynamics for U(1) fermions; area law at all rates in 1D, and a transition in d > 1 (Question 1). [APS](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.6.043246)
- *One-line context, N-flavor Majorana chains:* Fava, Piroli, Swann, Bernard, Nahum, "Nonlinear sigma models for monitored dynamics of free fermions", arXiv:2302.12820 (Feb 2023; v2 Dec 2023 adds the phase diagram for more than one flavor); Phys. Rev. X 13, 041045 (2023) [V]. Multiflavor Majorana chain; NLSM for an orthogonal N×N matrix in the replica limit N → 1; (log L)² entanglement in the nontrivial phase; disentangled phases from proliferating topological defects. [arXiv](https://arxiv.org/abs/2302.12820); [DOAJ](https://doaj.org/article/27e83e747e3349c8a9713f413fd1b5e9)
- **Feng, Skinner, Nahum, "Measurement-induced phase transitions on dynamical quantum trees"**: arXiv:2210.07264 (Oct 2022); PRX Quantum 4, 030333 (2023) [V].
  - Model: qubit processes with **real** measurements on tree-like spacetime graphs (collapsing or expanding).
  - Results: a transition at a nontrivial measurement strength for both real and forced measurements. The real-measurement case lies at the boundary between two types of critical scaling and has a smaller entangling phase. The work is motivated by the lack of exact solutions for Born-rule MIPTs, even on trees.
  - An expansion-process experimental protocol is suggested.
  - Sources: [arXiv](https://arxiv.org/abs/2210.07264v1); [DOAJ record](https://libmast.utm.my/Record/doaj.art-8ffab0ad52d745c5b26e31b09faf7b2f)
- Predecessor on trees with forced measurements, all-to-all circuits and Landau–Ginsburg theory: Nahum, Roy, Skinner, Ruhman, PRX Quantum 2, 010352 (2021) [journal reference V via DOI; **arXiv ID (2009.11311) from memory, NOT verified this session**]. [DOI](https://doi.org/10.1103/PRXQuantum.2.010352)
- All-to-all connectivity appears experimentally only in qubit (not free-fermion) cross-entropy experiments (Kamakari et al., Question 7).

### Inferences
Where the problem is controlled:
1. **N = ∞ Brownian SYK₂ (O(2) Goldstone saddle):** exact, but at N = ∞ it gives a stable 1D log phase.
2. **Large flavor number with fluctuations (Fava et al. 2024):** shows that for U(1) fermions in 1D the log phase eventually gives way to an area law at all rates. So the 1D large-N "log ↔ area" transition of Jian et al. / Zhang et al. is best read as an artifact of N = ∞ for the U(1) class. This is my inference from combining the two sources; neither states it in the retrieved text.
3. **d > 1:** large-N and d = 1+ε control are consistent with a genuine transition.
4. **Majorana classes:** the metallic phase survives at finite N even in 1D ((log L)²), so the class matters.

Trees are exactly analysable only for qubits so far.

### Gaps
- **No free-fermion (Gaussian/matchgate) version of the quantum-tree MIPT was found** (one targeted search).
- I found no finite-N numerics for monitored SYK₂ clusters in all-to-all geometry. arXiv:2301.05195, "A numerical study of measurement-induced phase transitions in the Sachdev–Ye–Kitaev model", was seen as a title only: [arXiv](https://arxiv.org/html/2301.05195.pdf).
- The large-N critical exponent values of Jian et al. were not retrieved.
- Liu, Zhang, Chen, "Non-unitary dynamics of Sachdev–Ye–Kitaev chain" (SciPost Phys. 2021, arXiv:2008.11955) was not verified this session.

---

## 6. Adaptive / feedback / steering protocols, absorbing-state transitions, and monitored open fermions (dephasing, loss)

### Takeaway
For free fermions, adaptive (feedback) protocols produce **steering/absorbing-state transitions**: parity-conserving (PC) universality in Ravindranath–Yang–Chen 2025. These are distinct from the entanglement transition and are visible at the level of the averaged density matrix. Feedback that creates a skin effect instead **kills** the entanglement transition, even with long-range hopping. Open-system effects (dephasing, losses) add an effective temperature or new non-Hermitian structure. In a 1D lossy chain, a log to area transition is reported; the paper's v1 called it a crossover.

### Cited Findings
- **Ravindranath, Yang, Chen, "Free fermions under adaptive quantum dynamics"**: arXiv:2306.16595 (Jun 2023); Quantum 9, 1685 (3 Apr 2025) [V].
  - Protocol: free-fermion unitaries, projective measurements and corrective unitaries, with a classical flag on each site.
  - Results: an entanglement transition and a **steering transition** into a charge-density-wave state, which are fundamentally distinct. The steering transition is in the **parity-conserving universality class**, attributed to fermion parity combined with the classical labelling. The study is numerical.
  - Sources: [Quantum](https://quantum-journal.org/papers/q-2025-04-03-1685/); [arXiv HTML](https://arxiv.org/html/2306.16595v3)
- **Ravindranath, Han, Yang, Chen, "Entanglement steering in adaptive circuits with feedback"**: arXiv:2211.05162 (Nov 2022); Phys. Rev. B 108, L041103 (12 Jul 2023) [V]. Not described as free-fermion in the sources seen.
  - Adaptive random circuits with feedback show transitions both for trajectories and for the averaged channel.
  - The PC critical line satisfies p × r ≈ 0.55, and the trajectory entanglement transition is bounded above by p_c^EE ≈ 0.45; the two transitions are distinct.
  - Sources: [APS accepted](https://link.aps.org/accepted/10.1103/PhysRevB.108.L041103); [arXiv](https://arxiv.org/pdf/2211.05162)
- **Buchhold, Müller, Diehl, "Revealing measurement-induced phase transitions by pre-selection"**: arXiv:2208.10506 (Aug 2022) [V-arXiv; **no journal reference found**].
  - Proposal: fix the outcomes in advance ("pre-selection"), via feedback, to create a unique dark/absorbing state with macroscopic order. This replaces the maximally mixed state of the unconditioned ensemble.
  - The MIPT then becomes visible in the averaged density matrix through a local order parameter, which links it to quantum absorbing-state transitions.
  - Sources: [arXiv](https://arxiv.org/pdf/2208.10506); [ar5iv](https://ar5iv.labs.arxiv.org/html/2208.10506)
- The 2022 PRL of Müller, Diehl, Buchhold frames the long-range transitions as dark-state transitions (Question 2). [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.010605)
- arXiv:2309.01315 (2023), on steering-induced transitions in measurement-only circuits, uses pre-selection to make such phases more accessible experimentally [T]. [arXiv HTML](https://arxiv.org/html/2309.01315v3)
- **Wang, Fang, Ren, "Absence of (measurement-induced) entanglement transition due to feedback-induced skin effect"** (the title changed between versions): arXiv:2209.11241 (Sep 2022) [V-arXiv; **journal reference not verified**].
  - "Generalized monitoring" (measurement plus conditional feedback) with open boundaries makes particles pile up at an edge, similar to the non-Hermitian skin effect.
  - The result is short-range entanglement and no transition.
  - The authors propose trapped ions; no post-selection is needed, and the physics is observable at the level of particle number.
  - Sources: [arXiv v3](https://arxiv.org/abs/2209.11241v3); [arXiv PDF](https://arxiv.org/pdf/2209.11241)
- Long-range extension: arXiv:2212.08090 (Dec 2022), "Absence of logarithmic and algebraic scaling entanglement phases due to skin effect". Power-law hopping does not restore a transition [T + secondary; the authors, "Feng et al.", were not verified]. [arXiv](https://arxiv.org/pdf/2212.08090)
- **Ladewig, Diehl, Buchhold, "Monitored open fermion dynamics: Exploring the interplay of measurement, decoherence, and free Hamiltonian evolution"**: arXiv:2203.00027 (Mar 2022); Phys. Rev. Research 4, 033001 (2022) [V].
  - Setting: dephasing drives the state into a non-Gaussian mixed state.
  - Results: at strong dephasing a new scale (an effective temperature) emerges with increased mixedness, yet density–density correlations and subsystem parity remain scale invariant.
  - Methods: exact small-system solution, correlation-matrix trajectories, and a bosonic replica field theory with perturbative RG.
  - Sources: [arXiv](https://arxiv.org/abs/2203.00027); [Cologne repository](https://kups.ub.uni-koeln.de/65996)
- **Soares, Le Gal, Schirò, "Entanglement transition due to particle losses in a monitored fermionic chain"**: arXiv:2408.03700 (Aug 2024); Phys. Rev. B 111, 064313 (2025) [V-sec: the DOI comes from the arXiv and IP-Paris listings; the APS page was not seen].
  - Model: a 1D free-fermion chain with hopping and pairing, where local losses act as quantum jumps. The pairing interpolates between number- and parity-conserving drives.
  - Result: a **log to area** entanglement transition, in qualitative agreement with the no-click (non-Hermitian) prediction over most of the phase diagram.
  - Method: the Faber-polynomial method plus Monte Carlo wave-function sampling.
  - **Version change:** v1 described a crossover; the published version describes a transition.
  - Sources: [arXiv](https://www.arxiv.org/abs/2408.03700); [IP-Paris](https://researchportal.ip-paris.fr/en/publications/entanglement-transition-due-to-particle-losses-in-a-monitored-fer/)
- No-click (non-Hermitian) reference: Le Gal, Turkeshi, Schirò, "Volume-to-area law entanglement transition in a non-Hermitian free fermionic chain", arXiv:2210.11937; SciPost Phys. 14, 138 (2023) [V-sec: the journal reference was seen only in a reference list]. [arXiv](https://arxiv.org/pdf/2210.11937)
- arXiv:2606.09306 (Jun 2026), "Range-controlled entanglement in Lindbladian skin states of monitored fermions" [T]. [arXiv](https://arxiv.org/pdf/2606.09306)

### Inferences
- In free fermions, adaptive protocols move the "transition" into the averaged dynamics: an absorbing-state transition with PC or directed-percolation-like character. This is observable without post-selection, but it is a **different** transition from the trajectory entanglement transition (stated by Ravindranath et al. in both papers).
- **Evidence strength:**
  - Steering universality: good numerical evidence.
  - Lossy 1D log to area transition: moderate, given the v1 → v2 change from crossover to transition and the general 1D finite-size caveats.

### Gaps
- No d ≥ 2 adaptive free-fermion study was found.
- Journal references for Wang–Fang–Ren and for the pre-selection paper were not verified.
- The authors of arXiv:2212.08090 were not verified.
- Starchl et al., "Generalized Zeno effect and entanglement dynamics induced by fermion counting" (arXiv:2406.07673, 1D), was seen as a title only [T]: [arXiv](https://arxiv.org/pdf/2406.07673).

---

## 7. Experiments (2022–2026) on monitored free-fermion / matchgate dynamics, the post-selection problem, and workarounds

### Takeaway
I found **no experiment (through Nov 2026) that realizes a monitored free-fermion or matchgate MIPT**; two targeted searches returned only theory. The closest hardware work is the qubit cross-entropy demonstration on IBM devices (Kamakari et al., PRL 2025; circuit family not confirmed as free-fermion). Free-fermion-relevant workarounds for post-selection are mostly proposals:
- pre-selection / feedback to dark states;
- particle-number-level observables in feedback protocols;
- ultracold-atom and trapped-ion proposals;
- cross-entropy benchmarks, which should be applicable because Gaussian circuits are classically simulable (my inference).

### Cited Findings
- A search for hardware realizations of monitored free-fermion (Gaussian/matchgate) MIPTs on trapped-ion, superconducting or ultracold platforms found none; the results were theory papers and reference lists (searches run Nov 2026). Representative theory hits: [arXiv:2503.21427](https://arxiv.org/html/2503.21427) and [PRX 13, 041046](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.13.041046).
- **Li, Zou, Glorioso, Altman, Fisher, "Cross entropy benchmark for measurement-induced phase transitions"**: arXiv:2209.00609 (Sep 2022); Phys. Rev. Lett. 130, 220404 (2023) [V].
  - Method: the linear cross-entropy between measurement-record distributions of two circuits that share a bulk but differ in initial state, which avoids post-selection.
  - Result: χ = 1 in the volume-law phase (thermodynamic limit) and χ < 1 in the area-law phase.
  - Sources: [arXiv](https://arxiv.org/pdf/2209.00609); [APS accepted](https://link.aps.org/accepted/10.1103/PhysRevLett.130.220404)
- **Kamakari, Sun, Li, Thio, Gujarati, Fisher, Motta, Minnich, "Experimental demonstration of scalable cross-entropy benchmarking to detect measurement-induced phase transitions on a superconducting quantum processor"**: arXiv:2403.00938 (Mar 2024; v2 May 2025); Phys. Rev. Lett. 134, 120401 (28 Mar 2025) [V].
  - IBM devices with 1D and **all-to-all** connectivity, up to 22 qubits, in under 8 device-hours.
  - The data collapse with critical exponents that agree with theory within uncertainty.
  - The gate set (Clifford, free-fermion or other) was not confirmed in the retrieved text.
  - Sources: [arXiv](https://arxiv.org/abs/2403.00938v1); [Nature Index record](https://www.nature.com/nature-index/article/10.1103/physrevlett.134.120401); [IBM Research](https://research.ibm.com/publications/experimental-demonstration-of-scalable-cross-entropy-benchmarking-to-detect-measurement-induced-phase-transitions-on-a-superconducting-quantum-processor)
- Free-fermion-specific proposals and workarounds:
  - Pre-selection to dark/absorbing states with a local order parameter (Buchhold–Müller–Diehl, arXiv:2208.10506). [arXiv](https://arxiv.org/pdf/2208.10506)
  - Feedback-induced skin effect: trapped-ion proposal with no post-selection and particle-number-level readout (Wang–Fang–Ren, arXiv:2209.11241). [arXiv](https://arxiv.org/abs/2209.11241v3)
  - Ultracold atoms under continuous measurement (Yamamoto–Hamazaki, PRB 107, L220201). [arXiv](https://arxiv.org/pdf/2301.07290)
  - Steering transitions visible in the averaged channel (Ravindranath et al., PRB 2023; Quantum 2025). [APS accepted](https://link.aps.org/accepted/10.1103/PhysRevB.108.L041103)
  - A post-selection-free analytic description of typical trajectories in the Zeno-regime quasiperiodic chain (arXiv:2602.12340). This is a theory tool, not an experimental protocol. [arXiv HTML](https://arxiv.org/html/2602.12340)
- Further related titles, not verified beyond the title [T]:
  - "Theory of free fermions dynamics under partial post-selected monitoring", arXiv:2312.14022: [arXiv HTML](https://arxiv.org/html/2312.14022v1)
  - "Measurement-induced phases of matter require feedback", arXiv:2210.07256: [arXiv](https://arxiv.org/pdf/2210.07256)
  - Two 2026 preprints whose nature (experiment or theory) could not be determined: "Measurement-enhanced entanglement in a monitored superconducting chain", arXiv:2604.04375 ([arXiv](https://arxiv.org/pdf/2604.04375)), and "Super-logarithmic entanglement scaling in a monitored superconducting chain", arXiv:2607.07835 ([arXiv](https://arxiv.org/pdf/2607.07835)).

### Inferences
- The Li et al. volume vs area-law criterion does not transfer directly to free fermions, which have no volume-law phase. For free fermions a cross-entropy or purification-type diagnostic would have to resolve the metallic (ℓ^{d−1} ln ℓ) vs area-law distinction, presumably through purification timescales. Its classical side is efficient for Gaussian circuits. This is my inference; I found no paper demonstrating it.
- The most realistic routes to a free-fermion MIPT experiment, as of Nov 2026, appear to be:
  - 2D (d ≥ 2) monitored U(1) or matchgate circuits, where the transition genuinely exists, using particle-number second cumulants (Poboiko et al.'s observable) combined with cross-entropy or decoders;
  - feedback/steering variants whose transitions appear in averaged observables.

  This is an inference.

### Gaps
- Because the search budget was exhausted, I could not verify the details (arXiv IDs, volumes) of the standard non-free-fermion MIPT experiments mentioned in search summaries:
  - Koh, Sun, Motta, Minnich (superconducting, mid-circuit readout, Nature Physics 2023);
  - Google Quantum AI (Nature 2023);
  - Noel et al. (trapped ions);
  - Agrawal et al. (learnability transition);
  - the trapped-ion implementation of the Feng–Skinner–Nahum tree protocol mentioned in one summary.
- No neural-network or decoder workaround specific to free fermions was found.
- The gate set of the Kamakari et al. experiment was not confirmed.
