# Monitored free fermions without U(1): symmetry classes, Born-rule vs forced measurements, Majorana loop models, and the non-Hermitian / no-click limit (notes as of 22 Nov 2026)

*How these were verified (read first).* WebFetch could not reach arxiv.org, export.arxiv.org, web3.arxiv.org or api.semanticscholar.org from this session (DNS not allowed). Every arXiv ID and journal reference below was therefore checked through WebSearch results that rendered the arXiv abstract page, the publisher page (APS, SciPost, Quantum, Springer), DOAJ, or institutional repositories (NSF-PAR, SISSA IRIS, IP Paris, Lancaster). I had abstract-level access only, never full texts. The shared web-search budget ran out before the last checks, so some items are flagged **[journal unverified]** or **[unverified]**. Status legend: **[V]** means the arXiv ID and journal reference were both confirmed. **[V-arXiv]** means the arXiv ID was confirmed but no journal version was found or confirmed. Dates give the arXiv posting month and the journal year.

U(1) context (one line only; another researcher covers it). For number-conserving 1D free fermions the current field-theory and numerics picture is area law at any finite monitoring rate, with no genuine transition at finite rate. The SU(N) NLSM of Fava–Piroli–Bernard–Nahum is in [arXiv:2407.08045](https://arxiv.org/abs/2407.08045) [V-arXiv]. GPU numerics showed that L ~ 10^4 is needed to see this ([arXiv:2508.18468](https://arxiv.org/html/2508.18468)). The Majorana / no-U(1) case below is the qualitatively different one.

## 1. Classification and NLSM: Jian–Bauer–Keselman–Ludwig (PRB 2022) and Fava–Piroli–Swann–Bernard–Nahum (PRX 2023)

### Takeaway
Generic Gaussian Majorana dynamics in 1+1D has only fermion parity as a symmetry. It maps to Altland–Zirnbauer class DIII, and its effective theory is the principal chiral NLSM on SO(N). The Born rule gives the replica limit N→1 and forced (or Haar tensor-network) ensembles give N→0. In both limits there is a stable "metallic" / nontrivially entangled phase with S ~ (ln L)^2, which is asymptotically exact, plus a stable area-law phase. They are separated by a genuine transition driven by vortex (Z2-defect) unbinding.

### Cited Findings
- **JBKL, arXiv:2012.04666 (Dec 2020); Phys. Rev. B 106, 134206 (2022) [V].**
  - The most general Haar-random Gaussian fermionic tensor-network (non-unitary circuit) ensemble corresponds to a *unitary* non-interacting fermion problem with static disorder in AZ class DIII, in D = d+1 dimensions. Class DIII hosts a stable critical metallic phase in D = 2 and D = 3 ([arXiv](https://arxiv.org/abs/2012.04666); journal ref via [Cornell listing](https://www.lassp.cornell.edu/bibcite/reference/1174)).
  - The other nine AZ classes can be reached from DIII by Clifford-algebra extensions, i.e. by imposing extra symmetries ([arXiv](https://arxiv.org/abs/2012.04666)).
  - Numerics: boundary states of D = 2 and D = 3 tensor networks ([arXiv](https://arxiv.org/abs/2012.04666)).
- JBKL's contents list sections on "the critical entanglement phase as the symmetry-class-DIII disordered metallic phase" and on "logarithmic corrections to scaling of the 2nd moment of the fermion correlation function in the DIII metallic phase" ([ar5iv](https://ar5iv.labs.arxiv.org/html/2012.04666)).
- **Fava, Piroli, Swann, Bernard, Nahum, arXiv:2302.12820 (24 Feb 2023; v2 8 Dec 2023); Phys. Rev. X 13, 041045 (Dec 2023) [V].**
  - Model and method: a multi-flavor Majorana chain with Hamiltonian evolution and continuous monitoring of local fermion-parity operators. The replica trick maps the dynamics to imaginary-time evolution of an effective spin chain, with the number of flavors as the large control parameter ([arXiv](https://arxiv.org/abs/2302.12820); [PRX](https://link.aps.org/doi/10.1103/PhysRevX.13.041045)).
  - Effective theory: an NLSM for an orthogonal N×N matrix (SO(N)) in the replica limit N→1 ([PRX](https://link.aps.org/doi/10.1103/PhysRevX.13.041045)).
  - The required replica limit is N→1, not the N→0 of familiar disordered systems. The extra replica encodes the Born-rule probability of a trajectory ([ar5iv](https://ar5iv.labs.arxiv.org/html/2302.12820)).
  - Phases and transitions: "together with known results for the RG beta function", the derivation establishes stable nontrivially entangled and disentangled phases at N→1. Transitions out of the nontrivial phase map to vortex-unbinding transitions of the sigma model. The phase diagram also has separate critical points on the measurement-only axis ([PRX](https://link.aps.org/doi/10.1103/PhysRevX.13.041045); [DOAJ record](https://doaj.org/article/27e83e747e3349c8a9713f413fd1b5e9)).
  - Scaling: in the nontrivial phase an asymptotically exact calculation gives a bipartite entanglement S ~ (ln L)^2, "in contrast to findings in previously-studied models". Rényi entropies in that phase take a universal form ([PRX](https://link.aps.org/doi/10.1103/PhysRevX.13.041045)).
  - Version note: the accepted v2 (Dec 2023) adds a discussion of the phase diagram for more than one flavor ([arXiv record](https://arxiv.org/abs/2302.12820)).
- **Jian–Shapourian–Bauer–Ludwig (JSBL), arXiv:2302.09094 (Feb 2023) [V-arXiv].** The transition is described by the principal chiral NLSM with an SO(n) field, with n→0 for forced and n→1 for Born-rule measurements. "Both transitions are driven by proliferating Z2 topological defects" ([ar5iv](https://ar5iv.labs.arxiv.org/html/2302.09094)).
- Contrast with U(1): the U(1) follow-up uses an SU(N) NLSM rather than SO(N). It finds an area law for all monitoring rates, with nontrivial scaling when the correlation length is large ([arXiv:2407.08045](https://arxiv.org/abs/2407.08045)).
- 2025 GPU paper, background statement: the NLSM predicts a MIPT at *finite* monitoring rate for Majorana fermions under homodyne parity measurements, and none for Dirac (U(1)) fermions ([arXiv:2508.18468](https://arxiv.org/html/2508.18468); [PRB DOI 10.1103/x7vx-tqr4](https://journals.aps.org/prb/abstract/10.1103/x7vx-tqr4)).

### Inferences
- **Target-space dictionary.** JBKL's class-DIII statement, Fava's SO(N) at N→1, and JSBL's SO(n) at n→0 or n→1 describe the same principal chiral model. They differ only in the replica limit. This is my synthesis; for DIII, the fermionic-replica NLSM target is the orthogonal group, which is consistent with all three abstracts.
- **Sign of the beta function.** This is standard background, not re-verified here. The one-loop beta function of the SO(N) PCM is proportional to the dual Coxeter number, N−2. For N > 2 the model is asymptotically free and flows to the disordered, strong-coupling regime. For N < 2, which covers both N→1 and N→0, the coupling flows to weak coupling and the stiffness grows like ln L. This is weak antilocalization in class DIII.
  - Since S ~ (stiffness) × ln L, the stiffness growth gives S ~ (ln L)^2.
  - So the replica limit, not the group structure, makes the "metal" stable.
- **Transition mechanism.** π1(SO(N)) = Z2 for N ≥ 3. Fava's "vortex unbinding" and JSBL's "Z2 topological defects" therefore name the same mechanism. The area-law phase is the defect-proliferated phase.
- **Naming.** The "metallic" phase (JBKL, Anderson language), the "nontrivially entangled" phase (Fava) and the "Goldstone" phase (loop models, §3) all have (ln L)^2 entanglement.
- **Slow approach to (ln L)^2.** The coupling flows only logarithmically, so numerics at moderate L can look like S ≈ c_eff ln L with a slowly drifting c_eff. This plausibly explains the earlier "log-law critical phase" claims (§4).
- **Answer to (i), symmetry classes in 1D:**
  - (a) Generic Majorana / no U(1), class DIII: a stable critical (ln L)^2 phase and a genuine transition. This holds with Born-rule and with forced measurements.
  - (b) U(1)-conserving: no stable critical phase asymptotically (one line, covered elsewhere).
  - (c) Measurement-only Gaussian circuits: a stable (ln L)^2 Goldstone phase only if loops are non-orientable. If loops are orientable there are only critical points (percolation/Potts) between area-law phases (§3).
  - (d) Deterministic no-click (non-Hermitian) dynamics: exact volume, log or area-law phases are possible, but they are not generic under the Born rule (§5).

### Gaps
- No full texts were read. Still missing:
  - the explicit beta-function coefficients and the universal prefactor of (ln L)^2 in Fava et al.;
  - the system sizes used in the Fava et al. and JBKL numerics;
  - whether JBKL explicitly state the replica limit n→0 for their Haar ensemble (it is implied, but not verified).
- I could not verify whether the single-flavor nearest-neighbour Majorana chain (the monitored Ising/Kitaev chain) has an extra symmetry that changes its class. Fava's control parameter is a large number of flavors, so the single-flavor case relies on a universality argument.

## 2. Born rule vs forced (postselection-free) measurements

### Takeaway
For generic 1D Majorana circuits (class DIII), Born-rule and forced measurements share the AZ class and both have the (ln L)^2 phase and the transition. The transitions lie in *different universality classes*: forced is the DIII Anderson-localization class, while Born-rule is a new class (numerics from JSBL, 2023). The *existence* of a metallic phase can depend on the replica limit in other classes, e.g. class D in decoherence settings. Full postselection (no-click) is a third, deterministic limit with different phases.

### Cited Findings
- **JSBL, arXiv:2302.09094 (Feb 2023). No journal publication was found as of Nov 2026, but the search budget ran out before this was exhaustive [V-arXiv].**
  - Model: a generic circuit on a 1D Majorana chain with only fermion parity as a symmetry.
  - Result: several numerically extracted critical exponents give clear evidence that the Born-rule and forced-measurement transitions lie in different universality classes. The forced transition is in an Anderson-localization universality class; the Born-rule transition is in "a novel universality class beyond Anderson localization".
  - Both share the same AZ class (DIII), so symmetry class alone does not explain the difference ([arXiv](https://arxiv.org/abs/2302.09094); [ar5iv](https://ar5iv.labs.arxiv.org/html/2302.09094)).
  - Methods: a method to extract the correlation-length exponent ν, and a Monte Carlo scheme to sample Born-rule trajectories. The paper cites the result that monitored Gaussian circuits with forced measurements are equivalent to Anderson-localization problems ([ar5iv](https://ar5iv.labs.arxiv.org/html/2302.09094)).
- Fava et al. (PRX 2023) make the general point that the Born rule requires N→1 rather than the N→0 of quenched disorder ([ar5iv](https://ar5iv.labs.arxiv.org/html/2302.12820)).
- **Counter-example where existence depends on the limit (class D, decoherence context).** A 2025 paper on decoherence-induced self-dual criticality in topological states (arXiv:2502.14034, Feb 2025) [journal unverified] states: "the N→0 limit, corresponding to generic quenched disorder, is dramatically different and moreover exhibits a Majorana metallic phase which is absent in the N→1 limit" ([arXiv](https://arxiv.org/html/2502.14034v2)).
- U(1) one-liner: Poboiko–Pöpperl–Gornyi–Mirlin (arXiv:2304.03138) [journal unverified in this session] note that the Born rule needs a different replica limit from forced measurements. Their analysis argues against a MIPT for U(1) free fermions ([arXiv](https://pdf.arxiv.org/pdf/2304.03138), as summarized in search results).
- Measured critical states, not dynamics: Yang–Mao–Jian (arXiv:2301.08255, Jan 2023) [journal unverified] measured the critical TFIM ground state. With an optimally *biased* forced measurement, the effective central charge returns to its unmeasured value independent of measurement strength, a feature shared with Born-rule averages ([arXiv](https://arxiv.org/pdf/2301.08255v1)).
- **Fully postselected, exactly solvable case: Granet, Zhang, Dreyer, arXiv:2212.10584 (Dec 2022); Phys. Rev. Lett. 130, 230401 (2023) [V].** These are translation-invariant, time-periodic Gaussian circuits with post-selected weak measurements. The critical exponent ν is computed *exactly* for the Hartley, von Neumann and Rényi entropies. There are volume-to-area and log-to-area transitions ([arXiv](https://arxiv.org/abs/2212.10584)). Correction: this paper is sometimes misattributed to "Granet–Zhang–Dalmonte". The arXiv record lists Dreyer.
- In partially postselected monitored Kitaev chains (arXiv:2407.13837, Jul 2024) [authors, title and journal unverified], imperfect postselection destroys the no-click critical phase. Only the fully postselected limit can be critical ([arXiv](https://arxiv.org/html/2407.13837v2)).

### Inferences
- In DIII, both replica limits lie in N < 2, so the phase structure is robust to Born vs forced: one (ln L)^2 phase and one area-law phase, with a Z2-defect-driven transition. Only the critical point changes. The forced (N→0) transition is the 2D class-DIII Anderson metal–insulator transition; the Born (N→1) transition is new.
- In classes where the leading beta-function coefficient changes sign between N = 0 and N = 1, the *existence* of the metal can depend on the measurement protocol. The 2025 class-D statement is such an example.
  - Lessons from Anderson localization therefore transfer reliably only to forced-measurement ensembles.
  - This is my inference; the class-D coefficient was not verified here.
- Three different ensembles must not be conflated:
  1. Born-rule trajectories (N→1);
  2. forced / uniform-weight ensembles (N→0, Anderson);
  3. the single postselected no-click trajectory. This is deterministic, can show volume law, and is fragile to admixture of jumps (§5).

### Gaps
- Numerical values of ν and the other exponents, system sizes, and the journal status of JSBL were not retrieved.
- I found no paper that computes the Born-rule (N→1) critical exponents of the DIII PCM analytically.

## 3. Measurement-only and Clifford Majorana circuits, and loop models

### Takeaway
Gaussian measurement-only Majorana circuits map *exactly* to 2D loop models (Klocke–Buchhold, PRX 2023), and entanglement is counted by the loop strands (arcs) that cross the cut.
- **Orientable loops:** Potts-type physics (percolation, Ising, coupled Potts). There are only critical points, which separate distinct (e.g. topological vs trivial) area-law phases, with ln L entanglement at criticality. Examples are the projective TFIM, Lavasani–Alavirad–Barkeshli and Sang–Hsieh.
- **Non-orientable loops (crossings):** the CPLC universality class (Nahum et al., PRB 2013), with a stable Goldstone phase where S ~ log^2 L. Its continuum description is the RP^{n−1} sigma model at n = 1.
- **Weak (non-projective) measurements:** these destabilize the percolation universality.

### Cited Findings
- **Klocke & Buchhold, arXiv:2305.18559 (29 May 2023); Phys. Rev. X 13, 041028 (2023), DOI 10.1103/PhysRevX.13.041028, published 10 Nov 2023 [V].**
  - Mapping: 1D measurement-only circuits map to 2D loop models. Gaussian Majorana circuits admit a *microscopic* mapping; for generic Clifford circuits the mapping may emerge only on a coarse-grained scale ([arXiv PDF](https://www.arxiv.org/pdf/2305.18559); [bibliographic record](https://www.bohrium.com/paper-details/majorana-loop-models-for-measurement-only-quantum-circuits/929934851790864619-479)).
  - Orientable case: depending on the measured operators, the circuits realize the 1-state Potts model (percolation), the 2-state Potts model (Ising) and coupled Potts models.
  - Non-orientable case: when orientability is broken, the circuit falls into the CPLC class and has a Goldstone phase with "a peculiar, universal log²(L)" entanglement ([arXiv PDF](https://www.arxiv.org/pdf/2305.18559)).
- **Nahum, Serna, Somoza, Ortuño, "Loop models with crossings", arXiv:1303.2342 (Mar 2013); Phys. Rev. B 87, 184204 (2013) [V].**
  - Model and phases: the completely packed loop model with crossings (CPLC) on the square lattice. Continuous transitions of a new kind separate a Goldstone phase from short-loop phases.
  - Continuum theory: a replica limit of the sigma model on real projective space RP^{n−1}.
  - Methods: RG compared with Monte Carlo up to L = 10^6. One point of the phase diagram corresponds to self-avoiding trails at their Θ-point ([arXiv](https://www.arxiv.org/abs/1303.2342)).
- Jacobsen–Read–Saleur (cond-mat/0205033) [journal unverified] argued that the Goldstone phase is generic for loop models with crossings when the loop weight lies between −2 and 2 ([arXiv](https://arxiv.org/abs/cond-mat/0205033)).
- Lorentz-mirror connection, indirect support only. arXiv:1502.05305 (OSp(3|2) superspin chain) [journal unverified] describes a lattice Lorentz gas with mirrors tilted by ±π/4. It states that the associated loop model with fugacity below 2 should show the Goldstone-phase behaviour of the sigma model. I could not confirm that Nahum et al. (2013) discuss the mirror model ([arXiv](https://arxiv.org/pdf/1502.05305)).
- **Nahum & Skinner, arXiv:1911.11169 (Nov 2019). Published in Phys. Rev. Research (2020); the volume/article number "2, 023288" is commonly cited but I could NOT verify it [V-arXiv, journal partially verified].**
  - Model: 1D domain-wall defects carrying Majorana modes diffuse because of coupling to a bath. The fermion parity of a pair is nonlocal and is "measured" only when the two defects meet.
  - The late-time density for γ+γ→∅ is exactly twice the classical A+A→∅ value.
  - In one regime only two-body entanglement is generated, and the number of arcs crossing the midpoint grows as log L ([arXiv](https://arxiv.org/abs/1911.11169v1); [ar5iv](https://ar5iv.labs.arxiv.org/html/1911.11169); [DOAJ-indexed record](https://openaccess.library.uitm.edu.my/find/Record/doaj-art-809dcb41d47b4d5fa2af9586de0e9dbc)).
- **Lang & Büchler, "Entanglement transition in the projective transverse field Ising model", arXiv:2006.09748 (Jun 2020) [V-arXiv; journal ref not verified].**
  - Setup: two noncommuting projective measurements and no unitaries.
  - The transition is between two *area-law* steady states, and a classical model reproduces the dynamics.
  - The mapping to bond percolation puts the critical point at p = 0.5, where a CFT gives the universal entanglement scaling ([arXiv](https://arxiv.org/pdf/2006.09748)).
  - Follow-ups:
    - arXiv:2012.03857 finds percolation surface exponents for entanglement clusters; Clifford models deviate ([arXiv](https://www.arxiv.org/pdf/2012.03857)) [authors and journal unverified].
    - arXiv:2511.17370 (Nov 2025), attributed in search results to Roser, Springer, Büchler and Lang, proposes robust experimental protocols to detect or bound the transition ([arXiv](https://arxiv.org/pdf/2511.17370)) [unverified].
- **Sang & Hsieh, "Measurement-protected quantum phases", arXiv:2004.09509 (Apr 2020); Phys. Rev. Research 3, 023200 (2021), DOI 10.1103/PhysRevResearch.3.023200, published 22 Jun 2021 [V].**
  - These are hybrid random-unitary plus projective-measurement circuits with long-range order inside the area-law phase.
  - Glass and trivial phases are perturbatively stable to competing Ising-symmetric measurements.
  - There is a spin-glass area-law to paramagnetic volume-law transition and a critical regime ([arXiv](https://arxiv.org/abs/2004.09509v3); [DOAJ](https://doaj.org/article/db71ade01f054298a8b32d3d218472a6)).
- **Lavasani, Alavirad, Barkeshli, arXiv:2004.07243 (Apr 2020); Nature Physics 17, 342–347 (2021), DOI 10.1038/s41567-020-01112-z [V].**
  - Symmetric random circuits with two competing measurement types plus unitaries. Without unitaries, a purely measurement-induced critical point maps exactly to *two copies of 2D percolation*.
  - The phases are SPT, trivial and volume-law.
  - The transitions appear only in trajectory-averaged entanglement, not in operator expectation values ([Springer](https://link.springer.com/article/10.1038/s41567-020-01112-z); [arXiv](https://arxiv.org/abs/2004.07243)).
- **Merritt & Fidkowski, "Entanglement transitions with free fermions", arXiv:2210.05681 (Oct 2022); Phys. Rev. B 107, 064303 (2023), published 16 Feb 2023 [V].** They use Majorana operators to study 1D entanglement dynamics under random free-fermion unitaries and projective measurements ([arXiv](https://arxiv.org/pdf/2210.05681); [NSF-PAR](https://par.nsf.gov/servlets/purl/10468306)). I did not retrieve their specific phase-diagram claims.
- **Kells, Meidan, Romito, arXiv:2112.09787 (Dec 2021; originally titled "Topological transitions with continuously monitored free fermions"); SciPost Phys. 14, 031 (2023), published 13 Mar 2023 [V].**
  - Two non-commuting sets of *non-projective* measurements stabilize area-law phases of distinct topological order.
  - The transition belongs to a different universality class from stroboscopic projective circuits.
  - With unitaries added, a region of sub-volume entanglement separates the two topological phases.
  - The transitions are diagnosed by combining the bipartite and topological entanglement entropies ([SciPost](https://scipost.org/SciPostPhys.14.3.031); [arXiv](https://arxiv.org/abs/2112.09787v1)).
- **Nehra, Romito, Meidan, arXiv:2404.07918 (Apr 2024); Quantum 9, 1697 (2025), published 8 Apr 2025 [V].** The percolation universality of the projective-limit free-fermion transition is unstable to reducing the system–detector coupling ([Quantum](https://quantum-journal.org/papers/q-2025-04-08-1697/); [arXiv](https://arxiv.org/abs/2404.07918v4)).
- **Klocke, Simm, Zhu, Trebst, Buchhold, arXiv:2409.02171 (Sep 2024); Phys. Rev. B 111, 224301 (2025) [V; journal ref from search-rendered records].**
  - Measurement-only adaptations of the Kitaev and Yao–Kivelson models, in classes BDI and D. "Majorana loop models" are presented as a general framework in arbitrary dimension.
  - The two entanglement phases are separated by a continuous transition with quantum Lifshitz scaling, with exponents in two distinct universality classes. These are linked to NLSMs (CP^{n−1} vs RP^{n−1}) and give efficient numerics.
  - A summary table contrasts circuits (loop fugacity n = 1) with Hamiltonians (n = √2) ([arXiv](https://www.arxiv.org/abs/2409.02171); [Trebst-group PDF](https://www.thp.uni-koeln.de/trebst/pubs/KitaevCircuitsLoopModels.pdf)).

### Inferences
- **Entanglement from arcs (standard fact; no verified source quote).**
  - In Gaussian measurement-only dynamics from a Majorana-dimer state, the state at all times is a Majorana pairing.
  - Each pair with one partner in A and one in B contributes (ln 2)/2, so S_A = (ln 2/2) × (number of loop strands, i.e. arcs, crossing the cut).
  - Expected scaling of the arc count: O(1) in short-loop phases, ~ln L at percolation-type critical points (cf. Nahum–Skinner's arc counting), and ~(ln L)^2 in the Goldstone phase.
- **Orientability.** Orientable loops correspond to CP^{n−1}-type theories at n = 1, which have no Goldstone phase, only percolation/Potts critical points. Non-orientable loops correspond to RP^{n−1} at n = 1, which has a Goldstone phase.
  - This is the standard loop-model assignment.
  - The search tool flagged that the column alignment in the Klocke et al. table was unclear, so check it in the paper.
- **Physical origin of crossings.** Crossings typically come from measuring or swapping Majoranas beyond a bipartite nearest-neighbour pattern; this breaks orientability. Purely nearest-neighbour parity measurements (projective TFIM, Kitaev-chain-type) are orientable. They therefore give percolation criticality between two area-law phases, topological vs trivial, rather than a critical phase.
- **Lorentz mirror model.** The CPLC at n = 1 is equivalent to the Lorentz mirror model: mirrors of two orientations are the two non-crossing node types and empty sites are crossings. Majorana measurement-only circuits with crossings are thus quantum realizations of mirror-model trajectory statistics. This is supported only indirectly by arXiv:1502.05305.
- **Fine tuning.** Loop-model (projective/Clifford) criticality is fine-tuned within generic Gaussian dynamics:
  - weak measurements change the universality class (Kells et al.) and destabilize percolation (Nehra et al.);
  - this is consistent with the "separate critical points on the measurement-only axis" in Fava et al.

### Gaps
- Not obtained:
  - the system sizes and the (ln L)^2 prefactor in Klocke–Buchhold;
  - whether Nahum–Skinner themselves discuss a (ln L)^2 regime (with swaps or braiding). The search tool attributed (ln L)^2 only to Fava et al.;
  - the journal references for Lang–Büchler and Jacobsen–Read–Saleur;
  - whether the Goldstone phase of the RP^0 loop model and the SO(N→1) PCM metal are the same phase once unitaries are added to a non-orientable measurement-only circuit;
  - the detailed results of Merritt–Fidkowski.
- Related Kitaev-circuit work was seen in search results, but the ID–journal pairing is unverified:
  - Lavasani–Luo–Vijay, "Monitored quantum dynamics and the Kitaev spin liquid" ([arXiv:2207.02877](https://pdf.arxiv.org/pdf/2207.02877)), likely Phys. Rev. B 108, 115135 (2023) ([APS accepted MS](https://link.aps.org/accepted/10.1103/PhysRevB.108.115135));
  - Zhu–Tantivasadakarn–Trebst on an *interacting* monitored Majorana spin liquid ([arXiv:2303.17627](https://ar5iv.labs.arxiv.org/html/2303.17627)).

## 4. Monitored quantum Ising / Kitaev chains: is there a true transition?

### Takeaway
The 2021–2023 numerics on the monitored TFIM, using QSD and quantum-jump unravelings at moderate L, reported a log-law "critical phase" with c_eff and a transition or crossover to area law. The results depend on the unraveling.

Since 2023, field theory has reinterpreted this, and 2026 work on a monitored superconducting (pairing) chain supports it: for Born-rule monitored Majorana chains *without* U(1), the nontrivial phase is asymptotically (ln L)^2, with a genuine transition to area law. The "log phase" claims are therefore superseded in interpretation, not in data.

The measurement-only TFIM has a genuine transition, but it is a percolation point between two area-law phases.

### Cited Findings
- **Turkeshi, Biella, Fazio, Dalmonte, Schirò, arXiv:2103.09138 (Mar 2021); Phys. Rev. B 103, 224210 (2021) [V].** They compare quantum-state diffusion ("infinite clicks") with the no-click limit, which is postselected and described by a non-Hermitian Hamiltonian. Both show similar phenomenology: a critical phase goes to an area-law phase as γ increases ([arXiv](https://www.arxiv.org/pdf/2103.09138); [LPTMS/HAL notice](https://www.lptms.universite-paris-saclay.fr/blog/measurement-induced-entanglement-transitions-in-the-quantum-ising-chain-from-infinite-to-zero-clicks-archive-ouverte-hal-4/)). **[Interpretation superseded by Fava et al. 2023 for Born-rule QSD; the no-click part remains exact for that trajectory.]**
- **Piccitto, Russomanno, Rossini, arXiv:2111.11300 (Nov 2021); Phys. Rev. B 105, 064305 (2022). v4 (Dec 2022) corrects Eq. D5 and Figs. 5–6 [V].**
  - The setup is a TFIM with dephasing, under two Gaussianity-preserving unravelings of the *same* Lindbladian.
  - QSD gives a crossover between area-law and log-law scaling, for which they draw a phase diagram. Their own word is a crossover.
  - The QJ unraveling gives only log scaling.
  - The non-Hermitian approximation conflicts with the full results ([arXiv](https://arxiv.org/abs/2111.11300v3); [UniPi record](https://arpi.unipi.it/handle/11568/1164710)).
- **Turkeshi, Dalmonte, Fazio, Schirò, arXiv:2111.03500 (Nov 2021); Phys. Rev. B 105, L241114 (2022); erratum Phys. Rev. B 107, 079901 (2023) [V; erratum details not seen].**
  - Picture: entanglement is carried by ballistic non-Hermitian quasiparticles that measurements stochastically *reset*.
  - For the QJ-monitored Ising chain it predicts a log critical phase, an area-law phase, and a continuous transition with c_eff vanishing as a square root. This agrees with exact numerics ([APS](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.L241114); [arXiv](https://arxiv.org/abs/2111.03500)).
  - Apparent conflict with Piccitto's "QJ ⇒ only log": the two use different jump protocols (monitored observable and Gaussian construction), so the comparison is not like-for-like.
- **Paviglianiti & Silva, arXiv:2302.06477 (Feb 2023); Phys. Rev. B 108, 184302 (2023) [V].** Using the quantum Fisher information, multipartite entanglement in the no-click trajectory shows the same phase boundaries as the entropy. With quantum jumps the transition remains, but a possible third phase appears with logarithmic entropy and bounded multipartiteness ([arXiv](https://www.arxiv.org/abs/2302.06477); [SISSA IRIS](https://iris.sissa.it/handle/20.500.11767/139292)).
- **Di Fresco, Spagnolo, Valenti, Carollo, arXiv:2302.10132; Quantum 8, 1326 (2024) [V].** In the no-click Ising chain, the QFI signals the criticality through non-analytic behaviour ([Quantum](https://quantum-journal.org/papers/q-2024-04-30-1326/)).
- **Tirrito, Santini, Fazio, Collura, arXiv:2212.09405 (Dec 2022); SciPost Phys. 15, 096 (2023), published 14 Sep 2023 [V].**
  - Setup: an Ising chain with random projective measurements of the transverse magnetization, starting from a GHZ state.
  - Instead of entanglement they use the full counting statistics of a local order parameter (the magnetization) over trajectories. The distributions distinguish paramagnetic from ferromagnetic order relaxation as a function of measurement rate.
  - In the referee exchange, the authors clarify that this is FCS over trajectories, not standard quantum FCS ([SciPost](https://www.scipost.org/SciPostPhys.15.3.096); [arXiv](https://arxiv.org/pdf/2212.09405)).
- **Zerba & Silva, arXiv:2301.07383 (Jan 2023); SciPost Phys. Core 6, 051 (2023), accepted 8 Jun 2023 [V].** For zero-click trajectories of the TFIM under local transverse-field measurements, and for a long-range Kitaev chain, the stationary-state entanglement transitions have the same nature as the transition of the non-Hermitian Hamiltonian's vacuum. They suggest generalizing the area-law theorem to non-Hermitian Hamiltonians ([SciPost](https://scipost.org/SciPostPhysCore.6.3.051); [arXiv](https://arxiv.org/abs/2301.07383v1)).
- **Le Gal, Turkeshi, Schirò, PRX Quantum 5, 030329 (2024), arXiv:2312.13419 [V].** They report a sub-volume logarithmic phase in the weak-monitoring regime of the Ising chain. Jumps strongly renormalize the effective non-Hermitian dynamics ([DOAJ](https://doaj.org/article/a1633e763ee04d76a53c6fb7c950f3c8); [arXiv](https://arxiv.org/abs/2312.13419v3)).
- **Paviglianiti, Turkeshi, Schirò, Silva, arXiv:2310.02686 (Oct 2023); Quantum 8, 1576 (2024) [V].** For measurement-*altered* Ising ground states rather than dynamics, forced on-site measurements can *increase* entanglement ([Quantum](https://quantum-journal.org/papers/q-2024-12-23-1576/)).
- Fava et al. (PRX 2023) present (ln L)^2 as "in contrast to findings in previously-studied models" ([PRX](https://link.aps.org/doi/10.1103/PhysRevX.13.041045)).
- **Guo, Chen, Wei, "Measurement-enhanced entanglement in a monitored superconducting chain", arXiv:2604.04375 (Apr 2026) [V-arXiv; no journal found].**
  - Using an NLSM calculation and free-fermion simulations, they find that for pairing Δ > 0 and small but finite γ, the steady-state entanglement is S ~ ln² L.
  - The "measurement-enhanced entanglement" therefore does *not* persist in the thermodynamic limit. The effective theory is an SO(R) NLSM for the remaining gapless modes ([arXiv](https://arxiv.org/abs/2604.04375)).
  - The companion derivation is Guo & Wei, "Super-logarithmic entanglement scaling in a monitored superconducting chain", arXiv:2607.07835 (Jul 2026) [V-arXiv] ([arXiv](https://arxiv.org/abs/2607.07835)).
- arXiv:2611.16912 (Nov 2026), "Symmetry-Controlled Criticality and Zeno Localization in a Monitored Quasiperiodic Superconductor" [title only; content and authors unverified] ([arXiv](https://arxiv.org/pdf/2611.16912)).

### Inferences
- **Current understanding (Nov 2026).**
  - For Born-rule continuously monitored Majorana chains without U(1), such as the Kitaev/BCS chain with Δ ≠ 0 and very likely the TFIM, the theory predicts a genuine transition between a (ln L)^2 phase and an area-law phase (Fava 2023; Guo–Chen–Wei 2026).
  - The 2021–2022 "log-law phase with c_eff" (Turkeshi 2021/2022; Piccitto 2022) is best read as the intermediate-L appearance of the (ln L)^2 phase. The square-root vanishing of c_eff comes from a phenomenological quasiparticle-resetting theory, not a controlled RG.
  - The Δ = 0 (U(1)) limit is different: it is area law asymptotically.
- **Unraveling dependence.** QSD vs QJ (Piccitto) is expected to change non-universal numbers and finite-size behaviour. The NLSM picture suggests the asymptotic phase structure is shared by Born-rule unravelings without U(1). This is my inference and it is not tested at L ~ 10^4 for the TFIM.
- **Measurement-only TFIM.** It has a genuine transition (percolation at p = 1/2), but it separates two area-law phases. There is no critical phase.

### Gaps
- I found no 2024–2026 study testing ln L vs (ln L)^2 at very large L (~10^4) *specifically* for the single-flavor monitored TFIM.
- I found no Born-vs-forced comparison for the TFIM dynamics, and no verified statement on whether the TFIM's nearest-neighbour Majorana structure adds a symmetry beyond class DIII.
- Precise critical rates and exponents for the TFIM under the Born rule are not available from the sources I could reach.

## 5. Non-Hermitian / no-click and jump-resolved dynamics

### Takeaway
The no-click trajectory is deterministic quadratic non-Hermitian evolution, and it is *exactly* solvable. Thermodynamic-limit quasiparticle pictures and translation-invariant circuits with exact ν show volume-to-area and log-to-area transitions governed by the imaginary (decay) spectrum of the quasiparticles, not by PT breaking. When jumps and Born weights are restored, these results are not generally robust:
- jumps strongly renormalize the dynamics;
- volume-law phases are not expected to survive;
- some area-law transitions survive with modified c_eff;
- partial postselection destroys the no-click critical phase.

### Cited Findings
- **Le Gal, Turkeshi, Schirò, arXiv:2210.11937 (Oct 2022); SciPost Phys. 14, 138 (2023), published 31 May 2023 [V].**
  - Model: a non-Hermitian SSH chain, the no-click limit of a continuously monitored free-fermion chain.
  - Method: entanglement dynamics computed in the thermodynamic limit.
  - There is a volume-to-area-law transition. It does not coincide with PT-symmetry breaking; it happens once the *full decay spectrum of the quasiparticles becomes gapped* ([SciPost](https://www.scipost.org/SciPostPhys.14.5.138); [arXiv](https://arxiv.org/abs/2210.11937v2)).
- **Turkeshi & Schirò, arXiv:2201.09895 (Jan 2022); Phys. Rev. B 107, L020403 (2023) [V].**
  - Setup: a quench in the anisotropic XY chain with a complex transverse field, solved exactly in the thermodynamic limit via Jordan–Wigner.
  - There are two dynamical phases. In one, entanglement and correlations saturate fast. In the other, entanglement grows logarithmically and correlations spread *faster* than in the unitary case, with collapses and revivals that form a modulated light cone ([arXiv](https://arxiv.org/abs/2201.09895v1); [IP Paris record](https://researchportal.ip-paris.fr/en/publications/entanglement-and-correlation-spreading-in-non-hermitian-spin-chai/)).
- Zerba & Silva (2023) identify the no-click stationary state with the non-Hermitian vacuum, and the transitions have the same nature ([SciPost](https://scipost.org/SciPostPhysCore.6.3.051)).
- Granet–Zhang–Dreyer (PRL 130, 230401, 2023) give an exact ν for postselected periodic Gaussian circuits, with volume-to-area and log-to-area transitions ([arXiv](https://arxiv.org/abs/2212.10584)).
- arXiv:2507.03768 (Jul 2025), "Entanglement transitions in structured and random nonunitary Gaussian circuits" [authors and journal unverified], recovers the periodic-circuit transition algebraically. The associated Möbius transformations are elliptic or parabolic for SU(2), so only volume or sub-volume-law phases arise in the compact case ([arXiv](https://arxiv.org/pdf/2507.03768)).
- **Kawabata, Numasawa, Ryu, arXiv:2206.05384 (Jun 2022); Phys. Rev. X 13, 021007 (2023) [V].**
  - The non-Hermitian skin effect suppresses entanglement propagation, giving an area-law steady state.
  - Competition between unitary dynamics and the skin effect induces an entanglement transition with no disorder, interactions or randomness ([NSF-PAR](https://par.nsf.gov/servlets/purl/10429292); [arXiv](https://arxiv.org/pdf/2206.05384)).
  - Follow-up: "Disorder-Induced Entanglement Phase Transitions in Non-Hermitian Systems with Skin Effects", Quantum, published 7 Aug 2026 ([Quantum](https://quantum-journal.org/papers/q-2026-08-07-2186/)) [arXiv ID and authors unverified].
- Turkeshi–Dalmonte–Fazio–Schirò (PRB 2022) model quantum jumps as stochastic resetting of non-Hermitian quasiparticles. This phenomenological theory matches exact numerics ([arXiv](https://arxiv.org/abs/2111.03500)).
- **Le Gal, Turkeshi, Schirò, "Entanglement dynamics in monitored systems and the role of quantum jumps", arXiv:2312.13419 (Dec 2023; v3 Jun 2024); PRX Quantum 5, 030329 (Aug 2024) [V].**
  - New diagnostic: the statistics of entanglement changes *after and between* jumps, studied in the Ising chain, the monitored SSH chain and a U(1) model.
  - Jumps strongly renormalize the effective non-Hermitian dynamics, giving large deviations from the no-click limit.
  - There is a sub-volume logarithmic phase in the weak-monitoring Ising regime ([arXiv](https://arxiv.org/abs/2312.13419v3); [DOAJ](https://doaj.org/article/a1633e763ee04d76a53c6fb7c950f3c8)).
- **Soares, Le Gal, Schirò, "Entanglement transition due to particle losses in a monitored fermionic chain", arXiv:2408.03700 (Aug 2024) [V-arXiv; journal unverified].**
  - The chain has hopping and pairing, so no U(1).
  - A loss-driven transition appears both in QJ trajectories and in the no-click limit.
  - The QJ effective central charge is *higher* than in the no-click limit, so jumps can increase entanglement.
  - The no-click transition is tied to the opening of a gap in the imaginary part of the quasiparticle spectrum at γ/J = 2 ([arXiv](https://arxiv.org/pdf/2408.03700)).
- arXiv:2407.13837 (2024) [unverified details] finds that partial postselection destroys the no-click critical phase in monitored Kitaev chains ([arXiv](https://arxiv.org/html/2407.13837v2)).
- With jumps, Paviglianiti–Silva find the transition persists together with a possible extra phase ([arXiv](https://www.arxiv.org/abs/2302.06477)). Piccitto et al. report that the non-Hermitian approximation conflicts with the full unravelings ([arXiv](https://arxiv.org/abs/2111.11300v3)).
- arXiv:2503.02418 (Mar 2025), "Measurement-Induced Crossover of Quantum Jump Statistics in Postselection-Free Many-Body Dynamics" [title only] ([arXiv](https://arxiv.org/pdf/2503.02418)).

### Inferences
- **What is exact.**
  - Exact: no-click dynamics of quadratic non-Hermitian Hamiltonians, i.e. momentum-decoupled quasiparticles with complex energies in the thermodynamic limit (Le Gal 2023; Turkeshi–Schirò 2023; Zerba–Silva 2023); and periodic postselected circuits, where ν is exact (Granet–Zhang–Dreyer 2023).
  - Not exact: the quasiparticle-resetting theory for jumps is phenomenological, though benchmarked numerically.
- **Mechanism.** In no-click dynamics the transitions are spectral: entanglement is carried by quasiparticles with the smallest decay rates. An imaginary gap kills long-range quasiparticle pairs and gives an area law. Volume law needs a finite density of non-decaying modes.
- **Survival with jumps.**
  - Born-rule Gaussian dynamics without U(1) is bounded asymptotically by (ln L)^2 (§1). So no-click *volume-law* phases cannot survive once jumps are restored (inference). This is consistent with the "large deviations" and "sub-volume log phase" in Le Gal 2024.
  - *Area-law* transitions tied to an imaginary gap can survive with jumps, with modified c_eff or phase boundaries (Soares 2024; Paviglianiti–Silva 2023).
  - The no-click critical phase does not survive partial postselection (2407.13837).
- **Toeplitz / Szegő.** For translation-invariant Gaussian states, block entanglement reduces to (block-)Toeplitz determinants. Szegő/Widom asymptotics then naturally give area, volume or log behaviour (Fisher–Hartwig for log) from the smoothness and winding of the symbol. This is a methodological inference; see Gaps.

### Gaps
- I could not verify that any of these papers explicitly uses Szegő or Toeplitz asymptotics. The searches returned "exact in the thermodynamic limit" and "quasiparticle" wording only.
- The search budget ran out before I could check:
  - Jian–Yang–Bi–Chen, "Yang–Lee edge singularity triggered entanglement transition";
  - Gopalakrishnan–Gullans on non-Hermitian purification;
  - Bao–Choi–Altman, "Symmetry enriched phases of quantum circuits".
- Whether Kawabata et al.'s models conserve U(1) was not verified; I believe they are Hatano–Nelson-type hopping models.
- I found no study of whether no-click volume-law phases survive at small jump rates at very large L.

## 6. 2024–2026 developments (field theory, numerics, experiments)

### Takeaway
After 2023 the Majorana (DIII) picture consolidated:
- loop-model and symmetry classifications extended to Kitaev-type circuits (PRB 2025);
- weak-measurement (detector-coupling) control of universality (Quantum 2025);
- interacting class-DIII extensions conjecturing the same critical point (PRB, 2025/26);
- monitored superconducting chains confirming ln² L, which overturns "measurement-enhanced entanglement" (2026 preprints);
- large-scale GPU numerics contrasting Majorana and U(1) (PRB 2025).

I found no experiment observing free-fermion (Majorana) monitored criticality or loop-model entanglement on hardware.

### Cited Findings
- **Klocke et al., Phys. Rev. B 111, 224301 (2025), arXiv:2409.02171 [V].** Covers classes BDI and D and Majorana loop models in any dimension, with a quantum Lifshitz transition in two universality classes ([arXiv](https://www.arxiv.org/abs/2409.02171)).
- **Nehra–Romito–Meidan, Quantum 9, 1697 (2025) [V].** Percolation is unstable to weaker detector coupling ([Quantum](https://quantum-journal.org/papers/q-2025-04-08-1697/)).
- **"Free-Fermion Measurement-Induced Volume- to Area-Law Entanglement Transition in the Presence of Fermion Interactions", arXiv:2510.23706 (Oct 2025); Phys. Rev. B, DOI 10.1103/874y-tk4t, an Editors' Suggestion. Only the first author's given names ("Matthew S.") were visible; full author list and volume unverified.**
  - Scope: 1D monitored *interacting* Majoranas in class DIII.
  - They conjecture that the MIPT critical point of the weakly interacting class-DIII circuit is identical to the noninteracting one.
  - The interaction-induced volume-law phase arises through a dangerously irrelevant "mass" tied to the arrest of quantum chaos.
  - The argument does not apply when there are extra conservation laws such as U(1) ([arXiv](https://arxiv.org/html/2510.23706v1); [APS](https://journals.aps.org/prb/abstract/10.1103/874y-tk4t)).
- **Guo–Chen–Wei, arXiv:2604.04375 (Apr 2026), and Guo–Wei, arXiv:2607.07835 (Jul 2026) [V-arXiv].** In monitored BCS/Kitaev chains, S ~ ln² L at Δ > 0 and small γ, from an SO(R) NLSM; "measurement-enhanced entanglement" does not survive in the thermodynamic limit ([arXiv](https://arxiv.org/abs/2604.04375); [arXiv](https://arxiv.org/abs/2607.07835)).
- **GPU numerics, arXiv:2508.18468 (Aug 2025); Phys. Rev. B, DOI 10.1103/x7vx-tqr4 [V-arXiv + DOI; volume unverified].**
  - Sizes: L = 16384 in 1D and 160×160 in 2D.
  - Confirming the NLSM-predicted *absence* of a MIPT for U(1) fermions, under projective and homodyne measurements, requires L ~ 10^4. A MIPT at finite rate for Majorana parity monitoring is the stated contrast ([arXiv](https://arxiv.org/html/2508.18468); [APS](https://journals.aps.org/prb/abstract/10.1103/x7vx-tqr4)).
  - Follow-up arXiv:2605.10758 (May 2026) treats U(1) fermions with disordered or quasiperiodic potentials up to L ≤ 18000. The critical monitoring strength is consistent with zero, and it reiterates that 1D Majorana parity monitoring has a finite-rate MIPT ([arXiv](https://arxiv.org/html/2605.10758v1)).
- Decoherence context, arXiv:2502.14034 (2025): in class D, the N→0 limit hosts a Majorana metal that is absent at N→1 ([arXiv](https://arxiv.org/html/2502.14034v2)).
- **Topology in monitored Majorana circuits** [authors and journals unverified]:
  - "Topology and Spectrum in Measurement-Induced Phase Transitions", arXiv:2412.11097 (Dec 2024; v4), attributed to Fuji and collaborators, proposes a bulk topological invariant and gapless edge modes for MIPTs in monitored Majorana circuits ([arXiv](https://arxiv.org/html/2412.11097v4)).
  - "Topological Modes in Monitored Quantum Dynamics", arXiv:2411.04191 (Nov 2024), reports a critical window with logarithmic entanglement between two area-law phases in a noninteracting model ([arXiv](https://arxiv.org/html/2411.04191v2)).
- Non-Hermitian and jumps (2024–2026): Le Gal–Turkeshi–Schirò (PRX Quantum 2024), Soares–Le Gal–Schirò (2024), partial postselection (2024), structured/random non-unitary Gaussian circuits (2025), and the skin-effect disorder paper (Quantum, Aug 2026) are all in §5.
- **Majorana defects, 2026 conference abstract.** "Hybrid dynamics of diffusing Majorana defects: amplitudes, parity regimes, and braiding-induced connectivity" (APS March Meeting 2026, Skinner co-author). The late-time decay stays in the A+A→∅ class, but the amplitude depends on topology ([APS abstract](https://meetings-archive.aps.org/smt/2026/mar-m37/7/)).
- **Experiments.** I found no experimental paper observing loop-model or (ln L)^2 entanglement in a monitored or measurement-only Majorana free-fermion circuit. Adjacent hardware work:
  - Majorana *Hamiltonian* dynamics on an IBM processor (arXiv:2510.20034) ([arXiv](https://arxiv.org/pdf/2510.20034));
  - braiding Majorana modes on programmable simulators (arXiv:2203.15083) ([arXiv](https://arxiv.org/pdf/2203.15083));
  - an experimental-detection *protocol* for the projective TFIM transition (arXiv:2511.17370, 2025) ([arXiv](https://arxiv.org/pdf/2511.17370)).
- Seminar activity: "Nonlinear sigma model description of monitored free Majorana fermions", UZH FKP seminar, spring 2024 ([UZH](https://www.physik.uzh.ch/en/seminars/fkpseminar/FS2024/Nonlinear-sigma-model-description-of-monitored-free-Majorana-fermions.html)).
- **Title-only 2025–2026 preprints** surfaced by searches; relevance not checked, so do not cite for content:
  - arXiv:2501.03851, "Gapless SPT states in measurement-only circuits" ([link](https://arxiv.org/pdf/2501.03851));
  - arXiv:2509.15634, "Measurement-driven transitions between area law phases" ([link](https://arxiv.org/html/2509.15634v1));
  - arXiv:2509.08584, "Spectral transitions of the entanglement Hamiltonian in monitored free fermions" ([link](https://arxiv.org/html/2509.08584v1));
  - arXiv:2507.07179, "Magic phase transitions in monitored Gaussian fermions" ([link](https://arxiv.org/pdf/2507.07179));
  - arXiv:2509.16758, "Error stabilized logical qubits in qudit generalizations of the monitored Kitaev model" ([link](https://arxiv.org/html/2509.16758v1));
  - arXiv:2503.21427, "Measurement-induced entanglement phase transition in free fermion systems" ([link](https://arxiv.org/html/2503.21427v3));
  - arXiv:2604.20828, "Arrow of time as an indicator of MIPTs" ([link](https://arxiv.org/pdf/2604.20828)).

### Inferences
- **Consensus as of Nov 2026, from the sources reached.**
  - The DIII (Majorana, no U(1)) case is the robust example of a genuine 1D free-fermion MIPT, between an asymptotically (ln L)^2 phase and an area-law phase.
  - Theory, the 2026 pairing-chain work, and repeated statements in U(1)-focused numerics papers all support this.
  - Interactions are conjectured not to change the DIII critical point (2510.23706).
- **Open directions:**
  - analytic Born-rule (N→1) critical exponents;
  - very-large-L numerics on the single-flavor TFIM;
  - experiments, of which none were found.

### Gaps
- The final web-search budget was exhausted, so several 2024–2026 items are title-only and some journal references remain unverified (see the ledger).
- A 2025–2026 review dedicated to monitored free fermions without U(1) was not located.

### Reference ledger (verification status)

| Reference | arXiv | Journal ref | Status |
|---|---|---|---|
| Jian, Bauer, Keselman, Ludwig | 2012.04666 | PRB 106, 134206 (2022) | V |
| Fava, Piroli, Swann, Bernard, Nahum | 2302.12820 | PRX 13, 041045 (2023) | V |
| Jian, Shapourian, Bauer, Ludwig | 2302.09094 | none found | V-arXiv |
| Nahum, Skinner | 1911.11169 | Phys. Rev. Research (2020); "2, 023288" NOT verified | partial |
| Klocke, Buchhold | 2305.18559 | PRX 13, 041028 (2023) | V |
| Klocke, Simm, Zhu, Trebst, Buchhold | 2409.02171 | PRB 111, 224301 (2025) | V (search-rendered) |
| Nahum, Serna, Somoza, Ortuño | 1303.2342 | PRB 87, 184204 (2013) | V |
| Jacobsen, Read, Saleur | cond-mat/0205033 | not verified | V-arXiv |
| Sang, Hsieh | 2004.09509 | PRR 3, 023200 (2021) | V |
| Lavasani, Alavirad, Barkeshli | 2004.07243 | Nat. Phys. 17, 342 (2021) | V |
| Merritt, Fidkowski | 2210.05681 | PRB 107, 064303 (2023) | V |
| Kells, Meidan, Romito | 2112.09787 | SciPost Phys. 14, 031 (2023) | V |
| Nehra, Romito, Meidan | 2404.07918 | Quantum 9, 1697 (2025) | V |
| Lang, Büchler | 2006.09748 | not verified | V-arXiv |
| Turkeshi, Biella, Fazio, Dalmonte, Schirò | 2103.09138 | PRB 103, 224210 (2021) | V |
| Piccitto, Russomanno, Rossini | 2111.11300 | PRB 105, 064305 (2022) | V |
| Turkeshi, Dalmonte, Fazio, Schirò | 2111.03500 | PRB 105, L241114 (2022); Erratum PRB 107, 079901 (2023) | V |
| Paviglianiti, Silva | 2302.06477 | PRB 108, 184302 (2023) | V |
| Di Fresco, Spagnolo, Valenti, Carollo | 2302.10132 | Quantum 8, 1326 (2024) | V |
| Paviglianiti, Turkeshi, Schirò, Silva | 2310.02686 | Quantum 8, 1576 (2024) | V |
| Tirrito, Santini, Fazio, Collura | 2212.09405 | SciPost Phys. 15, 096 (2023) | V |
| Zerba, Silva | 2301.07383 | SciPost Phys. Core 6, 051 (2023) | V |
| Le Gal, Turkeshi, Schirò (non-Hermitian SSH) | 2210.11937 | SciPost Phys. 14, 138 (2023) | V |
| Turkeshi, Schirò | 2201.09895 | PRB 107, L020403 (2023) | V |
| Le Gal, Turkeshi, Schirò (quantum jumps) | 2312.13419 | PRX Quantum 5, 030329 (2024) | V |
| Kawabata, Numasawa, Ryu | 2206.05384 | PRX 13, 021007 (2023) | V |
| Granet, Zhang, Dreyer | 2212.10584 | PRL 130, 230401 (2023) | V |
| Soares, Le Gal, Schirò | 2408.03700 | not verified | V-arXiv |
| Guo, Chen, Wei / Guo, Wei | 2604.04375 / 2607.07835 | none found | V-arXiv |
| (first author "Matthew S."), interacting DIII | 2510.23706 | PRB, DOI 10.1103/874y-tk4t (volume unverified) | partial |
| GPU monitored fermions | 2508.18468 | PRB, DOI 10.1103/x7vx-tqr4 (volume unverified) | partial |
| Fava, Piroli, Bernard, Nahum (U(1), context only) | 2407.08045 | not verified | V-arXiv |
| Poboiko, Pöpperl, Gornyi, Mirlin (U(1), context only) | 2304.03138 | not verified in this session | V-arXiv |
| Yang, Mao, Jian | 2301.08255 | not verified | V-arXiv |
| Decoherence self-dual criticality review | 2502.14034 | not verified | V-arXiv |
| Partially postselected Kitaev chains | 2407.13837 | not verified; authors unverified | partial |
| Skin effect plus disorder | unknown | Quantum (7 Aug 2026), q-2026-08-07-2186 | partial |
