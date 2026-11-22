# AI-driven and automated discovery of quantum error-correcting codes (code constructions and encodings), 2023 – Nov 2026

**Read first: verification status.**
- WebFetch could not reach arxiv.org, huggingface.co or pasqal.com from this session, so no full texts were read. Every finding below comes from WebSearch summaries of primary pages: arXiv abstract/HTML listings, journal pages, GitHub repositories and project blogs. The numbers are therefore abstract-level and should be spot-checked against the PDFs. arXiv IDs are given as search returned them. Some 2026 IDs have unusually high 5-digit suffixes (e.g. 2609.37056, 2610.41183), and none were cross-checked on arxiv.org.
- The session-wide WebSearch budget (200 calls per turn, shared by all agents) ran out partway through this assignment. As a result, Q3 (bosonic/autonomous QEC), most of Q4 (the fermion-to-qubit mapping literature), Q5 (subsystem/Floquet) and the pre-2023 seminal works could NOT be verified in this session. They appear only in the "Gaps" subsections, marked **[UNVERIFIED lead]**, and must not be cited as verified.
- Corrections to the brief (verified):
  - The HGP-optimization paper (arXiv:2501.09622) is by Bruno C. A. Freire, Nicolas Delfosse and Anthony Leverrier. No "Le Régent" appears on it ([arXiv](https://arxiv.org/abs/2501.09622); [EC Zoo HGP page](https://errorcorrectionzoo.org/c/hypergraph_product)).
  - "He & Liu" are Austin Yubo He and Zi-Wen Liu ([arXiv:2502.14372](https://arxiv.org/abs/2502.14372v1)).
  - Su et al. appeared as "Discovery of optimal quantum codes via reinforcement learning", Phys. Rev. Applied 23, 034048 (2025) ([APS](https://link.aps.org/doi/10.1103/PhysRevApplied.23.034048)).
- Per a coordinator scope update, these notes put the newest works first (2025 to Nov 2026), with emphasis on LLM-agent, AlphaEvolve/OpenEvolve/FunSearch-style and other automated search. "Gaps and opportunities" covers the literature in general and is not tailored to any group.

**Newest-first index of the 2025–2026 works found (preprints unless a venue is given)**

| Date | Work | Search method | Target | Headline claim (as reported by the authors) | Verifier / score |
|---|---|---|---|---|---|
| 2026-11 | Jain, Raveendran, Albarghouthi, QECLean, [arXiv:2611.21236](https://arxiv.org/abs/2611.21236) | Lean formalization, with AI agents used to formalize BB codes | verification infrastructure | formally verifies properties of prominent and recent codes | Lean |
| 2026-11 | Cheng, Guo, Pan, tropidist, [arXiv:2611.08247](https://arxiv.org/abs/2611.08247) | exact distance certification | verification infrastructure | audit of 26 papers: 10 overstated at least one circuit- or code-distance estimate | tropical/rank-based DP with certified bounds |
| 2026-11 | Batista, Falcão, Machado, GFlowQEC, [arXiv:2611.06819](https://arxiv.org/abs/2611.06819) | GFlowNet | lifted-product, n=180 | >2× more new top codes per paid evaluation than tuned parallel tempering (top-10 ranking) | circuit-level simulation scores |
| 2026-11 | Du et al., OmniQ, [arXiv:2611.02833](https://arxiv.org/abs/2611.02833) | reproduction benchmark for AI scientists | 60 computational-quantum papers | AI scientists often fail to reproduce published results | expert-audited tasks |
| 2026-10 | Baek, Won, Kim, [arXiv:2610.41183](https://arxiv.org/abs/2610.41183) | AlphaEvolve, researcher in the loop, Lean | entanglement-assisted (EA) codes | an infinite family that closes 54 gaps in the EAQECC tables | exact (k,c,d) evaluator, nonexistence certificates, Lean |
| 2026-10 | LiftDiff, [arXiv:2610.25213](https://arxiv.org/pdf/2610.25213) | guided discrete diffusion | lifted-product | first discrete-diffusion generator of QECCs | not retrieved |
| 2026-10 | Schneider, Nachmani, BB-Net, [arXiv:2610.18130](https://arxiv.org/abs/2610.18130) | learned distance predictor | BB codes | 96.39% vs 94.86% for heuristics; a third party could not reproduce this | heuristic upper bounds, certified subset |
| 2026-10 | Cruz-Benito, [arXiv:2610.06623](https://arxiv.org/abs/2610.06623) | LLM program evolution | weight-5 BB/PBB | 1,142 proposals; 816 of 936 CSS proposals certified d≥5 | certified distance, random controls |
| 2026-10 | An, Yang, [arXiv:2610.03214](https://arxiv.org/abs/2610.03214) | Lean-checked certificates | 11 families, up to 1872 qubits | certifies the [[144,12,12]] distance | Lean kernel |
| 2026-09-30 | Unitary Foundation QEC Challenge, [blog](https://unitary.foundation/posts/2026_qec_challenge/) | open leaderboard for human and AI submissions | qLDPC/stabilizer | about 2,610 verified codes, mostly machine-found | automated verifier, Pareto frontier per class |
| 2026-09 | Seddas, Shashkov, Viazovska, Abbe, LinCodeEvolve, [arXiv:2609.37056](https://arxiv.org/abs/2609.37056v1) | LLM program evolution (EvoTune/ShinkaEvolve) | classical binary linear codes | 7 new record codes including [200,21,77]; 22 table entries raised | exact minimum distance, exhaustive enumeration, review by the table maintainer |
| 2026-09 | Lu, Yang, Guo, [arXiv:2609.06572](https://arxiv.org/abs/2609.06572) | method not determined | weight-8 BB | [[144,6,d≥15]]; census of 53 codes at n=144 | BP-OSD, circuit-level |
| 2026-08 | Sella, Pinsler, Heinis, [arXiv:2608.15754](https://arxiv.org/abs/2608.15754) | PPO reinforcement learning | BB codes with a fixed decoder | results not retrieved | area under the logical-vs-physical error curve |
| 2026-08 | Qian, Li, [arXiv:2608.08996](https://arxiv.org/abs/2608.08996) | multi-agent LLM | coset-orbit balanced-product codes, n≤400, w≤10 | [[288,16,18]] at w=7, [[288,18,18]] at w=9, [[234,28,18]] at w=10 | BP-OSD, code capacity |
| 2026-07 | Ge Yan et al., OmniQEC, [arXiv:2607.25865](https://arxiv.org/abs/2607.25865v1) | LLM "AI scientist" (Claude, GPT and DeepSeek backends) | 4 qLDPC families | outperforms the [[72,12,6]] and [[144,12,12]] BB codes at budgets of 98 and 240 physical qubits | circuit-level LER |
| 2026-06 | Heitritter, Brown, Hardikar (qBraid), [arXiv:2606.25870](https://arxiv.org/abs/2606.25870) | AlphaEvolve plus an open-source replication | fermion-to-qubit GSE encodings | first GSE encodings beyond d=3 for dense molecular Hamiltonians (d=5; d=6 on one instance) | fixed, exact verifier |
| 2026-06 | Liu, Marquardt, SCE, [arXiv:2606.24808](https://arxiv.org/abs/2606.24808) | LLM with an algebraic mutation grammar | lifted-product | beats published LP benchmarks in some regimes | kd²/n or circuit-level LER; exact distances for the top codes |
| 2026-06 | Cruz-Benito, Cross, Kremer, Faro (IBM), [arXiv:2606.02418](https://arxiv.org/abs/2606.02418) | OpenEvolve (MAP-Elites) | BB and non-CSS PBB codes, n≤360 | 465 distinct codes, including an indecomposable [[288,16,12]] | staged checks: parameters, distance, equivalence class |
| 2026-05 | Lean-QEC, [arXiv:2605.16523](https://arxiv.org/html/2605.16523) | formalization | verification infrastructure | notes that certifying a distance lower bound is NP-hard in general | Lean |
| 2026-01 | Y. Chengyu, Meister, Carty, Lin, Bondesan, [arXiv:2601.18562](https://arxiv.org/html/2601.18562) | Bayesian optimization (Gaussian-process surrogate) | qLDPC, including HGP | [[144,36]] (3× the gross code's k) and [[144,26]] | Monte Carlo LER at code capacity |
| 2026-01 | Yanay, [arXiv:2601.08014](https://arxiv.org/abs/2601.08014) | RL with quantum hardware in the loop (Quantum Lego) | small stabilizer codes | qubit error reduced by 97% and 85% on two devices | hardware evaluation |
| 2025-06 | [arXiv:2506.11552](https://arxiv.org/pdf/2506.11552) | variational QEC (learned encodings) | encodings | not retrieved | not retrieved |
| 2025-04 | [arXiv:2504.00613](https://arxiv.org/abs/2504.00613) | FunSearch adaptation | classical deletion codes | rediscovered the VT code (with proof); larger short codes for multiple deletions | exact |
| 2025-02 | He, Liu, [arXiv:2502.14372](https://arxiv.org/abs/2502.14372v1) | RL weight reduction | stabilizer/qLDPC | d up to about 35; 1–2 orders of magnitude less overhead at weight 6 | distance constraint |
| 2025-01 | Freire, Delfosse, Leverrier, [arXiv:2501.09622](https://arxiv.org/abs/2501.09622) (ISIT 2025) | random walk, simulated annealing, RL | HGP codes, erasure channel | beats PEG; the optimized codes have girth 4 | erasure-channel performance |

---

## Q1. Qubit stabilizer codes: what did RL and search achieve for code discovery and noise-adapted codes, and what newer 2025–2026 works exist?

### Takeaway
All five RL works named in the brief exist (with the author correction above), as does the HGP paper. RL built from scratch over small stabilizer codes reached n of about 13–25 with d≤5. It mostly matched known optima, for example saturating the CSS linear-programming bound at n=13. The clearly new results are noise-tailored codes (biased noise, device noise, erasure). Since 2025 the work has moved to RL that modifies large existing codes (weight reduction up to d≈35), to RL with hardware in the loop (2026), and to RL co-design of codes and decoders within the BB family (2026). I found no paper in which AlphaZero or MCTS discovers codes.

### Cited Findings (newest first)
- **Sella, Pinsler, Heinis (Aug 2026), "Designing Quantum Error Correcting Codes to fit decoders via Reinforcement Learning", arXiv:2608.15754.**
  - *Method:* the decoder architecture is fixed and a PPO agent is trained over bivariate-bicycle (BB) codes to maximize decoder performance under a depolarizing model.
  - *Reward:* the area under the fitted physical-to-logical error-rate curve, a single scalar. A "warm start" uses a predictive model.
  - *Scale, results and limitations:* not retrievable from the snippets.
  - Sources: [arXiv](https://arxiv.org/abs/2608.15754); [HTML](https://arxiv.org/html/2608.15754v1)
- **Yanay (Jan 2026), "Learning Better Error Correction Codes with Hybrid Quantum-Assisted Machine Learning", arXiv:2601.08014.** The arXiv listing names Yariv Yanay; news coverage says "Yanay and colleagues".
  - *Method:* classical RL over Quantum-Lego constructions. Candidate codes are evaluated on two commercial quantum devices, tailoring them to device-specific errors and to induced photon loss.
  - *Result:* the learned codes reduce qubit error by 97% and 85% in the two cases. In both, the best correction uses 7-qubit stabilizers.
  - *Limitation (as summarized):* the devices do not run error correction themselves; the gains come from evaluating candidate codes on hardware.
  - Sources: [arXiv](https://arxiv.org/abs/2601.08014); [Quantum Zeitgeist](https://quantumzeitgeist.com/quantum-machine-learning-error-correction-codes-systems-hybrid-assisted-achieves/)
- **He & Liu (Feb 2025), "Discovering highly efficient low-weight quantum error-correcting codes with reinforcement learning", arXiv:2502.14372** (quant-ph, cross-listed in cs.AI, cs.IT and cs.LG).
  - *Method:* RL for stabilizer-code weight reduction.
  - *Result:* for weight-6 codes, physical-qubit overhead is 1–2 orders of magnitude lower than in existing results. Distances reach d≈35. The authors describe the new low-weight codes as substantially outperforming the state of the art in practically relevant regimes.
  - *Stated insight:* "decreasing weight with distance constraints is a significantly more approachable problem compared to increasing distance with weight constraints, especially for learning methods."
  - Sources: [arXiv](https://arxiv.org/abs/2502.14372v1)
- **Freire, Delfosse, Leverrier (Jan 2025, v3 Jun 2025; accepted at IEEE ISIT 2025), "Optimizing hypergraph product codes with random walks, simulated annealing and reinforcement learning", arXiv:2501.09622.**
  - *Methods:* random-walk exploration, simulated annealing and RL. The dimension and distance of an HGP code depend only on its classical seed codes, but the performance against a given channel does not, so performance is what is optimized.
  - *Objective:* performance on the quantum erasure channel. The baseline is Progressive Edge-Growth (PEG).
  - *Result:* simulated annealing and RL find HGP codes that outperform PEG codes. All PEG codes have girth 6 and all optimized codes have girth 4, which the authors read as "maximizing girth is not a good heuristic."
  - Sources: [arXiv](https://arxiv.org/abs/2501.09622); [HTML v1](https://arxiv.org/html/2501.09622v1)
- **Su, Cao, Hu, Yanay, Tahan, Swingle (arXiv May 2023; Phys. Rev. Applied 23, 034048, 2025), arXiv:2305.06378.**
  - *Method:* the Quantum-Lego construction is turned into a game with three actions (add a lego, contract legs, end the game), and an RL agent maximizes an objective.
  - *Results:*
    - The agent increases distance beyond naive concatenation and saturates the linear-programming bound for CSS codes on 13 qubits.
    - For biased Pauli noise it finds the best known CSS code for the task at ≲20 qubits.
    - It finds a [[17,1,3]] code that has lower adversarial distance but protects better than the surface, XZZX and 2D color codes.
  - Sources: [arXiv](https://arxiv.org/abs/2305.06378); [Phys. Rev. Applied](https://link.aps.org/doi/10.1103/PhysRevApplied.23.034048)
- **Olle, Zen, Puviani, Marquardt (arXiv Nov 2023; npj Quantum Information 10, 126, 2024), "Simultaneous discovery of quantum error correction codes and encoders with a noise-aware reinforcement learning agent", arXiv:2311.04750.**
  - *Method:* an RL agent discovers both a code and its encoding circuit from scratch, for a given gate set, connectivity and error model.
  - *Reward:* based on the Knill–Laflamme conditions and computed with a vectorized Clifford simulator.
  - *Noise awareness:* one agent learns to switch strategy as a noise parameter varies.
  - *Scale:* the arXiv v3 abstract says 20 physical qubits and d=5. The journal abstract says up to 25 physical qubits and d=5, with a "roadmap" to 100 qubits and d=10.
  - Sources: [npj QI](https://www.nature.com/articles/s41534-024-00920-y); [arXiv v3](https://arxiv.org/abs/2311.04750v3); [APS March Meeting 2024 talk](https://meetings-archive.aps.org/mar/2024/a49/6)
- **Mauron, Farrelly, Stace (arXiv May 2023; New J. Phys. 26(2), 023024, 2024; DOI 10.1088/1367-2630/ad23a6), "Optimization of tensor network codes with reinforcement learning", arXiv:2305.11470.**
  - *Method:* an RL agent in the projective-simulation framework chooses how seed codes are contracted into a tensor-network code.
  - *Result:* it consistently finds the best possible codes for the given environment and allowed actions, including codes with k>1. After 1000 trials it finds an optimal code 10% of the time, against a theoretical 0.16% for random search (a factor of 65).
  - *Caveat:* optimality is relative to the restricted action set.
  - Sources: [arXiv](https://arxiv.org/abs/2305.11470); [DOAJ record](https://libmast.utm.my/Record/doaj.art-3b02e8090624479380cda0befa609512)
- **"Engineering Quantum Error Correction Codes Using Evolutionary Algorithms" (Sep 2024), arXiv:2409.13017** (authors not verified in this session).
  - *Method:* an evolutionary algorithm searches for an optimal stabilizer code for a given error model, n and k. Codes are represented as binary strings, so valid codes can be generated at random and then mutated and crossed over.
  - Scale and results were not retrieved.
  - Sources: [arXiv](https://arxiv.org/abs/2409.13017)
- **Variational QEC (Jun 2025), "Learning Encodings by Maximizing State Distinguishability: Variational Quantum Error Correction", arXiv:2506.11552.** Only the title was retrieved. Sources: [arXiv](https://arxiv.org/pdf/2506.11552)
- **Earlier (2022): Tandeitnik & Guerreiro, "Evolving Quantum Circuits", arXiv:2210.05058.** With suitable fitness functions, a genetic algorithm repeatedly outperformed random search. It evolved codes equivalent to Shor's 9-qubit code and the 5-qubit perfect code, so these were rediscoveries. Sources: [arXiv](https://arxiv.org/pdf/2210.05058)
- **Earlier (2019): Nautrup, Delfosse, Dunjko, Briegel, Friis, "Optimizing Quantum Error Correction Codes with Reinforcement Learning", Quantum 3, 215 (16 Dec 2019), DOI 10.22331/q-2019-12-16-215, arXiv:1812.08451.**
  - An RL framework for optimizing and fault-tolerantly adapting QEC codes, simulated on surface codes.
  - One arXiv abstract version reports near-optimal solutions, in terms of physical-qubit count, for a surface-code memory of about 70 physical qubits.
  - It demonstrates transfer learning: agents train on fast simulations and then adapt under real conditions.
  - Sources: [arXiv](https://arxiv.org/abs/1812.08451v2); [Quantum](https://quantum-journal.org/papers/q-2019-12-16-215/)
- **AlphaZero/MCTS:** I found no AlphaZero or MCTS paper that discovers new codes. MCTS does appear in two out-of-scope settings:
  - scheduling syndrome-measurement circuits for existing commuting-stabilizer codes (AlphaSyndrome, arXiv:2601.12509, ACM DOI 10.1145/3779212.3790123);
  - stabilizer-state preparation (arXiv:2603.17743).
  - Sources: [AlphaSyndrome arXiv](https://arxiv.org/pdf/2601.12509); [ACM](https://dl.acm.org/doi/10.1145/3779212.3790123); [state preparation](https://arxiv.org/pdf/2603.17743)

### Inferences
- RL from scratch over codes and encoders has not gone beyond about 25 qubits at d=5 (Olle et al.), and Quantum-Lego RL stays at ≲20 qubits. The 2025–26 gains in scale came from shrinking the action space, either to transformations of already-good codes (He & Liu) or to algebraic families (BB, Sella et al.). The design of the search space, more than the learner, seems to be what made scaling possible.
- "New" is least ambiguous for noise-adapted codes: the biased-noise [[17,1,3]] code, the device-tailored 7-qubit stabilizers and the erasure-optimized HGP codes. Standard tables index only adversarial distance, so these codes are judged on a metric the tables do not record.
- Two results show that proxy metrics can mislead channel-specific optimization: the optimized HGP codes have girth 4, below the girth-6 PEG codes they beat (Freire et al.), and Su et al.'s best biased-noise code has the lower adversarial distance.

### Gaps
- Not retrieved: Sella et al.'s results; the scale and authors of arXiv:2409.13017; which devices Yanay used; He & Liu's exact RL state/action design and their distance-verification method.
- **[UNVERIFIED lead]** Cao & Lackey's "Quantum Lego" framework (PRX Quantum, 2022), on which Su et al. and Yanay build.
- **[UNVERIFIED lead]** Earlier variational or learned encodings: QVECTOR (Johnson et al., 2017) and quantum-autoencoder QEC (Locher, Cardarelli, Müller, Quantum 2023).
- **[UNVERIFIED lead]** Zen et al., RL for fault-tolerant logical state preparation (2024). This is an FT-gadget circuit, so out of scope.
- "Scaling the Automated Discovery of Quantum Circuits via Reinforcement Learning with Gadgets" ([arXiv:2503.11638](https://arxiv.org/pdf/2503.11638)) appeared in search but is a circuit-discovery paper (out of scope).

---

## Q2. qLDPC codes: has ML, LLM, evolutionary or RL search been used for BB, lifted-product, HGP or other qLDPC families, or for weight reduction and layout-aware codes? How were candidates evaluated?

### Takeaway
Yes. By November 2026 qLDPC families are the main target of automated code discovery. The methods used include:
- RL for weight reduction, for HGP optimization, and for BB code–decoder co-design;
- Bayesian optimization, GFlowNets and discrete diffusion;
- most visibly in 2026, LLM-guided program evolution and multi-agent "AI scientists": IBM's OpenEvolve campaigns on BB and PBB codes, Liu–Marquardt's structured concept evolution on lifted-product codes, Qian–Li's multi-agent search over balanced-product codes, and OmniQEC.

Candidates are evaluated in two ways:
- **Exact or certified distance:** exhaustive or rank-based methods and Lean certificates, with QDistRnd, MIP (HiGHS) and SAT serving as baselines.
- **Monte Carlo logical error rate:** under code-capacity noise with BP-OSD, or under circuit-level noise.

Since Sept/Oct 2026, a public, verifier-backed leaderboard (Unitary Foundation) collects mostly machine-found codes. Layout-aware discovery remains mostly non-AI.

### Cited Findings (newest first)

**LLM-guided and agentic search**
- **Cruz-Benito (Oct 2026), "Large Language Model-Guided Discovery of Weight-Five Bivariate Bicycle Codes", arXiv:2610.06623.** Any co-authors were not verified.
  - *Scope:* extends the IBM workflow to weight-5 BB and PBB codes, a more constrained family than weight ≥6.
  - *Catalogue:* 1,142 distinct proposals, of which 1,081 are non-baseline proposals generated by LLM programs. A basis-independent analysis flags 409 entries as disconnected.
  - *Distance:* of 936 distinct CSS proposals with a positive recorded distance, 816 (87.18%) are certified at d≥5.
  - *Comparison with random sampling:* the LLM campaign has a higher certified yield than the random-sampling controls, but the random controls cover more connected classes.
  - *Leaderboard:* entries from Table 3 ([[420,8,14]], [[432,8,15]], [[420,6,16]], [[420,8,13]]) were submitted to the Unitary Foundation leaderboard by GitHub user "cbjuan". The same handle hosts the cbjuan/qcode-discovery mirror, so these are author submissions that passed the leaderboard's verifier, not independent checks.
  - Sources: [arXiv](https://arxiv.org/abs/2610.06623); [PR 2879](https://github.com/unitaryfoundation/qldpc-challenge/pull/2879); [PR 2877](https://github.com/unitaryfoundation/qldpc-challenge/pull/2877); [PR 2875](https://github.com/unitaryfoundation/qldpc-challenge/pull/2875); [PR 2878](https://github.com/unitaryfoundation/qldpc-challenge/pull/2878); [repo PR 9](https://github.com/qiskit-community/qcode-discovery/pull/9)
- **Qian & Li (10 Aug 2026), "Multi-agent discovery of practical quantum LDPC codes", arXiv:2608.08996** (Dongheng Qian, Fudan; Tianyi Li, UW–Madison).
  - *Method:* a multi-agent framework whose programs instantiate coset-orbit balanced-product codes. This search space includes bicycle and lifted-product constructions as well as non-normal subgroup actions.
  - *Constraints:* binary CSS codes, n≤400, overall weight w≤10.
  - *Results:* rate–distance performance that is leading or competitive in each weight class, for example [[288,16,18]] at w=7, [[288,18,18]] at w=9 and [[234,28,18]] at w=10. Logical failure rates are low under code-capacity depolarizing noise with BP-OSD.
  - Whether these distances are exact or upper bounds was not retrieved.
  - Sources: [arXiv](https://arxiv.org/abs/2608.08996)
- **Ge Yan et al. (28 Jul 2026), "OmniQEC: discovering practical quantum error-correcting codes by an AI scientist", arXiv:2607.25865.** There are eight authors. A CatalyzeX listing names Min-Hsiu Hsieh and Yuxuan Du among them; this was not confirmed on the arXiv page.
  - *Method:* LLM reasoning plus agentic tool use. It optimizes the circuit-level logical error rate (LER) directly, under objectives and hardware budgets the user specifies.
  - *Scale:* 4 qLDPC construction families × 3 LLM backends (Claude, GPT, DeepSeek) × 14 total-physical-qubit budgets per backend.
  - *Result:* the discovered codes outperform the BB codes [[72,12,6]] and [[144,12,12]] under "complete-implementation" budgets of 98 and 240 physical qubits respectively. How those budgets are counted relative to n was not retrieved.
  - *Robustness:* 7 frontier codes per backend (21 in total) were re-run across full physical-error-rate sweeps, and the advantage persisted beyond the single error rate used during search.
  - *Metric:* −log10(p_LER), per round per logical qubit.
  - Sources: [arXiv](https://arxiv.org/abs/2607.25865v1); [HTML](https://arxiv.org/html/2607.25865v1); [CatalyzeX](https://www.catalyzex.com/author/Min-Hsiu%20Hsieh); [project GitHub](https://github.com/OmniQ-QID-Lab/Omniq-QiD-Lab.github.io)
- **Liu & Marquardt (Jun 2026), "Large-Language-Model Discovery of Quantum LDPC Codes through Structured Concept Evolution", arXiv:2606.24808** (Max Planck Institute for the Science of Light).
  - *Method:* structured concept evolution (SCE). The LLM edits construction programs within an algebraic mutation grammar, acting as an "intelligent mutation operator" in an evolutionary search over lifted-product CSS families.
  - *Constraints:* fewer than 600 data qubits and check weight ≤12.
  - *Score:* kd²/n, or circuit-level LER.
  - *Results:*
    - It beats published lifted-product benchmarks in some regimes, and the top codes have exactly certified distances.
    - At equal budget, the best LLM run outperformed random-mutation and hill-climbing baselines on both objectives.
    - Optimizing at circuit level often gives lower LER even when the distances are smaller.
  - *Inconsistency:* a separate search snippet mentions lifted-product codes "at n≈1500" trading rate against suppression, which conflicts with the 600-qubit limit. This is unresolved.
  - Sources: [arXiv](https://arxiv.org/abs/2606.24808); [alphaXiv](https://www.alphaxiv.org/abs/2606.24808)
- **Cruz-Benito, Cross, Kremer, Faro (IBM, 1 Jun 2026), "Evolutionary Discovery of Bivariate Bicycle Codes with LLM-Guided Search", arXiv:2606.02418.**
  - *Method:* OpenEvolve, an open-source MAP-Elites-style implementation of AlphaEvolve and FunSearch ideas. LLMs mutate Python programs that generate BB code ansätze and a new non-CSS "perturbed BB" (PBB) ansatz. PBB codes add two perturbation polynomials to the BB stabilizers, and they are valid codes only if AC^T + BD^T (mod 2) is symmetric.
  - *Verifier:* a staged pipeline that checks parameters, then distance, then equivalence class.
  - *Scale and cost:* 5 campaigns, about 1,650 iterations, about 2×10^5 candidates screened, about US$400 of LLM inference, and n≤360.
  - *Result:* 465 distinct codes (97 CSS BB and 368 non-CSS PBB). The search recovered known high-performing codes and found new ones, including an indecomposable [[288,16,12]]. Per a third-party summary, that code's kd²/n figure of merit (FOM) is 8.0; the campaign's best FOM, 12.0, belongs to a different code.
  - *Ablation:* eight arms on eight Stage-2 lattices, comprising random trinomials at two budgets, an exponent-tuple genetic algorithm, and a genetic algorithm over the LLM's ansatz representation.
    - Random search shows diminishing returns: a 10× larger budget gives only a 1.7× improvement.
    - On the lattices where the LLM fails, the exponent-tuple GA finds codes with d=4–6 and FOM ≤4.0, about 3× below the LLM's 12.0.
    - A generator-level GA produces more valid samples, but all have d≤2.
  - *Limitations:* IBM presents the outputs as candidates whose practical applicability needs further verification. Coverage also reports that the evolution first converged on codes with arbitrarily high k and no meaningful error correction, a reward-hacking failure. Which secondary outlet reported this was not pinned down.
  - Sources: [arXiv](https://arxiv.org/abs/2606.02418); [IBM Research blog](https://research.ibm.com/blog/ai-for-qec); [code](https://github.com/qiskit-community/qcode-discovery); [Quantum Computing Report](https://quantumcomputingreport.com/ibm-is-using-ai-to-help-identify-new-quantum-error-correction-codes/); [NAND Research](https://nand-research.com/ibms-llm-guided-evolutionary-framework-for-quantum-error-correction-code-discovery/); [OpenEvolve write-up "Does the LLM actually matter?"](https://algorithmicsuperintelligence.ai/blog/openevolve-qec/); [third-party summary](https://github.com/henryksobon300-code/qcode-discovery)

**Generative models, surrogates and Bayesian optimization**
- **Batista, Falcão, Machado (4 Nov 2026), "GFlowQEC: Generative Flow Networks that Learn to Sample Quantum LDPC Codes", arXiv:2611.06819** (16 pages).
  - *Method:* described as the first GFlowNet for quantum code design. It builds a code one group element at a time and proposes codes with probability tied to a power of their score. It is trained once on stored scores and then frozen.
  - *Bottleneck it addresses:* circuit-level scoring is so slow that "only a few hundred can be run".
  - *Scale:* lifted-product codes with 180 physical qubits.
  - *Result:* ranking ten of its draws per paid evaluation, it finds more than twice as many new top-scoring codes (inequivalent to any training code) as tuned parallel tempering. When both methods rank a hundred draws, a supervised screen beats GFlowQEC, and the screen does best when it ranks GFlowQEC's draws.
  - Sources: [arXiv](https://arxiv.org/abs/2611.06819)
- **LiftDiff (Oct 2026), "LiftDiff: Guided Discrete Diffusion for Quantum Code Construction", arXiv:2610.25213.** Authors were not retrieved. It claims to be the first discrete-diffusion model for generating QECCs, specifically lifted-product codes, and to yield codes with competitive rates, distances and low LER. Sources: [arXiv](https://arxiv.org/pdf/2610.25213)
- **Schneider & Nachmani (Oct 2026), "Distance Prediction for Quantum Codes", arXiv:2610.18130.**
  - *Method:* BB-Net, an attention-based network that predicts the distance of BB codes. It is trained on a dataset whose labels are distance upper bounds from high-budget heuristic searches.
  - *Claims:* 96.39% accuracy against 94.86% for state-of-the-art heuristics on the full benchmark, and 98.94% on the subset with certified exact distances.
  - *Independent reproduction:* 97.9% on the certified subset (the paper reports 98.9%), but 90.3% ±2 on the full benchmark (the paper reports 96.4%). QDistRnd reaches 92.3%, so the claim to be the most accurate method overall did not hold.
  - Sources: [arXiv](https://arxiv.org/abs/2610.18130); [official code](https://github.com/NadavSc/BBNet); [reproduction](https://github.com/jykim-UC/bbnet-reproduction)
- **Y. Chengyu, R. Meister, C. Carty, S.-K. Lin, R. Bondesan (Jan 2026), "Bayesian Optimization for Quantum Error-Correcting Code Discovery", arXiv:2601.18562.** Author names are as cited in later papers' reference lists.
  - *Method:* Bayesian optimization with a Gaussian-process surrogate. Each candidate's LER is estimated by Monte Carlo under depolarizing code-capacity noise.
  - *Claim:* "the first machine learning model that can accurately predict the logical error of quantum LDPC codes under code capacity noise", tested on HGP codes.
  - *Results (snippet-level and partly inconsistent):* a [[144,36]] code that prioritizes rate (three times the gross code's rate), and a [[144,26]] code with competitive per-qubit LER compared with the gross code. These may be two separate codes.
  - Sources: [arXiv HTML](https://arxiv.org/html/2601.18562); [PDF](https://arxiv.org/pdf/2601.18562)

**RL on qLDPC codes:** see Q1 for Sella et al. (BB codes, 2026), He & Liu (weight reduction to d≈35, 2025) and Freire et al. (HGP codes on the erasure channel, 2025).

**Non-AI structured searches (baselines and context)**
- **Lu, Yang, Guo (Sep 2026), "Quantum bivariate bicycle codes with weight-8 checks surpassing the BB benchmark", arXiv:2609.06572.** The search method was not determined from the snippets.
  - A census of 53 weight-8 codes at n=144. The best, [[144,6,d≥15]], beats the benchmark distance of 12, and [[144,10,12]] matches it.
  - At n=72, [[72,14,8]] has more than twice the kd²/n of the BB code of the same length.
  - The trade-off is a lower circuit-level pseudo-threshold, about 0.1% against about 0.4%, which the authors attribute to deeper circuits.
  - Sources: [arXiv](https://arxiv.org/abs/2609.06572)
- **"Sequences of Bivariate Bicycle Codes from Covering Graphs", arXiv:2511.13560.** New weight-8 BB codes, including [[64,14,8]] and [[144,14,14]]. Sources: [arXiv](https://arxiv.org/html/2511.13560)
- **"Self-dual bivariate bicycle codes with transversal Clifford gates", arXiv:2510.05211**, later published in npj QI (2026). Sources: [arXiv](https://arxiv.org/abs/2510.05211); [npj QI](https://www.nature.com/articles/s41534-026-01303-1)
- **Layout-aware but non-AI: "Coprime Bivariate Bicycle Codes and Their Layouts on Cold Atoms"**, Quantum (23 Feb 2026). Sources: [Quantum](https://quantum-journal.org/papers/q-2026-02-23-2009/)
- **The benchmark code: Bravyi et al., "High-threshold and low-overhead fault-tolerant quantum memory", arXiv:2308.07915.** This is the source of the [[144,12,12]] "gross" BB code that OmniQEC, the BO paper and the weight-8 paper all use as their reference. Sources: [arXiv](https://arxiv.org/pdf/2308.07915); [EC Zoo BB code](https://errorcorrectionzoo.org/c/qcga)

**Community verification and leaderboard**
- **Unitary Foundation "QEC Challenge"** (blog post by Farrokh Labib, 30 Sep 2026): "An open leaderboard for human- and AI-discovered quantum codes".
  - The post argues that QEC suits automated research because good codes "can be judged by numbers".
  - It states that the challenge does not yet measure fault-tolerant logical operations, only how well codes store qubits.
  - It invites new codes, better metrics, faster verifiers and refutations of existing entries.
  - The challenge is part of Unitary Foundation's AutoQEC project with NVIDIA, funded under Phase 1 of the US Department of Energy's Genesis Mission.
  - Sources: [blog](https://unitary.foundation/posts/2026_qec_challenge/); [blog PR](https://github.com/unitaryfoundation/unitary.foundation/pull/1092)
- **How the leaderboard works.**
  - Submissions must pass a verifier (`verify/qldpc_verify.py`). CSS and general stabilizer codes have separate boards.
  - Each cell (locality class × check weight) is ranked as a Pareto frontier, not by a single winner.
  - Most entries were found by machines, and their provenance names the model that produced them.
  - Snapshots returned by search show "1500+ verified codes from 40+ contributors" at one point and about 2,610 codes later. The top entries include a [[410,48,30]] generalized-bicycle code credited to a Gemini model. I did not check that entry's check weight or whether its d is exact.
  - Sources: [GitHub](https://github.com/unitaryfoundation/qldpc-challenge); [leaderboard](https://unitaryfoundation.github.io/qldpc-challenge/); [qecchallenge.com](http://qecchallenge.com/)
- **Example entry: [[360,12,20]], a non-CSS PBB code (kd²/n = 13.333).** None of its 24 checks was proven globally optimal within the time budget, so d=20 is an upper bound. Two independent methods reached d≤20, while a BP-OSD estimate gave 44. Sources: [entry page](https://qecchallenge.com/codes/360-12-20.html)
- **Another example:** a [[672,28,24]] weight-6 two-block code over C28×Dic3 with exact d (submitted by TimoKellerMath). Sources: [PR 5348](https://github.com/unitaryfoundation/qldpc-challenge/pull/5348)
- **Launch coverage** (14 Oct 2026) framed the challenge as an extension of the "ecdsa.fail" community-leaderboard model. A related Unitary Foundation post, "Agentic leaderboards for quantum information research", exists, but its content was not retrieved. Sources: [Quantum Zeitgeist](https://quantumzeitgeist.com/google-quantum-unitary-foundation-launched-qec/); [Unitary post](https://unitary.foundation/posts/2026_agentic_leaderboards/)

**How candidates were evaluated (cited)**
- Exact or certified distance was used by:
  - IBM's staged pipeline ([arXiv:2606.02418](https://arxiv.org/abs/2606.02418));
  - SCE's top codes ([arXiv:2606.24808](https://arxiv.org/abs/2606.24808));
  - the weight-5 campaign, which certifies d≥5 ([arXiv:2610.06623](https://arxiv.org/abs/2610.06623));
  - Lean certificates ([arXiv:2610.03214](https://arxiv.org/abs/2610.03214));
  - tropidist ([arXiv:2611.08247](https://arxiv.org/abs/2611.08247)).
- Baseline distance tools:
  - QDistRnd appears as a baseline in the BB-Net reproduction ([reproduction](https://github.com/jykim-UC/bbnet-reproduction)).
  - On a 74-instance benchmark, HiGHS (MIP) certified 13 instances and plain SAT 15, against 70 for tropidist ([Unitary issue #4716](https://github.com/unitaryfoundation/qldpc-challenge/issues/4716)).
- Monte Carlo LER was used at three levels:
  - code capacity: the BO paper ([arXiv](https://arxiv.org/html/2601.18562)), and Qian & Li with BP-OSD ([arXiv](https://arxiv.org/abs/2608.08996));
  - circuit level: OmniQEC ([arXiv](https://arxiv.org/abs/2607.25865v1)), SCE ([arXiv](https://arxiv.org/abs/2606.24808)) and GFlowQEC ([arXiv](https://arxiv.org/abs/2611.06819));
  - depolarizing noise with an area-under-curve reward: Sella et al. ([arXiv](https://arxiv.org/abs/2608.15754)).

### Inferences
- **The LLM edits constructions, not codes.** Every successful LLM-guided qLDPC search has the LLM edit *programs that instantiate algebraic families*: BB/PBB polynomials, lifts, balanced products, group elements. Validity (commutation) is then guaranteed or cheaply checked. The LLM's contribution is structured mutation, and the family choice bounds what can be found. In IBM's ablation the LLM's edge was finding algebraic families with d≥6 and FOM up to 12, not producing more valid samples.
- **kd²/n can be a misleading figure of merit.**
  - SCE and OmniQEC both report that optimizing at circuit level picks different winners.
  - The weight-8 BB codes have about 2× the kd²/n but about 4× lower circuit-level pseudo-threshold.
  - So kd²/n-only "beats the gross code" claims, and leaderboard FOMs, are weak evidence of practical advantage unless weight, connectivity and syndrome-circuit cost are matched.
- **Weight class decides dominance.** Qian & Li's [[288,16,18]] at w=7 dominates IBM's indecomposable [[288,16,12]] at w=6 in distance, at the cost of one more unit of check weight. This illustrates why the leaderboard ranks within weight and locality classes. Similarly, the Gemini-credited [[410,48,30]] has FOM = 48·900/410 ≈ 105, against 12 for the gross code. That comparison is meaningful only within its (unretrieved) weight class, and d may be only an upper bound.
- **Layout and connectivity.** No AI-driven search found in this session explicitly optimizes 2D or biplanar layouts or long-range-coupler constraints. OmniQEC budgets only the total number of physical qubits, and the cold-atom layout work is non-AI.

### Gaps
- Not confirmed from snippets: which simulators were used (stim/sinter) and which decoders, apart from BP-OSD in Qian & Li; Sella et al.'s results; the full IBM ablation table; OmniQEC's code parameters; LiftDiff's authors and numbers; whether Qian & Li's distances are exact.
- **[UNVERIFIED lead]** How Bravyi et al. (Nature, 2024) searched for BB codes and computed their distances.
- **[UNVERIFIED lead]** Lin & Pryadko, quantum two-block group-algebra codes, searched with QDistRnd distances (2023–24).
- **[UNVERIFIED lead]** Voss et al., multivariate bicycle codes (2024).
- **[UNVERIFIED lead]** Scruby, Hillmann, Roffe, trivariate/3D single-shot codes (2024).
- **[UNVERIFIED lead]** Sabo et al. and Hastings, weight-reduction baselines (the comparison points for He & Liu).
- **[UNVERIFIED lead]** Microsoft's 4D geometric codes (2025), reportedly found by lattice search.
- **[UNVERIFIED lead]** Pecorari et al., La-cross codes for long-range neutral atoms (2024).
- **[UNVERIFIED lead]** Berthusen et al., 2D-local qLDPC implementations (2024).

---

## Q3. Bosonic and autonomous QEC codes: automated or learned discovery and later developments

### Takeaway
This question could not be verified in this session, because the search budget ran out before targeted searches were run. The LLM- and agent-focused searches that were run (AlphaEvolve, OpenEvolve, FunSearch, "agentic" and "AI scientist" combined with QEC) returned no 2025–2026 work on discovering bosonic or autonomous-QEC codes. They were not aimed at bosonic codes, however, so this absence is weak evidence.

### Cited Findings
- The result sets of the agentic and LLM-guided searches contained only qubit qLDPC, EAQECC, fermionic-encoding and FT-circuit works, with no bosonic or autonomous-QEC entries. Examples: [Quantum Computing Report on IBM](https://quantumcomputingreport.com/ibm-is-using-ai-to-help-identify-new-quantum-error-correction-codes/); [qBraid AlphaEvolve blog](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction); [Lahiri & Marquardt FT circuits](https://arxiv.org/html/2611.18504); [OmniQEC](https://arxiv.org/abs/2607.25865v1)

### Inferences
- If the unverified leads below are accurate, bosonic and autonomous-QEC discovery up to 2023 relied on gradient-based and RL optimization over differentiable master-equation simulations. The 2026 wave of LLM-guided program evolution does not appear to have reached bosonic codes yet. This is a plausible gap, since fidelity under loss or dephasing can be computed exactly at small Fock cutoffs, which makes for a cheap, exact verifier.

### Gaps
- **[UNVERIFIED lead]** Wang, Rajabzadeh, Lee, Safavi-Naeini, "Automated discovery of autonomous quantum error correction schemes", PRX Quantum 3, 020302 (2022), believed to be arXiv:2108.02766. Recalled as gradient-based optimization of Hamiltonians and engineered dissipation (an "AutoQEC" framework) that rediscovers known schemes and proposes new small-Hilbert-space codes.
- **[UNVERIFIED lead]** Zeng, Zhou, Rinaldi, Gneiting, Nori, "Approximate autonomous quantum error correction with reinforcement learning", PRL 131, 050601 (2023).
- **[UNVERIFIED lead]** Sivak et al., "Model-free quantum control with reinforcement learning", PRX (2022), and "Real-time quantum error correction beyond break-even", Nature (2023). The latter is a GKP qubit whose QEC-protocol parameters were optimized by RL.
- **[UNVERIFIED lead]** Earlier numerically optimized or channel-adapted codes:
  - Reimpell & Werner, iterative code optimization (PRL 2005);
  - Fletcher, Shor, Win, SDP-optimal recovery (PRA 2007);
  - Albert et al., "Performance and structure of single-mode bosonic codes" (PRA 2018), with numerically optimized codes for amplitude damping;
  - Fösel, Tighineanu, Weiss, Marquardt, RL for quantum feedback (PRX 2018).
- Possible name collision: the Unitary Foundation's 2026 "AutoQEC project" (verified, Q2) is distinct from the "AutoQEC" framework recalled for Wang et al. 2022 (unverified). Do not conflate them.
- Any 2024–2026 bosonic work, such as learned GKP lattices, multimode codes or dissipative cat-code variants, remains entirely unchecked.

---

## Q4. Fermionic codes and fermion-to-qubit encodings: discovery by search, SAT, optimization, RL or LLM agents

### Takeaway
One agentic result is verified: qBraid's AlphaEvolve-driven search over error-correcting Generalized Superfast Encodings (GSE), June 2026. It reaches exact distance 5, and distance 6 on one 20-mode instance, for dense molecular Hamiltonians. The SAT, enumeration and optimization literature on fermion-to-qubit mappings (Chiew–Strelchuk, Fermihedral, ternary trees, Bonsai, treespilation) and on error-detecting and error-correcting encodings could not be re-verified in this session. No RL or LLM work on Majorana fermion codes surfaced.

### Cited Findings
- **Heitritter, Brown, Hardikar (qBraid, Jun 2026), "Evolving Quantum Error-Correcting Encodings for Molecular Simulation", arXiv:2606.25870.** qBraid's research page lists the paper as arXiv:2605.06792; this ID discrepancy is unresolved.
  - *Method:* LLM-guided evolutionary search using early access to Google's AlphaEvolve, plus an open-source replication.
  - *Target:* GSE fermion-to-qubit encodings.
  - *Verifier:* fixed while the search ran. AlphaEvolve returned readable Python that was scored by a verifier it was never allowed to modify.
  - *Results:*
    - The first GSE/superfast encodings beyond distance 3 for dense molecular Hamiltonians: exact d=5, and d=6 on one 20-mode instance.
    - A fixed-distance d=5 run evaluated 1,411 candidates and compressed the 272-qubit seed panel to 168 qubits, the "38% qubit cut" in the commentary's title.
    - qBraid reports logical error rates 3.4–7.9× lower under exact decoding.
  - A commentary contrasts the 168 qubits with about 800 for a conventional surface-code route; this secondary claim was not verified.
  - The work was scheduled as a talk at the 2nd "AI for Quantum" workshop, organized by NVIDIA and Quantinuum.
  - Sources: [arXiv](https://arxiv.org/abs/2606.25870); [HTML](https://arxiv.org/html/2606.25870v1); [qBraid blog](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction); [qBraid research](https://qbraid.com/research); [Medium commentary](https://medium.com/@adnanmasood/quantum-sundays-76-the-38-qubit-cut-what-alphaevolve-actually-found-inside-a-quantum-code-582b891cdb38); [Quantum Zeitgeist](https://quantumzeitgeist.com/brown-colleagues-designs-distance-5-gse-encodings-molecular-simulation/)
- **Context (titles only; not known to be AI-search papers):** "Efficient and Noise-Resilient Molecular Quantum Simulation with the Generalized Superfast Encoding" (arXiv:2511.09322) and "First fault-tolerant quantum memory demonstration for a generalized superfast encoding" (arXiv:2609.08957). Sources: [2511.09322](https://arxiv.org/pdf/2511.09322); [2609.08957](https://arxiv.org/pdf/2609.08957)
- **Access to the method:** AlphaEvolve became generally available through Google Cloud, which matters for reproducing AlphaEvolve-based searches. Sources: [Google Cloud blog](https://cloud.google.com/blog/products/ai-machine-learning/alphaevolve-is-available-for-everyone)

### Inferences
- Fermionic encodings suit agentic search well:
  - The verifier is exact and cheap at tens of modes: the distance of the encoding's stabilizer group and the Pauli weights of the mapped Hamiltonian terms.
  - The objective is naturally multi-objective (qubit count, distance, operator weight), which matches the quality-diversity (MAP-Elites) loop of AlphaEvolve and OpenEvolve.
- There are no public tables for fermionic encodings, so the qBraid result's novelty is judged against prior GSE constructions (d≤3), not against codetables.

### Gaps
- **[UNVERIFIED lead]** Enumeration, SAT and optimization of mappings:
  - Chiew & Strelchuk, "Discovering optimal fermion-qubit mappings through algorithmic enumeration", Quantum (2023);
  - Liu, Che, Zhou, Shi, Li, "Fermihedral: On the optimal compilation for fermion-to-qubit encoding", ASPLOS 2024 (SAT-based);
  - Jiang, Kalev, Mruczkiewicz, Neven, ternary-tree mappings, Quantum (2020; earlier);
  - Miller, Zimborás, Knecht, Maniscalco, García-Pérez, "Bonsai algorithm", PRX Quantum (2023);
  - Miller, Glos, Zimborás, "Treespilation" (2024);
  - "HATT", Hamiltonian-adaptive ternary trees (Liu et al., 2024);
  - Chien & Klassen, encodings optimized for both Hamiltonian and hardware (2022).
- **[UNVERIFIED lead]** Error-detecting and error-correcting encodings:
  - Bravyi, Terhal, Leemhuis, "Majorana fermion codes", NJP (2010; earlier);
  - Jiang, McClean, Babbush, Neven, Majorana loop stabilizer codes, Phys. Rev. Applied (2019);
  - Setia, Bravyi, Mezzacapo, Whitfield, superfast encodings, PRR (2019);
  - Derby, Klassen (with Bausch and Cubitt), compact encodings, PRB (2021), and low-weight high-distance fermionic encodings (2021);
  - Chen, Gorshkov, Xu, error-correcting codes for fermionic quantum simulation, SciPost Physics (2024).
- I could not check whether any RL or LLM-agent approach, other than qBraid's, has targeted Majorana codes or mapping enumeration.

---

## Q5. Subsystem, dynamical and Floquet codes: any automated or ML search for codes or measurement schedules?

### Takeaway
I found no automated, ML or LLM search for subsystem, Floquet or dynamical code definitions, or for the measurement schedules that define them. Targeted searches could not be run because of the budget. The nearest verified neighbors are out-of-scope circuit works: an LLM that evolves dynamic Clifford circuits under a fault-tolerance verifier, and MCTS syndrome scheduling. A verifier that certifies schedule-dependent (circuit) distances at scale also appeared in November 2026.

### Cited Findings
- **Lahiri & Marquardt (Nov 2026), "Agentic AI for Discovering Fault-Tolerant Quantum Circuits", arXiv:2611.18504.** An LLM writes and evolves Python programs that produce dynamic Clifford circuits for several quantum-computing tasks, and a verifier checks fault tolerance. It is out of scope (FT gadgets), but its program space, dynamic circuits with mid-circuit measurements, is the same object that defines Floquet and dynamical codes. The claimed depth reductions come from secondary coverage and are unverified. Sources: [arXiv HTML](https://arxiv.org/html/2611.18504); [CASRAI](https://casrai.org/news/agentic-ai-discovering-fault-tolerant-quantum-circuits-2026); a [Quantum Zeitgeist piece](https://quantumzeitgeist.com/fault-tolerant-circuits-llm-designs-error/) on LLM-designed fault-tolerant circuits, possibly about the same work (not confirmed)
- **AlphaSyndrome (arXiv:2601.12509; ACM DOI 10.1145/3779212.3790123).** MCTS scheduling of syndrome-measurement circuits for commuting-stabilizer codes, balancing error propagation against the decoder's correctable region. It is out of scope because it schedules existing codes. Sources: [arXiv](https://arxiv.org/pdf/2601.12509); [ACM](https://dl.acm.org/doi/10.1145/3779212.3790123)
- **tropidist.** It certifies circuit distances: it confirms a conjectured circuit distance of 11 for a three-color schedule of the gross code, and it handles memory circuits with up to 8×10^5 fault mechanisms. This is the kind of verifier a schedule search for dynamical codes would need. Sources: [arXiv](https://arxiv.org/abs/2611.08247); [GitHub](https://github.com/Fanerst/tropidist)
- **Unitary Foundation leaderboard.** It scores storage only, not fault-tolerant operations. Sources: [blog](https://unitary.foundation/posts/2026_qec_challenge/)

### Inferences
- For Floquet and dynamical codes, the quantity that must be verified is the circuit (spacetime) distance of the schedule. Certifying it at scale only became practical in late 2026 (tropidist), and the leaderboards have no category for schedule-defined codes, so search in this area appears to have been limited by the verifier more than by the searcher.

### Gaps
- **[UNVERIFIED lead]** Hastings & Haah's honeycomb/Floquet code (2021; earlier context).
- **[UNVERIFIED lead]** Delfosse & Paetznick, "Spacetime codes of Clifford circuits" (2023).
- **[UNVERIFIED lead]** Rodatz, Poór, Kissinger, "Floquetifying stabiliser codes with distance-preserving rewrites" (ZX-calculus, 2024). Its rewrite rules are a natural action space for search, but I found no ML search over them.
- **[UNVERIFIED lead]** Fu & Gottesman, "Error correction in dynamical codes" (2024).
- **[UNVERIFIED lead]** McEwen, Bacon, Gidney, "Relaxing hardware requirements for surface code circuits using time-dynamics" (2023; circuit level).
- No information on automated search for subsystem codes was gathered.

---

## Q6. LLM-guided discovery of codes (classical or quantum): FunSearch/AlphaEvolve-style results on codes, bounds or related combinatorial objects (2024–2026)

### Takeaway
Yes, and it accelerated sharply in 2026. The verified results are:
- **Classical deletion codes (2025):** a FunSearch adaptation rediscovered the VT code and found larger short codes for multiple deletions.
- **Classical binary linear codes (Sept 2026):** LinCodeEvolve found 7 new record codes, including [200,21,77], and raised 22 table entries. These were checked by exhaustive enumeration and reviewed by the table maintainer.
- **Entanglement-assisted quantum codes (Oct 2026):** AlphaEvolve plus Lean produced an infinite family that closes 54 table gaps.
- **qLDPC codes:** IBM's OpenEvolve campaigns, the weight-5 BB campaign, SCE, the multi-agent framework and OmniQEC (Q2).
- **Fermionic GSE encodings:** qBraid with AlphaEvolve (Q4).

The strongest "beyond the tables" claims come where public tables and an exact verifier exist, namely classical linear codes and EAQECCs.

### Cited Findings (newest first)
- **Baek, Won, Kim (30 Oct 2026), "Researcher-in-the-Loop Evolutionary Discovery of Entanglement-Assisted Quantum Code Families", arXiv:2610.41183** (cs.AI and cs.IT; 25 pages).
  - *Pipeline:* LLM-guided program evolution (AlphaEvolve is the evolution engine, per the repository), researcher direction, exact evaluation of (k,c,d), generalization, nonexistence certificates, and formal proof checking in Lean.
  - *Result:* early AlphaEvolve runs produced finite-qubit codes whose structure pointed to an infinite family. The authors generalized it as Theorem 1, including a prime-power construction, and checked the algebraic construction in Lean. Together with entanglement lifting, the family closes 54 listed gaps in archived EAQECC tables.
  - *Iterative versus independent program generation:* success rates on longer codes were 46.1% against 19.5% at n=13, and 31.6% against 3.1% at n=15.
  - Sources: [arXiv](https://arxiv.org/abs/2610.41183); [code, data and Lean proofs](https://github.com/sgmlc1234/llm_guided_EAQECC_search)
- **Seddas, Shashkov, Viazovska, Abbe (EPFL, 29 Sep 2026), "Evolving Towards Better Codes: LLM-Guided Search for High-Distance Binary Linear Codes", arXiv:2609.37056.**
  - *Method:* LinCodeEvolve, built on the EvoTune framework and the ShinkaEvolve codebase. It evolves code-construction programs against an exact minimum-distance evaluator. When progress stalls, a strategy loop that mixes diversity-driven search with expert input switches strategies.
  - *Results:* seven new record codes, including [200,21,77]. Six have short quasi-cyclic descriptions. Standard code-modification techniques then raise 22 entries in the code tables.
  - *Verification:* each code was checked by exhaustive enumeration, and the table maintainer has reviewed them for inclusion.
  - The authors conclude that LLM-guided search can supplement existing methods in coding theory.
  - Sources: [arXiv](https://arxiv.org/abs/2609.37056v1)
- **qLDPC campaigns (Jun–Nov 2026):** IBM's OpenEvolve on BB/PBB codes, the weight-5 BB extension, Liu–Marquardt SCE, the Qian–Li multi-agent framework and OmniQEC. Details are under Q2. Sources: [IBM blog](https://research.ibm.com/blog/ai-for-qec); [arXiv:2606.24808](https://arxiv.org/abs/2606.24808); [arXiv:2608.08996](https://arxiv.org/abs/2608.08996); [arXiv:2607.25865](https://arxiv.org/abs/2607.25865v1)
- **Fermionic GSE encodings with AlphaEvolve (qBraid, Jun 2026):** details under Q4. Sources: [arXiv](https://arxiv.org/abs/2606.25870)
- **LLM-Guided Search for Deletion-Correcting Codes (2025), arXiv:2504.00613** (authors not verified in this session).
  - *Method:* FunSearch adapted to find deletion-correcting codes at short lengths.
  - *Single deletion:* the search finds a function that the authors prove constructs the conjectured-optimal Varshamov–Tenengolts code, a rediscovery.
  - *Multiple binary deletions and quaternary edit correction:* larger short codes than prior search-based and explicit constructions, though these results are "largely empirical".
  - Sources: [arXiv](https://arxiv.org/abs/2504.00613)
- **Formalization with agents:** QECLean (Nov 2026) used AI agents extensively to formalize and verify recent codes such as BB codes in Lean. Sources: [arXiv](https://arxiv.org/abs/2611.21236)
- **Reproducibility of AI scientists:** Du, Zhao, Wang, Lu, Zhang, Tao (Nov 2026), "Reimagining computational quantum research with AI scientists", arXiv:2611.02833.
  - It introduces OmniQ, an expert-audited environment that turns 60 published computational-quantum studies into reproduction tasks.
  - Across five LLM-powered AI scientists and six third-party agent harnesses, AI scientists often failed to reproduce published results.
  - The authors propose "Residual Agents", which learn from verified workflow repairs.
  - Sources: [arXiv](https://arxiv.org/abs/2611.02833); [Quantum Zeitgeist](https://quantumzeitgeist.com/omniq-platform-verifies-quantum-papers/)
- **Adjacent general-agent work (not code-specific as far as retrieved):** Krenn et al., "Towards autonomous quantum physics research using LLM agents with access to intelligent tools" (arXiv:2511.11752); and OpenEvolve applied to error-mitigation software (arXiv:2609.38883; out of scope). Sources: [arXiv:2511.11752](https://arxiv.org/pdf/2511.11752); [arXiv:2609.38883](https://arxiv.org/pdf/2609.38883)

### Inferences
- The successful LLM-guided code searches share four features:
  1. They evolve *programs* (constructions), not code instances.
  2. They score them with exact or certified evaluators that the LLM cannot modify (stated explicitly by qBraid; implied by the exact evaluators in LinCodeEvolve and EAQECC).
  3. They restrict the search to algebraic families in which validity holds automatically.
  4. Increasingly, they add a human-in-the-loop or proof stage that turns finite finds into families: the EAQECC Lean pipeline and LinCodeEvolve's expert strategy loop.
- Iterative program refinement matters. The EAQECC comparison (46.1% vs 19.5%; 31.6% vs 3.1%) and IBM's ablation (where the LLM's edge was finding structured families) both suggest the LLM adds value through structured, cumulative edits rather than through raw sampling.

### Gaps
- **[UNVERIFIED lead]** FunSearch (Romera-Paredes et al., Nature 2024), including cap-set constructions.
- **[UNVERIFIED lead]** The AlphaEvolve white paper (Google DeepMind, 2025), including a recalled improvement of the kissing-number lower bound in 11 dimensions, which is a spherical-code problem.
- **[UNVERIFIED lead]** Georgiev, Gómez-Serrano, Tao, Wagner, "Mathematical exploration and discovery at scale" (2025).
- **[UNVERIFIED lead]** Nagda, Raghavan, Thakurta, AlphaEvolve-generated combinatorial structures such as MAX-k-CUT gadgets and Ramanujan graphs (2025).
- **[UNVERIFIED lead]** PatternBoost (Charton, Ellenberg, Wagner, Williamson, 2024) and Wagner's cross-entropy RL constructions (2021).
- I recall the authors of arXiv:2504.00613 as F. Weindel and R. Heckel (TUM), but this was not verified.
- I could not determine which code tables LinCodeEvolve's 22 entries belong to (presumably codetables.de) or how many record attempts failed.

---

## Q7. Across all of the above: rediscovery fraction, computational bottlenecks, and which verifiers made the search work

### Takeaway
Few works quantify rediscovery. Where numbers exist, a large share of the outputs are known, trivial, disconnected or only bounded from above. The bottlenecks are distance evaluation and circuit-level LER simulation. For distance, the hard part is certifying lower bounds, which is NP-hard in general; circuit-level simulation is costly enough that only a few hundred evaluations are affordable per search. The searches became credible through exact enumeration, rank-based (tropical) dynamic programming, Lean certificates, deduplication by equivalence class, and verifiers fixed outside the agent's reach. BP-OSD-based distance estimates and learned surrogates have proven unreliable as the final arbiter.

### Cited Findings
- **Rediscovery and triviality**
  - Su et al. saturate the CSS LP bound at n=13, which means reaching optimal known parameters ([arXiv](https://arxiv.org/abs/2305.06378)).
  - Mauron et al. find "the best possible codes" within a restricted action set ([arXiv](https://arxiv.org/abs/2305.11470)).
  - Tandeitnik & Guerreiro's genetic algorithm rediscovered the Shor and 5-qubit codes ([arXiv](https://arxiv.org/pdf/2210.05058)).
  - The FunSearch adaptation for deletion codes rediscovered the VT code ([arXiv](https://arxiv.org/abs/2504.00613)).
  - IBM's search "recovered known high-performing codes and found new ones". Its evolution first converged on high-k codes with no protection, and its generator-level GA produced only d≤2 codes ([IBM blog](https://research.ibm.com/blog/ai-for-qec); [arXiv](https://arxiv.org/abs/2606.02418); [summary](https://github.com/henryksobon300-code/qcode-discovery)).
  - In the weight-5 campaign, 409 of 1,142 entries are flagged as disconnected, and 816 of 936 CSS proposals are certified at d≥5 ([arXiv](https://arxiv.org/abs/2610.06623)).
  - GFlowQEC counts only codes inequivalent to its training set as new ([arXiv](https://arxiv.org/abs/2611.06819)).
- **Genuinely beyond the tables**
  - LinCodeEvolve: 7 records plus 22 raised entries, maintainer-reviewed ([arXiv](https://arxiv.org/abs/2609.37056v1)).
  - EAQECC: 54 table gaps closed ([arXiv](https://arxiv.org/abs/2610.41183)).
  - GSE encodings with d=5 and d=6, beyond the prior d≤3 ([arXiv](https://arxiv.org/abs/2606.25870)).
- **Bottlenecks**
  - Certifying a distance lower bound is NP-hard in general (Lean-QEC) ([arXiv](https://arxiv.org/html/2605.16523)).
  - Exact distance computation and heuristic upper-bound searches "suffer from severe computational bottlenecks" ([BB-Net](https://arxiv.org/abs/2610.18130)).
  - Circuit-level scoring allows only "a few hundred" evaluations ([GFlowQEC](https://arxiv.org/abs/2611.06819)).
  - Monte Carlo LER evaluation is costly enough to motivate a Gaussian-process surrogate ([BO](https://arxiv.org/html/2601.18562)).
  - Random search gives diminishing returns (10× budget for 1.7× gain), and the LLM cost was about $400 for about 2×10^5 candidates ([arXiv:2606.02418](https://arxiv.org/abs/2606.02418)).
- **Reliability of verifiers and claims**
  - **tropidist audit of 26 papers.** Code distances reported as exact held up, but in 10 papers at least one circuit-distance or code-distance estimate was too high; the repository reports a worst-case overstatement of at least 8. The tool certifies a code distance of 24 for a modular hyperbolic [[31200,1954]] code, and a coverage article reports a certified distance of 16 for IBM's two-gross circuit. Sources: [arXiv](https://arxiv.org/abs/2611.08247); [GitHub](https://github.com/Fanerst/tropidist); [Quantum Zeitgeist](https://quantumzeitgeist.com/ibm-quantum-code-distances-rank-based/)
  - **Benchmark comparison:** on 74 instances, tropidist certified 70, against 13 for HiGHS (MIP) and 15 for plain SAT ([Unitary issue #4716](https://github.com/unitaryfoundation/qldpc-challenge/issues/4716)).
  - **Leaderboard [[360,12,20]]:** d is only an upper bound, while a BP-OSD estimate gave 44 ([entry](https://qecchallenge.com/codes/360-12-20.html)).
  - **BB-Net reproduction:** 90.3% overall against the claimed 96.4%, with QDistRnd at 92.3% ([reproduction](https://github.com/jykim-UC/bbnet-reproduction)).
  - **Lean certificates (An & Yang, arXiv:2610.03214):**
    - The kernel checks certificates by computation. In enumeration form, certificate size is provably polynomial in n at fixed distance, and smaller for product families.
    - Coverage: eleven families and 39 parameter sets, the widest at 1,872 qubits.
    - On an 18-qubit toric code, the whole-file check drops from 42 s to 9 s.
    - The [[144,12,12]] distance is certified.
    - The arXiv listing's title omits "error-correcting".
    - Sources: [arXiv](https://arxiv.org/abs/2610.03214); [QECCertificates](https://github.com/QCL-SUAT/QECCertificates)
  - **Further verification tools:** QECLean ([arXiv](https://arxiv.org/abs/2611.21236)) and PC-QEC, which produces symbolic certificates independent of distance for syndrome-extraction circuits ([arXiv](https://arxiv.org/html/2611.12985)).
  - **AI scientists' reproducibility:** they often fail to reproduce published computational-quantum results ([OmniQ](https://arxiv.org/abs/2611.02833)).
- **Verifiers used inside the search loops**
  - IBM: parameters, then distance, then equivalence classes ([arXiv](https://arxiv.org/abs/2606.02418)).
  - qBraid: a fixed verifier the LLM could not modify, plus exact decoding ([blog](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)).
  - LinCodeEvolve: exact minimum distance plus exhaustive enumeration ([arXiv](https://arxiv.org/abs/2609.37056v1)).
  - EAQECC: exact (k,c,d), nonexistence certificates and Lean ([arXiv](https://arxiv.org/abs/2610.41183)).
  - Olle et al.: a Knill–Laflamme reward computed with a vectorized Clifford simulator ([npj QI](https://www.nature.com/articles/s41534-024-00920-y)).
  - Freire et al.: erasure-channel performance ([arXiv](https://arxiv.org/abs/2501.09622)).

### Inferences
- **Computed from the weight-5 numbers:**
  - 409/1,142 ≈ 35.8% of catalogue entries are disconnected, which presumably means decomposable into smaller codes and so not new as codes.
  - 816/936 ≈ 87.2% of CSS proposals with a recorded distance are certified at d≥5, so about 12.8% had d<5 or were uncertified.
- Novelty is easiest to establish where public tables and exact verifiers exist (classical linear codes, EAQECCs). For qLDPC codes at fixed weight and locality there was no codetables-like reference until the 2026 leaderboard, so "new" often meant new within a family or new at a given n.
- The reward-hacking and false-novelty modes observed so far are:
  - high k with no distance (IBM);
  - codes that are valid but have d≤2 (GA baselines);
  - disconnected or decomposable codes (weight-5);
  - distances that are only upper bounds (the leaderboard);
  - inflated circuit-distance claims (the tropidist audit);
  - winners that look good at code capacity but not at circuit level (SCE; the weight-8 BB codes).

### Gaps
- No paper reports a systematic rediscovery rate against codetables.de, the Error Correction Zoo or the leaderboard. IBM's split between known and new codes was not retrieved.
- Which works used stim/sinter, Magma or QDistRnd inside the loop could not be confirmed from snippets, beyond the cases cited above.

---

## Q8. Gaps and opportunities

### Takeaway
The 2025–2026 boom is concentrated on qubit CSS qLDPC memories within algebraic families (BB, lifted-product, balanced-product, generalized-bicycle) and is scored by kd²/n or by memory LER. Several areas remain essentially untouched by agentic search, as far as found:
- Floquet, dynamical and subsystem codes;
- bosonic and autonomous QEC;
- most fermionic and Majorana codes (beyond the qBraid GSE result);
- codes optimized for logical gates or magic;
- independent re-verification and proof of agent-found results.

### Cited Findings (sources that state or demonstrate the gaps)
- The leaderboard does not yet measure fault-tolerant logical operations, only storage ([Unitary blog](https://unitary.foundation/posts/2026_qec_challenge/)).
- IBM presents its codes as candidates needing further verification ([IBM blog](https://research.ibm.com/blog/ai-for-qec)).
- Circuit- and code-distance overstatements appear in 10 of 26 audited papers ([tropidist](https://arxiv.org/abs/2611.08247)).
- Increasing distance at fixed weight is the harder direction for learning methods ([He & Liu](https://arxiv.org/abs/2502.14372v1)).
- Optimizing at circuit level selects different codes from optimizing kd²/n ([SCE](https://arxiv.org/abs/2606.24808)), and weight-8 codes with higher kd²/n have a 4× lower pseudo-threshold ([weight-8 BB](https://arxiv.org/abs/2609.06572)).
- Only "a few hundred" circuit-level evaluations are affordable, and supervised screening beats the generative sampler at large rank budgets ([GFlowQEC](https://arxiv.org/abs/2611.06819)).
- The learned distance predictor's headline claim was not reproduced ([reproduction](https://github.com/jykim-UC/bbnet-reproduction)).
- AI scientists often fail to reproduce published results ([OmniQ](https://arxiv.org/abs/2611.02833)).
- From-scratch RL is still at 25 qubits and d=5, with 100 qubits and d=10 only a stated "roadmap" ([Olle et al.](https://www.nature.com/articles/s41534-024-00920-y)).
- Proof-backed generalization from finite finds to infinite families has been demonstrated for EAQECCs ([Baek et al.](https://arxiv.org/abs/2610.41183)), but no comparable proof stage was reported for any of the qLDPC searches in Q2.

### Inferences: concrete unexplored problems (general, not tailored to any group)
1. **Agentic search for Floquet, dynamical and subsystem codes.** The action space could be measurement schedules (sequences of low-weight checks) or ZX rewrite rules. The verifier would be the circuit (spacetime) distance from detector-error models, which tropidist-class tools can now certify, plus circuit-level LER. No such work was found.
2. **LLM program evolution for bosonic and autonomous QEC.** The search could run over Hamiltonian and dissipator constructions, codeword families, or GKP and multimode lattices, with exact master-equation fidelities at small cutoffs as a cheap, exact verifier. No 2025–26 agentic work was found.
3. **Fermionic codes beyond GSE.** Three directions are open:
   - extending AlphaEvolve-style search to Majorana fermion codes and error-correcting compact or lattice encodings, with explicit operator-weight bounds;
   - hardware-native (connectivity-constrained) encodings that also carry a distance;
   - proofs that finite finds generalize to families.
4. **Raising distance at fixed weight or locality**, the harder direction according to He & Liu. Any claim here should be benchmarked within leaderboard weight classes.
5. **Certified distances inside the loop.** Exact lower-bound certification (tropidist, Lean certificates) could be used as part of the reward, for example with budgeted certification as shaping, rather than applied after the fact. Leaderboards should flag or penalize entries whose distance is only an upper bound.
6. **Calibrated surrogates.** Distance and LER predictors need uncertainty estimates and multi-fidelity schemes (code capacity, then phenomenological, then circuit-level noise). The BB-Net reproduction failure and the GP-surrogate claims show that this is unsettled.
7. **Canonical forms and equivalence testing** for non-CSS and group-algebra codes, so that rediscovery rates can be measured against codetables, the Error Correction Zoo and leaderboards. IBM's equivalence classes and the weight-5 connectivity analysis are only first steps.
8. **Objectives beyond memory.** Targets could include transversal, fold-transversal or automorphism gates, codes relevant to distillation (for example triorthogonal-like codes), single-shot properties, and the cost of syndrome circuits. None of these are scored by the current leaderboard.
9. **Beyond binary CSS qubit codes.** Qudit qLDPC codes, non-CSS families at large n (PBB is a start), and non-stabilizer codes have seen little or no agentic work in what was found.
10. **From finite discoveries to theorems for qLDPC codes.** The researcher-in-the-loop Lean pipeline used for EAQECCs could be applied to LLM-found qLDPC families, to get asymptotic distance and rate guarantees, and QECLean could be used to formally verify agent-found codes.
11. **Scaling hardware-in-the-loop code tailoring** (Yanay) beyond 7-qubit stabilizers, including re-optimization as device noise drifts.
12. **Compute-matched, open baselines.** The field needs a standard protocol that reports the number of expensive evaluations, LLM cost and wall time, and includes random, GA, parallel-tempering and hill-climbing baselines. IBM, SCE and GFlowQEC each ran ad hoc versions of such comparisons.
13. **Independent re-verification as standard practice.** Given the tropidist audit and OmniQ's reproduction failures, independent re-checks of agent-found claims could become a routine role for community leaderboards.
14. **Layout- and connectivity-aware AI search** for qLDPC codes, with explicit 2D, biplanar or long-range-coupler embeddings as objectives. Current layout work is non-AI, and OmniQEC constrains only the total qubit budget.

### Gaps
- Because the search budget ran out, I could not check whether some of these problems, especially items 1–3, are addressed in works this session did not reach.
- The Unitary Foundation post "Agentic leaderboards for quantum information research" may describe further plans for leaderboards, but its content was not retrieved ([post](https://unitary.foundation/posts/2026_agentic_leaderboards/)).
