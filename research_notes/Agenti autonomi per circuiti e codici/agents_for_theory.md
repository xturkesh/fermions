# LLMs and LLM agents for theory in quantum information, QEC and quantum many-body physics (focus: 2024 – Nov 2026)

Scope and method. This is a literature map, newest work first within each theme. Every arXiv ID was confirmed through the arXiv API (export.arxiv.org) on 2026-11-22. AI-contribution statements are quoted or paraphrased from each paper's abstract or its "AI use / disclosure" section, read from the arXiv HTML or PDF. Model names (e.g. "GPT-5.6 Sol", "GPT-6 Astra", "Claude Fable 5", "Gemini 4") appear exactly as the authors wrote them. Many of them postdate my training data, and I could not check them independently. Pages outside arXiv (openai.com, anthropic.com, scottaaronson.blog, henryyuen.net, github.com) could not be fetched from this environment. Facts from those pages come only from web-search summaries and are marked "(search summary, not fetched)". Out of scope: code and circuit search, except where it produced theorems or conjectures.

---

## 1. AI-assisted proofs and results in quantum information, quantum complexity, QEC bounds and entanglement theory (plus physics exemplars)

### Takeaway
Between September 2025 and November 2026 the AI's role changed. In Aaronson–Witteveen (GPT-5-Thinking, Sept 2025) it suggested one technical step. By mid-2026 frontier models were producing the core argument for open problems, which humans then verified, rewrote and, more and more often, formalized in Lean. Examples are I3322 non-attainment, the optimal constant for quantum ordered search, the FGG QAOA conjecture, superactivation of private capacity, an explicit Tsirelson counterexample, a quantum-Hedetniemi counterexample, and OpenAI's computer-assisted Haldane-gap proof. OpenAI also released whole batches of AI-generated results, including quantum parallel repetition (Aug 2026), and a Haldane gap, a 2D area law and N(6)=3 MUBs (Oct 2026). Those batches are only partly formalized and are still being audited.

### Cited Findings

**2025 precursors**
- **Aaronson & Witteveen (later with Phillip Harris), "Limits to black-box amplification in QMA", arXiv:2509.21131.** v1 25 Sep 2025; v4 1 Oct 2026; journal ref. Quantum 10, 2227 (2026).
  - *Result:* an oracle relative to which no polynomial-resource QMA verification gets completeness closer to 1 than doubly exponential, or soundness that is super-exponentially small. This matches Jeffery–Witteveen amplification and makes the Aaronson (2008) QMA vs QMA₁ separation quantitative. — [arXiv abs](https://arxiv.org/abs/2509.21131)
  - *AI vs human:* the v1 footnote says the idea of studying the rational function r(z)=Tr[(I−P(z))⁻¹] and bounding the growth of rational functions "was suggested to us by GPT-5-Thinking". — [v1 HTML](https://arxiv.org/html/2509.21131v1)
  - *Later simplification:* v4's "Author contributions and use of AI" says the GPT-5 function required the no-instances to be an interval. After posting, PH (added as coauthor) "suggested to replace Tr[P_rej(θ)⁻¹] by the simpler det[P_rej(θ)]", which gave the "simplified argument". — [v4 HTML](https://arxiv.org/html/2509.21131v4)
  - *Aaronson's blog account (search summary, not fetched):* about 30 minutes of back-and-forth; the first answer was confident but wrong. He added that he and Witteveen could have solved it themselves "given a week or two". — [Aaronson, "The QMA Singularity"](https://scottaaronson.blog/?p=9183)
  - *Verification:* human proof plus peer review. *Limitation:* the AI supplied one key idea, and a simpler human idea later replaced it.
- **Bubeck et al. (14 authors), "Early science acceleration experiments with GPT-5", arXiv:2511.16072 (20 Nov 2025, 89 pp.).**
  - Case studies across mathematics, physics, astronomy, CS, biology and materials, including "four new results in mathematics (carefully verified by the human authors)". — [arXiv abs](https://arxiv.org/abs/2511.16072)
  - *Physics chapters:*
    - Lupsasca: GPT-5 Pro "(re)derived nontrivial Lie point symmetries—including an SL(2,R) algebra" of the stationary, axisymmetric wave equation (black holes). The chapter is filed under "independent rediscovery of known results".
    - Scherrer: power spectra of gravitational radiation from cosmic strings.
    - Spears: reduced-physics modelling of thermonuclear burn. — [HTML](https://arxiv.org/html/2511.16072)
  - I found no quantum-information case study in the parts of the table of contents I inspected.

**Physics exemplar: the GPT-5.2 gluon-amplitude conjecture**
- **Guevara, Lupsasca, Skinner, Strominger, Weil ("on behalf of OpenAI"), "Single-minus gluon tree amplitudes are nonzero", arXiv:2602.12176** (12 Feb 2026; v2 9 Mar 2026).
  - *Result:* single-minus tree n-gluon amplitudes are nonzero on "half-collinear" configurations in Klein space or for complexified momenta. The paper gives a piecewise-constant closed form for one minus-helicity gluon decaying into n−1 plus-helicity gluons. — [arXiv abs](https://arxiv.org/abs/2602.12176)
  - *AI vs human:* "The key formula (39) … was first conjectured by GPT-5.2 Pro and then proved by a new internal OpenAI model."
  - *Verification:* "checked by hand using the Berends–Giele recursion", and shown to "nontrivially obey the soft theorem, cyclicity, Kleiss–Kuijf, and decoupling identities". — [HTML](https://arxiv.org/html/2602.12176v1)
  - *Limitations:* a special kinematic region; no formal proof; the proving model is internal and not public.

**Lab-scale releases (2026)**
- **OpenAI, "Ten Advances in Mathematics and Theoretical Computer Science" (1 Aug 2026).**
  - *Process:* per an independent audit, "an internal model produced the arguments, humans prepared the manuscripts with model assistance, and the model subsequently formalized the results in Lean".
  - *Lean status:* the repo snapshot (commit 94bc0fe, 2 Aug 2026) specifies Lean 4.32.0 and lists sorry_count = 0 for the principal results.
  - *Topics* include quantum information. Chapter 6 is quantum parallel repetition. — [Sienicki & Sienicki, arXiv:2608.14673](https://arxiv.org/abs/2608.14673)
  - *Claimed result (search summary, not fetched):* "Exponential parallel repetition for arbitrary finite, two-player quantum games" (Lean file QuantumParallelRepetition.lean), from an internal version of OpenAI's "Astra" model. — [ten-proofs README](https://raw.githubusercontent.com/openai/ten-proofs/main/README.md); [Forklog](https://forklog.com/en/news/openai-says-internal-astra-model-produced-10-results-in-mathematics-and-theoretical-computer-science)
  - *Exponent (search summary, unverified):* the gap exponent in OpenAI's bound is reported as 13. A follow-up preprint by Zhao Song claims exponent 3. — [IACR ePrint 2026/1702](https://eprint.iacr.org/2026/1702)
  - *Expert reaction (search summary, not fetched):* Henry Yuen presumes the proof correct given the Lean check, but writes that it "doesn't give me any understanding". — [Yuen blog](https://www.henryyuen.net/posts/on-openai-and-quantum-parallel-repetition/)
  - *Audit outcome:* an early claim of a "polarity error" in Chapter 6's greedy-conditioning lemma was withdrawn. PDF extraction had dropped an overbar from the formula. No confirmed substantive error remains in the assessments examined. — [arXiv:2608.14673 HTML](https://arxiv.org/html/2608.14673)
- **OpenAI math release, 6 Oct 2026 (github.com/openai/math):** 722 AI-written manuscripts in 372 problem families. A survey of 38 families states: "No Lean development was compiled or kernel-checked for this survey." — [Xie, arXiv:2610.17090](https://arxiv.org/abs/2610.17090). Entries relevant here, all from [arXiv:2610.17090 HTML](https://arxiv.org/html/2610.17090):
  - *Family 268, "The spin-one Haldane gap"* — marked "not formalized; computer-assisted".
    - Claims a unique ground state for even L ≥ 60 and a uniform gap γ_L > (4/105)·log(80/79) > 4.79×10⁻⁴ (L ≥ 60), and > 3.8×10⁻³ (L ≥ 2304). The numerical gap is 0.4105.
    - Odd rings are not covered. A companion paper discharges Tasaki's Assumption 2 for boundary-field chains.
    - The new tool is a "coupled purity bootstrap" relating Gibbs purity and spatial purity. It does not need frustration-freeness.
  - *Family 265, "Two-dimensional area law from a global gap"* — "not formalized". Claims S(A) ≤ C(q,R,J,Δ)|∂A| for square-lattice Hamiltonians with a unique gapped ground state; frustration is allowed.
  - *Family 266, "Three mutually unbiased bases in dimension six"* — claims N(6)=3 (Zauner's conjecture) as a computer-assisted theorem, "conditional on arithmetic and compiler conditions".
    - The headline computation is C++ in binary64 with hand-proved rounding margins (4,470 s) and uses "neither interval arithmetic nor exact algebra".
    - A companion proves N(6) ≤ 5 exactly with a rational verifier. The Lean formalization covers only the companion.
  - *Family 221* — a Mézard–Parisi formula for diluted spin glasses (Lean, 405 files).
  - OpenAI's preprint "The periodic spin-one Haldane gap" is cited as dated 24 Sept 2026. — [Grover, arXiv:2611.07400 HTML](https://arxiv.org/html/2611.07400)
- **Verification and extensions of the Haldane-gap proof**
  - **Sirui Lu, "Bridging Gaps with Tensor Networks: From the Haldane Gap and Beyond", arXiv:2611.11049 (8 Nov 2026).**
    - *Checks:* independently reproduced "all of OpenAI's finite trace and matrix-product trial contractions in exact arithmetic". The finite seeds use certified eigenvalue enclosures (rational/interval arithmetic, Ostrowski–Weyl bounds) on symmetry-resolved blocks. — [PDF](https://arxiv.org/pdf/2611.11049)
    - *Restatement:* recast the proof as tensor-network diagrams and identified reflection positivity as the general condition behind it.
    - *Stronger result:* sharpened OpenAI's bound (Δ_L > log20/784 ≈ 0.003821 for even L ≥ 2304) by a factor of 19, to Δ_L > 0.0729 for all even L ≥ 36. Also certified gaps for seven other chains and ladders.
    - *Formalization:* the general criterion is formalized in Lean 4, with kernel-checked seed arithmetic (github.com/LionSR/CertifiedGaps). — [arXiv abs](https://arxiv.org/abs/2611.11049)
  - **Tarun Grover, "Extensions of the Haldane gap proof: higher dimensions, near degeneracy, fermions and tensor-network states", arXiv:2611.07400 (5 Nov 2026).** The proof idea is an inequality between the purity of the Gibbs state and of a "dual Gibbs state" obtained by space-time rotation.
    - *AI vs human:* ChatGPT Astra and Codex were used "extensively". The author developed the main ideas for higher dimensions, twists and fermions; AI refined them and worked through the details.
    - For PEPS, "AI developed much of the bootstrap argument". — [HTML](https://arxiv.org/html/2611.07400)
  - **Hal Tasaki, arXiv:2610.34330 (25 Oct 2026).** Uniqueness of the translation-invariant ground state of the infinite S=1 Heisenberg chain and of the J–D model at small |D|, built on OpenAI's machinery and numerical inputs.
    - "The author formulated the problems and proposed the core strategies. ChatGPT Astra made substantial contributions to developing the proofs."
    - Uniqueness among non-translation-invariant ground states is not established. — [arXiv abs](https://arxiv.org/abs/2610.34330)
  - **Hal Tasaki, arXiv:2609.34119 (v2 Oct 2026).** Proves that Ogata's and Tasaki's SPT indices are equivalent for U(1)⋊Z₂ chains.
    - Combined with OpenAI's boundary-field gap theorem, this gives a locally unique gapped ground state of the infinite S=1 chain with both indices equal to −1.
    - ChatGPT Astra made "substantial contributions to the development of the proof". — [arXiv abs](https://arxiv.org/abs/2609.34119)
- **Google DeepMind, AlphaProof Nexus → quantum optics theorem**
  - **Tsoukalas et al., "Advancing Mathematics Research with AI-Driven Formal Proof Search", arXiv:2605.22763 (21 May 2026).** A Lean proof-search agent that "autonomously resolved 9 of 353 open Erdős problems" and "proved 44/492 OEIS conjectures". It is "being deployed in … quantum optics research". — [arXiv abs](https://arxiv.org/abs/2605.22763)
  - **Krenn, Firsching, Tsoukalas, Gajjala, Gu, Chaudhuri, arXiv:2611.26300 (19 Nov 2026).**
    - *Result:* no complex-weighted monochromatic quantum graph exists for even N ≥ 4 when D ≥ N. So, in this regime, N-photon D-dimensional GHZ states cannot be generated by that probabilistic pair-source model.
    - *AI vs human:* the proof was "discovered through AI-assisted formal proof search in Lean using AlphaProof Nexus"; the authors supply the physical interpretation.
    - *Limitation:* only the D ≥ N regime. — [arXiv abs](https://arxiv.org/abs/2611.26300)
  - *Caveat on Lean-certified results (search summary, not fetched):* Lean verifies a proof only against the statement as written in Lean. — [Science Media Centre](https://sciencemediacentre.es/en/alphaproof-nexus-new-ai-tool-designed-tackle-mathematical-proofs)

**Researcher-led 2026 results with explicit AI disclosures**

*Bell nonlocality: three near-simultaneous resolutions of the Pál–Vértesi I3322 conjecture*
- **Jef Pauwels, arXiv:2608.29734 (30 Aug 2026).**
  - *Result:* the I3322 quantum supremum is not attained in finite dimension, so C_q(3,3,2,2) is not closed and approaching the supremum needs unbounded dimension. — [arXiv abs](https://arxiv.org/abs/2608.29734)
  - *How the proof arose:* on 25 July 2026 the conjecture was given to "OpenAI Codex (GPT-5.6-sol)". After a nudge ("Make a very strong push"), "roughly three hours later, the first end-to-end proof appeared". It was written in non-standard terms ("classical shadow", "positive cocycle").
  - *Verification:* Claude served "as translator and hostile referee"; the core argument was formalized in Lean 4.
  - The author could not verify Seth Douglas's concurrent AI-assisted repository. — [HTML](https://arxiv.org/html/2608.29734)
- **Seth Douglas, arXiv:2609.05555 (3 Sep 2026; v2 11 Sep).**
  - *Result:* no finite-dimensional strategy attains S, while a "spatial" strategy on ℓ²(Z)⊗ℓ²(Z) does; the dimension needed is D(ε)=Θ(log 1/ε).
  - *Correction:* v2 "replaces an unsupported step in the earlier lower argument".
  - *Verification scope:* exact arithmetic and Lean check specified finite facts only; "the complete analytic proof is not formalized". — [arXiv abs](https://arxiv.org/abs/2609.05555)
  - *AI use:* LLMs were used for "proof discovery, implementation, verification engineering, adversarial auditing", across multiple model families. Auditing sessions were instructed to refute. — [HTML](https://arxiv.org/html/2609.05555)
- **Andrea Coladangelo, "The I3322 Bell inequality requires infinite dimensions", arXiv:2609.06038 (5 Sep 2026, 79 pp.).**
  - *AI vs human:* "An approximate proof, containing all the main ideas, was produced by GPT 5.5 Pro after various rounds of interaction." That first proof was "much shorter, and arguably quite difficult to parse by a human". The author revised it for correctness and clarity, with help from GPT 5.6 Sol.
  - *Technique:* a Bellman-function certificate. — [HTML](https://arxiv.org/html/2609.06038)

*Quantum query complexity*
- **Carolan & Childs, "Optimal Quantum Algorithms for Ordered Search", arXiv:2609.34293 (28 Sep 2026).**
  - *Result:* closes the longstanding constant-factor question with (1/π)·ln n + o(log n) queries.
  - *AI vs human:* the zero-error algorithm was "discovered by Claude Fable 5" from a continuum relaxation. The exact algorithm was "discovered by GPT-5.6-Sol (informed by Claude's zero-error algorithm)" and solves the Farhi–Goldstone–Gutmann–Sipser polynomial program. — [arXiv abs](https://arxiv.org/abs/2609.34293)
  - *Verification:* "The core technical results in this paper were discovered autonomously by AI systems. All proofs have been verified and refined by the authors." — [HTML](https://arxiv.org/html/2609.34293)

*QAOA*
- **Kol, Ben-Shahar, Sulimany, Englund, "A Machine-Verified Proof of a Quantum-Optimization Conjecture", arXiv:2606.29687 (29 Jun 2026).**
  - *Result:* proves the FGG conjecture that depth-p QAOA on the ring of disagrees attains ratio (2p+1)/(2p+2) exactly.
  - *Human part:* formalized QAOA and the known parts of the problem, and reduced the conjecture to one open Lean statement.
  - *AI part:* Claude Fable 5 constructed the proof and "uncovered a hidden dynamical symmetry".
  - *Verification:* Lean 4 accepts the proof with Mathlib's standard classical axioms. Humans check only that the formal statement is faithful. — [arXiv abs](https://arxiv.org/abs/2606.29687); [HTML](https://arxiv.org/html/2606.29687)

*Nonlocal games, operator algebras, quantum graph theory*
- **Gao, Yu, Zhi, "An Explicit Counterexample to Tsirelson's Problem via a Linear System Game", arXiv:2610.10248 (7 Oct 2026).**
  - *Result:* a binary linear-system game that separates C_qa from C_qc. It has 1,417,152 equations in 1,889,684 variables, and its classical value is computed exactly.
  - *Verification:* the full system is specified in Lean 4, and the separation is formalized with Mathlib. — [arXiv abs](https://arxiv.org/abs/2610.10248)
  - *AI use:* "GPT 6 Astra and Fable 5.1" helped with exploration, the construction, the Lean code and drafting. — [HTML](https://arxiv.org/html/2610.10248)
- **Julius A. Zeiss, "A counterexample to the quantum Hedetniemi conjecture", arXiv:2609.20690 (17 Sep 2026).**
  - *Result:* χ(G×H) ≤ 1538 < 1539 = min(χ_q(G), χ_q(H)). The counterexample also holds for the C*-algebraic variants.
  - *Verification:* Lean 4 formalization plus exact integer certificates. — [arXiv abs](https://arxiv.org/abs/2609.20690)
  - *AI vs human:* "the proof idea underlying the main theorems and the construction of the counterexample is the result of prompting a generative AI system"; the author verified and simplified it. — [HTML](https://arxiv.org/html/2609.20690)
- **Zeiss, arXiv:2609.15757 (14 Sep 2026).** Determines the optimal entanglement-assisted zero-error source-coding cost under a balanced-difference promise, and resolves the spectral part of Cao et al.'s Conjecture 6.3.
  - *Verification:* every lemma is checked in Lean 4.19.0 with Mathlib.
  - *AI use:* OpenAI Codex helped with exploration, proofs and Lean. — [arXiv abs](https://arxiv.org/abs/2609.15757); [HTML](https://arxiv.org/html/2609.15757)
- **Zhi Li & Xiaofei Shi, arXiv:2608.06235 (6 Aug 2026).**
  - *Result:* the sharp bound ‖z‖₁ ≤ √2·min{n,m}·‖z‖_ε, via an L₁ noncommutative Khintchine inequality with Haar-unitary coefficients. Gives an improved data-hiding bound. The upper bound is Lean-checked. — [arXiv abs](https://arxiv.org/abs/2608.06235)
  - *AI use:* ChatGPT and Claude drafted the manuscript and developed the Lean code. The authors "manually verified that the Lean statements match". — [HTML](https://arxiv.org/html/2608.06235)

*Channel capacities and QKD*
- **Chengkai Zhu & Xin Wang, arXiv:2609.10520 (9 Sep 2026).**
  - *Result:* resolves whether two zero-private-capacity channels can jointly communicate privately. A four-level channel plus a qubit erasure channel with erasure probability 1/2 give more than 0.0001903 private bits per product use.
  - *AI vs human:* "The initial activation example was identified through interactions with large language models" (QudeLeap's AI Quantum Scientist).
  - *Verification:* formalized in Lean 4. — [arXiv abs](https://arxiv.org/abs/2609.10520); [HTML](https://arxiv.org/html/2609.10520)
- **Ben Shahar, Morag, Englund, Sulimany, "Quantum Cryptography and Autoformalization in Lean", arXiv:2610.17286 (12 Oct 2026).**
  - *Result:* an end-to-end Lean theorem for composable security of entanglement-based BB84 against coherent attacks (postselection technique, finite block size), written in an LOCC-enforcing language. The same pipeline was then reused "as a discovery machine, reducing the finite-size penalty on the key rate". — [arXiv abs](https://arxiv.org/abs/2610.17286)
  - *AI use:* Claude and GPT agents did the Lean translation. A local GLM-5.3-Flash model, "under the direction of Claude Opus 5", wrote several formal results. — [HTML](https://arxiv.org/html/2610.17286)

*Entanglement measures*
- **Holevo & Shirokov, arXiv:2608.27363 (27 Aug 2026).** Builds a state at which the entanglement of formation has no global supporting affine functional. This refutes an assumption used in several articles. Found "using Wootters' formula and the help of Claude Fable 5". — [arXiv abs](https://arxiv.org/abs/2608.27363)
- **Ma, Shirokov, Zhang, arXiv:2610.18816 (13 Oct 2026).** The discrete and continuous definitions of entanglement of formation coincide in infinite dimensions; "the main idea of the proof (realized with the AI-assistance)". — [arXiv abs](https://arxiv.org/abs/2610.18816)
- **Holevo & Shirokov, arXiv:2609.12667 (11 Sep 2026).** Log-Sobolev-based continuity bounds; the authors thank "GPT 5.6 for valuable discussion and technical help". — [arXiv abs](https://arxiv.org/abs/2609.12667)

*Quantum cryptography and complexity*
- **Seyoon Ragavan, arXiv:2607.21811 (23 Jul 2026).** The first plain-model, one-time, information-theoretically secure and efficient unclonable encryption of one bit. "GPT-5.6 Sol Ultra found this proof in an extended conversation with the author and drafted a preliminary version." — [arXiv abs](https://arxiv.org/abs/2607.21811)
- **Coladangelo, Qipeng Liu, Ziyi Xie, arXiv:2608.17629 (18 Aug 2026).** A simultaneous Goldreich–Levin reduction for two entangled parties, giving unclonable encryption from BB84 states. "Discovered by GPT-5.6 Ultra after a few interactions." — [arXiv abs](https://arxiv.org/abs/2608.17629)
- **Longcheng Li, Qian Li, Xingjian Li, Qipeng Liu, arXiv:2608.03824 (4 Aug 2026).** Perfectly complete quantum key agreement cannot be built from quantumly secure one-way functions in a black-box way (QROM). "GPT-5.6 Sol Ultra found this proof in a one-shot conversation and drafted a preliminary version." — [arXiv abs](https://arxiv.org/abs/2608.03824)
- **Aparna Gupte, arXiv:2611.18857 (13 Nov 2026).** A quasipolynomial quantum algorithm for the vector dihedral coset problem with smooth modulus. It gives no new attack on LWE. "The main algorithm was obtained in a GPT-5.6 Sol conversation." — [arXiv abs](https://arxiv.org/abs/2611.18857)
- **Yanqiao Wang, Jin-Peng Liu, Peng Li, Yang Liu, arXiv:2606.24899 (12 Jun 2026).** Case study of human–AI co-discovery, using the AIM system, of sign-embedding quantum algorithms. "The decisive scientific judgments remained human", for example rejecting a Cayley-trapezoidal approximation whose validity "required a hidden condition". — [arXiv abs](https://arxiv.org/abs/2606.24899)

*Mathematical many-body physics and open systems*
- **Ruizhi Liu & Jinmin Yi, arXiv:2610.24154 (16 Oct 2026).** Translation-invariant pure quasifree fermionic states, and through Jordan–Wigner–Araki also spin states, whose half-chain factors are of type II_∞ or III_λ. "The counterexample was discovered with AI assistance" (GPT-5.6 Sol, Claude Opus 5.5). — [arXiv abs](https://arxiv.org/abs/2610.24154); [HTML](https://arxiv.org/html/2610.24154)
- **Kamijima, Nakajima, Noda, Sasaya, arXiv:2610.18440 (13 Oct 2026).** A central limit theorem for log Z_N(β) of the Gaussian q-body SYK model at small β. A case study of human–AI proof refinement followed by Lean 4 formalization of the main arguments (github.com/TokyoProbabilityGroup/SYK_CLT). — [arXiv abs](https://arxiv.org/abs/2610.18440)
- **Becker & Zworski, arXiv:2609.13121 (11 Sep 2026).** Optimal trace-norm relaxation γ_h = λ₁(h)/(2h) for Witten Lindbladians. "An appendix by Chat GPT 6 presents a stronger result specialised to positive initial data." — [arXiv abs](https://arxiv.org/abs/2609.13121)
- **Projansky, Wang, Reiss, Whitfield, arXiv:2611.26315 (19 Nov 2026).** All non-contextual Pauli Hamiltonians are "free-fermion in disguise". Claude Opus 5.5 "generate[d] the Lean code that formalizes Theorems 2 and 4". — [HTML](https://arxiv.org/html/2611.26315)
- **Hongru Zhao, arXiv:2610.00210 (21 Sep 2026).** Minimax trace-norm tomography rates with periodic Clifford measurements. There are "Lean formalizations of all labeled theoretical results", and GPT-5.6 Sol / GPT-6 Astra helped develop the proofs. — [arXiv abs](https://arxiv.org/abs/2610.00210); [HTML](https://arxiv.org/html/2610.00210)
- **Valentin Link, arXiv:2610.37012 (27 Oct 2026).** The TEMPO tensor network can be made independent of system size. "This structure was discovered by Claude Opus 5.5 during guided algorithm development." — [arXiv abs](https://arxiv.org/abs/2610.37012)

### Inferences
- **AI share of the work:** the typical 2026 pattern in quantum-information theory is that the AI proposes the core construction or proof, and humans verify, simplify, rewrite and often formalize it. Human judgment still decides problem selection, statement fidelity and exposition (Carolan–Childs, Coladangelo, Pauwels, Wang et al.).
- **Lean as the default check:** single-author or small-team papers now routinely ship Lean formalizations (Zeiss, Li–Shi, Zhao, Gao–Yu–Zhi, Pauwels, Zhu–Wang). AI writes most of that Lean, so formalization cost is no longer the bottleneck.
- **Credit and priority:** the I3322 problem was solved three times within one week (30 Aug – 5 Sep 2026), each time with heavy AI assistance. Problems that frontier models can now reach may get solved in parallel by many groups.
- **What breaks most often:**
  - unformalized analytic or numerical steps (Douglas v2's unsupported step; OpenAI's binary64 MUB computation; the Haldane and 2D area-law results are not formalized);
  - AI arguments that need "translation" before humans can read them (Pauwels, Coladangelo).

### Gaps
- I could not read the OpenAI, Anthropic, Aaronson or Yuen pages directly. The parallel-repetition exponent (13), Song's exponent-3 claim, and the model attributions ("Astra") rest on search summaries.
- I found no independent peer review yet of OpenAI's quantum parallel-repetition, Haldane-gap, 2D area-law or N(6)=3 claims. The Haldane gap is independently re-verified only for its finite inputs plus the general criterion (Lu); the 2D area law has no independent check that I could find.
- AI involvement is not stated in the abstracts of several Lean-verified QEC papers I did not examine in full, e.g. Bharti, Haug, Tao, Ye, arXiv:2610.02621.

---

## 2. Formalization infrastructure for quantum theory (Lean and other proof assistants) and autoformalization agents

### Takeaway
By late 2026 there is a large set of Lean 4 libraries for quantum theory:
- quantum-information theory: Lean-QuantumInfo (now merged into Physlib), Lean-Quantum, Lean-QIT, the MIP*=RE components, TNLean;
- QEC: QECLean, Lean-QEC, QECert, kernel-size distance certificates;
- quantum programs: the QHL-based libraries QEDIA and DeepMind's, plus Lean-QuantumAlg;
- physics: Physlib/PhysLean, and the Osterwalder–Schrader / Glimm–Jaffe free-field QFT construction.

Most of these libraries were built largely by LLM agents. The repeated lesson is that the kernel check is cheap and trustworthy. The hard part is making sure the formal statement means the intended physics.

### Cited Findings

**Quantum-information-theory libraries**
- **Meiburg, Lessa, Soldati, "A Formalization of the Generalized Quantum Stein's Lemma in Lean", arXiv:2510.08672 (9 Oct 2025).** Formalizes the Hayashi–Yamasaki proof (the original proof had a gap).
  - The authors call it "the most technically demanding theorem in physics with a computer-verified proof to date".
  - Formalization "rectified minor imprecisions" in that proof. Part of the Lean-QuantumInfo library. — [arXiv abs](https://arxiv.org/abs/2510.08672)
  - *Merger:* Physlib is described as "formerly PhysLean and Lean-QuantumInfo". — [Tooby-Smith, arXiv:2603.08139](https://arxiv.org/abs/2603.08139)
- **Kasaura … Yamasaki, "Lean-Quantum: Toward AI-Assisted Formalization of Quantum Information", arXiv:2607.05492 (6 Jul 2026).**
  - *Content:* the data-processing inequality for sandwiched Rényi relative entropy, with strong subadditivity as a corollary.
  - *Infrastructure:* operator monotonicity and convexity, Lieb–Ando trace inequalities, Choi, Kraus and Stinespring representations, Haar measure.
  - It provides "the last missing component needed to complete the Lean formalization of the generalized quantum Stein's lemma". — [arXiv abs](https://arxiv.org/abs/2607.05492)
- **Zhu … Wang, "Lean-QIT", arXiv:2607.09632 (10 Jul 2026).** Formalizes Schumacher compression, the HSW classical capacity, and entanglement-assisted capacity with its strong converse. Operational definitions are kept separate from entropic characterizations. — [arXiv abs](https://arxiv.org/abs/2607.09632)
- **Lu, Tjoa, Cirac, "Multi-agent Autoformalization of Tensor Network Theory" (TNLean), arXiv:2607.07857 (8 Jul 2026).**
  - *Content:* the fundamental theorem of matrix-product states, extended toward 1D SPT phases. "For some statements, the agents were able to explore new proof routes".
  - *Main bottleneck:* "enforcing mathematical intent". — [arXiv abs](https://arxiv.org/abs/2607.07857)
  - *AI share:* "The Lean 4 code and the blueprint have been fully generated by the multi-agent AI system … under our supervision." — [HTML](https://arxiv.org/html/2607.07857)
- **Lu, Deng, Zhu, Ji, "Long-horizon autoformalization of a core theorem underlying MIP* = RE" (FormalFlow), arXiv:2609.19814 (17 Sep 2026).**
  - *Result:* quantum soundness of the classical low individual-degree test, done in 63 days; 126,367 lines of Lean, "all generated by agents".
  - *Errors found:* the work "corrects side conditions and intermediate error bounds in the published proof and improves the error exponent". — [arXiv abs](https://arxiv.org/abs/2609.19814)
  - *Tools:* TeXRA, Claude Code, OpenCode and Codex. A follow-up formalized the Pauli basis test using "orchestrating GPT-6 Astra sessions". — [HTML](https://arxiv.org/html/2609.19814)
- **Tianrun Zhao & Nengkun Yu, "Formalizing CHSH Rigidity in Lean 4", arXiv:2604.03884 (4 Apr 2026).** Found a gap in McKague–Yang–Scarani (arXiv:1203.2976). — [arXiv abs](https://arxiv.org/abs/2604.03884)
  - *The gap:* the definition X′_B = (B₀+B₁)/|B₀+B₁| fails when B₀±B₁ is not invertible; there is a concrete counterexample in dimension 2.
  - *LLM failure mode:* coding assistants "would quietly add hypotheses to a lemma, and the proof became trivial". — [HTML](https://arxiv.org/html/2604.03884)
- **Palsberg, Huang, Almanei, arXiv:2610.03364 (2 Oct 2026).** Clifford+T universality for any distance measure satisfying four axioms, formalized in Lean. — [arXiv abs](https://arxiv.org/abs/2610.03364)
- **Yin, Taylor, Englund, Koppens, "AxQM", arXiv:2609.05157 (4 Sep 2026).**
  - 1,019 kernel-checkable proof tasks over 479 items from Nielsen & Chuang, in a custom finite-dimensional QM library; "the largest proof-synthesis benchmark in physics by a factor of four".
  - Grading by the Lean kernel: no `sorry` and no new axioms. — [arXiv abs](https://arxiv.org/abs/2609.05157)
  - The paper documents a semantic defect in which "the circuit's action on an arbitrary eigenvector is taken as a hypothesis", so the task "is solved trivially". It notes that "this semantic deficiency again evades the Lean kernel". — [HTML](https://arxiv.org/html/2609.05157)

**Physics libraries**
- **Joseph Tooby-Smith's Physlib (formerly HepLean/PhysLean):**
  - HepLean — [arXiv:2405.08863](https://arxiv.org/abs/2405.08863) (May 2024)
  - index notation — [arXiv:2411.07667](https://arxiv.org/abs/2411.07667) (Nov 2024)
  - Wick's theorem (static and normal-ordered versions) — [arXiv:2505.07939](https://arxiv.org/abs/2505.07939) (May 2025)
  - SU(5) model building, "from scans to theorems" — [arXiv:2603.28406](https://arxiv.org/abs/2603.28406) (Mar 2026)
- **Douglas, Hoback, Mei, Nissim, "Formalization of QFT", arXiv:2603.15770 (16 Mar 2026).** The free bosonic Euclidean QFT in 4D satisfies the Glimm–Jaffe (OS-type) axioms.
  - Intended as proof that "existing AI tools" can carry extended mathematical-physics arguments into machine-checked proofs.
  - The initially assumed results (Minlos, nuclearity of Schwartz space, Goursat) have since been proven or avoided. — [arXiv abs](https://arxiv.org/abs/2603.15770)

**QEC formalization and distance certificates**
- **Jain, Raveendran, Albarghouthi, "Specifying & Verifying Quantum Error-Correcting Codes" (QECLean), arXiv:2611.21236 (16 Nov 2026).** "Code combinators" prove distances analytically.
  - It shows that IBM's [[360,12,≤24]] code has distance exactly 24, previously known only as an upper bound, and formalizes the surface-code family.
  - No AI involvement is stated in the abstract. — [arXiv abs](https://arxiv.org/abs/2611.21236)
- **Ehatamm, Lee, Wu, Tao, "End-to-End Formalization of Quantum Error Correction" (Lean-QEC), arXiv:2605.16523 (15 May 2026).**
  - Uses a verified reduction of distance to SAT, giving Lean-checked distances for [[90,8,10]] and [[70,6,9]] bivariate-bicycle codes.
  - The approach scales to 144 qubits only when run outside the kernel. — [arXiv abs](https://arxiv.org/abs/2605.16523)
- **Shuoming An & Fusheng Yang, arXiv:2610.03214 (2 Oct 2026).** Short certificates whose correctness the Lean kernel decides by computation, with a size theorem.
  - Certifies [[144,12,12]]. It settles the affine-permutation codes printed as [[1152,580,≤12]] and [[2304,1156,≤14]] at distances 12 and 14.
  - Covers 11 families and 39 parameter sets, trusting only the three standard axioms. — [arXiv abs](https://arxiv.org/abs/2610.03214)
  - "No solver in the trusted base". The largest instance is n=1872, d=7. — [HTML](https://arxiv.org/html/2610.03214)
- **He, Chen, Xiong, Ni, Zhang, Wang, "QECert", arXiv:2611.05047 (3 Nov 2026).** Formalizes hypergraph, balanced and lifted products, quantum Tanner codes (asymptotically good qLDPC codes), Bravyi–Poulin–Terhal bounds, constant-overhead magic-state distillation, and energy barriers. — [arXiv abs](https://arxiv.org/abs/2611.05047)
- **Bharti, Haug, Tao, Ye, arXiv:2610.02621 (2 Oct 2026).** Subsystem codes with dressed distance proven without numerical search; Lean snippets throughout the paper. — [arXiv abs](https://arxiv.org/abs/2610.02621)

**Quantum programs and algorithms (quantum Hoare logic in Lean)**
- **Shirobokov, Jaderberg, … Albarghouthi, Jain, … Babbush, Chaudhuri, Kohli, "Discovering formally verified quantum programs with AI", arXiv:2611.09822 (6 Nov 2026).**
  - *Setup:* a Lean quantum-program library in which "Gemini 4 agents synthesize novel programs against formal problem specifications and verify their correctness and cost".
  - *Reported gains:* parallel search up to 2.57×, quantum arithmetic 290×, fusion 12×, photochemistry 8×, drug metabolism 2.5×. — [arXiv abs](https://arxiv.org/abs/2611.09822)
  - *Workflow:* researchers write a Lean quantum-Hoare-logic (QHL) specification of correctness and cost (T gates). Agent teams iterate on Lean-kernel feedback, share results through a wiki and "conferences", and "an independent verifier audits every proof against a frozen specification".
  - *Problems:* arcsin, parallel Grover search, cytochrome P450, fusion α-particle stopping power, non-Born–Oppenheimer dynamics. The point is to certify programs too large to simulate. — [HTML](https://arxiv.org/html/2611.09822)
- **Ren, Zhang, Liu, Li, Qi, Ying, "QEDIA", arXiv:2610.17708 (13 Oct 2026).** Formalizes QHL+ (with classical variables) in Lean and proves it sound. The library and the program proofs (QFT, phase estimation, repeat-until-success, Hamiltonian simulation) were built "under human supervision, using an LLM-based prover agent". — [arXiv abs](https://arxiv.org/abs/2610.17708)
- **Lei Zhang, Jing, Zhao, Liu, Wang, "Lean-QuantumAlg", arXiv:2610.30216 (21 Oct 2026).**
  - An "Algorithm-Solves-Problem contract" ties correctness, success probability and resource bounds to a single execution.
  - Nine AI-formalized algorithms, from Grover up to attacks on RSA-2048 and P-256. Fewer than 1% of imported source code needs statement-level human review once the shared framework has been audited. — [arXiv abs](https://arxiv.org/abs/2610.30216)
  - Related work: Shor in Lean [arXiv:2607.14082](https://arxiv.org/abs/2607.14082); QNN certification [arXiv:2607.12981](https://arxiv.org/abs/2607.12981).
- **Liyi Li, Tagba, Potanin, arXiv:2610.34723 (26 Oct 2026).** Translates the Qafny quantum proof system into classical Hoare logic with array theory; the completeness direction is mechanized in Rocq. — [arXiv abs](https://arxiv.org/abs/2610.34723)
- **Older context (Coq):**
  - CoqQ — [arXiv:2207.11350](https://arxiv.org/abs/2207.11350)
  - the VOQC verified circuit optimizer — [arXiv:1912.02250](https://arxiv.org/abs/1912.02250)
  - survey "Formal Verification of Quantum Programs: Theory, Tools and Challenges" — [arXiv:2110.01320](https://arxiv.org/abs/2110.01320)
  - Mingsheng Ying on symbolic specification of quantum data and operations — [arXiv:2512.22383](https://arxiv.org/abs/2512.22383)

**Autoformalization agents applied to quantum theorems**
- **Ren, Li, Qi, "MerLean", arXiv:2602.16554 (18 Feb 2026).** Fully automated LaTeX → Lean 4 → LaTeX round trip, so the result can be reviewed in readable form. Run on three theoretical quantum-computing papers: 2,050 Lean declarations from 114 statements, all three formalized end-to-end. — [arXiv abs](https://arxiv.org/abs/2602.16554). The agent is based on Claude Code. — [HTML](https://arxiv.org/html/2602.16554)
- **Xi He, Sirui Lu, Bei Zeng, arXiv:2510.20728 (23 Oct 2025).** TeXRA extended with an independent Lean 4 verification layer, used to certify code constructions and no-go results (details in §4). — [arXiv abs](https://arxiv.org/abs/2510.20728)
- **Ax-Prover (Breen et al.), arXiv:2510.12787.** A multi-agent Lean prover with a quantum-physics benchmark. The arXiv admin withdrew it because it "does not meet arXiv's research content quality standards". — [arXiv abs](https://arxiv.org/abs/2510.12787)
- **Benchmarks for these agents:** Lean-QuantumAlg-Bench (36 tasks) and Lean-QIT-Bench (40 tasks), arXiv:2607.21533 (23 Jul 2026).
  - Best difficulty-weighted scores are 60.4 and 59.6 out of 100 (GPT-5.5, Kimi K3, DeepSeek V4-Pro, MiniMax M3).
  - Giving agents access to a verified library ("library-augmented deduction") improved all eight comparisons, by up to 15.9 points.
  - Weak areas: "quantum simulation, quantum learning, quantum information measures, and entanglement theory". — [arXiv abs](https://arxiv.org/abs/2607.21533)
- **Position paper: Krippendorf, Sarra, Tooby-Smith, arXiv:2610.14059 (6 Oct 2026).** Argues that AI "remains weaker at choosing definitions and designing libraries than at writing proofs", and that "a checked proof does not establish that its definitions faithfully capture the intended physics". — [arXiv abs](https://arxiv.org/abs/2610.14059)

### Inferences
- The infrastructure is now dense enough to retrieve against. Library access improves agent proving scores, and DeepMind's quantum-program agents use a shared QHL library as their action space.
- The main residual risk is semantic: added hypotheses, statement drift, and definitions that misencode the physics. Practical defences have converged on four things:
  - frozen specifications checked by an independent verifier (DeepMind);
  - explicit human-reviewed contracts (Lean-QuantumAlg);
  - LaTeX back-translation for review (MerLean);
  - blueprints with periodic human review (TNLean, FormalFlow).
- For QEC, distance certification now has three distinct trust models: SAT replay inside the kernel, kernel-sized certificates with a size theorem, and analytic code combinators.

### Gaps
- I found no published head-to-head comparison of these Lean QI libraries, and no statement of which are being upstreamed into Mathlib or Physlib beyond Lean-QuantumInfo's merger.
- AxQM baseline scores were not found in the parts I read.
- I did not find Isabelle-based quantum work on arXiv in this session; Isabelle Marries Dirac returned no arXiv hit.

---

## 3. LLM agents doing analytic many-body and field-theory calculations (Hartree–Fock, phase diagrams, stat-mech mappings, Bethe ansatz, perturbation theory)

### Takeaway
Agents now do real analytic many-body work, but with tight external checks. Examples:
- Hartree–Fock derivations from papers (Kim group);
- coordinate Bethe-ansatz solutions of new integrable chains, checked by independent exact diagonalization;
- extensions of OpenAI's Haldane-gap proof (Grover, Tasaki);
- autonomous DMRG phase diagrams cross-checked against ED.

Where checks are only numerical, agents can be "right for the wrong reason", for example passing partition-function checks while misidentifying the tractable structure.

### Cited Findings
- **Pozsgay & Vona, "Bethe Ansatz with a Large Language Model", arXiv:2603.29932 (31 Mar 2026).**
  - *Result:* the LLM "semi-autonomously solved the task in all cases": coordinate Bethe-ansatz solutions for three integrable Hamiltonians, two of them new. One is PT-symmetric but breaks left-right invariance. Another has a nested Bethe ansatz whose nesting level is free-fermionic without U(1) invariance — "found by the LLM".
  - *Errors:* "a few mistakes along the way", corrected after the humans spotted them.
  - *Verification:* results "checked against exact diagonalization (performed by separate programs)", and derivations checked by the authors.
  - *Models:* ChatGPT 5.2 Pro and 5.4 Pro. — [arXiv abs](https://arxiv.org/abs/2603.29932)
- **Pan, Mudur, Taranto, Tikhanovskaya, Venugopalan, Bahri, Brenner, Eun-Ah Kim, "Quantum Many-Body Physics Calculations with Large Language Models", arXiv:2403.03154 (Mar 2024; Commun. Phys. 8, 49 (2025)).**
  - Multi-step prompt templates for Hartree–Fock derivations on 15 papers. "With correction of intermediate steps", GPT-4 derived the correct final HF Hamiltonian in 13 cases and made minor errors in 2. Average step score 87.5/100.
  - LLMs were also used to extract template information from papers and to score the steps. — [arXiv abs](https://arxiv.org/abs/2403.03154)
- **Pan, Roggeveen, … Brenner, Eun-Ah Kim, "CMT-Benchmark", arXiv:2510.05228 (Oct 2025; ICLR 2026).** 50 expert-written condensed-matter-theory problems: HF, ED, QMC/VMC, DMRG, statistical mechanics, model building.
  - *Grading:* programmatic, including "symbolic handling of non-commuting operators via normal ordering".
  - *Results:* the best model (GPT-5) solves 30%. The 17-model average is 11.4 ± 2.1%. 18 problems are solved by no model. — [arXiv abs](https://arxiv.org/abs/2510.05228)
  - The arXiv entry notes that CMT-Benchmark was referenced in the Gemini 3 Deep Think release (Feb 2026). — [arXiv abs](https://arxiv.org/abs/2510.05228)
- **Nys, Jreissaty, Carrasquilla, "QuantumMatterExplorer", arXiv:2610.32001 (22 Oct 2026).**
  - An LLM agent that takes a Hamiltonian and parameter ranges and returns a phase diagram, with DMRG "cross-checked against independent exact diagonalization".
  - On a frustrated bilinear-biquadratic spin-1 chain with single-ion anisotropy, whose phase diagram was published only recently, it found Néel, Haldane and dimerized phases plus incommensurate correlations, within an 8-hour compute budget. — [arXiv abs](https://arxiv.org/abs/2610.32001)
- **Wanyu Zhao & Wanbing Zhao, "Can AI Agents Discover Statistical Mechanical Mappings?", arXiv:2607.26367 (29 Jul 2026).**
  - StatMechBench-v0: six Ising-type problems covering transfer matrices, gauge-removable disorder, and planar/Pfaffian structure.
  - Numerical feedback helps agents repair code, but agents "can also pass the numerical checks while misidentifying the underlying tractable class". The authors call for symbolic checks and structural invariants. — [arXiv abs](https://arxiv.org/abs/2607.26367)
- **Deng, … Di Luo, "Towards Verifiable and Self-Correcting AI Physicists for Quantum Many-Body Simulations", arXiv:2604.00149 (31 Mar 2026).** QMP-Bench: 100 tasks drawn from 21 journals. PhysVEC adds programming and "principle-based physical validity" verifiers and outperforms LLM baselines. — [arXiv abs](https://arxiv.org/abs/2604.00149)
- **Yi Zhou, "From Paper to Program", arXiv:2604.04089 (Apr 2026; v7 Sep 2026).** Having LLMs write specifications first raised success rates:
  - DMRG: 16/16 specification–implementation pairings succeeded, against 6/13 direct attempts.
  - Pfaffian HFB→MPS conversion: four self-paired workflow attempts passed and no direct attempt did. — [arXiv abs](https://arxiv.org/abs/2604.04089)
- **AI-assisted rigorous many-body mathematics** (details in §1): Haldane-gap extensions (Grover 2611.07400; Tasaki 2610.34330 and 2609.34119); exotic half-chain factors (Liu–Yi 2610.24154); the SYK central limit theorem (2610.18440); Witten Lindbladians (2609.13121); free fermions in disguise (2611.26315). — see the links in §1.
- **Field theory and amplitudes**
  - The gluon amplitude formula was checked with Berends–Giele recursion, soft and Kleiss–Kuijf identities. — [arXiv:2602.12176 HTML](https://arxiv.org/html/2602.12176v1)
  - QFT-Bench: 75 graduate-to-research problems with symbolic verification. The strongest models exceed 95%, but show "reasoning unfaithfulness" — correct answers whose reasoning does not support them. — [arXiv:2610.33134](https://arxiv.org/abs/2610.33134)
  - SCALAR (actor–critic–judge on QFT and string-theory problems): multi-turn critique helps, but depends on the actor–critic pairing. — [arXiv:2605.06772](https://arxiv.org/abs/2605.06772)
  - DualityCert: a symbolic verifier for Seiberg-duality claims ('t Hooft anomaly matching, R-charges, central charges). Verifier-gated retries improve repair success by 7–8 percentage points. A certificate means "no tested inconsistency", not proof. — [arXiv:2607.23614](https://arxiv.org/abs/2607.23614)
  - Solver Agent: a ledger-tracked multi-agent system with independent verifier agents, used for F-theory uplift theorems. — [arXiv:2609.35958](https://arxiv.org/abs/2609.35958)
- **Earlier context**
  - TPBench: 57 problems in high-energy theory and cosmology; research-level problems "mostly unsolved" by early-2025 models. — [arXiv:2502.15815](https://arxiv.org/abs/2502.15815)
  - A follow-up on test-time scaling. — [arXiv:2506.20729](https://arxiv.org/abs/2506.20729)
  - Lu et al. argue that physics needs "physics-aware verification tools". — [arXiv:2506.06214](https://arxiv.org/abs/2506.06214)

### Inferences
- How results were checked:
  - an independent ED code (Bethe ansatz, DMRG);
  - consistency identities (amplitudes);
  - CAS with operator normal ordering (CMT-Benchmark);
  - exact or interval arithmetic plus Lean (Haldane gap).
- The most credible pipelines separate the agent from the checker: separate programs, separate verifier agents, frozen specifications.
- Mean-field and HF derivations were the first analytic many-body task automated (2024). Research-level CMT problems are still hard even for 2025 frontier models (30% best). The jump in 2026 came mainly in rigorous mathematical physics rather than in routine derivations.

### Gaps
- I found no LLM-agent paper deriving QEC thresholds through stat-mech mappings, such as random-bond Ising / Nishimori-line mappings. A targeted arXiv search returned only StatMechBench-v0.
- I found no systematic agent work on Schrieffer–Wolff or other high-order perturbative effective Hamiltonians.
- I found no 2025–2026 Kim-group follow-up specifically on autonomous HF agents beyond CMT-Benchmark.
- I could not determine whether the CritPt and CMT-Benchmark scores have been updated for 2026 models.

---

## 4. Conjecture generation and AlphaEvolve-style search for counterexamples or extremal constructions in quantum information

### Takeaway
Most AI "discoveries" in quantum information in 2026 are counterexamples, constructions or sharp constants. They were found in chat-style interaction with frontier models or by LLM-guided program evolution, and then made rigorous with exact certificates (Lean, SAT/DRAT, exact arithmetic). DeepMind's broad AlphaEvolve campaigns contained no quantum-information problems in my keyword scan. Formal proof search (AlphaProof Nexus) did produce a quantum-optics no-go theorem.

### Cited Findings
- **AlphaEvolve itself:**
  - white paper — [arXiv:2506.13131](https://arxiv.org/abs/2506.13131) (Jun 2025)
  - Georgiev, Gómez-Serrano, Tao, Wagner, "Mathematical exploration and discovery at scale": 67 problems in analysis, combinatorics, geometry and number theory — [arXiv:2511.02864](https://arxiv.org/abs/2511.02864) (Nov 2025)
  - My keyword scan of the latter found no quantum-information problems. Its "uncertainty principle" problem is Fourier-analytic. — [HTML](https://arxiv.org/html/2511.02864)
- **Krenn et al., GHZ no-go via AlphaProof Nexus, arXiv:2611.26300.** Formal proof search found an algebraic argument: a weighted sum of amplitudes whose coefficients factor vertex by vertex, combined with a row–column imbalance. This rules out monochromatic quantum graphs for even N ≥ 4 and D ≥ N. — [arXiv abs](https://arxiv.org/abs/2611.26300)
- **Baek, Won, Jon-Lark Kim, arXiv:2610.41183 (30 Oct 2026).**
  - *Discovery:* LLM-guided program evolution produced finite entanglement-assisted code constructions. Researchers generalized them to the family [[n,1,n−1;n−q−1]]_q.
  - *Verification:* Lean checks the algebraic certificate; ten nonexistence results come with checked DRAT proofs.
  - *Outcome:* closes 54 gaps in code tables.
  - *Generalization test:* validation-selected programs succeeded 46.1% and 31.6% of the time at unseen n = 13 and 15, against 19.5% and 3.1% for independently generated programs. — [arXiv abs](https://arxiv.org/abs/2610.41183)
- **He, Lu, Zeng, arXiv:2510.20728.**
  - A Lean-certified catalogue of 14,116 nonadditive codes with transversal diagonal gates, from which closed-form infinite families were extracted.
  - A ((6,4,2)) code implementing diag(1,1,1,i).
  - For ((7,2,3)) transversal T in the BD₁₆ setting: 10 of 12 candidates realized, 2 excluded by no-go proofs. — [arXiv abs](https://arxiv.org/abs/2510.20728)
- **Counterexamples found in chat with frontier models** (details in §1): quantum Hedetniemi (2609.20690), the EoF supporting functional (2608.27363), type II_∞/III_λ half-chain factors (2610.24154), the private-capacity superactivation example (2609.10520), and the explicit Tsirelson game (2610.10248).
- **Sharp constants and bounds:**
  - √2·min(n,m) for the trace norm against the injective norm (2608.06235);
  - (1/π)·ln n for ordered search (2609.34293);
  - (2p+1)/(2p+2) for QAOA on the ring (2606.29687);
  - reduced finite-size penalty in BB84 key rates via the Lean "discovery machine" (2610.17286). — see §1 links.
- **Entropy inequalities (classical):** Wong, Liu, Song, Gohari, Cheuk Ting Li, arXiv:2606.05729 (4 Jun 2026). Fine-tuned 0.6B models with guided beam search prove 85% of 60 held-out Shannon-type inequalities with 10–15 variables. For comparison, GPT-5.5 zero-shot solves 1.7% and Psitip 33.3%. — [arXiv abs](https://arxiv.org/abs/2606.05729)
- **OpenAI MUB claim N(6)=3** (Oct 2026). Computer-assisted, binary64, not interval arithmetic. For verification purposes it is closer to a strongly supported conjecture than to a certified theorem. — [arXiv:2610.17090 HTML](https://arxiv.org/html/2610.17090)
- **Search-based items, listed for completeness (out of scope):**
  - LLM agents plus symmetry to restrict syndrome search for distance-≥3 magic-state factories, "followed by deterministic solving and independent verification" — [arXiv:2610.06535](https://arxiv.org/abs/2610.06535)
  - LLM-guided bivariate-bicycle code search — [arXiv:2606.02418](https://arxiv.org/abs/2606.02418), [arXiv:2610.06623](https://arxiv.org/abs/2610.06623)

### Inferences
- The productive loop is: an LLM proposes candidates or programs, an exact checker accepts or rejects them, and humans generalize finite patterns into families and proofs. Validation-selected program evolution generalizes better than one-shot generation (Baek et al.).
- Finite-dimensional extremal problems are the natural targets: chromatic numbers, nonlocal-game values, small channels, code families. Their certificates are integers, rationals or SAT proofs that Lean or DRAT can check.

### Gaps
- Targeted arXiv searches found no LLM- or AlphaEvolve-driven discovery of quantum entropy inequalities (beyond SSA/DPI formalization), holographic entropy-cone facets, magic monotones or stabilizer-rank bounds, or Lieb–Robinson-type bounds. The searches combined "entropy inequalities", "magic"/"stabilizer Rényi"/"stabilizer rank" and "Lieb-Robinson" with LLM, GPT, Claude, Gemini or AlphaEvolve terms.

---

## 5. Agents that reproduce, check or find errors in published theory; research-level benchmarks (what frontier agents get right and wrong)

### Takeaway
Formalization, increasingly done by agents, has become the main way errors are found in quantum-information and physics papers:
- CHSH rigidity;
- the 2HDM stability theorem, which was invalidated;
- MIP*=RE side conditions;
- the generalized quantum Stein's lemma.

By contrast, LLMs acting as direct error-finders or end-to-end reproducers remain weak. SPOT reached at most 21% recall; PRBench's best score is 34%; CritPt's research-scale challenges were solved about 6–10% of the time by 2025 models. Benchmarks that grade with machine verification show agents are strong on well-posed calculations but brittle on intent, structure and long workflows.

### Cited Findings

**Errors found through formalization**
- **Tooby-Smith, arXiv:2603.08139 (9 Mar 2026).** Formalizing the 2006 Maniatis–von Manteuffel–Nachtmann–Nagel paper on 2HDM potential stability "revealed an error … invalidating their main theorem". Described as "the first non-trivial error in a physics paper found through formalization". — [arXiv abs](https://arxiv.org/abs/2603.08139)
- **Zhao & Yu, arXiv:2604.03884:** a gap in a CHSH-rigidity proof, with a dimension-2 counterexample. — [HTML](https://arxiv.org/html/2604.03884)
- **FormalFlow, arXiv:2609.19814:** corrected the side conditions and error bounds of the published low individual-degree-test proof (MIP*=RE). — [arXiv abs](https://arxiv.org/abs/2609.19814)
- **Meiburg et al., arXiv:2510.08672:** "rectified minor imprecisions" in the Hayashi–Yamasaki proof. — [arXiv abs](https://arxiv.org/abs/2510.08672)

**Errors in AI-generated proofs, and in audits of them**
- **Douglas, arXiv:2609.05555 v2:** replaced an unsupported step. — [arXiv abs](https://arxiv.org/abs/2609.05555)
- **Sienicki & Sienicki audit of OpenAI's ten results, arXiv:2608.14673.**
  - Chapter 8 remains the main reservation; one specialist review called for "major revision of compressed analytic passages".
  - The Chapter 6 (quantum parallel repetition) objection was an artifact of PDF extraction.
  - The audit cites OpenAI's February 2026 "First Proof" exercise, in which "one attempt initially considered likely correct was reassessed after expert and community criticism". — [HTML](https://arxiv.org/html/2608.14673)
- **Xie survey, arXiv:2610.17090:** its source comparisons "correct mathematical statements, hypotheses and counts" in OpenAI's release. — [arXiv abs](https://arxiv.org/abs/2610.17090)

**LLMs as direct verifiers or reviewers**
- **SPOT (Son et al.), arXiv:2505.11855 (May 2025).** 83 published papers containing 91 errors serious enough for errata or retraction. No model surpassed 21.1% recall or 6.1% precision (o3 best). Models rarely rediscovered the same errors across runs and made "student-level misconceptions". — [arXiv abs](https://arxiv.org/abs/2505.11855)
- **Zubeldia, Bolliet, Villaescusa-Navarro, Villanueva-Domingo, Cranmer, arXiv:2611.17052 (12 Nov 2026; NeurIPS 2026).**
  - 357 injected errors across 30 arXiv papers (astrophysics and ML).
  - With the same prompt, three frontier LLMs differ by about 50 percentage points in detection.
  - Their reviewer "Skepthical" checks mathematics "with computer algebra" and numerics "with executed code". Claude Code improves over a single Claude call on maths and numerics. — [arXiv abs](https://arxiv.org/abs/2611.17052)

**Reproduction of theory and computational papers**
- **Yuxuan Du et al., arXiv:2611.02833 (2 Nov 2026).** 59 computational quantum papers; five agentic systems and six agent harnesses evaluated. "Current AI scientists often fail to reproduce published results."
  - Their "Residual Agents" raise DeepSeek V4.1 Flash from 68% to 85%, "matching Claude Code with Opus 5 at a third of its cost". — [arXiv abs](https://arxiv.org/abs/2611.02833)
- **PRBench, arXiv:2603.27646 (Mar 2026).** 30 end-to-end physics paper-reproduction tasks; the best agent (Codex with GPT-5.3-Codex) scores 34%. — [arXiv abs](https://arxiv.org/abs/2603.27646)
- **PHYRD, arXiv:2610.32575 (Oct 2026).** 63 "Reproduce" problems (deriving published results) and 43 "Discover" problems. Whole-problem accuracy is 49–68% on Reproduce. — [arXiv abs](https://arxiv.org/abs/2610.32575)
- **PRL-Bench, arXiv:2604.15411 (Apr 2026).** 100 PRL papers since Aug 2025, including quantum information; the best overall score is below 50. — [arXiv abs](https://arxiv.org/abs/2604.15411)
- **MerLean** formally reproduced three quantum-computing papers in Lean. — [arXiv:2602.16554](https://arxiv.org/abs/2602.16554)
- **The Kim group** reproduced HF derivations from 15 papers. — [arXiv:2403.03154](https://arxiv.org/abs/2403.03154)

**Research-level benchmarks**
- **CritPt, arXiv:2509.26574 (Sep 2025; v4 May 2026).** 71 composite research challenges and 190 checkpoint tasks written by more than 50 researchers, with machine-verifiable answers. Best base-model accuracy is 5.7% (GPT-5 high), about 10% with coding tools. — [arXiv abs](https://arxiv.org/abs/2509.26574)
- **TPBench** — [arXiv:2502.15815](https://arxiv.org/abs/2502.15815)
- **CMT-Benchmark** (best 30%) — [arXiv:2510.05228](https://arxiv.org/abs/2510.05228)
- **QFT-Bench** (more than 95% for the strongest models, with unfaithful reasoning) — [arXiv:2610.33134](https://arxiv.org/abs/2610.33134)
- **FrontierScience (OpenAI):** Olympiad and PhD-level "Research" tracks, with rubric grading — [arXiv:2601.21165](https://arxiv.org/abs/2601.21165)
- **PhysicsEval** — [arXiv:2508.00079](https://arxiv.org/abs/2508.00079)
- **StatMechBench-v0** — [arXiv:2607.26367](https://arxiv.org/abs/2607.26367)
- **Lean-QuantumAlg/QIT-Bench** (about 60/100) — [arXiv:2607.21533](https://arxiv.org/abs/2607.21533)
- **AxQM** — [arXiv:2609.05157](https://arxiv.org/abs/2609.05157)

### Inferences
- **Frontier agents get right:**
  - closed-form QFT and textbook-to-graduate calculations (QFT-Bench);
  - HF derivations once intermediate steps are corrected;
  - implementations, when given specifications (DMRG 16/16);
  - long formal proofs under a blueprint (126k lines for FormalFlow; TNLean);
  - creative constructions and counterexamples in finite-dimensional quantum information (§1).
- **Frontier agents get wrong:**
  - full research-scale workflows (CritPt, PRBench, PRL-Bench, Du et al.);
  - recognizing structure, as opposed to matching numbers (StatMechBench);
  - faithful reasoning (QFT-Bench);
  - stating theorems faithfully: quietly added hypotheses (CHSH), semantically vacuous tasks that the kernel accepts (AxQM), statement drift (FormalFlow);
  - reliably detecting errors in papers (SPOT).
- Error-finding works best when it is a side effect of formalization or CAS/code re-execution, not when an LLM reads the text and judges.

### Gaps
- I found no benchmark specific to quantum-information theory that measures agents finding errors in published QI papers. SPOT and the reviewer audit cover general and astrophysics/ML papers.
- I found no systematic, field-wide formalization audit of a corpus of QI papers that reports error rates. Tooby-Smith raises the question but gives no estimate.

---

## 6. Which theory tasks in quantum information and many-body physics have cheaply checkable outputs (suited to autonomous agents)?

### Takeaway
The work that agents could do autonomously in 2026 is the work with an independent, inexpensive checker:
- Lean-kernel-checkable finite certificates (codes, games, graphs, chromatic numbers, small channels);
- exact or interval-arithmetic diagonalization of small systems that seeds rigorous bounds (spectral gaps);
- an independent ED code (Bethe ansatz, DMRG phase diagrams);
- CAS or recursion identities (amplitudes, operator algebra with normal ordering);
- frozen Lean QHL specifications (quantum-algorithm correctness and cost).

Purely numerical agreement is not a sufficient check on its own.

### Cited Findings
- **Lean kernel plus finite certificates:**
  - code distances: An & Yang, kernel-sized certificates — [arXiv:2610.03214](https://arxiv.org/abs/2610.03214); QECLean — [arXiv:2611.21236](https://arxiv.org/abs/2611.21236); Lean-QEC SAT replay — [arXiv:2605.16523](https://arxiv.org/abs/2605.16523)
  - EAQECC certificates and DRAT nonexistence proofs — [arXiv:2610.41183](https://arxiv.org/abs/2610.41183)
  - the explicit Tsirelson game (Lean/Mathlib) — [arXiv:2610.10248](https://arxiv.org/abs/2610.10248)
  - Hedetniemi exact integer certificates — [arXiv:2609.20690](https://arxiv.org/abs/2609.20690)
- **Certified small-system numerics feeding rigorous bounds:** Haldane-gap seeds use rational or interval eigenvalue enclosures (Ostrowski–Weyl), then a doubling argument extends them to all L. — [arXiv:2611.11049 PDF](https://arxiv.org/pdf/2611.11049)
- **Weaker computational standard:** OpenAI's MUB computation used binary64 with hand-proved margins, and "neither interval arithmetic nor exact algebra". — [arXiv:2610.17090 HTML](https://arxiv.org/html/2610.17090)
- **Independent ED:**
  - Bethe-ansatz results "checked against exact diagonalization (performed by separate programs)" — [arXiv:2603.29932](https://arxiv.org/abs/2603.29932)
  - DMRG cross-checked against ED — [arXiv:2610.32001](https://arxiv.org/abs/2610.32001)
- **CAS / symbolic checks:**
  - normal-ordered operator grading — [arXiv:2510.05228](https://arxiv.org/abs/2510.05228)
  - Berends–Giele, soft and Kleiss–Kuijf checks — [arXiv:2602.12176 HTML](https://arxiv.org/html/2602.12176v1)
  - anomaly matching in DualityCert — [arXiv:2607.23614](https://arxiv.org/abs/2607.23614)
  - CAS-verified review (Skepthical) — [arXiv:2611.17052](https://arxiv.org/abs/2611.17052)
- **Formal specifications for algorithms:**
  - DeepMind's frozen QHL specifications, with an independent verifier — [arXiv:2611.09822 HTML](https://arxiv.org/html/2611.09822)
  - Lean-QuantumAlg contracts (less than 1% of code needs statement review) — [arXiv:2610.30216](https://arxiv.org/abs/2610.30216)
  - reduction of the FGG conjecture to a single Lean statement — [arXiv:2606.29687](https://arxiv.org/abs/2606.29687)
- **Failure modes when the checker is weak:**
  - agents pass numerical checks while misidentifying structure — [arXiv:2607.26367](https://arxiv.org/abs/2607.26367)
  - Lean accepts semantically vacuous statements — [arXiv:2609.05157 HTML](https://arxiv.org/html/2609.05157)
  - LLMs add hypotheses to make proofs trivial — [arXiv:2604.03884 HTML](https://arxiv.org/html/2604.03884)

### Inferences
Tasks ranked from cheapest to most expensive to check:

| Task | Checker | Notes |
|---|---|---|
| (a) Finite combinatorial or algebraic objects: code parameters, nonlocal-game values with exact classical value, quantum or classical chromatic numbers, small channel examples | integer or rational certificates checked in Lean, SAT/DRAT | Seconds to hours |
| (b) Statements reducible to finite-dimensional operator inequalities with closed-form test cases, e.g. EoF via Wootters' formula, entropy inequalities | random-state numerics (counterexample search), then a Lean proof using Lean-Quantum's operator-monotonicity infrastructure | |
| (c) Spectral-gap and correlation-decay bounds from finite-size criteria | ED with interval arithmetic for seeds, plus a formal proof of the general criterion | The Haldane template |
| (d) Integrability claims: Bethe ansatz, free fermions in disguise | ED spectra of small chains, plus Yang–Baxter / commutation identities by CAS | |
| (e) Mean-field / HF derivations and effective Hamiltonians | CAS with normal ordering, plus numerical self-consistency against ED on small clusters | |
| (f) Quantum-algorithm correctness and resources | Lean QHL specifications | |
| (g) Asymptotic analytic theorems, e.g. parallel repetition, area laws | Full formalization still needed | Expensive; human review of definitions remains the bottleneck |

### Gaps
- I found no study that measures the cost (compute or human hours) of each checking modality across these task classes. The ranking above is inferred from the cited examples.

---

## Gaps and opportunities

### Takeaway
The highest-value openings combine (i) a problem class where 2026 frontier models have already shown they can produce the core argument, with (ii) a cheap independent checker that now exists:
- Lean QI and QEC libraries;
- certified-ED gap criteria;
- SAT/DRAT certificate pipelines;
- QHL specifications;
- CAS identity checks.

### Cited Findings
These are reusable checkers and templates that agents could build on:
- **CertifiedGaps** (exact-arithmetic seeds + Lean criterion; already applied to seven other chains and ladders) — [arXiv:2611.11049](https://arxiv.org/abs/2611.11049)
- **QI libraries:** Lean-Quantum (DPI, SSA, operator monotonicity) — [arXiv:2607.05492](https://arxiv.org/abs/2607.05492); Lean-QIT (capacities) — [arXiv:2607.09632](https://arxiv.org/abs/2607.09632)
- **QEC libraries:** QECert (qLDPC, BPT bounds, energy barriers) — [arXiv:2611.05047](https://arxiv.org/abs/2611.05047); QECLean (code combinators) — [arXiv:2611.21236](https://arxiv.org/abs/2611.21236)
- **Quantum programs:** DeepMind QHL library — [arXiv:2611.09822](https://arxiv.org/abs/2611.09822)
- **Formalization workflows:** MerLean — [arXiv:2602.16554](https://arxiv.org/abs/2602.16554); FormalFlow — [arXiv:2609.19814](https://arxiv.org/abs/2609.19814)
- **Open verification debts flagged in the literature:**
  - Haldane gap and 2D area law not formalized; N(6)=3 MUB result only in binary64 — [arXiv:2610.17090 HTML](https://arxiv.org/html/2610.17090)
  - Douglas's I3322 analytic proof not formalized — [arXiv:2609.05555](https://arxiv.org/abs/2609.05555)
  - Tasaki's uniqueness result does not cover non-translation-invariant states — [arXiv:2610.34330](https://arxiv.org/abs/2610.34330)
  - QECLean-style analytic distances are still needed beyond the SAT scale — [arXiv:2605.16523](https://arxiv.org/abs/2605.16523)
- **Weak spots of current proof agents:** quantum simulation, quantum learning, QI measures, entanglement theory — [arXiv:2607.21533](https://arxiv.org/abs/2607.21533)

### Inferences (opportunities, with how outputs would be verified)
1. **Haldane-type certified gaps for new models.** Spin-1 XXZ and single-ion anisotropy, bilinear-biquadratic chains away from AKLT, ladders, odd rings, and the 2D extensions Grover sketches.
   - Agents generate MPS trial states and transfer-matrix seeds and run the purity bootstrap.
   - *Verification:* interval or rational ED enclosures plus a Lean-checked criterion (CertifiedGaps pattern).
2. **Close the OpenAI verification debts.**
   - Re-run the N(6)=3 MUB computation in interval or exact arithmetic and formalize it.
   - Formalize the Haldane-gap criterion end to end, and the 2D area-law proof.
   - Formalize the parallel-repetition improvements (exponent 13 → 3 claimed).
   - *Verification:* Lean kernel, plus independent human statement review.
3. **Agentic formalization audit of high-use QI lemmas** — self-testing and rigidity, postselection and de Finetti, decoupling, continuity bounds. These are the kinds of places where gaps were found (CHSH, MIP*=RE, GQSL, 2HDM).
   - *Verification:* Lean, with a LaTeX back-translation for review (MerLean) and blueprint review.
   - *Expected output:* errata, sharpened constants, improved exponents.
4. **Discovering quantum entropy inequalities and monotones.**
   - LLM or program-evolution proposals for new constrained inequalities, holographic-cone facets, magic-monotone bounds, and continuity constants.
   - *Verification:* random-state and SDP counterexample search for falsification, then Lean proofs built on Lean-Quantum's operator infrastructure. No LLM work was found here (§4 gaps).
5. **Extremal finite constructions in nonlocality and quantum graph theory** — Bell functionals with non-closure, quantum chromatic gaps, state-discrimination profiles.
   - *Verification:* exact integer or rational certificates plus Lean, following the Tsirelson-game and Hedetniemi patterns. NPA-type upper bounds would need rationally rounded SDP certificates; I found no examples of this in the cited works.
6. **Integrable and free-fermion model discovery.**
   - Agents propose Hamiltonians, then derive Bethe-ansatz or free-fermion solutions and frustration-graph criteria.
   - *Verification:* independent ED spectra, CAS checks of Yang–Baxter and commutation relations, and Lean for the graph-theoretic criteria (Pozsgay–Vona; Projansky et al.).
7. **Mapping QEC thresholds onto stat-mech models**, e.g. random-bond Ising, Nishimori line, gauge-removable disorder.
   - *Verification:* exact partition functions or transfer matrices on small lattices, symbolic gauge-invariance checks, and structural invariants rather than numbers alone (StatMechBench lesson).
   - No agent work was found here (§3 gaps).
8. **Mean-field, HF and effective-Hamiltonian derivations at scale** — moiré HF, high-order Schrieffer–Wolff.
   - *Verification:* CAS with normal ordering (CMT-Benchmark grading), plus comparison against ED on small clusters.
9. **Analytic QEC distance and bound proofs beyond the SAT scale** — code families via combinators, BPT-type locality bounds, energy barriers.
   - *Verification:* QECLean/QECert-style Lean proofs, with kernel-sized certificates (An & Yang) for specific instances.
10. **Verified quantum-algorithm resource reductions** — Hamiltonian simulation and arithmetic, following DeepMind's QHL workflow.
    - *Verification:* frozen Lean specifications of correctness and T-count, with an independent comparator.
11. **Reproduction plus error-finding sweeps over QI theory papers**, with computational claims re-executed (Du et al.'s platform) and mathematical claims CAS- or Lean-checked (Skepthical / MerLean).
    - *Verification:* the original numbers or figures, plus formal statement checks.
    - This would fill the missing QI-specific error-detection benchmark.

### Gaps
- Whether these opportunities are already being pursued privately (e.g. inside OpenAI, DeepMind or Anthropic) is unknowable from public sources.
- Anthropic's only primary science report I located is on mathematics, not quantum physics. "Learning more about Claude's mathematical capabilities" (Aug 2026) reports a Riemann-zeta critical-line proportion bound raised from 41.6% to 67.2%, about 31M output tokens, and a Lean formalization checked by comparator (search summary, not fetched). — [Anthropic](https://www.anthropic.com/research/riemann-zeta)
- Google's "Accelerating Scientific Research with Gemini" case studies contained no quantum-information case in a keyword scan. — [arXiv:2602.03837](https://arxiv.org/abs/2602.03837)
- Aletheia (DeepMind) targets mathematics only. — [arXiv:2602.10177](https://arxiv.org/abs/2602.10177)
- I found no public cost-per-result data for quantum-information results, except AlphaProof Nexus's "few hundred dollars" per Erdős problem. — [arXiv:2605.22763](https://arxiv.org/abs/2605.22763)
