# Agentic parent-Hamiltonian reconstruction: status check and literature map (as of 2026-11-22)

*Method note. Phase 1 used WebSearch (standard and extended). In phase 2, after the network was opened, I checked metadata (authors, dates, journal references, DOIs, abstracts) for about 130 arXiv IDs through the arXiv API with curl. I also ran six boolean arXiv API searches: newest "parent Hamiltonian" papers; Hamiltonian learning/reconstruction combined with LLM, agent or AI-scientist terms; AlphaEvolve/FunSearch/LLM-guided evolution combined with quantum terms; single-eigenstate and covariance-matrix follow-ups; and title checks. I read the **full texts** (arXiv HTML) of the two closest works, SciExplorer (arXiv:2509.24978v6) and Kleidermacher et al. (arXiv:2611.00620v1). For all other items, details come from abstracts and API metadata, plus search snippets where marked. WebFetch could not reach arxiv.org, alphaxiv.org or pith.science (hosts not in the session's network allowlist), so all arXiv reading was done with curl. "Title only" means the existence and title were confirmed but the content was not checked.*

---

## Q1. Has "agentic" parent-Hamiltonian reconstruction been done (2023 to Nov 2026)? Who, at what system sizes, with which verifier?

### Takeaway
**Partially.** No paper, preprint or announcement was found in which an LLM agent, AI-scientist pipeline or LLM-guided program/evolutionary search (FunSearch/AlphaEvolve-style) starts from a *given target state or eigenstate* (or its correlation matrix or entanglement data) and outputs a local parent Hamiltonian certified by a parent-Hamiltonian verifier (energy variance ⟨H²⟩−⟨H⟩², null space of the correlation matrix, spectral gap).

Two LLM-agent works come closest:
- **SciExplorer** (Nägele & Marquardt, *PRX* 16, 031002 (2026)). A GPT-5-based agent infers spin-1/2 Hamiltonians of *unknown* systems from ground-state expectation values of operators it defines itself. Systems are 9–10 spins (up to about N=12 allowed). It is scored by per-spin fidelity with the true ground state, and it checks hypotheses internally with an ED forward model it writes itself.
- **Kleidermacher et al.** (arXiv:2611.00620, Oct 2026). A Claude Code (Sonnet 4.6) agent proposes and ranks Hamiltonian and Lindbladian models for *experimental* few-body cavity-QED data, using BIC/AICc, nested-sampling Bayesian evidence and Bayes factors.

Neither paper mentions parent Hamiltonians, the covariance-matrix method or Qi–Ranard. In LLM work, parent Hamiltonians appear only as benchmark generators (Farias et al., arXiv:2611.17246).

### Cited Findings
- **SciExplorer: the closest match.** Maximilian Nägele and Florian Marquardt, "Agentic Exploration of Physics Models", arXiv:2509.24978 (29 Sep 2025; v6 8 Jul 2026), *Phys. Rev. X* **16**, 031002 (2026), DOI 10.1103/xnqc-q6nt. — [arXiv abs](https://arxiv.org/abs/2509.24978); [PRX](https://journals.aps.org/prx/abstract/10.1103/xnqc-q6nt)
  - *Architecture.* A general LLM tool-use agent with no domain-specific blueprints. Tools let it run experiments, code, plot and save results, and the code is public. The default LLM backbone is **GPT-5**. — [full text](https://arxiv.org/html/2509.24978v6); [GitHub](https://github.com/MaxNaeg/SciExplorer)
  - *Ground-state task.* The agent defines Hermitian operators and receives their expectation values in the ground state of a hidden Hamiltonian. It is told only that the system is 1D or 2D and made of spins. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Systems in the ground-state setting.*
    - 1D transverse-field Ising, 10 spins
    - 2D Heisenberg, 9 spins
    - 1D cluster-Ising with three-body terms, 10 spins
    - an arbitrary translationally invariant Hamiltonian on 10 spins
    - variants with tunable field or coupling and a choosable N ("do not choose values larger than about N=12")
    - — [full text](https://arxiv.org/html/2509.24978v6)
  - *Systems in the dynamics setting.* An arbitrary Hamiltonian on 3 spins, the 10-spin Heisenberg and transverse-field Ising chains, and a partial-observation variant in which only spins 0 and 1 of a 10-spin chain can be initialized and observed. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Benchmark verifiers.* Ground-state tasks use fidelity per spin, |⟨ψ_agent|ψ_true⟩|^{2/N}. Dynamics tasks use the normalized Hilbert–Schmidt overlap of the traceless parts of H_true and H_agent, which penalizes rescaling. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Agent-internal verification.* The agent first probes symmetry and structure: single-spin values, nearest-neighbour correlators, and chirality-sensitive correlators such as ⟨S^x_j S^y_{j+1} − S^y_j S^x_{j+1}⟩. It deduces open boundaries, translational invariance and which coupling families are absent. It then simulates candidate Hamiltonians to accept or reject them, while keeping unrefuted options such as longer-range couplings. In the cluster-Ising run it wrote sparse ED for N=10 and concluded that "the data uniquely select a cluster–Ising transverse-field model with an A-linear cluster term on an open chain of N=10 spins". — [full text](https://arxiv.org/html/2509.24978v6)
  - *Noise.* With finite measurement shots, results become less consistent, but the agent recovers a highly accurate model in at least one of six attempts per noise level. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Failure modes.* (i) Too few or too similar experiments: in the ground-state setting the agent "may fail to consider enough expectation values, leading to an incorrect model that nevertheless fits all observed measurement results". (ii) Correct structure with wrong parameters, often sign errors or missing factors of two, without considering an overall scale or sign. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Dependence on the LLM.* GPT-5 "dramatically outperforms" Gemini 2.5 Pro, which fails even with the full framework. GPT-5 nano always fails, GPT-5 mini succeeds in 1 of 5 attempts, and open-source models struggle even on easier mechanical tasks. A GPT-5 exploration costs about 2 USD. — [full text](https://arxiv.org/html/2509.24978v6)
  - *Not present.* A full-text search finds no mention of parent Hamiltonians, covariance or correlation matrices, energy variance, null spaces, Qi–Ranard, Chertkov–Clark or Bairey et al. — [full text](https://arxiv.org/html/2509.24978v6)
- **Hannah C. Kleidermacher, Bennet Windt, Miguel Bello, James Zou and Jelena Vučković, "Discovering Models for Experimental Quantum Systems with Agentic AI", arXiv:2611.00620 (31 Oct 2026, v1).** — [arXiv abs](https://arxiv.org/abs/2611.00620)
  - *Architecture.* The agent is implemented in Claude Code (Claude Sonnet 4.6) with four Claude Agent Skills: experimental-background, literature-search (fetches arxiv.org/html full texts), explore-data and validate-phase (a formal gate between phases).
  - *Workflow.* Three stages: (1) literature search and a research plan listing candidate Hamiltonian/Lindbladian models, evaluation methods and a ranking strategy, followed by self-revision; (2) generation of a Python tool library; (3) phased analysis execution with monitoring and autonomous code fixes. An optional follow-up loop is available. The only experiment-specific input is one document written by a human. — [full text](https://arxiv.org/html/2611.00620v1)
  - *Case 1.* Two SnV⁻ centers in a diamond photonic-crystal cavity, with 110 lifetime traces.
    - The agent found two emitters (the human analysis had found one dominant emitter, a fact hidden from the agent).
    - It derived an upper bound N ≤ 4 on the number of emitters from photoluminescence peaks.
    - "Drawing on prior work in quantum Hamiltonian learning (Wang 2017; Gentile 2021; Poteshman 2025)", it screened models by BIC and then ranked them by nested-sampling evidence and Bayes factors, filtered by validity checks.
    - **The three- and four-emitter models had *higher* Bayesian evidence, but the agent disqualified them because the fitted lifetimes were shorter than the instrument response time.**
    - — [full text](https://arxiv.org/html/2611.00620v1)
  - *Case 2.* A single V_Si center in a SiC resonator. The agent tested 32 candidate models and ran an "effect-importance" analysis based on ΔAICc. Its top model has two optical transitions plus Lorentzian spectral diffusion, with g/2π = 171 ± 1 MHz. It found that mode backscattering is degenerate with detuning and so unsupported. — [full text](https://arxiv.org/html/2611.00620v1)
  - *Case 3.* About 10 V_Si emitters in a whispering-gallery-mode disk. The agent fitted disordered open two-mode Tavis–Cummings Lindblad models with N = 3–6 emitters, computing g⁽²⁾ through the quantum regression theorem and fitting with CMA-ES. It found BIC_{N=5} < BIC_{N=4} < BIC_{N=3} and two comparably good N=4 solutions. The run used more than two days of compute, and the agent wrongly dismissed spectral diffusion on limited evidence. — [full text](https://arxiv.org/html/2611.00620v1)
  - *Reproducibility of plans.* Over 30 independently generated research plans, the mean pairwise F1 for candidate models was 0.64 (SnV⁻) and 0.84 (V_Si) with Sonnet 4.6, and 0.61 and 0.86 with Claude Fable 5. For analysis strategies it was 0.58/0.66 and 0.66/0.68. Agreement on retrieved literature was F1 ≈ 0.30–0.41. — [full text](https://arxiv.org/html/2611.00620v1)
  - *Reproducibility of parameters.* Over 35 analysis runs from one fixed plan, the winning-model parameters had relative interquartile ranges below 1.8% (SnV⁻) and below 0.06% (V_Si). — [full text](https://arxiv.org/html/2611.00620v1)
  - *How the authors frame it.* Autonomous Hamiltonian/Lindbladian model discovery from quantum experiments "remains largely unexplored". They cite SciExplorer and non-LLM structure-learning work such as Gentile et al. 2021, Bakshi et al. 2024, Wallace et al. 2025, Poteshman et al. 2025 and Ivashkov et al. 2026. — [full text](https://arxiv.org/html/2611.00620v1)
- **Parent Hamiltonians as an LLM benchmark (the reverse direction).** Tiago de Souza Farias, Caio Moraes Porto, Matheus da Silva Fonseca, Nicolas Armando Cabrera Carpi, Anderson de Araujo Oliveira and Celso Jorge Villas-Boas, arXiv:2611.17246 (12 Nov 2026).
  - *Setup.* A frozen LLM writes complete parameter vectors for variational circuits from symbolic descriptions of the circuit and Hamiltonian.
  - *Benchmarks.* General random Hamiltonians, "structured parent Hamiltonians" and molecules.
  - *Results.* For general Hamiltonians the errors are larger than with Adam, and feedback helps little. — [arXiv abs](https://arxiv.org/abs/2611.17246) A search snippet adds that the LLM's direct energy predictions can fall below the spectral lower bound. — [arXiv HTML](https://arxiv.org/html/2611.17246)
- **Systematic arXiv checks on 2026-11-22.**
  - (i) A boolean query combining Hamiltonian learning/reconstruction/inference/discovery, parent Hamiltonian, model discovery, Hamiltonian identification and effective Hamiltonian with LLM(s), agentic, AI agent or AI scientist returned 50 hits. The only LLM-based quantum-physics papers among them were Kleidermacher et al. and Farias et al.; the rest were LLM model-discovery work in other fields or agents not built on LLMs. SciExplorer's abstract does not use these phrases. — [2611.00620](https://arxiv.org/abs/2611.00620); [2611.17246](https://arxiv.org/abs/2611.17246)
  - (ii) Among the ~70 newest of 264 "parent Hamiltonian" papers (Mar 2025 to Nov 2026), the only one involving an LLM was Farias et al. — [arXiv API query](https://export.arxiv.org/api/query?search_query=abs:%22parent%20Hamiltonian%22&sortBy=submittedDate&sortOrder=descending&max_results=70)
  - (iii) A query for AlphaEvolve, FunSearch, OpenEvolve, ShinkaEvolve or LLM-guided program evolution combined with quantum terms returned 9 hits. All concern codes, circuits, encodings, compilation or Trotter circuits, and none concerns Hamiltonian reconstruction (Q2d). — [arXiv abs 2606.02418](https://arxiv.org/abs/2606.02418); [arXiv abs 2510.19550](https://arxiv.org/abs/2510.19550)
- **LLM-driven symbolic regression has not been applied to quantum Hamiltonians.** LLM-SR (ICLR 2025) and LLM-SRBench (ICML 2025; 239 problems, best symbolic accuracy 31.5%) contain no quantum-Hamiltonian tasks according to search. The closest quantum symbolic-regression work, QDisc (arXiv:2604.16015), uses symbolic regression guided by a VAE and no LLM. — [LLM-SR](https://arxiv.org/abs/2404.18400); [LLM-SRBench PMLR](https://proceedings.mlr.press/v267/shojaee25a.html); [QDisc](https://arxiv.org/abs/2604.16015)

### Inferences
- SciExplorer is, in substance, Hamiltonian learning from ground-state data. That is the parent-Hamiltonian inverse problem posed as an interactive game, but solved by hypothesize-simulate-compare with ED rather than by covariance-matrix linear algebra, and scored against hidden ground truth. Two changes would turn it into a parent-Hamiltonian reconstruction agent usable on states with no known Hamiltonian: (a) give it the correlation-matrix null-space solver as a tool, and (b) replace fidelity with the hidden ground state by intrinsic certificates (zero variance, null-space dimension, DMRG/ED ground-state check and gap).
- The two closest works use opposite verifier philosophies. SciExplorer uses ground truth from benchmarks; Kleidermacher et al. use statistical model selection on real data. Parent-Hamiltonian reconstruction offers a third kind: cheap, exact algebraic certificates. That makes it unusually well suited to LLM search with the verifier in the loop.
- "Partially" is the right label. "No" would ignore SciExplorer and Kleidermacher et al. "Yes" would overstate them: neither starts from a given target state, and neither uses parent-Hamiltonian certificates.

### Gaps
- Non-arXiv venues (workshop papers, talks, company blogs) were not exhaustively searched. KITP's 2026 "AI for Quantum Matter" program shows community interest; its contents were not checked. — [KITP aiqmatter26](https://www.kitp.ucsb.edu/activities/aiqmatter26)
- SciExplorer reports per-system success distributions as figures; exact success rates for each quantum task were not transcribed here.

---

## Q2. Map (1): agentic and LLM approaches to Hamiltonian learning and model discovery (newest first)

### Takeaway
LLM agents for quantum physics multiplied in 2025–2026, but most tackle forward simulation, reproducing published results, laboratory automation, circuits or codes. Only SciExplorer and Kleidermacher et al. infer Hamiltonian or Lindbladian models from data; a third, the NV-centre agent of Isogawa et al., forms Hamiltonian-level hypotheses during experiments. LLM-guided evolutionary search in quantum (AlphaEvolve and similar) is now common for codes, circuits and encodings, never for Hamiltonian reconstruction. Generic LLM model-discovery frameworks outside quantum (Bayesian experiment design with LLM proposals; benchmarks) provide reusable templates. Across all of these, success tracks the presence of a hard, automatic verifier.

### Cited Findings
**(a) LLM agents that infer Hamiltonians or models from quantum data**
- Kleidermacher et al. 2026 (arXiv:2611.00620) and Nägele & Marquardt (SciExplorer, *PRX* 16, 031002), detailed in Q1. — [2611.00620](https://arxiv.org/abs/2611.00620); [2509.24978](https://arxiv.org/abs/2509.24978)
- Takuya Isogawa, Ryotaro Okabe, Nutdech Phadetsuwannukun, Mingda Li and Paola Cappellaro, "Agentic AI for Scientific Reasoning in Autonomous Quantum Sensing Experiments", arXiv:2607.25145 (27 Jul 2026).
  - *Autonomous run.* An LLM agent ran NV-centre experiments: it selected an NV, calibrated the resonance, measured T₂* by Ramsey, and added CPMG to check a weak feature possibly due to a nearby ¹³C.
  - *Offline benchmarks.* GPT-5.4, GPT-5.5 and GPT-5.6 Sol were evaluated. Requiring an expected-signal calculation kept false positives low; higher reasoning effort without it produced *more* false-positive resonance judgments.
  - *Design lesson.* The agent forms hypotheses while deterministic code controls the hardware. — [arXiv](https://arxiv.org/abs/2607.25145)

**(b) LLM agents for forward many-body simulation (candidate verifier back-ends)**
- Jannes Nys, Andrew Jreissaty and Juan Carrasquilla, "QuantumMatterExplorer: An Autonomous Agent for the Discovery of Quantum Many-Body Phases of Matter", arXiv:2610.32001 (22 Oct 2026).
  - *Architecture.* Input is only a Hamiltonian and its parameter ranges; output is a phase diagram. An LLM reasons over a long horizon with a persistent written plan, and DMRG results are cross-checked against independent ED. The paper includes a 13-page research note generated by the agent.
  - *Details from search snippets.* 222 DMRG runs over 158 (J₂, D, N) points, Néel, Haldane and dimerized phases found, and "Qwen3.6-27B" reported as the best of four backends. — [arXiv](https://arxiv.org/abs/2610.32001); [GitHub](https://github.com/jwnys/quantum_matter_explorer)
- Ken Deng, Xiangfei Wang, Guijing Duan, Chen Mo, Junkun Huang, Runqing Zhang, Ling Qian, Zhiguo Huang, Jize Han and Di Luo, "Towards Verifiable and Self-Correcting AI Physicists for Quantum Many-Body Simulations", arXiv:2604.00149 (31 Mar 2026; v2 10 May 2026). QMP-Bench has 100 research-level tasks taken from 21 journal articles. PhysVEC is a multi-agent framework with self-verification and error correction. — [arXiv](https://arxiv.org/abs/2604.00149)
- Weitang Li, Jiajun Ren, Lixue Cheng and Cunxi Gong, "Autonomous Quantum Simulation through Large Language Model Agents", arXiv:2601.10194 (15 Jan 2026). LLM agents run tensor-network simulations with about 90% success. Three configurations are compared: baseline, single agent with in-context learning, and multi-agent with in-context learning. — [arXiv](https://arxiv.org/abs/2601.10194)
- Ignacio Gustin, Luis Mantilla Calderón, Juan B. Pérez-Sánchez, Jérôme F. Gonthier, Yuma Nakamura, Karthik Panicker, Manav Ramprasad, Zijian Zhang, Yunheng Zou, Varinia Bernales and Alán Aspuru-Guzik, "El Agente Cuántico: Automating quantum simulations", arXiv:2512.18847 (Dec 2025; v3 Mar 2026; published in *Rep. Prog. Phys.*). A multi-agent system turns natural-language intent into validated computations across quantum-software frameworks. — [arXiv](https://arxiv.org/abs/2512.18847); [IOP](https://iopscience.iop.org/article/10.1088/1361-6633/ae8933)
- Yuxuan Du, Hanqing Zhao, Xinbiao Wang, Yiming Lu, Kaining Zhang and Dacheng Tao, "Reimagining computational quantum research with AI scientists", arXiv:2611.02833 (2 Nov 2026). OmniQ is an expert-audited environment that converts 59 papers into reproduction tasks, spanning classical simulation, quantum algorithms and "quantum-property characterization". AI scientists often fail to reproduce the papers. — [arXiv](https://arxiv.org/abs/2611.02833)
- Jue Xu, Yue Wang and Qi Zhao, "QSim-Arena: an evolving developer for quantum simulation advantage", arXiv:2611.15000 (11 Nov 2026). Agents run a closed verification-and-benchmarking loop over specifications of Hamiltonians, initial states and accuracy targets. — [arXiv](https://arxiv.org/abs/2611.15000)

**(c) LLM agents for experiments and ideas (context)**
- Sören Arlt, Xuemei Gu and Mario Krenn, AI-Mandel, arXiv:2511.11752 (13 Nov 2025). Agents generate and implement quantum-optics research ideas using PyTheus. — [arXiv](https://arxiv.org/abs/2511.11752)
- Shuxiang Cao et al. (Aspuru-Guzik), k-agents, arXiv:2412.07978; published in *Patterns*, DOI 10.1016/j.patter.2025.101372. LLM agents automate quantum-computing laboratory experiments. — [arXiv](https://arxiv.org/abs/2412.07978)
- Niko R. Reed et al. (Walsworth group), "Principles and Techniques for Agentic Research in Experimental Physics", arXiv:2610.17193. Ten requirements for collaborative agents. — [arXiv](https://arxiv.org/abs/2610.17193)
- Marco Knipfer, Alexander Roman, Konstantin T. Matchev, Katia Matcheva and Sergei Gleyzer, "AI Agents for Variational Quantum Circuit Design", arXiv:2602.19387. — [arXiv](https://arxiv.org/abs/2602.19387)
- S. K. Rithvik, "Anubuddhi", arXiv:2512.15736. — [arXiv](https://arxiv.org/abs/2512.15736)
- Survey: Jiaqi Wei et al., "From AI for Science to Agentic Science", arXiv:2508.14111. — [arXiv](https://arxiv.org/abs/2508.14111)
- Not LLM-based: Shiyu Zhou, Yuxuan Zhang, Sebastian Wetzel, Roger Melko and Xiu-Zhe Luo, "Autonomous phase discovery", arXiv:2609.33802. Differentiable programming plus unsupervised learning search Hamiltonian space adaptively. — [arXiv](https://arxiv.org/abs/2609.33802)

**(d) LLM-guided program and evolutionary search in quantum (FunSearch/AlphaEvolve-style)**
- *AlphaEvolve inside Google Quantum AI.* In "Quantum computation of molecular geometry via many-body nuclear spin echoes", arXiv:2510.19550 (C. Zhang, …, T. E. O'Brien; 22 Oct 2025), AlphaEvolve evolved Python functions that generate circuits approximating Hamiltonian evolution, according to a search snippet. Google Cloud reports circuits with "10x lower error" for molecular simulation on Willow. — [arXiv](https://arxiv.org/abs/2510.19550); [Google Cloud blog](https://cloud.google.com/blog/products/ai-machine-learning/alphaevolve-is-available-for-everyone)
- Aritra Lahiri and Florian Marquardt, "Agentic AI for Discovering Fault-Tolerant Quantum Circuits", arXiv:2611.18504 (13 Nov 2026). Agentic evolutionary program synthesis of Clifford circuits for logical state preparation, logical gates and code switching, guided by a circuit verifier that checks faults exhaustively. All discovered circuits are fault-tolerant, and several match or beat benchmarks. — [arXiv](https://arxiv.org/abs/2611.18504)
- Ge Yan, Shanchuan Li, Pengyue Ma, Qixin Zhang, Pingchuan Ma, Jianping Wang, Min-Hsiu Hsieh and Yuxuan Du, "OmniQEC", arXiv:2607.25865 (28 Jul 2026). An LLM orchestrator coordinates code generation, code-level screening, syndrome-extraction synthesis and decoder-based circuit evaluation. Search snippets add that the discovered codes beat bivariate-bicycle codes at 98 and 240 physical qubits in simulation. — [arXiv](https://arxiv.org/abs/2607.25865)
- Juan Cruz-Benito, Andrew W. Cross, David Kremer and Ismael Faro (IBM), "Evolutionary Discovery of Bivariate Bicycle Codes with LLM-Guided Search", arXiv:2606.02418 (Jun 2026). LLMs mutate Python programs that generate code ansätze, with about 1,650 evolutionary iterations over five campaigns and certified code parameters. — [arXiv](https://arxiv.org/abs/2606.02418)
- Kenny Heitritter, James Brown and Tarini Hardikar, "Evolving Quantum Error-Correcting Encodings for Molecular Simulation", arXiv:2606.25870 (Jun 2026). An LLM edits a program, an external verifier scores it, and high scorers are kept and re-mutated; applied to Generalized Superfast fermion-to-qubit encodings. — [arXiv](https://arxiv.org/abs/2606.25870)
- Daniil Fisher, Valentin Khrulkov, Mikhail Saygin, Ivan Oseledets and Stanislav Straupe, "LLM-Guided Evolutionary Search for Algebraic T-Count Optimization", arXiv:2603.29894. — [arXiv](https://arxiv.org/abs/2603.29894)
- Jae-Hyun Baek, Byung-Sun Won and Jon-Lark Kim, "Researcher-in-the-Loop Evolutionary Discovery of Entanglement-Assisted Quantum Code Families", arXiv:2610.41183. LLM-guided program evolution plus exact evaluation, nonexistence checks and Lean proofs; the result generalizes to the family [[n,1,n−1;n−q−1]]_q. — [arXiv](https://arxiv.org/abs/2610.41183)
- Lily Jiaxin Wan, Deming Chen, Klara Nahrstedt and Bo Chen, "Toward Quantum Software Automation: A Quantum-Aware Harness for LLM-Guided Evolution", arXiv:2609.38883. — [arXiv](https://arxiv.org/abs/2609.38883)
- Tanuj Khattar, Nicholas C. Rubin, Craig Gidney, Michael Broughton, Guang Hao Low and Ryan Babbush, "Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI", arXiv:2610.23832. AI agents compile subroutines whose structure is classically verifiable, certified by fuzz testing. — [arXiv](https://arxiv.org/abs/2610.23832)
- Andres Paz, Christian Tarta, Cordelia Yuqiao Li, Mayee Sun, Sarju Patel and Sylvie Lausier, "StabilizerBench", arXiv:2604.21287. 192 stabilizer codes from 12 families (4–196 qubits, distances 2–21) and three tasks: state preparation, constrained optimization and fault-tolerant synthesis. A search snippet says frontier LLMs were "largely unsuccessful" at fault-tolerant synthesis. — [arXiv](https://arxiv.org/abs/2604.21287)
- Debarshi Kundu, Swaroop Ghosh and Vasant G. Honavar, "Weight-Space Programming of Variational Quantum Circuits", arXiv:2611.11537. A frontier LLM writes circuit angles without feedback over 139 requests (VQE for spin chains and H₂, QAOA, classification). — [arXiv](https://arxiv.org/abs/2611.11537)

**(e) LLM-driven symbolic regression and model discovery (templates from outside quantum)**
- Parshin Shojaee, Kazem Meidani, Shashank Gupta, Amir Barati Farimani and Chandan K. Reddy, "LLM-SR", arXiv:2404.18400, ICLR 2025 (oral). Equations are programs: the LLM proposes skeletons and evolutionary search refines them, with constants fitted to data. — [arXiv](https://arxiv.org/abs/2404.18400)
- Parshin Shojaee, Ngoc-Hieu Nguyen, Kazem Meidani, Amir Barati Farimani, Khoa D. Doan and Chandan K. Reddy, "LLM-SRBench", arXiv:2504.10415, ICML 2025 (oral). 239 problems designed to resist memorization; best symbolic accuracy 31.5% (from the PMLR page via search). — [arXiv](https://arxiv.org/abs/2504.10415); [PMLR](https://proceedings.mlr.press/v267/shojaee25a.html)
- Jianke Yang, Ohm Venkatachalam, Mohammad Kianezhad, Sharvaree Vadgama and Rose Yu, "Think like a Scientist: Physics-guided LLM Agent for Equation Discovery" (KeplerAgent), arXiv:2602.12259. The agent infers symmetries first and uses them as priors to restrict candidate equations. — [arXiv](https://arxiv.org/abs/2602.12259)
- Kevin Murphy, "Model Discovery Agent: LLM-assisted Bayesian experiment design for data-efficient discovery of mechanistic world models", arXiv:2608.09696 (v5 Sep 2026). LLM proposals for open-ended model classes, experiment design by value of information, and approximate Bayesian inference over structures and parameters, including BIC. — [arXiv](https://arxiv.org/abs/2608.09696)
- Further templates:
  - JooKyung Yuh et al., "Bayesian Model Discovery with Language Models", arXiv:2610.25747
  - Kanishk Gandhi et al., "BoxingGym", arXiv:2501.01540
  - Tianshi Zheng et al., "NewtonBench", arXiv:2510.07172 (ICLR 2026)
  - Marco Virgolin et al., "automodel: A Skill for Agentic Model Discovery" for coding agents such as Claude Code and Codex, arXiv:2611.24347
  - Ze-Yu Peng et al., "DeepInflation", arXiv:2601.14288, combining LLMs, a symbolic-regression engine and RAG for inflationary potentials
  - Michael Y. Li, Emily B. Fox and Noah D. Goodman, ICML 2024, arXiv:2402.17879
  - — [2610.25747](https://arxiv.org/abs/2610.25747); [2501.01540](https://arxiv.org/abs/2501.01540); [2510.07172](https://arxiv.org/abs/2510.07172); [2611.24347](https://arxiv.org/abs/2611.24347); [2601.14288](https://arxiv.org/abs/2601.14288); [2402.17879](https://arxiv.org/abs/2402.17879)
- Paulin de Schoulepnikoff, Hendrik Poulsen Nautrup, Hans J. Briegel and Gorka Muñoz-Gil, "Discovering quantum phenomena with Interpretable Machine Learning" (QDisc), arXiv:2604.16015. VAEs plus symbolic methods find order-parameter descriptors from raw quantum data. — [arXiv](https://arxiv.org/abs/2604.16015)

**(f) LLMs deriving Hamiltonians analytically**
- Haining Pan, Nayantara Mudur, Will Taranto, Maria Tikhanovskaya, Subhashini Venugopalan, Yasaman Bahri, Michael P. Brenner and Eun-Ah Kim, "Quantum Many-Body Physics Calculations with Large Language Models", arXiv:2403.03154, *Commun. Phys.* 8, 49 (2025). GPT-4 with multi-step prompt templates derives Hartree–Fock Hamiltonians. A search snippet adds that the final Hamiltonian was correct in 13 of 15 papers *with corrections to intermediate steps*. — [arXiv](https://arxiv.org/abs/2403.03154); [Google Research](https://research.google/pubs/quantum-many-body-physics-calculations-with-large-language-models/)

### Inferences
- **The pattern.** Reliable agentic results in quantum couple LLM proposals with a cheap, hard verifier:
  - exhaustive fault checks (Lahiri & Marquardt)
  - certified code parameters (IBM bivariate-bicycle search)
  - Lean proofs (EAQECC)
  - fuzz testing (Khattar et al.)
  - ED cross-checks (QuantumMatterExplorer)
  - Bayesian evidence (Kleidermacher et al.)
  - an expected-signal calculation (Isogawa et al.)

  Parent-Hamiltonian reconstruction has an exact native verifier, zero variance of H in the target state, plus checkable gap and ground-state conditions. It is therefore a natural next target for this pattern, yet nobody has published it.
- LLM-SR's "skeleton plus fitted constants" maps onto "operator-basis skeleton plus null-space coefficients", and KeplerAgent's "infer symmetries first" maps onto symmetry-adapted operator bases.
- Murphy's Model Discovery Agent (LLM proposals, value-of-information experiment design, BIC/SMC evidence) is a ready-made architecture for the *measurement-data* version of the problem.

### Gaps
- Some LLM model names and figures (QuantumMatterExplorer's "Qwen3.6-27B"; OmniQEC's backends) come from search snippets, not full texts.
- I found no LLM agent applied to fitting spin Hamiltonians to neutron or THz spectroscopy data. Non-LLM machine-learning versions exist (Q4c).

---

## Q3. Map (2): non-LLM, agent-like precursors

### Takeaway
The architectural ancestor is the Bristol/Microsoft QMLA line. Gentile et al. (*Nat. Phys.* 2021) and QMLA (*NJP* 2022) propose candidate models, fit each with quantum Hamiltonian learning, and compare them; the open-system extension (oQMLA) appeared in *PR Applied* in 2025. Bayesian model selection over model *structure* has matured in parallel: trans-dimensional MCMC (Poteshman et al., *Quantum* 2026) and Lindblad-model ranking for quantum emitters (Wallace et al., *PRR* 2025). Adaptive and RL experimental design (Granade 2012 → Hincks 2018 → Fiderer 2021 → SymQNet 2026) supplies the ability to choose the next measurement. Kleidermacher et al. (2026) cite exactly this lineage and reuse its Bayesian model selection inside an LLM agent.

### Cited Findings
- **Antonio A. Gentile, Brian Flynn, Sebastian Knauer, Nathan Wiebe, Stefano Paesani, Christopher E. Granade, John G. Rarity, Raffaele Santagati and Anthony Laing, "Learning models of quantum systems from experiments", arXiv:2002.06169, *Nature Physics* (2021), DOI 10.1038/s41567-021-01201-7.** An agent built on unsupervised machine learning (no LLM) learns models of an NV-centre experiment. With experimental data, 74% of protocol instances retrieve plausible models (from a search snippet of the Nature page). — [arXiv](https://arxiv.org/abs/2002.06169); [Nature Physics](https://www.nature.com/articles/s41567-021-01201-7)
- **Brian Flynn, Antonio Andreas Gentile, Nathan Wiebe, Raffaele Santagati and Anthony Laing, "Quantum Model Learning Agent: characterisation of quantum systems through machine learning", arXiv:2112.08409, *New J. Phys.* (2022), DOI 10.1088/1367-2630/ac68ff.** QMLA reverse-engineers Hamiltonian descriptions and is tested on simulated experiments. — [arXiv](https://arxiv.org/abs/2112.08409)
- **Lorenzo Fioroni, Ivan Rojkov and Florentin Reiter, "Learning agent-based approach to the characterization of open quantum systems" (oQMLA), arXiv:2501.05350, *Phys. Rev. Applied* 24, 034011 (2025), DOI 10.1103/h3mk-g5bv.** QMLA learns a Hamiltonian "via adaptive generations of informative experiments and their validation against simulated models"; oQMLA adds jump operators to capture incoherent noise. — [arXiv](https://arxiv.org/abs/2501.05350)
- **Stewart Wallace, Yoann Altmann, Brian D. Gerardot, Erik M. Gauger and Cristian Bonato, "Learning the dynamics of Markovian open quantum systems from experimental data", arXiv:2410.17942, *Phys. Rev. Research* 7, 033217 (2025).** A Bayesian MCMC algorithm outputs a *ranked list of interpretable master-equation models* matching measured traces, benchmarked on single emitters and emitter pairs. This is the non-LLM analogue of Kleidermacher et al. — [arXiv](https://arxiv.org/abs/2410.17942)
- **Abigail N. Poteshman, Jiwon Yun, Tim H. Taminiau and Giulia Galli, "Trans-dimensional Hamiltonian model selection and parameter estimation from sparse, noisy data", arXiv:2506.18802, *Quantum* 10, 2055 (2026).** Hybrid MCMC estimates both parameters *and model dimension* in ill-posed regimes. — [arXiv](https://arxiv.org/abs/2506.18802)
- **Christopher E. Granade, Christopher Ferrie, Nathan Wiebe and D. G. Cory, "Robust Online Hamiltonian Learning", arXiv:1207.1655, *New J. Phys.* 14, 103013 (2012).** Sequential Monte Carlo combined with Bayesian experimental design. — [arXiv](https://arxiv.org/abs/1207.1655)
- **Ian Hincks, Thomas Alexander, Michal Kononenko, Benjamin Soloway and David G. Cory, "Hamiltonian Learning with Online Bayesian Experiment Design in Practice", arXiv:1806.02427.** Online Bayesian design on an NV centre. — [arXiv](https://arxiv.org/abs/1806.02427)
- **Lukas J. Fiderer, Jonas Schuff and Daniel Braun, "Neural-Network Heuristics for Adaptive Bayesian Quantum Estimation", arXiv:2003.02183, *PRX Quantum* 2, 020303 (2021).** — [arXiv](https://arxiv.org/abs/2003.02183)
- **Vincent Lim, Ellen Novoseller, Jeffrey Ichnowski, Huang Huang and Ken Goldberg, "Policy-Based Bayesian Experimental Design for Non-Differentiable Implicit Models" (RL-DAD), arXiv:2203.04272.** — [arXiv](https://arxiv.org/abs/2203.04272)
- **Yash Vardhan Tomar and Dheeraj Peddireddy, "SymQNet: Amortized Acquisition for Low-Latency Adaptive Hamiltonian Learning", arXiv:2606.12808 (IEEE QCE 2026 workshop).** An RL acquisition policy, conditioned on the posterior, keeps Bayesian feedback; benchmarked on transverse-field Ising models. — [arXiv](https://arxiv.org/abs/2606.12808)
- **Jawaher Kaldari and Saif Al-Kuwari, "Adaptive Hamiltonian Identification via Quantum Reinforcement Learning under Non-Stationary Dynamics", arXiv:2611.26241 (19 Nov 2026).** — [arXiv](https://arxiv.org/abs/2611.26241)
- **Roy Liu et al. (incl. Matthew B. Stone, Zhantao Chen), "Observation geometry for uncertainty-aware Hamiltonian inference and experimental design in quantum magnets", arXiv:2608.10350.** Neural surrogates conditioned on the Hamiltonian, Bayesian inference and "observation geometry" quantify the *identifiability* of microscopic interactions and drive adaptive multimodal design, applied to inelastic neutron scattering on NiPS₃. — [arXiv](https://arxiv.org/abs/2608.10350)
- **Arielle Sanford, Andrew T. Kamen, Frederic T. Chong and Andy J. Goldschmidt, "Data-Driven Hamiltonian Reduction for Superconducting Qubits via Meta-Learning" (HAML), arXiv:2604.24912.** — [arXiv](https://arxiv.org/abs/2604.24912)
- **Algorithmic, non-agentic searches over Hamiltonian space.** Chertkov–Clark EHC and SHC, Pakrouski's "Automatic design of Hamiltonians", and Turkeshi et al.'s entanglement-guided search (all in Q4). — [Pakrouski](https://arxiv.org/abs/1907.05898)

### Inferences
- QMLA's propose → fit with quantum Hamiltonian learning → compare by Bayesian evidence loop, and its structure-selecting descendants (trans-dimensional MCMC, ranked master-equation models), are what an LLM agent should *call as tools*. The LLM's job is generating hypotheses from the literature, which is what Kleidermacher et al. did.
- SciExplorer's under-measurement failure is exactly what information-gain policies address (SymQNet, Fiderer et al., Murphy's Model Discovery Agent, observation geometry). Pairing an LLM proposer with such a policy for parent-Hamiltonian reconstruction is untested.

### Gaps
- QMLA's specific exploration strategies (for example, genetic-algorithm search) were not re-checked in full text this session.

---

## Q4. Map (3): classical and ML baselines, including provable theory

### Takeaway
The baselines an agent would build on and be compared against:
- correlation-matrix (covariance) reconstruction from one eigenstate (Qi & Ranard; Chertkov & Clark; Greiter, Schnells & Thomale), with known sensitivity to the operator basis (Brito et al. 2024)
- local-patch learning (Bairey, Arad & Lindner)
- entanglement-guided search with a Bisognano–Wichmann ansatz, and entanglement-Hamiltonian tomography (Turkeshi et al.; Kokail et al.; Joshi et al.)
- neural-network approaches, now including measures of eigenstate "learnability" and inverse design
- provable learning from Gibbs states (Anshu et al.; Haah, Kothari & Tang; Bakshi et al.) and from dynamics, including structure-learning and ansatz-free Heisenberg-limited protocols (Huang et al.; Bakshi et al.; Zhao; Hu et al.; Zhou & Gong; Gong et al. 2026)

New in 2025–2026: 2D Gutzwiller reconstruction (negative result), conserved-operator extraction from iPEPS, differentiable inverse design for Moore–Read, and conserved-operator classifications. **No provable sample-complexity guarantee was found for learning a Hamiltonian from ground or eigenstate copies.**

### Cited Findings
**(a) Single-eigenstate and correlation-matrix (covariance) methods**
- **Xiao-Liang Qi and Daniel Ranard, "Determining a local Hamiltonian from a single eigenstate", arXiv:1712.01850, *Quantum* 3, 159 (2019).** Two-point equal-time correlators of local observables in one eigenstate generically determine a range-k local Hamiltonian, whether the eigenstate is ground or excited. The "k-correlation spectrum" bounds sensitivity, and numerics cover random Hamiltonians on up to 12 qubits (search snippet). — [arXiv](https://arxiv.org/abs/1712.01850); [Quantum](https://quantum-journal.org/papers/q-2019-07-08-159/pdf/)
- **James R. Garrison and Tarun Grover, "Does a single eigenstate encode the full Hamiltonian?", arXiv:1503.00729, *PRX* 8, 021026 (2018).** Connects the question to which operators satisfy the eigenstate thermalization hypothesis. — [arXiv](https://arxiv.org/abs/1503.00729)
- **Eli Chertkov and Bryan K. Clark, "Computational inverse method for constructing spaces of quantum models from wave functions", arXiv:1802.01590, *PRX* 8, 031029 (2018).** The eigenstate-to-Hamiltonian construction uses the quantum covariance matrix to return the *space* of Hamiltonians with a given eigenstate. — [arXiv](https://arxiv.org/abs/1802.01590)
- **Eli Chertkov, Benjamin Villalonga and Bryan K. Clark, "Engineering Topological Models with a General-Purpose Symmetry-to-Hamiltonian Approach", arXiv:1910.10165, *PRR* 2, 023348 (2020).** The symmetric Hamiltonian construction (SHC) outputs all local Hamiltonians consistent with given symmetries; the open-source code is qosy. — [arXiv](https://arxiv.org/abs/1910.10165)
- **Martin Greiter, Vera Schnells and Ronny Thomale, "Method to identify parent Hamiltonians for trial states", arXiv:1802.07827, *PRB* 98, 081113 (2018).** Finds exact local parents directly or through annihilation operators, and approximate parents otherwise. — [arXiv](https://arxiv.org/abs/1802.07827)
- **Eyal Bairey, Itai Arad and Netanel H. Lindner, "Learning a local Hamiltonian from local measurements", arXiv:1807.04564, *PRL* 122, 020504 (2019).** Local constraints from patches give costs linear in system size.
  - Follow-up: Bairey, Guo, Poletti, Lindner and Arad, Lindbladians from steady states, *NJP* 22, 032001 (2020).
  - Follow-up: Zubida, Yitzhaki, Lindner and Bairey, "Optimal short-time measurements for Hamiltonian learning", arXiv:2108.08824.
  - — [arXiv](https://arxiv.org/abs/1807.04564); [1907.11154](https://arxiv.org/abs/1907.11154); [2108.08824](https://arxiv.org/abs/2108.08824)
- **Shi-Yao Hou, Ningping Cao, Sirui Lu, Yi Shen, Yiu-Tung Poon and Bei Zeng, "Determining system Hamiltonian from eigenstate measurements without correlation functions", arXiv:1903.06569, *NJP* 22, 083088 (2020).** Companion: Chenfeng Cao, Shi-Yao Hou, Ningping Cao and Bei Zeng, "Supervised learning in Hamiltonian reconstruction from local measurements on eigenstates", arXiv:2007.05962, *J. Phys.: Condens. Matter* 33, 064002. — [1903.06569](https://arxiv.org/abs/1903.06569); [2007.05962](https://arxiv.org/abs/2007.05962)
- **Maxime Dupont, Nicolas Macé and Nicolas Laflorencie, "From eigenstate to Hamiltonian: Prospects for ergodicity and localization", arXiv:1907.12124, *PRB* 100, 134201 (2019).** Searches for possibly *multiple* parent Hamiltonians of one state. — [arXiv](https://arxiv.org/abs/1907.12124)
- **Xhek Turkeshi and Marcello Dalmonte, "Parent Hamiltonian Reconstruction of Jastrow-Gutzwiller Wavefunctions", arXiv:1909.11327, *SciPost Phys.* 8, 042 (2020).** Approximate spin-chain parents. Exact local parents appear only in the free-fermion case, and the method certifies the ground state, not only an eigenstate (search snippet). — [arXiv](https://arxiv.org/abs/1909.11327); [SciPost](https://scipost.org/SciPostPhys.8.3.042)
- **Kevin Zhang, Samuel Lederer, Kenny Choo, Titus Neupert, Giuseppe Carleo and Eun-Ah Kim, "Hamiltonian reconstruction as metric for variational studies", arXiv:2102.00019, *SciPost Phys.* 13, 063 (2022).** Follow-up: Leo Joon Il Moon, …, Eun-Ah Kim and Peter L. McMahon, "Hamiltonian-reconstruction distance as a success metric for the Variational Quantum Eigensolver", arXiv:2403.11995, *Quantum Sci. Technol.* (2026). — [2102.00019](https://arxiv.org/abs/2102.00019); [2403.11995](https://arxiv.org/abs/2403.11995)
- **J. Alexander Jacoby and J. B. Marston, "Reconstructing Hamiltonians from Correlations and Entanglement", arXiv:2112.05107.** Compares reconstruction from correlations with reconstruction from entanglement spectra. — [arXiv](https://arxiv.org/abs/2112.05107)
- **Lucas Z. Brito, Stephen Carr, J. Alexander Jacoby and J. B. Marston, "Hamiltonian Reconstruction: the Correlation Matrix and Incomplete Operator Bases", arXiv:2311.09302, *PRB* 110, 125118 (2024).** Studies robustness to *undercomplete and overcomplete operator bases* and proposes an approximation scheme for undercomplete ones. — [arXiv](https://arxiv.org/abs/2311.09302)
- **Lucas Z. Brito and J. B. Marston, "Reconstructing Spin Hamiltonians of 2D Gutzwiller-Projected Wavefunctions", arXiv:2510.15034 (16 Oct 2025).** Correlation-matrix reconstruction on finite square and triangular lattices finds that "no spin Hamiltonian with simple local interaction terms stabilizes such states". — [arXiv](https://arxiv.org/abs/2510.15034)
- **Parent-Hamiltonian optimization by Rattacaso and coauthors.**
  - Davide Rattacaso, Gianluca Passarelli, Angelo Russomanno, Procolo Lucignano, Giuseppe E. Santoro and Rosario Fazio, "Parent Hamiltonian Reconstruction via Inverse Quantum Annealing", arXiv:2303.11200, *PRL* 132, 160401 (2024).
  - Rattacaso, Passarelli, Mezzacapo, Lucignano and Fazio, "Optimal parent Hamiltonians for time-dependent states", arXiv:2105.10187, *PRA* 104, 022611 (2021).
  - — [2303.11200](https://arxiv.org/abs/2303.11200); [2105.10187](https://arxiv.org/abs/2105.10187)
- **Non-Hermitian parent Hamiltonians.**
  - Yin Tang and W. Zhu, arXiv:2307.03107, *PRB* 108, 165114 (2023): a generalized covariance matrix for a biorthogonal pair of eigenstates.
  - Xu-Dan Xie, Zheng-Yuan Xue and Dan-Bo Zhang, arXiv:2408.15781, *PRB* 110, 235113 (2024): reconstruction from a single right or left eigenstate.
  - Ruohan Shen, Yuchen Guo and Shuo Yang, arXiv:2301.12448, *PRL* 130, 220401 (2023): non-Hermitian parents from MPS.
  - Experimental: Yuntao Xiao et al. (Peng Xue), arXiv:2607.28964, *PRL* 137, 110401 (2026), photonics.
  - — [2307.03107](https://arxiv.org/abs/2307.03107); [2408.15781](https://arxiv.org/abs/2408.15781); [2301.12448](https://arxiv.org/abs/2301.12448); [2607.28964](https://arxiv.org/abs/2607.28964)
- **Conserved-operator generalizations (2023–2026).**
  - Federico Petrovich, R. Rossignoli and M. Cerezo, "Classification of quantum states through conserved operators", arXiv:2610.20988. Conserved operators are obtained "from the kernel of a quantum covariance matrix" and generate parent Hamiltonians with the state as exact ground state.
  - J. Pawłowski, P. Łydżba and M. Mierzejewski, "Local integrals of motion encoded in a few eigenstates", arXiv:2603.01859.
  - Moudgalya and Motrunich, "Numerical Methods for Detecting Symmetries and Commutant Algebras", arXiv:2302.03028, *PRB* 107, 224312.
  - — [2610.20988](https://arxiv.org/abs/2610.20988); [2603.01859](https://arxiv.org/abs/2603.01859); [2302.03028](https://arxiv.org/abs/2302.03028)
- **Kiryl Pakrouski, "Automatic design of Hamiltonians", arXiv:1907.05898, *Quantum* 4, 315 (2020).** A variational Hamiltonian ansatz with a loss given by a weighted sum of desired properties, tested on the fractional quantum Hall effect. — [arXiv](https://arxiv.org/abs/1907.05898)

**(b) Entanglement-guided search and entanglement-Hamiltonian tomography**
- **X. Turkeshi, T. Mendes-Santos, G. Giudici and M. Dalmonte, "Entanglement guided search for parent Hamiltonians", arXiv:1807.06113, *PRL* 122, 150606 (2019).** An ansatz relating the reduced density matrix to the parent Hamiltonian (exact at field-theory level, Bisognano–Wichmann) plus minimization over relative entropies. — [arXiv](https://arxiv.org/abs/1807.06113)
- **Marcello Dalmonte, Benoît Vermersch and Peter Zoller, "Quantum Simulation and Spectroscopy of Entanglement Hamiltonians", arXiv:1707.04455, *Nature Physics* (2018), DOI 10.1038/s41567-018-0151-7.** — [arXiv](https://arxiv.org/abs/1707.04455)
- **Christian Kokail, Rick van Bijnen, Andreas Elben, Benoît Vermersch and Peter Zoller, "Entanglement Hamiltonian Tomography in Quantum Simulation", arXiv:2009.09000, *Nature Physics* 17, 936–942 (2021).** Applied to data from 10- and 20-ion simulators. — [arXiv](https://arxiv.org/abs/2009.09000)
- **Kokail, Sundar, Zache, Elben, Vermersch, Dalmonte, van Bijnen and Zoller, "Quantum Variational Learning of the Entanglement Hamiltonian", arXiv:2105.04317, *PRL* 127, 170501 (2021).** Spatial deformations of the many-body Hamiltonian, realized on the device, serve as the ansatz. — [arXiv](https://arxiv.org/abs/2105.04317)
- **Joshi, Kokail, van Bijnen, Kranzl, Zache, Blatt, Roos and Zoller, "Exploring Large-Scale Entanglement in Quantum Simulation", arXiv:2306.00057, *Nature* 624, 539–544 (2023).** On 51 ions, the entanglement Hamiltonian was learned for subsystems of up to 20 sites (search snippet). — [arXiv](https://arxiv.org/abs/2306.00057)
- **Siyi Yang, Yi-Ming Ding and Zheng Yan, "Exploring the limit of the Lattice-Bisognano-Wichmann form describing the Entanglement Hamiltonian: A quantum Monte Carlo study", arXiv:2511.00950, *PRB* (2026).** Tests where the Bisognano–Wichmann ansatz fails. — [arXiv](https://arxiv.org/abs/2511.00950)

**(c) Neural-network and other machine-learning approaches**
- Agnes Valenti, Evert van Nieuwenburg, Sebastian Huber and Eliska Greplova, "Hamiltonian Learning for Quantum Error Correction", arXiv:1907.02540, *PRR* 1, 033092 (2019). — [arXiv](https://arxiv.org/abs/1907.02540)
- Maksymilian Kliczkowski, Jarosław Pawłowski and Masudul Haque, "Information in Many-body Eigenstates: A Question of Learnability", arXiv:2605.03043 (May 2026). An encoder-decoder network reconstructs the Hamiltonian of a non-integrable chain. **Eigenstates at the spectral edges allow more accurate reconstruction from fewer eigenstates than mid-spectrum ones.** — [arXiv](https://arxiv.org/abs/2605.03043)
- Jie Liu and Xin Wang, iPINN-HL, arXiv:2506.10379, *PRR* 7, 043137 (2025). Kris Tucker et al., arXiv:2404.05526, *PR Applied* 22, 044080 (2024): learning from continuous weak measurements. — [2506.10379](https://arxiv.org/abs/2506.10379); [2404.05526](https://arxiv.org/abs/2404.05526)
- Ang-Kun Wu, David A. Dahlbom, Kipton Barros, Cristian D. Batista and Hao Zhang, "Learning Spin Hamiltonians from Interacting Magnon Spectra", arXiv:2611.25321 (Nov 2026). A network that learns an inverse map "inherits the approximations of the forward solver used to generate its training data", shown against linear spin-wave theory. — [arXiv](https://arxiv.org/abs/2611.25321)
- Martin Mootz, Chuankun Huang, Liang Luo, Jigang Wang and Yong-Xin Yao, "Learning Spin Hamiltonians from Terahertz Two-Dimensional Coherent Spectroscopy", arXiv:2608.14460. — [arXiv](https://arxiv.org/abs/2608.14460)
- Inverse design with ML or autodiff:
  - Alessandro Sinibaldi, Luciano Loris Viteritti, Riccardo Rende, Federico Becca and Giuseppe Carleo, "Designing Strongly Correlated Quantum Phases of Matter with Foundation Neural-Network Quantum States", arXiv:2609.40209
  - Yi Yang, Xin Wan and Zi-Xiang Hu, "Differentiable Inverse Design of Two-Body Interactions for the Moore-Read State", arXiv:2607.23695: autodiff plus ED, overlaps above 99% for up to 12 electrons and 91.1% at 14
  - Koji Inui and Yukitoshi Motome, arXiv:2402.15802, *PRR* 6, 033080 (2024)
  - — [2609.40209](https://arxiv.org/abs/2609.40209); [2607.23695](https://arxiv.org/abs/2607.23695); [2402.15802](https://arxiv.org/abs/2402.15802)

**(d) Provable Hamiltonian learning, including structure learning**
- **Anurag Anshu, Srinivasan Arunachalam, Tomotaka Kuwahara and Mehdi Soleimanifar, arXiv:2004.07266, *Nature Physics* (2021), DOI 10.1038/s41567-021-01232-0 (FOCS 2020).** First sample-efficient learning from Gibbs states. — [arXiv](https://arxiv.org/abs/2004.07266)
- **Jeongwan Haah, Robin Kothari and Ewin Tang, "Optimal learning of quantum Hamiltonians from high-temperature Gibbs states", arXiv:2108.04842, *Nat. Phys.* 20, 1027–1031 (2024).** — [arXiv](https://arxiv.org/abs/2108.04842)
- **Ainesh Bakshi, Allen Liu, Ankur Moitra and Ewin Tang, "Learning quantum Hamiltonians at any temperature in polynomial time", arXiv:2310.02243 (v2 May 2026).** — [arXiv](https://arxiv.org/abs/2310.02243)
- **Hsin-Yuan Huang, Yu Tong, Di Fang and Yuan Su, "Learning many-body Hamiltonians with Heisenberg-limited scaling", arXiv:2210.03030, *PRL* 130, 200403 (2023).** — [arXiv](https://arxiv.org/abs/2210.03030)
- **Bakshi, Liu, Moitra and Tang, "Structure learning of Hamiltonians from real-time evolution", arXiv:2405.00082, FOCS 2024 (v4 May 2026 strengthens the main theorem).** — [arXiv](https://arxiv.org/abs/2405.00082)
- **Andrew Zhao, "Learning the structure of any Hamiltonian from minimal assumptions", arXiv:2410.21635, STOC 2025.** — [arXiv](https://arxiv.org/abs/2410.21635)
- **Hong-Ye Hu, Muzhou Ma, Weiyuan Gong, Qi Ye, Yu Tong, Steven T. Flammia and Susanne F. Yelin, "Ansatz-free Hamiltonian learning with Heisenberg-limited scaling", arXiv:2502.11900, *PRX Quantum* 6, 040315 (2025).** — [arXiv](https://arxiv.org/abs/2502.11900)
- **2026 theory papers.**
  - Taiqi Zhou and Weiyuan Gong, "Optimal Ansatz-free Hamiltonian Learning In Situ", arXiv:2606.19486
  - Weiyuan Gong, Muzhou Ma, Sitan Chen, Jordan Cotler and Hsin-Yuan Huang, "Out-of-control Hamiltonian Learning", arXiv:2610.06709, which removes assumptions of fast arbitrary control
  - Constantin Cedillo Vayson de Pradenne, Jordan Cotler and Hsin-Yuan Huang, "Learning Hamiltonians at Long Times", arXiv:2606.05690
  - Francisco Escudero Gutiérrez, Junseo Lee and Sebastian Zur, "Hamiltonian locality testing and certification do not achieve the Heisenberg limit", arXiv:2610.03205
  - Bora Baran and Timothy Heightman, arXiv:2507.21374
  - Haoya Li, Yu Tong, Hongkang Ni, Tuvia Gefen and Lexing Ying (bosons), arXiv:2307.04690
  - — [2606.19486](https://arxiv.org/abs/2606.19486); [2610.06709](https://arxiv.org/abs/2610.06709); [2606.05690](https://arxiv.org/abs/2606.05690); [2610.03205](https://arxiv.org/abs/2610.03205)
- **Hamiltonian learning at scale and in experiments.**
  - Nathan Johnson and Eugene Dumitrescu, "Hamiltonian Learning at Scale", arXiv:2608.23707: tensor-network emulation of noisy learning workflows up to 300 qubits
  - S.-A. Guo et al. (L.-M. Duan), "Hamiltonian learning for 300 trapped ion qubits with long-range couplings", arXiv:2408.03801, *Sci. Adv.* 11, eadt4713 (2025)
  - — [2608.23707](https://arxiv.org/abs/2608.23707); [2408.03801](https://arxiv.org/abs/2408.03801)
- **Ground-state theory gap.** Phase-1 WebSearch and phase-2 arXiv queries found no provable sample-complexity result for learning a local Hamiltonian from *ground or eigenstate* copies. Related limits: Thomas Schuster, Dominik Kufel, Norman Y. Yao and Hsin-Yuan Huang, "Hardness of recognizing phases of matter", arXiv:2510.08503. — [arXiv](https://arxiv.org/abs/2510.08503)

**(e) Verification status of the example works named in the assignment** (all checked through arXiv API metadata in phase 2)
- Verified:
  - Qi & Ranard (*Quantum* 3, 159)
  - Chertkov & Clark (*PRX* 8, 031029)
  - Greiter, Schnells & Thomale (*PRB* 98, 081113)
  - Bairey, Arad & Lindner (*PRL* 122, 020504)
  - Turkeshi, Mendes-Santos, Giudici & Dalmonte (*PRL* 122, 150606)
  - Kokail et al. (*Nat. Phys.* 17, 936)
  - Dalmonte, Vermersch & Zoller (*Nat. Phys.* 2018)
  - Anshu et al. (*Nat. Phys.* 2021)
  - Haah, Kothari & Tang (*Nat. Phys.* 20, 1027 (2024), not only FOCS)
  - Huang, Tong, Fang & Su (*PRL* 130, 200403)
  - Bakshi, Liu, Moitra & Tang (structure learning, FOCS 2024; any temperature, arXiv:2310.02243)
  - Gentile et al. (*Nat. Phys.* 2021)
  - QMLA (Flynn, Gentile, Wiebe, Santagati, Laing; *NJP* 2022)
  - — [arXiv API listing](https://export.arxiv.org/api/query?id_list=1712.01850,1802.01590,1802.07827,1807.04564,1807.06113,2009.09000,1707.04455,2004.07266,2108.04842,2210.03030,2405.00082,2310.02243,2002.06169,2112.08409&max_results=20)
- Also verified: Garrison & Grover (*PRX* 8, 021026) and Chertkov, Villalonga & Clark (*PRR* 2, 023348). — [1503.00729](https://arxiv.org/abs/1503.00729); [1910.10165](https://arxiv.org/abs/1910.10165)

### Inferences
- **A verifier stack built from these baselines.** Write H = Σ c_a O_a. The variance of H in |ψ⟩ is cᵀMc, where M is the covariance matrix of the basis operators. The null space of M is therefore the eigenstate certificate (Qi–Ranard; Chertkov–Clark; Petrovich et al.). The smallest nonzero eigenvalue measures robustness. Completeness of the basis matters (Brito et al. 2024). Ground-state status and the gap need separate checks: DMRG/ED, injectivity or martingale bounds for tensor networks (Q5), or the ground-state-certifying method of Turkeshi–Dalmonte. None of these steps needs an LLM to be trusted, which is what makes them good verifiers.
- **Theory coverage.** Gibbs-state and dynamics access have strong guarantees, including structure learning without an ansatz. The ground-state and eigenstate setting, the true "parent Hamiltonian" setting, has uniqueness evidence, sensitivity bounds and empirical learnability results (Kliczkowski et al.) but no comparable provable theory. Agent-driven empirical exploration is most needed exactly there.

### Gaps
- The full texts of the 2026 theory papers were not read. Their exact bounds are not reproduced here.

---

## Q5. Map (4): parent Hamiltonians in QEC and circuit contexts

### Takeaway
For tensor-network, stabilizer and shallow-circuit states, parent Hamiltonians are constructive:
- injective PEPS are unique ground states of local parents
- MPS parents are always gapped (simplified proof in 2026)
- stabilizer states are fixed by their stabilizer group, learnable by Bell sampling
- shallow-circuit states are learnable through local inversions (STOC 2024/2025)

Activity in 2025–2026:
- generalized (eigenstate) parents of MPS, stabilizer states and W states as quantum many-body scars
- gapped and gapless parents of the *deformed toric code*
- ZX-calculus inverse design
- parents for bosonic cat-code resources
- certifying states on hardware with parent Hamiltonians
- learning Hamiltonians and Lindbladians *for* QEC hardware, with QEC as a learning primitive

Agentic LLM search in QEC is booming (codes, fault-tolerant circuits, encodings) but has not been applied to parent Hamiltonians of code states.

### Cited Findings
**(a) Tensor-network parent Hamiltonians and gap certificates**
- David Pérez-García, Frank Verstraete, J. Ignacio Cirac and Michael M. Wolf, "PEPS as unique ground states of local Hamiltonians", arXiv:0707.2260, *Quant. Inf. Comp.* 8, 650–663 (2008). Injectivity implies a unique ground state of the local parent (AKLT is an example). — [arXiv](https://arxiv.org/abs/0707.2260)
- Ignacio Cirac, David Pérez-García, Norbert Schuch and Frank Verstraete, review, arXiv:2011.12127, *Rev. Mod. Phys.* 93, 045003 (2021). — [arXiv](https://arxiv.org/abs/2011.12127)
- Carlos Fernández-González, Norbert Schuch, Michael M. Wolf, J. Ignacio Cirac and David Pérez-García, arXiv:1210.6613, *Commun. Math. Phys.* 333, 299 (2015). Gapless frustration-free "uncle" Hamiltonians for non-injective MPS. — [arXiv](https://arxiv.org/abs/1210.6613)
- Milán Ádám Rozmán, András Molnár and Norbert Schuch, "Lower Bounds on Spectral Gaps of Parent Hamiltonians via Tensor Networks", arXiv:2607.19078 (Jul 2026). An improved martingale method computes the overlap of local ground spaces "exactly and efficiently", outperforms other gap lower bounds, and gives a simplified proof that parents of any well-behaved MPS are gapped. Related: Michael J. Kastoryano, hexagonal AKLT gap, arXiv:2609.27993 (title only). — [2607.19078](https://arxiv.org/abs/2607.19078); [2609.27993](https://arxiv.org/abs/2609.27993)
- Norbert Schuch, Andras Molnar and David Perez-Garcia, "Simple Hamiltonians for Matrix Product State models", arXiv:2503.10767 (2025). — [arXiv](https://arxiv.org/abs/2503.10767)
- He-Ran Wang, Lorenzo Gotta and Sanjay Moudgalya, "Exhaustive Construction of Generalized Parent Hamiltonians for Matrix Product States", arXiv:2611.00203 (30 Oct 2026). — [arXiv](https://arxiv.org/abs/2611.00203)
- Lei Gioia, Sanjay Moudgalya and Olexei I. Motrunich, "Distinct Types of Parent Hamiltonians for Quantum States: Insights from the W State as a Quantum Many-Body Scar", arXiv:2510.24713. — [arXiv](https://arxiv.org/abs/2510.24713)
- Wen-Tao Xu, Miguel Frías Pérez and Mingru Yang, "Extracting conserved operators from a projected entangled pair state", arXiv:2511.20619 (Nov 2025; v2 Apr 2026).
  - *Method.* Extracts geometrically k-local operators, including Hamiltonians, that have a given infinite PEPS as an (approximate) eigenstate. Static structure factors are computed by differentiating generating functions, and conserved operators are the directions of vanishing fidelity susceptibility.
  - *Results.* Finds parents "beyond the standard construction" with better locality, including a 4-site plaquette Hamiltonian that approximately has the short-range RVB state as ground state, and a Hamiltonian for which the deformed toric code is an excited eigenstate (a scar). — [arXiv](https://arxiv.org/abs/2511.20619)
- Yuhan Liu, Alberto Ruiz-de-Alarcón, Georgios Styliaris, Xiao-Qi Sun, David Pérez-García and J. Ignacio Cirac, "Parent Lindbladians for Matrix Product Density Operators", arXiv:2501.10552, *PRR* 8, 013210 (2026). — [arXiv](https://arxiv.org/abs/2501.10552)

**(b) Code states and stabilizer states**
- Rahul Sahay, Curt von Keyserlingk, Ruben Verresen and Carolyn Zhang, "Enforced Gaplessness from States with Exponentially Decaying Correlations", arXiv:2503.01977 (2025). For the deformed toric code, any local parent must be gapless or have a degeneracy growing with system size. For local frustration-free Hamiltonians the finite-size gap is at most O(1/L³). — [arXiv](https://arxiv.org/abs/2503.01977)
- Nandagopal Manoj, Zack Weinstein and Jason Alicea, "Gapped Parent Hamiltonians for the Strongly Deformed Toric Code", arXiv:2607.28740 (Jul 2026). Rigorous *gapped* parents exist if one allows Wilson-loop sums whose coefficients decay exponentially, evading the Sahay et al. no-go by relaxing locality. — [arXiv](https://arxiv.org/abs/2607.28740)
- Shane Dooley, "Parent Hamiltonians for stabilizer quantum many-body scars", arXiv:2601.10805 (Jan 2026). Factorizing Pauli strings turns stabilizer elements into local few-body annihilators, so cluster, toric-code and "antipodal toric code" states become zero-energy mid-spectrum scars. ED confirms the construction. — [arXiv](https://arxiv.org/abs/2601.10805)
- Šimon Bräuer and Klaus Mølmer, "Symmetry-resolved parent Hamiltonians for entangled bosonic cat resources", arXiv:2607.02997 (2026). Parents for GHZ-, cluster- and W-type multimode cat states reduce to qubit parents based on correlations, stabilizers and exchange symmetry, and provide a common basis for preparation, dissipative stabilization and certification. — [arXiv](https://arxiv.org/abs/2607.02997)
- Ashley Montanaro, "Learning stabilizer states by Bell sampling", arXiv:1707.04012. O(n) copies suffice.
  - Sabee Grewal, Vishnu Iyer, William Kretschmer and Daniel Liang, arXiv:2304.13915 (STOC 2024).
  - Rebecca Chang, Matthias C. Caro, Martin Larocca and Maxwell West, "Sample-optimal learning of stabilizer states", arXiv:2609.10974.
  - Guglielmo Lami and Mario Collura, "Learning the stabilizer group of a Matrix Product State", arXiv:2401.16481.
  - — [1707.04012](https://arxiv.org/abs/1707.04012); [2304.13915](https://arxiv.org/abs/2304.13915); [2609.10974](https://arxiv.org/abs/2609.10974); [2401.16481](https://arxiv.org/abs/2401.16481)

**(c) Circuit-generated states**
- Fumiyoshi Kobayashi, Kosuke Mitarai and Keisuke Fujii, "Parent Hamiltonian as a benchmark problem for variational quantum eigensolvers", arXiv:2109.11759, *PRA* 105, 052415 (2022). — [arXiv](https://arxiv.org/abs/2109.11759)
- María García Díaz, Gael Sentís, Ramon Muñoz-Tapia and Anna Sanpera, "Volumes of parent Hamiltonians for benchmarking quantum simulators", arXiv:2111.13734. — [arXiv](https://arxiv.org/abs/2111.13734)
- Hsin-Yuan Huang, Yunchao Liu, Michael Broughton, Isaac Kim, Anurag Anshu, Zeph Landau and Jarrod R. McClean, "Learning shallow quantum circuits", arXiv:2401.10095 (STOC 2024).
  - Zeph Landau and Yunchao Liu, "Learning quantum states prepared by shallow circuits in polynomial time", arXiv:2410.23618 (STOC 2025): learns a constant-depth circuit preparing |ψ⟩ = U|0ⁿ⟩.
  - Ning Bao and Gun Suer, "Bootstrapping shallow circuits", arXiv:2403.14746.
  - — [2401.10095](https://arxiv.org/abs/2401.10095); [2410.23618](https://arxiv.org/abs/2410.23618); [2403.14746](https://arxiv.org/abs/2403.14746)
- Wucheng Zhang, "Sequential Generation of Two-dimensional Super-area-law States with Local Parent Hamiltonian", arXiv:2505.02914, *PRX Quantum* 7, 010311 (2026). — [arXiv](https://arxiv.org/abs/2505.02914)
- Marcin Szyniszewski, "Quantum many-body scars with tunable entanglement and Hamiltonian inverse design from ZX-calculus", arXiv:2607.22495. Parent-Hamiltonian methods combined with "local ZX identities that certify exact annihilation of the target state" produce frustration-free Hamiltonians written in the same diagrammatic language. — [arXiv](https://arxiv.org/abs/2607.22495)
- Guy-Philippe Nadon et al. (Jacob Biamonte), "Quantum State Certification via Effective Parent Hamiltonians from Local Measurement Data", arXiv:2603.04499 (Mar 2026). Fidelity is lower-bounded from local measurements of parent-Hamiltonian terms. On IBM hardware this certifies genuine multipartite entanglement of W_n for n ≤ 6 and positive fidelity bounds up to 13 qubits. — [arXiv](https://arxiv.org/abs/2603.04499)
- Farias et al. 2026 use parents of parameterized circuits to benchmark LLMs (Q1). — [arXiv](https://arxiv.org/abs/2611.17246)

**(d) Hamiltonian and Lindbladian learning for QEC hardware**
- Valenti et al., *PRR* 2019, with code (Q4c). — [arXiv](https://arxiv.org/abs/1907.02540)
- Dominik Hangleiter, Ingo Roth, Jonas Fuksa, Jens Eisert and Pedram Roushan, "Robustly learning the Hamiltonian dynamics of a superconducting quantum processor", arXiv:2108.08319, *Nat. Commun.* 15, 9595 (2024). — [arXiv](https://arxiv.org/abs/2108.08319)
- Nikita Romanov, Petr Ivashkov, Weiyuan Gong, Ishaan Kannan, Andi Gu, Hong-Ye Hu and Susanne F. Yelin, "Learning Arbitrary Lindbladians with Quantum Error Correction", arXiv:2606.18188. Ansatz-free Lindbladian learning: Hamiltonian components not obscured by dissipation are learnable at the Heisenberg limit. — [arXiv](https://arxiv.org/abs/2606.18188)
- Petr Ivashkov, Nikita Romanov, Weiyuan Gong, Andi Gu, Hong-Ye Hu and Susanne F. Yelin, "Ansatz-Free Learning of Lindbladian Dynamics In Situ", arXiv:2603.05492. — [arXiv](https://arxiv.org/abs/2603.05492)
- Marina Gordeychuk and Oleg Kiselev, "Code-Visible Liouvillian Modes and Coherent Logical Errors in Quantum Codes", arXiv:2609.18978. Coherent Hamiltonian drifts produce logical phase components that Pauli twirling removes. — [arXiv](https://arxiv.org/abs/2609.18978)
- Ewout van den Berg, Brad Mitchell, Ken Xuan Wei and Moein Malekakhlagh (IBM), "Large-scale Lindblad learning from time-series data", arXiv:2512.08165. — [arXiv](https://arxiv.org/abs/2512.08165)
- Further hardware Lindblad learning:
  - Johann Bock Severin et al. (Morten Kjaergaard), arXiv:2605.00626
  - Timothy Heightman et al., "Lindbladian Learning with Neural Differential Equations", arXiv:2603.07778
  - Tristan Kraft et al. (Peter Zoller), "Bounded-Error Quantum Simulation via Hamiltonian and Lindbladian Learning", arXiv:2511.23392, *PRX* 16, 031037 (2026)
  - — [2605.00626](https://arxiv.org/abs/2605.00626); [2603.07778](https://arxiv.org/abs/2603.07778); [2511.23392](https://arxiv.org/abs/2511.23392)

**(e) Agentic search in QEC and circuits (verifier lessons)**
- Lahiri & Marquardt (exhaustive fault checks), OmniQEC (fast proxy loop and slow circuit-level loop), IBM bivariate-bicycle LLM evolution (certified parameters), Baek et al. (Lean proofs), Heitritter et al. (external verifier for encodings), StabilizerBench (benchmark). See Q2d. — [2611.18504](https://arxiv.org/abs/2611.18504); [2607.25865](https://arxiv.org/abs/2607.25865); [2606.02418](https://arxiv.org/abs/2606.02418); [2604.21287](https://arxiv.org/abs/2604.21287)

### Inferences
- **Textbook constructions (standard facts, not from the sources above).**
  - For a stabilizer or code state, H = −Σ g_i is a commuting, frustration-free parent with gap 2.
  - For |ψ⟩ = U|0ⁿ⟩ with U of depth d, H = −Σ_i U Z_i U† is a commuting parent with unique ground state and gap 2, whose terms are supported on light cones.

  Reconstructing a parent for these states therefore reduces to learning the circuit or stabilizer group, which has provable algorithms. **An agent adds value through *constrained* parents:**
  - few-body, native or 2-local terms (Dooley's factorization of Pauli strings is a hand-made instance)
  - quasi-local versus strictly local trade-offs (deformed toric code: Sahay et al. versus Manoj et al.)
  - certified gaps (Rozmán et al. martingale)
  - symbolic certificates (Szyniszewski's ZX identities)
  - scar-type generalized parents (Wang–Gotta–Moudgalya; Gioia et al.; Xu et al.)
- Certifying states with parent Hamiltonians on hardware (Nadon et al.) means a reconstructed parent is directly useful as a *witness* for code-state preparation. That could serve as an agent's objective in a QEC setting.
- Hamiltonian learning for QEC hardware is moving toward syndromes as data (Romanov et al.) and coherent errors visible to the code (Gordeychuk & Kiselev). An LLM agent running QMLA-style model selection over coherent-error models of a code patch would join maps (2) and (4); none was found.

### Gaps
- No paper was found that reconstructs non-trivial (non-stabilizer-sum) parent Hamiltonians of code states with search or agents, or that applies an LLM agent to Hamiltonian learning for QEC hardware.

---

## Q6. If the exact thing has not been done, what is the closest existing work?

### Takeaway
Closest first:
1. SciExplorer: an LLM agent learning spin Hamiltonians on 9–10 spins from ground-state expectation values it chooses.
2. Kleidermacher et al.: an LLM agent selecting Hamiltonian and Lindbladian models of experimental cavity-QED data by Bayesian methods.
3. QMLA, oQMLA and Gentile et al., plus the Bayesian structure selectors (Wallace et al.; Poteshman et al.): propose-and-compare without an LLM.
4. Algorithmic inverse searches: Chertkov–Clark EHC/SHC, Turkeshi et al., Pakrouski, Xu et al. (iPEPS) and Yang et al. (differentiable).
5. LLM evolutionary search with exact verifiers, aimed at QEC circuits and codes (Lahiri–Marquardt; OmniQEC; IBM).

### Cited Findings
- (1) SciExplorer. GPT-5 agent; ground-state systems are the 1D transverse-field Ising (10 spins), 2D Heisenberg (9), cluster-Ising (10) and an arbitrary translationally invariant model (10), with N up to about 12. Benchmark verifier: fidelity per spin. Internal check: ED forward simulation. Stated failure: non-identifiability when too few expectation values are measured. — [full text](https://arxiv.org/html/2509.24978v6)
- (2) Kleidermacher et al. Claude Sonnet 4.6 agent in Claude Code; Tavis–Cummings/Lindblad models with up to N=6 emitters; BIC, AICc, nested-sampling log Z and Bayes factors; reproducibility quantified. — [full text](https://arxiv.org/html/2611.00620v1)
- (3) QMLA, oQMLA, Gentile et al., Wallace et al., Poteshman et al. — [2112.08409](https://arxiv.org/abs/2112.08409); [2501.05350](https://arxiv.org/abs/2501.05350); [2410.17942](https://arxiv.org/abs/2410.17942); [2506.18802](https://arxiv.org/abs/2506.18802)
- (4) Chertkov–Clark, Chertkov–Villalonga–Clark, Turkeshi et al., Pakrouski, Xu et al., Yang et al. — [1802.01590](https://arxiv.org/abs/1802.01590); [1910.10165](https://arxiv.org/abs/1910.10165); [1807.06113](https://arxiv.org/abs/1807.06113); [2511.20619](https://arxiv.org/abs/2511.20619); [2607.23695](https://arxiv.org/abs/2607.23695)
- (5) Lahiri–Marquardt; OmniQEC; IBM bivariate-bicycle evolution. — [2611.18504](https://arxiv.org/abs/2611.18504); [2607.25865](https://arxiv.org/abs/2607.25865); [2606.02418](https://arxiv.org/abs/2606.02418)
- Closest in the reverse direction: parent Hamiltonians used as LLM benchmarks (Farias et al.) and as VQE benchmarks (Kobayashi et al.). — [2611.17246](https://arxiv.org/abs/2611.17246); [2109.11759](https://arxiv.org/abs/2109.11759)

### Inferences
- A minimal "agentic parent Hamiltonian" paper could combine:
  - SciExplorer's harness, which is public
  - the correlation-matrix null-space solver as an exact inner tool, with Brito et al.'s handling of incomplete bases
  - SHC/qosy-style symmetry-adapted operator bases
  - DMRG/ED validation gates in the style of QuantumMatterExplorer
  - the martingale gap bounds of Rozmán et al. for MPS-like states

  Each piece exists; the combination has not been published, according to these searches.

### Gaps
- SciExplorer's prompts give the agent generic coding tools only. Whether an agent would *discover* the covariance-matrix method unprompted is untested.

---

## Gaps and opportunities (known hard parts where an agent could add value)

### Takeaway
The hard parts are well documented:
- non-uniqueness and identifiability, including multiple viable solutions and dependence on where in the spectrum the state sits
- choice and completeness of the locality ansatz and operator basis (exact local parents often do not exist, and some states force gaplessness unless quasi-locality is allowed)
- eigenstate versus ground state, and certifying a gap
- noise, and approximations inherited from forward solvers
- sample complexity, with no theory for ground states
- symmetry structure

LLM agents already infer symmetries and boundary conditions, choose measurements, write solvers and run model selection. They are unreliable without hard verifiers and depend strongly on the backbone LLM. The opportunity is an agent that proposes, runs certified linear-algebra and DMRG/ED verifiers, and reports negative results honestly, together with a public benchmark, which does not yet exist.

### Cited Findings
- **Non-uniqueness and identifiability.**
  - Single-eigenstate uniqueness holds generically within a range-k ansatz. — [Qi & Ranard](https://arxiv.org/abs/1712.01850)
  - The inverse method returns *spaces* of models; multiple parents can exist. — [Chertkov & Clark](https://arxiv.org/abs/1802.01590); [Dupont et al.](https://arxiv.org/abs/1907.12124)
  - Generalized-parent classes exist for MPS and W states. — [Wang–Gotta–Moudgalya](https://arxiv.org/abs/2611.00203); [Gioia et al.](https://arxiv.org/abs/2510.24713)
  - The SciExplorer agent found wrong models that fit all its measurements. — [SciExplorer](https://arxiv.org/html/2509.24978v6)
  - Kleidermacher et al.'s agent found two comparably good N=4 solutions. In the SnV⁻ case, the highest-evidence models (3–4 emitters) were unphysical and had to be rejected by an instrument-response check. Bayesian evidence alone is therefore not a sufficient verifier. — [Kleidermacher et al.](https://arxiv.org/html/2611.00620v1)
  - Spectral-edge eigenstates are more learnable than mid-spectrum ones. — [Kliczkowski et al.](https://arxiv.org/abs/2605.03043)
  - Identifiability can be quantified through observation geometry. — [Liu et al.](https://arxiv.org/abs/2608.10350)
- **Locality ansatz and operator basis.**
  - Reconstruction is sensitive to undercomplete or overcomplete bases. — [Brito et al. 2024](https://arxiv.org/abs/2311.09302)
  - No simple local spin Hamiltonian stabilizes 2D Gutzwiller-projected states. — [Brito & Marston 2025](https://arxiv.org/abs/2510.15034)
  - Jastrow–Gutzwiller states have only approximate parents. — [Turkeshi & Dalmonte](https://arxiv.org/abs/1909.11327)
  - Strictly local parents of the deformed toric code must be gapless, while quasi-local gapped ones exist. — [Sahay et al.](https://arxiv.org/abs/2503.01977); [Manoj et al.](https://arxiv.org/abs/2607.28740)
  - Learning theory without an ansatz exists for dynamics access. — [Hu et al.](https://arxiv.org/abs/2502.11900); [Bakshi et al.](https://arxiv.org/abs/2405.00082); [Zhao](https://arxiv.org/abs/2410.21635)
- **Eigenstate versus ground state, and the gap.**
  - Covariance methods enforce only zero variance. — [Turkeshi & Dalmonte (snippet)](https://scipost.org/SciPostPhys.8.3.042)
  - PEPS gaps in 2D are generally open. — [Pérez-García et al.](https://arxiv.org/abs/0707.2260)
  - New martingale gap bounds give a simplified proof that MPS parents are gapped. — [Rozmán et al.](https://arxiv.org/abs/2607.19078)
  - Injectivity certifies the gap for VQE parent benchmarks. — [Kobayashi et al.](https://arxiv.org/abs/2109.11759)
- **Noise and solver approximations.**
  - Neural inverse maps inherit the approximations of the forward solver. — [Wu et al.](https://arxiv.org/abs/2611.25321)
  - SciExplorer becomes less consistent with finite shots. — [SciExplorer](https://arxiv.org/html/2509.24978v6)
  - Learning can be made robust to SPAM errors. — [Hangleiter et al.](https://arxiv.org/abs/2108.08319)
  - Local patches keep costs linear. — [Bairey et al.](https://arxiv.org/abs/1807.04564)
- **Sample complexity.**
  - Optimal or polynomial results exist for Gibbs states. — [Anshu et al.](https://arxiv.org/abs/2004.07266); [Haah–Kothari–Tang](https://arxiv.org/abs/2108.04842); [Bakshi et al.](https://arxiv.org/abs/2310.02243)
  - Heisenberg-limited results exist for dynamics. — [Huang et al.](https://arxiv.org/abs/2210.03030)
  - Control-free settings have been studied. — [Gong et al.](https://arxiv.org/abs/2610.06709); [Zhou & Gong](https://arxiv.org/abs/2606.19486)
  - Nothing comparable was found for ground states (Q4).
- **Symmetry structure.**
  - SHC produces symmetric operator bases. — [Chertkov–Villalonga–Clark](https://arxiv.org/abs/1910.10165)
  - Commutant-algebra detection. — [Moudgalya–Motrunich](https://arxiv.org/abs/2302.03028)
  - The SciExplorer agent deduces translational invariance and boundary conditions on its own. — [SciExplorer](https://arxiv.org/html/2509.24978v6)
  - KeplerAgent infers symmetries first. — [KeplerAgent](https://arxiv.org/abs/2602.12259)
  - The Bisognano–Wichmann ansatz has known limits. — [Yang–Ding–Yan](https://arxiv.org/abs/2511.00950)
- **LLM reliability and cost.**
  - GPT-5 recovers models where Gemini 2.5 Pro and GPT-5 nano fail; open-source models struggle; about 2 USD per run. — [SciExplorer](https://arxiv.org/html/2509.24978v6)
  - Variance is concentrated in planning: research-plan F1 0.58–0.86, literature F1 about 0.3–0.4. — [Kleidermacher et al.](https://arxiv.org/html/2611.00620v1)
  - The agent drew a wrong conclusion about spectral diffusion from limited evidence. — [Kleidermacher et al.](https://arxiv.org/html/2611.00620v1)
  - More reasoning effort produced more false positives unless an expected-signal calculation was required. — [Isogawa et al.](https://arxiv.org/abs/2607.25145)
  - LLMs write worse variational parameters than Adam for general Hamiltonians. — [Farias et al.](https://arxiv.org/abs/2611.17246)
  - AI scientists often fail to reproduce published work. — [OmniQ](https://arxiv.org/abs/2611.02833)
  - LLM-SRBench was designed against memorization. — [LLM-SRBench](https://arxiv.org/abs/2504.10415)

### Inferences
Each opportunity is a proposal derived from the findings above, not from a source.
1. **Agentic covariance-matrix loop.**
   - *Proposals.* The LLM proposes symmetry-adapted operator bases (SHC/qosy, commutant detection), with range, sector and closed-form coupling profiles, decides which correlators to request, and simplifies the null-space output symbolically.
   - *Exact numerics.* Covariance matrix → null space → variance → DMRG/ED ground-state and gap check, or martingale bounds for MPS → robustness (smallest nonzero eigenvalue).
   - *Scoring.* Parsimony, with explicit "no local parent in a range-k basis" verdicts for 2D Gutzwiller-type and deformed-toric-code-type states.
2. **Symbolic parents through LLM-SR-style program evolution.** Evolve programs that output term sets and coupling functions J(r), scored by variance plus parsimony. Sanity checks: rediscover Haldane–Shastry from projected-Fermi-sea correlators (Brito & Marston), and AKLT or Majumdar–Ghosh. Use certificates such as ZX identities (Szyniszewski) where possible.
3. **Constrained parents for code and circuit states.** Few-body or native terms, certified gaps, trade-offs between quasi-locality and gap, and scar-type parents (Dooley; Wang–Gotta–Moudgalya; Xu et al.). Verifiers: variance, ED/DMRG gap, injectivity or martingale bounds, weight and term counts. Use reconstructed parents as hardware *certification witnesses* (Nadon et al.).
4. **Measurement-data route.** LLM proposals plus information-gain or RL design policies (Murphy's Model Discovery Agent, SymQNet, Fiderer et al.), local patches (Bairey et al.), and Bayesian structure selection (QMLA; Wallace et al.; Poteshman et al.; Kleidermacher et al.). This targets SciExplorer's under-measurement failure.
5. **Entanglement-Hamiltonian route.** Ansätze beyond Bisognano–Wichmann (Kokail et al. variational; Joshi et al. 51-ion data), verified by relative entropy or fidelity.
6. **QEC hardware.** Agentic selection of coherent-error models from syndrome data (Romanov et al.; Gordeychuk & Kiselev), with exhaustive or simulation-based verifiers.
7. **Benchmark gap.** No public benchmark exists for LLM-agent parent-Hamiltonian reconstruction. A tiered suite could include:
   - exact parents (AKLT, Majumdar–Ghosh, Haldane–Shastry, cluster, toric code)
   - random local Hamiltonians on 8–16 qubits (Qi–Ranard used up to 12)
   - circuit and stabilizer states
   - scars and generalized parents
   - non-Hermitian and Lindbladian cases
   - noisy snapshot data
   - negative instances where no local parent exists

   Verifiers: variance, gap, fidelity and parsimony. Repeated runs would quantify variance across LLM backbones, as Kleidermacher et al. and SciExplorer did.

### Gaps
- No provable sample-complexity result was found for learning a local Hamiltonian from ground-state or eigenstate copies.
- The task suites of OmniQ, QMP-Bench and QSim-Arena were not inspected item by item, so I cannot rule out that some contain Hamiltonian-reconstruction tasks.
- Full texts of most 2026 preprints other than SciExplorer and Kleidermacher et al. were not read; numbers quoted from them come from abstracts or snippets.
