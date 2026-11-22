# MBL fase o regime? Number entropy e lattice gauge theories: stato dell'arte e cosa possiamo dimostrare

*Note di lavoro, 22 novembre 2026.*

## Struttura della cartella

- `lit/`: cinque mappe della letteratura, in inglese:
  - `lit_mbl_debate.md`: il dibattito MBL fase/regime;
  - `lit_mbl_rigorous.md`: i risultati rigorosi;
  - `lit_number_entropy.md`: la number entropy;
  - `lit_lgt.md`: le lattice gauge theories;
  - `lit_math_toolkit.md`: Anderson–Bernoulli 1D/2D, random hopping, random flux.
- `theorems/`: due note con enunciati e dimostrazioni.
  - `number_entropy_theorems.md`;
  - `lgt_theorems.md`.
  - `referee_errata.md`: il report di un referee avversario indipendente, che ha ri-derivato le prove e ri-eseguito i check. **Dove contraddice le note, vale il referee.** Le correzioni sono già incorporate qui sotto.
- `code/`: script Python (numpy/scipy) per tutti i check numerici, con i relativi output `.out`.
  - I file `.npz` delle run l-bit (~30 MB) non sono inclusi.
  - Alcuni script scrivono su path assoluti della scratchpad della sessione.

## Affidabilità delle fonti

In questa sessione arxiv.org era bloccato dal proxy e il budget di ricerca web degli agenti si è esaurito presto. Di conseguenza:

- nessun paper è stato letto per intero;
- le affermazioni bibliografiche vengono da abstract e snippet dei motori di ricerca, oppure dalla memoria (nei file sono taggate `[snippet]`/`[memory]`);
- ho verificato a posteriori l'esistenza e il contenuto essenziale di alcuni lavori chiave: DHMP 2308.08462, DGHP 2408.04338, Aceituno Chávez et al. 2312.13420, KUFS Ann. Phys. 2012.12436, L. Li CMP 2022;
- **prima di citare qualsiasi cosa, va controllata sul testo.**

Le dimostrazioni nelle note invece sono autocontenute e sono state ricontrollate dal referee, quindi su quelle si può lavorare.

---

## 0. TL;DR

1. **MBL come fase in 1D (random, short-range) non è né dimostrato né confutato.** L'unico teorema "full MBL" è quello di Imbrie (2016), condizionato a un'assunzione non dimostrata sulla statistica dei livelli (LLA), per una catena senza U(1).
   - Il consenso 2026 distingue un **"MBL regime"** (prethermal), indiscutibile e descritto dagli l-bit su tutte le scale accessibili, da una **"MBL phase"** asintotica.
   - Se la fase esiste, sta a W ≳ 15–20 secondo gli argomenti di avalanche: ben oltre il W_c ≈ 3.7 dell'ED.
   - Le transizioni viste in ED sono crossover che derivano.
2. **Rigoroso e incondizionato**, per sistemi interagenti a densità di energia positiva, c'è solo:
   - quasi-localizzazione fino a tempi t* ~ exp[c (log W)^{2−ε}] (De Roeck–Huveneers–Meeus–Prośniak 2023);
   - conduttanza termica o(1/L) a ogni temperatura (De Roeck–Giacomin–Huveneers–Prośniak 2024/25), via segmenti non-risonanti rari: è il primo risultato di trasporto a tempo infinito senza LLA.

   A densità di energia nulla (bordo inferiore dello spettro della XXZ random) Elgart–Klein (2022–26) hanno localizzazione, quasi-località e un cono di luce logaritmico.
3. **Number entropy.** Il dato numerico (S_N ≈ (ν/2) ln ln t su tutte le finestre accessibili) non è contestato; l'interpretazione sì. Abbiamo:
   - un **teorema semplice ma pulito**: LIOM quasi-locali e U(1)-invarianti implicano S_N(t) limitata uniformemente in t e L, con costanti esplicite;
   - una **spiegazione esatta** del perché, in un modello l-bit localizzato per costruzione, il valore di saturazione di S_N cresce con L (Aceituno Chávez et al.): è un effetto di taglia finita che converge esponenzialmente a un limite indipendente da L;
   - un **bound incondizionato a tempi finiti**, S_N = O(1) fino a t* ≥ exp[c ln²(W/J)].

   Conclusione: **la number entropy non può decidere la questione con dati a taglia e tempo finiti.** Può farlo solo l'illimitatezza asintotica di S_N (o di Var N_A), e quella è fuori dalla portata della numerica.
4. **LGT.** Non esiste alcun teorema di localizzazione dinamica per una LGT. Ma:
   - la Z2 LGT 1D con materia libera (Smith–Knolle–Kovrizhin–Moessner) è **mappata esattamente** su fermioni liberi con potenziale Bernoulli i.i.d. mediato sui settori;
   - i teoremi di Anderson–Bernoulli 1D rendono quindi **rigorosa la disorder-free localization (DFL)**: localizzazione dinamica mediata, trasporto limitato, imbalance infinito positivo per ogni h ≠ 0, area law per l'entanglement condizionato, number entropy limitata. Questo è un paper concreto e fattibile;
   - **controesempio**: certi stati iniziali di gauge producono una legge di settore a dimeri correlati e la DFL fallisce (energie critiche, trasporto superdiffusivo). Quindi la DFL è una proprietà dello stato iniziale, non solo dell'Hamiltoniana;
   - nei **quantum link models U(1)** la legge di Gauss rende la number entropy limitata *cinematicamente* (S_N ≤ ln(2S+1) o ln(4S+1)) per qualunque dinamica: lì non è un diagnostico di localizzazione;
   - la DFL interagente in 1D **è** MBL con disordine binario, quindi rigorosamente si arriva solo a statement pretermici (ADHH).

---

## 1. MBL: fase termodinamica o regime?

### 1.1 Cosa viene affermato

**Campo "fase" (maggioritario).**
- Esiste una fase MBL in 1D a disordine forte, descritta da LIOM (l-bit) quasi-locali.
- La transizione è di tipo Kosterlitz–Thouless e guidata da avalanche: le RG di Goremykina–Vasseur–Serbyn, Dumitrescu et al., Morningstar–Huse(–Imbrie) danno ξ ~ exp(c/√|W−W_c|).
- Il W_c vero è molto più grande di quello ED: Morningstar–Colmenarez–Khemani–Luitz–Huse 2022 stimano W_c ≳ 20 per la catena di Heisenberg con h ∈ [−W, W].
- Riferimenti: Abanin et al., Ann. Phys. 2021; Sierant–Lewenstein–Zakrzewski PRL 2020; Long–Crowley–Khemani–Chandran PRL 2023.

**Campo "regime" (scettici).**
- Šuntajs–Bonča–Prosen–Vidmar (PRE 2020): t_Th ∝ L² e^{W/Ω}, quindi W*(L) ∝ L e nessuna MBL nel limite termodinamico.
- Sels–Polkovnikov (PRE 2021): il tempo di rilassamento cresce ~ e^{W} ma resta finito a ogni W studiato.
- Sels (PRB 2022): delocalizzazione indotta da bagno fino a W ~ 20.
- Evers–Bera e collaboratori; Prelovšek–Mierzejewski.
- Gruppo Sirker (number entropy, 2020–2026): sostengono lento trasporto di particelle e assenza di vera MBL.

**Cosa concedono entrambi (punto di convergenza, 2023–26).**
- W_c^ED ≈ 3.7 non è la transizione.
- Esiste un regime con tempi di rilassamento enormi ma (forse) finiti, ben descritto dagli l-bit.
- Avalanche e risonanze sono i meccanismi.
- Nessun dato numerico o sperimentale arriva al regime asintotico W ≳ 20.
- La review di Sierant–Lewenstein–Scardicchio–Vidmar–Zakrzewski (Rep. Prog. Phys. 2025), scritta da autori di entrambi i campi, tratta la fase come questione aperta.

### 1.2 Cosa è dimostrato

| Risultato | Cosa dice | Condizioni |
|---|---|---|
| Imbrie, J. Stat. Phys. 2016 | Esiste un unitario quasi-locale che diagonalizza H: LIOM, nessun trasporto per ogni t | Catena Ising con campo trasverso debole, **senza U(1)**, disordine continuo; **assunzione LLA non dimostrata** |
| Hamza–Sims–Stolz 2012; Abdul-Rahman–Nachtergaele–Sims–Stolz 2016; Abdul-Rahman–Stolz 2015 | Catene XY/XX disordinate (fermioni liberi): bound di Lieb–Robinson a velocità zero, entanglement e trasporto limitati per ogni t uniformemente in L, area law per tutti gli autostati | Incondizionato (localizzazione a 1 corpo) |
| Elgart–Klein–Stolz 2018; Beaud–Warzel; Elgart–Klein 2022–2026 | XXZ random in fase Ising: localizzazione, quasi-località, propagazione lenta, cono di luce logaritmico | **Solo finestre di energia al bordo inferiore** (densità di energia nulla), uniforme in N e L |
| De Roeck–Huveneers–Meeus–Prośniak, Physica A 2023 (2308.08462) | Insieme completo di quantità quasi-locali con dinamica trascurabile fino a t* ~ exp[c (log W)^{2−ε}] | Disordine forte, incondizionato, **tempi finiti** |
| De Roeck–Giacomin–Huveneers–Prośniak, 2408.04338 | Conduttanza termica di una catena di lunghezza L = o(1/L) a ogni T (assenza di conduzione normale) | Incondizionato, senza LLA; segmenti completamente MBL con probabilità ≥ e^{−cL} fanno da colli di bottiglia |
| Mastropietro 2015–2021 | Fermioni quasi-periodici interagenti: correlazioni dello stato fondamentale localizzate, Drude weight nullo | **Solo T = 0** |
| De Roeck–Hannani 2025; Yin–Nandkishore–Lucas 2024 | Modelli giocattolo (Quantum Sun; geometrie expander) con MBL e statistica di Poisson dimostrate | Non sono catene a corto raggio |
| ADHH (CMP 2017), Else–Bauer–Nayak, Yin–Lucas, Yin–Surace–Lucas | Prethermalizzazione: leggi di conservazione emergenti fino a tempi ~ exp(c U/J) | Senza disordine; tempi finiti |

**Nessun teorema di delocalizzazione** per sistemi disordinati interagenti. Avalanche, instabilità in d ≥ 2, interazioni a lungo raggio, continuo e SU(2) sono tutti argomenti euristici o numerici. Per SU(2) il nucleo algebrico è elementare: non esistono LIOM U(1)-tipo SU(2)-invarianti.

### 1.3 Argomenti plausibili

**Pro fase:**
- conteggio delle risonanze: ζ(W) ~ 1/ln(W/J) → 0, quindi prima o poi ζ < ζ_c e le avalanche si arrestano;
- l'RG KT con linea di punti fissi MBL è internamente coerente;
- Imbrie condizionato;
- nelle catene quasi-periodiche manca il seme raro delle avalanche.

**Contro:**
- tutti gli stimatori derivano con L, senza crossing convergente;
- il lower bound su W_c continua a salire (3.7 → 5.4 → ≳ 20);
- LLA non dimostrata;
- tutti i dati a W ≲ 10 sono pre-asintotici, anche per la stessa teoria delle risonanze.

**Ostacolo tecnico centrale per una prova:** la statistica dei livelli con pochi variabili random (~2^n livelli, n variabili). È lo stesso muro di Imbrie (LLA) e di Elgart–Klein (niente LIOM).

---

## 2. Number entropy

### 2.1 Il dibattito

- **KUFS** (Kiefer-Emmanouilidis, Unanyan, Fleischhauer, Sirker), PRL 124, 243601 (2020) e seguiti (PRB 2021, Ann. Phys. 2021, SciPost 2022):
  - S_N ≈ (ν/2) ln ln t, cioè Var N_A ~ (ln t)^ν, senza segni di saturazione;
  - lo interpretano come trasporto di particelle lentissimo;
  - con la Hartley number entropy vedono crescita tipica e non da regioni rare;
  - nel 2026 parlano di "path-entropy obstruction".
- **Luitz–Bar Lev**, PRB 2020: le distribuzioni dei valori saturati di S_N convergono in L, con cutoff ≈ una particella che fluttua attraverso il taglio. La crescita di S è puramente configurazionale.
- **Ghosh–Žnidarič**, PRB 2022: risonanze a due livelli; avvicinamento a legge di potenza alla saturazione, che non cresce con L.
- **Aceituno Chávez, Artiaco, Klein Kvorning, Herviou, Bardarson**, PRL 2024 (2312.13420). Il modello è un l-bit **localizzato per costruzione**: H = U H′ U†, con U un circuito di profondità finita di gate a due siti U(1)-invarianti e H′ diagonale con interazioni esponenzialmente decrescenti. In questo modello:
  - S_N cresce come ln ln t;
  - il valore di saturazione cresce con L fino a L = 28.

  Conclusione degli autori: queste osservazioni non bastano a escludere la MBL.
- **Altri:** Sierant–Zakrzewski (TDVP, L = 50); Sharifan–Kiefer-Emmanouilidis–Sirker 2026 (Aubry–André); Faulend–Buljan–Štrkalj 2026 (quantum beats); Weisse–Gerstner–Sirker 2025 (non convergenza dell'unitario LIOM). Tutti visti solo come snippet.

### 2.2 Cosa possiamo dire rigorosamente

Dettagli in `theorems/number_entropy_theorems.md`; tutte le prove sono verificate dal referee.

**Teorema A (T1): LIOM U(1) implicano S_N limitata.**

Ipotesi: H = U D U†, con D diagonale nella base delle occupazioni e [U, N] = 0. Sia δ_A := ‖U† N_A U − N_A‖. Allora:
- ‖N_A(t) − U N_A U†‖ = δ_A esattamente, per ogni t;
- ‖N_A(t) − N_A‖ ≤ 2δ_A;
- per uno stato di Fock (es. Néel), Var N_A(t) ≤ 4δ_A²;
- S_N(t) ≤ ½ ln[2πe(4δ_A² + 1/12)], uniformemente in t e L.

L'ultimo passo usa il bound di massima entropia per variabili intere, H(X) ≤ ½ ln 2πe(σ² + 1/12) (Massey; Cover–Thomas).

Costanti esplicite di δ_A:
- **circuito brickwork** di profondità d_u: δ_A ≤ d_u |∂A|. In questo caso p_t è supportata in n₀ ± 2 d_u |∂A|, quindi *tutte* le Rényi, Hartley inclusa, sono ≤ ln(4 d_u |∂A| + 1);
- **generatori quasi-locali** con code e^{−r/ξ}: δ_A ≤ J ξ² |∂A|;
- **prodotti random di rotazioni locali** (stile Imbrie, inclusi blocchi risonanti): bound in probabilità e in media, per cui S_N richiede solo un momento logaritmico della costante.

Estensioni:
- **T1′**: l'U(1)-invarianza di U non serve. Basta che i LIOM τ_x = U n_x U† commutino con N, il che è automatico se [H, N] = 0 e D ha spettro semplice. In quel caso il bound è in termini delle code di U† n_x U (non dei LIOM stessi: correzione del referee).
- **Code:** basta il primo momento Σ_j (j+1) η(j) < ∞; ad esempio code a legge di potenza j^{−β} con β > 2.

Novità: gli ingredienti sono standard. Il lemma "un unitario U(1) quasi-locale sposta carica solo al bordo" è nello stile di Bachmann–Bols–De Roeck–Fraas e Kapustin–Sopenko. L'enunciato è folklore euristico (KUFS stessi scrivono "since the conserved charges are local, S_N has to be bounded"), ma **non l'abbiamo trovato scritto come teorema**.

**Proposizione B: perché la saturazione cresce con L in un modello l-bit localizzato.**

Per circuiti l-bit con D a spettro semplice, sia p̄ la distribuzione defasata (media a tempo infinito di p_t) di N_A. Allora:
- p̄ è **indipendente da D**: dipende solo dagli autovettori U|c⟩;
- p̄ dipende solo dai gate entro 4 d_u − 2 siti dal taglio, quindi la sua legge è **esattamente indipendente da L** per L ≥ 8 d_u − 4;
- per ogni realizzazione, la media temporale Ŝ_N(L) di S_N(t) è ≤ H(p̄), per concavità.

Numerica in un modello esattamente nella classe T1:
- Ŝ_N(L) cresce con L, mentre E H(p̄) resta piatta;
- il gap E[H(p̄) − Ŝ_N] **si dimezza ogni due siti**;
- dopo il transiente iniziale S_N cresce lentamente per molte decadi.

Si riproducono quindi le due feature qualitative di Chávez et al. in un modello dove la limitatezza è un teorema. Interpretazione: le code di U fissano la *taglia* delle fluttuazioni, mentre le code di D fissano soltanto *l'orologio* e, a L finito, le fluttuazioni temporali residue.

Caveat del referee:
- il nostro D contiene termini a k corpi per ogni k, quindi le fasi sono razionalmente indipendenti; con D a pochi corpi (come forse in Chávez) va ricontrollato;
- non abbiamo potuto leggere le Eq. (8)–(10) di Chávez, quindi la loro appartenenza alla classe T1 è dedotta solo dall'abstract;
- i numeri 4.88 e 4.17 valgono per catena aperta e taglio a metà.

**Teorema C (T3): stati gaussiani U(1).**

N_A segue una distribuzione Poisson-binomiale con parametri ν_k (gli autovalori di C_A). Valgono:
- S ≥ 4 ln 2 · Var;
- 2 Var ≤ S^{(2)} ≤ 4 ln 2 · Var;
- ½ ln(8 Var/π) ≤ S_N ≤ ½ ln[2πe(Var + 1/12)];
- S_N = ½ ln(2πe Var) + o(1) per Var → ∞.

Quindi:
- **S_N = ½ ln S^{(2)} + O(1)** in modo bilatero (essenzialmente già in KUFS, SciPost 2020);
- **S_N ≤ ½ ln S + O(1)** solo come bound superiore, perché si possono avere S → ∞ e S_N → 0 (correzione del referee: "S_N ~ ½ ln S" non è un teorema bilatero);
- per stati interagenti U(1) l'unico bound universale è S_N ≤ S, ed è ottimale.

**Teorema D (T4): bound incondizionato a tempi finiti.**

Per la catena t-V/XXZ con campi random h_x ∈ [−W, W] si costruisce una carica vestita Ñ_n localizzata al taglio, con:
- ‖Ñ_n − N_A‖ ≤ 2|∂A| J/γ;
- ‖[H, Ñ_n]‖ ≤ |∂A| J n! (10 J/γ)^n;
- probabilità ≥ 1 − |∂A| 2^{8n−2} γ/W.

Ne segue S_N = O(1) fino a t* ≥ J⁻¹ exp[(ln x)²/(32 ln 2)(1 − o(1))], con x = pW/(10 |∂A| J), uniformemente in L, con probabilità 1 − p.

È un analogo elementare, su un solo taglio, di DHMP. Le costanti sono inutili a W fisici (servono W ≳ 10³ J). Il suo valore è dire esattamente cosa manca: o un controllo multiscala delle risonanze a tutti gli ordini vicino al taglio, oppure un cono di luce logaritmico a densità di energia positiva.

**Corollario condizionale.** Un cono di luce logaritmico implica rigorosamente S_N ≤ 2 ln ln t + O(1) (ln ln t con un cono di luce stretto). Quindi **ln ln t è il massimo compatibile con un log light cone**, e i dati osservati ci stanno dentro.

### 2.3 Logica: cosa si può e non si può concludere

- **Valido:** S_N illimitata (t → ∞ a L = ∞, oppure il valore saturato per L → ∞) esclude LIOM U(1) con costanti uniformi; per Var N_A serve il momento secondo.
- **Non valido:**
  - "S_N cresce a (L, t) finiti, quindi qualcosa di asintotico";
  - "S_N cresce, quindi conducibilità σ_dc > 0";
  - "S_N limitata, quindi localizzazione": è necessario ma non sufficiente; frammentazione e sistemi Stark la limitano comunque.
- **I valori osservati sono minuscoli.** S_N ≈ 0.3–1.0 richiede solo δ_A ≳ 0.08–0.3, molto sotto qualunque costante naturale. Dati a taglia finita non possono falsificare la classe T1.
- **Diagnostici migliori** (limitati uniformemente in L sotto T1):
  - (D1) max (o quantili) sugli autostati di Var_E(N_A) in un settore a N fissato: test statico, accessibile con ED completa;
  - (D2) varianza ed entropia della distribuzione defasata p̄_L (ensemble diagonale) al posto di S_N(t; L) tardivo;
  - (D3) il caso peggiore sugli stati iniziali, ‖N_A(t) − N_A‖.

---

## 3. Lattice gauge theories

### 3.1 Letteratura (sintesi)

- **DFL libera:**
  - Smith–Knolle–Kovrizhin–Moessner, PRL 118 (2017), arXiv:1701.04748: Z2 in 1D, fermioni liberi in un disordine binario "dinamico";
  - Smith et al. PRB 97 (2018) in 2D: Bernoulli e/o random π-flux, percolazione quantistica.
- **DFL interagente:**
  - Smith et al. PRL 119 (2017): XXZ in campo binario, "quantum disentangled liquid";
  - Brenes–Dalmonte–Heyl–Scardicchio PRL 2018 (Schwinger);
  - Karpov et al. PRL 2021 e Chakraborty et al. 2022/23 (QLM 2D, frammentazione/percolazione);
  - Osborne–McCulloch–Halimeh; Cataldi et al. 2026 (SU(2)).
- **Esperimento Google:** Gyawali et al., 2410.06557, pubblicato su Science 2026 secondo gli snippet.
- **Reinterpretazioni recenti:**
  - Jeyaretnam et al. 2025: la DFL di Schwinger viene da frammentazione approssimata;
  - Sala–Giudici–Halimeh 2024: nei QLM è un effetto classico da automa cellulare;
  - Tian et al. 2026: il disordine binario non sopprime le risonanze.
- **Rigoroso sulle LGT:**
  - **nessun teorema di localizzazione dinamica per una LGT** (trovato solo lavoro di equilibrio, es. Lieb flux phase, Goller–Porta 2025);
  - gauge protection rigorosa solo dove si applica ADHH (penalità a spettro intero), in appendici;
  - Yin–Lucas e Yin–Surace–Lucas per il decadimento del falso vuoto.

### 3.2 Cosa possiamo dimostrare

Dettagli in `theorems/lgt_theorems.md`. Mappature verificate con ED a ~1e−14.

**Teorema E: dualità esatta.**
- Modello: Z2 LGT 1D con termine elettrico a due link, −h Σ σ^x_{b_{j−1}} σ^x_{b_j} (è la forma che rende il modello libero; va confrontata con le equazioni di SKKM).
- È unitariamente equivalente a fermioni liberi accoppiati a cariche statiche q_j = ±1, con h_q = −J Δ + 2h diag(q).
- C'è un unitario esplicito, locale attraverso il taglio (un solo gate controllato).
- Con un link di bordo e gauge in |↑…↑⟩_z, le cariche sono **esattamente i.i.d. Bernoulli(1/2)**, per qualunque stato di materia.

**Teorema F: DFL rigorosa** (condizionata a (DL_N), cioè E_q sup_t |⟨x|e^{−i h_q t}|y⟩| ≤ C e^{−η|x−y|} uniformemente in N):
- (a) localizzazione dinamica mediata sui settori;
- (b) trasporto limitato attraverso ogni taglio;
- (c) **imbalance a tempo infinito strettamente positivo per ogni h ≠ 0**, uniformemente in N. Si basa su un'identità nuova, E_α a_α = 2h⟨φ_α, VΓφ_α⟩, più Khintchine–Szarek; nel bound va usato J + |h|;
- (d) area law uniforme nel tempo per l'entanglement *condizionato* S_G, e per l'entanglement completo materia+gauge per stati iniziali a carica definita;
- (e) number entropy limitata uniformemente nel tempo (costante del teorema ≈ 6.9).

Stato dell'input matematico:
- Ge–Zhao (2020) dà (DL) su Z per Bernoulli, a tutte le energie e per ogni h ≠ 0; per Bernoulli i.i.d. non ci sono energie critiche (teorema di Furstenberg);
- **gap (referee):** il nostro modello ha un bordo fisso, quindi serve la versione su semiretta o su intervallo finito, uniforme in N;
- candidato naturale: Macera–Sodin, AHP 2022 (2110.00097), correlatori di autofunzioni per operatori quasi-1D su intervalli finiti. Da verificare.

Precisazioni (referee):
- il bound di Lieb–Robinson a velocità zero vale solo **mediato sui settori** (dipende dallo stato). In norma operatoriale fallisce perché il settore pulito q ≡ +1 è balistico;
- S_G può essere < S_N (la Prop. 5(3) della nota è sbagliata per S_G).

**Proposizione G: la DFL può fallire.**
- Stato di gauge prodotto: link dispari in |↑⟩_z, link pari in autostati σ^x alternati.
- Produce coppie di cariche uguali (q_{2k−1}, q_{2k}) = (ε_k, ε_k), cioè il **random dimer model** con energie ±2h.
- Per h < J/2 ci sono energie critiche ±2h. Via Jitomirskaya–Schulz-Baldes–Stolz (le ipotesi di genericità vanno controllate) il trasporto è superdiffusivo; numericamente ⟨x²⟩ ~ t^{1.35–1.5}.
- Conclusione: **la DFL è una proprietà della legge di settore indotta dallo stato iniziale.**

**Entropie.** Struttura esatta delle entropie: entanglement condizionato S_G, completo S_full, materia mescolata S_mix e Holevo χ = S_mix − S_G. Questo rende esatto il punto "superposizione vs campionamento" dell'esperimento Google. Congetture con numerica:
- S_full ~ ln t (struttura l-bit esatta);
- S_mix a volume law.

**2D.**
- Mappatura esatta su Bernoulli 2D (gauge z-polarizzato) o su **random π-flux binario** (gauge x-polarizzato con CDW).
- Coupling forte (h > 2J): confinamento nei cluster di carica uguale fino a t ~ h/J².
- Oltre quella scala c'è un problema di percolazione quantistica supercritica, aperto.
- Risultati matematici esistenti: Ding–Smart (bordo dello spettro) e L. Li, CMP 393 (2022) (disordine grande, fuori da finitamente molte energie). Il π-flux binario è aperto.

**U(1).**
- **Quantum link models a spin S:** sommando la legge di Gauss su A si ottiene N_A = E_b − E_{a−1} + cost. Quindi, per **qualunque** stato in un settore di carica, S_N ≤ ln(2S+1) (taglio al bordo) oppure ≤ ln(4S+1) (intervallo nel bulk). Lì la number entropy satura per ragioni cinematiche.
- **Link di Wilson:** la media sui tagli di S_N è limitata uniformemente nel tempo dalla conservazione dell'energia (termine g²E²).
- **Pretermico (corollari di ADHH, enunciato da verificare):** conservazione emergente dell'energia elettrica o del numero di stringhe per g² ≫ w, e del numero di particelle per m ≫ w, fino a exp(c·scala/log³). Gauge protection con penalità quadratica.
- **Una carica nella Z2 standard** (termine a un link): è una particella in un potenziale random walk V(x) ~ √|x|.
  - Simon–Spencer escludono lo spettro assolutamente continuo;
  - la localizzazione è aperta: incrementi indipendenti ma valori correlati, degenerazioni esatte per incrementi interi.
- **DFL interagente 1D** = XXZ in campo binario: è MBL con disordine binario, nemmeno coperta dal framework di Imbrie (che richiede disordine continuo).

---

## 4. Dove possiamo fare qualcosa (ordinato per rapporto payoff/difficoltà)

1. **Paper "Rigorous disorder-free localization in a Z2 lattice gauge theory".** Difficoltà bassa–media, payoff alto.
   - Contenuto: Teoremi E, F, G, più la tassonomia delle entropie.
   - Lavoro mancante:
     - fissare l'input (DL_N) uniforme in N con bordi, via Macera–Sodin o rifacendo le large deviations su intervalli;
     - verificare Ge–Zhao e Jitomirskaya–Schulz-Baldes–Stolz sul testo;
     - confrontare il modello con SKKM.
   - Sarebbe il primo teorema di DFL.
2. **Nota rigorosa sulla number entropy.** Difficoltà bassa, payoff medio–alto per il dibattito.
   - Contenuto: Teorema A (con T1′), Proposizione B (indipendenza da D e località di p̄; Ŝ_N ≤ H(p̄); gap che si chiude in L), Teorema D e la logica del §2.3.
   - Proposta di diagnostici statici D1/D2 da applicare ai dati XXZ esistenti (ED fino a L ≈ 20–24).
   - Da fare prima: leggere davvero 2312.13420 (Eq. 8–10) e testare la Prop. B con D a pochi corpi.
3. **Entanglement logaritmico rigoroso in un sistema senza disordine:** S_full ≤ C ln t nella Z2 1D. Difficoltà media–alta, payoff alto.
   - La struttura l-bit è esatta e la stima di sensitività di Duhamel è già pronta.
   - Manca un argomento di decoupling telescopico per una sovrapposizione di 2^N settori.
4. **Assenza di catastrofe di ortogonalità** per un flip di due cariche nella catena Bernoulli localizzata, che dà memoria rigorosa delle osservabili di spin ⟨σ^z_b(t)⟩.
5. **ADHH per LGT, scritto bene:** gauge protection e conservazione emergente nei QLM e nella DFL a h forte, in qualunque dimensione e con interazioni. È formalizzazione, utile.
6. **Potenziale random walk** (carica singola nella Z2 standard, mesone U(1)): prima incrementi continui (Wegner per traslazione rigida + MSA), poi quelli interi.
7. **Free-fermion e catena del codice (random hopping J_i = u_i^δ, E = 0 critica).** Non esiste uno studio dedicato della number entropy. Stima RG: Var ~ ln ln t, quindi S_N ≲ ½ ln ln ln t, praticamente invisibile. Con il teorema C, S_N si ottiene gratis dagli autovalori ν_k già calcolati. Attenzione: per hopping puramente bipartito l'imbalance infinito da CDW è **identicamente zero** (chiralità), localizzato o no: serve un altro diagnostico.
8. **Aperti e molto difficili:** LLA per blocchi quasi-liberi (Anderson XX + piccola Δ); π-flux binario 2D; DFL interagente; transizione di percolazione quantistica.

---

## 5. Osservazioni sul codice Fortran (`src/`), senza modifiche

- `init_state` riempie solo `u(2i, i)` per *i pari*. Lo stato iniziale è quindi |0001 0001 …⟩: densità 1/4 e periodo 4, con L/4 particelle e non L/2; le colonne dispari di `umat0` sono nulle. Se si voleva il Néel, è un bug.
- `build_ham` con `i = ls` scrive `hmat(ls, ls+1)` e `hmat(ls+1, ls)`, **fuori dai limiti** dell'array (possibile corruzione di memoria). Il termine di bordo viene poi assegnato a parte per `bc = .true.`.
  - `-fcheck=all` è sintassi gfortran; con ifort serve `-check all`.
- Efficienza: `expm` (Padé) per ogni tempo costa O(L³) per punto. Conviene diagonalizzare una volta per realizzazione.
- Accuratezza: per hopping su molte decadi (J = u^δ), gli autovalori vicini a zero vanno calcolati con accuratezza *relativa*: matrice chirale, quindi SVD bidiagonale (Demmel–Kahan, dqds).
- Var N_A = Σ ν(1−ν) e S_N (convoluzione Poisson-binomiale) costano quasi nulla, partendo dagli autovalori già calcolati in `entropy`.

---

## 6. Fonti verificate in questa sessione (link)

- De Roeck, Huveneers, Meeus, Prośniak, *Rigorous and simple results on very slow thermalization…*: https://arxiv.org/abs/2308.08462
- De Roeck, Giacomin, Huveneers, Prośniak, *Absence of normal heat conduction…*: https://arxiv.org/abs/2408.04338 (commento JCCM: https://condmatjclub.org/uploads/2024/09/JCCM_September_2024_02.pdf)
- Aceituno Chávez et al., *Ultraslow growth of number entropy in an ℓ-bit model*: https://arxiv.org/abs/2312.13420
- Kiefer-Emmanouilidis et al., PRL 124, 243601: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.243601 ; Ann. Phys.: https://arxiv.org/abs/2012.12436 ; https://arxiv.org/abs/2010.00565
- Kim, Chandran, Abanin, *Local integrals of motion and the logarithmic lightcone*: https://arxiv.org/abs/1412.3073
- L. Li, *Anderson–Bernoulli localization at large disorder on the 2D lattice*, CMP 393, 151 (2022): https://arxiv.org/abs/2002.11580
- Macera, Sodin, *Anderson localisation for quasi-one-dimensional random operators*: https://arxiv.org/abs/2110.00097
- Macera, *Non-Lyapunov annealed decay for 1d Anderson eigenfunctions*: https://arxiv.org/abs/2307.07626
- Gyawali et al., *Observation of disorder-free localization using a (2+1)D LGT on a quantum processor*: https://arxiv.org/abs/2410.06557

Tutto il resto (centinaia di riferimenti) è nei file `lit/`, con il livello di verifica indicato per ciascuno.
