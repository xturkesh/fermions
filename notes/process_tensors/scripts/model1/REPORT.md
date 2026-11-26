# Modello 1 — collision model Haar con memoria e rumore sull'ambiente

Tutti i numeri sono riproducibili con gli script in questa cartella (`common.py`, `exact.py`, `t0_verify.py`, `t0b_identity.py`, `t1_noise.py`, `t1b_scaling.py`, `t2_quenched.py`, `t3_hs.py`, `t4_spectrum_exact.py`, `t4b_c1.py`). Ogni formula è stata verificata contro Monte Carlo esatto (process tensor completo) a taglie piccole.

## 1. Modello e notazione

- Sistema S (dim d), ambiente E = L qudit di dim q, D = q^L, ρ_E(0) = I/D.
- Step: input fresco di S (gamba di Choi), unitaria Haar U su S⊗E, poi canale di rumore N su ogni qudit di E (Kraus {K_i}).
- Υ_n = stato di Choi normalizzato sulle gambe (R_1 O_1 … R_n O_n); Υ_p, Υ_f = marginali su passato (n1 step) e futuro (n2 step).
- Media su k repliche: dopo ogni step Haar, l'operatore a k repliche su E sta in span{P_σ, σ∈S_k}; si propagano i momenti m_σ = Tr(P_σ^{-1} Y) con
  y_β = Σ_α Wg(α^{-1}β, dD) d^{#(α^{-1}π_in^{-1})} d^{#(π_out β)} m_α,  m'_γ = Σ_β G_γβ y_β,  G_γβ = Tr(P_γ^{-1} N^{⊗k}(P_β)) (fattore per qudit)^L.
  π_in, π_out sono le contrazioni delle gambe di S (I = traccia, F = swap, cicli per Tr Υ^k). La chiusura vale per **qualsiasi** canale (anche non unitale), perché ogni step contiene un Haar globale su S⊗E.

## 2. k = 2: forme chiuse (Task 1)

Costanti per qudit: g_I = Tr N(I)², g_F = Σ_ij |Tr K_i†K_j|² (verificate contro la matrice di permutazioni); G_I = g_I^L, G_F = g_F^L.

- depolarizzante: g_I = q, g_F = u q² + 1 − u (u = (1−p)²)
- dephasing (qubit): g_I = 2, g_F = 4[(1−p)² + p²]
- amplitude damping γ (non unitale): g_I = 2(1+γ²), g_F = (2−γ)² + γ²

Con a = (D²−1)/(d²D²−1), b = D(d²−1)/(d²D²−1), c_I = (D G_I − G_F)/(D(D²−1)), c_F = (D G_F − G_I)/(D(D²−1)):

T_FF = [[a, b], [c_I a, c_I b + c_F]] (gambe di S in swap),  T_II = [[1, 0], [c_I + c_F b, c_F a]] (gambe tracciate),  v0 = (1, Tr ρ_E(0)²).

- P(Υ_n) = e·T_FF^n v0;  P(Υ_p) = P(Υ_{n1}) (causalità);  P(Υ_f) = e·T_FF^{n2} T_II^{n1} v0.
- I2(n1,n2) = log2[P(Υ)/(P(Υ_p)P(Υ_f))] (Rényi-2 annealed).
- Senza rumore: P(n) = a^n + (1−a^n)/D², quindi una "curva di Page temporale": I2 → 2 log2 D (massimo permesso dalla bond dimension) per n1, n2 ≳ log_d D, e I2 ≈ (d^{2n1}−1)(d^{2n2}−1)/(D² ln 2) per D² ≫ d^{2n}. Esempi d=2, n1=n2=1,2,3,4,6,10: D=16 → 0.05, 0.84, 3.47, 6.03, 7.83, 8.00 bit.
- Tempi lunghi con rumore (forma chiusa): I2_∞ = −log2[(e·r_+)(l_+·w)], con r_+, l_+ autovettori dominanti di T_FF (l·r = 1) e w = (1, y*), y* = (c_I + c_F b)/(1 − c_F a).

### Transizione di Markovianizzazione

Due livelli "nudi": a ≈ 1/d² (settore Markoviano: l'ambiente agisce come bagno fresco) e c_F + c_I b ≈ (g_F/q²)^L (settore di memoria). Si incrociano a

L* ≈ ln(d²)/ln(q²/g_F)

(verificato: depol. p=0.1 → 9.01 vs 9.02; dephasing p=0.25 → 2.94 vs 2.95; amp. damping γ=0.3 → 4.70 vs 4.71). Interpretazione: la transizione avviene quando la produzione di entropia del settore di memoria, L·log2(q²/g_F) bit per step, eguaglia il tasso massimo di entropia del Choi, 2 log2 d.

- **Limite di scala** p → 0, L → ∞ a x = pL fisso: I2_∞/(2 log2 D) → Θ(x* − x), salto netto, con x* = ln(d²)/κ, κ = −∂_p ln(g_F/q²)|_0:
  depolarizzante x* = ln(d²)/(2(1−q⁻²)) = 0.924 (d=q=2); dephasing x* = ln 2; amplitude damping x* = ln 4.
  Tabella (depol.): a x = 0.90 → 0.092, 0.798, 0.936, 0.969 per p = 0.1, 0.03, 0.01, 0.005; a x = 0.95 → 0 per tutti i p.
- **Gap all'avoided crossing**: Δλ_min ≈ 2√(a b c_I) ≈ 2√(d²−1)/(d² D*), esponenzialmente piccolo in L* (p = 0.1: 1.447e-3 misurato, 1.446e-3 predetto). A L finito la caduta di I2 inizia qualche unità di L prima di L* (singolarità logaritmica 2 log2(1 − a/c_F)).

## 3. Quenched vs annealed (Task 2)

Propagazione esatta per campione delle tre purezze (d=q=2, n1=n2=8, 300 campioni per L≤2, 120 per L=3):

| L | p | I2 annealed | ⟨I2⟩ quenched | dev. std |
|---|---|---|---|---|
| 1 / 2 / 3 | 0.30 | 0.976 / 1.267 / 0.471 | 0.978 / 1.269 / 0.472 | 0.12 / 0.05 / 0.01 |
| 1 / 2 / 3 | 0.42 | 0.611 / 0.425 / 0.059 | 0.636 / 0.430 / 0.059 | 0.12 / 0.04 / 0.003 |
| 1 / 2 / 3 | 0.55 | 0.316 / 0.121 / 0.014 | 0.339 / 0.123 / 0.014 | 0.09 / 0.02 / 0.002 |

- Quenched ≈ annealed entro 0.03 bit, con fluttuazioni che calano con L (auto-media). Il collasso della memoria (andamento non monotono in L) c'è nei campioni tipici.
- Senza rumore, per ogni campione I2 = 2 log2 D esattamente (dev. std 0): a tempi lunghi lo stato ridotto di E⊗(riferimento) tende a I/D² campione per campione.
- von Neumann (taglie minime, n1=n2=2, Υ completo): stesso andamento. p = 0.2: L=1 0.90 bit, L=2 0.70; p = 0.42: 0.40 contro 0.16; senza rumore L=2 > L=1.
- Verdetto: per L ≤ 3 il fenomeno sopravvive nel quenched e in von Neumann. La nitidezza a L → ∞ per la quantità tipica (limite di replica k → 1) resta da dimostrare.

## 4. k = 3, 4: distanza dalla Markovianità (Task 3)

δ_HS = E‖Υ − Υ_p⊗Υ_f‖²₂ / E Tr Υ² richiede A = E Tr Υ² (k=2), B = E Tr[Υ(Υ_p⊗Υ_f)] (k=3: (12) sul passato, (13) sul futuro), C = E[Tr Υ_p² Tr Υ_f²] (k=4: (12) e (34)).

- **Identità (rumore unitale, ρ_E(0) = I/D): B = C**, esatta in aritmetica razionale.
  Schizzo: la replica col passato tracciato vede E = I/D al taglio; l'ultima unitaria del passato assorbe un twirl V su E, che sostituisce la contrazione col futuro con Tr(A_E)/D.
  Quindi **δ_HS = 1 − 2^{−I2}**: per rumore unitale k = 3, 4 non aggiungono nulla a k = 2. Fallisce per amplitude damping (D=2, n1=2, n2=1, γ=0.5: B = 0.10517, C = 0.10453, confermato da MC).
- Grande D: δ_HS ≈ (d^{2n1}−1)(d^{2n2}−1)/D² (n1=n2=1, D=256: δ D²/d^{2n} = 0.5624 contro 0.5625).
- Distanza di traccia (MC, n1=n2=1): T = 0.312, 0.158, 0.0795 per D = 4, 8, 16, cioè T ≈ 1.27/D. Operativamente quasi Markoviano per D ≫ d^n, coerente con la tipicità di Figueroa-Romero–Modi–Pollock. A tempi lunghi invece δ_HS → 1 − D⁻² (memoria massima).
- Con rumore δ_HS collassa oltre L* (p = 0.3, n1=n2=6, L = 1…5: 0.49, 0.57, 0.25, 0.023, 0.0018).

## 5. Spettro e free probability (Task 4)

Confronto di E Tr Υ^k (k = 2, 3, 4, aritmetica esatta) con la misura indotta (free Poisson / Marchenko–Pastur) con N = d^{2n}, M = D²: E Tr ρ^k = Σ_σ N^{#σ} M^{#(σ⁻¹γ)} / Π_{j<k}(NM + j).

- **Senza rumore: free Poisson a grande D, per ogni n.** Al punto più difficile c = N/M = 1 i rapporti (k=2, 3, 4) sono:
  D=16: 0.992, 0.982, 0.970;  D=64: 0.9993, 0.9984, 0.9974;  D=256: 0.99995, 0.99987, 0.99980 (scarto ~1/D²).
  A tempi lunghi (n = 8) il rapporto è 1.000 per ogni D.
- **Con rumore non è free Poisson.** Fissando M_eff dalla purezza, i momenti k = 3, 4 escono fino a 4× e 25× sopra la previsione (p = 0.2, L = 2, 3). Lo spettro diventa moltiplicativo, di tipo prodotto (nel limite Markoviano Υ ≈ ⊗_t J_t).

## 6. Note tecniche e problemi aperti

- Il trasferimento S_k in virgola mobile perde precisione per k ≥ 3 a D·n grandi (cancellazioni con Weingarten). Le tabelle k = 3, 4 a grande D usano `exact.py` (razionali; per il rumore solo depolarizzante con p razionale). La versione float è validata con MC a taglie piccole.
- Aperti:
  1. Limite di replica k → 1 (von Neumann, quenched) a L grande: la transizione resta netta? Servono numerici oltre L = 3 (Υ è un RMPS in tempo: MPS/campionamento).
  2. Haar globale su S⊗E è all-to-all: con ambiente locale (brickwork) x* = (pL)* compete con la profondità di scrambling.
  3. Rumore non unitale: correzione genuina a k = 3 da quantificare.
  4. Misure operative (traccia/diamante) analitiche, oltre il proxy HS.
  5. Statistica spettrale del caso rumoroso (log-normale / prodotto?).
