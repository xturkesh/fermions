# FermionsLean: prototipo Lean 4 + Mathlib

Prototipo di formalizzazione di fatti **rigorosi** sui fermioni liberi monitorati, cioè su stati
gaussiani che conservano il numero di particelle (determinanti di Slater). Ogni stato è descritto
dalla sua matrice di correlazione `C j k = ⟨c_j† c_k⟩`, un proiettore hermitiano.

**Stato (22 nov 2026).** Compilato con Lean `v4.37.0-rc1` e Mathlib al commit `a495318`
(21 nov 2026). Zero errori, zero warning, zero `sorry`. Tutti i teoremi elencati sotto dipendono
solo dagli assiomi standard `propext`, `Classical.choice` e `Quot.sound` (verificato con
`#print axioms`).

## Contenuto

### `FermionsLean/Measurement.lean`: misura proiettiva di `n_i`

| Teorema | Significato |
|---|---|
| `idem_update` | Identità algebrica astratta dietro la regola di aggiornamento. |
| `measOcc_mul_self`, `measEmp_mul_self` | Dopo la misura (esito 1 o 0), `C'` è ancora un proiettore: gaussianità e purezza sono preservate. |
| `measOcc_isHermitian`, `measEmp_isHermitian` | Anche l'hermiticità è preservata. |
| `measOcc_apply_self/row/col`, `measEmp_apply_self/row/col` | Il sito misurato ha occupazione definita ed è disaccoppiato dal resto. |
| `trace_measOcc`, `trace_measEmp` | Il numero di particelle `Tr C` si conserva. |
| `born_prob_nonneg`, `born_prob_le_one` | `0 ≤ C_ii ≤ 1`: le probabilità di Born sono probabilità vere. |
| `born_average` | La media di Born di `C'` è la `C` defasata `(1-P_i) C (1-P_i) + P_i C P_i`. |

Il significato di `born_average`: le funzioni a due punti mediate sulle traiettorie seguono
l'evoluzione lineare di defasamento (Lindblad). Una transizione indotta dalla misura può quindi
comparire solo in quantità non lineari nello stato, come entanglement e cumulanti.

Regole di aggiornamento, con `P_i = single i i 1`:

    esito 1 (prob. C_ii):     C' = C - C P_i C / C_ii + P_i
    esito 0 (prob. 1 - C_ii): C' = 1 - [(1-C) - (1-C) P_i (1-C) / (1-C_ii) + P_i]

### `FermionsLean/CutRank.lean`: entanglement attraverso un taglio `A|B`

| Teorema | Significato |
|---|---|
| `block_defect` | `C_AA - C_AA² = C_AB C_BA` per ogni `C` idempotente. |
| `rank_block_defect` | `rank(C_AA - C_AA²) = rank C_AB` per un proiettore hermitiano. |
| `rank_toBlocks₁₂_measure_le` | Una misura di occupazione (qualunque sito, qualunque esito) **non aumenta mai** `rank C_AB`. |
| `rank_toBlocks₁₂_blockDiag_le` | Le porte che agiscono solo dentro `A` o solo dentro `B` non aumentano `rank C_AB`. |

### `FermionsLean/EntropyBound.lean`: entropia e fluttuazioni di carica

| Teorema | Significato |
|---|---|
| `four_log_two_mul_le_binEntropy` | `4 ln2 · p(1-p) ≤ h(p)` su `[0,1]`. Non era in Mathlib. |
| `four_log_two_mul_var_le_entropy` | `S_A ≥ 4 ln2 · Var(N_A)`, espresso sullo spettro `ν_k` di `C_AA`. |
| `entropy_le_log_two_mul_card` | `S_A ≤ ln2 · #{k : ν_k ∉ {0,1}}`. |

### `FermionsLean/MinCut.lean`: bound di "min-cut" gaussiano

| Teorema | Significato |
|---|---|
| `rank_sub_mul_self_eq_card` | Per `M` hermitiana, `rank(M - M²)` conta gli autovalori di `M` diversi da 0 e 1. |
| `entropy_le_log_two_mul_rank` | **`S_A ≤ ln2 · rank C_AB`** per ogni determinante di Slater. |

Insieme a `rank_toBlocks₁₂_measure_le`, questo dice che l'entanglement attraverso il taglio può
crescere solo tramite porte che lo attraversano. Le misure e le porte locali non aumentano mai il
bound `ln2 · rank C_AB`.

## Cosa manca

Elenco in ordine di difficoltà crescente:

1. `ν_k ∈ [0,1]` per gli autovalori di `C_AA`, così che `S_A ≥ 4 ln2 · Tr C_AA(1-C_AA)` diventi un
   enunciato direttamente sulla matrice (circa 1 settimana).
2. Porte a cavallo del taglio: una porta a due siti aumenta `rank C_AB` al più di 2. Un sito misurato
   (disaccoppiato) si può spostare da `A` a `B` senza cambiare il rango. Da qui segue il **min-cut
   spazio-temporale** per circuiti gaussiani ibridi (qualche settimana).
3. La misura di Born sulle traiettorie (`ProbabilityTheory.Kernel.traj`) e un argomento di Peierls.
   Ne segue una **legge d'area uniforme in L** per probabilità di misura proiettiva sopra una soglia
   esplicita (circa 2–4 mesi). Il risultato dà un bound superiore rigoroso su qualunque punto di
   transizione.
4. Il ponte many-body: mostrare che queste regole matriciali coincidono con la meccanica quantistica
   nello spazio di Fock, e che `S_A` è l'entropia di von Neumann (circa 2–5 mesi-persona).
5. **Esistenza o assenza della transizione in 1D: oggi è fuori portata.** Non esiste una
   dimostrazione matematica da formalizzare; si veda il report in `../reports/`.

## Build

```bash
cd lean
lake exe cache get   # scarica la cache di Mathlib (serve cache.mathlib.org)
lake build
```

Nel container cloud in cui è stato sviluppato, il proxy bloccava `cache.mathlib.org` e
`release.lean-lang.org`. La toolchain è stata quindi scaricata dalle release GitHub di Lean e i
moduli di Mathlib necessari sono stati compilati da sorgente (circa un'ora su 4 core). Lake segnala
che i file potrebbero usare il nuovo `module` system; è solo un suggerimento e non influisce sulla
correttezza.
