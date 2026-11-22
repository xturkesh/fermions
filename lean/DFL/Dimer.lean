import DFL.Sector

/-!
# Disorder-free localization depends on the initial state (paper: Section `sec:dimer`)

A product gauge state with the even links polarized along `z` and the odd links in prescribed
`σˣ` eigenstates (labels alternating along the chain), together with density-wave matter, induces
the **random dimer** sector law: the charges come in equal pairs `q_{2k} = q_{2k+1}`, and the
`2^{N/2}` dimer configurations are equally likely. (The resulting one-body problem is the random
dimer model, which has critical energies `± 2h` for `0 < 2h < J`; the transport lower bounds of
Jitomirskaya–Schulz-Baldes–Stolz are infinite-volume statements and are not formalized here.)
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Prescribed label of the odd link `b = 2k+1`: `(k+1) mod 2`. -/
def dimerLabel (b : Fin N) : ZMod 2 := ((b.val / 2 + 1 : ℕ) : ZMod 2)

/-- The link configurations allowed by the dimer gauge state. -/
def dimerLinks (x : Cfg N) : Prop := ∀ b : Fin N, b.val % 2 = 1 → x b = dimerLabel b

instance (x : Cfg N) : Decidable (dimerLinks x) := by unfold dimerLinks; infer_instance

/-- The dimer gauge state in the `σˣ` basis (unnormalized amplitude `1` on allowed links). -/
def dimerGauge (x : Cfg N) : ℂ := if dimerLinks x then 1 else 0

/-- The full initial state: dimer gauge state ⊗ density-wave matter (unnormalized). -/
def dimerInit (N : ℕ) : Cfg N × Cfg N → ℂ :=
  fun p => (if p.1 = fockCfg (univ.filter (fun x : Fin N => Even x.val)) then 1 else 0) *
    dimerGauge p.2

/-- Dimer charge configurations: `q_{2k} = q_{2k+1}` for every pair inside the chain. -/
def IsDimer (q : Cfg N) : Prop :=
  ∀ k : Fin N, k.val % 2 = 0 → ∀ hk : k.val + 1 < N, q k = q ⟨k.val + 1, hk⟩

instance (q : Cfg N) : Decidable (IsDimer q) := by unfold IsDimer; infer_instance

/-- **The dimer sector law**: in the dual frame the dimer initial state is the uniform
superposition of density-wave matter over all dimer charge configurations. -/
theorem dimerInit_dual (hN : Even N) (n q : Cfg N) :
    dimerInit N (dualEquiv (n, q)) =
      if n = fockCfg (univ.filter (fun x : Fin N => Even x.val)) ∧ IsDimer q then 1 else 0 := by
  sorry

/-- The number of dimer charge configurations is `2^{N/2}`. -/
theorem card_isDimer (hN : Even N) :
    (univ.filter (fun q : Cfg N => IsDimer q)).card = 2 ^ (N / 2) := by
  sorry

/-- **Sector averaging over the dimer law**: every gauge-invariant observable evolves as the
uniform average over dimer charge configurations of free-fermion expectations. -/
theorem sector_average_dimer (hN : Even N) (J h : ℝ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ dimerInit N) =
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        expect (sectorBlock O q)
          (evo (Hsec J h q) t *ᵥ Pi.single (fockCfg (univ.filter (fun x : Fin N => Even x.val))) 1)
          := by
  sorry

end DFL
