# Aerothermodynamics Calculator — Freestream → Shock → Surface

A single-file, dependency-free web calculator for compressible & hypersonic aerothermodynamics.
Enter **freestream conditions** and it returns stagnation temperature/pressure, post-shock state,
shock stand-off distance, post-shock Mach number, surface pressure & static temperature versus wall
curvature, and stagnation-point heat flux — solved **in parallel for five gas models** so you can
read off the physics of each fidelity level side by side:

| Model | Physics |
|---|---|
| **Ideal gas** | Calorically perfect air, γ = 1.4, constant Cp. Closed-form isentropic + normal shock. |
| **Real gas** | Thermally-perfect air with vibrational excitation (variable Cp(T)); conservation-based shock. |
| **3-species** | Equilibrium air, O₂ ⇌ 2O (O₂, O, N₂). |
| **4-species** | Adds N₂ ⇌ 2N (+ N). |
| **5-species** | Adds NO ⇌ N + O (+ NO) — full neutral equilibrium air. |

Valid across **subsonic, transonic, supersonic and hypersonic** regimes (below M = 1 the bow shock
vanishes and only isentropic stagnation conditions are shown).

## Use it

Just open `index.html` in any browser — no build, no server, no internet. Or host it free on
**GitHub Pages** (Settings → Pages → deploy from branch → root) and use the published URL.

The page loads pre-filled with the M∞ ≈ 15, 50 km blunt-body case and computes immediately.

### Inputs
- Velocity by **Mach** or **m/s**
- Atmosphere by **altitude** (US Standard 1976, auto-fills T, p, ρ) or **manual T, p**
- Nose/leading-edge radius, body type (sphere or 2-D cylinder), wall temperature, wall catalycity

### Outputs
- Freestream summary (a∞, ρ∞, Re/m, total enthalpy, dynamic pressure, regime classification)
- Model-comparison table: T₀, T₂, p₂, p₀, ρ₂/ρ∞, M₂, stand-off Δ/Rₙ and Δ (mm), heat flux, species fractions, dissociation %
- Per-model equilibrium composition (mole fractions of O₂, N₂, NO, O, N)
- Surface distribution: Cₚ, pressure and **static temperature vs wall angle** (the "static temperature at different wall curvatures")

## Methods & references
- **Perfect gas / normal shock:** Anderson, *Modern Compressible Flow*.
- **Equilibrium air:** statistical-thermodynamic dissociation constants Kp(T) from molecular constants
  (characteristic vibration/rotation temperatures, ground electronic degeneracies) and bond energies
  derived from formation enthalpies, with element- and Dalton-pressure constraints
  (Vincenti & Kruger, *Introduction to Physical Gas Dynamics*; Anderson, *Hypersonic and
  High-Temperature Gas Dynamics*, Ch. 11).
- **Stagnation heat flux:** Fay & Riddell (1958); Sutton & Graves, NASA TR R-376 (1971).
- **Shock stand-off:** Billig (1967) correlation + density-ratio thin-shock-layer relation.
- **Surface pressure:** modified Newtonian, Cₚ = Cₚ,max·cos²θ.

## Validation
`verify.py` re-runs the physics core in Python against the M∞ ≈ 15 / 50 km case. It reproduces
Sutton–Graves ≈ 4.4 MW/m², ideal T₀ ≈ 12.7 kK collapsing to ≈ 5.4 kK with full chemistry, the
monotone ρ₂/ρ∞ rise across the model ladder, and exact element conservation (N/O = 3.762).

```
python verify.py
```

## Scope / limitations
- Single translational temperature — no two-temperature (T–T_v) thermal nonequilibrium, and no
  finite-rate kinetics. The equilibrium result is the **bounding** counterpart to nonequilibrium CFD,
  not a replacement for it.
- Chemistry is neutral, 5-species air; ionization and electronic excitation (significant above
  ~8000–9000 K) are not modelled.
- Sutherland viscosity under-predicts μ above ~2000 K.

For verification and teaching — not a substitute for converged CFD.

## License
MIT — see `LICENSE`.
