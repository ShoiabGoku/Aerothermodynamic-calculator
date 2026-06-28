# Compressible & Hypersonic Aerodynamics Calculator

A free, single-file, dependency-free web tool for **compressible-flow and aerothermodynamics studies** — usable by
anyone, for any project, in any browser (works offline; hosts free on GitHub Pages). It combines the classic
perfect-gas relations with a high-temperature **real-gas / equilibrium-air** module for hypersonics.

### 🔗 Live app: https://shoiabgoku.github.io/Aerothermodynamic-calculator/

## Tools (tabs)

| Tab | What it does |
|---|---|
| **Isentropic** | Solve from *any* input (M, p₀/p, T₀/T, ρ₀/ρ, A/A* sub/sup, Mach angle μ, Prandtl–Meyer angle ν) → all the others. |
| **Normal shock** | Solve from any of M₁, M₂, p₂/p₁, ρ₂/ρ₁, T₂/T₁, p₀₂/p₀₁; gives the full jump + pitot ratio. |
| **Oblique shock** | Wedge / compression ramp. Input deflection θ **or** wave angle β; weak & strong solutions; θ–β–M plot; attached/detached check. **Perfect gas *or* real-gas / 5-species equilibrium air.** |
| **Cone (Taylor–Maccoll)** | Sharp cone at zero incidence — shock angle, surface Mach, surface Cₚ, by integrating the Taylor–Maccoll ODE. |
| **Prandtl–Meyer** | Expansion fan: M₂, ν, μ, and p/T/ρ ratios through a given turn. |
| **Aerothermo (real gas)** | The hypersonic module: freestream → bow shock → surface, solved in parallel for **ideal · thermally-perfect · 3/4/5-species equilibrium air** — stagnation T₀/p₀, post-shock state, shock stand-off, surface p/T vs wall angle, stagnation heat flux (Fay–Riddell & Sutton–Graves). |
| **Atmosphere** | US Standard Atmosphere 1976 (to 86 km): T, p, ρ, speed of sound, viscosity. |
| **Trajectory sweep** | Sweep velocity, Mach, or altitude and **plot** stagnation T₀, heat flux and shock stand-off — e.g. along a re-entry trajectory. |

### Global features
- **Configurable γ** for all perfect-gas tools.
- **SI ⇄ US/Imperial** unit toggle for dimensional quantities.
- **CSV export** of the current results table and **Print / Save-as-PDF**.
- Inline SVG plots (θ–β–M diagram, trajectory curves) — no external libraries.

## Use it
Open `index.html` in any browser — no build, no server, no internet. Or use the live GitHub Pages link above.

## Methods & references
- **Perfect-gas relations** (isentropic, normal/oblique/conical shock, Prandtl–Meyer): Anderson, *Modern Compressible Flow*; NACA Report 1135.
- **Conical flow**: Taylor–Maccoll ODE, shot for the shock angle.
- **Real-gas oblique shock**: the shock-normal Mach is resolved into the equilibrium normal-shock solver (tangential velocity preserved); the kinematic relation tan(β−θ)/tanβ = ρ₁/ρ₂ closes it — reducing to θ–β–M for a perfect gas.
- **Equilibrium air**: statistical-thermodynamic dissociation constants Kp(T) from molecular partition functions and bond energies (Vincenti & Kruger; Anderson, *Hypersonic and High-Temperature Gas Dynamics*).
- **Stagnation heat flux**: Fay & Riddell (1958); Sutton & Graves, NASA TR R-376 (1971). **Stand-off**: Billig (1967). **Atmosphere**: US Standard 1976.

## Verification
The physics core is re-derived independently in Python and checked against textbook values:
- `verify.py` — equilibrium-air aerothermo (Sutton–Graves ≈ 4.4 MW/m²; ideal T₀ ≈ 12.7 kK → 5.4 kK with full chemistry; element conservation N/O = 3.762).
- `verify2.py` — perfect-gas tools (isentropic M=2: p₀/p=7.824, A/A*=1.6875, ν=26.38°; oblique M=3/θ=20° → β=37.76°/M₂=1.99; cone M=3/σ=15° → β=25.3°; all matching NACA 1135 / Anderson).

```
python verify.py
python verify2.py
```

## Scope / limitations
- Perfect-gas tools assume a calorically perfect gas (single γ).
- The real-gas / equilibrium module uses a single translational temperature (no two-temperature thermal nonequilibrium, no finite-rate kinetics) and neutral 5-species air (no ionization). It is the equilibrium **bounding** counterpart to nonequilibrium CFD, not a replacement.
- Sutherland viscosity under-predicts μ above ~2000 K.

For verification and education — not a substitute for converged CFD.

## License
MIT — see `LICENSE`.
