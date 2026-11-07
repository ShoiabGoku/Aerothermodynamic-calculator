# Compressible & Hypersonic Aerodynamics Calculator

A free, single-file, dependency-free web tool for **compressible-flow and aerothermodynamics studies** — usable by
anyone, for any project, in any browser (works offline; hosts free on GitHub Pages). It combines the classic
perfect-gas relations with a high-temperature **real-gas / equilibrium-air** module for hypersonics, and explains
every result in plain language so students and engineers from any branch (civil, mechanical, chemical, physics)
can use it — not only aerospace specialists.

### 🔗 Live app: https://shoiabgoku.github.io/Aerothermodynamic-calculator/

## Interface
- **Light "thermal" design.** One colour scale runs through the whole app, from cold navy to hot crimson, so colour
  always means temperature. Each tool has its own icon and accent colour.
- **Home page with a live flow model.** Air streams past a blunt body, coloured by temperature. Drag the speed from
  Mach 0.2 to 30 and watch the bow shock form above Mach 1, the nose heat up and the shock move closer to the body.
  The real-air/ideal-gas switch shows how much cooler real air stays once its molecules dissociate. The nose
  temperature comes from the same 5-species equilibrium solver as the Aerothermo tab; the shock shape is the Billig
  correlation.
- **"What do you want to find out?" cards.** Each tool is introduced by the question it answers, the fields that use it
  (civil, mechanical, chemical, aerospace) and a one-line description.
- **Plain-language help on every result.** Result labels have a dotted underline. Hover one for a definition, or click
  it to open the **Basics & glossary** drawer (34 searchable terms with everyday examples).
- **Everyday pictures.** The Isentropic and Normal-shock tabs restate the result for 15 °C sea-level air. The Aerothermo
  tab adds a **"How hot is it?" temperature ladder** that places each model's nose temperature next to familiar
  temperatures (boiling water, melting iron and tungsten, the Sun's surface).
- **Wind pressure for structures.** The Atmosphere tab turns a wind speed into the dynamic pressure q = ½ρV² at any
  altitude and at sea level (the basis of IS 875 / ASCE 7 wind loads).
- **Deep links.** Every tab has its own address, e.g. `…/#aero`, `…/#stube`, `…/#atm`, so you can link straight to a tool.

## Tools (tabs)

| Tab | What it does |
|---|---|
| **Home** | Live flow model, a guide to the tools, the three ideas behind compressible flow, and uses outside aerospace. |
| **Isentropic** | Solve from *any* input (M, p₀/p, T₀/T, ρ₀/ρ, A/A* sub/sup, Mach angle μ, Prandtl–Meyer angle ν) → all the others. |
| **Normal shock** | Solve from any of M₁, M₂, p₂/p₁, ρ₂/ρ₁, T₂/T₁, p₀₂/p₀₁; gives the full jump + pitot ratio. |
| **Oblique shock** | Wedge / compression ramp. Input deflection θ **or** wave angle β; weak & strong solutions; θ–β–M plot; attached/detached check. **Perfect gas *or* real-gas / 5-species equilibrium air.** |
| **Cone (Taylor–Maccoll)** | Sharp cone at zero incidence — shock angle, surface Mach, surface Cₚ, by integrating the Taylor–Maccoll ODE. |
| **Prandtl–Meyer** | Expansion fan: M₂, ν, μ, and p/T/ρ ratios through a given turn. |
| **Aerothermo (real gas)** | The hypersonic module: freestream → bow shock → surface, solved in parallel for **ideal · thermally-perfect · 3/4/5-species equilibrium air** — stagnation T₀/p₀, post-shock state, shock stand-off, surface p/T vs wall angle, stagnation heat flux (Fay–Riddell & Sutton–Graves). Includes a **nonequilibrium post-shock relaxation** card: finite-rate Park R1–R5 chemistry (single T, third-body efficiencies, equilibrium-constant backward rates) integrated behind a chemically-frozen shock — T(x) and Yᵢ(x) profiles plus 95% relaxation lengths, directly comparable to a Laminar Finite-Rate Fluent stagnation line — **and the Park two-temperature (T–Tᵥ) model**: vibrationally-frozen jump, Landau–Teller relaxation with Millikan–White + Park-correction times, dissociation at Tₐ=√(T·Tᵥ), overlaid as dashed T and Tᵥ profiles. |
| **Shock tube / tunnel** | Ideal shock tube and reflected-shock tunnel with 9 driver/driven gases: incident and reflected shocks, contact surface, expansion fan, reservoir (region 5) and optional nozzle test section. Includes an x–t wave diagram, property snapshots and the full equation sheet on the page. |
| **Atmosphere** | US Standard Atmosphere 1976 (to 86 km): T, p, ρ, speed of sound, viscosity. Works in reverse (pressure, density or temperature → altitude) and gives the wind pressure on structures. |
| **Trajectory sweep** | Sweep velocity, Mach, or altitude and **plot** stagnation T₀, heat flux and shock stand-off — e.g. along a re-entry trajectory. |

### Global features
- **Configurable γ** for all perfect-gas tools.
- **SI ⇄ US/Imperial** unit toggle for dimensional quantities.
- **CSV export** of the current results table and **Print / Save-as-PDF**.
- Inline SVG plots (θ–β–M diagram, trajectory curves) — no external libraries.

## Use it
Open `index.html` in any browser — no build, no server, no internet. Or use the live GitHub Pages link above.

### 📲 Install as an app (PWA)
The calculator is an installable Progressive Web App: visit the live link and click the **Install app** button in the
top bar (or the install icon in Chrome/Edge's address bar; on Android use "Add to Home screen", on iOS Safari
"Share → Add to Home Screen"). You get a desktop/home-screen icon, a standalone window, and the whole tool keeps
working **fully offline** — the service worker (`sw.js`) caches the app shell and refreshes it automatically whenever
a new version is pushed. Icons are generated by `make_icons.py`; the app metadata lives in `manifest.json`.

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

- `verify3.py` — shock tube (Ms=2 → T₅/T₁=2.500 exactly; He driver → T₅ ≈ 4863 K).
- `verify4.py` — nonequilibrium relaxation (finite-rate integration relaxes from the frozen shock state exactly onto the independently-computed equilibrium state; O₂-only relaxes in ~2 mm vs ~17 mm for 5-species at the reference condition).
- `verify5.py` — two-temperature Park model (vibrationally-frozen jump reproduces the calorically-perfect T₂ ≈ 12.3 kK with Tᵥ = T∞; T and Tᵥ merge in ~1 mm; end state lands on the 1-T equilibrium).

```
python verify.py
python verify2.py
python verify3.py
python verify4.py
```

## Scope / limitations
- Perfect-gas tools assume a calorically perfect gas (single γ).
- The equilibrium real-gas models use a single temperature and neutral 5-species air (no ionization); thermal nonequilibrium is covered only by the two-temperature relaxation card. The equilibrium solution and the finite-rate relaxation card together bound and profile what a single-temperature nonequilibrium CFD (e.g. Fluent Laminar Finite-Rate with the Park mechanism) should produce — they are a verification aid, not a CFD replacement.
- Sutherland viscosity under-predicts μ above ~2000 K.

For verification and education — not a substitute for converged CFD.

## License
MIT — see `LICENSE`.
