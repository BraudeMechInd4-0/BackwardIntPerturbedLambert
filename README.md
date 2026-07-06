# Perturbed Lambert Solvers — AMDR Benchmark

This repository contains the codebase accompanying the publication *"Backward Picard-Chebyshev Integration Improving Convergence and Efficiency in Perturbed Lambert Solvers"* by Elad Denenberg and Adham Salih.
_Full citation will be added upon publication._

It is a MATLAB implementation and benchmark suite for eight perturbed Lambert solvers under a
full J2–J6 + atmospheric drag + solar radiation pressure + lunar perturbation force
model. Includes a 75,000-case benchmark, convergence history analysis, and a
T\_switch/N\_switch sensitivity grid. All results and figures are reproducible from a
single pipeline script.

---

## Requirements

### MATLAB
- **MATLAB R2023b or later** (tested on R2023b / R2024a)

### Required Toolboxes
| Toolbox | Reason |
|---------|--------|
| Parallel Computing Toolbox | `parfor` in `run_benchmark_batch` (Stage 2) and `generate_benchmark_testcases` (Stage 1) |
| Aerospace Toolbox | `juliandate()`, used for the epoch in `run_full_pipeline.m`, `generate_benchmark_testcases.m`, `run_convergence_analysis.m`, `run_switch_sensitivity.m` |

> Both `generate_benchmark_testcases` and `run_benchmark_batch` auto-detect the Parallel
> Computing Toolbox via `license('test','Distrib_Computing_Toolbox')` and fall back to
> serial automatically. You can also force serial by passing `use_parallel=false` as the
> last argument to `generate_benchmark_testcases`, or replacing `parfor` with `for` in
> `run_benchmark_batch.m`.

**No Aerospace Toolbox? Use Vallado's `jday.m` instead.** All four call sites set an
epoch via `jdepoch = juliandate(2026, 1, 1, 12, 0, 0);`. Replace each with:
```matlab
[jd, jdfrac] = jday(2026, 1, 1, 12, 0, 0);   % vallado/software/matlab/jday.m
jdepoch = jd + jdfrac;
```
`jday.m` is part of the Vallado submodule (see Installation below), so no extra
toolbox is required once you make this substitution.

---

## Installation

This repo uses two git submodules — a Vallado astrodynamics fork (pre-patched, see
"Notes on the Vallado submodule" below) and matlab2tikz. Clone without `--recursive`
(the submodules are populated by the installer instead, to avoid an unnecessarily
large download — see Option B2's warning below for why):

```bash
git clone https://github.com/<your-org>/<repo>.git
cd <repo>
```

### Option A — Automatic (recommended)

```bash
./install.sh        # macOS/Linux/WSL/git-bash
.\install.ps1        # Windows PowerShell
```

This fetches `matlab2tikz` in full and the Vallado fork as a **lean, sparse checkout**
(only `software/matlab/`, ~2 MB — see below for why this matters), then writes a
`startup.m` that adds both to the MATLAB path.

### Option B — Manual

**B1 — Lean (matches what `install.sh` does):**

```bash
git submodule update --init -- matlab2tikz
git submodule update --init --no-checkout --filter=blob:none -- vallado
cd vallado
git sparse-checkout init --cone
git sparse-checkout set software/matlab
git checkout $(git -C .. ls-tree HEAD -- vallado | awk '{print $3}')
cd ..
```

**B2 — Simple, but large:**

```bash
git submodule update --init --recursive
```

⚠️ **This downloads the entire Vallado fork — approximately 237 MB** — because the
fork also carries `stk/` (132 MB of STK scenario files), `software/cpp` +
`software/python` (~58 MB combined), and `datalib/` (38 MB), none of which this
pipeline uses. Only `software/matlab/` (~2 MB) is actually needed. Use Option A or B1
unless you specifically want the rest of the Vallado library.

| Submodule | Path | Source |
|-----------|------|--------|
| Vallado astrodynamics library (pre-patched fork) | `vallado/` | https://github.com/BraudeMechInd4-0/fundamentals-of-astrodynamics |
| matlab2tikz | `matlab2tikz/` | https://github.com/matlab2tikz/matlab2tikz |

No manual file-patching step is needed either way — the Vallado fork already contains
the `lambhodograph.m`/`lambertb.m` fixes (see "Notes on the Vallado submodule" below)
and the vectorized `sun_v.m`/`moon_v.m` ephemerides the force model relies on.

---

## Running the full pipeline

```matlab
run_full_pipeline
```

Edit the flags at the top of `run_full_pipeline.m` to skip stages you have already
completed:

```matlab
run_generation    = true;   % Stage 1 — generate 75k test cases   (~5 min)
run_benchmark     = true;   % Stage 2 — run 8-solver benchmark     (~8-10 h, 8-core)
run_analysis      = true;   % Stage 3 — summary statistics table
run_select_cases  = true;   % Stage 4 — select convergence cases
run_conv_analysis = true;   % Stage 5 — convergence history (debug)
run_sensitivity   = true;   % Stage 6 — T/N_switch sensitivity grid
run_figures       = true;   % Stage 7 — paper figures (TikZ + PNG)
```

All outputs land under `benchmark_results/`.

---

## Pipeline stages

| Stage | Script called | Output |
|-------|--------------|--------|
| 1 | `generate_benchmark_testcases(75000, testcases_file)` | `testcases_75k.mat` |
| 2 | `run_benchmark_batch(...)` | `results_75k_consolidated.mat` |
| 3 | `analyze_benchmark_results(results, testcases, output_dir)` | summary `.tex` / `.csv` |
| 4 | `select_convergence_cases(testcases, cases_file, N_per_cell)` | `convergence_cases.mat` |
| 5 | `run_convergence_analysis(testcases, selected_indices, history_file, Cr, jdepoch)` | `convergence_history.mat` |
| 6 | `run_switch_sensitivity(testcases, output_dir, n_sample)` | `switch_sensitivity.mat` |
| 7a | `visualize_benchmark_paper(results, testcases, figures_dir, Tolr [, filter])` | `paper_figures/*.{tikz,png}` |
| 7b | `visualize_convergence_paper(convergence_history, conv_figures_dir)` | `convergence_analysis/*.{tikz,png}` |
| 7c | `visualize_switch_sensitivity(sensitivity, figures_dir [, filter])` | heatmaps + LaTeX tables |
| 7d | `visualize_direction_analysis(results, figures_dir [, filter])` | direction stats |

All stage functions also accept file paths as first arguments for standalone use
(e.g. `select_convergence_cases('benchmark_results/testcases_75k.mat', ...)`).

---

## Solvers

### 8-solver benchmark set

| Tag | File | Phase 1 | Phase 2 | Direction |
|-----|------|---------|---------|-----------|
| **T** | `prtlambertT.m` | Thompson | — | forward |
| **MPS** | `prtlambertMPS.m` | MPS | — | forward |
| **B** | `prtlambertB.m` | Broyden (FD J₀) | — | forward |
| **TB** | `prtlambertTB.m` | Thompson | Broyden | forward |
| **MPSB** | `prtlambertMPSB.m` | MPS | Broyden | forward |
| **TRB** | `prtlambertTRB.m` | Thompson | Broyden | smart preamble |
| **MPSRB** | `prtlambertMPSRB.m` | MPS | Broyden | smart preamble |
| **BRB** | `prtlambertBRB.m` | Broyden (FD J₀) | Broyden | smart preamble |

**Smart preamble** (TRB, MPSRB, BRB): propagate forward from r₁ and backward from r₂;
choose the direction whose residual is smaller. Phase 1 starts from the better end.

**FD J₀** (B, BRB): the initial Broyden Jacobian is computed via three finite-difference
perturbations before the iteration loop, matching the actual d(residual)/d(velocity)
sensitivity. Cost: 3 extra integrations at Phase 1 start.

### Extended 11-solver set (convergence analysis only)

Adds `prtlambertTR.m`, `prtlambertMPSR.m`, `prtlambertBR.m` — backward-only variants
(no Broyden phase, always use `odeMPCIrev`) — and `prtlambertBRB.m` from the benchmark
set. These four are used in `run_convergence_analysis` to isolate the contribution of
the backward propagation direction.

### Solver signatures

**Forward solvers with explicit ODE solver (T, MPS, B, MPSB):**
```matlab
[v1t, v2t, stats, history] = prtlambertXxx( ...
    forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, options)
```

**Forward solvers with Learning_rate (T, TB):**
```matlab
[v1t, v2t, stats, history] = prtlambertXxx( ...
    forcemodel, odesolver, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
```

**Reverse/hybrid solvers — no odesolver arg, use `odeMPCIrev` internally (TR, TRB, MPSR, MPSRB):**
```matlab
[v1t, v2t, stats, history] = prtlambertXxx( ...
    forcemodel, r1, r2, v1, dm, de, nrev, dtsec, Learning_rate, options)
```

**Broyden reverse/hybrid solvers — no odesolver, no Learning_rate (BR, BRB):**
```matlab
[v1t, v2t, stats, history] = prtlambertXxx( ...
    forcemodel, r1, r2, v1, dm, de, nrev, dtsec, options)
```

---

## Test-case generation

`generate_benchmark_testcases` draws from a stratified 5 × 5 grid:

- **5 revolution counts** (nrev = 0 … 4), 20% each
- **5 orbital regimes** per nrev, with these fractions:

| Regime | Fraction | Perigee altitude | Apogee altitude |
|--------|----------|-----------------|-----------------|
| LEO | 15% | 350 – 2 000 km | ≤ 2 000 km |
| MEO | 15% | 2 000 – 20 200 km | ≤ 20 200 km |
| GEO | 15% | 35 286 – 36 286 km | 35 286 – 36 286 km |
| HEO | 15% | 350 – 1 000 km | 20 000 – 42 000 km |
| hard | 40% | 200 – 350 km (sub-LEO) | up to 42 000 km |

Each test-case struct has fields:
`id, r1, v1, r2, v2, m, A, tm, dm, de, nrev, is_hard, regime`

---

## Key functions reference

### ODE / integration

| Function | Description |
|----------|-------------|
| `odeMPCI.m` | Chebyshev-Picard forward integrator (Bai & Junkins 2011) |
| `odeMPCIrev.m` | Chebyshev-Picard backward integrator (reverse time) |

### Force models

| Function | Model |
|----------|-------|
| `orbit_eq.m` | Two-body |
| `orbit_eq_J2_drag.m` | J2 + atmospheric drag |
| `orbit_eq_J6_drag.m` | J2–J6 + atmospheric drag (zonal terms computed inline) |
| `orbit_eq_J6_drag_SRP_moon.m` | J2–J6 + drag + SRP + lunar perturbation (used in benchmark) |
| `drag_accel.m` | Drag acceleration helper (embedded exponential density model) |
| `sun_v.m` | Vectorized Sun position (from JD); part of the `vallado` submodule (fork addition) |
| `moon_v.m` | Vectorized Moon position (from JD); part of the `vallado` submodule (fork addition) |

### Orbital utilities

| Function | Description |
|----------|-------------|
| `kep_elements.m` | r, v → Keplerian elements |
| `hitearth.m` | Earth-impact check during transfer |
| `posnvelos.m` | Position/velocity sampling helper (test-case generation) |

### Benchmark & analysis

| Function | Description |
|----------|-------------|
| `generate_benchmark_testcases.m` | Stratified test-case generator |
| `run_benchmark_batch.m` | Parallel 8-solver benchmark runner with checkpointing |
| `analyze_benchmark_results.m` | Statistical summary tables |
| `select_convergence_cases.m` | Stratified case selection for convergence analysis |
| `run_convergence_analysis.m` | 11-solver debug run with error history |
| `run_switch_sensitivity.m` | T\_switch × N\_switch grid sweep (TRB, MPSRB, BRB) |
| `visualize_benchmark_paper.m` | Publication benchmark figures |
| `visualize_convergence_paper.m` | Per-solver convergence graphs |
| `visualize_switch_sensitivity.m` | Sensitivity heatmaps and LaTeX tables |
| `visualize_direction_analysis.m` | Backward-direction selection analysis |

---

## Benchmark parameters (defaults)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Tolr` | 1 × 10⁻⁸ km | Position residual tolerance |
| `maxIter` | 200 | Max solver iterations |
| `T_switch` | `Tolr × 1000` | Phase-switch residual threshold |
| `N_switch` | `floor(maxIter × 0.05)` = 10 | Phase-switch iteration threshold |
| `Cd` | 2.2 | Drag coefficient |
| `Cr` | 1.2 | Reflectivity coefficient (SRP) |
| `jdepoch` | 2026-01-01 12:00 UTC | Atmosphere / SRP reference epoch |
| `N_per_cell` | 8 | Convergence cases per nrev × regime cell (200 total) |
| `switch_n_sample` | 8 000 | Cases sampled for sensitivity grid |

---

## Output directory layout

```
benchmark_results/
├── testcases_75k.mat
├── results_75k_consolidated.mat
├── convergence_cases.mat
├── convergence_history.mat
├── switch_sensitivity.mat
├── paper_figures/
│   └── *.{tikz,png}
└── convergence_analysis/
    └── *.{tikz,png}
```

---

## Third-party libraries

| Directory | Library | Mechanism | License |
|-----------|---------|-----------|---------|
| `vallado/` | Vallado *fundamentals-of-astrodynamics* (pre-patched fork) | git submodule | see upstream repo |
| `matlab2tikz/` | matlab2tikz (Schlömer et al.) | git submodule | BSD 2-Clause |

---

## Notes on the Vallado submodule

This repository points the `vallado` submodule at a **fork**
(`https://github.com/BraudeMechInd4-0/fundamentals-of-astrodynamics`) that carries two
fixes on top of upstream CelesTrak, applied to `vallado/software/matlab/`:

- **`lambhodograph.m`** — the `constastro;` workspace-injection call is replaced with
  inline `mu1`/`twopi` constants (and `mu` renamed to `mu1` throughout). `constastro`
  fails silently inside `parfor` workers and `mu` collides with MATLAB's built-in
  `mu()` function, returning `v1t = [0 0 0]`.
- **`lambertb.m`** — the `constmath`/`constastro` calls are commented out for the same
  reason.
- **`sun_v.m`, `moon_v.m`** (new files) — vectorized adaptations of Vallado's
  `sun.m`/`moon.m`, used by `orbit_eq_J6_drag_SRP_moon.m` for batch ephemeris
  evaluation.

If you point the submodule at a different commit or a different Vallado fork/mirror
(`git submodule update --remote vallado`, or repointing the URL), you must **re-apply**
these fixes yourself, otherwise the benchmark will silently return zero velocity
vectors in parallel runs.
