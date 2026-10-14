# AGENTS.md

Guidance for AI coding agents working in this repository.

## Repository layout

This is a DrWatson research project (Julia 1.12) with a nested Julia package:

- **Root project (`DifferentiableEvaporation`)**: the DrWatson "research" environment used by `scripts/`, `notebooks/` and plotting. It depends on `EvaporationModel` via a `path = "src/EvaporationModel"` dev entry in the root `Manifest.toml`. Scripts start with `using DrWatson; @quickactivate "DifferentiableEvaporation"` and use `datadir()`, etc. for paths. `src/site_data.jl` (`include(srcdir("site_data.jl"))`) loads site forcing (`site_forcing`, constant piecewise interpolation) and parameters (`site_parameters`, `model_parameters`) from `data/exp_pro`; it lives in the root project, not the package, because it depends on the `datadir()` layout. `src/dummy_src_file.jl` is DrWatson template leftover.
- **`src/EvaporationModel/`**: the actual model, a standalone package with its own `Project.toml`, `test/` env and `docs/` env (CI tests and docs only this package).
- **`scripts/data_download/`, `scripts/data_processing/`**: Python pipeline (conda env from `environment.yml`) that downloads and processes FluxDataKit eddy covariance, soil (SoilGrids, HiHydroSoil, ESDB), vegetation (MODIS VCF), ICOS soil moisture and land-cover translation data. Run from inside each folder (`python main.py`, paths are relative); outputs go to `data/exp_raw` and `data/exp_pro`. Some sources need manual steps or credentials (ESDB request form, GEE project id in `data_download/conf.py`, interactive ICOS login) — see the READMEs there.
- **`notebooks/`**: exploratory Julia scripts / Pluto notebooks (files starting with `### A Pluto.jl notebook ###`) and Quarto derivations (`.qmd`).
- `/data` and `/figures` are git-ignored.

## Commands

Tests (from repo root; mirrors CI `julia-runtest` on `src/EvaporationModel`):

```bash
julia --project=src/EvaporationModel -e 'using Pkg; Pkg.test()'
```

The full suite is slow (AD, sensitivity, allocation checks). To run a subset, activate the test env and include `common.jl` first, since it defines the shared inputs (`p_model`, `toy_model`, `treatments`, `constant_forcings`, ...) that every test file relies on:

```julia
# julia --project=src/EvaporationModel/test   (with EvaporationModel dev'ed into it)
# paste the `using ...` block from test/runtests.jl, then:
include("src/EvaporationModel/test/common.jl")
include("src/EvaporationModel/test/test_processes.jl")
```

Docs (Documenter; doctests also run as part of the test suite via `test_doctests.jl`):

```bash
cd src/EvaporationModel
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

Formatting: JuliaFormatter with `style = "blue"` (`.JuliaFormatter.toml` at root and in the package).

Root project has `ForwardDiff` `nansafe_mode = true` set as a preference in `Project.toml`.

## EvaporationModel architecture

A two-layer soil (force-restore, Noilhan & Planton / Mahfouf style) + canopy interception bucket land surface model, written to be differentiable (ForwardDiff and Enzyme) end to end.

- **State** `u = [w_1, w_2, w_r]`: surface soil moisture, root-zone soil moisture, canopy water storage. The ODE state appends one accumulator per `model.accumulated` diagnostic.
- **Parameters** `p`: a `ComponentArray` with named fields (`h, z_0ms, w_sat, a, p_soil, b, w_res, w_wp, w_fc, C_1sat, C_2ref, C_3, d_1, d_2, z_obs, kB⁻¹, g_d, r_smin, k_ext`), unpacked with `@unpack` in `compute_diagnostics`. Adding a parameter means updating that unpack and every `ComponentArray` that builds `p` (e.g. `test/common.jl`, scripts).
- **Forcings**: a `NamedTuple` of callables of time `t` (`P, T_a, u_a, p_a, VPD_a, SW_in, R_n, LAI`) plus scalars `lon`, `utc_offset`. Time is Unix seconds (`datetime2unix`). In tests, constant forcings use `Returns(x)` rather than closures over non-const globals, because the latter give silently wrong Enzyme derivatives.
- **`model.jl`**: `ProcessBasedModel{FT}` (mutable, `@kwdef`) → `initialize!` builds the RHS closure and `ODEProblem` → `solve!(model; alg=..., kwargs...)` solves and stores in `model.diagnostics` the means over the windows between the `saveat` times of the diagnostics in `model.accumulated` (default `(:λE_tot, :H_tot, :G)`), computed by `window_means` from accumulator states: one extra ODE state per accumulated diagnostic, appended after `[w_1, w_2, w_r]`, kept out of the error control by the per-state tolerances of `state_tolerances` (also stored in `model.prob`). Instantaneous diagnostics: `model.f_diagnostics(u, p, t)`. Optional `yaxarray_output=true` packs everything in a `YAXArray` (off by default, breaks AD). `compute_diagnostics` is the single place where all process functions are wired together; `compute_tendencies!` only forms `du` from its output.
- **Process files** (`canopy.jl`, `evaporation.jl`, `resistances.jl`, `soil.jl`, `soil_fluxes.jl`, `ground_heat_flux.jl`, `sensible_heat.jl`): pure scalar functions. Alternative formulations are selected by dispatch on singleton method types as first argument (e.g. `surface_resistance(JarvisStewart(), ...)`, `ground_heat_flux(SantanelloFriedl03(), ...)`, `soil_evaporation_efficiency(Pielke92(), ...)`). Many aerodynamic quantities come from `Bigleaf.jl` (pinned git rev via `[sources]`).
- **Thresholds (`thresholds.jl`)**: every discontinuity (floors, caps, clamps, steps) takes a `ThresholdTreatment` — `HardThresholds()` (exact) or `KavetskiSmoothing(; s_w, s_r, s_f, s_E)` (Kavetski & Kuczera 2007 smoothing). Use the overloaded `max(t, a, b, s)`, `min(t, ...)`, `clamp(t, ...)`, `step_function(t, ...)` with the matching scale (`moisture_scale`, `storage_scale`, `factor_scale`, `flux_scale`) instead of plain `max`/`min`/`if` so the model stays smoothable. The treatment is stored on the model and threaded through all process calls.
- **Type genericity / AD**: functions must stay generic over number types (Dual, Enzyme shadows) and allocation-free in the RHS (checked by `test_allocations.jl` with AllocCheck). Helpers `value_type` / `of_value_type` in `utils.jl` create constants of the right type; avoid branching on state values (e.g. the canopy condensation solve is always evaluated twice).
- **Solvers**: exported `Tsit5`, `Euler`, `Heun`, `ImplicitEuler`; implicit Jacobians via `AutoForwardDiff()` or `AutoEnzyme(; function_annotation=Enzyme.Duplicated)` (Duplicated needed since the RHS is a closure over forcings). Default `abstol = reltol = 1e-6`.
- **Tests** cover: individual processes, full model, allocations, AD of the RHS (vs FiniteDiff via DifferentiationInterface), AD through the solver, SciMLSensitivity, doctests — each looped over both threshold `treatments`.

`scripts/reference/` holds scripts generating reference runs/gradient checks against synthetic and site forcing (outputs under `datadir("reference")`).
