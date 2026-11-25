# Δt study: SciML backward Euler and Kavetski smoothing

*Updated 2026-10-14.* Smoothing gate of the Terrarium port (`2026-10-01-PLAN_terrarium_port_reference_trajectories.md`, PR 0). Goal: find the largest acceptable Δt per solver and threshold treatment, and decide which smoothings the Terrarium port needs.

Target: backward Euler with 8–16 sub-steps per 30-min window (SUMMA uses 8–32 per hour; Kavetski & Kuczera 2007 recommend fixed-step implicit Euler on a smoothed model).

## Setup

- **Accuracy**, against the truth of the same treatment. Each threshold keeps the time-stepping error well below the error of the observations the model will be calibrated against:
  - RMSE(w₁), RMSE(w₂) < 0.004 m³ m⁻³ at the window bounds: 10% of the 0.04 m³ m⁻³ accuracy target of satellite soil moisture (SMAP, Entekhabi et al. 2010).
  - RMSE(λE_tot) < 5 W m⁻² on the 30-min window means: small next to eddy-covariance errors of about 20–50 W m⁻². These errors combine random error that grows with the flux (Richardson et al. 2006) and the ~20% energy-balance closure gap (Wilson et al. 2002).
  - < 1% relative error in the cumulative latent heat ∫λE_tot dt, from the window means (cumulative ET in energy units). Poorly converged schemes can match states but bias the accumulated water fluxes (Clark & Kavetski 2010).

  Max errors are reported but don't count toward the pass/fail.
- **Truth**: `Kvaerno5()`, abstol = reltol = 1e-10, ForwardDiff Jacobian (`AutoFiniteDiff()` fallback).
- **Solvers**:
  - Fixed step:
    - Δt = 1800/N with N ∈ {1, 2, 4, 8, 16, 32} sub-steps per window (Δt ≥ 56.25 s, finer than SUMMA);
    - Heun;
    - **BE-N** = `ImplicitEuler(; nlsolve=fixed_step_newton(), autodiff=AutoForwardDiff())`, with `fixed_step_newton()` = `NLNewton(; max_iter=30, always_new=true, relax=BackTracking(), check_div=false)` (exported by EvaporationModel); `always_new=true` (a new Jacobian every step) because the default quasi-Newton, which reuses the Jacobian across steps, fails after forcing jumps, and `check_div=false` because the divergence check stops Newton after 2 iterations at the onset of heavy rain; `adaptive=false`, tol 1e-6.
  - Adaptive:
    - `Heun()` (embedded Euler error estimate) and `Tsit5()`, at 1e-6;
    - `ImplicitEuler()` and `KenCarp4()` at 1e-6, with `dtmax=1800`.
- **Treatments**: `HardThresholds()` and the default `KavetskiSmoothing()` (`s_w=1e-3, s_r=0.01, s_f=0.01, s_E=2`). Each run is compared with the truth of its own treatment, i.e. the numerical error against the exact solution of the same (smoothed) model equations, as Clark & Kavetski (2010, §6.1) do with models smoothed after Kavetski & Kuczera (2007). The smoothing bias (smoothed truth against hard truth) is reported separately.
- **Cost**: RHS evaluations per 30-min window (Clark & Kavetski 2010, §6.1.3): `sol.stats.nf`, which counts every RHS call, those of the ForwardDiff Jacobians included (OrdinaryDiffEq.jl issue #2537).
- **Case**: DE-Tha, the year 2018 (a drought year), after a 14-day spin-up that is left out of the metrics (as SUMMA-SUNDIALS evaluates long runs; Spiteri et al., 2024). The forcing windows are the 30-min `saveat` bounds (the default `tstops = nothing`): since PR #33 every evaluation within a window, also at its end, uses that window's forcing, so no `tstops` shifted 1 ms before each window edge are needed. Also reported: RRMSE (RMSE over the mean absolute truth) of λE, H, w₁, w₂ and the total storage ρ_w d₂ w₂ + w_r.

## Steps

### 1. Port `scripts/reference/` to the current API

- `setup.jl`:
  - `site_forcing` returns `(forcings, times, ds)`: use `times.span`, `times.bounds`.
  - Forcings are vectors of one value per window (PR #33): no `tstops` (the window bounds default to `saveat`), so no `tstops_before`.
- Diagnostics:
  - Set `accumulated = (:λE_tot, :H_tot, :G, :E_t, :E_i, :E_s, :Q_s, :K_2)`.
  - Rewrite `to_dataframe` as the states at the bounds plus the window means (it currently uses `saveval`).
  - Take the cumulative ET and runoff + drainage from the accumulators and drop `trapz`.
- `active_thresholds` needs `λE_pc`: either return `λE_pc1` from `compute_diagnostics` or recompute it in the script.
- `dt_study.jl`: add adaptive `Heun()` at 1e-6 to `ADAPTIVE`.
- Add `OrdinaryDiffEqSDIRK` to the root `Project.toml` (the Newton settings come from `fixed_step_newton`).
- Don't commit anything yet: run on the uncommitted working tree. `git_info()` only flags uncommitted changes in `src`, so the report's hash doesn't cover the uncommitted script changes.

### 2. Run the baseline

`julia --project=. scripts/reference/dt_study.jl`, serially (the job's 15 GB of memory doesn't fit many parallel Julia processes). It saves the results (one row per run: solver statistics, the accuracy criterion and the relative RMSE) with `produce_or_load` in `data/sims/reference/dt_study_*.jld2` (rerun with `rerun = true`) and as `data/sims/reference/dt_study.csv`, and plots `plots/reference/dt_study_rrmse.png` and `dt_study_rrmse_cost.png`. Present the results for review before going on to step 3.

### 3. Fix BE-N failures, if any

Re-run the study after each change. Suspects, in this order:

1. **Canopy**:
   - Bound w_r at 0 by limiting E_i through a `ThresholdTreatment` kernel, as `f_1` does for w₁.
   - Give hard f_wet a finite slope at w_r = 0.
2. **G switch**: smooth the switch in R_ns, if failures remain.
3. **Other suspects, only if needed**:
   - `beta_to_r_ss` has no guard against β → 0.
   - The C₂ pole and `w_geq` as w₂ → w_sat.
   - `surface_runoff` as w₂ → w_sat.
   - The Lee–Pielke `min(w_1, w_fc)`.

Goal: BE-N accurate at N ≤ 16, ideally at N = 8.

### 4. Decide and document

- Add the decision (required smoothings, solver, Δt) here and in `data/sims/reference/README.md`.
- Update the Terrarium plan, which still mentions the Vern9 truth, `gradient_check.md` and `generate_terrarium_references.jl`.
- Regenerating the reference CSVs is a follow-up.

## Verification

- `julia --project=src/EvaporationModel -e 'using Pkg; Pkg.test()'` passes after any `src` change, including the allocation tests.
- `dt_study.jl` runs end to end on the uncommitted working tree.
- The smoothing bias stays within the accuracy criterion, and BE-N at N ≤ 16 meets it on all cases, or the report says why not.

## References

- Clark & Kavetski (2010), *WRR* 46, W10510. https://doi.org/10.1029/2009WR008894
- Entekhabi et al. (2010), *Proc. IEEE* 98(5), 704–716. https://doi.org/10.1109/JPROC.2010.2043918
- Kavetski & Kuczera (2007), *WRR* 43, W03411. https://doi.org/10.1029/2006WR005195
- Richardson et al. (2006), *Agric. For. Meteorol.* 136, 1–18. https://doi.org/10.1016/j.agrformet.2006.01.007
- Spiteri et al. (2024), *JAMES*. https://doi.org/10.1029/2024MS004256
- Van Beusekom, Spiteri & Clark (2026), *JAMES*. https://doi.org/10.1029/2025MS005353
- Wilson et al. (2002), *Agric. For. Meteorol.* 113, 223–243. https://doi.org/10.1016/S0168-1923(02)00109-0
