# Jet RTL10000 validation: handoff notes (2026-10-30)

Branch: `pr_twolayer_mods`. These notes cover the validation sweep after the new two-layer
k-epsilon model (`kEpsilonTwoLayerFunctorMaterial` + blended TKED kernels) went in.
The next session will have 96 cores, so the jet runs below can be launched directly.

## Status of the whole validation sweep

| Case | Status |
|---|---|
| `ercoftac_030_bfs` (5 tests) | **Regolded, not committed, tests NOT rerun after the regold.** 35 CSVs in `gold/` copied from 8000-iteration runs (byte-identical to the outputs). `bfs_input.i` is back to `num_iterations = 8000` and has the two-layer FunctorMaterial; the `validation` spec has the new two-layer `cli_args`. Standard and two-layer moved closer to ERCOFTAC; the Realizable / quadratic / cubic diffs came only from the step-corner cell drifting with iteration count. Figures: `../../ercoftac_030_bfs/k_epsilon_*_old_vs_new.png` (untracked). |
| `ercoftac_032_channel_flow` (395, 590) | **Passed.** It uses the MOOSE k-epsilon kernels, so this only checks the `moose` submodule bump. |
| `single-stream-jet/RTL10000` | **Regolded at 3000 iterations, test passes (300/300 validation points), not committed.** See "Resolution" below. |
| `heat_transfer/flat_plate` | Not run. No k-epsilon; it only checks the MOOSE bump. Needs `-j20`. |
| `vortex_shedding/cylinder` (`strouhal`) | Not run. No k-epsilon; it only checks the MOOSE bump. Needs `-j16`. |

The `moose` submodule pointer was already modified before this session. Check whether it belongs in the PR.

## How to run the jet

From the repo root:

```bash
./run_tests -C validation/free_flow/isothermal/single-stream-jet/RTL10000/validation -j96 --sep-files
```

- `-j` is the slot count. The test has `min_parallel = 96`, so with fewer slots it is skipped ("insufficient slots").
- The positional path is rejected. Use `-C <spec file>`.
- `--sep-files` writes `jet_RTL10000.{tester,validation,runner_run}_out.txt`.
- If the iteration count goes up, add `--max-time <s>`. `max_time` is 1200 s, and 1500 iterations took about 8 min on 96 ranks (rough estimate from file timestamps).
- **Gotcha:** running the input by hand (outside `run_tests`) does **not** regenerate `combined_halfwidths.csv`. `jet_RTL10000_postprocess.py` produces it by calling `compute_halfwidths.py` in the test folder. Run that script yourself after a manual run, or the half-width will be stale.

## Uncommitted changes to `jet_RTL10000.i`

1. New two-layer syntax: `LinearFVTurbulentBlended{Advection,Diffusion}` on TKED with `blending_factor = two_layer_lambda`. `two_layer_lambda` / `two_layer_epsilon` on `TKED_source_sink`, and `two_layer_lambda` / `two_layer_mu_t` on `compute_mu_t`. A new `[FunctorMaterials/two_layer]` block (`kEpsilonTwoLayerFunctorMaterial`, Wolfstein, `Rey_star = 60`, `DeltaRey = 10`). `two_layer_flavor` was removed from the aux kernel.
2. `num_iterations` changed from 1500 to **3000** for the convergence test. **Set it to the final value before committing.** The spec comment and the gold both assume 1500.

## How the test checks results

- **Custom eval** (`jet_RTL10000_postprocess.py`): `combined_halfwidths.csv` must match `gold/` with rtol 3e-2. This is what fails first ("CUSTOM EVAL FAILED").
- **CSVDiff:** `jet_RTL10000_out_centerline_0002.csv`, rel_err 3e-2.
- **Validation** (`jet_RTL10000_vnv.py`): at each point, |exp - current| must lie within (1 ± 0.02)·|exp - gold| ± 1e-3. The band is tied to gold, so improvements fail too.

## Findings

### 1500 iterations (the official test run, 12:04)
- Custom eval failed. Validation: 198 of 300 points out of band (centerline 112: 56 better, 56 worse; half-width 86: 65 better, 21 worse).
- Up to z/D ≈ 15, the new result matches gold to within 0.8%. Beyond z/D ≈ 20 the jet spreads and decays faster: V_c is up to 6.4% lower and r½ up to 6.2% higher.
- RMS vs experiment: centerline 0.130 → 0.151 (worse), half-width 0.231 → 0.217 (better).
- At 1500 the far field is still evolving. Momentum residuals are still falling (~7e-5), pressure climbs back up from 1e-5 to 2e-5 after iteration 1200, and TKED is stuck at a noisy ~1e-6.

### 3000 iterations (manual run, 14:01)
- The 3000-iteration centerline, radial CSVs and `.e` are in the folder. **`combined_halfwidths.csv` and all `*_out.txt` in the folder are still from the 1500 run.** The 1500 centerline CSV was overwritten.
- Max relative change vs gold:

| z/D | r½, 1500 | r½, 3000 | V_c, 3000 |
|---|---|---|---|
| 0–15 | < 0.8% (2% at the first point) | same | < 0.8% |
| 15–20 | 1.2% | 4.7% | 5.2% |
| 20–25 | 4.5% | 9.0% | 8.8% |
| 25–30 | 6.2% | 9.2% | 8.8% |

- RMS vs experiment at 3000: centerline 0.215 (gold 0.130), half-width 0.205 (gold 0.231).
- **Conclusion:** the far field is not converged at 1500 iterations. Between 1500 and 3000 it moves by up to 5%, in the same direction as the old-vs-new change. The gold "converged after 1500 iterations" probably wasn't converged in the far field either.

### Hypotheses
1. **Convergence drift: confirmed** (above). The open question is where the far field settles.
2. **The two-layer blend may switch on outside the jet: unverified.** The only wall is `front` (`walls = 'front'`), so near the wall nothing should change, which matches the near field. But Re_y = sqrt(k)·y/ν may drop below `Rey_star = 60` in the low-k ambient air. That would give λ < 1 and switch the algebraic two-layer ε on at the jet edge, where entrainment controls spreading. The effect would grow downstream, matching the observed trend. To check it, output λ (see step 3 below). There is no `netCDF4` in `/opt/python`, so sample λ with a VectorPostprocessor to CSV rather than reading the Exodus file.

## Next steps

1. **Find the converged answer.** Run through `run_tests` at about 6000 iterations (set `num_iterations` in the input, raise `--max-time`). Check that the residuals plateau and that r½ / V_c stop moving between 3000 and 6000.
2. **Run the old model to the same convergence, for a fair comparison.** The current code rejects the old syntax (`two_layer_flavor` on the aux kernel was removed), so the old model needs an older build. The gold came from commit `ccce2f5` (2026-09-17, "Regold jet_RTL10000 validation data for deviatoric stress terms"). Build that commit in a git worktree (`make -j6`, per the build cap) and run the same iteration count.
3. **Output `two_layer_lambda`.** Add an AuxVariable plus a functor aux, and sample it along a radial line at a few z/D (for example 5, 15 and 25) with a `LineValueSampler`. If λ < 1 away from the `front` wall, that's a model issue (the blend region needs limiting to near-wall cells), not a reason to regold.
4. Only regold once steps 1 and 2 show the converged new model is the intended answer. If the converged centerline is really worse than the old model, raise it in the PR review instead of just regolding. When regolding, also update `num_iterations` and the "1500 iterations" comment in `validation`.
5. Then run `flat_plate` (`-j20`) and `cylinder` (`-j16`), rerun the BFS tests to confirm the regold, and commit.

## Files

- `jet_RTL10000_old_vs_new.png`: gold vs new at 1500 iterations (centerline, half-width, error vs experiment with the validation band, relative change, residual history).
- `jet_RTL10000_1500_vs_3000.png`: gold vs new at 1500 and 3000 iterations, and drift vs gold against the 3% tolerance.
- Analysis scripts, in `/tmp` (may not survive a reboot):
  - `/tmp/jet_compare/jet_compare.py <jet dir> <out dir>`: metrics plus the 6-panel figure. It reads the CSVs and `jet_RTL10000.runner_run_out.txt` from the jet folder.
  - `/tmp/jet_3000/scripts/cmp.py`, `plot3.py`: the 1500-vs-3000 comparison. They read copies of the outputs in `/tmp/jet_3000`.
  - `/tmp/bfs_compare/compare.py`: the BFS old-vs-new figures.
- Both PNGs and this file are untracked.

## Resolution (2026-10-30, afternoon)

- Old model (worktree `/tmp/op_base`, old syntax, on the bumped MOOSE) and new model both run to 6000 iterations.
  New vs old at 6000: < 0.9% everywhere (2.35% at the first half-width point). Old at 6000 vs gold: up to 9.4%.
  So the far-field drift came from the gold being unconverged at 1500 iterations, not from the two-layer model. Hypothesis 2 isn't needed.
- RMS vs experiment: centerline gold 0.130 / old 6000 0.210 / new 6000 0.211; half-width 0.231 / 0.212 / 0.205.
  Converged, both models decay the centerline too fast in the far field. Mention this in the PR.
- Figure: `jet_RTL10000_old_vs_new_6000.png` (script in `/tmp/jet_6000/scripts/compare6000.py`).
- Regold: `num_iterations = 3000` (within 1.1% of 6000), gold CSVs come from a `run_tests` run, and the spec comment is updated.
  Takes 875 to 910 s on 96 ranks, so it fits under `max_time = 1200`. `run_tests` here has no `--max-time` flag.
- Running on teton0019: `ssh teton0019 "bash -lc 'module load use.moose moose-dev-mpich/2026.10.19; export APPTAINER_SHARENS=false; moose-dev-exec ./run_tests ...'"`
  (or `mpiexec -n 96 moose-dev-exec <exe> ...`).
- Still to do: flat_plate (-j20), cylinder (-j16), rerun the BFS tests, check the moose submodule pointer, commit.
