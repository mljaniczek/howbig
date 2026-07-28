# DiNA Simulation Study

## Project Overview
Simulation study for the DiNA Tutorial manuscript comparing 5 differential network analysis methods (JGL, SSJGL, DiNGO, multiGGM, JRF) under controlled Gaussian graphical model settings.

## Adding a new method

See **`ADDING_A_METHOD.md`** at the project root for the full
checklist. Short version: every new method (including hyperparameter
variants) needs to be wired into ~17 places across 03/15/19/21/22/23
plus method_params() and the report label/color maps. Skip a step
and you'll get either `error: unknown method` (fit dispatch missed)
or silently NA inference columns (scoring dispatch missed).

## Simulation Design (Manuscript Section 5)

### Topology families

Three parallel topology families. Each setting label encodes its
family + dimensions + perturbation:

| Family | Prefix | Topology | m_pa | Within-block density | Mean degree |
|---|---|---|---:|---:|---:|
| **Legacy** | `SF-` `SF+` `SFC` | BA tree | 1 | ~8% | ~1.9 |
| **Current** | `SF2-` `SF2+` `SF2C` | BA scale-free | 2 | ~15.7% | ~3.8 |
| **AR(1)** | `AR1-` `AR1+` `AR1C` | Banded, NON-modular | n/a | (banded, ρ=0.4, bandwidth=1) | ~2 |

The legacy SF- series (BA tree) is retained for cross-comparison; all
*new* simulation runs use the SF2- or AR1- families. Results never
mix on disk because the prefix changes the setting label, and
gold-standard files live at `gold_standards/<setting>/rep_XXXX.rds`.

### Network Structure (SF2- family — current production)
- **p_block = 25 fixed**; B = ceil(p / p_block)
- Within-block topology: Barabási–Albert scale-free (`m_pa = 2`)
- Within-block density: **~15.7%** (47 edges per 25-node block)
- Between-block probability: `0.01 × (5/B)` — keeps the *expected*
  number of between-block edges roughly constant across p
- Edge weights: empirical (from `mxDist.csv`) by default; uniform
  U[-0.4, -0.1] ∪ [0.1, 0.4] available via `wt_method = "uniform"`
  (see `SF2-p250n1000-uniform`)

### AR(1) banded family
- Non-modular: Θ[i,j] = ρ^|i-j| for 1 ≤ |i-j| ≤ bandwidth, 0 otherwise
- ρ = 0.4, bandwidth = 1 (default; `ar1_bandwidth = 2` for "banded
  AR")
- Perturbation operates on a contiguous block of nodes 1..p_block_ar1
  (default 25), mirroring modular perturbation semantics for
  cross-topology comparability

### Precision Matrix Construction
- Row-sum normalization (Danaher 2014 style): divide each
  off-diagonal by `1.5 × sum(|off-diag in row|)`; diagonal set to 1
- AR(1) family bypasses row-sum normalization and uses direct
  diagonal-boost-to-PD on the banded Θ
- ε = 0.1 baseline; PD safety boost adds |λ_min| + 0.01 when needed

### Perturbation Settings
- **Setting A** (Structural / Deletion): zero out within-Block-1 edges
  in g2 (BA family) or zero out off-diagonals among nodes 1..25 (AR1)
- **Setting B** (Both / Deletion + Activation): Setting A + scale
  Block-2 (or nodes 26..50) off-diagonal magnitudes by 1.5×
- **Setting C** (Magnitude / Activation): scale Block-1 (or nodes
  1..25) off-diagonal magnitudes by 1.5×; no structural change

### Ground Truth (FINAL — magnitude, not XOR)
- **Δ_A (PRIMARY)**: magnitude-threshold differential —
  `|Θ_g1 − Θ_g2| > 0.01`. Captures both structural and magnitude
  perturbations. This is what `04_calc_metrics.R` scores against.
- **Δ_A_XOR** (secondary diagnostic): structural XOR of adjacencies.
  Kept for the `diff_XOR_*` columns; misses Setting C entirely
  (returns 0 differential edges since structure is unchanged).
- **ΔP(edge)**: partial correlation difference restricted to union
  of edges; **ΔP(global)**: same across all pairs.

### Canonical sweeps (post-refactor)
| Sweep | Vary | Hold fixed | Setting list (via `get_sweep_settings()`) |
|---|---|---|---|
| sample_size | n | p=250, BA scale-free, Setting A | SF2-p250n{250,500,1000,2000,5000} |
| feature_size | p | n=1000, BA scale-free, Setting A | SF2-p{50,100,200,250,300,500}n1000 |
| perturbation | A / B / C | p=250, n=1000, BA scale-free | SF2-, SF2+, SF2C p250n1000 |
| **topology (NEW)** | tree / scale-free / banded | p=250, n=1000, A | SF-p250n1000, SF2-p250n1000, AR1-p250n1000 |

50–100 replicates per setting in production; bump to 500 reps at one
canonical setting for tight calibration confidence intervals.

### Performance Metrics
- AUPR, AUC for differential ranking (rank-based, threshold-free)
- TPR @ α, Achieved FDR at the BH operating point for calibrated
  inference (FPR ~ 0 in our regime; omitted from main tables)
- Runtime per method (Section: Computational scalability)

## Pipeline Structure
```
00_helper_functions.R      → vBIC(), lambda_plots()
01_generate_gold_standard  → 7 modular functions for block-structured network generation
02_validate_gold_standard  → Structural validation checks
03_estimate.R              → JGL (BIC + CV), SSJGL, DiNGO, multiGGM, JRF wrappers
04_calc_metrics.R          → TPR/FPR computation
05_run_simulation.R        → Parallel simulation driver (future + furrr)
06_summarize_results.R     → Results aggregation + manuscript tables
```

## Key Data Generation Functions (01_generate_gold_standard_networks.R)
- `generate_block_graph()` → single block's igraph
- `assemble_block_network()` → full p×p adjacency with between-block edges
- `assign_edge_weights()` → empirical or uniform weights
- `build_precision_matrix()` → diagonal dominance construction
- `apply_perturbation()` → Setting A/B group 2 creation
- `compute_gold_standard_truths()` → ΔA, ΔP
- `generate_gold_standard()` → orchestrator returning full gold standard object
- `gen_dat()` → sample MVN data (center-only, no scaling)

## Key Tuning Decisions
- BIC tends to over-regularize at p≈n; **CV (5-fold) preferred** for lambda selection
- Data-adaptive lambda range: median(|off-diagonal S_pool|) × [1/50, 2]
- gen_dat() centers only (no scale) to preserve covariance structure

## Pipeline overview (post-refactor, session 89b4bcc)

The estimation pipeline is now split into **fit** and **score** phases
with raw outputs saved per rep. This means you can sweep scoring
options (PIP mask, diff threshold, FDR target, null correction) without
re-fitting anything expensive.

### Storage layout

```
results/raw_fits/
  <method>/
    <setting>/
      _config_<hash>/
        config.rds       -- the method's fit-time hyperparameters
        rep_0001.rds     -- one file per rep; contains $fit (raw output),
        rep_0002.rds        $params, $seed_*, $status, $runtime_sec, etc.
        ...
results/method_results/
  <method>/
    <setting>.rds        -- scored metrics table (one row per rep)
    <setting>.meta.rds   -- which config_hash + scoring options produced
                            the table (for provenance / reproducibility)
results/null_distributions/
  <setting>_<method>.rds -- permutation/bootstrap nulls (15_compute_null.R)
```

The `_config_<hash>` subdirectory name is an 8-char digest of the
method's fit-time hyperparameters (`method_params()` in
`00_helper_functions.R`). Different hyperparameter configurations produce
different subdirectories, so results never silently mix. To change a
method's hyperparameters, **add a new method label** in `method_params()`;
do not reuse an existing one.

### Typical workflow

```bash
# One-time: generate gold standards (resumable, deterministic seeds)
Rscript 05a_generate_gold_standards.R \
  --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"

# Phase 1: fit. Per-rep files; resumable via --target-reps.
# Start small to verify, then grow. Memory-conscious (future_pwalk).
Rscript 05b_run_methods.R --target-reps 5 \
  --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
  --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"
# ...later, resume to target 100 reps:
Rscript 05b_run_methods.R --target-reps 100 \
  --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
  --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"

# Phase 1.5 (optional): compute null distributions for inference columns.
# Uses fixed-GS mode; run AFTER 05a --gold-mode fixed.
Rscript 15_compute_null.R --n-perm 100 \
  --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
  --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"

# Phase 2: score. Fast; rerun whenever a scoring option changes.
Rscript 05d_score_results.R \
  --methods "JGL_BIC,JGL_BIC_GGL,SSJGL,SSJGL_GGL,multiGGM,JRF" \
  --settings "SF-p50n1000,SF-p125n100,SF-p250n1000,SF+p250n1000"

# Phase 3: merge metrics into one table for the reports.
Rscript 05c_merge_results.R

# Knit reports (09, 11, 12).
Rscript -e 'rmarkdown::render("09_simulation_report.Rmd")'
```

### 05d CLI highlights (sweep without refitting)

```bash
# Try a different PIP mask threshold for multiGGM
Rscript 05d_score_results.R --pip-mask-thr 0.10 --methods "multiGGM" ...

# Looser FDR target for the inference columns
Rscript 05d_score_results.R --fdr-target 0.10 --methods "SSJGL" ...

# Raw empirical p-values (no multiple-testing correction)
Rscript 05d_score_results.R --null-correction none ...

# Score a specific hyperparameter config (not the latest)
Rscript 05d_score_results.R --config 3bbfee91 ...
```

### Resumability

`05b_run_methods.R` checks which `rep_XXXX.rds` files already exist in
the target `_config_<hash>` directory and skips them. Start with
`--target-reps 5` to verify the pipeline end-to-end, then come back and
rerun with `--target-reps 100` to fill in the remaining 95 reps. If a
run crashes or is killed, just rerun — nothing is lost beyond any
in-flight rep. Atomic writes (`rep_XXXX.rds.tmp` -> `rename`) prevent
partial files from being treated as complete.

### Inference (pooled null)

`apply_cached_inference()` (in `04_calc_metrics.R`) supports three
multiple-testing corrections via `correction = "pooled" | "bh" | "none"`.
`"pooled"` (the default) concatenates all n_perm x n_pairs null values
into one pool and computes per-edge p-values against the pool; this
gives p-value resolution of `~1 / (n_perm * n_pairs)` — fine enough to
survive BH at FDR 0.05 on ~31k tests with n_perm as low as 30.
`"bh"` is the old per-edge approach (only works at small p).
`"none"` uses raw empirical p < alpha with no multiple-testing
correction. The default when you're not sure is `"pooled"`.

## Post-presentation followup (after April 2026 deadline)

Items intentionally deferred from the April presentation production
run because they would risk destabilizing the run mid-flight or were
out of scope for the time window. The user has months of revision
time after presentation feedback to incorporate any of these.

- [ ] **Switch R's BLAS to Apple's Accelerate framework** on both
  the main machine and the secondary (sweep) machine. R ships with
  a stock reference BLAS that's notably slower than Apple's
  hand-tuned vecLib, and on Apple Silicon the swap commonly cuts
  JGL / SSJGL / multiGGM fit times by 30–60%. The swap is a one-
  time per-machine setup (symlink `libRblas.dylib` and
  `libRlapack.dylib` to Apple's `libBLAS.dylib` /
  `libLAPACK.dylib` under
  `/System/Library/Frameworks/Accelerate.framework/.../vecLib.framework/...`).
  Defer because: changing BLAS mid-production-run risks numerical
  surprises that would be expensive to diagnose under deadline
  pressure. Do this *before* any subsequent large simulation
  campaign and re-run a single rep to verify outputs match the
  pre-swap run within numerical noise.
- [ ] **Limit `--n-workers` to performance-core count** on Apple
  Silicon machines. `parallel::detectCores()` returns
  `P_cores + E_cores`; the E (efficiency) cores run statistical
  workloads ~2-3× slower than P cores, so a batch with mixed-core
  workers is bottlenecked by the E-core stragglers. Use
  `sysctl hw.perflevel0.physicalcpu` to read the P-core count and
  set `--n-workers` to that. (Relevant for the next big run, not
  the current April one.)
- [ ] **GPU acceleration is not currently feasible.** None of the
  methods in this pipeline (JGL, SSJGL, multiGGM, JRF) have GPU
  bindings and at p ≤ 500 the matrix sizes are likely too small
  for GPU offload to overcome transfer overhead. Documented here
  so future-us doesn't waste time investigating Activity Monitor
  showing `%GPU = 0` — that's the correct outcome for this stack.

- [x] **Switch `m_pa = 1` → `m_pa = 2` (DONE 2026-05).** Introduced as
  a parallel SF2- family rather than an in-place change to legacy
  SF-, so the BA-tree results remain valid for cross-comparison.
  See `NOTES_network_density.md` for the analysis and
  `TODOS_FUTURE_SIMULATIONS.md` for what to bundle into the SF2-
  production run.

- [ ] **Investigate multiGGM MCMC config to reduce saving-phase
  memory peak.** Production multiGGM runs at `2 chains × (2000
  burnin + 1000 saved)` and the saving-phase memory peak (~4 GB
  per worker, dominated by `2 × nsave × p²` posterior storage) is
  the binding constraint on the 16 GB machine — only 1–2 multiGGM
  workers fit in green pressure, even though CPU has more
  capacity. Wall time would be unchanged because total iterations
  are the same; what changes is the per-worker memory peak.

  Candidate alternatives to benchmark on a small (≤ 10-rep)
  comparison run, ideally on `SF-p250n1000`:

  1. **1 chain, longer burnin**: `1 × (5000 + 1000) = 6000` iters.
     Saving peak: ~500 MB (half of current). Same wall time.
     Cost: lose Gelman–Rubin $\hat R$ inter-chain convergence
     diagnostic; more vulnerable to undetected local-mode trapping.
  2. **2 chains, fewer samples**: `2 × (2000 + 500) = 5000` iters.
     Saving peak: ~500 MB. Slightly faster wall time. PIP
     resolution drops from 1/2000 to 1/1000 (still adequate).
     Keeps the convergence diagnostic.
  3. **2 chains, longer burnin, fewer samples**: `2 × (4000 + 500)
     = 9000` iters. ~50% longer wall time but better convergence
     assurance + smaller memory peak. Probably the cleanest
     candidate.

  Validation criteria (compare against current 2-chain default):
  - Per-edge PIP correlation > 0.99
  - Differential AUC within ±0.01 (MC noise)
  - FDR-controlled differential selections at α = 0.05 overlap
    > 95% with the 2-chain baseline
  - Within-chain effective sample size > 100 per parameter (for
    the 1-chain candidate)

  If a candidate passes all four, adopt it as the production
  config and document the comparison in `NOTES_bayesian_fdr.md`
  under "MCMC config justification." This unlocks 4-worker
  multiGGM on 16 GB machines (currently capped at 1–2) and gives
  a defensible answer to reviewers asking why the chosen MCMC
  setup is appropriate.

  Defer because: the 2-chain default is the conservative,
  well-diagnosed choice and switching mid-production would risk
  introducing variance from a config change. Investigation is
  methodological work, not deadline-critical.

## Revisit Before Full Production

Deferred items to consider before committing to the big production run:

- [ ] **SSJGL PIP collapse at p=250.** At larger p, SSJGL's
  `problist1` / `problist2` both collapse to near-zero (absolute
  values too small to display at 4 decimals), so PIP-based masking
  or scoring doesn't work at production scale. At p=20 these PIPs
  are well-calibrated; at p=125 they're borderline; at p=250 they
  collapse. Possible fixes: longer v0 ladder (currently 3 steps),
  tuned lambda0/lambda1/lambda2. Investigate before relying on
  SSJGL PIPs in the report. For now, rely on magnitude-based
  scoring (the default in `score_method()`).

- [x] **CLAUDE.md simulation-design block (DONE 2026-05).** Rewritten
  to cover the three topology families (SF- BA-tree legacy, SF2-
  BA scale-free current, AR1- banded), the four canonical sweeps
  (sample / feature / perturbation / topology), the magnitude-based
  Delta_A primary ground truth, and the renamed
  `within_density_random` parameter.

## Open investigation: DiNGO behavior at p=250 and paper replication

Paused mid-stream during session `99625b1..` to pivot to SSJGL PIP
masking. Resume here when ready.

### What we observed
- At `SF-p250n1000` (p=250, n=1000, our own setting), DiNGO's primary
  differential metric is essentially random:
  `g1_AUC ≈ 0.94` (good per-group), `diff_AUC ≈ 0.56`, `diff_AUPR ≈ 0.005`.
- Wall-clock: one DiNGO call is ~19 minutes at p=250.
- Pilot-rho caching (compute BIC-optimal rho once per setting, reuse
  across reps) does NOT produce speedup: Step 1 (rho tuning via
  glasso grid) is only ~1 minute; the remaining ~18 minutes is
  Step 2 (regression of pooled-precision residuals on the group
  indicator) + post-processing. Rho is extremely stable across reps
  at p=250 (CV ≈ 0% on 2 reps).
- A fixed-rho DiNGO call was actually *slower* than a full-tuning
  call in the one-off benchmark (1995 s vs 1144 s). Most likely
  machine-load noise, but could be an unexpected code-path
  difference in `iDINGO::dingo()` — worth a short investigation.

### What we suspect
- Our `DINGO-paper` setting is misnamed: we replicate paper
  dimensions (p=100, n=75, K=2, scale-free, delete-a-block,
  eigenvalue-boost PD) but NOT the paper's data-generating process.
  The DiNGO paper samples data from DiNGO's own decomposition
  model `ρ_ij = P_ij + Q_ij · x`. We sample from two generic MVNs
  whose precisions differ by a block deletion. Methods typically
  look good on their own DGP and worse on foreign DGPs — this
  may explain the apparent p=250 failure independently of scale.
- Our production default is `B = 0` (no bootstrap). The paper
  uses B ≈ 100 to produce standardized Wald-type z-scores. Worth
  a single-rep comparison before drawing conclusions.

### Concrete next steps when resuming
1. **Verify against the paper.** Pull up Ha, Baladandayuthapani, Do
   (2015) Bioinformatics, Section "Simulation Study," and
   cross-check: (a) single vs modular topology, (b) whether data
   is sampled from the decomposition, (c) default B used in paper
   simulations.
2. **B=100 spot check.** Run one rep of `DINGO-paper` and one rep
   of `SF-p250n1000` with `run_dingo(..., B = 100)` and compare
   diff_AUC / diff_AUPR vs B=0. ~10 min extra per call at p=100;
   probably much more at p=250.
3. **Add a true DiNGO-model generator.** Sample data from
   `ρ_ij = P + Q · x` directly. Keep `DINGO-paper` (rename to
   `DINGO-paper-approx`?) alongside as a structural approximation.
4. **Profile where the 19 min at p=250 actually goes.** `Rprof()`
   around a single `iDINGO::dingo()` call. If Step 2 is dominated
   by an embarrassingly parallel per-edge regression, there may be
   an easy parallelism win.
5. **Re-examine the fixed-rho-slower anomaly.** One clean re-run
   with nothing else on the machine to confirm or rule out
   measurement noise.

Diagnostic artifacts kept from the session:
- Pilot rho at p=250 from `pilot_rho("SF-p250n1000")`: **0.0705** (identical across reps 1 and 2).
- Partial log: `/tmp/dingo_p250_benchmark.log` (may be cleaned up).
- Benchmark script (one-shot): `/tmp/dingo_p250_benchmark.R`.

Related commits in session:
- `bf5ee59` — multiGGM PIP-mask (where this thread started).
- `99625b1` — CLAUDE.md production checklist (prior entry).
