# Future Simulation TODOs

Items deferred from the May 2026 SF2- topology rerun, organized by
when they need attention. Items already landed in the SF2- rerun
(m_pa = 2 switch, AR(1) family, within_density rename, magnitude
ground truth, topology sweep) are not listed here — see CLAUDE.md
for the design-block updates.

## Tier A — Bundle with the SF2- production run if compute allows

Things that regenerate gold standards or require a new permutation
null. Cheapest to do them now while gold standards are being made
fresh anyway.

- [ ] **Run `31_topology_comparison.R` and verify output.** Script
  is written but hasn't successfully executed. Produces four
  figures for the topology-sweep panel of the poster
  (`topology_comparison_full`, `_g1_only`, `_diff_only`, `_theta`).
  Fastest path: run once, sanity-check the PDFs, and if the
  visuals need polish adjust colors / boundaries / zoom range
  inline. Only requires one gold-standard rep per family (fixed
  seed = 42 in the script), so it's cheap to iterate.

- [ ] **Topology sweep wired into reports/poster.** Setting list is
  defined via `get_sweep_settings("topology")` in
  `00_helper_functions.R`, but `06_results_tables.R` and
  `30_poster_figures.R` still split only on sample / feature /
  perturbation. Add a fourth sweep block to both that compares
  `SF-p250n1000` vs `SF2-p250n1000` vs `AR1-p250n1000` at
  Setting A.

- [ ] **SSJGL PIP collapse investigation at SF2- densities.**
  SSJGL's `problist1` / `problist2` collapse to ~0 at p ≥ 250
  under the legacy BA-tree topology. With denser SF2- blocks (~47
  edges/block vs 24), the PIP signal-to-noise should be different.
  Concrete action: try a 5-step `v0` ladder
  `c(0.001, 0.005, 0.01, 0.05, 0.1)` instead of the current 3-step.
  Add a new method label `SSJGL-extended` so the comparison is
  clean.

- [ ] **multiGGM MCMC config tuning before scale-up.** The two-chain
  memory peak (~4 GB / worker) is the binding constraint above
  p=300. Validate the `2 chains × (4000 burnin + 500 saved)`
  candidate from CLAUDE.md "Post-presentation followup" on
  `SF2-p250n1000` (10 reps) before running multiGGM on the full
  SF2- sweep. Validation criteria: per-edge PIP correlation > 0.99
  with the 2-chain baseline, diff AUC within ±0.01, FDR-controlled
  selection overlap > 95%.

- [ ] **DiNGO B=100 bootstrap + DGP investigation.** Two specific
  checks from the "Open investigation" block in CLAUDE.md:
  - Spot-check at `SF2-p250n1000` with `run_dingo(..., B = 100)`
    vs B=0 to see if DiNGO's random-differential result evaporates
    at production bootstrap depth.
  - Add a true DiNGO-DGP generator (sample data from
    `ρ_ij = P + Q · x` directly) and compare DiNGO performance on
    its native DGP vs the MVN-from-two-precision DGP. Should
    let us cleanly attribute the apparent p=250 failure to
    DGP-mismatch vs scale.

- [ ] **Calibration-CI anchor setting.** 500 reps at one canonical
  setting (suggest `SF2Cp250n1000` — the magnitude-only
  perturbation that exposed the most calibration breakage on
  legacy SF-) so the headline calibration claims have tight
  confidence intervals. The current 50–100 reps give point
  estimates only.

## Tier B — After SF2- rerun ships

Things that don't need to bundle with the rerun and can wait until
after the new tables/figures are stable.

- [ ] **R BLAS swap to Apple Accelerate** (existing item in
  CLAUDE.md post-presentation followup). 30–60% wall-time
  improvement on JGL / SSJGL / multiGGM fits; one-time per-machine
  setup. Do before any *subsequent* large campaign.

- [ ] **Limit `--n-workers` to performance-core count on Apple
  Silicon.** Existing item — `sysctl hw.perflevel0.physicalcpu`
  gives the P-core count.

- [ ] **Bayesian-FDR direct posterior procedure on PIPs.** Currently
  using strict-PIP threshold (≥ 0.95) as the calibrated FDR proxy
  for scale-multiGGM. The full Newton 2004 direct posterior
  procedure is implemented but sensitive to the differential-PIP
  definition (XOR vs magnitude). Resolve and switch to the direct
  procedure if the magnitude-based PIP definition produces
  well-calibrated FDR.

- [ ] **Investigate fixed-rho-slower DiNGO anomaly.** From CLAUDE.md
  open investigation. A one-off benchmark showed fixed-rho DiNGO
  slower than full-tuning (1995s vs 1144s at p=250). Likely
  measurement noise but worth a clean re-run before the manuscript
  freezes any claims about DiNGO scalability.

## Tier C — New experiments worth adding

Items that would expand the simulation's claim coverage. Lower
priority than getting the SF2- rerun done, but worth scheduling.

- [ ] **AR(1) sample-size and feature-size sweeps.** Currently
  only `AR1-p250n1000` and the +/C variants exist (topology-sweep
  members). To make AR(1) a first-class comparison topology, add
  `AR1-p250n{250,500,2000,5000}` and `AR1-p{50,100,200,300,500}n1000`
  mirroring the SF2- naming.

- [ ] **AR(1) wider bandwidth.** Currently bandwidth=1 (true AR(1)).
  Add `AR1b2-p250n1000` with `ar1_bandwidth = 2` to test how
  methods handle wider-banded precision.

- [ ] **Zero-inflated / non-monotonic transformation setting.**
  Sketched in `09_simulation_report.Rmd` TODO at line 2637. Add
  `SF2-zero-inflated` that post-processes MVN samples to mimic
  scRNA-seq drop-out:
    1. `y = round(exp(z))`           (log-normal counts)
    2. with prob `π(y) = exp(-λ y)`, replace y with 0
  Would test whether JRF (tree-based) degrades gracefully under
  non-monotonic distortion where Gaussian-likelihood methods fail.
  Adds complexity to the ground truth: what does "true edge" mean
  when much of the data is zero?

- [ ] **GGL/FGL `weights = "uniform"` comparison.** `SF2-p250n1000-uniform`
  is in the registry; need to actually run it and compare
  side-by-side with the empirical-weight `SF2-p250n1000` results
  to verify calibration findings hold under uniform weights.

- [ ] **JRF-paper setting as a distinct row in heatmaps.**
  `todos.R:41`. Currently treated as identical to the JGL setting
  and skipped in heatmaps. Visually duplicating it (labeled
  "JRF-paper", with a footnote that the setting is identical to
  JGL) improves heatmap readability.

## Tier D — Method-level investigations (don't need new sims)

These can use existing raw fits or require only re-scoring.

- [ ] **JGL_BIC vs JGL_CV supplement.** `todos.R:21`. Add a
  comparison subsection justifying the BIC default.

- [ ] **SSJGL v0 ladder mathematical justification.** `todos.R:16`.
  Currently 3 steps because more steps produce too-sparse
  estimates; the empirical observation needs a mathematical
  rationale before manuscript freeze.

- [ ] **SSJGL shared-PIP matrix interpretation.** `todos.R:23`.
  Document what `problist_shared` represents (P(edge in at least
  one group)? P(edge in both groups)?) and add to the manuscript
  glossary.

- [ ] **multiGGM hyperparameter exploration.** `todos.R:29-30`.
  Try `v0 = 0.01` (vs default `0.0004`) and `w = 0.5` (vs Shaddox
  default `0.9`) on a small subset to see if PIP collapse at
  p ≥ 250 lifts.

- [ ] **DiNGO partial-correlations → precision matrix back-derivation.**
  `todos.R:17`. If we can back-derive precision from DiNGO's
  partial correlations, we can score DiNGO with the same
  magnitude-based metrics as JGL/SSJGL/multiGGM.

- [ ] **JRF metric calculation methodology.** `todos.R:27`. Look up
  what metrics the JRF paper actually computes and align our
  scoring to it for the JRF-paper replication setting.

## Tier E — Documentation, software, tutorial

These are out of scope for the simulation pipeline but tracked here
so they don't get lost.

- [ ] **Gold-standard example plots for each setting.** `todos.R:36`.
  Plot adjacency / precision / partial correlation / sample data
  for one representative gold standard per setting. Supplementary
  material.

- [ ] **Cross-paper metrics-reported table.** `todos.R:33`. Which
  metrics (AUC / TPR / FPR / FDR / KL / ...) each source paper
  reports, side-by-side with what we report.

- [ ] **Install / startup automation.** `install_startup.R`. Make a
  proper package-loading script.

- [ ] **Symmetry / PD check for Sigma matrices.**
  `01_generate_gold_standard_networks.R:516`. Decide whether to
  validate covariance matrices in addition to precision matrices
  during gold-standard generation.

- [ ] **`02_validate_gold_standard.R` documentation.**
  `todos.R:11`. Explain when validation runs and what happens on
  failure.

- [ ] **Extend Kate Hoff Shutta's SpiderLearner to multi-network.**
  `todos.R:34`. Research collaboration, separate manuscript.

- [ ] **"nbnet" unified package.** `todos.R:24`. Long-term software
  goal — wrap all five methods behind one consistent interface.

## Out of scope for this file

The following are tracked elsewhere and intentionally not listed
here to avoid drift:
- Manuscript glossary expansions → manuscript repository
- Software-paper drafts for `spikeyglass` and `multiGGM` → those
  package repositories
- GPU acceleration → CLAUDE.md (documented as not feasible)
- DiNGO behavior at p=250 → CLAUDE.md "Open investigation" block
- multiGGM saving-phase memory peak → CLAUDE.md
  post-presentation followup
