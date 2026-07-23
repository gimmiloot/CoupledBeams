# Journal

Здесь ведётся рабочий журнал проекта: этапы, решения и важные исследовательские заметки.

## 2026-07-22

- Separated future frequency-map spectrum generation from figure rendering in
  a documentation-only policy. The future contract defines `fast_plot` as one
  sequential branch-informed beta path per case/model with sorted roots 1--10,
  root 11 as the mandatory K10 guard, periodic/event-driven global checks, and
  triggered strict recovery; root 12 and `full12_resolved` are not fast-mode
  requirements. `certified_audit` remains the research-grade path for
  counterexamples and independent validation, while `plot_only` must perform
  zero root, matrix, SVD, or cache calculations.

- Recorded that the existing S3_12/S3_14 dimensional-frequency PDFs are valid
  certified outputs. Their runtime reflects spectral certification rather
  than PDF rendering, so presentation changes must reuse saved CSV data. This
  documentation task did not recalculate roots, alter numerical results, or
  modify formulas, matrices, solver settings, FEM, article files, or the
  repository-root README.

- Completed the fixed-manifest research Step 3A targeted lower-envelope
  screen with `epsilon_lower_envelope_step3a_v1`. The runner validates all 28
  rows and their full-precision near/buffer epsilon provenance against the
  corrected `factorized_straight_spectrum_v2` thresholds before invoking the
  unchanged branch-informed EB/Timoshenko spectrum layer. It writes CSV-first
  prefix/mode/control/verification/cost products, six compact plots, separate
  primary and force-recompute verification caches, and an unexecuted paired
  Step-3B proposal. Added 53 targeted tests, including synthetic CLI and
  solver-free plot-only coverage.

- The full run resolved K10/root 11 for 28/28 geometries and 56/56 primary
  model spectra; 40/56 also resolved the optional root 12. All 18
  baseline-control/prefix comparisons passed the corrected factorized oracle.
  Independent verification was triggered for 27 geometries, including every
  provisional/near/quality case and the worst sample from each prefix group;
  all passed root and cluster agreement. `S3_12` confirmed a prefix-5
  violation of `1.73946990918e-2`, and `S3_14` confirmed a prefix-6 violation
  of `5.09348548033e-4`. No case was unresolved or numerically indeterminate,
  so the decision is `counterexample_found`.

- Wrote the 18-row paired near/buffer Step-3B proposal for prefixes 2--10 but
  did not execute it or refine any non-baseline epsilon threshold. The result
  is finite 28-case evidence, not a continuous-domain lower-envelope proof.
  No formula, matrix, determinant, coefficient ordering, shared solver
  default/tolerance, FEM/3D FEM, article workspace, or repository-root README
  was changed.

- Completed research step 2.5b with the versioned
  `branch_informed_continuation_v1` spectrum layer and gateway. Exact beta=0
  axial/bending parent blocks preserve the existing EB and Timoshenko unknown
  orderings. Isolated roots use adaptive projected windows; close roots use
  left/right null-subspace clusters and reduced candidates followed by
  unchanged full-6x6 stationary/SVD verification. Seeds cannot create root
  records directly. The global root-11 guard, triggered strict fallback,
  force-global comparison, local-independent refinement, and primary/force
  cache scopes are reported separately. The companion general helper is now
  `general_complete_svd_v2`, removing its former direct seed-acceptance path;
  no production solver default was changed.

- The full gateway resolved `K10_guard_resolved` for 122/122 audited
  model/geometries and `full12_resolved` for 103/122. R1--R3 at base and
  `epsilon +/- 1e-6`, B07/G01/G02/M02, straight-oracle comparisons, accepted
  clusters, local-independent refinements, root-11 guards, and the requested
  force-global samples all passed. The branch-informed pilot included 21/21
  geometries, with 42/42 model spectra resolved at K10 and no changed
  first-ten roots, `N_true`, or first-failure results in the pilot comparison.
  The decision is `ready_for_targeted_step3`.

- Wrote, but did not execute, the future-only 28-case manifest
  `scripts/analysis/thickness_mismatch/audits/data/eb_epsilon_lower_envelope_step3_cases.csv`.
  It uses only full-precision corrected `epsilon_near_n` and
  `epsilon_buffer_n` values, deduplicates prefixes 4/5 and 9/10, and selects a
  compact set of baseline, small-angle, 45/90-degree, high-mu, signed-eta, and
  mixed probes. No step-3 lower-envelope search, FEM/3D FEM, Gmsh, CalculiX,
  physical-model, determinant, coefficient-ordering, or article change was
  made.

- Implemented research step 2.5 as an offline general-spectrum completeness
  layer around the unchanged coupled EB and Timoshenko `6x6` matrices. The
  primary and independent verification configurations combine determinant
  brackets, shifted and half-step grids, normalized SVD valleys, adaptive
  refinement, continuation and cross-model seed windows, and the straight
  factorized oracle only where `beta=0`, `eta=0`. Every accepted root is
  checked in the row-normalized full matrix. Close-root deduplication requires
  Lambda, self-MAC, and compatible search history; exact nullity and coalesced
  continuation tracks remain distinct. Added an algorithm-versioned cache,
  operation counters, a stable audit entry point, synthetic regressions, and
  an explicit auto-spectrum option while leaving the historical pilot default
  on `legacy`.

- Ran the full first-12 audit for the 21 pilot geometries and the requested
  R1--R3 small-angle stresses. The strict general audit resolved both models
  for 17/21 pilot cases; B07, G01, G02, and M02 retain explicit unresolved
  rows. The straight comparison passed 431/432 oracle rows and recovered 65
  roots absent from raw sign-scan prefixes; G02 EB root 12 remains the sole
  oracle mismatch. The stress audit has 26 failing/unresolved rows and a
  minimum retained pair gap of `9.19871808003e-4`. The corrected auto pilot
  included 20/21 cases, excluded M02 without legacy fallback, and changed no
  first-ten roots or `N_true` values. False-safe geometry counts were unchanged
  across all 53 rule-comparison rows, although retention/loss summaries moved
  after M02 exclusion. The decision is `not_ready_for_step3`. No formula,
  determinant, shared solver/tolerance, FEM, article, or step-3 workflow was
  changed or run.

## 2026-07-21

- Corrected research step 2 after demonstrating that the general 6x6
  determinant sign scan can miss two close simple axial/bending roots inside
  one `Lambda=0.01` interval. The earlier prefix-5 through prefix-10 values are
  superseded but preserved with their cache under the named
  `legacy_pre_factorized_root_fix` directories. The corrected source
  `factorized_straight_spectrum_v2` uses the exact axial family and the exact
  4x4 bending block extracted from the unchanged Timoshenko matrix at
  `beta=0`, `eta=0`; it preserves cross-family multiplicity and keeps the raw
  general scan only as an independent completeness audit. No formula,
  determinant entry, shared root solver/tolerance, FEM workflow, or article
  file changed.

- Recomputed the full `epsilon=0.005..0.060`, `K=10`, first-12 baseline.
  Prefix 1 remains right-censored safe through `0.060`. Corrected conservative
  endpoints for prefixes 2--10 are `0.049140625`, `0.037009766`,
  `0.029705078` for prefixes 4--5, `0.024823242`, `0.021326172`,
  `0.018695312`, and `0.016643555` for prefixes 9--10. All nine first-loss
  brackets pass independent force-recompute verification. Prefixes 2--4 are
  unchanged within tolerance; prefixes 5--10 are classified as corrected due
  to a missing root. Presentation rounding is separate from conservative
  four/five-decimal floors.

- The corrected run has 1004/1004 resolved quality rows, 23592/23592 passing
  factorized-spectrum rows (11796 EB and 11796 Timoshenko), nine passing R1--R3
  plus/minus-epsilon regressions, and 720/720 passing first-12 mu-invariance
  rows for `mu=0,0.3,0.7,0.9`. The independent raw general scan misses 155
  factorized roots over all audit scopes (91 EB and 64 Timoshenko); all are
  confirmed by local full-6x6 SVD refinement. All 2754 required axial records
  pass their block/full-matrix checks, including 27 Timoshenko axial records
  absent from the raw sign scan. First-loss
  semantics remain separate from four safe and four unsafe re-entry events,
  53 family reorder events, and 721 points with late individual passes.
  Research step 3 was not implemented or run.

## 2026-07-20

- Completed the selected 21-geometry `K = 10` epsilon a-priori pilot for the
  safe-spectrum-prefix direction. The manifest-driven runner reused the
  existing EB/Timoshenko root, shape, predictor, MAC, cluster, cache, and local
  thickness helpers; all 21 points completed without root or candidate-boundary
  warnings. The CSV-only analysis compared baseline and fold-calibrated
  `epsilon_0`/`epsilon_max` rules with Rules A--D and Rule A-gap. On the 14-case
  baseline transfer, E0-ref retained 0.7863 of usable EB frequencies with zero
  observed false-safe, versus 0.4017 for Emax-ref; Rule A retained 0.8803. The
  two matched-`epsilon_max` triplets still spanned two and six `N_true` modes,
  and calibrated geometry-only rules produced false-safe cases across the
  broader repeated held-out folds. Reduced EB-rule searches also showed
  substantial 16-to-32 grid sensitivity. The pilot therefore supports further
  cascade testing, not replacement of the EB-based certificate or a universal
  guarantee.

- Adopted the `K = 10` safe-spectrum-prefix certification objective: solve the
  first ten sorted Euler--Bernoulli frequencies, use EB-only modal indicators
  to select a conservative prefix `N_hat`, and compute frequencies
  `N_hat + 1, ..., 10` with Timoshenko. The target is sorted-spectrum
  `N_true`, with homologous-mode MAC and cluster checks retained as quality
  diagnostics.
  Calibration is to maximize retained EB frequencies subject to zero observed
  false-safe on calibration geometries, followed by held-out complete-geometry
  transfer checks. The documentation plan is recorded in
  `docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md`; no new
  computation was implemented by that planning change itself.

- Implemented the first CSV-first `K = 10` safe-prefix postprocessor at
  `scripts/analysis/thickness_mismatch/postprocess/analyze_eb_safe_prefix_certification.py`.
  It uses sorted-spectrum targets, reconstructs only EB-derived indicators,
  calibrates Rules A--D on complete train geometries, evaluates deterministic
  cross-source and leave-one-parameter folds, and records false-safe,
  conservative-loss, overlap, exclusion, predictor-consistency, and primitive
  operation-count diagnostics. Legacy source defaults remain unchanged; K-aware
  fields preserve the exact meaning of existing `first8` columns. Complete
  production K=10 source grids were not run as part of this implementation.
