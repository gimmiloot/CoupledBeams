# Journal

Здесь ведётся рабочий журнал проекта: этапы, решения и важные исследовательские заметки.

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
