# Journal

Здесь ведётся рабочий журнал проекта: этапы, решения и важные исследовательские заметки.

## 2026-07-20

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
