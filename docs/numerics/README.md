# Numerical Policies

This directory contains project-wide numerical workflow policies. These
policies define how calculations are organized, checked, resumed, and
reported; they do not define the governing equations of individual physical
models.

- [FEM-3C numerical robustness and limited verification](nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3c-numerical-robustness-and-dissertation-verification) --
  four completed quarter-period time/mesh controls pass the predeclared guide;
  updated evolving-w model difference 7.76472% on its own scale. Full-period
  1D complete, all8 spatial PARTIAL; both medium full-period 3D jobs complete. Limited straight-rod verification
  complete with qualifications; full-period comparison remains illustrative.
  The [scientific summary](nlsp_straight_rod_3d_fem_verification_summary.md) keeps
  linear/static/dynamic evidence, fixed denominators and limits of the claim.
  No fitted coefficients/loads or universal/experimental validation; energy PARTIAL.

- [FEM-3B nonlinear-correction evolution](nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3b-nonlinear-correction-evolution-and-longer-horizon) --
  read-only old-data diagnostic, saved-state p64/p48 preliminary dynamics and one
  horizon chosen before two medium preload+dynamic jobs. Static initial offsets,
  nonlinear evolution and numerical/output qualifications are kept separate.
  Both .25T1 native jobs complete; measured evolving-w difference7.30261%
  and observed signal resolution are qualified by untested3D mesh/time accuracy,
  full1D spatial PARTIAL and native energy PARTIAL. No physical validation PASS.

- [FEM-3AR controlled short continuation](nlsp_nonlinear_dynamic_3d_fem_pilot.md#fem-3ar-controlled-continuation-after-native-elke-output-path-failure) --
  two corrected medium 3D jobs plus frozen-state p64 references complete .05T1;
  transfer/release/motion pass, native energy PARTIAL. Observed correction
  agreement mostly retains the static initial offset; no dynamic accuracy
  certification or automatically extended study. Historical FEM-3A remains failed.

- [FEM-3A short preload/release pilot](nlsp_nonlinear_dynamic_3d_fem_pilot.md) --
  source/protocol/input and1D saved-state preflight pass; first linear medium
  job aborts with native access violation before any accepted output can be
  verified. BLOCKED_BY_SOLVER, zero further jobs/trajectories/figures; actual
  transfer/release/energy comparison unverified. No model or source change.

- [FEM-2R controlled continuation](nlsp_nonlinear_static_3d_fem_validation.md#fem-2r-controlled-continuation-after-input-serialization-failure) --
  six actual static jobs on saved meshes complete; refined Delta w=-9.16268e-6,
  same sign as1D,3.19353% full-profile correction difference on the common scale.
  Signal exceeds observed mesh/recovery/rounding measures; diagnostic complete
  with qualifications, not universal V0/inertia validation. Historical failure
  stays preserved; no new mesh,1D solve, nonlinear dynamics or automatic FEM-3.

- [FEM-2 nonlinear static comparison](nlsp_nonlinear_static_3d_fem_validation.md) --
  frozen load and 1D static preflight PASS; first medium linear CCX input rejected
  before equilibrium. Overall PARTIAL; no 3D correction or mesh comparison,
  no automatic retry, new mesh or dynamic calculation.

- [FEM-1R: one additional .020 mesh](nlsp_linear_rectangular_3d_fem_validation.md#fem-1r-one-additional-mesh-refinement) --
  same frozen rod/reference,24 complete vectors,all8 identities and preset.1%
  mesh changes accepted; four-grid linear evidence, no nonlinear execution.
- [Full-family rectangular FEM-1](nlsp_linear_rectangular_3d_fem_validation.md) --
  three real C3D10 meshes, four-family shape identification PASS, numerical mesh
  convergence/all-family quantitative comparison PARTIAL; linear only.
- [NLSP local 3D FEM readiness](nlsp_3d_fem_environment_readiness.md) --
  verified Windows Gmsh/CCX paths, historical modal reuse and documented nonlinear
  capabilities; audit only, no FEM calculation or new mesh.
- [Frequency-map computation policy](frequency_map_computation_policy.md) --
  the canonical `frequency-map-v1` contract for ordinary maps, strict audits,
  and rendering from saved data.
- [`scripts/sweep_grid_policy.py`](../../scripts/sweep_grid_policy.py) -- the
  existing helper for primary parameter grids and explicitly requested local
  refinements.

Model-specific equations, root-quality thresholds, and physical assumptions
remain in their own theory, implementation, and validation documents.
