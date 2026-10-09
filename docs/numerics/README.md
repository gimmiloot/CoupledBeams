# Numerical Policies

This directory contains project-wide numerical workflow policies. These
policies define how calculations are organized, checked, resumed, and
reported; they do not define the governing equations of individual physical
models.

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
