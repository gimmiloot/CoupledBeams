# Numerical Policies

This directory contains shared numerical policies and model-specific numerical
reports. Their presence together does not make their physical models interchangeable.
For other directions use the [research index](../research_index.md#current-entry-points).

## NLSP

- **Модель:** nonlinear spatial seven-field rods — шесть движений сечения
  и независимое сокращение толщины c; принятое действие V0, кубические
  балансы и переменная инерция. Два плеча соединены идеальным жёстким
  common-DOF узлом с c1=c2, без переноса NSL/KV-податливости.
- **Каноническая теория:** [семиполевая модель](../theory/weakly_nonlinear_spatial_rod.md).
  [Отчёт руководителю](../reports/nlsp_beta15_parametric_instability/README.md)
  раскрывает механику и область результата; исходник/PDF доступен в Git.
- **Реализация:** [энергия и состояния](../../scripts/lib/weakly_nonlinear_spatial_rod.py),
  [пространственная динамика](../../scripts/lib/weakly_nonlinear_spatial_dynamics.py),
  [двухплечевая сборка](../../scripts/lib/weakly_nonlinear_coupled_spatial_dynamics.py).
- **Последнее состояние:** [NLSP-D32](../memory/decisions.md#nlsp-d32) /
  [NLSP-K33](../memory/knowledge.md#nlsp-k33),
  [завершающие проверки](nlsp_beta15_pre_supervisor_checks.md):
  OUT-рост типа II при β15/μ.1/A.005 подтверждён в p12 и воспроизведён в p14.
  Обычное OUT2-отклонение возбуждает рост, но физический прогиб за5T мал.
  Причина PL neutral residual объяснена; строгий gate остаётся PARTIAL.
- **Устранено адресно:** ошибки ввода FEM-2 и native output FEM-3A получили
  продолжения FEM-2R/FEM-3AR. [Сводка прямого стержня](nlsp_straight_rod_3d_fem_verification_summary.md)
  сохраняет историю и границы; исходные отказы не переписываются.
- **Остаётся:** полная пространственная сходимость, подготовленные strict
  критерии, [1D/3D нелинейное сравнение](nlsp_spatial_nonlinear_3d_fem_verification.md)
  и нативная энергия FEM квалифицированы; c_eff — proxy.
  Не закрыты DEFERRED_TWO_PARAMETER_ACTION_EXPANSION и физическая
  обоснованность идеального углового узла/c1=c2. Результат одной ветви
  не доказывает общей устойчивости, динамики типа I или насыщения.
- **Остановка:** обсуждение с руководителем. Новые p/временные уточнения,
  ODE/FEM/Floquet, параметры и физическая валидация не разрешены автоматически.
- **Воспроизведение:** параметры, команды и source hashes — в
  [техническом отчёте](nlsp_beta15_pre_supervisor_checks.md);
  массивы — локально в `results/nlsp_beta15_pre_supervisor_checks/d3f6375a27dba883/`.
  Читать научный итог можно без этих массивов; команды требуют своих входов.

## Historical stage navigation

The entries below describe their original bounded stages. Later continuations
do not retroactively change PARTIAL/FAIL or NOT_RUN records. Use the passport
above for the current stop; these historical links are optional evidence paths.

- [Beta15 final checks before supervisor discussion](nlsp_beta15_pre_supervisor_checks.md):
  neutral residual cause identified with the unchanged PL gate still PARTIAL;
  displacement-only OUT2 projection and one short p14 tangent control through5T.
  [Standalone Russian report and editable source](../reports/nlsp_beta15_parametric_instability/README.md).

- [Beta15 planar 1:3 internal resonance](nlsp_beta15_planar_1to3_internal_resonance.md):
  canonical direct/indirect coupling and controlled algebraic type-II candidates;
  type-I modulation remains separate, no OUT/Floquet or new ODE.

- [Beta15 nonlinear planar resonance](nlsp_beta15_nonlinear_planar_resonance.md):
  formal PL correction encounters a resolved near 1:3 planar denominator;
  isolated-orbit gate fails, M1/M2 are not available. No ODE or new geometry.

- [Out-of-plane instability feasibility](nlsp_out_of_plane_instability_feasibility.md):
  action/Jacobian agreement and canonical formal growth channel confirmed;
  beta15 actual direct pair remains oscillatory, true-orbit instability unproved.


- [Five-angle resonant coupling screening](nlsp_beta_resonant_coupling_screening.md):
  fixed mu=.1, no ODE; no direct pair-growth candidate, planar correspondence
  at beta15 unresolved. Mean terms are decisive for the closest self channel.


- [Resonant spatial coupling screening](nlsp_resonant_spatial_coupling_screening.md):
  full-Jacobian modal/harmonic diagnostic without new ODEs; C interactions
  resolved, no direct growth candidate, mu=0 correspondence unresolved.

- [Construction C modal sensitivity](nlsp_construction_C_modal_sensitivity.md):
  one modal amplitude, displacement-only OUT seeds, fixed-energy restricted
  gain and bounded spatial/time controls; no new FEM or parameter map.

- [Authorized two-arm sensitivity continuation](nlsp_coupled_planar_sensitivity_resume.md):
  mixed/sign consistency, four-input reference-energy gain, independent full
  nonlinear responses and bounded numerical qualifications. Historical results
  remain separate; no parameter map or stability certificate.

This directory contains project-wide numerical workflow policies. These
policies define how calculations are organized, checked, resumed, and
reported; they do not define the governing equations of individual physical
models.

- [Addressed torsional tangent diagnostic](nlsp_coupled_torsion_tangent_diagnostic.md) --
  historical failure preserved; independent short temporal controls identify
  the full-solver numerical floor and recover observed delta-squared response.

- [Two-arm planar sensitivity pilot](nlsp_coupled_planar_sensitivity_pilot.md) --
  endpoint lifting, common seven-DOF rigid joint, tangent/full nonlinear
  checks and restricted reference-energy initial-state gain. Two frozen
  geometries and bounded 2Tref programme; no new FEM or parameter map.

- [Seven-field spatial nonlinear verification](nlsp_spatial_nonlinear_3d_fem_verification.md) --
  six 1D and four actual medium/fine 3D trajectories to .25T1; overall
  NUMERICAL_PARTIAL. Full w/v agree in scale, but v correction differs
  substantially. Weak-form/planar checks, all14, mesh/recovery, dynamic twist
  and native-energy qualifications stay separate; no automatic follow-up.

- [Saved spatial-profile audit](nlsp_spatial_profile_audit.md) --
  limited postprocessing of u, theta, c and 3D effective contraction after FEM-3C.
  Raw/interpolated recovery, 21/41/81 slabs, actual medium/fine states,
  independent volume-strain averages and saved p48/p64 profiles are kept
  distinct; no new physical solve or altered historical status.

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

- [Construction C to one OUT period and cost audit](nlsp_construction_C_longer_sensitivity.md) — saved-state p12 continuation,Gdisp3=1.07027; conditionalp16 NOT_RUN bybudget,newinterval numericalPARTIAL; immutable output failure/authorizedrepeat,no maps.

- [β15 periodic orbit and spatial stability](nlsp_beta15_periodic_orbit_spatial_stability.md) — full p12 HB9/11, resolved OUT candidate atA=.005, sigma=.0034114101; three HB gates fail, no ODE/FEM or automatic dynamics.


- [β15 dynamic OUT growth validation](nlsp_beta15_dynamic_out_growth_validation.md) —
  A=.005 Type-II Ψ=0, p12: direct tangent and full seven-field dynamics confirm
  saved Floquet growth through5T; gain1.279292, max OUT deflection .735%A.
  Half-scale error-order qualification remains numerical-floor limited; no physical validation.


- [β15 OUT growth, p14 spatial check](nlsp_beta15_out_growth_p14_convergence.md) —
  same A=.005 Type-II Ψ=0 branch; σ=.0034113155, p12/p14 difference .00277%,
  physical mode matching and full81OUT controls pass. Zero-velocity OUT2
  excites the growing subspace; PL neutral control PARTIAL,0ODE/FEM/p16.
