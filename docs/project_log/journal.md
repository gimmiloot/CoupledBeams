# Journal

Здесь ведётся рабочий журнал проекта: этапы, решения и важные исследовательские заметки.

## 2026-10-09

- Completed explicitly authorized FEM-1R: one .020 C3D10 grid with20752 nodes/
  12687 elements,one24-mode CalculiX job,all8 shape matches stable. All inherited
  .1% mesh checks pass (.01649-.06880%); four-grid frequencies decrease with
  shrinking absolute changes. Updated finite-grid discrepancies.47666% axial,
  1.36648-1.74302% bending,4.40336-4.97569% twist retain continuum/BC/warping
  qualifications; no extrapolation or fitting. Parent94 artifact hashes/source
  data preserved,0 new1D solves/old-job repeats,32.84s/1200s,353.41MiB peak.
  Separate thin continuation CLI/config/tests,bundle63d44daae533389c,2 figures,
  appended canonical report and NLSP-D10/K11. Linear baseline ready for a bounded
  FEM-2 decision if separately authorized; no fifth grid/static/dynamic execution.

## 2026-10-08

- Completed the authorized full-family rectangular linear FEM-1: pre-FEM h=.10
  puts first MH acoustic at sorted8; generalized CT retained,8 roots complete in
  guarded omega3.4651 window. Three audited C3D10 meshes yield24 real modes each;
  all8 MAC identities/4 families PASS, zero extra in-window/ambiguous matches.
  Mesh convergence/all-family comparison PARTIAL:5/8 exceed preset.1% medium-fine
  tolerance. Apparent fine differences.49% axial,1.40-1.81% bending,4.46-5.04% twist
  retain mesh/end/warping qualifications. New scoped CLI/config/tests/report,
  strict FRD adapter (old reader/code untouched),3 PDF/PNG figures and D09/K10.
  Three sequential CCX calls, no extension/fourth mesh/nonlinear jobs;26.31s
  primary numerical work/3600s. Existing readiness edits and frozen files preserved.

- Audited local3D FEM readiness without solver jobs, new meshes, installation,
  archive expansion or code/model/result changes. Gmsh4.15.2 and cached CCX/MT2.22
  x64 version probes confirmed; expected CCX -v exit201 verified from local source.
  Historical linear solid evidence and reusable C3D10/DAT/FRD/convergence functions
  distinguished from angular joint-model/MAC issues. Local2.22 manual confirms
  STATIC/DYNAMIC with NLGEOM; new rectangular/static/transient workflows still
  need adaptation and execution tests, including constitutive/initial-state/load/
  energy contracts. Added one readiness report/navigation and one substantive
  NLSP-K09 infrastructure record; historical scientific statuses remain unchanged.

- Completed bounded physical sanity checks of the frozen planar quartic action:
  reused p64 .05 one-T1 history and exactlyONE new p64 .025 tight run toT1.
  Full-history amplitude maxima show2/4 patterns; normalized leading deviations
  reduce approximately4x. Exact sign symmetry and formal EA/(8L) bulk stretching
  checks pass; strain/reaction tables retain finite-p strong residuals and explicit
  degree4 rotational truncation accounting. Overall DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS;
  no physical/3D validation or independent small-amplitude p/time certification.
  New focused CLI/config/tests/note,3 PDF/PNG figures and NLSP-D08/K08; frozen
  model/RHS/Jacobian/IC profiles/BC/basis and historical strict/PARTIAL untouched.
  One ODE118.55s, primary numerical work237.72s conservative/600s; matching cache/
  report/plot perform zero new ODE/BVP/eigen/symbolic work. No next stage.

- Extended the unchanged prepared four-field IVP to one fixed LINEAR T1 with
  exact saved p48/p64 initial coordinates and tight/allowed-extra settings.
  Exactly3 exploratory runs complete; prefix/temporal/energy-mass PASS, spatial
  7/8 PARTIAL (theta_t max1.590e-4>1e-4). Absolute differences grow1.94–2.73x;
  full-horizon scale changes are qualified separately. No model/RHS/IC/BC/basis
  change or repeat projection/MP/BVP/eigen/symbolic audit. Added one config and
  explicit --one-T1 mode, optional memory-mapped history buffer, blockwise
  all8/time-window diagnostics,3 PDF/PNG figures and targeted tests. ODE150.29s,
  primary numerical work232.85s/1200s. Historical sources/strict XFAIL preserved;
  appended NLSP-D07/K07. Default short horizon remains0.1T1; no next stage.

- Completed a separate bounded prepared-state precision/feasibility continuation.
  Exact Gram/analytic moments and one initial-only endpoint-constrained L2 policy
  represent the same frozen target at p48/p64; initial1e-6 checks pass. Independent
  MP45/70 quadrature evidence explains stored-float strong/weak discrepancy; the
  unchanged relative2e-12 gate remains FAIL, strict verification PARTIAL.
  Explicit EXPLORATORY_NOT_CERTIFIED mode completes exactly3 short0.1T1 runs:
  temporal8/8 PASS, spatial7/8 PARTIAL (theta_t max1.405e-4>1e-4). Energy/safety/mass
  pass;17.34s integration,53.49s charged numerical work. Added explicit-q0 runner
  input/prefix metadata while preserving default numerical path/frozen action.
  Historical reports/results/D-K records retained; appended NLSP-D06/K06.
  Cache/report/plot perform zero new ODE/BVP/eigen/history evaluations. No full5T1,
  new p/amplitude, physical correction or next research stage.

- Prepared-state final verification:190 PASS,2 strict XFAIL,8 deselected;
  zero ODE/model eigensolves. The failed relative weak/action checks remain
  explicit qualifications; frozen files, staging and D/K append-only history
  are preserved. Cache/report/plot, link checks and git diff --check pass.

- Завершён [prepared initial-state gate](../theory/planar_prepared_initial_state.md)
  после явного задания пользователя. Исходный main HEAD0b2f31a77814340c32ce77eac85c2ef6ea3aa910,
  staging и все frozen scientific files/historical bundles сохранены.
  Stat/harm p64/p96 и endpoint jets PASS; общий O2 state и единственная
  degree<=5 Theta3 приw3=0 дают compatibility through cubic order.
  Finite epsilon4/5 residuals сохранены, periodic orbit не заявлена.
- Full four-field initial projection PARTIAL: p48 сохраняет c_ss endpoint
  error1.94e-5 и O2 force coefficient residual1.77e-6; p64 O2 u/c проходит,
  но theta_ss endpoint error3.03e-6 не проходит1e-6. Обе разрешённые пары
  заблокированы до trajectories. New short temporal/spatial NOT_RUN;
  NLSP_PREPARED_INITIAL_STATE_PILOT=PARTIAL. Улучшение nonlinear convergence
  не установлено; old eight-component norms восстановлены без ODE.
- Дополнительный strong/weak floating check проходитabsolute2e-12, но
  неrelative2e-12 приp48/p64 (1.14e-11/4.31e-11); energy power PASS.
  Эти проверки сохраняются strict XFAIL, без правки frozen RHS или gate.
  Finalbundle5ea8d41faf8ede54: primary3.857s,0 ODE/model eigensolves,
  two PDF/PNG figures, cache/report/plot zero evaluations. LONG CLOSED,
  EB/RLB-KV PAUSED, angular same-clamp reference UNAVAILABLE; новых этапов нет.

- Final verification:220 unique checks PASS (215combined+3Timoshenko+2source-cache);
  три старых ODE tests намеренно deselected. Finalcacheb3ea4eb6ac95d6e1
  имеет прямые audit/reference hashes; numerical AST/data не изменены.
  Links и protected file/source/staging preservation PASS.


- По отдельному заданию выполнен [контроль ведущего продольного отклика
  второго порядка](../theory/planar_second_order_axial_response.md).
  Initial main HEADb608aa118247819dfd118ddafb2a0a6155512c80, Version0.6.1;
  protected action/RHS/Jacobian/IC/BC и historical bundles сохранены.
  Это диагностическая asymptotic specialization исходной задачи, не новая
  production-модель и не сокращение четырёх независимых nonlinear fields.
- При epsilon_a=A/h0 общий continuous background W=h0*w_hat,
  Theta=h0*theta_hat использует прежние omega1/T1 при всех p. Forcing
  независимо извлечён из quartic action/cubic residuals: constant и2omega1,
  включая положительный inertia source. Все2(p−1) M-H координаты retained;
  exact-time matrix function и аналитические скорости/ускорения без ODE.
  Forced energy проверяется по power identity, не как постоянная энергия.
- Primary p16/24/32/48/64 и один conditional p96 выполнены в fixed1200s
  actual-cost scope. Ошибка relative-parent metadata после сохранения p96
  spectral state устранена адресным resume без повторной eigendecomposition;
  исходные primary bundle и checkpoint сохранены. Spatial convergence PARTIAL
  отделена от NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC=COMPLETE; eigen/forcing/
  quadrature/exact-time controls PASS. Неравномерное изменение differences
  с p не превращается в утверждение об ошибке физической модели.
- На refined p64→96 толькоu2 проходит обе нормы1e-3; c2,u2_t,c2_t
  остаются выше gate. Finite-dimensional exact-time control проходит,
  но velocity difficulty присутствует уже без time-integration error.
  Sampling799581→1599159: изменение sampled maxima2.71e-5, контрольPASS;
  это не continuous supremum certificate. Whole-time common-space evolution
  преобладает над хвостом, без phase-only или sole-cause утверждения.
  Charged numerical work373.599s/1200s;6primary M-H eigh,2cached restores,
  0ODE/root solves. Final bundleb3ea4eb6ac95d6e1; p>96 не рассчитывался.
- Historical nonlinear comparisons сохраняют actual timestamps: p24/p32
  large-amplitude до5T1, p32 small до5T1, p24 small только сохранённый prefix,
  nonlinear p48 только0.1T1. Общий continuous фон диагностического forcing
  не объявлен абсолютно идентичным projected semidiscrete background.
  Three PDF/PNG figures; response traces показывают разрешённый0...0.1T1,
  convergence norms покрывают5T1. Новых time integrations0.
- Начальная incompatibility не устранялась; leading response не заменяет
  полную cubic trajectory и не включает nonlinear feedback. Старые pilot/
  recovery PARTIAL сохранены, full nonlinear spatial convergence unresolved.
  LONG closed, EB/RLB-KV paused, angular same-clamp reference UNAVAILABLE.
  Новый basis, IC, physical model, amplitude или следующий этап не выбраны.

## 2026-10-07

- По отдельному явному заданию выполнена
  [адресная диагностика planar solver](../theory/weakly_nonlinear_planar_time_pilot.md).
  Исходный main HEAD1510d75c106a28a4da899c7eea1a337f11791ce6, Version0.6.0;
  старый bundle c97287772bc461ef и все его hashes/timestamps проверены,
  trajectories/manifest не пересчитывались и не менялись. Прежние L2/max
  разности воспроизведены; физическая L2-проекция показывает преобладание
  разности эволюции в общем пространстве, без доказательства phase-only причины.
- Из protected quartic action независимо подтверждён начальный axial trace
  u_tt=±1.10185433e-5 при A=.0025, порядок A². Для остальных трёх полей
  traces исчезают по continuous linear eigenpair identities. Это mismatch
  гладкой совместности у неподвижных границ, не недопустимость weak IVP,
  не доказанная ошибка программы и не единственная причина всех разностей.
  Initial fields, четыре independent fields, clamps и Shen basis сохранены.
- В existing planar helper разделены запросы energy/gradient/Hessian;
  variable mass и все inertia terms сохранены. Pointwise V/derivatives/RHS/
  Jacobian/energy/weak equivalence PASS. Short old/new p32 на0...0.1T1
  имеют нулевые разности всех восьми metrics; runtime9.02→7.13s.
- Всего ровно3 short controls,0 new full runs. p48 strict short13.81s;
  guarded full forecast863.33s превышает remaining859.59s. Charged total
  profiling/integration40.41s при fixed900s budget, без повышения лимита.
  NLSP_PLANAR_P48_SPATIAL_CHECK=REFINEMENT_DEFERRED_BY_BUDGET;
  NLSP_PLANAR_SOLVER_RECOVERY=PARTIAL. Полное p32→48 сравнение отсутствует,
  smaller-amplitude old PARTIAL control не продолжался.
- New evidence results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff/
  сохраняет original execution db47d4efb6941bed с прозрачным cache revision
  для plot lookups и двух guards в неисполненном full-run path: проверка
  PASS всех коротких controls и различение master/comparison timestamps.
  Исходные numerical data не менялись, повторных интегрирований нет.
  31 targeted и165 combined tests PASS;3 ODE tests намеренно deselected.
  Две diagnostic PDF/PNG figures;
  report/plot-only и matching cache compute выполняют0 integrations.
  LONG closed, EB/RLB-KV paused, angular same-clamp reference UNAVAILABLE
  сохранены. Новый научный этап или смена basis/model/IC не выбраны.

- Завершён bounded [первый planar time pilot](../theory/weakly_nonlinear_planar_time_pilot.md)
  по explicit user request f352f4d1-b074-4046-871a-313a14573314. Initial main
  HEAD7b3d667d5418cde247a6ca09fe9186945177537b,16modified/8untracked prior
  audit files сохранены вместе с index; old helpers/appendix/source/reference
  hashes не менялись. Новые helper/CLI/config/tests/note дискретизируют
  accepted quartic action, не exact untruncated model. Четыре independent
  поля, clamps по значениям, без slope/inextensibility/quasistatic-c условий.
- До main зафиксированы p16/24/32,nq2p+1, две A/h=.05/.025 и5T1,
  Radau analytic Jacobian, три time levels, numerical gates и smoke-based
  2400s scheduling/480s per-integration budget. В main6 full histories,
  включая обе finalp32 amplitudes, и один saved p24 small prefix до2.648T1.
  Root solves0, model derivation1;5,090,009 RHS,13 Jacobians,26 Newton LU,
  5,089,987 mass factorizations. Integration stopped at deadline; saving/
  postprocessing completed at2425.512s. Бюджет и допуски не повышались.
- Linear MH first3 и Tim first3, continuous initial pair/projections,
  exact-in-time semidiscrete controls, weak/action/Jacobian/RHS energy и
  source-polynomial saved-snapshot energy проходят. Last temporal pair
  PASS. Last spatial p24→32: q u/c .274%/3.058%, velocities u/theta/c
  2.200%/.1432%/5.750%, выше заданных gates. Поэтому общий
  NLSP_PLANAR_TIME_PILOT=PARTIAL; small-amplitude neighboring-p check
  incomplete, p48/extra tightening не запускались в исчерпанном бюджете.
- Finalp32 drift9.53e-11/5.66e-11; mass positive, all sampled safety measures
  bounded. Inducedu/c nonzero; normalized departures decrease with smallerA,
  но complete four-field nonlinear effect/physical validation не объявлены.
  Три PDF/PNG figures из full trajectories,247 unique targeted/regression
  tests PASS, cache compute/report/plot-only zero integration/root/derivation.
  README/CHANGELOG/navigation обновлены; append-only NLSP-D02/K02 фиксируют
  PARTIAL и остановку. LONG closed, EB/RLB-KV paused, angular same-clamp
  out-of-plane reference UNAVAILABLE. Нет periodic/Floquet/threshold/maps,
  3D truth, source/model fitting или автоматического следующего этапа.

- По новому явному заданию завершён отдельный
  [семиполевой nonlinear spatial action audit](../theory/weakly_nonlinear_spatial_rod.md).
  Исходный main HEAD7b3d667d5418cde247a6ca09fe9186945177537b, checkout clean;
  до редактирования сохранены документы,658 tracked hashes, index hash и104
  protected old-result hashes. Supplied MD сохранён побайтово SHA87d7ce3c...063f.
  Новая модель не переоткрывает LONG и не переносит EB/RLB spring/KV.
- Принятые V0/kinematics вынесены в isolated helper/note, frozen equations и
  oldsolvers/variants/tolerances не менялись. Independent quartic-action A и
  matrix-exponential-balance B:21/21 exact zero differences,21/21 supplied MATCH.
  Exact Fraction ring (SymPy unavailable, no installs); boundary/energy/
  reflection/dimension/axial/independent planar/rigid/split identities PASS.
  Numerical evaluator uses full Rodrigues/Jr with analytic derivatives.
- Fixed manufactured3profiles,2samples,epsilon .04→.0025: aggregate error
  6.76753e-8→1.03264e-12, order≈4; final small rotation/c components below
  predeclared reporting floor flagged. No PDE integration or trajectory claim.
  G20 CT from isotropic generalized rectangular torsion, correct book-axis swap;
  section-rotation clamp distinguished from historical book_slope_clamp.
- Bounded direct first6+guard7 and3 MH family pairs match independently
  established old blocks. Direct old min-max/scalar count9/9 belowomega2.7;
 10 direct profiles, two split profiles/roots including c/R PASS. Rank14/
  duality0/45/90; saved planar first3 at45/90 pass. Out-of-plane angular
  same-clamp spectrum reference unavailable; no forced source comparison.
- Final bundle results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/:
  model/identity/coefficients,polynomials/differences,profiles/brackets/counts,
  manufactured errors,joint diagnostics/generated appendix. Runtime2.07s,
 715 boundary evaluations, cache/report-only zero derivations/roots.
  Earlier developmental bundles retained. New57 tests plus189 unchanged
  relevant MH/Bishop/Yartsev/rectangular Tim regressions:246unique PASS.
  Initial test metadata-key mismatch corrected without physics/tolerance change;
  review added explicit frozen-reference cache hashes and strict parser guards.
- Registered Crespo1988 journal scan (10pages,SHA2949bede...4e128), checked
  relevant page images. No-shear/u-order/EA-scaling not imported. README,
  CHANGELOG/navigation/assumptions/scripts guides updated; append-only NLSP-D01/
  NLSP-K01 added. AllNLSP audit statuses PASS, supplied MATCH; c continuity
  remains reduced variational closure. No physical nonlinear validation,
  nonlinear time integration, Floquet/periodic/modal reduction, critical amplitude,
  maps, forcing/damping/FEM or fitting. Afteraudit stop; no automatic next stage.

- Завершены [production Lambda(beta) large checks](../theory/mindlin_herrmann_timoshenko_lambda_beta_large_checks.md)
  по отдельному implementation-only заданию. Initial main/HEAD
  0f961acb227b504439408d3484b13a72a7acf3be:16modified/21untracked.
  Diff/status/source и old-bundle hashes сохранены в OS temp; изменения
  пользователя сохранены. Physics helpers/API/source fixtures не менялись.
- Canonical Lambda подтверждена equations.tex/project_rules/main_note:
  Lambda^4=rho*A_ref*omega^2*l_ref^4/(E*I_ref),fixed .5/.20/.05 reference,
  factor300. RLB-2B использует локально square parameter, RLB-2C явно
  восстанавливает canonical mapping; старые notation/formulas не переписаны.
- Один новый CLI/config использует два existing arm kernels и тот же
  invariant joint operator. Length mu0/.25/.5 и contrast delta_h0/.2/.4,
  beta0:2.5:90. Shared baseline даёт185 unique main cases,222 logical CSV
  cases плюс20 sparse negative-parameter symmetry controls. Total mass .01.
  Independent stepped beta0 reference uses literal positive-X state continuity.
- Все10 statuses PASS. Length direct collapse6.01e-12 relative frequency;
  stepped-reference1.12e-15. Swap length6.15e-12/section6.06e-12; full
  kinematic/resultant L2 <=4.49e-11. Baseline frequency3.99e-12,
  Lambda1.992e-12. Every first12+guard13 count/quality passes.
- 5 full seed scans и200 sorted-frequency bracket predictors, no fallback,
  74 bounded frequency subdivisions;0 ceiling expansion. Determinant32196,
  count4879,catalog7862,reference583,independent-flag975 evaluations;
  final compute137.27s. Это numerical acceleration,не modal continuation.
  68 curvature flags retained/reviewed at same beta by independent state QR,
  max frequency difference7.21e-12; no smoothing or beta refinement.
- Max clamp6.79e-12,PDE1.48e-11,force1.95e-11,moment6.70e-13,c2.30e-14,
  R1.25e-14,Gram5.17e-11. Nonzero SVD condition2.848--10639.1, unchanged
  tolerances.36 new +233 prior +4 rectangular Tim=273PASS,7 unrelated
  deselected. Two3-panel PDF/PNG figures saved in canonical bundle
  results/mindlin_herrmann_timoshenko_lambda_beta_large_checks/356c3be4953268f2/.
  Prototype retained; final version adds Gram/independent anomaly checks.
- README/CHANGELOG и navigation updated. Literature/BibTeX/other memory,
  old solvers/tests/results preserved. Independent sorted positions only,
  no energy classification/across-beta MAC/branch identity, sensitivity or
  applicability analysis, coefficient fit, damping/FEM/nonlinear work.
  c1=c2 remains reduced variational closure,не direct finite3D joint derivation.
  Эти geometry checks достаточны для implementation gate в заявленном
  диапазоне; дальнейший scientific parameter study автоматически не запускался.

- По новому заданию завершён [bounded thickness screening](../theory/coupled_longitudinal_theory_thickness_screening.md):
  пять s_h=1/1.25/1.5/1.75/2, beta0/45/90, три прежние theories.
  Initial main/HEAD0f961acb227b504439408d3484b13a72a7acf3be:16modified и
  17untracked prior-stage files. Snapshot diff/status/hashes сохранён в
  OS temp; прежние edits/source/solver/test/results сохранены, Git mutations нет.
- Exact Fraction/scaling audit PASS: A,m,C,S~h; I,J,j,H,B,r~h^3;
  ratios contain I/A=h^2/12, common b cancels. Static MH layer sqrt(H/EA)~h;
  full spectral h^2 law не утверждается. Масса растёт, κ=5/6 не меняется.
  Начальный ceiling3.75 сохранён; unchanged bounded basis требует prefix
  certificates below Tim optical cutoff. Реальные windows для s_h1.75/2:
  3.2942708/2.8010653, всё ещё покрывают13 roots. Above-cutoff tail не считался.
- Все45 inventories PASS, без ceiling expansion/catalog retry/failed intervals.
  15 direct beta0 checks PASS по first13 frequencies/full kinematic forms:
  max relative frequency6.35e-12,L2 profile8.41e-12. Artificial joint remains
  transparent, c/R variational closure не переинтерпретирован как3Delasticity.
- Global E/MH0.2708143% (beta0,s_h2,k4), RL/MH0.4792888%
  (beta0,s_h2,k12), против0.14103/0.17514% baseline. Рост maxima имеет
  geometry/position dependence и не равен s_h^2. Максимальное отличие RL от
  MH больше отличия E во всех15 cases; все best overlaps diagonal. Max D_c0.1124794
  (beta90,s_h2,k5), без energy fractions/classes и без across-case tracking.
- Close pair beta90,s_h1.25,k9->10: MH gap3.67266e-5, nonzero SVD
  condition52194.6. Residuals и Gram проходят unchanged tolerances;
  это diagnostic context, не crossing/veering/resonance conclusion.
  Max clamp5.06e-12, force6.39e-12, joint rotation6.14e-12, Gram1.48e-10.
- 35 new +198 prior hierarchy/MH/Bishop +4 rectangular Tim checks =237PASS;
  7 unrelated tests deselected. Added only CLI/config/tests/theory note;
  source index/BibTeX и branch memory unchanged. README/CHANGELOG/navigation
  обновлены. Canonical local bundle7be0fce968fd2b35 под
  results/coupled_longitudinal_theory_thickness_screening/ содержит scaling,
  roots/counts/brackets/profiles/overlaps/D_c/direct checks и3 figures.
  Нет fitting/power-law fit/applicability threshold, h>2, beta refinement,
  unequal-arm/section/material sweeps, FEM/3D truth или nonlinear equations.
  Дальнейший thickness range автоматически не расширялся.

- Завершён [bounded hierarchy screening](../theory/coupled_longitudinal_theory_hierarchy_screening.md)
  по актуальному прикреплённому заданию. Initial main/HEAD
  0f961acb227b504439408d3484b13a72a7acf3be:16 modified tracked +12 untracked
  prior-stage files. Snapshot initial diff/status/hashes сохранён в OS temp;
  пользовательские изменения сохранены, Git mutations не выполнялись.
- Переиспользованы Bishop H=0 axial kernel, неизменённые MH/Tim arms и
  geometry/joint maps. Новый comparator-layer имеет только u,w,theta,N,Q,M;
  planar Love J=nu^2 rho Iy и variational N=(EA-J omega^2)u_x, без c/R и
  дополнительных slope constraints. Common centroid/bending clamps не
  отождествляются с independently resolved MH contraction clamp.
- Оба comparator gates PASS до screening: beta0 direct CC recovery на
  splits .5/.35, duality/rank6, beta45 swap/reflection eigenpair remapping.
  Relative frequency recovery <=7.00e-12; kinematic remap <=2.83e-14.
  Energy/Schur count корректно сохраняет global roots, совпадающие с arm
  poles; effective count queries записаны отдельно, частоты не сдвигались.
- Один G20 equal-arm case, beta0/5/15/30/45/60/75/90. Все24 inventories
  certifiable12+guard13; MH5/45/90/beta0 reference roots переиспользованы
  read-only, MH15/30/60/75 computed неизменённым solver. Нет ceiling
  expansion, unresolved/failed numerical intervals. Max force residual
  3.01e-11, clamp1.49e-12, nonzero singular condition1060; old tolerances.
- Max spectral difference E/MH0.1410296% (beta0,k5), RL/MH0.1751407%
  (beta0,k11). Все geometric row/column best matches diagonal; нет
  position-exchange cases. Minimum gap около1.2637% (MH,beta75,k9->10)
  только diagnostic context. MH reference не объявлялся truth.
- D_c с тремя отдельными norms, без modal classification:0.04162--0.08575
  для defined fields;9 beta0 positions NOT_DEFINED_SMALL_FIELD. Theta
  overlaps используют machine-scale SMALL_NORM. Gauss200->300 checks
  beta0/45/90 pass. Нет energy fractions/types, across-beta tracking,
  автоматической перестановки roots или arbitrary overlap threshold.
- 32 new +166 unchanged +4 rectangular Tim tests =202PASS;7 unrelated
  deselected. Source index/BibTeX/fixtures, previous solvers/tests/bundles,
  other-branch memory unchanged. README/CHANGELOG и navigation обновлены.
  Canonical local bundle9a35ff23c3f43c66 под
  results/coupled_longitudinal_theory_hierarchy_screening/: full tables,
  brackets/counts, profiles, overlaps, D_c, diagnostics и3 figures.
  Earlier smoke/draft attempts retained; corrected dictionary-key and
  direct-catalog ceiling bookkeeping errors без изменения physics/tolerance.
  Нет beta refinement, applicability/safe-prefix conclusions, mu/tau/thickness
  scans, коэффициентного fit, FEM/nonlinear work. Следующий этап не запускался.

## 2026-10-06

- General-beta continuation принят по отдельному заданию после beta0PASS:
  [canonical joint note §§9–16](../theory/mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract).
  Current main/HEAD0f961acb227b504439408d3484b13a72a7acf3be содержит16modified
  tracked/10untracked prior-stage files; initial diff/status/hashes сохранены
  в OS temp, прежние изменения сохранены. Projectβ восстановлен из
  RLB physical coordinate contract/helper, не введён новый angle definition.
- Existing joint-helper расширен общей frame_boundary_matrix; old beta0 API
  guard делегирует той же сборке. Source/single-arm physics,κ,variants,
  frozen equations, old CLI/tests и beta0bundle3059d70b1b50ea2e не изменены.
  c остаётся scalar reduced DOF;θ/M меняют знак при physical reflection,
  оставаясь invariant при proper rotations. Force map выведен из dual work.
- All8 structural gates PASS: rank8/JJᵀ=2I, work max6.67e−16, beta90
  axis error6.13e−17, frozen beta0 matrix/frequency/profile differences0
  на3splits/54 modes. Small-angle first3 frequency relative differences:
  6.67e−15/3.76e−11/3.76e−7 для1e−6/1e−4/1e−2deg. Swap/reflection
  6 low eigenpair controls доf*=.7 — same geometry remapping, не tracking.
- После structuralPASS выполнен один equal-arm G20 pilotβ5/45/90,totalL1.
  Derived exact energy/Schur count использует independently verified4MH/9Tim
  fixed-arm pole catalogs и negative nodal inertia. Complete inventories
  23/24/24 доfixed ceilingf*=3.75 покрывают12+guard13; subdivisions0/0/2,
  failed intervals0. Все71 normalized full modes/resultants сохранены.
  Max joint residual1.37e−11, mass Gram3.64e−11, nonzero condition461.22;
  unchanged tolerances соблюдены. Коэффициенты/знаки не подбирались.
- MHTIM_GENERAL_BETA_JOINT_GATE=PASS,37 new +129 unchanged MH/Bishop/beta0
  +4 rectangularTimo regressions=170PASS. Expanded canonical note,
  assumptions/research/results/scripts navigation,README/CHANGELOG updated.
  Source index/BibTeX/fixtures/memory не менялись; source closure остаётся
  variational reduced1D,не3Delasticity welded-region proof. Нет scientific
  angle sweep, hierarchy/applicability study, across-beta MAC/tracking,
  nonlinear/FEM work. Existing reduced-limit checks только regression tests.
  Result root results/mindlin_herrmann_timoshenko_general_beta_joint/;
  здесь остановка, следующий parameter study не запускался.

- По текущему заданию принят [published reduced rigid-joint closure](../theory/mindlin_herrmann_timoshenko_rigid_joint.md)
  с общими d,theta,c и dual nodal balances. Rucka PDF4–5 (14),(15),(26)–(28)
  повторно сверены по изображениям: c continuity inferred from common scalar
  DOF assembly, не отдельная напечатанная формула. Jang PDF3,5 (12),(28)–(30),
  (36) подтверждает end DOF/natural pair. Это variational reduced1D closure,
  не direct3D elasticity proof сварного finite joint; не названо ошибочным.
- Initial main/HEAD 0f961acb227b504439408d3484b13a72a7acf3be содержал16 modified
  tracked и6 untracked previous-stage files. Initial diff/status/hashes
  сохранены во временной папке ОС; старые edits не удалены. Production preset,
  single/source solvers, source fixtures/PDFs, old tests, baseline equations,
  source index/BibTeX и memory другой ветви сохранены побайтово.
- Endpoint signs независимо получены из energy [p delta q]_0^L. Both local
  x run positive outer clamp->joint; beta0 t1=(1,0),n1=(0,-1),t2=(-1,0),
  n2=(0,1),joint signs++. Global R transparency соответствует R1+R2=0
  в этих координатах. 8 conditions/rank8 и arbitrary-state dual work PASS.
- Unchanged direct single bundle342ce44bff81c36f переиспользован с hash
  checks; reference roots не пересчитывались. Bounded independent joint
  matrices на splits .5,.35,.65 дают7MH/11Tim, saturated min-max counts,
  first12+guard13. Max relative frequency difference4.56e-13; c L2 error
  4.81e-13, min massMAC>.9999999999999995, max joint residual1.78e-12.
  Reflection .35<->.65 и independent segmented expm/QR passes. Retried/failed
  intervals0; никаких signed fitting/tolerance/physics changes.
- Все7 gates включая MHTIM_BETA0_JOINT_GATE=PASS. 31 new tests +98 unchanged
  single/source/Bishop +4 existing rectangular Timo regressions=133PASS;
  import/smoke/provenance/cache checks проходят. One new joint helper/CLI
  and canonical note, navigation/README/CHANGELOG updated. Generated bundle
  в results/mindlin_herrmann_timoshenko_beta0_joint/. Старые результаты не
  изменены. beta!=0 spectrum/closure не валидированы; коэффициенты не
  подбирались, nonlinear equations/FEM truth не вводились. Здесь остановка.

- По явному выбору пользователя принят отдельный project closure
  [Jang bare isotropic reduced M-H + Timoshenko](../theory/mindlin_herrmann_timoshenko_single_rod.md#14-production-formulation-decision):
  PRODUCTION_MHTIM_FORMULATION_SELECTED=PASS. Project κ=5/6 найден в
  frozen G20 rectangular contract и существующем comparator/test; его
  GI/GA convention проверена, κ один для MH и Tim, j=r=rhoI.
  PRODUCTION_MHTIM_KAPPA_RESOLVED. Это не установленный source κ Jang.
- Настоящий Fernandes PDF energies-15-07725-v2.pdf зарегистрирован:
  publisher article, 26 pp., DOI10.3390/en15207725, SHA25627cb6de0...c35a2fc.
  PDF6–7 (19)–(21) подтверждает повторное использование Ng factors, но
  Lamé normal block не заменяет выбранный reduced C,nu*C. Source I не
  превращён в polar Ip или выдуманный rectangular specimen. Ref32 Doyle
  зарегистрирован как citation chain, книга не проверена. Предыдущий
  unresolved mapping audit сохранён, 12/pi² не принят в production.
- Один finite straight G20 rod E=rho=1,nu=.3,b=.2,h=.05,L=L_ref=1:
  source (1),(12),(28) подтверждает essential u=c=w=theta=0 как full clamp
  within planar 4-field kinematics. Получены variation/PDE/resultants/state,
  bounded exact cos/sin+anchored exponentials и независимый short-step
  state expm/QR. Min-max upper counts насыщены: 7 MH и11 bending roots с
  guards, без нового WW/map/transfer-products framework. Finite spectrum
  PASS. First6 MH f*=.500706144,1.001229213,1.501381967,2.000968673,
  2.499780436,2.997590005; contraction wave cutoff11.558994422 выше low
  inventory. Это normalized control, не experimental/source-table fit.
- Numerical hierarchy elementary/planar Rayleigh–Love/MH passes, bending
  один и тот же. Overall HIERARCHY_SINGLE_ROD=PARTIAL_PASS: reduced axial
  theories не имеют independent c и не описывают её clamp boundary layer;
  u′=0 к second-order equation не добавлено. Planar RL J=nu²rhoIy,H=0,
  не polar standalone Bishop. Angular-joint contraction condition открыт.
- Добавлены named preset/finite functions в existing MH helper, один finite
  CLI/config и19 tests. Итого98 targeted +4 existing rectangular Timo
  regressions=102 PASS; source-check, smoke и hash-validated zero-root reuse
  проходят. Rucka/Jang outputs/records/cases/tolerances точно сохранены,
  source Jang numeric κ остаётся null. Initial diff/status main/HEAD
  0f961acb227b504439408d3484b13a72a7acf3be сохранены во временной папке;
  13 protected files byte-identical. Все предыдущие user changes сохранены.
  README/CHANGELOG/navigation обновлены для нового finite command.
  Memory, frozen baseline, Bishop, anisotropic/viscous/nonlinear workflows
  не менялись. Новая two-beam MH system, joint BC, nonlinear derivation и
  3D FEM не выполнялись.

- Выполнен [rectangular M-H prescription gate](../theory/mindlin_herrmann_timoshenko_single_rod.md#13-production-rectangular-m-h-correction-prescription).
  Начальный main/HEAD 0f961acb227b504439408d3484b13a72a7acf3be: tracked diff
  пуст, два untracked PDF сохранены. hdl_85788.pdf — Ng 2014 accepted
  version; ssrn-5985611.pdf — Elishakoff–Tharu review, not peer reviewed,
  не ожидаемый Fernandes et al. 2022. Последний не найден среди 92 PDF
  проекта; по заданию ему не присвоены выдуманные metadata/key/hash.
- Ng PDF 10 (2) прямо содержит S1=12/pi², S2=S1[(1+nu)/(.87+1.12nu)]²;
  H=S1*GI и j=S2*rhoI — прямой factor-role mapping. Normal block источника
  использует напечатанный 3D Lamé lambda и отличается от текущего reduced
  C,nu*C. Exact stationary reduction даёт EA/(1−nu²), не EA. Circular
  squared corrections препринта не служат вторым rectangular подтверждением.
  Doyle/Graff зарегистрированы только как citation chain, книги не проверены.
- Итог RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED и
  PRODUCTION_MH_COEFFICIENTS_UNRESOLVED. Fixtures разделяют source variants
  и unadopted candidate rectangular_literature_default; helper/preset не
  добавлен. При nu=.3 только arithmetic candidate output:
  1.2158542037080533 / 1.4127769143961029. Production model/cutoff run не
  выполнен. Старые Rucka/Jang records/cases/tolerances и fixed-point outputs
  сохранены; техническое tuple/list различие snapshot устранено одинаковой
  JSON representation, без numeric tolerance. 79 targeted tests + 3 existing
  rectangular Timoshenko regressions PASS; source hash-check PASS.
- Source index/BibTeX, существующие literature/theory notes, assumptions,
  research/results navigation, journal и CHANGELOG обновлены. Local audit
  хранит page snapshots/provenance; прежние results не пересчитывались.
  README и script guides оставлены: пользовательский CLI/API не меняется,
  unresolved status остаётся актуальным. Memory другой ветви, Bishop,
  Timoshenko kappa, source CLI и frozen baseline сохранены. Coupled model,
  angular-joint BC, nonlinear equations, FEM и coefficient fit не выполнялись.

- Выполнен [single-rod M-H + Timoshenko source audit](../theory/mindlin_herrmann_timoshenko_single_rod.md).
  Начальный main/HEAD 7ef59a5d8340d46735363d789ec557cca5144422 содержит
  незакоммиченные изменения предыдущей литературной регистрации и другого
  circular EB reviewer этапа. Их initial diff/hashes сохранены, staged
  содержимое и пользовательские изменения сохранены. Memory не менялась.
- По локальным Rucka и Jang восстановлены поля, разные reduced constitutive
  contributions Jang, corrected energies, PDE, boundary quantities, размеры,
  локальное разделение, acoustic/optical dispersion и low-k limits.
  Точное family mapping: K₁ᴹᴴ=K₁ᵀⁱᵐ=κ_b, K₂ᴹᴴ=K₂ᵀⁱᵐ=1. Fitted Rucka
  variant ему не соответствует: MH_SOURCE_VARIANTS_NOT_EQUIVALENT для
  опубликованных prescriptions. Source-specific факторы не стали defaults.
- Rucka Fig.4: square 6×6 mm, steel source properties/factors, 0–500 kHz;
  counts longitudinal/flexural по одному на всём [100,120] kHz подтверждены.
  Cutoffs contraction/shear 345680.639263 / 262945.871880 Hz. Jang Fig.9(a)
  восстановлен условно при **явном** κ_b=5/6: contraction cutoff
  1476250.985206 Hz независимо от κ_b, shear 779995.291328 Hz для этого
  input. Численное κ_b авторов не установлено; никакого fit по графику.
  Один common comparison использует geometry/material Rucka.
- Добавлены один helper, один CLI, source fixtures, targeted tests и
  канонический отчёт. Source reproduction MHTIM_VARIANT_DEPENDENT;
  PRODUCTION_MH_COEFFICIENTS_UNRESOLVED. Numerical/HF/full-operator checks
  проходят, графики qualitative/conditional qualitative. Первые tests нашли
  one-ulp strict comparison и ошибку renderer; исправлены с неизменными
  numeric contract/factors. Старые namespace bundles сохранены, current.json
  указывает актуальный. Plot-only и matching-cache reuse не вычисляют roots.
- README/CHANGELOG и theory/research/results/script navigation обновлены.
  Frozen equations, helpers, Bishop benchmarks и closed kinematic gate
  сохранены. Coupled rods, angular-joint M-H BC, L-joint production,
  nonlinear equations, anisotropic/viscous branches и 3D FEM не выполнялись.
- Финальные проверки: 74 tests (26 новых + 48 Bishop), затем 3 существующих
  rectangular Timoshenko regression checks — PASS. Source hashes, smoke
  compute/plot-only, matching-cache reuse с zero root evaluations и
  `git diff --check` проходят. Numeric gates заданы в fixtures до расчёта;
  display/source precision не подменены machine tolerances.
  Source-check дополнительно сверяет SI parameters/factors с печатными
  строками, чтобы изменённая geometry не сохраняла source-reproduction label.

## 2026-10-05

- Документационно завершён [circular EB reviewer diagnostic](../laminated_beams/circular_eb_rotational_spring_rigid_limit.md):
  сохранён итог 6/6 локальных descendants, ограничения target prefix,
  прежних direct MAC и немонотонности seed 06. Добавлены D25/K26,
  обновлены текущий контекст и навигация, зарегистрированы локальные
  junction sources. Расчёты и модели не менялись. Следующий выбранный
  вопрос — transmission/equilibrium и thin-joint asymptotics;
  FEM-калибровка реального узла отложена пользователем.

- Выполнена [литературная регистрация M-H + Timoshenko](../literature/mindlin_herrmann_timoshenko_sources.md).
  Начальное состояние: cwd/git root D:/PHD/CoupledBeams/CoupledBeams,
  main, HEAD 7ef59a5d8340d46735363d789ec557cca5144422; tracked diff пуст,
  три пользовательских untracked PDF сохранены под исходными именами.
  Найдены Rucka 2010, Jang–Park–Lee 2014, Banerjee–Ananthapuvirajah 2019
  и Liu et al. 2021; Banerjee сохраняет прежний citation key. Source index
  и bibliography синхронизированы, добавлена одна сравнительная заметка.
- Прочитаны постановки и source validations в пределах указанных разделов.
  Rucka: contraction — source ψ; K₁ᴹᴴ=1.1, K₂ᴹᴴ=2.1 и K₁ᵀⁱᵐ=.95 fitted
  по скоростям при 100/120 kHz; K₂ᵀⁱᵐ выбран по Lamb cutoff. Все отмечены
  SOURCE-SPECIFIC / NOT A PROJECT DEFAULT. Jang: разные reduced laws
  для M-H и Timoshenko, предупреждения к (A3), (A7), (A10). Banerjee —
  modular Rayleigh–Love/Timoshenko precedent, не M-H. Liu 2021 упоминает
  Bishop в обзоре, но не использует его в demonstration/Appendix.
- Оригиналы Mindlin–Herrmann и Martin–Gopalakrishnan–Doyle не найдены в
  локальном каталоге; сохранены как cited / full text unavailable, без
  выдуманных canonical records. Адресная внешняя проверка ограничена
  недостающими метаданными трёх работ; широкого поиска и OCR не было.
- Текущий кандидат combined in-plane model — M-H axial + Timoshenko;
  проектная реализация и валидация не выполнены. Bishop standalone,
  benchmarks и прежний COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED сохранены.
  README обновлён для новой навигации, CHANGELOG — для регистрации этапа.
  Assumptions и scoped memory другой ветви не менялись: новых проектных
  физических предпосылок не принято. Solvers, tests, results, article files
  и другие исследовательские ветви не менялись; расчёты и FEM не запускались.
- Документационные проверки: строгий синтаксический разбор всего BibTeX,
  69 уникальных ключей с полным соответствием source index, 56 DOI без
  дубликатов; 66 canonical PDF-путей существуют. SHA256 четырёх основных
  работ подтверждены, все 76 PDF побайтно сохранены; новые локальные ссылки
  разрешаются, `git diff --check` без ошибок. Сам BibTeX engine не запустился:
  незавершённая установка MiKTeX и недоступная запись его user configuration.
  Установка/настройка окружения не выполнялась; синтаксис проверен отдельным
  временным parser, без нового проектного скрипта.

- Выполнен [аудит общей кинематики Тимошенко—Бишопа одного прямоугольного
  стержня](../theory/timoshenko_bishop_single_rod.md). Начальный checkout
  чистый, main, HEAD a15e3af0e2ccc271a48df484e5d9865dc1cf4596.
  Знаки и оси восстановлены по текущему rectangular helper и его theory
  note; Marais и Banerjee проверены по локальным PDF. Yucel остаётся
  отсутствующим локальным текстом, новых источников не искали.
- Из двух displacement fields выведены деформации, все смешанные члены,
  точные моменты и вариационный оператор минимального поля. При центре
  cross coefficients равны нулю; смещённая ось даёт ненулевые члены.
  Но полный закон Гука для centroid-only contraction возвращает
  C11*I_y и G*A вместо E*I_y и kappa*G*A. Полная Poisson contraction
  добавляет изгибную градиентную инерцию и энергию. Итог
  COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED: гибридное замыкание не выбрано,
  combined spectrum/иерархия A/B/C остановлены на hard gate.
- Добавлены один algebra-only CLI, один файл тестов и канонический отчёт;
  README/CHANGELOG и навигация обновлены. SymPy/Lean недоступны, использованы
  ручной общий вывод и точная Fraction-алгебра без установки зависимостей.
  Исходные Bishop/Timoshenko helpers, литературные fixtures/benchmarks,
  equations.tex, source index/bibliography и scoped memory не менялись.
  Два стержня, угловой узел, нелинейные уравнения и FEM не выполнялись.
- Проверки нового этапа: 48 targeted tests пройдены (24 новых и 24
  неизменённых литературных), 4.37 s; `git diff --check` без ошибок.

- Выполнен отдельный узкий [литературный этап Рэлея—Бишопа](../theory/bishop_literature_reproduction.md):
  один диагностический модуль/CLI, транскрипция локальных PDF с SHA256,
  fixtures и целевые тесты. C, F и U=P=0 разделены; пониженные модели H=0,
  rigid mode, знак Γ′=−N, интерфейс, энергия и обобщённая массовая
  ортогональность проверены. Базовая теория/API не изменены.
- Marais §4: пять корней и конечная проверка полноты, независимый expm
  при 50/70 разрядах. (23),(24) однозначно реконструированы по (11);
  четыре печатные частоты не совпадают, только вторая проходит обе
  интерпретации четырёхзначной печати. Подбора параметров не было.
- Попов–Садовский: сохранены 30 округлённых отношений, включая `10.998`;
  минимум (6) ν=.336842781730, округляется до .337. Раздельно вычислены
  F-F, C-C и UP; (13),(15) согласованы. Утверждение с. 278 не подтверждается
  для (11); отдельные столбцы рис. 5 визуально отличаются. Исходных
  неокруглённых измерений нет; научные оговорки записаны в tracked отчёте,
  а не только в ignored results. Начальная попытка счётчика 30/31 и её
  ограниченный повтор сохранены, окончательный поиск хранит 31-й guard.
- Исходные изменения регистрации литературы сохранены; README и CHANGELOG
  дополнены из-за нового запуска и научного результата. Bibliography
  использована без новых изменений. Scoped memory другой ветви EB/RLB
  не перезаписана. Угловой узел, прямоугольная конструкция, Timoshenko +
  Bishop, FEM и нелинейная модель не начинались; этап остановлен.
- Завершение проверки: 24 целевых теста пройдены; все четыре CLI-действия,
  повторное чтение cache без решения и `git diff --check` проверены.
  Обнаруженное при дополнительном тесте усиление почти нулевого столбца
  пониженной F-F системы исправлено в численном масштабировании без
  изменения уравнений/допусков; история сохранена в локальном source audit.

- Зарегистрированы 11 новых локальных публикаций: шесть по продольным
  моделям (включая главу книги) и пять по нелинейной пространственной
  динамике. Исходное состояние: `main`, HEAD
  `a610b3405be655896c829dfc353737b0a5c6dd8c`, tracked diff пустой,
  11 untracked PDF. В `docs/literature/pdf/` проверены все 73 PDF:
  62 tracked, 11 untracked, ignored PDF в этой папке нет. Проверены также
  остальные локальные PDF; новых литературных источников вне неё не найдено.
- Синхронизированы [source index](../literature/source_index.md) и
  [bibliography](../literature/bibliography.bib); добавлены карты
  [продольных моделей](../literature/longitudinal_rod_models_sources.md) и
  [пространственной нелинейной динамики](../literature/nonlinear_inplane_outofplane_sources.md).
  В них разделены проверка метаданных, доступность PDF, прочитанные места
  и ещё не выполненные проверки формул/расчётов. Дубликатов новых публикаций
  не найдено; версии и переводы связаны, PDF не переименованы и не изменены.
  Из ожидаемых самостоятельных работ 12 полных текстов не найдены;
  они не скачивались и не объявлены прочитанными.
- Воспроизведение опубликованных задач Рэлея—Бишопа остаётся отдельным
  планируемым этапом; в checkout его результаты не обнаружены. Вывод
  нелинейных уравнений отложен до обсуждения с руководителем, условия
  нашего углового узла не выбраны. Научные расчёты, модели, прежние
  исследовательские статусы и scoped memory не изменялись.

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
