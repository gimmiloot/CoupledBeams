# Assumptions

For project-wide rules on source-of-truth priority, branch identity,
diagnostic workflow, and model-extension checks, see `../project_rules.md`.

Здесь фиксируются активные допущения и рабочие гипотезы проекта.

## Working Notes

## Separate prepared planar initial-data control (NLSP)

- The [prepared-state audit](planar_prepared_initial_state.md) defines the
  distinct case `prepared_axial_O2_with_cubic_endpoint_compatibility`;
  the historical zero-u/c pilot and its PARTIAL remain unchanged.
- epsilon_a=A/h0, W=h0*w_hat and Theta=h0*theta_hat retain the original
  continuous analytic pair and one common normalization. A single verified
  existing numerical profile represents U_star,C_star, the sum of constant
  and second-harmonic leading responses; it is not continuum truth.
- Initial u=epsilon_a^2 U_star, c=epsilon_a^2 C_star, w=epsilon_a W and
  theta=epsilon_a Theta+epsilon_a^3 Theta3 have zero velocities. The selected
  rule is w3=0 and the unique degree<=5 Hermite polynomial Theta3 in s/L,
  with endpoint value/first/second jets obtained from the frozen cubic action.
  This interpolation rule is neither a new boundary condition nor a claim
  of a unique physical preparation, nonlinear normal mode or periodic orbit.
- u,c remain independent during autonomous evolution. No forcing is added,
  no quasistatic c closure, phase/amplitude adjustment or energy matching.
  V0, quartic action, cubic residuals, material/kappa, basis and essential-only
  clamps are unchanged. Higher finite-amplitude residuals are retained.
- Preparation accuracy1e-6 is a separate numerical admission rule, declared
  before the audit; old trajectory1e-3/1e-4 and energy1e-6 gates are unchanged.
  Failed projection of either participant forbids a short convergence run.
  No new physical applicability threshold or scientific direction is selected.

## Four-field planar free-motion pilot (NLSP)

- The [numerical pilot](weakly_nonlinear_planar_time_pilot.md) restricts the
  accepted quartic action to v=Phi=psi=0, retaining independent u,w,theta,c.
  One homogeneous straight G20 rod has u=w=theta=c=0 at both ends, without
  slope constraints, an internal joint or book_slope_clamp.
- Shen polynomials P_n-P_(n+2), maximum degrees16/24/32, provide p-1 numerical
  unknowns per field. Gauss nq=2p+1 integrates quartic products exactly.
  The numerical note uses d_h only for the full spatial coefficient vector;
  it does not denote the geometric beta or a rotation-vector chart.
  This is spatial approximation of the distributed action, not a physical
  two-mode model. No quasistatic c, no-shear or inextensibility constraint.
- Initial bending is one unchanged continuous analytic Timoshenko pair,
  with one common max|w_hat|=1 normalization, A/h=.05/.025 and zero velocities,
  u=c=0. No static correction; the horizon is5 fixed linear periods.
- Variable theta mass and both conjugate inertial terms are retained.
  Constant-M0 whitening is an invertible numerical basis transformation.
  Radau uses a verified analytic Jacobian and solves the actual mass system;
  no inverse expansion, damping, filtering or energy projection.
- Declared numerical convergence gates: w/theta and velocities1e-4,
  u/c and velocities1e-3 relative to their own nonzero characteristic scales,
  quartic discrete energy drift1e-6. Small-neighborhood safety and bounded
  wall-time policies are numerical stopping rules, not physical applicability.
- A planar trajectory cannot establish out-of-plane stability. The angular
  same-clamp out-of-plane reference remains UNAVAILABLE. No parameter maps,
  periodic continuation, Floquet, critical amplitude or spring/KV transfer;
  V0 and prior algebra stay unchanged, LONG closed, EB/RLB-KV paused.

## Seven-field spatial nonlinear continuation (NLSP)

- The [isolated canonical note](weakly_nonlinear_spatial_rod.md) fixes
  q=(u,w,v,Phi,psi,theta,c), B=[t,n,k], a=(Phi,-psi,theta), R=exp([a]x).
  All fields and normalized jets have amplitude degree one; geometry and
  scales stay fixed. No no-shear, inextensible or quasistatic-c constraint.
- V0 is an explicitly adopted project nonlinear reduced constitutive law,
  with constant effective coefficients of the original section, not a
  published nonlinear Jang/Yartsev system or a full 3D Hooke-law reduction.
  Its objective measures, positive reference energy, unstressed origin and
  accepted linear limits motivate the choice; they do not validate a real
  joint's nonlinear threshold. Additional c/strain-curvature elastic products
  are omitted by assumption, not because all have higher amplitude order.
- J(c) is obtained from the declared affine section mass displacement and
  original mass. No extra (1+c) volume/mass factor is applied. Constant elastic
  coefficients and c-dependent inertia are distinct adopted reductions.
  C_T retains the existing rectangular generalized torsion reduction, not GIp;
  no dynamic warping/bimoment or monoclinic Sbar16 is introduced.
- Quartic action gives cubic coordinate residuals and their conjugate fluxes.
  The rotational covector includes P^T Jr^T; it is not the raw body moment.
  Common translation, rotation-vector and c DOFs give 7+7 dual joint rows;
  c1=c2 remains a reduced variational closure, not finite-joint 3D elasticity.
  Outer U=z=c=0 clamps do not impose centerline slopes and differ from
  Yartsev's book_slope_clamp. New notation is local; frozen equations stay intact.
- Exact A/B and supplied-draft comparison, manufactured-jet O(epsilon_a^4)
  residual checks and bounded linear/split controls audit this mathematical
  model. Plane invariance is not plane stability. No nonlinear trajectory,
  periodic orbit, modal reduction, Floquet multiplier or critical amplitude
  was calculated. LONG remains closed; EB/RLB-KV remains paused, not transferred.

## Production Lambda(beta) implementation checks

- The [large geometry check](mindlin_herrmann_timoshenko_lambda_beta_large_checks.md)
  uses canonical Lambda^4=rho*A_ref*omega^2*l_ref^4/(E*I_ref), with fixed
  l_ref=.5,b_ref=.20,h_ref=.05,E=rho=1 for every curve. The RLB-2B local
  frequency label is canonical Lambda^2; its old formula/data are untouched.
- Length mu0/.25/.5 preserves L1+L2=1 and identical section. Descriptive
  thickness_contrast delta_h0/.2/.4 gives h1/h2=h0(1∓delta_h) at equal
  lengths, preserving total mass. This local delta_h is not another eta/tau.
- Both arms retain project_jang_reduced_rectangular, kappa5/6 and all
  boundary/joint rules. Two existing arm bases/Dirichlet maps are composed
  with the same joint operator; no new wave law is introduced. Straight
  stepped reference uses literal positive-X state transmission.
-37 fixed beta values0:2.5:90 have independently sorted spectra. Prior
  frequencies seed brackets only; count/quality gates apply at each point.
  Same-system swap controls reflect the angle bisector: c/R invariant,
  w/theta/Q/M change sign. No across-beta forms/branch assignment.
- All map/reference/swap/root gates PASS. Curvature flags are numerical
  diagnostics independently checked at the same points, not physics or
  smoothness thresholds. c1=c2 remains reduced variational closure rather
  than finite-joint3D elasticity; no applicability/sensitivity/veering claim.

## Bounded homogeneous thickness screening

- The [2026-10-07 thickness diagnostic](coupled_longitudinal_theory_thickness_screening.md)
  changes h=h0*s_h at s_h=1/1.25/1.5/1.75/2, beta0/45/90 only; b=.20,
  h0=.05,E=rho=1,nu=.3,L1=L2=.5,kappa5/6 fixed. No mass preservation;
  m grows with h, I/J/j/r with h^3. Same Tim coefficients at each h for A/B/C.
- q_h=I/(A L_arm^2) and lambda_h=L_arm/h are derived reporting quantities;
  lambda_h is not the project's frequency Lambda. Normalized refined
  coefficient ratios contain h^2, while static MH layer length sqrt(H/EA)
  scales as h. Neither establishes a spectral power law. Common b cancels
  in this homogeneous planar 1D model, not necessarily in 3D elasticity.
- The unchanged finite basis supports sub-optical windows. Certificates
  cover first12+guard13 inside them; unsupported ceiling tails are not
  asserted complete. Same accepted residual/conditioning tolerances.
  Full beta0 direct-profile/interface checks pass at all five thicknesses.
- Increasing h reduces slenderness to L_arm/h=5,b/h=2. Only comparison
  among adopted 1D theories is claimed, not physical 2D/3D validation.
  MH clamp/contraction closure qualifications remain unchanged. Fixed-case
  overlaps and D_c give diagnostics; no energy classes, across-h/beta
  tracking, root reordering, threshold search or range extension is implied.

## Bounded coupled longitudinal-theory screening

- The [2026-10-07 screening](coupled_longitudinal_theory_hierarchy_screening.md)
  compares elementary, planar Rayleigh–Love and unchanged Jang-project MH
  for one equal-arm G20 rectangle on eight declared angles. Same material,
  Timoshenko basis, kappa=5/6 and physical centroid/bending clamps.
- Planar Love uses J=nu^2 rho Iy, H=0 and the variational harmonic force
  N=(EA-J omega^2)u_x. No polar Ip, independent c/R or u_x clamp is appended.
  A/B omit the resolved contraction clamp/joint constraint of MH; this
  boundary limitation is retained. MH c-continuity remains reduced-frame
  closure, not finite-joint 3D elasticity.
- Each theory/angle is independently sorted. MH is a reference model;
  fixed-beta geometric displacement/rotation overlaps only diagnose
  correspondence and never reorder roots. D_c compares c with−nu u_x,
  with machine-scale null-field handling. No modal energy labels, across-beta
  continuation, applicability thresholds or extrapolation are authorized.

## General-angle reduced M-H--Timoshenko joint

- The [general stage §§9--16](mindlin_herrmann_timoshenko_rigid_joint.md#9-general-angle-geometry-действующий-project-contract)
  extends the same reduced point-joint closure and unchanged production arms.
  Geometry is imported from the project physical t/n contract: positive beta
  turns the joint-to-right-clamp ray upward, both local x point clamp-to-joint.
  Proper frame rotations leave c/theta scalar DOFs unchanged; displacement
  and nodal-force transforms are dual. Endpoint signs are still minus at0,
  plus atL; all eight invariant joint rows have rank8.
- Under an improper EX reflection, keep the same signed planar rotation axis:
  t*=S t,n*=-S n, (u,c,w,theta,N,R,Q,M)*=(u,c,-w,-theta,N,R,-Q,-M).
  Thus c/R are scalar but theta/M are pseudoscalar for physical reflection.
  This is not a sign fit or a change of the local energy/PDE convention.
- MHTIM_GENERAL_BETA_JOINT_GATE=PASS: frozen beta0 regression, small-angle
  limits, rank/duality/right angle, same-geometry swap/reflection and one
  equal-arm G20 pilot5/45/90 pass. First12+guard13 are covered by exact
  energy/Schur count using certified fixed-arm poles. No beta map, tracking,
  applicability hierarchy or new coefficient prescription is adopted.
- c1=c2 remains published variational reduced-frame closure, not a direct
  finite3D welded-joint derivation. The historical beta0/single/source
  qualifications below remain attached to their original stages.

## Reduced M-H--Timoshenko rigid-joint beta0 gate

- The [joint note](mindlin_herrmann_timoshenko_rigid_joint.md) adopts common
  translations, theta and c from Rucka's reduced common-DOF frame assembly.
  Additional c compatibility is variational 1D closure, not a direct 3D
  elasticity derivation for the finite welded joint. The massless rigid point
  contributes no energy/inertia; all eight balances follow from endpoint work.
- Arm/source physics and the production Jang preset remain unchanged.
  Both positive local coordinates run outer clamp to joint; at beta0,
  t1=(1,0),n1=(0,-1),t2=(-1,0),n2=(0,1),theta rotation about k=-EZ.
  Outward endpoint sign is - at local0, + at localL; both joint signs are +.
  Hence c1=c2 and R1+R2=0; reflection to global X makes R continuous.
- Only homogeneous G20,total L=1,splits .5/.35/.65 are verified: energy-domain
  transparency, saturated min-max inventories7 MH/11 Timo, combined12+guard13,
  mass-overlap/component L2 including c, interface work and arm swap pass.
  MHTIM_BETA0_JOINT_GATE=PASS; beta!=0 spectrum/closure is not validated.
  No fitting, sweeps, nonlinear/FEM work or promotion of Ng coefficients.

## Isolated Mindlin--Herrmann / Timoshenko source audit

- Subsequent explicit project decision: PRODUCTION_MHTIM_FORMULATION_SELECTED
  selects JANG_BARE_ISOTROPIC_REDUCED_MH_TIMOSHENKO, independently of the
  historical source statuses below. project_jang_reduced_rectangular uses
  the accepted G20 rectangular K=5/6 in both H=K*GI and S=K*GA, j=r=rhoI,
  C=EA/(1-nu^2); PRODUCTION_MHTIM_KAPPA_RESOLVED. Numeric source Jang kappa
  remains unstated, source CLI still requires explicit input. Fernandes's
  newly available (19)--(21) confirms Ng-factor reuse only; neither the
  factors 12/pi^2 nor its Lame normal block are adopted for production.
- Finite control is one normalized G20 rod, E=rho=1, nu=.3, b=.2,h=.05,
  L=L_ref=1. Source fields and Jang (12) support u=c=w=theta=0 at both
  ends as full clamp within the planar four-field theory. c is independent,
  so c=0 does not impose u_x=0. That single-rod stage chose no joint closure;
  the subsequent qualified reduced closure/beta0 gate is recorded above.
  MHTIM_SINGLE_ROD_FINITE_SPECTRUM=PASS. Elementary and planar Rayleigh--Love
  use only u=0 per end and unchanged bending clamps; planar RL J=nu^2*rhoI,
  H=0 exactly. Their unresolved contraction boundary layer qualifies
  HIERARCHY_SINGLE_ROD=PARTIAL_PASS despite all numerical gates passing.
- Root completeness in this fixed case uses a Young lower quadratic form
  for MH and exact simply-supported lower spectrum for Timoshenko. Saturated
  min-max count bounds (7 MH,11 bending, including guards), not sign scans
  alone, certify the bounded inventories. Above-cutoff finite branches and
  any new geometry/parameter map are outside this gate.
- Rectangular prescription gate (2026-10-06): Ng PDF 10 (2) directly maps
  S1=12/pi^2 to H=S1*GI and S2=S1*((1+nu)/(.87+1.12nu))^2 to j=S2*rhoI.
  These are printed decimal constants, not exact rational physical constants.
  However its printed 3D Lame normal block D=(2mu+lambda)A,F=lambda*A is
  different from the current C,nu*C reduced block at the same E,nu.
  The second new PDF is a circular SSRN review; the expected Fernandes
  full text was not found. RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED and
  PRODUCTION_MH_COEFFICIENTS_UNRESOLVED are retained. The named fixture
  candidate rectangular_literature_default is unadopted; no helper/default,
  constitutive change or transfer assumption is introduced. Source Rucka/
  Jang, unstated Jang kappa, Bishop and project Timoshenko kappa are preserved.
- The [one-rectangle source audit](mindlin_herrmann_timoshenko_single_rod.md)
  uses planar, homogeneous isotropic, unstressed linear source energies.
  Local q=(u,c,w,theta); c is independent dimensionless contraction along z,
  theta is the legacy-sign Timoshenko rotation. b along y, h along z;
  source I is I_y=b*h^3/12, not I_p. This is not full two-direction lateral
  contraction or a new literal common 3D displacement/Hooke law.
- Jang's source closure sums axial/M-H and bending/Timoshenko reduced
  stress-work contributions. Axial suppresses sigma_yy,tau_xy,tau_yz;
  bending additionally suppresses sigma_zz and returns Q11=E. Corrected
  energy (5) / resultants (A9) fix the implemented shear prescription;
  uncorrected (A7) and printed (A3),(A10) are not copied literally.
- Rucka K_MH1=1.1, K_MH2=2.1, K_Tim1=.95 and K_Tim2=12*.95/pi^2 are
  source-specific fitted/selected factors, not project defaults. Jang uses
  K_MH1=K_Tim1=kappa_b, K_MH2=K_Tim2=1; numeric kappa_b is unstated in the
  audited text. The explicit 5/6 run is a conditional control only.
- Centered homogeneous geometry eliminates kinetic first-moment cross
  terms; local potential separation follows the source constitutive blocks.
  Acoustic/optical dispersion labels are continued from k=0, not descendant
  mode identities. Figure comparisons are qualitative without digital data.
- `MHTIM_VARIANT_DEPENDENT`, `MH_SOURCE_VARIANTS_NOT_EQUIVALENT` for the
  published prescriptions, and `PRODUCTION_MH_COEFFICIENTS_UNRESOLVED` are
  distinct statuses. Production promotion, finite-rod extra BC and angular
  joint conditions are not selected. The closed Bishop gate is preserved.

## Isolated single-rod Timoshenko--Bishop kinematic audit

- The [single-rod audit](timoshenko_bishop_single_rod.md) is restricted to
  a straight, unstressed, homogeneous isotropic centered rectangle, linear
  motion and one bending plane. Local legacy signs are Ux=u-z*psi,
  Q=kappa*G*A*(w'-psi); b is along y, h along z, I_y=b*h^3/12.
- Centroid-only Poisson contraction and a full-strain compatible transverse
  distortion are candidate fields, not accepted combined models. Exact
  axial--bending cross terms vanish for both centered fields, but the raw
  first field gives C11*I_y and G*A instead of E*I_y and kappa*G*A; the
  second adds bending gradient energy/inertia. A relaxed hybrid preserving
  the original self terms additionally needs an explicit reduction and a
  reflection-preserving shear closure. Neither is silently adopted.
- Scientific status: COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED. This finite
  audit stops before combined spectra, hierarchy acceptance or boundary
  selection. It does not establish physical linear coupling, select angular
  joint conditions, alter the circular literature controls, or replace the
  frozen equations/API. The Lamé coefficient lambda_L and J_B,H_B are local
  notation only; no project-wide renaming is implied.

## Isolated longitudinal Rayleigh–Bishop literature diagnostics

- The [literature reproduction](bishop_literature_reproduction.md) uses linear,
  local, isotropic circular rods only, with the exact Marais two-section and
  Popov–Sadovsky uniform-rod geometry. No angled-joint conditions, rectangular
  section, bending coupling, damping, FEM or nonlinear extension is implied.
- Local coefficients are m=rho*A, J=nu²*rho*Ip, H=nu²*G*Ip; Ip is the polar
  area moment. Source Marais eta=nu, mu=G and lambda=omega² are not project
  geometric eta/mu or Lambda. N=-Gamma' y; state order (U,U',N,P).
  C (U=U'=0), F (N=P=0), and UP (U=P=0) are distinct physical boundaries.
  H=0 is a separate second-order problem, with only U=0 or N=0 per end.
- Marais's thick/short circular sections are retained to reproduce its rod
  theory, without a 3D accuracy claim. In Popov, division by rho*A permits
  using c,nu,d alone; no E/rho values are invented. The same printed c is
  used for all models. c and fitted nu are experiment-derived, and the 30
  rounded ratios do not recover raw experimental precision.
- The homogeneous sin(n*pi*x/L) limit is checked against the boundary matrix;
  the elementary H=J=0 limit agrees with the longitudinal cos/sin block of
  equations.tex and analytic/formulas.py. Existing baseline formulas are frozen.

## Literature-only KV benchmarks

- Failla 2014 example 6.1 and Hong–Kim 1999 example 1 use their own straight-beam
  boundaries, material inputs and reference time. Translational supports exist
  only in these benchmark assemblies; production angled-joint physics is unchanged.
- Source variables and signed transformations are defined in the
  [transcription](../laminated_beams/inplane_kelvin_voigt_literature_benchmarks.md#source-transcription-and-mapping).
  Failla's deflection `psi`, moment `mu` and eigenvalue `omega` are not project
  section rotation, length mismatch or real elastic frequency. Hong's E,G,nu
  remain the independent printed inputs, without enforcing another G relation.

## Bounded in-plane Kelvin–Voigt joint pilot

- The arms remain linearly elastic, without distributed damping. One
  nonnegative rotational dashpot acts in parallel with the existing massless
  spring on relative section rotation, retaining translation compatibility
  and external clamps. RLB uses section rotation, not the centreline slope.
- The local convention is `exp(p*t)`, `p=-alpha+i*omega_d`, with fixed
  reference time and `z=p*t_ref`; old elastic symbols, state order and
  coordinate signs remain unchanged.
- This pilot uses two isolated descendants, identical H/L/L/H arms,
  beta0=5°, kappa_theta=1 and three small positive damping values. See the
  [imported theory](../laminated_beams/inplane_kelvin_voigt_joint_theory.tex)
  and [numerical scope](../laminated_beams/inplane_kelvin_voigt_pilot.md).

## Chapter-2 anisotropic rigid-joint pilot

- The first two-monoclinic-rod diagnostic uses an ideal point joint with no
  mass, rotary inertia, eccentricity, elastic compliance, extra joint torsion
  stiffness, independent warping coordinate, or bimoment. Generalized torsion
  remains condensed in the existing Chapter-2 `C_T`. Both external ends use
  the source-confirmed book slope clamp. These assumptions define only the
  small elastic pilot documented in
  `docs/anisotropic_rods/yartsev_ch2_rigid_angular_joint.md`; they do not
  define a stable or final coupled-rod model.
- Its canonical state is the book state
  `[w_i, psi_i, Phi_i, Q_i, M_i, M_{T,i}]^T`. The project bases remain
  `e_z`, `t_i`, `n_i=e_z x t_i`, and the translation to the old out-of-plane
  notation is isolated in
  `docs/anisotropic_rods/yartsev_ch2_notation_translation.md`.

## Chapter-2 rectangular orthotropic EB comparator

- The comparator is restricted to the real-elastic HMS/DX-209 endpoint
  `theta_1=theta_2=0`, where `Sbar16=0`. It retains the book state and the
  existing ideal point-joint conditions, removes transverse shear and bending
  rotary inertia, and retains torsional inertia `rho I_p`. It introduces no
  warping coordinate or bimoment.
- Rectangular torsion uses the existing generalized stiffness
  `C_SV=Cbar=C_T` at `theta=0`; `G I_p` is not a valid replacement. The
  independent validation FEM uses Hermite EB bending and linear Saint-Venant
  torsion with consistent mass. These assumptions define only the finite
  `PARTIAL_PASS` gate documented in
  `docs/anisotropic_rods/yartsev_ch2_rectangular_eb_validation.md`, not a
  general anisotropic EB model or production API.

## Out-of-plane EB + Saint-Venant torsion

- The out-of-plane subsystem is treated as a separate linear model for a
  planar two-beam system: Euler--Bernoulli out-of-plane bending plus
  Saint-Venant torsion. It is independent of the existing in-plane
  bending/axial subsystem in the ideal linear planar setting. Its coordinate
  convention is fixed in `docs/theory/out_of_plane_eb_torsion.md`: `e_z`
  points downward, `t_i` runs from clamp to joint, `n_i=e_z x t_i`, and
  tangent vectors must be written as `t_i` rather than `tau_i`.

## Diagnostic thickness mismatch

- The mass-preserving thickness-mismatch model with
  `eta=(r_2-r_1)/(r_1+r_2)` is diagnostic-only. It keeps the total length
  `2l`, the equal-radius base radius `r_0` at `eta=0`, and the total mass fixed
  through `tau_1=(1-eta)/sqrt(1+2 mu eta+eta^2)` and
  `tau_2=(1+eta)/sqrt(1+2 mu eta+eta^2)`. It does not replace the baseline
  equal-radius determinant, article model, or FEM model. See
  `docs/thickness_mismatch/README.md`.

- Базовые обозначения для theory-facing формул: `l_1`, `l_2` — длины плеч, `l=(l_1+l_2)/2` — базовая длина, `r` — радиус круглого сечения, `\Lambda` — безразмерная частота, `\beta` — угол сопряжения, `\varepsilon` — параметр толщины, как в `docs/literature/pdf/Статья-Дорофеев-2025.pdf`.
- Локальный параметр `\mu` используется как параметр несимметрии длины в текущем расширении задачи: `\mu=(l_2-l_1)/(l_1+l_2)`, `l_1=l(1-\mu)`, `l_2=l(1+\mu)`. Эквивалентно, в формулах с базовой длиной `L`: `L_1/L = 1-\mu`, `L_2/L = 1+\mu`. В опубликованной статье Дорофеева 2025 этого параметра нет, потому что там рассматриваются одинаковые стержни длины `L`.
- Случай `\mu=0` фиксирует симметричный эталон для рабочей нумерации ветвей. Нумерация относится к отслеживаемой ветви, а не к текущему месту этой ветви в отсортированном спектре.
- Пунктирные reference-линии на графиках `\Lambda(\mu)` соответствуют одиночным стержням с закреплением заделка--шарнир: `\Lambda_n^{(1)}(\mu)=\alpha_n/(1-\mu)` и `\Lambda_n^{(2)}(\mu)=\alpha_n/(1+\mu)`, где `\alpha_n` — корни `\tan\alpha=\tanh\alpha`. Эти линии используются только как ориентиры для интерпретации, а не как замена анализу форм колебаний.
- При использовании `docs/literature/pdf/2003JSVb.pdf` знаки в матрице `T` формулы `(39)` нужно сверять вручную с уравнениями непрерывности `(10)--(11)`; печатную запись этого определителя не стоит переносить в теорию или код без дополнительной проверки.
