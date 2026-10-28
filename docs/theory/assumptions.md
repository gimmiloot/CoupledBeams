# Assumptions

For project-wide rules on source-of-truth priority, branch identity,
diagnostic workflow, and model-extension checks, see `../project_rules.md`.

Здесь фиксируются активные допущения и рабочие гипотезы проекта.

## Working Notes

## LONG: one horizontal external spring, 2026-10-26

- [Canonical construction and verification](isotropic_two_springs_horizontal_spring.md):
  unchanged isotropic RL_FULL/Timoshenko rods and two exact nodal translations;
  first outer clamp, internal rotational spring, one global horizontal spring
  at the second outer end. No outer guide or rotational spring.
- Each angle has its own unstressed equilibrium and a horizontal spring at
  natural length. The linearized spring Hessian is rank one; no transverse
  geometric stiffness, spring inertia, damping or preload is included.
- Exact horizontal RIGID constrains only `cos(beta)*u2+sin(beta)*w2`.
  Physical sorted positions include exact zero mechanisms. Static nullity
  is checked before positive roots; no artificial regularization is allowed.
- Local beta-radian asymptotics require an isolated simple eigenvalue and
  a nonresonant axial problem. The finite-joint coefficient uses the TOTAL
  angular derivative of outer horizontal displacement, including axial response.
  Bounded verification does not authorize parameter maps or geometry optimization.

## NSL-3A-ANGLE-C: two complex targets and modal envelopes, 2026-10-20

- [The targeted confirmation](../nonsymmetric_laminated_beams/nsl3a_angle_complex_confirmation.md)
  uses the unchanged passive NSL-3A pencil KJ=k_theta+c_theta*p at beta15,
  kappa1,d=.001, only MIRROR5/SAME_LOCAL5 from the validated ANGLE seeds.
  No constitutive, internal damping, support or joint assumption changes.
- Figure1 uses elastic predictors only; anglewise maxima among positions1–6
  are not one tracked branch. Figure2 uses saved mass-normalized elastic
  profiles with signed mass-overlap alignment and shared ordinate ranges.
- Figure3 is the envelope of ONE complex eigenmode in each construction:
  exp(-2*pi*alpha*N/omega_d), N=omega_d*t/(2*pi). Equal N is not equal time.
  Accumulated separation over200 cycles is not evidence of strong per-cycle
  damping, a forced response or behavior with other loss mechanisms.
- Three new exploratory figures are not approved article illustrations.
  Previously completed stages, their qualifications and bundles stay unchanged.

## NSL-3A-ANGLE: elastic angle screening, 2026-10-20

- [The separate screening](../nonsymmetric_laminated_beams/nsl3a_elastic_angle_screening.md)
  keeps the NSL-2 PROJECT_REDDY_COUPLED model, coefficients, clamps and point
  spring unchanged; c_theta=0 in every spectral operator, kappa_theta=1.
- Only beta15/30/45/60/75 and MIRROR/SAME_LOCAL are compared. The first six
  sorted positions are observations, not branch labels. The seventh is a
  completeness guard. Beta30 forms/frequencies are reused unchanged.
- P, s, G, alpha_pred and zeta_pred are elastic energy indicators.
  Illustrative d_theta=.001 enters postprocessing only; no new damped
  eigenvalue or finite-viscosity accuracy is inferred.
- Exact zeros require MIRROR reflection class eta=-1. Small SAME_LOCAL
  values are not rounded to structural zeros. Comparisons at different
  angles use global physical displacements/rotations at common arm coordinates;
  the diagnostic common mass metric is justified here by equal densities,
  equal m,j2 and j1=0. MAC does not imply continuous branch continuation.
- Maxima are restricted to the fixed grid/first six positions. No physical
  assumption or historical NSL/RLB/LONG/NLSP status is revised.

## NSL-3A: one passive rotational Kelvin–Voigt element, 2026-10-20

- [NSL-3A note](../nonsymmetric_laminated_beams/nsl3a_kelvin_voigt_pilot.md)
  adds only R_J=c_theta*(Delta psi_dot)^2/2 to the accepted NSL-2 model.
  The node remains massless, concentric and translationally rigid; rods
  remain linearly elastic PROJECT_REDDY_COUPLED FSDT with full B_e,j1.
- The time convention is exp(p*t), K_J=k_theta+c_theta*p, including the
  real part of p. No material or distributed damping is added.
  Rayleigh dissipation here denotes a single joint loss function, not
  proportional mass/stiffness damping.
- Complete Hermitian modal energies give p^2*M+p*C+K=0 and exactly
  alpha=C/(2*M) for an oscillatory eigenpair. Elastic-form prediction is
  only first order in c. omega_d=Im(p), zeta=alpha/abs(p).
- MIRROR eta=-1 is structurally protected for both k and c; SAME_LOCAL
  cannot use that reduction. Only MIRROR4/6 and the active MIRROR5,
  SAME_LOCAL4/6 are evaluated, at beta30,kappa1,d=.001 and saved NSL-2 units.
  NSL-2R column equilibration is required; fixed scales preserve analyticity.
- This limited diagnostic is PASS. It does not model a finite physical
  joint or authorize NSL-3B, other viscosities, maps, veering, FEM or FRF.
  Historical assumptions and qualifications below remain unchanged.

## NSL-2: separate elastic angular joint, 2026-10-20

- [NSL-2 theory](../nonsymmetric_laminated_beams/coupled_joint_theory.md)
  extends the accepted PROJECT_REDDY_COUPLED NSL-1 state only. Both arm
  coordinates run from section-rotation clamps to the massless point joint.
  RLB physical maps d=u*t+w*n, rotation=-psi*k, F=N*t+Q*n, moment=-M*k
  are retained. The local NSL material z in Ux=u+z*psi is embedded along n;
  stack order is increasing local z. The source-coordinate conversion
  z_P=-z_R is not applied a second time to accepted NSL coefficients.
- RIGID imposes common translation/rotation and force/moment equilibrium.
  SPRING keeps common translation and independent rotations with
  U_J=k_theta*(psi1-psi2)^2/2, no extra mass, eccentricity, prestress,
  translational spring or damping. k_theta=0 is the hinge; RIGID is exact.
- Reference-axis shifts are representation changes only: u*=u+z0*psi,
  M*=M-z0*N and the transformed B,D,j1,j2 are accompanied by the
  pullback of the full joint and clamp rows to the ORIGINAL physical point.
  Neutral-axis B*=0 does not authorize dropping the induced j1*.
- For equal arms, physical reflection exchanges states by Tr=diag(1,-1,-1,1,-1,-1).
  MIRROR pairs reflected stacks; SAME_LOCAL does not in general preserve
  reflection. Eta=-1 has Delta_psi=0 and is spring-independent.
  A continuous straight HHLL rod requires local HHLL/LLHH on opposite axes.
- Only the fixed H/L pilot at beta30,kappa1 and prescribed hinge/rigid
  and beta0 controls are performed. D_ref=.20*.05^3/12,l_ref=1,m_ref=.01
  remain common. beta90 is structural only. No Kelvin–Voigt NSL-3 work.
- Overall PARTIAL_PASS preserves strict endpoint failures in the straight
  control, despite matching frequencies. Main pilot prefixes pass.
  Earlier NSL-1, RLB, LONG and NLSP assumptions/statuses are not revised.

## NSL-1: separate nonsymmetric single-beam FSDT direction, 2026-10-20

- [NSL theory](../nonsymmetric_laminated_beams/single_beam_theory.md)
  is a new direction, not a change to symmetric RLB or the MH/Timoshenko
  line. It uses Ux=u+z*psi, Uz=w, gamma_xz=w'+psi; state
  (u,w,psi,N,Q,M), project z increasing bottom-to-top. The Reddy source
  coordinate is z_R=-z_project; Bachoo Fig.1 uses the upward coordinate
  but lists its [0/90]_n stack top-to-bottom.
- SOURCE_BACHOO directly integrates b*Qbar11 times 1,z,z² and K*b*Qbar55.
  PROJECT_REDDY_COUPLED retains ABD strain coordinates (xx membrane,
  xx curvature) and jointly eliminates the other four strains under zero
  Ny,Nxy,My,Mxy. No independent inertia is assigned to those eliminated
  strains; kinetic energy is the declared three-field FSDT energy with
  m=b*I0, j1=b*I1, j2=b*I2. This is a chosen reduced model, not a unique
  deduction from 3D elasticity. Source and project coefficients are labeled.
- Only coaxial/cross-ply 0/90 layers, perfect bonds, initially straight
  unstressed constant sections and positive constitutive/mass forms are
  accepted. Off-axis spatial couplings are explicitly rejected. K is an
  explicit input, applied once; 5/6 is printed by Bachoo and inherited
  from the previous H/L diagnostic for the project control.
- NSL-1 has one CC single beam, u=w=psi=0 at x=0,L, with no angular joint,
  spring, damping, forced response or FEM. The only spectral case is
  PROJECT_REDDY_COUPLED H/H/L/L, chi=.4,L=1,b=.20,h=.05.
- Joint Schur reduction, a future KV law and future frame results are not
  attributed to Bachoo. Printed Eq.(11d) differs from the determinant of
  Eq.(9) by the missing density power; the determinant governs the local
  model. Eq.(4b) also conflicts with orthotropic reciprocity, although its
  Q12 contribution vanishes in Qbar11 for the admitted 0/90 plies.
  See [source comparison](../nonsymmetric_laminated_beams/source_comparison.md).
- B_e,j1 are local NSL coupling notation. At B_e=j1=0 the project operator
  recovers the unchanged symmetric RLB operator exactly. Historical
  equations.tex and src/my_project/analytic/ have no NSL block and remain
  unchanged. Earlier RLB, LONG, NLSP statuses retain their original scopes.

## Separate linear one-axis spatial KV joint, 2026-10-19

- The explicitly authorized [diagnostic](uniaxial_kelvin_voigt_split.md)
  extends the planar MH-FREE/KV idealization to initially planar, unstrained,
  isotropic rectangular rods with aligned principal section axes. It uses
  only the quadratic part of the accepted seven-field action. This is a new
  stated joint idealization, not an inferred property of a real 3D joint.
- The full center displacement is continuous. For the physical rotation
  `vartheta=Phi*t-psi*n+theta*k`, both relative components in the initial
  plane are constrained, while only `theta1-theta2` carries the existing
  spring/damper. The joint is massless and has no eccentricity, translational
  compliance or other damping. Contraction uses only `R1=R2=0`.
- In this diagnostic `P=I-k*k.T` is the planar projector; the historical
  named-rotation sign matrix is called `D_z=diag(1,-1,1)` locally. This
  explicit distinction does not rename frozen historical notation/code.
  Physical moments are `MT*t-Mpsi*n+Mtheta*k`, as fixed by virtual work.
- External clamps prescribe independent section rotations, in particular
  `v=psi=Phi=0`, not `book_slope_clamp`. Rectangular torsion retains the
  existing generalized C_T, not G*I_p. The same theory and BC are used on
  both sides of the rigid/KV comparison.
- The exact linear OUT independence relies on these assumptions and the
  absence of mixed quadratic terms. It makes no nonlinear, real-joint or
  3D-accuracy claim. Historical rigid `c1=c2`, NLSP and EB/RLB assumptions
  remain unchanged. Only matrices at beta0/45/90 and one beta45 spectrum
  were authorized; this completed diagnostic does not open a wider study.

## Fixed longitudinal/Kelvin–Voigt comparison, 2026-10-19

- The separately authorized [comparison](longitudinal_kelvin_voigt_comparison.md)
  uses only G20 beta45, kappa_theta1, d_theta.001 and the existing physical
  k_theta,c_theta,l_ref,t_ref. Elementary, planar RL and MH-FREE share the
  unchanged Tim block, translation/force conditions and Jang rotation signs.
  External MH c=0 and reduced-model absence of c remain different boundary
  representations; MH is not 3D ground truth. No new constitutive assumption.
- At exp(p*t), the energy-derived RL natural force is N=(EA+J*p^2)*U',
  J=nu^2*rho*Iy. Each theory retains its own full modal mass; common u,w,theta
  projection is solely a fixed-case geometric correspondence diagnostic.
  The algebraic factor M_common/M_full does not introduce a new shape or
  reference normalization. Complex zeta means alpha/abs(p); its elastic
  predictor is alpha_pred/omega0.
- Only saved positions1–6 enter elastic predictions; only positions1/6 enter
  local complex comparison. Accepted elastic/complex MH bundles are reused,
  not recomputed. New E/RL local matrices preserve the state and analytic
  characteristic zeros with fixed nonsingular numerical scales and pivots.
  Positive finite k_theta,c_theta represent a passive point rotational element.
  No other c closure, geometry, viscosity, spectral search or parameter map.

## Separate contraction-closure sensitivity, 2026-10-19

- Explicitly authorized [bounded diagnostic](contraction_joint_sensitivity.md):
  unchanged Jang-project MH/Tim arms and full external u=c=w=theta=0 clamps;
  massless point joint with common centroid displacement, independent rotations
  and U_J=k_theta(theta1-theta2)^2/2. G20 equal halves, beta0/45,
  kappa_theta=1 only. Old rigid production closure is preserved.
- Alternative c trace spaces are FREE (R1=R2=0), CONTINUOUS (c1=c2,
  R1+R2=0), BLOCKED (c1=c2=0, reactive R allowed). All retain the other
  six scalar conditions. R=Hc' is conjugate to contraction, not ordinary
  bending moment. Its effective shear-stress first-moment interpretation
  does not assert a pointwise 3D stress solution or necessary point-joint transfer.
- U_c=k_c(c1-c2)^2/2 is described analytically only: zero/infinite k_c
  gives FREE/CONTINUOUS, never BLOCKED without a common-mode constraint.
  No k_c sweep. Strict potential positivity follows from positive bulk
  coefficients and external clamps; min-max uses nested MH domains.
- Elementary/planar RL comparators retain the same Tim block and spring,
  but omit independent c and its external clamp. Their differences from
  MH include boundary representations; MH is not an exact reference theory.
- Weak KV uses exp(p t), k_theta+c_theta p and full MH mass including j c^2;
  t_ref=l_ref^2 sqrt(rho A/B_ref), l_ref=.5 and d_theta=.001. Only two
  accepted fixed-case modal correspondences at beta45 may receive complex
  confirmation after elastic gates. No production EB/RLB changes or new maps.


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

## NSL-VIS-1: scope of the finite-dimensional spectral maps

The [NSL-VIS-1 maps](../nonsymmetric_laminated_beams/comparative_spectral_curves.md)
use the unchanged PROJECT_REDDY_COUPLED physical model, homogeneous BASE
and the HHLL/LLHH family with common A,D,S,m,j2. BASE is not H/L/L/H.
For these figures the explicit local frequency coordinate is
\(\Lambda=\sqrt{\Omega}=\sqrt{\omega t_{\rm ref}}\), with the common NSL-2
reference values. This does not redefine another model's determinant parameter.
Finite spring stiffness adds only its elastic quadratic form; RIGID imposes
equal section rotations as a separate essential constraint.

`frequency-map-v1 / fast_plot / sorted_positions` is implemented by the
existing Ritz method at N36/N44. Seven finite-dimensional eigenvalues are
retained; the guard is not a proof of continuum completeness. Saved analytic
references and two-order convergence support the research figures within
their declared accuracy. Lines through calculated points do not establish
modal branch identity, crossing or veering. Modal P/G are elastic quantities,
with exact eta=-1 zeros distinguished from numerical residues; no positive
viscosity or forced response is computed in this stage.

## NLSP: separate six-field linear candidate (2026-10-21)

The explicit user-defined [six/seven-field comparison](../numerics/nlsp_six_vs_seven_linear_comparison.md)
uses q₆=(u,w,v,Φ,ψ,θ), all at ordinary first amplitude order. Its axial
stiffness is EA₀, with physical ν=.3 retained in G,S,C_T. There is no
independent c field, its inertia or gradient energy, and no static/dynamic
condensation of the seven-field model. This is a separate linear idealization,
not the c=0 restriction with axial coefficient C=EA₀/(1−ν²).

The fixed isotropic rectangular arms, section-rotation clamps, coordinate
orientation and ideal rigid massless joint for common global translations
and rotations remain those of NLSP. Removing c also removes its external
and joint conditions; it adds no condition on u_s. Only J(0) is used.
No nonlinear law or nonlinear-stability transfer is assumed. The preserved
seven-field model remains a comparison model, not a three-dimensional exact
reference. This assumption is local to NLSP and changes no baseline,
LONG, NSL/RLB, anisotropic or viscoelastic-joint contract.

## NLSP: six-field quartic action and planar orbit (2026-10-21)

The subsequent explicit user contract extends the separate six-field candidate
through the [quartic action](../numerics/nlsp_six_vs_seven_planar_orbit.md).
Use the same exponential coordinates, Gamma and body curvature/angular
velocity as NLSP, with D6=diag(EA0,S,S) and fixed
J0=diag(j_parallel+j_perp,j_perp,j_parallel). All six fields, including u,
remain first amplitude order. Integrate and vary the action truncated at total
amplitude degree four; the production equations have degree at most three.
Geometry, section properties, rigid joint and clamps remain those of the
linear six-field contract above. No deformed-section stiffness law is added.

This action equals the seven-field restriction c=0, with the complete
(C-EA0)*Gamma_1^2/2 correction through degrees two, three and four. Only six
fields are varied. Physical nu=.3 is retained in G,S and torsion. This is a
new prescribed idealization, not exact elimination of the seventh equation.
At nu=0 only the restricted common-field actions coincide. Planar mass is
constant and planar nonlinear inertia vanishes; full spatial geometric
inertia remains even with fixed J0. These facts were checked independently.

The bounded p14/A=.005 calculation establishes a reversible Psi=0 orbit
with PL2 fundamental and PL4 third-harmonic participation. Its relation to
the historical Type-II orbit is supported locally by the physical seed,
phase and composition; no global branch continuation is asserted. Passing
HB11/HB13 tolerances neither proves nonlinear spatial convergence nor
transfers the historical seven-field OUT instability to the new model.

## 2026-10-22: isolated full Rayleigh–Love/Timoshenko two-spring verification

The separately authorized [two-spring diagnostic](isotropic_two_springs_rl_timoshenko.md)
uses homogeneous isotropic rectangles, small planar motion, fixed unstressed
geometry at each beta, common material/section, no damping, eccentricity,
prestress or joint mass. Both local coordinates run from outer supports to joint.
RL_FULL has J=rho*nu²*Ip, Ip=17I for b/h=4; the old planar J=rho*nu²*I
and Elementary J=0 are explicitly labelled comparison variants, with identical
Timoshenko bending. They do not replace the selected MH production model.
An external rotation spring is attached specifically to arm 2; finite-support
arm interchange requires moving the support types as well as the geometry.
RIGID is an exact constraint. Frequency Lambda is sqrt(omega*t_ref), with
l=(L1+L2)/2. The accepted range is the bounded low-frequency verification;
no 3D validation, plateau-angle rule or two-stiffness map follows automatically.

## LONG: направленное продольно-податливое сопряжение (2026-10-24)

Отдельная диагностическая ветвь, явно разрешённая заданием пользователя;
каноническая RL_FULL/Timoshenko-модель точного сопряжения не изменяется.
Локальные координаты и состояния прежние. В фиксированном направлении оси
первого стержня допускается относительное смещение c1=u1+cosβ*u2+sinβ*w2
с энергией k_parallel*c1²/2; поперечная связь c2=w1−sinβ*u2+cosβ*w2=0 точная.
Узел идеализирован как безмассовая направляющая с линейной продольной пружиной;
нет эксцентриситета, начального напряжения, вязкости или изменения масс.
χ=k_parallel*l/(EA)≥0; RIGID означает точное c1=0, не штрафную жёсткость.
При конечном χ освобождена ровно одна поступательная степень свободы;
обе конечные поворотные пружины сохранены. Каждая β — отдельная ненапряжённая
конфигурация. Направление направляющей фиксировано вдоль оси первого стержня,
поэтому перестановка стержней не является симметрией нового узла.
Варьируется только χ при ε=.005, μ=.01, b/h=4, ν=.3, k_s=5/6, l=E=ρ=1.
Сравниваются отсортированные частоты. При простой собственной паре и знаке
Kz+C2ᵀ*r2=ω²Mz производная на радиан равна
(k_parallel*c1*partial_beta(c1)+r2*partial_beta(c2))/(2V), r2=−Q1.
[Постановка, проверки и результаты](isotropic_two_springs_longitudinal_slip.md).

## Формальная малоугловая асимптотика RL_FULL + Timoshenko с точным узлом

2026-10-24, [LONG-D21](../memory/decisions.md#long-d21).
Это аналитический режим канонической двухпружинной модели, а не замена её
физического модуля и не ветвь продольно-податливого узла.
ε→0, β_rad=εξ при ограниченном ξ; фиксированы μ=.01, b/h=4, ν=.3, k_s=5/6,
равные сечения двух стержней и конечные положительные безразмерные κ_joint/end.
Исследуются низкие собственные значения ζ=Λ⁴=O(1); λ_i=1∓μ обозначают длины.
При масштабе u=aεU, w=aW, ψ=(a/l)Ψ ограниченная энергия обосновывает малость
продольного поля в норме, а не поточечное отношение u/w. Ведущая EB-кинематика
требует корректного вариационного предела сдвиговой формы; ограниченность
энергии сама по себе не даёт малости всей сдвиговой энергии.
В ведущей задаче продольные поля минимизируются и дают ξ²W2(λ2)²/4
в потенциальной энергии. Поперечная масса остаётся ведущей, продольная и
вращательная имеют порядок ε², добавочная RL_FULL-масса — ε⁴ при этих полях.
Простота пары и сходимость производных необходимы для переноса кривизны;
масштабирование β_pl дополнительно требует устойчивого порога, ненулевой
опорной чувствительности и отдельного контроля хвоста до 90°.
Ни равномерность при ξ→∞, ни универсальность подобия не предполагаются доказанными.
[Полный вывод и ограничения](isotropic_two_springs_angle_asymptotics.md).
Локальные обозначения ζ, ξ и L_Σ не меняют нотацию базового `equations.tex`.

## Разные сечения двухпружинной RL_FULL/Тимошенко-конструкции при постоянном объёме

2026-10-25, [LONG-D26](../memory/decisions.md#long-d26).
Отдельное расширение геометрии по явному заданию: μ,τ∈{−.1,0,.1},
τ=(h2−h1)/(h1+h2), b_i=4h_i, одинаковые E,ρ,ν, k_s=5/6.
При d=1+τ²+2μτ положены h1=h0(1−τ)/sqrt(d), h2=h0(1+τ)/sqrt(d).
Суммарные длина 2l, объём 2A0l и масса постоянны; сечения каждого стержня
постоянны вдоль его оси. Внешняя пружина остаётся на стержне 2, обе
поступательные связи точные. Коэффициенты EA,EI,k_sGA,ρA,ρI,J_RL=ρν²I_P
индивидуальны; I_P=17I. Размерные пружины и Λ определяются базовыми A0,I0,l,
ε0=.005, а не текущим сечением какого-либо стержня. Новых элементов узла нет.

Формальная ведущая задача относится к ε0→0 при ограниченном ξ=β_rad/ε0
и фиксированных μ,τ,κ. Её изгибная и массовая формы взвешены I_i/I0 и A_i/A0;
продольная минимизация даёт ξ²P/S_ax. Для простой пары C=P/(2S_ax ζ0 M0).
Проверены алгебраические тождества и численное решение; общая теорема
асимптотической сходимости неоднородной полной модели не предполагается.
Разложение производной по τ локально при τ=0 и не задаёт единственного
аддитивного разложения конечного изменения при |τ|=.1.
Обозначения a_i=A_i/A0 и j_i=I_i/I0 локальны для
[отчёта](isotropic_two_springs_thickness_asymmetry.md), где приведены энергии,
знаки, вывод коэффициентов и ограничения. Базовые `equations.tex` не изменены.

## Продольная пружина на внешнем конце второго стержня

2026-10-25, [LONG-D28](../memory/decisions.md#long-d28).
Отдельная постановка при прежних одинаковых сечениях, mu=.01, epsilon=.005,
b/h=4 и двух прежних конечных положительных поворотных пружинах.
Две поступательные связи узла точные. Новая энергия ku*u2(0)^2/2 заменяет
только продольное ограничение внешнего конца второго стержня; w2(0)=0
и его поворотная пружина сохраняются. Локальное обозначение chi_u=ku*l/(EA).
Для каждой ненапряжённой геометрии продольная пружина и направляющая
ориентированы вдоль оси второго стержня; фиксированное глобальное направление
не предполагается. Сила N=(EA−J_RL*omega²)u′ сохраняет инерционный член.
Естественное условие −N2(0)+ku*u2(0)=0 получено из энергии; RIGID означает
точное u2(0)=0, без штрафной жёсткости. Дополнительной массы опоры нет.

При фиксированной chi_u и ограниченном xi=beta_rad/epsilon ведущая
продольная податливость равна l/(EA)*(2+1/chi_u). При chi_u=0 ведущая
добавка нулевая, при точном закреплении равна xi²P/2. Общая с узловой
пружиной ведущая задача не утверждает равенства полных динамических моделей.
Вывод, численные проверки и ограничения относятся только к
[новому отчёту](isotropic_two_springs_external_axial_spring.md);
поперечная податливость, вязкость и замена поворотной пружины не вводятся.

## Поперечная пружина на внешнем конце второго стержня

2026-10-25, [LONG-D29](../memory/decisions.md#long-d29).
Отдельная модель с одинаковыми сечениями и прежними mu=.01, epsilon=.005,
b/h=4. Обе поворотные пружины и две точные поступательные связи узла
сохранены. Продольное условие u2(0)=0 точное; продольной внешней пружины нет.
Только w2(0)=0 заменяется энергией kw*w2(0)^2/2 и естественным условием
−Q2(0)+kw*w2(0)=0. Локальное обозначение kappa_w=kw*l³/(EI).
Пружина и направляющая ориентированы поперёк локальной недеформированной
оси второго стержня для каждой ненапряжённой конфигурации. Масса опоры,
преднапряжение, геометрическая жёсткость и фиксированное глобальное
направление не вводятся. RIGID означает точное w2(0)=0, без штрафа.

При фиксированных kappa_w,κ,mu и ограниченном xi=beta_rad/epsilon
обе внешние продольные фиксации сохраняют суммарную податливость 2l/(EA)
и ведущую добавку xi²P/2. Изгибное пространство допускает свободное W2(0)
при конечной kappa_w; в энергию включено kappa_w*W2(0)^2/2.
Естественное условие ведущей задачи W2'''(0)+kappa_w*W2(0)=0 согласовано
с Q≈−EI*W'''/l³. Для простой первой пары C=P/(4ζ0 M0).
Это локальное ведущее приближение, не теорема сходимости собственных пар.
Его полный численный контроль остановлен по сильной невязке;
[новый отчёт](isotropic_two_springs_external_transverse_spring.md)
разделяет аналитический вывод, принятую полную модель и численный отказ.

## Отдельная модель со свободным внешним поворотом, 2026-10-25

По новому заданию [замены внешней поворотной пружины](isotropic_two_springs_external_rotational_replacement.md)
полагается κ_e=0 при κ_j=.1 или 10. Внешнее продольное закрепление u₂(0)=0
и обе поступательные связи узла точны. Для конечной κ_w внешние условия
−Q₂(0)+k_w w₂(0)=0 и M₂(0)=0; поворот ψ₂(0) не закреплён.
При RIGID точно задано w₂(0)=0, при κ_w=0 отсутствует только поперечная
пружина. Геометрия, RL_FULL/Тимошенко, локальное направление направляющей
и ненапряжённость конфигураций сохранены. Суммарная продольная податливость
остаётся 2l/(EA), ведущая добавка — ξ²P/2. Это самостоятельное исследование,
не замена условий прежних моделей и не изменение исторических статусов.

## HS-L1: отдельный ведущий предел горизонтальной опоры, 2026-10-27

Для [HS-L1](isotropic_two_springs_horizontal_leading.md) фиксируются κ_j>0,
μ и χ_s=kl/(EA)=O(1) при ε→0, β_rad=εξ и ограниченной ξ. В H5 χ_s=0.5,
поэтому κ_s=χ_s/ε²; в H6 горизонтальная связь точная. Предполагаются низшие
изгибные значения ζ=O(1), поля u=aεU, w=aW, ψ=(a/l)Ψ. Ведущая сдвиговая
связь Ψ=W' не заменяет полную модель Тимошенко в сохранённых расчётах.
После продольной минимизации добавка к жёсткостной форме равна
χ_s ξ²[W₂(0)−W₂(λ₂)]²/(1+2χ_s). Второй внешний поворот и поперечное
перемещение свободны. Точный предел насыщения задаётся равенством двух
перемещений W₂(0)=W₂(λ₂), с реакцией, исключаемой из силовых условий.
Вывод формальный на ограниченных ξ; равномерная теорема спектральной
сходимости и перенос насыщения на полную геометрию 90° не утверждаются.
Это допущения новой приведённой задачи, не изменение исторических моделей.
