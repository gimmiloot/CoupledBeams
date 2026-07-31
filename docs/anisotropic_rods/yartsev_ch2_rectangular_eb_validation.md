# Rectangular orthotropic EB validation

## 1. Scope

This note records a project comparator, not a source reproduction. It is a
finite elastic validation gate for two rectangular HMS/DX-209 rods at
`theta_1=theta_2=0`, with the existing ideal rigid angular joint. It checks an
Euler--Bernoulli plus generalized Saint-Venant limit against the unchanged
Chapter-2 Timoshenko model, exact straight-rod families, and an independent 1D
FEM.

Overall status: `PARTIAL_PASS`.

Original fixed-elements-per-arm accuracy gate: `PARTIAL_PASS`.

Targeted length-proportional refinement: `FAIL_CONVERGENCE_ORDER`.

## 2. Orthotropic endpoint theta=0

The existing material-rotation machinery is used without copying material
constants. At `theta=0`, the numerical endpoint gives

\[
\bar S_{16}=0,
\]

to machine precision. Hence bending and torsion separate. All calculations
use the real-elastic HMS/DX-209 material; no off-axis case is inferred from
this endpoint.

## 3. Rectangular geometry and section quantities

The base cross-section is `a=0.005 m`, `b=0.020 m`, so

\[
A=ab,\qquad I_y=\frac{a^3b}{12},\qquad
I_p=\frac{ab(a^2+b^2)}{12}.
\]

`A` has units `m^2`; `I_y` and `I_p` have units `m^4`. `I_p` is retained in
the torsional inertia `rho I_p`, but it is not used as the rectangular
Saint-Venant torsion constant. The unequal-length check keeps
`L_1+L_2=0.4 m`; the slender-limit check scales both `a` and `b`, preserving
`a/b=1/4`.

## 4. EB bending equations

The accepted book state and ordering are retained:

\[
y_i^{EB}=[w_i,\psi_i,\Phi_i,Q_i,M_i,M_{T,i}]^T.
\]

The bending part is

\[
w_i'=\psi_i,\qquad
\psi_i'=\frac{M_i}{E_{x,i}I_{y,i}},\qquad
Q_i'=-\rho_iA_i\omega^2w_i,\qquad
M_i'=-Q_i.
\]

There is no transverse-shear deformation and no bending rotary inertia in
this comparator. Dimensions are consistent: `M/(E_x I_y)` has units `1/m`,
and `rho A omega^2 w` has units `N/m`.

## 5. Generalized Saint-Venant torsion

The torsion equations are

\[
\Phi_i'=\frac{M_{T,i}}{C_{SV,i}},\qquad
M_{T,i}'=-\rho_iI_{p,i}\omega^2\Phi_i.
\]

The constitutive stiffness is the already implemented generalized rectangular
value

\[
C_{SV}=\bar C.
\]

At `theta=0`, the calculation gives `C_T=Cbar=C_SV` with zero relative
difference at double precision. Thus Timoshenko and EB use the same torsional
constitutive law; their endpoint differences arise only from transverse shear
and bending rotary inertia.

## 6. Why C_SV is not G I_p for a rectangle

`I_p` is the polar second moment used by the rotary kinetic energy. A
noncircular Saint-Venant torsion problem has a nonuniform warping/shear field,
so its torsional rigidity is not obtained by multiplying a shear modulus by
`I_p`. This comparator therefore reuses the converged rectangular generalized
torsion series for `Cbar`; it never substitutes `G I_p`.

## 7. External EB clamps

Because `w'=psi` in this EB theory, slope and section-rotation clamps coincide:

\[
w_i(0)=\psi_i(0)=\Phi_i(0)=0.
\]

With physical reaction unknowns `r_i=[Q_i(0),M_i(0),M_{T,i}(0)]^T`, the clamp
map is

\[
B_i^{EB}=\begin{bmatrix}0_{3\times3}\\I_{3\times3}\end{bmatrix}.
\]

No additional clamp variant is introduced.

## 8. Reused rigid-joint conditions

The project bases remain `e_z`, `t_i`, `n_i=e_z x t_i`, with physical vectors

\[
\boldsymbol\vartheta_i=\Phi_i\,\mathbf t_i-\psi_i\,\mathbf n_i,\qquad
\mathbf m_i=M_{T,i}\,\mathbf t_i-M_i\,\mathbf n_i.
\]

The existing `J_book(beta)` is imported from the rigid-joint helper. Its signs,
rows, unknown ordering, and physical vector convention are not copied or
changed.

## 9. EB coupled boundary matrix

For each arm,

\[
H_i^{EB}(\omega)=\exp(A_i^{EB}(\omega)L_i)B_i^{EB},
\]

and the coupled matrix is

\[
D_{EB}(\omega,\beta)=J_{book}(\beta)
\operatorname{blockdiag}(H_1^{EB},H_2^{EB}).
\]

It is `6 x 6`. The implementation retains physical raw matrices and a
positively row/column-equilibrated matrix for root search. CSV evidence keeps
both raw and scaled determinant/SVD diagnostics; no expanded determinant is
formed.

## 10. Direct straight EB reference

The independent straight reference of length `L_total` is

\[
D_{EB}^{straight}=C_{right}^{EB}\exp(A_{EB}L_{total})B_{left}^{EB},
\]

where `C_right^EB` selects `[w,psi,Phi]`. This `3 x 3` matrix does not call the
coupled determinant and is the reference for artificial split invariance.

## 11. Exact fixed-fixed bending and torsion limits

For bending,

\[
\lambda_b^4=\frac{\rho A\omega^2}{E_xI_y},\qquad
\cosh(\lambda_bL)\cos(\lambda_bL)-1=0.
\]

For torsion,

\[
\omega_n=\frac{n\pi}{L}\sqrt{\frac{C_{SV}}{\rho I_p}}.
\]

The first three transfer roots agree with the respective exact families to
`9.44e-14` and `1.58e-13` maximum relative error. Family labels are kept
separate because the full sorted spectrum interleaves them.

## 12. Unequal-length invariance

The splits `(0.20,0.20)`, `(0.15,0.25)`, `(0.10,0.30)`, and `(0.05,0.35) m`
were compared with independent direct rods of length `0.4 m` at `beta=0`.
For the first seven roots, the maximum relative differences are
`5.19e-12` for Timoshenko and `1.23e-12` for EB, both below `1e-8`. The
physical end-map recovery and arm-specific transfer scaling therefore pass
this artificial-joint-position check.

## 13. Timoshenko vs EB comparison

The first six sorted roots plus a seventh guard were compared for equal
`0.2 m` arms at `beta=0,30,90 deg`. Neighbor gaps showed no near-degenerate
cluster under the declared `1e-3` relative-gap diagnostic, so sorted matching
was used. No branch identity or accuracy ranking is claimed: `delta_f` and
`delta_lambda` measure only the difference between the two continuum models.

## 14. Slender-limit check

At `beta=0`, the first three bending families were isolated before comparing
the models. For section scales `1`, `0.5`, and `0.25`, the Timoshenko--EB
relative differences decrease monotonically for all three families. For
example, mode 1 decreases from `6.18e-3` to `1.55e-3` and `3.87e-4`. This is a
mathematical limit check, not a thickness parameter study.

## 15. Independent 1D FEM formulation

Each node has local DOFs `[w,psi,Phi]`. Bending uses the standard two-node
Hermite EB stiffness `EI/L_e^3` matrix and consistent mass
`rho A L_e/420` matrix. Torsion uses

\[
K_t=\frac{C_{SV}}{L_e}\begin{bmatrix}1&-1\\-1&1\end{bmatrix},\qquad
M_t=\frac{\rho I_pL_e}{6}\begin{bmatrix}2&1\\1&2\end{bmatrix}.
\]

The two outer nodes are fixed in all three DOFs. The generalized symmetric
eigenproblem is solved for seven positive roots. This code does not use the
old circular/isotropic FEM and does not implement a Timoshenko FEM.

## 16. FEM rigid-joint implementation

The FEM joint is assembled independently with
`q_J=[w_J,theta_t,theta_n]^T`. Arm 1 uses
`[w,psi,Phi]=[w_J,-theta_n,theta_t]`; arm 2 uses
`w=w_J`, `psi=sin(beta) theta_t+cos(beta) theta_n`, and
`Phi=-cos(beta) theta_t+sin(beta) theta_n`. Local endpoint DOFs are mapped to
these common DOFs by congruence. `J_book` is not used to impose FEM
constraints; it is used only afterward to check the kinematic residual.
Post-assembly generalized joint reactions satisfy equilibrium within the
declared numerical residual gate.

## 17. FEM convergence and comparison

Meshes `4,8,16,32,64` elements per arm were run for equal arms at
`beta=0,30,90 deg` and for `(0.1,0.3) m` at `beta=0`. All reduced stiffness and
mass matrices are symmetric, all reduced mass matrices are positive definite,
no zero mode remains after the outer clamps, and seven roots are present. The
first six errors converge without systematic divergence and mesh 64 meets the
`5e-4` modes-1--6 gate (`3.71e-4` worst case).

The strict first-three `1e-5` gate is not met: the worst case is `5.18e-5`
(the equal `beta=0` case is `2.51e-5`). This isolated shortfall is retained as
`PARTIAL_PASS`. It is consistent with the second-order dispersion of the
specified linear torsion element with consistent rotary-inertia mass; no
element coefficient, stiffness, mass, or matching rule was tuned to remove
it.

## Targeted length-proportional FEM refinement

The original unequal-length mesh used the same element count on both arms.
For `(L_1,L_2)=(0.1,0.3) m` and `(N_1,N_2)=(64,64)`, this means
`h_1=0.0015625 m` but `h_2=0.0046875 m`. That original raw result remains
`PARTIAL_PASS`, with a maximum first-three relative error of
`5.181993e-5`.

The targeted run retained the same analytic model, FEM element stiffness and
consistent mass, joint reduction, dense generalized eigensolver, and
acceptance thresholds. Only the two arm element counts were allowed to
differ. The prescribed sequence gave equal physical element lengths:

| N1 | N2 | h1 (m) | h2 (m) | reduced size | max error, roots 1--3 | max error, roots 1--6 |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 24 | 0.0125 | 0.0125 | 93 | `4.016435e-4` | `3.618169e-3` |
| 16 | 48 | 0.00625 | 0.00625 | 189 | `1.004021e-4` | `9.038326e-4` |
| 32 | 96 | 0.003125 | 0.003125 | 381 | `2.510018e-5` | `2.259131e-4` |
| 64 | 192 | 0.0015625 | 0.0015625 | 765 | `6.180685e-6` | `5.648621e-5` |

For the first torsional root, mode 2, the relative-error sequence is
`4.016435e-4`, `1.004021e-4`, `2.510018e-5`, and `6.180685e-6`.
The three empirical orders are `2.000125`, `2.000020`, and `2.021859`.
The diagnostic Richardson frequency is `766.861317901724 Hz`, coinciding
with the accepted analytic value at the saved precision. It is diagnostic
only; all gates use the raw FEM roots.

At `(64,192)`, the analytic frequencies for roots 1--7 are
`353.18108306`, `766.86131790`, `973.55725545`, `1533.72263580`,
`1908.56148220`, `2300.58395371`, and `3067.44527161 Hz`. The raw FEM
relative errors are `6.416670e-7`, `6.180685e-6`, `1.237610e-8`,
`2.513662e-5`, `1.541698e-8`, `5.648621e-5`, and `1.004084e-4`.
Thus the unchanged first-three and first-six accuracy gates both pass, and
the seventh positive root is present.

The equal-physical-step control `(L_1,L_2)=(0.2,0.2) m`,
`(N_1,N_2)=(128,128)`, uses the same `h=0.0015625 m`. Its first-three
maximum error is also `6.180685e-6`, and its seven FEM frequencies agree
root by root with the finest unequal-split case to the saved precision. The
unequal-length joint assembly therefore introduces no separate anomaly.

All symmetry, positive-mass, zero-mode, root-count, joint-kinematic, and
joint-equilibrium checks pass with their prior tolerances. However, the
unchanged refinement gate requires every one of roots 1--6 to be no worse on
the finest mesh than on the preceding level, with the prior small numerical
allowance. Mode 1 increases from `2.659182e-8` to `6.416670e-7` between
`(32,96)` and `(64,192)`. At the finest level the raw matrix condition
numbers are `cond(K)=2.467e12` and `cond(M)=8.601e7`, consistent with a dense
generalized-eigenvalue conditioning floor, but the failure is not waived.

The targeted status is therefore `FAIL_CONVERGENCE_ORDER`. The original
accuracy shortfall is closed by the raw finest result without changing the
model or threshold, but the overall rectangular EB validation remains
`PARTIAL_PASS` because the targeted refinement did not pass every declared
gate.

## 18. Limitations

This finite gate validates only the declared orthotropic endpoint, lengths,
angles, first-seven sorted spectra, and 1D EB continuum FEM. It does not make
the coupled formulation a stable baseline, prove root completeness beyond the
guard root, establish branch identity, or validate a physical three-dimensional
joint. Generated CSV/report files remain ignored local evidence.

## 19. Explicit exclusions

- no `theta != 0`;
- no material bending--torsion coupling in the EB comparator;
- no complex roots or damping;
- no unequal thickness;
- no Timoshenko FEM or 3D FEM;
- no coupled mode shapes or MAC;
- no broad beta, length-ratio, or parameter maps;
- no production anisotropic API.
