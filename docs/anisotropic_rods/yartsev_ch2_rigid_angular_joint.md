# Yartsev Chapter-2 rigid angular joint

## 1. Scope and assumptions

This note derives the first diagnostic model of two Chapter-2 monoclinic
Timoshenko rods connected at an ideal rigid angular point joint. Each arm keeps
the accepted `state_corrected` one-rod equations without alteration.

The joint has no mass, rotary inertia, eccentricity, elastic compliance,
additional torsional spring, or independent warping coordinate. It has no
bimoment. Section warping remains condensed in the existing generalized
torsional stiffness `C_T`. Both outer ends use the source-confirmed book slope
clamp. These assumptions define an ideal point connection, not a finite joint
geometry or a validated production model.

## 2. Book notation

For arm `i`, the canonical state is

```text
y_i = [w_i, psi_i, Phi_i, Q_i, M_i, M_{T,i}]^T.
```

Here `w_i` is out-of-plane displacement, `psi_i` is the independent
Timoshenko bending rotation, `Phi_i` is twist, `Q_i` is shear force, `M_i` is
bending moment, and `M_{T,i}` is torque. Python uses `psi_b` and `Phi_t` for
the two rotation variables.

## 3. Existing project bases

The rods lie in one plane. The normal `e_z` is directed downward, each `t_i`
points from the external clamp of arm `i` to the joint, and

```text
n_i = e_z x t_i,
t_i x n_i = e_z.
```

The symbol `tau_i` is not used for tangents because it is reserved for
thickness factors.

## 4. Notation translation

The exact translation to the old isotropic out-of-plane state is recorded
separately in [the notation-translation note](yartsev_ch2_notation_translation.md).
The new derivation does not change the old model's notation.

## 5. Geometry of the two arms

With `c = cos(beta)` and `s = sin(beta)`, the unchanged project convention is

```text
t_2 = -c t_1 + s n_1,
n_2 = -s t_1 - c n_1.
```

Therefore

```text
beta = 0:       t_2 = -t_1,  n_2 = -n_1,
beta = 90 deg:  t_2 =  n_1,  n_2 = -t_1.
```

## 6. Physical rotation and moment vectors

In book notation the physical vectors are

```text
vartheta_i = Phi_i t_i - psi_i n_i,
m_i        = M_{T,i} t_i - M_i n_i.
```

The signs follow directly from the notation translation. Orthonormality gives
the energy pairing

```text
m_i . delta(vartheta_i)
  = (M_{T,i} t_i - M_i n_i)
    . (delta(Phi_i) t_i - delta(psi_i) n_i)
  = M_{T,i} delta(Phi_i) + M_i delta(psi_i).
```

This is an exact symbolic reduction using `t_i.t_i=n_i.n_i=1` and
`t_i.n_i=0`; deterministic numerical checks independently exercise it.

## 7. Invariant joint conditions

The ideal rigid point joint satisfies

```text
w_1 = w_2,
vartheta_1 = vartheta_2,
Q_1 + Q_2 = 0,
m_1 + m_2 = 0.
```

The rotation condition is equality of physical section-rotation vectors. A
condition such as `w_1'=w_2'` is not used as a replacement.

## 8. Scalar joint conditions

Expanding arm 2 in the arm-1 basis gives

```text
vartheta_2
  = Phi_2(-c t_1 + s n_1) - psi_2(-s t_1 - c n_1)
  = (-c Phi_2 + s psi_2) t_1
    + (s Phi_2 + c psi_2) n_1.
```

Equating this with `vartheta_1 = Phi_1 t_1 - psi_1 n_1` gives, without a
sign guess,

```text
Phi_1 + c Phi_2 - s psi_2 = 0,
psi_1 + c psi_2 + s Phi_2 = 0.
```

Similarly,

```text
m_2
  = M_{T,2}(-c t_1 + s n_1) - M_2(-s t_1 - c n_1)
  = (-c M_{T,2} + s M_2) t_1
    + (s M_{T,2} + c M_2) n_1.
```

The `t_1` projection of `m_1+m_2=0` is the fifth row below. Multiplying the
`n_1` projection by `-1` gives the sixth row. The complete six scalar
conditions are

```text
w_1 - w_2 = 0,
Phi_1 + Phi_2 cos(beta) - psi_2 sin(beta) = 0,
psi_1 + psi_2 cos(beta) + Phi_2 sin(beta) = 0,
Q_1 + Q_2 = 0,
M_{T,1} - M_{T,2} cos(beta) + M_2 sin(beta) = 0,
M_1 - M_2 cos(beta) - M_{T,2} sin(beta) = 0.
```

## 9. Full `J_book(beta)` matrix

For

```text
Y_J = [
  w_1, psi_1, Phi_1, Q_1, M_1, M_{T,1},
  w_2, psi_2, Phi_2, Q_2, M_2, M_{T,2}
]^T,
```

the joint equations are `J_book(beta) Y_J=0`, where

```text
J_book = [
 [1, 0, 0, 0, 0, 0,  -1,  0,  0, 0,  0,  0],
 [0, 0, 1, 0, 0, 0,   0, -s,  c, 0,  0,  0],
 [0, 1, 0, 0, 0, 0,   0,  c,  s, 0,  0,  0],
 [0, 0, 0, 1, 0, 0,   0,  0,  0, 1,  0,  0],
 [0, 0, 0, 0, 0, 1,   0,  0,  0, 0,  s, -c],
 [0, 0, 0, 0, 1, 0,   0,  0,  0, 0, -c, -s]
].
```

It has size `6 x 12`. With
`P=diag(1,-1,1,1,-1,1)` and `P_2=block_diag(P,P)`, the implemented sign gate
is

```text
J_book(beta) = J_old(beta) P_2.
```

`J_old` is formed from the conditions in
[`out_of_plane_eb_torsion.md`](../theory/out_of_plane_eb_torsion.md). Its old
bending-rotation and bending-moment rows are multiplied by `-1`, which leaves
the homogeneous conditions unchanged and fixes their row orientation so the
displayed equality is exact. The old note and implementation are not changed.

## 10. Virtual-work check

For compatible virtual motion

```text
delta w_1 = delta w_2 = delta w_J,
delta vartheta_1 = delta vartheta_2 = delta vartheta_J,
```

and equilibrated resultants `Q_2=-Q_1`, `m_2=-m_1`,

```text
delta W_J
  = (Q_1+Q_2) delta w_J
    + (m_1+m_2) . delta vartheta_J
  = 0.
```

This is the exact symbolic check. The test and pilot also use at least 100
fixed-seed random cases for `beta=0,30,90 deg`.

## 11. External source clamps

Both external ends use the confirmed book condition

```text
w_i(0)=0,
w_i'(0)=0,
Phi_i(0)=0,
w_i' = psi_i + Q_i/K_{s,i},
K_{s,i} = k G_{xz,i} A_i.
```

For physical reactions `r_i=[Q_i(0),M_i(0),M_{T,i}(0)]^T`,

```text
B_i = [
 [0,       0, 0],
 [-1/K_s,  0, 0],
 [0,       0, 0],
 [1,       0, 0],
 [0,       1, 0],
 [0,       0, 1]
].
```

The implementation reuses the existing `book_slope_clamp` matrix.

## 12. Coupled boundary matrix construction

Let `T_i^{phys}(omega)` be the physical transfer matrix recovered from the
existing scaled `state_corrected` propagation without duplicating its private
state scales. Then

```text
H_i^{phys}(omega) = T_i^{phys}(omega) B_i,
y_i(L_i) = H_i^{phys}(omega) r_i.
```

For `r=[r_1,r_2]^T`,

```text
D_coupled(omega,beta)
  = J_book(beta) block_diag(H_1^{phys},H_2^{phys}),
det D_coupled = 0.
```

The physical matrix is `6 x 6`. Numerical root searches apply only positive
row and column equilibration factors; these cannot move determinant zeros.
Raw and scaled quality measures remain separately reportable.

The independent straight reference uses a single corrected transfer matrix
over `L_total=L_1+L_2`:

```text
D_straight(omega) = C_right T^{phys}(omega,L_total) B_left,
```

where `C_right` selects `w`, `psi+Q/K_s`, and `Phi`. It does not use the
coupled matrix.

## 13. `beta=0`, `beta=90 deg`, and orthotropic limits

At `beta=0`, the scalar conditions reduce to

```text
w_1=w_2,  Phi_1+Phi_2=0,  psi_1+psi_2=0,
Q_1+Q_2=0,  M_{T,1}-M_{T,2}=0,  M_1-M_2=0.
```

At `beta=90 deg`, they reduce to

```text
Phi_1-psi_2=0,  psi_1+Phi_2=0,
M_{T,1}+M_2=0,  M_1-M_{T,2}=0.
```

At the orthotropic endpoint `theta_1=theta_2=0`, `Sbar16_i=0`, so bending and
torsion decouple within each arm. At `beta=0`, row/column permutation must then
separate a bending `4 x 4` block and a torsion `2 x 2` block.

## 14. Warping boundary

No independent warping coordinate and no bimoment are introduced. Generalized
torsion remains exactly the existing Chapter-2 model condensed into `C_T`.

## 15. Pilot verification status

Status: `PASS` for the declared rigid-joint and small elastic pilot gate. The
sign transformation, ranks, deterministic random states, virtual work,
`beta=0`/`90 deg` limits, orthotropic block split, equal-arm exchange, first
seven `beta=0` roots against the independent straight rod, all 21 pilot-root
quality checks, and all three seventh-root guards passed. Generated evidence
is under `results/anisotropic_rods/yartsev_ch2_coupled_joint_pilot/`. This is
not a stable baseline, final coupled-rod model, or parameter study.

## 16. Explicit exclusions

- no Euler--Bernoulli model;
- no Saint-Venant replacement or comparison;
- no complex roots or damping;
- no FEM;
- no parameter maps, shapes, MAC, or localization study.
