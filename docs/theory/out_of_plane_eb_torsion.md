# Out-of-Plane EB + Saint-Venant Torsion Subsystem

This note fixes the sign and coordinate conventions for a new linear
out-of-plane subsystem of the planar two-beam model. It is a theory and
lightweight determinant note only: it is not an article task, plotting
workflow, FEM workflow, or replacement for the existing in-plane determinant.

## 1. Scope

The model describes linear out-of-plane oscillations of a planar two-beam
system. Each rod has Euler--Bernoulli out-of-plane bending and Saint-Venant
torsion.

The displacement/rotation variables are:

```text
z_i(s_i)      out-of-plane displacement
phi_i(s_i)    bending rotation component
psi_i(s_i)    torsional rotation
```

The force/moment variables are:

```text
Q_i            out-of-plane shear force
M_i            bending moment
T_i            torsional moment
```

The in-plane problem uses bending plus axial motion, with variables such as
`w_i` and `u_i`. The out-of-plane problem uses out-of-plane bending plus
torsion, with variables `z_i` and `psi_i`. In the ideal linear planar model,
these two subsystems are independent.

## 2. Local Coordinate Convention

The coordinate convention is part of the model and must not be changed without
a separate theory/code consistency audit.

- The undeformed two-beam system lies on a horizontal table.
- `e_z` is directed downward, normal to the table.
- `t_1` is directed from the clamp of rod 1 to the coupling joint.
- `t_2` is directed from the clamp of rod 2 to the coupling joint.
- `n_i` is defined by

```text
n_i = e_z x t_i
```

Therefore:

```text
t_i x n_i = e_z
```

The coupling angle `beta` is measured so that:

```text
t_2 = -t_1 cos(beta) + n_1 sin(beta)
n_2 = -t_1 sin(beta) - n_1 cos(beta)
```

Do not denote the tangent vectors by `tau_i`: `tau_i` is already used for
thickness factors. Use `t_i` for tangents and `tau_i` only for thickness
factors.

Checks:

```text
beta = 0:
    t_2 = -t_1
    n_2 = -n_1

beta = 90 deg:
    t_2 = n_1
    n_2 = -t_1
```

## 3. Rotation and Moment Vector Convention

The joint rotation vector on rod `i` is decomposed as:

```text
theta_i = psi_i t_i + phi_i n_i
```

where `psi_i` is torsional rotation about `t_i`, and `phi_i` is the bending
rotation component about `n_i`.

The moment vector on rod `i` is decomposed as:

```text
m_i = T_i t_i + M_i n_i
```

where `T_i` is torsional moment, and `M_i` is bending moment.

This choice is coupled to the sign convention for `phi_i` and `M_i`. The
virtual-work pairing is:

```text
m_i . delta theta_i = T_i delta psi_i + M_i delta phi_i
```

This is the basic consistency check for the component signs.

## 4. Joint Conditions in Force Variables

At the joint, the force-variable conditions are:

```text
z_1 = z_2

psi_1 + psi_2 cos(beta) + phi_2 sin(beta) = 0
phi_1 - psi_2 sin(beta) + phi_2 cos(beta) = 0

Q_1 + Q_2 = 0

T_1 - T_2 cos(beta) - M_2 sin(beta) = 0
M_1 - M_2 cos(beta) + T_2 sin(beta) = 0
```

Equivalently:

```text
psi_1 = -psi_2 cos(beta) - phi_2 sin(beta)
phi_1 =  psi_2 sin(beta) - phi_2 cos(beta)

T_1 = T_2 cos(beta) + M_2 sin(beta)
M_1 = M_2 cos(beta) - T_2 sin(beta)
```

The variables `phi` and `psi` are not balanced like independent scalar
forces: they are components of one and the same joint rotation vector.
Similarly, `M` and `T` are components of one and the same moment vector.

## 5. Relation to Derivatives of z and psi

Use dimensional arc coordinate `s_i` from clamp to joint.

With the chosen convention:

```text
phi_i = - d z_i / d s_i
```

because:

```text
theta_i x t_i = z_{i,s} e_z
n_i x t_i = -e_z
```

Therefore `phi_i = -z_{i,s}`. The associated bending and torsion resultants
are:

```text
M_i = E J_i d phi_i / d s_i = -E J_i d^2 z_i / d s_i^2
Q_i = d M_i / d s_i = -E J_i d^3 z_i / d s_i^3
T_i = G J_{p,i} d psi_i / d s_i
```

These signs follow from the coordinate convention, not from an arbitrary
choice. Classical Euler--Bernoulli bending theory uses displacement, rotation,
bending moment, and shear force as boundary quantities; this convention is one
consistent orientation of those quantities.

## 6. Dimensionless Variables and Solution Ansatz

Use:

```text
x_i = s_i / l_i,      0 <= x_i <= 1
```

with `x_i=0` at the clamp and `x_i=1` at the joint.

For out-of-plane bending, use the same clamped-end basis as the existing
Euler--Bernoulli determinant:

```text
z_i'''' - alpha_i^4 z_i = 0

z_i(x_i) =
    A_i (cos(alpha_i x_i) - cosh(alpha_i x_i))
  + B_i (sin(alpha_i x_i) - sinh(alpha_i x_i))
```

This basis satisfies:

```text
z_i(0) = 0
z_i'(0) = 0
```

For torsion:

```text
psi_i'' + gamma_i^2 psi_i = 0
psi_i(0) = 0
psi_i(x_i) = P_i sin(gamma_i x_i)
```

The unknown vector is:

```text
X = (A_1, B_1, A_2, B_2, P_1, P_2)^T
```

## 7. Arguments alpha_i and gamma_i

For bending:

```text
L_1 = 1 - mu
L_2 = 1 + mu

alpha_i = Lambda L_i / sqrt(tau_i)
```

where `tau_i` are the same thickness factors as in the existing thickness
mismatch model.

For torsion:

```text
gamma_i = sqrt(E/G) epsilon Lambda^2 L_i
```

For isotropic material:

```text
G = E / (2(1 + nu))
gamma_i = sqrt(2(1 + nu)) epsilon Lambda^2 L_i
```

For circular rods, `gamma_i` does not contain `tau_i`, because torsional
stiffness and polar inertia both scale as `r_i^4`.

## 8. Thickness Factors

Reuse the existing mass-preserving thickness-mismatch factors:

```text
tau_1 = (1 - eta) / sqrt(1 + 2 mu eta + eta^2)
tau_2 = (1 + eta) / sqrt(1 + 2 mu eta + eta^2)
```

For bending:

```text
a_i = tau_i^(-1/2)      rotation scale
b_i = tau_i^3           bending moment scale
c_i = tau_i^(5/2)       shear force scale
```

For torsion:

```text
e_i = tau_i^4           torsional moment scale
```

The exponent in `e_i = tau_i^4` comes from the circular-section polar moment:
`J_p ~ r^4`.

Define the torsion-to-bending moment scaling coefficient:

```text
chi_T = (J_p0 / J_0) epsilon sqrt(G/E)
```

For a circular section:

```text
J_p0 = 2 J_0
chi_T = 2 epsilon sqrt(G/E)
      = 2 epsilon / sqrt(2(1 + nu))
```

## 9. Compact Determinant Form

Define joint values:

```text
Z_i   = z_i(1)
R_i   = z_i'(1) / alpha_i
K_i   = z_i''(1) / alpha_i^2
V_i   = z_i'''(1) / alpha_i^3

Psi_i = psi_i(1)
D_i   = psi_i'(1) / gamma_i
```

For the ansatz above:

```text
Z_i = A_i C_i^- + B_i S_i^-
R_i = -A_i S_i^+ + B_i C_i^-
K_i = -A_i C_i^+ - B_i S_i^+
V_i =  A_i S_i^- - B_i C_i^+
```

where:

```text
C_i^- = cos(alpha_i) - cosh(alpha_i)
C_i^+ = cos(alpha_i) + cosh(alpha_i)
S_i^- = sin(alpha_i) - sinh(alpha_i)
S_i^+ = sin(alpha_i) + sinh(alpha_i)
```

and:

```text
Psi_i = P_i sin(gamma_i)
D_i   = P_i cos(gamma_i)
```

The six determinant equations are:

```text
Z_1 - Z_2 = 0

Psi_1 + Psi_2 cos(beta) - Lambda a_2 R_2 sin(beta) = 0

-Lambda a_1 R_1 - Psi_2 sin(beta) - Lambda a_2 R_2 cos(beta) = 0

c_1 V_1 + c_2 V_2 = 0

chi_T e_1 D_1 - chi_T e_2 D_2 cos(beta) + b_2 K_2 sin(beta) = 0

-b_1 K_1 + b_2 K_2 cos(beta) + chi_T e_2 D_2 sin(beta) = 0
```

In compact matrix form:

```text
M_perp(Lambda, beta, mu, eta, epsilon, nu) X = 0
det M_perp = 0
```

The full expanded `6x6` matrix is intentionally not repeated here. The compact
block definitions above are the preferred source for avoiding sign mistakes.

## 10. beta=0 Checks

At `beta=0`, the system splits into a bending block and a torsion block. The
bending equations, up to nonzero row scalings used in determinant assembly,
are:

```text
Z_1 - Z_2 = 0
a_1 R_1 + a_2 R_2 = 0
c_1 V_1 + c_2 V_2 = 0
b_1 K_1 - b_2 K_2 = 0
```

The torsion block is:

```text
Psi_1 + Psi_2 = 0
e_1 D_1 - e_2 D_2 = 0
```

After row permutation, the matrix is block diagonal with:

```text
bending variables: A_1, B_1, A_2, B_2
torsion variables: P_1, P_2
```

For `eta=0` and `beta=0`, bending roots must match the clamped-clamped
straight beam of total length `2`:

```text
cosh(2 Lambda) cos(2 Lambda) - 1 = 0
```

Torsion roots must satisfy:

```text
sin(gamma_1 + gamma_2) = 0
```

Since `L_1 + L_2 = 2`:

```text
2 epsilon sqrt(E/G) Lambda^2 = n pi

Lambda_n^torsion =
    sqrt(n pi / (2 epsilon sqrt(E/G)))
```

For isotropic material:

```text
sqrt(E/G) = sqrt(2(1 + nu))
```

These are sanity checks, not numerical results for the full coupled
`beta>0` problem.
