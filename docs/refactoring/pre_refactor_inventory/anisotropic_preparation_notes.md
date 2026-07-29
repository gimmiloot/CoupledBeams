# Anisotropic-rods preparation notes

This document is preparation only. It does not implement an anisotropic model,
approve an API, change the isotropic determinant, or create any proposed future
directory.

## 1. Current isotropic assumptions embedded in code

The current project does not have one neutral cross-section/material data
layer. Isotropic circular-rod assumptions are distributed across stable source
and diagnostic code.

| Assumption | Representative locations | Current expression/role | Refactor concern |
| --- | --- | --- | --- |
| scalar `E` and `rho` | `src/my_project/analytic/formulas.py::BeamParams` | dimensional frequency scale and circular section | scalar properties cannot express directional stiffness |
| scalar `E`, `rho`, `A`, `I` | `src/my_project/fem/python_fem.py` | `EI`, `EA`, `rhoA`; baseline nondimensional frame | production baseline must remain unchanged and provide isotropic-limit oracle |
| scalar Poisson ratio `nu` | `scripts/lib/variable_length_timoshenko.py` | global `NU`; used in `G` and circular `kappa` | anisotropy cannot generally derive shear response from one `nu` |
| `G=E/(2(1+nu))` | `variable_length_timoshenko.py`, `out_of_plane_fem_1d.py`, several comparison/FEM scripts | isotropic shear modulus | invalid as a general anisotropic constitutive rule |
| circular radius | `formulas.py::BeamParams.r`, thickness mismatch factors, Timoshenko `Section` | one scalar defines section geometry | noncircular sections need axes, orientation, multiple inertias, warping assumptions |
| `A=pi*r^2` | `formulas.py`, `python_fem.py`, `variable_length_timoshenko.py`, solid-FEM diagnostics | area and axial/mass scaling | currently duplicated in stable and diagnostic layers |
| `I=pi*r^4/4` | same locations | one bending inertia, equal in both principal directions | hides directional bending stiffness and material/geometry axes |
| polar inertia `J=2I` for circle | `out_of_plane_fem_1d.py` and out-of-plane theory/code | torsional stiffness and rotary mass | Saint-Venant torsion for noncircular/anisotropic sections needs reviewed torsion constant/warping model |
| circular shear coefficient | `variable_length_timoshenko.py::circular_shear_coefficient` and copied diagnostic implementations | `kappa=(6+12nu+6nu^2)/(7+12nu+4nu^2)` | geometry/material-specific modeling choice, not universal |
| single scalar `EI` | baseline analytic/FEM and Timoshenko diagnostics | one in-plane bending channel per arm | anisotropy may require `EI_y`, `EI_z`, coupling, or a matrix |
| single scalar `EA` | baseline analytic/FEM and axial rows | axial force law | extension–bending coupling cannot be represented |
| single scalar `GJ` | out-of-plane EB + torsion subsystem | Saint-Venant torsion | bending–torsion coupling and material-axis rotation are absent |
| diagonal energy partitions | Timoshenko energy helpers | bending, shear, axial terms summed separately | coupled constitutive energy needs cross terms and positive-definiteness checks |
| identical material on both arms | many defaults and geometry constructors | one `E`, `nu`, `rho` reused | different arm materials are not first-class data |
| circular mass-preserving eta scaling | `ThicknessMismatchFactors`, `TauFactors` | area scales `tau^2`, inertia `tau^4` | this is a physical circular-radius model, not generic section scaling |

The duplication inventory identifies the central candidates, but the first
anisotropic step should not be “deduplicate all formulas.” The physical model
and constitutive variables must be decided before a common layer is designed.

## 2. Candidate common data model

The following names are proposals, not an approved API.

### `RodEffectiveProperties`

Candidate immutable value object for section-level effective properties after
the beam theory is selected. Possible fields:

```text
mass_per_length
rotary_inertia_tensor_or_components
axial_stiffness
bending_stiffness_y
bending_stiffness_z
torsional_stiffness
shear_stiffness_y
shear_stiffness_z
extension_bending_coupling
bending_torsion_coupling
reference_axes
units_and_normalization
provenance
```

For a generally coupled theory, the scalar stiffness fields should be views or
validated reductions of a constitutive matrix, not independent duplicated
sources of truth.

### `CoupledRodGeometry`

Candidate geometry/topology value object:

```text
arm_lengths
centerline_tangents
joint_angle_beta
cross_section_descriptors_per_arm
local_frames
right_arm_local_to_global_transform
reference_length
joint_and_boundary_condition_metadata
```

It must preserve the production convention
`q_global = T @ q_local` and the existing right-arm axial/transverse mappings.
Geometry must not silently include constitutive assumptions.

### `MaterialOrientation`

Candidate representation of material axes relative to each rod/section frame:

```text
orientation_matrix_or_minimal_angles
material_frame_definition
section_frame_definition
handedness
variation_along_s
swap_operation
```

Any orientation matrix requires orthogonality/determinant checks and an
explicit local-to-global direction. A single angle is sufficient only after a
planar orthotropic scope is approved.

### `ConstitutiveMatrix`

Candidate generalized force–strain relation, with ordering frozen explicitly,
for example a reviewed subset of

```text
[N, M_y, M_z, T, Q_y, Q_z]^T = C [epsilon, kappa_y, kappa_z, twist, gamma_y, gamma_z]^T
```

The ordering above is illustrative, not approved. The chosen matrix must state
beam theory, symmetry, units, positive-definiteness/stability requirements,
allowed coupling terms, and rotation law between material and rod frames.

## 3. Scientific questions to resolve first

1. Is the target material orthotropic, transversely isotropic, or generally
   anisotropic?
2. Is the first model strictly planar, or must it support full spatial bending
   and torsion?
3. Are material axes fixed to each cross-section, and how do they rotate on
   the right arm?
4. Is extension–bending coupling physically allowed by material/section
   symmetry?
5. Is bending–torsion coupling required, and is warping restrained or free?
6. Can the two arms use different materials as well as different sections?
7. What is the exact rod-swap transformation for geometry, material axes, and
   constitutive matrices, and should it preserve the sorted spectrum?
8. What isotropic limit must reproduce the existing in-plane determinant,
   out-of-plane EB+torsion model, and FEM baseline?
9. Which beam theory is required: Euler–Bernoulli, Timoshenko, a coupled
   anisotropic beam, Vlasov/warping torsion, or another reduction?
10. What independent FEM reference is credible: anisotropic 1D beam elements,
    shell/solid material orientations, or both?

Physical meaning gates for every new assumption:

- dimensional consistency and energy positivity;
- intended boundary and rigid-joint conditions;
- symmetric-arm and material-axis-aligned limits;
- isotropic circular limit;
- `beta=0`, `mu=0`, and arm-swap limits;
- consistency with the production local-to-global transform;
- independence of the reference FEM from analytic implementation details.

## 4. Required verification before implementation

1. Write and approve a theory note with generalized variables, force/strain
   ordering, frames, and constitutive symmetry.
2. Symbolically verify rotation of the constitutive law and the isotropic
   reduction.
3. Define unit tests for zero coupling, orthotropic aligned axes, swapped arms,
   and circular isotropic limits.
4. Choose benchmark problems with known analytic or independent FEM results.
5. Decide how branch identity is seeded when material orientation adds new
   crossings or degeneracies.
6. Decide whether current `eta` remains a circular-radius diagnostic only or
   is replaced by a more general section-parameter model; do not overload it.
7. Keep current determinant and solvers frozen until the new model is a
   separate explicitly selected path.

## 5. Proposed future directories

Only after the scientific/API review, a later task may consider:

```text
docs/anisotropic_rods/
src/my_project/analytic/anisotropic_rods/
scripts/analysis/anisotropic_rods/
```

None of these directories was created in this task.

