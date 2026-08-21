# Banerjee and Williams 1996: Composite Timoshenko Dynamic Stiffness

- PDF: `docs/literature/pdf/banerjee_williams_1996_composite_timoshenko_dynamic_stiffness.pdf`
- Citation key: `banerjee_williams_1996_composite_timoshenko_dynamic_stiffness`
- Bibliographic form: J. R. Banerjee and F. W. Williams, "Exact dynamic stiffness matrix for composite Timoshenko beams with applications", Journal of Sound and Vibration 194(4), 573--585 (1996), doi:10.1006/jsvi.1996.0378.
- Role for CoupledBeams: methodological source for exact dynamic stiffness and composite Timoshenko eigenfrequency calculations.

## What Matters For The Project

- The dynamic stiffness matrix includes transverse shear deformation and section rotary inertia.
- The governing theory includes material bending--torsion coupling caused by fibrous-composite anisotropy.
- Explicit matrix elements are used with the Wittrick--Williams algorithm to calculate free-vibration frequencies.
- The paper demonstrates how one uniform member can be used in an assembly, but it does not supply the CoupledBeams rigid-joint law.

## Terminology And Notation

- Important terms: exact dynamic stiffness matrix, composite Timoshenko beam, material bending--torsion coupling, shear deformation, rotatory inertia and Wittrick--Williams algorithm.
- The first-page model uses bending displacement, section rotation and torsional rotation for one beam; warping stiffness is neglected in the displayed governing equations.

## Critical Places In The Paper

- pp. 573--574 and Sec. 1: scope, prior composite-beam work and motivation for including shear and rotary inertia.
- Sec. 2: governing equations and derivation of the exact dynamic stiffness terms.
- Sec. 3: assembly/use of the dynamic stiffness matrix and the Wittrick--Williams count.
- Sec. 4: numerical comparisons with and without the refined effects.
- Sec. 5: conclusions and stated computational benefits.

## Mismatch With The Current Model

- The paper's member matrix is not the Yartsev Chapter-2 state transfer and must not replace it silently.
- The paper does not derive the physical compatibility and equilibrium conditions for the current rigid angular connection.
- Its treatment of torsion and warping is not identical to the project's generalized-torsion source line.

## Use As Source

Use as a method reference for exact composite Timoshenko dynamic stiffness and eigenfrequency calculation. Do not use it as the derivation of the submitted article's determinant or joint conditions.
