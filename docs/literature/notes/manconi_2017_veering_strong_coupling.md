# Manconi and Mace 2017: Veering And Strong Coupling Effects In Structural Dynamics

- PDF: `docs/literature/pdf/vib_139_02_021009.pdf`
- Citation key: `manconi_2017_veering_strong_coupling`
- Bibliographic form: E. Manconi and B. Mace, "Veering and Strong Coupling Effects in Structural Dynamics", Journal of Vibration and Acoustics, 139(2), article 021009, 2017, doi:10.1115/1.4035109.
- Role for CoupledBeams: key general theoretical source for the veering line; primary by mechanism, not by geometry.

## What Matters For The Project

- The paper gives the cleanest general distinction needed here: rapid veering is associated with weak coupling, while strong coupling can give slow evolution of eigenvalue loci rather than rapid veering.
- It defines an `uncoupled-blocked system`, its `skeleton`, and `critical points` where skeleton branches intersect.
- It links rapid eigenvalue changes near critical points with rapid eigenvector rotation.
- It gives a conservative language for the current `mu` question: if our branches evolve slowly or do not show a small paired gap, Manconi and Mace support avoiding a strict veering claim.

## Terminology And Notation

- Important terms: `mode veering`, `weak coupling`, `strong coupling`, `uncoupled-blocked system`, `skeleton`, `critical point`, eigenvector rotation.
- Their small parameter is a generic coupling order, not the project parameter `mu`.
- The `skeleton` construction is useful conceptually, but should not be imported as project notation unless a corresponding blocked reference system is explicitly defined for CoupledBeams.

## Critical Places In The Paper

- p. 021009-1: abstract and introduction state the weak-coupling / strong-coupling distinction and define the general purpose of the paper.
- Sec. 2, especially Eq. (3): weak coupling in terms of small off-diagonal mass/stiffness/gyroscopic coupling relative to diagonal terms.
- Sec. 2.2 and Eqs. (17)--(19): perturbation analysis near a critical point, with rapid eigenvalue and eigenvector changes in a small parameter neighbourhood.
- Figs. 1, 3, and 5: skeleton/coupled-locus picture and eigenvector rotation for weak coupling.
- Fig. 6 and Sec. 5.1.2: strong coupling example with gradual eigenvalue evolution instead of rapid veering.
- Sec. 5.3 and Fig. 15: continuous-system examples, including a plate-beam case, useful by mechanism but not as a direct geometry match.
- Conclusion around pp. 021009-9--021009-10: summary of when veering occurs and why weak coupling matters.

## Mismatch With CoupledBeams

- The examples are general discrete/continuous structural systems, not the verified two-rod coupled-beam determinant of this project.
- The paper's blocked coordinates and coupling order are model-dependent; they do not automatically define a `mu = 0` or `mu`-sweep skeleton for CoupledBeams.
- It supports interpretation, not a proof that any selected CoupledBeams branch has strict veering.

## Use As Source

Use as a main theoretical source for the veering/slow-evolution distinction, especially when discussing whether the `mu` sweep shows strict veering or only gradual modal-character reorganization. It should be cited before weaker or geometry-only sources when the point is `rapid veering requires weak coupling`.
