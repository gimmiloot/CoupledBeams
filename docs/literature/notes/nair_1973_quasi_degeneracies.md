# Nair and Durvasula 1973: Quasi-Degeneracies In Plate Vibration Problems

- PDF: `docs/literature/pdf/nair1973.pdf`
- Citation key: `nair_1973_quasi_degeneracies`
- Bibliographic form: P. S. Nair and S. Durvasula, "On quasi-degeneracies in plate vibration problems", International Journal of Mechanical Sciences, 15, 975--986, 1973.
- Role for CoupledBeams: auxiliary source for spectral-reorganization vocabulary and criteria; not a direct geometry analogue.

## What Matters For The Project

- The paper separates true frequency crossings from cases where two frequency curves approach, come close, and veer away.
- It recommends the term `quasi-degeneracy` for what older plate literature often called a `transition`.
- It ties rapid changes of nodal patterns and modal character to the neighbourhood of quasi-degeneracies.
- It gives a useful warning for the current project: rapid modal reorganization is not enough, by itself, to claim strict veering.

## Terminology And Notation

- Important terms: `frequency crossing`, `transition`, `quasi-degeneracy`, `symmetry group`, `nodal pattern`.
- Plate-specific notation includes plate dimensions `a`, `b`, stiffness `D`, operators `H`, `H'`, perturbation operator, eigenvalue parameters, and mode shapes/eigenfunctions.
- The notation is not imported into CoupledBeams; only the distinction between crossing, quasi-degeneracy, and rapid modal-pattern change is reused.

## Critical Places In The Paper

- pp. 975--976: title-page metadata, summary, terminology, and the statement that `transitions` are preferably called `quasi-degeneracies`.
- pp. 975--976: notation block for the perturbation/symmetry discussion.
- Analytical discussion around the symmetry cases: curves from different symmetry groups can cross; curves from the same symmetry group cannot cross.
- Conclusion around pp. 985--986: in unsymmetric plate configurations, crossings are absent and only quasi-degeneracies occur.
- Representative rectangular/skew-plate examples and figures: useful for the link between close frequency curves and rapid nodal-pattern changes.

## Mismatch With CoupledBeams

- The paper is about plate vibration problems, not two rigidly connected Euler--Bernoulli rods.
- Its symmetry-group crossing criteria are plate-specific and cannot be copied as a proof for the coupled-rod system.
- The project may use the paper for mechanism-level vocabulary, but any CoupledBeams veering claim still requires tracked branch pairs, a small paired gap, and pairwise mode-shape evidence from local data.

## Use As Source

Use as a supporting source for `quasi-degeneracy`, true-crossing caution, and rapid modal-character change. Do not use it as the main source for the CoupledBeams geometry or as standalone evidence for strict veering in the `beta = 15 deg`, `r = 5 mm` data.
