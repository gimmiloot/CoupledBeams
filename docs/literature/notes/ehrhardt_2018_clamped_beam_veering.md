# Ehrhardt et al. 2018: Clamped-Beam Veering In Bending And Torsion

- PDF: `docs/literature/pdf/ehrhardt2018.pdf`
- Citation key: `ehrhardt_2018_clamped_beam_veering`
- Bibliographic form: D. A. Ehrhardt, T. L. Hill, S. A. Neild, and J. E. Cooper, "Veering and nonlinear interactions of a clamped beam in bending and torsion", Journal of Sound and Vibration, 416, 1--16, 2018, doi:10.1016/j.jsv.2017.11.045.
- Role for CoupledBeams: strongest close beam-like analogue by mechanism; not a direct geometry source.

## What Matters For The Project

- The paper gives a beam assembly where symmetry-preserving tuning produces crossing, while symmetry-breaking tuning produces veering.
- It uses mode-shape correlation / self-MAC style evidence to distinguish independent crossing modes from mixed veering modes.
- The physical system is a clamped-clamped beam with a perpendicular cross-beam and movable tip masses, so it is close in beam mechanics even though the geometry differs from the two rigidly connected rods.
- It is a useful analogue for interpreting mode-shape mixing and modal-character exchange, especially when frequency curves alone are ambiguous.

## Terminology And Notation

- Important terms: `linear normal mode veering`, `nonlinear normal modes`, `crossing`, `veering`, `bending`, `torsion`, `self-MAC`, mode-shape correlation.
- The tuning variables are movable-mass positions and asymmetry measures, not the project parameter `mu`.
- The paper distinguishes symmetric configurations, where crossing can occur, from asymmetric configurations, where mixed bending/torsion veering is observed.

## Critical Places In The Paper

- pp. 1--2: introduction explains close eigenvalues, perturbations in symmetry, eigenvector sensitivity, self-MAC/correlation, and crossing versus veering.
- Sec. 2 and Fig. 1: physical clamped-clamped cross-beam assembly and finite-element/reduced-order modelling context.
- Sec. 3.1 and Fig. 3: key linear result; symmetric tuning crosses, asymmetric tuning veers, with self-MAC/correlation used as modal evidence.
- Sec. 3.2 and Figs. 4--5: nonlinear normal mode crossing/veering analogue and deformation-shape mixing.
- Secs. 4--5 and Figs. 6--8: forced-response and experimental comparisons; useful if the project later needs nonlinear or experimental language, secondary for the present linear `mu` question.
- Sec. 6: concise conclusion that the tunable beam demonstrates linear crossing/veering and related nonlinear behaviour.

## Mismatch With CoupledBeams

- The mechanism is close, but the structure is a clamped beam plus cross-beam with movable masses, not two rods coupled at an angle through the verified CoupledBeams determinant.
- Bending-torsion interaction is not the same as the current axial/transverse coupling and arm-order evolution.
- The nonlinear normal mode material should not be used as evidence for the current strict linear `mu` claim.

## Use As Source

Use as a strong auxiliary source and the best beam-like mechanism analogue among the four new papers. It is especially useful for the crossing-versus-veering distinction under symmetry breaking and for the need to check mode-shape correlation rather than relying only on frequency curves.
