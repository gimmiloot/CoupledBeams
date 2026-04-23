# Fontanela et al. 2021: Nonlinear Localisation In Two Coupled Beams

- PDF: `docs/literature/pdf/s11071-020-05760-x.pdf`
- Citation key: `fontanela_2021_nonlinear_localisation_coupled_beams`
- Bibliographic form: F. Fontanela, A. Vizzaccaro, J. Auvray, B. Niedergesäß, A. Grolet, L. Salles, and N. Hoffmann, "Nonlinear vibration localisation in a symmetric system of two coupled beams", Nonlinear Dynamics, 103(4), 3417--3428, 2021, doi:10.1007/s11071-020-05760-x.
- Role for CoupledBeams: source for localization in a two-beam system; not a primary source for strict linear veering.

## What Matters For The Project

- The paper gives a symmetric two-beam experimental system where nonlinear effects produce localized vibration on one beam.
- It is useful for the project's localization vocabulary and for the warning that localization can arise through nonlinear bifurcation, not only through linear veering.
- The experimental setup is geometrically closer than many plate or abstract examples because it uses two weakly coupled cantilever beams, but the mechanism is contact/clearance nonlinearity.
- It supports localization discussion, not a strict linear avoided-crossing claim.

## Terminology And Notation

- Important terms: `vibration localisation`, `symmetry breaking bifurcation`, `clearance nonlinearity`, `piecewise linear stiffness`, `in-phase mode`, `out-of-phase mode`, `localised state`.
- The model uses two masses with grounding stiffness/damping, coupling stiffness `k_c`, and nonlinear stopper/contact force `f_nl`.
- The spelling `localisation` follows the paper; project text may use `localization` consistently while noting the source title.

## Critical Places In The Paper

- pp. 3417--3418: abstract and introduction state the nonlinear-localisation claim and the weakly coupled two-beam experimental validation.
- Sec. 2.1 and Eqs. (1)--(6): two-degree-of-freedom model, coupling spring, nonlinear force, and the linear in-phase/out-of-phase baseline.
- Figs. 3--4: backbone curves and bifurcating localized branches.
- Sec. 3.1 and Fig. 7: experimental two-beam test setup, weak coupling through the slender connection.
- Figs. 9--11: measured linear/nonlinear responses and the localized state on either beam.
- Sec. 4: summary that nonlinear localized states can be stable and triggered in the driven system.

## Mismatch With CoupledBeams

- The decisive mechanism is clearance/contact nonlinearity and finite-amplitude response, not the linear eigenvalue problem in the current CoupledBeams `mu` sweep.
- The beams are weakly coupled cantilevers represented by a two-degree-of-freedom model, not rods rigidly coupled at an angle.
- It should not be used as a main source for strict veering or for the determinant/root structure.

## Use As Source

Use as a strong auxiliary source for localization in two-beam systems and as background for arm-wise localization metrics. Use only as background, not primary evidence, when the claim is strict linear veering.
