# Ryabov and Yartsev 2021: Monoclinic Composite Strip

- PDF: `docs/literature/pdf/ryabov_yartsev_2021_monoclinic_strip.pdf` (Russian original)
- Citation key: `ryabov_yartsev_2021_monoclinic_strip`
- Bibliographic form: V. M. Ryabov and B. A. Yartsev, "Nonclassical vibrations of a monoclinic composite strip", Vestnik St. Petersburg University, Mathematics 54(4), 437--446 (2021), doi:10.1134/S1063454121040166. Russian original doi:10.21638/spbu01.2021.415.
- Role for CoupledBeams: closest journal source for the rectangular monoclinic single-rod material model and terminology used by the submitted article.

## What Matters For The Project

- The paper models damped flexural--torsional vibration of a constant rectangular monoclinic composite strip.
- Its basis is refined Timoshenko bending, Voigt--Lekhnitskii generalized torsion and the elastic--viscoelastic correspondence principle with complex moduli.
- It defines the directional moduli `E_x(theta)`, `G_{xy}(theta)` and `G_{xz}(theta)` and mutual-influence coefficients that couple bending and torsion.
- At `theta=0°` and `90°` the mutual-influence coefficients vanish and the equations reduce to independent orthotropic Timoshenko bending and torsion equations.
- The model and numerical method are checked against measured natural frequencies and loss factors for free--free specimens.

## Terminology And Notation

- Main variables: transverse displacement `w`, twist angle `Phi`, `E_x`, `G_{xy}`, `G_{xz}`, mutual-influence coefficients, shear coefficient `k`, area `A`, moments `I_y` and `I_p`, density and generalized torsional rigidity.
- The local Russian article writes `k=5/6` for the rectangular section used in its model.
- The official English translation and Russian original are one scientific source and share the single citation key above.

## Critical Places In The Paper

- Russian pp. 695--696: scope, physical motivation and the rectangular strip geometry.
- Russian pp. 696--697, equations (1)--(5): coupled equations, definitions of `E_x`, `G_{xy}`, `G_{xz}`, boundary conditions and orthotropic limits.
- Russian pp. 697--698: Laplace-transform frequency solution, determinant condition and complex-frequency iteration.
- Russian pp. 698--700 and Fig. 2: comparison of calculated and experimental frequencies and loss factors.
- Later numerical sections: influence of reinforcing-fibre angle and strip length on coupled modes.

## Mismatch With The Current Model

- The paper treats one free--free or cantilever strip and does not derive a rigid joint between two arms.
- Its damped/complex-frequency workflow is broader than the real-elastic frequency calculations used for the submitted article's comparisons.
- The project must continue to use the verified local Chapter-2 formulas and recorded corrections; the journal paper is corroborating literature, not permission to rewrite them.

## Use As Source

Use as the closest journal source for the monoclinic rectangular-rod equations, modulus definitions, bending--torsion terminology, orthotropic limits and experimental validation. Pair it with the Yartsev 2024 monograph for the current Chapter-2 source line and with the separately derived project joint conditions for the two-arm geometry.
