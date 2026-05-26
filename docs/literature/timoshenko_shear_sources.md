# Timoshenko Beam Theory and Shear-Coefficient Sources

## Purpose

This note collects local PDF sources for the planned thickness-mismatch and
diagnostic Timoshenko-correction work. It focuses on Timoshenko beam theory,
the shear correction coefficient, circular/cylindrical rods, critical-frequency
limitations, and frame-oriented extensions that may later matter for coupled
beams. The bibliographic data below were extracted from the local PDFs in
`docs/literature/pdf/`; no internet lookup was used.

## Source table

| ID | PDF file | Bibliographic entry | Role in project | Notes |
|---|---|---|---|---|
| A2-1 | `A2-1.pdf` | Thomas Kramer and Michael Helmut Gfrerer, "A numerical assembly technique for free vibration of plane frames using a shifted and deflated Newton method," *Journal of Sound and Vibration*, 583, 118435, 2024. DOI: `10.1016/j.jsv.2024.118435`. | Frame-oriented Timoshenko-Ehrenfest beam source. Useful for future coupled-beam/frame assembly with axial and transverse dynamics. | Uses a rectangular-section shear coefficient `k = 5/6`; this is not a source for the circular-rod coefficient. |
| A2-2 | `A2-2.pdf` | W. P. Howson and F. W. Williams, "Natural frequencies of frames with axially loaded Timoshenko members," *Journal of Sound and Vibration*, 26(4), 503--515, 1973. DOI: not found in local PDF. | Classic dynamic-stiffness frame source using Timoshenko members, axial load, rotary inertia, and shear deflection. Useful as a frame benchmark and as evidence that Timoshenko members were used in plane-frame frequency calculations. | Local PDF is scanned; title-page data and pages were read from the first-page image. |
| A2-3 | `A2-3.pdf` | G. M. L. Gladwell, "The vibration of frames," *Journal of Sound and Vibration*, 1(4), 402--425, 1964. DOI: not found in local PDF. | General frame-vibration background source. Useful for historical frame methods and Rayleigh--Ritz style frame modeling, but not a Timoshenko or shear-coefficient source. | Local PDF is scanned; metadata contains PII `0022-460X(64)90056-2`. |
| A2-4 | `A2-4.pdf` | Alejandro R. Ratazzi, Diana V. Bambill, and Carlos A. Rossit, "Free Vibrations of Beam System Structures with Elastic Boundary Conditions and an Internal Elastic Hinge," *Chinese Journal of Engineering*, 2013, Article ID 624658, 9 pages, 2013. DOI: `10.1155/2013/624658`. | Related coupled-beam/frame source with an internal elastic hinge and elastic boundary conditions. Useful for future joint-flexibility planning and for comparison against FEM/measurements, but it uses Bernoulli--Euler theory rather than Timoshenko theory. | Same bibliographic source as the older local `работа2.pdf`; not a source for circular-rod shear coefficient. |
| T1 | `Т1.pdf` | A. Diaz-de-Anda, J. Flores, L. Gutierrez, R. A. Mendez-Sanchez, G. Monsivais, and A. Morales, "Experimental study of the Timoshenko beam theory predictions," *Journal of Sound and Vibration*, 331, 5732--5744, 2012. DOI: `10.1016/j.jsv.2012.07.041`. | Main experimental/theoretical validation source for Timoshenko predictions, cylindrical rods, rectangular beams, 3-D FEM comparison, critical frequency, and the second Timoshenko spectrum. | The paper explicitly discusses cylindrical rods and rectangular beams under free-free boundary conditions. |
| T2 | `Т2.pdf` | A. Diaz-de-Anda, A. Pimentel, J. Flores, A. Morales, L. Gutierrez, and R. A. Mendez-Sanchez, "Locally periodic Timoshenko rod: Experiment and theory," *Journal of the Acoustical Society of America*, 117(5), 2814--2819, 2005. DOI: `10.1121/1.1880732`. | Primary local source for circular/cylindrical aluminum rods modeled with Timoshenko beam theory and the transfer matrix method. Also provides the working circular-section shear coefficient used below. | The PDF states `k = (6 + 12*nu + 6*nu^2)/(7 + 12*nu + 4*nu^2)` and uses `k = 0.925` for `nu = 0.3`. |
| T3 | `Т3.pdf` | J. A. Franco-Villafane and R. A. Mendez-Sanchez, "On the accuracy of the Timoshenko beam theory above the critical frequency: best shear coefficient," arXiv:1405.4885v2 `[physics.class-ph]`, 2014; preprint submitted to Elsevier, May 26, 2014. DOI: not found in local PDF. | Main source in this local set for non-uniqueness and best-fit shear coefficients, including one-coefficient and two-coefficient Timoshenko variants and below/above-critical-frequency fits. | The text cites Kaneko, Cowper, Hutchinson, and related coefficient formulas. |
| T4 | `Т4.pdf` | N. G. Stephen, "On 'A check on the accuracy of Timoshenko's beam theory'," *Journal of Sound and Vibration*, 257(4), 809--812, 2002. DOI: not found in local PDF. | Short critique/source note on shear-coefficient accuracy, second-spectrum behavior, and why coefficient choice is not unique. Useful as a caution against treating any one `k` as universal. | Mentions Cowper, Hutchinson, and two-coefficient theory; not a circular-rod coefficient source. |

## Source roles

### 1. Timoshenko theory for circular / cylindrical rods

- **T2, Diaz-de-Anda et al. 2005.** This is the primary local source for
  circular/cylindrical rods. It studies a locally periodic Timoshenko rod from
  experiment and theory, uses aluminum rods of circular cross-section, and
  computes normal-mode frequencies and amplitudes with Timoshenko beam theory
  and the transfer matrix method. For this project it is the cleanest source
  for a diagnostic circular-rod Timoshenko model.
- **T1, Diaz-de-Anda et al. 2012.** This source studies cylindrical rods and
  rectangular beams experimentally and theoretically, with FEM comparison. It
  is important for checking how Timoshenko predictions behave near and above
  the critical frequency, not only in a low-frequency circular-rod regime.

### 2. Experimental validation of Timoshenko beam theory

- **T1, Diaz-de-Anda et al. 2012.** This is the main validation source because
  it compares measured resonances and wave amplitudes with Timoshenko theory
  and 3-D FEM. It supports using TBT as a diagnostic model below the critical
  frequency and documents what becomes delicate near the second spectrum.
- **T2, Diaz-de-Anda et al. 2005.** This paper validates the transfer-matrix
  Timoshenko calculation against EMAT measurements for locally periodic rods.
  It is especially useful because the tested specimens are aluminum rods and
  the shear coefficient is explicitly chosen for circular cross-sections.
- **A2-2, Howson and Williams 1973.** This source compares theory and
  experiment for an H-shaped frame with Timoshenko members. It is useful for
  future frame validation, but it should not be used to choose the circular-rod
  shear coefficient.

### 3. Shear correction coefficient

- **T2, Diaz-de-Anda et al. 2005.** The source explicitly gives the circular
  coefficient
  `k = (6 + 12*nu + 6*nu^2)/(7 + 12*nu + 4*nu^2)` and uses `k = 0.925` for
  `nu = 0.3`, a typical aluminum Poisson ratio. This is the selected baseline
  coefficient for future diagnostic Timoshenko models of circular rods unless
  the project later records a different choice.
- **A2-1, Kramer and Gfrerer 2024.** This paper uses `k = 5/6` as the
  rectangular-section shear coefficient in a Timoshenko-Ehrenfest frame
  setting. It is useful for frame formulation, but the coefficient belongs to
  a rectangular cross-section and must not be transferred to circular rods.
- **T3 and T4.** These sources are useful because they show that the shear
  coefficient is a modeling parameter rather than a universal constant. They
  compare several formulas and discuss one-coefficient versus two-coefficient
  versions.

### 4. Best shear coefficient / non-uniqueness

- **T3, Franco-Villafane and Mendez-Sanchez 2014.** This is the clearest local
  source for the "best coefficient" issue. It compares Timoshenko predictions
  with experimental data and finds that different coefficients can be optimal
  below and above the critical frequency, including a two-coefficient version.
- **T4, Stephen 2002.** This note explicitly treats the accuracy of different
  Timoshenko coefficient choices and argues against a naive universal
  interpretation. It cites Cowper and Hutchinson and is a useful warning source
  when selecting `kappa`.
- **Cowper / Hutchinson / Kaneko.** Separate PDFs for these works were not
  among the eight files, but they are cited and discussed inside T3/T4. For
  this project, `kappa` should be recorded as a modeling parameter; different
  formulas and even two-coefficient theories exist, and the T2 circular-rod
  value is only the current baseline diagnostic choice. No Cowper-like
  circular formula has been promoted here to a project diagnostic option yet.

### 5. Critical frequency / second Timoshenko spectrum / valid frequency range

- **T1, Diaz-de-Anda et al. 2012.** This is the main local experimental source
  for the critical-frequency and second-spectrum discussion. It reports
  measurements below and above the critical frequency and compares TBT, FEM,
  and experiment.
- **T3, Franco-Villafane and Mendez-Sanchez 2014.** This source explicitly
  fits coefficients below and above the critical frequency. It is useful for
  planning any future diagnostic that asks whether a single `k` remains
  adequate outside the low-frequency regime.
- **T4, Stephen 2002.** This note is a compact cautionary source about the
  second spectrum and the accuracy of TBT at short wavelengths. It should be
  read before making strong claims about the valid frequency range.
- **Current diagnostic convention.** The local diagnostic scripts record the
  Timoshenko cut-off
  `Omega_c = sqrt(kappa*G*A/(rho*I))`. For the coupled equal-rod normalization
  `Omega = epsilon*Lambda^2`, `r = 2*epsilon`, and `L_segment = 1`, this gives
  `Lambda_c = (kappa/(2*(1 + nu)))**0.25/epsilon`.

### 6. Timoshenko theory for frames

- **A2-1, Kramer and Gfrerer 2024.** This is the most direct modern
  frame-oriented Timoshenko source in the set. It gives a numerical assembly
  technique for plane frames based on Timoshenko-Ehrenfest beam theory and is
  relevant to future coupled-beam/frame formulations.
- **A2-2, Howson and Williams 1973.** This is the classic frame-oriented
  source with axially loaded Timoshenko members and dynamic stiffness assembly.
  It is useful for future frame benchmarks, but it is not a circular-rod
  coefficient source.
- **A2-3 and A2-4.** These are frame/coupled-beam background sources rather
  than Timoshenko-frame sources. They are still useful for historical frame
  methods, internal hinges, elastic boundary conditions, and comparison against
  exact/FEM/experimental results.
