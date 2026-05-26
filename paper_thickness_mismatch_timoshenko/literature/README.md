# Literature

Article-specific literature notes for the planned thickness-mismatch /
Timoshenko work belong here. Keep source-specific warnings in the repository
literature notes when they affect the whole project.

## Candidate Timoshenko and shear-coefficient sources

The project-wide note for these local PDFs is
`docs/literature/timoshenko_shear_sources.md`.

- Circular rods: `Т2.pdf` (Diaz-de-Anda et al. 2005) is the primary source for
  locally periodic circular aluminum rods modeled with Timoshenko beam theory
  and the transfer matrix method. `Т1.pdf` (Diaz-de-Anda et al. 2012) extends
  the experimental/TBT/FEM discussion to cylindrical rods and rectangular
  beams.
- Shear coefficient: `Т2.pdf` records the circular-section working choice
  `k = (6 + 12*nu + 6*nu^2)/(7 + 12*nu + 4*nu^2)`, with `k = 0.925` for
  `nu = 0.3`. `A2-1.pdf` uses `k = 5/6` for a rectangular frame member and is
  not a circular-rod coefficient source.
- Experimental validation: `Т1.pdf` and `Т2.pdf` compare Timoshenko predictions
  with measurements; `A2-2.pdf` compares a Timoshenko-member frame model with
  H-frame experiments.
- Frames: `A2-1.pdf` and `A2-2.pdf` are the main Timoshenko/frame sources.
  `A2-3.pdf` and `A2-4.pdf` are related frame/coupled-beam background sources
  but do not set the Timoshenko circular-rod coefficient.
- Critical frequency / second spectrum: `Т1.pdf`, `Т3.pdf`, and `Т4.pdf` are
  the main local sources for the critical frequency, second Timoshenko
  spectrum, valid range, and non-unique/best-fit coefficient discussion.
