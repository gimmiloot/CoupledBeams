# Yartsev 2024 — Source Note and Local Scan Map

## 1. Bibliographic record

Б. А. Ярцев. *Связанные колебания композитных конструкций*: монография.
Санкт-Петербург: ФГУП «Крыловский государственный научный центр», 2024.
216 с.: ил. ISBN `978-5-6048511-4-2`. УДК `534.12:678.067`. ББК `35.719`.
Язык: русский.

The data were checked directly against the first bibliographic page of the
local Chapter-1 fragment. The book prints the place as the abbreviation
`СПб.`; the bibliography normalizes it to `Санкт-Петербург`. No web metadata
were used. The project citation key is
`yartsev_2024_coupled_composite_structures`.

## 2. Role in the project

This monograph is the primary literature source for the `anisotropic_rods`
direction. It supplies the material-transformation background and the selected
Chapter-2 single-rod model: one rectangular monoclinic Timoshenko rod with
bending--torsion coupling, checked in free-free and cantilever configurations.
It does not replace the verified isotropic theory, determinants, root solvers,
or FEM baseline.

All seven local PDFs listed below are fragments of this one book. They have one
BibTeX record and one citation key. The PDFs are currently local untracked
source material and may be absent from a public clone.

## 3. Local PDF-fragment map

| Local relative path | Size | PDF pages | Printed book pages | Contents | Status and overlap | SHA256 |
| --- | ---: | ---: | --- | --- | --- | --- |
| `docs/literature/pdf/Глава 1_compressed.pdf` | 7,159,712 B (6.828 MiB) | 23 | 2--5, 11--51 | Bibliographic page; contents; start of the preface; Chapter 1 on elasticity, viscoelasticity and transformed material properties; opening of Chapter 2 | `current; incomplete; overlap`: pages 6--10 are absent; pages 50--51 repeat in the Chapter-2 fragment | `30d0a703c444f65e5b462158fdfac9ecc06636fa50a86d38fc9366dfbc4f4733` |
| `docs/literature/pdf/Глава 2_compressed.pdf` | 7,123,779 B (6.794 MiB) | 17 | 50--83 | End of Chapter 1; complete Chapter 2; opening of Chapter 3 | `current; complete fragment; overlap`: 50--51 repeat from Chapter 1 and 82--83 repeat in Chapter 3 part 1 | `c3aff7b82f5d25a92447fad0658f256c68582849aea8627c8894744ab35bdd83` |
| `docs/literature/pdf/Глава 3 Часть 1.pdf` | 8,353,196 B (7.966 MiB) | 20 | 82--121 | End of Chapter 2; Chapter 3 through page 121, including multilayer-plate formulation and numerical experiments | `current; complete fragment; overlap`: 82--83 repeat from Chapter 2; continuation starts at 122 without a gap | `5d15d307069220a492d2e1e4966c1f50a3e499c496eafe4036b08e98dc64c7c6` |
| `docs/literature/pdf/Глава 3 Часть 2_compressed.pdf` | 3,397,002 B (3.240 MiB) | 9 | 122--137 | Continuation and end of Chapter 3; opening of Chapter 4 | `current; complete fragment; overlap`: no gap after 121; 136--137 repeat in the Chapter-4 fragment | `a5a82879f3d01322abfa52ba50b986319136628eea22ddc0e9ee61ba7e86ebcb` |
| `docs/literature/pdf/Глава 4_compressed.pdf` | 7,338,240 B (6.998 MiB) | 20 | 136--175 | End of Chapter 3; complete Chapter 4 on closed-profile thin-walled rods; opening of Chapter 5 | `current; complete fragment; overlap`: 136--137 repeat from Chapter 3 and 174--175 repeat in the applications fragment | `e96e5157ecd2e9e4efcc13bbb40c48c29903630f891b549ed023a44ad900c50b` |
| `docs/literature/pdf/Применения_compressed.pdf` | 6,342,448 B (6.049 MiB) | 12 | 174--197 | End of Chapter 4; complete Chapter 5, “Применение в технике”; opening of the bibliography | `current; complete fragment; overlap`: 174--175 repeat from Chapter 4 and 196--197 repeat in the bibliography fragment | `20b9063b52d4f572fba6bec5e344fc6bbb2900728080aa80136b895be64b5d3e` |
| `docs/literature/pdf/Литература_compressed.pdf` | 5,198,924 B (4.958 MiB) | 10 | 197--215 | Complete bibliography on 197--212; subject index on 213--215 | `current; incomplete; overlap`: page 197 repeats from Chapter 5; declared volume page 216 is absent | `5220c09f839d68a6570536e2cce917a4b13d77b207ede15e6a9b93e55f5341f8` |

The local set therefore covers printed pages 2--5 and 11--215. It has no gaps
from Chapter 1 page 11 through the end of the subject index on page 215. The
known gaps are front-matter pages 6--10 and the declared final volume page 216.
The latter may be an unnumbered or publication page, but its content cannot be
inferred from the available scan. No superseded Chapter-3 part-1 duplicate was
present in this inventory.

All 111 PDF pages opened and rendered successfully. Low-resolution review of
every rendered page, supplemented by pixel-statistics screening, found no
fully white or visibly damaged page. This is a scan-integrity check, not OCR or
a proof that every small printed symbol is legible.

## 4. Chapter map and applicability to CoupledBeams

| Chapter | Content | Project applicability |
| --- | --- | --- |
| 1. Элементы механики деформируемого твёрдого тела | Linear elasticity and viscoelasticity; complex moduli; transformation of a unidirectional layer; mutual-influence coefficients; material parameters in Table 1.2 | Material and damping definitions used by the current Chapter-2 source line. Critical pages: 24--25 for (1.32), (1.34); 30--31 for (1.41), (1.42); 40--46 for (1.50)--(1.56) and Table 1.2. |
| 2. Моноклинный стержень | Rectangular monoclinic rod; Timoshenko bending plus generalized torsion; material bending--torsion interaction; free-free and cantilever rods; experimental identification | Current source-reproduction model. Critical pages: 52--55 for (2.1)--(2.18), 56--57 for the free-free specimen and Figure 2.2, and 64--68 for the cantilever study and Figure 2.8. |
| 3. Многослойные пластины | Multilayer-plate theories and numerical experiments | Relevant to the broader coupled-vibration theory of composites, but not the source model for the first rod stage. |
| 4. Тонкостенные стержни замкнутого профиля | Spatial theory with longitudinal motion, two bending directions, shear, torsion and warping; begins on page 137 | Possible more general future stage. It is not a direct replacement for Chapter 2 and is not the current model. |
| 5. Применение в технике | Adaptive and vibration-absorbing structures; engineering use and control of coupling | Engineering context and future design motivation, not a current computational baseline. |

## 5. Current selected source line

The selected line is deliberately narrow:

- Chapter 2;
- elastic specialization for equation and boundary gates, and the book's
  complex material description for source reproduction;
- Timoshenko bending plus generalized torsion;
- one rectangular monoclinic rod in free-free and cantilever configurations;
- corrected internally consistent equation variants;
- no coupled rods, coupled-joint conditions, or production anisotropic API.

Euler--Bernoulli comparison, Saint-Venant replacement, coupled rods, and FEM
validation have not been implemented for this direction.

## 6. Internal-consistency check of the formulas

The source-reproduction stage kept the literally printed and internally
consistent variants separate.

- Equation (2.1) is printed without `I_y` in the rotary-inertia term. The
  `state_corrected` form restores the dimensionally required
  `rho * I_y * psi_tt`.
- The signs printed for `d0` and `f0` after (2.16) do not recover an independent
  positive torsional spectrum at the orthotropic endpoints
  `theta = 0°` and `90°`.
- `eliminated_corrected` uses positive `d0` and `f0`. Its first eight positive
  elastic roots agree with the corrected state system to a maximum relative
  difference of `1.408351e-09` (about `1.4e-9`).

These are recorded source-specific corrections, not silent changes to the
book, the verified isotropic theory, or any production model.

## 7. Reproduction evidence

The existing diagnostic evidence gives
`PASS_WITHIN_GRAPH_RESOLUTION` for the calculated solid curves in Figure 2.2.
The accepted pair is `state_corrected` with restored `I_y` and
`eliminated_corrected` with positive `d0`, `f0`. The comparison used a separate
rounded graphical reading with declared uncertainties of `±0.06 kHz` and
`±0.06` in `eta*1e2`; all 98 calculated comparisons were within those bounds.
No material parameter was fitted.

The book supplies no digital curve table. Experimental markers were not
digitized and are not covered by this status.

The subsequent cantilever source gate used the same corrected Chapter-2 model.
The frequency panel of Figure 2.8 alone was ambiguous at the scan resolution,
but the saved frequency and loss-factor curves together gave
`BOOK_SLOPE_CLAMP_CONFIRMED`. Thus the calculated Figure-2.8 curves use
`w(0)=0`, `w'(0)=0`, `Phi(0)=0`. Existing complex roots were reused in this
decision: no determinant evaluation, root solve, continuation, or scientific
matrix exponential was run, and no parameter was fitted.

## 8. Known source warnings

1. Do not use the literal (2.1) inertia term as a physical baseline without
   restoring `I_y`.
2. Do not use the printed `d0`/`f0` signs as the current physical baseline;
   that variant remains diagnostic only.
3. Do not describe `PASS_WITHIN_GRAPH_RESOLUTION` as exact numerical agreement
   or experimental validation.
4. Do not present Chapter 4 as a direct replacement for Chapter 2.
5. The Figure-2.8 source baseline is `book_slope_clamp`. The alternative
   `timoshenko_section_clamp`, which sets section rotation rather than
   centerline slope to zero, remains physically meaningful but is not the
   source baseline.
6. Do not transfer the external source clamp automatically to a future
   internal rigid joint. Coupled-rod compatibility and vector equilibrium must
   be derived and verified separately.
7. Do not promote this one-rod diagnostic helper to a coupled-rod or production
   anisotropic API without a separate model, limiting-case contract, and
   verification stage.

## 9. Scan and research limitations

- The scan has no text layer; no whole-book OCR was performed.
- Front-matter pages 6--10 and declared volume page 216 are absent.
- Several neighboring fragments intentionally overlap by one two-page spread.
- The local PDFs are untracked and may be absent from the public repository.
- The scan inventory registers the monograph only. Its 239-item bibliography
  was not imported into the project bibliography.
- Current evidence covers the Chapter-2 free-free Figure-2.2 benchmark and the
  cantilever Figure-2.8 source-boundary decision. It does not cover coupled
  rods or internal-joint conditions.

## 10. Project links

- [BibTeX bibliography](../literature/bibliography.bib) — entry
  `yartsev_2024_coupled_composite_structures`.
- [Source-index entry](../literature/source_index.md#yartsev_2024_coupled_composite_structures).
- [Chapter-2 single-rod reproduction note](yartsev_ch2_single_rod_reproduction.md).
- [Chapter-2 cantilever reproduction note](yartsev_ch2_cantilever_reproduction.md).
- [Free-free diagnostic entry point](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py).
- [Cantilever diagnostic entry point](../../scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py).
- Local generated [result report](../../results/anisotropic_rods/yartsev_ch2_free_free/single_rod_reproduction_report.md) and directory `results/anisotropic_rods/yartsev_ch2_free_free/`; generated evidence is Git-ignored and may be absent from a public clone.
- Local cantilever evidence is under `results/anisotropic_rods/yartsev_ch2_cantilever/`, `results/anisotropic_rods/yartsev_ch2_cantilever_quick_gate/`, and `results/anisotropic_rods/yartsev_ch2_cantilever_boundary_source_check/`; it may likewise be absent from a public clone.
