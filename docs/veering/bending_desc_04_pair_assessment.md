# Bending Desc 04 Pair Assessment

## Question

Can the existing `beta = 15 deg`, `r = 5 mm` results identify a genuinely
close competing tracked branch for `bending_desc_04` as `mu` varies?

## Data Used

- `docs/veering/data/paired_branch_candidates_beta15_r5.csv`
- `docs/veering/data/paired_gap_scan_beta15_r5.csv`
- `docs/veering/data/pair_mac_beta15_r5.csv`
- `docs/veering/data/target_descendants_beta15_r5_analysis.csv`
- `docs/veering/data/flat_mu_bending_energy_cases.csv`
- `results/branchwise_fem_grid.csv`
- `results/branchwise_fem_near_crossings.csv`
- `results/branchwise_fem_min_gap_map.csv`

No new FE solves, eigenvalue solves, source CSV edits, or solver-code changes
were made.

## Candidate Competing Tracked Pairs

The available tracked-grid scan gives `bending_desc_03` as the closest tracked
bending neighbour of `bending_desc_04`, with
`gap_min_rel = 0.0803359182206152` at `mu = 0`. This is a boundary minimum:
the gap grows immediately as `mu` increases, so it does not define a convincing
avoided-crossing window.

The previously used conservative pair, `bending_desc_04` vs `bending_desc_01`,
has an interior minimum around `mu = 0.27`, with
`gap_min_rel = 0.123295265709942`. This window is more relevant to the staged
local-order reorganization of `bending_desc_04`, but the gap is still large.

The next closest tracked bending branch found in the existing grid is
`bending_desc_06`, with `gap_min_rel = 0.270372762660338`. This is not a close
competitor.

## Gap Comparison

The summary table is
`docs/veering/data/bending_desc_04_pair_search_beta15_r5.csv`.

The best numerical gap is with `bending_desc_03`, but it is not an internal
minimum. The best reorganization-window candidate is `bending_desc_01`, but its
minimum relative gap is about 12 percent in the current `Lambda` scale. The
existing `branchwise_fem_near_crossings.csv` does not list a low-branch close
pair involving `bending_desc_04` at `beta = 15 deg`.

## Why The Best Candidate Is Not Convincing

`bending_desc_03` is close only at the left boundary of the sampled interval.
That is not enough for a strict veering claim because there is no approach,
local minimum, and separation around an interior interaction region.

`bending_desc_01` gives a real interior minimum, but the gap is not small and
the existing data do not include pairwise MAC / cross-MAC evidence showing
exchange between the two tracked branches.

The observed changes in `bending_desc_04` are therefore better supported as
modal-character reorganization along one tracked branch, with localization
tendency from the existing arm-energy data, than as a confirmed paired
avoided-crossing event.

## Can Strict Veering Be Claimed From Existing Data?

No. The current evidence lacks the required combination of:

- a close competing tracked pair;
- a small interior paired-gap minimum;
- pairwise MAC or equivalent cross-branch mode-shape exchange evidence.

The absence of pairwise MAC is a direct limitation of the existing stored data.
No surrogate MAC was introduced.

## Most Defensible Current Label

For `bending_desc_04`, the most defensible current label remains:

`exchange of modal character without strict veering`

With the terminology from `nair_1973_quasi_degeneracies`, it is also reasonable
to describe the broader situation as quasi-degeneracy-like rapid modal
reorganization, but not as strict veering.
