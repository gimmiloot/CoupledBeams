# Strict Veering Assessment

## Scope

This assessment covers the existing `CoupledBeams` results for:

- `beta = 15 deg`
- `r = 5 mm`
- `0 <= mu <= 0.9`
- candidate branches `bending_desc_01`, `bending_desc_02`, and
  `bending_desc_04`

No new FE solves, eigenvalue solves, solver-code changes, or article edits
were made. The outputs in this folder are post-processing artifacts built from
existing CSV/PNG results.

## Data Used

- `docs/veering/data/target_descendants_beta15_r5_analysis.csv`
  - tracked histories for `bending_desc_01`, `bending_desc_02`, and
    `bending_desc_04`
  - includes `Lambda`, current sorted position, local half-wave counts,
    nearest single-rod reference families, and along-track `tracking_mac`
- `docs/veering/data/target_descendants_beta15_r5_selected_points.csv`
  - selected qualitative states for the three target branches
- `docs/veering/data/flat_mu_bending_energy_cases.csv`
  - existing arm-energy split; among the current target branches it directly
    covers `bending_desc_04`
- `results/branchwise_fem_grid.csv`
  - existing tracked-branch FEM grid used to check candidate competing
    branches and sorted-spectrum positions
- `results/branchwise_fem_near_crossings.csv`
  - existing near-crossing scan used to check whether the target low branches
    already have a documented close tracked neighbor
- `results/branchwise_fem_min_gap_map.csv`
  - existing minimum-gap map used as context for where small gaps occur
- Existing figures copied into `docs/veering/figures/`

No stored eigenvector arrays or branch-pair mode-shape histories were found in
the existing result files. Therefore true pairwise MAC exchange evidence could
not be computed without running new FE solves.

## Methodological Rule: Tracked Branches Vs Sorted Spectrum

All veering decisions here are made at the level of tracked branch identities.
The current position of a branch in the sorted spectrum is recorded as a
descriptor only. A change in sorted position is not treated as evidence of
veering by itself.

For `bending_desc_01` and `bending_desc_04`, the paired comparisons are
therefore conservative gap-scan candidates, not confirmed competing branches.
For `bending_desc_02`, the paired comparison is a control comparison.

## Criteria For Strict Veering

The label `strict veering` is allowed only if all of the following are present:

- a concrete pair of tracked branches;
- a clear local minimum of the gap between those tracked branches;
- no evidence of a true crossing;
- pairwise evidence of exchange between the two branches;
- the modal-character change is tied to the small-gap region.

The current data do not satisfy these criteria for any of the three candidate
branches.

## Candidate Assessment

### `bending_desc_01`

`bending_desc_01` has a real modal-character change along the tracked branch:
the local order changes from `(1,1)` to `(1,2)`, and the branch has a
nonmonotonic `Lambda(mu)` curve.

The nearest available low tracked comparison used here is `bending_desc_04`
over `mu = 0.20..0.35`. The minimum gap in that scan is not small:

- `gap_min_abs = 0.721047903846701`
- `gap_min_rel = 0.123295265709942`

No pairwise MAC exchange evidence is available from existing stored data.
Therefore this branch should not be called strict veering. The conservative
label is:

`exchange of modal character without strict veering`

### `bending_desc_02`

`bending_desc_02` is the control case. Its local modal order remains `(1,1)`
over the tracked range, and the available descriptors do not show modal
exchange.

The control comparison against `bending_desc_01` gives a large minimum gap:

- `gap_min_abs = 0.766082111455375`
- `gap_min_rel = 0.177760104621130`

No pairwise MAC is available, but the existing descriptors already support the
negative conclusion. The conservative label is:

`smooth evolution / no veering`

### `bending_desc_04`

`bending_desc_04` has the strongest modal reorganization among the current
candidate branches. Its local order changes from `(2,2)` through `(1,2)` to
`(1,3)`, and existing energy data for `r = 5 mm` show a strong right-arm
bending-energy bias at large `mu`.

The nearest available low tracked comparison used here is `bending_desc_01`
over `mu = 0.20..0.35`. The minimum gap in that scan is the same pairwise gap
seen from the other side and is still not small:

- `gap_min_abs = 0.721047903846701`
- `gap_min_rel = 0.123295265709942`

The energy split supports localization / arm-bias in this branch, but it is
available only at `mu = 0`, `0.5`, and `0.9`, so it cannot establish
coincidence with a small paired gap. Pairwise MAC exchange evidence is also
missing. The conservative label is:

`exchange of modal character without strict veering`

## What Is Established Confidently

- The three target branches are tracked branches, not merely sorted-spectrum
  slots.
- `bending_desc_02` behaves as a useful no-veering control in the current
  evidence set.
- `bending_desc_01` and `bending_desc_04` undergo modal-character changes
  along their tracked histories.
- `bending_desc_04` has existing evidence of arm-energy redistribution /
  localization tendency.
- No current target branch has enough paired evidence for `strict veering`.

## What Remains Unproven

- A confirmed competing tracked branch for `bending_desc_01` or
  `bending_desc_04`.
- A small avoided-crossing gap tied to the modal-character change.
- Pairwise exchange of modal character between two tracked branches.
- Coincidence between localization and a small paired gap.

## Final Labels

- `bending_desc_01`: `exchange of modal character without strict veering`
- `bending_desc_02`: `smooth evolution / no veering`
- `bending_desc_04`: `exchange of modal character without strict veering`

Overall, the present evidence supports modal-character evolution and, for
`bending_desc_04`, localization tendency. It does not support the term
`strict veering`.

## Missing Data Required For Any Stronger Claim

- Stored mode vectors or branch-pair mode-shape descriptors for the relevant
  windows.
- Pairwise MAC / cross-MAC tables across candidate interaction windows.
- Minimum-gap tables for confirmed competing tracked branch pairs.
- Localization metrics for both branches in each suspected pair.
- Evidence that modal-character exchange occurs between two branches, not only
  along one tracked branch.
