# Dataset Index

Curated copies in this folder were copied from existing project outputs. No
new calculations were run.

## Tables

- `target_descendants_beta15_r5_analysis.csv`
  - Source: `results/target_descendants_beta15_r5_analysis.csv`.
  - Scope: `beta = 15 deg`, `r = 5 mm`, `0 <= mu <= 0.9`.
  - Coverage: `bending_desc_01`, `bending_desc_02`, `bending_desc_04`.
  - Contents: tracked `Lambda`, current sorted position, local half-wave
    counts on left/right arms, nearest clamped-supported reference family,
    tracking MAC, and branch-behavior labels.

- `target_descendants_beta15_r5_selected_points.csv`
  - Source: `results/target_descendants_beta15_r5_selected_points.csv`.
  - Scope: selected qualitative states from the same `beta = 15 deg`,
    `r = 5 mm` analysis.
  - Coverage: three selected points for each of `bending_desc_01`,
    `bending_desc_02`, and `bending_desc_04`.
  - Contents: compact state table for figure interpretation, including
    `Lambda`, sorted position, local half-wave counts, nearest reference
    families, and qualitative labels.

- `flat_mu_bending_energy_cases.csv`
  - Source: `results/flat_mu_bending_energy_cases.csv`.
  - Scope: existing arm-energy cases for flat-`mu` bending descendants.
  - Coverage in this veering line: directly covers `bending_desc_04` among
    the three current candidate branches; the same table also includes
    `bending_desc_05` and `bending_desc_06`.
  - Contents: `beta`, branch id, radius, `mu`, current sorted index,
    frequency in Hz, left/right bending-energy fractions, and axial fraction.

## Candidate Coverage

- `bending_desc_01`: full tracked table and selected states are present; no
  copied arm-energy table specific to this branch is currently available.
- `bending_desc_02`: full tracked table and selected states are present; no
  copied arm-energy table specific to this branch is currently available.
- `bending_desc_04`: full tracked table, selected states, and copied
  arm-energy cases are present.

## Missing For Strict Veering Check

- Paired tracked branch comparison for each suspected interaction.
- Minimum-gap tables between candidate and competing branches.
- Explicit localization metrics for competing branch pairs.
- Branch-pair MAC / modal-correlation tables across suspected exchange
  windows.
