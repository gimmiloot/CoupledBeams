# Mu Sweep Slow Evolution Assessment

## Question

Does the existing `beta = 15 deg`, `r = 5 mm` `mu` sweep support strict
veering, or is the safer interpretation that the observed branch
reorganization is slow evolution / quasi-degeneracy-like behaviour / exchange
of modal character?

The main target is `bending_desc_04`. `bending_desc_01` and
`bending_desc_02` are used as comparison cases.

No new FE solves, eigenvalue solves, solver-code changes, source-data edits,
figure edits, or manuscript edits were made.

## Data Used

- `docs/veering/README.md`
- `docs/veering/literature_assessment.md`
- `docs/veering/terminology.md`
- `docs/veering/candidate_cases.md`
- `docs/veering/strict_veering_assessment.md`
- `docs/veering/bending_desc_04_pair_assessment.md`
- `docs/veering/data/target_descendants_beta15_r5_analysis.csv`
- `docs/veering/data/target_descendants_beta15_r5_selected_points.csv`
- `docs/veering/data/flat_mu_bending_energy_cases.csv`
- `docs/veering/data/paired_branch_candidates_beta15_r5.csv`
- `docs/veering/data/paired_gap_scan_beta15_r5.csv`
- `docs/veering/data/strict_veering_decision_table_beta15_r5.csv`
- `docs/veering/data/bending_desc_04_pair_search_beta15_r5.csv`
- `results/branchwise_fem_grid.csv`
- `results/branchwise_fem_near_crossings.csv`
- `results/branchwise_fem_min_gap_map.csv`
- `docs/literature/notes/nair_1973_quasi_degeneracies.md`
- `docs/literature/notes/manconi_2017_veering_strong_coupling.md`
- `docs/literature/notes/ehrhardt_2018_clamped_beam_veering.md`
- `docs/literature/notes/lacarbonara_2005_imperfect_beams_veering.md`
- `docs/literature/notes/fontanela_2021_nonlinear_localisation_coupled_beams.md`

The derived decision table for this assessment is:

- `docs/veering/data/mu_sweep_interpretation_beta15_r5.csv`

## Operational Criteria

### `strict veering`

Use this label only when all of the following are present:

- a concrete competing pair of tracked branches;
- a genuinely small paired gap;
- an interior local gap minimum, not only a boundary minimum;
- evidence that the interaction is localized in a narrow `mu` window;
- pairwise exchange evidence, such as pairwise MAC / cross-MAC or an
  equivalent two-branch mode-shape exchange diagnostic;
- the modal-character change is tied to that branch-pair interaction, not only
  observed along one tracked branch.

### `quasi-degeneracy-like`

Use this label when the data show notable modal reorganization or close
spectral context, but the strict-veering pair evidence is missing, weak, or
not tied to a narrow avoided-crossing window. This follows the caution in
`nair_1973_quasi_degeneracies`: rapid or notable modal change is not enough by
itself to prove strict veering.

### `exchange of modal character without strict veering`

Use this label when a tracked branch preserves its continuation identity but
changes local modal descriptors, such as arm-wise half-wave order or
localization, and no branch-pair exchange has been demonstrated. This label is
tracked-branch based and does not rely on sorted-spectrum relabeling.

### `slow evolution compatible with comparatively stronger coupling`

Use this label when the available evidence fits broad or staged evolution
better than rapid veering:

- the paired gap is not especially small;
- no convincing narrow interaction window is present;
- modal descriptors change gradually or in staged intervals over `mu`;
- pairwise exchange evidence is unavailable;
- localization, if present, is not shown to coincide with a small paired gap.

This is only a compatibility statement with the Manconi--Mace caution that
rapid veering is associated with weak coupling while stronger coupling can give
slow eigenvalue-locus evolution. It does not prove strong coupling for
CoupledBeams, because no explicit Manconi-style coupling/skeleton model has
been constructed here.

## `bending_desc_04`

`bending_desc_04` remains the strongest current candidate for
modal-character reorganization, but not for strict veering.

The tracked local modal order changes in stages:

- `(2,2)` for `0 <= mu <= 0.21375`;
- `(1,2)` for `0.21625 <= mu <= 0.28375`;
- `(1,3)` for `0.28625 <= mu <= 0.9`.

The minimum stored tracking MAC along this branch is about
`0.986810080707937`, so the tracked identity is preserved in the available
continuation data.

The best available competing-pair evidence is not sufficient for strict
veering:

- the closest numerical tracked bending neighbour in
  `bending_desc_04_pair_search_beta15_r5.csv` is `bending_desc_03`, with
  `gap_min_abs = 0.591623872049765` and
  `gap_min_rel = 0.0803359182206152`, but the minimum occurs at the boundary
  `mu = 0` and the gap grows immediately;
- the more relevant reorganization-window candidate is `bending_desc_01`,
  with an interior minimum near `mu = 0.26625`, but the relative gap is still
  large: `gap_min_rel = 0.123295265709942`;
- no pairwise MAC / cross-MAC or equivalent branch-pair exchange evidence is
  available in the stored data.

The localization evidence supports modal-character reorganization, but not a
strict pair interaction. For `r = 5 mm`, the stored bending energy split is:

- `mu = 0`: left/right bending energy approximately `0.50 / 0.50`;
- `mu = 0.5`: left/right approximately `0.112 / 0.888`;
- `mu = 0.9`: left/right approximately `0.048 / 0.952`.

This is a strong right-arm localization tendency, but it is sampled only at
three `mu` values and is not tied to a small paired-gap minimum.

Conclusion for `bending_desc_04`: strict veering is not supported. The most
defensible current description is staged exchange of modal character without
strict veering, with quasi-degeneracy-like reorganization and slow evolution
compatible with comparatively stronger coupling. The last phrase is
interpretive caution, not proof of strong coupling.

## Comparison With `bending_desc_01` And `bending_desc_02`

`bending_desc_01` is a real modal-character-change case, but weaker than
`bending_desc_04` for the slow-evolution / quasi-degeneracy interpretation.
Its local order changes from `(1,1)` to `(1,2)` between the sampled points
`mu = 0.04625` and `mu = 0.04875`, while its nonmonotonic `Lambda(mu)` curve
has a maximum near `mu = 0.285`. The available pair scan with
`bending_desc_04` has an interior minimum near `mu = 0.26625`, but the gap is
large: `gap_min_rel = 0.123295265709942`. The early descriptor change is not
shown to be a narrow branch-pair exchange. The conservative label remains
`exchange of modal character without strict veering`.

`bending_desc_02` behaves as the control case. Its local order remains `(1,1)`
over the full `0 <= mu <= 0.9` range, and the comparison against
`bending_desc_01` has a large boundary minimum:
`gap_min_rel = 0.17776010462112982` at `mu = 0`. Existing descriptors support
smooth evolution / no veering.

Relative to these controls, `bending_desc_04` is still the strongest
quasi-degeneracy-like candidate because it has staged local-order
reorganization and stored arm-energy localization evidence, even though it
lacks the paired evidence required for strict veering.

## Strict Veering: What Is Missing

The current evidence lacks:

- a genuinely close competing tracked pair for `bending_desc_04`;
- a small interior paired-gap minimum;
- a narrow `mu` interaction window tied to the modal-character change;
- pairwise MAC / cross-MAC or equivalent two-branch exchange evidence;
- localization metrics for both branches in a suspected pair;
- an explicit uncoupled-blocked/skeleton construction that would justify a
  direct Manconi-style coupling classification.

## Is Slow-Evolution / Quasi-Degeneracy Interpretation More Consistent?

Yes, for the current evidence set and with conservative wording.

For `bending_desc_04`, the key observations are broad/staged modal-order
change, right-arm localization tendency at larger `mu`, no small interior
paired gap, and no branch-pair exchange evidence. That combination is more
consistent with slow evolution / quasi-degeneracy-like reorganization than
with strict veering.

The Manconi--Mace literature note supports this cautious framing because rapid
veering is expected in weak-coupling settings, while stronger coupling can
lead to slower eigenvalue-locus evolution. However, the project has not
defined or measured the corresponding coupling parameter, so the correct
wording is:

`consistent with slow evolution under comparatively stronger coupling`

not:

`strong coupling is proved`.

## Most Defensible Current Wording For The Project

For `bending_desc_04`:

`The existing beta = 15 deg, r = 5 mm, mu-sweep evidence does not support a
strict veering claim. The tracked branch shows staged exchange of modal
character and a right-arm localization tendency. Because no close competing
tracked pair, small interior paired gap, narrow interaction window, or pairwise
mode-shape exchange evidence is available, the behaviour is better described
as quasi-degeneracy-like modal-character reorganization / slow evolution,
consistent with the Manconi--Mace caution that rapid veering requires weak
coupling. This does not prove strong coupling in the present model.`

For `bending_desc_01`:

`exchange of modal character without strict veering`

For `bending_desc_02`:

`smooth evolution / no veering control`

## What Additional Data Would Be Needed To Overturn This Conclusion

To support strict veering in a future revision, the project would need:

- stored mode vectors or mode-shape descriptors across the suspected
  interaction windows;
- pairwise MAC / cross-MAC tables for the candidate branch pairs;
- a confirmed competing tracked branch with a genuinely small interior gap;
- evidence that modal-character exchange occurs between the two branches, not
  only along one tracked branch;
- localization or arm-energy metrics for both branches in the suspected pair;
- a narrow-window analysis showing approach, minimum gap, and separation;
- if using Manconi--Mace language strongly, an explicit
  uncoupled-blocked/skeleton or coupling-order model for CoupledBeams.
