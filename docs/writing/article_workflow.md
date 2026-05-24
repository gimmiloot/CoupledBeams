# Article Workflow

This document separates diagnostic work from article-ready material.

## Diagnostic Results

Diagnostic outputs are exploratory. They may show useful patterns, warnings,
failed assignments, or checks that still need refinement. They should remain in
`results/` unless a task explicitly promotes them to article material.

Diagnostic plots must not overwrite article figures. Diagnostic conclusions
should not be copied into article prose until the checks below are complete.

## Promotion Checklist

Before a diagnostic result becomes article material:

- confirm that the plotted branches are descendant branches or explicitly label
  the plot as sorted frequencies;
- check whether any low-MAC, low-margin, or suspicious assignment is present;
- check thin-rod applicability with the project criterion `2*r_i/l_i <= 0.1`
  for each circular rod;
- verify the relevant limiting case, such as eta=0 for the
  thickness-mismatch extension;
- state the normalization used for any FEM comparison;
- verify that reference curves use the intended boundary conditions;
- confirm that the figure is reproducible from a documented command;
- move or regenerate the selected figure into the article figure workflow only
  after the article task explicitly requests it.

## Figure Selection Workflow

Candidate figures should be listed in the article notes first. Each candidate
needs:

- source script and output path;
- parameter set;
- branch identity convention;
- applicability warnings;
- FEM or shape-evidence status, when relevant;
- decision status: diagnostic-only, candidate, accepted, or rejected.

## TODO Placeholder Policy

New article skeletons may contain TODO markers for title, authors, abstract,
section content, bibliography, and final figure selection. Do not fill these
with invented scientific claims. Replace TODO markers only after the relevant
diagnostic result has passed the promotion checklist.
