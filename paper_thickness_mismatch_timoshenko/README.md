# Thickness-Mismatch / Timoshenko Article Skeleton

Status: skeleton / planned article.

This folder is reserved for a future article on the influence of thickness
mismatch in coupled circular rods and on regimes where Timoshenko-type theory
or validation may be required.

The current mass-preserving eta model based on Euler-Bernoulli thin-rod theory
is treated here as a baseline diagnostic model. Regions where
`2*r_i/l_i > 0.1` for either rod require extra care, validation, or a
Timoshenko treatment before they can support article-level claims.

Branch identity follows `../docs/project_rules.md`: a branch is a descendant of
a mode shape, and sorted position is diagnostic metadata. Article-ready figures
must be separated from diagnostic outputs in `results/` and promoted only after
the checks in `../docs/writing/article_workflow.md`.

TODO:

- title;
- authors;
- journal;
- UDK/MSC;
- abstract;
- bibliography;
- final figure set;
- validated Timoshenko/FEM comparison plan.

## Structure

- `notes/` -- planning notes, outline, and candidate figures.
- `figures/` -- future article-ready figure files only.
- `results/` -- article-specific processed result tables, not raw diagnostics.
- `tex/` -- LaTeX skeleton and future article source.
- `literature/` -- article-specific literature notes.
