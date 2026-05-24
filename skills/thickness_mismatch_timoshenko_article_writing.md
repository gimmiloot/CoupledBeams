# Thickness-Mismatch Timoshenko Article Writing

Use this skill for the planned article about how rod-thickness mismatch affects
the coupled-beam spectrum and where Timoshenko theory or validation is needed.

## Scope

- The article concerns two coupled circular rods with different radii.
- The current Euler-Bernoulli / thin-rod eta model is a baseline diagnostic
  model, not the final theory for thick or questionable regimes.
- Questionable regimes are identified by the project applicability rule
  `2*r_i/l_i > 0.1` for at least one rod.
- Planned article-ready evidence should include theory/FEM/Timoshenko
  comparison where the thin-rod model is questionable.

## Branch And Mode Language

- Describe branch identity as a descendant branch: the continuation of a mode
  shape selected at the tracking start point.
- Do not describe a branch as merely the k-th sorted frequency at each
  parameter value.
- Treat `sorted_position` as diagnostic metadata only.
- Do not accept low-MAC assignments or large sorted-position jumps as physical
  conclusions without local refined shape analysis.

## Evidence Standards

- Do not transfer diagnostic claims into article text without the required
  checks recorded in `docs/project_rules.md` and `docs/writing/article_workflow.md`.
- Isolated-rod grids are interpretation aids, not proof of modal identity or
  modal exchange.
- FEM comparison requires a stated normalization and boundary-condition check.
- Claims about veering or modal exchange require tracked descendants, small
  spectral gap evidence, and mode-shape evidence.

## Style

- Use academic article prose, not internal repository-note language.
- Stay close in tone to the current article style when appropriate, but do not
  mechanically copy its structure.
- Keep diagnostic caveats out of the main narrative unless they are promoted to
  article-ready limitations or validation statements.
