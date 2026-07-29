# Pre-refactor technical snapshot

Captured at local backup timestamp `20260729_094757` (Europe/Moscow). The
working tree was already treated as user-owned even though the initial
tracked/untracked status happened to be clean.

## Git state

| Field | Value |
| --- | --- |
| Repository root | `<repository-root>` |
| Working directory | `<repository-root>` |
| Branch | `main` |
| HEAD | `df701f36569723444e7b131e7cfdad8c894889db` |
| Last commit | `df701f3 Version 0.3.0` |
| Remote | `origin https://github.com/gimmiloot/CoupledBeams` (fetch/push) |
| Tracked files | 300 |
| Modified tracked files | 0 |
| Staged files | 0 |
| Untracked files | 0 |
| Ignored files | 9,459 |
| Pre-existing worktree files, excluding `.git` | 9,759 |

Initial `git status --short` was empty:

```text
(no entries)
```

Initial `git status --ignored --short` reported these top-level ignored
groups. The full output is preserved externally in
`git_status_ignored_short.txt`.

```text
!! .pytest_cache/
!! private/
!! results/_smoke/
!! results/eb_epsilon_apriori_pilot/
!! results/eb_epsilon_apriori_pilot_branch_continuation_v1/
!! results/eb_epsilon_apriori_pilot_complete_spectrum_v1/
!! results/eb_epsilon_baseline_thresholds/
!! results/eb_epsilon_lower_envelope_step3a/
!! results/eb_rule_ab_exact_pareto/
!! results/eb_rule_s_cost_break_even/
!! results/eb_timo_branch_continuation_gateway/
!! results/eb_timo_clean_mode_shapes/
!! results/eb_timo_counterexample_dimensional_frequency_beta/
!! results/eb_timo_general_spectrum_completeness/
!! results/eb_timo_mode_shapes_eps0p03_beta45_eta0_modes4_6/
!! results/eb_validity_fixed_epsilon_geometry_scan/
!! results/eb_validity_vs_timoshenko_stage1/
!! results/eb_vs_timoshenko_3d_validation/
!! results/eb_vs_timoshenko_lambda_beta_cases/
!! results/eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan/
!! results/eb_vs_timoshenko_lambda_mu_cases/
!! results/eb_vs_timoshenko_longitudinal_suspect_modes/
!! results/timoshenko_mode_shape_diagnostics/
!! results/timoshenko_shape_bug_audit/
!! results/timoshenko_shape_construction_audit/
!! results/variable_length_timoshenko_limits_audit.csv
!! results/variable_length_timoshenko_limits_audit.md
!! scripts/**/__pycache__/
!! src/my_project/**/__pycache__/
!! tests/__pycache__/
```

The worktree overlay contains zero files because there were no modified,
staged, or untracked paths. Deleted and renamed path counts were also zero.

## Environment

No package or environment change was made.

| Component | Version |
| --- | --- |
| OS | `Windows-10-10.0.19045-SP0` |
| Python | `3.12.4` (MSC v.1940, 64-bit) |
| pip | `25.0.1` |
| numpy | `2.1.3` |
| scipy | `1.15.2` |
| matplotlib | `3.9.2` |
| mpmath | `1.4.1` |
| pytest | `8.3.4` |

The full environment freeze is stored only in the external backup as
`environment_pip_freeze.txt`.

## Scientific source-of-truth policy

This snapshot records, but does not redefine, the existing priority:

1. `docs/theory/equations.tex`;
2. `docs/theory/assumptions.md`;
3. `docs/literature/source_index.md`;
4. `src/my_project/`;
5. diagnostic scripts and generated results.

The detailed policy is in `AGENTS.md` and `docs/project_rules.md`. No mismatch
was patched during this audit.

## Frozen mathematical and scientific files

All requested paths exist in this checkout.

| Path | SHA256 |
| --- | --- |
| `docs/theory/equations.tex` | `ca9a0d66056784b2da74378b05a84ab84e7cb23d67744ac3a9d0ad6b452abe7a` |
| `docs/theory/assumptions.md` | `ab866d874c4a7c15cef2b0899a34358b9d5e2bdb4f129523161064e4c55da3c8` |
| `docs/project_rules.md` | `99b954c7d1a02570c819b161a61246114fb7214fc5b5deec4f16c9db63f9aac4` |
| `src/my_project/analytic/formulas.py` | `299aa9f0648d1963a1f860cb9750340cb248d670358ad65f5cbc628efcff1215` |
| `src/my_project/analytic/solvers.py` | `70bf81f438e5fb5f25a443f4fab7a7f3e995f01e6c08e8e260c2e6a4129110dd` |
| `src/my_project/analytic/formulas_thickness_mismatch.py` | `ec7febd4875167f999bbc6e78e55c6a75024bf87a8679e29e4ab4857e14ba8dd` |
| `scripts/lib/analytic_branch_tracking.py` | `9db73707b8363b7c777bb8636b06b700791534733420cff6ea2da97fb723df6f` |
| `scripts/lib/analytic_coupled_rods_shapes.py` | `a5a832524646da281c3d8e78030bf96dc7b9a2c741b4e2f5960a193d30ce4549` |
| `scripts/lib/variable_length_timoshenko.py` | `ea1aa0ddd88bbd7b5493db24ab18cc9d5f6de7fe255a3c2e0614e52b193a356e` |
| `scripts/lib/branch_informed_spectrum_continuation.py` | `9ebc4bb5cd2e347205a486764fae6875c13630cbfac9b3a06a7f212b36a5a8d5` |
| `scripts/lib/general_spectrum_completeness.py` | `c9f3d18063b6cc4892fd45353b4d1f137d23017b309ee7effd9be5fe1c6d8679` |

The same data are in `manifests/key_scientific_files_sha256.csv`; all 300
tracked files are in `manifests/tracked_files_sha256.csv`.

## Scientific anchors read from existing canonical outputs

No root, spectrum, mode reconstruction, or FEM calculation was performed to
produce this table.

| Quantity | Value | Canonical source file | Source field or section |
| --- | --- | --- | --- |
| Baseline determinant unknown ordering | `Z=(A1,B1,A2,B2,P1,P2)` | `docs/theory/equations.tex` | equal-length determinant, column-order statement near line 138 |
| Baseline length convention | `l=(l1+l2)/2`, `l1=l(1-mu)`, `l2=l(1+mu)` | `docs/theory/assumptions.md` | base notation and local `mu` definition |
| Base thickness convention | `epsilon=r0/(2l)` | `docs/thickness_mismatch/README.md` | Parameters Held Fixed |
| `eta=0` limit | `tau1=tau2=1`; eta determinant reduces directly to equal-radius baseline | `docs/thickness_mismatch/README.md` | Determinant Convention |
| Mass-preserving factors | `tau1=(1-eta)/sqrt(1+2 mu eta+eta^2)`, `tau2=(1+eta)/sqrt(1+2 mu eta+eta^2)`; `(1-mu)tau1^2+(1+mu)tau2^2=2` | `docs/thickness_mismatch/README.md` | Radius-Ratio Parameter |
| Spectrum target | first `K=10` sorted frequencies with a 10% dimensional-frequency error criterion | `docs/thickness_mismatch/eb_safe_prefix_stage_closure.md` | sections 2–3 |
| Root 11 role | mandatory right-hand completeness guard for roots 1–10; root 12 is optional diagnostic margin | `docs/thickness_mismatch/eb_safe_prefix_stage_closure.md` | section 3 |
| Step-3A status | `counterexample_found`; 28/28 geometries passed K10/root-11; 2 confirmed counterexamples | `results/eb_epsilon_lower_envelope_step3a/eb_epsilon_lower_envelope_step3a_report.md` | Spectrum Quality, Counterexample Assessment, Decision |
| `S3_12` | `epsilon_0=0.029408510742187498`, `beta=90 deg`, `mu=0.7`, `eta=0`, `N_true=4`, `delta_f,5=0.11739469909177262`, `V_5=0.0173946990918` | `docs/thickness_mismatch/eb_safe_prefix_stage_closure.md` | section 4.3 |
| `S3_14` | `epsilon_0=0.024798906738281248`, `beta=45 deg`, `mu=0.5`, `eta=-0.1`, `N_true=4`, `delta_f,5=0.10050934854803291`, `V_6=0.000509348548` | `docs/thickness_mismatch/eb_safe_prefix_stage_closure.md` | section 4.3 and Step-3A report |
| Rule A | `T_A=0.20310844707256814`; 147/155 development objective; `rule_A_no_false_safe_on_checked_sets_benchmark_only` | `results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md` | Exact calibration, Decisions |
| Rule B | `T_s=0.16762413001084248`, `T_r=0.046719283392029604`; 153/155; Phase-I status `rule_B_safety_survives_cost_test_required` | `results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md` | Exact calibration, Decisions |
| Rule S | `T_s=0.16762413001084248`; 153/155; finite checked safety survived, but production-selector path closed | `docs/thickness_mismatch/eb_safe_spectrum_prefix_research_plan.md` | Exact experiment and Stage closure |
| Rule-S / Rule-B equivalence | identical `N_hat` for all 49 checked geometries; rotary threshold nonbinding | `results/eb_rule_ab_exact_pareto/eb_rule_ab_exact_pareto_report.md` | Rule-S extraction |
| Rule-S cost benchmark | `rule_S_cost_not_beneficial`; 4/4 nonempty suffixes required full K10 fallback; median 1.319 s direct vs 1.411 s hybrid | `results/eb_rule_s_cost_break_even/rule_S_cost_break_even_report.md` | Result and operation accounting |
| Out-of-plane EB + torsion | separate diagnostic subsystem; determinant sanity tests and 1D FEM validation workflow exist; no full 3D-FEM validation claim | `docs/theory/out_of_plane_eb_torsion.md` | sections 11–12 |

The out-of-plane generated directories named in documentation are not present
in this checkout. That is recorded as a results/provenance gap, not interpreted
as a mathematical failure.

## Backup anchors

| Artifact | Value |
| --- | --- |
| External directory | `<external-local-backup>` |
| Git bundle | `CoupledBeams_all_refs.bundle` |
| Bundle SHA256 | `fd6ed898056066d8806ad80ce8b6596bec50611f4561051dc0f6638fed4ef51c` |
| Bundle verification | exit 0; complete history; 8 refs; SHA-1 object format |
| Initial full-manifest SHA256 | `d2f341b313f54304f6586479587afce8766bd125881ff8cdf130cd5bea2cf043` |
| Ignored snapshot files | 9,368 (`results`: 9,316; `private`: 52) |
| Ignored snapshot bytes | 1,369,126,508 |
| Ignored snapshot hash mismatches | 0 |
| Worktree overlay files | 0 |

## Lightweight verification

Bytecode generation was disabled with `PYTHONDONTWRITEBYTECODE=1`. Matplotlib
configuration/cache, where imports required it, was redirected to a local OS
temporary directory outside the repository.

The syntax/import-structure audit parsed every tracked Python file under
`src/`, `scripts/`, and `tests/` using the standard-library AST:

```text
command class: ast.parse(Path(...).read_text()) for git-tracked *.py
files parsed: 181
parse errors: 0
runtime: 2.2 s
```

The selected pytest command used `-q -p no:cacheprovider` and covered only
formula/matrix identities, pure metrics, geometry transforms, synthetic helper
logic, documented module imports, filename formatting, and out-of-plane
determinant identities. Results:

```text
collected: 40
passed: 40
failed: 0
skipped: 0
runtime reported by pytest: 6.34 s
```

The full suite was intentionally not run. Static review found tests that call
real EB/Timoshenko root searches, branch continuation, mode reconstruction,
and FEM eigenvalue solves (including analytic smoke, completeness/gateway,
out-of-plane FEM, branch regression, and transform-residual tests). Running
those would violate this task's zero-calculation constraint. No snapshot or
golden file was updated.

```text
new EB root calculations = 0
new Timoshenko root calculations = 0
new FEM calculations = 0
```
