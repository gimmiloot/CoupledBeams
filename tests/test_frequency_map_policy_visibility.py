from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "docs/numerics/frequency_map_computation_policy.md"
NUMERICS_INDEX = ROOT / "docs/numerics/README.md"
REDIRECT = ROOT / "docs/thickness_mismatch/frequency_plot_computation_policy.md"
CANONICAL_REPO_PATH = "docs/numerics/frequency_map_computation_policy.md"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _assert_markers(text: str, markers: tuple[str, ...]) -> None:
    missing = [marker for marker in markers if marker not in text]
    assert not missing, f"Missing semantic markers: {missing}"


def test_canonical_policy_and_numerics_index_exist() -> None:
    assert CANONICAL.is_file()
    assert NUMERICS_INDEX.is_file()
    _assert_markers(
        _read(NUMERICS_INDEX),
        ("frequency_map_computation_policy.md", "scripts/sweep_grid_policy.py"),
    )


def test_canonical_policy_contains_required_contract_markers() -> None:
    text = _read(CANONICAL)
    _assert_markers(
        text,
        (
            "frequency-map-v1",
            "fast_plot",
            "certified_audit",
            "plot_only",
            "sorted_positions",
            "tracked_branches",
            "frequency_map_policy",
            "calculation_mode",
            "spectrum_semantics",
            "sweep_parameter",
            "parameter_grid",
            "K_plot",
            "K_guard",
            "guard_root_role",
            "neighbour_audit",
            "local_repair_policy",
            "strict_audit_default",
            "policy_override_reason",
        ),
    )


def test_old_policy_path_is_a_short_compatibility_redirect() -> None:
    text = _read(REDIRECT)
    _assert_markers(text, (CANONICAL_REPO_PATH, "older links"))
    assert len(text) < 1200
    assert len([line for line in text.splitlines() if line.strip()]) <= 8
    for duplicated_section in (
        "## `fast_plot`",
        "## `certified_audit`",
        "## `plot_only`",
        "Required local policy instance",
        "frequency_map_policy: frequency-map-v1",
    ):
        assert duplicated_section not in text


def test_root_governance_files_expose_the_policy() -> None:
    agents = _read(ROOT / "AGENTS.md")
    _assert_markers(agents, (CANONICAL_REPO_PATH, "fast_plot", "plot_only"))

    project_rules = _read(ROOT / "docs/project_rules.md")
    _assert_markers(
        project_rules,
        (
            CANONICAL_REPO_PATH,
            "frequency-map-v1",
            "fast_plot",
            "policy_override_reason",
            "inherits",
        ),
    )


def test_project_navigation_exposes_stable_policy() -> None:
    root_readme = _read(ROOT / "README.md")
    _assert_markers(root_readme, ("docs/numerics/", CANONICAL_REPO_PATH))

    research_index = _read(ROOT / "docs/research_index.md")
    row = next(
        line
        for line in research_index.splitlines()
        if "Frequency-map computation policy" in line
    )
    _assert_markers(
        row,
        (
            "stable-baseline",
            "numerics/frequency_map_computation_policy.md",
            "fast_plot",
            "certified_audit",
        ),
    )

    scripts_readme = _read(ROOT / "scripts/README.md")
    _assert_markers(
        scripts_readme,
        (
            CANONICAL_REPO_PATH,
            "frequency-map-v1",
            "plot_only",
            "Script Proliferation Control",
        ),
    )


def test_active_model_directions_inherit_project_policy() -> None:
    direction_readmes = (
        ROOT / "docs/thickness_mismatch/README.md",
        ROOT / "docs/anisotropic_rods/README.md",
        ROOT / "docs/laminated_beams/README.md",
    )
    for path in direction_readmes:
        _assert_markers(
            _read(path),
            ("numerics/frequency_map_computation_policy.md", "frequency-map-v1"),
        )


def test_rlb_2d_note_declares_future_local_instance_without_status_promotion() -> None:
    note = _read(
        ROOT
        / "docs/laminated_beams"
        / "rectangular_weakly_orthotropic_mu_beta_graphs_note.md"
    )
    _assert_markers(
        note,
        (
            "frequency_map_policy: frequency-map-v1",
            "future_default_calculation_mode: fast_plot",
            "spectrum_semantics: sorted_positions",
            "K_plot: 8",
            "K_guard: 9",
            "guard_root_role: completeness_only",
            "branch_tracking: false",
            "mac: false",
            "PARTIAL_PASS_PLOTTED_RANGE_QUALIFICATIONS",
        ),
    )


def test_active_navigation_does_not_use_superseded_policy_path() -> None:
    active_navigation = (
        ROOT / "README.md",
        ROOT / "docs/project_rules.md",
        ROOT / "docs/research_index.md",
        ROOT / "docs/results_index.md",
        ROOT / "scripts/README.md",
        ROOT / "scripts/STATUS.md",
        ROOT / "scripts/analysis/thickness_mismatch/README.md",
        ROOT / "docs/thickness_mismatch/README.md",
        ROOT / "docs/anisotropic_rods/README.md",
        ROOT / "docs/laminated_beams/README.md",
    )
    old_path = "docs/thickness_mismatch/frequency_plot_computation_policy.md"
    for path in active_navigation:
        assert old_path not in _read(path), path
