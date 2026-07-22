from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts.analysis.thickness_mismatch.maps import (
    plot_counterexample_dimensional_frequency_beta as target,
)
from scripts.lib import branch_informed_spectrum_continuation as branch


def fake_branch(index: int, beta_deg: float, *, model: str = target.MODEL_EB) -> branch.ContinuedBranch:
    return branch.ContinuedBranch(
        model=model,
        branch_id=f"B{index:03d}",
        parent_family="bending" if index % 3 else "axial",
        beta_deg=beta_deg,
        predicted_Lambda=float(index),
        Lambda=float(index),
        right_vector=(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        left_vector=(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        sigma_1=1.0e-9,
        sigma_ratio=1.0e-7,
        self_MAC=1.0,
        nullity=1,
        cluster_id="",
        refinement_status=branch.SEED_REFINED_TO_NEW_ROOT,
        detection_source="test",
    )


def fake_cases() -> tuple[target.CaseSpec, ...]:
    return (
        target.CaseSpec("S3_12", 0.029, 0.7, 0.0, 90.0),
        target.CaseSpec("S3_14", 0.024, 0.5, -0.1, 45.0),
    )


def fake_points(*, unresolved: bool = False) -> dict[tuple[str, str], tuple[target.SpectrumPoint, ...]]:
    output: dict[tuple[str, str], tuple[target.SpectrumPoint, ...]] = {}
    for case in fake_cases():
        for model in target.MODELS:
            points = []
            for beta in (0.0, 45.0, 90.0):
                k10 = not (unresolved and beta == 45.0 and model == target.MODEL_EB)
                points.append(
                    target.SpectrumPoint(
                        case_id=case.case_id,
                        model=model,
                        beta_deg=beta,
                        branches=tuple(fake_branch(index, beta, model=model) for index in range(1, 13)),
                        k10_guard_resolved=k10,
                        full12_resolved=k10,
                        guard_status="guard_resolved" if k10 else "guard_unresolved",
                        strict_fallback_used=False,
                        strict_fallback_count=0,
                        cache_status="test",
                        warnings=() if k10 else ("test_unresolved",),
                    )
                )
            output[(case.case_id, model)] = tuple(points)
    return output


def synthetic_plot_rows() -> list[dict[str, object]]:
    return target.build_frequency_rows(fake_cases(), fake_points(), 10)


def test_loads_full_precision_counterexample_epsilons_from_step3a_outputs() -> None:
    cases = target.load_step3a_cases(target.DEFAULT_STEP3A_DIR)
    values = {case.case_id: case.epsilon_0 for case in cases}
    assert values == {
        "S3_12": float("2.9408510742187498e-02"),
        "S3_14": float("2.4798906738281248e-02"),
    }


def test_entrypoint_does_not_embed_rounded_prompt_epsilons() -> None:
    source = target.SCRIPT_PATH.read_text(encoding="utf-8")
    assert "0.0294085107" not in source
    assert "0.0247989067" not in source


def test_canonical_dimensional_scaling_is_common_and_quadratic() -> None:
    epsilon = float("2.9408510742187498e-02")
    params = target.beam_params_from_epsilon(epsilon)
    assert params.E == target.CANONICAL_E_PA
    assert params.rho == target.CANONICAL_RHO_KG_M3
    assert params.L_total == target.CANONICAL_L_TOTAL_M
    assert params.eps == pytest.approx(epsilon, abs=2.0e-16)
    frequencies = target.dimensional_frequencies((1.0, 2.0, 3.0), epsilon)
    scale = target.frequency_scale(params)
    np.testing.assert_allclose(frequencies, scale * np.asarray((1.0, 4.0, 9.0)))
    scale_rows = target.build_scale_rows((target.load_step3a_cases(target.DEFAULT_STEP3A_DIR)[0],))
    row = scale_rows[0]
    assert row["scaling_formula"] == target.FREQUENCY_FORMULA
    assert row["source_helper_file"] == target.SCALING_SOURCE
    assert row["EB_scale_factor"] == row["Timoshenko_scale_factor"]
    assert row["EB_Timoshenko_scale_equality_check"] is True


def test_beta_grid_contains_exact_counterexample_angles() -> None:
    grid = target.beta_grid(0.0, 90.0, 0.5)
    assert len(grid) == 181
    assert 45.0 in grid
    assert 90.0 in grid


def test_frequency_rows_include_exactly_target_ten_and_not_guard_roots() -> None:
    rows = target.build_frequency_rows(fake_cases(), fake_points(), 10)
    assert {int(row["sorted_index"]) for row in rows} == set(range(1, 11))
    assert len(rows) == 2 * 2 * 3 * 10
    assert all(int(row["sorted_index"]) not in {11, 12} for row in rows)
    quality = target.build_quality_rows(fake_cases(), fake_points())
    assert all(float(row["root11"]) == 11.0 for row in quality)
    assert all(row["full12_resolved"] is True for row in quality)


def test_unresolved_k10_is_nan_and_never_substituted() -> None:
    rows = target.build_frequency_rows(fake_cases(), fake_points(unresolved=True), 10)
    unresolved = [
        row
        for row in rows
        if row["case_id"] == "S3_12"
        and row["model"] == target.MODEL_EB
        and row["beta_deg"] == 45.0
    ]
    assert len(unresolved) == 10
    assert all(math.isnan(float(row["Lambda"])) for row in unresolved)
    assert all(math.isnan(float(row["dimensional_frequency"])) for row in unresolved)
    assert all(row["branch_id"] == "" for row in unresolved)


def test_figure_has_exact_line_and_text_contract() -> None:
    fig, ax = target.create_case_figure(synthetic_plot_rows(), "S3_12")
    try:
        assert len(ax.lines) == 20
        for index in range(10):
            eb = ax.lines[2 * index]
            timo = ax.lines[2 * index + 1]
            assert eb.get_linestyle() == "--"
            assert timo.get_linestyle() == "-"
            assert eb.get_color() == timo.get_color()
        assert ax.get_legend() is None
        assert ax.get_title() == ""
        assert fig._suptitle is None
        assert len(ax.texts) == 0
        assert len(fig.texts) == 0
        assert ax.get_xlabel() == r"$\beta,\ ^\circ$"
        assert ax.get_ylabel() == r"$f,\ \mathrm{Hz}$"
        assert not any(line.get_marker() not in {"None", None, ""} for line in ax.lines)
        assert not any(gridline.get_visible() for gridline in ax.get_xgridlines())
        assert not any(gridline.get_visible() for gridline in ax.get_ygridlines())
    finally:
        plt.close(fig)


def test_exactly_two_graphical_outputs_and_csv_round_trip(tmp_path: Path) -> None:
    rows = synthetic_plot_rows()
    target.write_csv(tmp_path / target.MAIN_CSV, rows, target.MAIN_FIELDS)
    target.create_plots(rows, tmp_path)
    graphics = sorted(path.name for path in tmp_path.iterdir() if path.suffix.lower() in {".pdf", ".png", ".svg"})
    assert graphics == sorted(target.PDF_NAMES.values())

    def forbidden_solver(*_args, **_kwargs):
        raise AssertionError("plot-only attempted a root calculation")

    result = target.main(
        ["--output-dir", str(tmp_path), "--plot-only"],
        sweep_solver=forbidden_solver,
    )
    assert result["root_calculations"] == 0
    assert all((tmp_path / name).exists() for name in target.PDF_NAMES.values())


def test_plot_only_requires_no_step3a_or_gateway_inputs(tmp_path: Path) -> None:
    rows = synthetic_plot_rows()
    target.write_csv(tmp_path / target.MAIN_CSV, rows, target.MAIN_FIELDS)
    missing = tmp_path / "does_not_exist"
    args = target.parse_args(
        [
            "--step3a-dir",
            str(missing),
            "--gateway-dir",
            str(missing),
            "--output-dir",
            str(tmp_path),
            "--plot-only",
        ]
    )
    result = target.execute(args)
    assert result["root_calculations"] == 0


def test_parser_disables_abbreviated_options() -> None:
    with pytest.raises(SystemExit):
        target.parse_args(["--beta-s", "1"])
