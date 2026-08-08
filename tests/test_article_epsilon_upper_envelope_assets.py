from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

import pytest
from PIL import Image

from scripts.lib import article_epsilon_upper_envelope_assets as assets


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "coarse_grid_v1"
    / "compact_finalization_epsilon_005_resolved"
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    output = tmp_path_factory.mktemp("article_epsilon_assets") / "article_ready_v1"
    source_paths = assets._source_paths(SOURCE_DIR)
    source_before = {name: _sha(path) for name, path in source_paths.items()}
    first = assets.generate_article_assets(
        SOURCE_DIR,
        output,
        language="ru",
        dpi=600,
        overwrite=True,
        generator_arguments=("--test",),
    )
    deterministic_paths = assets._deterministic_asset_paths(output)
    first_hashes = {path.relative_to(output).as_posix(): _sha(path) for path in deterministic_paths}
    second = assets.generate_article_assets(
        SOURCE_DIR,
        output,
        language="ru",
        dpi=600,
        overwrite=True,
        generator_arguments=("--test",),
    )
    second_hashes = {
        path.relative_to(output).as_posix(): _sha(path)
        for path in assets._deterministic_asset_paths(output)
    }
    source_after = {name: _sha(path) for name, path in source_paths.items()}
    return {
        "output": output,
        "first": first,
        "second": second,
        "first_hashes": first_hashes,
        "second_hashes": second_hashes,
        "source_before": source_before,
        "source_after": source_after,
    }


def test_scientific_scope_is_strictly_isotropic_circular() -> None:
    assert assets.SCIENTIFIC_SCOPE == "isotropic_circular_coupled_rods_eb_timoshenko"


def test_anisotropic_modules_are_not_imported() -> None:
    assert not any("anisotropic_rods" in name for name in sys.modules)
    source = Path(assets.__file__).read_text(encoding="utf-8")
    assert "scripts.analysis.anisotropic_rods" not in source


def test_main_epsilon_set_contains_exactly_eight_levels(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    assert [float(row["epsilon_0"]) for row in rows] == list(assets.BASE_EPSILON_LEVELS)


def test_s3_epsilons_are_excluded_from_main_envelope(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    epsilons = {float(row["epsilon_0"]) for row in rows}
    assert all(epsilon not in epsilons for _, epsilon, _, _ in assets.REGRESSION_CONTROLS)


def test_base_case_count_is_1552(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    assert sum(int(row["geometry_case_count"]) for row in rows) == 1552


def test_regression_case_count_is_two(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_regression_controls.csv")
    assert [row["regression_id"] for row in rows] == ["S3_14", "S3_12"]


def test_every_base_epsilon_has_194_cases(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    assert {int(row["geometry_case_count"]) for row in rows} == {194}


def test_raw_envelope_has_expected_values(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    assert [int(row["N_up_h"]) for row in rows] == [10, 10, 10, 10, 9, 6, 4, 3]


def test_suffix_max_equals_raw_envelope(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    assert [row["N_up_h"] for row in rows] == [row["suffix_max_N_up_h"] for row in rows]


def test_all_main_statuses_are_article_valid(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    allowed = {
        "exact_on_adopted_finite_grid",
        "exact_by_K10_saturation_on_adopted_finite_grid",
    }
    assert {row["envelope_status"] for row in rows} <= allowed


def test_0010_and_0020_have_saturation_rationale(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    selected = {row["epsilon_0"]: row for row in rows}
    for epsilon in ("0.010", "0.020"):
        assert selected[epsilon]["envelope_status"].startswith("exact_by_K10_saturation")
        assert "K=10" in selected[epsilon]["exactness_reason"]


def test_0050_has_exact_maximum_four(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_main_envelope.csv")
    row = next(row for row in rows if row["epsilon_0"] == "0.050")
    assert row["N_up_h"] == "4"
    assert row["envelope_status"] == "exact_on_adopted_finite_grid"


def test_0050_has_23_maximizers(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_maximizing_geometries.csv")
    assert sum(row["epsilon_0"] == "0.050" for row in rows) == 23


def test_maximizing_geometries_equal_level_maximum(generated: dict[str, object]) -> None:
    output = generated["output"]
    maxima = {
        row["epsilon_0"]: int(row["N_up_h"])
        for row in _csv(output / "data" / "article_main_envelope.csv")
    }
    for row in _csv(output / "data" / "article_maximizing_geometries.csv"):
        assert int(row["N_true"]) == maxima[row["epsilon_0"]]


def test_unresolved_cases_are_not_distributed(generated: dict[str, object]) -> None:
    output = generated["output"]
    distribution = _csv(output / "data" / "article_N_true_distribution.csv")
    main = {row["epsilon_0"]: row for row in _csv(output / "data" / "article_main_envelope.csv")}
    for epsilon, level in main.items():
        rows = [row for row in distribution if row["epsilon_0"] == epsilon]
        assert sum(int(row["resolved_geometry_count"]) for row in rows) == int(
            level["resolved_case_count"]
        )


def test_distribution_counts_reconcile_to_194(generated: dict[str, object]) -> None:
    output = generated["output"]
    distribution = _csv(output / "data" / "article_N_true_distribution.csv")
    for epsilon in {_row["epsilon_0"] for _row in distribution}:
        rows = [row for row in distribution if row["epsilon_0"] == epsilon]
        unresolved = {int(row["unresolved_geometry_count_at_epsilon"]) for row in rows}
        assert len(unresolved) == 1
        assert sum(int(row["resolved_geometry_count"]) for row in rows) + unresolved.pop() == 194


def test_s3_controls_are_not_in_main_envelope(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_regression_controls.csv")
    assert all(row["included_in_main_envelope"] == "false" for row in rows)


def test_transition_brackets_are_correct(generated: dict[str, object]) -> None:
    rows = _csv(generated["output"] / "data" / "article_transition_brackets.csv")
    actual = [
        (int(row["N_left"]), int(row["N_right"]), row["epsilon_left"], row["epsilon_right"])
        for row in rows
    ]
    assert actual == [
        (10, 9, "0.025", "0.030"),
        (9, 6, "0.030", "0.040"),
        (6, 4, "0.040", "0.050"),
        (4, 3, "0.050", "0.060"),
    ]
    assert {row["status"] for row in rows} == {"finite_grid_transition_bracket"}


def test_main_figure_source_contains_only_eight_points(generated: dict[str, object]) -> None:
    rows = _csv(
        generated["output"] / "figures" / "article_finite_grid_upper_envelope_source.csv"
    )
    assert len(rows) == 8
    assert [int(row["N_up_h"]) for row in rows] == list(assets.EXPECTED_ENVELOPE)


def test_figures_have_requested_dpi(generated: dict[str, object]) -> None:
    output = generated["output"]
    for name in (
        "article_finite_grid_upper_envelope.png",
        "article_N_true_distribution_heatmap.png",
    ):
        with Image.open(output / "figures" / name) as image:
            dpi = image.info["dpi"]
        assert dpi[0] == pytest.approx(600, abs=0.1)
        assert dpi[1] == pytest.approx(600, abs=0.1)


def test_article_text_uses_finite_grid_qualification(generated: dict[str, object]) -> None:
    output = generated["output"]
    text = (output / "text" / "article_results_section_ru.md").read_text(encoding="utf-8")
    required = (
        "эмпирическая верхняя огибающая на принятой конечной сетке",
        "максимальное наблюдаемое число допустимых частот",
        "сопоставлялись одинаковые позиции в упорядоченных спектрах",
        "результат не является доказанной границей на непрерывном множестве геометрических параметров",
        "соединяющие линии на рисунке служат только для визуального восприятия",
    )
    assert all(phrase in text.lower() for phrase in required)


def test_article_text_avoids_overstated_claims_and_short_length(
    generated: dict[str, object],
) -> None:
    output = generated["output"]
    texts = [
        (output / "text" / name).read_text(encoding="utf-8").lower()
        for name in ("article_results_section_ru.md", "article_results_short_ru.md")
    ]
    forbidden = (
        "для любой геометрии",
        "строго доказано для всех параметров",
        "априорная оценка",
        "одинаковые физические моды",
        "получена непрерывная огибающая",
    )
    assert not any(phrase in text for phrase in forbidden for text in texts)
    russian_words = re.findall(r"[А-Яа-яЁё]+(?:[-—][А-Яа-яЁё]+)*", texts[1])
    assert 180 <= len(russian_words) <= 250


def test_scientific_operation_counters_are_zero(generated: dict[str, object]) -> None:
    manifest = json.loads(
        (generated["output"] / "metadata" / "article_asset_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    counter_keys = [key for key in manifest if key.endswith("_calls")]
    assert counter_keys
    assert all(manifest[key] == 0 for key in counter_keys)


def test_generated_hashes_are_deterministic(generated: dict[str, object]) -> None:
    assert generated["first_hashes"] == generated["second_hashes"]
    assert generated["second"]["deterministic_repeat_match"] is True
    assert generated["second"]["gates"]["article_assets_readiness_gate"] is True


def test_source_files_are_unchanged(generated: dict[str, object]) -> None:
    assert generated["source_before"] == generated["source_after"]


def test_raw_point_caches_are_not_read_or_deleted(generated: dict[str, object]) -> None:
    manifest = json.loads(
        (generated["output"] / "metadata" / "article_asset_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["raw_cache_files_read"] == 0
    assert manifest["raw_cache_files_deleted"] == 0


def test_maximizer_summary_matches_detail(generated: dict[str, object]) -> None:
    output = generated["output"]
    detail = _csv(output / "data" / "article_maximizing_geometries.csv")
    summary = _csv(output / "data" / "article_maximizer_summary.csv")
    counts = Counter(row["epsilon_0"] for row in detail)
    thin_counts = Counter(
        row["epsilon_0"] for row in detail if row["thinness_flag"] == "true"
    )
    for row in summary:
        assert int(row["maximizing_geometry_count"]) == counts[row["epsilon_0"]]
        assert int(row["thinness_flagged_count"]) == thin_counts[row["epsilon_0"]]
