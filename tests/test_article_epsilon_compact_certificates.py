from __future__ import annotations

import gc
import gzip
import inspect
import json
import math
import sys
import tracemalloc
from pathlib import Path

import pytest

from scripts.analysis.thickness_mismatch.audits import (
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.lib import article_epsilon_compact_certificates as compact


def _manifest(case_id: str = "AUE_00000000000000000000") -> dict[str, str]:
    return {
        "case_id": case_id, "epsilon_0": "0.01", "beta_deg": "30",
        "mu": "0.5", "eta": "-0.1", "s_max": "0.05",
    }


def _payload(*, n_true: int | None = 2, status: str = "resolved_prefix_early_stop") -> dict[str, object]:
    eb = [1.0, 2.0, 3.2, 4.1, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    timo = [1.0, 2.0, 3.0, 4.0, 4.9, 5.9, 6.9, 7.9, 8.9, 9.9, 10.9]
    deltas = compact.calculate_deltas(eb, timo)
    return {
        "execution_status": status, "N_true": n_true,
        "first_failed_mode": 3 if n_true == 2 else None,
        "prefix_guard_resolved_through": 3 if n_true == 2 else 10,
        "deltas": deltas,
        "models": {
            compact.MODEL_EB: {
                "resolved_roots": eb, "clusters": [],
                "brackets_and_continuation_metadata": [],
                "verification_evaluator": {"dense": [1.0] * 1000},
            },
            compact.MODEL_TIMO: {
                "resolved_roots": timo, "clusters": [],
                "brackets_and_continuation_metadata": [],
                "verification_evaluator": {"dense": [1.0] * 1000},
            },
        },
        "identity": {
            "algorithm_version": "test", "prefix_strategy": "paired",
            "strict_policy": "main-pass",
            "scientific_model_configuration": {
                "eb_matrix_evaluator_version": "eb-test",
                "timoshenko_basis_evaluator_version": "timo-test",
                "primary_settings": {"root_match_tol": 2e-4},
                "strict_settings": {},
            },
        },
        "force_strict_requested": 0, "force_strict_executed": 0,
    }


def _certificate(payload: dict[str, object] | None = None) -> dict[str, object]:
    return compact.extract_certificate(
        payload or _payload(), _manifest(), manifest_sha="manifest",
        source_commit="commit", source_path="cache/source.json.gz",
        source_sha="source", source_size=123,
    )


def test_compact_cli_is_explicit_and_default_is_unchanged() -> None:
    assert not runner.parse_args([]).compact_only
    args = runner.parse_args([
        "--compact-only", "--build-compact-point-certificates",
        "--no-new-point-solves",
    ])
    assert args.compact_only and args.no_new_point_solves
    with pytest.raises(SystemExit):
        runner.parse_args(["--compact-only", "--build-compact-point-certificates"])
    post = runner.parse_args([
        "--family-post-stage-only", "--use-compact-point-certificates",
        "--no-new-point-solves", "--defer-expensive-strict",
    ])
    assert post.family_post_stage_only


def test_production_family_stage_uses_compact_layer() -> None:
    source = inspect.getsource(runner.main)
    local_branch = source[source.index('if args.family_inventory_policy == "local-repair":'):]
    assert "build_compact_certificates" in local_branch
    assert "run_compact_family_poststage" in local_branch
    assert "family_integration.run_shadow" not in local_branch


def test_full_payload_converts_without_large_trace_fields() -> None:
    cert = _certificate()
    assert cert["scientific_scope"] == compact.SCIENTIFIC_SCOPE
    encoded = json.dumps(cert)
    assert "verification_evaluator" not in encoded
    assert "dense" not in encoded
    assert cert["schema_version"] == compact.SCHEMA_VERSION


def test_compact_certificate_reproduces_n_failure_guard_and_deltas() -> None:
    cert = _certificate()
    result = cert["result"]
    spectra = cert["spectra"]
    n_true, failure, guard, deltas = compact.reproduce_prefix(
        spectra[compact.MODEL_EB]["roots"], spectra[compact.MODEL_TIMO]["roots"]
    )
    assert (n_true, failure, guard) == (2, 3, 4)
    assert result["N_true"] == n_true
    assert result["first_failed_mode"] == failure
    assert result["required_guard"] == guard
    assert deltas == spectra["delta_f"]


def test_n_true_ten_requires_root_eleven() -> None:
    payload = _payload(n_true=10, status="resolved_full_K10")
    payload["models"][compact.MODEL_EB]["resolved_roots"] = [float(k) for k in range(1, 12)]
    payload["models"][compact.MODEL_TIMO]["resolved_roots"] = [float(k) for k in range(1, 12)]
    payload["deltas"] = {str(k): 0.0 for k in range(1, 11)}
    cert = _certificate(payload)
    assert cert["result"]["N_true"] == 10
    assert cert["result"]["required_guard"] == 11
    assert cert["result"]["required_guard_confirmed"] is True


def test_unresolved_and_scientific_zero_are_distinct() -> None:
    unresolved = _payload(n_true=None, status="attempted_unresolved")
    unresolved["first_failed_mode"] = None
    unresolved_cert = _certificate(unresolved)
    assert unresolved_cert["result"]["N_true"] is None
    zero = _payload(n_true=0)
    zero["models"][compact.MODEL_EB]["resolved_roots"][0] = 1.2
    zero["deltas"] = compact.calculate_deltas(
        zero["models"][compact.MODEL_EB]["resolved_roots"],
        zero["models"][compact.MODEL_TIMO]["resolved_roots"],
    )
    zero["first_failed_mode"] = 1
    zero_cert = _certificate(zero)
    assert zero_cert["result"]["N_true"] == 0
    assert zero_cert["result"]["required_guard"] == 2


def test_corrupt_scientific_result_is_rejected() -> None:
    payload = _payload()
    payload["N_true"] = 7
    with pytest.raises(ValueError, match="not reproduced"):
        _certificate(payload)


def test_atomic_compact_write_is_deterministic(tmp_path: Path) -> None:
    path = tmp_path / "case.json.gz"
    cert = _certificate()
    compact.atomic_write_gzip_json(path, cert)
    first = path.read_bytes()
    compact.atomic_write_gzip_json(path, cert)
    assert path.read_bytes() == first
    assert compact.load_certificate(path) == cert
    assert not list(tmp_path.glob("*.tmp"))


def test_stream_iterator_holds_only_current_payload(tmp_path: Path) -> None:
    rows = []
    paths = {}
    for index in range(80):
        case_id = f"AUE_{index:020x}"
        row = _manifest(case_id)
        rows.append(row)
        path = tmp_path / f"{case_id}.json.gz"
        with gzip.open(path, "wt", encoding="utf-8") as stream:
            json.dump({"blob": "x" * 200_000, "case_id": case_id}, stream)
        paths[case_id] = path
    tracemalloc.start()
    baseline = tracemalloc.get_traced_memory()[0]
    seen = 0
    for row, path, payload in compact.iter_full_payloads(rows, paths):
        assert payload["case_id"] == row["case_id"]
        seen += 1
    gc.collect()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert seen == 80
    assert peak - baseline < 8_000_000
    assert current - baseline < 2_000_000


def test_compact_pseudo_payload_contains_only_detector_inputs() -> None:
    pseudo = compact.compact_pseudo_payload(_certificate())
    assert set(pseudo["models"]) == set(compact.MODELS)
    assert "verification_evaluator" not in json.dumps(pseudo)


def test_scope_isolation_and_no_anisotropic_import() -> None:
    assert compact.SCIENTIFIC_SCOPE == "isotropic_circular_coupled_rods_eb_timoshenko"
    assert not any("anisotropic_rods" in name for name in sys.modules)


def test_source_roots_preserve_confirmed_multiplicity_slots() -> None:
    payload = _payload()
    roots = payload["models"][compact.MODEL_EB]["resolved_roots"]
    roots[1] = roots[0]
    payload["models"][compact.MODEL_EB]["clusters"] = [
        {"cluster_id": "c", "cluster_member_index": 1, "cluster_size": 2, "sorted_index": 1},
        {"cluster_id": "c", "cluster_member_index": 2, "cluster_size": 2, "sorted_index": 2},
    ]
    payload["N_true"] = None
    payload["execution_status"] = "attempted_unresolved"
    payload["first_failed_mode"] = None
    payload["deltas"] = {}
    cert = _certificate(payload)
    saved = cert["spectra"][compact.MODEL_EB]
    assert saved["roots"][0] == saved["roots"][1]
    assert [item["multiplicity_slot"] for item in saved["quality_by_rank"][:2]] == [1, 2]
