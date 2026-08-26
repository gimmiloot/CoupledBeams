"""Postprocess the frozen RLB-1C-ISO determinant inventories.

This script-local helper is deliberately downstream of the two independent
root workers.  It never writes either inventory, never supplies roots to a
determinant search, and contains no Ritz, Euler--Bernoulli, or FEM model.  Its
only public entry point is :func:`run_postprocessing`.

The comparison layer may import both frozen RLB modules and the independent
rectangular Timoshenko comparator because the caller must already have frozen
and hashed both inventories.  The helper verifies that contract before doing
any scientific comparison and verifies the inventory bytes again on return.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize_scalar

FloatArray = NDArray[np.float64]
MatrixProvider = Callable[[float], FloatArray]

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

# These modules are intentionally loaded only inside run_postprocessing(),
# after both inventory manifests and all six immutable CSV hashes pass.  Thus
# importing this helper from the main CLI cannot violate worker isolation.
legacy: Any = None
geometry: Any = None
coupled: Any = None
rlb: Any = None

INVENTORY_PREFIXES = ("rlb", "legacy_rectangular")
INVENTORY_SUFFIXES = (
    "root_candidates.csv",
    "rejected_candidates.csv",
    "roots.csv",
)
REQUIRED_CASE_IDS = tuple(f"ISO-{index:02d}" for index in range(1, 9))
SELECTED_POSITIONS = {"ISO-03": (1, 2, 3, 6), "ISO-04": (1, 2, 3)}
EXPECTED_CASE_CONTRACT = {
    "ISO-01": ("G20", 1.0, 1.0, 0.0),
    "ISO-02": ("G20", 0.7, 1.3, 0.0),
    "ISO-03": ("G20", 1.0, 1.0, 30.0),
    "ISO-04": ("G20", 1.0, 1.0, 90.0),
    "ISO-05": ("G20", 1.0, 1.0, -30.0),
    "ISO-06": ("G20", 0.7, 1.3, 30.0),
    "ISO-07": ("G20", 1.3, 0.7, 30.0),
    "ISO-08": ("G10", 1.0, 1.0, 30.0),
}
S1_OLD_TO_RLB = np.diag([1.0, 1.0, -1.0, 1.0, 1.0, -1.0])
S2_OLD_TO_RLB = np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])

EXPECTED_RLB_INVENTORY_SHA256 = (
    "4E0DDE518A134F5C58A8A36F25B6C40A9750CFBE361288A0212F9664360EA837"
)
EXPECTED_LEGACY_INVENTORY_SHA256 = (
    "6F2A06D39EC66D2E49DE8D519D201BFBF42150C1BA9D6E493EB2EFEB44F534D6"
)
SCIENTIFIC_COMPARISON_KIND = "RLB_VS_LEGACY_RECTANGULAR"
DIRECT_BETA0_AUXILIARY_KINDS = (
    "BETA0_DIRECT_PRIMARY_VERIFICATION",
    "BETA0_BENDING_PRIMARY_VERIFICATION",
    "BETA0_RLB_VS_DIRECT_FIXED_FIXED",
    "BETA0_LEGACY_COUPLED_VS_DIRECT_FIXED_FIXED",
    "BETA0_DIRECT_FIXED_FIXED_VS_AXIAL_BENDING_UNION",
    "BETA0_RLB_VS_AXIAL_BENDING_UNION",
)
DIRECT_BETA0_QUALIFICATION_REASON = (
    "CLOSED_FORM_MIXED_REGIME_BASIS_CONDITIONING"
)
SCIENTIFIC_RESULT_STATEMENT = (
    "FOUR_PLY_ISOTROPIC_LIMIT_VALIDATED_FOR_DECLARED_FINITE_CASES_"
    "AGAINST_INDEPENDENT_RECTANGULAR_TIMOSHENKO_DETERMINANT"
)
SCIENTIFIC_STATUSES = {
    "RLB-ISO-4PLY-CONSTITUTIVE": "PASS",
    "RLB-ISO-SECTION-REDUCTION": "PASS",
    "RLB-ISO-LEGACY-RECTANGULAR-ADAPTER": "PASS",
    "RLB-ISO-LOCAL-ARM-EQUIVALENCE": "PASS",
    "RLB-ISO-COUPLED-SPECTRUM": "PASS",
    "RLB-ISO-MODE-SHAPES": "PASS",
}
DIRECT_BETA0_AUXILIARY_STATUS = (
    "INCONCLUSIVE_NUMERICAL_BASIS_CONDITIONING"
)
ARM_EXCHANGE_PASS_STATUS = "PASS_AFTER_BOUNDED_LOCAL_REFINEMENT"
ARM_EXCHANGE_PARTIAL_STATUS = "PARTIAL_PASS_ROOT_LOCALIZATION"
SCIENTIFIC_OVERALL_STATUS = "PASS_WITH_AUXILIARY_NUMERICAL_QUALIFICATIONS"

DIRECT_NUMERIC_FIELDS = (
    "left_slot_start",
    "left_slot_end",
    "right_slot_start",
    "right_slot_end",
    "left_omega",
    "right_omega",
    "relative_difference",
    "left_neighbor_gap_before",
    "left_neighbor_gap_after",
    "right_neighbor_gap_before",
    "right_neighbor_gap_after",
    "left_root_count",
    "right_root_count",
    "left_multiplicity",
    "right_multiplicity",
    "left_nullity",
    "right_nullity",
    "tolerance",
    "alignment_slot",
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest().upper()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0.0 else "-Infinity"
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    data = [dict(row) for row in rows]
    fields: list[str] = []
    for row in data:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in data:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fields})


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _inventory_paths(output_directory: Path) -> tuple[Path, ...]:
    return tuple(
        output_directory / f"{prefix}_{suffix}"
        for prefix in INVENTORY_PREFIXES
        for suffix in INVENTORY_SUFFIXES
    )


def _validate_finite_case_contract(contract: Mapping[str, Any]) -> None:
    cases = list(contract["cases"])
    case_ids = tuple(str(case["case_id"]) for case in cases)
    if case_ids != REQUIRED_CASE_IDS:
        raise RuntimeError(f"postprocess requires exactly {REQUIRED_CASE_IDS}, got {case_ids}")
    for case in cases:
        case_id = str(case["case_id"])
        expected_geometry, expected_l1, expected_l2, expected_beta = (
            EXPECTED_CASE_CONTRACT[case_id]
        )
        actual = (
            str(case["geometry"]),
            float(case["L1"]),
            float(case["L2"]),
            float(case["beta_deg"]),
        )
        expected = (expected_geometry, expected_l1, expected_l2, expected_beta)
        if actual[0] != expected[0] or any(
            not math.isclose(left, right, rel_tol=0.0, abs_tol=8.0 * np.finfo(float).eps)
            for left, right in zip(actual[1:], expected[1:], strict=True)
        ):
            raise RuntimeError(
                f"{case_id} differs from the frozen finite case contract: "
                f"expected {expected}, got {actual}"
            )
    selected = {
        str(case_id): tuple(int(position) for position in positions)
        for case_id, positions in contract["mode_positions"].items()
    }
    if selected != SELECTED_POSITIONS:
        raise RuntimeError(
            f"mode positions differ from the frozen selection: {selected}"
        )


def _verify_inventory_files(
    contract_path: Path, output_directory: Path
) -> tuple[dict[str, str], dict[str, Any], dict[str, Any]]:
    """Verify worker manifests, semantic hashes, and every frozen CSV byte."""

    contract_hash = _sha256_file(contract_path)
    snapshot_path = output_directory / "case_contract.json"
    if not snapshot_path.is_file():
        raise RuntimeError("case_contract.json is missing")
    source_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    snapshot_contract = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if snapshot_contract != source_contract:
        raise RuntimeError("case_contract.json differs semantically from the worker input")
    preliminary_path = output_directory / "preliminary_gate.json"
    if not preliminary_path.is_file():
        raise RuntimeError("preliminary_gate.json is missing")
    preliminary = json.loads(preliminary_path.read_text(encoding="utf-8"))
    if preliminary.get("preliminary_status") != "PASS":
        raise RuntimeError("preliminary gate is not PASS")
    if preliminary.get("contract_sha256") != contract_hash:
        raise RuntimeError("preliminary gate and case contract hashes differ")

    manifests: dict[str, dict[str, Any]] = {}
    snapshots: dict[str, str] = {}
    for prefix in INVENTORY_PREFIXES:
        manifest_path = output_directory / f"{prefix}_inventory_manifest.json"
        if not manifest_path.is_file():
            raise RuntimeError(f"missing frozen inventory manifest: {manifest_path.name}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifests[prefix] = manifest
        if manifest.get("contract_sha256") != contract_hash:
            raise RuntimeError(f"{prefix} inventory contract hash differs")
        if manifest.get("comparison_performed") is not False:
            raise RuntimeError(f"{prefix} inventory was not frozen before comparison")
        if manifest.get("cross_model_roots_read") is not False:
            raise RuntimeError(f"{prefix} worker read cross-model roots")
        output_hashes = manifest.get("output_sha256", {})
        for suffix in INVENTORY_SUFFIXES:
            path = output_directory / f"{prefix}_{suffix}"
            if not path.is_file():
                raise RuntimeError(f"missing frozen inventory file: {path.name}")
            actual = _sha256_file(path)
            expected = output_hashes.get(path.name)
            if actual != expected:
                raise RuntimeError(f"frozen inventory hash mismatch: {path.name}")
            snapshots[path.name] = actual
    return snapshots, manifests["rlb"], manifests["legacy_rectangular"]


def _load_scientific_models_after_inventory_freeze() -> None:
    global legacy, geometry, coupled, rlb
    if all(module is not None for module in (legacy, geometry, coupled, rlb)):
        return
    legacy = importlib.import_module(
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams"
    )
    geometry = importlib.import_module("scripts.lib.reddy_inplane_geometry")
    coupled = importlib.import_module("scripts.lib.reddy_symmetric_coupled_beams")
    rlb = importlib.import_module("scripts.lib.reddy_symmetric_laminated_beam")


@dataclass(frozen=True)
class RootEventRecord:
    case_id: str
    event_id: str
    Omega: float
    omega: float
    multiplicity: int
    nullity: int
    cluster_id: str
    cluster_semantics: str
    cluster_multiplicity: int
    cluster_total_nullity: int
    cluster_center_Omega: float
    slots: tuple[int, ...]
    roles: tuple[str, ...]


@dataclass(frozen=True)
class SpectrumBlock:
    case_id: str
    block_id: str
    events: tuple[RootEventRecord, ...]
    slot_start: int
    slot_end: int
    root_count: int
    multiplicity: int
    nullity: int
    centre_omega: float
    centre_Omega: float
    clustered: bool
    cluster_semantics: str
    roles: tuple[str, ...]
    gap_before_omega: float = math.inf
    gap_after_omega: float = math.inf


def _parse_root_events(rows: Sequence[Mapping[str, str]]) -> dict[str, list[RootEventRecord]]:
    """Parse a slot inventory and reject incomplete cluster metadata."""

    required_fields = (
        "case_id",
        "sorted_slot",
        "role",
        "repeated_root_slot",
        "event_id",
        "omega",
        "Omega",
        "multiplicity",
        "detected_nullity",
        "cluster_id",
        "cluster_semantics",
        "cluster_multiplicity",
        "cluster_total_nullity",
        "cluster_center_Omega",
    )
    event_constant_fields = (
        "case_id",
        "event_id",
        "omega",
        "Omega",
        "multiplicity",
        "detected_nullity",
        "cluster_id",
        "cluster_semantics",
        "cluster_multiplicity",
        "cluster_total_nullity",
        "cluster_center_Omega",
    )
    grouped: dict[tuple[str, str], list[Mapping[str, str]]] = {}
    slots_by_case: dict[str, set[int]] = {case_id: set() for case_id in REQUIRED_CASE_IDS}
    for row in rows:
        missing = [field for field in required_fields if field not in row]
        if missing:
            raise RuntimeError(f"root inventory row is missing fields: {missing}")
        case_id = str(row["case_id"])
        if case_id not in REQUIRED_CASE_IDS:
            raise RuntimeError(f"unexpected root-inventory case: {case_id}")
        event_id = str(row["event_id"])
        if not event_id:
            raise RuntimeError(f"empty event_id in {case_id}")
        slot = int(row["sorted_slot"])
        if slot <= 0 or slot in slots_by_case[case_id]:
            raise RuntimeError(f"invalid or duplicate sorted slot {slot} in {case_id}")
        slots_by_case[case_id].add(slot)
        grouped.setdefault((case_id, event_id), []).append(row)
    by_case: dict[str, list[RootEventRecord]] = {case_id: [] for case_id in REQUIRED_CASE_IDS}
    for (case_id, event_id), event_rows in grouped.items():
        first = event_rows[0]
        for field in event_constant_fields:
            values = {str(row[field]) for row in event_rows}
            if len(values) != 1:
                raise RuntimeError(
                    f"inconsistent event field {field} for {case_id}/{event_id}"
                )
        multiplicity = int(first["multiplicity"])
        nullity = int(first["detected_nullity"])
        slots = tuple(sorted(int(row["sorted_slot"]) for row in event_rows))
        repeated_slots = tuple(
            sorted(int(row["repeated_root_slot"]) for row in event_rows)
        )
        if multiplicity <= 0 or nullity <= 0 or multiplicity != nullity:
            raise RuntimeError(
                f"invalid multiplicity/nullity for {case_id}/{event_id}: "
                f"{multiplicity}/{nullity}"
            )
        if len(event_rows) != multiplicity:
            raise RuntimeError(
                f"truncated repeated slots for {case_id}/{event_id}: "
                f"observed {len(event_rows)}, declared {multiplicity}"
            )
        if repeated_slots != tuple(range(1, multiplicity + 1)):
            raise RuntimeError(
                f"non-contiguous repeated_root_slot for {case_id}/{event_id}"
            )
        if slots != tuple(range(min(slots), min(slots) + multiplicity)):
            raise RuntimeError(f"non-contiguous event slots for {case_id}/{event_id}")
        Omega = float(first["Omega"])
        omega = float(first["omega"])
        cluster_center = float(first["cluster_center_Omega"])
        if not all(
            math.isfinite(value) and value > 0.0
            for value in (Omega, omega, cluster_center)
        ):
            raise RuntimeError(f"non-positive or non-finite root for {case_id}/{event_id}")
        event = RootEventRecord(
            case_id=case_id,
            event_id=event_id,
            Omega=Omega,
            omega=omega,
            multiplicity=multiplicity,
            nullity=nullity,
            cluster_id=str(first.get("cluster_id", "")),
            cluster_semantics=str(first.get("cluster_semantics", "ISOLATED")),
            cluster_multiplicity=int(first.get("cluster_multiplicity", "1") or 1),
            cluster_total_nullity=int(first.get("cluster_total_nullity", "1") or 1),
            cluster_center_Omega=cluster_center,
            slots=slots,
            roles=tuple(str(row["role"]) for row in sorted(event_rows, key=lambda item: int(item["sorted_slot"]))),
        )
        by_case.setdefault(case_id, []).append(event)
    for case_id, events in by_case.items():
        events.sort(key=lambda event: (min(event.slots), event.Omega, event.event_id))
        observed_slots = sorted(slots_by_case[case_id])
        if observed_slots != list(range(1, len(observed_slots) + 1)):
            raise RuntimeError(f"non-contiguous sorted slots in {case_id}")
        if len(observed_slots) < 13:
            raise RuntimeError(f"{case_id} inventory does not contain the root-13 guard")
        slot_roles = {
            slot: role for event in events for slot, role in zip(event.slots, event.roles)
        }
        for slot, role in slot_roles.items():
            expected_role = (
                "FIRST_12"
                if slot <= 12
                else "ROOT_13_GUARD" if slot == 13 else "GUARD_CLUSTER_COMPLETION"
            )
            if role != expected_role:
                raise RuntimeError(
                    f"invalid role {role!r} for {case_id} sorted slot {slot}"
                )

        clusters: dict[str, list[RootEventRecord]] = {}
        for event in events:
            if not event.cluster_id:
                if not (
                    event.cluster_semantics == "ISOLATED"
                    and event.multiplicity == 1
                    and event.nullity == 1
                    and event.cluster_multiplicity == 1
                    and event.cluster_total_nullity == 1
                    and event.cluster_center_Omega == event.Omega
                ):
                    raise RuntimeError(
                        f"inconsistent isolated-root metadata for {case_id}/{event.event_id}"
                    )
                continue
            clusters.setdefault(event.cluster_id, []).append(event)

        for cluster_id, members in clusters.items():
            semantics = {event.cluster_semantics for event in members}
            declared_multiplicity = {event.cluster_multiplicity for event in members}
            declared_nullity = {event.cluster_total_nullity for event in members}
            declared_centres = {event.cluster_center_Omega for event in members}
            if not (
                len(semantics) == 1
                and len(declared_multiplicity) == 1
                and len(declared_nullity) == 1
                and len(declared_centres) == 1
            ):
                raise RuntimeError(f"inconsistent metadata in {case_id}/{cluster_id}")
            observed_multiplicity = sum(event.multiplicity for event in members)
            observed_nullity = sum(event.nullity for event in members)
            if observed_multiplicity != next(iter(declared_multiplicity)):
                raise RuntimeError(f"truncated multiplicity in {case_id}/{cluster_id}")
            if observed_nullity != next(iter(declared_nullity)):
                raise RuntimeError(f"truncated nullity in {case_id}/{cluster_id}")
            centre = math.fsum(
                event.Omega * event.multiplicity for event in members
            ) / observed_multiplicity
            declared_centre = next(iter(declared_centres))
            if not math.isclose(
                centre,
                declared_centre,
                rel_tol=8.0 * np.finfo(float).eps,
                abs_tol=8.0 * np.finfo(float).eps * max(1.0, abs(centre)),
            ):
                raise RuntimeError(f"cluster centre mismatch in {case_id}/{cluster_id}")
            semantics_value = next(iter(semantics))
            exact = len(members) == 1 and members[0].multiplicity > 1
            if semantics_value == "EXACT_DEGENERATE_SUBSPACE" and not exact:
                raise RuntimeError(f"invalid exact-cluster semantics in {case_id}/{cluster_id}")
            if semantics_value == "NEAR_DEGENERATE_CLUSTER" and len(members) < 2:
                raise RuntimeError(f"invalid near-cluster semantics in {case_id}/{cluster_id}")
            if semantics_value not in {
                "EXACT_DEGENERATE_SUBSPACE",
                "NEAR_DEGENERATE_CLUSTER",
            }:
                raise RuntimeError(f"unknown cluster semantics in {case_id}/{cluster_id}")

        if len(observed_slots) > 13:
            guard_event = next(event for event in events if 13 in event.slots)
            completion_events = [
                event for event in events if any(slot > 13 for slot in event.slots)
            ]
            if not guard_event.cluster_id or any(
                event.cluster_id != guard_event.cluster_id for event in completion_events
            ):
                raise RuntimeError(f"invalid guard-cluster completion in {case_id}")
    return by_case


def _semantic_inventory_hash(
    rows: Sequence[Mapping[str, str]], manifest: Mapping[str, Any]
) -> str:
    by_case: dict[str, list[Mapping[str, str]]] = {}
    for row in rows:
        by_case.setdefault(str(row["case_id"]), []).append(row)
    manifest_cases = list(manifest.get("cases", []))
    if tuple(str(case["case_id"]) for case in manifest_cases) != REQUIRED_CASE_IDS:
        raise RuntimeError("worker manifest case order differs from the frozen contract")
    cases_payload: list[dict[str, Any]] = []
    for case_metadata in manifest_cases:
        case_id = str(case_metadata["case_id"])
        case_rows = sorted(by_case.get(case_id, []), key=lambda row: int(row["sorted_slot"]))
        if not case_rows:
            raise RuntimeError(f"frozen inventory has no rows for {case_id}")
        case_hashes = {str(row["case_inventory_sha256"]) for row in case_rows}
        if len(case_hashes) != 1:
            raise RuntimeError(f"inconsistent case inventory hashes for {case_id}")
        case_hash = next(iter(case_hashes))
        if case_hash != str(case_metadata["inventory_sha256"]):
            raise RuntimeError(f"root CSV and manifest case hashes differ for {case_id}")
        cases_payload.append(
            {
                "case_id": case_id,
                "case_inventory_sha256": case_hash,
                "status": str(case_metadata["inventory_status"]),
                "primary_slot_count": int(
                    case_metadata["primary"][
                        "slot_count_including_guard_cluster_completion"
                    ]
                ),
                "verification_slot_count": int(
                    case_metadata["verification"][
                        "slot_count_including_guard_cluster_completion"
                    ]
                ),
            }
        )
    payload = {
        "contract_sha256": str(manifest["contract_sha256"]),
        "builder_id": str(manifest["builder_id"]),
        "root_policy": manifest["root_policy"],
        "cases": cases_payload,
    }
    return _sha256_bytes(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def _spectrum_blocks(events: Sequence[RootEventRecord]) -> list[SpectrumBlock]:
    raw_groups: list[list[RootEventRecord]] = []
    cluster_group: dict[str, list[RootEventRecord]] = {}
    for event in events:
        if event.cluster_id:
            cluster_group.setdefault(event.cluster_id, []).append(event)
        else:
            raw_groups.append([event])
    raw_groups.extend(cluster_group.values())
    raw_groups.sort(key=lambda group: min(min(event.slots) for event in group))
    blocks: list[SpectrumBlock] = []
    for index, group in enumerate(raw_groups, start=1):
        ordered = tuple(sorted(group, key=lambda event: (event.Omega, event.event_id)))
        slots = sorted(slot for event in ordered for slot in event.slots)
        multiplicity = sum(event.multiplicity for event in ordered)
        nullity = sum(event.nullity for event in ordered)
        centre_omega = math.fsum(event.omega * event.multiplicity for event in ordered) / multiplicity
        centre_Omega = math.fsum(event.Omega * event.multiplicity for event in ordered) / multiplicity
        clustered = multiplicity > 1 or len(ordered) > 1
        blocks.append(
            SpectrumBlock(
                case_id=ordered[0].case_id,
                block_id=f"block_{index:04d}",
                events=ordered,
                slot_start=min(slots),
                slot_end=max(slots),
                root_count=len(slots),
                multiplicity=multiplicity,
                nullity=nullity,
                centre_omega=centre_omega,
                centre_Omega=centre_Omega,
                clustered=clustered,
                cluster_semantics=(
                    ordered[0].cluster_semantics if clustered else "ISOLATED"
                ),
                roles=tuple(role for event in ordered for role in event.roles),
            )
        )
    updated: list[SpectrumBlock] = []
    for index, block in enumerate(blocks):
        before = (
            block.centre_omega - blocks[index - 1].centre_omega
            if index > 0
            else math.inf
        )
        after = (
            blocks[index + 1].centre_omega - block.centre_omega
            if index + 1 < len(blocks)
            else math.inf
        )
        updated.append(replace(block, gap_before_omega=before, gap_after_omega=after))
    return updated


def _block_row(
    comparison_kind: str,
    case_id: str,
    left_label: str,
    right_label: str,
    left: SpectrumBlock | None,
    right: SpectrumBlock | None,
    *,
    isolated_tolerance: float,
    cluster_tolerance: float,
) -> dict[str, Any]:
    if left is None or right is None:
        return {
            "comparison_kind": comparison_kind,
            "case_id": case_id,
            "left_model": left_label,
            "right_model": right_label,
            "left_block": "" if left is None else left.block_id,
            "right_block": "" if right is None else right.block_id,
            "status": "FAIL",
            "reason": "UNMATCHED_SPECTRUM_BLOCK",
        }
    clustered = left.clustered or right.clustered
    tolerance = cluster_tolerance if clustered else isolated_tolerance
    relative = _relative_difference(left.centre_omega, right.centre_omega)
    structure_pass = bool(
        left.root_count == right.root_count
        and left.multiplicity == right.multiplicity
        and left.nullity == right.nullity
        and left.clustered == right.clustered
        and left.cluster_semantics == right.cluster_semantics
    )
    passed = structure_pass and relative <= tolerance
    return {
        "comparison_kind": comparison_kind,
        "case_id": case_id,
        "left_model": left_label,
        "right_model": right_label,
        "left_block": left.block_id,
        "right_block": right.block_id,
        "left_slot_start": left.slot_start,
        "left_slot_end": left.slot_end,
        "right_slot_start": right.slot_start,
        "right_slot_end": right.slot_end,
        "left_omega": left.centre_omega,
        "right_omega": right.centre_omega,
        "relative_difference": relative,
        "left_neighbor_gap_before": left.gap_before_omega,
        "left_neighbor_gap_after": left.gap_after_omega,
        "right_neighbor_gap_before": right.gap_before_omega,
        "right_neighbor_gap_after": right.gap_after_omega,
        "left_root_count": left.root_count,
        "right_root_count": right.root_count,
        "left_multiplicity": left.multiplicity,
        "right_multiplicity": right.multiplicity,
        "left_nullity": left.nullity,
        "right_nullity": right.nullity,
        "left_cluster_semantics": left.cluster_semantics,
        "right_cluster_semantics": right.cluster_semantics,
        "guard_or_completion": bool(
            any(role != "FIRST_12" for role in (*left.roles, *right.roles))
        ),
        "tolerance": tolerance,
        "status": "PASS" if passed else "FAIL",
        "reason": "" if passed else (
            "SPECTRUM_STRUCTURE_MISMATCH" if not structure_pass else "FREQUENCY_GATE_FAIL"
        ),
    }


def _compare_block_lists(
    comparison_kind: str,
    case_id: str,
    left_label: str,
    right_label: str,
    left_blocks: Sequence[SpectrumBlock],
    right_blocks: Sequence[SpectrumBlock],
    *,
    isolated_tolerance: float,
    cluster_tolerance: float,
) -> list[dict[str, Any]]:
    # Align declared multiplicity slots, not block-array indices.  If one side
    # groups slots 4--5 as a cluster while the other reports two isolated
    # blocks, this emits a localized structural mismatch and resumes at slot 6
    # instead of shifting every subsequent comparison.
    maximum_slot = max(
        max((block.slot_end for block in left_blocks), default=0),
        max((block.slot_end for block in right_blocks), default=0),
    )
    rows: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for slot in range(1, maximum_slot + 1):
        left = next(
            (block for block in left_blocks if block.slot_start <= slot <= block.slot_end),
            None,
        )
        right = next(
            (block for block in right_blocks if block.slot_start <= slot <= block.slot_end),
            None,
        )
        if left is None and right is None:
            continue
        pair = (
            "" if left is None else left.block_id,
            "" if right is None else right.block_id,
        )
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        row = _block_row(
            comparison_kind,
            case_id,
            left_label,
            right_label,
            left,
            right,
            isolated_tolerance=isolated_tolerance,
            cluster_tolerance=cluster_tolerance,
        )
        row["alignment_slot"] = slot
        rows.append(row)
    return rows


@dataclass(frozen=True)
class SmallRootEvent:
    Omega: float
    omega: float
    multiplicity: int
    nullity: int
    sources: tuple[str, ...]


def _positive_scale(matrix: FloatArray) -> tuple[FloatArray, FloatArray, FloatArray]:
    raw = np.asarray(matrix, dtype=float)
    row_norms = np.linalg.norm(raw, axis=1)
    row_reference = float(np.max(row_norms)) if row_norms.size else 0.0
    row_factors = np.ones(raw.shape[0])
    if row_reference > 0.0:
        row_factors = 1.0 / np.maximum(
            row_norms, float(np.finfo(float).eps ** 0.25) * row_reference
        )
    row_scaled = row_factors[:, None] * raw
    column_factors = 1.0 / np.maximum(np.linalg.norm(row_scaled, axis=0), 1.0)
    return row_scaled * column_factors[None, :], row_factors, column_factors


def _small_diagnostics(
    Omega: float, provider: MatrixProvider, frequency_scale: float, nullity_threshold: float
) -> tuple[float, float, int]:
    matrix = np.asarray(provider(float(Omega) / frequency_scale), dtype=float)
    scaled, _rows, _columns = _positive_scale(matrix)
    sign, logabs = np.linalg.slogdet(scaled)
    determinant = (
        float(sign * math.exp(float(logabs) / scaled.shape[0]))
        if sign != 0.0 and math.isfinite(float(logabs))
        else 0.0
    )
    singular = np.linalg.svd(scaled, compute_uv=False)
    sigma_ratio = float(singular[-1] / singular[0]) if singular[0] > 0.0 else 0.0
    nullity = int(np.count_nonzero(singular <= nullity_threshold * singular[0]))
    return determinant, sigma_ratio, nullity


def _direct_candidate_is_canonical_root(
    sources: Iterable[str], sigma_ratio: float, nullity: int, threshold: float
) -> bool:
    """Accept direct-control roots without duplicating a simple root basin.

    A simple zero of a real square determinant changes determinant sign.  The
    sigma-minimum detector is retained to find even/multiple roots, but a
    sigma-only nullity-one point near an already bracketed simple zero is not a
    second eigenfrequency.  Close *distinct* simple roots remain supported by
    their separate local determinant brackets.
    """

    determinant_supported = bool(
        {
            "determinant_sign_bracket",
            "determinant_zero_grid_certified",
            "determinant_close_pair_sign_bracket",
            "determinant_close_pair_zero_grid_certified",
        }
        & {str(source) for source in sources}
    )
    return bool(
        math.isfinite(float(sigma_ratio))
        and float(sigma_ratio) <= float(threshold)
        and int(nullity) >= 1
        and (determinant_supported or int(nullity) >= 2)
    )


def _certified_zero_grid_candidates(
    grid: FloatArray,
    determinants: FloatArray,
    *,
    source_prefix: str,
) -> tuple[list[tuple[float, str]], int]:
    """Return only isolated grid zeros certified by opposite neighbouring signs.

    A run of two or more numerical zeros is a zero plateau, not evidence of a
    simple determinant crossing.  Repeated physical roots remain discoverable
    through the independent sigma/nullity route.
    """

    values = np.asarray(grid, dtype=float)
    signs = np.asarray(determinants, dtype=float)
    if values.ndim != 1 or signs.shape != values.shape:
        raise ValueError("grid and determinant samples must be one-dimensional peers")
    candidates: list[tuple[float, str]] = []
    rejected_plateaus = 0
    index = 0
    while index < values.size:
        if signs[index] != 0.0:
            index += 1
            continue
        start = index
        while index + 1 < values.size and signs[index + 1] == 0.0:
            index += 1
        stop = index
        left = start - 1
        right = stop + 1
        isolated = start == stop
        bracketed = bool(
            isolated
            and left >= 0
            and right < values.size
            and math.isfinite(float(signs[left]))
            and math.isfinite(float(signs[right]))
            and signs[left] != 0.0
            and signs[right] != 0.0
            and signs[left] * signs[right] < 0.0
        )
        if bracketed:
            candidates.append(
                (float(values[start]), f"{source_prefix}_zero_grid_certified")
            )
        else:
            rejected_plateaus += 1
        index += 1
    return candidates, rejected_plateaus


def _canonicalize_small_root_basins(
    events: Sequence[SmallRootEvent],
    diagnostic: Callable[[float], tuple[float, float, int]],
    singular_threshold: float,
) -> list[SmallRootEvent]:
    """Collapse numerical sign jitter inside one connected singular basin.

    Candidate points are merged only when the complete sampled interval
    between them remains inside the already frozen ``sigma_min/sigma_max``
    root gate.  Two genuinely distinct close roots are preserved when the
    matrix becomes nonsingular between their zeros.  The representative is the
    member with the smallest singular ratio; multiplicity is never summed
    across numerical duplicates.
    """

    ordered = sorted(events, key=lambda event: event.Omega)
    if not ordered:
        return []
    groups: list[list[SmallRootEvent]] = [[ordered[0]]]
    for event in ordered[1:]:
        previous = groups[-1][-1]
        probes = np.linspace(previous.Omega, event.Omega, 5, dtype=float)[1:-1]
        connected = all(
            diagnostic(float(probe))[1] <= float(singular_threshold)
            for probe in probes
        )
        if connected:
            groups[-1].append(event)
        else:
            groups.append([event])

    canonical: list[SmallRootEvent] = []
    for group in groups:
        representative = min(group, key=lambda event: diagnostic(event.Omega)[1])
        maximum_nullity = max(event.nullity for event in group)
        canonical.append(
            SmallRootEvent(
                Omega=representative.Omega,
                omega=representative.omega,
                multiplicity=maximum_nullity,
                nullity=maximum_nullity,
                sources=tuple(
                    sorted({source for event in group for source in event.sources})
                ),
            )
        )
    return canonical


def _scan_small_matrix(
    provider: MatrixProvider,
    frequency_scale: float,
    contract: Mapping[str, Any],
    *,
    points: int,
    phase: float,
    audit: dict[str, Any] | None = None,
) -> list[SmallRootEvent]:
    policy = contract["root_policy"]
    thresholds = contract["thresholds"]
    lower = float(policy["Omega_min"])
    upper = float(policy["Omega_max"])
    step = (upper - lower) / (points - 1)
    start = lower + phase * step
    count = int(math.floor((upper - start) / step)) + 1
    grid = start + step * np.arange(count)
    cache: dict[float, tuple[float, float, int]] = {}

    def diagnostic(value: float) -> tuple[float, float, int]:
        key = float(value)
        if key not in cache:
            cache[key] = _small_diagnostics(
                key,
                provider,
                frequency_scale,
                float(policy["nullity_relative_threshold"]),
            )
        return cache[key]

    determinants = np.empty(grid.size)
    sigmas = np.empty(grid.size)
    for index, value in enumerate(grid):
        determinants[index], sigmas[index], _nullity = diagnostic(float(value))
    candidates, rejected_global_zero_plateaus = _certified_zero_grid_candidates(
        grid,
        determinants,
        source_prefix="determinant",
    )
    rejected_local_zero_plateaus = 0
    for index in range(grid.size - 1):
        left, right = float(grid[index]), float(grid[index + 1])
        f_left, f_right = determinants[index], determinants[index + 1]
        if (
            f_left != 0.0
            and f_right != 0.0
            and math.isfinite(f_left)
            and math.isfinite(f_right)
            and f_left * f_right < 0.0
        ):
            try:
                root = brentq(
                    lambda value: diagnostic(float(value))[0],
                    left,
                    right,
                    xtol=float(policy["root_xtol_Omega"]),
                    rtol=8.0 * np.finfo(float).eps,
                    maxiter=180,
                )
                candidates.append((float(root), "determinant_sign_bracket"))
            except (ValueError, RuntimeError):
                pass
    determinant_locations = sorted(value for value, _source in candidates)
    required_slots = 13
    if len(determinant_locations) >= required_slots:
        provisional_guard = determinant_locations[required_slots - 1]
        refinement_cutoff = min(
            upper,
            provisional_guard
            + float(policy["post_guard_tail_Omega"])
            + 2.0 * step,
        )
    else:
        provisional_guard = math.nan
        refinement_cutoff = upper
    prefilter = float(policy["sigma_prefilter"])
    for index in range(1, grid.size - 1):
        if float(grid[index - 1]) > refinement_cutoff:
            break
        if (
            sigmas[index] <= sigmas[index - 1]
            and sigmas[index] <= sigmas[index + 1]
            and sigmas[index] <= prefilter
        ):
            left, right = float(grid[index - 1]), float(grid[index + 1])
            guard_subintervals = int(
                policy.get("local_close_pair_guard_subintervals", 0)
            )
            if guard_subintervals > 0:
                if guard_subintervals < 2:
                    raise RuntimeError(
                        "local_close_pair_guard_subintervals must be zero or >= 2"
                    )
                local_grid = np.linspace(
                    left, right, guard_subintervals + 1, dtype=float
                )
                local_determinants = np.array(
                    [diagnostic(float(value))[0] for value in local_grid],
                    dtype=float,
                )
                local_sigmas = np.array(
                    [diagnostic(float(value))[1] for value in local_grid],
                    dtype=float,
                )
                certified_local_zeros, local_plateaus = (
                    _certified_zero_grid_candidates(
                        local_grid,
                        local_determinants,
                        source_prefix="determinant_close_pair",
                    )
                )
                candidates.extend(certified_local_zeros)
                rejected_local_zero_plateaus += local_plateaus
                for local_index in range(local_grid.size - 1):
                    local_left = float(local_grid[local_index])
                    local_right = float(local_grid[local_index + 1])
                    local_f_left = float(local_determinants[local_index])
                    local_f_right = float(local_determinants[local_index + 1])
                    if (
                        local_f_left != 0.0
                        and local_f_right != 0.0
                        and math.isfinite(local_f_left)
                        and math.isfinite(local_f_right)
                        and local_f_left * local_f_right < 0.0
                    ):
                        try:
                            local_root = brentq(
                                lambda value: diagnostic(float(value))[0],
                                local_left,
                                local_right,
                                xtol=float(policy["root_xtol_Omega"]),
                                rtol=8.0 * np.finfo(float).eps,
                                maxiter=180,
                            )
                            candidates.append(
                                (
                                    float(local_root),
                                    "determinant_close_pair_sign_bracket",
                                )
                            )
                        except (ValueError, RuntimeError):
                            pass
                for local_index in range(1, local_grid.size - 1):
                    if (
                        local_sigmas[local_index]
                        <= local_sigmas[local_index - 1]
                        and local_sigmas[local_index]
                        <= local_sigmas[local_index + 1]
                        and local_sigmas[local_index] <= prefilter
                    ):
                        local_result = minimize_scalar(
                            lambda value: diagnostic(float(value))[1] ** 2,
                            bounds=(
                                float(local_grid[local_index - 1]),
                                float(local_grid[local_index + 1]),
                            ),
                            method="bounded",
                            options={
                                "xatol": float(policy["root_xtol_Omega"]),
                                "maxiter": 180,
                            },
                        )
                        if local_result.success:
                            candidates.append(
                                (
                                    float(local_result.x),
                                    "sigma_close_pair_minimum",
                                )
                            )
            result = minimize_scalar(
                lambda value: diagnostic(float(value))[1] ** 2,
                bounds=(left, right),
                method="bounded",
                options={"xatol": float(policy["root_xtol_Omega"]), "maxiter": 180},
            )
            if result.success:
                candidates.append((float(result.x), "sigma_ratio_minimum"))
    candidates.sort(key=lambda item: item[0])
    deduplicated: list[tuple[float, set[str]]] = []
    for value, source in candidates:
        if deduplicated:
            previous = deduplicated[-1][0]
            tolerance = float(policy["dedup_atol_Omega"]) + float(policy["dedup_rtol"]) * max(
                abs(previous), abs(value)
            )
            if abs(value - previous) <= tolerance:
                old_value, sources = deduplicated[-1]
                if diagnostic(value)[1] < diagnostic(old_value)[1]:
                    deduplicated[-1] = (value, sources | {source})
                else:
                    sources.add(source)
                continue
        deduplicated.append((value, {source}))
    accepted: list[SmallRootEvent] = []
    for value, sources in deduplicated:
        _determinant, sigma, nullity = diagnostic(value)
        if _direct_candidate_is_canonical_root(
            sources,
            sigma,
            nullity,
            float(thresholds["root_singular_ratio"]),
        ):
            accepted.append(
                SmallRootEvent(
                    Omega=value,
                    omega=value / frequency_scale,
                    multiplicity=nullity,
                    nullity=nullity,
                    sources=tuple(sorted(sources)),
                )
            )
    canonical = _canonicalize_small_root_basins(
        accepted,
        diagnostic,
        float(thresholds["root_singular_ratio"]),
    )
    if audit is not None:
        audit.update(
            {
                "sample_count": int(grid.size),
                "raw_candidate_count": len(candidates),
                "deduplicated_candidate_count": len(deduplicated),
                "accepted_before_basin_canonicalization": len(accepted),
                "canonical_event_count": len(canonical),
                "rejected_global_zero_plateau_count": rejected_global_zero_plateaus,
                "rejected_local_zero_plateau_count": rejected_local_zero_plateaus,
                "determinant_provisional_guard_Omega": provisional_guard,
                "sigma_refinement_cutoff_Omega": refinement_cutoff,
            }
        )
    return canonical


def _small_events_to_blocks(
    case_id: str, events: Sequence[SmallRootEvent], contract: Mapping[str, Any]
) -> list[SpectrumBlock]:
    policy = contract["root_policy"]
    ordered = sorted(events, key=lambda event: event.Omega)
    groups: list[list[SmallRootEvent]] = []
    for event in ordered:
        if not groups:
            groups.append([event])
            continue
        previous = groups[-1][-1]
        tolerance = float(policy["cluster_atol_Omega"]) + float(policy["cluster_rtol"]) * max(
            abs(previous.Omega), abs(event.Omega)
        )
        if event.Omega - previous.Omega <= tolerance:
            groups[-1].append(event)
        else:
            groups.append([event])
    root_events: list[RootEventRecord] = []
    slot = 0
    for index, group in enumerate(groups, start=1):
        total = sum(event.multiplicity for event in group)
        cluster_id = f"direct_cluster_{index:04d}" if total > 1 else ""
        spread = max(event.Omega for event in group) - min(
            event.Omega for event in group
        )
        exact_resolution = float(policy["root_xtol_Omega"])
        exact_subspace = bool(
            total > 1 and (len(group) == 1 or spread <= exact_resolution)
        )
        semantics = (
            "EXACT_DEGENERATE_SUBSPACE"
            if exact_subspace
            else "NEAR_DEGENERATE_CLUSTER" if total > 1 else "ISOLATED"
        )
        for event_index, event in enumerate(group, start=1):
            slots = tuple(range(slot + 1, slot + event.multiplicity + 1))
            slot += event.multiplicity
            roles = tuple("FIRST_12" if number <= 12 else "ROOT_13_GUARD" for number in slots)
            root_events.append(
                RootEventRecord(
                    case_id=case_id,
                    event_id=f"direct_event_{index:04d}_{event_index:02d}",
                    Omega=event.Omega,
                    omega=event.omega,
                    multiplicity=event.multiplicity,
                    nullity=event.nullity,
                    cluster_id=cluster_id,
                    cluster_semantics=semantics,
                    cluster_multiplicity=total,
                    cluster_total_nullity=sum(item.nullity for item in group),
                    cluster_center_Omega=math.fsum(
                        item.Omega * item.multiplicity for item in group
                    ) / total,
                    slots=slots,
                    roles=roles,
                )
            )
        if slot >= 13 and (not cluster_id or group is groups[-1]):
            # Continue building; truncation below keeps a complete guard cluster.
            pass
    blocks = _spectrum_blocks(root_events)
    retained: list[SpectrumBlock] = []
    for block in blocks:
        if block.slot_start <= 13:
            retained.append(block)
        elif retained and retained[-1].clustered and block.block_id == retained[-1].block_id:
            retained.append(block)
    return retained


def _verified_small_spectrum(
    provider: MatrixProvider,
    frequency_scale: float,
    contract: Mapping[str, Any],
    *,
    case_id: str,
) -> tuple[list[SpectrumBlock], dict[str, Any]]:
    policy = contract["root_policy"]
    primary_scan_audit: dict[str, Any] = {}
    verification_scan_audit: dict[str, Any] = {}
    primary_events = _scan_small_matrix(
        provider,
        frequency_scale,
        contract,
        points=int(policy["primary_scan_points"]),
        phase=float(policy["primary_phase"]),
        audit=primary_scan_audit,
    )
    verification_events = _scan_small_matrix(
        provider,
        frequency_scale,
        contract,
        points=int(policy["verification_scan_points"]),
        phase=float(policy["verification_phase"]),
        audit=verification_scan_audit,
    )
    primary = _small_events_to_blocks(case_id, primary_events, contract)
    verification = _small_events_to_blocks(case_id, verification_events, contract)
    rows = _compare_block_lists(
        "DIRECT_PRIMARY_VS_VERIFICATION",
        case_id,
        "direct_primary",
        "direct_verification",
        primary,
        verification,
        isolated_tolerance=1.0e-8,
        cluster_tolerance=1.0e-8,
    )
    primary_slot_count = sum(block.root_count for block in primary)
    guard_not_at_boundary = bool(
        primary_slot_count >= 13
        and primary[-1].centre_Omega
        <= float(policy["Omega_max"]) - float(policy["post_guard_tail_Omega"])
    )
    status = bool(
        len(primary) > 0
        and primary_slot_count >= 13
        and guard_not_at_boundary
        and all(row["status"] == "PASS" for row in rows)
    )
    return primary, {
        "status": "PASS" if status else "FAIL",
        "primary_verification_rows": rows,
        "primary_event_count": len(primary_events),
        "verification_event_count": len(verification_events),
        "primary_scan_audit": primary_scan_audit,
        "verification_scan_audit": verification_scan_audit,
        "primary_slot_count_through_guard": primary_slot_count,
        "guard_not_at_scan_boundary": guard_not_at_boundary,
    }


def _frequency_scale(contract: Mapping[str, Any], geometry_id: str) -> float:
    material = contract["material"]
    geometry_data = contract["geometries"][geometry_id]
    area = float(geometry_data["width"]) * float(geometry_data["thickness"])
    inertia = float(geometry_data["width"]) * float(geometry_data["thickness"]) ** 3 / 12.0
    return float(contract["lengths"]["L_ref"]) ** 2 * math.sqrt(
        float(material["rho"]) * area / (float(material["E"]) * inertia)
    )


def _rectangular_section(contract: Mapping[str, Any], geometry_id: str) -> legacy.SectionProperties:
    material = contract["material"]
    geometry_data = contract["geometries"][geometry_id]
    return legacy.rectangular_section(
        E=float(material["E"]),
        nu=float(material["nu"]),
        rho=float(material["rho"]),
        width=float(geometry_data["width"]),
        thickness=float(geometry_data["thickness"]),
        K=float(material["K"]),
    )


def _direct_fixed_fixed_provider(
    section: legacy.SectionProperties, length: float
) -> MatrixProvider:
    def provider(omega: float) -> FloatArray:
        columns = legacy.clamped_endpoint_columns(length, omega, section)
        return np.vstack([columns["u"], columns["w"], columns["psi"]])

    return provider


def _direct_bending_provider(
    section: legacy.SectionProperties, length: float
) -> MatrixProvider:
    def provider(omega: float) -> FloatArray:
        columns = legacy.clamped_endpoint_columns(length, omega, section)
        return np.vstack([columns["w"][:2], columns["psi"][:2]])

    return provider


def _union_family_blocks(
    case_id: str,
    section: legacy.SectionProperties,
    length: float,
    frequency_scale: float,
    bending_blocks: Sequence[SpectrumBlock],
    contract: Mapping[str, Any],
) -> list[SpectrumBlock]:
    events: list[SmallRootEvent] = []
    for block in bending_blocks:
        for event in block.events:
            events.append(
                SmallRootEvent(
                    Omega=event.Omega,
                    omega=event.omega,
                    multiplicity=event.multiplicity,
                    nullity=event.nullity,
                    sources=("bending_fixed_fixed",),
                )
            )
    axial_step = math.pi / length * math.sqrt(section.EA / section.rhoA)
    n = 1
    while frequency_scale * axial_step * n <= float(contract["root_policy"]["Omega_max"]):
        omega = axial_step * n
        Omega = frequency_scale * omega
        if Omega >= float(contract["root_policy"]["Omega_min"]):
            events.append(
                SmallRootEvent(
                    Omega=Omega,
                    omega=omega,
                    multiplicity=1,
                    nullity=1,
                    sources=(f"axial_fixed_fixed_n={n}",),
                )
            )
        n += 1
    return _small_events_to_blocks(case_id, events, contract)


def _rlb_properties(contract: Mapping[str, Any], geometry_id: str) -> rlb.BeamProperties:
    material_data = contract["material"]
    geometry_data = contract["geometries"][geometry_id]
    elastic = float(material_data["E"])
    poisson = float(material_data["nu"])
    shear = elastic / (2.0 * (1.0 + poisson))
    material = rlb.LaminaMaterial(
        E1=elastic,
        E2=elastic,
        nu12=poisson,
        G12=shear,
        G13=shear,
        G23=shear,
        rho=float(material_data["rho"]),
        name="RLB-1C-ISO postprocess isotropic lamina",
    )
    thickness = float(geometry_data["thickness"])
    angles = tuple(float(value) for value in contract["spectral_stack_deg"])
    if len(angles) != 4:
        raise RuntimeError("postprocess contract must retain exactly four plies")
    plies = tuple(
        rlb.Ply(material, angle, thickness / 4.0, label=f"ply-{index}")
        for index, angle in enumerate(angles, start=1)
    )
    laminate = rlb.integrate_laminate(plies)
    return rlb.reduce_to_beam_properties(
        laminate,
        width=float(geometry_data["width"]),
        K=float(material_data["K"]),
        symmetry_tolerance=1.0e-12,
        reduction_tolerance=1.0e-12,
    )


@dataclass(frozen=True)
class ModeSeed:
    method: str
    case_id: str
    event_id: str
    omega: float
    vector_index: int
    coefficient: FloatArray
    boundary_null_residual: float
    scaled_sigma_ratio: float


@dataclass(frozen=True)
class PhysicalMode:
    method: str
    case_id: str
    event_id: str
    omega: float
    vector_index: int
    x1: FloatArray
    x2: FloatArray
    state1: FloatArray
    state2: FloatArray
    derivative1: FloatArray
    derivative2: FloatArray
    boundary_null_residual: float
    scaled_sigma_ratio: float
    normalization_factor: float = 1.0
    deterministic_sign: float = 1.0
    sign_pivot_index: int = -1


def _raw_nullspace_basis(
    matrix: FloatArray, expected_nullity: int, nullity_threshold: float
) -> tuple[FloatArray, float, float]:
    raw = np.asarray(matrix, dtype=float)
    scaled, _row_factors, column_factors = _positive_scale(raw)
    _left, singular, vh = np.linalg.svd(scaled, full_matrices=False)
    detected = int(np.count_nonzero(singular <= nullity_threshold * singular[0]))
    if detected != expected_nullity:
        raise RuntimeError(
            f"mode nullity changed after inventory freeze: expected {expected_nullity}, "
            f"recomputed {detected}"
        )
    scaled_basis = vh[-detected:, :].T
    basis = column_factors[:, None] * scaled_basis
    for index in range(basis.shape[1]):
        norm = float(np.linalg.norm(basis[:, index]))
        if norm == 0.0 or not math.isfinite(norm):
            raise RuntimeError("recovered determinant null vector is invalid")
        basis[:, index] /= norm
        pivot = int(np.argmax(np.abs(basis[:, index])))
        if basis[pivot, index] < 0.0:
            basis[:, index] *= -1.0
    raw_norm = float(np.linalg.norm(raw, ord=2))
    residual = max(
        float(np.linalg.norm(raw @ basis[:, index]))
        / max(raw_norm * float(np.linalg.norm(basis[:, index])), np.finfo(float).tiny)
        for index in range(basis.shape[1])
    )
    sigma_ratio = float(singular[-1] / singular[0]) if singular[0] > 0.0 else 0.0
    return basis, residual, sigma_ratio


class _ModeReconstructor:
    def __init__(self, contract: Mapping[str, Any]) -> None:
        self.contract = contract
        self.case_by_id = {str(case["case_id"]): case for case in contract["cases"]}
        self.rlb_properties = {
            geometry_id: _rlb_properties(contract, geometry_id)
            for geometry_id in contract["geometries"]
        }
        self.sections = {
            geometry_id: _rectangular_section(contract, geometry_id)
            for geometry_id in contract["geometries"]
        }
        self.seed_cache: dict[tuple[str, str, str], tuple[ModeSeed, ...]] = {}

    def seeds(self, method: str, event: RootEventRecord) -> tuple[ModeSeed, ...]:
        key = (method, event.case_id, event.event_id)
        if key in self.seed_cache:
            return self.seed_cache[key]
        case = self.case_by_id[event.case_id]
        length_1 = float(case["L1"])
        length_2 = float(case["L2"])
        beta_deg = float(case["beta_deg"])
        if method == "RLB":
            properties = self.rlb_properties[str(case["geometry"])]
            matrix = coupled.coupled_boundary_matrix(
                event.omega,
                math.radians(beta_deg),
                length_1,
                properties,
                length_2,
                properties,
            )
        elif method == "LEGACY_RECTANGULAR":
            section = self.sections[str(case["geometry"])]
            matrix = legacy.legacy_coupled_boundary_matrix_raw(
                event.omega,
                section,
                length_1,
                section,
                length_2,
                beta_deg=beta_deg,
            )
        else:
            raise ValueError(f"unknown mode reconstruction method: {method}")
        basis, residual, sigma_ratio = _raw_nullspace_basis(
            np.asarray(matrix, dtype=float),
            event.nullity,
            float(self.contract["root_policy"]["nullity_relative_threshold"]),
        )
        seeds = tuple(
            ModeSeed(
                method=method,
                case_id=event.case_id,
                event_id=event.event_id,
                omega=event.omega,
                vector_index=index + 1,
                coefficient=np.asarray(basis[:, index], dtype=float),
                boundary_null_residual=residual,
                scaled_sigma_ratio=sigma_ratio,
            )
            for index in range(basis.shape[1])
        )
        self.seed_cache[key] = seeds
        return seeds

    def evaluate(self, seed: ModeSeed, n_points: int) -> PhysicalMode:
        case = self.case_by_id[seed.case_id]
        length_1 = float(case["L1"])
        length_2 = float(case["L2"])
        beta_deg = float(case["beta_deg"])
        points = int(n_points)
        x1 = np.linspace(0.0, length_1, points)
        x2 = np.linspace(0.0, length_2, points)
        if seed.method == "RLB":
            properties = self.rlb_properties[str(case["geometry"])]
            initial_1 = np.asarray(coupled.CLAMP_BASIS @ seed.coefficient[:3], dtype=float)
            initial_2 = np.asarray(coupled.CLAMP_BASIS @ seed.coefficient[3:], dtype=float)
            state1 = np.vstack(
                [
                    coupled.arm_transfer_matrix(seed.omega, float(value), properties)
                    @ initial_1
                    for value in x1
                ]
            )
            state2 = np.vstack(
                [
                    coupled.arm_transfer_matrix(seed.omega, float(value), properties)
                    @ initial_2
                    for value in x2
                ]
            )
            state_matrix = coupled.arm_state_matrix(seed.omega, properties)
            derivative1 = state1 @ state_matrix.T
            derivative2 = state2 @ state_matrix.T
        else:
            section = self.sections[str(case["geometry"])]
            assembly = legacy.assemble_legacy_coupled_boundary(
                seed.omega,
                section,
                length_1,
                section,
                length_2,
                beta_deg=beta_deg,
            )
            old_x2 = -x2
            arm1 = legacy.evaluate_clamped_arm_fields(
                seed.omega, section, assembly.basis_1, seed.coefficient[:3], x1
            )
            arm2 = legacy.evaluate_clamped_arm_fields(
                seed.omega, section, assembly.basis_2, seed.coefficient[3:], old_x2
            )
            old_state1 = np.column_stack(
                [arm1.u, arm1.w, arm1.psi, arm1.N, arm1.Q, arm1.M]
            )
            old_state2 = np.column_stack(
                [arm2.u, arm2.w, arm2.psi, arm2.N, arm2.Q, arm2.M]
            )
            old_derivative1 = np.column_stack(
                [
                    arm1.u_prime,
                    arm1.w_prime,
                    arm1.psi_prime,
                    arm1.N_prime,
                    arm1.Q_prime,
                    arm1.M_prime,
                ]
            )
            old_derivative2 = np.column_stack(
                [
                    arm2.u_prime,
                    arm2.w_prime,
                    arm2.psi_prime,
                    arm2.N_prime,
                    arm2.Q_prime,
                    arm2.M_prime,
                ]
            )
            state1 = old_state1 @ S1_OLD_TO_RLB.T
            state2 = old_state2 @ S2_OLD_TO_RLB.T
            derivative1 = old_derivative1 @ S1_OLD_TO_RLB.T
            derivative2 = old_derivative2 @ (-S2_OLD_TO_RLB).T
        return PhysicalMode(
            method=seed.method,
            case_id=seed.case_id,
            event_id=seed.event_id,
            omega=seed.omega,
            vector_index=seed.vector_index,
            x1=x1,
            x2=x2,
            state1=np.asarray(state1, dtype=float),
            state2=np.asarray(state2, dtype=float),
            derivative1=np.asarray(derivative1, dtype=float),
            derivative2=np.asarray(derivative2, dtype=float),
            boundary_null_residual=seed.boundary_null_residual,
            scaled_sigma_ratio=seed.scaled_sigma_ratio,
        )


def _simpson_weights(length: float, points: int) -> FloatArray:
    if points < 3 or points % 2 == 0:
        raise ValueError("composite Simpson integration needs an odd point count >= 3")
    weights = np.ones(points, dtype=float)
    weights[1:-1:2] = 4.0
    weights[2:-1:2] = 2.0
    weights *= float(length) / (3.0 * (points - 1))
    return weights


def _mass_embedding(mode: PhysicalMode, mass: float, rotary: float) -> FloatArray:
    blocks: list[FloatArray] = []
    for coordinates, states in ((mode.x1, mode.state1), (mode.x2, mode.state2)):
        weights = _simpson_weights(float(coordinates[-1]), coordinates.size)
        blocks.extend(
            [
                np.sqrt(weights * mass) * states[:, 0],
                np.sqrt(weights * mass) * states[:, 1],
                np.sqrt(weights * rotary) * states[:, 2],
            ]
        )
    return np.concatenate(blocks)


def _modal_mass(mode: PhysicalMode, mass: float, rotary: float) -> float:
    vector = _mass_embedding(mode, mass, rotary)
    return float(vector @ vector)


def _scale_mode(
    mode: PhysicalMode, factor: float, *, deterministic_sign: float, pivot: int
) -> PhysicalMode:
    total = float(factor) * float(deterministic_sign)
    return replace(
        mode,
        state1=mode.state1 * total,
        state2=mode.state2 * total,
        derivative1=mode.derivative1 * total,
        derivative2=mode.derivative2 * total,
        normalization_factor=float(factor),
        deterministic_sign=float(deterministic_sign),
        sign_pivot_index=int(pivot),
    )


def _normalize_mode(mode: PhysicalMode, mass: float, rotary: float) -> PhysicalMode:
    modal_mass = _modal_mass(mode, mass, rotary)
    if modal_mass <= 0.0 or not math.isfinite(modal_mass):
        raise RuntimeError("mode has nonpositive or nonfinite mass norm")
    factor = 1.0 / math.sqrt(modal_mass)
    trial = _scale_mode(mode, factor, deterministic_sign=1.0, pivot=-1)
    embedding = _mass_embedding(trial, mass, rotary)
    pivot = int(np.argmax(np.abs(embedding)))
    sign = -1.0 if embedding[pivot] < 0.0 else 1.0
    return _scale_mode(mode, factor, deterministic_sign=sign, pivot=pivot)


def _block_for_slot(blocks: Sequence[SpectrumBlock], slot: int) -> SpectrumBlock | None:
    return next((block for block in blocks if block.slot_start <= slot <= block.slot_end), None)


def _seeds_for_block(
    reconstructor: _ModeReconstructor, method: str, block: SpectrumBlock
) -> tuple[ModeSeed, ...]:
    return tuple(
        seed
        for event in block.events
        for seed in reconstructor.seeds(method, event)
    )


def _kinematic_amplitude(mode: PhysicalMode, reference_length: float) -> float:
    return max(
        float(
            np.max(
                np.sqrt(
                    states[:, 0] ** 2
                    + states[:, 1] ** 2
                    + (reference_length * states[:, 2]) ** 2
                )
            )
        )
        for states in (mode.state1, mode.state2)
    )


def _mode_residual_row(
    mode: PhysicalMode,
    case: Mapping[str, Any],
    properties: rlb.BeamProperties,
    contract: Mapping[str, Any],
    *,
    requested_position: int,
    effective_position: int,
    grid_mass_convergence: float,
) -> dict[str, Any]:
    reference_length = float(contract["lengths"]["L_ref"])
    beta_rad = math.radians(float(case["beta_deg"]))
    joint = coupled.physical_joint_residuals(
        beta_rad, mode.state1[-1], mode.state2[-1]
    )
    kinematic_scale = max(_kinematic_amplitude(mode, reference_length), np.finfo(float).tiny)
    displacement_scale = max(
        float(np.max(np.hypot(mode.state1[:, 0], mode.state1[:, 1]))),
        float(np.max(np.hypot(mode.state2[:, 0], mode.state2[:, 1]))),
        np.finfo(float).tiny,
    )
    rotation_scale = max(
        float(np.max(np.abs(mode.state1[:, 2]))),
        float(np.max(np.abs(mode.state2[:, 2]))),
        kinematic_scale / reference_length,
        np.finfo(float).tiny,
    )
    outer_absolute = max(
        float(
            np.linalg.norm(
                [mode.state1[0, 0], mode.state1[0, 1], reference_length * mode.state1[0, 2]]
            )
        ),
        float(
            np.linalg.norm(
                [mode.state2[0, 0], mode.state2[0, 1], reference_length * mode.state2[0, 2]]
            )
        ),
    )
    force_scale = max(
        float(np.max(np.hypot(mode.state1[:, 3], mode.state1[:, 4]))),
        float(np.max(np.hypot(mode.state2[:, 3], mode.state2[:, 4]))),
        float(np.max(np.abs(mode.state1[:, 5]))) / reference_length,
        float(np.max(np.abs(mode.state2[:, 5]))) / reference_length,
        np.finfo(float).tiny,
    )
    moment_scale = max(
        float(np.max(np.abs(mode.state1[:, 5]))),
        float(np.max(np.abs(mode.state2[:, 5]))),
        reference_length * force_scale,
        np.finfo(float).tiny,
    )
    displacement_absolute = float(np.linalg.norm(joint.displacement))
    rotation_absolute = float(np.linalg.norm(joint.rotation))
    force_absolute = float(np.linalg.norm(joint.force))
    moment_absolute = float(np.linalg.norm(joint.moment))

    omega2 = mode.omega**2
    differential_components = np.zeros(6, dtype=float)
    differential_scales = np.zeros(6, dtype=float)
    stiffness_integral = 0.0
    modal_mass = 0.0
    for coordinates, states, derivatives in (
        (mode.x1, mode.state1, mode.derivative1),
        (mode.x2, mode.state2, mode.derivative2),
    ):
        residuals = np.column_stack(
            [
                derivatives[:, 0] - states[:, 3] / properties.A,
                derivatives[:, 1] - states[:, 4] / properties.S + states[:, 2],
                derivatives[:, 2] - states[:, 5] / properties.D,
                derivatives[:, 3] + properties.m * omega2 * states[:, 0],
                derivatives[:, 4] + properties.m * omega2 * states[:, 1],
                derivatives[:, 5] - states[:, 4] + properties.J * omega2 * states[:, 2],
            ]
        )
        term_sets = (
            (derivatives[:, 0], states[:, 3] / properties.A),
            (derivatives[:, 1], states[:, 4] / properties.S, states[:, 2]),
            (derivatives[:, 2], states[:, 5] / properties.D),
            (derivatives[:, 3], properties.m * omega2 * states[:, 0]),
            (derivatives[:, 4], properties.m * omega2 * states[:, 1]),
            (derivatives[:, 5], states[:, 4], properties.J * omega2 * states[:, 2]),
        )
        differential_components = np.maximum(
            differential_components, np.max(np.abs(residuals), axis=0)
        )
        differential_scales = np.maximum(
            differential_scales,
            np.array(
                [max(float(np.max(np.abs(term))) for term in terms) for terms in term_sets]
            ),
        )
        weights = _simpson_weights(float(coordinates[-1]), coordinates.size)
        modal_mass += float(
            weights
            @ (
                properties.m * (states[:, 0] ** 2 + states[:, 1] ** 2)
                + properties.J * states[:, 2] ** 2
            )
        )
        stiffness_integral += float(
            weights
            @ (
                states[:, 3] ** 2 / properties.A
                + states[:, 4] ** 2 / properties.S
                + states[:, 5] ** 2 / properties.D
            )
        )
    differential_normalized = differential_components / np.maximum(
        differential_scales, np.finfo(float).tiny
    )
    inertia_integral = omega2 * modal_mass
    energy_residual = abs(stiffness_integral - inertia_integral) / max(
        abs(stiffness_integral), abs(inertia_integral), np.finfo(float).tiny
    )
    thresholds = contract["thresholds"]
    displacement_normalized = displacement_absolute / displacement_scale
    rotation_normalized = rotation_absolute / rotation_scale
    force_normalized = force_absolute / force_scale
    moment_normalized = moment_absolute / moment_scale
    joint_compatibility = max(displacement_normalized, rotation_normalized)
    joint_equilibrium = max(force_normalized, moment_normalized)
    outer_normalized = outer_absolute / kinematic_scale
    passed = bool(
        outer_normalized <= float(thresholds["outer_clamp_residual"])
        and joint_compatibility <= float(thresholds["joint_compatibility"])
        and joint_equilibrium <= float(thresholds["joint_equilibrium"])
        and mode.boundary_null_residual <= float(thresholds["boundary_null_residual"])
        and energy_residual <= float(thresholds["energy_identity"])
        and grid_mass_convergence <= float(thresholds["grid_convergence"])
        and np.all(np.isfinite(differential_normalized))
    )
    return {
        "method": mode.method,
        "case_id": mode.case_id,
        "requested_position": requested_position,
        "effective_position": effective_position,
        "event_id": mode.event_id,
        "basis_vector": mode.vector_index,
        "grid_points": mode.x1.size,
        "omega": mode.omega,
        "outer_clamp_absolute": outer_absolute,
        "outer_clamp_normalized": outer_normalized,
        "joint_displacement_components": joint.displacement,
        "joint_rotation_components": joint.rotation,
        "joint_force_components": joint.force,
        "joint_moment_components": joint.moment,
        "joint_displacement_absolute": displacement_absolute,
        "joint_rotation_absolute": rotation_absolute,
        "joint_force_absolute": force_absolute,
        "joint_moment_absolute": moment_absolute,
        "joint_displacement_normalized": displacement_normalized,
        "joint_rotation_normalized": rotation_normalized,
        "joint_force_normalized": force_normalized,
        "joint_moment_normalized": moment_normalized,
        "joint_compatibility_normalized": joint_compatibility,
        "joint_equilibrium_normalized": joint_equilibrium,
        "differential_component_absolute": differential_components,
        "differential_component_normalized": differential_normalized,
        "maximum_differential_normalized": float(np.max(differential_normalized)),
        "differential_gate": "NON_GATING_DIAGNOSTIC_NO_FROZEN_THRESHOLD",
        "boundary_null_residual": mode.boundary_null_residual,
        "scaled_sigma_ratio": mode.scaled_sigma_ratio,
        "modal_mass_after_normalization": modal_mass,
        "stiffness_integral": stiffness_integral,
        "inertia_integral": inertia_integral,
        "energy_identity_relative": energy_residual,
        "grid_mass_convergence_relative": grid_mass_convergence,
        "normalization_factor": mode.normalization_factor,
        "deterministic_sign": mode.deterministic_sign,
        "sign_pivot_index": mode.sign_pivot_index,
        "status": "PASS" if passed else "FAIL",
    }


def _stiffness_embedding(mode: PhysicalMode, properties: rlb.BeamProperties) -> FloatArray:
    blocks: list[FloatArray] = []
    for coordinates, states in ((mode.x1, mode.state1), (mode.x2, mode.state2)):
        weights = _simpson_weights(float(coordinates[-1]), coordinates.size)
        blocks.extend(
            [
                np.sqrt(weights / properties.A) * states[:, 3],
                np.sqrt(weights / properties.S) * states[:, 4],
                np.sqrt(weights / properties.D) * states[:, 5],
            ]
        )
    return np.concatenate(blocks)


def _global_displacements(mode: PhysicalMode, beta_deg: float) -> tuple[FloatArray, FloatArray]:
    bases = geometry.reddy_inplane_geometry(beta_deg)
    displacement_1 = (
        mode.state1[:, 0, None] * bases.arm1.t[None, :]
        + mode.state1[:, 1, None] * bases.arm1.n[None, :]
    )
    displacement_2 = (
        mode.state2[:, 0, None] * bases.arm2.t[None, :]
        + mode.state2[:, 1, None] * bases.arm2.n[None, :]
    )
    return displacement_1, displacement_2


def _global_displacement_difference(
    left: PhysicalMode, right: PhysicalMode, mass: float, beta_deg: float
) -> float:
    left_displacements = _global_displacements(left, beta_deg)
    right_displacements = _global_displacements(right, beta_deg)
    numerator = 0.0
    denominator_left = 0.0
    denominator_right = 0.0
    for coordinates, left_values, right_values in (
        (left.x1, left_displacements[0], right_displacements[0]),
        (left.x2, left_displacements[1], right_displacements[1]),
    ):
        weights = _simpson_weights(float(coordinates[-1]), coordinates.size)
        numerator += float(weights @ (mass * np.sum((left_values - right_values) ** 2, axis=1)))
        denominator_left += float(weights @ (mass * np.sum(left_values**2, axis=1)))
        denominator_right += float(weights @ (mass * np.sum(right_values**2, axis=1)))
    return math.sqrt(numerator / max(denominator_left, denominator_right, np.finfo(float).tiny))


def _component_relative_l2(left: PhysicalMode, right: PhysicalMode) -> FloatArray:
    numerator = np.zeros(6, dtype=float)
    denominator_left = np.zeros(6, dtype=float)
    denominator_right = np.zeros(6, dtype=float)
    for coordinates, left_states, right_states in (
        (left.x1, left.state1, right.state1),
        (left.x2, left.state2, right.state2),
    ):
        weights = _simpson_weights(float(coordinates[-1]), coordinates.size)
        numerator += weights @ ((left_states - right_states) ** 2)
        denominator_left += weights @ (left_states**2)
        denominator_right += weights @ (right_states**2)
    return np.sqrt(
        numerator
        / np.maximum(
            np.maximum(denominator_left, denominator_right), np.finfo(float).tiny
        )
    )


def _mode_shape_rows(
    mode: PhysicalMode,
    case: Mapping[str, Any],
    *,
    requested_position: int,
    effective_position: int,
    comparison_semantics: str,
) -> list[dict[str, Any]]:
    beta_deg = float(case["beta_deg"])
    bases = geometry.reddy_inplane_geometry(beta_deg)
    displacements = _global_displacements(mode, beta_deg)
    rows: list[dict[str, Any]] = []
    for arm_index, (coordinates, states, global_values, basis) in enumerate(
        (
            (mode.x1, mode.state1, displacements[0], bases.arm1),
            (mode.x2, mode.state2, displacements[1], bases.arm2),
        ),
        start=1,
    ):
        length = float(coordinates[-1])
        reference_origin = -length * basis.t
        for node, (coordinate, state, displacement) in enumerate(
            zip(coordinates, states, global_values, strict=True)
        ):
            position = reference_origin + float(coordinate) * basis.t
            rows.append(
                {
                    "method": mode.method,
                    "case_id": mode.case_id,
                    "requested_position": requested_position,
                    "effective_position": effective_position,
                    "comparison_semantics": comparison_semantics,
                    "event_id": mode.event_id,
                    "basis_vector": mode.vector_index,
                    "grid_points": coordinates.size,
                    "arm": arm_index,
                    "node": node,
                    "xi": float(coordinate) / length,
                    "x_local": coordinate,
                    "reference_X": position[0],
                    "reference_Y": position[1],
                    "u": state[0],
                    "w": state[1],
                    "psi_physical_EZ": state[2],
                    "N": state[3],
                    "Q": state[4],
                    "M_physical_EZ": state[5],
                    "displacement_X": displacement[0],
                    "displacement_Y": displacement[1],
                    "displacement_Z": displacement[2],
                }
            )
    return rows


def _mass_orthonormal_embedding(
    modes: Sequence[PhysicalMode], mass: float, rotary: float, rank_threshold: float
) -> tuple[FloatArray, int, float]:
    matrix = np.column_stack([_mass_embedding(mode, mass, rotary) for mode in modes])
    gram = 0.5 * (matrix.T @ matrix + (matrix.T @ matrix).T)
    values, vectors = np.linalg.eigh(gram)
    maximum = float(np.max(values)) if values.size else 0.0
    rank = int(np.count_nonzero(values > rank_threshold * maximum)) if maximum > 0.0 else 0
    if rank != matrix.shape[1]:
        raise RuntimeError(
            f"physical cluster field basis lost rank: {rank} != {matrix.shape[1]}"
        )
    inverse_root = vectors @ np.diag(1.0 / np.sqrt(values)) @ vectors.T
    orthonormal = matrix @ inverse_root
    condition = float(maximum / np.min(values))
    return orthonormal, rank, condition


def _evaluate_seed_pair(
    reconstructor: _ModeReconstructor,
    left_seed: ModeSeed,
    right_seed: ModeSeed,
    properties: rlb.BeamProperties,
    points: int,
) -> tuple[PhysicalMode, PhysicalMode, float, float]:
    left_raw = reconstructor.evaluate(left_seed, points)
    right_raw = reconstructor.evaluate(right_seed, points)
    left = _normalize_mode(left_raw, properties.m, properties.J)
    right = _normalize_mode(right_raw, properties.m, properties.J)
    return left, right, _modal_mass(left_raw, properties.m, properties.J), _modal_mass(
        right_raw, properties.m, properties.J
    )


def _isolated_mode_comparison(
    reconstructor: _ModeReconstructor,
    left_block: SpectrumBlock,
    right_block: SpectrumBlock,
    case: Mapping[str, Any],
    contract: Mapping[str, Any],
    *,
    requested_position: int,
    effective_position: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    left_seeds = _seeds_for_block(reconstructor, "RLB", left_block)
    right_seeds = _seeds_for_block(reconstructor, "LEGACY_RECTANGULAR", right_block)
    if len(left_seeds) != 1 or len(right_seeds) != 1:
        raise RuntimeError("isolated mode comparison received a multidimensional block")
    properties = reconstructor.rlb_properties[str(case["geometry"])]
    comparisons: list[dict[str, Any]] = []
    residual_rows: list[dict[str, Any]] = []
    shape_rows: list[dict[str, Any]] = []
    raw_masses: dict[tuple[str, int], float] = {}
    evaluated: dict[tuple[str, int], PhysicalMode] = {}
    for points in (401, 801):
        left, right, left_mass, right_mass = _evaluate_seed_pair(
            reconstructor, left_seeds[0], right_seeds[0], properties, points
        )
        raw_masses[("RLB", points)] = left_mass
        raw_masses[("LEGACY_RECTANGULAR", points)] = right_mass
        left_vector = _mass_embedding(left, properties.m, properties.J)
        right_vector = _mass_embedding(right, properties.m, properties.J)
        inner = float(left_vector @ right_vector)
        alignment = -1.0 if inner < 0.0 else 1.0
        # Only sign-flip the already normalized comparison copy.  Calling the
        # normalization helper here would apply its amplitude factor twice.
        if alignment < 0.0:
            aligned_right = replace(
                right,
                state1=-right.state1,
                state2=-right.state2,
                derivative1=-right.derivative1,
                derivative2=-right.derivative2,
                deterministic_sign=-right.deterministic_sign,
            )
        else:
            aligned_right = right
        mac = min(1.0, max(0.0, inner**2))
        full_difference = float(
            np.linalg.norm(left_vector - alignment * right_vector)
        )
        stiffness_left = _stiffness_embedding(left, properties)
        stiffness_right = _stiffness_embedding(aligned_right, properties)
        resultant_difference = float(
            np.linalg.norm(stiffness_left - stiffness_right)
            / max(np.linalg.norm(stiffness_left), np.linalg.norm(stiffness_right), np.finfo(float).tiny)
        )
        global_difference = _global_displacement_difference(
            left, aligned_right, properties.m, float(case["beta_deg"])
        )
        component_differences = _component_relative_l2(left, aligned_right)
        passed = 1.0 - mac <= float(contract["thresholds"]["isolated_one_minus_MAC"])
        comparisons.append(
            {
                "case_id": left_block.case_id,
                "requested_position": requested_position,
                "effective_position": effective_position,
                "replacement_reason": (
                    "" if requested_position == effective_position else "NEXT_MATCHED_ISOLATED_AFTER_CLUSTER"
                ),
                "grid_points": points,
                "rlb_event_id": left.event_id,
                "legacy_event_id": right.event_id,
                "rlb_omega": left.omega,
                "legacy_omega": right.omega,
                "frequency_relative_difference": _relative_difference(left.omega, right.omega),
                "mass_MAC": mac,
                "one_minus_MAC": 1.0 - mac,
                "comparison_alignment_sign": alignment,
                "full_mass_norm_difference": full_difference,
                "global_displacement_relative_L2": global_difference,
                "resultant_energy_norm_relative_difference": resultant_difference,
                "u_relative_L2": component_differences[0],
                "w_relative_L2": component_differences[1],
                "psi_physical_relative_L2": component_differences[2],
                "N_relative_L2": component_differences[3],
                "Q_relative_L2": component_differences[4],
                "M_physical_relative_L2": component_differences[5],
                "rlb_sign_pivot": left.sign_pivot_index,
                "legacy_sign_pivot": right.sign_pivot_index,
                "status": "PASS" if passed else "FAIL",
            }
        )
        evaluated[("RLB", points)] = left
        evaluated[("LEGACY_RECTANGULAR", points)] = right
    grid_convergence = {
        method: _relative_difference(raw_masses[(method, 401)], raw_masses[(method, 801)])
        for method in ("RLB", "LEGACY_RECTANGULAR")
    }
    sign_consistency = {
        method: bool(
            evaluated[(method, 401)].deterministic_sign
            == evaluated[(method, 801)].deterministic_sign
        )
        for method in ("RLB", "LEGACY_RECTANGULAR")
    }
    for row in comparisons:
        row["rlb_grid_mass_convergence_relative"] = grid_convergence["RLB"]
        row["legacy_grid_mass_convergence_relative"] = grid_convergence[
            "LEGACY_RECTANGULAR"
        ]
        row["rlb_grid_sign_consistent"] = sign_consistency["RLB"]
        row["legacy_grid_sign_consistent"] = sign_consistency[
            "LEGACY_RECTANGULAR"
        ]
        if not (
            row["status"] == "PASS"
            and max(grid_convergence.values())
            <= float(contract["thresholds"]["grid_convergence"])
            and all(sign_consistency.values())
        ):
            row["status"] = "FAIL"
    for method in ("RLB", "LEGACY_RECTANGULAR"):
        for points in (401, 801):
            mode = evaluated[(method, points)]
            residual_rows.append(
                _mode_residual_row(
                    mode,
                    case,
                    properties,
                    contract,
                    requested_position=requested_position,
                    effective_position=effective_position,
                    grid_mass_convergence=grid_convergence[method],
                )
            )
            shape_rows.extend(
                _mode_shape_rows(
                    mode,
                    case,
                    requested_position=requested_position,
                    effective_position=effective_position,
                    comparison_semantics="ISOLATED_POINTWISE",
                )
            )
    return comparisons, residual_rows, shape_rows


def _cluster_mode_comparison(
    reconstructor: _ModeReconstructor,
    left_block: SpectrumBlock,
    right_block: SpectrumBlock,
    case: Mapping[str, Any],
    contract: Mapping[str, Any],
    *,
    requested_position: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    left_seeds = _seeds_for_block(reconstructor, "RLB", left_block)
    right_seeds = _seeds_for_block(reconstructor, "LEGACY_RECTANGULAR", right_block)
    properties = reconstructor.rlb_properties[str(case["geometry"])]
    subspace_rows: list[dict[str, Any]] = []
    residual_rows: list[dict[str, Any]] = []
    shape_rows: list[dict[str, Any]] = []
    metrics: dict[int, float] = {}
    raw_masses: dict[tuple[str, int, int], float] = {}
    all_modes: dict[tuple[str, int], list[PhysicalMode]] = {}
    for points in (401, 801):
        left_raw = [reconstructor.evaluate(seed, points) for seed in left_seeds]
        right_raw = [reconstructor.evaluate(seed, points) for seed in right_seeds]
        for method, modes in (("RLB", left_raw), ("LEGACY_RECTANGULAR", right_raw)):
            all_modes[(method, points)] = [
                _normalize_mode(mode, properties.m, properties.J) for mode in modes
            ]
            for index, mode in enumerate(modes):
                raw_masses[(method, points, index)] = _modal_mass(
                    mode, properties.m, properties.J
                )
        left_basis, left_rank, left_condition = _mass_orthonormal_embedding(
            left_raw,
            properties.m,
            properties.J,
            float(contract["root_policy"]["nullity_relative_threshold"]),
        )
        right_basis, right_rank, right_condition = _mass_orthonormal_embedding(
            right_raw,
            properties.m,
            properties.J,
            float(contract["root_policy"]["nullity_relative_threshold"]),
        )
        if left_rank != right_rank:
            singular = np.array([], dtype=float)
            metric = math.inf
            status = "FAIL"
        else:
            singular = np.linalg.svd(left_basis.T @ right_basis, compute_uv=False)
            singular = np.clip(singular, 0.0, 1.0)
            metric = 1.0 - float(np.min(singular)) ** 2
            status = (
                "PASS"
                if metric <= float(contract["thresholds"]["isolated_one_minus_MAC"])
                else "FAIL"
            )
        metrics[points] = metric
        angles = np.degrees(np.arccos(singular)) if singular.size else np.array([])
        subspace_rows.append(
            {
                "case_id": left_block.case_id,
                "requested_position": requested_position,
                "rlb_block": left_block.block_id,
                "legacy_block": right_block.block_id,
                "grid_points": points,
                "rlb_dimension": left_rank,
                "legacy_dimension": right_rank,
                "rlb_gram_condition": left_condition,
                "legacy_gram_condition": right_condition,
                "overlap_singular_values": singular,
                "principal_angles_deg": angles,
                "minimum_subspace_singular_value": (
                    float(np.min(singular)) if singular.size else math.nan
                ),
                "one_minus_minimum_overlap_squared": metric,
                "threshold_interpretation": "isolated_MAC_threshold_applied_to_1-s_min^2",
                "status": status,
            }
        )
    grid_metric_change = _relative_difference(metrics[401], metrics[801]) if all(
        math.isfinite(metrics[value]) for value in (401, 801)
    ) else math.inf
    for method in ("RLB", "LEGACY_RECTANGULAR"):
        for points in (401, 801):
            for index, mode in enumerate(all_modes[(method, points)]):
                grid_convergence = _relative_difference(
                    raw_masses[(method, 401, index)], raw_masses[(method, 801, index)]
                )
                residual_rows.append(
                    _mode_residual_row(
                        mode,
                        case,
                        properties,
                        contract,
                        requested_position=requested_position,
                        effective_position=requested_position,
                        grid_mass_convergence=grid_convergence,
                    )
                )
                shape_rows.extend(
                    _mode_shape_rows(
                        mode,
                        case,
                        requested_position=requested_position,
                        effective_position=requested_position,
                        comparison_semantics="CLUSTER_BASIS_NO_NATIVE_VECTOR_COMPARISON",
                    )
                )
    for row in subspace_rows:
        row["grid_metric_relative_change_401_801"] = grid_metric_change
    return subspace_rows, residual_rows, shape_rows


def _matched_isolated(
    left: SpectrumBlock | None,
    right: SpectrumBlock | None,
    tolerance: float,
) -> bool:
    return bool(
        left is not None
        and right is not None
        and not left.clustered
        and not right.clustered
        and left.root_count == right.root_count == 1
        and left.nullity == right.nullity == 1
        and _relative_difference(left.centre_omega, right.centre_omega) <= tolerance
    )


def _process_selected_modes(
    contract: Mapping[str, Any],
    blocks_by_method: Mapping[str, Mapping[str, Sequence[SpectrumBlock]]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    reconstructor = _ModeReconstructor(contract)
    cases = {str(case["case_id"]): case for case in contract["cases"]}
    isolated_rows: list[dict[str, Any]] = []
    subspace_rows: list[dict[str, Any]] = []
    residual_rows: list[dict[str, Any]] = []
    shape_rows: list[dict[str, Any]] = []
    processed_clusters: set[tuple[str, str, str]] = set()
    processed_replacements: set[tuple[str, int]] = set()
    frequency_tolerance = float(contract["thresholds"]["isolated_frequency_comparison"])
    for case_id, requested_positions in SELECTED_POSITIONS.items():
        left_blocks = blocks_by_method["RLB"][case_id]
        right_blocks = blocks_by_method["LEGACY_RECTANGULAR"][case_id]
        case = cases[case_id]
        for requested in requested_positions:
            left = _block_for_slot(left_blocks, requested)
            right = _block_for_slot(right_blocks, requested)
            if left is None or right is None:
                isolated_rows.append(
                    {
                        "case_id": case_id,
                        "requested_position": requested,
                        "status": "FAIL",
                        "reason": "REQUESTED_POSITION_ABSENT_FROM_FROZEN_INVENTORY",
                    }
                )
                continue
            clustered = left.clustered or right.clustered
            if clustered:
                cluster_key = (case_id, left.block_id, right.block_id)
                if cluster_key not in processed_clusters:
                    if (
                        left.root_count != right.root_count
                        or left.multiplicity != right.multiplicity
                        or left.nullity != right.nullity
                        or left.clustered != right.clustered
                        or left.cluster_semantics != right.cluster_semantics
                        or _relative_difference(left.centre_omega, right.centre_omega)
                        > float(contract["thresholds"]["cluster_center"])
                    ):
                        subspace_rows.append(
                            {
                                "case_id": case_id,
                                "requested_position": requested,
                                "rlb_block": left.block_id,
                                "legacy_block": right.block_id,
                                "status": "FAIL",
                                "reason": "CLUSTER_SPECTRUM_NOT_MATCHED",
                            }
                        )
                    else:
                        rows, residuals, shapes = _cluster_mode_comparison(
                            reconstructor,
                            left,
                            right,
                            case,
                            contract,
                            requested_position=requested,
                        )
                        subspace_rows.extend(rows)
                        residual_rows.extend(residuals)
                        shape_rows.extend(shapes)
                    processed_clusters.add(cluster_key)
                effective: int | None = None
                start = max(left.slot_end, right.slot_end) + 1
                for candidate in range(start, 13):
                    candidate_left = _block_for_slot(left_blocks, candidate)
                    candidate_right = _block_for_slot(right_blocks, candidate)
                    if _matched_isolated(candidate_left, candidate_right, frequency_tolerance):
                        effective = candidate
                        left = candidate_left
                        right = candidate_right
                        break
                if effective is None:
                    isolated_rows.append(
                        {
                            "case_id": case_id,
                            "requested_position": requested,
                            "status": "NOT_APPLICABLE",
                            "reason": "NO_SUBSEQUENT_MATCHED_ISOLATED_POSITION_BEFORE_GUARD",
                        }
                    )
                    continue
            else:
                effective = requested
            if clustered:
                replacement_key = (case_id, int(effective))
                if replacement_key in processed_replacements:
                    isolated_rows.append(
                        {
                            "case_id": case_id,
                            "requested_position": requested,
                            "effective_position": effective,
                            "status": "NOT_APPLICABLE",
                            "reason": "SHARED_CLUSTER_REPLACEMENT_ALREADY_EMITTED",
                        }
                    )
                    continue
                processed_replacements.add(replacement_key)
            if not _matched_isolated(left, right, frequency_tolerance):
                isolated_rows.append(
                    {
                        "case_id": case_id,
                        "requested_position": requested,
                        "effective_position": effective,
                        "status": "FAIL",
                        "reason": "ISOLATED_SPECTRUM_NOT_MATCHED",
                    }
                )
                continue
            comparisons, residuals, shapes = _isolated_mode_comparison(
                reconstructor,
                left,
                right,
                case,
                contract,
                requested_position=requested,
                effective_position=int(effective),
            )
            isolated_rows.extend(comparisons)
            residual_rows.extend(residuals)
            shape_rows.extend(shapes)
    if not subspace_rows:
        subspace_rows.append(
            {
                "status": "NOT_APPLICABLE",
                "reason": "NO_PRESELECTED_POSITION_BELONGED_TO_A_CLUSTER",
            }
        )
    return isolated_rows, subspace_rows, residual_rows, shape_rows


def _inventory_blocks(
    output_directory: Path,
    rlb_manifest: Mapping[str, Any],
    legacy_manifest: Mapping[str, Any],
) -> tuple[
    dict[str, dict[str, list[SpectrumBlock]]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    rlb_rows = _read_csv(output_directory / "rlb_roots.csv")
    legacy_rows = _read_csv(output_directory / "legacy_rectangular_roots.csv")
    if _semantic_inventory_hash(rlb_rows, rlb_manifest) != rlb_manifest.get(
        "semantic_inventory_sha256_before_comparison"
    ):
        raise RuntimeError("RLB semantic inventory hash mismatch")
    if _semantic_inventory_hash(legacy_rows, legacy_manifest) != legacy_manifest.get(
        "semantic_inventory_sha256_before_comparison"
    ):
        raise RuntimeError("legacy semantic inventory hash mismatch")
    rlb_events = _parse_root_events(rlb_rows)
    legacy_events = _parse_root_events(legacy_rows)
    missing_rlb = [case_id for case_id in REQUIRED_CASE_IDS if not rlb_events.get(case_id)]
    missing_legacy = [case_id for case_id in REQUIRED_CASE_IDS if not legacy_events.get(case_id)]
    if missing_rlb or missing_legacy:
        raise RuntimeError(
            f"incomplete frozen inventories: RLB={missing_rlb}, legacy={missing_legacy}"
        )
    blocks = {
        "RLB": {
            case_id: _spectrum_blocks(rlb_events[case_id]) for case_id in REQUIRED_CASE_IDS
        },
        "LEGACY_RECTANGULAR": {
            case_id: _spectrum_blocks(legacy_events[case_id])
            for case_id in REQUIRED_CASE_IDS
        },
    }
    return blocks, rlb_rows, legacy_rows


def _cross_method_spectrum_rows(
    contract: Mapping[str, Any],
    blocks: Mapping[str, Mapping[str, Sequence[SpectrumBlock]]],
) -> list[dict[str, Any]]:
    threshold = float(contract["thresholds"]["isolated_frequency_comparison"])
    cluster_threshold = float(contract["thresholds"]["cluster_center"])
    rows: list[dict[str, Any]] = []
    for case_id in REQUIRED_CASE_IDS:
        rows.extend(
            _compare_block_lists(
                "RLB_VS_LEGACY_RECTANGULAR",
                case_id,
                "RLB_FOUR_PLY_TRANSFER",
                "LEGACY_RECTANGULAR_CLOSED_FORM",
                blocks["RLB"][case_id],
                blocks["LEGACY_RECTANGULAR"][case_id],
                isolated_tolerance=threshold,
                cluster_tolerance=cluster_threshold,
            )
        )
    return rows


def _beta0_direct_rows(
    contract: Mapping[str, Any],
    blocks: Mapping[str, Mapping[str, Sequence[SpectrumBlock]]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    cases = {str(case["case_id"]): case for case in contract["cases"]}
    isolated_threshold = float(contract["thresholds"]["isolated_frequency_comparison"])
    cluster_threshold = float(contract["thresholds"]["cluster_center"])
    cache: dict[tuple[str, float], tuple[list[SpectrumBlock], list[SpectrumBlock], dict[str, Any]]] = {}
    rows: list[dict[str, Any]] = []
    details: dict[str, Any] = {}
    for case_id in ("ISO-01", "ISO-02"):
        case = cases[case_id]
        if float(case["beta_deg"]) != 0.0:
            raise RuntimeError(f"{case_id} is not the frozen beta=0 direct-control case")
        if float(case["L1"]) <= 0.0 or float(case["L2"]) <= 0.0:
            raise RuntimeError(f"{case_id} has a non-positive arm length")
        if not math.isclose(
            float(case["L1"]) + float(case["L2"]),
            float(contract["lengths"]["L_total"]),
            rel_tol=0.0,
            abs_tol=8.0 * np.finfo(float).eps,
        ):
            raise RuntimeError(f"{case_id} does not preserve the frozen total length")
        if "joint_mass_or_compliance" not in set(contract["explicit_exclusions"]):
            raise RuntimeError("direct beta=0 control requires the frozen ideal-joint exclusion")
        geometry_id = str(case["geometry"])
        total_length = float(case["L1"]) + float(case["L2"])
        key = (geometry_id, total_length)
        if key not in cache:
            section = _rectangular_section(contract, geometry_id)
            scale = _frequency_scale(contract, geometry_id)
            direct_blocks, direct_audit = _verified_small_spectrum(
                _direct_fixed_fixed_provider(section, total_length),
                scale,
                contract,
                case_id="DIRECT-FIXED-FIXED",
            )
            bending_blocks, bending_audit = _verified_small_spectrum(
                _direct_bending_provider(section, total_length),
                scale,
                contract,
                case_id="DIRECT-BENDING-FAMILY",
            )
            union_blocks = _union_family_blocks(
                "DIRECT-AXIAL-BENDING-UNION",
                section,
                total_length,
                scale,
                bending_blocks,
                contract,
            )
            cache[key] = (
                direct_blocks,
                union_blocks,
                {
                    "direct_primary_verification": direct_audit,
                    "bending_primary_verification": bending_audit,
                },
            )
        direct_blocks, union_blocks, audit = cache[key]
        details[case_id] = audit
        for audit_name in ("direct_primary_verification", "bending_primary_verification"):
            for audit_row in audit[audit_name]["primary_verification_rows"]:
                copied = dict(audit_row)
                copied["case_id"] = case_id
                copied["comparison_kind"] = f"BETA0_{audit_name.upper()}"
                rows.append(copied)
        for comparison_kind, left_label, left_blocks, right_label, right_blocks in (
            (
                "BETA0_RLB_VS_DIRECT_FIXED_FIXED",
                "RLB_FOUR_PLY_TRANSFER",
                blocks["RLB"][case_id],
                "DIRECT_RECTANGULAR_FIXED_FIXED",
                direct_blocks,
            ),
            (
                "BETA0_LEGACY_COUPLED_VS_DIRECT_FIXED_FIXED",
                "LEGACY_RECTANGULAR_CLOSED_FORM",
                blocks["LEGACY_RECTANGULAR"][case_id],
                "DIRECT_RECTANGULAR_FIXED_FIXED",
                direct_blocks,
            ),
            (
                "BETA0_DIRECT_FIXED_FIXED_VS_AXIAL_BENDING_UNION",
                "DIRECT_RECTANGULAR_FIXED_FIXED",
                direct_blocks,
                "INDEPENDENT_AXIAL_BENDING_UNION",
                union_blocks,
            ),
            (
                "BETA0_RLB_VS_AXIAL_BENDING_UNION",
                "RLB_FOUR_PLY_TRANSFER",
                blocks["RLB"][case_id],
                "INDEPENDENT_AXIAL_BENDING_UNION",
                union_blocks,
            ),
        ):
            rows.extend(
                _compare_block_lists(
                    comparison_kind,
                    case_id,
                    left_label,
                    right_label,
                    left_blocks,
                    right_blocks,
                    isolated_tolerance=isolated_threshold,
                    cluster_tolerance=cluster_threshold,
                )
            )
    return rows, details


def _symmetry_rows(
    contract: Mapping[str, Any],
    blocks: Mapping[str, Mapping[str, Sequence[SpectrumBlock]]],
) -> list[dict[str, Any]]:
    tolerance = float(contract["thresholds"]["symmetry_spectrum"])
    rows: list[dict[str, Any]] = []
    for method in ("RLB", "LEGACY_RECTANGULAR"):
        for relation, left_case, right_case in (
            ("ANGLE_REFLECTION_PLUS30_MINUS30", "ISO-03", "ISO-05"),
            ("ARM_EXCHANGE_0P7_1P3", "ISO-06", "ISO-07"),
        ):
            relation_rows = _compare_block_lists(
                relation,
                f"{left_case}_vs_{right_case}",
                f"{method}:{left_case}",
                f"{method}:{right_case}",
                blocks[method][left_case],
                blocks[method][right_case],
                isolated_tolerance=tolerance,
                cluster_tolerance=tolerance,
            )
            for row in relation_rows:
                row["method"] = method
                row["left_case"] = left_case
                row["right_case"] = right_case
            rows.extend(relation_rows)
    return rows


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )


def _direct_numeric_projection_sha256(
    rows: Sequence[Mapping[str, Any]],
) -> str:
    """Hash the ordered numerical payload of the superseded direct rows."""

    payload = [
        {
            "comparison_kind": str(row.get("comparison_kind", "")),
            "case_id": str(row.get("case_id", "")),
            **{field: str(row.get(field, "")) for field in DIRECT_NUMERIC_FIELDS},
        }
        for row in rows
        if str(row.get("comparison_kind", "")) in DIRECT_BETA0_AUXILIARY_KINDS
    ]
    encoded = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return _sha256_bytes(encoded)


def classify_existing_spectrum_rows(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Separate frozen scientific evidence from superseded beta=0 diagnostics.

    The function changes no numerical field.  The former status and reason of
    every auxiliary row are retained explicitly before the closing
    classification is applied.
    """

    classified: list[dict[str, Any]] = []
    counts = {kind: 0 for kind in DIRECT_BETA0_AUXILIARY_KINDS}
    scientific_count = 0
    for source in rows:
        row = dict(source)
        kind = str(row.get("comparison_kind", ""))
        if kind == SCIENTIFIC_COMPARISON_KIND:
            scientific_count += 1
            row["scientific_role"] = "SCIENTIFIC_EVIDENCE"
            row["used_for_scientific_status"] = "true"
            row.setdefault("original_status", "")
            row.setdefault("original_reason", "")
        elif kind in DIRECT_BETA0_AUXILIARY_KINDS:
            counts[kind] += 1
            row["scientific_role"] = "AUXILIARY_DIAGNOSTIC"
            row["original_status"] = row.get("original_status") or row.get("status", "")
            row["original_reason"] = row.get("original_reason") or row.get("reason", "")
            row["status"] = "SUPERSEDED_DIAGNOSTIC_ARTIFACT"
            row["used_for_scientific_status"] = "false"
            row["reason"] = DIRECT_BETA0_QUALIFICATION_REASON
        else:
            raise RuntimeError(f"unexpected spectrum comparison kind: {kind!r}")
        classified.append(row)

    if scientific_count != 104:
        raise RuntimeError(
            f"expected 104 scientific spectrum rows, got {scientific_count}"
        )
    if sum(counts.values()) != 164 or any(value <= 0 for value in counts.values()):
        raise RuntimeError(f"unexpected direct-beta0 auxiliary inventory: {counts}")
    return classified, {
        "scientific_row_count": scientific_count,
        "auxiliary_row_count": sum(counts.values()),
        "auxiliary_counts_by_kind": counts,
    }


def _validate_primary_spectrum_rows(
    rows: Sequence[Mapping[str, Any]], contract: Mapping[str, Any]
) -> dict[str, Any]:
    primary = [
        row
        for row in rows
        if str(row.get("comparison_kind", "")) == SCIENTIFIC_COMPARISON_KIND
    ]
    if len(primary) != 104:
        raise RuntimeError(f"expected 104 cross-model comparisons, got {len(primary)}")
    threshold = float(contract["thresholds"]["isolated_frequency_comparison"])
    by_case: dict[str, list[Mapping[str, Any]]] = {
        case_id: [] for case_id in REQUIRED_CASE_IDS
    }
    for row in primary:
        case_id = str(row.get("case_id", ""))
        if case_id not in by_case:
            raise RuntimeError(f"unknown scientific comparison case: {case_id}")
        by_case[case_id].append(row)
        if str(row.get("status", "")) != "PASS":
            raise RuntimeError(f"scientific spectrum row is not PASS: {case_id}")
        if float(row["relative_difference"]) > threshold:
            raise RuntimeError(f"scientific spectrum tolerance failed: {case_id}")
    for case_id, case_rows in by_case.items():
        slots = sorted(int(row["alignment_slot"]) for row in case_rows)
        if slots != list(range(1, 14)):
            raise RuntimeError(f"{case_id} does not contain slots 1--13: {slots}")
    worst = max(primary, key=lambda row: float(row["relative_difference"]))
    return {
        "comparison_count": len(primary),
        "case_count": len(by_case),
        "slots_per_case": 13,
        "requested_roots_per_case": 12,
        "guard_roots_per_case": 1,
        "maximum_relative_difference": float(worst["relative_difference"]),
        "maximum_case_id": str(worst["case_id"]),
        "maximum_alignment_slot": int(worst["alignment_slot"]),
        "threshold": threshold,
        "status": "PASS",
    }


def _validate_preliminary_evidence(output_directory: Path) -> dict[str, Any]:
    definitions = {
        "isotropic_qbar_checks.csv": (
            "Qbar_relative_residual",
            "shear_Qbar_relative_residual",
        ),
        "four_ply_laminate_checks.csv": ("maximum_residual",),
        "section_reduction_checks.csv": ("maximum_residual",),
        "legacy_circular_backcompat.csv": ("maximum_residual", "root_relative_residual"),
        "local_solution_space_checks.csv": ("maximum_residual",),
    }
    result: dict[str, Any] = {}
    for name, metric_fields in definitions.items():
        rows = _read_csv(output_directory / name)
        if not rows or any(str(row.get("status", "")) != "PASS" for row in rows):
            raise RuntimeError(f"preliminary evidence is not entirely PASS: {name}")
        values = [
            float(row[field])
            for row in rows
            for field in metric_fields
            if str(row.get(field, "")).strip() not in {"", "nan", "NaN"}
        ]
        result[name] = {
            "row_count": len(rows),
            "maximum_recorded_residual": max(values, default=0.0),
            "status": "PASS",
        }
    return result


def _validate_mode_evidence(
    output_directory: Path, contract: Mapping[str, Any]
) -> dict[str, Any]:
    comparisons = _read_csv(output_directory / "mode_shape_comparison.csv")
    residuals = _read_csv(output_directory / "mode_joint_residuals.csv")
    subspaces = _read_csv(output_directory / "mode_subspace_comparison.csv")
    if len(comparisons) != 14 or any(row.get("status") != "PASS" for row in comparisons):
        raise RuntimeError("expected 14 PASS pointwise mode comparisons")
    if len(residuals) != 28 or any(row.get("status") != "PASS" for row in residuals):
        raise RuntimeError("expected 28 PASS physical mode-residual rows")
    if any(
        row.get("status") not in {"PASS", "NOT_APPLICABLE"} for row in subspaces
    ):
        raise RuntimeError("mode subspace evidence contains a failed row")

    metrics = {
        "maximum_one_minus_MAC": max(float(row["one_minus_MAC"]) for row in comparisons),
        "maximum_outer_clamp_normalized": max(
            float(row["outer_clamp_normalized"]) for row in residuals
        ),
        "maximum_joint_compatibility_normalized": max(
            float(row["joint_compatibility_normalized"]) for row in residuals
        ),
        "maximum_joint_equilibrium_normalized": max(
            float(row["joint_equilibrium_normalized"]) for row in residuals
        ),
        "maximum_joint_force_normalized": max(
            float(row["joint_force_normalized"]) for row in residuals
        ),
        "maximum_joint_moment_normalized": max(
            float(row["joint_moment_normalized"]) for row in residuals
        ),
        "maximum_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in residuals
        ),
        "maximum_grid_mass_convergence_relative": max(
            float(row["grid_mass_convergence_relative"]) for row in residuals
        ),
        "maximum_energy_identity_relative": max(
            float(row["energy_identity_relative"]) for row in residuals
        ),
        "maximum_differential_normalized": max(
            float(row["maximum_differential_normalized"]) for row in residuals
        ),
    }
    thresholds = contract["thresholds"]
    required = {
        "maximum_one_minus_MAC": float(thresholds["isolated_one_minus_MAC"]),
        "maximum_outer_clamp_normalized": float(thresholds["outer_clamp_residual"]),
        "maximum_joint_compatibility_normalized": float(
            thresholds["joint_compatibility"]
        ),
        "maximum_joint_equilibrium_normalized": float(
            thresholds["joint_equilibrium"]
        ),
        "maximum_joint_force_normalized": float(thresholds["joint_equilibrium"]),
        "maximum_joint_moment_normalized": float(thresholds["joint_equilibrium"]),
        "maximum_boundary_null_residual": float(
            thresholds["boundary_null_residual"]
        ),
        "maximum_grid_mass_convergence_relative": float(
            thresholds["grid_convergence"]
        ),
        "maximum_energy_identity_relative": float(thresholds["energy_identity"]),
    }
    failed = {
        key: (metrics[key], limit)
        for key, limit in required.items()
        if metrics[key] > limit
    }
    if failed:
        raise RuntimeError(f"mode evidence gate failed: {failed}")
    return {
        "pointwise_comparison_rows": len(comparisons),
        "physical_residual_rows": len(residuals),
        "subspace_rows": len(subspaces),
        **metrics,
        "differential_residual_role": "NON_GATING_DIAGNOSTIC_NO_FROZEN_THRESHOLD",
        "status": "PASS",
    }


def _verify_frozen_reference_hashes(
    contract_file: Path, output_directory: Path
) -> dict[str, Any]:
    manifest = json.loads(
        (output_directory / "model_manifest.json").read_text(encoding="utf-8")
    )
    repository_root = contract_file.parents[2]
    verified: dict[str, str] = {}
    for record in manifest.get("frozen_reference_hashes", []):
        relative = str(record["path"])
        path = repository_root / relative
        expected = str(record["sha256"])
        if not path.is_file() or _sha256_file(path) != expected:
            raise RuntimeError(f"frozen reference changed: {relative}")
        verified[relative] = expected
    if len(verified) != 29:
        raise RuntimeError(f"expected 29 frozen references, got {len(verified)}")
    return {"verified_count": len(verified), "sha256": verified, "status": "PASS"}


def _verify_iso04_close_pair(output_directory: Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for method, name in (
        ("RLB", "rlb_root_candidates.csv"),
        ("LEGACY_RECTANGULAR", "legacy_rectangular_root_candidates.csv"),
    ):
        rows = _read_csv(output_directory / name)
        method_result: dict[str, Any] = {}
        for scan_id in ("primary", "verification"):
            values = sorted(
                float(row["Omega"])
                for row in rows
                if row.get("case_id") == "ISO-04"
                and row.get("scan_id") == scan_id
                and row.get("accepted") == "True"
                and row.get("canonical") == "True"
                and 309.0 < float(row["Omega"]) < 309.1
            )
            if len(values) != 2 or values[1] <= values[0]:
                raise RuntimeError(f"ISO-04 close pair was not preserved in {name}/{scan_id}")
            method_result[scan_id] = {
                "Omega": values,
                "separation_Omega": values[1] - values[0],
            }
        result[method] = method_result
    result["inventory_note"] = (
        "FIRST_ROOT_IS_SLOT_13_GUARD_SECOND_IS_ACCEPTED_CANONICAL_POST_GUARD_TAIL"
    )
    result["status"] = "PASS"
    return result


def _frequency_scale(contract: Mapping[str, Any], geometry_id: str) -> float:
    material = contract["material"]
    geometry = contract["geometries"][geometry_id]
    width = float(geometry["width"])
    thickness = float(geometry["thickness"])
    area = width * thickness
    inertia = width * thickness**3 / 12.0
    length_ref = float(contract["lengths"]["L_ref"])
    return length_ref**2 * math.sqrt(
        float(material["rho"]) * area / (float(material["E"]) * inertia)
    )


def _validate_saved_symmetry_evidence(
    output_directory: Path, contract: Mapping[str, Any]
) -> dict[str, Any]:
    rows = _read_csv(output_directory / "symmetry_checks.csv")
    expected_groups = {
        ("RLB", "ANGLE_REFLECTION_PLUS30_MINUS30"),
        ("RLB", "ARM_EXCHANGE_0P7_1P3"),
        ("LEGACY_RECTANGULAR", "ANGLE_REFLECTION_PLUS30_MINUS30"),
        ("LEGACY_RECTANGULAR", "ARM_EXCHANGE_0P7_1P3"),
    }
    groups: dict[tuple[str, str], list[Mapping[str, Any]]] = {
        key: [] for key in expected_groups
    }
    for row in rows:
        key = (str(row.get("method", "")), str(row.get("comparison_kind", "")))
        if key not in groups:
            raise RuntimeError(f"unexpected saved symmetry group: {key}")
        groups[key].append(row)
    if len(rows) != 52:
        raise RuntimeError(f"expected 52 saved symmetry rows, got {len(rows)}")
    for key, group_rows in groups.items():
        slots = sorted(int(row["alignment_slot"]) for row in group_rows)
        if slots != list(range(1, 14)):
            raise RuntimeError(f"saved symmetry slots differ for {key}: {slots}")

    expected_exception = ("LEGACY_RECTANGULAR", "ARM_EXCHANGE_0P7_1P3")
    unexpected_failures = [
        row
        for key, group_rows in groups.items()
        if key != expected_exception
        for row in group_rows
        if row.get("status") != "PASS"
    ]
    if unexpected_failures:
        raise RuntimeError("saved symmetry evidence has an unexpected failure")
    legacy_arm_rows = groups[expected_exception]
    legacy_arm_failures = [row for row in legacy_arm_rows if row.get("status") != "PASS"]
    if len(legacy_arm_failures) != 1 or int(
        legacy_arm_failures[0]["alignment_slot"]
    ) != 12:
        raise RuntimeError("saved legacy arm exchange differs from the sole slot-12 exception")
    if any(
        row.get("status") != "PASS" and row is not legacy_arm_failures[0]
        for row in legacy_arm_rows
    ):
        raise RuntimeError("saved legacy arm exchange contains another failure")
    tolerance = float(contract["thresholds"]["symmetry_spectrum"])
    if any(float(row["tolerance"]) != tolerance for row in rows):
        raise RuntimeError("saved symmetry threshold differs from the frozen contract")
    return {
        "row_count": len(rows),
        "group_count": len(groups),
        "rows_per_group": 13,
        "pre_refinement_pass_rows": sum(row.get("status") == "PASS" for row in rows),
        "sole_pre_refinement_exception": {
            "method": "LEGACY_RECTANGULAR",
            "comparison_kind": "ARM_EXCHANGE_0P7_1P3",
            "alignment_slot": 12,
            "relative_difference": float(
                legacy_arm_failures[0]["relative_difference"]
            ),
        },
        "threshold": tolerance,
        "status": "PASS_WITH_DECLARED_BOUNDED_REFINEMENT",
    }


def bounded_sign_change_refinement(
    *,
    old_Omega: float,
    bracket_left_Omega: float,
    bracket_right_Omega: float,
    frequency_scale: float,
    determinant_Omega: Callable[[float], float],
    diagnostics_Omega: Callable[[float], Mapping[str, float]],
) -> dict[str, Any]:
    """Refine one root only inside its own saved sign-changing bracket."""

    left = float(bracket_left_Omega)
    right = float(bracket_right_Omega)
    value_left = float(determinant_Omega(left))
    value_right = float(determinant_Omega(right))
    if not (
        left < float(old_Omega) < right
        and math.isfinite(value_left)
        and math.isfinite(value_right)
        and value_left * value_right < 0.0
    ):
        raise ValueError("frozen bracket is not a finite sign-changing bracket")
    refined = float(
        brentq(
            determinant_Omega,
            left,
            right,
            xtol=float(np.nextafter(0.0, 1.0)),
            rtol=4.0 * np.finfo(float).eps,
            maxiter=256,
        )
    )
    diagnostic = dict(diagnostics_Omega(refined))
    return {
        "old_Omega": float(old_Omega),
        "refined_Omega": refined,
        "old_omega": float(old_Omega) / float(frequency_scale),
        "refined_omega": refined / float(frequency_scale),
        "bracket_left_Omega": left,
        "bracket_right_Omega": right,
        "bracket_left_scaled_determinant": value_left,
        "bracket_right_scaled_determinant": value_right,
        **diagnostic,
    }


def _bounded_legacy_arm_exchange_refinement(
    contract: Mapping[str, Any], output_directory: Path
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    symmetry_rows = _read_csv(output_directory / "symmetry_checks.csv")
    relevant = [
        row
        for row in symmetry_rows
        if row.get("comparison_kind") == "ARM_EXCHANGE_0P7_1P3"
        and row.get("method") == "LEGACY_RECTANGULAR"
        and row.get("left_case") == "ISO-06"
        and row.get("right_case") == "ISO-07"
    ]
    if not relevant:
        raise RuntimeError("legacy ISO-06/ISO-07 symmetry evidence is missing")
    offender = max(relevant, key=lambda row: float(row["relative_difference"]))
    slot = int(offender["alignment_slot"])
    if slot != 12:
        raise RuntimeError(f"unexpected legacy arm-exchange offender slot: {slot}")

    root_rows = _read_csv(output_directory / "legacy_rectangular_roots.csv")
    roots = {
        row["case_id"]: row
        for row in root_rows
        if row.get("case_id") in {"ISO-06", "ISO-07"}
        and int(row["sorted_slot"]) == slot
    }
    if set(roots) != {"ISO-06", "ISO-07"}:
        raise RuntimeError("frozen legacy roots do not contain both offender rows")

    comparator = importlib.import_module(
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams"
    )
    case_by_id = {str(case["case_id"]): case for case in contract["cases"]}
    refined_by_case: dict[str, dict[str, Any]] = {}
    errors: dict[str, str] = {}
    for case_id in ("ISO-06", "ISO-07"):
        case = case_by_id[case_id]
        geometry_id = str(case["geometry"])
        geometry_values = contract["geometries"][geometry_id]
        material = contract["material"]
        section = comparator.rectangular_section(
            E=float(material["E"]),
            nu=float(material["nu"]),
            rho=float(material["rho"]),
            width=float(geometry_values["width"]),
            thickness=float(geometry_values["thickness"]),
            K=float(material["K"]),
        )
        length_1 = float(case["L1"])
        length_2 = float(case["L2"])
        beta_deg = float(case["beta_deg"])
        scale = _frequency_scale(contract, geometry_id)

        def determinant_Omega(Omega: float) -> float:
            return float(
                comparator.legacy_coupled_characteristic_determinant(
                    float(Omega) / scale,
                    section,
                    length_1,
                    section,
                    length_2,
                    beta_deg=beta_deg,
                    scaled=True,
                )
            )

        def diagnostics_Omega(Omega: float) -> Mapping[str, float]:
            diagnostic = comparator.evaluate_legacy_coupled_boundary(
                float(Omega) / scale,
                section,
                length_1,
                section,
                length_2,
                beta_deg=beta_deg,
            ).diagnostics
            return {
                "raw_determinant": float(diagnostic.raw_determinant),
                "raw_determinant_residual": abs(float(diagnostic.raw_determinant)),
                "scaled_determinant": float(diagnostic.scaled_determinant),
                "scaled_determinant_residual": abs(
                    float(diagnostic.scaled_determinant)
                ),
                "raw_sigma_ratio": float(diagnostic.raw_sigma_ratio),
                "scaled_sigma_ratio": float(diagnostic.scaled_sigma_ratio),
                "raw_boundary_null_residual": float(
                    diagnostic.raw_boundary_null_residual
                ),
                "scaled_boundary_null_residual": float(
                    diagnostic.scaled_null_residual
                ),
            }

        frozen = roots[case_id]
        try:
            refined_by_case[case_id] = bounded_sign_change_refinement(
                old_Omega=float(frozen["Omega"]),
                bracket_left_Omega=float(frozen["bracket_left_Omega"]),
                bracket_right_Omega=float(frozen["bracket_right_Omega"]),
                frequency_scale=scale,
                determinant_Omega=determinant_Omega,
                diagnostics_Omega=diagnostics_Omega,
            )
        except (ArithmeticError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            errors[case_id] = f"{type(error).__name__}: {error}"

    old_difference = float(offender["relative_difference"])
    threshold = float(contract["thresholds"]["symmetry_spectrum"])
    if not errors and len(refined_by_case) == 2:
        refined_difference = _relative_difference(
            refined_by_case["ISO-06"]["refined_omega"],
            refined_by_case["ISO-07"]["refined_omega"],
        )
        status = (
            ARM_EXCHANGE_PASS_STATUS
            if refined_difference <= threshold
            else ARM_EXCHANGE_PARTIAL_STATUS
        )
    else:
        refined_difference = math.nan
        status = ARM_EXCHANGE_PARTIAL_STATUS

    rows: list[dict[str, Any]] = []
    for case_id in ("ISO-06", "ISO-07"):
        frozen = roots[case_id]
        detail = refined_by_case.get(case_id, {})
        rows.append(
            {
                "comparison_kind": "LEGACY_ARM_EXCHANGE_BOUNDED_LOCAL_REFINEMENT",
                "case_id": case_id,
                "paired_case_id": "ISO-07" if case_id == "ISO-06" else "ISO-06",
                "sorted_slot": slot,
                "old_Omega": float(frozen["Omega"]),
                "old_omega": float(frozen["omega"]),
                "refined_Omega": detail.get("refined_Omega", ""),
                "refined_omega": detail.get("refined_omega", ""),
                "bracket_left_Omega": float(frozen["bracket_left_Omega"]),
                "bracket_right_Omega": float(frozen["bracket_right_Omega"]),
                "bracket_left_scaled_determinant": detail.get(
                    "bracket_left_scaled_determinant", ""
                ),
                "bracket_right_scaled_determinant": detail.get(
                    "bracket_right_scaled_determinant", ""
                ),
                "raw_determinant": detail.get("raw_determinant", ""),
                "raw_determinant_residual": detail.get(
                    "raw_determinant_residual", ""
                ),
                "scaled_determinant": detail.get("scaled_determinant", ""),
                "scaled_determinant_residual": detail.get(
                    "scaled_determinant_residual", ""
                ),
                "raw_sigma_ratio": detail.get("raw_sigma_ratio", ""),
                "scaled_sigma_ratio": detail.get("scaled_sigma_ratio", ""),
                "raw_boundary_null_residual": detail.get(
                    "raw_boundary_null_residual", ""
                ),
                "scaled_boundary_null_residual": detail.get(
                    "scaled_boundary_null_residual", ""
                ),
                "old_symmetry_relative_difference": old_difference,
                "refined_symmetry_relative_difference": refined_difference,
                "threshold": threshold,
                "own_frozen_bracket_used": True,
                "paired_case_root_used_as_seed": False,
                "global_scan_run": False,
                "refinement_error": errors.get(case_id, ""),
                "scientific_role": "AUXILIARY_DIAGNOSTIC",
                "used_for_scientific_status": False,
                "status": status,
            }
        )
    return rows, {
        "offender_slot": slot,
        "old_relative_difference": old_difference,
        "refined_relative_difference": refined_difference,
        "threshold": threshold,
        "refined_case_count": len(refined_by_case),
        "global_scan_run": False,
        "paired_case_root_used_as_seed": False,
        "errors": errors,
        "status": status,
    }


def _render_closing_report(summary: Mapping[str, Any]) -> str:
    spectrum = summary["spectrum_evidence"]
    modes = summary["mode_evidence"]
    refinement = summary["arm_exchange_refinement"]
    direct = summary["direct_beta0_auxiliary"]
    statuses = summary["scientific_statuses"]
    status_lines = "\n".join(f"- `{name}: {value}`" for name, value in statuses.items())
    auxiliary_kinds = "\n".join(
        f"- `{kind}`: {count} rows"
        for kind, count in direct["counts_by_kind"].items()
    )
    return f"""# RLB-1C-ISO: закрывающий отчёт

## Результат

Для восьми заранее объявленных случаев подтверждён четырёхслойный
изотропный предел Reddy-beam относительно независимого прямоугольного
closed-form determinant Тимошенко. Численный contract использует ровно четыре
равных слоя `[0/90/90/0]` и общий для двух моделей коэффициент `K=5/6`.

Разрешённая итоговая формулировка:

`{SCIENTIFIC_RESULT_STATEMENT}`

## Аналитический предел

Из изотропных свойств каждого слоя и симметричного четырёхслойного
интегрирования следуют

`EA = E*b*h`, `EI = E*b*h^3/12`, `KGA = (5/6)*G*b*h`,
`m = rho*b*h`, `J = rho*b*h^3/12`.

Проверки угловой инвариантности `Qbar` и transverse shear, стеков
`[0/90/90/0]`, `[0/0/0/0]`, `[17/-38/-38/17]`, а также условий `B=0` и
`I1=0` имеют статус `PASS`. Независимый comparator не импортирует RLB
transfer-модули и сохраняет прежний circular benchmark в below-cutoff,
cutoff и above-cutoff regimes. Максимальный residual локального пространства
clamped-arm states равен `{summary['preliminary_evidence']['local_solution_space_checks.csv']['maximum_recorded_residual']:.16e}`.

## Frozen coupled inventories

- RLB semantic hash: `{summary['rlb_inventory_sha256']}`;
- legacy rectangular semantic hash: `{summary['legacy_inventory_sha256']}`;
- cases: `{spectrum['case_count']}`;
- roots каждого метода: `{spectrum['comparison_count']}` = 8 x 13;
- first roots: 12 на случай; root 13 сохранён как guard;
- межмодельные сопоставления: `{spectrum['comparison_count']}/{spectrum['comparison_count']} PASS`;
- максимальная относительная разность: `{spectrum['maximum_relative_difference']:.16e}`
  (`{spectrum['maximum_case_id']}`, slot {spectrum['maximum_alignment_slot']},
  threshold `{spectrum['threshold']:.1e}`).

Близкие кандидаты ISO-04 около `Omega=309.0489` и `309.0784` остаются двумя
различными accepted canonical roots. Первый является slot-13 guard, второй
сохранён в post-guard candidate tail. Ни один inventory worker не читал roots
другого метода и не использовал их как seeds.

## Физические формы

Все 14 сохранённых pointwise comparisons и 28 строк физических residuals
имеют статус `PASS`. Максимумы: `1-MAC = {modes['maximum_one_minus_MAC']:.16e}`,
outer clamp `{modes['maximum_outer_clamp_normalized']:.16e}`, joint
compatibility `{modes['maximum_joint_compatibility_normalized']:.16e}`, force
equilibrium `{modes['maximum_joint_force_normalized']:.16e}`, moment
equilibrium `{modes['maximum_joint_moment_normalized']:.16e}`, energy identity
`{modes['maximum_energy_identity_relative']:.16e}` и 401/801 grid convergence
`{modes['maximum_grid_mass_convergence_relative']:.16e}`.

## Вспомогательные численные квалификации

Direct beta=0 fixed--fixed 3x3 representation не входит в научный gate. Его
164 старые строки сохранены с прежними числовыми полями и классифицированы как
`AUXILIARY_DIAGNOSTIC / SUPERSEDED_DIAGNOSTIC_ARTIFACT`; причина —
`{DIRECT_BETA0_QUALIFICATION_REASON}`. Состав:

{auxiliary_kinds}

Статус: `AUX-LEGACY-DIRECT-BETA0-FIXED-FIXED:
{DIRECT_BETA0_AUXILIARY_STATUS}`.

Единственное превышение auxiliary threshold для legacy arm exchange относилось
к ISO-06/ISO-07, slot 12. Каждый root уточнён независимо только внутри своего
frozen sign-changing bracket; зеркальный root не служил seed. Разность
уменьшилась с `{refinement['old_relative_difference']:.16e}` до
`{refinement['refined_relative_difference']:.16e}` при неизменном пороге
`{refinement['threshold']:.1e}`. Статус:
`AUX-LEGACY-ARM-EXCHANGE: {refinement['status']}`. Эта квалификация описывает
точность локализации root, а не установленную ошибку модели.

## Статусы

{status_lines}

- `AUX-LEGACY-DIRECT-BETA0-FIXED-FIXED: {DIRECT_BETA0_AUXILIARY_STATUS}`
- `AUX-LEGACY-ARM-EXCHANGE: {refinement['status']}`
- `ORIGINAL-CONTRACT-OVERALL: PARTIAL_PASS`
- `SCIENTIFIC_OVERALL: {SCIENTIFIC_OVERALL_STATUS}`

## Ограничения

Проверен только объявленный конечный набор ISO-01--ISO-08, а не произвольный
угол, геометрия или материал. Исторический Ritz FAIL предыдущего этапа не
входит в этот scientific overall. В closing stage не выполнялись global root
search, Ritz, Euler--Bernoulli, FEM, beta sweep, parameter study, torsion,
damping, создание нового solver, commit или push. Frozen scientific modules,
root inventories и historical outputs не изменялись.
"""


def run_postprocessing(
    contract_path: str | Path,
    output_directory: str | Path,
) -> dict[str, Any]:
    """Close the existing RLB-1C-ISO evidence without any new global search.

    Both determinant inventories, spectrum comparison, and selected mode data
    must already exist.  The only determinant solves performed here are the
    two explicitly bounded ISO-06/ISO-07 refinements inside each root's own
    frozen sign-changing bracket.
    """

    contract_file = Path(contract_path).resolve()
    output_dir = Path(output_directory).resolve()
    contract = json.loads(contract_file.read_text(encoding="utf-8"))
    _validate_finite_case_contract(contract)
    initial_inventory_hashes, rlb_manifest, legacy_manifest = _verify_inventory_files(
        contract_file, output_dir
    )
    _blocks, rlb_rows, legacy_rows = _inventory_blocks(
        output_dir, rlb_manifest, legacy_manifest
    )
    rlb_semantic_hash = _semantic_inventory_hash(rlb_rows, rlb_manifest)
    legacy_semantic_hash = _semantic_inventory_hash(legacy_rows, legacy_manifest)
    if rlb_semantic_hash != EXPECTED_RLB_INVENTORY_SHA256:
        raise RuntimeError("RLB inventory differs from the declared frozen semantic hash")
    if legacy_semantic_hash != EXPECTED_LEGACY_INVENTORY_SHA256:
        raise RuntimeError("legacy inventory differs from the declared frozen semantic hash")
    if len(rlb_rows) != 104 or len(legacy_rows) != 104:
        raise RuntimeError("each frozen inventory must contain exactly 104 root slots")

    frozen_references = _verify_frozen_reference_hashes(contract_file, output_dir)
    preliminary_evidence = _validate_preliminary_evidence(output_dir)
    iso04_close_pair = _verify_iso04_close_pair(output_dir)

    spectrum_path = output_dir / "spectrum_comparison.csv"
    original_spectrum_rows = _read_csv(spectrum_path)
    direct_numeric_hash_before = _direct_numeric_projection_sha256(
        original_spectrum_rows
    )
    classified_spectrum_rows, classification = classify_existing_spectrum_rows(
        original_spectrum_rows
    )
    spectrum_evidence = _validate_primary_spectrum_rows(
        classified_spectrum_rows, contract
    )
    _write_csv(spectrum_path, classified_spectrum_rows)
    persisted_spectrum_rows = _read_csv(spectrum_path)
    direct_numeric_hash_after = _direct_numeric_projection_sha256(
        persisted_spectrum_rows
    )
    if direct_numeric_hash_after != direct_numeric_hash_before:
        raise RuntimeError("direct-beta0 numerical payload changed during classification")

    mode_evidence = _validate_mode_evidence(output_dir, contract)
    symmetry_evidence = _validate_saved_symmetry_evidence(output_dir, contract)
    refinement_rows, refinement = _bounded_legacy_arm_exchange_refinement(
        contract, output_dir
    )
    refinement_name = "legacy_arm_exchange_local_refinement.csv"
    _write_csv(output_dir / refinement_name, refinement_rows)

    preliminary_manifest = json.loads(
        (output_dir / "preliminary_gate.json").read_text(encoding="utf-8")
    )
    if preliminary_manifest.get("statuses") != {
        key: SCIENTIFIC_STATUSES[key]
        for key in (
            "RLB-ISO-4PLY-CONSTITUTIVE",
            "RLB-ISO-LEGACY-RECTANGULAR-ADAPTER",
            "RLB-ISO-LOCAL-ARM-EQUIVALENCE",
            "RLB-ISO-SECTION-REDUCTION",
        )
    }:
        raise RuntimeError("preliminary scientific statuses differ from the frozen PASS gate")

    auxiliary_statuses = {
        "AUX-LEGACY-DIRECT-BETA0-FIXED-FIXED": DIRECT_BETA0_AUXILIARY_STATUS,
        "AUX-LEGACY-ARM-EXCHANGE": refinement["status"],
    }
    final_status = {
        "stage": "RLB-1C-ISO",
        "scientific_result_statement": SCIENTIFIC_RESULT_STATEMENT,
        "scientific_statuses": SCIENTIFIC_STATUSES,
        "auxiliary_statuses": auxiliary_statuses,
        "SCIENTIFIC_OVERALL": SCIENTIFIC_OVERALL_STATUS,
        "ORIGINAL-CONTRACT-OVERALL": "PARTIAL_PASS",
        "finite_case_scope": list(REQUIRED_CASE_IDS),
        "used_exactly_four_equal_plies": True,
        "spectral_stack_deg": [0.0, 90.0, 90.0, 0.0],
        "status": SCIENTIFIC_OVERALL_STATUS,
    }
    _write_json(output_dir / "final_status.json", final_status)

    closing_summary = {
        "rlb_inventory_sha256": rlb_semantic_hash,
        "legacy_inventory_sha256": legacy_semantic_hash,
        "preliminary_evidence": preliminary_evidence,
        "spectrum_evidence": spectrum_evidence,
        "mode_evidence": mode_evidence,
        "symmetry_evidence": symmetry_evidence,
        "arm_exchange_refinement": refinement,
        "direct_beta0_auxiliary": {
            "row_count": classification["auxiliary_row_count"],
            "counts_by_kind": classification["auxiliary_counts_by_kind"],
            "numeric_projection_sha256_before": direct_numeric_hash_before,
            "numeric_projection_sha256_after": direct_numeric_hash_after,
            "scientific_role": "AUXILIARY_DIAGNOSTIC",
            "status": "SUPERSEDED_DIAGNOSTIC_ARTIFACT",
            "used_for_scientific_status": False,
            "reason": DIRECT_BETA0_QUALIFICATION_REASON,
        },
        "scientific_statuses": SCIENTIFIC_STATUSES,
        "auxiliary_statuses": auxiliary_statuses,
        "scientific_overall": SCIENTIFIC_OVERALL_STATUS,
    }
    (output_dir / "report.md").write_text(
        _render_closing_report(closing_summary), encoding="utf-8"
    )

    old_run_manifest_path = output_dir / "run_manifest.json"
    old_run_manifest = json.loads(old_run_manifest_path.read_text(encoding="utf-8"))
    if old_run_manifest.get("workflow_state") == "CLOSING_AUDIT_COMPLETE":
        inventory_generation_provenance = old_run_manifest.get(
            "inventory_generation_provenance", {}
        )
    else:
        inventory_generation_provenance = old_run_manifest

    evidence_names = (
        "preliminary_gate.json",
        "isotropic_qbar_checks.csv",
        "four_ply_laminate_checks.csv",
        "section_reduction_checks.csv",
        "legacy_circular_backcompat.csv",
        "local_solution_space_checks.csv",
        "spectrum_comparison.csv",
        "symmetry_checks.csv",
        "mode_shape_comparison.csv",
        "mode_subspace_comparison.csv",
        "mode_joint_residuals.csv",
        "mode_shapes.csv",
        refinement_name,
        "final_status.json",
        "report.md",
    )
    final_inventory_hashes = {
        path.name: _sha256_file(path) for path in _inventory_paths(output_dir)
    }
    if final_inventory_hashes != initial_inventory_hashes:
        raise RuntimeError("closing audit changed frozen inventory bytes")

    run_manifest = {
        "stage": "RLB-1C-ISO",
        "run_manifest_schema_version": 2,
        "generated_utc": _utc_now(),
        "workflow_state": "CLOSING_AUDIT_COMPLETE",
        "contract_sha256": _sha256_file(contract_file),
        "inventory_generation_provenance": inventory_generation_provenance,
        "frozen_semantic_hashes": {
            "RLB": rlb_semantic_hash,
            "LEGACY_RECTANGULAR": legacy_semantic_hash,
        },
        "frozen_inventory_file_sha256": initial_inventory_hashes,
        "frozen_inventory_files_preserved_byte_for_byte": True,
        "frozen_reference_audit": frozen_references,
        "iso04_close_pair": iso04_close_pair,
        "scientific_evidence": {
            "preliminary": preliminary_evidence,
            "spectrum": spectrum_evidence,
            "modes": mode_evidence,
            "symmetry": symmetry_evidence,
        },
        "scientific_evidence_file_sha256": {
            name: _sha256_file(output_dir / name) for name in evidence_names
        },
        "superseded_auxiliary_artifacts": {
            "container": "spectrum_comparison.csv",
            "comparison_kinds": list(DIRECT_BETA0_AUXILIARY_KINDS),
            "row_count": classification["auxiliary_row_count"],
            "classification": "SUPERSEDED_DIAGNOSTIC_ARTIFACT",
            "scientific_role": "AUXILIARY_DIAGNOSTIC",
            "used_for_scientific_status": False,
            "reason": DIRECT_BETA0_QUALIFICATION_REASON,
            "numeric_projection_sha256_before": direct_numeric_hash_before,
            "numeric_projection_sha256_after": direct_numeric_hash_after,
        },
        "legacy_arm_exchange_local_refinement": refinement,
        "thresholds": contract["thresholds"],
        "scientific_statuses": SCIENTIFIC_STATUSES,
        "auxiliary_statuses": auxiliary_statuses,
        "SCIENTIFIC_OVERALL": SCIENTIFIC_OVERALL_STATUS,
        "ORIGINAL-CONTRACT-OVERALL": "PARTIAL_PASS",
        "scientific_result_statement": SCIENTIFIC_RESULT_STATEMENT,
        "global_root_searches_run_in_closing_stage": 0,
        "ritz_solves_run_in_closing_stage": 0,
        "direct_fixed_fixed_solves_run_in_closing_stage": 0,
        "bounded_local_root_refinements_run_in_closing_stage": len(
            refinement_rows
        ),
        "new_solver_modules_created_in_closing_stage": 0,
        "figures_created_in_closing_stage": 0,
        "explicit_exclusions_confirmed": [
            "global_root_scan",
            "Rayleigh_Ritz",
            "Euler_Bernoulli",
            "FEM",
            "beta_sweep",
            "parameter_study",
            "torsion",
            "damping",
            "new_solver",
            "article_work",
            "commit",
            "push",
        ],
    }
    _write_json(old_run_manifest_path, run_manifest)
    return _json_value(run_manifest)


DEFAULT_CONTRACT_PATH = (
    REPOSITORY_ROOT / "tests" / "data" / "reddy_four_ply_isotropic_limit_cases.json"
)
DEFAULT_OUTPUT_DIRECTORY = (
    REPOSITORY_ROOT
    / "results"
    / "laminated_beams"
    / "reddy_four_ply_isotropic_limit_validation"
)


def _build_closing_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--close-existing-results",
        action="store_true",
        help="close already frozen outputs without launching any global search",
    )
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _build_closing_parser().parse_args(argv)
    if not arguments.close_existing_results:
        raise SystemExit(
            "Refusing implicit work: pass --close-existing-results for the frozen closing audit."
        )
    summary = run_postprocessing(arguments.contract, arguments.output)
    printable = {
        "workflow_state": summary["workflow_state"],
        "scientific_overall": summary["SCIENTIFIC_OVERALL"],
        "scientific_result_statement": summary["scientific_result_statement"],
        "global_root_searches_run_in_closing_stage": summary[
            "global_root_searches_run_in_closing_stage"
        ],
        "ritz_solves_run_in_closing_stage": summary[
            "ritz_solves_run_in_closing_stage"
        ],
    }
    print(json.dumps(printable, ensure_ascii=False, sort_keys=True))
    return 0


__all__ = [
    "bounded_sign_change_refinement",
    "classify_existing_spectrum_rows",
    "run_postprocessing",
]


if __name__ == "__main__":
    raise SystemExit(main())
