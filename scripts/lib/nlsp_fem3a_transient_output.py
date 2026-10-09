"""Scoped, streaming CalculiX 2.22 static-to-transient output readers.

FRD PSTEP names output counter, actual increment and actual step. Installed
nonlingeo.c passes ttime+time to frd: the 100CL number is total time. DAT
printout.f uses the same total time. A static preload is never relabelled a
dynamic t=0 frame. Only the historical fixed-width numeric decoder is reused;
this module neither runs a solver nor interprets modal results.
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import _static_frd_record

_FRD_COMPONENTS = {"DISP": 3, "VELO": 3, "FORC": 3,
                   "STRESS": 6, "TOSTRAIN": 6, "ENER": 1}
_DAT_NAMES = {"displacements": "DISP", "velocities": "VELO", "forces": "FORC"}
_DAT_HEADER = re.compile(
    r"^\s*(displacements|velocities|forces)\s+\([^)]*\)\s+for set\s+(\S+)\s+and time\s+(\S+)", re.I)
_ENERGY_HEADER = re.compile(
    r"^\s*total (internal|kinetic) energy for set\s+(\S+)\s+and time\s+(\S+)", re.I)
_ANY_DAT_HEADER = re.compile(r"for set\s+\S+\s+and time\s+", re.I)


def _finite_number(text: str) -> float:
    number = float(text.replace("D", "E").replace("d", "E"))
    if not math.isfinite(number):
        raise ValueError("Nonfinite transient output number")
    return number


def _expected_ids(ids: Iterable[int]) -> np.ndarray:
    array = np.asarray(list(ids), dtype=np.int64)
    if array.ndim != 1 or not len(array) or len(set(array.tolist())) != len(array):
        raise ValueError("Expected node IDs must be nonempty and unique")
    return array


def read_transient_sta(path: str | Path, *, dynamic_step: int = 2) -> dict:
    """Read actual accepted increments; never pad a missing time prefix.

    Native .sta records describe accepted increments. 'ATT' reports how many
    attempts preceded acceptance, not separate rejected-increment trajectories.
    Duplicate or decreasing actual step/increment/time records are rejected.
    """
    rows = []
    if not Path(path).exists():
        return {"status": "NOT_PRESENT", "accepted_increments": [],
                "step_offsets": {}, "actual_dynamic_end": None}
    with Path(path).open(encoding="utf8", errors="strict") as stream:
        for raw in stream:
            tokens = raw.split()
            if len(tokens) != 7 or not all(re.fullmatch(r"\d+", s) for s in tokens[:4]):
                continue
            row = dict(zip(("step", "increment", "attempt", "iterations"),
                           map(int, tokens[:4])))
            row.update(zip(("total_time", "step_time", "increment_time"),
                           map(_finite_number, tokens[4:])))
            if row["step"] < 1 or row["increment"] < 1 or row["attempt"] < 1:
                raise ValueError("Invalid transient STA step/increment/attempt")
            if row["step_time"] <= 0 or row["increment_time"] <= 0:
                raise ValueError("Nonpositive accepted increment time")
            if rows:
                last = rows[-1]
                if ((row["step"], row["increment"]) <= (last["step"], last["increment"])
                        or row["total_time"] <= last["total_time"]):
                    raise ValueError("Duplicate/decreasing transient STA increment")
                if row["step"] == last["step"] and row["step_time"] <= last["step_time"]:
                    raise ValueError("Decreasing transient STA step time")
            rows.append(row)
    offsets = {}
    for row in rows:
        offset = row["total_time"] - row["step_time"]
        old = offsets.setdefault(row["step"], offset)
        # Native STA E13.6 timestamps carry only seven significant digits.
        if abs(offset - old) > 2e-6 * max(1., abs(row["total_time"])):
            raise ValueError("Inconsistent transient STA step offset")
    dynamic = [r for r in rows if r["step"] == dynamic_step]
    return {"status": "PARSED" if rows else "NO_ACCEPTED_INCREMENTS",
            "accepted_increments": rows, "step_offsets": offsets,
            "actual_dynamic_end": dynamic[-1]["step_time"] if dynamic else None,
            "actual_total_end": rows[-1]["total_time"] if rows else None,
            "reported_attempts": sum(r["attempt"] for r in rows),
            "reported_cutbacks": sum(r["attempt"] - 1 for r in rows),
            "timestamp_precision": "native STA E13.6; rounded printed timestamps"}


def _metadata(step: int, increment: int, total_time: float, static_end_time: float,
              dynamic_step: int, increments: list[dict] | None = None) -> dict:
    if total_time < 0:
        raise ValueError("Negative transient total time")
    if not math.isfinite(static_end_time) or static_end_time < 0:
        raise ValueError("Invalid static-end offset")
    if step not in (1, dynamic_step):
        raise ValueError(f"Unexpected step {step} in two-step transient output")
    meta = {"step": step, "increment": increment, "total_time": total_time,
            "step_time": total_time if step == 1 else total_time - static_end_time,
            "dynamic_time": total_time - static_end_time if step == dynamic_step else None,
            "time_source": "native total time; documented preload offset"}
    if increments is not None:
        matched = [r for r in increments if r["step"] == step and r["increment"] == increment]
        if len(matched) != 1:
            raise ValueError("FRD frame lacks a unique actual STA increment")
        row = matched[0]
        if abs(row["total_time"] - total_time) > 2e-6 * max(1., abs(total_time)):
            raise ValueError("FRD/STA total-time disagreement")
        meta["sta_total_time"] = row["total_time"]
        meta["sta_step_time"] = row["step_time"]
    if step == dynamic_step and meta["dynamic_time"] < -1e-10:
        raise ValueError("Dynamic frame predates the static preload end")
    return meta


def iter_transient_frd_blocks(path: str | Path, expected_node_ids: Iterable[int], *,
                              static_end_time: float = 1., dynamic_step: int = 2,
                              increments: list[dict] | None = None) -> Iterator[dict]:
    """Yield one validated nodal dataset, retaining at most one block in memory."""
    ids = _expected_ids(expected_node_ids)
    expected = set(ids.tolist())
    pstep = None
    total_time = None
    declared_count = None
    field = None
    nodes = None
    with Path(path).open(encoding="utf8", errors="strict") as stream:
        for number, raw in enumerate(stream, 1):
            s = raw.strip()
            if s.startswith("1PMODE") or "MODAL" in raw[:75]:
                raise ValueError("Modal FRD supplied to transient reader")
            if s.startswith("1PSTEP"):
                if nodes is not None:
                    raise ValueError("Unterminated transient FRD field")
                try:
                    pstep = tuple(int(t) for t in s.split()[-3:])
                except ValueError as exc:
                    raise ValueError(f"Invalid FRD PSTEP at line {number}") from exc
                if len(pstep) != 3 or pstep[0] < 1 or pstep[1] < 0:
                    raise ValueError("Invalid FRD step/increment metadata")
                total_time = None
            elif s.startswith("100C"):
                if nodes is not None:
                    raise ValueError("Unterminated transient FRD field")
                total_time = _finite_number(raw[12:24])
                declared_count = int(raw[24:36])
                # A binary dataset cannot be interpreted as ASCII.
                if len(raw.rstrip("\r\n")) > 74 and raw[74] == "2":
                    raise ValueError("Binary FRD requires a distinct decoder")
            elif raw[:3].strip() == "-4":
                if nodes is not None:
                    raise ValueError("Unterminated transient FRD field")
                tokens = s.split()
                field = tokens[1]
                if field in _FRD_COMPONENTS:
                    if pstep is None or total_time is None:
                        raise ValueError("FRD field lacks actual step/time metadata")
                    if declared_count != len(ids):
                        raise ValueError("FRD declared nodal count differs from source mesh")
                    nodes = {}
            elif nodes is not None and raw[:3].strip() == "-1":
                node, values = _static_frd_record(raw.rstrip("\r\n"), _FRD_COMPONENTS[field])
                if node in nodes:
                    raise ValueError(f"Duplicate FRD node {node} in {field}")
                nodes[node] = values
            elif nodes is not None and raw[:3].strip() == "-3":
                if set(nodes) != expected:
                    raise ValueError(f"Incomplete {field} nodes: expected {len(ids)}, got {len(nodes)}")
                meta = _metadata(pstep[2], pstep[1], total_time, static_end_time,
                                 dynamic_step, increments)
                yield {**meta, "dataset_counter": pstep[0], "name": field,
                       "node_ids": ids, "values": np.asarray([nodes[int(n)] for n in ids]),
                       "precision": "native FRD float32/E12.5 values"}
                nodes = None
                field = None
    if nodes is not None:
        raise ValueError("Partial trailing transient FRD field")


def iter_transient_frd(path: str | Path, expected_node_ids: Iterable[int], *,
                       static_end_time: float = 1., dynamic_step: int = 2,
                       fixed_node_ids: Iterable[int] = (),
                       required_dynamic_fields: tuple[str, ...] = ("DISP", "VELO"),
                       increments: list[dict] | None = None) -> Iterator[dict]:
    """Yield chronological complete frames, with static and dynamic labels intact."""
    ids = _expected_ids(expected_node_ids)
    fixed = set(map(int, fixed_node_ids))
    if not fixed.issubset(set(ids.tolist())):
        raise ValueError("Fixed-face node is absent from source mesh")
    fixed_rows = np.asarray([i for i, node in enumerate(ids) if int(node) in fixed], dtype=int)
    frame = None
    previous_key = None
    previous_counter = 0

    def finish(item):
        required = required_dynamic_fields if item["step"] == dynamic_step else ("DISP",)
        missing = set(required) - set(item["fields"])
        if missing:
            raise ValueError(f"Missing transient frame fields {sorted(missing)}")
        for name, label in (("DISP", "fixed_displacement_max"), ("VELO", "fixed_velocity_max")):
            item[label] = (float(np.max(np.abs(item["fields"][name][fixed_rows])))
                           if len(fixed_rows) and name in item["fields"] else None)
        return item

    for block in iter_transient_frd_blocks(path, ids, static_end_time=static_end_time,
                                           dynamic_step=dynamic_step, increments=increments):
        key = (block["step"], block["increment"], block["total_time"])
        if block["dataset_counter"] <= previous_counter:
            raise ValueError("Duplicate/decreasing FRD dataset counter")
        previous_counter = block["dataset_counter"]
        if previous_key is not None and key < previous_key:
            raise ValueError("Decreasing transient FRD frame order")
        if previous_key is not None and key[:2] == previous_key[:2] and key != previous_key:
            raise ValueError("Conflicting times for the same actual FRD increment")
        if frame is None or key != previous_key:
            if frame is not None:
                yield finish(frame)
            frame = {k: v for k, v in block.items() if k not in ("name", "values", "precision")}
            frame["fields"] = {}
            frame["dataset_counters"] = {}
        if block["name"] in frame["fields"]:
            raise ValueError("Duplicate field in transient FRD frame")
        frame["fields"][block["name"]] = block["values"]
        frame["dataset_counters"][block["name"]] = block["dataset_counter"]
        previous_key = key
    if frame is not None:
        yield finish(frame)


def _dat_time_metadata(total_time: float, static_end_time: float,
                       dynamic_step: int, increments: list[dict] | None) -> dict:
    if increments is not None:
        tolerance = 6e-8 * max(1., abs(total_time))  # DAT E14.7 timestamps
        rows = [r for r in increments if abs(r["total_time"] - total_time) <= tolerance]
        # STA timestamps have fewer digits than DAT: allow their own rounding.
        if not rows:
            rows = [r for r in increments if abs(r["total_time"] - total_time)
                    <= 6e-7 * max(1., abs(total_time))]
        if len(rows) != 1:
            raise ValueError("DAT timestamp lacks a unique actual STA increment")
        row = rows[0]
        return _metadata(row["step"], row["increment"], total_time, static_end_time,
                         dynamic_step, increments)
    step = 1 if total_time <= static_end_time else dynamic_step
    item = _metadata(step, None, total_time, static_end_time, dynamic_step)
    item["time_source"] = "DAT total-time offset only; increment not independently available"
    return item


def iter_transient_dat(path: str | Path, expected_node_sets: dict[str, Iterable[int]], *,
                       static_end_time: float = 1., dynamic_step: int = 2,
                       increments: list[dict] | None = None) -> Iterator[dict]:
    """Stream DAT nodal blocks (DAT U/RF; V only if actually present).

    The 2.22 manual does not promise structural V on NODE PRINT; FRD VELO is
    the primary velocity output. Unknown DAT sets/quantities are skipped.
    """
    sets = {name.upper(): _expected_ids(ids) for name, ids in expected_node_sets.items()}
    header = None
    values = {}
    seen = set()

    def finish():
        if header is None:
            return None
        name, setname, time = header
        if setname not in sets:
            return None
        ids = sets[setname]
        if set(values) != set(ids.tolist()):
            raise ValueError(f"Incomplete DAT {name} {setname}")
        key = (name, setname, time)
        if key in seen:
            raise ValueError("Duplicate transient DAT nodal block")
        seen.add(key)
        return {**_dat_time_metadata(time, static_end_time, dynamic_step, increments),
                "name": name, "set": setname, "node_ids": ids,
                "values": np.asarray([values[int(n)] for n in ids]),
                "precision": "native DAT E13.6 nodal values"}

    with Path(path).open(encoding="utf8", errors="strict") as stream:
        for raw in stream:
            found = _DAT_HEADER.match(raw)
            if found or _ANY_DAT_HEADER.search(raw):
                block = finish()
                if block is not None:
                    yield block
                values = {}
                header = ((_DAT_NAMES[found[1].lower()], found[2].upper(), _finite_number(found[3]))
                          if found else None)
            elif header is not None and raw.strip():
                tokens = raw.split()
                if len(tokens) != 4 or not re.fullmatch(r"\d+", tokens[0]):
                    raise ValueError("Invalid transient DAT nodal row")
                node = int(tokens[0])
                if node in values:
                    raise ValueError("Duplicate transient DAT node")
                values[node] = np.asarray([_finite_number(v) for v in tokens[1:]])
    block = finish()
    if block is not None:
        yield block


def parse_transient_dat_energies(path: str | Path, *, static_end_time: float = 1.,
                                 dynamic_step: int = 2, increments: list[dict] | None = None,
                                 element_set: str | None = None) -> dict:
    """Read native *EL PRINT,TOTALS=ONLY ELSE/ELKE; do not invent missing energy."""
    pending = None
    records = {}
    with Path(path).open(encoding="utf8", errors="strict") as stream:
        for raw in stream:
            found = _ENERGY_HEADER.match(raw)
            if found:
                if pending is not None:
                    raise ValueError("Missing native DAT total-energy scalar")
                pending = (found[1].lower(), found[2].upper(), _finite_number(found[3]))
            elif pending is not None and raw.strip():
                if len(raw.split()) != 1:
                    raise ValueError("Invalid native DAT total-energy scalar")
                value = _finite_number(raw.strip())
                kind, setname, time = pending
                pending = None
                if element_set is not None and setname != element_set.upper():
                    continue
                record = records.setdefault((setname, time), {
                    **_dat_time_metadata(time, static_end_time, dynamic_step, increments),
                    "element_set": setname})
                key = "internal_energy" if kind == "internal" else "kinetic_energy"
                if key in record:
                    raise ValueError("Duplicate native DAT total energy")
                record[key] = value
    if pending is not None:
        raise ValueError("Partial trailing native DAT total-energy block")
    ordered = sorted(records.values(), key=lambda r: (r["total_time"], r["element_set"]))
    complete = all("internal_energy" in r and "kinetic_energy" in r for r in ordered)
    for record in ordered:
        if "internal_energy" in record and "kinetic_energy" in record:
            record["mechanical_energy"] = record["internal_energy"] + record["kinetic_energy"]
    return {"status": "PARSED" if ordered and complete else "PARTIAL" if ordered else "NOT_PRESENT",
            "records": ordered, "energy_definition": "native internal + kinetic; no removed GRAV potential",
            "precision": "native total energy E13.6; output rounding retained"}


def parse_transient_stdout_energies(path: str | Path, *, static_end_time: float = 1.,
                                    dynamic_step: int = 2,
                                    increments: list[dict] | None = None) -> dict:
    """Read native printenergy.c diagnostics, including external/damping work.

    Native 'energy balance (relative)' has a history-dependent denominator;
    it is retained by its original name and never called drift relative to E0.
    Rounded stdout scalars are separate from DAT energy output. Records without
    an actual timestamp remain unassigned, rather than acquiring invented time.
    """
    names = {
        "initial energy (at start of step)": "initial_step_energy",
        "external work": "external_work",
        "work performed by the damping forces": "damping_work",
        "netto work": "net_work",
        "internal energy": "internal_energy",
        "kinetic energy": "kinetic_energy",
        "elastic contact energy": "elastic_contact_energy",
        "energy lost due to friction": "friction_energy",
        "total energy": "native_total_energy",
        "energy increase": "native_energy_increase",
        "energy balance (absolute)": "native_energy_balance_absolute",
        "energy balance (relative)": "native_energy_balance_relative_percent",
    }
    current_time = None
    current_step_time = None
    current_increment = None
    current_attempt = None
    record = None
    records = []

    def finish():
        if record is None:
            return
        if "internal_energy" in record and "kinetic_energy" in record:
            record["mechanical_energy"] = record["internal_energy"] + record["kinetic_energy"]
            e0 = record.get("initial_step_energy")
            if e0 is not None and e0 != 0:
                record["relative_change_to_initial_step_energy"] = (record["mechanical_energy"] - e0) / abs(e0)
        records.append(record)

    with Path(path).open(encoding="utf8", errors="strict") as stream:
        for raw in stream:
            found = re.match(r"\s*increment\s+(\d+)\s+attempt\s+(\d+)", raw, re.I)
            if found:
                current_increment, current_attempt = map(int, found.groups())
            found = re.match(r"\s*actual (total|step) time\s*=\s*(\S+)", raw, re.I)
            if found:
                value = _finite_number(found[2])
                if found[1].lower() == "total":
                    current_time = value
                else:
                    current_step_time = value
            if "=" not in raw:
                continue
            key, value = raw.split("=", 1)
            name = names.get(" ".join(key.lower().split()))
            if name is None:
                continue
            tokens = value.strip().split()
            if not tokens:
                raise ValueError("Missing native stdout energy scalar")
            scalar = _finite_number(tokens[0])
            if name == "initial_step_energy":
                finish()
                if current_time is None:
                    record = {"step": None, "increment": current_increment,
                              "total_time": None, "step_time": current_step_time,
                              "dynamic_time": None, "time_source": "NOT_PRESENT"}
                else:
                    record = _dat_time_metadata(current_time, static_end_time,
                                                dynamic_step, increments)
                record["stdout_increment"] = current_increment
                record["stdout_attempt"] = current_attempt
                record["stdout_step_time"] = current_step_time
                record["acceptance"] = ("STA_ACCEPTED_INCREMENT" if increments is not None
                                        and current_time is not None else "NOT_INDEPENDENTLY_VERIFIED")
            if record is None:
                # Unattached lines do not define a complete printenergy record.
                continue
            if name in record:
                raise ValueError("Duplicate scalar in native stdout energy record")
            record[name] = scalar
    finish()
    required = {"initial_step_energy", "external_work", "damping_work", "internal_energy", "kinetic_energy"}
    complete = all(required.issubset(r) for r in records)
    return {"status": "PARSED" if records and complete else "PARTIAL" if records else "NOT_PRESENT",
            "records": records, "precision": "native printf %e energy/work; seven significant digits",
            "native_relative_balance_qualification": "history-dependent native denominator, not E0 drift",
            "energy_definition": "internal + kinetic; removed GRAV potential excluded"}
