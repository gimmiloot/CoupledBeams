"""FEM-3AR parser precision regression: synthetic/native saved data, zero solves."""
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import nlsp_fem3a_transient_output as io

ROOT = Path(__file__).resolve().parents[1]
ACTUAL = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893/cases/linear"


@pytest.mark.parametrize("token, quantum", [
    ("0.151891E+01", 1e-5), (".518914E+00", 1e-6),
    ("0.1002595E+01", 1e-6), ("1.002594570", 1e-9),
    ("1.002595e+00", 1e-6), (".259457D-02", 1e-8)])
def test_native_lexical_timestamp_quantum(token, quantum):
    item = io._native_time_precision(token)
    assert item["token"] == token
    assert item["quantum"] == pytest.approx(quantum, rel=0, abs=1e-20)
    assert item["rounding_bound"] == pytest.approx(quantum / 2, rel=0, abs=1e-20)


def test_sta_actual_leading_zero_format_preserves_offset(tmp_path):
    p = tmp_path / "native.sta"
    p.write_text("1 1 1 1 0.100000E+01 0.100000E+01 0.100000E+01\n"
                 "2 1 1 2 0.100259E+01 0.259457E-02 0.259457E-02\n"
                 "2 2 1 2 0.151891E+01 0.518914E+00 0.129729E-02\n")
    result = io.read_transient_sta(p)
    row = result["accepted_increments"][-1]
    assert row["total_time"] == 1.51891
    assert row["native_time_tokens"]["total_time"] == "0.151891E+01"
    assert row["time_quantums"]["total_time"] == 1e-5
    low, high = result["step_offset_rounding_intervals"][2]
    assert low <= 1 <= high
    assert result["actual_dynamic_end"] == .518914


def test_sta_inconsistent_offset_beyond_printed_intervals_rejected(tmp_path):
    p = tmp_path / "native.sta"
    p.write_text("2 1 1 2 0.100259E+01 0.259457E-02 0.259457E-02\n"
                 "2 2 1 2 0.152000E+01 0.518914E+00 0.129729E-02\n")
    with pytest.raises(ValueError, match="Inconsistent transient STA step offset"):
        io.read_transient_sta(p)


def test_sta_all_offset_intervals_intersect_not_only_first_pair(tmp_path):
    p = tmp_path / "native.sta"
    p.write_text("2 1 1 1 1.01 0.01 0.01\n"
                 "2 2 1 1 1.0210 0.0200 0.0100\n"
                 "2 3 1 1 1.0290 0.0300 0.0100\n")
    # Pairwise intervals touch the first, but the final intersection is empty.
    with pytest.raises(ValueError, match="Inconsistent transient STA step offset"):
        io.read_transient_sta(p)


def test_dat_maps_actual_coarse_sta_using_native_rounding(tmp_path):
    sta = tmp_path / "native.sta"
    sta.write_text("2 1 1 2 0.100259E+01 0.259457E-02 0.259457E-02\n")
    rows = io.read_transient_sta(sta)["accepted_increments"]
    dat = tmp_path / "native.dat"
    dat.write_text(" displacements (vx,vy,vz) for set ALL_NODES and time 0.1002595E+01\n\n"
                   " 1 0.000000E+00 -0.004000E+00 0.000000E+00\n\n"
                   " total force (fx,fy,fz) for set OTHER and time 0.1002595E+01\n\n"
                   " 0.000000E+00 0.000000E+00 0.000000E+00\n")
    blocks = list(io.iter_transient_dat(dat, {"ALL_NODES": [1]}, increments=rows))
    assert len(blocks) == 1
    assert blocks[0]["increment"] == 1 and blocks[0]["step"] == 2
    assert blocks[0]["total_time"] == 1.002595
    assert blocks[0]["total_time_quantum"] == 1e-6
    assert blocks[0]["native_total_time_token"] == "0.1002595E+01"
    assert blocks[0]["dynamic_time"] == pytest.approx(.002595)


def test_time_match_outside_combined_rounding_is_rejected(tmp_path):
    row = {"step": 2, "increment": 1, "total_time": 1.002590,
           "step_time": .00259457, "time_rounding_bounds": {"total_time": 5e-6}}
    precision = io._native_time_precision("1.002596")
    with pytest.raises(ValueError, match="unique actual STA"):
        io._dat_time_metadata(precision["value"], 1., 2, [row], precision)


def test_rounding_match_ambiguity_not_resolved_by_nearest_time():
    rows = [{"step": 2, "increment": i, "total_time": 1.002590,
             "step_time": .00259457, "time_rounding_bounds": {"total_time": 5e-6}}
            for i in [1, 2]]
    precision = io._native_time_precision("1.002594570")
    with pytest.raises(ValueError, match="unique actual STA"):
        io._dat_time_metadata(precision["value"], 1., 2, rows, precision)


def test_frd_explicit_increment_native_rounding_contract():
    row = {"step": 2, "increment": 102, "total_time": 1.51891,
           "step_time": .518914, "time_rounding_bounds": {"total_time": 5e-6}}
    precision = io._native_time_precision("1.518914080")
    meta = io._metadata(2, 102, precision["value"], 1., 2, [row], precision)
    assert meta["total_time"] == 1.51891408
    assert meta["sta_total_time"] == 1.51891
    assert meta["sta_step_time"] == .518914
    assert meta["native_total_time_token"] == "1.518914080"
    assert meta["dynamic_time"] == pytest.approx(.51891408)


def test_saved_linear_sta_reparse_without_solver_call():
    result = io.read_transient_sta(ACTUAL / "motion.sta")
    rows = result["accepted_increments"]
    assert len(rows) == 103
    assert sum(r["step"] == 1 for r in rows) == 1
    assert sum(r["step"] == 2 for r in rows) == 102
    assert result["reported_cutbacks"] == 0
    assert result["actual_dynamic_end"] == .518914
    low, high = result["step_offset_rounding_intervals"][2]
    assert low <= 1 <= high
    assert all(r["time_rounding_bounds"]["total_time"] == 5e-6 for r in rows)
