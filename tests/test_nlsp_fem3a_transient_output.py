"""Synthetic and saved-output reader checks; no FEM/ODE/BVP/eigen execution."""
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import nlsp_fem3a_transient_output as output

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/nlsp_nonlinear_static_3d_fem_resume/210b74b8b166997c/cases/medium/nonlinear"


def frd_block(counter, increment, step, time, field, ids=(11, 29), values=None):
    count = 6 if field in ("STRESS", "TOSTRAIN") else 3
    values = np.zeros((len(ids), count)) if values is None else values
    header = "    1PSTEP" + " " * 14 + f"{counter:12d}{increment:12d}{step:12d} "
    time_header = "  100CL  101" + f"{time:12.5E}{len(ids):12d}" + " " * 38 + "1"
    lines = [header, time_header, f" -4  {field}        4    1", " -5  D1          1"]
    lines += [f" -1{node:10d}" + "".join(f"{v:12.5E}" for v in row)
              for node, row in zip(ids, values)]
    return "\n".join(lines + [" -3"]) + "\n"


def two_frames():
    return (frd_block(1, 10, 1, 1., "DISP")
            + frd_block(2, 1, 2, 1.01, "DISP", values=[[0, 0, 0], [0, -.005, 0]])
            + frd_block(3, 1, 2, 1.01, "VELO", values=[[0, 0, 0], [0, .001, 0]]))


def write(path, text):
    path.write_text(text, encoding="utf8")
    return path


def test_frd_static_dynamic_time_offset_and_real_ids(tmp_path):
    frames = list(output.iter_transient_frd(write(tmp_path / "test.frd", two_frames()), [29, 11],
                                            fixed_node_ids=[11]))
    assert len(frames) == 2
    assert frames[0]["step"] == 1 and frames[0]["dynamic_time"] is None
    assert frames[1]["step"] == 2 and frames[1]["increment"] == 1
    assert frames[1]["total_time"] == 1.01
    assert frames[1]["step_time"] == pytest.approx(.01)
    assert frames[1]["dynamic_time"] == pytest.approx(.01)
    assert frames[1]["fields"]["DISP"][0, 1] == -.005
    assert frames[1]["fields"]["VELO"][0, 1] == .001
    assert frames[1]["fixed_displacement_max"] == 0
    assert frames[1]["fixed_velocity_max"] == 0


def test_fixed_face_nonzero_is_visible_not_relabelled(tmp_path):
    text = two_frames().replace("1.00000E-03", "2.00000E-03")
    frames = list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29],
                                            fixed_node_ids=[29]))
    assert frames[-1]["fixed_velocity_max"] == .002
    assert frames[-1]["fixed_displacement_max"] == .005


def test_frd_streaming_does_not_read_all_text(tmp_path, monkeypatch):
    path = write(tmp_path / "test.frd", two_frames())
    monkeypatch.setattr(Path, "read_text", lambda *a, **k: pytest.fail("Nonstreaming read"))
    assert len(list(output.iter_transient_frd(path, [11, 29]))) == 2


@pytest.mark.parametrize("field", ["DISP", "VELO", "FORC", "STRESS", "TOSTRAIN"])
def test_frd_supported_actual_components(tmp_path, field):
    blocks = list(output.iter_transient_frd_blocks(
        write(tmp_path / "test.frd", frd_block(1, 1, 2, 1.01, field)), [11, 29]))
    assert blocks[0]["name"] == field
    assert blocks[0]["values"].shape == (2, 6 if field in ("STRESS", "TOSTRAIN") else 3)


def test_missing_velocity_is_not_invented(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 2, 1.01, "DISP"))
    with pytest.raises(ValueError, match="Missing transient frame fields"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_missing_displacement_is_not_invented(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 2, 1.01, "VELO"))
    with pytest.raises(ValueError, match="Missing transient frame fields"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_missing_nodes_rejected(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 1, 1., "DISP", ids=[11]))
    with pytest.raises(ValueError, match="declared nodal count"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_duplicate_nodes_rejected(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 1, 1., "DISP", ids=[11, 11]))
    with pytest.raises(ValueError, match="Duplicate FRD node"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_wrong_ids_rejected(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 1, 1., "DISP", ids=[11, 30]))
    with pytest.raises(ValueError, match="Incomplete DISP nodes"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_duplicate_frame_field_rejected(tmp_path):
    text = frd_block(1, 1, 1, 1., "DISP") + frd_block(2, 1, 1, 1., "DISP")
    with pytest.raises(ValueError, match="Duplicate field"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29]))


def test_duplicate_counter_rejected(tmp_path):
    text = frd_block(1, 1, 1, 1., "DISP") + frd_block(1, 1, 2, 1.01, "DISP")
    with pytest.raises(ValueError, match="dataset counter"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29]))


def test_partial_trailing_block_rejected(tmp_path):
    text = frd_block(1, 1, 1, 1., "DISP").removesuffix(" -3\n")
    with pytest.raises(ValueError, match="Partial trailing"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29]))


def test_unterminated_frame_rejected(tmp_path):
    text = frd_block(1, 1, 1, 1., "DISP").removesuffix(" -3\n") + frd_block(2, 1, 2, 1.01, "VELO")
    with pytest.raises(ValueError, match="Unterminated"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29]))


def test_modal_frd_rejected(tmp_path):
    path = write(tmp_path / "test.frd", "    1PMODE        1\n" + two_frames())
    with pytest.raises(ValueError, match="Modal FRD"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_nonfinite_fields_rejected(tmp_path):
    path = write(tmp_path / "test.frd", two_frames().replace("1.00000E-03", "1.00000E+999"))
    with pytest.raises(ValueError, match="Nonfinite"):
        list(output.iter_transient_frd(path, [11, 29]))


def test_dynamic_frame_before_release_rejected(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 2, .5, "DISP"))
    with pytest.raises(ValueError, match="predates"):
        list(output.iter_transient_frd_blocks(path, [11, 29]))


def test_unknown_fixed_node_rejected(tmp_path):
    with pytest.raises(ValueError, match="Fixed-face"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", two_frames()), [11, 29], fixed_node_ids=[30]))


def test_sta_two_step_actual_prefix(tmp_path):
    text = ("STEP INC ATT ITRS TOT TIME STEP TIME INC TIME\n"
            "1 1 1 1 1.000000E+00 1.000000E+00 1.000000E+00\n"
            "2 1 1 2 1.010000E+00 1.000000E-02 1.000000E-02\n"
            "2 2 2 3 1.015000E+00 1.500000E-02 5.000000E-03\n")
    sta = output.read_transient_sta(write(tmp_path / "test.sta", text))
    assert sta["status"] == "PARSED"
    assert sta["actual_dynamic_end"] == .015
    assert sta["step_offsets"][2] == pytest.approx(1.)
    assert sta["reported_attempts"] == 4 and sta["reported_cutbacks"] == 1
    assert len(sta["accepted_increments"]) == 3


@pytest.mark.parametrize("row", ["2 1 1 2 1.0100 .01 .01", "2 2 1 2 1.0050 .005 .005"])
def test_sta_duplicate_decreasing_rejected(tmp_path, row):
    text = "2 1 1 2 1.0100 .01 .01\n" + row + "\n"
    with pytest.raises(ValueError, match="Duplicate/decreasing"):
        output.read_transient_sta(write(tmp_path / "test.sta", text))


def test_sta_missing_is_visible(tmp_path):
    assert output.read_transient_sta(tmp_path / "missing.sta")["status"] == "NOT_PRESENT"


def test_frd_sta_increment_crosscheck(tmp_path):
    increments = [{"step": 1, "increment": 10, "total_time": 1., "step_time": 1.},
                  {"step": 2, "increment": 1, "total_time": 1.01, "step_time": .01}]
    frames = list(output.iter_transient_frd(write(tmp_path / "test.frd", two_frames()), [11, 29], increments=increments))
    assert frames[-1]["sta_step_time"] == .01
    increments[-1]["total_time"] = 1.02
    with pytest.raises(ValueError, match="time disagreement"):
        list(output.iter_transient_frd(tmp_path / "test.frd", [11, 29], increments=increments))


def dat_block(name, setname, time, ids=(11, 29), values=None):
    values = np.zeros((len(ids), 3)) if values is None else values
    lines = [f" {name} (vx,vy,vz) for set {setname} and time {time:.7E}", ""]
    lines += [f"{node} " + " ".join(f"{v:.6E}" for v in row) for node, row in zip(ids, values)]
    return "\n".join(lines) + "\n\n"


def test_dat_stream_fields_with_sta_real_increment(tmp_path):
    text = dat_block("displacements", "ALL_NODES", 1.) + dat_block("displacements", "ALL_NODES", 1.01)
    text += dat_block("forces", "LEFT_FIXED", 1.01, ids=[11])
    increments = [{"step": 1, "increment": 10, "total_time": 1., "step_time": 1.},
                  {"step": 2, "increment": 1, "total_time": 1.01, "step_time": .01}]
    blocks = list(output.iter_transient_dat(write(tmp_path / "test.dat", text), {"ALL_NODES": [29, 11], "LEFT_FIXED": [11]}, increments=increments))
    assert len(blocks) == 3
    assert blocks[0]["dynamic_time"] is None
    assert blocks[1]["step"] == 2 and blocks[1]["increment"] == 1
    assert blocks[2]["name"] == "FORC"


def test_dat_duplicate_block_rejected(tmp_path):
    text = dat_block("displacements", "ALL_NODES", 1.) * 2
    with pytest.raises(ValueError, match="Duplicate transient DAT"):
        list(output.iter_transient_dat(write(tmp_path / "test.dat", text), {"ALL_NODES": [11, 29]}))


def test_dat_incomplete_nodes_rejected(tmp_path):
    text = dat_block("displacements", "ALL_NODES", 1., ids=[11])
    with pytest.raises(ValueError, match="Incomplete DAT"):
        list(output.iter_transient_dat(write(tmp_path / "test.dat", text), {"ALL_NODES": [11, 29]}))


def test_dat_unknown_quantity_ends_nodal_block(tmp_path):
    text = dat_block("displacements", "ALL_NODES", 1.)
    text += " total internal energy for set SOLID and time 1.0000000E+00\n\n 2.000000E-07\n"
    assert len(list(output.iter_transient_dat(write(tmp_path / "test.dat", text), {"ALL_NODES": [11, 29]}))) == 1


def test_energy_native_totals_dynamic_time_and_definition(tmp_path):
    text = (" total internal energy for set SOLID and time 1.0000000E+00\n\n 2.000000E-07\n"
            " total kinetic energy for set SOLID and time 1.0000000E+00\n\n 0.000000E+00\n"
            " total internal energy for set SOLID and time 1.0100000E+00\n\n 1.900000E-07\n"
            " total kinetic energy for set SOLID and time 1.0100000E+00\n\n 1.000000E-08\n")
    result = output.parse_transient_dat_energies(write(tmp_path / "test.dat", text), element_set="SOLID")
    assert result["status"] == "PARSED"
    assert result["records"][0]["dynamic_time"] is None
    assert result["records"][1]["dynamic_time"] == pytest.approx(.01)
    assert result["records"][1]["mechanical_energy"] == pytest.approx(2e-7)
    assert "no removed GRAV" in result["energy_definition"]


def test_energy_missing_is_not_invented(tmp_path):
    result = output.parse_transient_dat_energies(write(tmp_path / "test.dat", "displacements\n"))
    assert result["status"] == "NOT_PRESENT" and result["records"] == []


def test_energy_single_type_is_partial(tmp_path):
    text = " total internal energy for set SOLID and time 1.01\n 1.2e-7\n"
    result = output.parse_transient_dat_energies(write(tmp_path / "test.dat", text))
    assert result["status"] == "PARTIAL"
    assert "mechanical_energy" not in result["records"][0]


def test_energy_duplicate_is_rejected(tmp_path):
    text = " total internal energy for set SOLID and time 1.01\n 1.2e-7\n" * 2
    with pytest.raises(ValueError, match="Duplicate native"):
        output.parse_transient_dat_energies(write(tmp_path / "test.dat", text))


def test_energy_incomplete_trailing_scalar_rejected(tmp_path):
    text = " total internal energy for set SOLID and time 1.01\n"
    with pytest.raises(ValueError, match="Partial trailing"):
        output.parse_transient_dat_energies(write(tmp_path / "test.dat", text))


def test_saved_parent_static_frd_readonly():
    arrays = np.load(SOURCE / "static_nodal_results.npz")
    frames = list(output.iter_transient_frd(SOURCE / "static.frd", arrays["node_ids"]))
    assert len(frames) == 1
    assert frames[0]["step"] == 1 and frames[0]["increment"] == 10
    assert frames[0]["dynamic_time"] is None
    assert set(frames[0]["fields"]) >= {"DISP", "FORC", "STRESS", "TOSTRAIN"}
    assert np.array_equal(frames[0]["fields"]["DISP"], arrays["U_FRD"])


def test_saved_parent_sta_readonly():
    sta = output.read_transient_sta(SOURCE / "static.sta")
    assert len(sta["accepted_increments"]) == 10
    assert sta["actual_dynamic_end"] is None
    assert sta["actual_total_end"] == 1.


def test_saved_parent_dat_readonly():
    arrays = np.load(SOURCE / "static_nodal_results.npz")
    # Only ALL_NODES is required; both support-force sets are safely skipped.
    blocks = list(output.iter_transient_dat(SOURCE / "static.dat", {"ALL_NODES": arrays["node_ids"]}))
    assert len(blocks) == 1 and blocks[0]["name"] == "DISP"
    assert np.array_equal(blocks[0]["values"], arrays["U"])


def test_reader_has_no_scientific_execution_entrypoints():
    source = Path(output.__file__).read_text(encoding="utf8")
    for forbidden in ("subprocess", "run_job(", "solve_ivp(", "generate_mesh(", "fem2_static_newton("):
        assert forbidden not in source


def test_same_increment_conflicting_time_rejected(tmp_path):
    text = frd_block(1, 1, 1, .5, "DISP") + frd_block(2, 1, 1, .6, "DISP")
    with pytest.raises(ValueError, match="Conflicting times"):
        list(output.iter_transient_frd(write(tmp_path / "test.frd", text), [11, 29]))


def test_partial_prefix_is_not_padded(tmp_path):
    path = write(tmp_path / "test.frd", two_frames())
    frames = list(output.iter_transient_frd(path, [11, 29]))
    assert frames[-1]["dynamic_time"] == pytest.approx(.01)
    assert len([f for f in frames if f["step"] == 2]) == 1


def test_negative_total_time_rejected(tmp_path):
    path = write(tmp_path / "test.frd", frd_block(1, 1, 1, -.1, "DISP"))
    with pytest.raises(ValueError, match="Negative transient"):
        list(output.iter_transient_frd(path, [11, 29]))


def native_stdout(time=1.01):
    return (f" increment 1 attempt 1\n actual step time={time - 1:.6e}\n actual total time={time:.6e}\n"
            " initial energy (at start of step) = 2.000000e-07\n\n"
            " since start of the step: \n external work = 0.000000e+00\n"
            " work performed by the damping forces = 0.000000e+00\n netto work = 0.000000e+00\n\n"
            " actual energy: \n internal energy = 1.900000e-07\n kinetic energy = 1.000000e-08\n"
            " elastic contact energy = 0.000000e+00\n energy lost due to friction = 0.000000e+00\n"
            " total energy  = 2.000000e-07\n energy increase = 0.000000e+00\n"
            " energy balance (absolute) = 0.000000e+00 \n energy balance (relative) = 0.000000 % \n")


def test_native_stdout_work_damping_time_and_denominator(tmp_path):
    result = output.parse_transient_stdout_energies(write(tmp_path / "test.stdout", native_stdout()))
    record = result["records"][0]
    assert result["status"] == "PARSED"
    assert record["dynamic_time"] == pytest.approx(.01)
    assert record["external_work"] == 0 and record["damping_work"] == 0
    assert record["mechanical_energy"] == pytest.approx(2e-7)
    assert record["relative_change_to_initial_step_energy"] == pytest.approx(0)
    assert "native_energy_balance_relative_percent" in record
    assert "history-dependent" in result["native_relative_balance_qualification"]


def test_native_stdout_sta_acceptance_mapping(tmp_path):
    row = {"step": 2, "increment": 1, "total_time": 1.01, "step_time": .01}
    result = output.parse_transient_stdout_energies(write(tmp_path / "test.stdout", native_stdout()), increments=[row])
    assert result["records"][0]["step"] == 2 and result["records"][0]["increment"] == 1
    assert result["records"][0]["acceptance"] == "STA_ACCEPTED_INCREMENT"


def test_native_stdout_no_time_does_not_invent_it(tmp_path):
    text = native_stdout().split(" initial energy", 1)[1]
    result = output.parse_transient_stdout_energies(write(tmp_path / "test.stdout", " initial energy" + text))
    assert result["records"][0]["total_time"] is None
    assert result["records"][0]["dynamic_time"] is None


def test_native_stdout_partial_energy_remains_partial(tmp_path):
    text = "actual total time=1.01\n initial energy (at start of step) = 2.000000e-7\n"
    result = output.parse_transient_stdout_energies(write(tmp_path / "test.stdout", text))
    assert result["status"] == "PARTIAL"
    assert "mechanical_energy" not in result["records"][0]


def test_saved_static_stdout_has_no_invented_energy():
    result = output.parse_transient_stdout_energies(SOURCE / "static.stdout.txt")
    assert result["status"] == "NOT_PRESENT" and result["records"] == []
