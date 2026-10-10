"""Saved endpoint arithmetic and source guards; no FEM/ODE jobs."""
from pathlib import Path
import inspect

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_recovery_endpoints as endpoints


def test_same_grid_max_signed_difference_and_physical_L2():
    x = np.linspace(0., 1., 41); first = np.zeros((41,7)); second = first.copy()
    second[:,1] = -2e-9
    result = endpoints.difference_metrics(first,second,x)
    assert result["w"]["absolute_max"] == 2e-9
    assert result["w"]["signed_81_minus_41_at_max"] == -2e-9
    assert result["w"]["L2_on_41_point_trapezoidal_grid"] == pytest.approx(2e-9)


def test_full_delta_evolution_sensitivity_does_not_mix_static_and_dynamic():
    x = np.linspace(0.,1.,41)
    names = ("linear_static","linear_final","nonlinear_static","nonlinear_final")
    first = {name:np.zeros((41,7)) for name in names}
    second = {name:value.copy() for name,value in first.items()}
    second["linear_static"][:,2] = 1e-9
    second["nonlinear_static"][:,2] = 4e-9
    second["linear_final"][:,2] = 2e-9
    second["nonlinear_final"][:,2] = 10e-9
    metrics,arrays = endpoints.compare_endpoint_fields(first,second,x)
    assert metrics["correction"]["static"]["v"]["absolute_max"] == pytest.approx(3e-9)
    assert metrics["correction"]["final"]["v"]["absolute_max"] == pytest.approx(8e-9)
    assert metrics["evolution"]["v"]["absolute_max"] == pytest.approx(5e-9)
    np.testing.assert_allclose(arrays["evolution81"][:,2],5e-9)
    assert "physical_orientation" in metrics["full"]["linear_static"]
    assert "physical_orientation" not in metrics["correction"]["static"]
    assert "physical_orientation" not in metrics["evolution"]


def test_physical_orientation_uses_complete_canonical_rotation():
    x = np.linspace(0.,1.,41); first = np.zeros((41,7)); second = first.copy()
    second[:,3:6] = [.003,-.004,.012]
    result = endpoints.difference_metrics(first,second,x)
    assert result["physical_orientation"]["maximum_principal_angle"] == pytest.approx(.013)


def test_cached_diagnostic_checks_source_and_has_no_recovery_calls(tmp_path,monkeypatch):
    output = tmp_path/"recovery_sensitivity_endpoints/medium";output.mkdir(parents=True)
    source = tmp_path/"saved.npy";source.write_bytes(b"immutable source")
    metrics = {"recovery_calls_81":4,"scientific_calls":{"CCX":0,"Radau":0}}
    endpoints.protocol.static.write_json(output/"metrics.json",metrics)
    endpoints.protocol.static.write_json(output/"manifest.json",{"source_hashes":{str(source):endpoints.native.sha(source)},
        "artifact_hashes":{"metrics.json":endpoints.native.sha(output/"metrics.json")}})
    monkeypatch.setattr(endpoints.protocol,"recover_spatial_sections",lambda *a,**k:pytest.fail("Cached endpoint rerecovered"))
    assert endpoints.audit_recovery_endpoints(tmp_path)==metrics
    source.write_bytes(b"changed")
    with pytest.raises(ValueError,match="source changed"):
        endpoints.audit_recovery_endpoints(tmp_path)


def test_no_native_ode_or_eigen_execution_routes():
    source = inspect.getsource(endpoints)
    assert "run_job(" not in source and "run_attempt(" not in source
    assert "solve_ivp" not in source and "Radau(" not in source
    assert "linear_eigenpairs(" not in source and "static_newton(" not in source
    assert "section_count=81" in source


def test_only_existing_medium_fine_scope(tmp_path):
    with pytest.raises(ValueError,match="two existing"):
        endpoints.audit_recovery_endpoints(tmp_path,"refined")


def test_only_derived_model_context_may_resolve_to_exact_first_processing_archive(tmp_path):
    comparison = tmp_path/"comparison_medium";comparison.mkdir()
    source = comparison/"one_d_three_d_comparison.json";source.write_bytes(b"original model context")
    digest = endpoints.native.sha(source)
    archive = comparison/"first_processing_evidence";archive.mkdir()
    (archive/source.name).write_bytes(source.read_bytes())
    source.write_bytes(b"postprocessing added orientation diagnostics")
    endpoints._verify_source(str(source),digest)
    (archive/source.name).write_bytes(b"changed original")
    with pytest.raises(ValueError,match="source changed"):
        endpoints._verify_source(str(source),digest)
