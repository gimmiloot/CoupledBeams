"""Bounded one-T1 continuation checks; saved data and orchestration only.

No new ODE, multiprecision projection, BVP or eigensystem is evaluated here.
The source's failed strict check is preserved independently of feasibility.
"""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest

from scripts.analysis import prepare_planar_initial_state as cli
from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
from scripts.lib import planar_prepared_initial_state as prep
from scripts.lib import planar_second_order_axial_response as leading
from scripts.lib import weakly_nonlinear_planar_dynamics as planar

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/planar_prepared_feasibility/284a4039177391d1"
BASELINE_HEAD="a39c196127eae0cb1f44ade5a924be8ea5acc784"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(autouse=True)
def no_new_scientific_runs(monkeypatch):
    import scipy.integrate
    import scipy.linalg
    import mpmath
    def forbidden(*args,**kwargs):
        raise AssertionError("One-T1 targeted tests may not run ODE/MP/BVP/eigensolvers")
    for module,names in (
            (scipy.integrate,("solve_ivp","Radau")),
            (scipy.linalg,("eigh","eig")),
            (np.linalg,("eigh","eig","eigvals")),
            (leading,("eigh",)),(planar,("eigh",)),
            (prep,("stable_initial_projection","project_saved_legendre","high_precision_source_jets")),
            (mpmath,("workdps",)),
            (runner,("integrate_case",))):
        for name in names:
            monkeypatch.setattr(module,name,forbidden)


@pytest.fixture(scope="module")
def config():
    return read(ROOT/"data/input/planar_prepared_one_T1.json")


@pytest.fixture(scope="module")
def restored(config):
    if not (SOURCE/"manifest.json").is_file():
        pytest.skip("Immutable short source missing; no reproduction is launched")
    return cli.load_one_T1_source(config)


def test_source_settings_and_actual_horizon_are_verified_against_own_manifest(restored):
    manifest=read(SOURCE/"manifest.json")
    for name,digest in manifest["artifact_hashes"].items():
        assert sha(SOURCE/name)==digest
    assert restored["source"].resolve()==SOURCE.resolve()
    assert restored["manifest"]==manifest
    assert restored["summary"]==read(SOURCE/"summary.json")
    t=restored["short_time"]
    assert t[0]==0. and np.all(np.diff(t)>0.)
    assert t[-1]==restored["summary"]["sampling"]["horizon"]
    assert t[-1]==.1*restored["summary"]["background"]["T1"]


@pytest.mark.parametrize("p",(48,64))
def test_restored_full_coordinates_equal_source_without_reprojection(restored,p):
    with np.load(SOURCE/"initial_projection"/f"p{p}.npz",allow_pickle=False) as saved:
        np.testing.assert_array_equal(restored["q0"][p],saved["q"])
        np.testing.assert_array_equal(restored["v0"][p],saved["velocity"])
    assert restored["q0"][p].shape==restored["v0"][p].shape==(4*(p-1),)
    assert restored["q0"][p].dtype==np.float64
    assert np.all(np.isfinite(restored["q0"][p]))
    np.testing.assert_array_equal(restored["v0"][p],0.)
    assert np.linalg.norm(restored["q0"][p][:p-1])>0.
    assert np.linalg.norm(restored["q0"][p][3*(p-1):])>0.


def test_source_loader_keeps_strict_failure_and_exploratory_qualification(restored):
    summary=restored["summary"]
    assert summary["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert summary["state_admitted_flag"] is False
    assert summary["statuses"]["NLSP_STRICT_INITIAL_VERIFICATION"]=="PARTIAL"
    assert any(not row["pass"] for row in summary["strict_table"])
    assert summary["new_ODE_integrations"]==3  # inherited source, not these tests
    assert summary["new_eigendecompositions"]==summary["new_BVP_solves"]==0


@pytest.mark.parametrize("name",("p48_tight","p64_tight","p64_allowed_extra"))
def test_restored_radau_prescription_is_exact_saved_source_vector(restored,name):
    expected=read(SOURCE/"short_controls"/name/"case.json")
    assert restored["settings"][name]==expected
    p=expected["p"]
    assert len(expected["atol"])==8*(p-1)
    assert all(value>0. for value in expected["atol"])
    assert expected["rtol"]==(2e-11 if name.endswith("allowed_extra") else 1e-10)
    assert expected["max_step"]>0.


def test_one_period_configuration_keeps_exact_three_cases_and_fixed_physical_state(config):
    assert config["schema"]=="nlsp-prepared-one-T1-v1"
    assert config["source_bundle"]==str(SOURCE.relative_to(ROOT)).replace(chr(92),"/")
    assert config["periods"]==1.
    assert config["degrees"]==[48,64]
    assert config["cases"]==[[48,"tight"],[64,"tight"],[64,"allowed_extra"]]
    assert config["amplitude_over_h"]==.05
    assert config["projection_policy"]==prep.CONSTRAINED_PROJECTION
    assert config["state_admitted"] is False
    assert config["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert config["budget"]=={"numerical_wall_seconds":1200,"maximum_integrations":3}
    assert config["storage"]["policy"]=="one_memmapped_float64_state_per_case"
    assert config["storage"]["block_rows"]==256


@pytest.mark.parametrize("key",("recompute_projection","derivative_BC","phase_alignment","period_fit",
                              "energy_matching","energy_classification","out_of_plane_stability","continuum_truth"))
def test_one_period_does_not_authorize_other_physics_or_alignment(config,key):
    assert config["semantics"][key] is False


def test_old_feasibility_default_horizon_and_prescriptions_are_not_promoted_to_one_period():
    current=read(ROOT/"data/input/planar_prepared_feasibility.json")
    source=read(SOURCE/"summary.json")["config"]
    assert current==source
    assert current["short_periods"]==.1
    assert current["cases"]==[[48,"tight"],[64,"tight"],[64,"allowed_extra"]]


def test_one_period_grid_preserves_real_short_timestamps_and_frequency_bound(restored):
    end=restored["summary"]["background"]["T1"]
    times,sampling=cli.one_T1_time_grid(restored,end)
    assert times.dtype==np.float64
    assert times[0]==0. and times[-1]==end
    assert np.all(np.diff(times)>0.)
    assert np.all(np.isin(restored["short_time"],times))
    omega=restored["summary"]["sampling"]["omega_upper_bound"]
    bound_spacing=2*np.pi/(16*omega)
    assert np.max(np.diff(times))<=bound_spacing*(1+2e-12)
    assert len(times)>len(restored["short_time"])
    for fraction in (0.,.1,.25,.5,.75,1.):
        assert fraction*end in times
    assert isinstance(sampling,dict)


def test_new_model_and_rhs_files_match_frozen_source_execution_hashes(restored):
    protected=("scripts/lib/weakly_nonlinear_spatial_rod.py",
               "scripts/lib/weakly_nonlinear_planar_dynamics.py",
               "scripts/lib/planar_prepared_initial_state.py")
    expected=restored["manifest"]["identity"]["code_hashes"]
    for path in protected:
        assert sha(ROOT/path)==expected[path]


def test_runner_storage_extension_preserves_entire_previous_numerical_path():
    # Read immutable Git object only; neither checkout nor a legacy integration.
    before=subprocess.run(["git","show",BASELINE_HEAD+
        ":scripts/analysis/simulate_weakly_nonlinear_planar_rod.py"],
        cwd=ROOT,check=True,stdout=subprocess.PIPE).stdout.decode("utf8")
    old=ast.parse(before)
    current=ast.parse((ROOT/"scripts/analysis/simulate_weakly_nonlinear_planar_rod.py").read_text(encoding="utf8"))
    old_functions={x.name:x for x in old.body if isinstance(x,ast.FunctionDef)}
    new_functions={x.name:x for x in current.body if isinstance(x,ast.FunctionDef)}
    assert old_functions.keys()==new_functions.keys()
    for name in old_functions.keys()-{"integrate_case"}:
        assert ast.dump(old_functions[name])==ast.dump(new_functions[name]),name
    original,new=old_functions["integrate_case"],new_functions["integrate_case"]
    assert [x.arg for x in new.args.kwonlyargs]==["initial_coordinates","history_buffer"]
    assert all(ast.dump(x)==ast.dump(ast.Constant(None)) for x in new.args.kw_defaults)
    branch=next(x for x in new.body if isinstance(x,ast.If)
                and ast.dump(x.test)==ast.dump(ast.parse("history_buffer is None").body[0].value))
    assignment=next(x for x in original.body if isinstance(x,ast.Assign)
                    and any(isinstance(t,ast.Name) and t.id=="history" for t in x.targets))
    assert ast.dump(branch.body[0])==ast.dump(assignment)
    index=new.body.index(branch);new.body[index:index+1]=branch.body
    new.args=original.args
    assert ast.dump(new)==ast.dump(original)


def test_test_module_never_invokes_integrators_precision_or_eigenfunctions():
    tree=ast.parse(Path(__file__).read_text(encoding="utf8"))
    forbidden={"integrate_case","solve_ivp","Radau","eigh","eig","eigvals",
               "stable_initial_projection","high_precision_source_jets","project_saved_legendre","workdps"}
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.id if isinstance(node.func,ast.Name) else (
                node.func.attr if isinstance(node.func,ast.Attribute) else None)
            assert name not in forbidden



@pytest.mark.parametrize("key,value",(
    ("degrees",[64]),("cases",[[48,"tight"]]*4),("periods",5.),
    ("amplitude_over_h",.025),("projection_policy",prep.UNCONSTRAINED_PROJECTION),
    ("state_admitted",True),("execution_mode","STRICT_ADMITTED"),
    ("budget",{"numerical_wall_seconds":1201,"maximum_integrations":3}),
    ("budget",{"numerical_wall_seconds":1200,"maximum_integrations":4}),
))
def test_source_authorization_rejects_any_unrequested_extension_before_loading(
        config,monkeypatch,key,value):
    altered=json.loads(json.dumps(config));altered[key]=value
    def forbidden(*args,**kwargs):
        raise AssertionError("An unrequested extension must stop before source loading")
    monkeypatch.setattr(cli,"validate_cache",forbidden)
    with pytest.raises(ValueError):
        cli.load_one_T1_source(altered)


@pytest.mark.parametrize("relative",("initial_projection/p48.npz","summary.json"))
def test_source_loader_refuses_changed_own_artifact_hash(config,monkeypatch,relative):
    original=cli.sha
    target=SOURCE/relative
    def changed(path):
        return "0"*64 if Path(path).resolve()==target.resolve() else original(path)
    monkeypatch.setattr(cli,"sha",changed)
    with pytest.raises(ValueError,match="artifact hash mismatch"):
        cli.load_one_T1_source(config)


def test_source_loader_refuses_changed_physics_or_time_input(config,monkeypatch):
    original=cli.sha
    target=ROOT/config["pilot_config"]
    monkeypatch.setattr(cli,"sha",lambda path:
                        "0"*64 if Path(path).resolve()==target.resolve() else original(path))
    with pytest.raises(ValueError,match="Frozen physics/time input mismatch"):
        cli.load_one_T1_source(config)


@pytest.mark.parametrize("changed",("horizon","source","sampling","budget"))
def test_cache_identity_includes_source_horizon_and_declared_numerical_contract(
        config,tmp_path,changed):
    original=tmp_path/"config.json"
    cli.write_json(original,config)
    before,item=cli.one_T1_identity(original)
    assert item["source"]["manifest_sha256"]==sha(SOURCE/"manifest.json")
    assert item["source"]["artifact_hashes"]==read(SOURCE/"manifest.json")["artifact_hashes"]
    altered=json.loads(json.dumps(config))
    if changed=="horizon":altered["periods"]=.5
    elif changed=="sampling":altered["sampling"]["samples_per_upper_bound_period"]=24
    elif changed=="budget":altered["budget"]["numerical_wall_seconds"]=1100
    else:
        different=tmp_path/"source"
        different.mkdir()
        manifest=read(SOURCE/"manifest.json")
        manifest["synthetic_identity_change"]=True
        cli.write_json(different/"manifest.json",manifest)
        altered["source_bundle"]=str(different)
    cli.write_json(original,altered)
    after,_=cli.one_T1_identity(original)
    assert before!=after  # identity test only; changed config is never executed


def test_mmap_history_loader_reads_actual_prefix_and_never_unwritten_future_rows(tmp_path):
    folder=tmp_path/"cases"/"partial"
    folder.mkdir(parents=True)
    times=np.array([0.,.1,.2])
    state=np.full((8,8),np.nan)
    state[:3]=np.arange(24).reshape(3,8)
    np.save(folder/"time.npy",times)
    np.save(folder/"state.npy",state)
    case={"actual_valid_rows":3,"storage_allocated_rows":8,"ndof":4,"time_end":.2,"status":"PARTIAL"}
    cli.write_json(folder/"case.json",case)
    result=cli.one_T1_history(tmp_path,"partial")
    assert result["q"].shape==result["velocity"].shape==(3,4)
    np.testing.assert_array_equal(result["time"],times)
    np.testing.assert_array_equal(result["q"],state[:3,:4])
    np.testing.assert_array_equal(result["velocity"],state[:3,4:])
    assert np.all(np.isfinite(result["q"])) and np.all(np.isfinite(result["velocity"]))
    assert result["case"]["status"]=="PARTIAL"
    assert isinstance(result["q"],np.memmap)


@pytest.mark.parametrize("corruption",("rows","end","width"))
def test_mmap_history_loader_rejects_inconsistent_real_prefix_metadata(tmp_path,corruption):
    folder=tmp_path/"cases"/"partial"
    folder.mkdir(parents=True)
    times=np.array([0.,.1,.2])
    state=np.zeros((8,8))
    case={"actual_valid_rows":3,"storage_allocated_rows":8,"ndof":4,"time_end":.2}
    if corruption=="rows":case["actual_valid_rows"]=4
    elif corruption=="end":case["time_end"]=.3
    else:state=state[:,:7]
    np.save(folder/"time.npy",times);np.save(folder/"state.npy",state)
    cli.write_json(folder/"case.json",case)
    with pytest.raises(ValueError,match="Actual prefix metadata mismatch"):
        cli.one_T1_history(tmp_path,"partial")


def synthetic_cache(folder,item):
    folder.mkdir(parents=True)
    result={"schema":"nlsp-prepared-one-T1-v1","synthetic_fixture":True,
            "execution_mode":"EXPLORATORY_NOT_CERTIFIED","state_admitted_flag":False,
            "statuses":{"NLSP_PREPARED_ONE_T1_FEASIBILITY":"PARTIAL"},
            "new_ODE_integrations":3,"new_eigendecompositions":0}
    cli.write_json(folder/"summary.json",result)
    cli.write_json(folder/"manifest.json",cli.manifest_for(folder,item))
    return result


@pytest.mark.parametrize("action",("compute","report-only","plot-only"))
def test_cached_one_period_modes_do_zero_extra_runs_or_preparation(
        tmp_path,monkeypatch,capsys,action):
    item={"synthetic_cache":"one-T1-no-evaluation"}
    bundle=tmp_path/"output"/"fixture"
    expected=synthetic_cache(bundle,item)
    def forbidden(*args,**kwargs):
        raise AssertionError("Cached one-period mode may only read saved data")
    for name in ("run_one_T1","run_feasibility","run_compute","load_one_T1_source",
                 "one_T1_time_grid","one_T1_compare","one_T1_case_diagnostics"):
        monkeypatch.setattr(cli,name,forbidden)
    plots=[]
    monkeypatch.setattr(cli,"plot_only",lambda path:plots.append(Path(path)))
    monkeypatch.setattr(cli,"one_T1_identity",lambda *args:("fixture",item))
    args=(["--compute","--one-T1","--output-dir",str(bundle.parent)]
          if action=="compute" else ["--"+action,str(bundle)])
    returned=cli.main(args)
    printed=json.loads(capsys.readouterr().out)
    assert returned==expected
    assert all(value==0 for value in printed["this_run_counters"].values())
    assert returned["new_ODE_integrations"]==3  # saved history, not current work
    assert plots==([bundle] if action=="plot-only" else [])


def test_mutually_exclusive_horizon_flags_cannot_create_an_extra_or_ambiguous_run():
    with pytest.raises(SystemExit) as captured:
        cli.main(["--compute","--one-T1","--feasibility"])
    assert captured.value.code==2



RESULT=ROOT/"results/planar_prepared_one_T1/795dcb14d3cd3a55"


@pytest.fixture(scope="module")
def completed_result():
    if not (RESULT/"manifest.json").is_file():
        pytest.skip("Actual one-T1 result absent; tests never recreate trajectories")
    summary=cli.validate_cache(RESULT)
    histories={name:cli.one_T1_history(RESULT,name) for name in summary["cases"]}
    return summary,histories


def test_actual_execution_is_exactly_three_exploratory_runs_with_zero_new_preparation(completed_result,restored):
    summary,_=completed_result
    assert set(summary["cases"])=={"p48_tight","p64_tight","p64_allowed_extra"}
    assert summary["new_ODE_integrations"]==summary["runtime"]["ODE_integrations"]==3
    assert summary["new_projection_MP_BVP_eigen_symbolic_calls"]==0
    assert summary["runtime"]["projection_MP_BVP_eigen_symbolic_calls"]==0
    assert summary["runtime"]["numerical_wall_seconds"]<=1200.
    assert summary["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert summary["state_admitted_flag"] is False
    assert summary["source_strict_status"]=="PARTIAL"
    assert summary["source_strict_table"]==restored["summary"]["strict_table"]
    assert summary["background"]==restored["summary"]["background"]
    assert summary["coefficients"]==restored["summary"]["coefficients"]
    assert summary["no_new_dynamic_constraints"] is True


@pytest.mark.parametrize("name",("p48_tight","p64_tight","p64_allowed_extra"))
def test_actual_full_time_arrays_settings_and_q0_match_declared_source(
        completed_result,restored,name):
    summary,histories=completed_result
    hist=histories[name];case=hist["case"];p=case["p"]
    t=hist["time"];end=summary["background"]["T1"]
    assert t[0]==0. and t[-1]==end
    assert np.all(np.diff(t)>0.)
    assert len(t)==case["samples"]==case["actual_valid_rows"]==case["storage_allocated_rows"]
    assert len(t)==summary["sampling"]["samples"]
    assert hist["q"].shape==hist["velocity"].shape==(len(t),4*(p-1))
    assert case["time_end"]==case["target_time_end"]==end
    assert case["status"]=="PASS" and case["periods"]==1.
    assert case["state_admitted"] is False
    assert case["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    np.testing.assert_array_equal(hist["q"][0],restored["q0"][p])
    np.testing.assert_array_equal(hist["velocity"][0],restored["v0"][p])
    old=restored["settings"][name]
    for key in ("rtol","max_step","atol","atol_coordinate_scales","velocity_scale_multiplier","nq","ndof"):
        assert case[key]==old[key]


def test_actual_three_histories_share_full_grid_preserve_source_and_snapshot_times(completed_result,restored):
    summary,histories=completed_result
    common=None;T1=summary["background"]["T1"]
    for name,hist in histories.items():
        if common is None:common=hist["time"]
        else:np.testing.assert_array_equal(hist["time"],common)
        assert np.all(np.isin(restored["short_time"],hist["time"]))
        with np.load(RESULT/"cases"/name/"snapshots.npz",allow_pickle=False) as data:
            np.testing.assert_array_equal(data["time"],T1*np.array([0.,.1,.25,.5,.75,1.]))
            assert data["fields"].shape==data["velocities"].shape==(6,501,4)
            np.testing.assert_array_equal(data["fields"][:,[0,-1],:],0.)
            np.testing.assert_array_equal(data["velocities"][:,[0,-1],:],0.)
    np.testing.assert_array_equal(histories["p64_tight"]["q"][0],histories["p64_allowed_extra"]["q"][0])


def test_actual_copy_of_continuous_target_is_byte_identical_and_not_refit():
    if not (RESULT/"manifest.json").is_file():
        pytest.skip("Actual result unavailable")
    assert sha(RESULT/"common_initial_state.npz")==sha(SOURCE/"common_initial_state.npz")


def test_actual_source_prefix_regressions_cover_entire_saved_point_one_period(completed_result,restored):
    summary,_=completed_result
    for name,record in summary["prefix_regression"].items():
        assert record["source_interval_complete"] is True
        assert record["initial_q_v_exact"] is True
        assert record["time_settings_identical"] is True
        assert record["time_end"]==restored["short_time"][-1]
        assert record["samples"]==len(restored["short_time"])
        assert record["horizon_complete"] is True
        assert record["required_end"]==restored["short_time"][-1]
        assert record["component_gates_pass"] is True
        assert record["pass"] is True and record["status"]=="PASS"
        assert len(record["fields"])==8
        assert record["relative_energy_difference_max"]<=restored["pilot"]["gates"]["energy_relative_drift"]
    assert summary["statuses"]["NLSP_PREPARED_ONE_T1_PREFIX_REGRESSION"]=="PASS"


@pytest.mark.parametrize("comparison",("spatial","temporal"))
def test_actual_comparison_all_eight_metrics_use_unchanged_full_horizon_reference(
        completed_result,restored,comparison):
    summary,_=completed_result
    report=summary[comparison+"_comparison"];T1=summary["background"]["T1"]
    assert report["time_end"]==T1 and report["samples"]==summary["sampling"]["samples"]
    assert report["required_end"]==T1 and report["horizon_complete"] is True
    assert set(report["fields"])=={part+"_"+f for part in ("q","velocity") for f in ("u","w","theta","c")}
    assert "full-comparison-horizon" in report["normalization"]
    assert "no phase alignment" in report["normalization"]
    assert "not continuum supremum" in report["maximum_semantics"]
    for name,row in report["fields"].items():
        field=name.split("_")[-1]
        expected=restored["pilot"]["gates"]["w_theta_relative" if field in ("w","theta") else "u_c_relative"]
        assert row["tolerance"]==expected
        assert row["relative_L2"]==row["max_time_L2_difference"]/max(row["reference_max_time_L2"],row["numerical_floor"])
        assert row["relative_max"]==row["max_space_time_difference"]/max(row["reference_max_space_time"],row["numerical_floor"])
        assert row["pass"]==(row["relative_L2"]<=expected and row["relative_max"]<=expected)
        assert 0.<=row["L2_peak_time_tau"]<=1.
        assert 0.<=row["max_peak_time_tau"]<=1.
        assert 0.<=row["max_peak_s_over_L"]<=1.
    assert report["status"]==("PASS" if all(r["pass"] for r in report["fields"].values()) else "PARTIAL")


@pytest.mark.parametrize("comparison",("spatial","temporal"))
def test_actual_saved_error_curves_and_cumulative_windows_are_full_scale_without_exclusions(
        completed_result,comparison):
    summary,_=completed_result
    report=summary[comparison+"_comparison"];T1=summary["background"]["T1"]
    with np.load(RESULT/(comparison+"_differences.npz"),allow_pickle=False) as data:
        times=data["time"];cuml2=data["cumulative_L2"];cummax=data["cumulative_max"]
        assert data["d_L2"].shape==data["d_max"].shape==(len(times),8)
        np.testing.assert_array_equal(cuml2,np.maximum.accumulate(data["d_L2"],axis=0))
        np.testing.assert_array_equal(cummax,np.maximum.accumulate(data["d_max"],axis=0))
        assert times[-1]==T1
        fields=[part+"_"+f for part in ("q","velocity") for f in ("u","w","theta","c")]
        for k,name in enumerate(fields):
            row=report["fields"][name]
            assert cuml2[-1,k]==row["max_time_L2_difference"]
            assert cummax[-1,k]==row["max_space_time_difference"]
        assert len(report["windows"])==8*5
        for row in report["windows"]:
            k=fields.index(row["component"])
            at=np.searchsorted(times,row["end_tau"]*T1,side="right")-1
            assert times[at]==row["end_tau"]*T1
            assert row["absolute_cumulative_L2"]==cuml2[at,k]
            assert row["absolute_cumulative_max"]==cummax[at,k]
            assert row["relative_L2_full_scale"]==cuml2[at,k]/data["normalization_L2"][k]
            assert row["relative_max_full_scale"]==cummax[at,k]/data["normalization_max"][k]


def test_actual_energy_safety_and_mass_bounds_pass_without_eigen_recomputation(completed_result,restored):
    summary,histories=completed_result
    for name,hist in histories.items():
        case=hist["case"]
        energy=np.load(RESULT/"cases"/name/"energy.npy",mmap_mode="r")
        assert len(energy)==len(hist["time"])
        assert np.all(np.isfinite(energy)) and np.all(energy>0.)
        assert np.max(abs((energy-energy[0])/energy[0]))==case["relative_energy_drift_max"]
        assert case["relative_energy_drift_max"]<=restored["pilot"]["gates"]["energy_relative_drift"]
        assert case["mass_lower_bound_min"]>=restored["pilot"]["safety"]["min_relative_mass_eigenvalue"]
        assert case["energy_and_mass_pass"] and case["safety_pass"]
        assert "no eigensolves" in case["mass_bound_method"]
    assert summary["statuses"]["NLSP_PREPARED_ONE_T1_ENERGY_AND_MASS"]=="PASS"


def test_actual_completed_execution_does_not_hide_spatial_theta_velocity_failure(completed_result):
    summary,_=completed_result
    spatial,temporal=summary["spatial_comparison"],summary["temporal_comparison"]
    assert temporal["status"]=="PASS"
    assert spatial["status"]=="PARTIAL"
    assert [name for name,row in spatial["fields"].items() if not row["pass"]]==["velocity_theta"]
    assert spatial["fields"]["velocity_theta"]["relative_max"]>1e-4
    assert spatial["fields"]["velocity_theta"]["relative_L2"]<=1e-4
    assert summary["statuses"]["NLSP_PREPARED_ONE_T1_SPATIAL_CHECK"]=="PARTIAL"
    assert summary["statuses"]["NLSP_PREPARED_ONE_T1_TEMPORAL_CHECK"]=="PASS"
    assert summary["statuses"]["NLSP_PREPARED_ONE_T1_FEASIBILITY"]=="COMPLETED_EXPLORATORY_NOT_CERTIFIED"
    assert "no periodic orbit or continuum truth claim" in summary["qualification"]


@pytest.mark.parametrize("action",("compute","report-only"))
def test_actual_cached_compute_and_report_have_zero_numerical_work(
        completed_result,monkeypatch,capsys,action):
    summary,_=completed_result
    def forbidden(*args,**kwargs):
        raise AssertionError("Actual completed one-T1 data must never be recomputed")
    for name in ("run_one_T1","load_one_T1_source","one_T1_time_grid","one_T1_compare","one_T1_case_diagnostics"):
        monkeypatch.setattr(cli,name,forbidden)
    item=read(RESULT/"manifest.json")["identity"]
    monkeypatch.setattr(cli,"one_T1_identity",lambda *args:(RESULT.name,item))
    args=(["--compute","--one-T1","--output-dir",str(RESULT.parent)]
          if action=="compute" else ["--report-only",str(RESULT)])
    result=cli.main(args);output=json.loads(capsys.readouterr().out)
    assert result==summary
    assert all(value==0 for value in output["this_run_counters"].values())


def test_partial_horizon_cannot_pass_full_T1_but_explicit_point_one_prefix_can(tmp_path):
    from types import SimpleNamespace
    # Software coverage fixture only; no Galerkin model or trajectory is solved.
    t=np.linspace(0.,.1,11)
    q=np.ones((len(t),4))*1e-5
    data={"time":t,"q":q,"velocity":np.zeros_like(q)}
    disc=SimpleNamespace(length=1.,reconstruct_series=lambda a,x:
                         a[:,None,:]*(x*(1-x))[None,:,None])
    pilot={"spatial":{"comparison_quadrature":12},"gates":{
        "relative_numerical_floor":1e-10,"w_theta_relative":1e-4,"u_c_relative":1e-3}}
    full=cli.one_T1_compare(data,data,disc,disc,pilot,1.,np.ones(8),tmp_path/"full.npz")
    assert full["component_gates_pass"] is True
    assert all(row["pass"] for row in full["fields"].values())
    assert full["time_end"]==.1 and full["required_end"]==1.
    assert full["horizon_complete"] is False and full["status"]=="PARTIAL"
    prefix=cli.one_T1_compare(data,data,disc,disc,pilot,1.,np.ones(8),
                              tmp_path/"prefix.npz",required_end=.1)
    assert prefix["time_end"]==prefix["required_end"]==.1
    assert prefix["horizon_complete"] and prefix["component_gates_pass"]
    assert prefix["status"]=="PASS"
    assert prefix["fields"]==full["fields"]  # only coverage changes, never norms
