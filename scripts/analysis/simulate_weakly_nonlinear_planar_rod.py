"""Bounded four-field free-motion pilot of the accepted quartic action.

This new action/time-history contract is separate from the frozen spatial
algebra audit. It is not a frequency map or nonlinear modal reduction.
Cached compute and plot/report-only never construct the symbolic model.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
VERSION = "nlsp-planar-time-pilot-v1"
CONFIG = ROOT / "data/input/weakly_nonlinear_planar_time_pilot.json"
OUTPUT = ROOT / "results/weakly_nonlinear_planar_time_pilot"
FIELDS = ("u", "w", "theta", "c")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf8")
    temp.replace(path)


def identity(config_path=CONFIG):
    config_path = Path(config_path)
    config = json.loads(config_path.read_text(encoding="utf8"))
    if config["schema"] != VERSION or config["fields"] != list(FIELDS):
        raise ValueError("Unsupported planar pilot contract")
    paths = [config_path, Path(__file__), ROOT/"scripts/lib/weakly_nonlinear_planar_dynamics.py",
             ROOT/"scripts/lib/weakly_nonlinear_spatial_rod.py",
             ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
             ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"]
    for key in ("audit_bundle", "linear_reference_bundle"):
        bundle = ROOT/config[key]
        manifest = json.loads((bundle/"manifest.json").read_text(encoding="utf8"))
        saved_hashes = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
        if not saved_hashes:
            raise ValueError(f"Reference has no recorded artifact hashes: {bundle}")
        for name, digest in saved_hashes.items():
            if sha(bundle/name) != digest:
                raise ValueError(f"Frozen reference artifact changed: {bundle/name}")
        paths += [bundle/"manifest.json", bundle/"result.json"]
    item = {"version": VERSION, "config": config,
            "hashes": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p) for p in paths},
            "python": sys.version, "dependencies": {n: importlib.metadata.version(n) for n in ("numpy", "scipy", "matplotlib")}}
    key = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()[:16]
    return key, item


def validate_bundle(bundle, expected=None):
    bundle = Path(bundle)
    manifest = json.loads((bundle/"manifest.json").read_text(encoding="utf8"))
    if expected is not None and manifest["identity"] != expected:
        raise ValueError("Cache identity differs from current inputs/model/environment")
    for name, digest in manifest["artifact_hashes"].items():
        if sha(bundle/name) != digest:
            raise ValueError(f"Cached artifact hash mismatch: {name}")
    return json.loads((bundle/"summary.json").read_text(encoding="utf8"))


def load_runtime():
    global np, eigh, solve_ivp, Radau, dynamics, rod, mh, rectangular_section
    import numpy as np
    from scipy.linalg import eigh
    from scipy.integrate import solve_ivp, Radau
    from scripts.lib import weakly_nonlinear_planar_dynamics as dynamics
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.lib import mindlin_herrmann_longitudinal as mh
    from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section


def setup(config):
    accepted = json.loads((ROOT/config["audit_bundle"]/"result.json").read_text(encoding="utf8"))
    coefficients = rod.RodCoefficients(**accepted["coefficients"])
    reference = json.loads((ROOT/config["linear_reference_bundle"]/"result.json").read_text(encoding="utf8"))
    omega = reference["timoshenko"]["roots"][0]["omega"]
    g = config["material_geometry"]
    section = rectangular_section(E=g["E"], rho=g["rho"], nu=g["nu"], width=g["b"], thickness=g["h"], K=5/6)
    source = mh.project_jang_reduced_rectangular(section)
    mode = mh.finite_mode(source, g["L"], omega, "timoshenko")
    peak = float((mh.finite_state_basis(source, g["L"], omega, [g["L"]/2], "timoshenko")@mode["coefficients"])[0,0])
    normalized = mode["coefficients"]/peak

    def shape(points, derivative=0):
        values = mh.finite_state_basis(source, g["L"], omega, np.asarray(points), "timoshenko", derivative)@normalized
        fields = np.zeros((len(values), 4))
        fields[:,1:3] = values[:,:2]
        return fields

    grid = np.linspace(0, g["L"], 501)
    if abs(np.max(np.abs(shape(grid)[:,1]))-1) > 2e-12 or abs(shape([g["L"]/2],1)[0,1]) > 2e-10:
        raise ArithmeticError("Continuous first-mode common normalization failed")
    info = {"omega": omega, "T1": 2*math.pi/omega, "peak_location": g["L"]/2,
            "mass_normalized_signed_w_peak": peak, "analytic_coefficients": normalized.tolist(),
            "normalization": "ALL analytic state coefficients divided by signed w(L/2), validated stationary/global maximum",
            "source_bundle": config["linear_reference_bundle"], "root_solves": 0}
    return coefficients, shape, info, reference


def make_discretization(coefficients, p, config, model=None):
    return dynamics.PlanarGalerkin(coefficients, p, length=config["material_geometry"]["L"], model=model)


def time_settings(discretization, amplitude, level, config):
    settings = config["time_levels"][level]
    p, length, n = discretization.coefficients, discretization.length, discretization.n
    cutoff = math.sqrt(p.C/p.jp)
    field_scales = np.array([amplitude, amplitude, amplitude/length, amplitude/length])
    masses = np.array([p.m, p.m, p.jp, p.jp])
    coordinate_scales = field_scales*np.sqrt(masses*length/n)
    atol_q = np.repeat(coordinate_scales, n)*settings["atol_relative"]
    atol = np.concatenate((atol_q, atol_q*cutoff))
    return {"rtol": settings["rtol"], "atol": atol,
            "max_step": 2*math.pi/cutoff*settings["max_step_cutoff_period_fraction"],
            "coordinate_scales": coordinate_scales.tolist(), "velocity_scale_multiplier": cutoff,
            "atol_relative": settings["atol_relative"]}


def preflight(config, coefficients, shape, initial, reference):
    started = time.perf_counter()
    rows = []
    symbolic = rod.derive_polynomials()
    for degree in config["spatial"]["degrees"]:
        disc = make_discretization(coefficients, degree, config, symbolic)
        eigen = {block: disc.linear_eigenpairs(block)["omega"][:3] for block in ("mh", "timoshenko")}
        old = {block: np.array([row["omega"] for row in reference[block]["roots"][:3]]) for block in eigen}
        random = np.random.default_rng(410+degree)
        scales = np.repeat([1e-6,1e-6,1e-8,1e-8], disc.n)
        q, v = random.normal(size=disc.ndof)*scales, random.normal(size=disc.ndof)*scales
        acceleration = disc.acceleration(q,v)
        weak = disc.weak_residual(q,v,acceleration)
        action = disc.mass_matrix(q)@acceleration+disc.inertial_terms(q,v)+disc.potential(q)["gradient"]
        power_scale = np.linalg.norm(v)*(np.linalg.norm(disc.potential(q)["gradient"])+np.linalg.norm(disc.mass_matrix(q)@acceleration))
        q0 = disc.project(shape(disc.x))
        projection = disc.reconstruct(q0)-shape(disc.x)
        matrices = disc.mass_matrix(q)
        high = dynamics.PlanarGalerkin(coefficients, degree, disc.length, nq=3*degree+7, model=symbolic)
        qh = high.from_raw_coefficients(disc.raw_coefficients(q))
        vh = high.from_raw_coefficients(disc.raw_coefficients(v))
        rows.append({"p":degree, "ndof":disc.ndof, "nq":disc.nq,
            "omega": {key:value.tolist() for key,value in eigen.items()},
            "linear_relative_errors": {key:(np.abs(eigen[key]/old[key]-1)).tolist() for key in eigen},
            "projection_L2": np.sqrt(disc.weights@(projection**2)).tolist(),
            "weak_action_max_abs":float(np.max(np.abs(weak-action))),
            "energy_rhs_relative":abs(disc.energy_rate(q,v,acceleration))/max(power_scale,1e-30),
            "M0_min_eigenvalue":float(np.linalg.eigvalsh(disc.M0)[0]),
            "mass_symmetry_abs":float(np.max(np.abs(matrices-matrices.T))),
            "quadrature_energy_difference":abs(disc.energy(q,v)-high.energy(qh,vh)),
            "highest_retained_omega":float(disc.linear_eigenpairs()["omega"][-1])})
    final = make_discretization(coefficients, config["spatial"]["degrees"][-1], config, symbolic)
    q0=final.project(.0025*shape(final.x)); v0=np.zeros(final.ndof)
    settings=time_settings(final,.0025,"tight",config)
    ndof=final.ndof
    linear_matrix=np.block([[np.zeros((ndof,ndof)),np.eye(ndof)],[-np.linalg.solve(final.M0,final.K),np.zeros((ndof,ndof))]])
    times=np.linspace(0, initial["T1"]/10,101)
    controls=[]
    for name, coordinates in (("first_bending",q0),("pure_axial",.00001*final.linear_eigenpairs("mh")["vectors"][:,0])):
        solved=solve_ivp(lambda t,y: linear_matrix@y,(0,times[-1]),np.r_[coordinates,v0],method="Radau",jac=linear_matrix,
                         rtol=settings["rtol"],atol=settings["atol"],max_step=settings["max_step"],t_eval=times)
        exact=final.linear_reference(coordinates,v0,times)
        controls.append({"name":name,"success":bool(solved.success),"time_end":float(times[-1]),
            "relative_q_error":float(np.max(np.abs(solved.y[:ndof].T-exact["q"]))/max(np.max(np.abs(exact["q"])),1e-30)),
            "relative_velocity_error":float(np.max(np.abs(solved.y[ndof:].T-exact["velocity"]))/max(np.max(np.abs(exact["velocity"])),1e-30)),
            "nfev":solved.nfev,"njev":solved.njev,"nlu":solved.nlu})
    passed=(max(rows[-1]["linear_relative_errors"]["mh"])<config["gates"]["linear_final_frequency_relative"] and
            all(row["energy_rhs_relative"]<config["gates"]["identity_scaled"] and row["weak_action_max_abs"]<2e-12 and row["M0_min_eigenvalue"]>0 for row in rows) and
            all(row["success"] and max(row["relative_q_error"],row["relative_velocity_error"])<1e-7 for row in controls))
    return {"status":"PASS" if passed else "FAIL","spatial_linear_controls":rows,"linear_time_controls":controls,
            "initial_eigenpair":initial,"zero_rhs_max":float(np.max(np.abs(final.rhs(0,np.zeros(2*ndof))))),
            "runtime_seconds":time.perf_counter()-started,"symbolic_derivations":1,"root_solves":0}


def safety_check(disc, q, policy):
    values=disc.reconstruct(q)
    gradients=disc.reconstruct(q,derivative=1)
    conditions={"min_one_plus_c":float(np.min(1+values[:,3])),"max_abs_c":float(np.max(np.abs(values[:,3]))),
                "max_abs_theta":float(np.max(np.abs(values[:,2]))),"max_abs_axial_gradient":float(np.max(np.abs(gradients[:,0]))),
                "max_abs_transverse_gradient":float(np.max(np.abs(gradients[:,1]))),
                "max_L_abs_curvature":float(disc.length*np.max(np.abs(gradients[:,2])))}
    if conditions["min_one_plus_c"]<=policy["min_one_plus_c"] or any(conditions[key]>policy[key] for key in conditions if key!="min_one_plus_c"):
        raise ArithmeticError(f"Declared small-neighborhood safety gate: {conditions}")


def integrate_case(disc, shape, initial, config, amplitude_ratio, level, times, deadline, *, initial_coordinates=None):
    amplitude=amplitude_ratio*config["material_geometry"]["h"]
    if initial_coordinates is None:
        q0=disc.project(amplitude*shape(disc.x))
    else:
        q0=np.asarray(initial_coordinates,dtype=float)
        if q0.shape!=(disc.ndof,) or not np.all(np.isfinite(q0)):
            raise ValueError("Explicit initial coordinates must be a finite full-sized vector")
        q0=q0.copy()
    v0=np.zeros(disc.ndof)
    settings=time_settings(disc,amplitude,level,config)
    history=np.empty((len(times),2*disc.ndof));history[0]=np.r_[q0,v0]
    cursor=1
    started=time.perf_counter(); disc.reset_counters()
    def rhs(t,y):
        if not np.all(np.isfinite(y)):
            raise ArithmeticError("NONFINITE_STATE")
        safety_check(disc,y[:disc.ndof],config["safety"])
        return disc.rhs(t,y)
    solver=Radau(rhs,0,np.r_[q0,v0],float(times[-1]),jac=disc.jacobian,
                 rtol=settings["rtol"],atol=settings["atol"],max_step=settings["max_step"])
    steps=[]; failure=None
    while solver.status=="running":
        if time.perf_counter()>deadline:
            failure="PREDECLARED_COMPUTATIONAL_BUDGET_EXHAUSTED"
            break
        old=solver.t
        try:
            message=solver.step()
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
            failure=type(error).__name__+": "+str(error)
            break
        if solver.status=="failed":
            failure=message or "RADAU_FAILED";break
        steps.append(solver.t-old)
        end=int(np.searchsorted(times,solver.t,side="right"))
        if end>cursor:
            history[cursor:end]=solver.dense_output()(times[cursor:end]).T
            cursor=end
    elapsed=time.perf_counter()-started
    stats={"status":"PASS" if cursor==len(times) and failure is None else "PARTIAL",
           "failure":failure,"p":disc.p,"ndof":disc.ndof,"nq":disc.nq,
           "amplitude_over_h":amplitude_ratio,"amplitude":amplitude,"time_level":level,
           "rtol":settings["rtol"],"atol":settings["atol"].tolist(),"max_step":settings["max_step"],
           "atol_coordinate_scales":settings["coordinate_scales"],"velocity_scale_multiplier":settings["velocity_scale_multiplier"],
           "time_end":float(times[cursor-1]),"target_time_end":float(times[-1]),"samples":cursor,
           "accepted_internal_steps":len(steps),"internal_time_steps":steps,"min_internal_step":min(steps,default=0),"max_internal_step":max(steps,default=0),
           "nfev":solver.nfev,"njev":solver.njev,"nlu":solver.nlu,"counters":disc.counters(),"integration_seconds":elapsed}
    return history[:cursor],stats


def series_measure(disc, q, v, config):
    maxima=np.zeros(4); velocity_maxima=np.zeros(4); energy=[]
    norm_rows=[]; velocity_norm_rows=[]
    mass_lower=1.;mass_upper=1.; min_scale=1.; max_theta=0.;max_c=0.;max_us=0.;max_ws=0.;max_curvature=0.;max_g1=0.;max_g2=0.
    explicit=[]
    for index in range(len(q)):
        values=disc.reconstruct(q[index]); speeds=disc.reconstruct(v[index]);grad=disc.reconstruct(q[index],derivative=1)
        maxima=np.maximum(maxima,np.max(np.abs(values),axis=0));velocity_maxima=np.maximum(velocity_maxima,np.max(np.abs(speeds),axis=0))
        norm_rows.append(np.sqrt(disc.weights@(values**2)))
        velocity_norm_rows.append(np.sqrt(disc.weights@(speeds**2)))
        energy.append(disc.energy(q[index],v[index]))
        scale=(1+values[:,3])**2
        mass_lower=min(mass_lower,float(np.min(scale)));mass_upper=max(mass_upper,float(np.max(scale)))
        min_scale=min(min_scale,float(np.min(1+values[:,3])))
        max_c=max(max_c,float(np.max(np.abs(values[:,3]))));max_theta=max(max_theta,float(np.max(np.abs(values[:,2]))))
        max_us=max(max_us,float(np.max(np.abs(grad[:,0]))));max_ws=max(max_ws,float(np.max(np.abs(grad[:,1]))))
        max_curvature=max(max_curvature,disc.length*float(np.max(np.abs(grad[:,2]))))
        th=values[:,2];us=grad[:,0];ws=grad[:,1]
        gamma1=us+ws*th-th**2/2-us*th**2/2-ws*th**3/6+th**4/24
        gamma2=ws-th-us*th-ws*th**2/2+th**3/6+us*th**3/6
        max_g1=max(max_g1,float(np.max(np.abs(gamma1))));max_g2=max(max_g2,float(np.max(np.abs(gamma2))))
        if index in (0,len(q)//4,len(q)//2,3*len(q)//4,len(q)-1):
            explicit.append({"sample":index,**disc.diagnostics(q[index])})
    energy=np.asarray(energy);drift=(energy-energy[0])/energy[0]
    return {"max_abs_fields":dict(zip(FIELDS,maxima.tolist())),"max_abs_velocities":dict(zip(FIELDS,velocity_maxima.tolist())),
            "energy_initial":float(energy[0]),"max_relative_energy_drift":float(np.max(np.abs(drift))),
            "relative_mass_eigenvalue_lower_bound":mass_lower,"relative_mass_eigenvalue_upper_bound":mass_upper,
            "relative_mass_condition_upper_bound":mass_upper/mass_lower,"mass_bound_method":"Loewner bounds of weighted theta Gram; sparse explicit generalized eigensolves",
            "min_one_plus_c":min_scale,"max_abs_c":max_c,"max_abs_theta":max_theta,"max_abs_u_s":max_us,"max_abs_w_s":max_ws,
            "max_L_abs_theta_s":max_curvature,"max_quartic_Gamma1":max_g1,"max_quartic_Gamma2":max_g2,"explicit_mass_spots":explicit},energy,drift,np.asarray(norm_rows),np.asarray(velocity_norm_rows)


def compare_histories(disc_a, data_a, disc_b, data_b, config):
    nodes,weights=np.polynomial.legendre.leggauss(config["spatial"]["comparison_quadrature"])
    points=(nodes+1)*disc_a.length/2;weights*=disc_a.length/2
    rows={}
    for part in ("q","velocity"):
        a,b=data_a[part],data_b[part]
        if len(a)!=len(b):
            return {"status":"PARTIAL","reason":"Unequal/unfinished time horizons"}
        abs_l2=np.zeros(4);abs_max=np.zeros(4);scale_l2=np.zeros(4);scale_max=np.zeros(4)
        for start in range(0,len(a),256):
            aa=disc_a.reconstruct_series(a[start:start+256],points)
            bb=disc_b.reconstruct_series(b[start:start+256],points)
            difference=aa-bb
            abs_l2=np.maximum(abs_l2,np.max(np.sqrt(np.einsum("tif,i,tif->tf",difference,weights,difference)),axis=0))
            abs_max=np.maximum(abs_max,np.max(np.abs(difference),axis=(0,1)))
            scale_l2=np.maximum(scale_l2,np.max(np.sqrt(np.einsum("tif,i,tif->tf",bb,weights,bb)),axis=0))
            scale_max=np.maximum(scale_max,np.max(np.abs(bb),axis=(0,1)))
        for i,field in enumerate(FIELDS):
            threshold=config["gates"]["w_theta_relative" if field in ("w","theta") else "u_c_relative"]
            floor=config["gates"]["relative_numerical_floor"]*max(np.max(scale_l2),1e-30)
            relative=float(abs_l2[i]/max(scale_l2[i],floor))
            rows[part+"_"+field]={"max_time_L2_difference":float(abs_l2[i]),"max_space_time_difference":float(abs_max[i]),
                                  "reference_max_time_L2":float(scale_l2[i]),"reference_max_space_time":float(scale_max[i]),
                                  "relative_L2":relative,"relative_max":float(abs_max[i]/max(scale_max[i],floor)),
                                  "numerical_floor":float(floor),"floor_limited":bool(scale_l2[i]<=floor),"tolerance":threshold,
                                  "pass":bool(relative<=threshold and abs_max[i]/max(scale_max[i],floor)<=threshold)}
    return {"status":"PASS" if all(row["pass"] for row in rows.values()) else "PARTIAL","fields":rows}


def case_name(p, amplitude_ratio, level):
    return f"p{p}_Aoverh{str(amplitude_ratio).replace('.','p')}_{level}"


def save_case(bundle, name, disc, history, stats, times, shape, initial, config):
    folder=bundle/"cases"/name;folder.mkdir(parents=True,exist_ok=True)
    ndof=disc.ndof;q=history[:,:ndof];velocity=history[:,ndof:];times=times[:len(q)]
    started=time.perf_counter()
    diagnostic,energy,drift,field_norms,velocity_norms=series_measure(disc,q,velocity,config)
    observations=disc.reconstruct_series(q,np.array([disc.length/4,disc.length/2]))
    snap_ids=[int(np.argmin(abs(times-target*initial["T1"]))) for target in config["sampling"]["snapshots_periods"]]
    snapshot_points=np.linspace(0,disc.length,201)
    raw_q=np.array([disc.raw_coefficients(row) for row in q[snap_ids]])
    linear=disc.linear_reference(q[0],velocity[0],times)
    compared=compare_histories(disc,{"q":q,"velocity":velocity},disc,linear,config)
    continuous_errors={}
    sample_shape=shape(disc.x); reference_scale=np.sqrt(disc.weights@(sample_shape**2))*stats["amplitude"]
    for part in ("q","velocity"):
        errors=np.zeros(4)
        for start in range(0,len(times),256):
            physical=disc.reconstruct_series(q[start:start+256] if part=="q" else velocity[start:start+256])
            factor=np.cos(initial["omega"]*times[start:start+256]) if part=="q" else -initial["omega"]*np.sin(initial["omega"]*times[start:start+256])
            difference=physical-stats["amplitude"]*factor[:,None,None]*sample_shape[None,:,:]
            errors=np.maximum(errors,np.max(np.sqrt(np.einsum("tif,i,tif->tf",difference,disc.weights,difference)),axis=0))
        continuous_errors[part]={field:{"max_time_L2_difference":float(errors[i]),"difference_over_A":float(errors[i]/stats["amplitude"]),
                    "relative_to_linear_characteristic_L2":float(errors[i]/(reference_scale[i]*(initial["omega"] if part=="velocity" else 1))) if reference_scale[i]>0 else None} for i,field in enumerate(FIELDS)}
    weak_spots=[]
    for index in snap_ids:
        acceleration=disc.acceleration(q[index],velocity[index])
        force=disc.potential(q[index])["gradient"]
        scale=np.linalg.norm(force)+np.linalg.norm(disc.mass_matrix(q[index])@acceleration)+np.linalg.norm(disc.inertial_terms(q[index],velocity[index]))
        weak_spots.append({"sample":index,"scaled_L2_coefficient_weak_residual":float(np.linalg.norm(disc.weak_residual(q[index],velocity[index],acceleration))/max(scale,1e-30))})
    stats.update({"diagnostics":diagnostic,"linear_semidiscrete_difference":compared,
                  "continuous_linear_difference":continuous_errors,"weak_equation_spots":weak_spots,
                  "postprocessing_seconds":time.perf_counter()-started})
    if diagnostic["max_relative_energy_drift"]>config["gates"]["energy_relative_drift"]:
        stats["quality_status"]="FAIL"
    else:
        stats["quality_status"]="PASS" if stats["status"]=="PASS" else "PARTIAL"
    target=folder/"trajectory.npz";temp=folder/"trajectory.tmp.npz"
    np.savez_compressed(temp,time=times,q=q,velocity=velocity,observations=observations,energy=energy,energy_drift=drift,
                        field_L2_norms=field_norms,velocity_L2_norms=velocity_norms,
                        snapshot_indices=snap_ids,snapshot_points=snapshot_points,
                        snapshots=disc.reconstruct_series(q[snap_ids],snapshot_points),
                        snapshot_velocities=disc.reconstruct_series(velocity[snap_ids],snapshot_points),
                        raw_snapshot_coefficients=raw_q,initial_projection=q[0],
                        initial_shape=shape(snapshot_points),analytic_linear_q=shape(snapshot_points)[None,:,:]*stats["amplitude"]*np.cos(initial["omega"]*times[snap_ids])[:,None,None])
    temp.replace(target)
    stats["artifact_hashes"]={"trajectory.npz":sha(target)}
    write_json(folder/"case.json",stats)
    return stats


def run_compute(bundle, config, item):
    if config["budget"]["status"]!="FIXED_AFTER_SMOKE":
        raise ValueError("Run --smoke; record budget in config before --compute")
    started=time.perf_counter(); deadline=started+config["budget"]["total_wall_seconds"]
    coefficients,shape,initial,reference=setup(config)
    pre=preflight(config,coefficients,shape,initial,reference)
    write_json(bundle/"preflight.json",pre)
    if pre["status"]!="PASS":
        raise ArithmeticError("Planar discretization/linear hard gate failed")
    symbolic=rod.derive_polynomials()
    discs={p:make_discretization(coefficients,p,config,symbolic) for p in config["spatial"]["degrees"]}
    final_p=config["spatial"]["degrees"][-1]
    fastest=float(discs[final_p].linear_eigenpairs()["omega"][-1])
    output_step=2*math.pi/fastest/config["sampling"]["steps_per_fastest_retained_linear_period"]
    end=config["periods"]*initial["T1"]
    times=np.linspace(0,end,int(math.ceil(end/output_step))+1)
    plan=[(p,.05,"tight") for p in config["spatial"]["degrees"]]+[(final_p,.05,"coarse"),(final_p,.05,"medium"),(final_p,.025,"tight"),(config["spatial"]["degrees"][-2],.025,"tight")]
    completed={}; failures=[]; integrations=0
    write_json(bundle/"plan.json",{"plan":plan,"budget":config["budget"],"output_samples":len(times),"output_dt":float(times[1]),"highest_linear_omega":fastest,"locked_before_main":True})
    for p,ratio,level in plan:
        name=case_name(p,ratio,level);folder=bundle/"cases"/name
        if (folder/"case.json").exists():
            cached=json.loads((folder/"case.json").read_text(encoding="utf8"))
            if cached["status"]=="PASS" and sha(folder/"trajectory.npz")==cached["artifact_hashes"]["trajectory.npz"]:
                completed[name]=cached;print(json.dumps({"case":name,"cache":"reused"}),flush=True);continue
        if time.perf_counter()>deadline:
            failures.append({"case":name,"reason":"TOTAL_PREDECLARED_BUDGET_EXHAUSTED"});break
        print(json.dumps({"case":name,"action":"integrating","remaining_budget_seconds":deadline-time.perf_counter()}),flush=True)
        try:
            case_deadline=min(deadline,time.perf_counter()+config["budget"]["per_case_wall_seconds"])
            history,stats=integrate_case(discs[p],shape,initial,config,ratio,level,times,case_deadline);integrations+=1
            completed[name]=save_case(bundle,name,discs[p],history,stats,times,shape,initial,config)
            print(json.dumps({"case":name,"status":stats["status"],"integration_seconds":stats["integration_seconds"],"energy_drift":stats["diagnostics"]["max_relative_energy_drift"]}),flush=True)
            if stats["status"]!="PASS":
                failures.append({"case":name,"reason":stats["failure"]});break
        except (ArithmeticError,ValueError) as error:
            failures.append({"case":name,"reason":str(error)});write_json(bundle/"failure.json",failures);break
    comparisons={}
    def comparison(label,p_a,ratio_a,level_a,p_b,ratio_b,level_b):
        a,b=case_name(p_a,ratio_a,level_a),case_name(p_b,ratio_b,level_b)
        if a in completed and b in completed and completed[a]["status"]==completed[b]["status"]=="PASS":
            with np.load(bundle/"cases"/a/"trajectory.npz") as za,np.load(bundle/"cases"/b/"trajectory.npz") as zb:
                comparisons[label]=compare_histories(discs[p_a],za,discs[p_b],zb,config)
        else:
            comparisons[label]={"status":"PARTIAL","reason":"Required completed trajectory unavailable"}
    comparison("spatial_p16_p24",16,.05,"tight",24,.05,"tight")
    comparison("spatial_p24_p32",24,.05,"tight",32,.05,"tight")
    comparison("temporal_coarse_medium",final_p,.05,"coarse",final_p,.05,"medium")
    comparison("temporal_medium_tight",final_p,.05,"medium",final_p,.05,"tight")
    comparison("small_amplitude_spatial",24,.025,"tight",32,.025,"tight")
    statuses={"NLSP_PLANAR_DISCRETIZATION":pre["status"],"NLSP_PLANAR_LINEAR_TIME_REFERENCE":pre["status"],
              "NLSP_PLANAR_TIME_INTEGRATION":"PASS" if len(completed)==len(plan) and not failures else "PARTIAL",
              "NLSP_PLANAR_SPATIAL_CONVERGENCE":comparisons["spatial_p24_p32"]["status"],
              "NLSP_PLANAR_TEMPORAL_CONVERGENCE":comparisons["temporal_medium_tight"]["status"],
              "NLSP_PLANAR_SMALL_AMPLITUDE_LIMIT":"PARTIAL",
              "NLSP_PLANAR_ENERGY_AND_MASS":"PASS" if completed and all(r["quality_status"]=="PASS" for r in completed.values()) else "PARTIAL"}
    large=completed.get(case_name(final_p,.05,"tight"));small=completed.get(case_name(final_p,.025,"tight"))
    amplitude_check={}
    if large and small and large["status"]==small["status"]=="PASS":
        for field in FIELDS:
            key="q_"+field
            amplitude_check[field]={"large_normalized_difference":large["linear_semidiscrete_difference"]["fields"][key]["max_time_L2_difference"]/large["amplitude"],
                                    "small_normalized_difference":small["linear_semidiscrete_difference"]["fields"][key]["max_time_L2_difference"]/small["amplitude"]}
        statuses["NLSP_PLANAR_SMALL_AMPLITUDE_LIMIT"]="PASS" if all(r["small_normalized_difference"]<r["large_normalized_difference"] for r in amplitude_check.values()) else "PARTIAL"
    statuses["NLSP_PLANAR_TIME_PILOT"]="PASS" if all(value=="PASS" for value in statuses.values()) and comparisons["small_amplitude_spatial"]["status"]=="PASS" else "PARTIAL"
    effect_resolution={}
    if large and small:
        for field in FIELDS:
            key="q_"+field
            spatial=comparisons["spatial_p24_p32"].get("fields",{}).get(key,{})
            temporal=comparisons["temporal_medium_tight"].get("fields",{}).get(key,{})
            uncertainty=max(spatial.get("max_time_L2_difference",float("inf")),temporal.get("max_time_L2_difference",float("inf")))
            effect=large["continuous_linear_difference"]["q"][field]["max_time_L2_difference"]
            effect_resolution[field]={"large_amplitude_effect_L2":effect,"numerical_uncertainty_L2":uncertainty if math.isfinite(uncertainty) else None,
                                      "status":"RESOLVED" if statuses["NLSP_PLANAR_TIME_PILOT"]=="PASS" and effect>uncertainty else "EFFECT_NOT_RESOLVED"}
    summary={"statuses":statuses,"initial_eigenpair":initial,"coefficients":coefficients.values(),"config":config,"cases":completed,"comparisons":comparisons,
             "amplitude_check":amplitude_check,"failures":failures,"time_integrations_this_run":integrations,"root_solves":0,"symbolic_derivations_this_run":1,
             "effect_resolution":effect_resolution,"effect_status":"RESOLVED" if effect_resolution and all(r["status"]=="RESOLVED" for r in effect_resolution.values()) else "EFFECT_NOT_RESOLVED",
             "compute_wall_seconds":time.perf_counter()-started,"discretization_order":list(FIELDS),"not_a_periodic_orbit":True}
    write_json(bundle/"summary.json",summary)
    artifacts={str(p.relative_to(bundle)):sha(p) for p in bundle.rglob("*") if p.is_file() and p.name!="manifest.json" and p.suffix!=".tmp"}
    git={key:subprocess.run(["git"]+command,cwd=ROOT,capture_output=True,text=True).stdout.strip() for key,command in
         (("head",["rev-parse","HEAD"]),("branch",["branch","--show-current"]),("status",["status","--short"]))}
    write_json(bundle/"manifest.json",{"identity":item,"artifact_hashes":artifacts,"git":git,"command":sys.argv,"integration_seconds":sum(c["integration_seconds"] for c in completed.values())})
    return summary


def plot_bundle(bundle):
    load_runtime()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    summary=validate_bundle(bundle);bundle=Path(bundle)
    config=summary["config"];p=config["spatial"]["degrees"][-1];T=summary["initial_eigenpair"]["T1"]
    cases=[]
    for ratio in config["amplitude_over_h"]:
        name=case_name(p,ratio,"tight")
        if name in summary["cases"]:
            cases.append((ratio,summary["cases"][name],np.load(bundle/"cases"/name/"trajectory.npz")))
    if not cases:
        raise ValueError("No final-degree trajectories available to plot")
    plt.rcParams.update({"font.family":"DejaVu Serif","font.size":10,"axes.grid":True,"grid.alpha":.2})
    figures=[]
    fig,ax=plt.subplots(figsize=(7,3.4),layout="constrained")
    for ratio,case,z in cases:
        ax.plot(z["time"]/T,z["observations"][:,1,1]/case["amplitude"],label=f"A/h={ratio:g}",lw=1)
    t=cases[0][2]["time"];ax.plot(t/T,np.cos(summary["initial_eigenpair"]["omega"]*t),"k--",label="linear",lw=.8)
    ax.set(xlabel=r"$t/T_1$",ylabel=r"$w(L/2,t)/A$");ax.legend();figures.append((fig,"transverse_motion"))
    fig,axes=plt.subplots(2,1,figsize=(7,5),sharex=True,layout="constrained")
    for ratio,case,z in cases:
        axes[0].plot(z["time"]/T,z["observations"][:,0,0],label=f"A/h={ratio:g}",lw=.65)
        axes[1].plot(z["time"]/T,z["observations"][:,0,3],label=f"A/h={ratio:g}",lw=.65)
    axes[0].set(ylabel=r"$u(L/4,t)$");axes[1].set(xlabel=r"$t/T_1$",ylabel=r"$c(L/4,t)$");axes[0].legend();figures.append((fig,"generated_axial_contraction"))
    fig,ax=plt.subplots(figsize=(7,3.4),layout="constrained")
    for ratio,case,z in cases:
        ax.plot(z["time"]/T,z["energy_drift"],label=f"A/h={ratio:g}",lw=.8)
    ax.set(xlabel=r"$t/T_1$",ylabel=r"$[E_h(t)-E_h(0)]/E_h(0)$");ax.legend();figures.append((fig,"energy_drift"))
    folder=bundle/"figures";folder.mkdir(exist_ok=True)
    for fig,name in figures:
        fig.savefig(folder/(name+".pdf"));fig.savefig(folder/(name+".png"),dpi=220);plt.close(fig)
    for _,_,z in cases:z.close()
    return [str(folder/(name+".pdf")) for _,name in figures]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--check",action="store_true")
    actions.add_argument("--smoke",action="store_true")
    actions.add_argument("--compute",action="store_true")
    actions.add_argument("--report-only",type=Path)
    actions.add_argument("--plot-only",type=Path)
    parser.add_argument("--config",type=Path,default=CONFIG)
    parser.add_argument("--output-dir",type=Path,default=OUTPUT)
    args=parser.parse_args()
    if args.report_only or args.plot_only:
        bundle=args.report_only or args.plot_only
        summary=validate_bundle(bundle)
        figures=plot_bundle(bundle) if args.plot_only else []
        print(json.dumps({"bundle":str(bundle),"statuses":summary["statuses"],"figures":figures,
                          "time_integrations":0,"root_solves":0,"symbolic_derivations":0},ensure_ascii=False));return
    key,item=identity(args.config);config=item["config"];bundle=args.output_dir/key
    if args.compute and (bundle/"manifest.json").exists():
        summary=validate_bundle(bundle,item)
        print(json.dumps({"bundle":str(bundle),"cache":"reused","statuses":summary["statuses"],
                          "time_integrations":0,"root_solves":0,"symbolic_derivations":0},ensure_ascii=False));return
    for variable in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"):
        os.environ[variable]=str(config["integrator"]["blas_threads"])
    load_runtime()
    if args.check or args.smoke:
        coefficients,shape,initial,reference=setup(config)
        if args.check:
            result=preflight(config,coefficients,shape,initial,reference)
        else:
            disc=make_discretization(coefficients,16,config)
            duration=config["budget"]["smoke_duration"]
            history,result=integrate_case(disc,shape,initial,config,.05,"medium",np.linspace(0,duration,101),time.perf_counter()+120)
            result["estimated_5T1_seconds"]=result["integration_seconds"]*5*initial["T1"]/duration
        args.output_dir.mkdir(parents=True,exist_ok=True)
        write_json(args.output_dir/("preflight.json" if args.check else "smoke.json"),result)
        print(json.dumps(result,ensure_ascii=False));return
    bundle.mkdir(parents=True,exist_ok=True)
    summary=run_compute(bundle,config,item)
    write_json(args.output_dir/"current.json",{"fingerprint":key,"bundle":str(bundle.relative_to(ROOT)) if bundle.is_relative_to(ROOT) else str(bundle)})
    print(json.dumps({"bundle":str(bundle),"statuses":summary["statuses"],"failures":summary["failures"]},ensure_ascii=False))


if __name__=="__main__":
    main()
