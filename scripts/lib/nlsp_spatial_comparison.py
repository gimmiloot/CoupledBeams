"""Read-only spatial seven-field comparisons and pre-FEM signal diagnostics.

Every trajectory is loaded from saved scientific attempts. This module does
not construct an action, solve equilibrium, integrate time, or run native
FEM. Physical finite rotation matrices are separate diagnostics of the
quartic coordinate trajectory, never substituted into its equations.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_triangular

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.nlsp_spatial_1d_program import cached_case, evaluate_dense_records

VERSION = "seven-field-saved-spatial-comparison-v2-independent-mixed-p-control"
FIELDS = rod.FIELD_ORDER
DEFAULT_T1 = 10.37828159055014
HISTORICAL_OUTPUT_INDICATORS = (5e-9, 7.78842e-9)
HISTORICAL_PLANAR_PHI_FLOOR = 1.44790e-5
HISTORICAL_FINE_PLANAR_PHI_FLOOR = 1.10627e-5


def _write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf8")


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _physical_scales(length, h, cutoff, velocity=False):
    result = np.array((h, h, h, 1., 1., 1., 1.))
    return result*cutoff if velocity else result


def field_metrics(first, second, x, times, *, fixed_scale=1., numerical_floor=1e-10,
                  reference_floor=None, tolerance=None):
    """Full-horizon sampled metrics; no phase or amplitude adjustment."""
    first, second, x, times = map(np.asarray, (first, second, x, times))
    if first.shape != second.shape or first.shape != (len(times), len(x)):
        raise ValueError("Physical field/time shapes disagree")
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("Nonfinite saved field")
    difference = first-second
    l2 = np.sqrt(np.trapezoid(difference**2, x=x, axis=1))
    first_l2 = np.sqrt(np.trapezoid(first**2, x=x, axis=1))
    second_l2 = np.sqrt(np.trapezoid(second**2, x=x, axis=1))
    index = np.unravel_index(np.argmax(abs(difference)), difference.shape)
    scale_l2 = float(second_l2.max()); scale_max = float(np.max(abs(second)))
    floor = float(numerical_floor*max(scale_l2, 1e-30) if reference_floor is None else reference_floor)
    result = {"absolute_max": float(np.max(abs(difference))), "max_time_L2": float(l2.max()),
              "signed_difference_at_absolute_max": float(difference[index]),
              "maximum_time": float(times[index[0]]), "maximum_x": float(x[index[1]]),
              "reference_max_L2": scale_l2, "reference_max_abs": scale_max,
              "common_pair_max_L2": float(max(first_l2.max(), second_l2.max())),
              "common_pair_max_abs": float(max(np.max(abs(first)), np.max(abs(second)))),
              "numerical_floor": floor, "floor_limited": bool(scale_l2 <= floor or scale_max <= floor),
              "relative_L2": float(l2.max()/max(scale_l2, floor)),
              "relative_max": float(np.max(abs(difference))/max(scale_max, floor)),
              "fixed_physical_scale": float(fixed_scale),
              "relative_max_fixed": float(np.max(abs(difference))/fixed_scale),
              "relative_L2_fixed": float(l2.max()/(fixed_scale*np.sqrt(x[-1]-x[0]))),
              "time_range": [float(times[0]), float(times[-1])],
              "sampled_maxima_not_continuous_suprema": True}
    if tolerance is not None:
        result["tolerance"] = tolerance
        result["status"] = "PASS" if max(result["relative_L2"], result["relative_max"]) <= tolerance else "PARTIAL"
    return result


def physical_rotations(fields):
    """R=exp(hat(Phi,-psi,theta)), exact diagnostic of accepted coordinates."""
    fields = np.asarray(fields)
    if fields.shape[-1] != 7 or not np.isfinite(fields).all():
        raise ValueError("Canonical seven finite physical fields required")
    vectors = fields[..., (3, 4, 5)]*np.array((1., -1., 1.))
    rotations = np.empty(vectors.shape[:-1]+(3, 3))
    for index in np.ndindex(vectors.shape[:-1]):
        rotations[index] = rod.rotation_and_right_jacobian(vectors[index])[0]
    return rotations


def retained_curvatures(fields, gradients):
    """Canonical accepted chi=a_s−a×a_s/2+a×(a×a_s)/6; no new derivation."""
    fields, gradients = map(np.asarray, (fields, gradients))
    if fields.shape != gradients.shape or fields.shape[-1] != 7:
        raise ValueError("Field/gradient shapes differ")
    a = fields[..., (3, 4, 5)]*np.array((1., -1., 1.))
    ass = gradients[..., (3, 4, 5)]*np.array((1., -1., 1.))
    return ass-np.cross(a, ass)/2+np.cross(a, np.cross(a, ass))/6


def orientation_difference(first, second):
    """Principal angle of R_first.T R_second, stable at small rotations."""
    first, second = map(np.asarray, (first, second))
    if first.shape != second.shape or first.shape[-2:] != (3, 3):
        raise ValueError("Physical orientation matrices disagree")
    relative = np.swapaxes(first, -1, -2)@second
    vector = np.stack((relative[..., 2, 1]-relative[..., 1, 2],
                       relative[..., 0, 2]-relative[..., 2, 0],
                       relative[..., 1, 0]-relative[..., 0, 1]), axis=-1)/2
    sine = np.linalg.norm(vector, axis=-1)
    cosine = np.clip((np.trace(relative, axis1=-2, axis2=-1)-1)/2, -1., 1.)
    return np.arctan2(sine, cosine)


def exact_diagnostic_curvatures(fields, gradients):
    """J_r(a) a_s, separate from the retained cubic curvature of the solver."""
    fields, gradients = map(np.asarray, (fields, gradients))
    if fields.shape != gradients.shape or fields.shape[-1] != 7:
        raise ValueError("Canonical field/gradient shapes differ")
    a = fields[..., (3, 4, 5)]*np.array((1., -1., 1.))
    ass = gradients[..., (3, 4, 5)]*np.array((1., -1., 1.))
    result = np.empty_like(a)
    for index in np.ndindex(a.shape[:-1]):
        result[index] = rod.rotation_and_right_jacobian(a[index])[1]@ass[index]
    return result


def physical_curvature_from_rotations(rotations, x):
    """Two explicit derivative diagnostics of saved SO(3) section orientations.

    The first uses a nonuniform-coordinate finite difference of R and takes
    vee(skew(R.T R_s)). The second takes the principal log of each adjacent
    relative orientation divided by its actual material interval. Neither
    quantity is a native solid strain output or identical to Phi. Recovery
    and differentiation sensitivity must remain part of its qualification.
    """
    from scipy.spatial.transform import Rotation
    rotations, x = np.asarray(rotations), np.asarray(x)
    if (rotations.ndim != 4 or rotations.shape[1:] != (len(x), 3, 3)
            or len(x) < 3 or np.any(np.diff(x) <= 0) or not np.isfinite(rotations).all()):
        raise ValueError("Saved orientation/material-x dimensions are invalid")
    orthogonal = rotations@np.swapaxes(rotations, -1, -2)
    if np.max(abs(orthogonal-np.eye(3))) > 2e-8 or np.min(np.linalg.det(rotations)) <= 0:
        raise ValueError("Physical curvature requires proper saved orientation matrices")
    derivative = np.gradient(rotations, x, axis=1, edge_order=2)
    local = np.swapaxes(rotations, -1, -2)@derivative
    chi = np.stack((local[..., 2, 1]-local[..., 1, 2],
                    local[..., 0, 2]-local[..., 2, 0],
                    local[..., 1, 0]-local[..., 0, 1]), axis=-1)/2
    relative = np.swapaxes(rotations[:, :-1], -1, -2)@rotations[:, 1:]
    logarithm = Rotation.from_matrix(relative.reshape(-1, 3, 3)).as_rotvec().reshape(relative.shape[:-2]+(3,))
    geodesic = logarithm/np.diff(x)[None, :, None]
    return {"chi_gradient": chi, "interval_chi_geodesic": geodesic,
            "interval_x": (x[:-1]+x[1:])/2,
            "orthogonality_max_absolute": float(np.max(abs(orthogonal-np.eye(3)))),
            "gradient_symmetric_defect_max": float(np.max(abs((local+np.swapaxes(local, -1, -2))/2))),
            "definition": "vee(skew(R.T R_s)) or log(R_i.T R_(i+1))/dx; first component is section-orientation twist-curvature proxy",
            "Phi_is_not_chi1": True,
            "native_solid_strain_or_torsional_stress_not_reconstructed": True}


def axis_plane_diagnostic(fields, x):
    """Instantaneous best-plane residual, distinct from nonlinear coupling."""
    fields, x = np.asarray(fields), np.asarray(x)
    if fields.ndim != 3 or fields.shape[1:] != (len(x), 7):
        raise ValueError("Expected (nt,nx,7) physical axis fields")
    weights = np.empty(len(x))
    weights[1:-1] = (x[2:]-x[:-2])/2
    weights[0], weights[-1] = (x[1]-x[0])/2, (x[-1]-x[-2])/2
    mass = weights.sum()
    singular, normals, maximum, rms = [], [], [], []
    for row in fields:
        axis = np.column_stack((x+row[:, 0], row[:, 1], row[:, 2]))
        centered = axis-(weights@axis)/mass
        _, values, right = np.linalg.svd(centered*np.sqrt(weights[:, None]/mass), full_matrices=False)
        normal = right[-1]
        if normals and np.dot(normal, normals[0]) < 0:
            normal = -normal
        distance = centered@normal
        singular.append(values); normals.append(normal)
        maximum.append(np.max(abs(distance))); rms.append(np.sqrt(weights@(distance**2)/mass))
    normals = np.asarray(normals)
    return {"singular_values": np.asarray(singular), "normal": normals,
            "max_distance": np.asarray(maximum), "RMS_distance": np.asarray(rms),
            "plane_change_from_initial_radians": np.arccos(np.clip(abs(normals@normals[0]), 0., 1.)),
            "definition": "original material x, current axis, trapezoid-weighted centered instantaneous plane SVD",
            "nonplanarity_is_not_by_itself_nonlinear_coupling": True}


def planning_signal(signal, p_difference, *, indicators=HISTORICAL_OUTPUT_INDICATORS, multiplier=10.,
                    change_label="observed_p48_p64_change"):
    values = [float(signal), float(p_difference)]+list(map(float, indicators))
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Invalid nonnegative planning signal/indicator")
    indicator = max(values[1:]); threshold = multiplier*indicator
    return {"signal": values[0], change_label: values[1],
            "historical_output_recovery_indicators": values[2:], "planning_multiplier": multiplier,
            "planning_threshold": threshold, "signal_to_largest_indicator": signal/indicator if indicator else None,
            "resolved_for_pre_FEM_planning": signal > threshold,
            "not_physical_accuracy_gate_or_strict_error_bound": True}


def _load_case(path):
    path = Path(path); metadata = cached_case(path)
    if metadata is None:
        raise ValueError("Missing saved scientific case: "+str(path))
    with np.load(path/"trajectory.npz", allow_pickle=False) as z:
        nonlinear = {k: z[k].copy() for k in ("times", "q", "velocity")}
    with np.load(path/"linear_trajectory.npz", allow_pickle=False) as z:
        linear = {k: z[k].copy() for k in ("times", "q", "velocity")}
    if not np.array_equal(nonlinear["times"], linear["times"]):
        raise ValueError("Own linear/nonlinear actual timestamps differ")
    return metadata, nonlinear, linear


def _summary_curve(values, x, times, midpoint=None):
    l2 = np.sqrt(np.trapezoid(values**2, x=x, axis=1))
    maximum = np.max(abs(values), axis=1)
    index = np.unravel_index(np.argmax(abs(values)), values.shape)
    mid = int(np.argmin(abs(x-(x[-1]+x[0])/2))) if midpoint is None else midpoint
    return {"absolute_max": float(maximum.max()), "max_time_L2": float(l2.max()),
            "signed_maximum": float(values[index]), "maximum_time": float(times[index[0]]),
            "maximum_x": float(x[index[1]]), "midspan_absolute_max": float(np.max(abs(values[:, mid]))),
            "initial_midspan": float(values[0, mid]), "final_midspan": float(values[-1, mid])}, l2, maximum, values[:, mid]


def analyze_stage_b(cases, discs, output, *, T1=DEFAULT_T1, h=.1, length=1., chunk_size=64):
    """Saved four/six-case diagnostic; reconstruct in blocks, never integrate.

    Optional isolated_w_p48/isolated_v_p48 must occur as a pair. They enable
    an actual p comparison of the mixed response, rather than treating joint
    sensitivity as the uncertainty of three separately evolved solutions.
    """
    required = ("joint_p48", "joint_p64", "isolated_w_p64", "isolated_v_p64")
    if any(name not in cases for name in required):
        raise ValueError("Joint p48/p64 and exactly two isolated p64 controls required")
    optional = ("isolated_w_p48", "isolated_v_p48")
    if any(name in cases for name in optional) and not all(name in cases for name in optional):
        raise ValueError("Both isolated p48 controls are required for an independent mixed p comparison")
    mixed_p_available = all(name in cases for name in optional)
    required = required+optional if mixed_p_available else required
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    loaded = {name: _load_case(cases[name]) for name in required}
    times = loaded["joint_p64"][1]["times"]
    for name, (_, nonlinear, linear) in loaded.items():
        if not np.array_equal(times, nonlinear["times"]):
            raise ValueError("Stage B comparisons require the same exact saved output timestamps: "+name)
    if times[0] != 0 or times[-1] > .25*T1*(1+1e-12):
        raise ValueError("Saved Stage B horizon exceeds the frozen quarter period")
    x = np.linspace(0., length, 801)
    gauss_x, gauss_weights = np.polynomial.legendre.leggauss(100)
    gauss_x, gauss_weights = (gauss_x+1)*length/2, gauss_weights*length/2
    snapshots = np.unique([int(np.argmin(abs(times-t))) for t in np.linspace(0., times[-1], 5)])
    positions = (.25*length, .5*length, .5*length, .5*length, .25*length, .25*length, .25*length)
    observation_indices = np.array([int(np.argmin(abs(x-p))) for p in positions])
    indices = np.arange(7)
    fields, velocities, gradients, gauss_fields, gauss_velocities = {}, {}, {}, {}, {}
    # Small time-series and representative profiles are retained, rather than
    # all four dense physical space/time reconstructions.
    data = {"times": times, "x": x, "snapshot_indices": snapshots,
            "snapshot_times": times[snapshots], "observations_x": np.array(positions)}
    parts = ("q", "velocity", "linear_q", "linear_velocity", "correction", "evolution",
             "velocity_correction", "velocity_evolution")+( ("mixed", "mixed_evolution") if mixed_p_available else () )
    spatial = {part: [] for part in parts}
    accum = {}
    for part in spatial:
        count = 3 if part in ("mixed", "mixed_evolution") else 7
        accum[part] = {key: np.zeros(count) for key in
            ("abs_max", "abs_l2", "reference_l2", "reference_max", "common_l2", "common_max", "max_time", "max_x", "signed_max")}
    metric_curves = {part: {"L2": [], "max": []} for part in spatial}
    response_curves = {part: [] for part in ("delta", "evolution", "mix", "mix_evolution")}
    response_l2 = {part: [] for part in response_curves}
    response_max = {part: [] for part in response_curves}
    response_extrema = {part: {"max": np.zeros(7 if part in ("delta", "evolution") else 3),
        "signed": np.zeros(7 if part in ("delta", "evolution") else 3),
        "time": np.zeros(7 if part in ("delta", "evolution") else 3),
        "x": np.zeros(7 if part in ("delta", "evolution") else 3)} for part in response_curves}
    initial_delta, initial_mix, initial_gauss_delta, initial_gauss_mix = None, None, None, None
    initial_velocity_delta = initial_gauss_velocity_delta = initial_mix48 = initial_gauss_mix48 = None
    plane_linear, plane_nonlinear = [], []
    curvature_spatial_max = np.zeros(3)
    torsion_max = torsion_sensitivity = 0.
    for start in range(0, len(times), chunk_size):
        stop = min(start+chunk_size, len(times))
        block = slice(start, stop)
        for name in required:
            degree = 48 if name.endswith("p48") else 64
            _, nonlinear, linear = loaded[name]; disc = discs[degree]
            fields[name] = disc.reconstruct_series(nonlinear["q"][block], x)
            gauss_fields[name] = disc.reconstruct_series(nonlinear["q"][block], gauss_x)
            if name.startswith("joint"):
                fields[name+"_linear"] = disc.reconstruct_series(linear["q"][block], x)
                gauss_fields[name+"_linear"] = disc.reconstruct_series(linear["q"][block], gauss_x)
                velocities[name] = disc.reconstruct_series(nonlinear["velocity"][block], x)
                gauss_velocities[name] = disc.reconstruct_series(nonlinear["velocity"][block], gauss_x)
                velocities[name+"_linear"] = disc.reconstruct_series(linear["velocity"][block], x)
                gauss_velocities[name+"_linear"] = disc.reconstruct_series(linear["velocity"][block], gauss_x)
                gradients[name] = disc.reconstruct_series(nonlinear["q"][block], x, derivative=1)
                data.setdefault(name+"_linear_observations", []).append(fields[name+"_linear"][:, observation_indices, indices])
            observations = fields[name][:, observation_indices, indices]
            data.setdefault(name+"_observations", []).append(observations)
        first, second = fields["joint_p48"], fields["joint_p64"]
        delta48 = first-fields["joint_p48_linear"]
        delta64 = second-fields["joint_p64_linear"]
        delta_gauss48 = gauss_fields["joint_p48"]-gauss_fields["joint_p48_linear"]
        delta_gauss64 = gauss_fields["joint_p64"]-gauss_fields["joint_p64_linear"]
        mix = second[..., :3]-fields["isolated_w_p64"][..., :3]-fields["isolated_v_p64"][..., :3]
        gauss_mix = gauss_fields["joint_p64"][..., :3]-gauss_fields["isolated_w_p64"][..., :3]-gauss_fields["isolated_v_p64"][..., :3]
        velocity_delta48 = velocities["joint_p48"]-velocities["joint_p48_linear"]
        velocity_delta64 = velocities["joint_p64"]-velocities["joint_p64_linear"]
        gauss_velocity_delta48 = gauss_velocities["joint_p48"]-gauss_velocities["joint_p48_linear"]
        gauss_velocity_delta64 = gauss_velocities["joint_p64"]-gauss_velocities["joint_p64_linear"]
        if mixed_p_available:
            mix48 = first[..., :3]-fields["isolated_w_p48"][..., :3]-fields["isolated_v_p48"][..., :3]
            gauss_mix48 = gauss_fields["joint_p48"][..., :3]-gauss_fields["isolated_w_p48"][..., :3]-gauss_fields["isolated_v_p48"][..., :3]
        if initial_delta is None:
            initial_delta = {48: delta48[0].copy(), 64: delta64[0].copy()}
            initial_mix = mix[0].copy()
            initial_gauss_delta = {48: delta_gauss48[0].copy(), 64: delta_gauss64[0].copy()}
            initial_gauss_mix = gauss_mix[0].copy()
            initial_velocity_delta = {48: velocity_delta48[0].copy(), 64: velocity_delta64[0].copy()}
            initial_gauss_velocity_delta = {48: gauss_velocity_delta48[0].copy(), 64: gauss_velocity_delta64[0].copy()}
            if mixed_p_available:
                initial_mix48, initial_gauss_mix48 = mix48[0].copy(), gauss_mix48[0].copy()
        evolution48, evolution64 = delta48-initial_delta[48], delta64-initial_delta[64]
        pairs = {"q": (first, second), "velocity": (velocities["joint_p48"], velocities["joint_p64"]),
                 "linear_q": (fields["joint_p48_linear"], fields["joint_p64_linear"]),
                 "linear_velocity": (velocities["joint_p48_linear"], velocities["joint_p64_linear"]),
                 "correction": (delta48, delta64), "evolution": (evolution48, evolution64),
                 "velocity_correction": (velocity_delta48, velocity_delta64),
                 "velocity_evolution": (velocity_delta48-initial_velocity_delta[48], velocity_delta64-initial_velocity_delta[64])}
        gauss_pairs = {"q": (gauss_fields["joint_p48"], gauss_fields["joint_p64"]),
                       "velocity": (gauss_velocities["joint_p48"], gauss_velocities["joint_p64"]),
                       "linear_q": (gauss_fields["joint_p48_linear"], gauss_fields["joint_p64_linear"]),
                       "linear_velocity": (gauss_velocities["joint_p48_linear"], gauss_velocities["joint_p64_linear"]),
                       "correction": (delta_gauss48, delta_gauss64),
                       "evolution": (delta_gauss48-initial_gauss_delta[48], delta_gauss64-initial_gauss_delta[64]),
                       "velocity_correction": (gauss_velocity_delta48, gauss_velocity_delta64),
                       "velocity_evolution": (gauss_velocity_delta48-initial_gauss_velocity_delta[48], gauss_velocity_delta64-initial_gauss_velocity_delta[64])}
        if mixed_p_available:
            pairs.update(mixed=(mix48, mix), mixed_evolution=(mix48-initial_mix48, mix-initial_mix))
            gauss_pairs.update(mixed=(gauss_mix48, gauss_mix),
                               mixed_evolution=(gauss_mix48-initial_gauss_mix48, gauss_mix-initial_gauss_mix))
        for part, (aa, bb) in pairs.items():
            difference = aa-bb
            ga, gb = gauss_pairs[part]
            d_l2 = np.sqrt(np.einsum("tij,i,tij->tj", ga-gb, gauss_weights, ga-gb))
            d_max = np.max(abs(difference), axis=1)
            local = accum[part]
            local["abs_l2"] = np.maximum(local["abs_l2"], d_l2.max(axis=0))
            local["reference_l2"] = np.maximum(local["reference_l2"], np.sqrt(np.einsum("tij,i,tij->tj", gb, gauss_weights, gb)).max(axis=0))
            local["reference_max"] = np.maximum(local["reference_max"], np.max(abs(bb), axis=(0, 1)))
            local["common_l2"] = np.maximum(local["common_l2"], np.maximum(
                np.sqrt(np.einsum("tij,i,tij->tj", ga, gauss_weights, ga)).max(axis=0),
                np.sqrt(np.einsum("tij,i,tij->tj", gb, gauss_weights, gb)).max(axis=0)))
            local["common_max"] = np.maximum(local["common_max"], np.maximum(np.max(abs(aa), axis=(0, 1)), np.max(abs(bb), axis=(0, 1))))
            for j in range(difference.shape[-1]):
                ij = np.unravel_index(np.argmax(abs(difference[..., j])), difference[..., j].shape)
                maximum = abs(difference[ij[0], ij[1], j])
                if maximum >= local["abs_max"][j]:
                    local["abs_max"][j], local["max_time"][j], local["max_x"][j] = maximum, times[start+ij[0]], x[ij[1]]
                    local["signed_max"][j] = difference[ij[0], ij[1], j]
            metric_curves[part]["L2"].append(d_l2); metric_curves[part]["max"].append(d_max)
        response = {"delta": delta64, "evolution": evolution64, "mix": mix, "mix_evolution": mix-initial_mix}
        data.setdefault("mix_evolution_observations", []).append((mix-initial_mix)[:, observation_indices[:3], np.arange(3)])
        gauss_response = {"delta": delta_gauss64, "evolution": delta_gauss64-initial_gauss_delta[64],
                          "mix": gauss_mix, "mix_evolution": gauss_mix-initial_gauss_mix}
        for part, values in response.items():
            gv = gauss_response[part]
            response_l2[part].append(np.sqrt(np.einsum("tij,i,tij->tj", gv, gauss_weights, gv)))
            response_max[part].append(np.max(abs(values), axis=1))
            response_curves[part].append(values[:, len(x)//2])
            extrema = response_extrema[part]
            for j in range(values.shape[-1]):
                ij = np.unravel_index(np.argmax(abs(values[..., j])), values[..., j].shape)
                if abs(values[ij[0], ij[1], j]) >= extrema["max"][j]:
                    extrema["max"][j] = abs(values[ij[0], ij[1], j])
                    extrema["signed"][j] = values[ij[0], ij[1], j]
                    extrema["time"][j], extrema["x"][j] = times[start+ij[0]], x[ij[1]]
        chi48, chi64 = retained_curvatures(first, gradients["joint_p48"]), retained_curvatures(second, gradients["joint_p64"])
        curvature_spatial_max = np.maximum(curvature_spatial_max, np.max(abs(chi48-chi64), axis=(0, 1)))
        torsion_max = max(torsion_max, float(np.max(abs(chi64[..., 0]))))
        torsion_sensitivity = max(torsion_sensitivity, float(np.max(abs(chi64[..., 0]-chi48[..., 0]))))
        plane_linear.append(axis_plane_diagnostic(fields["joint_p64_linear"], x))
        plane_nonlinear.append(axis_plane_diagnostic(second, x))
        for global_index in snapshots[(snapshots >= start) & (snapshots < stop)]:
            k = global_index-start
            for name in ("joint_p48", "joint_p64", "joint_p64_linear"):
                data.setdefault(name+"_snapshots", []).append(fields[name][k])
            data.setdefault("physical_rotation_snapshots", []).append(physical_rotations(second[k]))
            data.setdefault("retained_curvature_snapshots", []).append(chi64[k])
            data.setdefault("mixed_translation_snapshots", []).append(mix[k])
    cutoff = math.sqrt(discs[64].coefficients.C/discs[64].coefficients.jp)
    metrics = {}
    for part, local in accum.items():
        rows = {}
        floor = 1e-10*max(float(local["reference_l2"].max()), 1e-30)
        fixed = _physical_scales(length, h, cutoff, "velocity" in part)
        for j, name in enumerate(FIELDS[:len(local["abs_max"])]):
            tolerance = 1e-3 if name in ("u", "c") else 1e-4
            relative_l2 = float(local["abs_l2"][j]/max(local["reference_l2"][j], floor))
            relative_max = float(local["abs_max"][j]/max(local["reference_max"][j], floor))
            rows[name] = {"absolute_max": float(local["abs_max"][j]), "max_time_L2": float(local["abs_l2"][j]),
                "reference_max_L2": float(local["reference_l2"][j]), "reference_max_abs": float(local["reference_max"][j]),
                "common_pair_max_L2": float(local["common_l2"][j]), "common_pair_max_abs": float(local["common_max"][j]),
                "relative_L2": relative_l2, "relative_max": relative_max, "numerical_floor": floor,
                "floor_limited": bool(local["reference_l2"][j] <= floor or local["reference_max"][j] <= floor),
                "fixed_physical_scale": float(fixed[j]), "relative_max_fixed": float(local["abs_max"][j]/fixed[j]),
                "relative_L2_fixed": float(local["abs_l2"][j]/(fixed[j]*math.sqrt(length))),
                "maximum_time": float(local["max_time"][j]), "maximum_x": float(local["max_x"][j]),
                "signed_difference_at_absolute_max": float(local["signed_max"][j]), "tolerance": tolerance,
                "status": "PASS" if max(relative_l2, relative_max) <= tolerance else "PARTIAL"}
        metrics[part] = {"status": "PASS" if all(v["status"] == "PASS" for v in rows.values()) else "PARTIAL", "fields": rows}
        for quantity in ("L2", "max"):
            data["spatial_"+part+"_"+quantity] = np.concatenate(metric_curves[part][quantity])
    for name in required:
        data[name+"_observations"] = np.concatenate(data[name+"_observations"])
        if name.startswith("joint"):
            data[name+"_linear_observations"] = np.concatenate(data[name+"_linear_observations"])
    data["mix_evolution_observations"] = np.concatenate(data["mix_evolution_observations"])
    for name in list(data):
        if isinstance(data[name], list):
            data[name] = np.asarray(data[name])
    responses = {}
    for part in response_curves:
        mid, l2, maximum = (np.concatenate(mapping[part]) for mapping in (response_curves, response_l2, response_max))
        data[part+"_midspan"], data[part+"_L2"], data[part+"_max"] = mid, l2, maximum
        responses[part] = {FIELDS[j]: {"absolute_max": float(maximum[:, j].max()), "max_time_L2": float(l2[:, j].max()),
            "midspan_absolute_max": float(np.max(abs(mid[:, j]))), "initial_midspan": float(mid[0, j]),
            "final_midspan": float(mid[-1, j]), "signed_global_maximum": float(response_extrema[part]["signed"][j]),
            "maximum_time": float(response_extrema[part]["time"][j]),
            "maximum_x": float(response_extrema[part]["x"][j])} for j in range(mid.shape[1])}
    plane_summary = {}
    for name, blocks in (("linear", plane_linear), ("nonlinear", plane_nonlinear)):
        for key in ("singular_values", "normal", "max_distance", "RMS_distance"):
            data["axis_"+name+"_"+key] = np.concatenate([block[key] for block in blocks])
        normals = data["axis_"+name+"_normal"]
        data["axis_"+name+"_plane_change"] = np.arccos(np.clip(abs(normals@normals[0]), 0., 1.))
        plane_summary[name] = {"max_plane_distance": float(data["axis_"+name+"_max_distance"].max()),
            "max_RMS_plane_distance": float(data["axis_"+name+"_RMS_distance"].max()),
            "max_plane_change_radians": float(data["axis_"+name+"_plane_change"].max()),
            "instantaneous_nonplanarity_is_not_a_nonlinear_coupling_proof": True}
    bending_planning = {name: planning_signal(responses["evolution"][name]["absolute_max"],
        metrics["evolution"]["fields"][name]["absolute_max"]) for name in ("w", "v")}
    mixed_planning = {name: planning_signal(responses["mix_evolution"][name]["absolute_max"],
        metrics["mixed_evolution" if mixed_p_available else "evolution"]["fields"][name]["absolute_max"])
        for name in ("u", "w", "v")}
    phi = metrics["q"]["fields"]["Phi"]
    torsion = {"Phi_max": phi["reference_max_abs"], "Phi_p48_p64_max_difference": phi["absolute_max"],
        "retained_chi1_max": torsion_max, "retained_chi1_p48_p64_difference": torsion_sensitivity,
        "historical_planar_medium_Phi_proxy_floor": HISTORICAL_PLANAR_PHI_FLOOR,
        "historical_planar_fine_Phi_proxy_floor": HISTORICAL_FINE_PLANAR_PHI_FLOOR,
        "Phi_pre_FEM_planning": planning_signal(phi["reference_max_abs"], phi["absolute_max"],
            indicators=(HISTORICAL_PLANAR_PHI_FLOOR,), multiplier=10.),
        "chi1_does_not_equal_Phi_or_global_surface_rotation": True,
        "physical_torsion_verification": "NOT_RESOLVED_BEFORE_NATIVE_3D_RECOVERY"}
    primary_spatial = all(metrics["q"]["fields"][field]["status"] == "PASS" for field in ("w", "v"))
    evolution_resolved = all(row["resolved_for_pre_FEM_planning"] for row in bending_planning.values())
    mixed_resolved = any(mixed_planning[field]["resolved_for_pre_FEM_planning"] for field in ("w", "v"))
    own_cases_completed = all(row[0]["status"] == "PASS" for row in loaded.values())
    reached_target = math.isclose(times[-1], .25*T1, rel_tol=1e-12, abs_tol=0.)
    all14_passed = metrics["q"]["status"] == metrics["velocity"]["status"] == "PASS"
    stage_c_allowed = (primary_spatial and evolution_resolved and mixed_resolved and own_cases_completed
                       and reached_target and mixed_p_available)
    summary = {"version": VERSION,
        "status": "STAGE_B_SPATIAL_SIGNAL_RESOLVED_WITH_QUALIFICATIONS" if stage_c_allowed and all14_passed else "NUMERICAL_PARTIAL",
        "scientific_calls": 0, "stage_c_allowed": stage_c_allowed,
        "resolved_spatial_bending": stage_c_allowed,
        "stage_c_decision": {"all_requested_cases_completed": own_cases_completed, "reached_frozen_horizon": reached_target,
            "primary_w_v_displacement_spatial_gates": primary_spatial,
            "full_velocity_spatial_status": metrics["velocity"]["status"],
            "independent_mixed_p_controls_available": mixed_p_available,
            "both_w_v_evolving_corrections_exceed_planning_indicators": evolution_resolved,
            "at_least_one_transverse_evolving_mixed_response_exceeds_planning_indicators": mixed_resolved,
            "rule": "all requested own cases complete; primary w/v displacement gates pass; both evolving bend corrections and at least one transverse evolving mixed response exceed ten times actual matching p-change/output indicators; full velocities retain their separate status",
            "isolated_control_uncertainty_not_independently_certified": not mixed_p_available,
            "this_decision_is_not_a_physical_accuracy_or_complete_seven_field_PASS": True},
        "actual_time_end": float(times[-1]),
        "target_time_end": .25*T1, "T1": T1,
        "sampling": "saved exact common timestamps; 100 Gauss nodes for L2 (old policy), 801 material-x sampled maxima",
        "all14_spatial": {"status": "PASS" if all14_passed else "PARTIAL",
                          "q": metrics["q"], "velocity": metrics["velocity"]},
        "spatial_correction": metrics["correction"], "spatial_evolution": metrics["evolution"],
        "spatial_linear_displacement": metrics["linear_q"], "spatial_linear_velocity": metrics["linear_velocity"],
        "spatial_velocity_correction": metrics["velocity_correction"], "spatial_velocity_evolution": metrics["velocity_evolution"],
        "spatial_mixed": metrics.get("mixed", {"status": "NOT_RUN"}),
        "spatial_mixed_evolution": metrics.get("mixed_evolution", {"status": "NOT_RUN"}),
        "responses": responses, "bending_signal_planning": bending_planning, "mixed_signal_planning": mixed_planning,
        "isolated_p64_controls_have_no_independent_p_control": not mixed_p_available,
        "mixed_response_p_sensitivity_uses_joint_control_only": not mixed_p_available,
        "independent_mixed_p_control_available": mixed_p_available,
        "torsion": torsion, "retained_curvature_p_changes": curvature_spatial_max.tolist(),
        "axis_plane": plane_summary, "full_horizon_fixed_denominators": True,
        "historical_p64_reference_gate_denominators_preserved": True,
        "common_pair_characteristic_scales_are_additional_not_gate_replacements": True,
        "physical_rotation_definition": "exact exp(hat(Phi,-psi,theta)) diagnostic; quartic solver unchanged",
        "source_cases": {name: {"path": str(Path(cases[name])), "case_sha256": _sha(Path(cases[name])/"case.json"),
            "status": loaded[name][0]["status"]} for name in required},
        "temporal_1D": "PARTIAL_SINGLE_TIGHT_LEVEL", "sampled_maxima_not_continuous_suprema": True}
    np.savez_compressed(output/"stage_b_comparison.npz", **data)
    _write(output/"stage_b_comparison.json", summary)
    with (output/"all14_spatial.csv").open("w", newline="", encoding="utf8") as stream:
        columns = ("component", "absolute_max", "max_time_L2", "relative_max", "relative_L2", "tolerance", "status")
        writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader()
        for part in ("q", "velocity"):
            for name, row in metrics[part]["fields"].items():
                writer.writerow({"component": name+("_t" if part == "velocity" else ""), **{key: row[key] for key in columns[1:]}})
    return summary


def evaluate_saved_one_d(case, disc, times):
    """Exact accepted nonlinear dense polynomials and saved linear modes."""
    case = Path(case); cached_case(case)
    times = np.asarray(times, dtype=float)
    states = evaluate_dense_records(case/"accepted_dense.npz", times)
    with np.load(case/"linear_modes.npz", allow_pickle=False) as z:
        omega, vectors, mass = z["omega"], z["vectors"], z["M0"]
    with np.load(case/"static_states.npz", allow_pickle=False) as z:
        initial = z["q_linear"]
    if vectors.shape != (disc.ndof, disc.ndof) or initial.shape != (disc.ndof,):
        raise ValueError("Saved full linear basis identity/dimension mismatch")
    amplitudes = vectors.T@(mass@initial)
    angles = times[:, None]*omega
    q_linear = (np.cos(angles)*amplitudes)@vectors.T
    v_linear = (-np.sin(angles)*(omega*amplitudes))@vectors.T
    if len(times) and times[0] == 0:
        q_linear[0], v_linear[0] = initial, 0.
    return {"times": times, "q_nonlinear": states[:, :disc.ndof], "v_nonlinear": states[:, disc.ndof:],
            "q_linear": q_linear, "v_linear": v_linear, "scientific_calls": 0}


def read_fem_trajectory(case):
    """Read saved native sections and separately confirmed static preload."""
    case = Path(case)
    with np.load(case/"section_history.npz", allow_pickle=False) as z:
        history = {name: z[name].copy() for name in
                   ("time", "x", "fields", "raw_rotation_x", "raw_rotation_matrices")}
    with np.load(case/"initial_sections.npz", allow_pickle=False) as z:
        initial = {name: z[name].copy() for name in
                   ("x", "fields", "raw_rotation_x", "raw_rotation_matrices")}
    if (history["time"].ndim != 1 or np.any(history["time"] <= 0)
            or np.any(np.diff(history["time"]) <= 0)):
        raise ValueError("Saved native dynamic frames must have positive increasing physical times")
    if not np.array_equal(history["x"], initial["x"]):
        raise ValueError("Material section coordinates changed between preload and dynamics")
    if not np.array_equal(history["raw_rotation_x"], initial["raw_rotation_x"]):
        raise ValueError("Raw material rotation centers changed between preload and dynamics")
    history["time"] = np.r_[0., history["time"]]
    history["fields"] = np.concatenate((initial["fields"][None, ...], history["fields"]))
    history["raw_rotation_matrices"] = np.concatenate((initial["raw_rotation_matrices"][None, ...], history["raw_rotation_matrices"]))
    history["origin"] = "confirmed_static_preload_at_0_then_actual_native_dynamic_samples"
    history["zero_frame_is_native_dynamic"] = False
    return history


def _interpolate(values, source_times, target_times, method):
    values, source_times, target_times = map(np.asarray, (values, source_times, target_times))
    if (source_times.ndim != 1 or np.any(np.diff(source_times) <= 0) or target_times.ndim != 1
            or not np.isfinite(target_times).all() or np.any(np.diff(target_times) <= 0)
            or target_times[0] < source_times[0] or target_times[-1] > source_times[-1]):
        raise ValueError("Interpolation must remain inside actual source coverage")
    if method == "pchip":
        return PchipInterpolator(source_times, values, axis=0, extrapolate=False)(target_times)
    if method != "linear":
        raise ValueError("Only predeclared linear and shape-preserving cubic interpolation")
    flat = values.reshape(len(source_times), -1)
    return np.column_stack([np.interp(target_times, source_times, row) for row in flat.T]).reshape(
        (len(target_times),)+values.shape[1:])


def _project_so3(matrices):
    left, _, right = np.linalg.svd(np.asarray(matrices))
    sign = np.linalg.det(left@right)
    left[..., :, -1] *= sign[..., None]
    return left@right


def _raw_rotations_on_sections(history):
    """Linear space transfer of full R with known face orientation R=I."""
    x, raw_x = np.asarray(history["x"]), np.asarray(history["raw_rotation_x"])
    rotations = np.asarray(history["raw_rotation_matrices"])
    if rotations.shape != (len(history["time"]), len(raw_x), 3, 3):
        raise ValueError("Saved full section orientation dimensions disagree")
    source_x = np.r_[x[0], raw_x, x[-1]]
    identities = np.broadcast_to(np.eye(3), (len(history["time"]), 1, 3, 3))
    values = np.concatenate((identities, rotations, identities), axis=1)
    transferred = _interpolate(np.swapaxes(values, 0, 1), source_x, x, "linear")
    return _project_so3(np.swapaxes(transferred, 0, 1))


def compare_fem_pair(linear, nonlinear, case_1d, disc, output, *, T1=DEFAULT_T1,
                     common_times=None, h=.1, stage_b=None):
    """Read-only 1D/3D comparison on the predeclared 201-point time grid.

    linear/nonlinear may be saved case paths or dictionaries produced by
    read_fem_trajectory. Full finite orientations are compared as matrices;
    coordinate differences never replace a physical rotation composition.
    """
    if not isinstance(linear, dict):
        linear = read_fem_trajectory(linear)
    if not isinstance(nonlinear, dict):
        nonlinear = read_fem_trajectory(nonlinear)
    if not np.array_equal(linear["x"], nonlinear["x"]):
        raise ValueError("Linear/nonlinear recovered material grids differ")
    x = np.asarray(linear["x"])
    grid = np.linspace(0., .25*T1, 201) if common_times is None else np.asarray(common_times)
    if len(grid) != 201 or grid[0] != 0. or np.any(np.diff(grid) <= 0):
        raise ValueError("The predeclared common 201-point physical time grid is required")
    if grid[-1] > min(linear["time"][-1], nonlinear["time"][-1]):
        raise ValueError("FEM actual overlap does not reach the requested comparison grid")
    values = evaluate_saved_one_d(case_1d, disc, grid)
    one_l = disc.reconstruct_series(values["q_linear"], x)
    one_n = disc.reconstruct_series(values["q_nonlinear"], x)
    one_delta = one_n-one_l; one_evolution = one_delta-one_delta[0]
    one_R_l, one_R_n = physical_rotations(one_l), physical_rotations(one_n)
    one_gradient_l = disc.reconstruct_series(values["q_linear"], x, derivative=1)
    one_gradient_n = disc.reconstruct_series(values["q_nonlinear"], x, derivative=1)
    one_chi_l, one_chi_n = exact_diagnostic_curvatures(one_l, one_gradient_l), exact_diagnostic_curvatures(one_n, one_gradient_n)
    raw_R_l, raw_R_n = _raw_rotations_on_sections(linear), _raw_rotations_on_sections(nonlinear)
    reconstructed = {}
    for method in ("linear", "pchip"):
        lin = _interpolate(linear["fields"], linear["time"], grid, method)
        non = _interpolate(nonlinear["fields"], nonlinear["time"], grid, method)
        delta = non-lin
        reconstructed[method] = {"linear": lin, "nonlinear": non, "correction": delta,
            "evolution": delta-delta[0],
            "R_linear": _project_so3(_interpolate(raw_R_l, linear["time"], grid, method)),
            "R_nonlinear": _project_so3(_interpolate(raw_R_n, nonlinear["time"], grid, method))}
    primary, alternative = reconstructed["linear"], reconstructed["pchip"]
    one_values = {"linear": one_l, "nonlinear": one_n, "correction": one_delta, "evolution": one_evolution}
    metrics = {}
    interpolation = {}
    for kind in one_values:
        metrics[kind], interpolation[kind] = {}, {}
        floor = 1e-10*max(float(np.sqrt(np.trapezoid(primary[kind]**2, x=x, axis=1)).max()), 1e-30)
        for i, field in enumerate(FIELDS):
            fixed = h if i < 3 else 1.
            row = field_metrics(one_values[kind][..., i], primary[kind][..., i], x, grid,
                                fixed_scale=fixed, reference_floor=floor)
            row["relative_max_common"] = row["absolute_max"]/max(row["common_pair_max_abs"], floor)
            row["relative_L2_common"] = row["max_time_L2"]/max(row["common_pair_max_L2"], floor)
            if field == "c":
                row["qualification"] = "3D effective contraction proxy; not identical M-H generalized coordinate"
            if field in ("Phi", "psi", "theta"):
                row["qualification"] = "common rotation-vector chart component; primary finite orientation comparison uses full R"
            metrics[kind][field] = row
            interpolation[kind][field] = field_metrics(primary[kind][..., i], alternative[kind][..., i], x, grid,
                                                       fixed_scale=fixed)
    one_correction_R = np.swapaxes(one_R_l, -1, -2)@one_R_n
    fem_correction_R = np.swapaxes(primary["R_linear"], -1, -2)@primary["R_nonlinear"]
    one_evolution_R = np.swapaxes(one_correction_R[0], -1, -2)[None, ...]@one_correction_R
    fem_evolution_R = np.swapaxes(fem_correction_R[0], -1, -2)[None, ...]@fem_correction_R
    orientation = {"linear": orientation_difference(one_R_l, primary["R_linear"]),
                   "nonlinear": orientation_difference(one_R_n, primary["R_nonlinear"]),
                   "NL_minus_L_change_difference": orientation_difference(
                       one_correction_R, fem_correction_R),
                   "evolving_NL_minus_L_change_difference": orientation_difference(one_evolution_R, fem_evolution_R),
                   "interpolation_linear": orientation_difference(primary["R_linear"], alternative["R_linear"]),
                   "interpolation_nonlinear": orientation_difference(primary["R_nonlinear"], alternative["R_nonlinear"])}
    orientation_summary = {kind: {"maximum_principal_angle_radians": float(value.max()),
        "max_time_L2_angle": float(np.sqrt(np.trapezoid(value**2, x=x, axis=1)).max())} for kind, value in orientation.items()}
    planning = {}
    for field in ("w", "v"):
        i = FIELDS.index(field)
        signal = float(np.max(abs(primary["evolution"][..., i])))
        p_difference = 0. if stage_b is None else stage_b["spatial_evolution"]["fields"][field]["absolute_max"]
        interpolation_difference = interpolation["evolution"][field]["absolute_max"]
        planning[field] = planning_signal(signal, p_difference,
            indicators=HISTORICAL_OUTPUT_INDICATORS+(interpolation_difference,))
        planning[field]["single_mesh_single_timestep_not_independent_certification"] = True
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    payload = {"time": grid, "x": x, "one_d_linear": one_l, "one_d_nonlinear": one_n,
               "one_d_correction": one_delta, "one_d_evolution": one_evolution,
               "one_d_R_linear": one_R_l, "one_d_R_nonlinear": one_R_n}
    payload.update({"one_d_gradients_linear": one_gradient_l, "one_d_gradients_nonlinear": one_gradient_n,
        "one_d_chi_exact_linear": one_chi_l, "one_d_chi_exact_nonlinear": one_chi_n,
        "one_d_chi_retained_linear": retained_curvatures(one_l, one_gradient_l),
        "one_d_chi_retained_nonlinear": retained_curvatures(one_n, one_gradient_n)})
    curvature_summary = {}
    for kind in ("linear", "nonlinear"):
        curvature = physical_curvature_from_rotations(primary["R_"+kind], x)
        alternative_curvature = physical_curvature_from_rotations(alternative["R_"+kind], x)
        payload["three_d_chi_orientation_gradient_"+kind] = curvature["chi_gradient"]
        payload["three_d_chi_geodesic_interval_"+kind] = curvature["interval_chi_geodesic"]
        payload["three_d_pchip_chi_orientation_gradient_"+kind] = alternative_curvature["chi_gradient"]
        payload["curvature_interval_x"] = curvature["interval_x"]
        exact = one_chi_l if kind == "linear" else one_chi_n
        curvature_summary[kind] = {"one_d_exact_vs_3D_orientation_twist_proxy": field_metrics(
                exact[..., 0], curvature["chi_gradient"][..., 0], x, grid, fixed_scale=1./disc.length),
            "orientation_derivative_method_change": field_metrics(
                (curvature["chi_gradient"][:, :-1, 0]+curvature["chi_gradient"][:, 1:, 0])/2,
                curvature["interval_chi_geodesic"][..., 0], curvature["interval_x"], grid, fixed_scale=1./disc.length),
            "time_interpolation_twist_proxy_change": field_metrics(
                curvature["chi_gradient"][..., 0], alternative_curvature["chi_gradient"][..., 0], x, grid, fixed_scale=1./disc.length),
            "orthogonality_max_absolute": curvature["orthogonality_max_absolute"],
            "gradient_symmetric_defect_max": curvature["gradient_symmetric_defect_max"],
            "definition": curvature["definition"], "Phi_is_not_chi1": True,
            "recovery_and_orientation_differentiation_remain_qualified": True}
    payload.update({"one_d_R_correction": one_correction_R, "three_d_R_correction": fem_correction_R,
                    "one_d_R_evolution": one_evolution_R, "three_d_R_evolution": fem_evolution_R})
    payload.update({"three_d_"+kind: value for kind, value in primary.items()})
    payload.update({"three_d_pchip_"+kind: value for kind, value in alternative.items()})
    payload.update({"orientation_difference_"+kind: value for kind, value in orientation.items()})
    np.savez_compressed(output/"one_d_three_d_comparison.npz", **payload)
    summary = {"version": VERSION, "scientific_calls": 0, "time_grid": {"count": 201,
        "start": float(grid[0]), "end": float(grid[-1]), "normalized_end": float(grid[-1]/T1)},
        "metrics": metrics, "interpolation": interpolation, "physical_orientation": orientation_summary,
        "section_orientation_curvature": curvature_summary,
        "signal_planning": planning, "primary_time_interpolation": "linear", "alternative": "shape-preserving cubic PCHIP",
        "native_samples_and_interpolated_values_are_distinct": True,
        "zero_state_source": "confirmed static preload, not a native dynamic output frame",
        "physical_orientation_space_transfer": "full raw R, face R=I, linear material-x transfer then nearest SO(3)",
        "physical_orientation_time_transfer": "component interpolation of full R followed by nearest SO(3); no angle addition",
        "orientation_correction_definition": "C(t)=R_L(t).T@R_NL(t); evolving C=C(0).T@C(t); group diagnostic, not vector subtraction",
        "physical_model_accuracy_status": "NOT_AUTOMATICALLY_ASSIGNED",
        "three_d_temporal_certification": "PARTIAL_SINGLE_TIME_LEVEL",
        "effective_contraction": "QUALIFIED_PROXY_RECOVERY_SENSITIVITY_REMAINS",
        "no_phase_amplitude_or_time_fitting": True}
    _write(output/"one_d_three_d_comparison.json", summary)
    return summary


def analyze_mesh_comparison(bundle):
    """Read-only medium/fine comparison of the two actual dynamic pairs.

    The common output grid must already match exactly. This function never
    performs an additional transfer, calls a solver, or substitutes mesh
    changes for a continuum error bound. Orientation-derived twist curvature
    is retained as a qualified recovery/differentiation diagnostic.
    """
    bundle = Path(bundle)
    paths = {level: bundle/("comparison_"+level)/"one_d_three_d_comparison.npz" for level in ("medium", "fine")}
    if not all(path.exists() for path in paths.values()):
        return {"status": "NOT_RUN", "reason": "Both actual medium and fine comparison arrays are required",
                "scientific_calls": 0, "available_levels": [name for name, path in paths.items() if path.exists()]}
    arrays = {}
    for level, path in paths.items():
        with np.load(path, allow_pickle=False) as z:
            arrays[level] = {name: z[name].copy() for name in z.files}
    medium, fine = arrays["medium"], arrays["fine"]
    if not np.array_equal(medium["time"], fine["time"]) or not np.array_equal(medium["x"], fine["x"]):
        raise ValueError("Mesh comparison requires the same frozen physical time/material-x grid")
    times, x = fine["time"], fine["x"]
    metrics, interpolation, orientation = {}, {}, {}
    payload = {"time": times, "x": x}
    for kind in ("linear", "nonlinear", "correction", "evolution"):
        key = "three_d_"+kind
        difference = medium[key]-fine[key]
        payload[kind+"_difference"] = difference
        payload[kind+"_difference_L2"] = np.sqrt(np.trapezoid(difference**2, x=x, axis=1))
        payload[kind+"_difference_max"] = np.max(abs(difference), axis=1)
        metrics[kind], interpolation[kind] = {}, {}
        common_floor = 1e-10*max(float(np.sqrt(np.trapezoid(fine[key]**2, x=x, axis=1)).max()), 1e-30)
        for i, field in enumerate(FIELDS):
            fixed = .1 if i < 3 else 1.
            row = field_metrics(medium[key][..., i], fine[key][..., i], x, times,
                                fixed_scale=fixed, reference_floor=common_floor)
            if field == "c":
                row["qualification"] = "effective contraction proxy; not a direct M-H coordinate convergence test"
            if field in ("Phi", "psi", "theta"):
                row["qualification"] = "rotation-vector chart diagnostic; full physical R compared separately"
            metrics[kind][field] = row
            interpolation[kind][field] = {}
            for level, current in arrays.items():
                interpolation[kind][field][level] = field_metrics(current[key][..., i],
                    current["three_d_pchip_"+kind][..., i], x, times, fixed_scale=fixed)
    for kind in ("linear", "nonlinear", "correction", "evolution"):
        key = "three_d_R_"+kind
        angle = orientation_difference(medium[key], fine[key])
        payload["orientation_mesh_angle_"+kind] = angle
        orientation[kind] = {"maximum_principal_angle_radians": float(angle.max()),
            "max_time_L2_angle": float(np.sqrt(np.trapezoid(angle**2, x=x, axis=1)).max())}
    curvature, torsion = {}, {}
    for kind in ("linear", "nonlinear"):
        recovered = {level: physical_curvature_from_rotations(current["three_d_R_"+kind], x)
                     for level, current in arrays.items()}
        for level, current in recovered.items():
            payload[level+"_chi_gradient_"+kind] = current["chi_gradient"]
            payload[level+"_chi_geodesic_"+kind] = current["interval_chi_geodesic"]
        curvature[kind] = {"medium_fine_twist_proxy_change": field_metrics(
            recovered["medium"]["chi_gradient"][..., 0], recovered["fine"]["chi_gradient"][..., 0],
            x, times, fixed_scale=1./(x[-1]-x[0])), "levels": {}}
        for level, current in recovered.items():
            curvature[kind]["levels"][level] = {"twist_proxy_max": float(np.max(abs(current["chi_gradient"][..., 0]))),
                "geodesic_interval_twist_proxy_max": float(np.max(abs(current["interval_chi_geodesic"][..., 0]))),
                "orientation_derivative_method_change": field_metrics(
                    (current["chi_gradient"][:, :-1, 0]+current["chi_gradient"][:, 1:, 0])/2,
                    current["interval_chi_geodesic"][..., 0], current["interval_x"], times,
                    fixed_scale=1./(x[-1]-x[0])),
                "gradient_symmetric_defect_max": current["gradient_symmetric_defect_max"]}
            exact_key = "one_d_chi_exact_"+kind
            if exact_key in arrays[level]:
                curvature[kind]["levels"][level]["one_d_exact_vs_orientation_twist_proxy"] = field_metrics(
                    arrays[level][exact_key][..., 0], current["chi_gradient"][..., 0], x, times,
                    fixed_scale=1./(x[-1]-x[0]))
            else:
                curvature[kind]["levels"][level]["one_d_exact_vs_orientation_twist_proxy"] = {"status": "NOT_AVAILABLE"}
    for level, current in arrays.items():
        nonlinear_phi = current["three_d_nonlinear"][..., 3]
        linear_phi = current["three_d_linear"][..., 3]
        delta_phi = nonlinear_phi-linear_phi
        evolving_phi = delta_phi-delta_phi[0]
        floor = HISTORICAL_PLANAR_PHI_FLOOR if level == "medium" else HISTORICAL_FINE_PLANAR_PHI_FLOOR
        signal = float(np.max(abs(evolving_phi)))
        interpolation_change = interpolation["evolution"]["Phi"][level]["absolute_max"]
        p_change = metrics["evolution"]["Phi"]["absolute_max"]
        torsion[level] = {"nonlinear_rotation_vector_Phi_max": float(np.max(abs(nonlinear_phi))),
            "linear_rotation_vector_Phi_max": float(np.max(abs(linear_phi))),
            "nonlinear_minus_linear_Phi_max": float(np.max(abs(delta_phi))),
            "evolving_Phi_max": signal, "historical_planar_proxy_indicator": floor,
            "Phi_evolution_planning_diagnostic": planning_signal(signal, p_change,
                indicators=(floor, interpolation_change), change_label="observed_medium_fine_change"),
            "Phi_is_not_twist_curvature": True,
            "raw_3D_section_orientation_and_native_strain_not_interchangeable": True}
    # Separate full motion, nonlinear corrections, and evolving corrections.
    # No single universal PASS is inferred from small mesh changes.
    result = {"version": VERSION, "status": "MESH_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS",
        "scientific_calls": 0, "actual_common_time_range": [float(times[0]), float(times[-1])],
        "mesh_metrics": metrics, "interpolation_metrics": interpolation, "physical_orientation_mesh_changes": orientation,
        "section_orientation_curvature": curvature, "torsion_proxy_resolution": torsion,
        "mesh_change_is_not_a_continuum_error_bound": True,
        "single_3D_time_level_temporal_certification": "PARTIAL",
        "torsion_requires_both_finite_orientation_and_recovery_differentiation_qualification": True,
        "no_phase_amplitude_time_or_coefficient_fitting": True,
        "source_comparison_hashes": {level: _sha(path) for level, path in paths.items()}}
    np.savez_compressed(bundle/"mesh_comparison.npz", **payload)
    _write(bundle/"mesh_comparison.json", result)
    return result


class SavedShenReconstruction:
    """Basis-only adapter for saved coordinates; no action or solver assembly."""

    def __init__(self, case):
        from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
        metadata = cached_case(case)
        if metadata is None or metadata.get("fields") != list(FIELDS):
            raise ValueError("Canonical saved seven-field case required")
        request = metadata["request"]
        self.p, self.n, self.ndof = request["p"], request["p"]-1, 7*(request["p"]-1)
        self.length = request["length"]
        self.coefficients = rod.RodCoefficients(**request["coefficients"])
        self.slices = {field: slice(i*self.n, (i+1)*self.n) for i, field in enumerate(FIELDS)}
        self._raw_basis = lambda points, derivative: PlanarGalerkin._raw_basis(self, points, derivative)
        xi, weights = np.polynomial.legendre.leggauss(request["nq"])
        raw = self._raw_basis((xi+1)*self.length/2, 0)
        gram = raw.T@((weights*self.length/2)[:, None]*raw)
        p = self.coefficients
        masses = (p.m, p.m, p.m, p.jp+p.jb, p.jb, p.jp, p.jp)
        self.transforms = []
        for mass in masses:
            lower = np.linalg.cholesky(mass*gram)
            self.transforms.append(solve_triangular(lower.T, np.eye(self.n), lower=False)
                if request["whiten"] else np.eye(self.n))
        with np.load(Path(case)/"static_states.npz", allow_pickle=False) as z:
            for kind in ("linear", "nonlinear"):
                reproduced = np.concatenate([a@z["q_"+kind][self.slices[field]] for field, a in zip(FIELDS, self.transforms)])
                if np.max(abs(reproduced-z["raw_"+kind])) > 5e-12:
                    raise ValueError("Saved coefficient basis/scaling reproduction failed")

    def reconstruct_series(self, rows, points, derivative=0):
        rows = np.asarray(rows)
        if rows.ndim != 2 or rows.shape[1] != self.ndof:
            raise ValueError("Saved coordinate rows have the wrong seven-field dimension")
        raw = self._raw_basis(points, derivative)
        return np.stack([rows[:, self.slices[field]]@(raw@transform).T
                         for field, transform in zip(FIELDS, self.transforms)], axis=-1)


def _linear_saved_initial(case, initial, times):
    """Full stored spectral operator, zero initial speed, no new eigensolve."""
    with np.load(Path(case)/"linear_modes.npz", allow_pickle=False) as z:
        omega, vectors, mass = z["omega"], z["vectors"], z["M0"]
    if vectors.shape != (len(initial), len(initial)) or mass.shape != vectors.shape:
        raise ValueError("Full saved linear operator dimensions disagree")
    amplitudes = vectors.T@(mass@initial)
    angles = np.asarray(times)[:, None]*omega
    q = (np.cos(angles)*amplitudes)@vectors.T
    v = (-np.sin(angles)*(omega*amplitudes))@vectors.T
    if times[0] == 0:
        q[0], v[0] = initial, 0.
    return q, v


def write_initial_state_decomposition(cases, output, *, discs=None, T1=DEFAULT_T1, common_times=None):
    """Saved-operator decomposition; this does not revise a pre-FEM decision.

    N(qNL0)−L(qL0) = [L(qNL0)−L(qL0)] + [N(qNL0)−L(qNL0)].
    Mixed translations are decomposed using the same identity separately for
    joint and both isolated nonlinear initial states. Rotational coordinate
    components are additive chart diagnostics, never sums of physical R.
    """
    required = ("joint_p48", "joint_p64", "isolated_w_p64", "isolated_v_p64")
    if not all(name in cases for name in required):
        raise ValueError("Joint p48/p64 and isolated p64 saved cases required")
    optional = ("isolated_w_p48", "isolated_v_p48")
    if any(name in cases for name in optional) and not all(name in cases for name in optional):
        raise ValueError("Optional isolated p48 controls must be supplied together")
    names = required+optional if all(name in cases for name in optional) else required
    times = np.linspace(0., T1/4, 201) if common_times is None else np.asarray(common_times)
    if len(times) != 201 or times[0] != 0 or times[-1] > T1/4*(1+1e-12) or np.any(np.diff(times) <= 0):
        raise ValueError("The predeclared 201-point quarter-period grid is required")
    reconstructions = {} if discs is None else dict(discs)
    source, coordinates = {}, {}
    payload = {"time": times}
    for name in names:
        case = Path(cases[name]); metadata = cached_case(case)
        if metadata is None or metadata["status"] != "PASS":
            raise ValueError("Completed immutable saved case required: "+name)
        degree = 48 if name.endswith("p48") else 64
        if degree not in reconstructions:
            reconstructions[degree] = SavedShenReconstruction(case)
        disc = reconstructions[degree]
        if disc.ndof != metadata["ndof"]:
            raise ValueError("Saved case and reconstruction dimension disagree")
        nonlinear = evaluate_dense_records(case/"accepted_dense.npz", times)
        with np.load(case/"static_states.npz", allow_pickle=False) as z:
            initial_l, initial_n = z["q_linear"], z["q_nonlinear"]
        from_l = _linear_saved_initial(case, initial_l, times)
        from_n = _linear_saved_initial(case, initial_n, times)
        coordinates[name] = {"nonlinear_q": nonlinear[:, :disc.ndof], "nonlinear_v": nonlinear[:, disc.ndof:],
                             "linear_from_linear_q": from_l[0], "linear_from_linear_v": from_l[1],
                             "linear_from_nonlinear_q": from_n[0], "linear_from_nonlinear_v": from_n[1]}
        source[name] = {"path": str(case), "case_sha256": _sha(case/"case.json"),
                       "saved_modes_sha256": _sha(case/"linear_modes.npz"), "new_eigen_ODE_static_calls": 0}
    decomposed, summaries = {}, {}
    for degree in (48, 64):
        joint = coordinates["joint_p"+str(degree)]
        total = joint["nonlinear_q"]-joint["linear_from_linear_q"]
        initial = joint["linear_from_nonlinear_q"]-joint["linear_from_linear_q"]
        same = joint["nonlinear_q"]-joint["linear_from_nonlinear_q"]
        decomposed["joint_p"+str(degree)] = (total, initial, same)
        payload["joint_p"+str(degree)+"_total_q"] = total
        payload["joint_p"+str(degree)+"_propagated_initial_q"] = initial
        payload["joint_p"+str(degree)+"_same_initial_nonlinear_q"] = same
        payload["joint_p"+str(degree)+"_linear_from_nonlinear_q"] = joint["linear_from_nonlinear_q"]
        payload["joint_p"+str(degree)+"_same_initial_nonlinear_v"] = joint["nonlinear_v"]-joint["linear_from_nonlinear_v"]
        if all("isolated_"+field+"_p"+str(degree) in coordinates for field in ("w", "v")):
            w, v = (coordinates["isolated_"+field+"_p"+str(degree)] for field in ("w", "v"))
            mixed = joint["nonlinear_q"]-w["nonlinear_q"]-v["nonlinear_q"]
            mixed_initial = joint["linear_from_nonlinear_q"]-w["linear_from_nonlinear_q"]-v["linear_from_nonlinear_q"]
            mixed_same = mixed-mixed_initial
            decomposed["mixed_p"+str(degree)] = (mixed, mixed_initial, mixed_same)
            for label, rows in zip(("total_q", "propagated_initial_q", "same_initial_nonlinear_q"), decomposed["mixed_p"+str(degree)]):
                payload["mixed_p"+str(degree)+"_"+label] = rows
    for name, components in decomposed.items():
        degree = int(name[-2:]); disc = reconstructions[degree]
        x = np.linspace(0., disc.length, 801)
        gauss, weights = np.polynomial.legendre.leggauss(100)
        gauss, weights = (gauss+1)*disc.length/2, weights*disc.length/2
        count = 3 if name.startswith("mixed") else 7
        payload.setdefault("x", x)
        accum = {label: {"max": np.zeros(count), "L2": np.zeros(count), "obs": [], "snapshot": []}
                 for label in ("total", "propagated_initial", "same_initial_nonlinear")}
        observation_x = np.array((.25, .5, .5, .5, .25, .25, .25))*disc.length
        observation_index = np.array([int(np.argmin(abs(x-p))) for p in observation_x])[:count]
        identity_max = 0.
        for start in range(0, len(times), 64):
            stop = min(start+64, len(times))
            values = []
            for label, rows in zip(accum, components):
                physical = disc.reconstruct_series(rows[start:stop], x)[..., :count]
                physical_gauss = disc.reconstruct_series(rows[start:stop], gauss)[..., :count]
                values.append(physical)
                accum[label]["max"] = np.maximum(accum[label]["max"], np.max(abs(physical), axis=(0, 1)))
                accum[label]["L2"] = np.maximum(accum[label]["L2"], np.sqrt(np.einsum("tij,i,tij->tj", physical_gauss, weights, physical_gauss)).max(axis=0))
                accum[label]["obs"].append(physical[:, observation_index, np.arange(count)])
                for index in (0, 50, 100, 150, 200):
                    if start <= index < stop:
                        accum[label]["snapshot"].append(physical[index-start])
            identity_max = max(identity_max, float(np.max(abs(values[0]-values[1]-values[2]))))
        rows = {}
        for i, field in enumerate(FIELDS[:count]):
            char = max(float(accum["total"]["max"][i]), 1e-10*float(accum["total"]["max"].max()), 1e-30)
            rows[field] = {"total_max": float(accum["total"]["max"][i]),
                "total_max_time_L2": float(accum["total"]["L2"][i]),
                "propagated_initial_max": float(accum["propagated_initial"]["max"][i]),
                "same_initial_nonlinear_max": float(accum["same_initial_nonlinear"]["max"][i]),
                "propagated_initial_max_over_total_characteristic": float(accum["propagated_initial"]["max"][i]/char),
                "same_initial_nonlinear_max_over_total_characteristic": float(accum["same_initial_nonlinear"]["max"][i]/char),
                "characteristic_denominator": char,
                "max_ratios_are_not_additive_percentages": True}
        for label, current in accum.items():
            observations = np.concatenate(current["obs"])
            payload[name+"_"+label+"_observations"] = observations
            payload[name+"_"+label+"_snapshots"] = np.asarray(current["snapshot"])
            for i, field in enumerate(FIELDS[:count]):
                rows[field]["final_"+label+"_observation"] = float(observations[-1, i])
        summaries[name] = {"fields": rows, "physical_identity_max_absolute": identity_max,
            "coefficient_identity_max_absolute": float(np.max(abs(components[0]-components[1]-components[2]))),
            "observation_x": observation_x[:count].tolist()}
    result = {"version": VERSION, "status": "READ_ONLY_DECOMPOSITION_COMPLETE", "scientific_calls": 0,
        "operation_counts": {"saved_full_linear_operator_evaluations": 2*len(names),
            "basis_only_reconstruction_adapters": 0 if discs is not None else len(reconstructions),
            "new_eigen_ODE_static_calls": 0},
        "grid": {"count": 201, "start": float(times[0]), "end": float(times[-1])},
        "components": summaries, "sources": source,
        "new_excitation_or_initial_state_scientific_case": False,
        "exact_saved_full_linear_operators_no_modal_truncation": True,
        "evolving_NL_minus_L_is_not_same_initial_nonlinear_evolution": True,
        "rotational_components_are_generalized_coordinate_charts_not_additive_physical_orientations": True,
        "pre_FEM_decision_unchanged": True, "diagnostic_not_a_new_stage_C_gate": True,
        "p48_p64_comparison_of_components_requires_own_absolute_scales": True}
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/"initial_state_decomposition.npz", **payload)
    _write(output/"initial_state_decomposition.json", result)
    return result


def transfer_moment_origin(moment, force, old_origin, new_origin):
    """Same resultant, M_new=M_old+(old_origin−new_origin)×F."""
    moment, force, old_origin, new_origin = map(np.asarray, (moment, force, old_origin, new_origin))
    if any(a.shape != (3,) or not np.isfinite(a).all() for a in (moment, force, old_origin, new_origin)):
        raise ValueError("Finite global three-vectors required for moment origin transfer")
    return moment+np.cross(old_origin-new_origin, force)


def audit_static_torsional_supports(bundle, cases, *, discs=None, action_result=None):
    """Qualified independent static Mx evidence; no dynamic RF interpretation."""
    bundle = Path(bundle)
    root = Path(__file__).resolve().parents[2]
    action_result = (root/"results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json"
                     if action_result is None else Path(action_result))
    manifest = json.loads((action_result.parent/"manifest.json").read_text(encoding="utf8"))
    hashes = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
    if hashes.get(action_result.name) != _sha(action_result):
        raise ValueError("Frozen action artifact hash mismatch in static torsional audit")
    polynomial = rod.Polynomial.deserialize(json.loads(action_result.read_text(encoding="utf8"))["polynomials"]["V4"])
    endpoint_flux = polynomial.derivative("Phi_s").substitute({field: 0 for field in FIELDS})
    expected = rod.Polynomial.symbol("CT")*rod.Polynomial.symbol("Phi_s")
    if endpoint_flux != expected:
        raise ValueError("Accepted clamped endpoint torsional flux does not reduce to CT*Phi_s")
    source_1d, one = {}, {}
    for degree in (48, 64):
        case = Path(cases["joint_p"+str(degree)])
        metadata = cached_case(case)
        disc = SavedShenReconstruction(case) if discs is None else discs[degree]
        with np.load(case/"static_states.npz", allow_pickle=False) as z:
            states = {kind: z["q_"+kind].copy() for kind in ("linear", "nonlinear")}
        ct = metadata["request"]["coefficients"]["CT"]
        one[str(degree)] = {}
        for kind, coordinate in states.items():
            fields = disc.reconstruct_series(coordinate[None, :], np.array((0., disc.length)))[0]
            gradients = disc.reconstruct_series(coordinate[None, :], np.array((0., disc.length)), derivative=1)[0]
            if np.max(abs(fields)) > 1e-12:
                raise ValueError("Saved static state does not retain seven essential endpoint values")
            flux = np.array([endpoint_flux.evaluate({**disc.coefficients.values(), "Phi_s": row[3]}) for row in gradients])
            support = flux*np.array((-1., 1.))
            gauss, weights = np.polynomial.legendre.leggauss(100)
            gauss, weights = (gauss+1)*disc.length/2, weights*disc.length/2
            values = disc.reconstruct_series(coordinate[None, :], gauss)[0]
            qw, qv = metadata["loads"]
            body_mx = float(weights@(values[:, 1]*qv-values[:, 2]*qw))
            one[str(degree)][kind] = {"CT": ct, "endpoint_Phi_s": gradients[:, 3].tolist(),
                "LEFT_FIXED_Mx": float(support[0]), "RIGHT_FIXED_Mx": float(support[1]),
                "summed_support_Mx": float(support.sum()), "deformed_axis_body_Mx_diagnostic": body_mx,
                "summed_support_plus_body_Mx_diagnostic": float(support.sum()+body_mx),
                "body_Mx_definition": "integral (w*q_v−v*q_w) ds in the common global frame",
                "body_axis_torque_is_separate_diagnostic_not_equilibrium_gate": True,
                "linear_geometry_balance_is_not_tested_on_deformed_axis": kind == "linear",
                "endpoint_flux_source": "exact derivative of saved V4, all field endpoint values zero; no slope constraints"}
        source_1d[str(degree)] = {"case_sha256": _sha(case/"case.json"), "static_states_sha256": _sha(case/"static_states.npz")}
    comparisons = {}
    native_sources = {}
    config_path = bundle/"config.json"
    config = json.loads(config_path.read_text(encoding="utf8")) if config_path.exists() else {}
    fem1_bundle = root/config.get("source_bundles", {}).get("FEM1", {}).get(
        "path", "results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad")
    mesh_manifest = json.loads((fem1_bundle/"manifest.json").read_text(encoding="utf8"))
    mesh_hashes = mesh_manifest.get("artifact_hashes", mesh_manifest.get("artifacts", {}))
    for level in ("medium", "fine"):
        comparisons[level] = {}
        geometry = None
        for kind in ("linear", "nonlinear"):
            path = bundle/"FEM"/(level+"_"+kind)/"preload_equilibrium.json"
            if not path.exists():
                comparisons[level][kind] = {"status": "NOT_AVAILABLE"}
                continue
            native = json.loads(path.read_text(encoding="utf8"))
            if geometry is None:
                from scripts.analysis.solid_fem_single_rod_fixed_fixed import read_gmsh_inp_mesh_data
                mesh_path, audit_path = (fem1_bundle/"meshes"/level/name for name in ("solid_mesh.inp", "mesh_audit.json"))
                for source in (mesh_path, audit_path):
                    if mesh_hashes.get(source.relative_to(fem1_bundle).as_posix()) != _sha(source):
                        raise ValueError("Frozen source mesh/audit hash mismatch in support origin transfer")
                mesh = read_gmsh_inp_mesh_data(mesh_path)
                audit = json.loads(audit_path.read_text(encoding="utf8"))
                geometry = {face: np.mean([mesh.nodes[int(node)] for node in audit[key]], axis=0)
                            for face, key in (("LEFT_FIXED", "fixed_left_ids"), ("RIGHT_FIXED", "fixed_right_ids"))}
                native_sources[level+"_mesh"] = {"mesh_path": str(mesh_path), "mesh_sha256": _sha(mesh_path),
                                                 "audit_sha256": _sha(audit_path)}
            centers = {"LEFT_FIXED": np.array((0., 0., 0.)), "RIGHT_FIXED": np.array((disc.length, 0., 0.))}
            old_moments, corrections, new_moments = {}, {}, {}
            for face in ("LEFT_FIXED", "RIGHT_FIXED"):
                old = np.asarray(native["supports"][face]["moment_about_face_centroid"])
                force = np.asarray(native["supports"][face]["force"])
                corrected = transfer_moment_origin(old, force, geometry[face], centers[face])
                old_moments[face], new_moments[face] = old.tolist(), corrected.tolist()
                corrections[face] = (corrected-old).tolist()
            values = np.array([new_moments[face][0] for face in ("LEFT_FIXED", "RIGHT_FIXED")])
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite saved independent support moment")
            expected_values = np.array([one["64"][kind]["LEFT_FIXED_Mx"], one["64"][kind]["RIGHT_FIXED_Mx"]])
            scale = max(float(np.max(abs(values))), float(np.max(abs(expected_values))), 1e-30)
            comparisons[level][kind] = {"status": "READ_ONLY_COMPARISON", "LEFT_FIXED_Mx": float(values[0]),
                "RIGHT_FIXED_Mx": float(values[1]), "signed_1D_minus_3D": (expected_values-values).tolist(),
                "max_absolute_1D_3D_difference": float(np.max(abs(expected_values-values))),
                "nonzero_common_scale": scale,
                "relative_max_on_common_scale": float(np.max(abs(expected_values-values))/scale),
                "relative_difference_near_zero_linear_moment_is_not_validation": kind == "linear",
                "original_native_moment_vectors": old_moments, "origin_transfer_correction": corrections,
                "moments_about_geometric_face_centers": new_moments,
                "original_arithmetic_fixed_node_centroids": {face: point.tolist() for face, point in geometry.items()},
                "new_geometric_face_centers": {face: point.tolist() for face, point in centers.items()},
                "source": "DAT RF minus independently integrated consistent body load; original arithmetic fixed-node centroid moment transferred to exact axis/face center",
                "native_moment_about_face_centroid_label_means_arithmetic_fixed_node_mean": True,
                "origin_transfer_formula": "M_axis=M_node_mean+(node_mean−axis_center)×F_support",
                "not_inferred_from_global_balance": True}
            native_sources[level+"_"+kind] = {"path": str(path), "sha256": _sha(path)}
    result = {"version": VERSION, "status": "QUALIFIED_STATIC_TORSIONAL_SUPPORT_COMPARISON",
        "scientific_calls": 0, "endpoint_flux_exact_polynomial_identity": "PASS",
        "action_result_sha256": _sha(action_result), "action_manifest_sha256": _sha(action_result.parent/"manifest.json"),
        "one_d": one, "native_preload_comparisons": comparisons,
        "source_1D": source_1d, "source_native": native_sources,
        "dynamic_RF_as_independent_support_torque": "NOT_TESTED",
        "one_d_endpoint_flux_is_retained_action_generalized_moment": True,
        "moment_units": "unchanged normalized force times length",
        "no_explicit_distributed_local_torque_load": True,
        "dead_transverse_force_on_deflected_axis_can_have_global_axial_moment": True,
        "does_not_identify_or_validate_all_dynamic_torsion_or_warping_terms": True}
    _write(bundle/"static_torsional_support_audit.json", result)
    return result


def render_figures(bundle):
    """At most four main PDF/PNG figures from already saved comparison data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    bundle = Path(bundle)
    if not (bundle/"stage_b_comparison.npz").exists():
        return []
    meta = json.loads((bundle/"stage_b_comparison.json").read_text(encoding="utf8"))
    with np.load(bundle/"stage_b_comparison.npz", allow_pickle=False) as z:
        stage_b = {name: z[name].copy() for name in z.files}
    stage_c = None
    level = None
    for name in ("fine", "medium"):
        path = bundle/("comparison_"+name)/"one_d_three_d_comparison.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as z:
                stage_c = {key: z[key].copy() for key in z.files}
            level = name
            break
    T1 = meta["T1"]
    colors = plt.cm.viridis(np.linspace(.12, .88, 5))
    plots = []

    def finish(fig, name):
        fig.tight_layout()
        for suffix in ("pdf", "png"):
            path = bundle/(name+"."+suffix)
            fig.savefig(path, dpi=180, bbox_inches="tight")
            plots.append(path.name)
        plt.close(fig)

    def axes_style(axes):
        for ax in np.ravel(np.asarray(axes, dtype=object)):
            ax.grid(alpha=.22)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)

    if stage_c is None:
        x = stage_b["x"]; snapshot_times = stage_b["snapshot_times"]
        profiles_1d = stage_b["joint_p64_snapshots"]; profiles_3d = None
        source_label = "1D only; no independent 3D trajectory is available"
    else:
        x = stage_c["x"]; selected = np.linspace(0, len(stage_c["time"])-1, 5).astype(int)
        snapshot_times = stage_c["time"][selected]
        profiles_1d = stage_c["one_d_nonlinear"][selected]
        profiles_3d = stage_c["three_d_nonlinear"][selected]
        source_label = "1D and "+level+" 3D; displayed FEM profiles use labelled postprocessing time transfer"
    length = x[-1]-x[0]
    fig = plt.figure(figsize=(11.5, 7.8))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 1.))
    spatial = fig.add_subplot(grid[0, :], projection="3d")
    yview, zview = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    for j, (time_value, one) in enumerate(zip(snapshot_times, profiles_1d)):
        label = r"$t/T_1=$"+f"{time_value/T1:.3f}"
        spatial.plot((x+one[:, 0])/length, -one[:, 1]/length, -one[:, 2]/length, color=colors[j], lw=1.2)
        yview.plot((x+one[:, 0])/length, -one[:, 1]/length, color=colors[j], label=label)
        zview.plot((x+one[:, 0])/length, -one[:, 2]/length, color=colors[j])
        if profiles_3d is not None:
            three = profiles_3d[j]
            spatial.plot((x+three[:, 0])/length, -three[:, 1]/length, -three[:, 2]/length, color=colors[j], ls="--", lw=1.)
            yview.plot((x+three[:, 0])/length, -three[:, 1]/length, color=colors[j], ls="--")
            zview.plot((x+three[:, 0])/length, -three[:, 2]/length, color=colors[j], ls="--")
    # Equal physical unit scale in the three-dimensional view; transverse
    # projections below disclose their own axes instead of magnifying a rod.
    spatial.set(xlim=(0., 1.), ylim=(-.5, .5), zlim=(-.5, .5), xlabel="Global X/L", ylabel="Global Y/L", zlabel="Global Z/L")
    spatial.set_box_aspect((1., 1., 1.)); spatial.view_init(elev=19, azim=-61)
    spatial.set_title("Current rod axis, equal physical scale; no displacement magnification", fontsize=10)
    yview.set(xlabel="Global X/L", ylabel="Global Y/L", title="Projection on X–Y")
    zview.set(xlabel="Global X/L", ylabel="Global Z/L", title="Projection on X–Z")
    yview.legend(fontsize=8, loc="best")
    if profiles_3d is not None:
        zview.legend(handles=(Line2D([], [], color="k", label="1D", ls="-"),
                              Line2D([], [], color="k", label="3D "+level, ls="--")), fontsize=9)
    axes_style((yview, zview))
    fig.suptitle("Spatial free motion on 0…0.25T₁\n"+source_label, fontsize=11)
    finish(fig, "spatial_free_motion")

    fig, axes = plt.subplots(2, 2, figsize=(11.3, 7.5))
    if stage_c is None:
        times = stage_b["times"]
        one_linear = stage_b["joint_p64_linear_observations"]
        one_nonlinear = stage_b["joint_p64_observations"]
        three_linear = three_nonlinear = None
    else:
        times = stage_c["time"]
        observation_indices = np.array([int(np.argmin(abs(stage_c["x"]-p))) for p in (.25, .5, .5, .5, .25, .25, .25)])
        component_indices = np.arange(7)
        one_linear = stage_c["one_d_linear"][:, observation_indices, component_indices]
        one_nonlinear = stage_c["one_d_nonlinear"][:, observation_indices, component_indices]
        three_linear = stage_c["three_d_linear"][:, observation_indices, component_indices]
        three_nonlinear = stage_c["three_d_nonlinear"][:, observation_indices, component_indices]
    for ax, field, position in zip(axes.ravel(), (1, 2, 5, 4), ("L/2", "L/2", "L/4", "L/4")):
        ax.plot(times/T1, one_linear[:, field], color="tab:blue", label="1D linear")
        ax.plot(times/T1, one_nonlinear[:, field], color="tab:orange", label="1D nonlinear")
        if three_linear is not None:
            ax.plot(times/T1, three_linear[:, field], color="tab:blue", ls="--", label="3D linear "+level)
            ax.plot(times/T1, three_nonlinear[:, field], color="tab:orange", ls="--", label="3D nonlinear "+level)
        name = ("u", "w", "v", "Φ", "ψ", "θ", "c")[field]
        ax.set(xlabel="t/T₁", ylabel=name+("/L" if field < 3 else " [rad]"), title=name+"("+position+",t)")
    axes[0, 0].legend(fontsize=8)
    axes_style(axes)
    fig.suptitle("Two bending components and section rotation coordinates\n"+source_label, fontsize=11)
    finish(fig, "two_bending_components")

    fig, axes = plt.subplots(2, 2, figsize=(11.3, 7.5))
    if stage_c is None:
        delta, evolution = stage_b["delta_midspan"], stage_b["evolution_midspan"]
        fem_delta = fem_evolution = None
    else:
        midpoint = int(np.argmin(abs(stage_c["x"]-.5*length)))
        delta, evolution = stage_c["one_d_correction"][:, midpoint], stage_c["one_d_evolution"][:, midpoint]
        fem_delta, fem_evolution = stage_c["three_d_correction"][:, midpoint], stage_c["three_d_evolution"][:, midpoint]
    for j, field in enumerate((1, 2)):
        for row, one, three, label in ((0, delta, fem_delta, "NL−L"), (1, evolution, fem_evolution, "(NL−L)(t)−(NL−L)(0)")):
            ax = axes[row, j]
            ax.plot(times/T1, one[:, field], label="1D", color="tab:blue")
            if three is not None:
                ax.plot(times/T1, three[:, field], label="3D "+level, color="tab:orange", ls="--")
            ax.axhline(0., color="k", lw=.5)
            ax.set(xlabel="t/T₁", ylabel=("Δ" if row == 0 else "δₑᵥₒₗ ")+FIELDS[field]+"/L", title=label+", "+FIELDS[field]+"(L/2)")
    axes[0, 0].legend(fontsize=9); axes_style(axes)
    fig.suptitle("Nonlinear correction and its evolution; own static initial states retained", fontsize=11)
    finish(fig, "nonlinear_spatial_response")

    fig, axes = plt.subplots(2, 2, figsize=(12., 8.))
    bt = stage_b["times"]/T1
    for i, field in enumerate(("u", "w", "v")):
        axes[0, 0].plot(bt, stage_b["mix_evolution_observations"][:, i], label=field+("(L/4)" if i == 0 else "(L/2)"))
    axes[0, 0].set(xlabel="t/T₁", ylabel="Mixed evolving displacement/L", title="Joint−isolated w−isolated v, minus its initial value")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(bt, stage_b["joint_p64_observations"][:, 3], label="1D Φ coordinate, L/2", color="tab:blue")
    if stage_c is not None:
        midpoint = int(np.argmin(abs(stage_c["x"]-.5*length)))
        axes[0, 1].plot(times/T1, stage_c["three_d_nonlinear"][:, midpoint, 3], label="3D rotation-vector proxy, L/2", color="tab:orange", ls="--")
    floor = HISTORICAL_PLANAR_PHI_FLOOR
    axes[0, 1].axhspan(-floor, floor, color="gray", alpha=.17,
                       label="Historical planar proxy level; indicator, not error bound")
    axes[0, 1].set(xlabel="t/T₁", ylabel="Φ [rad]", title="Twist-related coordinate diagnostic; resolution requires separate evidence")
    axes[0, 1].legend(fontsize=7)
    labels, ratios, colors_status = [], [], []
    for part in ("q", "velocity"):
        for name, row in meta["all14_spatial"][part]["fields"].items():
            labels.append(name+("ₜ" if part == "velocity" else ""))
            ratios.append(max(row["relative_max"], row["relative_L2"])/row["tolerance"])
            colors_status.append("tab:blue" if row["status"] == "PASS" else "tab:red")
    axes[1, 0].bar(np.arange(14), np.maximum(ratios, 1e-6), color=colors_status)
    axes[1, 0].axhline(1., color="k", ls="--", lw=.8)
    axes[1, 0].set_yscale("log"); axes[1, 0].set_xticks(np.arange(14), labels, rotation=45, fontsize=8)
    axes[1, 0].set(ylabel="Observed relative change / inherited tolerance", title="All fourteen p48→p64 spatial checks")
    axes[1, 1].plot(bt, stage_b["axis_linear_RMS_distance"], label="linear axis")
    axes[1, 1].plot(bt, stage_b["axis_nonlinear_RMS_distance"], label="nonlinear axis")
    axes[1, 1].set(xlabel="t/T₁", ylabel="Best-plane RMS distance/L", title="Instantaneous nonplanarity; distinct from nonlinear coupling")
    axes[1, 1].legend(fontsize=8)
    axes_style((axes[0, 0], axes[0, 1], axes[1, 1]))
    axes[1, 0].grid(axis="y", alpha=.22)
    control_label = ("Isolated controls have no independent p check" if meta["isolated_p64_controls_have_no_independent_p_control"]
                     else "Matching isolated p48/p64 mixed-response controls completed")
    fig.suptitle("Spatial coupling and numerical qualifications\n"+control_label+"; 3D temporal certification remains partial", fontsize=11)
    finish(fig, "spatial_coupling_diagnostics")
    _write(bundle/"figure_data_provenance.json", {"figures": plots, "scientific_calls": 0,
        "source_stage_b_comparison_sha256": _sha(bundle/"stage_b_comparison.npz"),
        "3D_level": level, "no_displacement_magnification_or_independent_curve_normalization": True,
        "three_dimensional_view_uses_equal_physical_scale": True,
        "source_3D_samples_are_explicitly_postprocessed": stage_c is not None})
    return plots
