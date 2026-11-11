"""Publication figures from accepted reflection-class CSV and crossing fields only.

No basis evaluations or spectral searches. Optional branch continuation uses
only archived physical fields and the assignment helpers of the existing tracker.
The input contract differs from the sorted-map publication renderer: a verified
class spectrum plus a separately verified double-eigenvalue field archive.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.lib import frequency_geometry_style as line_geometry
MAPS = ROOT / "results/rl_full_rr_maps"
ORDER = ["U", "W", "Psi", "N", "Q", "M"]
COLORS = {"+": "#1f77b4", "-": "#ff7f0e"}
BRANCH_COLORS = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                 "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")


def tracking_helpers():
    # The library also exports baseline spectral routines: none are called here.
    for path in (ROOT, ROOT / "src"):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    from scripts.lib import analytic_branch_tracking
    return analytic_branch_tracking


def mass_features(arms, beta_deg, contract):
    """Euclidean embedding of the FULL inertia product on saved Gauss nodes.

    At different angles, compare global translations, the common normal rotation
    and axial derivatives on corresponding local material coordinates. The four
    legacy component keys are merely storage slots for weighted feature blocks;
    they are NOT the EB displacement metric or reconstructed EB components.
    """
    h = np.sqrt(12)*contract["l"]*contract["epsilon"]
    area, inertia = h*h, h**4/12
    rho, extra = contract["rho"], contract["rho"]*contract["nu"]**2*2*inertia
    b = np.deg2rad(beta_deg)
    tangents = (np.array([1., 0.]), np.array([-np.cos(b), -np.sin(b)]))
    normals = (np.array([0., -1.]), np.array([-np.sin(b), np.cos(b)]))
    components = {}
    for a, t, n, side in zip(arms, tangents, normals, ("left", "right")):
        y, weights = a["y"], a["weights"]
        if np.any(weights <= 0):
            raise ValueError("Nonpositive saved quadrature weights")
        d = y[:, 0, None]*t + y[:, 1, None]*n
        components["u_"+side] = (np.sqrt(weights*rho*area)[:, None]*d).ravel()
        components["w_"+side] = np.r_[np.sqrt(weights*rho*inertia)*y[:, 2],
                                      np.sqrt(weights*extra)*a["dy"][:, 0]]
    norm2 = sum(float(np.dot(v, v)) for v in components.values())
    if not np.isfinite(norm2) or norm2 <= 0:
        raise ValueError("Nonpositive physical mode mass")
    np.testing.assert_allclose(norm2, mass_inner(arms, arms, contract), rtol=3e-15)
    return components


def continue_class(previous, candidates, sorted_lambdas, step):
    """Existing one-to-one assignment, restricted to ONE reflection class."""
    tracker = tracking_helpers()
    diagnostics = tracker.assignment_diagnostics(previous_points=previous,
        current_states=candidates, sorted_lambdas=sorted_lambdas,
        freq_weight=0., shape_metric="full")
    # Frequency is not the identity criterion. Check insensitivity to the legacy
    # optional frequency penalty without accepting any different assignment.
    for penalty in (.03, .1):
        other = tracker.assignment_diagnostics(previous_points=previous,
            current_states=candidates, sorted_lambdas=sorted_lambdas,
            freq_weight=penalty, shape_metric="full")
        if not np.array_equal(other.assignment, diagnostics.assignment):
            raise ValueError("Branch assignment depends on frequency penalty")
    margins = []
    for row, col in enumerate(diagnostics.assignment):
        selected = diagnostics.mac[row, col]
        margin = selected - np.max(np.delete(diagnostics.mac[row], col))
        # Existing tracker warning threshold; the actual data and sensitivity
        # distribution are reported, not interpreted as a universal tolerance.
        if not np.isfinite(selected) or selected < tracker.DEFAULT_MAC_WARNING_THRESHOLD or margin <= 0:
            raise ValueError("Ambiguous or low-overlap saved branch continuation")
        margins.append(float(margin))
    points = tracker.build_next_points_from_assignment(source_points=previous,
        current_states=candidates, sorted_lambdas=sorted_lambdas, diagnostics=diagnostics,
        step_type="saved_beta", step_index=step,
        mac_warning_threshold=tracker.DEFAULT_MAC_WARNING_THRESHOLD, shape_metric="full")
    return points, margins


def archived_branches(rows, contract, verified, spectrum_rows):
    """Start with the lowest eight at beta=mu=0; never re-sort their identities."""
    tracker = tracking_helpers()
    angles = sorted({float(r["beta_deg"]) for r in rows})
    if len(angles) != 197 or (angles[0], angles[-1]) != (0., 90.):
        raise ValueError("Expected the complete symmetric angular archive")
    source = {(float(r["beta_deg"]), int(r["rank"])): r for r in spectrum_rows}
    previous, curves, initial, grids = {}, {i: [] for i in range(1, 9)}, {}, None
    minimum_gap = {c: float("inf") for c in ("+", "-")}
    for step, beta in enumerate(angles):
        group = sorted([r for r in rows if float(r["beta_deg"]) == beta],
                       key=lambda r: int(r["original_sorted_position"]))
        if [int(r["original_sorted_position"]) for r in group] != list(range(1, 10)):
            raise ValueError("Classification does not reconstruct the nine-position prefix")
        paths = {r["source_npz"] for r in group}
        if len(paths) != 1 or next(iter(paths)) not in verified:
            raise ValueError("Unverified classified mode archive")
        archive = ROOT / next(iter(paths))
        if sha(archive) != verified[next(iter(paths))]:
            raise ValueError("Classified NPZ changed")
        states = {"+": [], "-": []}
        by_position = {}
        with np.load(archive, allow_pickle=False) as data:
            for row in group:
                cls, position = row["symmetry_class"], int(row["original_sorted_position"])
                if (cls not in states or row["classification_status"] != {"+": "PURE_PLUS", "-": "PURE_MINUS"}[cls]
                        or row["prefix_complete_below"] != "True" or row["shape_key"] != f"mode{position}"):
                    raise ValueError("Unconfirmed classified mode")
                original = source[(beta, position)]
                if (float(row["Lambda"]) != float(original["Lambda"])
                        or sha(ROOT / original["record"]) != original["record_sha256"]):
                    raise ValueError("Classified frequency does not match the source spectrum")
                arms = [{key: data[f'{row["shape_key"]}_arm{i}_{key}']
                         for key in ("x", "weights", "y", "dy")} for i in (1, 2)]
                current_grids = [a[key] for a in arms for key in ("x", "weights")]
                if grids is None:
                    grids = [v.copy() for v in current_grids]
                if any(not np.array_equal(a, b) for a, b in zip(grids, current_grids)):
                    raise ValueError("Different material/quadrature grids: no interpolation allowed")
                components = mass_features(arms, beta, contract)
                state = tracker.AnalyticModeState(epsilon=contract["epsilon"], beta=np.deg2rad(beta), mu=0.,
                    current_sorted_index=position, Lambda=float(row["Lambda"]), coeff=np.empty(0),
                    components=components,
                    shape_vector=tracker.shape_vector_from_components(components, shape_metric="full"),
                    smallest_singular_value=np.nan, singular_value_ratio=np.nan)
                states[cls].append(state)
                by_position[position] = row
        sorted_lambdas = [float(r["Lambda"]) for r in group]
        for cls in ("+", "-"):
            minimum_gap[cls] = min(minimum_gap[cls], float(np.diff([s.Lambda for s in states[cls]]).min()))
            if step == 0:
                chosen = [s for s in states[cls] if s.current_sorted_index <= 8]
                points = [tracker.make_branch_point(state=s,
                    branch_id=tracker.branch_id_from_base_sorted_index(s.current_sorted_index, prefix="reflection"),
                    base_sorted_index=s.current_sorted_index, mac_to_previous=1., relative_lambda_jump=0.,
                    sorted_lambdas=sorted_lambdas, step_type="base", step_index=0,
                    mac_warning_threshold=tracker.DEFAULT_MAC_WARNING_THRESHOLD) for s in chosen]
                margins = [None]*len(points)
            else:
                points, margins = continue_class(previous[cls], states[cls], sorted_lambdas, step)
            for point, margin in zip(points, margins):
                row = by_position[point.current_sorted_index]
                branch = point.base_sorted_index
                if step == 0:
                    initial[branch] = dict(symmetry_class=cls, class_rank=int(row["class_rank"]))
                curves[branch].append(dict(row, branch_id=point.branch_id, base_sorted_index=branch,
                    mass_overlap=point.mac_to_previous, assignment_margin=margin,
                    plot_source="class_spectrum.csv"))
            previous[cls] = points
    overlaps = [r["mass_overlap"] for rs in curves.values() for r in rs[1:]]
    margins = [r["assignment_margin"] for rs in curves.values() for r in rs[1:]]
    return curves, dict(semantics="tracked_branches", base=dict(beta_deg=0., mu=0., positions=list(range(1, 9))),
        initial_classes=initial, modes_used=1773, tracked_nodes=1576, intervals_checked=1568,
        minimum_mass_overlap=min(overlaps), minimum_MAC_squared=min(overlaps)**2,
        minimum_assignment_margin=min(margins), minimum_sampled_within_class_gap=minimum_gap,
        overlap_quantiles=np.quantile(overlaps, [0, .01, .5, .99, 1]).tolist(),
        thresholds_all_accepted={str(t): bool(min(overlaps) >= t) for t in (.8, .95, .99)},
        frequency_penalty_assignments_unchanged=[0., .03, .1],
        class_rank_changes=sum(int(r["class_rank"]) != initial[k]["class_rank"] for k, rs in curves.items() for r in rs),
        tracker="scripts/lib/analytic_branch_tracking.py: assignment_diagnostics/build_next_points_from_assignment",
        inertia="global translations + section rotation + full Rayleigh-Love axial derivative; saved Gauss weights",
        spectral_calls=False, missing_nodes=0)


def insert_branch_local_roots(curves, local, solution):
    """Add only previously computed first-class roots; preserve branch colors/IDs."""
    result = {k: [dict(r) for r in rs] for k, rs in curves.items()}
    if not solution["first_class_roots_certified"]:
        raise ValueError("Uncertified crossing")
    for branch, cls in ((1, "+"), (2, "-")):
        template = result[branch][0]
        if template["symmetry_class"] != cls or int(template["class_rank"]) != 1:
            raise ValueError("The first pair is not the expected class pair")
        existing = {float(r["beta_deg"]): r for r in result[branch]}
        for row in local[cls]:
            beta = float(row["beta_deg"])
            if beta in existing or beta == solution["beta_deg"]:
                expected = existing[beta]["Lambda"] if beta in existing else solution["Lambda"]
                np.testing.assert_allclose(float(row["Lambda"]), float(expected), rtol=3e-12, atol=0)
                continue
            result[branch].append(dict(beta_deg=beta, Lambda=row["Lambda"], symmetry_class=cls,
                class_rank=1, branch_id=template["branch_id"], base_sorted_index=branch,
                plot_source="block_roots.csv"))
        result[branch].append(dict(beta_deg=solution["beta_deg"], Lambda=solution["Lambda"],
            symmetry_class=cls, class_rank=1, branch_id=template["branch_id"], base_sorted_index=branch,
            plot_source="crossing_solution.json"))
        result[branch].sort(key=lambda r: float(r["beta_deg"]))
    return result


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def verify_manifest(folder, root=ROOT):
    """Check every recorded data artifact; report historical writing-skill drift."""
    manifest = read(folder / "manifest.json")
    if manifest["status"] != "COMPLETE":
        raise ValueError("Incomplete source archive")
    checked, editorial_drift = {}, {}
    pairs = [(root / p, h) for p, h in manifest["source_sha256"].items()]
    pairs += [(folder / p, h) for p, h in manifest["artifacts"].items()]
    pairs += [(root / p, h) for p, h in manifest.get("implementation_sha256", {}).items()]
    for path, expected in pairs:
        actual = sha(path)
        name = path.resolve().relative_to(root.resolve()).as_posix()
        if actual != expected:
            if name.startswith(".agents/skills/dissertation-writing-ru/"):
                editorial_drift[name] = dict(recorded=expected, current=actual)
                continue
            raise ValueError(f"SHA256 mismatch: {name}")
        checked[name] = actual
    return manifest, dict(verified_sha256=checked, historical_editorial_drift=editorial_drift)


def first_classes(rows, *, local=False):
    result = {}
    for cls in ("+", "-"):
        selected = [r for r in rows if r["symmetry_class"] == cls and int(r["class_rank"]) == 1]
        selected.sort(key=lambda r: float(r["beta_deg"]))
        angles = [float(r["beta_deg"]) for r in selected]
        if not selected or len(set(angles)) != len(angles):
            raise ValueError("Missing or duplicate first class roots")
        for r in selected:
            if local:
                valid = r["first_class_certified"] == "True" and int(r["nullity"]) == 1
            else:
                valid = (r["classification_status"] == ("PURE_PLUS" if cls == "+" else "PURE_MINUS")
                         and r["prefix_complete_below"] == "True")
            if not valid or not np.isfinite(float(r["Lambda"])):
                raise ValueError("Unconfirmed first class root")
        result[cls] = selected
    if [r["beta_deg"] for r in result["+"]] != [r["beta_deg"] for r in result["-"]]:
        raise ValueError("Class angular grids differ")
    return result


def validate_guard_record(record, rows, beta):
    """Read the accepted count/quality certificate; never recompute a count."""
    c = record["contract"]
    expected = dict(epsilon=.005, nu=.3, aspect=1., E=1., rho=1., l=1.,
                    ks=5/6, variant="RL_FULL", mu=0., tau=0.)
    if c["model"] != expected or c["K_plot"] != 8 or c["K_guard"] != 9:
        raise ValueError("Wrong physical contract or guard prefix")
    case = c["case"]
    if case["k_joint"] != "RIGID" or case["k_end"] != "RIGID" or case["mu"] != 0 or case["tau"] != 0:
        raise ValueError("Wrong joint, clamp or symmetric geometry")
    if (record["beta_deg"] != beta or record["mu"] != 0 or record["tau"] != 0
            or c["spectrum_semantics"] != "sorted_positions"
            or c["coordinates"] != "outer clamp to joint"
            or c["state_order"] != ["u", "w", "psi", "N", "Q", "M"]):
        raise ValueError("Wrong point coordinates or spectrum semantics")
    if case["beta"] != beta:  # The production case contract stores degrees.
        raise ValueError("Case angle differs from the saved degree coordinate")
    np.testing.assert_allclose(c["time_reference"], 200., atol=6e-14, rtol=0)
    roots = record["roots"]
    if [r["sorted_position"] for r in roots] != list(range(1, 10)):
        raise ValueError("Nine saved frequency positions are required")
    search = record["search"]
    if (record["status"] != "PASS" or search["status"] != "PASS" or search["lower_count"] != 0
            or search["requested_guard"] != 9 or search["computed_ranks"] != 9
            or search["guard_right_count"] != 9 or search["unresolved_intervals"]):
        raise ValueError("Unresolved saved guard/count certificate")
    if (record["quality"]["status"] != "PASS"
            or len(record["quality"]["gates"]) != 9
            or not all(all(g.values()) for g in record["quality"]["gates"])):
        raise ValueError("Saved physical quality checks failed")
    if sorted(int(r["rank"]) for r in rows) != list(range(1, 10)):
        raise ValueError("Missing or duplicate CSV frequency position")
    for row, root in zip(sorted(rows, key=lambda r: int(r["rank"])), roots):
        if row["status"] != "PASS" or any(float(row[k]) != root[k] for k in ("Lambda", "omega", "Omega")):
            raise ValueError("CSV/point frequency mismatch")
        if root["role"] != ("plot" if root["sorted_position"] <= 8 else "guard"):
            raise ValueError("Incorrect saved rank role")
    if np.any(np.diff([r["Lambda"] for r in roots]) < 0):
        raise ValueError("Unsorted frequency positions")
    if search["guard_right_omega"] <= roots[-1]["omega"]:
        raise ValueError("Guard endpoint does not enclose the ninth position")


def sorted_spectrum(folder, map_name):
    """Eight curves from spectrum.csv, with the ninth position checked in JSON/NPZ."""
    manifest = read(folder / "manifest.json")
    if manifest["status"] != "PASS" or sha(folder / "spectrum.csv") != manifest["artifacts"]["spectrum.csv"]:
        raise ValueError("Unverified spectrum.csv")
    spec = next(m for m in manifest["identity"]["maps"] if m["name"] == map_name)
    if spec["mode"] != "BETA" or spec["fixed"] != dict(mu=0., tau=0.):
        raise ValueError("A symmetric angular map is required")
    policy = manifest["identity"]["local_policies"][map_name]
    if (policy["K_plot"], policy["K_guard"], policy["spectrum_semantics"]) != (8, 9, "sorted_positions"):
        raise ValueError("Unexpected map policy")
    rows = [r for r in table(folder / "spectrum.csv") if r["map"] == map_name]
    points = [p for p in manifest["points"] if p["map"] == map_name]
    if (len(points) != 197 or len(rows) != 197*9
            or sorted(p["beta_deg"] for p in points) != spec["grid"]):
        raise ValueError("Incomplete symmetric angular map")
    verified, guard_checks = {}, []
    for point in points:
        beta = point["beta_deg"]
        subset = [r for r in rows if float(r["beta_deg"]) == beta]
        for row in subset:
            if (row["record"] != point["record"] or row["record_sha256"] != point["record_sha256"]
                    or float(row["mu"]) != 0 or float(row["tau"]) != 0 or float(row["aspect"]) != 1):
                raise ValueError("CSV/manifest point mismatch")
        path = ROOT / point["record"]
        if point["status"] != "PASS" or sha(path) != point["record_sha256"]:
            raise ValueError("Point JSON checksum/status mismatch")
        verified[point["record"]] = sha(path)
        record = read(path)
        validate_guard_record(record, subset, beta)
        shape = path.parent.parent / record["shape_file"]
        if sha(shape) != record["shape_sha256"]:
            raise ValueError("Point NPZ checksum mismatch")
        verified[shape.relative_to(ROOT).as_posix()] = sha(shape)
        with np.load(shape, allow_pickle=False) as data:
            for rank in range(1, 10):
                for arm in (1, 2):
                    for field in ("y", "dy", "ends"):
                        value = data[f"mode{rank}_arm{arm}_{field}"]
                        if value.shape[1] != 6 or not np.all(np.isfinite(value)):
                            raise ValueError("Missing/nonfinite physical fields")
        guard_checks.append(dict(beta_deg=beta, positions=len(record["roots"]),
                                 lower_count=record["search"]["lower_count"],
                                 guard_right_count=record["search"]["guard_right_count"],
                                 quality=record["quality"]["status"]))
    curves = {rank: sorted([r for r in rows if int(r["rank"]) == rank], key=lambda r: float(r["beta_deg"]))
              for rank in range(1, 9)}
    return curves, dict(source=folder.relative_to(ROOT).as_posix(), map=map_name,
                         verified_sha256=verified, point_guard_checks=guard_checks,
                         raw_csv_positions=len(rows), plotted_csv_positions=197*8,
                         semantics="sorted_positions", count_recomputed=False)


def insert_saved_crossing(curves, solution):
    """Add two already certified coincident roots, never interpolate missing ranks."""
    if not solution["first_class_roots_certified"]:
        raise ValueError("The crossing is not certified as the first pair")
    result = {rank: [dict(r, plot_source="spectrum.csv") for r in rows] for rank, rows in curves.items()}
    beta = solution["beta_deg"]
    for rank in (1, 2):
        if any(float(r["beta_deg"]) == beta for r in result[rank]):
            raise ValueError("A full-spectrum point already exists at the crossing")
        result[rank].append(dict(beta_deg=beta, Lambda=solution["Lambda"], omega=solution["omega"],
                                 rank=rank, plot_source="crossing_solution.json"))
        result[rank].sort(key=lambda r: float(r["beta_deg"]))
    return result


def branch_behaviour(classes):
    result = {}
    for cls, rows in classes.items():
        values = np.array([float(r["Lambda"]) for r in rows])
        result[cls] = dict(minimum=float(values.min()), maximum=float(values.max()),
                           relative_span_percent=float((values.max()-values.min())/values.min()*100),
                           minimum_increment=float(np.diff(values).min()),
                           maximum_increment=float(np.diff(values).max()))
    if result["+"]["minimum_increment"] <= 0 or result["-"]["relative_span_percent"] >= .02:
        raise ValueError("Saved class spectra do not support the requested panel descriptions")
    return result


def mass_inner(left, right, contract):
    """Saved Gauss weights, including rotation and full Rayleigh--Love inertia."""
    h = np.sqrt(12) * contract["l"] * contract["epsilon"]
    area, inertia = h**2, h**4 / 12
    rho, nu = contract["rho"], contract["nu"]
    extra = rho * nu**2 * 2 * inertia
    total = 0.0
    for a, b in zip(left, right):
        if not np.array_equal(a["x"], b["x"]) or not np.array_equal(a["weights"], b["weights"]):
            raise ValueError("Incompatible saved quadrature grids")
        y, z = a["y"], b["y"]
        total += np.dot(a["weights"], rho * area * (y[:, 0]*z[:, 0] + y[:, 1]*z[:, 1])
                        + rho * inertia * y[:, 2]*z[:, 2] + extra * a["dy"][:, 0]*b["dy"][:, 0])
    return float(total)


def coordinates(arms, beta_deg, length):
    """Local-to-global convention from the accepted audit, using saved endpoints."""
    b = np.deg2rad(beta_deg)
    tangents = (np.array([1., 0.]), np.array([-np.cos(b), -np.sin(b)]))
    normals = (np.array([0., -1.]), np.array([-np.sin(b), np.cos(b)]))
    scenes = []
    for arm, tangent, normal in zip(arms, tangents, normals):
        x = np.r_[0., arm["x"], length]
        if np.any(np.diff(x) <= 0):
            raise ValueError("Saved interior coordinates must be strictly ordered")
        y = np.vstack((arm["ends"][0], arm["y"], arm["ends"][1]))
        scenes.append(((x - length)[:, None]*tangent, y[:, 0, None]*tangent + y[:, 1, None]*normal))
    return scenes


def presentation_scenes(modes, beta_deg, length, amplitude):
    """One sign per mode, ONE displacement multiplier for both unit-mass modes."""
    scenes, signs = {}, {}
    for cls, arms in modes.items():
        raw = coordinates(arms, beta_deg, length)
        vertical = raw[0][1][:, 1]
        signs[cls] = 1 if vertical[np.argmax(abs(vertical))] >= 0 else -1
        scenes[cls] = [(x, signs[cls]*d) for x, d in raw]
    maximum = max(float(np.linalg.norm(d, axis=1).max()) for scene in scenes.values() for _, d in scene)
    if maximum <= 0:
        raise ValueError("Zero displacement fields")
    return scenes, dict(signs=signs, maximum_saved_displacement=maximum,
                        common_displacement_multiplier=amplitude*length/maximum,
                        displayed_maximum_over_length=amplitude,
                        per_mode_amplitude_normalization=False)


def load_fields(folder, contract, solution):
    checks = read(folder / "crossing_form_checks.json")
    if checks["block_nullities"] != [1, 1] or checks["full_eigenspace_dimension"] != 2:
        raise ValueError("Unexpected saved crossing multiplicity")
    modes = {}
    with np.load(folder / "crossing_forms.npz", allow_pickle=False) as data:
        if data["state_order"].tolist() != ORDER:
            raise ValueError("Unexpected physical state ordering")
        for name in ("beta_deg", "omega", "Lambda"):
            if float(data[name]) != solution[name]:
                raise ValueError(f"Crossing field/solution mismatch: {name}")
        for cls, prefix in (("+", "plus"), ("-", "minus")):
            modes[cls] = [{k: data[f"{prefix}_arm{i}_{k}"].copy()
                           for k in ("x", "weights", "y", "dy", "ends")} for i in (1, 2)]
    gram = [[mass_inner(modes[a], modes[b], contract) for b in ("+", "-")] for a in ("+", "-")]
    # Roundoff consistency with accepted diagnostics, not a new solver tolerance.
    tolerance = 256*np.finfo(float).eps
    np.testing.assert_allclose(gram, checks["mass_gram"], atol=tolerance, rtol=0)
    np.testing.assert_allclose(gram, np.eye(2), atol=tolerance, rtol=0)
    for cls, arms in modes.items():
        if not all(np.all(np.isfinite(a[k])) for a in arms for k in a):
            raise ValueError("Nonfinite saved fields")
        sigma = 1 if cls == "+" else -1
        for field in ("y", "dy", "ends"):
            # Verify the accepted archive relation, without defining a new operator.
            np.testing.assert_array_equal(arms[1][field], sigma*arms[0][field]*[1, -1, -1, 1, -1, -1])
        if any(np.any(a["weights"] <= 0) for a in arms):
            raise ValueError("Nonpositive quadrature weights")
    return modes, dict(recomputed_full_mass_gram=gram, source_checks=checks,
                       archive_roundoff_comparison_tolerance=tolerance,
                       exact_saved_class_relations=True)


def draw(output, full, local, solution, scenes, presentation, dpi, font_size, sorted_curves=None, branches=None, contract=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    material = {k: v for k, v in (contract or dict(epsilon=.005, aspect=1.)).items()
                if k in ("epsilon", "aspect", "E", "rho", "l", "nu", "ks", "mu", "tau")}
    geometries, captions = {}, {}
    def frequency_line(ax, name, rows, color):
        geometry = line_geometry.FrequencyGeometry([float(r["beta_deg"]) for r in rows],
            axis="BETA", fixed={}, material=material)
        geometries[name] = geometry
        line_geometry.plot_frequency_curve(ax, geometry, [float(r["Lambda"]) for r in rows],
            color=color, linewidth=1.05)

    style = {"font.family": "serif", "font.serif": ["Times New Roman"], "font.size": font_size,
             "mathtext.fontset": "stix", "pdf.fonttype": 42, "axes.grid": False,
             "figure.facecolor": "white", "axes.facecolor": "white", "path.simplify": False}
    with plt.rc_context(style):
        if sorted_curves is None and branches is None:
            fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.8))
            fig.subplots_adjust(left=.075, right=.98, bottom=.24, top=.97, wspace=.34)
            for ax, data, panel in zip(axes, (full, local), ("а", "б")):
                for cls in ("+", "-"):
                    frequency_line(ax, panel+cls, data[cls], COLORS[cls])
                ax.set_xlabel(r"$\beta,\ ^\circ$", labelpad=1)
                ax.set_ylabel(r"$\Lambda$", labelpad=2)
                ax.tick_params(direction="in", top=True, right=True, width=.6, length=2.5)
                for spine in ax.spines.values():
                    spine.set_linewidth(.6)
                formatter = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
                ax.xaxis.set_major_formatter(formatter)
                ax.yaxis.set_major_formatter(formatter)
                ax.text(.5, -.25, f"({panel})", ha="center", va="top", transform=ax.transAxes)
            axes[0].set(xlim=(0, 90), ylim=(2.25, 4.85), xticks=[0, 30, 60, 90], yticks=[2.5, 3, 3.5, 4, 4.5])
            axes[1].set(xlim=(6.1, 6.2), ylim=(3.917, 3.939), xticks=[6.1, 6.15, 6.2], yticks=[3.92, 3.93, 3.94])
            axes[1].plot(solution["beta_deg"], solution["Lambda"], "o", color="black", markersize=2.5)
            captions["class_frequencies"] = line_geometry.add_geometry_caption(fig)
            for ext in ("pdf", "png"):
                fig.savefig(output / f"class_frequencies.{ext}", dpi=dpi, bbox_inches="tight",
                            metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
            plt.close(fig)

        else:
            fig, ax = plt.subplots(figsize=(6.3, 3.7))
            fig.subplots_adjust(left=.10, right=.98, bottom=.15, top=.97)
            curves = branches if branches is not None else sorted_curves
            for rank in range(1, 9):
                rows = curves[rank]
                frequency_line(ax, str(rank), rows, BRANCH_COLORS[rank-1])
            ax.set(xlim=(0, 90), ylim=(2, 15), xticks=[0, 15, 30, 45, 60, 75, 90],
                   yticks=[2, 4, 6, 8, 10, 12, 14], xlabel=r"$\beta,\ ^\circ$", ylabel=r"$\Lambda$")
            ax.tick_params(direction="in", top=True, right=True, width=.6, length=2.5)
            for spine in ax.spines.values():
                spine.set_linewidth(.6)
            name = "branch_frequencies" if branches is not None else "sorted_frequencies"
            captions[name] = line_geometry.add_geometry_caption(fig)
            for ext in ("pdf", "png"):
                fig.savefig(output / f"{name}.{ext}", dpi=dpi, bbox_inches="tight",
                            metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
            plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(6.3, 1.45))
        fig.subplots_adjust(left=.02, right=.98, bottom=.24, top=.96, wspace=.10)
        scale = presentation["common_displacement_multiplier"]
        positions = np.vstack([x + scale*d for scene in scenes.values() for x, d in scene]
                              + [x for scene in scenes.values() for x, _ in scene])
        lower, upper = positions.min(axis=0), positions.max(axis=0)
        margin = .06*(upper[0]-lower[0])
        for ax, cls, panel in zip(axes, ("+", "-"), ("а", "б")):
            for x, d in scenes[cls]:
                ax.plot(x[:, 0], x[:, 1], color=".65", linewidth=.55)
                z = x + scale*d
                ax.plot(z[:, 0], z[:, 1], color="black", linewidth=1.5)
            ax.set(xlim=(lower[0]-margin, upper[0]+margin), ylim=(lower[1]-margin, upper[1]+margin))
            ax.set_aspect("equal")
            ax.axis("off")
            ax.text(.5, -.12, f"({panel})", ha="center", va="top", transform=ax.transAxes)
        for ext in ("pdf", "png"):
            fig.savefig(output / f"crossing_modes.{ext}", dpi=dpi,
                        metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
        plt.close(fig)

    diagnostics = line_geometry.write_geometry_diagnostics(output, geometries)
    return dict(geometry_diagnostics=diagnostics, geometry_captions=captions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--classification", type=Path, default=MAPS / "_audit_reflection_classification_20261107")
    parser.add_argument("--intersection", type=Path, default=MAPS / "_audit_reflection_intersection_20261107")
    parser.add_argument("--output-parent", type=Path, default=MAPS / "figures")
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--font-size", type=float, default=10)
    parser.add_argument("--displacement-fraction", type=float, default=.28)
    parser.add_argument("--frequency-view", choices=("classes", "sorted", "branches"), default="classes")
    parser.add_argument("--spectrum-bundle", type=Path, default=MAPS / "b775d0980d7bad77")
    parser.add_argument("--map-name", default="square_beta_mu0_tau0_refined")
    args = parser.parse_args()
    if args.dpi <= 0 or args.font_size <= 0 or args.displacement_fraction <= 0:
        raise ValueError("Presentation settings must be positive")
    cm, cv = verify_manifest(args.classification)
    im, iv = verify_manifest(args.intersection)
    class_path = (args.classification / "class_spectrum.csv").resolve().relative_to(ROOT).as_posix()
    if im["source_sha256"].get(class_path) != sha(ROOT / class_path):
        raise ValueError("Crossing archive does not certify the selected class spectrum")
    contract = im["physical_contract"]
    expected = dict(model="Rayleigh-Love full transverse inertia / Timoshenko", aspect=1., mu=0., tau=0.,
                    epsilon=.005, E=1., rho=1., l=1., nu=.3, ks=5/6, joint="RIGID", external="two clamped",
                    time_reference=200., state_order=ORDER)
    if contract != expected:
        raise ValueError("The accepted symmetric square contract is required")
    solution = read(args.intersection / "crossing_solution.json")
    if not solution["first_class_roots_certified"]:
        raise ValueError("First roots not certified in source")
    np.testing.assert_allclose(solution["omega"], solution["Lambda"]**2/contract["time_reference"], rtol=2e-15)
    full = first_classes(table(args.classification / "class_spectrum.csv"))
    if any(len(rows) != 197 or float(rows[0]["beta_deg"]) != 0 or float(rows[-1]["beta_deg"]) != 90
           for rows in full.values()):
        raise ValueError("Expected 197 archived angles covering 0--90 degrees")
    local = first_classes(table(args.intersection / "block_roots.csv"), local=True)
    behaviour = branch_behaviour(full)
    sorted_curves, sorted_checks, branches, branch_checks = None, None, None, None
    if args.frequency_view in ("sorted", "branches"):
        curves, sorted_checks = sorted_spectrum(args.spectrum_bundle, args.map_name)
        if args.frequency_view == "sorted":
            sorted_curves = insert_saved_crossing(curves, solution)
        else:
            source_rows = [r for r in table(args.spectrum_bundle / "spectrum.csv") if r["map"] == args.map_name]
            branches, branch_checks = archived_branches(table(args.classification / "class_spectrum.csv"),
                contract, sorted_checks["verified_sha256"], source_rows)
            branches = insert_branch_local_roots(branches, local, solution)
    modes, checks = load_fields(args.intersection, contract, solution)
    scenes, presentation = presentation_scenes(modes, solution["beta_deg"], contract["l"], args.displacement_fraction)
    # Curves retain the independent local roots verbatim; the saved simultaneous
    # solution is a separate black point, never a replacement/interpolated value.
    from matplotlib.font_manager import FontProperties, findfont
    import matplotlib
    font = Path(findfont(FontProperties(family="Times New Roman"), fallback_to_default=False))
    sources = [args.classification / n for n in ("manifest.json", "class_spectrum.csv")]
    sources += [args.intersection / n for n in ("manifest.json", "block_roots.csv", "crossing_solution.json",
                                               "crossing_forms.npz", "crossing_form_checks.json")]
    if sorted_checks is not None:
        sources += [args.spectrum_bundle / n for n in ("manifest.json", "spectrum.csv")]
    if branches is not None:
        sources += [ROOT / "scripts/lib/analytic_branch_tracking.py",
                    ROOT / "scripts/lib/analytic_coupled_rods_shapes.py"]
    identity = dict(script_sha256=sha(__file__), inputs={p.relative_to(ROOT).as_posix(): sha(p) for p in sources},
                    geometry_sources=line_geometry.provenance(),
                    font_sha256=sha(font), numpy=np.__version__, matplotlib=matplotlib.__version__,
                    dpi=args.dpi, font_size=args.font_size, displacement_fraction=args.displacement_fraction,
                    frequency_view=args.frequency_view, map_name=args.map_name if sorted_checks else None)
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:12]
    output = args.output_parent / ("section1_3_1_" + fingerprint)
    output.mkdir(parents=True, exist_ok=False)
    geometry_report = draw(output, full, local, solution, scenes, presentation, args.dpi, args.font_size, sorted_curves, branches, contract)
    if branches is not None:
        with (output / "plotted_frequencies.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["branch_id", "base_sorted_index", "symmetry_class",
                "class_rank", "beta_deg", "Lambda", "original_sorted_position", "mass_overlap",
                "assignment_margin", "shape_key", "source_npz", "plot_source"], extrasaction="ignore")
            writer.writeheader()
            for rows in branches.values():
                writer.writerows(rows)
        branch_checks["colors"] = {rows[0]["branch_id"]: BRANCH_COLORS[k-1] for k, rows in branches.items()}
        first_values = {c: {float(r["beta_deg"]): float(r["Lambda"]) for r in rows}
                        for c, rows in full.items()}
        branch_checks["crossing"] = dict(beta_deg=solution["beta_deg"], Lambda=solution["Lambda"],
            first_pair_branches=[branches[i][0]["branch_id"] for i in (1, 2)],
            gap_before=first_values["+"][6.1]-first_values["-"][6.1],
            gap_after=first_values["+"][6.2]-first_values["-"][6.2])
        if not branch_checks["crossing"]["gap_before"] < 0 < branch_checks["crossing"]["gap_after"]:
            raise ValueError("First branch pair does not cross")
        branch_checks["plotted_values"] = sum(map(len, branches.values()))
    elif sorted_curves is None:
        with (output / "plotted_frequencies.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["panel", "source", "symmetry_class", "beta_deg", "Lambda"])
            writer.writeheader()
            for panel, source, curves in (("a", "class_spectrum.csv", full), ("b", "block_roots.csv", local)):
                for cls, rows in curves.items():
                    for r in rows:
                        writer.writerow(dict(panel=panel, source=source, symmetry_class=cls,
                                             beta_deg=r["beta_deg"], Lambda=r["Lambda"]))
            writer.writerow(dict(panel="b", source="crossing_solution.json", symmetry_class="both",
                                 beta_deg=solution["beta_deg"], Lambda=solution["Lambda"]))
    else:
        with (output / "plotted_frequencies.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["source", "sorted_position", "beta_deg", "Lambda"])
            writer.writeheader()
            for rank, rows in sorted_curves.items():
                for row in rows:
                    writer.writerow(dict(source=row["plot_source"], sorted_position=rank,
                                         beta_deg=row["beta_deg"], Lambda=row["Lambda"]))
    with (output / "plotted_mode_axes.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["symmetry_class", "arm", "node", "x_local", "X0", "Y0",
                         "dX_signed", "dY_signed", "X_display", "Y_display"])
        scale = presentation["common_displacement_multiplier"]
        for cls, scene in scenes.items():
            for i, (x, d) in enumerate(scene):
                local_x = np.r_[0., modes[cls][i]["x"], contract["l"]]
                for j, (xj, dj) in enumerate(zip(x, d)):
                    node = "clamp" if j == 0 else "joint" if j == len(x)-1 else "saved_gauss_node"
                    writer.writerow([cls, i+1, node, local_x[j], *xj, *dj, *(xj+scale*dj)])
    write_json(output / "checks.json", dict(field_checks=checks, presentation=presentation,
                full_nodes_per_class={k: len(v) for k, v in full.items()},
                local_nodes_per_class={k: len(v) for k, v in local.items()},
                mode_panel_branch_behaviour=behaviour, sorted_spectrum_checks=sorted_checks,
                branch_continuation=branch_checks,
                frequency_interpolation=False, smoothing=False, spectral_calculations=False))
    style_sources = {
        "temp_codex_backup/chapter1_manuscript_sources_20261103/perm_2026.pdf": "Fig. 5, PDF page 10 / printed 67: frequency panels",
        "temp_codex_backup/chapter1_manuscript_sources_20261103/coupled_2025.pdf": "Fig. 5, PDF page 6 / printed 78: black modes, thin undeformed axes",
        "private/article/lambda_beta_mu0_eps0p0025_eta_0.png": "Local source figure accompanying private/article/main.tex; styling reference only",
    }
    write_json(output / "manifest.json", dict(status="COMPLETE", calculation_mode="saved_data_plot_only",
        **geometry_report,
        identity=identity, physical_contract=contract, classification_verification=cv, intersection_verification=iv,
        source_figure_style={p: dict(sha256=sha(ROOT/p), use=note) for p, note in style_sources.items()},
        style="Times New Roman / STIX; white; geometric solid/dashed criterion; no grid, titles or legends; panel letters below",
        artifacts={p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()},
        interpretation=("Eight branches selected at beta=mu=0 and continued by full-inertia shape assignment "
                        "within each reflection class, using analytic_branch_tracking helpers. Colors follow "
                        "base branch IDs, never sorted positions. Saved local roots and the simultaneous "
                        "crossing augment only the first branch of each class. No new spectral calculations."
                        if branches is not None else
                        "Frequency figure: sorted positions 1--8 from spectrum.csv. The saved coincident roots "
                        "are added only to positions 1 and 2; no frequencies are interpolated. "
                        "Mode panels: increasing class + and weakly varying class -, at one crossing."
                        if sorted_curves else "Two reflection classes; class rank is not a sorted frequency position. "
                        "Line segments join saved nodes only; the black point is the saved simultaneous solution.")))
    print(output.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()
