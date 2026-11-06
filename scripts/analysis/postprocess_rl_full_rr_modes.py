"""Saved-field mechanical audit; no basis evaluation, root search or Ritz.

Script Proliferation Control: a stable read-only postprocessor has a different
input/output contract from the immutable K8 verifier. It reuses its frozen
energy integrator, coordinate frames and full-mass inner product; solver source
and fingerprint stay unchanged. Diagnostic-only, not an article figure runner.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack, contextmanager
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.lib import isotropic_two_springs_rl_timoshenko as p
from scripts.lib import isotropic_two_springs_thickness_asymmetry as section_physics

SOURCE = ROOT/"results/rl_full_rr_prefix/ec9e338a2985a308"
OUTPUT = ROOT/"results/rl_full_rr_modal_audit"
ANGLES = [0., 5., 6.1194, 7., 30., 90.]
FLOOR = 1e-12  # relative energy/mass reporting floor, not an acceptance tolerance
AMPLITUDE_FLOOR = 1e-8  # node/scale diagnostic floor relative to peak displacement
PAIRS = [(0., 5.), (5., 6.1194), (6.1194, 7.), (5., 7.), (7., 30.), (30., 90.), (0., 30.), (0., 90.)]
FREQUENCY_POLICY = dict(frequency_map_policy="frequency-map-v1", calculation_mode="plot_only",
    spectrum_semantics="sorted_positions", sweep_parameter="beta_deg", parameter_grid=ANGLES,
    K_plot=8, K_guard=9, guard_root_role="completeness_and_quality_only",
    neighbour_audit="saved_fields_only; user_requested_global_mass_MAC_not_branch_tracking",
    local_repair_policy="none; no_new_spectral_points", strict_audit_default=False,
    source_calculation_mode="certified_audit_A_B")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def write_json(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def write_csv(path, rows):
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@contextmanager
def no_solves():
    def forbidden(*args, **kwargs):
        raise AssertionError("Saved-field postprocessing cannot evaluate bases or solve spectra")
    with ExitStack() as stack:
        for module, names in [(p, ("prefix", "count_at", "ritz", "mode", "basis", "boundary", "arm_dtn")),
                              (p.tim, ("finite_roots", "finite_mode", "finite_state_basis")),
                              (p.axial, ("basis", "explicit_frequencies"))]:
            for name in names:
                stack.enter_context(patch.object(module, name, forbidden))
        yield


def load_source(source):
    manifest = json.loads((source/"manifest.json").read_text())
    if manifest["status"] != "PASS" or manifest["stages"] != {"A": "PASS", "B": "PASS"}:
        raise ValueError("Accepted A-B archive required")
    for name, expected in manifest["artifacts"].items():
        if sha(source/name) != expected:
            raise ValueError("Changed source artifact: "+name)
    for name, expected in manifest["identity"]["sources"].items():
        if sha(ROOT/name) != expected:
            raise ValueError("Changed source implementation: "+name)
    records, modes = {}, {}
    for name in manifest["point_files"]:
        record = json.loads((source/name).read_text())
        if record["stage"] != "B":
            continue
        beta = record["beta_deg"]
        if record["status"] != "PASS" or len(record["roots"]) != 9:
            raise ValueError("Nine accepted ranks required")
        records[beta] = dict(record, source_point=name)
        with np.load(source/record["shape_file"], allow_pickle=False) as arrays:
            modes[beta] = [dict(arms=[{key: arrays[f"mode{k}_arm{arm}_{key}"]
                    for key in ("x", "weights", "y", "dy", "ends")} for arm in (1, 2)]) for k in range(1, 10)]
    if sorted(records) != ANGLES:
        raise ValueError("Only the six accepted angles are supported")
    with (source/"spectrum.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["stage"] == "B":
                root = records[float(row["beta_deg"])]["roots"][int(row["rank"])-1]
                if any(float(row[key]) != root[key] for key in ("omega", "Omega", "Lambda")):
                    raise ValueError("CSV/JSON frequency disagreement")
    model = p.Model(**records[0.]["contract"]["model"])
    for record in records.values():
        if (record["contract"]["model"] != records[0.]["contract"]["model"] or
            record["contract"]["tau"] != 0 or record["contract"]["case"]["mu"] != .01 or
            record["contract"]["case"]["k_joint"] != "RIGID" or record["contract"]["case"]["k_end"] != "RIGID"):
            raise ValueError("Shared homogeneous RR geometry required for cross-angle mass metric")
    return manifest, model, records, modes


def global_fields(model, case, mode):
    return [dict(position=(a["x"]-L)[:, None]*f.t, displacement=a["y"][:, :2]@f.translation.T,
                 joint_displacement=a["ends"][1, :2]@f.translation.T,
                 joint_force=a["ends"][1, 3:5]@f.translation.T)
            for a, L, f in zip(mode["arms"], case.lengths(model), p.geometry.frames(case.beta))]


def node_brackets(x, value, peak):
    """Sign-changing component-zero brackets on stored samples, no root finder.

    Fixed-end zeros excluded. Touching zeros not certified by this diagnostic.
    Below-floor fields are not assigned a node count.
    """
    threshold = AMPLITUDE_FLOOR*peak
    if max(abs(value)) <= threshold:
        return dict(status="BELOW_AMPLITUDE_FLOOR", count=None, brackets=[])
    nonzero = np.flatnonzero(abs(value) > threshold)
    brackets = [[float(x[i]), float(x[j])] for i, j in zip(nonzero[:-1], nonzero[1:]) if value[i]*value[j] < 0]
    return dict(status="RESOLVED_SIGN_BRACKETS", count=len(brackets), brackets=brackets)


def analyze(model, records, modes):
    adapter = SimpleNamespace(arms=(model, model), section=model.section, l=model.l)
    rows, scales, nodes, checks, parity = [], [], [], [], []
    for beta in ANGLES:
        record = records[beta]; case = p.Case(**record["contract"]["case"])
        for k, (root, mode) in enumerate(zip(record["roots"], modes[beta]), 1):
            mass, stiffness, springs = section_physics.energies(adapter, case, mode["arms"])
            M, K = float(mass.sum()), float(stiffness.sum())
            mf, kf = mass.sum(axis=0)/M, stiffness.sum(axis=0)/K
            energy_error = abs(K/M/root["omega"]**2-1)
            inner_error = abs(M-p.inner(model, mode, mode))
            check = dict(beta_deg=beta, rank=k, role="research" if k <= 8 else "guard", mass=M,
                stiffness=K, energy_error=energy_error, inner_error=inner_error,
                fractions_error=max(abs(mf.sum()-1), abs(kf.sum()-1), abs(mass.sum(axis=1).sum()/M-1), abs(stiffness.sum(axis=1).sum()/K-1)))
            policy = record["contract"]["policy"]
            if (energy_error > policy["energy_tol"] or abs(M-1) > policy["mass_gram_tol"] or
                check["fractions_error"] > 1e-12 or inner_error > 1e-12 or any(springs)):
                raise ArithmeticError("Saved-field energy/mass closure failed: "+str(check))
            checks.append(check)
            if k == 9:
                continue
            g = global_fields(model, case, mode)
            peak = max(np.linalg.norm(a["displacement"], axis=1).max() for a in g)
            joint = sum(a["joint_displacement"] for a in g)/2
            rms = np.sqrt(sum(a["weights"]@(a["y"][:, 0]**2+a["y"][:, 1]**2) for a in mode["arms"])/sum(case.lengths(model)))
            end_u = [a["ends"][1, 0] for a in mode["arms"]]
            chord = sum(model.section.EA*u*u/L for u, L in zip(end_u, case.lengths(model)))
            forces = np.array([a["ends"][1, 3:] for a in mode["arms"]])
            row = dict(beta_deg=beta, rank=k, omega=root["omega"], Lambda=root["Lambda"],
                F_omega=root["omega"]/records[0.]["roots"][k-1]["omega"],
                K=K, M=M, K_axial=float(kf[0]), K_bending=float(kf[1]), K_shear=float(kf[2]),
                M_axial=float(mf[0]), M_transverse=float(mf[1]), M_rotation=float(mf[2]), M_RL=float(mf[3]),
                K_arm1=float(stiffness[0].sum()/K), K_arm2=float(stiffness[1].sum()/K),
                M_arm1=float(mass[0].sum()/M), M_arm2=float(mass[1].sum()/M),
                joint_dx_over_peak=float(abs(joint[0])/peak), joint_dy_over_peak=float(abs(joint[1])/peak),
                joint_norm_over_rms=float(np.linalg.norm(joint)/rms),
                joint_displacement_mismatch=float(np.linalg.norm(g[0]["joint_displacement"]-g[1]["joint_displacement"])/peak),
                axial_chord_K=chord, axial_chord_share=float(chord/stiffness[:, 0].sum()) if kf[0] >= FLOOR else None,
                axial_chord_status="RESOLVED" if kf[0] >= FLOOR else "BELOW_ENERGY_FLOOR",
                N1=float(forces[0, 0]), N2=float(forces[1, 0]), Q1=float(forces[0, 1]), Q2=float(forces[1, 1]),
                Moment1=float(forces[0, 2]), Moment2=float(forces[1, 2]),
                joint_force_balance_norm=float(np.linalg.norm(g[0]["joint_force"]+g[1]["joint_force"])),
                energy_error=energy_error, source_point=record["source_point"], source_npz=record["shape_file"])
            for arm in range(2):
                for j, name in enumerate(("axial", "bending", "shear")):
                    row[f"K_{name}_arm{arm+1}"] = float(stiffness[arm, j]/K)
                for j, name in enumerate(("axial", "transverse", "rotation", "RL")):
                    row[f"M_{name}_arm{arm+1}"] = float(mass[arm, j]/M)
            rows.append(row)
            for arm, (a, L) in enumerate(zip(mode["arms"], case.lengths(model)), 1):
                for j, component in enumerate(("u", "w")):
                    data = np.r_[a["ends"][0, j], a["y"][:, j], a["ends"][1, j]]
                    x = np.r_[0., a["x"], L]
                    entry = node_brackets(x, data, peak)
                    entry.update(beta_deg=beta, rank=k, arm=arm, component=component, length=L,
                        fractions=[[lo/L, hi/L] for lo, hi in entry["brackets"]],
                        straight_s_over_2l_brackets=([[(lo if arm == 1 else 2*model.l-hi)/(2*model.l),
                                                    (hi if arm == 1 else 2*model.l-lo)/(2*model.l)]
                                                   for lo, hi in entry["brackets"]] if beta == 0 else None))
                    nodes.append(entry)
                    norm = float(a["weights"]@a["y"][:, j]**2)
                    gradient = float(a["weights"]@a["dy"][:, j]**2)
                    meaningful = max(abs(data)) > AMPLITUDE_FLOOR*peak and norm > 0
                    keff = np.sqrt(gradient/norm) if meaningful else None
                    scales.append(dict(beta_deg=beta, rank=k, arm=arm, component=component,
                        status="DEFINED_GRADIENT_INDICATOR" if meaningful else "BELOW_AMPLITUDE_FLOOR",
                        norm_integral=norm, gradient_integral=gradient, k_eff=keff,
                        k_eff_h=keff*model.section.thickness if meaningful else None,
                        k_eff_b=keff*model.section.width if meaningful else None,
                        pi_over_k=math_pi/keff if meaningful and keff > 0 else None,
                        sign_nodes=entry["count"], sign_lobes=entry["count"]+1 if entry["count"] is not None else None))
            if beta == 0:
                moments = np.zeros(8)
                for arm, (a, glob) in enumerate(zip(mode["arms"], g)):
                    coordinate = a["x"] if arm == 0 else 2*model.l-a["x"]
                    z = (coordinate-model.l)/model.l
                    for degree in range(8):
                        moments[degree] += a["weights"]@(glob["displacement"][:, 1]*z**degree)/(2*model.l*peak)
                even, odd = np.linalg.norm(moments[::2]), np.linalg.norm(moments[1::2])
                parity.append(dict(rank=k, normalized_moments=moments.tolist(),
                    indicated_parity="symmetric" if odd < even else "antisymmetric",
                    forbidden_moment_ratio=float(min(even, odd)/max(even, odd))))
    return rows, scales, nodes, checks, parity


math_pi = float(np.pi)


def comparisons(model, modes):
    matrices, best, packets = [], [], []
    for beta1, beta2 in PAIRS:
        first, second = modes[beta1], modes[beta2]
        n1 = np.array([p.inner(model, x, x) for x in first]); n2 = np.array([p.inner(model, x, x) for x in second])
        overlap = np.array([[p.inner(model, a, b, beta1, beta2) for b in second] for a in first])/np.sqrt(n1[:, None]*n2[None, :])
        mac = overlap**2
        if mac.max() > 1+1e-10:
            raise ArithmeticError("Mass MAC exceeds one")
        for i in range(9):
            for j in range(9):
                matrices.append(dict(beta_from=beta1, beta_to=beta2, rank_from=i+1, rank_to=j+1,
                                     overlap=float(overlap[i, j]), MAC=float(mac[i, j]),
                                     includes_guard=i == 8 or j == 8))
            if i < 8:
                ordered = np.argsort(mac[i])[::-1]
                best.append(dict(beta_from=beta1, beta_to=beta2, rank_from=i+1,
                    best_rank=int(ordered[0])+1, best_MAC=float(mac[i, ordered[0]]), second_rank=int(ordered[1])+1,
                    second_MAC=float(mac[i, ordered[1]]), diagonal_MAC=float(mac[i, i]),
                    captured_mass_projection_1to9=float(mac[i].sum()), status="DISCRETE_MATCH_ONLY_NOT_BRANCH_ID"))
        for pair in ((1, 2), (3, 4), (5, 6), (7, 8)):
            ids = np.array(pair)-1
            singular = np.linalg.svd(overlap[np.ix_(ids, ids)], compute_uv=False)
            packets.append(dict(beta_from=beta1, beta_to=beta2, ranks=list(pair),
                singular_values=singular.tolist(), min_projection=float(min(singular)**2),
                mean_projection=float(np.mean(singular**2))))
    return matrices, best, packets


def frequency_intervals(rows):
    by = {(r["beta_deg"], r["rank"]): r for r in rows}
    intervals = [(0., 5.), (5., 6.1194), (6.1194, 7.), (5., 7.), (7., 30.), (30., 90.)]
    return [dict(rank=k, beta_from=a, beta_to=b,
        omega_ratio=by[b, k]["omega"]/by[a, k]["omega"],
        increment_percent=100*(by[b, k]["omega"]/by[a, k]["omega"]-1),
        discrete_percent_per_degree=100*(by[b, k]["omega"]/by[a, k]["omega"]-1)/(b-a))
        for k in range(1, 9) for a, b in intervals]


def figures(folder, model, records, modes, rows, matrices):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = plt.cm.tab10(np.arange(8))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for k in range(1, 9):
        series = [r for r in rows if r["rank"] == k]
        for ax, subset in zip(axes, (series, series[:4])):
            ax.scatter([r["beta_deg"] for r in subset], [100*(r["F_omega"]-1) for r in subset],
                       color=colors[k-1], marker="o" if k % 2 else "s", label=str(k), s=32)
    for ax in axes:
        ax.set(xlabel="beta (degrees)", ylabel="100 [omega_k(beta)/omega_k(0) - 1] (%)")
        ax.grid(alpha=.2)
    axes[0].set_title("Sorted positions 1-8: saved angles only")
    axes[1].set_title("Small angles; no interpolation or threshold fit")
    axes[0].legend(title="Rank", ncol=4, fontsize=8)
    fig.savefig(folder/"frequency_changes.png", dpi=170); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), constrained_layout=True)
    for ax, key, title in zip(axes, ("K_axial", "K_arm1", "M_arm1"),
                             ("Axial stiffness share (%)", "Stiffness in arm 1 (%)", "Mass form in arm 1 (%)")):
        values = np.array([[100*next(r[key] for r in rows if r["rank"] == k and r["beta_deg"] == b) for b in ANGLES] for k in range(1, 9)])
        image = ax.imshow(values, aspect="auto", vmin=0, vmax=100, cmap="viridis")
        for i in range(8):
            for j in range(6):
                ax.text(j, i, f"{values[i,j]:.1f}", ha="center", va="center", color="white" if values[i,j] < 55 else "black", fontsize=7)
        ax.set(xticks=range(6), xticklabels=[f"{b:g}" for b in ANGLES], yticks=range(8), yticklabels=range(1, 9),
               xlabel="beta (degrees; discrete columns)", ylabel="Sorted rank", title=title)
        ax.tick_params(axis="x", rotation=45)
    fig.colorbar(image, ax=axes, shrink=.8)
    fig.savefig(folder/"energy_and_localization.png", dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    for ax, (a, b) in zip(axes.flat, [(5., 6.1194), (6.1194, 7.), (7., 30.), (30., 90.)]):
        block = np.array([[next(r["MAC"] for r in matrices if r["beta_from"] == a and r["beta_to"] == b and r["rank_from"] == i and r["rank_to"] == j) for j in range(1, 9)] for i in range(1, 9)])
        im = ax.imshow(block, vmin=0, vmax=1, cmap="magma")
        for i, j in zip(*np.where(block >= .08)):
            ax.text(j, i, f"{block[i,j]:.2f}", ha="center", va="center", color="black" if block[i,j] > .65 else "white", fontsize=8)
        ax.set(xticks=range(8), yticks=range(8), xticklabels=range(1, 9), yticklabels=range(1, 9),
               xlabel=f"Rank at {b:g} deg", ylabel=f"Rank at {a:g} deg", title=f"Global full-mass MAC: {a:g} -> {b:g}")
    fig.colorbar(im, ax=axes, shrink=.8)
    fig.savefig(folder/"global_mass_MAC.png", dpi=170); plt.close(fig)

    selected = [(0., 1), (0., 4), (0., 6), (0., 8), (6.1194, 1), (6.1194, 2), (30., 6), (90., 8)]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), constrained_layout=True)
    for ax, (beta, k) in zip(axes.flat, selected):
        case = p.Case(**records[beta]["contract"]["case"]); mode = modes[beta][k-1]
        fields = global_fields(model, case, mode)
        peak = max(np.linalg.norm(a["displacement"], axis=1).max() for a in fields)
        all_y = np.concatenate([a["displacement"][:, 1] for a in fields])
        sign = 1 if all_y[np.argmax(abs(all_y))] >= 0 else -1
        scale = .16/peak
        for arm, L, frame in zip(mode["arms"], case.lengths(model), p.geometry.frames(beta)):
            coord = np.r_[0., arm["x"], L]
            values = np.vstack((arm["ends"][0], arm["y"], arm["ends"][1]))
            pos = (coord-L)[:, None]*frame.t
            disp = values[:, :2]@frame.translation.T
            ax.plot(pos[:, 0], pos[:, 1], "k--", lw=.7)
            moved = pos+sign*scale*disp
            ax.plot(moved[:, 0], moved[:, 1], lw=1.4)
        joint = fields[0]["joint_displacement"]*sign*scale
        ax.plot(joint[0], joint[1], "ko", ms=3)
        ax.set(title=f"beta={beta:g} deg, rank {k}", xlabel="X/l", ylabel="Y/l")
        ax.axis("equal"); ax.grid(alpha=.15)
    fig.suptitle("Global displacement shapes; each plot has its own amplitude and arbitrary sign")
    fig.savefig(folder/"selected_global_shapes.png", dpi=170); plt.close(fig)


def run(source=SOURCE, output=OUTPUT):
    if output.resolve().is_relative_to(source.resolve()):
        raise ValueError("Postprocessing output must be outside the source bundle")
    with no_solves():
        import matplotlib
        manifest, model, records, modes = load_source(source)
        identity = dict(version="saved-rl-full-rr-modal-audit-v1", source=source.relative_to(ROOT).as_posix(),
            source_manifest_sha256=sha(source/"manifest.json"), source_artifacts=manifest["artifacts"],
            code_sha256=sha(Path(__file__)), energy_floor=FLOOR, amplitude_floor=AMPLITUDE_FLOOR,
            frequency_policy=FREQUENCY_POLICY,
            source_implementation_hashes=manifest["identity"]["sources"],
            versions=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__),
            method="stored NPZ quadrature only; global translation + material psi/u-prime full-mass metric; no root/basis calls")
        key = digest(identity)[:16]; folder = output/key
        if (folder/"manifest.json").exists():
            existing = json.loads((folder/"manifest.json").read_text())
            if not all(sha(folder/n) == h for n, h in existing["artifacts"].items()):
                raise ValueError("Changed postprocessing artifacts")
            return folder, existing
        folder.mkdir(parents=True, exist_ok=False)
        rows, scales, nodes, checks, parity = analyze(model, records, modes)
        matrices, best, packets = comparisons(model, modes)
        write_json(folder/"identity.json", identity)
        write_csv(folder/"modal_components.csv", rows)
        write_csv(folder/"frequency_intervals.csv", frequency_intervals(rows))
        write_csv(folder/"spatial_scales.csv", scales)
        write_csv(folder/"global_mass_MAC.csv", matrices)
        write_csv(folder/"discrete_matches.csv", best)
        write_json(folder/"node_brackets.json", nodes)
        write_json(folder/"straight_parity_moments.json", parity)
        write_json(folder/"pair_subspaces.json", packets)
        write_json(folder/"quality_checks.json", dict(status="PASS", modes=checks,
            maximum_energy_error=max(c["energy_error"] for c in checks),
            maximum_fraction_error=max(c["fractions_error"] for c in checks),
            spectral_solves=0, basis_evaluations=0, angles=ANGLES,
            reporting_floor=FLOOR, numerical_tolerances=records[0.]["contract"]["policy"]))
        figures(folder, model, records, modes, rows, matrices)
        # Original manifest and every source file must still match after plotting.
        load_source(source)
        result = dict(status="PASS", identity=identity, primary_modes=48, guard_quality_checks=6,
            spectral_solves=0, basis_evaluations=0,
            artifacts={x.name: sha(x) for x in sorted(folder.iterdir()) if x.is_file()})
        write_json(folder/"manifest.json", result)
        return folder, result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    folder, result = run(args.source.resolve(), args.output.resolve())
    print(json.dumps(dict(folder=str(folder), status=result["status"], primary_modes=result["primary_modes"],
                          spectral_solves=0, basis_evaluations=0), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
