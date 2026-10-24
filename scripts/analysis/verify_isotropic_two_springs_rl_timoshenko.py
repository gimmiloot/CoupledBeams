"""Bounded two-spring verification. One diagnostic CLI; no production changes.

Script Proliferation Control: external rotational compliance changes the trace
space and input/output contract. Reuses LONG/Tim arm bases, geometry and pole
catalogue; separate from immutable G20 and symmetric-laminate workflows.
"""
from pathlib import Path
import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
import math
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.linalg import block_diag, null_space
from scipy.optimize import brentq
from scripts.lib import isotropic_two_springs_rl_timoshenko as physics
from scripts.lib import coupled_longitudinal_comparators as legacy
from scripts.lib import contraction_joint_sensitivity as spring
from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as rectangular
from scripts.analysis.screen_coupled_longitudinal_theory_hierarchy import tim_catalog

CONFIG = ROOT/"data/input/isotropic_two_springs_rl_timoshenko.json"
OUTPUT = ROOT/"results/isotropic_two_springs_rl_timoshenko"
LIMIT_CONFIG = ROOT/"data/input/isotropic_two_springs_limit_angles.json"
FINITE_CONFIG = ROOT/"data/input/isotropic_two_springs_finite_angles.json"
SATURATION_CONFIG = ROOT/"data/input/isotropic_two_springs_saturation.json"
PLATEAU_CONFIG = ROOT/"data/input/isotropic_two_springs_plateau_criterion.json"
GRID_CONFIG = ROOT/"data/input/isotropic_two_springs_plateau_grid.json"
EARLY_CONFIG = ROOT/"data/input/isotropic_two_springs_early_saturation.json"
GROWTH_CONFIG = ROOT/"data/input/isotropic_two_springs_growth_transition.json"
REFINE_CONFIG = ROOT/"data/input/isotropic_two_springs_growth_refinement.json"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, value):
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(p.suffix+".tmp")
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else x.tolist()), encoding="utf8")
    tmp.replace(p)


def csvwrite(p, rows):
    if not rows:
        return
    with p.open("w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)


def identity(config):
    names = [CONFIG, Path(__file__), Path(physics.__file__), Path(legacy.__file__),
             Path(spring.__file__), Path(rectangular.__file__), Path(physics.tim.__file__),
             Path(physics.axial.__file__), Path(physics.geometry.__file__),
             ROOT/"scripts/lib/reddy_inplane_geometry.py",
             ROOT/"scripts/lib/reddy_symmetric_coupled_beams_ritz.py",
             ROOT/"scripts/analysis/verify_mindlin_herrmann_timoshenko_general_beta_joint.py",
             ROOT/"scripts/analysis/screen_coupled_longitudinal_theory_hierarchy.py"]
    ref = ROOT/config["reference"]
    frozen = json.loads((ref/"SHA256.json").read_text())
    if any(sha(ref/f) != h for f, h in frozen.items()):
        raise ValueError("Independent reference bundle changed")
    record = {"config": config, "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in names},
              "reference_hashes": frozen, "versions": {"python": platform.python_version(),
              "numpy": np.__version__, "scipy": scipy.__version__}}
    key = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16]
    return key, record


class Run:
    def __init__(self, config, folder):
        self.config, self.folder, self.policy = config, folder, config["policy"]
        self.model = physics.Model(**config["material"])
        self.catalogs = {}
        self.new_points = self.reused_points = self.new_roots = self.old_roots = self.ritz_solves = 0
        self.pole_roots = 0

    def cats(self, model, case):
        result = {}
        for L in case.lengths(model):
            if L not in self.catalogs:
                p = self.policy
                base = {"policy": {"scan_intervals": p["pole_scan_intervals"],
                    "root_xtol": p["root_xtol"], "root_rtol": p["root_rtol"]}}
                self.catalogs[L] = tim_catalog(model.carrier, L,
                                    model.omega(p["Lambda_upper"]), base)
                self.pole_roots += len(self.catalogs[L]["roots"])
            result[L] = self.catalogs[L]
        return result

    def solve(self, case, variant="RL_FULL", neighbour=None):
        model = replace(self.model, variant=variant)
        contract = {**asdict(case), "variant": variant}
        key = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:16]
        path = self.folder/"points"/(key+".json")
        if path.exists():
            rec = json.loads(path.read_text(encoding="utf8"))
            if rec["case"] != contract or sha(path.with_suffix('.npz')) != rec["shape_sha256"]:
                raise ValueError("Case cache/provenance mismatch")
            self.reused_points += 1
            return rec
        roots, search = physics.prefix(model, case, self.cats(model, case), self.policy, neighbour)
        modes = [physics.mode(model, case, r["omega"], self.policy["quadrature_order"]) for r in roots]
        gram = np.array([[physics.inner(model, a, b) for b in modes] for a in modes])
        err = float(np.max(abs(gram-np.eye(3))))
        for i, (r, mode) in enumerate(zip(roots, modes), 1):
            d = mode["diagnostics"]
            gates = {"boundary": max(np.max(d["joint_scaled"]), np.max(d["external_scaled"])) <= self.policy["boundary_tol"],
                "ode": d["ode_error"] <= self.policy["ode_tol"],
                "energy": d["energy_error"] <= self.policy["energy_tol"],
                "singularity": d["singular_ratio"] <= self.policy["singular_tol"],
                "simple": d["nonzero_condition"] <= self.policy["nonzero_condition_max"],
                "mass_gram": err <= self.policy["mass_gram_tol"]}
            r.update(sorted_position=i, Omega=r["omega"]*model.time,
                Lambda=float(model.Lambda(r["omega"])), diagnostics=d,
                coefficients=mode["coefficients"], gates=gates)
        rec = {"case": contract, "roots": roots, "count": search, "mass_gram_error": err,
               "status": "PASS" if all(all(r["gates"].values()) for r in roots) else "FAIL"}
        if rec["status"] != "PASS":
            write(self.folder/"failures"/(key+'.json'), rec)
            raise ArithmeticError(f"Physical gate failed {contract}")
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path.with_suffix('.npz'),
            **{f"mode{i+1}_arm{j+1}_{field}": a[field]
               for i, mode in enumerate(modes) for j, a in enumerate(mode["arms"])
               for field in ("x", "weights", "y", "dy", "ends")})
        rec["shape_sha256"] = sha(path.with_suffix('.npz'))
        rec["shape_file"] = path.with_suffix('.npz').relative_to(self.folder).as_posix()
        write(path, rec)
        self.new_points += 1; self.new_roots += 3
        return rec

    def shapes(self, record):
        data = np.load(self.folder/record["shape_file"])
        return [{"arms": [{f: data[f"mode{i+1}_arm{j+1}_{f}"]
                for f in ("x", "weights", "y", "dy", "ends")} for j in range(2)]}
                for i in range(3)]

    def ritz_check(self, rec):
        case = physics.Case(**{k: rec["case"][k] for k in ("beta", "mu", "k_joint", "k_end")})
        model = replace(self.model, variant=rec["case"]["variant"])
        values = []
        for n in self.policy["ritz_orders"]:
            w, coef, diag = physics.ritz(model, case, n)
            self.ritz_solves += 1
            rel = float(max(abs(w/np.array([r["omega"] for r in rec["roots"]])-1)))
            values.append({"order": n, "omega": w.tolist(), "relative_frequency_error": rel, **diag})
            if rel > self.policy["ritz_frequency_tol"]:
                raise ArithmeticError("Ritz comparison gate failed")
            if diag['mass_gram_error'] > self.policy['mass_gram_tol'] or diag['energy_error'] > self.policy['energy_tol']:
                raise ArithmeticError("Ritz physical mass/energy gate failed")
        return {"case": rec["case"], "orders": values}


def structural(run):
    m = run.model; rng = np.random.default_rng(20261022); checks = []
    for beta in (0, 30, 90):
        for kj in (0, 1, "RIGID"):
            for ke in (0, 1, "RIGID"):
                c = physics.Case(beta, k_joint=kj, k_end=ke)
                op = np.zeros((12, 24))
                op[:3, :6] = physics.external_matrix(m, c, 0)
                op[3:6, 12:18] = physics.external_matrix(m, c, 1)
                op[6:, 6:12] = physics.joint_matrix(m, c)[:, :6]
                op[6:, 18:] = physics.joint_matrix(m, c)[:, 6:]
                scaled = op/np.linalg.norm(op, axis=1)[:, None]
                kernel = null_space(scaled)
                states = (kernel@rng.normal(size=12)).reshape(4, 6)
                p, sp = physics.trace_map(m, c)
                z = rng.normal(size=p.shape[1]); delta = (p@z).reshape(4, 3)
                q = states[:, :3].flatten()
                coord = np.linalg.lstsq(p, q, rcond=None)[0]
                work = sum(sign*state[3:]@v for sign, state, v in zip((-1, 1, -1, 1), states, delta))+coord@sp@z
                rank = np.linalg.matrix_rank(scaled)
                dual = max(np.max(abs(f.translation.T@f.translation-np.eye(2))) for f in physics.geometry.frames(beta))
                if rank != 12 or abs(work) > 1e-12 or dual > 1e-14:
                    raise ArithmeticError("Structural work/rank gate failed")
                checks.append({"beta": beta, "k_joint": kj, "k_end": ke,
                    "rank": int(rank), "trace_dimension": p.shape[1], "work_error": float(abs(work)), "duality": float(dual)})
    # Entire local basis satisfies independently assembled first-order balances.
    ode = []
    for variant in physics.VARIANTS:
        mod = replace(m, variant=variant)
        w = float(mod.omega(4.))
        b = physics.basis(mod, .99, w, [.1, .5, .9])
        d = physics.basis(mod, .99, w, [.1, .5, .9], 1)
        error = np.max(abs(d-np.einsum('ij,njk->nik', physics.state_matrix(mod, w), b)))/np.max(abs(d))
        ode.append({"variant": variant, "relative_ode": float(error)})
        if error > run.policy["ode_tol"]:
            raise ArithmeticError("Local ODE gate failed")
    return {"status": "PASS", "rank_work": checks, "basis_ode": ode,
            "positivity": "Positive rod energy and mass; clamp1 excludes all zero-energy rigid motions, including the double hinge"}


def reference_stage(run):
    cfg = run.config; ref = ROOT/cfg["reference"]
    table = list(csv.DictReader((ref/'results/spectrum.csv').open()))
    summary = json.loads((ref/'results/summary.json').read_text())
    rows, ritz = [], []
    for beta in cfg["control_angles"]+[summary["gap_minimum"]["beta_deg"]]:
        rec = run.solve(physics.Case(beta, cfg["mu"]))
        matched = [r for r in table if abs(float(r['beta_deg'])-beta) < 1e-10]
        for n, root in enumerate(rec['roots'], 1):
            val = float(matched[0][f'Lambda{n}' if n < 3 else 'Lambda3_guard']) if matched else None
            diff = root['Lambda']/val-1 if val else None
            rows.append({"case_id": f'A_beta{beta:.12g}', "model": "RL_FULL", "beta_deg": beta,
                "mu": cfg['mu'], "k_joint": 'RIGID', "k_end": 'RIGID', "sorted_position": n,
                "omega": root['omega'], "Lambda": root['Lambda'], "reference_model": 'CHATGPT_RL_FULL',
                "reference_Lambda": val, "relative_difference": diff,
                "comparison_kind": 'CHATGPT_REFERENCE_REPRODUCTION',
                "status": 'UNAVAILABLE_EXACT_ANGLE' if val is None else 'PASS' if abs(diff) <= run.policy['reference_frequency_tol'] else 'FAIL'})
        ritz.append(run.ritz_check(rec))
    if any(r['status'] == 'FAIL' for r in rows):
        write(run.folder/'reference_failure.json', rows)
        raise ArithmeticError('ChatGPT reproduction failed')
    csvwrite(run.folder/'reference_comparison.csv', rows)
    write(run.folder/'reference_ritz.json', ritz)
    return {"status": "PASS", "points": len(ritz), "max_relative_difference": max(abs(r['relative_difference']) for r in rows if r['relative_difference'] is not None),
            "rounded_angle_qualification": "6.1194 differs from the archive minimum angle 6.1194280284140214; both evaluated, no rounded-answer fitting"}


def map_stage(run):
    rows = list(csv.DictReader((ROOT/run.config['reference']/'results/spectrum.csv').open()))
    results, records, previous = [], [], None
    for row in rows:
        beta = float(row['beta_deg'])
        rec = run.solve(physics.Case(beta, run.config['mu']), neighbour=previous)
        previous = [r['omega'] for r in rec['roots']]
        records.append(rec)
        values = rec['roots']
        results.append({"beta_deg": beta, **{f'omega{i+1}': r['omega'] for i, r in enumerate(values)},
            **{f'Lambda{i+1}': r['Lambda'] for i, r in enumerate(values)},
            'gap_Omega': values[1]['Omega']-values[0]['Omega'],
            'gap_Lambda': values[1]['Lambda']-values[0]['Lambda'],
            'reference_relative_max': max(abs(r['Lambda']/float(row[f'Lambda{i+1}' if i < 2 else 'Lambda3_guard'])-1) for i, r in enumerate(values)),
            'count_status': rec['count']['status']})
    csvwrite(run.folder/'spectrum.csv', results)
    if max(r['reference_relative_max'] for r in results) > run.policy['reference_frequency_tol']:
        raise ArithmeticError('Map reference comparison failed; no interpretation allowed')
    local = [r for r in records if 5-1e-10 <= r['case']['beta'] <= 7+1e-10]
    shapes = [run.shapes(r) for r in local]
    macs = []; margins = []; adjacent = []
    for ra, rb, a, b in zip(local[:-1], local[1:], shapes[:-1], shapes[1:]):
        ba, bb = ra['case']['beta'], rb['case']['beta']
        mac = np.array([[physics.inner(run.model, aa, bbb, ba, bb)**2 for bbb in b[:2]] for aa in a[:2]])
        macs.extend(np.diag(mac)); margins.extend([mac[0, 0]-mac[0, 1], mac[1, 1]-mac[1, 0]])
        adjacent.append({'beta0': ba, 'beta1': bb, 'physical_global_MAC': mac.tolist()})
    a, b = shapes[0], shapes[-1]
    far = [[physics.inner(run.model, aa, bb, 5., 7.)**2 for bb in b[:2]] for aa in a[:2]]
    far_local = [[physics.inner(run.model, aa, bb)**2 for bb in b[:2]] for aa in a[:2]]
    minimum = min(results, key=lambda r:r['gap_Omega'])
    flat = [r for r in results if r['beta_deg'] >= 10-1e-10]
    first = np.array([r['Lambda1'] for r in flat]); omegas = np.array([r['omega1'] for r in flat])
    fractions = []
    for rec, ss in ((local[0], shapes[0]), (local[-1], shapes[-1])):
        for i, sh in enumerate(ss[:2], 1):
            arm = sh['arms'][1]; s = run.model.section
            mm = arm['weights']@(s.rhoA*(arm['y'][:, 0]**2+arm['y'][:, 1]**2)+s.rhoI*arm['y'][:, 2]**2+run.model.J*arm['dy'][:, 0]**2)
            fractions.append({'beta': rec['case']['beta'], 'position': i, 'long_arm_mass_fraction': float(mm)})
    jumps = []
    for a, b in zip(results[:-1], results[1:]):
        if max(abs(b[f'omega{k}']/a[f'omega{k}']-1) for k in (1, 2)) > run.policy['neighbour_relative_jump']:
            jumps.append([a['beta_deg'], b['beta_deg']])
    veering = (min(macs) >= run.policy['mac_min'] and min(margins) >= run.policy['assignment_margin_min']
               and far[0][1] > far[0][0] and far[1][0] > far[1][1] and minimum['gap_Omega'] > 0)
    result = {'status': 'PASS', 'point_count': len(records), 'sampled_minimum': minimum,
        'flatness': {'Lambda_relative_range': float(np.ptp(first)/min(first)), 'omega_relative_range': float(np.ptp(omegas)/min(omegas)),
                     'sampled_maximum_beta': flat[int(np.argmax(first))]['beta_deg']},
        'local_pair': {'status': 'NUMERICALLY_CONFIRMED_LOCAL_VEERING' if veering else 'VEERING_CANDIDATE',
             'seed': [5., 1, 2], 'adjacent_min_MAC': min(macs), 'minimum_margin': min(margins),
             'endpoint_global_MAC': far, 'endpoint_local_pullback_MAC': far_local,
             'long_arm_mass_fractions': fractions, 'adjacent': adjacent},
        'neighbour_jump_flags': jumps, 'global_reference_relative_max': max(r['reference_relative_max'] for r in results)}
    write(run.folder/'map_diagnostics.json', result)
    return result


def independent_scan(matrix, model, policy):
    roots = []; step = policy['independent_scan_Lambda_step']
    def det(lam):
        return float(np.linalg.det(physics.balance(matrix(float(model.omega(lam))))[0]))
    a = policy['Lambda_lower']; fa = det(a)
    for b in np.arange(a+step, policy['Lambda_upper']+step/2, step):
        fb = det(b)
        if fa*fb < 0:
            lam = brentq(det, a, b, xtol=2e-13, rtol=2e-14)
            roots.append(float(model.omega(lam)))
            if len(roots) == 3:
                return roots
        a, fa = b, fb
    raise ArithmeticError('Independent legacy scan did not resolve three roots')


def regression_stage(run):
    rows, matrices, ritz, legacy_shapes = [], [], [], []
    m = run.model; p = run.policy
    jobs = [('B1', 'ELEMENTARY_LEGACY', b, 'RIGID') for b in run.config['control_angles']]
    jobs += [('B2', 'RL_PLANAR_LEGACY', b, 'RIGID') for b in run.config['control_angles']]
    jobs += [('B3', 'RL_PLANAR_LEGACY', b, k) for b in run.config['spring_angles'] for k in (0, 1, 100)]
    for label, variant, beta, kj in jobs:
        mod = replace(m, variant=variant); case = physics.Case(beta, run.config['mu'], kj)
        new = run.solve(case, variant); lengths = case.lengths(mod)
        if label == 'B1':
            old = lambda w: rectangular.legacy_coupled_boundary_matrix_raw(w, mod.section, lengths[0], mod.section, lengths[1], beta_deg=beta)
        elif label == 'B2':
            old = lambda w: legacy.boundary(mod.carrier, lengths, w, 'rayleigh_love_planar', physics.geometry.frames(beta))
        else:
            sc = spring.Case(beta, variant='rayleigh_love_planar', rotation='SPRING')
            old = lambda w: spring.boundary(mod.carrier, lengths, w, sc, kj*mod.section.EI/mod.l, scaled=False)
        reference = independent_scan(old, mod, p); run.old_roots += 3
        # Independent scan has no new-root seeds. Its prefix is qualified by the
        # proved same operator and the independent DtN count of the new assembly.
        if label != 'B1':
            w = float(mod.omega(4.2)); newb = physics.boundary(mod, case, w)
            oldb = old(w)
            if label == 'B2':
                newb = newb[[0, 1, 2, 3, 4, 5, 6, 7, 10, 8, 9, 11]]
            er = float(np.max(abs(newb-oldb))/np.max(abs(oldb)))
            matrices.append({'case': new['case'], 'relative_matrix_error': er})
            if er > 1e-12:
                raise ArithmeticError('Legacy operator mapping failed')
        for i, (root, wr) in enumerate(zip(new['roots'], reference), 1):
            rel = root['omega']/wr-1
            rows.append({'case_id': f'{label}_b{beta}_k{kj}', 'model': variant, 'beta_deg': beta,
                'mu': case.mu, 'k_joint': kj, 'k_end': 'RIGID', 'sorted_position': i,
                'omega': root['omega'], 'Lambda': root['Lambda'], 'reference_model': label,
                'reference_Lambda': float(mod.Lambda(wr)), 'relative_difference': rel,
                'comparison_kind': 'SAME_PHYSICS_REGRESSION', 'status': 'PASS' if abs(rel) <= p['same_physics_frequency_tol'] else 'FAIL'})
        if label == 'B1' and beta in (0, 6.1194, 90):
            shapes = run.shapes(new)
            for i, w in enumerate(reference):
                coef = rectangular.coupled_mode_coefficients(w, mod.section, lengths[0], mod.section, lengths[1], beta_deg=beta).coefficients
                sh = {'arms': []}
                for arm, (L, a) in enumerate(zip(lengths, shapes[i]['arms'])):
                    sign = 1 if arm == 0 else -1
                    f = rectangular.evaluate_clamped_arm_fields(w, mod.section, rectangular.timoshenko_spatial_basis(w, mod.section), coef[3*arm:3*arm+3], sign*a['x'])
                    y = np.column_stack((sign*f.u, sign*f.w, f.psi, f.N, f.Q, sign*f.M))
                    dy = np.column_stack((f.u_prime, f.w_prime, sign*f.psi_prime, sign*f.N_prime, sign*f.Q_prime, f.M_prime))
                    sh['arms'].append({'y': y, 'dy': dy, 'weights': a['weights']})
                mac = physics.inner(mod, shapes[i], sh)**2/physics.inner(mod, sh, sh)
                legacy_shapes.append({'beta': beta, 'position': i+1, 'physical_mass_MAC': mac})
                if mac < .999999:
                    raise ArithmeticError('Legacy field coordinate conversion failed')
    if any(r['status'] != 'PASS' for r in rows):
        write(run.folder/'regression_failure.json', rows)
        raise ArithmeticError('Same-physics regression failed')
    # Signed differences in A1 and the finite-joint B3 controls, never a map.
    theory_cases = [physics.Case(b) for b in run.config['control_angles']]
    theory_cases += [physics.Case(b, k_joint=k) for b in run.config['spring_angles'] for k in (0, 1, 100)]
    for case in theory_cases:
        beta = case.beta
        full = run.solve(case)
        for variant in ('RL_PLANAR_LEGACY', 'ELEMENTARY_LEGACY'):
            rec = run.solve(case, variant)
            for i, (r, ref) in enumerate(zip(rec['roots'], full['roots']), 1):
                rows.append({'case_id': f'B4_b{beta}', 'model': variant, 'beta_deg': beta,
                    'mu': case.mu, 'k_joint': case.k_joint, 'k_end': case.k_end, 'sorted_position': i,
                    'omega': r['omega'], 'Lambda': r['Lambda'], 'reference_model': 'RL_FULL',
                    'reference_Lambda': ref['Lambda'], 'relative_difference': r['omega']/ref['omega']-1,
                    'comparison_kind': 'DIFFERENT_THEORY_COMPARISON', 'status': 'COMPUTED'})
    # Finite-end checks, with an independently assembled variational pencil.
    for kwargs in run.config['end_checks']:
        rec = run.solve(physics.Case(**kwargs)); ritz.append(run.ritz_check(rec))
    monotonic = []
    for beta in run.config['spring_angles']:
        spectra = [run.solve(physics.Case(beta, k_joint=k), 'RL_PLANAR_LEGACY') for k in run.config['spring_values']]
        values = np.array([[r['omega'] for r in s['roots']] for s in spectra])
        minimum = float(np.min(np.diff(values, axis=0)))
        monotonic.append({'beta': beta, 'minimum_increment': minimum})
        if minimum < -1e-10:
            raise ArithmeticError('Internal-spring min-max gate failed')
    direct = []
    # Single beam CC bending, independent finite-boundary determinant. Axial
    # first root is analytically above this prefix; no fictitious cut modes.
    def single_matrix(w):
        return physics.tim.finite_boundary_matrix(m.carrier, 2*m.l, w, 'timoshenko')
    single_roots = independent_scan(single_matrix, m, p); run.old_roots += 3
    for mu in (0., .01, -.01):
        rec = run.solve(physics.Case(0, mu))
        relative = max(abs(np.array([r['omega'] for r in rec['roots']])/single_roots-1))
        direct.append({'mu': mu, 'relative_frequency_error': float(relative)})
        if relative > p['same_physics_frequency_tol']:
            raise ArithmeticError('Straight split transparency failed')
    csvwrite(run.folder/'regression_comparison.csv', rows)
    result = {'status': 'PASS', 'matrix_checks': matrices, 'legacy_shape_checks': legacy_shapes,
              'end_spring_ritz': ritz, 'joint_monotonicity': monotonic,
              'straight_split': direct, 'single_rod_omega': single_roots,
              'legacy_count_qualification': 'Unmodified legacy boundary operators with independent sign scans; NOT the G20-only historical count API. Prefix equivalence follows from operator maps plus independently certified new count.'}
    write(run.folder/'regression_diagnostics.json', result)
    return result


def supplemental_checks(run):
    """Read accepted roots only: direct forms, end-limit order, source Ritz."""
    m = run.model
    checks = []
    for mu in (0., .01, -.01):
        rec = run.solve(physics.Case(0, mu))
        modes = run.shapes(rec)
        for root, mode in zip(rec['roots'], modes):
            w = root['omega']; L = 2*m.l
            b = physics.tim.finite_boundary_matrix(m.carrier, L, w, 'timoshenko')
            a, _, col = physics.balance(b)
            _, _, vh = np.linalg.svd(a); coef = col*vh[-1]
            direct = {'arms': []}
            for i, arm in enumerate(mode['arms']):
                x = arm['x'] if i == 0 else L-arm['x']
                vals = physics.tim.finite_state_basis(m.carrier, L, w, x, 'timoshenko')@coef
                grad = physics.tim.finite_state_basis(m.carrier, L, w, x, 'timoshenko', 1)@coef
                y, dy = np.zeros((len(x), 6)), np.zeros((len(x), 6))
                y[:, [1, 2, 4, 5]], dy[:, [1, 2, 4, 5]] = vals, grad
                if i == 1:
                    transform = np.array([-1, -1, 1, 1, 1, -1])
                    y *= transform; dy *= -transform
                direct['arms'].append({'y': y, 'dy': dy, 'weights': arm['weights']})
            mac = physics.inner(m, mode, direct)**2/physics.inner(m, direct, direct)
            checks.append({'mu': mu, 'position': root['sorted_position'], 'direct_split_MAC': mac})
            if mac < 1-1e-8:
                raise ArithmeticError('Direct physical shape transparency failed')
    end = [run.solve(physics.Case(30, k_joint=1, k_end=k)) for k in (0, 1, 100, 'RIGID')]
    inc = np.diff([[r['omega'] for r in rec['roots']] for rec in end], axis=0)
    if inc.min() < -1e-10:
        raise ArithmeticError('External-spring min-max gate failed')
    # The archive is imported without invoking its map/plot entry point.
    spec = importlib.util.spec_from_file_location('independent_reference_ritz',
            ROOT/run.config['reference']/'independent_ritz.py')
    archived = importlib.util.module_from_spec(spec); spec.loader.exec_module(archived)
    source_checks = []
    for beta in (0, 6.1194, 90):
        rec = run.solve(physics.Case(beta)); source = archived.solve(.01, beta, degree=20, n_modes=3)
        run.ritz_solves += 1
        error = float(max(abs(source['Lambda']/[r['Lambda'] for r in rec['roots']]-1)))
        source_checks.append({'beta': beta, 'relative_Lambda_error': error,
                              'diagnostics': source['diagnostics']})
        if error > run.policy['ritz_frequency_tol']:
            raise ArithmeticError('Archived independent Ritz failed')
    result = {'direct_shape_checks': checks, 'external_spring_minimum_increment': float(inc.min()),
              'archived_independent_Ritz': source_checks}
    write(run.folder/'supplemental_checks.json', result)
    return result


def reuse_points(run, previous):
    """Explicit development restart: unchanged physics/config, rechecked kernels.

    Not silent cross-version reuse. Only the orchestration file may differ;
    all physics/dependency bytes and the frozen configuration must match.
    Existing shapes and root brackets are retained, not re-solved.
    """
    old = json.loads((previous/'contract.json').read_text(encoding='utf8'))
    _, current = identity(run.config)
    own = Path(__file__).relative_to(ROOT).as_posix()
    if old['config'] != current['config'] or old['reference_hashes'] != current['reference_hashes'] or old['versions'] != current['versions']:
        raise ValueError('Restart contract mismatch')
    if any(v != current['files'][k] for k, v in old['files'].items() if k != own):
        raise ValueError('Restart physics/dependency changed')
    evidence = []
    for path in sorted((previous/'points').glob('*.json')):
        rec = json.loads(path.read_text()); c = rec['case']
        model = replace(run.model, variant=c['variant'])
        case = physics.Case(**{k: c[k] for k in ('beta', 'mu', 'k_joint', 'k_end')})
        if sha(previous/rec['shape_file']) != rec['shape_sha256']:
            raise ValueError('Restart shape hash mismatch')
        for r in rec['roots']:
            a = physics.balance(physics.boundary(model, case, r['omega']))[0]
            sv = np.linalg.svd(a, compute_uv=False)
            if sv[-1]/sv[0] > run.policy['singular_tol'] or not all(r['gates'].values()):
                raise ValueError('Restart root rejected by current operator')
        dest = run.folder/'points'/path.name; dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)
        shutil.copyfile(previous/rec['shape_file'], run.folder/rec['shape_file'])
        evidence.append({'file': path.name, 'sha256': sha(path), 'shape_sha256': rec['shape_sha256']})
    write(run.folder/'reused_points.json', {'source': str(previous), 'source_contract_sha256': sha(previous/'contract.json'),
                                           'validation': 'unchanged physical dependencies + current kernel checks', 'points': evidence})
    complete = previous/'manifest.json'
    if complete.exists():
        manifest = json.loads(complete.read_text(encoding='utf8'))
        if any(sha(previous/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Completed predecessor artifacts changed')
        for name in ('reference_comparison.csv', 'reference_ritz.json', 'spectrum.csv',
                     'map_diagnostics.json', 'regression_comparison.csv',
                     'regression_diagnostics.json', 'supplemental_checks.json',
                     'stages.json', 'pole_catalogs.json'):
            shutil.copyfile(previous/name, run.folder/name)
        write(run.folder/'reused_checks.json', {'source_manifest': str(complete),
            'source_manifest_sha256': sha(complete), 'source_counts': {k: manifest[k]
            for k in ('new_points', 'new_roots', 'legacy_roots', 'pole_roots', 'ritz_solves')},
            'reason': 'explicit restart; numerical physics/config unchanged; checkpoint kernels and all previous artifacts verified'})
        return True
    return False


def plot(folder):
    # Pure saved-table rendering: no assembly, SVD, roots or cache mutation.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader((folder/'spectrum.csv').open()))
    beta = np.array([float(r['beta_deg']) for r in rows])
    l1, l2 = (np.array([float(r[k]) for r in rows]) for k in ('Lambda1', 'Lambda2'))
    for name, zoom in (('frequency_reference', False), ('first_frequency_zoom', True)):
        fig, ax = plt.subplots(figsize=(7, 4), layout='constrained')
        mask = beta >= 10 if zoom else np.ones(len(beta), dtype=bool)
        ax.plot(beta[mask], l1[mask], lw=1.5, label=r'$\Lambda_1$')
        if not zoom:
            ax.plot(beta, l2, lw=1.5, label=r'$\Lambda_2$'); ax.legend()
        ax.set(xlabel=r'$\beta$ (deg)', ylabel=r'$\Lambda$')
        ax.ticklabel_format(useOffset=False, axis='y'); ax.grid(alpha=.2)
        for ext in ('png', 'pdf'):
            # Deterministic PDFs permit plot-only without invalidating the cache.
            metadata = {'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None
            fig.savefig(folder/f'{name}.{ext}', dpi=200, metadata=metadata)
        plt.close(fig)


def limit_identity(cfg):
    """Explicit reuse of frozen verification, permitting this CLI extension only."""
    source = ROOT/cfg['source_bundle']
    old = json.loads((source/'manifest.json').read_text(encoding='utf8'))
    if sha(source/'manifest.json') != cfg['source_manifest_sha256']:
        raise ValueError('Verification manifest changed')
    if any(sha(source/f) != h for f, h in old['artifacts'].items()):
        raise ValueError('Verification artifacts changed')
    _, ident = identity(cfg)
    own = Path(__file__).relative_to(ROOT).as_posix()
    if old['identity']['files'][own] != cfg['source_cli_sha256']:
        raise ValueError('Unexpected predecessor CLI')
    if any(ident['files'][f] != h for f, h in old['identity']['files'].items() if f != own):
        raise ValueError('Frozen physical implementation/configuration changed')
    for field in ('material', 'mu', 'reference'):
        if cfg[field] != old['identity']['config'][field]:
            raise ValueError('Different physical contract')
    for field, value in old['identity']['config']['policy'].items():
        if field != 'Lambda_lower' and cfg['policy'][field] != value:
            raise ValueError('Acceptance/search policy changed')
    if ident['versions'] != old['identity']['versions']:
        raise ValueError('Environment changed; explicit review required')
    ident['files'][LIMIT_CONFIG.relative_to(ROOT).as_posix()] = sha(LIMIT_CONFIG)
    ident['source_manifest_sha256'] = sha(source/'manifest.json')
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    return key, ident


def seed_limit_points(run):
    """Copy only identical accepted cases; retain source bytes/hashes as provenance."""
    source = ROOT/run.config['source_bundle']; evidence = []
    available = [(f, json.loads(f.read_text())) for f in sorted((source/'points').glob('*.json'))]
    for label, (joint, end) in run.config['cases'].items():
        for beta in run.config['map']['parameter_grid']:
            case = physics.Case(beta, run.config['mu'], joint, end)
            contract = {**asdict(case), 'variant': 'RL_FULL'}
            matches = [(f, r) for f, r in available if all(r['case'][k] == v for k, v in contract.items())]
            if not matches:
                continue
            src, record = matches[0]
            if record['status'] != 'PASS' or not all(all(r['gates'].values()) for r in record['roots']):
                raise ValueError('Rejected predecessor point')
            key = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:16]
            path = run.folder/'points'/(key+'.json')
            proof = {'ID': label, 'beta': beta, 'source': src.relative_to(ROOT).as_posix(),
                     'sha256': sha(src), 'shape_sha256': record['shape_sha256']}
            evidence.append(proof)
            if path.exists():
                continue
            record['case'] = contract
            record['reused_from'] = proof
            record['shape_file'] = 'points/'+key+'.npz'
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src.with_suffix('.npz'), path.with_suffix('.npz'))
            write(path, record)
    write(run.folder/'source_reuse.json', evidence)
    run.catalogs = {float(k): v for k, v in json.loads((source/'pole_catalogs.json').read_text()).items()}


def limit_preflight(run, labels=('HR', 'RH', 'HH')):
    out = []; m = run.model
    for label in labels:
        joint, end = run.config['cases'][label]
        for beta in run.config['preflight_angles']:
            case = physics.Case(beta, run.config['mu'], joint, end)
            probes = []
            for lam in run.config['low_count_Lambda']:
                _, count, diag = physics.count_at(m, case, float(m.omega(lam)), run.cats(m, case), run.policy)
                if count != 0:
                    raise ArithmeticError('Low-frequency prefix not empty')
                probes.append({'Lambda': lam, 'count': count, **diag})
            record = run.solve(case)
            out.append({'ID': label, 'beta': beta, 'low_count': probes,
                        'Lambda': [r['Lambda'] for r in record['roots']], 'status': record['status']})
            write(run.folder/'preflight.json', {'status': 'IN_PROGRESS', 'points': out})
            print(f'preflight {label} beta={beta}: PASS', flush=True)
    result = {'status': 'PASS', 'points': out,
              'zero_mechanism': 'Excluded by positive rod energy and clamp1; fixed joint translation then removes arm2 rigid rotation.'}
    write(run.folder/'preflight.json', result)
    return result


def tail_ranges(rows):
    """Discrete omega range, never a Lambda range or interpolated envelope."""
    result = []
    for row in rows:
        tail = [r for r in rows if r['beta_deg'] >= row['beta_deg']]
        lo = min(tail, key=lambda r: r['omega1']); hi = max(tail, key=lambda r: r['omega1'])
        result.append({'ID': row['ID'], 'beta_deg': row['beta_deg'],
            'R': (hi['omega1']-lo['omega1'])/lo['omega1'],
            'beta_min': lo['beta_deg'], 'beta_max': hi['beta_deg'],
            'omega_min': lo['omega1'], 'omega_max': hi['omega1'], 'nodes': len(tail)})
    return result


def limit_map(run, finite=False, summarize=None):
    rows, checks, neighbours, counts = [], [], [], []; m = run.model
    for label, (joint, end) in run.config['cases'].items():
        previous = None
        for beta in run.config['map']['parameter_grid']:
            case = physics.Case(beta, run.config['mu'], joint, end)
            record = run.solve(case, neighbour=previous)
            omega = np.array([r['omega'] for r in record['roots']])
            # Independent count below the first and between adjacent roots.
            # Reused RR already has accepted brackets; evaluate only new limit cases.
            if label != 'RR':
                queries = [float(m.omega(run.policy['Lambda_lower'])),
                           *(.5*(omega[:-1]+omega[1:])).tolist()]
                ns = []
                for query, expected in zip(queries, (0, 1, 2)):
                    _, n, d = physics.count_at(m, case, query, run.cats(m, case), run.policy)
                    if n != expected:
                        raise ArithmeticError('Neighbour/count mismatch')
                    ns.append({'count': n, **d})
                counts.append({'ID': label, 'beta': beta, 'queries': ns})
            if previous is not None:
                jump = float(max(abs(omega/previous-1)))
                neighbours.append({'ID': label, 'beta': beta, 'relative_jump': jump,
                    'trigger': jump > run.policy['neighbour_relative_jump'],
                    'resolution': 'accepted count-isolated roots, mid-gap counts and physical gates; no interpolation'})
            previous = omega
            row = {'ID': label, 'beta_deg': beta}
            for i, r in enumerate(record['roots'], 1):
                row.update({f'omega{i}': r['omega'], f'Omega{i}': r['Omega'], f'Lambda{i}': r['Lambda']})
                d = r['diagnostics']
                checks.append({'ID': label, 'beta_deg': beta, 'position': i,
                    'physical': max(np.max(d['joint_scaled']), np.max(d['external_scaled'])),
                    'joint_absolute': float(np.max(abs(np.array(d['joint_absolute'])))),
                    'external_absolute': float(np.max(abs(np.array(d['external_absolute'])))),
                    'ode': d['ode_error'], 'energy': d['energy_error'], 'singular': d['singular_ratio'],
                    'condition': d['nonzero_condition'], 'mass_gram': record['mass_gram_error'],
                    'mass_error': abs(d['mass']-1), 'status': record['status']})
            row['provenance'] = 'verification_bundle' if 'reused_from' in record else ('new_finite_map' if finite else 'new_limit_map')
            row['point_file'] = next(p.name for p in (run.folder/'points').glob('*.json')
                if p.stem == hashlib.sha256(json.dumps(record['case'], sort_keys=True).encode()).hexdigest()[:16])
            rows.append(row)
        print(f'map {label}: {len(previous)} roots x 91 points accepted', flush=True)
    csvwrite(run.folder/'spectrum.csv', rows); csvwrite(run.folder/'checks.csv', checks)
    write(run.folder/'neighbour_audit.json', neighbours); write(run.folder/'midgap_counts.json', counts)
    if summarize is not None:
        return summarize(run, rows, checks, neighbours)
    if finite:
        return finite_summary(run, rows, checks, neighbours)
    tails = sum((tail_ranges([r for r in rows if r['ID'] == label]) for label in run.config['cases']), [])
    csvwrite(run.folder/'tail_ranges.csv', tails)
    source = ROOT/run.config['source_bundle']
    dense = [{'ID': 'RR_archive', **{k: float(v) for k, v in r.items() if k in ('beta_deg', 'omega1')}}
             for r in csv.DictReader((source/'spectrum.csv').open())]
    dense_tail = tail_ranges(dense); csvwrite(run.folder/'rr_archive_tail_ranges.csv', dense_tail)
    controls = [r for r in tails+dense_tail if r['beta_deg'] in run.config['tail_intervals']]
    csvwrite(run.folder/'control_intervals.csv', controls)
    order = []
    for beta in run.config['map']['parameter_grid']:
        by = {r['ID']: r for r in rows if r['beta_deg'] == beta}
        for lo, hi in [('HH', 'HR'), ('HR', 'RR'), ('HH', 'RH'), ('RH', 'RR')]:
            for i in (1, 2, 3):
                diff = by[hi][f'omega{i}']/by[lo][f'omega{i}']-1
                order.append({'beta': beta, 'position': i, 'lower': lo, 'upper': hi, 'relative_margin': diff})
                if diff < -run.config['comparison_relative_tolerance']:
                    raise ArithmeticError('Min-max violated')
    csvwrite(run.folder/'minmax.csv', order)
    extrema = []; pairs = []
    for label in run.config['cases']:
        part = [r for r in rows if r['ID'] == label]
        for a, b, c in zip(part[:-2], part[1:-1], part[2:]):
            if (b['omega1']-a['omega1'])*(c['omega1']-b['omega1']) < 0:
                extrema.append({'ID': label, 'beta': b['beta_deg'], 'omega1': b['omega1'],
                                'neighbours': [a['beta_deg'], c['beta_deg']]})
        small = min(part, key=lambda r: r['omega2']/r['omega1']-1)
        pairs.append({'ID': label, 'beta': small['beta_deg'],
            'relative_gap': small['omega2']/small['omega1']-1,
            'gap_Omega': small['Omega2']-small['Omega1']})
    summary = {'status': 'PASS', 'points': len(rows), 'interior_extrema': extrema, 'first_pair': pairs,
        'physical_maxima': {k: max(r[k] for r in checks) for k in ('physical', 'joint_absolute',
            'external_absolute', 'ode', 'energy', 'singular', 'condition', 'mass_gram', 'mass_error')},
        'minmax_minimum_margin': min(r['relative_margin'] for r in order),
        'source_points': sum(r['provenance'] == 'verification_bundle' for r in rows),
        'neighbour_triggers': [r for r in neighbours if r['trigger']],
        'range_semantics': 'discrete confirmed nodes; omega proportional to Omega, not Lambda; R(90)=0 is trivial'}
    write(run.folder/'summary.json', summary)
    return summary


def plot_limits(folder):
    """Saved CSV only. Neither eigenvalues nor state/boundary matrices evaluated."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = list(csv.DictReader((folder/'spectrum.csv').open()))
    tails = list(csv.DictReader((folder/'tail_ranges.csv').open()))
    colors = dict(zip(('RR', 'HR', 'RH', 'HH'), ('#245a81', '#d37e14', '#32916c', '#a84e9b')))
    for name, table, field, ylabel in [('first_frequency_limits', rows, 'Lambda1', r'$\Lambda_1$'),
            ('remaining_range_limits', tails, 'R', r'$100\,\mathcal{R}(\beta)$, %')]:
        fig, ax = plt.subplots(figsize=(7, 4.2), layout='constrained')
        for label, color in colors.items():
            part = [r for r in table if r['ID'] == label]
            factor = 100 if field == 'R' else 1
            ax.plot([float(r['beta_deg']) for r in part], [factor*float(r[field]) for r in part],
                    lw=1.5, label=label, color=color)
        ax.set(xlabel=r'$\beta$, deg', ylabel=ylabel, xlim=(0, 90)); ax.legend(); ax.grid(alpha=.2)
        if field == 'R':
            # Display-only scale: retain the exact zero at 90 degrees without
            # hiding the small, nonzero tail ranges. Not a plateau criterion.
            ax.set_yscale('symlog', linthresh=.001)
            ax.set_ylim(bottom=0)
        for ext in ('png', 'pdf'):
            fig.savefig(folder/f'{name}.{ext}', dpi=220, bbox_inches='tight',
                        metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
        plt.close(fig)
    # Separate y scales expose actual variation; descriptive interval, no plateau criterion.
    fig, axs = plt.subplots(2, 2, figsize=(9, 6), layout='constrained')
    for ax, (label, color) in zip(axs.flat, colors.items()):
        part = [r for r in rows if r['ID'] == label and float(r['beta_deg']) >= 10]
        ax.plot([float(r['beta_deg']) for r in part], [float(r['Lambda1']) for r in part], color=color)
        ax.set(title=label, xlabel=r'$\beta$, deg', ylabel=r'$\Lambda_1$'); ax.grid(alpha=.2)
        ax.ticklabel_format(useOffset=False, axis='y')
    for ext in ('png', 'pdf'):
        fig.savefig(folder/f'first_frequency_limits_zoom.{ext}', dpi=220, bbox_inches='tight',
                    metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
    plt.close(fig)


def saved_limit_pair_tables(folder):
    """Read accepted tables/records and hashes only; no physical operator calls."""
    def checked_bundle(path):
        manifest = json.loads((path/'manifest.json').read_text(encoding='utf8'))
        if manifest['status'] != 'PASS':
            raise ValueError(f'Unaccepted source: {path}')
        if any(sha(path/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError(f'Source artifact mismatch: {path}')
        return manifest
    manifest = checked_bundle(folder)
    cfg = manifest['identity']['config']
    expected_material = dict(E=1., rho=1., nu=.3, epsilon=.005, aspect=4., l=1., ks=5/6)
    cases = {'RR': ['RIGID', 'RIGID'], 'HR': [0, 'RIGID'],
             'RH': ['RIGID', 0], 'HH': [0, 0]}
    if (cfg['material'] != expected_material or cfg['mu'] != .01 or cfg['cases'] != cases
            or cfg['map']['spectrum_semantics'] != 'sorted_positions'):
        raise ValueError('Unexpected limit-map physical contract')
    rows = list(csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')))
    tables = {label: [] for label in cases}
    for row in rows:
        record = json.loads((folder/'points'/row['point_file']).read_text(encoding='utf8'))
        label = row['ID']; beta = float(row['beta_deg'])
        if (label not in cases or record['status'] != 'PASS' or record['count']['status'] != 'PASS'
                or record['case'] != dict(beta=beta, mu=.01, k_joint=cases[label][0],
                                         k_end=cases[label][1], variant='RL_FULL')):
            raise ValueError('Unaccepted or incompatible saved point')
        for position, root in enumerate(record['roots'], 1):
            if root['sorted_position'] != position or not all(root['gates'].values()):
                raise ValueError('Saved modal gate failed')
            for field in ('omega', 'Omega', 'Lambda'):
                if float(row[f'{field}{position}']) != root[field]:
                    raise ValueError('CSV/accepted-root mismatch')
            if not math.isclose(root['Lambda']**2, root['Omega'], rel_tol=2e-14):
                raise ValueError('Frequency normalization mismatch')
        values = [float(row[f'Lambda{i}']) for i in (1, 2, 3)]
        if not all(math.isfinite(x) for x in values) or not 0 < values[0] < values[1] < values[2]:
            raise ValueError('Invalid sorted prefix')
        tables[label].append({k: float(row[k]) for k in
                              ('beta_deg', 'Lambda1', 'Lambda2', 'omega1', 'omega2')})
    for table in tables.values():
        table.sort(key=lambda r: r['beta_deg'])
        if [r['beta_deg'] for r in table] != list(range(91)):
            raise ValueError('Missing/duplicate grid point; no interpolation or solve allowed')
    archive = ROOT/cfg['source_bundle']
    old = checked_bundle(archive)
    if (sha(archive/'manifest.json') != cfg['source_manifest_sha256']
            or old['identity']['config']['material'] != expected_material
            or old['identity']['config']['mu'] != .01):
        raise ValueError('RR archive contract mismatch')
    dense = list(csv.DictReader((archive/'spectrum.csv').open(encoding='utf8')))
    if any(r['count_status'] != 'PASS' for r in dense):
        raise ValueError('Unaccepted dense RR point')
    zoom_rr = [{k: float(r[k]) for k in ('beta_deg', 'Lambda1', 'Lambda2', 'omega1', 'omega2')}
               for r in dense if 0 <= float(r['beta_deg']) <= 15]
    zoom_rr.sort(key=lambda r: r['beta_deg'])
    if not zoom_rr or zoom_rr[0]['beta_deg'] != 0 or zoom_rr[-1]['beta_deg'] != 15:
        raise ValueError('Incomplete RR zoom interval')
    for row in zoom_rr:
        if not 0 < row['Lambda1'] < row['Lambda2'] or not all(math.isfinite(v) for v in row.values()):
            raise ValueError('Invalid dense RR frequencies')
    return tables, zoom_rr, [(folder, manifest), (archive, old)]


def plot_limit_pairs(folder, russian=False):
    """Presentation-only follow-up. Preserve original figures, CSV, NPZ, manifests."""
    folder = Path(folder).resolve()
    tables, zoom_rr, sources = saved_limit_pair_tables(folder)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = folder/'visualizations'/('lambda12_ru' if russian else 'lambda12')
    output.mkdir(parents=True, exist_ok=True)
    colors = dict(zip(('RR', 'HR', 'RH', 'HH'), ('#245a81', '#d37e14', '#32916c', '#a84e9b')))
    names = {'RR': 'Жёсткий узел — заделка', 'HR': 'Шарнир — заделка',
             'RH': 'Жёсткий узел — шарнир', 'HH': 'Шарнир — шарнир'}
    images = []; started = time.perf_counter()
    def finish(fig, name):
        for ext in ('png', 'pdf'):
            path = output/f'{name}.{ext}'
            fig.savefig(path, dpi=300, metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
            images.append(path.name)
        plt.close(fig)
    def axes_style(ax, zoom=False):
        ax.set(xlabel='β, °', ylabel='Λ', xlim=(0, 15 if zoom else 90), ylim=(1.35, 4.9))
        ax.set_xticks(range(0, 16, 3) if zoom else range(0, 91, 15))
        ax.grid(alpha=.18, linewidth=.6)
        ax.spines[['top', 'right']].set_visible(False)
    def pair(ax, table):
        beta = [r['beta_deg'] for r in table]
        ax.plot(beta, [r['Lambda1'] for r in table], color=colors['RR'], lw=1.7, label='Λ₁')
        ax.plot(beta, [r['Lambda2'] for r in table], color=colors['HR'], lw=1.7, ls='--', label='Λ₂')
        ax.legend(frameon=False)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.labelsize': 12,
                         'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42}):
        for label, table in tables.items():
            fig, ax = plt.subplots(figsize=(6.4, 4.2), layout='constrained')
            axes_style(ax); pair(ax, table)
            if russian:
                ax.set_title(names[label], fontsize=12)
            finish(fig, 'lambda12_'+label)
        fig, ax = plt.subplots(figsize=(6.4, 4.2), layout='constrained')
        axes_style(ax); ax.set_ylabel('Λ₂'); ax.set_ylim(3.25, 4.9)
        for (label, table), style in zip(tables.items(), ('-', '--', '-.', ':')):
            ax.plot([r['beta_deg'] for r in table], [r['Lambda2'] for r in table],
                    color=colors[label], ls=style, lw=1.7, label=names[label].lower() if russian else label)
        ax.legend(frameon=False); finish(fig, 'lambda2_comparison')
        fig, axs = plt.subplots(2, 2, figsize=(10, 7), layout='constrained', sharey=True)
        for ax, (label, table) in zip(axs.flat, tables.items()):
            axes_style(ax, zoom=True); ax.set_title(names[label] if russian else label, fontsize=12)
            pair(ax, zoom_rr if label == 'RR' else [r for r in table if r['beta_deg'] <= 15])
        finish(fig, 'lambda12_small_angles')
    observations = {}
    for label, table in tables.items():
        closest = min(table, key=lambda r: r['Lambda2']-r['Lambda1'])
        tail = [r for r in table if r['beta_deg'] >= 10]
        observations[label] = {'minimum_grid_Lambda_gap_beta': closest['beta_deg'],
            'minimum_grid_Lambda_gap': closest['Lambda2']-closest['Lambda1'],
            'Lambda2_at_10': tail[0]['Lambda2'], 'Lambda2_at_90': tail[-1]['Lambda2'],
            'Lambda2_relative_range_10_90': max(r['Lambda2'] for r in tail)/min(r['Lambda2'] for r in tail)-1}
    # Hash verification, rather than a new physical check, establishes preservation.
    for path, manifest in sources:
        if any(sha(path/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Source changed during rendering')
    record = {'status': 'PASS', 'calculation_mode': 'plot_only', 'spectrum_semantics': 'sorted_positions',
        'new_eigenvalue_solves': 0, 'matrix_assemblies': 0, 'SVD_calls': 0, 'shape_recoveries': 0,
        'sources': [{'bundle': path.relative_to(ROOT).as_posix(), 'manifest_sha256': sha(path/'manifest.json'),
                     'spectrum_sha256': sha(path/'spectrum.csv'), 'verified_artifacts': len(m['artifacts'])}
                    for path, m in sources],
        'material': sources[0][1]['identity']['config']['material'], 'mu': .01, 'variant': 'RL_FULL',
        'normalization': 'Lambda^4=rho*A*omega^2*l^4/(E*I); Omega=Lambda^2',
        'common_grid': {'beta_degrees': [0, 90], 'step': 1, 'accepted_points_per_case': 91},
        'small_angle_RR_archive_points': len(zoom_rr),
        'case_colors': colors, 'individual_line_styles': {'Lambda1': 'solid blue', 'Lambda2': 'dashed orange'},
        'PNG_dpi': 300, 'PDF': 'vector', 'render_seconds': time.perf_counter()-started,
        'preserved_first_frequency_comparison': {f: sha(folder/f) for f in
            ('first_frequency_limits.png', 'first_frequency_limits.pdf')},
        'images': {f: sha(output/f) for f in images}, 'observations_from_saved_values': observations,
        'CLI_sha256': sha(Path(__file__)),
        'command': 'python scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py '
                   '--preset limit-angles --plot-only '+folder.relative_to(ROOT).as_posix()+' --plot-set '+('lambda12-ru' if russian else 'lambda12'),
        'interpretation': 'No new MAC, tracking, veering proof, root refinement or plateau threshold.'}
    if russian:
        record['presentation_change'] = {'case_names': names, 'original_images_preserved': True,
            'original_manifest_sha256': sha(folder/'visualizations/lambda12/manifest.json'),
            'scope': 'Russian case titles/legend only; identical frequencies, styles, axes and normalization.'}
    write(output/'manifest.json', record)
    print(output, flush=True)


def limit_extremum_check(run):
    """Triggered RH maximum only: two midpoints and one checked parabola vertex.

    The original integer-grid range remains the main statistic. No interpolated
    value is reported as a computed eigenfrequency, and no angular optimizer runs.
    """
    grid = list(csv.DictReader((run.folder/'spectrum.csv').open()))
    rh = [r for r in grid if r['ID'] == 'RH']
    maximum = max(rh, key=lambda r: float(r['omega1']))
    beta = float(maximum['beta_deg'])
    if beta != 84.:
        raise ValueError('Targeted diagnostic trigger changed; review required')
    neighbours = [float(maximum[f'omega{i}']) for i in (1, 2, 3)]
    records = []
    for angle in (83.5, 84.5):
        records.append(run.solve(physics.Case(angle, run.config['mu'], 'RIGID', 0), neighbour=neighbours))
    left, center, right = records[0]['roots'][0]['omega'], neighbours[0], records[1]['roots'][0]['omega']
    offset = .25*(left-right)/(left-2*center+right)
    if not abs(offset) < .5:
        raise ArithmeticError('Local extremum not bracketed by diagnostic points')
    vertex = run.solve(physics.Case(84.+offset, run.config['mu'], 'RIGID', 0), neighbour=neighbours)
    records.append(vertex)
    best = max([center]+[r['roots'][0]['omega'] for r in records])
    changes = []
    for start in run.config['tail_intervals']:
        tail = [float(r['omega1']) for r in rh if float(r['beta_deg']) >= start]
        changes.append({'beta_start': start, 'R_grid': (max(tail)-min(tail))/min(tail),
            'R_with_checked_points': (best-min(tail))/min(tail),
            'absolute_R_increment': (best-max(tail))/min(tail)})
    result = {'status': 'PASS', 'trigger': 'RH interior maximum on integer grid at 84 degrees',
        'fixed_protocol': '83.5,84.5 then one actual three-root solve at local quadratic vertex; no further refinement',
        'new_angles': [r['case']['beta'] for r in records],
        'omega1': [r['roots'][0]['omega'] for r in records],
        'checked_vertex_beta': vertex['case']['beta'],
        'relative_increase_over_grid_max': best/center-1, 'ranges': changes,
        'qualification': 'checked sampled maximum; not a certified continuous extremum or plateau angle'}
    write(run.folder/'rh_extremum_check.json', result)
    return result


def limits_main(args):
    cfg = json.loads(LIMIT_CONFIG.read_text()); key, ident = limit_identity(cfg)
    folder = OUTPUT/'limit_angle_maps'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = args.verify_cache or folder
        manifest = json.loads((target/'manifest.json').read_text())
        if manifest['identity'] != ident or any(sha(target/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Limit map cache mismatch; no overwrite')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True); write(folder/'contract.json', ident)
    run = Run(cfg, folder); seed_limit_points(run); start = time.perf_counter()
    inherited = {'new_points': 0, 'new_roots': 0}
    if args.reuse_points_from:
        source = args.reuse_points_from
        old = json.loads((source/'manifest.json').read_text())
        own = Path(__file__).relative_to(ROOT).as_posix()
        expected = {k: v for k, v in ident.items() if k != 'files'}
        previous = {k: v for k, v in old['identity'].items() if k != 'files'}
        if previous != expected or any(ident['files'][f] != h for f, h in old['identity']['files'].items() if f != own):
            raise ValueError('Limit checkpoint physical contract changed')
        if any(sha(source/f) != h for f, h in old['artifacts'].items()):
            raise ValueError('Limit checkpoint artifacts changed')
        for path in (source/'points').glob('*'):
            shutil.copyfile(path, folder/'points'/path.name)
        shutil.copyfile(source/'preflight.json', folder/'preflight.json')
        inherited = {k: old[k] for k in ('new_points', 'new_roots')}
        write(folder/'reused_limit_checkpoint.json', {'source': str(source),
            'manifest_sha256': sha(source/'manifest.json'), 'inherited_computations': inherited,
            'reason': 'Only CLI adds triggered RH extremum diagnostics; all physical sources, configuration, saved gates and artifacts match.'})
    logpath = folder/'run_ledger.json'; ledger = json.loads(logpath.read_text()) if logpath.exists() else []
    if not args.reuse_points_from and (folder/'reused_limit_checkpoint.json').exists():
        inherited = json.loads((folder/'reused_limit_checkpoint.json').read_text())['inherited_computations']
    try:
        preflight = folder/'preflight.json'
        if not preflight.exists() or json.loads(preflight.read_text())['status'] != 'PASS':
            limit_preflight(run)
        if args.stage != 'preflight':
            summary = limit_map(run)
            summary['RH_extremum'] = limit_extremum_check(run)
            write(folder/'summary.json', summary); plot_limits(folder)
    except Exception as exc:
        write(folder/'failure.json', {'error': repr(exc), 'new_points': run.new_points, 'new_roots': run.new_roots})
        raise
    finally:
        ledger.append({'command': sys.argv, 'new_points': run.new_points, 'new_roots': run.new_roots,
            'reused_point_calls': run.reused_points, 'seconds': time.perf_counter()-start})
        write(logpath, ledger)
    if args.stage != 'preflight':
        manifest = {'status': summary['status'], 'identity': ident, 'git_head': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True).strip(), 'new_points': inherited['new_points']+sum(r['new_points'] for r in ledger),
            'new_roots': inherited['new_roots']+sum(r['new_roots'] for r in ledger), 'source_points': summary['source_points'],
            'pole_roots': 0, 'ritz_solves': 0, 'complex_solves': 0,
            'artifacts': {p.relative_to(folder).as_posix(): sha(p) for p in folder.rglob('*')
                          if p.is_file() and p.name != 'manifest.json'}}
        write(folder/'manifest.json', manifest)
    print(folder, flush=True)


def finite_identity(cfg):
    """Same arm physics, tolerances and source data; a six-case parameter preset."""
    _, ident = limit_identity(cfg)
    source = ROOT/cfg['limit_bundle']
    old = json.loads((source/'manifest.json').read_text(encoding='utf8'))
    if sha(source/'manifest.json') != cfg['limit_manifest_sha256'] or old['status'] != 'PASS':
        raise ValueError('Limit source manifest mismatch')
    if any(sha(source/f) != h for f, h in old['artifacts'].items()):
        raise ValueError('Limit source artifact mismatch')
    for field in ('material', 'mu', 'policy', 'map'):
        if cfg[field] != old['identity']['config'][field]:
            raise ValueError('Finite map changed frozen physical/numerical contract')
    own = Path(__file__).relative_to(ROOT).as_posix()
    if any(ident['files'][f] != h for f, h in old['identity']['files'].items() if f != own):
        raise ValueError('Limit-map physical dependencies changed')
    ident['files'][FINITE_CONFIG.relative_to(ROOT).as_posix()] = sha(FINITE_CONFIG)
    return hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16], ident


def finite_summary(run, rows, checks, neighbours):
    """Compare accepted sorted positions, never infer mode descendants."""
    limits = list(csv.DictReader((ROOT/run.config['limit_bundle']/'spectrum.csv').open()))
    lookup = {(r['ID'], float(r['beta_deg'])): r for r in limits+rows}
    order = []
    for series, chain in run.config['series'].items():
        for beta in run.config['map']['parameter_grid']:
            for low, high in zip(chain[:-1], chain[1:]):
                for i in (1, 2, 3):
                    margin = float(lookup[high, beta][f'omega{i}'])/float(lookup[low, beta][f'omega{i}'])-1
                    order.append(dict(series=series, beta_deg=beta, position=i, lower=low,
                                      upper=high, relative_margin=margin))
                    if margin < -run.config['comparison_relative_tolerance']:
                        raise ArithmeticError('Finite-spring min-max ordering failed')
    csvwrite(run.folder/'minmax.csv', order)
    tails = sum((tail_ranges([r for r in rows if r['ID'] == label]) for label in run.config['cases']), [])
    csvwrite(run.folder/'control_intervals.csv', [r for r in tails if r['beta_deg'] in run.config['tail_intervals']])
    comparisons = []
    for label in run.config['cases']:
        part = [r for r in rows if r['ID'] == label]
        closest = min(part, key=lambda r: r['Lambda2']-r['Lambda1'])
        gaps = {i: [1-r[f'omega{i}']/float(lookup['RR', r['beta_deg']][f'omega{i}']) for r in part] for i in (1, 2)}
        comparisons.append({'ID': label, 'minimum_sampled_Lambda_gap_beta': closest['beta_deg'],
            'minimum_sampled_Lambda_gap': closest['Lambda2']-closest['Lambda1'],
            'relative_omega_gap_at_that_point': closest['omega2']/closest['omega1']-1,
            'Lambda1_at_0_10_20_90': [next(r['Lambda1'] for r in part if r['beta_deg'] == b) for b in (0, 10, 20, 90)],
            'Lambda2_at_0_10_20_90': [next(r['Lambda2'] for r in part if r['beta_deg'] == b) for b in (0, 10, 20, 90)],
            'max_relative_omega_deficit_to_RR': {str(i): max(v) for i, v in gaps.items()},
            'relative_omega_deficit_to_RR_at_90': {str(i): v[-1] for i, v in gaps.items()}})
    summary = {'status': 'PASS', 'points': len(rows), 'first_pair_and_limits': comparisons,
        'source_points': sum(r['provenance'] == 'verification_bundle' for r in rows),
        'physical_maxima': {k: max(r[k] for r in checks) for k in ('physical', 'joint_absolute',
            'external_absolute', 'ode', 'energy', 'singular', 'condition', 'mass_gram', 'mass_error')},
        'minmax_minimum_margin': min(r['relative_margin'] for r in order),
        'neighbour_triggers': [r for r in neighbours if r['trigger']],
        'qualification': 'Finite 1-degree grid; sorted positions; no modal tracking, plateau angle or continuous gap minimum.'}
    write(run.folder/'summary.json', summary)
    return summary


def plot_finite_maps(folder):
    """Read-only spectral input: no solver, matrix, SVD or shape recovery calls."""
    folder = Path(folder).resolve()
    manifest_path = folder/'manifest.json'
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding='utf8'))
        if manifest['status'] != 'PASS' or any(sha(folder/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Finite plot source failed provenance check')
    cfg = json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    rows = list(csv.DictReader((folder/'spectrum.csv').open()))
    limits = ROOT/cfg['limit_bundle']
    old = json.loads((limits/'manifest.json').read_text(encoding='utf8'))
    if sha(limits/'manifest.json') != cfg['limit_manifest_sha256'] or any(sha(limits/f) != h for f, h in old['artifacts'].items()):
        raise ValueError('Plot limit provenance failed')
    allrows = rows+list(csv.DictReader((limits/'spectrum.csv').open()))
    tables = {}
    for label in list(cfg['cases'])+['HR', 'RH', 'RR']:
        part = [{k: float(r[k]) for k in ('beta_deg', 'Lambda1', 'Lambda2', 'Lambda3')}
                for r in allrows if r['ID'] == label]
        part.sort(key=lambda r: r['beta_deg'])
        if [r['beta_deg'] for r in part] != cfg['map']['parameter_grid']:
            raise ValueError('Incomplete plot grid; do not interpolate')
        if any(not 0 < r['Lambda1'] < r['Lambda2'] < r['Lambda3'] for r in part):
            raise ValueError('Invalid sorted data')
        tables[label] = part
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = folder/'figures'; output.mkdir(exist_ok=True)
    images = []
    def finish(fig, name):
        for ext in ('png', 'pdf'):
            path = output/f'{name}.{ext}'
            fig.savefig(path, dpi=300, metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
            images.append(path.name)
        plt.close(fig)
    def style(ax, zoom=False):
        ax.set(xlabel='β, °', ylabel='Λ', xlim=(0, 20 if zoom else 90), ylim=(1.35, 4.9))
        ax.set_xticks(range(0, 21, 5) if zoom else range(0, 91, 15))
        ax.grid(alpha=.18, linewidth=.6); ax.spines[['top', 'right']].set_visible(False)
    def name(label):
        joint, end = cfg['cases'][label]
        return (rf'Внутренняя пружина: $\kappa_{{joint}}={joint:g}$; заделка' if end == 'RIGID'
                else rf'Жёсткий узел; внешняя пружина: $\kappa_{{end}}={end:g}$')
    def pair(ax, label, zoom=False):
        part = [r for r in tables[label] if not zoom or r['beta_deg'] <= 20]
        for i, color, ls in ((1, '#245a81', '-'), (2, '#d37e14', '--')):
            ax.plot([r['beta_deg'] for r in part], [r[f'Lambda{i}'] for r in part],
                    color=color, ls=ls, lw=1.7, label='Λ₁' if i == 1 else 'Λ₂')
        ax.legend(frameon=False); ax.set_title(name(label), fontsize=11)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.labelsize': 12,
                         'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42}):
        for series, chain in cfg['series'].items():
            fig, ax = plt.subplots(figsize=(7, 4.6), layout='constrained'); style(ax)
            ax.set_ylabel('Λ₁')
            for label, color, ls in zip(chain, ('#777777', '#245a81', '#d37e14', '#32916c', '#222222'), (':', '-', '-', '-', '--')):
                part = tables[label]
                legend = {'HR': 'шарнир — заделка', 'RH': 'жёсткий узел — шарнир', 'RR': 'жёсткий узел — заделка'}.get(label)
                if legend is None:
                    joint, end = cfg['cases'][label]
                    legend = rf'$\kappa_{{joint}}={joint:g}$' if series == 'A' else rf'$\kappa_{{end}}={end:g}$'
                ax.plot([r['beta_deg'] for r in part], [r['Lambda1'] for r in part],
                        color=color, ls=ls, lw=1.6, label=legend)
            ax.legend(frameon=False, fontsize=10); finish(fig, 'lambda1_series_'+series)
            fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), layout='constrained', sharey=True)
            for ax, label in zip(axs, chain[1:-1]):
                style(ax, True); pair(ax, label, True)
            finish(fig, 'lambda12_small_angles_'+series)
        for label in cfg['cases']:
            fig, ax = plt.subplots(figsize=(7, 4.6), layout='constrained')
            style(ax); pair(ax, label); finish(fig, 'lambda12_'+label)
    write(folder/'plot_manifest.json', {'calculation_mode': 'plot_only', 'new_eigenvalue_solves': 0,
        'spectrum_sha256': sha(folder/'spectrum.csv'), 'limit_spectrum_sha256': sha(limits/'spectrum.csv'),
        'PNG_dpi': 300, 'PDF': 'vector', 'spectrum_semantics': 'sorted_positions',
        'images': {f: sha(output/f) for f in images}})


def finite_main(args):
    cfg = json.loads(FINITE_CONFIG.read_text()); key, ident = finite_identity(cfg)
    folder = OUTPUT/'finite_one_parameter_angle_maps'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = (args.verify_cache or folder).resolve()
        old = json.loads((target/'manifest.json').read_text())
        if old['identity'] != ident or any(sha(target/f) != h for f, h in old['artifacts'].items()):
            raise ValueError('Finite cache mismatch; no overwrite')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True)
    write(folder/'contract.json', ident)
    run = Run(cfg, folder); seed_limit_points(run); start = time.perf_counter()
    logpath = folder/'run_ledger.json'; ledger = json.loads(logpath.read_text()) if logpath.exists() else []
    try:
        preflight = folder/'preflight.json'
        if not preflight.exists() or json.loads(preflight.read_text())['status'] != 'PASS':
            limit_preflight(run, tuple(cfg['cases']))
        if args.stage != 'preflight':
            summary = limit_map(run, finite=True)
            plot_finite_maps(folder)
    except Exception as exc:
        write(folder/'failure.json', {'error': repr(exc), 'new_points': run.new_points, 'new_roots': run.new_roots})
        raise
    finally:
        ledger.append({'command': sys.argv, 'new_points': run.new_points, 'new_roots': run.new_roots,
            'reused_point_calls': run.reused_points, 'seconds': time.perf_counter()-start})
        write(logpath, ledger)
    if args.stage != 'preflight':
        _, after = finite_identity(cfg)
        if after != ident:
            raise ValueError('Sources changed during calculation')
        manifest = {'status': summary['status'], 'identity': ident, 'git_head': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True).strip(),
            'git_status': subprocess.check_output(['git', 'status', '--short'], text=True),
            'new_points': sum(r['new_points'] for r in ledger), 'new_roots': sum(r['new_roots'] for r in ledger),
            'source_points': summary['source_points'], 'points': summary['points'],
            'pole_roots': run.pole_roots, 'ritz_solves': 0, 'complex_solves': 0,
            'artifacts': {p.relative_to(folder).as_posix(): sha(p) for p in folder.rglob('*')
                          if p.is_file() and p.name != 'manifest.json'}}
        write(folder/'manifest.json', manifest)
    print(folder, flush=True)


def saturation_sources(cfg):
    """Saved scalar spectra only. Hashes/gates establish authority, no new physics."""
    tables, contracts = {}, []
    expected = dict(E=1., rho=1., nu=.3, epsilon=.005, aspect=4., l=1., ks=5/6)
    for source in cfg['sources']:
        folder = ROOT/source['bundle']
        manifest = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if sha(folder/'manifest.json') != source['manifest_sha256'] or manifest['status'] != 'PASS':
            raise ValueError('Saturation source manifest mismatch')
        if any(sha(folder/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Saturation source artifact mismatch')
        old = manifest['identity']['config']; contracts.append(old)
        if old['material'] != expected or old['mu'] != .01 or old['map']['spectrum_semantics'] != 'sorted_positions':
            raise ValueError('Saturation physical contract mismatch')
        for row in csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')):
            record = json.loads((folder/'points'/row['point_file']).read_text(encoding='utf8'))
            if record['status'] != 'PASS' or record['count']['status'] != 'PASS':
                raise ValueError('Unaccepted source point')
            case = record['case']; label = row['ID']; beta = float(row['beta_deg'])
            if case != dict(beta=beta, mu=.01, k_joint=old['cases'][label][0], k_end=old['cases'][label][1], variant='RL_FULL'):
                raise ValueError('Source point belongs to different physical case')
            for i, r in enumerate(record['roots'], 1):
                if not all(r['gates'].values()) or any(float(row[f'{f}{i}']) != r[f] for f in ('omega', 'Omega', 'Lambda')):
                    raise ValueError('CSV differs from accepted prefix')
                if not (math.isclose(r['Lambda']**2, r['Omega'], rel_tol=2e-14)
                        and math.isclose(r['omega']*200., r['Omega'], rel_tol=2e-14)):
                    raise ValueError('Source normalization mismatch')
            tables.setdefault(label, []).append({k: float(row[k]) for k in ('beta_deg', 'omega1', 'Omega1', 'Lambda1')})
    for label, table in tables.items():
        table.sort(key=lambda r: r['beta_deg'])
        if [r['beta_deg'] for r in table] != cfg['beta_deg']:
            raise ValueError('Missing source nodes; no interpolation or spectrum regeneration')
    if len(tables) != 10 or any(contracts[0][key] != contracts[1][key] for key in ('material', 'mu', 'policy', 'map')):
        raise ValueError('Incompatible source maps')
    return tables


def saturation_measures(label, table):
    """Differences in degrees, unmodified nonmonotone data and exact discrete tails."""
    beta = np.array([r['beta_deg'] for r in table]); w = np.array([r['omega1'] for r in table])
    if not np.array_equal(beta, np.arange(91)) or np.any(w <= 0) or max(w) <= w[0]:
        raise ValueError('Expected positive growing spectrum on integer 0..90 grid')
    f = (w-w[0])/(max(w)-w[0]); logw = np.log(w)
    s1 = np.gradient(logw, beta, edge_order=2)
    tails = tail_ranges([{'ID': label, **r} for r in table])
    rows = []
    for i, (point, tail) in enumerate(zip(table, tails)):
        row = {'ID': label, **point, 'F': float(f[i]), 'R': tail['R'],
               'tail_beta_min': tail['beta_min'], 'tail_beta_max': tail['beta_max'],
               'tail_omega_min': tail['omega_min'], 'tail_omega_max': tail['omega_max'],
               'S_h1_per_degree': float(s1[i]), 'abs_S_h1': float(abs(s1[i]))}
        for h in (2, 3):
            row[f'S_h{h}_per_degree'] = float((logw[i+h]-logw[i-h])/(2*h)) if h <= i <= 90-h else None
            row[f'abs_S_h{h}'] = abs(row[f'S_h{h}_per_degree']) if row[f'S_h{h}_per_degree'] is not None else None
        rows.append(row)
    return rows


def saturation_transitions(rows, window):
    """Grid knot search for descriptive geometry; no frequency solver or smoother."""
    part = [r for r in rows if r['beta_deg'] <= window]
    beta = np.array([r['beta_deg'] for r in part]); x = beta/window
    y = np.array([r['F'] for r in part]); label = part[0]['ID']
    chord = y[0]+(y[-1]-y[0])*x
    distance = (y-chord)/math.sqrt(1+(y[-1]-y[0])**2)
    geometric = int(np.argmax(abs(distance)))
    profiles, best = [], None
    for knot in range(1, window):
        design = np.column_stack((np.ones(len(x)), x, np.maximum(0., x-knot/window)))
        q, r = np.linalg.qr(design, mode='reduced')
        coef = np.linalg.solve(r, q.T@y)
        prediction = design@coef
        sse = float(np.sum((prediction-y)**2))
        profile = {'ID': label, 'window_deg': window, 'knot_deg': knot, 'SSE_F': sse,
                   'RMSE_F': math.sqrt(sse/len(x)), 'intercept': float(coef[0]),
                   'slope_left_F_per_degree': float(coef[1]/window),
                   'slope_right_F_per_degree': float((coef[1]+coef[2])/window)}
        profiles.append(profile)
        if best is None or sse < best[0]['SSE_F']:
            best = profile, prediction
    best_profile, prediction = best
    transitions = []
    for method, i in [('chord', geometric), ('segmented', best_profile['knot_deg'])]:
        point = part[i]
        transitions.append({'ID': label, 'window_deg': window, 'method': method,
            'beta_deg': int(beta[i]), 'F': point['F'], 'Lambda1': point['Lambda1'],
            'Omega1': point['Omega1'], 'omega1': point['omega1'], 'R_to_90': point['R'],
            'S_h1_per_degree': point['S_h1_per_degree'],
            'S_h2_per_degree': point['S_h2_per_degree'], 'S_h3_per_degree': point['S_h3_per_degree'],
            'tail_beta_min': point['tail_beta_min'], 'tail_beta_max': point['tail_beta_max'],
            'geometric_signed_distance': float(distance[i]) if method == 'chord' else None,
            'fit_RMSE_F': best_profile['RMSE_F'] if method == 'segmented' else None})
    curves = [{'ID': label, 'window_deg': window, 'beta_deg': int(b), 'F': float(v),
               'chord_F': float(c), 'signed_distance': float(d), 'segmented_F': float(pred)}
              for b, v, c, d, pred in zip(beta, y, chord, distance, prediction)]
    return transitions, profiles, curves


def plot_saturation(folder):
    """Render already derived CSVs only; even regressions are excluded here."""
    folder = Path(folder).resolve()
    manifest_file = folder/'manifest.json'
    if manifest_file.exists():
        m = json.loads(manifest_file.read_text(encoding='utf8'))
        if any(sha(folder/f) != h for f, h in m['artifacts'].items()):
            raise ValueError('Postprocessing cache corrupted')
    rows = list(csv.DictReader((folder/'diagnostics.csv').open(encoding='utf8')))
    transitions = list(csv.DictReader((folder/'transitions.csv').open(encoding='utf8')))
    fits = list(csv.DictReader((folder/'fit_curves.csv').open(encoding='utf8')))
    cfg = json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    groups = [('Пределы', ['RR', 'HR', 'RH', 'HH']), ('Внутренняя пружина', ['J0.1_ER', 'J1_ER', 'J10_ER']),
              ('Внешняя пружина', ['JR_E0.1', 'JR_E1', 'JR_E10'])]
    names = {'RR': 'жёсткий узел — заделка', 'HR': 'шарнир — заделка',
             'RH': 'жёсткий узел — шарнир', 'HH': 'шарнир — шарнир'}
    colors = dict(zip(names, ('#245a81', '#d37e14', '#32916c', '#a84e9b')))
    for _, ids in groups[1:]:
        for label, k, color in zip(ids, (.1, 1, 10), ('#245a81', '#d37e14', '#32916c')):
            names[label] = rf'$\kappa_{{joint}}={k:g}$' if label.endswith('_ER') else rf'$\kappa_{{end}}={k:g}$'
            colors[label] = color
    tables = {label: [r for r in rows if r['ID'] == label] for label in names}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = folder/'figures'; output.mkdir(exist_ok=True); images = []
    def style(ax, ylabel, xmax=90):
        ax.set(xlabel='β, °', ylabel=ylabel, xlim=(0, xmax))
        ax.grid(alpha=.18, linewidth=.6); ax.spines[['top', 'right']].set_visible(False)
    def lines(ax, ids, field, factor=1., low=0, high=90):
        for label in ids:
            table = [r for r in tables[label] if low <= float(r['beta_deg']) <= high]
            ax.plot([float(r['beta_deg']) for r in table], [factor*float(r[field]) for r in table],
                    color=colors[label], lw=1.6, label=names[label])
    def save(fig, name):
        for ext in ('png', 'pdf'):
            path = output/f'{name}.{ext}'
            fig.savefig(path, dpi=300, metadata={'CreationDate': None, 'ModDate': None} if ext == 'pdf' else None)
            images.append(path.name)
        plt.close(fig)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.labelsize': 11,
                         'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42}):
        for (title, ids), suffix in zip(groups, ('limits', 'internal', 'external')):
            fig, ax = plt.subplots(figsize=(7, 4.6), layout='constrained')
            style(ax, 'F'); lines(ax, ids, 'F'); ax.set_ylim(-.02, 1.04)
            ax.legend(frameon=False, loc='lower right'); save(fig, 'gain_fraction_'+suffix)
        fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), layout='constrained', sharey=True)
        for ax, (title, ids) in zip(axs, groups):
            style(ax, 'R, %'); lines(ax, ids, 'R', 100)
            ax.set_yscale('symlog', linthresh=cfg['display_only_scales']['R_percent_symlog_linthresh'])
            ax.set_title(title); ax.legend(frameon=False, fontsize=8)
        save(fig, 'remaining_range')
        fig, axs = plt.subplots(2, 3, figsize=(14, 8), layout='constrained')
        for j, (title, ids) in enumerate(groups):
            style(axs[0, j], 'S, град⁻¹', 25); lines(axs[0, j], ids, 'S_h1_per_degree', low=1, high=25)
            axs[0, j].set_title(title); axs[0, j].axhline(0, color='.55', lw=.6)
            axs[0, j].legend(frameon=False, fontsize=8)
            style(axs[1, j], 'S, град⁻¹'); lines(axs[1, j], ids, 'S_h1_per_degree', low=10, high=89)
            axs[1, j].set_xlim(10, 90)
            axs[1, j].set_yscale('symlog', linthresh=cfg['display_only_scales']['S_symlog_linthresh'])
            axs[1, j].axhline(0, color='.55', lw=.6)
        save(fig, 'signed_sensitivity')
        fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), layout='constrained', sharey=True)
        for ax, (title, ids) in zip(axs, groups):
            style(ax, '|S|, град⁻¹'); lines(ax, ids, 'abs_S_h1', low=1, high=89)
            ax.set_yscale('log'); ax.set_title(title); ax.legend(frameon=False, fontsize=8)
        save(fig, 'absolute_sensitivity')
        fig, axs = plt.subplots(2, 2, figsize=(11, 7.8), layout='constrained', sharex=True, sharey=True)
        for i, label in enumerate(cfg['illustration_cases']):
            part = [r for r in fits if r['ID'] == label and int(r['window_deg']) == 90]
            for j, (method, field, caption) in enumerate((('chord', 'chord_F', 'Хорда'), ('segmented', 'segmented_F', 'Два участка'))):
                ax = axs[i, j]; style(ax, 'F')
                ax.plot([float(r['beta_deg']) for r in part], [float(r['F']) for r in part], color=colors[label], lw=1.7, label='Сохранённые частоты')
                ax.plot([float(r['beta_deg']) for r in part], [float(r[field]) for r in part], color='.35', ls='--', label=caption)
                t = next(r for r in transitions if r['ID'] == label and int(r['window_deg']) == 90 and r['method'] == method)
                b = int(t['beta_deg']); ax.axvline(b, color='.25', ls=':', lw=1)
                ax.scatter([b], [float(t['F'])], color=colors[label], s=24, zorder=4)
                ax.set_title(names[label].capitalize()+f'\nХарактерный переход: {b}°')
                ax.legend(frameon=False, fontsize=8, loc='lower right')
        save(fig, 'transition_methods')
    write(folder/'plot_manifest.json', {'mode': 'plot_only', 'new_eigenvalue_solves': 0,
        'regression_fits': 0, 'PNG_dpi': 300, 'PDF': 'vector',
        'images': {f: sha(output/f) for f in images},
        'note': 'Symlog widths are display settings, never plateau criteria; R(90)=0 is a singleton identity. S figures show centered nodes 1..89 only; one-sided endpoint estimates remain in CSV.'})


def saturation_main(args):
    cfg = json.loads(SATURATION_CONFIG.read_text(encoding='utf8'))
    tables = saturation_sources(cfg)
    ident = {'config': cfg, 'files': {Path(__file__).relative_to(ROOT).as_posix(): sha(Path(__file__)),
        SATURATION_CONFIG.relative_to(ROOT).as_posix(): sha(SATURATION_CONFIG)},
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    folder = OUTPUT/'saturation_diagnostics'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = (args.verify_cache or folder).resolve()
        m = json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if m['identity'] != ident or any(sha(target/f) != h for f, h in m['artifacts'].items()):
            raise ValueError('Saturation cache mismatch')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True); write(folder/'contract.json', ident)
    diagnostics, transitions, profiles, curves, comparison, derivative = [], [], [], [], [], []
    for label, table in tables.items():
        rows = saturation_measures(label, table); diagnostics.extend(rows)
        current = []
        for window in cfg['windows_deg']:
            ts, ps, cs = saturation_transitions(rows, window)
            current.extend(ts); profiles.extend(ps); curves.extend(cs)
        transitions.extend(current)
        full = {t['method']: t for t in current if t['window_deg'] == 90}
        comparison.append({'ID': label,
            **{f'{field}_{method}': full[method][field] for method in ('chord', 'segmented')
               for field in ('beta_deg', 'F', 'Lambda1', 'omega1', 'R_to_90', 'S_h1_per_degree')},
            'angle_difference_chord_minus_segmented_deg': full['chord']['beta_deg']-full['segmented']['beta_deg'],
            **{f'R_at_{b}': rows[b]['R'] for b in cfg['control_angles_deg']}})
        changes = {f'h{h}': [[a['beta_deg'], b['beta_deg']] for a, b in zip(rows[:-1], rows[1:])
                  if a[f'S_h{h}_per_degree'] is not None and b[f'S_h{h}_per_degree'] is not None
                  and a[f'S_h{h}_per_degree']*b[f'S_h{h}_per_degree'] < 0] for h in (1, 2, 3)}
        derivative.append({'ID': label, 'sign_change_intervals_deg': changes,
            'grid_maximum_omega_beta': max(rows, key=lambda r: r['omega1'])['beta_deg'],
            'F_at_90': rows[-1]['F'],
            'max_abs_difference_h1_h3_first_25': max(abs(r['S_h1_per_degree']-r['S_h3_per_degree']) for r in rows if 3 <= r['beta_deg'] <= 25),
            'max_abs_difference_h1_h3_after_10': max(abs(r['S_h1_per_degree']-r['S_h3_per_degree']) for r in rows if 10 <= r['beta_deg'] <= 87)})
    for name, rows in [('diagnostics', diagnostics), ('transitions', transitions), ('fit_objectives', profiles),
                       ('fit_curves', curves), ('comparison', comparison)]:
        csvwrite(folder/(name+'.csv'), rows)
    write(folder/'derivative_diagnostics.json', derivative)
    plot_saturation(folder)
    saturation_sources(cfg)  # Verify preservation again, without numerical reevaluation.
    write(folder/'manifest.json', {'status': 'PASS', 'scope': 'Descriptive saved-spectrum diagnostics; no accepted plateau angle',
        'identity': ident, 'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'new_eigenvalue_solves': 0, 'matrix_assemblies': 0, 'shape_recoveries': 0,
        'source_cases': 10, 'source_angle_points': 910, 'regression_fits': len(profiles),
        'artifacts': {p.relative_to(folder).as_posix(): sha(p) for p in folder.rglob('*') if p.is_file() and p.name != 'manifest.json'}})
    print(folder, flush=True)


def plateau_tails(label, table):
    """Positive saved omega data; no modal reconstruction or curve fitting."""
    if [r['beta_deg'] for r in table] != list(range(91)) or any(r['omega1'] <= 0 for r in table):
        raise ValueError('Plateau analysis requires the accepted positive 0:1:90 map')
    tails = tail_ranges([{'ID': label, **r} for r in table])
    rows = [{**point, **tail, 'R_percent': 100*tail['R']} for point, tail in zip(table, tails)]
    if any(r['R'] < 0 for r in rows) or any(b['R'] > a['R'] for a, b in zip(rows, rows[1:])) or rows[-1]['R'] != 0:
        raise ValueError('Nested-tail range monotonicity failed')
    return rows


def plateau_at(rows, delta_percent, minimum_length):
    """First qualifying integer node, or an explicit missing value, never a clamp."""
    if delta_percent <= 0 or not 0 < minimum_length <= 90:
        raise ValueError('Positive tolerance and a nonzero admissible tail required')
    # Compare percentages directly: saved event thresholds are exactly included.
    hit = next(r for r in rows if r['R_percent'] <= delta_percent)
    j = int(hit['beta_deg']); previous = rows[j-1] if j else None
    eligible = j <= 90-minimum_length
    if previous and previous['R_percent'] <= delta_percent:
        raise ValueError('First-hit condition failed')
    cutoff = rows[int(90-minimum_length)]
    return {'ID': rows[0]['ID'], 'delta_percent': float(delta_percent), 'L_min_deg': minimum_length,
        'status': 'QUALIFYING_PLATEAU' if eligible else 'NO_QUALIFYING_PLATEAU',
        'beta_pl_deg': j if eligible else None,
        'first_tolerance_hit_deg': j, 'cutoff_beta_deg': cutoff['beta_deg'],
        'existence_threshold_percent': cutoff['R_percent'],
        'remaining_length_deg': 90-j if eligible else None,
        'Lambda_pl': hit['Lambda1'] if eligible else None,
        'Omega_pl': hit['Omega1'] if eligible else None,
        'omega_pl': hit['omega1'] if eligible else None,
        'q_pl': hit['omega1']/max(r['omega1'] for r in rows) if eligible else None,
        'R_pl_percent': hit['R_percent'] if eligible else None,
        'previous_beta_deg': j-1 if eligible and previous else None,
        'previous_R_percent': previous['R_percent'] if eligible and previous else None,
        'tail_beta_min': hit['beta_min'] if eligible else None,
        'tail_beta_max': hit['beta_max'] if eligible else None}


def plateau_order_intervals(tables, pairs, lower, upper, minimum_length):
    """Exact order intervals of discrete inverses, without ranking missing values."""
    result = []
    for a, b in pairs:
        events = sorted({lower, upper} | {r['R_percent'] for label in (a, b) for r in tables[label]
            if r['beta_deg'] <= 90-minimum_length and lower < r['R_percent'] < upper})
        intervals = []
        for i, d in enumerate(events):
            pa = plateau_at(tables[a], d, minimum_length)['beta_pl_deg']
            pb = plateau_at(tables[b], d, minimum_length)['beta_pl_deg']
            relation = ('BOTH_UNAVAILABLE' if pa is None and pb is None else
                'A_UNAVAILABLE' if pa is None else 'B_UNAVAILABLE' if pb is None else
                'EARLIER' if pa < pb else 'LATER' if pa > pb else 'TIED')
            end = events[i+1] if i+1 < len(events) else upper
            if intervals and intervals[-1]['relation_A_to_B'] == relation:
                intervals[-1]['delta_upper_percent'] = end
                intervals[-1]['upper_inclusive'] = i == len(events)-1
            else:
                intervals.append({'A': a, 'B': b, 'L_min_deg': minimum_length,
                    'delta_lower_percent': d, 'delta_upper_percent': end,
                    'lower_inclusive': True, 'upper_inclusive': i == len(events)-1,
                    'relation_A_to_B': relation})
        result.extend(intervals)
    return result


def plot_plateau(folder):
    """Only render saved exact threshold events; never evaluate a spectrum."""
    folder = Path(folder).resolve()
    if (folder/'manifest.json').exists():
        m = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f) != h for f, h in m['artifacts'].items()):
            raise ValueError('Plateau cache corrupted')
    rows = list(csv.DictReader((folder/'plot_steps.csv').open(encoding='utf8')))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = {'RR': 'жёсткий узел — заделка', 'HR': 'шарнир — заделка',
             'RH': 'жёсткий узел — шарнир', 'HH': 'шарнир — шарнир'}
    colors = dict(zip(names, ('#245a81', '#d37e14', '#32916c', '#a84e9b')))
    for prefix in ('J', 'E'):
        for k, color in zip((.1, 1, 10), ('#245a81', '#d37e14', '#32916c')):
            label = f'J{k:g}_ER' if prefix == 'J' else f'JR_E{k:g}'
            names[label] = rf'$\kappa_{{joint}}={k:g}$' if prefix == 'J' else rf'$\kappa_{{end}}={k:g}$'
            colors[label] = color
    groups = [('limits', ['RR','HR','RH','HH']),
              ('internal', ['HR','J0.1_ER','J1_ER','J10_ER','RR']),
              ('external', ['RH','JR_E0.1','JR_E1','JR_E10','RR'])]
    out = folder/'figures'; out.mkdir(exist_ok=True); images = {}
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42}):
        for name, labels in groups:
            fig, ax = plt.subplots(figsize=(8, 5.2), layout='constrained')
            for label in labels:
                part = [r for r in rows if r['ID'] == label]
                x = [float(r['delta_percent']) for r in part]
                y = [float(r['beta_pl_deg']) if r['beta_pl_deg'] else np.nan for r in part]
                ls = '--' if name != 'limits' and label == 'RR' else ':' if name != 'limits' and label in ('HR','RH') else '-'
                ax.step(x, y, where='post', lw=1.65, color=colors[label], ls=ls, label=names[label])
            ax.set(xscale='log', xlim=(.01,5), ylim=(0,63), xlabel='δ, %',
                   ylabel=r'$\beta_{\mathrm{pl}}^{\mathrm{disc}}$, °', title='Минимальная протяжённость 30°')
            ticks = [.01,.02,.05,.1,.2,.5,1,2,5]
            ax.set_xticks(ticks, [f'{v:g}' for v in ticks]); ax.set_yticks(range(0,61,10))
            ax.grid(alpha=.18); ax.spines[['top','right']].set_visible(False)
            ax.legend(frameon=False, fontsize=9, loc='upper right')
            for ext in ('png','pdf'):
                p = out/f'beta_plateau_{name}.{ext}'
                fig.savefig(p, dpi=300, metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
                images[p.name] = sha(p)
            plt.close(fig)
    write(folder/'plot_manifest.json', {'new_eigenvalue_solves':0, 'postprocessing_evaluations':0,
        'mode':'plot_only', 'PNG_dpi':300, 'PDF':'vector', 'images':images,
        'semantics':'Right-continuous exact steps of the saved discrete criterion; missing values break curves'})


def plateau_main(args):
    cfg = json.loads(PLATEAU_CONFIG.read_text(encoding='utf8'))
    source_tables = saturation_sources(cfg)
    old = ROOT/cfg['previous_diagnostics']['bundle']
    old_manifest = json.loads((old/'manifest.json').read_text(encoding='utf8'))
    if sha(old/'manifest.json') != cfg['previous_diagnostics']['manifest_sha256'] or any(
            sha(old/p) != h for p, h in old_manifest['artifacts'].items()):
        raise ValueError('Previous range diagnostics provenance failed')
    ident = {'config':cfg, 'files':{Path(__file__).relative_to(ROOT).as_posix():sha(__file__),
        PLATEAU_CONFIG.relative_to(ROOT).as_posix():sha(PLATEAU_CONFIG)},
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}}
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    folder = OUTPUT/'plateau_criterion_step_A'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = (args.verify_cache or folder).resolve()
        m = json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if m['identity'] != ident or any(sha(target/p) != h for p,h in m['artifacts'].items()):
            raise ValueError('Plateau cache mismatch')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    started = time.perf_counter(); folder.mkdir(parents=True, exist_ok=True)
    write(folder/'contract.json', ident)
    tables = {label:plateau_tails(label, table) for label,table in source_tables.items()}
    previous = {(r['ID'],float(r['beta_deg'])):float(r['R']) for r in csv.DictReader((old/'diagnostics.csv').open(encoding='utf8'))}
    all_tails = [r for table in tables.values() for r in table]
    if any(r['R'] != previous[r['ID'],r['beta_deg']] for r in all_tails):
        raise ValueError('Range differs from previous diagnostics')
    main = []; existence = []; steps = []
    for label, table in tables.items():
        for length in cfg['L_min_deg']:
            current = [plateau_at(table,d,length) for d in cfg['delta_percent']]
            main.extend(current)
            values = [r['beta_pl_deg'] for r in current if r['beta_pl_deg'] is not None]
            if any(b>a for a,b in zip(values,values[1:])):
                raise ValueError('Tolerance inverse monotonicity failed')
            existence.append({'ID':label,'L_min_deg':length, 'cutoff_beta_deg':90-length,
                'minimum_delta_percent':table[90-length]['R_percent']})
        events = sorted(set(cfg['delta_percent']) | {r['R_percent'] for r in table
            if r['beta_deg'] <= 60 and min(cfg['delta_percent']) <= r['R_percent'] <= max(cfg['delta_percent'])})
        steps.extend(plateau_at(table,d,30) for d in events)
    statuses = []
    for label in tables:
        for d in cfg['delta_percent']:
            selected = [r for r in main if r['ID']==label and r['delta_percent']==d]
            valid = [r['beta_pl_deg'] for r in selected if r['beta_pl_deg'] is not None]
            if len(set(valid))>1 or any(a['beta_pl_deg'] is None and b['beta_pl_deg'] is not None for a,b in zip(selected,selected[1:])):
                raise ValueError('Minimum-length restriction changed first-hit angle')
            statuses.append({'ID':label,'delta_percent':d, **{f'status_L{r["L_min_deg"]}':r['status'] for r in selected},
                             **{f'beta_L{r["L_min_deg"]}':r['beta_pl_deg'] for r in selected}})
    orders = plateau_order_intervals(tables,cfg['comparison_pairs'],min(cfg['delta_percent']),max(cfg['delta_percent']),30)
    for name, rows in [('tail_ranges',all_tails),('plateau_criteria',main),('existence_thresholds',existence),
        ('existence_status',statuses),('plot_steps',steps),('comparison_intervals',orders),
        ('control_ranges',[r for r in all_tails if r['beta_deg'] in (10,20,40)])]:
        csvwrite(folder/(name+'.csv'), rows)
    arithmetic_seconds = time.perf_counter()-started
    plot_start = time.perf_counter(); plot_plateau(folder); plot_seconds = time.perf_counter()-plot_start
    saturation_sources(cfg)
    write(folder/'diagnostics.json', {'status':'PASS','range_replay_max_abs_error':0.,
        'positive_frequency_nested_range_checks':910,'criterion_combinations':len(main),
        'qualifying_combinations':sum(r['beta_pl_deg'] is not None for r in main),
        'no_qualifying_combinations':sum(r['beta_pl_deg'] is None for r in main),
        'monotonicity_and_length_restriction':'PASS','comparison_intervals':len(orders),
        'normalization':'Omega=200*omega=Lambda**2; source-reader relative check 2e-14',
        'discrete_bracket':'previous_R > delta >= R_at_first_hit, no certified continuous enclosure'})
    write(folder/'manifest.json', {'status':'PASS','scope':'Discrete criterion sensitivity, no final engineering tolerance chosen',
        'identity':ident,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'new_eigenvalue_solves':0,'spectral_matrix_assemblies':0,'shape_recoveries':0,
        'reused_spectral_points':910,'postprocessing':{'tail_ranges':910,'requested_criteria':len(main),'plot_events':len(steps),
            'comparison_intervals':len(orders),'arithmetic_seconds':arithmetic_seconds,'render_seconds':plot_seconds},
        'artifacts':{p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}})
    print(folder, flush=True)


def grid_identity(cfg):
    """Fixed nine-pair extension; no change to the accepted arm or spring physics."""
    _, ident = finite_identity(cfg)
    for prefix in ('finite', 'step_A'):
        folder = ROOT/cfg[prefix+'_bundle']
        old = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if sha(folder/'manifest.json') != cfg[prefix+'_manifest_sha256'] or old['status'] != 'PASS':
            raise ValueError('Grid predecessor manifest mismatch')
        if any(sha(folder/f) != h for f, h in old['artifacts'].items()):
            raise ValueError('Grid predecessor artifact mismatch')
        if prefix == 'finite':
            for field in ('material', 'mu', 'policy', 'map'):
                if cfg[field] != old['identity']['config'][field]:
                    raise ValueError('Grid changed accepted physical/numerical contract')
        ident[prefix+'_manifest_sha256'] = sha(folder/'manifest.json')
    ident['files'][GRID_CONFIG.relative_to(ROOT).as_posix()] = sha(GRID_CONFIG)
    return hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16], ident


def grid_order(rows, boundaries, cfg):
    """Min-max tests for sorted frequencies, not for the plateau angle."""
    lookup = {(r['ID'], float(r['beta_deg'])): r for r in boundaries+rows}
    labels = {tuple(v): k for k, v in cfg['cases'].items()}; grid = cfg['stiffness_grid']
    pairs = []
    for fixed in grid:
        for a, b in zip(grid[:-1], grid[1:]):
            pairs += [('joint', labels[a, fixed], labels[b, fixed]),
                      ('end', labels[fixed, a], labels[fixed, b])]
    for label, (joint, end) in cfg['cases'].items():
        pairs += [('HH', 'HH', label), ('RR', label, 'RR'),
                  ('joint_RIGID', label, f'JR_E{end:g}'),
                  ('end_RIGID', label, f'J{joint:g}_ER')]
    result = []
    for beta in cfg['map']['parameter_grid']:
        for kind, low, high in pairs:
            for i in (1, 2, 3):
                margin = float(lookup[high, beta][f'omega{i}'])/float(lookup[low, beta][f'omega{i}'])-1
                if margin < -cfg['comparison_relative_tolerance']:
                    raise ArithmeticError(f'Grid min-max failed: {low}, {high}, beta={beta}, k={i}')
                result.append(dict(kind=kind, beta_deg=beta, position=i, lower=low,
                                   upper=high, relative_margin=margin))
    return result


def grid_summary(run, rows, checks, neighbours):
    cfg, folder = run.config, run.folder
    boundaries = sum((list(csv.DictReader((ROOT/cfg[name]/'spectrum.csv').open(encoding='utf8')))
                      for name in ('limit_bundle', 'finite_bundle')), [])
    order = grid_order(rows, boundaries, cfg); csvwrite(folder/'minmax.csv', order)
    tails, plateaus, sensitivity, curves = [], [], [], []
    rr = {float(r['beta_deg']): float(r['omega1']) for r in boundaries if r['ID'] == 'RR'}
    for label, (joint, end) in cfg['cases'].items():
        part = [r for r in rows if r['ID'] == label]
        if len(part) != 91:
            raise ValueError('Incomplete grid case')
        for r in part:
            for i in (1, 2, 3):
                if not (math.isclose(r[f'omega{i}']*run.model.time, r[f'Omega{i}'], rel_tol=2e-14)
                        and math.isclose(r[f'Lambda{i}']**2, r[f'Omega{i}'], rel_tol=2e-14)):
                    raise ValueError('Grid normalization failed')
        tail = plateau_tails(label, part); tails.extend(tail)
        local = []
        for length in cfg['minimum_length_deg']:
            angles = []
            for delta in cfg['delta_percent']:
                p = plateau_at(tail, delta, length)
                j = p['beta_pl_deg']
                p.update(kappa_joint=joint, kappa_end=end,
                    omega_relative_to_RR_at_plateau=None if j is None else p['omega_pl']/rr[j],
                    tolerance_margin_percent=None if j is None else delta-p['R_pl_percent'],
                    previous_excess_percent=None if j is None or not j else p['previous_R_percent']-delta)
                local.append(p)
                if j is not None:
                    angles.append(j)
            if any(b > a for a, b in zip(angles, angles[1:])):
                raise ValueError('Plateau tolerance order failed')
        for delta in cfg['delta_percent']:
            a, b = [p for p in local if p['delta_percent'] == delta]
            if b['beta_pl_deg'] is not None and b['beta_pl_deg'] != a['beta_pl_deg']:
                raise ValueError('Minimum length incorrectly moved first hit')
        plateaus.extend(local)
        chosen = [p for p in local if p['L_min_deg'] == cfg['main_minimum_length_deg']]
        angles = [p['beta_pl_deg'] for p in chosen]
        sensitivity.append(dict(ID=label, kappa_joint=joint, kappa_end=end,
            beta_01=angles[0], beta_02=angles[1], beta_05=angles[2],
            change_01_to_02=None if None in angles[:2] else angles[1]-angles[0],
            change_02_to_05=None if None in angles[1:] else angles[2]-angles[1]))
        extrema = [b['beta_deg'] for a, b, c in zip(part[:-2], part[1:-1], part[2:])
                   if (b['omega1']-a['omega1'])*(c['omega1']-b['omega1']) < 0]
        closest = min(part, key=lambda r: r['omega2']/r['omega1']-1)
        curves.append(dict(ID=label, kappa_joint=joint, kappa_end=end,
            interior_extrema_beta=extrema, maximum_beta=max(part, key=lambda r: r['omega1'])['beta_deg'],
            Lambda1_0=part[0]['Lambda1'], Lambda1_90=part[-1]['Lambda1'],
            minimum_sampled_relative_gap=closest['omega2']/closest['omega1']-1,
            minimum_sampled_gap_beta=closest['beta_deg']))
    csvwrite(folder/'tail_ranges.csv', tails); csvwrite(folder/'plateaus.csv', plateaus)
    csvwrite(folder/'criterion_sensitivity.csv', sensitivity)
    comparisons = []
    for i, left in enumerate(cfg['cases']):
        for right in list(cfg['cases'])[i+1:]:
            for delta in cfg['delta_percent']:
                pair = [next(p for p in plateaus if p['ID'] == label and p['delta_percent'] == delta
                              and p['L_min_deg'] == cfg['main_minimum_length_deg']) for label in (left, right)]
                a, b = [p['beta_pl_deg'] for p in pair]
                comparisons.append(dict(left=left, right=right, delta_percent=delta,
                    left_beta=a, right_beta=b, difference=None if a is None or b is None else a-b,
                    relation='UNDEFINED' if a is None or b is None else 'EARLIER' if a < b else 'LATER' if a > b else 'EQUAL'))
    csvwrite(folder/'pair_order.csv', comparisons)
    control_plateaus = []
    for label in sorted({r['ID'] for r in boundaries}):
        part = [{k: float(r[k]) for k in ('beta_deg', 'omega1', 'Omega1', 'Lambda1')}
                for r in boundaries if r['ID'] == label]
        tail = plateau_tails(label, part)
        control_plateaus += [plateau_at(tail, d, L) for L in cfg['minimum_length_deg'] for d in cfg['delta_percent']]
    csvwrite(folder/'boundary_plateaus.csv', control_plateaus)
    summary = dict(status='PASS', points=len(rows), grid_pairs=9,
        source_points=sum(r['provenance'] == 'verification_bundle' for r in rows),
        physical_maxima={k: max(r[k] for r in checks) for k in ('physical', 'joint_absolute',
            'external_absolute', 'ode', 'energy', 'singular', 'condition', 'mass_gram', 'mass_error')},
        minmax_tests=len(order), minmax_minimum_margin=min(r['relative_margin'] for r in order),
        neighbour_triggers=[r for r in neighbours if r['trigger']], curves=curves,
        plateau_records=len(plateaus), missing_plateaus=[p for p in plateaus if p['beta_pl_deg'] is None],
        qualification='Nine actual stiffness pairs and 1-degree beta grid; sorted positions; discrete tail criterion, no continuous surface certification.')
    write(folder/'summary.json', summary)
    return summary


def plot_grid(folder):
    """Only saved scalar CSV; no spectral assembly, reconstruction or plateau solve."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    folder = Path(folder)
    if (folder/'manifest.json').exists():
        manifest = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        for name in ('plateaus.csv', 'spectrum.csv', 'contract.json'):
            if sha(folder/name) != manifest['artifacts'][name]:
                raise ValueError('Grid plot source changed')
    cfg = json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    data = list(csv.DictReader((folder/'plateaus.csv').open(encoding='utf8')))
    data = [r for r in data if float(r['L_min_deg']) == cfg['main_minimum_length_deg']]
    spec = list(csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')))
    figures = folder/'figures'; figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
    values = [float(r['beta_pl_deg']) for r in data if r['beta_pl_deg']]
    norm = Normalize(min(values) if values else 0, max(values) if values else 60)
    cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('#eeeeee')
    def save(fig, name):
        fig.savefig(figures/(name+'.png'), dpi=300, bbox_inches='tight')
        fig.savefig(figures/(name+'.pdf'), bbox_inches='tight', metadata={'CreationDate': None, 'ModDate': None})
        plt.close(fig)
    matrices = []
    for delta in cfg['delta_percent']:
        z = np.full((3, 3), np.nan)
        for r in data:
            if float(r['delta_percent']) == delta and r['beta_pl_deg']:
                x = cfg['stiffness_grid'].index(float(r['kappa_joint']))
                y = cfg['stiffness_grid'].index(float(r['kappa_end']))
                z[y, x] = float(r['beta_pl_deg'])
        matrices.append(z)
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
        fig.subplots_adjust(left=.14, right=.88, bottom=.16, top=.87)
        art = ax.pcolormesh([-1.5, -.5, .5, 1.5], [-1.5, -.5, .5, 1.5], z,
                           norm=norm, cmap=cmap, edgecolors='white', linewidth=2)
        for y in range(3):
            for x in range(3):
                value = z[y, x]
                ax.text(x-1, y-1, '—' if np.isnan(value) else f'{value:.0f}°',
                        ha='center', va='center', fontsize=15,
                        color='black' if np.isnan(value) or norm(value) > .55 else 'white')
        ax.set(xticks=[-1, 0, 1], xticklabels=['0.1', '1', '10'],
               yticks=[-1, 0, 1], yticklabels=['0.1', '1', '10'],
               xlabel=r'$\kappa_{\rm joint}$', ylabel=r'$\kappa_{\rm end}$',
               title=rf'$\delta={delta:g}\%$, $L_{{\min}}=30^\circ$')
        fig.colorbar(art, ax=ax, label=r'$\beta_{\rm pl}^{\rm disc}$, °')
        save(fig, 'plateau_'+str(delta).replace('.', 'p'))
    fig = plt.figure(figsize=(12, 4.2))
    fig.subplots_adjust(left=.02, right=.96, bottom=.16, top=.88, wspace=.16)
    x, y = np.meshgrid([-1, 0, 1], [-1, 0, 1])
    for i, (delta, z) in enumerate(zip(cfg['delta_percent'], matrices), 1):
        ax = fig.add_subplot(1, 3, i, projection='3d')
        ax.scatter(x, y, z, c=z.ravel(), cmap=cmap, norm=norm, s=40, depthshade=False)
        ax.plot_wireframe(x, y, z, color='0.5', linewidth=.6)
        ax.set(xticks=[-1, 0, 1], xticklabels=['0.1', '1', '10'],
               yticks=[-1, 0, 1], yticklabels=['0.1', '1', '10'],
               xlabel=r'$\kappa_{\rm joint}$', ylabel=r'$\kappa_{\rm end}$',
               zlabel=r'$\beta_{\rm pl}$, °', zlim=(0, 60), title=rf'$\delta={delta:g}\%$')
        ax.view_init(elev=24, azim=-55)
    save(fig, 'plateau_nodes_3d')
    for xmax, name in ((90, 'lambda12_grid'), (20, 'lambda12_grid_small_angles')):
        fig, axes = plt.subplots(3, 3, figsize=(10, 8), sharex=True, sharey=True)
        fig.subplots_adjust(left=.07, right=.99, bottom=.08, top=.94, wspace=.12, hspace=.28)
        for label, (joint, end) in cfg['cases'].items():
            ax = axes[cfg['stiffness_grid'].index(end), cfg['stiffness_grid'].index(joint)]
            part = sorted([r for r in spec if r['ID'] == label and float(r['beta_deg']) <= xmax],
                          key=lambda r: float(r['beta_deg']))
            for k, color, style in ((1, '#1f77b4', '-'), (2, '#d62728', '--')):
                ax.plot([float(r['beta_deg']) for r in part], [float(r[f'Lambda{k}']) for r in part],
                        color=color, linestyle=style, linewidth=1.4, label=rf'$\Lambda_{k}$')
            ax.set(title=rf'$\kappa_{{\rm joint}}={joint:g},\ \kappa_{{\rm end}}={end:g}$',
                   xlim=(0, xmax), xticks=[0, 30, 60, 90] if xmax == 90 else [0, 5, 10, 15, 20])
            ax.grid(alpha=.2)
        for ax in axes[-1]: ax.set_xlabel('β, °')
        for ax in axes[:, 0]: ax.set_ylabel('Λ')
        axes[0, 0].legend(loc='lower right')
        save(fig, name)


def grid_main(args):
    cfg = json.loads(GRID_CONFIG.read_text()); key, ident = grid_identity(cfg)
    folder = OUTPUT/'plateau_grid_step_B'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = (args.verify_cache or folder).resolve()
        old = json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if old['identity'] != ident or any(sha(target/f) != h for f, h in old['artifacts'].items()):
            raise ValueError('Grid cache mismatch; no overwrite')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True)
    write(folder/'contract.json', ident)
    run = Run(cfg, folder); seed_limit_points(run); start = time.perf_counter()
    logpath = folder/'run_ledger.json'; ledger = json.loads(logpath.read_text()) if logpath.exists() else []
    try:
        preflight = folder/'preflight.json'
        if not preflight.exists() or json.loads(preflight.read_text())['status'] != 'PASS':
            limit_preflight(run, cfg['preflight_labels'])
        ritzpath = folder/'ritz_checks.json'
        ritz = json.loads(ritzpath.read_text()) if ritzpath.exists() else []
        for label, beta in cfg['ritz_points']:
            if any(r['ID'] == label and r['beta'] == beta for r in ritz): continue
            record = run.solve(physics.Case(beta, cfg['mu'], *cfg['cases'][label]))
            check = run.ritz_check(record)
            if any(v['constraint_residual'] > cfg['ritz_constraint_tol'] for v in check['orders']):
                raise ArithmeticError('Ritz kinematic constraint failed')
            ritz.append(dict(ID=label, beta=beta, **check)); write(ritzpath, ritz)
            print(f'Ritz {label} beta={beta}: PASS', flush=True)
        if args.stage != 'preflight':
            summary = limit_map(run, finite=True, summarize=grid_summary)
            plot_grid(folder)
    except Exception as exc:
        write(folder/'failure.json', dict(error=repr(exc), new_points=run.new_points, new_roots=run.new_roots))
        raise
    finally:
        ledger.append(dict(command=sys.argv, new_points=run.new_points, new_roots=run.new_roots,
            reused_point_calls=run.reused_points, ritz_solves=run.ritz_solves,
            pole_roots=run.pole_roots, seconds=time.perf_counter()-start))
        write(logpath, ledger)
    if args.stage != 'preflight':
        _, after = grid_identity(cfg)
        if after != ident: raise ValueError('Sources changed during grid calculation')
        manifest = dict(status=summary['status'], identity=ident,
            git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            git_status=subprocess.check_output(['git', 'status', '--short'], text=True),
            points=summary['points'], source_points=summary['source_points'], complex_solves=0,
            **{k: sum(r[k] for r in ledger) for k in ('new_points', 'new_roots', 'ritz_solves', 'pole_roots')},
            artifacts={p.relative_to(folder).as_posix(): sha(p) for p in folder.rglob('*')
                       if p.is_file() and p.name != 'manifest.json'})
        write(folder/'manifest.json', manifest)
    print(folder, flush=True)


def early_sources(cfg):
    """Scalar saved-data reader; never assembles a boundary matrix or reads modes."""
    source = ROOT/cfg['source_bundle']
    manifest = json.loads((source/'manifest.json').read_text(encoding='utf8'))
    if sha(source/'manifest.json') != cfg['source_manifest_sha256'] or manifest['status'] != 'PASS':
        raise ValueError('Early-analysis source manifest mismatch')
    if any(sha(source/f) != h for f, h in manifest['artifacts'].items()):
        raise ValueError('Early-analysis source artifact mismatch')
    old = manifest['identity']['config']
    for key in ('material', 'mu', 'cases'):
        if old[key] != cfg[key]: raise ValueError('Early-analysis physical contract mismatch')
    if old['map']['parameter_grid'] != cfg['beta_deg'] or old['map']['spectrum_semantics'] != 'sorted_positions':
        raise ValueError('Source grid or semantics changed')
    tables = {label: [] for label in cfg['cases']}
    for row in csv.DictReader((source/'spectrum.csv').open(encoding='utf8')):
        label = row['ID']; beta = int(row['beta_deg'])
        point = json.loads((source/'points'/row['point_file']).read_text(encoding='utf8'))
        if point['status'] != 'PASS' or point['count']['status'] != 'PASS' or len(point['roots']) != 3:
            raise ValueError('Source prefix not accepted')
        if point['case'] != dict(beta=beta, mu=cfg['mu'], k_joint=cfg['cases'][label][0],
                                 k_end=cfg['cases'][label][1], variant=cfg['variant']):
            raise ValueError('Source point physical provenance failed')
        for i, r in enumerate(point['roots'], 1):
            if not all(r['gates'].values()) or any(float(row[f'{f}{i}']) != r[f] for f in ('omega', 'Omega', 'Lambda')):
                raise ValueError('CSV and accepted roots disagree')
            tol = cfg['numerical_checks']['normalization_relative']
            if not (math.isclose(r['omega']*200, r['Omega'], rel_tol=tol)
                    and math.isclose(r['Lambda']**2, r['Omega'], rel_tol=tol)):
                raise ValueError('Source frequency normalization failed')
        tables[label].append({k: float(row[k]) for k in ('beta_deg', 'omega1', 'Omega1', 'Lambda1')})
    for part in tables.values():
        part.sort(key=lambda r: r['beta_deg'])
        if [r['beta_deg'] for r in part] != cfg['beta_deg']:
            raise ValueError('Missing or duplicate source nodes; no interpolation allowed')
    return tables


def early_measures(label, table, intervals):
    """Signed finite-interval log changes and full-tail ranges, not derivatives."""
    tails = plateau_tails(label, table)
    w0, w90 = table[0]['omega1'], table[-1]['omega1']
    maximum = max(r['omega1'] for r in table)
    for r in tails:
        w = r['omega1']; attained = w-w0
        r.update(G90_percent=100*(w90-w)/w,
                 gain_from_0_percent=100*(w/w0-1),
                 F=None if maximum == w0 else attained/(maximum-w0),
                 remaining_span_over_attained_gain_percent=None if attained == 0 else
                    100*(r['omega_max']-r['omega_min'])/abs(attained),
                 omega_change_from_15_percent=100*(w/table[15]['omega1']-1))
        if abs(r['G90_percent']) > r['R_percent']+2e-13:
            raise ValueError('Endpoint change exceeds enclosing tail range')
    averages = []
    for a, b in intervals:
        x, y = table[a], table[b]
        value = 100*math.log(y['omega1']/x['omega1'])/(b-a)
        averages.append(dict(ID=label, beta_from=a, beta_to=b,
            mean_log_percent_per_degree=value,
            signed_interval_change_percent=100*(y['omega1']/x['omega1']-1),
            omega_from=x['omega1'], omega_to=y['omega1']))
    reference = abs(averages[0]['mean_log_percent_per_degree'])
    for r in averages:
        value = abs(r['mean_log_percent_per_degree'])
        r.update(magnitude_relative_to_initial=None if reference == 0 else value/reference,
                 reduction_factor=None if value == 0 else reference/value)
    return tails, averages


def plot_early(folder):
    """Four saved-data figures; no criterion estimation or spectral functions."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    folder = Path(folder)
    names = ('contract.json', 'profiles.csv', 'interval_sensitivity.csv', 'comparison.csv')
    if (folder/'manifest.json').exists():
        m = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f) != m['artifacts'][f] for f in names):
            raise ValueError('Early plot inputs changed')
    cfg = json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    read = lambda name: list(csv.DictReader((folder/(name+'.csv')).open(encoding='utf8')))
    profile, intervals, comp = read('profiles'), read('interval_sensitivity'), read('comparison')
    figures = folder/'figures'; figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'pdf.fonttype': 42,
                        'axes.spines.top': False, 'axes.spines.right': False})
    def save(fig, name):
        fig.savefig(figures/(name+'.png'), dpi=300, bbox_inches='tight')
        fig.savefig(figures/(name+'.pdf'), bbox_inches='tight', metadata={'CreationDate': None, 'ModDate': None})
        plt.close(fig)
    number = lambda rows, name: np.array([float(r[name]) for r in rows])
    labels = list(cfg['cases'])
    # Shared main scale; insets resolve the weak change using a common percent scale.
    early = [r for r in profile if float(r['beta_deg']) <= 25]
    zoom = [r for r in early if float(r['beta_deg']) >= 10]
    yl = (min(number(early, 'Lambda1'))-.08, max(number(early, 'Lambda1'))+.12)
    zl = (min(number(zoom, 'omega_change_from_15_percent'))-.05,
          max(number(zoom, 'omega_change_from_15_percent'))+.05)
    fig, axes = plt.subplots(3, 3, figsize=(11, 9), sharex=True, sharey=True)
    fig.subplots_adjust(left=.07, right=.98, bottom=.07, top=.94, wspace=.15, hspace=.3)
    for label in labels:
        j, e = cfg['cases'][label]
        ax = axes[[.1, 1., 10.].index(e), [.1, 1., 10.].index(j)]
        part = [r for r in early if r['ID'] == label]
        ax.plot(number(part, 'beta_deg'), number(part, 'Lambda1'), color='#1f77b4', lw=1.4)
        ax.set(title=rf'{label}: $\kappa_{{\rm joint}}={j:g},\ \kappa_{{\rm end}}={e:g}$',
               xlim=(0, 25), ylim=yl, xticks=[0, 5, 10, 15, 20, 25])
        ax.grid(alpha=.18)
        inset = ax.inset_axes([.44, .13, .52, .45])
        small = [r for r in part if float(r['beta_deg']) >= 10]
        inset.plot(number(small, 'beta_deg'), number(small, 'omega_change_from_15_percent'),
                   color='#1f77b4', marker='.', markersize=2, lw=.9)
        inset.axhline(0, color='0.6', lw=.5)
        inset.set(xlim=(10, 25), ylim=zl, xticks=[10, 15, 25], title=r'$100(\omega/\omega_{15}-1),\ \%$')
        inset.title.set_fontsize(7); inset.tick_params(labelsize=7); inset.grid(alpha=.15)
    for ax in axes[-1]: ax.set_xlabel('β, °')
    for ax in axes[:, 0]: ax.set_ylabel('Λ₁')
    save(fig, '01_lambda1_early')
    # Signed common symlog color scales retain possible zeros and negative values.
    values = np.array([[float(r['mean_log_percent_per_degree']) for r in intervals if r['ID'] == label] for label in labels])
    ratios = np.array([[100*float(r['magnitude_relative_to_initial']) for r in intervals if r['ID'] == label] for label in labels])
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.8))
    fig.subplots_adjust(left=.15, right=.96, bottom=.12, top=.9, wspace=.4)
    for n, (ax, matrix, title, colorlabel) in enumerate(zip(axes, (values, ratios),
            ('Средняя логарифмическая чувствительность', 'Модуль относительно интервала 0–5°'),
            ('%/град', '% исходного модуля'))):
        vmax = max(abs(matrix.min()), abs(matrix.max()))
        norm = SymLogNorm(linthresh=.01, vmin=-vmax if matrix.min() < 0 else 0, vmax=vmax)
        image = ax.pcolormesh(np.arange(6)-.5, np.arange(10)-.5, matrix,
            cmap='coolwarm' if matrix.min() < 0 else 'viridis', norm=norm, rasterized=False)
        ax.set_xlim(-.5, 4.5); ax.set_ylim(8.5, -.5)
        ax.set(xticks=range(5), xticklabels=[f'{a}–{b}°' for a,b in cfg['intervals_deg']],
               yticks=range(9), yticklabels=[f'{k} ({cfg["cases"][k][0]:g}; {cfg["cases"][k][1]:g})' for k in labels] if n == 0 else labels,
               title=title)
        for y in range(9):
            for x in range(5):
                value = matrix[y, x]
                ax.text(x, y, f'{value:.3g}', ha='center', va='center', fontsize=8,
                        color='black' if norm(value) > .55 else 'white')
        colorbar = fig.colorbar(image, ax=ax, shrink=.8, label=colorlabel)
        colorbar.solids.set_rasterized(False)
    save(fig, '02_interval_sensitivity')
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    fig.subplots_adjust(left=.08, right=.99, bottom=.2, top=.94)
    x = np.arange(9)
    for i, (beta, color) in enumerate(zip(cfg['control_angles_deg'], ('#1f77b4', '#ff7f0e', '#2ca02c'))):
        ax.bar(x+(i-1)*.24, [float(r[f'R_{beta}_percent']) for r in comp], width=.23, color=color, label=f'{beta}°')
    ax.set(xticks=x, xticklabels=[f'{k}\n({cfg["cases"][k][0]:g}; {cfg["cases"][k][1]:g})' for k in labels],
           ylabel='Остаточный размах R, %', xlabel='Случай (κ_joint; κ_end)')
    ax.grid(axis='y', alpha=.2); ax.legend(title='Начало хвоста', ncols=3)
    save(fig, '03_remaining_range')
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    fig.subplots_adjust(left=.075, right=.99, bottom=.08, top=.92, wspace=.3, hspace=.45)
    for col, label in enumerate(cfg['example_ids']):
        part = [r for r in profile if r['ID'] == label]
        row = next(r for r in comp if r['ID'] == label)
        tol = float(row['beta_tol_05_deg']); ax, lower = axes[:, col]
        early = [r for r in part if float(r['beta_deg']) <= 30]
        ax.plot(number(early, 'beta_deg'), number(early, 'Lambda1'), color='#1f77b4', lw=1.4)
        ax.axvspan(*cfg['illustrative_transition_window_deg'], color='0.85', alpha=.55)
        for a in (ax, lower):
            a.axvline(15, color='0.4', ls=':', lw=1)
            a.axvline(tol, color='#d62728', ls='--', lw=1)
            a.grid(alpha=.18)
        j,e = cfg['cases'][label]
        ax.set(xlim=(0,30), xlabel='β, °', ylabel='Λ₁',
               title=rf'{label}: $\kappa_{{\rm joint}}={j:g},\ \kappa_{{\rm end}}={e:g}$')
        lower.plot(number(part, 'beta_deg'), number(part, 'R_percent'), color='#1f77b4', lw=1.4)
        lower.axhline(.5, color='#d62728', ls='--', lw=.8)
        lower.scatter([15], [float(row['R_15_percent'])], color='#ff7f0e', zorder=3)
        lower.set(xlim=(5,90), ylim=(0, max(3.,float(row['R_10_percent'])*1.15)),
                  xlabel='β, °', ylabel='R, %', title=f'R(15°)={float(row["R_15_percent"]):.3f}%; β_tol={tol:.0f}°')
    save(fig, '04_transition_vs_tolerance')


def early_main(args):
    cfg = json.loads(EARLY_CONFIG.read_text(encoding='utf8')); tables = early_sources(cfg)
    ident = dict(config=cfg, files={Path(__file__).relative_to(ROOT).as_posix(): sha(__file__),
        EARLY_CONFIG.relative_to(ROOT).as_posix(): sha(EARLY_CONFIG)},
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    folder = OUTPUT/'early_saturation_analysis'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target = (args.verify_cache or folder).resolve(); m = json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if m['identity'] != ident or any(sha(target/f) != h for f,h in m['artifacts'].items()):
            raise ValueError('Early-analysis cache mismatch')
        print(f'REUSED {target}: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True); write(folder/'contract.json', ident)
    source = ROOT/cfg['source_bundle']
    old_ranges = {(r['ID'],float(r['beta_deg'])):float(r['R']) for r in csv.DictReader((source/'tail_ranges.csv').open(encoding='utf8'))}
    old_angles = [r for r in csv.DictReader((source/'plateaus.csv').open(encoding='utf8'))
                  if float(r['delta_percent']) == cfg['displayed_tolerance_percent'] and float(r['L_min_deg']) == cfg['displayed_minimum_length_deg']]
    profiles, averages, controls, comp = [], [], [], []
    started = time.perf_counter()
    for label, table in tables.items():
        tail, intervals = early_measures(label, table, cfg['intervals_deg'])
        if any(r['R'] != old_ranges[label,r['beta_deg']] for r in tail):
            raise ValueError('Tail range differs from Step B')
        for r in intervals:
            a,b = r['beta_from'],r['beta_to']; value = r['mean_log_percent_per_degree']
            alternate = [100*math.log(table[b]['Omega1']/table[a]['Omega1'])/(b-a),
                         200*math.log(table[b]['Lambda1']/table[a]['Lambda1'])/(b-a)]
            if max(abs(v-value) for v in alternate) > cfg['numerical_checks']['log_identity_absolute']:
                raise ValueError('Log-change normalization failed')
        total = sum(r['mean_log_percent_per_degree']*(r['beta_to']-r['beta_from']) for r in intervals)
        if abs(total-100*math.log(table[30]['omega1']/table[0]['omega1'])) > cfg['numerical_checks']['log_identity_absolute']:
            raise ValueError('Interval log changes do not telescope')
        accepted = next(r for r in old_angles if r['ID'] == label)
        if accepted['status'] != 'QUALIFYING_PLATEAU': raise ValueError('Missing saved tolerance angle')
        beta_tol = int(accepted['beta_pl_deg'])
        first = next(r['beta_deg'] for r in tail if r['R_percent'] <= .5 and r['beta_deg'] <= 60)
        if first != beta_tol: raise ValueError('Saved tolerance angle replay failed')
        j,e = cfg['cases'][label]
        profiles.extend(tail); averages.extend(intervals)
        controls.extend([r for r in tail if r['beta_deg'] in cfg['control_angles_deg']])
        comp.append(dict(ID=label,kappa_joint=j,kappa_end=e,
            **{f'S_{r["beta_from"]}_{r["beta_to"]}_percent_per_degree':r['mean_log_percent_per_degree'] for r in intervals},
            **{f'reduction_{r["beta_from"]}_{r["beta_to"]}':r['reduction_factor'] for r in intervals},
            **{f'R_{b}_percent':tail[b]['R_percent'] for b in cfg['control_angles_deg']},
            **{f'G90_{b}_percent':tail[b]['G90_percent'] for b in cfg['control_angles_deg']},
            F15=tail[15]['F'],gain_0_15_percent=tail[15]['gain_from_0_percent'],
            remaining_span_over_gain_0_15_percent=tail[15]['remaining_span_over_attained_gain_percent'],
            beta_tol_05_deg=beta_tol,beta_sat_deg=None,
            full_grid_maximum_beta=max(tail,key=lambda r:r['omega1'])['beta_deg']))
    for name, rows in [('profiles',profiles),('interval_sensitivity',averages),('control_ranges',controls),('comparison',comp)]:
        csvwrite(folder/(name+'.csv'),rows)
    arithmetic_seconds=time.perf_counter()-started
    plot_start=time.perf_counter();plot_early(folder);plot_seconds=time.perf_counter()-plot_start
    early_sources(cfg)
    write(folder/'diagnostics.json',dict(status='PASS',range_replay_max_abs_error=0,
        interval_normalization_and_telescoping='PASS',source_prefixes=819,
        saved_tolerance_angles_reproduced=9,beta_sat_assigned=False))
    write(folder/'manifest.json',dict(status='PASS',scope='Saved-spectrum early-growth diagnostics; no beta_sat estimator',
        identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        new_eigenvalue_solves=0,spectral_matrix_assemblies=0,shape_recoveries=0,
        reused_spectral_points=819,postprocessing=dict(intervals=len(averages),ranges=len(profiles),
            requested_tail_controls=len(controls),saved_tolerance_controls=9,
            arithmetic_seconds=arithmetic_seconds,render_seconds=plot_seconds),
        artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def growth_slopes(log_values, h):
    """Degree-grid differences; h=1/2, second-order one-sided boundary layers."""
    y = np.asarray(log_values, dtype=float)
    if y.shape != (91,) or h not in (1, 2) or not np.all(np.isfinite(y)):
        raise ValueError('Expected finite 0:1:90 log frequencies and h=1 or 2')
    out = np.empty(91)
    for i in range(91):
        if i < h:
            out[i] = (-3*y[i]+4*y[i+h]-y[i+2*h])/(2*h)
        elif i > 90-h:
            out[i] = (3*y[i]-4*y[i-h]+y[i-2*h])/(2*h)
        else:
            out[i] = (y[i+h]-y[i-h])/(2*h)
    return out


def growth_pick(q, level, minimum_length=30):
    """Persistent discrete slope threshold; no endpoint substitution or slack."""
    q = np.asarray(q, dtype=float)
    if q.shape != (91,) or np.any(q < 0) or not np.all(np.isfinite(q)) or level <= 0:
        raise ValueError('Invalid normalized sensitivity')
    if not 0 < minimum_length <= 90:
        raise ValueError('Invalid remaining interval')
    envelope = np.maximum.accumulate(q[::-1])[::-1]
    hits = np.flatnonzero((envelope <= level) & (np.arange(91) <= 90-minimum_length))
    first = int(hits[0]) if len(hits) else None
    result = dict(status='QUALIFYING_TRANSITION' if first is not None else 'NO_QUALIFYING_TRANSITION',
        beta_pl_deg=first, q_star=level, L_min_deg=minimum_length,
        Q_at_transition=None if first is None else float(envelope[first]),
        previous_Q=None if first in (None, 0) else float(envelope[first-1]),
        previous_beta_deg=None if first in (None, 0) else first-1,
        tail_peak_beta=None if first is None else first+int(np.argmax(q[first:])))
    if first is not None and (np.any(q[first:] > level) or np.any(envelope[:first] <= level)):
        raise ArithmeticError('First persistent crossing failed')
    return result, envelope


def plot_growth(folder):
    """Six figures from derived CSVs only; vector meshes, no refitting."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    folder = Path(folder); names = ('contract.json', 'profiles.csv', 'transitions.csv')
    if (folder/'manifest.json').exists():
        m = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f) != m['artifacts'][f] for f in names):
            raise ValueError('Growth plot input hashes changed')
    cfg = json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    read = lambda f: list(csv.DictReader((folder/f).open(encoding='utf8')))
    profiles, transitions = read('profiles.csv'), read('transitions.csv')
    primary = [r for r in transitions if int(r['h_deg']) == 1 and int(r['reference_end_deg']) == 5]
    figures = folder/'figures'; figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    def save(fig, stem):
        fig.savefig(figures/(stem+'.png'), dpi=300, bbox_inches='tight')
        fig.savefig(figures/(stem+'.pdf'), bbox_inches='tight', metadata={'CreationDate':None,'ModDate':None})
        plt.close(fig)
    def part(label, h=1):
        return [r for r in profiles if r['ID']==label and int(r['h_deg'])==h and int(r['reference_end_deg'])==5]
    arr = lambda rows, key: np.array([float(r[key]) for r in rows])
    colors = {0.20:'#2ca02c',0.15:'#ff7f0e',0.10:'#d62728'}
    accepted = [int(r['beta_pl_deg']) for r in primary if r['beta_pl_deg']]
    lo, hi = (min(accepted), max(accepted)) if accepted else (0,60)
    edges = np.arange(lo-.5,hi+1.5); norm=BoundaryNorm(edges,256)
    for level in cfg['q_levels']:
        fig,ax=plt.subplots(figsize=(5.2,4.2));fig.subplots_adjust(left=.15,right=.88,bottom=.14,top=.90)
        matrix=np.full((3,3),np.nan)
        for r in primary:
            if float(r['q_star'])!=level:continue
            j,e=cfg['cases'][r['ID']]; y=[.1,1,10].index(e);x=[.1,1,10].index(j)
            if r['beta_pl_deg']:matrix[y,x]=int(r['beta_pl_deg'])
        cmap=plt.colormaps['viridis'].copy();cmap.set_bad('0.85')
        mesh=ax.pcolormesh(np.arange(4)-1.5,np.arange(4)-1.5,np.ma.masked_invalid(matrix),cmap=cmap,norm=norm,rasterized=False)
        for y in range(3):
            for x in range(3):
                value=matrix[y,x];txt='—' if np.isnan(value) else f'{value:.0f}°'
                ax.text(x-1,y-1,txt,ha='center',va='center',fontsize=14,
                    color='black' if np.isnan(value) or value >= (lo+hi)/2 else 'white')
        ax.set(xticks=[-1,0,1],xticklabels=['0.1','1','10'],yticks=[-1,0,1],yticklabels=['0.1','1','10'],
            xlabel='κ_joint',ylabel='κ_end',title=rf'$q_*={level:.2f}$; $h=1°$, опора 0–5°')
        cb=fig.colorbar(mesh,ax=ax,ticks=range(lo,hi+1),label='β_pl, °');cb.solids.set_rasterized(False)
        save(fig,f'01_map_q{level:.2f}')
    fig,axes=plt.subplots(3,3,figsize=(11,8.5),sharex=True,sharey=True)
    fig.subplots_adjust(left=.07,right=.98,bottom=.08,top=.92,wspace=.15,hspace=.28)
    for label,(j,e) in cfg['cases'].items():
        ax=axes[[.1,1,10].index(e),[.1,1,10].index(j)]
        for h,color,ls in [(1,'#1f77b4','-'),(2,'#9467bd','--')]:
            rows=[r for r in part(label,h) if float(r['beta_deg'])<=25]
            ax.plot(arr(rows,'beta_deg'),arr(rows,'q'),color=color,ls=ls,lw=1.2,label=f'h={h}°')
        for level,color in colors.items():ax.axhline(level,color=color,lw=.8,ls=':')
        ax.set(yscale='symlog',xlim=(0,25),title=f'{label}: ({j:g}; {e:g})',xticks=[0,5,10,15,20,25]);ax.set_yscale('symlog',linthresh=.01)
        ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=8);fig.text(.5,.97,'q(β): уровни 0.20 / 0.15 / 0.10',ha='center')
    for ax in axes[-1]:ax.set_xlabel('β, °')
    for ax in axes[:,0]:ax.set_ylabel('q')
    save(fig,'02_q_early')
    fig,axes=plt.subplots(2,3,figsize=(11.5,7));fig.subplots_adjust(left=.08,right=.99,bottom=.08,top=.92,hspace=.42,wspace=.3)
    for col,label in enumerate(cfg['example_ids']):
        rows=part(label);wide=part(label,2)
        for ax in axes[:,col]:
            ax.plot(arr(rows,'beta_deg'),arr(rows,'q'),color='0.65',lw=1,label='q, h=1°')
            ax.plot(arr(rows,'beta_deg'),arr(rows,'Q_tail'),color='#1f77b4',lw=1.5,label='Q_tail, h=1°')
            ax.plot(arr(wide,'beta_deg'),arr(wide,'Q_tail'),color='#9467bd',ls='--',lw=1,label='Q_tail, h=2°')
            ax.grid(alpha=.18);ax.set_xlabel('β, °');ax.set_ylabel('q, Q_tail')
        upper,lower=axes[:,col];upper.set(xlim=(0,25),ylim=(0,.6),title=label)
        late=[r for r in rows+wide if float(r['beta_deg'])>=25]
        lower.set(xlim=(25,90),ylim=(0,1.15*max(arr(late,'Q_tail'))))
        lower.set_yscale('symlog',linthresh=1e-6)
        for r in primary:
            if r['ID']!=label:continue
            lev=float(r['q_star']);upper.axhline(lev,color=colors[lev],ls=':',lw=.8)
            if r['beta_pl_deg']:upper.scatter([int(r['beta_pl_deg'])],[float(r['Q_at_transition'])],color=colors[lev],s=25,zorder=5)
    axes[0,0].legend(fontsize=7);save(fig,'03_tail_envelope')
    fig,axes=plt.subplots(1,3,figsize=(11.5,3.7));fig.subplots_adjust(left=.07,right=.99,bottom=.17,top=.87,wspace=.28)
    for ax,label in zip(axes,cfg['example_ids']):
        rows=[r for r in part(label) if float(r['beta_deg'])<=25]
        ax.plot(arr(rows,'beta_deg'),arr(rows,'Lambda1'),color='#1f77b4',lw=1.4)
        for r in primary:
            if r['ID']==label and r['beta_pl_deg']:
                lev=float(r['q_star']);ax.axvline(int(r['beta_pl_deg']),color=colors[lev],ls='--',lw=1,
                    label=f'q*={lev:.2f}: {r["beta_pl_deg"]}°')
        ax.set(xlim=(0,25),xlabel='β, °',ylabel='Λ₁',title=label);ax.grid(alpha=.2);ax.legend(fontsize=7)
    save(fig,'04_frequency_transitions')


def growth_main(args):
    cfg=json.loads(GROWTH_CONFIG.read_text(encoding='utf8'));tables=early_sources(cfg)
    geo=ROOT/cfg['geometry_source_bundle'];gm=json.loads((geo/'manifest.json').read_text(encoding='utf8'))
    if sha(geo/'manifest.json')!=cfg['geometry_manifest_sha256'] or any(sha(geo/f)!=v for f,v in gm['artifacts'].items()):
        raise ValueError('Geometry diagnostic provenance mismatch')
    ident=dict(config=cfg,files={Path(__file__).relative_to(ROOT).as_posix():sha(__file__),
        GROWTH_CONFIG.relative_to(ROOT).as_posix():sha(GROWTH_CONFIG)},
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    folder=OUTPUT/'growth_transition_sensitivity'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target=(args.verify_cache or folder).resolve();m=json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if m['identity']!=ident or any(sha(target/f)!=v for f,v in m['artifacts'].items()):
            raise ValueError('Growth cache identity or artifacts mismatch')
        print(f'REUSED {target}: 0 eigenvalue solves');return
    started=time.perf_counter();folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident)
    profiles=[];transitions=[];references=[];q_error=0.;derivative_error=0.
    for label,table in tables.items():
        tails=plateau_tails(label,table);y=np.log([r['omega1'] for r in table]);z=np.log([r['Lambda1'] for r in table])
        for ref in cfg['reference_intervals_end_deg']:
            s0=float((y[ref]-y[0])/ref);l0=float((z[ref]-z[0])/ref)
            if s0<=0 or l0<=0:raise ValueError('Nonpositive initial reference slope')
            references.append(dict(ID=label,reference_end_deg=ref,S0_per_degree=s0,Lambda_reference_slope=l0))
            for h in cfg['derivative_halfwidth_deg']:
                signed=growth_slopes(y,h);alternate=growth_slopes(z,h);q=abs(signed)/s0;ql=abs(alternate)/l0
                q_error=max(q_error,float(np.max(abs(q-ql))));derivative_error=max(derivative_error,float(np.max(abs(signed-2*alternate))))
                if q_error>cfg['numerical_checks']['q_identity_absolute'] or derivative_error>cfg['numerical_checks']['derivative_identity_absolute']:
                    raise ArithmeticError('Log normalization identity failed')
                last=None
                for level in sorted(cfg['q_levels'],reverse=True):
                    result,envelope=growth_pick(q,level,cfg['minimum_length_deg']);beta=result['beta_pl_deg']
                    if last is not None and beta is not None and beta<last:raise ArithmeticError('Threshold order failed')
                    last=beta
                    first_local=next((i for i in range(5,91) if q[i]<=level),None)
                    transitions.append(dict(ID=label,h_deg=h,reference_end_deg=ref,**result,
                        Lambda1=None if beta is None else table[beta]['Lambda1'],
                        R_percent=None if beta is None else tails[beta]['R_percent'],
                        q_at_transition=None if beta is None else float(q[beta]),
                        first_local_below_after_5_deg=first_local,
                        persistent_delay_deg=None if beta is None or first_local is None else beta-first_local,
                        late_q_max_after_25=float(max(q[25:])),late_q_peak_beta=25+int(np.argmax(q[25:]))))
                if np.any(np.diff(envelope)>0):raise ArithmeticError('Tail envelope must not increase')
                for i,r in enumerate(table):
                    profiles.append(dict(ID=label,h_deg=h,reference_end_deg=ref,**r,
                        signed_S_per_degree=float(signed[i]),S0_per_degree=s0,q=float(q[i]),
                        Q_tail=float(envelope[i]),Q_peak_beta=i+int(np.argmax(q[i:])),R_percent=tails[i]['R_percent'],
                        stencil='centered' if h<=i<=90-h else 'one_sided_second_order'))
    lookup={(r['ID'],r['q_star'],r['h_deg'],r['reference_end_deg']):r for r in transitions}
    sensitivity=[]
    for label in cfg['cases']:
        for level in cfg['q_levels']:
            one=lookup[label,level,1,5];two=lookup[label,level,2,5]
            b1,b2=one['beta_pl_deg'],two['beta_pl_deg'];shift=None if None in (b1,b2) else b2-b1
            sensitivity.append(dict(ID=label,q_star=level,beta_h1=b1,beta_h2=b2,h2_minus_h1_deg=shift,
                beta_ref4=lookup[label,level,1,4]['beta_pl_deg'],beta_ref6=lookup[label,level,1,6]['beta_pl_deg'],
                qualification='SCHEME_SENSITIVE' if shift is None or abs(shift)>1 else 'WITHIN_ONE_GRID_STEP'))
    contrast=[dict(ID=label,beta_q020=lookup[label,.20,1,5]['beta_pl_deg'],beta_q010=lookup[label,.10,1,5]['beta_pl_deg'],
        delta_deg=None if any(lookup[label,q,1,5]['beta_pl_deg'] is None for q in (.20,.10)) else
            lookup[label,.10,1,5]['beta_pl_deg']-lookup[label,.20,1,5]['beta_pl_deg']) for label in cfg['cases']]
    order=[];labels=list(cfg['cases'])
    def relation(a,b):return 'UNAVAILABLE' if None in (a,b) else ('EARLIER' if a<b else 'LATER' if a>b else 'TIED')
    for h in cfg['derivative_halfwidth_deg']:
        for ref in cfg['reference_intervals_end_deg']:
            for level in cfg['q_levels']:
                for i,a in enumerate(labels):
                    for b in labels[i+1:]:
                        order.append(dict(ID_A=a,ID_B=b,h_deg=h,reference_end_deg=ref,q_star=level,
                            relation=relation(lookup[a,level,h,ref]['beta_pl_deg'],lookup[b,level,h,ref]['beta_pl_deg']),
                            baseline_same_q=relation(lookup[a,level,1,5]['beta_pl_deg'],lookup[b,level,1,5]['beta_pl_deg']),
                            baseline_q020=relation(lookup[a,.20,1,5]['beta_pl_deg'],lookup[b,.20,1,5]['beta_pl_deg'])))
    # Recover previous stiffness contracts, never identify cases by frequency proximity.
    prior={}
    for src in gm['identity']['config']['sources']:
        sm=json.loads((ROOT/src['bundle']/'manifest.json').read_text(encoding='utf8'))
        if sha(ROOT/src['bundle']/'manifest.json')!=src['manifest_sha256']:raise ValueError('Prior geometry source changed')
        prior.update(sm['identity']['config']['cases'])
    comparison=[dict(ID=label,status='NO_COMPARABLE_STIFFNESS_PAIR',matching_old_IDs=';'.join(k for k,v in prior.items() if v==pair))
        for label,pair in cfg['cases'].items()]
    if any(r['matching_old_IDs'] for r in comparison):raise ValueError('Comparable geometry case requires explicit saved-fit comparison')
    for name,rows in [('profiles',profiles),('transitions',transitions),('reference_slopes',references),
            ('sensitivity',sensitivity),('criterion_difference',contrast),('pair_order',order),('geometry_comparison',comparison)]:
        csvwrite(folder/(name+'.csv'),rows)
    seconds=time.perf_counter()-started;started=time.perf_counter();plot_growth(folder);render=time.perf_counter()-started
    early_sources(cfg)
    write(folder/'diagnostics.json',dict(status='PASS',q_normalization_max_abs_error=q_error,
        derivative_normalization_max_abs_error=derivative_error,accepted_transitions=sum(r['beta_pl_deg'] is not None for r in transitions),
        missing_transitions=sum(r['beta_pl_deg'] is None for r in transitions),
        scheme_sensitive_primary_cases=sum(r['qualification']=='SCHEME_SENSITIVE' for r in sensitivity),
        exact_geometry_overlap=0))
    write(folder/'manifest.json',dict(status='PASS',scope='Discrete normalized-slope criterion; no final q_star selected',identity=ident,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),new_eigenvalue_solves=0,
        spectral_matrix_assemblies=0,shape_recoveries=0,reused_spectral_points=819,
        postprocessing=dict(profiles=len(profiles),transitions=len(transitions),arithmetic_seconds=seconds,render_seconds=render),
        artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def refined_slopes(beta, log_values, archived_slopes, halfwidth):
    """Explicit stencils on the two grids; never differentiate the mixed spacing."""
    values = dict(zip(beta, log_values)); out = []; stencils = []; comparable = []
    if halfwidth not in (.25, .5) or len(values) != 118:
        raise ValueError('Expected 91 archived plus 27 fractional angles')
    for b in beta:
        both = 3.5 <= b <= 11.5
        if 3 < b < 12:
            h = halfwidth if halfwidth == .25 or both else .25
            if b-h not in values or b+h not in values:
                raise ValueError('Missing explicit derivative endpoints')
            out.append((values[b+h]-values[b-h])/(2*h))
            stencils.append(f'centered_{h:g}' if h == halfwidth else 'primary_edge_fallback_0.25')
        else:
            out.append(archived_slopes[int(b)]); stencils.append('archived_h1')
        comparable.append(both)
    return np.array(out), stencils, comparable


def refined_pick(beta, q, level, length=30):
    """First persistent crossing on saved nodes, with an explicit predecessor."""
    beta = np.asarray(beta); q = np.asarray(q)
    if beta.shape != q.shape or np.any(np.diff(beta) <= 0) or not np.all(np.isfinite(q)) or np.any(q < 0):
        raise ValueError('Invalid refined sensitivity grid')
    envelope = np.maximum.accumulate(q[::-1])[::-1]
    valid = np.flatnonzero((beta <= 90-length) & (envelope <= level))
    i = int(valid[0]) if len(valid) else None
    r = dict(status='PASS' if i is not None else 'NO_QUALIFYING_TRANSITION',
        beta_pl_deg=None if i is None else float(beta[i]),
        previous_beta_deg=None if i in (None, 0) else float(beta[i-1]),
        q_at_transition=None if i is None else float(q[i]),
        previous_q=None if i in (None, 0) else float(q[i-1]),
        Q_at_transition=None if i is None else float(envelope[i]),
        previous_Q=None if i in (None, 0) else float(envelope[i-1]),
        tail_peak_beta=None if i is None else float(beta[i+np.argmax(q[i:])]))
    if np.any(np.diff(envelope) > 0) or (i is not None and
            (np.any(q[i:] > level) or np.any(envelope[:i] <= level))):
        raise ArithmeticError('Persistent crossing or envelope failed')
    return r, envelope


def refinement_identity(cfg):
    tables = early_sources(cfg)
    _, ident = identity(cfg)
    ident['files'][REFINE_CONFIG.relative_to(ROOT).as_posix()] = sha(REFINE_CONFIG)
    manifests = {}
    for field in ('source', 'old_growth', 'catalog'):
        source = ROOT/cfg[field+'_bundle']
        m = json.loads((source/'manifest.json').read_text(encoding='utf8'))
        if sha(source/'manifest.json') != cfg[field+'_manifest_sha256'] or m['status'] != 'PASS':
            raise ValueError('Refinement source manifest mismatch')
        if any(sha(source/f) != h for f, h in m['artifacts'].items()):
            raise ValueError('Refinement source artifact mismatch')
        manifests[field] = m
    source = manifests['source']['identity']
    for field in ('material', 'mu', 'cases', 'policy'):
        if cfg[field] != source['config'][field]:
            raise ValueError('Refinement physics or quality gates changed')
    for name, h in source['files'].items():
        if name != Path(__file__).relative_to(ROOT).as_posix() and sha(ROOT/name) != h:
            raise ValueError('Verified spectral implementation changed: '+name)
    if cfg['map']['parameter_grid'] != [3+i*.25 for i in range(37)]:
        raise ValueError('Refinement grid outside authorized scope')
    old = manifests['old_growth']['identity']['config']
    if any(old[f] != cfg[f] for f in ('material', 'mu', 'cases', 'source_manifest_sha256')):
        raise ValueError('Old derivative contract mismatch')
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    return key, ident, tables


def refinement_spectra(run, preflight_only=False):
    """Bounded prefix calls; integers and pole catalogs are strictly read-only reuse."""
    cfg = run.config; source = ROOT/cfg['source_bundle']
    saved = {(r['ID'], float(r['beta_deg'])): r for r in csv.DictReader((source/'spectrum.csv').open())}
    catalogs = ROOT/cfg['catalog_bundle']/'pole_catalogs.json'
    run.catalogs = {float(k): v for k, v in json.loads(catalogs.read_text()).items()}
    audits_path = run.folder/'new_count_checks.json'
    audits = json.loads(audits_path.read_text()) if audits_path.exists() else {}
    def point(label, b, previous=None):
        if b not in cfg['map']['parameter_grid'] or label not in cfg['cases']:
            raise ValueError('Unauthorized spectral point')
        if b.is_integer():
            path = source/'points'/saved[label, b]['point_file']
            return json.loads(path.read_text()), path, 'archive'
        case = physics.Case(b, cfg['mu'], *cfg['cases'][label])
        casekey = hashlib.sha256(json.dumps({**asdict(case), 'variant':'RL_FULL'}, sort_keys=True).encode()).hexdigest()[:16]
        path = run.folder/'points'/(casekey+'.json')
        if not path.exists() and len(list((run.folder/'points').glob('*.json'))) >= cfg['maximum_new_points']:
            raise ArithmeticError('New spectral point budget exceeded')
        if previous is None:
            previous = np.array([float(saved[label, math.floor(b)][f'omega{i}']) for i in (1,2,3)])
        rec = run.solve(case, neighbour=previous)
        if rec['status'] != 'PASS' or rec['count']['status'] != 'PASS' or not all(all(r['gates'].values()) for r in rec['roots']):
            raise ArithmeticError('Cached or new prefix is not accepted')
        if casekey not in audits:
            w = np.array([r['omega'] for r in rec['roots']]); counts = []
            for query, expected in zip([float(run.model.omega(run.policy['Lambda_lower'])), *(.5*(w[:-1]+w[1:]))], (0,1,2)):
                _, n, d = physics.count_at(run.model, case, float(query), run.cats(run.model, case), run.policy)
                if n != expected: raise ArithmeticError('Refinement low-prefix count failed')
                counts.append(dict(query=float(query), count=n, expected=expected, **d))
            audits[casekey] = dict(ID=label, beta_deg=b, status='PASS', queries=counts)
            write(audits_path, audits)
        return rec, path, 'new_fractional'
    preflight = []
    for label in cfg['preflight_labels']:
        for b in cfg['preflight_angles']:
            rec, path, origin = point(label, b)
            w = np.array([r['omega'] for r in rec['roots']])
            jumps = [float(max(abs(w/np.array([float(saved[label, side][f'omega{i}']) for i in (1,2,3)])-1)))
                     for side in (math.floor(b), math.ceil(b))]
            if max(jumps) > run.policy['neighbour_relative_jump']:
                raise ArithmeticError('Preflight archive-neighbour jump needs local review')
            preflight.append(dict(ID=label, beta_deg=b, Lambda=[r['Lambda'] for r in rec['roots']],
                archive_neighbour_relative_jump=max(jumps), status='PASS'))
    write(run.folder/'preflight.json', dict(status='PASS', points=preflight))
    print('Refinement preflight: 12 accepted fractional points', flush=True)
    if preflight_only: return None
    rows = []; checks = []; neighbours = []; reuse = []
    for label in cfg['cases']:
        previous = None
        for b in cfg['map']['parameter_grid']:
            rec, path, origin = point(label, b, previous)
            w = np.array([r['omega'] for r in rec['roots']])
            if np.any(np.diff(w) <= 0): raise ArithmeticError('Unordered prefix')
            if previous is not None:
                jump = float(max(abs(w/previous-1)))
                neighbours.append(dict(ID=label, beta_deg=b, relative_jump=jump,
                    trigger=jump > run.policy['neighbour_relative_jump'], gap12_relative=float(w[1]/w[0]-1),
                    resolution='count-isolated prefix, independent midgap counts and physical gates'))
                if jump > run.policy['neighbour_relative_jump']:
                    raise ArithmeticError('Refined neighbour trigger requires local review')
            previous = w
            row = dict(ID=label, beta_deg=b)
            for i, r in enumerate(rec['roots'], 1):
                row.update({f'{f}{i}': r[f] for f in ('omega','Omega','Lambda')})
                if not (math.isclose(r['Omega'],200*r['omega'],rel_tol=2e-14) and math.isclose(r['Omega'],r['Lambda']**2,rel_tol=2e-14)):
                    raise ArithmeticError('New frequency normalization failed')
                d = r['diagnostics']
                checks.append(dict(ID=label,beta_deg=b,position=i,provenance=origin,
                    physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),
                    ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],
                    condition=d['nonzero_condition'],mass_gram=rec['mass_gram_error'],mass_error=abs(d['mass']-1),status=rec['status']))
            row.update(provenance=origin, point_file=path.relative_to(ROOT).as_posix(), point_sha256=sha(path))
            rows.append(row)
            reuse.append(dict(ID=label,beta_deg=b,origin=origin,point_file=row['point_file'],sha256=sha(path),shape_sha256=rec['shape_sha256']))
        print(f'Refined {label}: 37 accepted nodes (10 archived, 27 fractional)',flush=True)
    if len(rows)!=333 or sum(r['origin']=='archive' for r in reuse)!=90 or len(audits)!=243:
        raise ArithmeticError('Unexpected local-grid accounting')
    # Min-max for each independent spring, including the guard position.
    order=[]; by={(r['ID'],r['beta_deg']):r for r in rows}
    for a, pair in cfg['cases'].items():
        for b, upper in cfg['cases'].items():
            if sum(x==y for x,y in zip(pair,upper))!=1 or not all(x<=y for x,y in zip(pair,upper)):continue
            for angle in cfg['map']['parameter_grid']:
                for i in (1,2,3):
                    margin=by[b,angle][f'omega{i}']/by[a,angle][f'omega{i}']-1
                    if margin < -cfg['comparison_relative_tolerance']:raise ArithmeticError('Refined frequency min-max failed')
                    order.append(dict(lower=a,upper=b,beta_deg=angle,position=i,relative_margin=margin))
    for name, data in [('spectrum',rows),('checks',checks),('node_provenance',reuse),('minmax',order),('neighbour_audit',neighbours)]:
        csvwrite(run.folder/(name+'.csv'),data)
    return rows,checks,neighbours


def refinement_analysis(folder, cfg, tables, spectra):
    oldfolder=ROOT/cfg['old_growth_bundle']
    read=lambda p:list(csv.DictReader(p.open(encoding='utf8')))
    oldprofiles=read(oldfolder/'profiles.csv');oldtrans=read(oldfolder/'transitions.csv')
    oldangles={(r['ID'],float(r['q_star'])):float(r['beta_pl_deg']) for r in oldtrans
               if float(r['h_deg'])==1 and int(r['reference_end_deg'])==5 and r['beta_pl_deg']}
    profiles=[];transitions=[];references=[];normalization_error=0.;derivative_error=0.;stitch=[]
    for label, archived in tables.items():
        merged={r['beta_deg']:dict(r,ID=label) for r in archived}
        merged.update({r['beta_deg']:r for r in spectra if r['ID']==label})
        table=[merged[b] for b in sorted(merged)];beta=np.array(sorted(merged));tails=tail_ranges(table)
        y=np.log([r['omega1'] for r in table]);z=np.log([r['Lambda1'] for r in table])
        oldpart=[r for r in oldprofiles if r['ID']==label and float(r['h_deg'])==1 and int(r['reference_end_deg'])==5]
        archived_d=np.array([float(r['signed_S_per_degree']) for r in oldpart])
        dy=growth_slopes(np.log([r['omega1'] for r in archived]),1)
        if not np.array_equal(dy,archived_d):raise ArithmeticError('Archived derivative was not reproduced')
        dz=growth_slopes(np.log([r['Lambda1'] for r in archived]),1)
        for ref in cfg['reference_intervals_end_deg']:
            s0=float((y[np.where(beta==ref)[0][0]]-y[0])/ref)
            l0=float((z[np.where(beta==ref)[0][0]]-z[0])/ref)
            if min(s0,l0)<=0:raise ArithmeticError('Nonpositive reference slope')
            references.append(dict(ID=label,reference_end_deg=ref,S0=s0))
            for h in cfg['derivative_halfwidth_deg']:
                signed,stencils,valid=refined_slopes(beta,y,archived_d,h)
                alternate,_,_=refined_slopes(beta,z,dz,h);q=abs(signed)/s0
                normalization_error=max(normalization_error,float(max(abs(q-abs(alternate)/l0))))
                derivative_error=max(derivative_error,float(max(abs(signed-2*alternate))))
                last=None
                for level in sorted(cfg['q_levels'],reverse=True):
                    result,envelope=refined_pick(beta,q,level,cfg['minimum_length_deg']);b=result['beta_pl_deg']
                    if b is not None and last is not None and b<last:raise ArithmeticError('Threshold monotonicity failed')
                    last=b;idx=None if b is None else int(np.where(beta==b)[0][0])
                    transitions.append(dict(ID=label,h_deg=h,reference_end_deg=ref,q_star=level,**result,
                        Lambda1=None if idx is None else table[idx]['Lambda1'],
                        R_percent=None if idx is None else 100*tails[idx]['R'],
                        late_q_max_after12=float(max(q[beta>=12])),
                        beta_old=oldangles.get((label,level)),
                        difference_deg=None if b is None else b-oldangles[label,level]))
                for i,r in enumerate(table):
                    profiles.append(dict(ID=label,h_deg=h,reference_end_deg=ref,beta_deg=float(beta[i]),
                        omega1=r['omega1'],Lambda1=r['Lambda1'],signed_S=float(signed[i]),q=float(q[i]),
                        Q_tail=float(envelope[i]),Q_peak_beta=float(beta[i+np.argmax(q[i:])]),
                        old_q=None if not beta[i].is_integer() else abs(float(archived_d[int(beta[i])]))/s0,
                        R_percent=100*tails[i]['R'],stencil=stencils[i],stencil_comparison_valid=valid[i]))
                if ref==5:
                    for b in (3,3.25,3.5,4,11,11.5,11.75,12,13):
                        i=int(np.where(beta==b)[0][0]);stitch.append(dict(ID=label,h_deg=h,beta_deg=b,
                            q=float(q[i]),stencil=stencils[i],all_thresholds_pass=bool(q[i]<=min(cfg['q_levels']))))
    if normalization_error>cfg['numerical_checks']['q_identity_absolute'] or derivative_error>cfg['numerical_checks']['derivative_identity_absolute']:
        raise ArithmeticError('Refined normalization identity failed')
    lookup={(r['ID'],r['q_star'],r['h_deg'],r['reference_end_deg']):r for r in transitions}
    primary=[r for r in transitions if r['h_deg']==.25 and r['reference_end_deg']==5]
    sensitivity=[];contrast=[];order=[]
    for label in cfg['cases']:
        for level in cfg['q_levels']:
            base=lookup[label,level,.25,5]['beta_pl_deg']
            for h in cfg['derivative_halfwidth_deg']:
                for ref in cfg['reference_intervals_end_deg']:
                    b=lookup[label,level,h,ref]['beta_pl_deg']
                    sensitivity.append(dict(ID=label,q_star=level,h_deg=h,reference_end_deg=ref,beta_pl_deg=b,
                        delta_from_primary_deg=None if None in (b,base) else b-base))
        vals={q:lookup[label,q,.25,5]['beta_pl_deg'] for q in cfg['q_levels']}
        contrast.append(dict(ID=label,beta_q010=vals[.1],beta_q015=vals[.15],beta_q020=vals[.2],
            q010_minus_q020=None if None in vals.values() else vals[.1]-vals[.2]))
    for a,pair in cfg['cases'].items():
        for b,upper in cfg['cases'].items():
            if a>=b:continue
            direction='joint' if pair[1]==upper[1] else 'end' if pair[0]==upper[0] else 'diagonal'
            for level in cfg['q_levels']:
                oldmargin=oldangles[b,level]-oldangles[a,level]
                for h in cfg['derivative_halfwidth_deg']:
                    for ref in cfg['reference_intervals_end_deg']:
                        av=lookup[a,level,h,ref]['beta_pl_deg'];bv=lookup[b,level,h,ref]['beta_pl_deg']
                        margin=None if None in (av,bv) else bv-av
                        order.append(dict(ID_A=a,ID_B=b,direction=direction,q_star=level,h_deg=h,reference_end_deg=ref,
                            beta_B_minus_A=margin,old_beta_B_minus_A=oldmargin,
                            strict_reversal=margin is not None and margin*oldmargin<0))
    for name,data in [('profiles',profiles),('transitions',transitions),('old_new_comparison',primary),
        ('reference_slopes',references),('sensitivity',sensitivity),('criterion_difference',contrast),('pair_order',order),('stitch_checks',stitch)]:
        csvwrite(folder/(name+'.csv'),data)
    return dict(q_normalization_max_abs_error=normalization_error,derivative_identity_max_abs_error=derivative_error,
        accepted_transitions=sum(r['beta_pl_deg'] is not None for r in transitions),
        missing_transitions=sum(r['beta_pl_deg'] is None for r in transitions))


def plot_refinement(folder):
    """Read only derived CSVs: six figures, no model assembly or postprocessing."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    folder=Path(folder);names=['contract.json','profiles.csv','old_new_comparison.csv']
    if (folder/'manifest.json').exists():
        m=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f)!=m['artifacts'][f] for f in names):raise ValueError('Refinement plot sources changed')
    cfg=json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    read=lambda f:list(csv.DictReader((folder/f).open(encoding='utf8')))
    profiles=read('profiles.csv');trans=read('old_new_comparison.csv')
    figures=folder/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,stem):
        fig.savefig(figures/(stem+'.png'),dpi=300,bbox_inches='tight')
        fig.savefig(figures/(stem+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    def part(label,h=.25):
        return [r for r in profiles if r['ID']==label and float(r['h_deg'])==h and int(r['reference_end_deg'])==5]
    arr=lambda rows,k:np.array([float(r[k]) for r in rows])
    colors={.1:'#d62728',.15:'#ff7f0e',.2:'#2ca02c'}
    accepted=[float(r['beta_pl_deg']) for r in trans if r['beta_pl_deg']];lo=min(accepted) if accepted else 0;hi=max(accepted) if accepted else 60
    for level in cfg['q_levels']:
        matrix=np.full((3,3),np.nan)
        for r in trans:
            if float(r['q_star'])==level and r['beta_pl_deg']:
                j,e=cfg['cases'][r['ID']];matrix[[.1,1,10].index(e),[.1,1,10].index(j)]=float(r['beta_pl_deg'])
        fig,ax=plt.subplots(figsize=(5.4,4.4),layout='constrained');cmap=plt.colormaps['viridis'].copy();cmap.set_bad('0.85')
        mesh=ax.pcolormesh(10.**np.arange(-1.5,2),10.**np.arange(-1.5,2),np.ma.masked_invalid(matrix),cmap=cmap,norm=Normalize(lo,hi),rasterized=False)
        for y,e in enumerate((.1,1,10)):
            for x,j in enumerate((.1,1,10)):
                v=matrix[y,x];ax.text(j,e,'—' if np.isnan(v) else f'{v:.2f}°',ha='center',va='center',fontsize=14,color='white' if v<(lo+hi)/2 else 'black')
        ax.set(xscale='log',yscale='log',xticks=[.1,1,10],yticks=[.1,1,10],xticklabels=['0.1','1','10'],yticklabels=['0.1','1','10'],xlabel='κ_joint',ylabel='κ_end',title=f'q* = {level:.2f}'+(' — рабочий вариант' if level==.15 else ''))
        ax.minorticks_off();cb=fig.colorbar(mesh,ax=ax,label='β_pl, °');cb.solids.set_rasterized(False);save(fig,f'01_map_q{level:.2f}')
    fig,axes=plt.subplots(3,3,figsize=(11,8),sharex=True,sharey=True,layout='constrained')
    for label,(j,e) in cfg['cases'].items():
        ax=axes[[.1,1,10].index(e),[.1,1,10].index(j)]
        for h,color,ls in [(.25,'#1f77b4','-'),(.5,'#9467bd','--')]:
            rows=[r for r in part(label,h) if 3<=float(r['beta_deg'])<=12]
            ax.plot(arr(rows,'beta_deg'),arr(rows,'q'),color=color,ls=ls,lw=1.3,label=f'h = {h:g}°')
        for q,color in colors.items():ax.axhline(q,color=color,ls=':',lw=.9,label=f'q* = {q:.2f}')
        ax.set(xlim=(3,12),ylim=(0,1.1*max(float(r['q']) for r in profiles
            if int(r['reference_end_deg'])==5 and 3<=float(r['beta_deg'])<=12)),
            title=f'{label}: ({j:g}; {e:g})',xticks=[3,5,7,9,12]);ax.grid(alpha=.15)
    axes[0,0].legend(fontsize=7,ncol=2)
    for ax in axes[-1]:ax.set_xlabel('β, °')
    for ax in axes[:,0]:ax.set_ylabel('q')
    save(fig,'02_q_stencils')
    fig,axes=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
    for ax,label in zip(axes,cfg['example_ids']):
        rows=[r for r in part(label) if 3<=float(r['beta_deg'])<=12];old=[r for r in rows if r['old_q']]
        ax.plot(arr(old,'beta_deg'),arr(old,'old_q'),'o--',color='0.5',lw=1,ms=3,label='Архив: h = 1°')
        ax.plot(arr(rows,'beta_deg'),arr(rows,'q'),color='#1f77b4',lw=1.4,label='Уточнение: h = 0.25°')
        for q,color in colors.items():ax.axhline(q,color=color,ls=':',lw=.8)
        ax.set(xlim=(3,12),ylim=(0,.45),xlabel='β, °',ylabel='q',title=label);ax.grid(alpha=.18);ax.legend(fontsize=7)
    save(fig,'03_resolution')
    fig,axes=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
    for ax,label in zip(axes,cfg['example_ids']):
        rows=[r for r in part(label) if float(r['beta_deg'])<=20]
        ax.plot(arr(rows,'beta_deg'),arr(rows,'Lambda1'),color='#1f77b4',lw=1.4)
        r=next(r for r in trans if r['ID']==label and float(r['q_star'])==.15)
        if r['beta_pl_deg']:ax.axvline(float(r['beta_pl_deg']),color=colors[.15],ls='--',label=f'β_pl = {float(r["beta_pl_deg"]):.2f}°')
        ax.set(xlim=(0,20),xlabel='β, °',ylabel='Λ₁',title=label);ax.grid(alpha=.18);ax.legend(fontsize=8)
    save(fig,'04_frequency')


def refinement_main(args):
    cfg=json.loads(REFINE_CONFIG.read_text(encoding='utf8'));key,ident,tables=refinement_identity(cfg)
    folder=OUTPUT/'growth_transition_refinement'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        target=(args.verify_cache or folder).resolve();m=json.loads((target/'manifest.json').read_text(encoding='utf8'))
        if m['identity']!=ident or any(sha(target/f)!=h for f,h in m['artifacts'].items()):
            raise ValueError('Refinement cache mismatch; no overwrite')
        print(f'REUSED {target}: 0 eigenvalue solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident)
    run=Run(cfg,folder);start=time.perf_counter();ledgerpath=folder/'run_ledger.json'
    ledger=json.loads(ledgerpath.read_text()) if ledgerpath.exists() else []
    try:
        result=refinement_spectra(run,args.stage=='preflight')
        if result is not None:
            rows,checks,neighbours=result;diagnostics=refinement_analysis(folder,cfg,tables,rows)
            diagnostics.update(status='PASS',maxima={f:max(r[f] for r in checks) for f in
                ('physical','ode','energy','singular','condition','mass_gram','mass_error')},
                neighbour_max_relative_jump=max(r['relative_jump'] for r in neighbours),
                minimum_gap12_relative=min(r['gap12_relative'] for r in neighbours))
            write(folder/'diagnostics.json',diagnostics);plot_refinement(folder)
    except Exception as exc:
        write(folder/'failure.json',dict(error=repr(exc),new_points=run.new_points,new_roots=run.new_roots));raise
    finally:
        ledger.append(dict(command=sys.argv,new_points=run.new_points,new_roots=run.new_roots,
            checkpoint_reuse_calls=run.reused_points,pole_roots=run.pole_roots,ritz_solves=run.ritz_solves,seconds=time.perf_counter()-start))
        write(ledgerpath,ledger)
    if result is not None:
        _,after,_=refinement_identity(cfg)
        if ident!=after:raise ValueError('Refinement sources changed during calculation')
        totals={k:sum(r[k] for r in ledger) for k in ('new_points','new_roots','pole_roots','ritz_solves')}
        if totals!=dict(new_points=243,new_roots=729,pole_roots=0,ritz_solves=0):raise ArithmeticError('Unexpected solve accounting')
        write(folder/'manifest.json',dict(status='PASS',identity=ident,
            git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            **totals,reused_local_points=90,reused_full_grid_points=819,local_points=333,merged_points=1062,
            complex_solves=0,postprocessing=dict(profiles=6372,transitions=162),
            artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['all', 'gates', 'reference', 'map', 'regression', 'preflight'], default='all')
    parser.add_argument('--preset', choices=['verification', 'limit-angles', 'finite-angle-maps', 'saturation-analysis', 'plateau-criterion', 'plateau-grid', 'early-saturation', 'growth-transition', 'growth-refinement'], default='verification')
    parser.add_argument('--plot-only', type=Path)
    parser.add_argument('--plot-set', choices=['original', 'lambda12', 'lambda12-ru'], default='original')
    parser.add_argument('--verify-cache', type=Path)
    parser.add_argument('--reuse-points-from', type=Path)
    args = parser.parse_args()
    if args.plot_set in ('lambda12', 'lambda12-ru') and (not args.plot_only or args.preset != 'limit-angles'):
        parser.error('lambda12 is available only with --preset limit-angles --plot-only PATH')
    if args.plot_only:
        if args.plot_set in ('lambda12', 'lambda12-ru'):
            plot_limit_pairs(args.plot_only, russian=args.plot_set == 'lambda12-ru')
        elif args.preset == 'finite-angle-maps':
            plot_finite_maps(args.plot_only)
        elif args.preset == 'saturation-analysis':
            plot_saturation(args.plot_only)
        elif args.preset == 'plateau-criterion':
            plot_plateau(args.plot_only)
        elif args.preset == 'plateau-grid':
            plot_grid(args.plot_only)
        elif args.preset == 'early-saturation':
            plot_early(args.plot_only)
        elif args.preset == 'growth-transition':
            plot_growth(args.plot_only)
        elif args.preset == 'growth-refinement':
            plot_refinement(args.plot_only)
        else:
            (plot_limits if args.preset == 'limit-angles' else plot)(args.plot_only)
        print('plot-only: 0 eigenvalue solves'); return
    if args.preset == 'limit-angles':
        limits_main(args); return
    if args.preset == 'finite-angle-maps':
        finite_main(args); return
    if args.preset == 'saturation-analysis':
        saturation_main(args); return
    if args.preset == 'plateau-criterion':
        plateau_main(args); return
    if args.preset == 'plateau-grid':
        grid_main(args); return
    if args.preset == 'early-saturation':
        early_main(args); return
    if args.preset == 'growth-transition':
        growth_main(args); return
    if args.preset == 'growth-refinement':
        refinement_main(args); return
    cfg = json.loads(CONFIG.read_text()); key, ident = identity(cfg)
    folder = OUTPUT/key
    if args.verify_cache:
        manifest = json.loads((args.verify_cache/'manifest.json').read_text())
        if manifest['identity'] != ident or any(sha(args.verify_cache/f) != h for f, h in manifest['artifacts'].items()):
            raise ValueError('Saved result identity/hash mismatch')
        print('cache verified: 0 eigenvalue solves'); return
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/'manifest.json').exists():
        manifest = json.loads((folder/'manifest.json').read_text())
        if manifest['identity'] == ident and all(sha(folder/f) == h for f, h in manifest['artifacts'].items()):
            print(f'REUSED {folder}: 0 eigenvalue solves'); return
        raise ValueError('Completed cache mismatch; no overwrite')
    write(folder/'contract.json', ident)
    run = Run(cfg, folder); start = time.perf_counter(); stages = {}
    reused_complete = False
    if args.reuse_points_from:
        reused_complete = reuse_points(run, args.reuse_points_from)
        if reused_complete:
            stages = json.loads((folder/'stages.json').read_text(encoding='utf8'))
    try:
        stages['gate1'] = structural(run); write(folder/'structural.json', stages['gate1'])
        if not reused_complete and args.stage in ('all', 'reference', 'map'):
            stages['gate2_reference'] = reference_stage(run)
            print('A1/A2 PASS', flush=True)
        if not reused_complete and args.stage in ('all', 'map'):
            stages['gate2_map'] = map_stage(run); print('A3/A4 complete', flush=True)
        if not reused_complete and args.stage in ('all', 'regression'):
            stages['gate3'] = regression_stage(run); print('B1-B4 complete', flush=True)
            stages['supplemental'] = supplemental_checks(run)
    except Exception as exc:
        write(folder/'failure.json', {'error': repr(exc), 'stages': stages,
                'new_points': run.new_points, 'new_roots': run.new_roots})
        raise
    write(folder/'stages.json', stages)
    if not reused_complete:
        write(folder/'pole_catalogs.json', {str(k): v for k, v in run.catalogs.items()})
    if args.stage == 'all':
        plot(folder)
        manifest = {'status': 'PASS', 'identity': ident, 'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            'git_status': subprocess.check_output(['git', 'status', '--short'], text=True),
            'command': sys.argv, 'new_points': run.new_points, 'reused_points': run.reused_points,
            'new_roots': run.new_roots, 'legacy_roots': run.old_roots,
            'pole_roots': run.pole_roots, 'ritz_solves': run.ritz_solves,
            'completed_checks_reused': reused_complete,
            'elapsed_seconds': time.perf_counter()-start}
        manifest['artifacts'] = {p.relative_to(folder).as_posix(): sha(p) for p in folder.rglob('*') if p.is_file() and p.name != 'manifest.json'}
        write(folder/'manifest.json', manifest)
    print(folder, flush=True)


if __name__ == '__main__':
    main()
