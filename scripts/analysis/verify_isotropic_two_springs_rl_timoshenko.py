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
STIFFNESS_CONFIG = ROOT/"data/input/isotropic_two_springs_growth_stiffness_grid.json"
MECHANICS_CONFIG = ROOT/"data/input/isotropic_two_springs_mechanical_interpretation.json"
VARIATION_CONFIG = ROOT/"data/input/isotropic_two_springs_angle_variation.json"
EPSILON_CONFIG = ROOT/"data/input/isotropic_two_springs_epsilon_scaling.json"
SLIP_CONFIG = ROOT/"data/input/isotropic_two_springs_longitudinal_slip.json"
THICKNESS_CONFIG = ROOT/"data/input/isotropic_two_springs_thickness_grid.json"


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
    def __init__(self, config, folder, engine=physics):
        self.config, self.folder, self.policy = config, folder, config["policy"]
        self.engine = engine
        self.model = engine.Model(**config["material"])
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
        roots, search = self.engine.prefix(model, case, self.cats(model, case), self.policy, neighbour)
        modes = [self.engine.mode(model, case, r["omega"], self.policy["quadrature_order"]) for r in roots]
        gram = np.array([[self.engine.inner(model, a, b) for b in modes] for a in modes])
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


def stiffness_slopes(beta, values):
    """Explicit, polynomial-exact stencils on the fixed 83-node composite grid."""
    expected=sorted(set(range(21)) | {3+i*.25 for i in range(37)} | set(range(22,91,2)))
    if list(beta)!=expected or len(values)!=83:raise ValueError('Unexpected composite grid')
    y=dict(zip(beta,values));out=[];stencils=[]
    for b in beta:
        if b==0:d=(-3*y[0]+4*y[1]-y[2])/2;name='left_second_order_h1'
        elif b==90:d=(3*y[90]-4*y[88]+y[86])/4;name='right_second_order_h2'
        elif b==20:d=-2*y[19]/3+y[20]/2+y[22]/6;name='nonuniform_19_20_22'
        else:
            h=.25 if 3<b<12 else 1 if b<20 else 2
            d=(y[b+h]-y[b-h])/(2*h);name=f'centered_h{h:g}'
        out.append(d);stencils.append(name)
    return np.array(out),stencils


def stiffness_sources(cfg):
    """Validate immutable archives; select exact physical pairs, never nearest modes."""
    _,ident=identity(cfg);ident['files'][STIFFNESS_CONFIG.relative_to(ROOT).as_posix()]=sha(STIFFNESS_CONFIG)
    sources={}
    for name in ('source','refined','limit','catalog'):
        folder=ROOT/cfg[name+'_bundle'];m=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if sha(folder/'manifest.json')!=cfg[name+'_manifest_sha256'] or m['status']!='PASS':raise ValueError('Stiffness source manifest changed')
        if any(sha(folder/f)!=v for f,v in m['artifacts'].items()):raise ValueError('Stiffness source artifact changed')
        sources[name]=m
    previous=sources['source']['identity']
    for name in ('material','mu','policy'):
        if previous['config'][name]!=cfg[name]:raise ValueError('Physical or quality contract changed')
    for name,h in previous['files'].items():
        if name!=Path(__file__).relative_to(ROOT).as_posix() and sha(ROOT/name)!=h:raise ValueError('Verified implementation changed: '+name)
    values=[10**x for x in (-1,-.5,0,.5,1)]
    if cfg['stiffness_values']!=values or cfg['cases']!={f'J{i}E{j}':[a,b] for i,a in enumerate(values) for j,b in enumerate(values)}:
        raise ValueError('Stiffness exponents must be evaluated without decimal truncation')
    oldcases=previous['config']['cases'];mapping={old:next(k for k,v in cfg['cases'].items() if v==pair) for old,pair in oldcases.items()}
    archives={}
    for field in ('source','refined'):
        folder=ROOT/cfg[field+'_bundle']
        for row in csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')):
            label=mapping[row['ID']];b=float(row['beta_deg'])
            if b not in cfg['map']['parameter_grid']:continue
            path=folder/'points'/row['point_file'] if field=='source' else ROOT/row['point_file']
            rec=json.loads(path.read_text());case=rec['case']
            if (case['beta'],case['mu'],case['k_joint'],case['k_end'],case['variant'])!=(b,cfg['mu'],*cfg['cases'][label],cfg['variant']):raise ValueError('Archive case mismatch')
            if rec['status']!='PASS' or rec['count']['status']!='PASS' or len(rec['roots'])!=3:raise ValueError('Archive prefix not accepted')
            for i,root in enumerate(rec['roots'],1):
                if not all(root['gates'].values()) or any(float(row[f'{f}{i}'])!=root[f] for f in ('omega','Omega','Lambda')):raise ValueError('Archive CSV/root mismatch')
            if sha(path.with_suffix('.npz'))!=rec['shape_sha256']:raise ValueError('Archive shapes changed')
            archives[label,b]=(rec,path)
    if len(archives)!=747:raise ValueError('Missing archived composite nodes')
    read=lambda field,name:list(csv.DictReader((ROOT/cfg[field+'_bundle']/name).open(encoding='utf8')))
    profiles=[{**r,'ID':mapping[r['ID']]} for r in read('refined','profiles.csv') if r['h_deg']=='0.25' and r['reference_end_deg']=='5']
    transitions=[{**r,'ID':mapping[r['ID']]} for r in read('refined','old_new_comparison.csv')]
    limits={(r['ID'],float(r['beta_deg'])):r for r in read('limit','spectrum.csv') if r['ID'] in ('HH','RR')}
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    return key,ident,archives,profiles,transitions,limits


def stiffness_thinning(cfg,profiles,transitions):
    beta=cfg['map']['parameter_grid'];out=[]
    for label in sorted({r['ID'] for r in profiles}):
        old={float(r['beta_deg']):r for r in profiles if r['ID']==label}
        y=np.log([float(old[b]['omega1']) for b in beta]);s0=(y[beta.index(5)]-y[0])/5
        slope,_=stiffness_slopes(beta,y);q=abs(slope)/s0
        omitted=max(float(r['q']) for b,r in old.items() if b not in beta)
        for level in cfg['q_levels']:
            got,_=refined_pick(beta,q,level,cfg['minimum_length_deg'])
            prior=next(r for r in transitions if r['ID']==label and float(r['q_star'])==level)
            passed=got['beta_pl_deg']==float(prior['beta_pl_deg']) and omitted<level
            out.append(dict(ID=label,q_star=level,beta_archive=float(prior['beta_pl_deg']),beta_thinned=got['beta_pl_deg'],
                omitted_q_max=omitted,tail20_q_max=float(max(v for b,v in zip(beta,q) if b>=20)),status='PASS' if passed else 'FAIL'))
            if not passed:raise ArithmeticError('Thinning changes transition or hides a late exceedance; review grid first')
    return out


def stiffness_spectra(run,archives,limits,preflight_only=False):
    cfg=run.config;beta=cfg['map']['parameter_grid'];catalog=ROOT/cfg['catalog_bundle']/'pole_catalogs.json'
    run.catalogs={float(k):v for k,v in json.loads(catalog.read_text()).items()}
    auditfile=run.folder/'count_checks.json';audits=json.loads(auditfile.read_text()) if auditfile.exists() else {}
    def enclosure(label,b,w):
        pair=cfg['cases'][label];low=[];high=[]
        for (a,angle),(rec,_) in archives.items():
            if angle!=b:continue
            other=cfg['cases'][a];freq=[r['omega'] for r in rec['roots']]
            if all(x<=y for x,y in zip(other,pair)):low.append(freq)
            if all(x>=y for x,y in zip(other,pair)):high.append(freq)
        margins=np.r_[w/np.max(low,axis=0)-1,np.min(high,axis=0)/w-1]
        if min(margins)<-cfg['comparison_relative_tolerance']:raise ArithmeticError('Archive stiffness enclosure failed')
        return float(min(margins))
    def point(label,b,previous=None):
        if (label,b) in archives:return (*archives[label,b],'archive')
        if label not in cfg['new_case_ids'] or b not in beta:raise ValueError('Unrequested spectral point')
        case=physics.Case(b,cfg['mu'],*cfg['cases'][label]);contract={**asdict(case),'variant':'RL_FULL'}
        key=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()[:16];path=run.folder/'points'/(key+'.json')
        if not path.exists() and len(list((run.folder/'points').glob('*.json')))>=cfg['expected_new_points']:raise ArithmeticError('Point budget exceeded')
        rec=run.solve(case,neighbour=previous)
        if rec['status']!='PASS' or rec['count']['status']!='PASS' or not all(all(r['gates'].values()) for r in rec['roots']):raise ArithmeticError('Unaccepted prefix')
        w=np.array([r['omega'] for r in rec['roots']]);margin=enclosure(label,b,w)
        if key not in audits:
            counts=[]
            for query,n0 in zip([float(run.model.omega(run.policy['Lambda_lower'])),*(.5*(w[:-1]+w[1:]))],(0,1,2)):
                _,n,d=physics.count_at(run.model,case,float(query),run.cats(run.model,case),run.policy)
                if n!=n0:raise ArithmeticError('Independent low/midgap count failed')
                counts.append(dict(query=float(query),expected=n0,count=n,**d))
            audits[key]=dict(ID=label,beta_deg=b,queries=counts,archive_enclosure_min_margin=margin,status='PASS');write(auditfile,audits)
        return rec,path,'new'
    pre=[]
    for label in cfg['preflight_case_ids']:
        previous=None
        for b in cfg['preflight_angles']:
            rec,_,_=point(label,b,previous);previous=np.array([r['omega'] for r in rec['roots']])
            pre.append(dict(ID=label,beta_deg=b,Lambda=[r['Lambda'] for r in rec['roots']],status='PASS'))
    write(run.folder/'preflight.json',dict(status='PASS',points=pre));print('5x5 preflight: 36 prefixes PASS',flush=True)
    if preflight_only:return None
    rows=[];checks=[];neighbours=[]
    for label in cfg['cases']:
        previous=None
        for b in beta:
            rec,path,origin=point(label,b,previous);w=np.array([r['omega'] for r in rec['roots']])
            if min(w)<=0 or np.any(np.diff(w)<=0):raise ArithmeticError('Unordered prefix')
            if previous is not None:
                jump=float(max(abs(w/previous-1)));neighbours.append(dict(ID=label,beta_deg=b,relative_jump=jump,
                    trigger=jump>run.policy['neighbour_relative_jump'],gap12_relative=float(w[1]/w[0]-1),
                    resolution='count-one intervals + independent low/midgap counts + physical gates + archived stiffness enclosure'))
            previous=w;row=dict(ID=label,beta_deg=b)
            for i,root in enumerate(rec['roots'],1):
                row.update({f'{f}{i}':root[f] for f in ('omega','Omega','Lambda')})
                if not (math.isclose(root['omega']*200,root['Omega'],rel_tol=2e-14) and math.isclose(root['Lambda']**2,root['Omega'],rel_tol=2e-14)):raise ArithmeticError('Normalization failed')
                d=root['diagnostics'];checks.append(dict(ID=label,beta_deg=b,position=i,origin=origin,
                    physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],
                    energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],
                    condition=d['nonzero_condition'],status=rec['status']))
            row.update(origin=origin,point_file=path.relative_to(ROOT).as_posix(),point_sha256=sha(path),shape_sha256=rec['shape_sha256']);rows.append(row)
        print(f'5x5 {label}: 83 prefixes '+('reused' if label not in cfg['new_case_ids'] else 'accepted'),flush=True)
    if len(rows)!=2075 or len(audits)!=1328:raise ArithmeticError('Grid accounting mismatch')
    by={(r['ID'],r['beta_deg']):r for r in rows};order=[]
    for j in range(5):
        for e in range(5):
            a=f'J{j}E{e}'
            for other,axis in [(f'J{j+1}E{e}','joint')]*(j<4)+[(f'J{j}E{e+1}','end')]*(e<4):
                for b in beta:
                    for k in (1,2,3):
                        margin=by[other,b][f'omega{k}']/by[a,b][f'omega{k}']-1
                        if margin < -cfg['comparison_relative_tolerance']:raise ArithmeticError('Spectral stiffness min-max failed')
                        order.append(dict(axis=axis,lower=a,upper=other,beta_deg=b,position=k,relative_margin=margin))
            for b in beta:
                if ('HH',b) not in limits:continue
                for k in (1,2,3):
                    w=by[a,b][f'omega{k}']
                    for lo,hi,margin in [('HH',a,w/float(limits['HH',b][f'omega{k}'])-1),(a,'RR',float(limits['RR',b][f'omega{k}'])/w-1)]:
                        if margin < -cfg['comparison_relative_tolerance']:raise ArithmeticError('HH/RR limit enclosure failed')
                        order.append(dict(axis='limit',lower=lo,upper=hi,beta_deg=b,position=k,relative_margin=margin))
    for name,data in [('spectrum',rows),('reused_nodes',[r for r in rows if r['origin']=='archive']),('checks',checks),('neighbour_audit',neighbours),('spectral_order',order)]:csvwrite(run.folder/(name+'.csv'),data)
    return rows,checks,neighbours


def stiffness_effects(transitions,frequencies):
    """Finite increments per half decade and mixed differences, not smooth fits."""
    angles={(r['ID'],r['q_star']):r['beta_pl_deg'] for r in transitions};levels=sorted({r['q_star'] for r in transitions})
    freq={r['ID']:r for r in frequencies};increments=[];mixed=[]
    for level in levels:
        for j in range(5):
            for e in range(5):
                a=f'J{j}E{e}'
                for b,axis in [(f'J{j+1}E{e}','joint')]*(j<4)+[(f'J{j}E{e+1}','end')]*(e<4):
                    av,bv=angles[a,level],angles[b,level]
                    increments.append(dict(q_star=level,axis=axis,lower=a,upper=b,log10_step=.5,
                        delta_beta_deg=None if None in (av,bv) else bv-av,
                        omega10_change_percent=100*((freq[b]['Lambda10']/freq[a]['Lambda10'])**2-1),
                        omega90_change_percent=100*((freq[b]['Lambda90']/freq[a]['Lambda90'])**2-1)))
                if j<4 and e<4:
                    ids=[a,f'J{j+1}E{e}',f'J{j}E{e+1}',f'J{j+1}E{e+1}'];v=[angles[x,level] for x in ids]
                    mixed.append(dict(q_star=level,lower_joint_index=j,lower_end_index=e,
                        mixed_difference_deg=None if None in v else v[3]-v[2]-v[1]+v[0]))
    return increments,mixed


def stiffness_analysis(folder,cfg,rows,oldtrans):
    profiles=[];transitions=[];frequencies=[];normalization_error=0.;stitches=[]
    for label,pair in cfg['cases'].items():
        table=[r for r in rows if r['ID']==label];beta=[r['beta_deg'] for r in table];by={r['beta_deg']:r for r in table}
        y=np.log([r['omega1'] for r in table]);z=np.log([r['Lambda1'] for r in table]);s0=float((y[beta.index(5)]-y[0])/5);l0=float((z[beta.index(5)]-z[0])/5)
        if s0<=0:raise ArithmeticError('Nonpositive initial slope')
        signed,stencils=stiffness_slopes(beta,y);other,_=stiffness_slopes(beta,z);q=abs(signed)/s0
        normalization_error=max(normalization_error,float(max(abs(q-abs(other)/l0))))
        if normalization_error>cfg['numerical_checks']['q_identity_absolute']:raise ArithmeticError('q normalization identity failed')
        tails=tail_ranges(table);last=None
        if max(v for b,v in zip(beta,q) if b>=20)>=min(cfg['q_levels']):raise ArithmeticError('Late sensitivity trigger; inspect locally before acceptance')
        for level in sorted(cfg['q_levels'],reverse=True):
            result,Q=refined_pick(beta,q,level,cfg['minimum_length_deg']);b=result['beta_pl_deg']
            first=next((b0 for b0,v in zip(beta,q) if b0>=5 and v<=level),None)
            if b is not None and first is not None and b>first:raise ArithmeticError('Late peak delays local transition; inspect locally')
            if b is not None and last is not None and b<last:raise ArithmeticError('Threshold order failed')
            last=b;idx=None if b is None else beta.index(b)
            prior=next((r for r in oldtrans if r['ID']==label and float(r['q_star'])==level),None)
            if prior is not None and b!=float(prior['beta_pl_deg']):raise ArithmeticError('Accepted archived transition changed')
            transitions.append(dict(ID=label,kappa_joint=pair[0],kappa_end=pair[1],q_star=level,**result,
                Lambda_pl=None if b is None else by[b]['Lambda1'],R_percent=None if idx is None else 100*tails[idx]['R'],
                transition_origin='archived_accepted' if prior is not None else 'new_composite_grid'))
            if level==.15:frequencies.append(dict(ID=label,kappa_joint=pair[0],kappa_end=pair[1],beta_pl_deg=b,
                Lambda_pl=None if b is None else by[b]['Lambda1'],Lambda0=by[0]['Lambda1'],Lambda10=by[10]['Lambda1'],Lambda90=by[90]['Lambda1'],S0=s0))
        for i,r in enumerate(table):
            row=dict(ID=label,beta_deg=beta[i],omega1=r['omega1'],Lambda1=r['Lambda1'],signed_S=float(signed[i]),S0=s0,
                q=float(q[i]),Q_tail=float(Q[i]),tail_peak_beta=beta[i+int(np.argmax(q[i:]))],stencil=stencils[i],R_percent=100*tails[i]['R'])
            profiles.append(row)
            if beta[i] in (3,3.25,11.75,12,13,19,20,22,24,88,90):stitches.append(row)
    increments,mixed=stiffness_effects(transitions,frequencies)
    for name,data in [('profiles',profiles),('transitions',transitions),('frequency_levels',frequencies),('increments',increments),('mixed_differences',mixed),('stitch_checks',stitches)]:csvwrite(folder/(name+'.csv'),data)
    # A small optional diagnostic reads already saved mode1 arrays; no recovery or tracking.
    energies=[];m=physics.Model(**cfg['material']);s=m.section;by={(r['ID'],r['beta_deg']):r for r in rows}
    for label in cfg['energy_case_ids']:
        f=next(r for r in frequencies if r['ID']==label)
        for b in sorted({0,10,f['beta_pl_deg']} - {None}):
            r=by[label,b];path=ROOT/r['point_file'];data=np.load(path.with_suffix('.npz'));mass=rod=0.;ends=[]
            for arm in (1,2):
                y=data[f'mode1_arm{arm}_y'];dy=data[f'mode1_arm{arm}_dy'];wt=data[f'mode1_arm{arm}_weights'];ends.append(data[f'mode1_arm{arm}_ends'])
                mass+=wt@(s.rhoA*(y[:,0]**2+y[:,1]**2)+s.rhoI*y[:,2]**2+m.J*dy[:,0]**2)
                rod+=wt@(s.EA*dy[:,0]**2+s.EI*dy[:,2]**2+s.KGA*(dy[:,1]-y[:,2])**2)
            delta=ends[0][1,2]-ends[1][1,2];rotation=ends[1][0,2];kj,ke=cfg['cases'][label]
            joint=kj*s.EI/m.l*delta**2;end=ke*s.EI/m.l*rotation**2;total=rod+joint+end;error=abs(total/(mass*r['omega1']**2)-1)
            if error>cfg['policy']['energy_tol'] or abs(mass-1)>cfg['policy']['mass_gram_tol']:raise ArithmeticError('Saved-field energy decomposition failed')
            energies.append(dict(ID=label,beta_deg=b,delta_psi=float(delta),psi_end=float(rotation),mass=float(mass),
                V_joint=float(joint/2),V_end=float(end/2),V_rods=float(rod/2),joint_fraction=float(joint/total),
                end_fraction=float(end/total),energy_error=float(error),source_sha256=sha(path.with_suffix('.npz'))))
            data.close()
    csvwrite(folder/'energy_samples.csv',energies)
    return dict(q_normalization_max_abs_error=normalization_error,accepted_transitions=sum(r['beta_pl_deg'] is not None for r in transitions),
        missing_transitions=sum(r['beta_pl_deg'] is None for r in transitions),tail20_q_max=max(r['q'] for r in profiles if r['beta_deg']>=20),
        energy_samples=len(energies),energy_decomposition_max_error=max(r['energy_error'] for r in energies))


def plot_stiffness(folder):
    """Five saved-table figures; 3D facets are explicitly graphical interpolation."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    folder=Path(folder);names=['contract.json','transitions.csv','frequency_levels.csv']
    if (folder/'manifest.json').exists():
        m=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f)!=m['artifacts'][f] for f in names):raise ValueError('Stiffness plot inputs changed')
    cfg=json.loads((folder/'contract.json').read_text(encoding='utf8'))['config'];v=cfg['stiffness_values']
    read=lambda n:list(csv.DictReader((folder/n).open(encoding='utf8')))
    t=read('transitions.csv');f=read('frequency_levels.csv');primary={r['ID']:r for r in t if float(r['q_star'])==.15}
    matrix=np.array([[float(primary[f'J{j}E{e}']['beta_pl_deg']) if primary[f'J{j}E{e}']['beta_pl_deg'] else np.nan for j in range(5)] for e in range(5)])
    figures=folder/'figures';figures.mkdir(exist_ok=True);labels=['0.1','0.316','1','3.162','10'];colors=plt.colormaps['viridis'](np.linspace(.05,.85,5));styles=['-','--','-.',':','-'];markers=['o','s','^','D','v']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.savefig(figures/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(figures/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.2,5.1),layout='constrained');edges=10.**np.arange(-1.25,1.51,.5)
    cmap=plt.colormaps['viridis'].copy();cmap.set_bad('0.85');mesh=ax.pcolormesh(edges,edges,np.ma.masked_invalid(matrix),cmap=cmap,vmin=np.nanmin(matrix),vmax=np.nanmax(matrix),rasterized=False)
    for e in range(5):
        for j in range(5):
            z=matrix[e,j];ax.text(v[j],v[e],'—' if np.isnan(z) else f'{z:.2f}°',ha='center',va='center',fontsize=11,color='white' if z<(np.nanmin(matrix)+np.nanmax(matrix))/2 else 'black')
    ax.set(xscale='log',yscale='log',xticks=v,yticks=v,xticklabels=labels,yticklabels=labels,xlabel='κ_joint',ylabel='κ_end',title='q* = 0.15');ax.minorticks_off();cb=fig.colorbar(mesh,ax=ax,label='β_pl, °');cb.solids.set_rasterized(False);save(fig,'01_transition_map')
    x,y=np.meshgrid(cfg['stiffness_exponents'],cfg['stiffness_exponents']);valid=np.isfinite(matrix.ravel())
    fig=plt.figure(figsize=(7.5,5.6));ax=fig.add_subplot(111,projection='3d')
    if np.all(valid):
        tri=Triangulation(x.ravel(),y.ravel());ax.plot_trisurf(tri,matrix.ravel(),cmap='viridis',alpha=.68,edgecolor='0.55',linewidth=.35)
    ax.scatter(x.ravel()[valid],y.ravel()[valid],matrix.ravel()[valid],color='black',s=22,depthshade=False)
    ax.set(xticks=cfg['stiffness_exponents'],yticks=cfg['stiffness_exponents'],xticklabels=labels,yticklabels=labels,xlabel='κ_joint (log)',ylabel='κ_end (log)',zlabel='β_pl, °')
    ax.set_zlabel('');ax.text2D(1.06,.58,'β_pl, °',transform=ax.transAxes,rotation=90,va='center')
    ax.view_init(elev=25,azim=-55);fig.text(.5,.04,'25 расчётных узлов; грани — графическая интерполяция',ha='center',fontsize=9);fig.subplots_adjust(bottom=.1,top=.98);save(fig,'02_surface')
    for axis,name in [('joint','03_joint_sections'),('end','04_end_sections')]:
        fig,ax=plt.subplots(figsize=(6.6,4.5),layout='constrained')
        for i in range(5):
            z=matrix[i,:] if axis=='joint' else matrix[:,i]
            ax.plot(v,z,color=colors[i],ls=styles[i],marker=markers[i],ms=5,lw=1.3,label=f'κ_{"end" if axis=="joint" else "joint"} = {labels[i]}')
        ax.set(xscale='log',xticks=v,xticklabels=labels,xlabel=f'κ_{axis}',ylabel='β_pl, °',ylim=(np.nanmin(matrix)-.15,np.nanmax(matrix)+.15));ax.minorticks_off();ax.grid(alpha=.18);ax.legend(fontsize=8,ncol=2);save(fig,name)
    # Frequency levels complement, rather than repeat, the transition maps.
    by={r['ID']:r for r in f};fig,axes=plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained')
    for ax,axis in zip(axes,('joint','end')):
        for fixed in cfg['plot_frequency_sections']:
            ids=[f'J{i}E{fixed}' if axis=='joint' else f'J{fixed}E{i}' for i in range(5)]
            ax.plot(v,[float(by[k]['Lambda90']) for k in ids],marker=markers[fixed],ms=4,color=colors[fixed],lw=1.3,label=f'κ_{"end" if axis=="joint" else "joint"} = {labels[fixed]}')
        ax.set(xscale='log',xticks=v,xticklabels=labels,xlabel=f'κ_{axis}',ylabel='Λ₁(90°)');ax.minorticks_off();ax.grid(alpha=.18);ax.legend(fontsize=8)
    save(fig,'05_frequency_levels')


def stiffness_main(args):
    cfg=json.loads(STIFFNESS_CONFIG.read_text());key,ident,archives,profiles,oldtrans,limits=stiffness_sources(cfg)
    folder=OUTPUT/'growth_stiffness_5x5'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        p=(args.verify_cache or folder).resolve();m=json.loads((p/'manifest.json').read_text(encoding='utf8'))
        if m['identity']!=ident or any(sha(p/f)!=h for f,h in m['artifacts'].items()):raise ValueError('Stiffness cache mismatch')
        print(f'REUSED {p}: 0 eigenvalue solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident)
    thin=stiffness_thinning(cfg,profiles,oldtrans);csvwrite(folder/'thinning_check.csv',thin)
    run=Run(cfg,folder);ledgerpath=folder/'run_ledger.json';ledger=json.loads(ledgerpath.read_text()) if ledgerpath.exists() else []
    start=time.perf_counter();spectral=postprocess=render=0.
    try:
        result=stiffness_spectra(run,archives,limits,args.stage=='preflight');spectral=time.perf_counter()-start
        if result is not None:
            rows,checks,neighbours=result;startpost=time.perf_counter();diag=stiffness_analysis(folder,cfg,rows,oldtrans);postprocess=time.perf_counter()-startpost
            diag.update(status='PASS',maxima={k:max(r[k] for r in checks) for k in ('physical','ode','energy','singular','mass_gram','condition')},
                neighbour_triggers=sum(r['trigger'] for r in neighbours),minimum_gap12_relative=min(r['gap12_relative'] for r in neighbours))
            write(folder/'diagnostics.json',diag);startplot=time.perf_counter();plot_stiffness(folder);render=time.perf_counter()-startplot
    except Exception as exc:
        write(folder/'failure.json',dict(error=repr(exc),new_points=run.new_points,new_roots=run.new_roots));raise
    finally:
        ledger.append(dict(command=sys.argv,new_points=run.new_points,new_roots=run.new_roots,checkpoint_reuse_calls=run.reused_points,
            pole_roots=run.pole_roots,ritz_solves=run.ritz_solves,spectral_seconds=spectral,postprocessing_seconds=postprocess,render_seconds=render,wall_seconds=time.perf_counter()-start));write(ledgerpath,ledger)
    if result is not None:
        _,after,*_=stiffness_sources(cfg)
        if ident!=after:raise ValueError('Sources changed during map')
        totals={k:sum(r[k] for r in ledger) for k in ('new_points','new_roots','pole_roots','ritz_solves')}
        if totals!=dict(new_points=1328,new_roots=3984,pole_roots=0,ritz_solves=0):raise ArithmeticError('Unexpected solve counters')
        write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            **totals,reused_points=747,archived_transitions_reused=27,total_parameter_points=2075,complex_solves=0,
            source_archive_points_for_thinning=1062,postprocessing=dict(profiles=2075,transitions=75),
            artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def mechanics_energy(model, case, arms):
    """Split the canonical energy of supplied fields; no new physical operator."""
    s=model.section;rod=[];mass=[]
    for a in arms:
        y,d,w=a['y'],a['dy'],a['weights']
        rod.append([float(w@(s.EA*d[:,0]**2)/2),float(w@(s.EI*d[:,2]**2)/2),
                    float(w@(s.KGA*(d[:,1]-y[:,2])**2)/2)])
        mass.append([float(w@(s.rhoA*y[:,0]**2)),float(w@(s.rhoA*y[:,1]**2)),
                     float(w@(s.rhoI*y[:,2]**2)),float(w@(model.J*d[:,0]**2))])
    jump=arms[0]['ends'][1,2]-arms[1]['ends'][1,2];end=arms[1]['ends'][0,2]
    energies=np.r_[np.sum(rod,axis=0),case.k_joint*s.EI/model.l*jump**2/2,case.k_end*s.EI/model.l*end**2/2]
    if min(energies)<0 or sum(energies)<=0:raise ArithmeticError('Nonpositive field energy')
    return dict(energies=energies,fractions=energies/sum(energies),total=float(sum(energies)),
                mass=float(np.sum(mass)),arm_energy=np.array(rod),mass_components=np.array(mass),
                delta_psi=float(jump),psi_end=float(end))


def mechanics_fields(model, case, omega, coefficients, order):
    """Re-evaluate an accepted analytic coefficient vector, without SVD/root search."""
    nodes,weights=np.polynomial.legendre.leggauss(order);arms=[]
    for L,c in zip(case.lengths(model),coefficients):
        x=(nodes+1)*L/2
        arms.append(dict(x=x,weights=weights*L/2,y=physics.basis(model,L,omega,x)@c,
            dy=physics.basis(model,L,omega,x,1)@c,ends=physics.basis(model,L,omega,[0.,L])@c))
    return arms


def mechanics_global(model, case, coefficients, omega, nodes):
    """Joint at the origin; both positive local coordinates point toward it."""
    fields=[]
    for L,c,f in zip(case.lengths(model),coefficients,physics.geometry.frames(case.beta)):
        x=np.linspace(0,L,nodes);y=physics.basis(model,L,omega,x)@c
        fields.append(dict(x=x,y=y,reference=(x-L)[:,None]*f.t,
            displacement=y[:,:2]@f.translation.T,axial=y[:,0,None]*f.t,transverse=y[:,1,None]*f.n))
    return fields


def mechanics_sources(cfg):
    _,ident=identity(cfg);ident['files'][MECHANICS_CONFIG.relative_to(ROOT).as_posix()]=sha(MECHANICS_CONFIG)
    source=ROOT/cfg['source_bundle'];m=json.loads((source/'manifest.json').read_text(encoding='utf8'))
    if sha(source/'manifest.json')!=cfg['source_manifest_sha256'] or m['status']!='PASS':raise ValueError('Mechanical source manifest changed')
    if any(sha(source/f)!=h for f,h in m['artifacts'].items()):raise ValueError('Mechanical source artifact changed')
    for k in ('material','mu','policy'):
        if m['identity']['config'][k]!=cfg[k]:raise ValueError('Mechanical physical/quality contract mismatch')
    for f,h in m['identity']['files'].items():
        if f!=Path(__file__).relative_to(ROOT).as_posix() and sha(ROOT/f)!=h:raise ValueError('Source implementation changed: '+f)
    catalog=ROOT/cfg['catalog_bundle'];cm=json.loads((catalog/'manifest.json').read_text(encoding='utf8'))
    if sha(catalog/'manifest.json')!=cfg['catalog_manifest_sha256'] or any(sha(catalog/f)!=h for f,h in cm['artifacts'].items()):raise ValueError('Pole catalog source changed')
    read=lambda n:list(csv.DictReader((source/n).open(encoding='utf8')))
    table={(r['ID'],float(r['beta_deg'])):r for r in read('spectrum.csv')};selected=[]
    transitions={r['ID']:r for r in read('transitions.csv') if float(r['q_star'])==.15}
    for label,pair in cfg['cases'].items():
        if pair!=m['identity']['config']['cases'][label] or float(transitions[label]['beta_pl_deg'])!=cfg['transition_angles'][label]:raise ValueError('Selected case or angle mismatch')
        for b in sorted(set(cfg['base_angles']+[cfg['transition_angles'][label]])):
            row=table[label,b];path=ROOT/row['point_file'];rec=json.loads(path.read_text(encoding='utf8'))
            if sha(path)!=row['point_sha256'] or sha(path.with_suffix('.npz'))!=row['shape_sha256']:raise ValueError('Selected point/shape changed')
            c=rec['case']
            if (c['beta'],c['mu'],c['k_joint'],c['k_end'],c['variant'])!=(b,cfg['mu'],*pair,'RL_FULL'):raise ValueError('Selected physical mismatch')
            if rec['status']!='PASS' or rec['count']['status']!='PASS':raise ValueError('Unaccepted source point')
            for k,r in enumerate(rec['roots'],1):
                if r['counts']!=[k-1,k] or not all(r['gates'].values()) or r['omega']!=float(row[f'omega{k}']):raise ValueError('Source prefix mismatch')
            selected.append((row,rec,path))
    if len(selected)!=29:raise ValueError('Exactly 29 accepted forms required')
    profiles=[r for r in read('profiles.csv') if r['ID'] in cfg['cases']]
    ident['source_manifest_sha256']=sha(source/'manifest.json')
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    return key,ident,selected,profiles


def mechanics_analysis(folder,cfg,selected,profiles):
    m=physics.Model(**cfg['material']);s=m.section;rows=[];checks=[];quads=[];provenance=[];plotmax=0.
    (folder/'forms').mkdir(exist_ok=True);components=['axial','bend','shear','joint','end'];massnames=['axial','transverse','rotation','RL_gradient']
    profile={(r['ID'],float(r['beta_deg'])):r for r in profiles}
    for source,rec,path in selected:
        label=source['ID'];b=float(source['beta_deg']);case=physics.Case(b,cfg['mu'],*cfg['cases'][label]);root=rec['roots'][0];w=root['omega'];coef=np.array(root['coefficients'])
        with np.load(path.with_suffix('.npz')) as data:
            arms=[{f:data[f'mode1_arm{i+1}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in range(2)]
        original=mechanics_energy(m,case,arms);mass_error=abs(original['mass']-1)
        factor=1/math.sqrt(original['mass']);coef*=factor
        for a in arms:
            for f in ('y','dy','ends'):a[f]*=factor
        # One fixed physical sign anchor, not nearest-frequency or modal tracking.
        anchor=-(physics.basis(m,case.lengths(m)[0],w,case.lengths(m)[0]/2)@coef[0])[1]
        scale=max(np.max(abs(a['y'][:,:2])) for a in arms)
        if abs(anchor)<1e-6*scale:raise ArithmeticError('Unstable shape sign anchor')
        sign=1 if anchor>0 else -1;coef*=sign
        for a in arms:
            for f in ('y','dy','ends'):a[f]*=sign
        en=mechanics_energy(m,case,arms);rebuilt=mechanics_fields(m,case,w,coef,120)
        reconstruction=max(np.max(abs(a[f]-z[f]))/max(np.max(abs(a[f])),1e-30) for a,z in zip(arms,rebuilt) for f in ('y','dy','ends'))
        fraction_spread=mass_spread=0.
        for order in cfg['quadrature_orders']:
            tested=mechanics_fields(m,case,w,coef,order);e=mechanics_energy(m,case,tested)
            df=float(max(abs(e['fractions']-en['fractions'])));dm=abs(e['mass']-en['mass']);fraction_spread=max(fraction_spread,df);mass_spread=max(mass_spread,dm)
            quads.append(dict(ID=label,beta_deg=b,order=order,mass=e['mass'],energy_error=abs(2*e['total']/(w*w*e['mass'])-1),fraction_max_abs_difference=df))
        qscale=max(np.max(abs(a['y'][:,:3]*[1,1,m.l])) for a in tested);fscale=max(np.max(abs(a['y'][:,3:]*[1,1,1/m.l])) for a in tested)
        joint=physics.joint_matrix(m,case)@np.concatenate([a['ends'][1] for a in tested]);physical=max(abs(joint)/[qscale,qscale,fscale,fscale,fscale*m.l,fscale*m.l])
        for i,a in enumerate(tested):
            er=physics.external_matrix(m,case,i)@a['ends'][0];physical=max(physical,max(abs(er)/[qscale,qscale,qscale/m.l if i==0 else fscale*m.l]))
        scales=np.array([1,1,m.l,qscale/fscale,qscale/fscale,qscale/(fscale*m.l)])
        ode=max(np.max(abs((a['dy']-a['y']@physics.state_matrix(m,w).T)*scales))/np.max(abs(a['dy']*scales)) for a in tested)
        # Alternate trace-matrix calculation of the SAME spring quadratic forms.
        frames=physics.geometry.frames(b);node=arms[0]['ends'][1,:2]@frames[0].translation.T
        z=np.r_[node,arms[0]['ends'][1,2],arms[1]['ends'][1,2],arms[1]['ends'][0,2]]
        trace_errors=[]
        for springname,index in [('joint',3),('end',4)]:
            only=replace(case,**{'k_end' if springname=='joint' else 'k_joint':0.})
            trace,stiff=physics.trace_map(m,only);via=float(z@stiff@z/2)
            trace_errors.append(abs(via-en['energies'][index])/max(en['total'],1e-30))
        global_fields=mechanics_global(m,case,coef,w,cfg['shape_plot_nodes'])
        coordinate_error=np.linalg.norm(global_fields[0]['displacement'][-1]-global_fields[1]['displacement'][-1])/scale
        energy_error=abs(2*en['total']/(w*w*en['mass'])-1);sum_error=abs(sum(en['fractions'])-1)
        gates=dict(physical=physical<=cfg['policy']['boundary_tol'],ode=ode<=cfg['policy']['ode_tol'],energy=energy_error<=cfg['policy']['energy_tol'],
            mass=mass_error<=cfg['mass_normalization_tol'],quadrature=fraction_spread<=cfg['quadrature_fraction_tol'] and mass_spread<=cfg['mass_normalization_tol'],
            reconstruction=reconstruction<=cfg['reconstruction_relative_tol'],fractions=sum_error<=cfg['fraction_sum_tol'],
            trace=max(trace_errors)<=cfg['trace_energy_relative_tol'],coordinates=coordinate_error<=cfg['coordinate_relative_tol'])
        if not all(gates.values()):raise ArithmeticError(f'Mechanical form gate failed {label}/{b}: {gates}')
        basename=f'{label}_beta{b:g}';arrays={'coefficients':coef,'omega':np.array(w)}
        for i,(a,g) in enumerate(zip(arms,global_fields),1):
            arrays.update({f'arm{i}_{k}':v for k,v in a.items()});arrays.update({f'plot_arm{i}_{k}':v for k,v in g.items()})
            plotmax=max(plotmax,float(np.max(np.linalg.norm(g['displacement'],axis=1))))
        np.savez_compressed(folder/'forms'/(basename+'.npz'),**arrays)
        r=dict(ID=label,kappa_joint=case.k_joint,kappa_end=case.k_end,beta_deg=b,beta_pl_deg=cfg['transition_angles'][label],
            omega1=w,omega2=rec['roots'][1]['omega'],omega3=rec['roots'][2]['omega'],Lambda1=root['Lambda'],gap12_relative=rec['roots'][1]['omega']/w-1,
            q=float(profile[label,b]['q']),Q_tail=float(profile[label,b]['Q_tail']),S0=float(profile[label,b]['S0']),
            mass=en['mass'],V=en['total'],delta_psi=en['delta_psi'],psi_end=en['psi_end'],D_joint=en['fractions'][3]/2,D_end=en['fractions'][4]/2,
            sign=sign,normalization_factor=factor,form_file=f'forms/{basename}.npz')
        for j,name in enumerate(components):r['V_'+name]=en['energies'][j];r['eta_'+name]=en['fractions'][j]
        for i in range(2):
            for j,name in enumerate(components[:3]):r[f'arm{i+1}_V_{name}']=en['arm_energy'][i,j]
            for j,name in enumerate(massnames):r[f'arm{i+1}_mass_{name}']=en['mass_components'][i,j]
        rows.append(r);checks.append(dict(ID=label,beta_deg=b,physical=float(physical),ode=float(ode),energy=energy_error,mass_error=mass_error,
            quadrature_fraction_spread=fraction_spread,quadrature_mass_spread=mass_spread,reconstruction=reconstruction,fraction_sum_error=sum_error,
            trace_energy_error=max(trace_errors),global_joint_error=float(coordinate_error),archived_singular=root['diagnostics']['singular_ratio'],
            archived_nonzero_condition=root['diagnostics']['nonzero_condition'],status='PASS'))
        provenance.append(dict(ID=label,beta_deg=b,point_file=source['point_file'],point_sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz')),frequency_origin='archived_sorted_position_1'))
    interval=[]
    for label in cfg['cases']:
        table={r['beta_deg']:r for r in rows if r['ID']==label}
        for a,b in zip(cfg['base_angles'][:-1],cfg['base_angles'][1:]):
            for name in components:
                delta=100*(table[b]['eta_'+name]-table[a]['eta_'+name])
                interval.append(dict(ID=label,beta_from=a,beta_to=b,component=name,delta_percentage_points=delta,mean_percentage_points_per_degree=delta/(b-a)))
    for filename,data in [('energies.csv',rows),('checks.csv',checks),('quadrature.csv',quads),('source_nodes.csv',provenance),('energy_intervals.csv',interval),('q_profiles.csv',profiles)]:csvwrite(folder/filename,data)
    write(folder/'display.json',dict(mass_normalization=1,displacement_multiplier=cfg['shape_max_display_displacement']*m.l/plotmax,
        maximum_plotted_displacement_over_l=cfg['shape_max_display_displacement'],sign_anchor=cfg['shape_sign'],source='same multiplier for all 29 mass-normalized forms'))
    return rows,checks


def mechanics_fd(run,cfg,rows):
    """Eight and only eight local prefixes; derivative check, not an angular map."""
    catalogs=json.loads((ROOT/cfg['catalog_bundle']/'pole_catalogs.json').read_text());run.catalogs={float(k):v for k,v in catalogs.items()}
    derivatives=[];h=cfg['fd_log_step']
    for label in cfg['fd_case_ids']:
        baseline=next(r for r in rows if r['ID']==label and r['beta_deg']==cfg['fd_beta']);pair=cfg['cases'][label]
        for axis in (0,1):
            found=[]
            for sign in (-1,1):
                k=list(pair);k[axis]*=math.exp(sign*h);case=physics.Case(cfg['fd_beta'],cfg['mu'],*k)
                if run.new_points+run.reused_points>=cfg['max_new_spectral_points']:raise ArithmeticError('Local derivative budget exhausted')
                rec=run.solve(case,neighbour=[baseline[f'omega{i}'] for i in (1,2,3)]);found.append(rec['roots'][0]['omega'])
            derivative=(math.log(found[1])-math.log(found[0]))/(2*h);expected=baseline['D_joint' if axis==0 else 'D_end'];error=abs(derivative-expected)
            if error>cfg['fd_absolute_tol']+cfg['fd_relative_tol']*abs(expected):raise ArithmeticError('Log-stiffness derivative check failed')
            derivatives.append(dict(ID=label,beta_deg=cfg['fd_beta'],axis='joint' if axis==0 else 'end',log_step=h,
                omega_minus=found[0],omega_plus=found[1],D_energy=expected,D_finite_difference=derivative,absolute_error=error,relative_error=error/abs(expected),status='PASS'))
    csvwrite(run.folder/'derivative_checks.csv',derivatives)


def plot_mechanics(folder):
    """Four figure sets from saved scalar tables and global coordinates only."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    folder=Path(folder);cfg=json.loads((folder/'contract.json').read_text(encoding='utf8'))['config']
    if (folder/'manifest.json').exists():
        manifest=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        names=['contract.json','energies.csv','q_profiles.csv','display.json']+[f for f in manifest['artifacts'] if f.startswith('forms/')]
        if any(sha(folder/f)!=manifest['artifacts'][f] for f in names):raise ValueError('Mechanical plot inputs changed')
    read=lambda n:list(csv.DictReader((folder/n).open(encoding='utf8')))
    rows=read('energies.csv');profiles=read('q_profiles.csv');display=json.loads((folder/'display.json').read_text(encoding='utf8'))
    components=['axial','bend','shear','joint','end'];colors=['#d98c24','#3177a8','#8c8c8c','#23996d','#8b5aa8'];labels=['Продольная','Изгиб','Сдвиг','Внутренняя пружина','Внешняя пружина']
    title=lambda k:rf'$\kappa_{{joint}}={cfg["cases"][k][0]:g},\ \kappa_{{end}}={cfg["cases"][k][1]:g}$'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    dest=folder/'figures';dest.mkdir(exist_ok=True)
    def save(fig,name):
        fig.savefig(dest/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(dest/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    fig,axes=plt.subplots(1,5,figsize=(13,3.5),sharey=True)
    for ax,label in zip(axes,cfg['cases']):
        selected=[r for r in rows if r['ID']==label and float(r['beta_deg']) in cfg['base_angles']];bottom=np.zeros(5)
        for name,color,legend in zip(components,colors,labels):
            values=100*np.array([float(r['eta_'+name]) for r in selected]);ax.bar(range(5),values,bottom=bottom,color=color,label=legend,width=.72);bottom+=values
        ax.set(xticks=range(5),xticklabels=cfg['base_angles'],xlabel='β, °',title=title(label),ylim=(0,100))
    axes[0].set_ylabel('Доля потенциальной энергии, %');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=5,frameon=False);fig.subplots_adjust(bottom=.27,wspace=.14);save(fig,'01_energy_partition')
    fig,axes=plt.subplots(1,5,figsize=(13,3.5),sharey=True)
    for ax,label in zip(axes,cfg['cases']):
        selected=[r for r in rows if r['ID']==label];beta=[float(r['beta_deg']) for r in selected]
        for name,color,style in [('joint',colors[3],'-o'),('end',colors[4],'--s')]:ax.plot(beta,[100*float(r['eta_'+name]) for r in selected],style,color=color,ms=3,lw=1.2,label='Внутренняя' if name=='joint' else 'Внешняя')
        ax.axvline(cfg['transition_angles'][label],color='.55',lw=.8,ls=':');ax.set(xlabel='β, °',title=title(label),xticks=[0,6,10,20]);ax.grid(alpha=.15)
    axes[0].set_ylabel('Энергия пружины / полная энергия, %');fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False);fig.subplots_adjust(bottom=.27,wspace=.14);save(fig,'02_spring_fractions')
    fig,axes=plt.subplots(len(cfg['shape_case_ids']),3,figsize=(11.5,7.9),sharex=True,sharey=True);factor=display['displacement_multiplier']
    for i,label in enumerate(cfg['shape_case_ids']):
        for j,b in enumerate((0,6,20)):
            ax=axes[i,j];r=next(r for r in rows if r['ID']==label and float(r['beta_deg'])==b)
            with np.load(folder/r['form_file']) as data:
                for arm in (1,2):
                    ref=data[f'plot_arm{arm}_reference'];disp=data[f'plot_arm{arm}_displacement'];axial=data[f'plot_arm{arm}_axial']
                    ax.plot(*ref.T,color='.65',lw=1);ax.plot(*(ref+factor*axial).T,color='#23996d',ls='--',lw=.9);ax.plot(*(ref+factor*disp).T,color='#3177a8',lw=1.4)
                    ax.plot(*ref[0],marker='s',ms=3,color='.2')
            ax.set_aspect('equal');ax.set(xlim=(-1.12,1.12),ylim=(-.25,.62));ax.grid(alpha=.12)
            if i==0:ax.set_title(f'β = {b}°')
            if j==0:ax.set_ylabel(title(label)+'\nY/l')
            if i==len(axes)-1:ax.set_xlabel('X/l')
    handles=[Line2D([],[],color='.65',label='Исходная ось'),Line2D([],[],color='#3177a8',label='Полное перемещение'),Line2D([],[],color='#23996d',ls='--',label='Только продольная составляющая')]
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.045));fig.text(.5,.02,f'Массовая норма M = 1; единый множитель перемещений ×{factor:.6g}; знак: dY первого плеча в середине > 0',ha='center',fontsize=8)
    fig.subplots_adjust(left=.12,right=.98,top=.95,bottom=.16,hspace=.15,wspace=.13);save(fig,'03_global_shapes')
    fig,axes=plt.subplots(2,3,figsize=(11,6),sharex=True)
    for j,label in enumerate(cfg['contrast_case_ids']):
        pp=[r for r in profiles if r['ID']==label and float(r['beta_deg'])<=20];rr=[r for r in rows if r['ID']==label]
        axes[0,j].plot([float(r['beta_deg']) for r in pp],[float(r['q']) for r in pp],color='.15',lw=1.2);axes[0,j].axhline(.15,color='#ba3d3d',ls='--',lw=1);axes[0,j].set(title=title(label),ylabel='q',ylim=(0,2.05))
        for name,color,legend in zip(components[:2],colors[:2],labels[:2]):axes[1,j].plot([float(r['beta_deg']) for r in rr],[100*float(r['eta_'+name]) for r in rr],'-o',color=color,lw=1.2,ms=3,label=legend)
        for ax in axes[:,j]:ax.axvline(cfg['transition_angles'][label],color='.45',ls=':',lw=1);ax.set(xticks=[0,3,6,10,20],xlim=(0,20));ax.grid(alpha=.15)
        axes[1,j].set(xlabel='β, °',ylabel='Доля энергии, %',ylim=(0,100))
    axes[1,1].legend(loc='center right',fontsize=8);fig.text(.5,.02,'Пунктир по вертикали: рабочий β_pl; энергетические линии соединяют только исследованные точки',ha='center',fontsize=8);fig.subplots_adjust(bottom=.12,wspace=.3,hspace=.25);save(fig,'04_slope_and_energy')


def mechanics_main(args):
    cfg=json.loads(MECHANICS_CONFIG.read_text(encoding='utf8'));key,ident,selected,profiles=mechanics_sources(cfg);folder=OUTPUT/'mechanical_interpretation'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        p=(args.verify_cache or folder).resolve();manifest=json.loads((p/'manifest.json').read_text(encoding='utf8'))
        if manifest['identity']!=ident or any(sha(p/f)!=h for f,h in manifest['artifacts'].items()):raise ValueError('Mechanical cache mismatch')
        print(f'REUSED {p}: 0 eigenvalue solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident);run=Run(cfg,folder);start=time.perf_counter()
    rows,checks=mechanics_analysis(folder,cfg,selected,profiles);analysis_time=time.perf_counter()-start;startfd=time.perf_counter()
    try:mechanics_fd(run,cfg,rows)
    except Exception as exc:
        write(folder/'failure.json',dict(error=repr(exc),new_points=run.new_points,new_roots=run.new_roots));raise
    fd_time=time.perf_counter()-startfd
    if run.new_points+run.reused_points!=8 or run.pole_roots or run.ritz_solves:raise ArithmeticError('Mechanical solve accounting mismatch')
    startplot=time.perf_counter();plot_mechanics(folder);render_time=time.perf_counter()-startplot
    _,after,_,_=mechanics_sources(cfg)
    if after!=ident:raise ValueError('Mechanical sources changed during analysis')
    write(folder/'diagnostics.json',dict(status='PASS',forms=29,maximum={k:max(r[k] for r in checks) for k in checks[0] if k not in ('ID','beta_deg','status')},minimum_gap12_relative=min(r['gap12_relative'] for r in rows)))
    write(folder/'run_ledger.json',dict(command=sys.argv,archived_points=29,archived_root_values=87,base_frequency_solves=0,
        archived_forms=29,quadrature_reconstructions=87,new_nullspace_recoveries_for_base_forms=0,
        finite_difference_points=8,new_points=run.new_points,checkpoint_reuse_calls=run.reused_points,new_roots=run.new_roots,
        pole_roots=run.pole_roots,ritz_solves=run.ritz_solves,analysis_seconds=analysis_time,finite_difference_seconds=fd_time,render_seconds=render_time))
    write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        base_frequency_solves=0,finite_difference_points=8,finite_difference_root_values=24,archived_points=29,complex_solves=0,
        artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def variation_constraints(beta_rad):
    """Endpoint order (u1,w1,u2,w2); fixed orthonormal basis (t1,n1)."""
    c,s=math.cos(beta_rad),math.sin(beta_rad)
    return np.array([[1.,0.,c,s],[0.,1.,-s,c]]),np.array([[0.,0.,-s,c],[0.,0.,-c,-s]])


def variation_work(reaction, geometric, omega, mass):
    """Signed log-frequency derivative PER RADIAN, for Kz+C.T r=lambda Mz."""
    work=float(reaction@geometric);rn=float(np.linalg.norm(reaction));gn=float(np.linalg.norm(geometric))
    return dict(W_beta=work,reaction_norm=rn,geometry_norm=gn,
        cos_gamma=work/(rn*gn) if rn*gn else None,s_rad=work/(omega*omega*mass))


def variation_fd(table,b,h):
    """Actual archived points only; endpoint one-sided stencils are explicit."""
    if b-h in table and b+h in table:points=[b-h,b,b+h];weights=np.array([-1.,0.,1.])/(2*h);scheme='central'
    elif b in (0,3) and b+2*h in table and b+h in table:points=[b,b+h,b+2*h];weights=np.array([-3.,4.,-1.])/(2*h);scheme='forward'
    elif b==12 and b-2*h in table and b-h in table:points=[b-2*h,b-h,b];weights=np.array([1.,-4.,3.])/(2*h);scheme='backward'
    else:return None
    return dict(h_deg=h,scheme=scheme,angles_deg=';'.join(map(str,points)),
                s_FD_deg=float(weights@np.log([float(table[x]['omega1']) for x in points])))


def variation_sources(cfg):
    _,ident=identity(cfg);ident['files'][VARIATION_CONFIG.relative_to(ROOT).as_posix()]=sha(VARIATION_CONFIG)
    for name in ('source','mechanics'):
        folder=ROOT/cfg[name+'_bundle'];manifest=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if sha(folder/'manifest.json')!=cfg[name+'_manifest_sha256'] or manifest['status']!='PASS':raise ValueError('Variation source manifest changed')
        if any(sha(folder/f)!=h for f,h in manifest['artifacts'].items()):raise ValueError('Variation source artifacts changed')
        for k in ('material','mu','policy'):
            if cfg[k]!=manifest['identity']['config'][k]:raise ValueError('Variation physical contract mismatch')
        for f,h in manifest['identity']['files'].items():
            if f!=Path(__file__).relative_to(ROOT).as_posix() and sha(ROOT/f)!=h:raise ValueError('Variation source implementation changed: '+f)
    source=ROOT/cfg['source_bundle'];read=lambda f:list(csv.DictReader((source/f).open(encoding='utf8')))
    spectra=[r for r in read('spectrum.csv') if r['ID'] in cfg['cases']]
    profiles=[r for r in read('profiles.csv') if r['ID'] in cfg['cases']]
    transitions=[r for r in read('transitions.csv') if r['ID'] in cfg['cases'] and float(r['q_star'])==.15]
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    return key,ident,spectra,profiles,transitions


def variation_analysis(folder,cfg,spectra,profiles,transitions):
    """Endpoint reactions from the accepted weak form; never call a root solver."""
    m=physics.Model(**cfg['material']);sec=m.section;rows=[];checks=[];differences=[];sources=[];gaps=[]
    previous={(r['ID'],float(r['beta_deg'])):r for r in csv.DictReader((ROOT/cfg['mechanics_bundle']/'energies.csv').open(encoding='utf8'))}
    prof={(r['ID'],float(r['beta_deg'])):r for r in profiles};transition={r['ID']:float(r['beta_pl_deg']) for r in transitions}
    for label,pair in cfg['cases'].items():
        table={float(r['beta_deg']):r for r in spectra if r['ID']==label}
        for b,r in sorted(table.items()):
            if b<=12:gaps.append(dict(ID=label,beta_deg=b,omega1=r['omega1'],omega2=r['omega2'],omega3=r['omega3'],gap12=float(r['omega2'])-float(r['omega1']),gap12_relative=float(r['omega2'])/float(r['omega1'])-1))
        for b in cfg['angles_deg']:
            source=table[b];path=ROOT/source['point_file'];rec=json.loads(path.read_text(encoding='utf8'))
            if sha(path)!=source['point_sha256'] or sha(path.with_suffix('.npz'))!=source['shape_sha256']:raise ValueError('Variation point changed')
            expected=(b,cfg['mu'],*pair,'RL_FULL');actual=tuple(rec['case'][k] for k in ('beta','mu','k_joint','k_end','variant'))
            if actual!=expected or rec['status']!='PASS' or rec['count']['status']!='PASS':raise ValueError('Variation point contract failed')
            for i,r in enumerate(rec['roots'],1):
                if r['counts']!=[i-1,i] or not all(r['gates'].values()) or r['omega']!=float(source[f'omega{i}']):raise ValueError('Variation prefix failed')
            case=physics.Case(b,cfg['mu'],*pair);root=rec['roots'][0];w=root['omega'];coef=np.array(root['coefficients'])
            with np.load(path.with_suffix('.npz')) as data:arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
            energy=mechanics_energy(m,case,arms);mass=energy['mass'];rebuilt=mechanics_fields(m,case,w,coef,cfg['quadrature_order'])
            reconstruction=max(np.max(abs(a[f]-d[f]))/max(np.max(abs(a[f])),1e-30) for a,d in zip(arms,rebuilt) for f in ('y','dy','ends'))
            ends=np.array([a['ends'][1] for a in arms]);z=ends[:,:2].ravel();rad=math.radians(b);C,Cb=variation_constraints(rad)
            reaction=-ends[0,3:5];geometric=Cb@z;work=variation_work(reaction,geometric,w,mass);sdeg=work['s_rad']*math.pi/180
            qscale=max(np.max(abs(a['y'][:,:3]*[1,1,m.l])) for a in arms);fscale=max(np.max(abs(a['y'][:,3:]*[1,1,1/m.l])) for a in arms)
            constraint=float(max(abs(C@z))/qscale);reaction_residual=max(abs(ends[:,3:5].ravel()+C.T@reaction))/fscale
            # Independent integration against admissible outer-end-zero test fields x/L.
            weak=[]
            for a,L in zip(arms,case.lengths(m)):
                y,d,wt,x=a['y'],a['dy'],a['weights'],a['x']
                weak.append([wt@((sec.EA-w*w*m.J)*d[:,0]/L-w*w*sec.rhoA*y[:,0]*x/L),
                             wt@(sec.KGA*(d[:,1]-y[:,2])/L-w*w*sec.rhoA*y[:,1]*x/L)])
            weak_error=float(np.max(abs(np.array(weak)-ends[:,3:5]))/fscale)
            state=np.concatenate(ends);B=np.diag([1.,-1.]);actual_C=physics.joint_matrix(m,case)[:2][:,[0,1,6,7]]
            operator_error=float(np.max(abs(C-B@actual_C)))
            h=cfg['operator_step_rad'];Cp=variation_constraints(rad+h)[0];Cm=variation_constraints(rad-h)[0]
            operator_fd=float(np.max(abs((Cp-Cm)/(2*h)-Cb)))
            # Existing frequency-dependent DtN reduction, evaluated at the SAVED omega.
            # Its beta derivative is a trace-functional check, not a new eigenproblem.
            P,Ks=physics.trace_map(m,case);D=block_diag(*(physics.arm_dtn(m,L,w)[0] for L in case.lengths(m)))
            Pb=np.zeros_like(P);Pb[9,:2]=[math.sin(rad),-math.cos(rad)];Pb[10,:2]=[-math.cos(rad),-math.sin(rad)]
            node=ends[0,:2]@physics.geometry.frames(b)[0].translation.T;t=np.r_[node,ends[0,2],ends[1,2],arms[1]['ends'][0,2]]
            Wred=float(t@(Pb.T@D@P+P.T@D@Pb)@t/2)
            Pp,Kp=physics.trace_map(m,replace(case,beta=b+math.degrees(h)));Pm,Km=physics.trace_map(m,replace(case,beta=b-math.degrees(h)))
            Wfd=float(t@((Pp.T@D@Pp+Kp)-(Pm.T@D@Pm+Km))@t/(4*h))
            denominator=w*w*mass;reduced_error=abs(Wred-work['W_beta'])/denominator;reduced_fd=abs(Wfd-work['W_beta'])/denominator
            trace_residual=float(np.linalg.norm((P.T@D@P+Ks)@t)/max(fscale,1e-30))
            scale_error=max(abs(variation_work(a*reaction,a*geometric,w,a*a*mass)['s_rad']-work['s_rad']) for a in (-3.,.2))
            # A fixed orthogonal basis change leaves both work and factor decomposition invariant.
            O=np.array([[.6,-.8],[.8,.6]]);rotated=variation_work(O@reaction,O@geometric,w,mass)
            rotation_error=abs(rotated['s_rad']-work['s_rad'])
            reused=(label,b) in previous;energy_reuse=0.
            if reused:energy_reuse=max(abs(energy['fractions'][i]-float(previous[label,b]['eta_'+name])) for i,name in enumerate(('axial','bend','shear','joint','end')))
            diag=root['diagnostics'];physical=max(max(diag['joint_scaled']),np.max(diag['external_scaled']));energy_error=abs(2*energy['total']/denominator-1)
            gate=dict(constraint=constraint<=cfg['policy']['boundary_tol'],reaction=reaction_residual<=cfg['policy']['boundary_tol'],weak=weak_error<=cfg['weak_reaction_tolerance'],
                operator=max(operator_error,operator_fd)<=cfg['operator_tolerance'],reduced=max(reduced_error,reduced_fd,trace_residual)<=cfg['operator_tolerance'],
                invariance=max(scale_error,rotation_error)<=cfg['operator_tolerance'],mass=abs(mass-1)<=cfg['mass_tolerance'],energy=energy_error<=cfg['policy']['energy_tol'],
                physical=physical<=cfg['policy']['boundary_tol'],ode=diag['ode_error']<=cfg['policy']['ode_tol'],reconstruction=reconstruction<=cfg['reconstruction_tolerance'],energy_reuse=energy_reuse<=cfg['energy_reuse_tolerance'])
            if not all(gate.values()):raise ArithmeticError(f'Variation gate {label}/{b}: {gate}')
            fd=[]
            for halfwidth in cfg['fd_halfwidths_deg']:
                diff=variation_fd(table,b,halfwidth)
                if diff is not None:
                    diff.update(ID=label,beta_deg=b,s_var_deg=sdeg,absolute_error=abs(diff['s_FD_deg']-sdeg),relative_error=abs(diff['s_FD_deg']-sdeg)/abs(sdeg) if abs(sdeg)>1e-12 else None);fd.append(diff)
            primary=fd[0];fd_status='EARLY_GRID_LIMITED'
            if b>=3:
                coarse=next(d for d in fd if d['h_deg']==.5);estimate=abs(coarse['s_FD_deg']-primary['s_FD_deg'])
                if primary['absolute_error']>estimate+1e-9 or primary['absolute_error']>coarse['absolute_error']:raise ArithmeticError('FD refinement does not confirm variation')
                fd_status='REFINEMENT_CONFIRMED'
            r=dict(ID=label,kappa_joint=pair[0],kappa_end=pair[1],beta_deg=b,beta_pl_deg=transition[label],omega1=w,Lambda1=root['Lambda'],omega2=rec['roots'][1]['omega'],omega3=rec['roots'][2]['omega'],
                gap12_relative=rec['roots'][1]['omega']/w-1,mass=mass,V=energy['total'],r1=reaction[0],r2=reaction[1],Cbeta_z1=geometric[0],Cbeta_z2=geometric[1],N1=ends[0,3],Q1=ends[0,4],u1=ends[0,0],w1=ends[0,1],
                **work,s_deg=sdeg,W_over_2V_rad=work['W_beta']/(2*energy['total']),S0=float(prof[label,b]['S0']),q_variational=abs(sdeg)/float(prof[label,b]['S0']),q_archived=float(prof[label,b]['q']),
                s_FD_primary=primary['s_FD_deg'],FD_h_deg=primary['h_deg'],FD_scheme=primary['scheme'],FD_status=fd_status,energy_origin='previous_mechanics' if reused else 'archived_fields_evaluated')
            for i,name in enumerate(('axial','bend','shear','joint','end')):r['eta_'+name]=float(previous[label,b]['eta_'+name]) if reused else energy['fractions'][i]
            rows.append(r);differences.extend(fd);checks.append(dict(ID=label,beta_deg=b,constraint=constraint,reaction_equilibrium=reaction_residual,weak_reaction=weak_error,operator=operator_error,operator_FD=operator_fd,
                reduced=reduced_error,reduced_FD=reduced_fd,trace_residual=trace_residual,scale_invariance=scale_error,basis_invariance=rotation_error,mass_error=abs(mass-1),energy_error=energy_error,physical=physical,ode=diag['ode_error'],reconstruction=reconstruction,energy_reuse=energy_reuse,status='PASS'))
            sources.append(dict(ID=label,beta_deg=b,point_file=source['point_file'],point_sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz')),source='archived_first_sorted_mode'))
    for name,data in [('variation.csv',rows),('checks.csv',checks),('finite_differences.csv',differences),('sources.csv',sources),('gap_profiles.csv',gaps)]:csvwrite(folder/name,data)
    return rows,checks


def variation_summary(folder):
    """Scalar diagnostics from saved rows, no new states or spectral samples."""
    rows=list(csv.DictReader((folder/'variation.csv').open(encoding='utf8')));gaps=list(csv.DictReader((folder/'gap_profiles.csv').open(encoding='utf8')))
    events=[];factors=[];components=[]
    for label in dict.fromkeys(r['ID'] for r in rows):
        rr=[r for r in rows if r['ID']==label];table={float(r['beta_deg']):r for r in rr};gg=[r for r in gaps if r['ID']==label]
        peak=max(rr,key=lambda r:float(r['s_deg']));peakW=max(rr,key=lambda r:float(r['W_beta']));gap=min(gg,key=lambda r:float(r['gap12_relative']));absolute=min(gg,key=lambda r:float(r['gap12']))
        slopes=[((float(b['s_deg'])-float(a['s_deg']))/(float(b['beta_deg'])-float(a['beta_deg'])),a['beta_deg'],b['beta_deg']) for a,b in zip(rr,rr[1:])];fast=min(slopes)
        events.append(dict(ID=label,sensitivity_sample_max_beta=float(peak['beta_deg']),W_sample_max_beta=float(peakW['beta_deg']),fastest_sample_decline_start=float(fast[1]),fastest_sample_decline_end=float(fast[2]),
            min_relative_gap_beta=float(gap['beta_deg']),min_relative_gap=float(gap['gap12_relative']),min_absolute_gap_beta=float(absolute['beta_deg']),min_absolute_gap=float(absolute['gap12']),beta_pl_deg=float(rr[0]['beta_pl_deg'])))
        for start,end in [(1.,3.),(2.,3.),(2.,10.),(3.,10.),(3.,8.)]:
            a,b=table[start],table[end];r=dict(ID=label,beta_from=start,beta_to=end)
            for k in ['reaction_norm','geometry_norm','cos_gamma','V','s_deg','W_beta']:r[k+'_ratio']=float(b[k])/float(a[k])
            r['s_ratio_reconstructed']=r['reaction_norm_ratio']*r['geometry_norm_ratio']*r['cos_gamma_ratio']/r['V_ratio']
            if abs(r['s_ratio_reconstructed']-r['s_deg_ratio'])>1e-10:raise ArithmeticError('Work factorization failed')
            factors.append(r)
        for r in rr:
            axial=float(r['N1'])*float(r['w1']);shear=-float(r['Q1'])*float(r['u1']);W=float(r['W_beta'])
            if abs(axial+shear-W)>1e-12:raise ArithmeticError('Fixed-frame work components failed')
            components.append(dict(ID=label,beta_deg=r['beta_deg'],N1_times_w1=axial,minus_Q1_times_u1=shear,W_beta=W,
                axial_term_fraction=axial/W if abs(W)>1e-12 else None,transverse_term_fraction=shear/W if abs(W)>1e-12 else None,
                geometry_norm_to_nodal_translation_ratio=float(r['geometry_norm'])/math.hypot(float(r['u1']),float(r['w1']))))
    for name,data in [('event_comparison.csv',events),('factor_changes.csv',factors),('work_components.csv',components)]:csvwrite(folder/name,data)


def plot_variation(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifest_path=folder/'manifest.json'
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text(encoding='utf8'))
        if any(sha(folder/f)!=h for f,h in manifest['artifacts'].items() if not f.startswith('figures/')):raise ValueError('Variation plot data changed')
    rows=list(csv.DictReader((folder/'variation.csv').open(encoding='utf8')));gaps=list(csv.DictReader((folder/'gap_profiles.csv').open(encoding='utf8')))
    labels=['J0E4','J4E0'];colors=['#3274a1','#e1812c'];titles=[r'$\kappa_{joint}=0.1,\ \kappa_{end}=10$',r'$\kappa_{joint}=10,\ \kappa_{end}=0.1$'];figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.16,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.savefig(figdir/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    def data(label):return [r for r in rows if r['ID']==label]
    def arr(rr,k):return np.array([float(r[k]) for r in rr])
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained',sharey=True)
    for ax,label,col,title in zip(axes,labels,colors,titles):
        rr=data(label);b=arr(rr,'beta_deg');ax.plot(b,arr(rr,'s_deg'),'o-',color=col,label='Вариационная');ax.plot(b,arr(rr,'s_FD_primary'),'s--',color='#555555',ms=4,label='Конечные разности');ax.axhline(0,color='gray',lw=.6);ax.set(title=title,xlabel='β, °',xlim=(0,12));ax.legend(fontsize=8)
    axes[0].set_ylabel('d ln ω₁ / dβ, 1/°');save(fig,'01_variation_and_FD')
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for label,col,title in zip(labels,colors,titles):
        rr=data(label);b=arr(rr,'beta_deg');axes[0].plot(b,arr(rr,'W_beta'),'o-',color=col,label=title);axes[1].plot(b,arr(rr,'W_over_2V_rad'),'o-',color=col,label=title)
    for ax in axes:ax.set(xlabel='β, °',xlim=(0,12));ax.legend(fontsize=8)
    axes[0].set_ylabel('Wβ при M=1');axes[1].set_ylabel('Wβ / (2V), 1/рад');save(fig,'02_virtual_work')
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for label,col,title in zip(labels,colors,titles):
        rr=data(label);b=arr(rr,'beta_deg')
        for ax,k,ylabel in zip(axes.flat,['reaction_norm','geometry_norm','cos_gamma','V'],['‖r‖ при M=1','‖Cβ z‖ при M=1','cos γ','2V при M=1']):ax.plot(b,arr(rr,k)*(2 if k=='V' else 1),'o-',color=col,label=title);ax.set(xlabel='β, °',ylabel=ylabel,xlim=(0,12))
    axes[0,0].legend(fontsize=8);save(fig,'03_work_factors')
    fig,axes=plt.subplots(3,2,figsize=(10,8),layout='constrained')
    for j,(label,col,title) in enumerate(zip(labels,colors,titles)):
        rr=data(label);b=arr(rr,'beta_deg');pl=float(rr[0]['beta_pl_deg']);gg=[r for r in gaps if r['ID']==label]
        axes[0,j].plot(b,arr(rr,'q_variational'),'o-',color=col);axes[0,j].axhline(.15,color='#aa4444',ls='--',lw=1);axes[0,j].set_title(title)
        axes[1,j].plot(b,100*arr(rr,'eta_axial'),'o-',color=col)
        axes[2,j].plot(arr(gg,'beta_deg'),arr(gg,'gap12_relative'),color=col)
        for ax in axes[:,j]:ax.axvline(pl,color='gray',ls=':',lw=1);ax.set(xlabel='β, °',xlim=(0,12))
    for ax,ylabel in zip(axes[:,0],['|s_var| / S₀','Продольная энергия, %','(ω₂ − ω₁) / ω₁']):ax.set_ylabel(ylabel)
    save(fig,'04_sensitivity_energy_gap')


def variation_main(args):
    cfg=json.loads(VARIATION_CONFIG.read_text(encoding='utf8'));key,ident,spectra,profiles,transitions=variation_sources(cfg);folder=OUTPUT/'angle_variation'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        p=(args.verify_cache or folder).resolve();manifest=json.loads((p/'manifest.json').read_text(encoding='utf8'))
        if manifest['identity']!=ident or any(sha(p/f)!=h for f,h in manifest['artifacts'].items()):raise ValueError('Variation cache mismatch')
        print(f'REUSED {p}: 0 eigenvalue solves, 0 field reconstructions');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident);start=time.perf_counter()
    try:rows,checks=variation_analysis(folder,cfg,spectra,profiles,transitions)
    except Exception as exc:write(folder/'failure.json',dict(error=repr(exc),new_spectral_solves=0));raise
    variation_summary(folder)
    analysis_time=time.perf_counter()-start;start=time.perf_counter();plot_variation(folder);plot_time=time.perf_counter()-start
    write(folder/'diagnostics.json',dict(status='PASS',maximum={k:max(r[k] for r in checks) for k in checks[0] if k not in ('ID','beta_deg','status')},fd_qualification='6 early points retain coarse-frequency-grid limitation; 14 points checked by refinement'))
    write(folder/'run_ledger.json',dict(command=sys.argv,archived_points=20,archived_frequency_values=60,archived_first_forms=20,new_spectral_solves=0,new_nullspace_recoveries=0,
        coefficient_field_evaluations=20,arm_DtN_linear_solves=40,energy_rows_reused=sum(r['energy_origin']=='previous_mechanics' for r in rows),energy_rows_evaluated=sum(r['energy_origin']!='previous_mechanics' for r in rows),analysis_seconds=analysis_time,render_seconds=plot_time))
    if variation_sources(cfg)[1]!=ident:raise ValueError('Variation inputs changed')
    write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),new_spectral_solves=0,
        artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}));print(folder,flush=True)


def epsilon_curvature(rows,windows,degrees):
    """Even log-frequency fits in radians, scaled design and zero intercept."""
    result=[]
    for ID in dict.fromkeys(r['ID'] for r in rows):
        rr=[r for r in rows if r['ID']==ID];w0=next(float(r['omega1']) for r in rr if float(r['beta_deg'])==0)
        for window in windows:
            fit=[r for r in rr if 0<float(r['beta_deg'])<=window];x=np.array([float(r['beta_deg'])/window for r in fit])**2;y=np.log([float(r['omega1'])/w0 for r in fit])
            for degree in degrees:
                n=degree//2
                if len(fit)<=n:continue
                A=np.column_stack([x**j for j in range(1,n+1)]);coef=np.linalg.lstsq(A,y,rcond=None)[0];residual=float(max(abs(A@coef-y)));a2=float(coef[0]/math.radians(window)**2)
                result.append(dict(ID=ID,epsilon=float(rr[0]['epsilon']),window_deg=window,degree=degree,positive_points=len(fit),a2_rad=a2,a2_epsilon_squared=a2*float(rr[0]['epsilon'])**2,max_log_residual=residual,relative_log_residual=residual/max(abs(y)),design_condition=float(np.linalg.cond(A))))
    return result


def epsilon_scalars(model,label,epsilon,rec,path):
    case=physics.Case(**{k:rec['case'][k] for k in ('beta','mu','k_joint','k_end')});root=rec['roots'][0];w=root['omega']
    with np.load(path.with_suffix('.npz')) as data:arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
    en=mechanics_energy(model,case,arms);a,b=[v['ends'][1] for v in arms];C,Cb=variation_constraints(math.radians(case.beta));z=np.r_[a[:2],b[:2]];reaction=-a[3:5];geom=Cb@z;work=variation_work(reaction,geom,w,en['mass']);WN=a[3]*a[1];WQ=-a[4]*a[0]
    if abs(WN+WQ-work['W_beta'])>1e-10*max(w*w*en['mass'],abs(work['W_beta'])):raise ArithmeticError('Epsilon virtual-work identity failed')
    if abs(en['mass']-1)>1e-10 or abs(2*en['total']/(w*w*en['mass'])-1)>5e-8:raise ArithmeticError('Epsilon energy/mass failed')
    r=dict(ID=f'E{epsilon:g}_{label}',pair=label,epsilon=epsilon,beta_deg=case.beta,kappa_joint=case.k_joint,kappa_end=case.k_end,
        omega1=w,omega2=rec['roots'][1]['omega'],omega3=rec['roots'][2]['omega'],Lambda1=root['Lambda'],Lambda2=rec['roots'][1]['Lambda'],Lambda3=rec['roots'][2]['Lambda'],
        gap12_relative=rec['roots'][1]['omega']/w-1,mass=en['mass'],V=en['total'],u1=a[0],w1=a[1],N1=a[3],Q1=a[4],W_N=WN,W_Q=WQ,**work,
        s_deg=work['s_rad']*math.pi/180,W_N_fraction=WN/work['W_beta'] if abs(work['W_beta'])>1e-12 else None,W_Q_fraction=WQ/work['W_beta'] if abs(work['W_beta'])>1e-12 else None,
        xi=math.sin(math.radians(case.beta))/epsilon,beta_rad_over_epsilon=math.radians(case.beta)/epsilon,
        point_file=path.relative_to(ROOT).as_posix(),point_sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz')))
    for i,name in enumerate(('axial','bend','shear','joint','end')):r['eta_'+name]=en['fractions'][i]
    return r


def epsilon_postprocess(folder,cfg,records):
    rows=[epsilon_scalars(physics.Model(**{**cfg['material'],'epsilon':e}),label,e,rec,path) for e,label,rec,path in records]
    fd=[];summary=[];similarity=[];factors=[];checks=[]
    for ID in dict.fromkeys(r['ID'] for r in rows):
        rr=sorted([r for r in rows if r['ID']==ID],key=lambda r:r['beta_deg']);table={r['beta_deg']:r for r in rr};w0=table[0]['omega1']
        for r in rr:
            r['F_omega']=r['omega1']/w0;r['G_omega']=r['F_omega']-1
            for h in (.125,.25,.5,1):
                b=r['beta_deg']
                if b-h in table and b+h in table:
                    val=(math.log(table[b+h]['omega1'])-math.log(table[b-h]['omega1']))/(2*h)
                    fd.append(dict(ID=ID,beta_deg=b,h_deg=h,s_FD_deg=val,s_var_deg=r['s_deg'],absolute_error=abs(val-r['s_deg']),relative_error=abs(val-r['s_deg'])/abs(r['s_deg']) if abs(r['s_deg'])>1e-12 else None))
        peak=max(rr,key=lambda r:r['s_deg']);gap=min(rr,key=lambda r:r['gap12_relative'])
        summary.append(dict(ID=ID,pair=rr[0]['pair'],epsilon=rr[0]['epsilon'],sample_peak_beta=peak['beta_deg'],sample_peak_s_deg=peak['s_deg'],gap_at_s_peak=peak['gap12_relative'],min_gap_beta=gap['beta_deg'],min_gap_relative=gap['gap12_relative'],F_1deg=table[1]['F_omega'],F_3deg=table[3]['F_omega'],F_10deg=table[10]['F_omega'],F_20deg=table[20]['F_omega'],growth_10_to_20=table[20]['omega1']/table[10]['omega1']-1))
        # Comparable scaled positions: beta=3/10 deg at epsilon=.005 and their actual available counterparts.
        a=table[3*rr[0]['epsilon']/.005];b=table[10*rr[0]['epsilon']/.005]
        factors.append(dict(ID=ID,beta_from=a['beta_deg'],beta_to=b['beta_deg'],reaction_ratio=b['reaction_norm']/a['reaction_norm'],node_motion_ratio=b['geometry_norm']/a['geometry_norm'],alignment_ratio=b['cos_gamma']/a['cos_gamma'],energy_ratio=b['V']/a['V'],s_ratio=b['s_deg']/a['s_deg']))
    for pair in cfg['cases']:
        tables=[{round(r['beta_deg']/e,8):r for r in rows if r['pair']==pair and r['epsilon']==e} for e in cfg['epsilon_values']]
        for key in sorted(set.intersection(*(set(t) for t in tables))):
            a,b,c=[t[key] for t in tables]
            similarity.append(dict(pair=pair,beta_rad_over_epsilon=a['beta_rad_over_epsilon'],beta_thin=a['beta_deg'],beta_reference=b['beta_deg'],beta_thick=c['beta_deg'],xi_thin=a['xi'],xi_reference=b['xi'],xi_thick=c['xi'],F_thin=a['F_omega'],F_reference=b['F_omega'],F_thick=c['F_omega'],relative_F_spread=(max(a['F_omega'],b['F_omega'],c['F_omega'])-min(a['F_omega'],b['F_omega'],c['F_omega']))/b['F_omega']))
    for e,label,rec,path in records:
        for k,r in enumerate(rec['roots'],1):
            if rec['status']!='PASS' or rec['count']['status']!='PASS' or r['counts']!=[k-1,k] or not all(r['gates'].values()):raise ArithmeticError('Epsilon archived/new gate failed')
            d=r['diagnostics'];checks.append(dict(epsilon=e,pair=label,beta_deg=rec['case']['beta'],position=k,boundary=max(np.max(d['joint_scaled']),np.max(d['external_scaled'])),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],condition=d['nonzero_condition'],mass_gram=rec['mass_gram_error'],status='PASS'))
    fits=epsilon_curvature(rows,cfg['curvature_windows_deg'],cfg['curvature_degrees'])
    for name,data in [('spectrum.csv',rows),('curvature.csv',fits),('finite_differences.csv',fd),('comparison.csv',summary),('similarity.csv',similarity),('mechanical_factors.csv',factors),('quality.csv',checks)]:csvwrite(folder/name,data)
    return rows,fits


def plot_epsilon(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if (folder/'manifest.json').exists():
        mm=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f)!=h for f,h in mm['artifacts'].items() if not f.startswith('figures/')):raise ValueError('Epsilon plot inputs changed')
    rows=list(csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')));figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.16,'axes.spines.top':False,'axes.spines.right':False})
    pairs=['J0E4','J4E0'];titles=[r'$\kappa_{joint}=0.1,\ \kappa_{end}=10$',r'$\kappa_{joint}=10,\ \kappa_{end}=0.1$'];eps=[.0025,.005,.01];cols=['#3274a1','#e1812c','#3a923a']
    def data(pair,e):return sorted([r for r in rows if r['pair']==pair and float(r['epsilon'])==e],key=lambda r:float(r['beta_deg']))
    def vals(rr,k):return np.array([float(r[k]) for r in rr])
    def save(fig,name):
        fig.savefig(figdir/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    for name,x,y,xlabel,ylabel,xmax in [('01_relative_frequency','beta_deg','F_omega','β, °','ω₁(β) / ω₁(0)',20),('02_variational_slope','beta_deg','s_deg','β, °','d ln ω₁ / dβ, 1/°',12),('03_geometric_scaling','xi','F_omega','ξ = sin β / ε','ω₁(β) / ω₁(0)',40)]:
        fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained',sharey=True)
        for ax,pair,title in zip(axes,pairs,titles):
            for e,col in zip(eps,cols):
                rr=data(pair,e);ax.plot(vals(rr,x),vals(rr,y),'.-',ms=3,color=col,label=f'ε = {e:g}')
            ax.set(title=title,xlabel=xlabel,xlim=(0,xmax));ax.legend(fontsize=8)
        axes[0].set_ylabel(ylabel);save(fig,name)
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for j,(pair,title) in enumerate(zip(pairs,titles)):
        for e,col in zip(eps,cols):
            chosen=[.5*e/.005,1*e/.005,3*e/.005,10*e/.005];rr=[r for r in data(pair,e) if float(r['beta_deg']) in chosen]
            axes[0,j].plot(vals(rr,'beta_deg'),vals(rr,'W_N'),'o-',color=col,label=f'WN, ε={e:g}');axes[0,j].plot(vals(rr,'beta_deg'),vals(rr,'W_beta'),'x--',color=col,label=f'Wβ, ε={e:g}')
            axes[1,j].plot(vals(rr,'beta_deg'),vals(rr,'W_Q'),'o-',color=col,label=f'ε={e:g}')
        axes[0,j].set_title(title);axes[0,j].legend(fontsize=7);axes[1,j].axhline(0,color='gray',lw=.7);axes[1,j].legend(fontsize=8)
        for ax in axes[:,j]:ax.set_xlabel('β, °')
    axes[0,0].set_ylabel('WN и Wβ при M=1');axes[1,0].set_ylabel('WQ при M=1');save(fig,'04_virtual_work_components')


def epsilon_main(args):
    cfg=json.loads(EPSILON_CONFIG.read_text(encoding='utf8'));_,ident,spectra,_,_=variation_sources(cfg)
    vf=ROOT/cfg['variation_bundle'];vm=json.loads((vf/'manifest.json').read_text(encoding='utf8'))
    if sha(vf/'manifest.json')!=cfg['variation_manifest_sha256'] or any(sha(vf/f)!=h for f,h in vm['artifacts'].items()):raise ValueError('Epsilon variation source changed')
    catalog_source=ROOT/cfg['catalog_bundle'];catalog_manifest=json.loads((catalog_source/'manifest.json').read_text())
    if sha(catalog_source/'manifest.json')!='cd83173000bc700b339156c2aafd94b990b4287ae518774b0fc9da2de8a18267' or sha(catalog_source/'pole_catalogs.json')!=catalog_manifest['artifacts']['pole_catalogs.json']:raise ValueError('Epsilon reference pole catalog changed')
    ident['files'][EPSILON_CONFIG.relative_to(ROOT).as_posix()]=sha(EPSILON_CONFIG);key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];folder=OUTPUT/'epsilon_scaling'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        saved=(args.verify_cache or folder).resolve();mm=json.loads((saved/'manifest.json').read_text(encoding='utf8'))
        if mm['identity']!=ident or any(sha(saved/f)!=h for f,h in mm['artifacts'].items()):raise ValueError('Epsilon cache changed')
        print(f'REUSED {saved}: 0 new spectral solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident)
    runs={};pointmap={};source={(r['ID'],float(r['beta_deg'])):r for r in spectra};section=[];execution=[];start=time.perf_counter()
    for e in cfg['epsilon_values']:
        local={**cfg,'material':{**cfg['material'],'epsilon':e}};dest=folder/f'eps_{e:g}';dest.mkdir(exist_ok=True);contract=dest/'material_contract.json'
        if contract.exists() and json.loads(contract.read_text())!=local['material']:raise ValueError('Epsilon folder physical identity mismatch')
        write(contract,local['material']);run=Run(local,dest);runs[e]=run
        catalog=dest/'pole_catalogs.json'
        if catalog.exists():run.catalogs={float(k):v for k,v in json.loads(catalog.read_text()).items()}
        elif e==.005:run.catalogs={float(k):v for k,v in json.loads((ROOT/cfg['catalog_bundle']/'pole_catalogs.json').read_text()).items()}
        m=run.model;s=m.section;h=math.sqrt(12)*m.l*e
        section.append(dict(epsilon=e,h=h,b=4*h,A=s.area,I=s.second_moment,Ip=17*s.second_moment,J_RL=m.J,time_ref=m.time,EA_l2_over_EI=s.EA*m.l*m.l/s.EI,k_for_kappa_01=.1*s.EI/m.l,k_for_kappa_10=10*s.EI/m.l))
    def get(e,label,b):
        key=(e,label,b)
        if key in pointmap:return pointmap[key]
        if e==.005 and (label,b) in source:
            s=source[label,b];path=ROOT/s['point_file'];rec=json.loads(path.read_text())
            if sha(path)!=s['point_sha256'] or sha(path.with_suffix('.npz'))!=s['shape_sha256']:raise ValueError('Epsilon reuse mismatch')
        else:
            run=runs[e];case=physics.Case(b,cfg['mu'],*cfg['cases'][label]);near=[r for (ee,ll,bb),(r,pp) in pointmap.items() if ee==e and ll==label]
            rec=run.solve(case,neighbour=[r['omega'] for r in near[-1]['roots']] if near else None);name=hashlib.sha256(json.dumps(rec['case'],sort_keys=True).encode()).hexdigest()[:16];path=run.folder/'points'/(name+'.json');write(run.folder/'pole_catalogs.json',{str(k):v for k,v in run.catalogs.items()})
        pointmap[key]=(rec,path);return rec,path
    def ledger(stage):
        entry=dict(stage=stage,new_points=sum(r.new_points for r in runs.values()),new_roots=sum(r.new_roots for r in runs.values()),checkpoint_reads=sum(r.reused_points for r in runs.values()),new_pole_roots=sum(r.pole_roots for r in runs.values()),elapsed_seconds=time.perf_counter()-start)
        p=folder/'executions.json';old=json.loads(p.read_text()) if p.exists() else [];write(p,old+[entry]);return entry
    try:
        for e in (.0025,.01):
            for label in cfg['cases']:
                for b in cfg['preflight_angles']:get(e,label,b)
        write(folder/'preflight.json',dict(status='PASS',points=40,physical_gates='unchanged',low_count=0,first_three_counts=[[0,1],[1,2],[2,3]]))
        if args.stage=='preflight':print(ledger('preflight'));print(folder);return
        for e in cfg['epsilon_values']:
            for label in cfg['cases']:
                for b in cfg['angles_deg']:get(e,label,b)
        records=lambda:[(e,label,rec,path) for (e,label,b),(rec,path) in sorted(pointmap.items())]
        rows,fits=epsilon_postprocess(folder,cfg,records());scalar_evaluations=len(rows);triggers=[]
        for label in cfg['cases']:
            ff=[r for r in fits if r['ID']==f'E0.0025_{label}' and r['degree']>=4 and r['window_deg']<=1.25];values=[r['a2_rad'] for r in ff];spread=(max(values)-min(values))/abs(np.median(values));trigger=spread>cfg['thin_refinement']['diagnostic_relative_spread']
            triggers.append(dict(pair=label,epsilon=.0025,initial_fit_relative_spread=spread,triggered=trigger,added_angles=cfg['thin_refinement']['allowed_added_angles'] if trigger else []))
            if trigger:
                for b in cfg['thin_refinement']['allowed_added_angles']:get(.0025,label,b)
        write(folder/'curvature_refinement.json',triggers)
        rows,fits=epsilon_postprocess(folder,cfg,records());scalar_evaluations+=len(rows);csvwrite(folder/'sections.csv',section)
        quadrature=[]
        for r in rows:
            e=r['epsilon'];b=r['beta_deg']
            if b not in [v*e/.005 for v in (.5,1,3,10)]:continue
            rec,path=pointmap[e,r['pair'],b];m=runs[e].model;case=physics.Case(b,cfg['mu'],*cfg['cases'][r['pair']]);root=rec['roots'][0]
            en=mechanics_energy(m,case,mechanics_fields(m,case,root['omega'],np.array(root['coefficients']),240))
            errors=dict(mass_error=abs(en['mass']-1),energy_error=abs(en['total']/r['V']-1),fraction_error=max(abs(en['fractions']-np.array([r['eta_'+n] for n in ('axial','bend','shear','joint','end')]))))
            if max(errors.values())>1e-10:raise ArithmeticError('Epsilon quadrature instability')
            quadrature.append(dict(ID=r['ID'],beta_deg=b,base_order=120,check_order=240,**errors,status='PASS'))
        csvwrite(folder/'quadrature.csv',quadrature)
        reuse=[dict(epsilon=e,pair=label,beta_deg=b,origin='archived' if e==.005 and (label,b) in source else 'new',point_file=path.relative_to(ROOT).as_posix(),point_sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz'))) for (e,label,b),(rec,path) in sorted(pointmap.items())];csvwrite(folder/'provenance.csv',reuse)
        # Adjacent jumps are retained with their already accepted count-one/gate evidence.
        audit=[]
        for ID in dict.fromkeys(r['ID'] for r in rows):
            rr=sorted([r for r in rows if r['ID']==ID],key=lambda r:r['beta_deg'])
            for a,b in zip(rr,rr[1:]):
                jump=max(abs(b[f'omega{k}']/a[f'omega{k}']-1) for k in (1,2))
                if jump>cfg['policy']['neighbour_relative_jump']:audit.append(dict(ID=ID,beta_from=a['beta_deg'],beta_to=b['beta_deg'],relative_jump=jump,resolution='accepted count-one intervals, low count zero, unchanged physical gates; genuine resolved slope'))
        csvwrite(folder/'neighbour_audit.csv',audit);entry=ledger('complete')
        nnew=sum(r['origin']=='new' for r in reuse);nold=len(reuse)-nnew
        write(folder/'run_ledger.json',dict(unique_new_points=nnew,new_root_values=3*nnew,archived_points=nold,archived_root_values=3*nold,new_pole_roots_total=sum(x['new_pole_roots'] for x in json.loads((folder/'executions.json').read_text())),first_mode_scalar_evaluations=scalar_evaluations,unique_first_mode_points=len(rows),saved_coefficient_quadrature_evaluations=len(quadrature),new_mode_reconstructions=3*nnew,ritz_solves=0,complex_solves=0,extra_spectral_repairs=0,neighbour_events=len(audit)))
        startplot=time.perf_counter();plot_epsilon(folder);write(folder/'render.json',dict(seconds=time.perf_counter()-startplot))
        write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}));print(folder,flush=True)
    except Exception as exc:ledger('failed');write(folder/'failure.json',dict(error=repr(exc)));raise


def slip_scalar(model, label, chi, rec, path):
    from scripts.lib import isotropic_two_springs_longitudinal_slip as slip
    c=slip.Case(**{k:rec['case'][k] for k in ('beta','mu','k_joint','k_end')},chi=chi)
    with np.load(path.with_suffix('.npz')) as data:
        arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
    en=mechanics_energy(model,c,arms);u,v=[a['ends'][1] for a in arms];C,Cb=variation_constraints(math.radians(c.beta));z=np.r_[u[:2],v[:2]];c1,c2=C@z;d1,d2=Cb@z
    k=chi*model.section.EA/model.l if chi!='RIGID' else None
    force=k*c1 if k is not None else -u[3];r2=-u[4];extra=k*c1*c1/2 if k is not None else 0.;V=en['total']+extra;w=rec['roots'][0]['omega']
    W1,W2=force*d1,r2*d2;energy_error=abs(2*V/(w*w*en['mass'])-1)
    if energy_error>5e-8 or abs(en['mass']-1)>1e-10:raise ArithmeticError('Slip energy/mass failed')
    if chi!='RIGID' and abs(force+u[3])>1e-9*max(abs(force),np.linalg.norm(u[3:5])):raise ArithmeticError('Slip force sign failed')
    check=(u[3]*u[1]-u[4]*u[0]+(u[4]*c1 if chi!='RIGID' else 0))
    if abs(W1+W2-check)>1e-9*w*w:raise ArithmeticError('Slip variational identity failed')
    row=dict(pair=label,chi=chi,beta_deg=c.beta,**{f'omega{i}':r['omega'] for i,r in enumerate(rec['roots'],1)},**{f'Lambda{i}':r['Lambda'] for i,r in enumerate(rec['roots'],1)},
        c1=c1,c2=c2,partial_c1=d1,partial_c2=d2,k_c1=force,r2=r2,W_parallel=W1,W_constraint=W2,W_total=W1+W2,V=V,mass=en['mass'],energy_error=energy_error,s_rad=(W1+W2)/(2*V),s_deg=(W1+W2)/(2*V)*math.pi/180,
        old_formula_s_deg=(u[3]*u[1]-u[4]*u[0])/(2*V)*math.pi/180,gap12_relative=rec['roots'][1]['omega']/w-1,point_file=path.relative_to(ROOT).as_posix(),point_sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz')))
    for name,value in zip(('axial','bend','shear','joint','end','parallel'),list(en['energies'])+[extra]):row['eta_'+name]=value/V
    return row


def slip_validate(folder, cfg, rows, run):
    """Targeted quadrature and triggered count checks, no additional roots."""
    slip=run.engine;m=run.model;checks=[];events=[];queried=set();scale=0.
    for r in rows:
        c=slip.Case(r['beta_deg'],cfg['mu'],*cfg['cases'][r['pair']],r['chi'])
        for factor in (-3.,.2):
            value=(factor*r['k_c1']*factor*r['partial_c1']+factor*r['r2']*factor*r['partial_c2'])/(2*factor**2*r['V'])
            scale=max(scale,abs(value-r['s_rad']))
        if r['chi'] in (0,.01,1,'RIGID') and r['beta_deg'] in (1,3,10):
            rec=json.loads((ROOT/r['point_file']).read_text());rt=rec['roots'][0]
            arms=mechanics_fields(m,c,rt['omega'],np.array(rt['coefficients']),240);en=mechanics_energy(m,c,arms)
            u,v=[a['ends'][1] for a in arms];C,_=variation_constraints(math.radians(c.beta));delta=(C@np.r_[u[:2],v[:2]])[0]
            extra=0 if c.chi=='RIGID' else c.chi*m.section.EA/m.l*delta**2/2
            err=abs((en['total']+extra)/r['V']-1);mass=abs(en['mass']-1)
            if max(err,mass)>1e-10:raise ArithmeticError('Slip quadrature check')
            checks.append(dict(pair=r['pair'],chi=c.chi,beta=c.beta,quadrature=240,energy_change=err,mass_error=mass))
    for pair in cfg['cases']:
        for chi in cfg['chi_values']:
            rr=sorted([r for r in rows if r['pair']==pair and r['chi']==chi],key=lambda r:r['beta_deg'])
            for a,b in zip(rr,rr[1:]):
                jump=max(abs(b[f'omega{k}']/a[f'omega{k}']-1) for k in (1,2))
                if jump<=cfg['policy']['neighbour_relative_jump']:continue
                for r in (a,b):
                    tag=(pair,chi,r['beta_deg'])
                    if tag in queried:continue
                    c=slip.Case(r['beta_deg'],cfg['mu'],*cfg['cases'][pair],chi)
                    for k in (1,2,3):
                        counts=[slip.count_at(m,c,r[f'omega{k}']*(1+t*1e-6),run.cats(m,c),cfg['policy'])[1] for t in (-1,1)]
                        if counts!=[k-1,k]:raise ArithmeticError('Slip triggered count check')
                    queried.add(tag)
                events.append(dict(pair=pair,chi=chi,beta_from=a['beta_deg'],beta_to=b['beta_deg'],jump=jump,status='RESOLVED_COUNT_ONE_RECHECK',additional_eigenvalue_solves=0))
        rr=[r for r in rows if r['pair']==pair and r['beta_deg']==3]
        if any(r['chi']==10. for r in rr):
            one,ten,rigid=[next(r for r in rr if r['chi']==c) for c in (1.,10.,'RIGID')]
            if any(not one[f'omega{k}']<=ten[f'omega{k}']<=rigid[f'omega{k}'] for k in (1,2,3)):raise ArithmeticError('Slip chi10 monotonicity')
    if scale>1e-10:raise ArithmeticError('Slip sign/normalization invariance')
    csvwrite(folder/'quadrature.csv',checks);csvwrite(folder/'neighbour_audit.csv',events)
    write(folder/'additional_validation.json',dict(status='PASS',normalization_sign_error_rad=scale,quadrature_evaluations=len(checks),neighbour_events=len(events),triggered_count_queries=len(queried)*6,additional_eigenvalue_solves=0,chi10_monotonicity='PASS',peak_status_note='chi=0 has no positive-growth peak; other maxima are sampled values, not continuous optimizers'))


def plot_slip(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    if (folder/'manifest.json').exists():
        mm=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/f)!=h for f,h in mm['artifacts'].items() if not f.startswith('figures/')):raise ValueError('Slip plot inputs changed')
    rows=list(csv.DictReader((folder/'spectrum.csv').open(encoding='utf8')));out=folder/'figures';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.18,'axes.spines.top':False,'axes.spines.right':False})
    chis=[0.,.0001,.001,.01,.1,1.,'RIGID'];colors=plt.cm.viridis(np.linspace(.05,.9,7));pairs=['J0E4','J4E0'];titles=[r'$\kappa_j=0.1,\ \kappa_e=10$',r'$\kappa_j=10,\ \kappa_e=0.1$']
    def get(pair,chi):return sorted([r for r in rows if r['pair']==pair and (r['chi']=='RIGID' if chi=='RIGID' else r['chi']!='RIGID' and float(r['chi'])==chi)],key=lambda r:float(r['beta_deg']))
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    for name,y,label in [('01_relative_frequency','F_omega','ω₁(β,χ) / ω₁(0,χ)'),('02_variational_slope','s_deg','d ln ω₁ / dβ, град⁻¹')]:
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for ax,pair,title in zip(axes,pairs,titles):
            for chi,col in zip(chis,colors):
                rr=get(pair,chi);ax.plot([float(r['beta_deg']) for r in rr],[float(r[y]) for r in rr],'.-',color=col,label='точная связь' if chi=='RIGID' else f'χ={chi:g}')
            ax.set(title=title,xlabel='β, °',ylabel=label,xlim=(0,20));ax.legend(fontsize=8)
        save(fig,name)
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for j,(pair,title) in enumerate(zip(pairs,titles)):
        for chi in (0.,.01,1.,'RIGID'):
            rr=[r for r in get(pair,chi) if float(r['beta_deg']) in (0,1,3,6,10,20)];col=colors[chis.index(chi)]
            for i,key in enumerate(('W_parallel','W_constraint')):axes[i,j].plot([float(r['beta_deg']) for r in rr],[float(r[key])/(2*float(r['V']))*math.pi/180 for r in rr],'o-',color=col,label='точная связь' if chi=='RIGID' else f'χ={chi:g}')
        axes[0,j].set_title(title)
        for i in (0,1):axes[i,j].set(xlabel='β, °',ylabel=('W∥' if i==0 else 'W₂')+' / (2V), град⁻¹');axes[i,j].legend(fontsize=8);axes[i,j].axhline(0,color='gray',lw=.6)
    save(fig,'03_work_decomposition')
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for ax,pair,title in zip(axes,pairs,titles):
        for b in (1,3,6):
            yy=[next(float(r['G_omega'])*100 for r in get(pair,c) if float(r['beta_deg'])==b) for c in chis]
            ax.plot(range(7),yy,'o-',label=f'β={b}°')
        ax.set_xticks(range(7),['0','10⁻⁴','10⁻³','10⁻²','10⁻¹','1','∞']);ax.set(title=title,xlabel='χ (категории; ∞ — точная связь)',ylabel='Относительный прирост ω₁, %');ax.legend()
    save(fig,'04_initial_gain')


def slip_main(args):
    from scripts.lib import isotropic_two_springs_longitudinal_slip as slip
    cfg=json.loads(SLIP_CONFIG.read_text());_,ident=identity(cfg)
    ident['files'].update({SLIP_CONFIG.relative_to(ROOT).as_posix():sha(SLIP_CONFIG),Path(slip.__file__).relative_to(ROOT).as_posix():sha(slip.__file__)})
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];folder=OUTPUT/'longitudinal_slip'/key
    source=ROOT/cfg['source_bundle'];sm=json.loads((source/'manifest.json').read_text())
    if sha(source/'manifest.json')!=cfg['source_manifest_sha256'] or any(sha(source/f)!=h for f,h in sm['artifacts'].items()):raise ValueError('Slip archive changed')
    growth=OUTPUT/'growth_stiffness_5x5/361e9a839ef5be30';gm=json.loads((growth/'manifest.json').read_text())
    if sha(growth/'manifest.json')!='edffda3081ab53afa6676b41368d7f585739fe466e7a1c448e134c24bb86ee42' or any(sha(growth/f)!=h for f,h in gm['artifacts'].items()):raise ValueError('Slip dense source changed')
    catfolder=OUTPUT/'51af749660af1df5';cm=json.loads((catfolder/'manifest.json').read_text())
    if sha(catfolder/'manifest.json')!='cd83173000bc700b339156c2aafd94b990b4287ae518774b0fc9da2de8a18267' or sha(catfolder/'pole_catalogs.json')!=cm['artifacts']['pole_catalogs.json']:raise ValueError('Slip pole catalog changed')
    if args.verify_cache or (folder/'manifest.json').exists():
        saved=args.verify_cache or folder;mm=json.loads((saved/'manifest.json').read_text())
        if mm['identity']!=ident or any(sha(saved/f)!=h for f,h in mm['artifacts'].items()):raise ValueError('Slip cache mismatch')
        print(f'REUSED {saved}: 0 spectral solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident);run=Run(cfg,folder,engine=slip)
    catalog=OUTPUT/'51af749660af1df5'/'pole_catalogs.json';run.catalogs={float(k):v for k,v in json.loads(catalog.read_text()).items()}
    old=list(csv.DictReader((source/'spectrum.csv').open()));old=[r for r in old if float(r['epsilon'])==.005]
    # Dense baseline nodes are accepted source points of the verified epsilon archive.
    growth=ROOT/'results/isotropic_two_springs_rl_timoshenko/growth_stiffness_5x5/361e9a839ef5be30'
    old+= [{**r,'pair':r['ID']} for r in csv.DictReader((growth/'spectrum.csv').open()) if r['ID'] in cfg['cases']]
    baseline={(r['pair'],float(r['beta_deg'])):r for r in old};points={}
    ritz=json.loads((folder/'preflight.json').read_text())['ritz'] if (folder/'preflight.json').exists() else []
    def get(pair,chi,b):
        key=(pair,chi,float(b))
        if key in points:return points[key]
        if chi=='RIGID' and (pair,float(b)) in baseline:
            ss=baseline[pair,float(b)];path=ROOT/ss['point_file'];rec=json.loads(path.read_text())
            if sha(path)!=ss['point_sha256'] or sha(path.with_suffix('.npz'))!=ss['shape_sha256']:raise ValueError('Slip rigid source mismatch')
        else:
            case=slip.Case(float(b),cfg['mu'],*cfg['cases'][pair],chi);rec=run.solve(case);name=hashlib.sha256(json.dumps(rec['case'],sort_keys=True).encode()).hexdigest()[:16];path=folder/'points'/(name+'.json')
        points[key]=(rec,path);return rec,path
    try:
        for pair in cfg['cases']:
            for chi in cfg['preflight_chi']:
                for b in cfg['angles_deg']:get(pair,chi,b)
                for b in (0,3):
                    rec,path=get(pair,chi,b)
                    for n in cfg['policy']['ritz_orders']:
                        if any(v['pair']==pair and v['chi']==chi and v['beta']==b and v['order']==n for v in ritz):continue
                        w,coef,diag=slip.ritz(run.model,slip.Case(b,cfg['mu'],*cfg['cases'][pair],chi),n);run.ritz_solves+=1
                        error=float(max(abs(w/np.array([r['omega'] for r in rec['roots']])-1)))
                        if error>cfg['policy']['ritz_frequency_tol']:raise ArithmeticError('Slip independent Ritz frequency gate')
                        ritz.append(dict(pair=pair,chi=chi,beta=b,order=n,relative_error=error,**diag))
        write(folder/'preflight.json',dict(status='PASS',points=78,ritz=ritz))
        if args.stage=='preflight':print(folder,run.new_points,'new points');return
        for pair in cfg['cases']:
            for chi in cfg['chi_values']:
                for b in cfg['angles_deg']:get(pair,chi,b)
        # Bounded independent FD verification near initial and post-peak angles.
        for pair in cfg['cases']:
            for chi in (0.,.01,1.):
                for b in cfg['derivative_check_angles']:
                    for h in cfg['derivative_check_h']:
                        for sign in (-1,1):get(pair,chi,b+sign*h)
        extra=[]
        for pair in cfg['cases']:
            one=get(pair,1.,3)[0]['roots'][0]['omega'];rigid=get(pair,'RIGID',3)[0]['roots'][0]['omega'];difference=1-one/rigid
            if difference>.01:get(pair,10.,3)
            extra.append(dict(pair=pair,beta=3,chi1_relative_difference=difference,chi10_triggered=difference>.01))
        write(folder/'chi10_trigger.json',extra)
        rows=[slip_scalar(run.model,p,c,rec,path) for (p,c,b),(rec,path) in points.items()];fd=[];summary=[];mono=[]
        for pair in cfg['cases']:
            for chi in cfg['chi_values']:
                rr=sorted([r for r in rows if r['pair']==pair and r['chi']==chi],key=lambda r:r['beta_deg']);table={r['beta_deg']:r for r in rr};w0=table[0]['omega1']
                for r in rr:
                    r['F_omega']=r['omega1']/w0;r['G_omega']=r['F_omega']-1
                    for h in (.125,.25,.5,1):
                        b=r['beta_deg']
                        if b-h in table and b+h in table:
                            val=(math.log(table[b+h]['omega1'])-math.log(table[b-h]['omega1']))/(2*h)
                            fd.append(dict(pair=pair,chi=chi,beta=b,h=h,s_var=r['s_deg'],s_FD=val,absolute_error=abs(val-r['s_deg']),relative_error=abs(val-r['s_deg'])/max(abs(r['s_deg']),1e-15)))
                peak=max(rr,key=lambda r:r['s_deg']);summary.append(dict(pair=pair,chi=chi,sample_peak_beta=peak['beta_deg'],sample_peak_s=peak['s_deg'],peak_resolution='boundary' if peak['beta_deg']==20 else 'discrete_sample_only',min_gap=min(r['gap12_relative'] for r in rr),**{f'F{b}':table[b]['F_omega'] for b in (1,3,6,10,20)}))
            for b in cfg['angles_deg']:
                rr=[next(r for r in rows if r['pair']==pair and r['chi']==c and r['beta_deg']==b) for c in cfg['chi_values']]
                for a,z in zip(rr,rr[1:]):
                    for k in (1,2,3):
                        change=z[f'omega{k}']/a[f'omega{k}']-1
                        if change < -1e-8:raise ArithmeticError('Slip spectral monotonicity failed')
                        mono.append(dict(pair=pair,beta=b,chi_from=a['chi'],chi_to=z['chi'],position=k,relative_change=change,status='PASS'))
        for r in rows:
            if r['chi']==10.:
                w0=next(x['omega1'] for x in rows if x['pair']==r['pair'] and x['chi']==1. and x['beta_deg']==0)
                # Only a check point: do not invent a chi=10 normalized curve.
                r['F_omega']=None;r['G_omega']=None
        csvwrite(folder/'spectrum.csv',rows);csvwrite(folder/'summary.csv',summary);csvwrite(folder/'finite_differences.csv',fd);csvwrite(folder/'monotonicity.csv',mono)
        quality=[];provenance=[]
        for (pair,chi,b),(rec,path) in points.items():
            origin='new' if path.is_relative_to(folder) else 'archived';provenance.append(dict(pair=pair,chi=chi,beta=b,origin=origin,path=path.relative_to(ROOT).as_posix(),sha256=sha(path),shape_sha256=sha(path.with_suffix('.npz'))))
            for k,r in enumerate(rec['roots'],1):
                if r['counts']!=[k-1,k] or not all(r['gates'].values()):raise ArithmeticError('Slip root gate failed')
                d=r['diagnostics'];quality.append(dict(pair=pair,chi=chi,beta=b,position=k,boundary=max(np.max(d['joint_scaled']),np.max(d['external_scaled'])),ode=d['ode_error'],energy=d['energy_error'],mass_gram=rec['mass_gram_error'],status='PASS'))
        csvwrite(folder/'quality.csv',quality);csvwrite(folder/'provenance.csv',provenance)
        slip_validate(folder,cfg,rows,run)
        nnew=sum(r['origin']=='new' for r in provenance);write(folder/'ledger.json',dict(new_unique_points=nnew,new_root_values=3*nnew,archive_points=len(points)-nnew,new_pole_roots=0,ritz_eigenproblems=len(ritz),ritz_new_this_execution=run.ritz_solves,new_mode_reconstructions=3*nnew,first_mode_postprocessing=len(rows),spectral_repairs=0,new_prefix_solves_this_execution=run.new_points,checkpoint_reads_this_execution=run.reused_points))
        plot_slip(folder);write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}));print(folder)
    except Exception as exc:write(folder/'failure.json',dict(error=repr(exc),new_points=run.new_points,ritz_solves=run.ritz_solves));raise


def thickness_slopes(beta, values):
    """Actual-node quadratic derivative; zero is fixed by reflection symmetry.

    At 24 degrees the weights use 23.5,24,26. At 90 they use 86,88,90.
    The wider diagnostic changes only nodes with both beta +/- 1 available.
    """
    b=np.asarray(beta,dtype=float);y=np.asarray(values,dtype=float)
    if b.shape!=y.shape or len(b)<3 or b[0]!=0 or np.any(np.diff(b)<=0):
        raise ValueError('Invalid thickness derivative grid')
    primary=np.gradient(y,b,edge_order=2);primary[0]=0.
    wider=primary.copy();lookup=dict(zip(b,y))
    for i,x in enumerate(b):
        if 1<=x<=23 and x-1 in lookup and x+1 in lookup:
            wider[i]=(lookup[x+1]-lookup[x-1])/2
    return primary,wider


def thickness_sources(cfg):
    """Pinned canonical exact-joint archives, with epsilon-specific provenance."""
    _,ident=identity(cfg);ident['files'][THICKNESS_CONFIG.relative_to(ROOT).as_posix()]=sha(THICKNESS_CONFIG)
    read=lambda p:list(csv.DictReader(p.open(encoding='utf8')))
    sources={}
    for name in ('old','epsilon'):
        folder=ROOT/cfg[name+'_bundle'];m=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if sha(folder/'manifest.json')!=cfg[name+'_manifest_sha256'] or m['status']!='PASS':
            raise ValueError('Thickness source manifest changed')
        if any(sha(folder/f)!=h for f,h in m['artifacts'].items()):raise ValueError('Thickness source artifact changed')
        for f,h in m['identity']['files'].items():
            if f!=Path(__file__).relative_to(ROOT).as_posix() and sha(ROOT/f)!=h:
                raise ValueError('Verified source implementation changed: '+f)
        previous=m['identity']['config']
        if previous['mu']!=cfg['mu'] or previous['policy']!=cfg['policy']:raise ValueError('Source policy mismatch')
        if {**previous['material'],'epsilon':.01}!=cfg['material']:raise ValueError('Material mismatch')
        sources[name]=read(folder/'spectrum.csv')
    ep=ROOT/cfg['epsilon_bundle']/'eps_0.01'
    if json.loads((ep/'material_contract.json').read_text())!=cfg['material']:raise ValueError('Epsilon archive contract mismatch')
    values=[10**x for x in (-1,-.5,0,.5,1)]
    if cfg['stiffness_values']!=values or cfg['cases']!={f'J{i}E{j}':[a,b] for i,a in enumerate(values) for j,b in enumerate(values)}:
        raise ValueError('Rounded stiffness grid')
    old=[];archives={}
    for name,rows in sources.items():
        for row in rows:
            if name=='epsilon' and (float(row['epsilon'])!=.01 or float(row['beta_deg']) not in cfg['map']['parameter_grid']):continue
            label=row['ID'] if name=='old' else row['pair'];b=float(row['beta_deg']);path=ROOT/row['point_file']
            rec=json.loads(path.read_text());case=rec['case']
            if (case['beta'],case['mu'],case['k_joint'],case['k_end'],case['variant'])!=(b,cfg['mu'],*cfg['cases'][label],'RL_FULL'):
                raise ValueError('Source case mismatch')
            if sha(path)!=row['point_sha256'] or sha(path.with_suffix('.npz'))!=rec['shape_sha256']:raise ValueError('Source point changed')
            if rec['status']!='PASS' or rec['count']['status']!='PASS' or len(rec['roots'])!=3:raise ValueError('Unaccepted source prefix')
            for k,r in enumerate(rec['roots'],1):
                if not all(r['gates'].values()) or any(float(row[f'{f}{k}'])!=r[f] for f in ('omega','Lambda')):raise ValueError('Source CSV mismatch')
            if name=='epsilon':archives[label,b]=(rec,path)
            else:old.append({**row,'beta_deg':b,**{f'{f}{k}':float(row[f'{f}{k}']) for f in ('omega','Lambda') for k in (1,2,3)}})
    if len(archives)!=cfg['expected_reused_points'] or len(old)!=2075:raise ValueError('Source accounting mismatch')
    transitions=[r for r in read(ROOT/cfg['old_bundle']/'transitions.csv') if float(r['q_star'])==.15]
    if len(transitions)!=25:raise ValueError('Missing old surface')
    key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    return key,ident,archives,old,transitions


def thickness_profiles(cfg,rows):
    profiles=[];transitions=[];error=0.
    for label in sorted({r['ID'] for r in rows}):
        table=sorted((r for r in rows if r['ID']==label),key=lambda r:r['beta_deg'])
        beta=[r['beta_deg'] for r in table];y=np.log([r['omega1'] for r in table]);z=np.log([r['Lambda1'] for r in table])
        slope,wide=thickness_slopes(beta,y);other,_=thickness_slopes(beta,z)
        for ref in cfg['reference_angles_deg']:
            s0=float((y[beta.index(ref)]-y[0])/ref);l0=(z[beta.index(ref)]-z[0])/ref
            if s0<=0:raise ArithmeticError('Nonpositive initial slope')
            q=abs(slope)/s0;error=max(error,float(max(abs(q-abs(other)/l0))))
            if error>cfg['numerical_checks']['q_identity_absolute']:raise ArithmeticError('Normalization identity failed')
            result,Q=refined_pick(beta,q,.15,cfg['minimum_length_deg'])
            check,_=refined_pick(beta,abs(wide)/s0,.15,cfg['minimum_length_deg'])
            a,b=result['beta_pl_deg'],check['beta_pl_deg'];delta=None if None in (a,b) else b-a
            tail=max(v for x,v in zip(beta,q) if x>=24)
            if tail>=.15 or (delta is not None and abs(delta)>.5):
                raise ArithmeticError(f'Local resolution/tail trigger {label}, reference {ref}: {delta}, {tail}; inspect before extension')
            transitions.append(dict(ID=label,reference_end_deg=ref,kappa_joint=cfg['cases'][label][0],kappa_end=cfg['cases'][label][1],
                **result,S0=s0,wide_beta_pl_deg=b,wide_minus_primary_deg=delta,late24_q_max=float(tail),
                Lambda_pl=None if a is None else table[beta.index(a)]['Lambda1']))
            for i,x in enumerate(beta):
                profiles.append(dict(ID=label,reference_end_deg=ref,beta_deg=x,signed_S=float(slope[i]),wide_signed_S=float(wide[i]),
                    q=float(q[i]),Q_tail=float(Q[i]),S0=s0,stencil='exact_symmetry' if x==0 else 'nonuniform_23.5_24_26' if x==24 else 'backward_86_88_90' if x==90 else 'central_actual_neighbours'))
    return profiles,transitions,error


def thickness_spectra(run,archives,preflight_only=False):
    cfg=run.config;beta=cfg['map']['parameter_grid']
    catalog=ROOT/cfg['epsilon_bundle']/'eps_0.01/pole_catalogs.json'
    run.catalogs={float(k):v for k,v in json.loads(catalog.read_text()).items()}
    write(run.folder/'pole_catalogs.json',run.catalogs)
    rows=[];checks=[];neighbours=[];auditfile=run.folder/'trigger_checks.json'
    audits=json.loads(auditfile.read_text()) if auditfile.exists() else {}
    pre=cfg['preflight_case_ids'];labels=pre+([] if preflight_only else [k for k in cfg['cases'] if k not in pre])
    for label in labels:
        previous=None
        for b in beta:
            case=physics.Case(b,cfg['mu'],*cfg['cases'][label]);contract={**asdict(case),'variant':'RL_FULL'}
            key=hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest()[:16]
            if (label,b) in archives:rec,path=archives[label,b];origin='archive_epsilon_scaling'
            else:
                rec=run.solve(case,neighbour=previous);path=run.folder/'points'/(key+'.json');origin='new'
            w=np.array([r['omega'] for r in rec['roots']])
            if np.any(np.diff(w)<=0) or rec['count']['status']!='PASS' or [r['counts'] for r in rec['roots']]!=[[0,1],[1,2],[2,3]]:
                raise ArithmeticError('Incomplete or unordered prefix')
            jump=0. if previous is None else float(max(abs(w/previous-1)))
            trigger=jump>run.policy['neighbour_relative_jump']
            if trigger and key not in audits:
                queries=[]
                for value,expected in zip([run.model.omega(run.policy['Lambda_lower']),*(.5*(w[:-1]+w[1:]))],(0,1,2)):
                    _,count,d=physics.count_at(run.model,case,float(value),run.cats(run.model,case),run.policy)
                    if count!=expected:raise ArithmeticError('Triggered count audit failed')
                    queries.append(dict(query=float(value),expected=expected,count=count,**d))
                audits[key]=dict(ID=label,beta_deg=b,queries=queries,status='PASS',resolution='count-one prefix and independent low/midgap counts; physical gates passed; retain actual steep growth')
                write(auditfile,audits)
            neighbours.append(dict(ID=label,beta_deg=b,relative_jump=jump,trigger=trigger,gap12_relative=float(w[1]/w[0]-1),status='PASS'))
            previous=w;row=dict(ID=label,beta_deg=b,kappa_joint=case.k_joint,kappa_end=case.k_end)
            for k,r in enumerate(rec['roots'],1):
                if not all(r['gates'].values()):raise ArithmeticError('Quality gate failed')
                row.update({f'{f}{k}':r[f] for f in ('omega','Omega','Lambda')})
                if not math.isclose(r['omega']*run.model.time,r['Lambda']**2,rel_tol=2e-14):raise ArithmeticError('Frequency normalization failed')
                d=r['diagnostics'];checks.append(dict(ID=label,beta_deg=b,position=k,origin=origin,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),
                    ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],status='PASS'))
            row.update(origin=origin,point_file=path.relative_to(ROOT).as_posix(),point_sha256=sha(path),shape_sha256=rec['shape_sha256']);rows.append(row)
        print(f'Thickness {label}: {len(beta)} prefixes PASS',flush=True)
        if label==pre[-1]:
            profiles,transitions,error=thickness_profiles(cfg,rows)
            write(run.folder/'preflight.json',dict(status='PASS',parameter_points=len(rows),transitions=transitions,q_identity_error=error))
            csvwrite(run.folder/'preflight_profiles.csv',profiles)
            print('Preflight: five complete angular curves and both tail criteria PASS',flush=True)
    if preflight_only:return None
    by={(r['ID'],r['beta_deg']):r for r in rows};order=[]
    for j in range(5):
        for e in range(5):
            a=f'J{j}E{e}'
            for other,axis in [(f'J{j+1}E{e}','joint')]*(j<4)+[(f'J{j}E{e+1}','end')]*(e<4):
                for b in beta:
                    for k in (1,2,3):
                        margin=by[other,b][f'omega{k}']/by[a,b][f'omega{k}']-1
                        if margin < -cfg['comparison_relative_tolerance']:raise ArithmeticError('Spectral stiffness monotonicity failed')
                        order.append(dict(axis=axis,lower=a,upper=other,beta_deg=b,position=k,relative_margin=margin,status='PASS'))
    if len(rows)!=2050 or sum(r['origin']!='new' for r in rows)!=cfg['expected_reused_points']:raise ArithmeticError('Grid accounting mismatch')
    for name,data in [('spectrum',rows),('reused_nodes',[r for r in rows if r['origin']!='new']),('quality',checks),('neighbour_audit',neighbours),('spectral_order',order)]:csvwrite(run.folder/(name+'.csv'),data)
    return rows,checks,neighbours,order


def thickness_analysis(folder,cfg,rows,old,oldtrans):
    profiles,transitions,error=thickness_profiles(cfg,rows);t={(r['ID'],r['reference_end_deg']):r for r in transitions}
    prior={r['ID']:float(r['beta_pl_deg']) for r in oldtrans};comparison=[]
    for label,pair in cfg['cases'].items():
        a=prior[label];b=t[label,5]['beta_pl_deg'];c=t[label,10]['beta_pl_deg']
        comparison.append(dict(ID=label,kappa_joint=pair[0],kappa_end=pair[1],beta_old_deg=a,beta_new5_deg=b,beta_scaled_deg=c,
            scaled_half_deg=None if c is None else c/2,D_deg=None if c is None else c/2-a,R_beta=None if c is None else c/(2*a),
            fixed_shift_deg=None if b is None else b-a,normalization_shift_deg=None if None in (b,c) else c-b,
            status5=t[label,5]['status'],status_scaled=t[label,10]['status']))
    effects=[];mixed=[];spans=[];frequency=[]
    for name,epsilon,angles in [('old',.005,prior),('new5',.01,{k:t[k,5]['beta_pl_deg'] for k in prior}),('scaled',.01,{k:t[k,10]['beta_pl_deg'] for k in prior})]:
        source=old if name=='old' else rows;by={(r['ID'],r['beta_deg']):r for r in source}
        freqs=[dict(ID=k,Lambda10=by[k,10]['Lambda1'],Lambda90=by[k,90]['Lambda1']) for k in prior]
        inc,mix=stiffness_effects([dict(ID=k,q_star=.15,beta_pl_deg=v) for k,v in angles.items()],freqs)
        effects.extend(dict(surface=name,epsilon=epsilon,**r,delta_beta_over_epsilon=None if r['delta_beta_deg'] is None else r['delta_beta_deg']/epsilon) for r in inc)
        mixed.extend(dict(surface=name,**r) for r in mix)
        for axis in ('joint','end'):
            for fixed in range(5):
                a,b=[f'J{i}E{fixed}' if axis=='joint' else f'J{fixed}E{i}' for i in (0,4)]
                v=None if None in (angles[a],angles[b]) else angles[b]-angles[a]
                spans.append(dict(surface=name,axis=axis,fixed_kappa=cfg['stiffness_values'][fixed],delta_beta_deg=v,delta_beta_over_epsilon=None if v is None else v/epsilon))
        for k in prior:
            b=angles[k];frequency.append(dict(surface=name,ID=k,beta_pl_deg=b,Lambda_pl=None if b is None else by[k,b]['Lambda1'],
                Lambda0=by[k,0]['Lambda1'],Lambda10=by[k,10]['Lambda1'],Lambda20=by[k,20]['Lambda1'],Lambda90=by[k,90]['Lambda1']))
    similarity=[];curves=[];oldby={(r['ID'],r['beta_deg']):r for r in old};newby={(r['ID'],r['beta_deg']):r for r in rows}
    for label in cfg['preflight_case_ids']:
        for source,eps in [(old,.005),(rows,.01)]:
            base=next(r['omega1'] for r in source if r['ID']==label and r['beta_deg']==0)
            for r in source:
                if r['ID']==label and r['beta_deg']<=12*eps/.005:
                    curves.append(dict(ID=label,epsilon=eps,beta_deg=r['beta_deg'],xi=math.radians(r['beta_deg'])/eps,F_omega=r['omega1']/base,interpolated=False))
        for label0,b in oldby:
            if label0==label and b<=12 and (label,2*b) in newby:
                f0=oldby[label,b]['omega1']/oldby[label,0]['omega1'];f1=newby[label,2*b]['omega1']/newby[label,0]['omega1']
                similarity.append(dict(ID=label,beta_old_deg=b,beta_new_deg=2*b,xi=math.radians(b)/.005,F_old=f0,F_new=f1,relative_difference=f1/f0-1,interpolated=False))
    for name,data in [('profiles',profiles),('transitions',transitions),('comparison',comparison),('increments',effects),('mixed_differences',mixed),('full_range_effects',spans),('frequency_levels',frequency),('similarity',similarity),('similarity_curves',curves)]:csvwrite(folder/(name+'.csv'),data)
    return dict(q_identity_max_error=error,accepted_transitions=sum(r['status']=='PASS' for r in transitions),missing_transitions=sum(r['status']!='PASS' for r in transitions),
        late24_q_max=max(r['late24_q_max'] for r in transitions),wide_scheme_max_angle_shift=max(abs(r['wide_minus_primary_deg'] or 0) for r in transitions),
        negative_angle_increments=sum(r['delta_beta_deg'] is not None and r['delta_beta_deg']<0 for r in effects),
        similarity_max_relative_error=max(abs(r['relative_difference']) for r in similarity))


def plot_thickness(folder):
    """Saved scalar tables only: no physical assembly, recovery, or root calls."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder);names=['comparison.csv','similarity_curves.csv','contract.json']
    if (folder/'manifest.json').exists():
        m=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        if any(sha(folder/n)!=m['artifacts'][n] for n in names):raise ValueError('Thickness plot inputs changed')
    rows=list(csv.DictReader((folder/'comparison.csv').open(encoding='utf8')));by={r['ID']:r for r in rows}
    matrix=lambda field:np.array([[float(by[f'J{j}E{e}'][field]) if by[f'J{j}E{e}'][field] else np.nan for j in range(5)] for e in range(5)])
    cfg=json.loads((folder/'contract.json').read_text())['config'];v=cfg['stiffness_values'];edges=10.**np.arange(-1.25,1.51,.5);labels=['0.1','0.316','1','3.162','10']
    figures=folder/'figures';figures.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    def save(fig,name):
        fig.savefig(figures/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(figures/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    def maps(fields,titles,name,difference=False):
        matrices=[matrix(f) for f in fields];lo=min(np.nanmin(a) for a in matrices);hi=max(np.nanmax(a) for a in matrices)
        if difference:lo,hi=-max(abs(lo),abs(hi),.25),max(abs(lo),abs(hi),.25)
        fig,axes=plt.subplots(1,len(fields),figsize=(5.7*len(fields),5),layout='constrained',squeeze=False)
        cmap=plt.colormaps['RdBu_r' if difference else 'viridis'].copy();cmap.set_bad('.85')
        for ax,a,title in zip(axes[0],matrices,titles):
            mesh=ax.pcolormesh(edges,edges,np.ma.masked_invalid(a),vmin=lo,vmax=hi,cmap=cmap,rasterized=False)
            for e in range(5):
                for j in range(5):
                    x=a[e,j];ax.text(v[j],v[e],'—' if np.isnan(x) else f'{x:.2f}',ha='center',va='center',color='black' if difference or x>(lo+hi)/2 else 'white')
            ax.set(xscale='log',yscale='log',xticks=v,yticks=v,xticklabels=labels,yticklabels=labels,xlabel='κ_joint',ylabel='κ_end',title=title);ax.minorticks_off()
        cb=fig.colorbar(mesh,ax=list(axes[0]),label='D, °' if difference else 'Характерный угол, °');cb.solids.set_rasterized(False);save(fig,name)
    maps(['beta_old_deg','beta_new5_deg'],['ε = 0.005; опора 0–5°','ε = 0.01; опора 0–5°'],'01_fixed_reference_maps')
    maps(['beta_old_deg','scaled_half_deg'],['ε = 0.005; β_pl','ε = 0.01; β_pl,scaled / 2'],'02_scaled_maps')
    maps(['D_deg'],['β_pl,scaled(0.01) / 2 − β_pl(0.005)'],'03_scaling_difference',True)
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),layout='constrained');colors=plt.colormaps['viridis'](np.linspace(.05,.85,5))
    old=matrix('beta_old_deg');new=matrix('scaled_half_deg')
    for ax,axis in zip(axes,('joint','end')):
        for i in range(5):
            a=old[i,:] if axis=='joint' else old[:,i];b=new[i,:] if axis=='joint' else new[:,i]
            ax.plot(v,a,'-o',color=colors[i],ms=3,label=f'κ_{"end" if axis=="joint" else "joint"}={labels[i]}')
            ax.plot(v,b,'--x',color=colors[i],ms=5)
        ax.set(xscale='log',xticks=v,xticklabels=labels,xlabel=f'κ_{axis}',ylabel='Приведённый угол, °');ax.minorticks_off();ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Сплошная: ε=0.005; штриховая: ε=0.01, масштабированная опора, угол / 2');save(fig,'04_spring_sections')
    curves=list(csv.DictReader((folder/'similarity_curves.csv').open(encoding='utf8')))
    fig,axes=plt.subplots(2,3,figsize=(11,6),layout='constrained')
    for ax,label in zip(axes.flat,cfg['preflight_case_ids']):
        for eps,style in [(.005,'-o'),(.01,'--x')]:
            s=sorted([r for r in curves if r['ID']==label and float(r['epsilon'])==eps],key=lambda r:float(r['xi']))
            ax.plot([float(r['xi']) for r in s],[float(r['F_omega']) for r in s],style,ms=2,lw=1,label=f'ε={eps}')
        ax.set(title=str(tuple(cfg['cases'][label])),xlabel='ξ = β_rad / ε',ylabel='F_ω');ax.grid(alpha=.2);ax.legend(fontsize=8)
    axes.flat[-1].axis('off');axes.flat[-1].text(.03,.65,'Только рассчитанные узлы.\nЛинии соединяют соседние точки.\nИнтерполяция для сравнения не нужна.',fontsize=10);save(fig,'05_frequency_similarity')


def thickness_main(args):
    cfg=json.loads(THICKNESS_CONFIG.read_text());key,ident,archives,old,oldtrans=thickness_sources(cfg)
    folder=OUTPUT/'growth_thickness_comparison'/key
    if args.verify_cache or (folder/'manifest.json').exists():
        p=(args.verify_cache or folder).resolve();m=json.loads((p/'manifest.json').read_text(encoding='utf8'))
        if m['identity']!=ident or any(sha(p/f)!=h for f,h in m['artifacts'].items()):raise ValueError('Thickness cache mismatch')
        print(f'REUSED {p}: 0 eigenvalue solves');return
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',ident)
    run=Run(cfg,folder);m=run.model;s=m.section;oldm=physics.Model(**{**cfg['material'],'epsilon':.005})
    h=math.sqrt(12)*m.l*m.epsilon;b=4*h
    if not (math.isclose(s.area,b*h) and math.isclose(s.second_moment,b*h**3/12) and math.isclose(m.J,m.rho*m.nu**2*17*s.second_moment)):
        raise ArithmeticError('Section geometry mismatch')
    write(folder/'physical_contract.json',dict(model='RL_FULL + Timoshenko; two exact translational constraints',material=cfg['material'],
        h=h,b=b,A=s.area,I=s.second_moment,Ip=17*s.second_moment,J_RL=m.J,time=m.time,area_ratio=s.area/oldm.section.area,
        I_ratio=s.second_moment/oldm.section.second_moment,dimensional_springs={k:[x*s.EI/m.l for x in pair] for k,pair in cfg['cases'].items()}))
    ledgerpath=folder/'run_ledger.json';ledger=json.loads(ledgerpath.read_text()) if ledgerpath.exists() else []
    start=time.perf_counter();spectral=postprocess=render=0.;result=None
    try:
        result=thickness_spectra(run,archives,args.stage=='preflight');spectral=time.perf_counter()-start
        if result is not None:
            rows,checks,neighbours,order=result;now=time.perf_counter();diag=thickness_analysis(folder,cfg,rows,old,oldtrans);postprocess=time.perf_counter()-now
            diag.update(status='PASS',maxima={k:max(r[k] for r in checks) for k in ('physical','ode','energy','singular','mass_gram','condition')},
                neighbour_triggers=sum(r['trigger'] for r in neighbours),minimum_gap12_relative=min(r['gap12_relative'] for r in neighbours),
                spectral_monotonicity_comparisons=len(order),minimum_spectral_margin=min(r['relative_margin'] for r in order))
            write(folder/'diagnostics.json',diag);now=time.perf_counter();plot_thickness(folder);render=time.perf_counter()-now
    except Exception as exc:
        write(folder/'failure.json',dict(error=repr(exc),new_points=run.new_points,new_roots=run.new_roots));raise
    finally:
        ledger.append(dict(command=sys.argv,new_points=run.new_points,new_roots=run.new_roots,checkpoint_reuse_calls=run.reused_points,
            pole_roots=run.pole_roots,ritz_solves=run.ritz_solves,spectral_seconds=spectral,postprocessing_seconds=postprocess,render_seconds=render,wall_seconds=time.perf_counter()-start));write(ledgerpath,ledger)
    if result is not None:
        if thickness_sources(cfg)[1]!=ident:raise ValueError('Sources changed during computation')
        totals={k:sum(r[k] for r in ledger) for k in ('new_points','new_roots','pole_roots','ritz_solves')}
        if totals!=dict(new_points=1996,new_roots=5988,pole_roots=0,ritz_solves=0):raise ArithmeticError('Unexpected solve accounting')
        write(folder/'manifest.json',dict(status='PASS',identity=ident,git_head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            **totals,total_parameter_points=2050,reused_points=54,archived_surface_points=25,archived_spectrum_points=2075,
            postprocessing=dict(profiles=4100,transitions=50,comparisons=25),
            artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    print(folder,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['all', 'gates', 'reference', 'map', 'regression', 'preflight'], default='all')
    parser.add_argument('--preset', choices=['verification', 'limit-angles', 'finite-angle-maps', 'saturation-analysis', 'plateau-criterion', 'plateau-grid', 'early-saturation', 'growth-transition', 'growth-refinement', 'growth-stiffness-grid', 'mechanical-interpretation', 'angle-variation', 'epsilon-scaling', 'longitudinal-slip', 'thickness-grid'], default='verification')
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
        elif args.preset == 'growth-stiffness-grid':
            plot_stiffness(args.plot_only)
        elif args.preset == 'mechanical-interpretation':
            plot_mechanics(args.plot_only)
        elif args.preset == 'angle-variation':
            plot_variation(args.plot_only)
        elif args.preset == 'epsilon-scaling':
            plot_epsilon(args.plot_only)
        elif args.preset == 'longitudinal-slip':
            plot_slip(args.plot_only)
        elif args.preset == 'thickness-grid':
            plot_thickness(args.plot_only)
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
    if args.preset == 'growth-stiffness-grid':
        stiffness_main(args); return
    if args.preset == 'mechanical-interpretation':
        mechanics_main(args); return
    if args.preset == 'angle-variation':
        variation_main(args); return
    if args.preset == 'epsilon-scaling':
        epsilon_main(args); return
    if args.preset == 'longitudinal-slip':
        slip_main(args); return
    if args.preset == 'thickness-grid':
        thickness_main(args); return
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
