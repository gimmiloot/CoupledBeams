"""Bounded seven-field action audit; no nonlinear evolution or parameter map.

This CLI has a distinct symbolic/action/jet contract and therefore cannot be
a preset of a linear spectrum CLI. Old arm solvers remain reference inputs.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import importlib.metadata
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from scripts.lib import weakly_nonlinear_spatial_rod as rod

VERSION = "nlsp-action-audit-v1"
CONFIG = ROOT / "data/input/weakly_nonlinear_spatial_rod.json"
OUTPUT = ROOT / "results/weakly_nonlinear_spatial_rod"
APPENDIX = ROOT / "docs/theory/weakly_nonlinear_spatial_rod_expansion_generated.md"
PDF = ROOT / "docs/literature/pdf/0020-76832990087-x.pdf"
P = np.diag([1., -1., 1.])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf8")
    temp.replace(path)


def identity(config_path=CONFIG):
    config_path = Path(config_path)
    config = json.loads(config_path.read_text(encoding="utf8"))
    if config["schema"] != "nlsp-cubic-audit-v1" or config["model"] != "seven_field_reference_section_V0":
        raise ValueError("Unsupported model/config contract")
    if config["field_order"] != list(rod.FIELD_ORDER):
        raise ValueError("Configured field ordering differs from the fixed model")
    expected = {"E":1.,"rho":1.,"nu":.3,"b":.20,"h":.05,"L":1.,"kappa":"5/6"}
    if config["material_geometry"] != expected:
        raise ValueError("This bounded reference audit supports only the declared G20 geometry and kappa=5/6")
    if config["length_scale"] <= 0 or config["time_scale"] != 1:
        raise ValueError("Positive reference length and the declared unit time scale are required")
    source = ROOT / config["supplied_expansion"]["path"]
    if sha(source) != config["supplied_expansion"]["sha256"]:
        raise ValueError("Supplied expansion SHA256 mismatch")
    if sha(ROOT/config["primary_source"]["path"]) != config["primary_source"]["sha256"]:
        raise ValueError("Primary source PDF SHA256 mismatch")
    paths = [config_path, source, Path(__file__), Path(rod.__file__), PDF,
             ROOT / "scripts/lib/mindlin_herrmann_longitudinal.py",
             ROOT / "scripts/lib/mindlin_herrmann_timoshenko_joint.py",
             ROOT / "scripts/lib/yartsev_ch2_monoclinic_rod.py",
             ROOT / "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"]
    # Validated saved references are scientific inputs, not just provenance.
    reference_manifests = []
    for frozen in (ROOT/"results/mindlin_herrmann_timoshenko_single_rod/342ce44bff81c36f",
                   ROOT/"results/mindlin_herrmann_timoshenko_general_beta_joint/a36e72f715cb526b"):
        manifest = frozen/"manifest.json"
        if manifest.exists():
            reference_manifests.append(manifest)
            reference = json.loads(manifest.read_text(encoding="utf8"))
            paths += [manifest] + [frozen/name for name in reference["artifact_hashes"]]
    hashes = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p) for p in paths}
    item = {"schema": config["schema"], "model_version": rod.MODEL_VERSION,
            "audit_version": VERSION, "hashes": hashes, "python": sys.version,
            "reference_manifests": [str(p.relative_to(ROOT)) for p in reference_manifests],
            "dependencies": {n: importlib.metadata.version(n) for n in ("numpy", "scipy")}}
    fingerprint = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()[:16]
    return config, item, fingerprint


class SuppliedPolynomialParser:
    """Strict adapter for the supplied 21 expanded sums, not a LaTeX parser.

    Recognizes only named model symbols, rational fractions, small integer
    powers, sums and products. Unknown syntax fails rather than being guessed.
    No eval, user instructions or general document commands are executed.
    """
    def __init__(self, expression):
        for old, new in ((r"\mathcal C_T", "CT"), (r"j_{\parallel}", "jp"),
                         (r"j_{\perp}", "jb"), (r"B_{\parallel}", "Bp"),
                         (r"B_{\perp}", "Bb"), (r"\Phi", "Phi"),
                         (r"\psi", "psi"), (r"\theta", "theta"), (r"\nu", "nu")):
            expression = expression.replace(old, new)
        expression = re.sub(r"(u|w|v|Phi|psi|theta|c)_\{(ss|st|tt|s|t)\}", r"\1_\2", expression)
        expression = expression.replace(r"\quad", " ").replace(r"\\", " ").replace("&", " ")
        self.tokens = re.findall(r"\\frac|[A-Za-z][A-Za-z_0-9]*|\d+|[{}()+\-*/^]", expression)
        compact = re.sub(r"\s+", "", expression)
        if "".join(self.tokens) != compact:
            raise ValueError("Unsupported supplied-expression syntax: " + expression)
        self.index = 0

    def peek(self):
        return self.tokens[self.index] if self.index < len(self.tokens) else None

    def take(self):
        token = self.peek()
        if token is None:
            raise ValueError("Unexpected end of supplied polynomial")
        self.index += 1
        return token

    def atom(self):
        token = self.take()
        if token == "+":
            return self.atom()
        if token == "-":
            return -self.atom()
        if token in ("{", "("):
            value = self.expression()
            if self.take() != ("}" if token == "{" else ")"):
                raise ValueError("Unbalanced polynomial group")
        elif token == r"\frac":
            value, denominator = self.atom(), self.atom()
            if any(key for key in denominator.terms) or not denominator:
                raise ValueError("Only rational constant denominators are allowed")
            value = value / denominator.terms[()]
        elif token.isdigit():
            value = rod.Polynomial(int(token))
        elif token in rod.SYMBOL_ORDER:
            value = rod.Polynomial.symbol(token)
        else:
            raise ValueError("Unknown supplied symbol: " + token)
        if self.peek() == "^":
            self.take()
            power = self.atom()
            if set(power.terms) != {()} or power.terms[()].denominator != 1:
                raise ValueError("Noninteger supplied power")
            if not 0 <= int(power.terms[()]) <= 3:
                raise ValueError("Unsupported supplied power")
            value = value ** int(power.terms[()])
        return value

    def product(self):
        value = self.atom()
        while self.peek() is not None and self.peek() not in ("+", "-", "}", ")"):
            if self.peek() == "*":
                self.take()
            other = self.atom()
            left_degree = max((sum(i < 42 for i in key) for key in value.terms), default=0)
            right_degree = max((sum(i < 42 for i in key) for key in other.terms), default=0)
            if left_degree + right_degree > 4:
                raise ValueError("Supplied expression exceeds the strict adapter degree")
            value = value * other
        return value

    def expression(self):
        value = self.product()
        while self.peek() in ("+", "-"):
            operation = self.take()
            term = self.product()
            value = value + term if operation == "+" else value - term
        return value

    def parse(self):
        value = self.expression()
        if self.peek() is not None:
            raise ValueError("Unconsumed supplied polynomial tokens")
        return value


def compare_supplied(model, path):
    text = Path(path).read_text(encoding="utf8")
    pattern = r"\\mathcal E_\{([^}]+)\}\^\{\[(\d)\]\}\s*&=(.*?)\\end\{aligned\}"
    blocks = re.findall(pattern, text, re.S)
    if len(blocks) != 21:
        raise ValueError(f"Expected exactly 21 supplied equations, found {len(blocks)}")
    seen, comparisons = set(), []
    for name, degree, body in blocks:
        name, degree = name.lstrip("\\"), int(degree)
        key = (name, degree)
        if key in seen or name not in rod.FIELD_ORDER or degree not in (1, 2, 3):
            raise ValueError("Duplicate or unsupported supplied equation")
        seen.add(key)
        supplied = SuppliedPolynomialParser(body).parse()
        if supplied != supplied.homogeneous(degree):
            raise ValueError("Supplied equation is not homogeneous at its claimed degree")
        difference = model.residual_a[rod.FIELD_ORDER.index(name)].homogeneous(degree) - supplied
        comparisons.append({"field": name, "degree": degree,
                            "status": "MATCH" if not difference else "DISCREPANCY",
                            "difference": difference.serialize(), "difference_text": str(difference)})
    return {"status": "MATCH" if all(c["status"] == "MATCH" for c in comparisons) else "DISCREPANCY",
            "sha256": sha(path), "comparisons": comparisons}


def restriction(fields):
    return {name + suffix: 0 for name in fields for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}


def exact_checks(model):
    p = model.symbols
    linear = (p["m"]*p["u_tt"]-p["C"]*p["u_ss"]-p["nu"]*p["C"]*p["c_s"],
              p["m"]*p["w_tt"]-p["S"]*(p["w_ss"]-p["theta_s"]),
              p["m"]*p["v_tt"]-p["S"]*(p["v_ss"]-p["psi_s"]),
              (p["jp"]+p["jb"])*p["Phi_tt"]-p["CT"]*p["Phi_ss"],
              p["jb"]*p["psi_tt"]-p["Bb"]*p["psi_ss"]-p["S"]*(p["v_s"]-p["psi"]),
              p["jp"]*p["theta_tt"]-p["Bp"]*p["theta_ss"]-p["S"]*(p["w_s"]-p["theta"]),
              p["jp"]*p["c_tt"]-p["H"]*p["c_ss"]+p["C"]*(p["c"]+p["nu"]*p["u_s"]))
    linear_differences = [a.homogeneous(1)-b for a, b in zip(model.residual_a, linear)]
    axial = restriction(("w", "v", "Phi", "psi", "theta"))
    axial_differences = [a.substitute(axial)-b.substitute(axial) for a, b in zip(model.residual_a, linear)]
    # Independent planar trigonometric balances, not restriction of A/B.
    theta, c = p["theta"], p["c"]
    cosine, sine = 1-theta**2/2+theta**4/24, theta-theta**3/6
    g1 = ((1+p["u_s"])*cosine+p["w_s"]*sine-1).truncate(3)
    g2 = (-(1+p["u_s"])*sine+p["w_s"]*cosine).truncate(3)
    N, Q = p["C"]*(g1+p["nu"]*c), p["S"]*g2
    planar = (p["m"]*p["u_tt"]-(N*cosine-Q*sine).truncate(3).total_derivative("s"),
              p["m"]*p["w_tt"]-(N*sine+Q*cosine).truncate(3).total_derivative("s"),
              rod.Polynomial(), rod.Polynomial(), rod.Polynomial(),
              (p["jp"]*(1+c)**2*p["theta_t"]).total_derivative("t")-p["Bp"]*p["theta_ss"]-((1+g1)*Q-g2*N).truncate(3),
              p["jp"]*p["c_tt"]-p["H"]*p["c_ss"]+p["C"]*(c+p["nu"]*g1)-p["jp"]*(1+c)*p["theta_t"]**2)
    planar_differences = [a.substitute(restriction(("v", "Phi", "psi")))-b.truncate(3) for a, b in zip(model.residual_a, planar)]
    signs = (1, 1, -1, -1, -1, 1, 1)
    reflection = {name+suffix: sign*p[name+suffix] for name, sign in zip(rod.FIELD_ORDER, signs)
                  for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}
    symmetry = [model.T4.substitute(reflection)-model.T4, model.V4.substitute(reflection)-model.V4]
    symmetry += [a.substitute(reflection)-sign*a for sign, a in zip(signs, model.residual_a)]
    symmetry += [a.substitute(reflection)-sign*a for sign, a in zip(signs, model.flux_a)]
    power = sum((p[name+"_t"]*f for name, f in zip(rod.FIELD_ORDER, model.flux_a)), rod.Polynomial())
    energy_identity = (model.T4+model.V4).total_derivative("t")-power.total_derivative("s")-sum(
        (p[name+"_t"]*e for name, e in zip(rod.FIELD_ORDER, model.residual_a)), rod.Polynomial())
    rigid = {"c": 0, "c_s": 0, **{name+"_s": 0 for name in ("Phi", "psi", "theta")}}
    # R e1-e1 directly from independently generated exponential series.
    rigid.update({name+"_s": model.rotation_b[i][0]-int(i == 0) for i, name in enumerate(("u", "w", "v"))})
    rigid_flux = [f.substitute(rigid).truncate(3) for f in model.flux_a]
    mass_symmetry = [model.T4.derivative(a+"_t").derivative(b+"_t")-
                     model.T4.derivative(b+"_t").derivative(a+"_t") for a in rod.FIELD_ORDER for b in rod.FIELD_ORDER]
    # Deliberately omit Jr transpose in a test copy, then omit c inertia.
    wrong_z = (model.body_b[0], -model.body_b[1], model.body_b[2])
    wrong_z_differences = [(a-b).truncate(3) for a, b in zip(model.residual_a[3:6], wrong_z)]
    omitted_c_inertia = p["jp"]*(1+c)*(model.omega_b[0]**2+model.omega_b[2]**2)
    # beta0 right arm to a single global coordinate: ds reverses orientation.
    reversal_signs = (-1, -1, 1, -1, -1, 1, 1)
    reversal = {name+suffix: sign*(-1 if suffix in ("_s", "_st") else 1)*p[name+suffix]
                for name, sign in zip(rod.FIELD_ORDER, reversal_signs)
                for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}
    reversal_energy = [model.T4.substitute(reversal)-model.T4, model.V4.substitute(reversal)-model.V4]
    reversal_flux = [f.substitute(reversal)+sign*f for sign, f in zip(reversal_signs, model.flux_a)]
    groups = {"linear": linear_differences, "pure_axial": axial_differences, "planar": planar_differences,
              "reflection": symmetry, "energy": [energy_identity], "rigid_flux": rigid_flux,
              "mass_symmetry": mass_symmetry, "straight_reversal_energy": reversal_energy,
              "straight_reversal_flux": reversal_flux,
              "boundary": [a-b for a, b in zip(model.flux_a, model.flux_b)]}
    result = {name: {"status": "PASS" if all(not p for p in polys) else "FAIL",
                     "coefficient_differences": [p.serialize() for p in polys]} for name, polys in groups.items()}
    result["negative_controls"] = {"status": "PASS" if all(wrong_z_differences) and omitted_c_inertia.truncate(3) else "FAIL",
                                   "omit_Jr_transpose": [p.serialize() for p in wrong_z_differences],
                                   "omit_c_inertia": omitted_c_inertia.truncate(3).serialize()}
    dimensions = {}
    for name in rod.SYMBOL_ORDER:
        if name in rod.COEFFICIENT_UNITS:
            dimensions[name] = rod.COEFFICIENT_UNITS[name]
        else:
            field, _, suffix = name.partition("_")
            dimensions[name] = (0, int(field in ("u","w","v"))-suffix.count("s"), -suffix.count("t"))
    def dimensional_differences(poly, expected):
        return [list(key) for key in poly.terms if tuple(sum(dimensions[rod.SYMBOL_ORDER[i]][axis] for i in key) for axis in range(3)) != expected]
    bad_dimensions = dimensional_differences(model.T4,(1,1,-2))+dimensional_differences(model.V4,(1,1,-2))
    for index, (residual,flux) in enumerate(zip(model.residual_a,model.flux_a)):
        bad_dimensions += dimensional_differences(residual,(1,0 if index<3 else 1,-2))
        bad_dimensions += dimensional_differences(flux,(1,1 if index<3 else 2,-2))
    result["dimensions"] = {"status":"PASS" if not bad_dimensions else "FAIL", "bad_monomials":bad_dimensions,
                            "units":"kg,m,s exponents; densities per original length"}
    return result


def manufactured_jet(profile, s, t):
    """Three predeclared smooth all-field test functions and analytic jets."""
    j = np.arange(1., 8.)
    Aj = (0.12+0.01*j)*((-1.)**(j+profile))
    Bj = 0.04+0.004*j
    ks, wt, ps, vt = 0.4*j*(profile+1), 0.3*(j+1), 0.17*(j+2), 0.23*(j+profile+1)
    phase = 0.1*j+0.2*profile
    a, b = ks*s+wt*t+phase, ps*s+vt*t-phase
    q = Aj*np.sin(a)+Bj*np.cos(b)
    qs = Aj*ks*np.cos(a)-Bj*ps*np.sin(b)
    qt = Aj*wt*np.cos(a)-Bj*vt*np.sin(b)
    qss = -Aj*ks**2*np.sin(a)-Bj*ps**2*np.cos(b)
    qst = -Aj*ks*wt*np.sin(a)-Bj*ps*vt*np.cos(b)
    qtt = -Aj*wt**2*np.sin(a)-Bj*vt**2*np.cos(b)
    return rod.FieldJet(q, qs, qt, qss, qst, qtt)


def amplitude_checks(config, coefficients, model):
    g = config["material_geometry"]
    EA, ell = g["E"]*g["b"]*g["h"], config["length_scale"]
    scales = np.array([EA/ell]*3+[EA]*4)
    flux_scales = np.array([EA]*3+[EA*ell]*4)
    rows, previous = [], None
    for epsilon in config["amplitude_epsilons"]:
        components, boundaries, per_profile = [], [], []
        for profile in range(config["manufactured_profiles"]):
            profile_errors = []
            for s, t in config["manufactured_samples"]:
                jet = manufactured_jet(profile, s, t).scaled(epsilon)
                full = rod.full_evaluate(jet, coefficients)
                cubic = rod.polynomial_evaluate(jet, coefficients, model)
                err = (full["residual"]-cubic["residual"])/scales
                components.append(err)
                boundaries.append((full["flux"]-cubic["flux"])/flux_scales)
                profile_errors.append(err)
            per_profile.append(float(np.linalg.norm(profile_errors)))
        norms = np.linalg.norm(np.array(components), axis=0)
        norm = float(np.linalg.norm(components))
        row = {"epsilon_a": epsilon, "component_errors": norms.tolist(), "aggregate_error": norm,
               "component_above_reporting_floor": (norms > config["policy"]["numerical_floor_scaled"]).tolist(),
               "profile_errors": per_profile, "boundary_error": float(np.linalg.norm(boundaries)),
               "ratios": None if previous is None else (np.divide(previous["components"], norms, out=np.zeros(7), where=norms != 0)).tolist(),
               "aggregate_ratio": None if previous is None else previous["norm"]/norm,
               "aggregate_order": None if previous is None else math.log2(previous["norm"]/norm),
               "component_orders": None if previous is None else [math.log2(x/y) if x>0 and y>0 else None for x,y in zip(previous["components"], norms)]}
        rows.append(row)
        previous = {"norm": norm, "components": norms}
    lower, upper = config["policy"]["amplitude_aggregate_final_order"]
    above_floor = rows[-1]["aggregate_error"] > config["policy"]["numerical_floor_scaled"]
    return {"status": "PASS" if above_floor and lower <= rows[-1]["aggregate_order"] <= upper else "PARTIAL",
            "residual_scales": scales.tolist(), "flux_scales": flux_scales.tolist(),
            "samples": config["manufactured_samples"], "rows": rows,
            "numerical_floor_scaled": config["policy"]["numerical_floor_scaled"],
            "floor_policy": "conservative predeclared reporting floor; flagged component slopes are retained but not used as acceptance evidence",
            "final_above_floor": above_floor, "interpretation": "manufactured jets, not PDE solutions or trajectory/frequency error"}


def frame_transforms(beta):
    from scripts.lib.mindlin_herrmann_timoshenko_joint import frames
    result = []
    for frame in frames(beta):
        B = np.array([[frame.t[0], frame.n[0], 0.], [frame.t[1], frame.n[1], 0.], [0., 0., -1.]])
        transform = np.zeros((7, 7))
        transform[:3, :3], transform[3:6, 3:6], transform[6, 6] = B, B@P, 1.
        result.append(transform)
    return result


def joint_operator(beta):
    T1, T2 = frame_transforms(beta)
    zero = np.zeros((7, 7))
    return np.block([[T1, zero, -T2, zero], [zero, T1, zero, T2]])


def joint_checks(config, coefficients, model):
    rng = np.random.default_rng(713)
    geometry, maxwork = [], 0.
    for beta in (0., 45., 90.):
        transforms = frame_transforms(beta)
        J = joint_operator(beta)
        for T in transforms:
            for _ in range(16):
                q, f, variation = rng.normal(size=(3, 7))
                error = abs(f@variation-(T@f)@(T@variation))
                maxwork = max(maxwork, error/max(1., np.linalg.norm(f)*np.linalg.norm(variation)))
        geometry.append({"beta_deg": beta, "transforms": [T.tolist() for T in transforms],
                         "rank": int(np.linalg.matrix_rank(J)), "JJt_error": float(np.max(abs(J@J.T-2*np.eye(14))))})
    signs = np.array([-1., -1., 1., -1., -1., 1., 1.])
    def reverse(jet):
        return rod.FieldJet(*(signs*getattr(jet, name)*(-1 if name in ("qs", "qst") else 1) for name in rod.JET_ORDER))
    original = manufactured_jet(1, .35, .17).scaled(.02)
    reversed_jet = reverse(original)
    f1 = rod.polynomial_evaluate(original, coefficients, model)["flux"]
    f2 = rod.polynomial_evaluate(reversed_jet, coefficients, model)["flux"]
    T1, T2 = frame_transforms(0)
    balance = T1@f1+T2@f2
    common_velocity = T1@original.qt
    power = float(common_velocity@balance)
    def integral(start, stop, right=False):
        nodes, weights = np.polynomial.legendre.leggauss(48)
        value = 0.
        for x, weight in zip(start+(nodes+1)*(stop-start)/2, weights):
            jet = manufactured_jet(1, x, .17).scaled(.02)
            if right:
                jet = reverse(jet)
            values = rod.polynomial_evaluate(jet, coefficients, model)
            value += weight*(values["T"]-values["V"])*(stop-start)/2
        return value
    whole = integral(0., 1.)
    splits = [{"split": split, "action_additivity_error": abs(whole-integral(0,split)-integral(split,1,True))} for split in (.5,.35)]
    maximum = max(maxwork, np.linalg.norm(balance), abs(power), *(r["action_additivity_error"] for r in splits))
    return {"status": "PASS" if maximum < config["policy"]["boundary_work_absolute"] and all(r["rank"]==14 for r in geometry) else "FAIL",
            "geometry": geometry, "max_duality_error": maxwork, "nonlinear_interface_balance": balance.tolist(),
            "nonlinear_interface_power": power, "action_splits": splits,
            "fixed_clamp_power": 0., "qualification": "c continuity is reduced common-DOF closure, not direct finite-joint 3D elasticity"}


def linear_state_from_action(coefficients, model):
    values = coefficients.values() | {name: 0. for name in rod.SYMBOL_ORDER if name not in rod.COEFFICIENT_ORDER}
    def hessian(poly, left, right):
        return np.array([[poly.derivative(a).derivative(b).evaluate(values) for b in right] for a in left])
    q, qs, qt = rod.FIELD_ORDER, [n+"_s" for n in rod.FIELD_ORDER], [n+"_t" for n in rod.FIELD_ORDER]
    M, K, A, D = hessian(model.T4,qt,qt), hessian(model.V4,qs,qs), hessian(model.V4,qs,q), hessian(model.V4,q,q)
    def state(omega):
        invK = np.linalg.inv(K)
        return np.block([[-invK@A, invK], [D-omega**2*M-A.T@invK@A, A.T@invK]])
    return state, {"M": M.tolist(), "K": K.tolist(), "A": A.tolist(), "D": D.tolist()}


def latex(poly):
    labels = {"Phi": r"\Phi", "psi": r"\psi", "theta": r"\theta", "nu": r"\nu",
              "jp": r"j_{\parallel}", "jb": r"j_{\perp}", "Bp": r"B_{\parallel}", "Bb": r"B_{\perp}", "CT": r"C_T"}
    terms = []
    for monomial, value in sorted(poly.terms.items()):
        factors = []
        for index in sorted(set(monomial)):
            name = rod.SYMBOL_ORDER[index]
            if "_" in name:
                field, suffix = name.split("_")
                label = labels.get(field,field)+"_{"+suffix+"}"
            else:
                label = labels.get(name,name)
            count = monomial.count(index)
            factors.append(label + ("^{"+str(count)+"}" if count>1 else ""))
        coefficient = abs(value)
        numerator = " ".join(factors)
        if coefficient.numerator != 1 or not numerator:
            numerator = str(coefficient.numerator)+(" "+numerator if numerator else "")
        if coefficient.denominator != 1:
            numerator = r"\frac{"+numerator+"}{"+str(coefficient.denominator)+"}"
        terms.append(("-" if value<0 else "+",numerator))
    if not terms:
        return "0"
    return "\n".join((sign if i or sign=="-" else "")+" "+term for i,(sign,term) in enumerate(terms))


def appendix_text(model, fingerprint):
    lines = ["# Семиполевая модель: generated quartic action / cubic expansion", "",
             "Generated by `verify_weakly_nonlinear_spatial_rod.py --compute`; do not hand-edit.",
             f"Model `{rod.MODEL_VERSION}`; fingerprint `{fingerprint}`.",
             "Exact rational coefficients; q=(u,w,v,Phi,psi,theta,c); all jet variables have amplitude degree one.",
             "Supplied analytical draft is stored separately and is not overwritten.", ""]
    for label, poly in (("T",model.T4),("V",model.V4)):
        for degree in (2,3,4):
            lines += [f"## {label}, degree {degree}", "", "\\[", "\\begin{aligned}",
                      latex(poly.homogeneous(degree)).replace("\n",r" \\"+"\n"), "\\end{aligned}", "\\]", ""]
    for name, expression in zip(rod.FIELD_ORDER, model.residual_a):
        lines += [f"## Residual {name}", ""]
        for degree in (1,2,3):
            lines += [f"### Degree {degree}", "", "\\[", "\\begin{aligned}",
                      latex(expression.homogeneous(degree)).replace("\n",r" \\"+"\n"), "\\end{aligned}", "\\]", ""]
    for name, expression in zip(rod.FIELD_ORDER, model.flux_a):
        lines += [f"## Boundary covector for {name}", "", "\\[", "\\begin{aligned}",
                  latex(expression).replace("\n",r" \\"+"\n"), "\\end{aligned}", "\\]", ""]
    return "\n".join(lines)


def run_audit(config, fingerprint):
    started = time.perf_counter()
    model = rod.derive_polynomials()
    comparisons = model.comparisons()
    if any(c["status"] != "PASS" for c in comparisons):
        raise RuntimeError("Path A/B promotion gate failed")
    # Reference module provides section reduction and independent old controls.
    from scripts.lib import yartsev_ch2_monoclinic_rod as yartsev
    g = config["material_geometry"]
    shear = g["E"]/(2*(1+g["nu"]))
    material = yartsev.BookMaterial(E1_real=g["E"],E2_real=g["E"],G12_real=shear,
                                    G13_real=shear,G23_real=shear,nu12=g["nu"],rho=g["rho"],
                                    eta1=0.,eta2=0.,eta12=0.,eta13=0.,eta23=0.)
    # This call is intentionally the existing rectangular generalized torsion.
    geometry = yartsev.Geometry(a=g["b"],b=g["h"],length=g["L"],shear_factor=5/6)
    reduction = yartsev.make_rod_point(0.,geometry=geometry,material=material,material_mode="elastic")
    CT = float(reduction.torsion.C_T.real)
    coefficients = rod.RodCoefficients.rectangular(g["E"],g["rho"],g["nu"],g["b"],g["h"],CT)
    checks = exact_checks(model)
    supplied = compare_supplied(model, ROOT/config["supplied_expansion"]["path"])
    amplitude = amplitude_checks(config, coefficients, model)
    joint = joint_checks(config, coefficients, model)
    state, matrices = linear_state_from_action(coefficients, model)
    linear = limited_linear_checks(state, matrices, config)
    mass_cases = []
    for c in (-.1,0.,.1):
        for signs in itertools.product((-1.,1.),repeat=3):
            q = np.zeros(7); q[3:6] = .1*np.array(signs)/math.sqrt(3); q[6]=c
            jet = rod.FieldJet(q, *(np.zeros(7) for _ in range(5)))
            mass = rod.quartic_mass_matrix(jet,coefficients,model)
            mass_cases.append(float(np.linalg.eigvalsh(mass).min()))
    mass = {"status": "PASS" if min(mass_cases)>0 else "FAIL", "sample_count": len(mass_cases),
            "minimum_eigenvalue": min(mass_cases), "scope": "sampled small neighborhood, not global truncated positivity"}
    result = {"fingerprint": fingerprint, "field_order": list(rod.FIELD_ORDER), "symbol_order": list(rod.SYMBOL_ORDER),
              "coefficients": coefficients.values(), "A_B_comparisons": comparisons,
              "exact_checks": checks, "supplied_comparison": supplied, "amplitude": amplitude,
              "joint": joint, "linear": linear, "linear_action_matrices": matrices, "mass": mass,
              "polynomials": {"T4": model.T4.serialize(), "V4": model.V4.serialize(),
                              "residuals_A": [p.serialize() for p in model.residual_a],
                              "residuals_B": [p.serialize() for p in model.residual_b],
                              "fluxes": [p.serialize() for p in model.flux_a]},
              "runtime_seconds": time.perf_counter()-started}
    statuses = {"NLSP_MODEL_SPECIFICATION":"PASS", "NLSP_VARIATIONAL_REDERIVATION":"PASS",
                "NLSP_SUPPLIED_EXPANSION_COMPARISON":supplied["status"],
                "NLSP_BOUNDARY_AND_JOINT_WORK":joint["status"],
                "NLSP_LINEAR_LIMIT":linear["status"],
                "NLSP_AXIAL_AND_PLANAR_LIMITS":"PASS" if all(checks[n]["status"]=="PASS" for n in ("pure_axial","planar")) else "FAIL",
                "NLSP_ENERGY_AND_SYMMETRY":"PASS" if all(v["status"]=="PASS" for v in checks.values()) and mass["status"]=="PASS" else "FAIL",
                "NLSP_AMPLITUDE_TRUNCATION_ORDER":amplitude["status"],
                "NLSP_STRAIGHT_SPLIT_CHECK":linear["split_status"]}
    statuses["NLSP_CUBIC_MODEL_AUDIT"] = "PASS" if all(v in ("PASS","MATCH") for v in statuses.values()) else "PARTIAL"
    result["statuses"] = statuses
    return result, appendix_text(model,fingerprint)


def linear_reference_controls(root, new_state_builder, new_joint_operator):
    """Callbacks: A14(omega), J14x28(beta_deg), q-first physical state.

    q=(u,w,v,Phi,psi,theta,c); efforts=(N,Qw,Qv,T,Mpsi,Mtheta,Rc).
    This function is a bounded diagnostic, not a new production solver.
    """
    from scripts.lib import mindlin_herrmann_longitudinal as mh
    from scripts.lib import yartsev_ch2_monoclinic_rod as book
    from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section
    root = Path(root)
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    section = rectangular_section(E=1., nu=.3, rho=1., width=.20, thickness=.05, K=5/6)
    model = mh.project_jang_reduced_rectangular(section)
    g = section.G
    material = book.BookMaterial(E1_real=1., E2_real=1., G12_real=g,
        G13_real=g, G23_real=g, nu12=.3, rho=1., eta1=0., eta2=0.,
        eta12=0., eta13=0., eta23=0.)
    # Book I_y=a**3*b/12 becomes new I_perp=h*b**3/12, not I_parallel.
    point = book.make_rod_point(0., geometry=book.Geometry(a=.20,b=.05,
        length=1.,shear_factor=5/6),material=material)
    ct = float(point.torsion.C_T.real)
    p = model.coefficients
    coefficients = {'m':p['m'],'jp':p['j'],'jb':point.material.rho*point.geometry.I_y,
        'C':p['C'],'H':p['H'],'S':p['S'],'Bp':p['B'],
        'Bb':float(point.properties.Ex.real)*point.geometry.I_y,'CT':ct,'nu':.3}
    embeddings = {'mh':[0,6,7,13],'inplane':[1,5,8,12],
                  'book_outplane_torsion':[2,4,3,9,11,10]}
    def old_matrix(omega):
        result = np.zeros((14,14))
        for key, a in (('mh',mh.harmonic_state_matrix(model,omega,'mh')),
                       ('inplane',mh.harmonic_state_matrix(model,omega,'timoshenko')),
                       ('book_outplane_torsion',book.state_matrix(omega,point).real)):
            ids=embeddings[key]
            result[np.ix_(ids,ids)]=a
        return result
    matrix_checks=[]
    for omega in (.317,2.,3.,9.):
        new=np.asarray(new_state_builder(omega),dtype=float)
        old=old_matrix(omega)
        difference=float(np.max(np.abs(new-old)/np.maximum(1.,np.abs(old))))
        if difference > 2e-13:
            raise ArithmeticError('New zero-Jacobian state matrix disagrees with old linear blocks')
        matrix_checks.append({'omega':omega,'max_scaled_coefficient_difference':difference})
    # Frozen saved inputs/profiles; no recalculation of the previous study.
    bundle=root/'results/mindlin_herrmann_timoshenko_single_rod/342ce44bff81c36f'
    if not bundle.is_dir():
        raise FileNotFoundError('Required saved single-rod reference unavailable; report PARTIAL')
    manifest=json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
    for name,digest in manifest['artifact_hashes'].items():
        if sha(bundle/name)!=digest:
            raise ValueError('Frozen single-rod reference artifact changed')
    direct=json.loads((bundle/'result.json').read_text(encoding='utf-8'))
    source_hashes={str((bundle/name).relative_to(root)):sha(bundle/name)
        for name in ('manifest.json','result.json','mh_modes.csv','timoshenko_modes.csv')}
    # Independent old Yartsev section-clamped fixed-fixed shooting, narrow interval.
    out_ids=[0,1,3,4]
    evaluations=0
    def book_out_bc(omega):
        nonlocal evaluations
        evaluations += 1
        a=book.scaled_state_matrix(omega,point)[np.ix_(out_ids,out_ids)]
        return expm(a)[:2,2:].real
    def out_det(omega):return float(np.linalg.det(book_out_bc(omega)))
    nodes=np.linspace(.001,2.70,181)
    samples=[out_det(w) for w in nodes]
    out_roots=[]
    for left,right,fl,fr in zip(nodes[:-1],nodes[1:],samples[:-1],samples[1:]):
        if fl*fr>=0:continue
        omega=brentq(out_det,left,right,xtol=1e-11,rtol=1e-12)
        sv=np.linalg.svd(book_out_bc(omega),compute_uv=False)
        out_roots.append({'omega':omega,'family':'outplane','bracket_omega':[float(left),float(right)],
                          'old_boundary_singular_ratio':float(sv[-1]/sv[0])})
    rotated=mh.project_jang_reduced_rectangular(rectangular_section(E=1.,nu=.3,rho=1.,
        width=.05,thickness=.20,K=5/6))
    out_count=mh.finite_count_upper_bound(rotated,1.,2.70,'timoshenko')
    # The lower relaxed-rotation spectrum bounds CC count. Below this interval
    # two old modes saturate the bound; no sign-scan completeness claim alone.
    if len(out_roots)!=out_count['upper_count']:
        raise ArithmeticError('Limited out-of-plane old reference count not certified')
    torsion_speed=math.sqrt(ct/(point.material.rho*point.geometry.I_p))
    torsion=[{'omega':n*math.pi*torsion_speed,'family':'torsion',
              'bracket_omega':[.98*n*math.pi*torsion_speed,1.02*n*math.pi*torsion_speed]}
             for n in range(1,4)]
    all_roots=[{'omega':r['omega'],'family':'mh','bracket_omega':r['bracket_omega']}
               for r in direct['mh']['roots'][:3]]
    all_roots += [{'omega':r['omega'],'family':'inplane','bracket_omega':r['bracket_omega']}
                  for r in direct['timoshenko']['roots'][:5]]
    all_roots+=out_roots+torsion
    prefix=sorted(all_roots,key=lambda r:r['omega'])[:7]
    ceiling=2.70
    in_count=mh.finite_count_upper_bound(model,1.,ceiling,'timoshenko')
    mh_count=mh.finite_count_upper_bound(model,1.,ceiling,'mh')
    scalar_count=math.floor(ceiling/(math.pi*torsion_speed))
    found=sum(r['omega']<=ceiling for r in all_roots)
    certificate=in_count['upper_count']+out_count['upper_count']+mh_count['upper_count']+scalar_count
    if found != certificate or found != 9:
        raise ArithmeticError('First6+guard7 full direct reference count not saturated')
    # New Jacobian shooting is independent of old bounded analytic finite basis.
    # Diagonal unit scaling followed by short-step expm / positive-diagonal QR.
    units=np.array([1.,1.,1.,1.,1.,1.,1.,p['C'],p['S'],p['S'],ct,
                    coefficients['Bb'],p['B'],p['H']])
    new_evaluations=0
    def end_frame(omega,length):
        nonlocal new_evaluations
        new_evaluations+=1
        a=np.asarray(new_state_builder(omega),dtype=float)
        a=a*units[None,:]/units[:,None]
        rate=float(np.max(np.abs(np.linalg.eigvals(a))))
        steps=max(1,math.ceil(rate*length))
        if steps>512:raise ArithmeticError('New linear shooting QR step budget exceeded')
        step=expm(a*(length/steps))
        frame=np.vstack((np.zeros((7,7)),np.eye(7)))
        for _ in range(steps):
            frame,r=np.linalg.qr(step@frame,mode='reduced')
            frame*=np.where(np.diag(r)>=0,1.,-1.)[None,:]
        return units[:,None]*frame
    def direct_bc(omega):
        result=end_frame(omega,1.)[:7]
        return result/units[:7,None]
    def split_bc(omega,split,beta=0.):
        basis=np.zeros((28,14))
        basis[:14,:7]=end_frame(omega,split)
        basis[14:,7:]=end_frame(omega,1-split)
        result=new_joint_operator(beta)@basis
        return result/units[:,None]
    def match(builder,reference):
        output=[]
        for r in reference:
            left,right=r['bracket_omega']
            # Isolate the saved root from other families in the full determinant.
            left=max(left, .9999*r['omega']);right=min(right,1.0001*r['omega'])
            f=lambda omega:float(np.linalg.det(builder(omega)))
            fl,fr=f(left),f(right)
            if fl*fr>0:raise ArithmeticError('New linear boundary matrix did not bracket old reference')
            omega=brentq(f,left,right,xtol=1e-11,rtol=1e-12)
            sv=np.linalg.svd(builder(omega),compute_uv=False)
            relative=abs(omega/r['omega']-1.)
            if relative>2e-8 or sv[-1]/sv[0]>1e-9:
                raise ArithmeticError(f'New linear spectrum/clamp failed: {r}, new={omega}, relative={relative}, sigma={sv[-1]/sv[0]}')
            output.append({'omega':omega,'frequency_hz':omega/(2*math.pi),
                'reference_omega':r['omega'],'family':r['family'],
                'relative_frequency_difference':relative,'boundary_singular_ratio':float(sv[-1]/sv[0]),
                'bracket_omega':[left,right]})
        return output
    direct_new=match(direct_bc,prefix)
    axial_new=match(direct_bc,[{'omega':r['omega'],'family':'mh',
        'bracket_omega':r['bracket_omega']} for r in direct['mh']['roots'][:3]])
    split_results=[]
    for split in (.50,.35):
        matches=match(lambda omega:split_bc(omega,split),prefix)
        axial_matches=match(lambda omega:split_bc(omega,split),[{'omega':r['omega'],'family':'mh','bracket_omega':r['bracket_omega']} for r in direct['mh']['roots'][:3]])
        split_results.append({'split':split,'roots':matches,'mh_family_first3':axial_matches,'status':'PASS'})
    # Kinematic and effort rows must be dual; theta/M is a vector component,
    # not the old planar scalar once the full 3D rotation is retained.
    structural=[]
    from scripts.lib.mindlin_herrmann_timoshenko_joint import frames
    P=np.diag([1.,-1.,1.])
    rng=np.random.default_rng(25741)
    for beta in (0.,45.,90.):
        transforms=[]
        for frame in frames(beta):
            t=np.r_[frame.t,0.];n=np.r_[frame.n,0.];k=np.cross(t,n)
            b=np.column_stack((t,n,k))
            t7=np.eye(7);t7[:3,:3]=b;t7[3:6,3:6]=b@P
            transforms.append(t7)
        joint=new_joint_operator(beta)
        expected=np.zeros((14,28))
        expected[:7,:7]=transforms[0];expected[:7,14:21]=-transforms[1]
        expected[7:,7:14]=transforms[0];expected[7:,21:28]=transforms[1]
        difference=float(np.max(np.abs(joint-expected)))
        rank=int(np.linalg.matrix_rank(joint))
        work=0.
        for t7 in transforms:
            for _ in range(12):
                f,d=rng.normal(size=(2,7))
                local=float(f@(t7.T@d));global_work=float((t7@f)@d)
                work=max(work,abs(local-global_work)/max(1.,abs(local),abs(global_work)))
        if rank!=14 or difference>2e-14 or work>2e-14:
            raise ArithmeticError('Full seven-field joint rank or force/moment duality failed')
        structural.append({'beta_deg':beta,'rank':rank,'matrix_difference':difference,'duality_max_error':work})
    # Only the in-plane sector has a saved same-section-rotation-clamp angular
    # reference. Historical Yartsev coupled frequencies use a different clamp.
    angular_bundle=root/'results/mindlin_herrmann_timoshenko_general_beta_joint/a36e72f715cb526b'
    angular=[]
    if angular_bundle.is_dir():
        angular_manifest=json.loads((angular_bundle/'manifest.json').read_text(encoding='utf-8'))
        for name,digest in angular_manifest['artifact_hashes'].items():
            if sha(angular_bundle/name)!=digest:
                raise ValueError('Frozen angular reference artifact changed')
        angular_result=json.loads((angular_bundle/'result.json').read_text(encoding='utf-8'))
        source_hashes[str((angular_bundle/'result.json').relative_to(root))]=sha(angular_bundle/'result.json')
        rows=[0,1,5,6,7,8,12,13]
        columns=[0,6,1,5,7,13,8,12]
        for beta in (45.,90.):
            case=next(c for c in angular_result['pilot'] if c['beta_deg']==beta)
            selected=[{'omega':r['omega'],'family':'coupled_inplane',
                       'bracket_omega':r['bracket_omega']} for r in case['roots'][:3]]
            matches=match(lambda omega:split_bc(omega,.5,beta)[np.ix_(rows,columns)],selected)
            angular.append({'beta_deg':beta,'sector':'inplane','roots':matches,'status':'PASS',
                            'reference':'frozen production general-beta first3 only; no angle sweep'})
    else:
        angular.append({'status':'UNAVAILABLE','reason':'Saved same-clamp angular reference absent'})
    # Profile checks use a spatial-eigenvector basis of the NEW Jacobian,
    # independently of the OLD PDE-based trigonometric/exponential basis.
    import csv
    profile_x=np.linspace(0.,1.,201)
    def spatial_basis(omega,ids,points,length=1.):
        sub=np.asarray(new_state_builder(omega))[np.ix_(ids,ids)]
        scale=units[ids]
        eig,vec=np.linalg.eig(sub*scale[None,:]/scale[:,None])
        columns=[]
        for i,r in enumerate(eig):
            physical=scale*vec[:,i]
            if r.imag>1e-8:
                value=np.exp(r*points[:,None])*physical[None,:]
                columns.extend((value.real,value.imag))
            elif abs(r.imag)<1e-8:
                anchor=length if r.real>0 else 0.
                value=np.exp(r.real*(points[:,None]-anchor))*physical.real[None,:]
                columns.append(value)
        if len(columns)!=len(ids):
            raise ArithmeticError('Unsupported bounded new linear spatial basis')
        return np.stack(columns,axis=2)
    def new_profile(omega,ids):
        displacement_count=len(ids)//2
        end=spatial_basis(omega,ids,np.array([0.,1.]))
        bc=np.concatenate((end[0,:displacement_count],end[1,:displacement_count]))
        bc/=np.linalg.norm(bc,axis=1)[:,None]
        _,_,right=np.linalg.svd(bc)
        return spatial_basis(omega,ids,profile_x)@right[-1]
    reversal=np.r_[[-1.,-1.,1.,-1.,-1.,1.,1.],[1.,1.,-1.,1.,1.,-1.,-1.]]
    def new_split_profile(omega,ids,split):
        n=len(ids);nq=n//2
        end1=spatial_basis(omega,ids,np.array([0.,split]),split)
        end2=spatial_basis(omega,ids,np.array([0.,1-split]),1-split)
        external=np.zeros((n,2*n))
        external[:nq,:n]=end1[0,:nq]
        external[nq:,n:]=end2[0,:nq]
        ends=np.zeros((28,2*n))
        ends[np.ix_(ids,range(n))]=end1[1]
        ends[np.ix_([i+14 for i in ids],range(n,2*n))]=end2[1]
        rows=list(ids[:nq])+[i+7 for i in ids[:nq]]
        j=new_joint_operator(0.)
        bc=np.concatenate((external,(j@ends)[rows]))
        bc/=np.linalg.norm(bc,axis=1)[:,None]
        _,_,right=np.linalg.svd(bc);coeff=right[-1]
        values=np.zeros((len(profile_x),n))
        left=profile_x<=split;right_mask=~left
        values[left]=spatial_basis(omega,ids,profile_x[left],split)@coeff[:n]
        values[right_mask]=(spatial_basis(omega,ids,1-profile_x[right_mask],1-split)@coeff[n:])*reversal[ids]
        endpoint_states=ends@coeff
        return values,(j@endpoint_states)[rows],(external@coeff)
    def saved_profile(block,index):
        filename='mh_modes.csv' if block=='mh' else 'timoshenko_modes.csv'
        with (bundle/filename).open(encoding='utf-8',newline='') as stream:
            rows=[r for r in csv.DictReader(stream) if int(r['family_index'])==index]
        if len(rows)!=201:
            raise ValueError('Expected saved 201-point physical reference profile')
        return np.array([[float(r[k]) for k in ('first_displacement','second_coordinate',
                        'first_resultant','second_resultant')] for r in rows])
    profile_checks=[]
    profiles=[]
    for family,indices,references in (
            ('mh',embeddings['mh'],direct['mh']['roots'][:3]),
            ('inplane',embeddings['inplane'],direct['timoshenko']['roots'][:3]),
            ('outplane',[2,4,9,11],out_roots),('torsion',[3,10],torsion[:2])):
        for index,r in enumerate(references,1):
            omega=r['omega']
            if family in ('mh','inplane'):
                reference=saved_profile(family,index)
            elif family=='outplane':
                bc=book_out_bc(omega)
                _,_,right=np.linalg.svd(bc)
                initial=np.r_[np.zeros(2),right[-1]]
                a=book.scaled_state_matrix(omega,point)[np.ix_(out_ids,out_ids)]
                physical_scale=np.array([1.,1.,coefficients['Bb'],coefficients['Bb']])
                reference=np.array([expm(a*x)@initial*physical_scale for x in profile_x]).real
            else:
                wave=index*math.pi
                reference=np.column_stack((np.sin(wave*profile_x),ct*wave*np.cos(wave*profile_x)))
            candidate=new_profile(omega,indices)
            dimensions=np.maximum(np.max(abs(reference),axis=0),1e-30)
            a=candidate/dimensions;b=reference/dimensions
            alpha=float(np.sum(a*b)/np.sum(a*a))
            candidate*=alpha
            component=np.linalg.norm(candidate-reference,axis=0)/np.maximum(np.linalg.norm(reference,axis=0),1e-30)
            displacement_count=len(indices)//2
            boundary=float(np.max(np.abs(candidate[[0,-1],:displacement_count])/dimensions[:displacement_count]))
            if max(component)>5e-8 or boundary>1e-9:
                raise ArithmeticError(f'New spatial-eigenbasis profile mismatch: {family}/{index}: {component}, bc={boundary}')
            for split_case in split_results:
                split=split_case['split']
                split_value,joint_value,external=new_split_profile(omega,indices,split)
                a=split_value/dimensions;b=reference/dimensions
                alpha=float(np.sum(a*b)/np.sum(a*a))
                split_value*=alpha
                component_split=np.linalg.norm(split_value-reference,axis=0)/np.maximum(np.linalg.norm(reference,axis=0),1e-30)
                joint_scaled=np.abs(alpha*joint_value)/dimensions
                external_scaled=np.abs(alpha*external)/np.tile(dimensions[:displacement_count],2)
                if max(component_split)>5e-8 or max(joint_scaled)>1e-9 or max(external_scaled)>1e-9:
                    raise ArithmeticError(f'Selected split profile/interface failed: {split}/{family}/{index}')
                split_case.setdefault('profile_comparison',[]).append({'family':family,'family_index':index,
                    'component_L2_relative':component_split.tolist(),'joint_scaled_residual':joint_scaled.tolist(),
                    'clamp_scaled_residual':float(max(external_scaled))})
                for x,value,old in zip(profile_x,split_value,reference):
                    split_case.setdefault('mode_profiles',[]).append({'family':family,'family_index':index,
                        'x':float(x),'global_direct_state':value.tolist(),'old_reference_state':old.tolist(),
                        'state_indices':indices})
            profile_checks.append({'family':family,'family_index':index,'omega':omega,
                'component_L2_relative':component.tolist(),'clamp_scaled_residual':boundary,
                'comparison':'new Jacobian spatial-eigenvector basis vs saved old PDE basis / old Yartsev expm / exact torsion'})
            for x,new,old in zip(profile_x,candidate,reference):
                profiles.append({'family':family,'family_index':index,'x':float(x),
                    'new_state':new.tolist(),'old_state':old.tolist(),'state_indices':indices})
    return {'status':'PASS','scope':'bounded linear limits, no nonlinear trajectory / no parameter sweep',
        'state_order':['u','w','v','Phi','psi','theta','c','N','Qw','Qv','T','Mpsi','Mtheta','Rc'],
        'clamp':'section_rotation_clamp: U=z=c=0; no centerline-slope constraint',
        'book_slope_clamp_reference_used':False,'coefficients':coefficients,
        'yartsev_mapping':{'a_book':.20,'b_book':.05,'Iy_book':'I_perp',
                          'I_perp':point.geometry.I_y,'Ip':point.geometry.I_p,
                          'Sbar16':float(point.properties.Sbar16.real),
                          'torsion_series_terms':point.torsion.terms_used,
                          'torsion_relative_tail_bound':point.torsion.estimated_relative_tail},
        'operator_comparison':matrix_checks,'saved_reference_hashes':source_hashes,
        'full_first6_guard7':direct_new,'mh_family_first3':axial_new,
        'outplane_reference':out_roots,'torsion_reference':torsion,
        'completeness':{'ceiling_omega':ceiling,'found_count':found,'upper_count':certificate,
            'inplane':in_count,'outplane':out_count,'mh':mh_count,'torsion_count':scalar_count},
        'straight_split':split_results,'split_status':'PASS','full_joint_structure':structural,
        'profile_comparison':profile_checks,'linear_profiles':profiles,'nonzero_angle_inplane':angular,
        'limited_old_reference_evaluations':evaluations,'new_linear_boundary_evaluations':new_evaluations,
        'nonzero_angle_old_same_clamp_reference':'UNAVAILABLE: historical Yartsev joint uses book_slope_clamp'}


def limited_linear_checks(state, matrices, config):
    result=linear_reference_controls(ROOT,state,joint_operator)
    result.setdefault("split_status","PASS")
    result["root_evaluations"]=result["limited_old_reference_evaluations"]+result["new_linear_boundary_evaluations"]
    return result


def compute(config_path=CONFIG, output=OUTPUT):
    config, inputs, fingerprint = identity(config_path)
    bundle = Path(output)/fingerprint
    manifest_path = bundle/"manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf8"))
        if manifest["identity"] != inputs:
            raise ValueError("Cache identity mismatch")
        for name, digest in manifest["artifacts"].items():
            if sha(bundle/name) != digest:
                raise ValueError("Cache artifact mismatch: "+name)
        return json.loads((bundle/"result.json").read_text(encoding="utf8")), bundle, {"cache_reused":True,"derivation_calls":0,"root_evaluations":0}
    bundle.mkdir(parents=True, exist_ok=True)
    try:
        result, appendix = run_audit(config,fingerprint)
        json_write(bundle/"result.json",result)
        (bundle/"expansion_generated.md").write_text(appendix,encoding="utf8")
        APPENDIX.write_text(appendix,encoding="utf8")
        manifest = {"identity":inputs, "config":config, "command":sys.argv,
                    "git":{name:subprocess.check_output(command,cwd=ROOT).decode().strip() for name,command in
                           (("HEAD",["git","rev-parse","HEAD"]),("branch",["git","branch","--show-current"]),
                            ("status",["git","status","--short"]),("diff_stat",["git","diff","--stat"]))},
                    "artifacts":{name:sha(bundle/name) for name in ("result.json","expansion_generated.md")},
                    "performance":{"runtime_seconds":result["runtime_seconds"],"derivation_calls":1,
                                   "root_evaluations":result["linear"].get("root_evaluations",0)}}
        json_write(manifest_path,manifest)
        json_write(Path(output)/"current.json",{"fingerprint":fingerprint})
        return result,bundle,{"cache_reused":False,"derivation_calls":1,"root_evaluations":manifest["performance"]["root_evaluations"]}
    except Exception as exc:
        json_write(bundle/"failure.json",{"type":type(exc).__name__,"reason":str(exc)})
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=CONFIG)
    parser.add_argument("--output-dir",type=Path,default=OUTPUT)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check-sources",action="store_true")
    modes.add_argument("--compute",action="store_true")
    modes.add_argument("--report-only",type=Path)
    args=parser.parse_args(argv)
    if args.check_sources:
        _,inputs,fingerprint=identity(args.config)
        print(json.dumps({"status":"PASS","fingerprint":fingerprint,"input_hashes":inputs["hashes"]},indent=2))
        return 0
    if args.report_only:
        config,inputs,_=identity(args.config)
        manifest=json.loads((args.report_only/"manifest.json").read_text(encoding="utf8"))
        if manifest["identity"]!=inputs:
            raise ValueError("Report-only input/code identity mismatch")
        for name,digest in manifest["artifacts"].items():
            if sha(args.report_only/name)!=digest:
                raise ValueError("Report-only artifact mismatch")
        result=json.loads((args.report_only/"result.json").read_text(encoding="utf8"))
        performance={"cache_reused":True,"root_evaluations":0,"derivation_calls":0}
        bundle=args.report_only
    else:
        result,bundle,performance=compute(args.config,args.output_dir)
    print(json.dumps({"bundle":str(bundle),"statuses":result["statuses"],"performance":performance},indent=2))
    return 0 if result["statuses"]["NLSP_CUBIC_MODEL_AUDIT"]=="PASS" else 2


if __name__=="__main__":
    raise SystemExit(main())
