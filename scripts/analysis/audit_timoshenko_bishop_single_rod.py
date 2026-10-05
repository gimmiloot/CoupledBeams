"""Exact, diagnostic-only audit of two candidate displacement fields.

This is a kinematics/constitutive compatibility workflow, not a spectral preset.
It reuses the existing rectangular section and Bishop Segment only for limits.
Local sparse section polynomials use Fraction; no CAS/dependency is installed.
No combined theory is selected and no roots or joint matrices are assembled.
See docs/theory/timoshenko_bishop_single_rod.md for the general manual proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
STATUS = "COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED"
# Existing G20 section, isotropic-limit note sections 6--8; exact algebra only.
PARAMETERS = dict(E=F(1), rho=F(1), nu=F(3, 10), b=F(1, 5), h=F(1, 20), kappa=F(5, 6))
SOURCES = {
    "docs/literature/pdf/s13370-014-0286-3.pdf": "19fec24bfb5dda35d27141a09661eaf68226e59fe08153a3f1a1c441551e5739",
    "docs/literature/pdf/Rayleigh_Love accepted version.pdf": "26db797736decf81428d3f990afac5a9b48a5ee60f06cfc31eb6b52b052fcd48",
}
REFERENCES = (
    "scripts/lib/bishop_longitudinal.py",
    "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
    "docs/theory/bishop_literature_reproduction.md",
    "docs/laminated_beams/reddy_four_ply_isotropic_limit_validation.md",
    "data/input/bishop_literature_sources.json",
    "tests/test_bishop_literature.py",
    "tests/test_timoshenko_bishop_single_rod.py",
    "docs/theory/timoshenko_bishop_single_rod.md",
)


# A section polynomial is {(power_y, power_z): exact coefficient}.
def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for powers, value in polynomial.items():
            result[powers] = result.get(powers, F(0)) + value
    return {key: value for key, value in result.items() if value}


def scale(polynomial, factor):
    return {key: value * factor for key, value in polynomial.items() if value * factor}


def multiply(a, b):
    result = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = (i + k, j + l)
            result[key] = result.get(key, F(0)) + x * y
    return {key: value for key, value in result.items() if value}


def section_derivative(polynomial, coordinate):
    result = {}
    for powers, value in polynomial.items():
        degree = powers[coordinate]
        if degree:
            target = list(powers)
            target[coordinate] -= 1
            result[tuple(target)] = degree * value
    return result


def moment(i, j, b, h, yc=F(0), zc=F(0)):
    """Exact rectangle integral; shifts are algebraic negative controls only."""
    return (((yc+b/2)**(i+1)-(yc-b/2)**(i+1))/(i+1)
            * ((zc+h/2)**(j+1)-(zc-h/2)**(j+1))/(j+1))


def integrate(polynomial, b, h, yc=F(0), zc=F(0)):
    return sum((value * moment(*powers, b, h, yc, zc)
                for powers, value in polynomial.items()), F(0))


# A linear field maps jets (name, x-derivative order, t-derivative order)
# to section polynomials. Only differentiation and quadratic energies are needed.
def field_add(*fields):
    keys = set().union(*(field.keys() for field in fields))
    return {key: value for key in keys
            if (value := add(*(field.get(key, {}) for field in fields)))}


def derivative(field, coordinate):
    if coordinate in ("x", "t"):
        return {(name, dx + (coordinate == "x"), dt + (coordinate == "t")): p
                for (name, dx, dt), p in field.items()}
    return {key: value for key, p in field.items()
            if (value := section_derivative(p, ("y", "z").index(coordinate)))}


def displacement(nu, b, h, full_poisson=False):
    """Candidates only: centroid contraction or compatible full-strain example."""
    one, y, z = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
    ux = {("u", 0, 0): one, ("psi", 0, 0): scale(z, -1)}
    uy = {("u", 1, 0): scale(y, -nu)}
    uz = {("w", 0, 0): one, ("u", 1, 0): scale(z, -nu)}
    if full_poisson:
        # r has zero section mean, so w remains the centroid displacement.
        c = (h*h-b*b)/12
        r = {(0, 2): F(1, 2), (2, 0): F(-1, 2), (0, 0): -c/2}
        uy[("psi", 1, 0)] = {(1, 1): nu}
        uz[("psi", 1, 0)] = scale(r, nu)
    return ux, uy, uz


def strains(fields):
    ux, uy, uz = fields
    return (
        derivative(ux, "x"), derivative(uy, "y"), derivative(uz, "z"),
        field_add(derivative(uy, "z"), derivative(uz, "y")),  # gamma_yz
        field_add(derivative(ux, "z"), derivative(uz, "x")),  # gamma_xz
        field_add(derivative(ux, "y"), derivative(uy, "x")),  # gamma_xy
    )


def constitutive(E, nu):
    G = E/(2*(1+nu))
    lame = E*nu/((1+nu)*(1-2*nu))
    return [[(lame + (2*G if i == j else 0)) if i < 3 and j < 3
             else (G if i == j else F(0)) for j in range(6)] for i in range(6)]


def quadratic_hessian(fields, metric):
    """H in density = 1/2 sum H_ij jet_i jet_j; retain section polynomials."""
    result = {}
    for i, left in enumerate(fields):
        for j, right in enumerate(fields):
            for a, pa in left.items():
                for b, pb in right.items():
                    key = a, b
                    result[key] = add(result.get(key, {}), scale(multiply(pa, pb), metric[i][j]))
    return {key: value for key, value in result.items() if value}


def candidate_hessians(E, rho, nu, b, h, full_poisson=False):
    fields = displacement(nu, b, h, full_poisson)
    velocities = [derivative(field, "t") for field in fields]
    mass = [[rho if i == j else F(0) for j in range(3)] for i in range(3)]
    return quadratic_hessian(velocities, mass), quadratic_hessian(strains(fields), constitutive(E, nu))


def integrated(hessian, b, h, yc=F(0), zc=F(0)):
    return {key: integrate(value, b, h, yc, zc) for key, value in hessian.items()}


def axial_bending(key):
    return (key[0][0] == "u") != (key[1][0] == "u")


def label(jet):
    name, dx, dt = jet
    return name + ("_" + "x"*dx + "t"*dt if dx+dt else "")


def records(hessian, b, h, zc=F(0)):
    return [{"left": label(a), "right": label(bjet),
             "density_hessian": [{"y_power": i, "z_power": j, "coefficient": str(value)}
                                 for (i, j), value in sorted(p.items())],
             "section_hessian": str(integrate(p, b, h, zc=zc)),
             "axial_bending": axial_bending((a, bjet))}
            for (a, bjet), p in sorted(hessian.items()) if a <= bjet]


def audit():
    p = PARAMETERS
    b, h = p["b"], p["h"]
    inputs = {key: p[key] for key in ("E", "rho", "nu", "b", "h")}
    cases = {}
    for name, full in (("centroid_contraction_raw_3d", False), ("full_poisson_example_raw_3d", True)):
        mass, stiffness = candidate_hessians(**inputs, full_poisson=full)
        cases[name] = {"mass": records(mass, b, h), "stiffness": records(stiffness, b, h)}
    raw = candidate_hessians(**inputs)
    shifted = {"z_center": "1/100", "interpretation": "coordinate negative control only",
               "mass": records(raw[0], b, h, F(1, 100)),
               "stiffness": records(raw[1], b, h, F(1, 100))}
    moments = {f"y{i}_z{j}": str(moment(i, j, b, h))
               for i, j in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2))}
    from scripts.lib import bishop_longitudinal as bishop
    from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as timo
    section = timo.rectangular_section(E=float(p["E"]), rho=float(p["rho"]), nu=float(p["nu"]),
                                      width=float(b), thickness=float(h), K=float(p["kappa"]))
    polar = moment(2, 0, b, h)+moment(0, 2, b, h)
    J, H = p["nu"]**2*p["rho"]*polar, p["nu"]**2*p["E"]/(2*(1+p["nu"]))*polar
    segment = bishop.Segment(1., section.EA, section.rhoA, float(H), float(J))
    ratio = (1-p["nu"])/((1+p["nu"])*(1-2*p["nu"]))
    return {
        "scientific_status": STATUS,
        "parameters_exact": {key: str(value) for key, value in p.items()},
        "parameter_origin": "G20 section, isotropic-limit note sections 6--8; no frequency run",
        "moments_exact": moments, "candidates": cases, "negative_control": shifted,
        "centered_mixed_coefficients_zero": all(
            row["section_hessian"] == "0" for case in cases.values()
            for matrix in case.values() for row in matrix if row["axial_bending"]),
        "raw_to_timoshenko_bending_rigidity_ratio": str(ratio),
        "raw_to_timoshenko_shear_rigidity_ratio": str(1/p["kappa"]),
        "pure_bending_gate": "MISMATCH_FOR_RAW_FIELD; relaxed hybrid not selected",
        "bishop_pure_axial_coefficients": {"EA": segment.EA, "m": segment.m, "H": segment.H, "J": segment.J},
        "timoshenko_reference_coefficients": {"EI": section.EI, "KGA": section.KGA, "rhoI": section.rhoI},
        "full_poisson_extra_section_moment": str(b*h*(b**4+5*b*b*h*h+h**4)/720),
        "spectrum_and_factorization": "NOT_RUN_KINEMATICS_HARD_GATE",
        "missing_local_source": "Yucel--Arpaci--Tufekci 2014: already listed as unavailable; not used",
    }


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true", required=True, help="exact algebra only; no spectra")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/timoshenko_bishop_single_rod")
    args = parser.parse_args(argv)
    for path, expected in SOURCES.items():
        if sha(ROOT/path) != expected:
            raise ValueError("Local source changed; repeat the source audit: " + path)
    result = audit()
    def git(*options):
        return subprocess.check_output(["git", *options], cwd=ROOT, text=True, encoding="utf-8").strip()
    identity = {"schema": "timoshenko-bishop-kinematics-audit-v1", "sources": SOURCES,
                "files": {name: sha(ROOT/name) for name in (*REFERENCES, str(Path(__file__).relative_to(ROOT)).replace("\\", "/"))},
                "head": git("rev-parse", "HEAD"), "python": platform.python_version(),
                "versions": {name: version(name) for name in ("numpy", "scipy")},
                "parameters": result["parameters_exact"]}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    directory = args.output_dir/fingerprint[:16]
    directory.mkdir(parents=True, exist_ok=True)
    report = directory/"audit.json"
    report.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    manifest = {"identity": identity, "fingerprint": fingerprint, "git_branch": git("branch", "--show-current"),
                "git_status": git("status", "--short"), "executable": sys.executable,
                "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
                "arithmetic": "fractions.Fraction; exact rational polynomial integration",
                "cache_policy": "always recompute; content-addressed output, never read a cache",
                "artifacts": {"audit.json": sha(report)}}
    (directory/"manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (args.output_dir/"current.json").write_text(json.dumps({"directory": str(directory.resolve()), "fingerprint": fingerprint})+"\n", encoding="utf-8")
    print(STATUS)
    print(directory)


if __name__ == "__main__":
    main()
