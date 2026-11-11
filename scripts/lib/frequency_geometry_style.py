"""Geometry-only line styling for saved RL/Timoshenko frequency curves.

Lengths and section dimensions come from the canonical model. No matrices,
spectra, mode fields or tracking are evaluated. Inserted drawing vertices are
linear points on saved line segments, never additional spectral observations.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import textwrap

THRESHOLD = 10.0
CAPTION = ("Сплошными линиями обозначены участки, на которых для обоих стержней "
           "выполняется условие $L_i/h_i\\geq10$; пунктирными — участки, на которых "
           "это условие нарушается хотя бы для одного стержня")
AXES = {"BETA": "beta_deg", "BETA_DEG": "beta_deg", "MU": "mu",
        "TAU": "tau", "EPS": "epsilon", "EPSILON": "epsilon"}
ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("scripts/lib/frequency_geometry_style.py", "scripts/lib/rl_full_rr_aspect.py",
           "scripts/lib/isotropic_two_springs_thickness_asymmetry.py",
           "scripts/lib/isotropic_two_springs_rl_timoshenko.py",
           "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py")


def provenance():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def classify(r1, r2):
    """Inclusive, unrounded comparison; branch identity is not an input."""
    r1, r2 = float(r1), float(r2)
    if not all(math.isfinite(r) and r > 0 for r in (r1, r2)):
        raise ValueError("Positive finite individual slenderness ratios required")
    chi = min(r1, r2)
    passed = chi >= THRESHOLD
    return dict(L1_over_h1=r1, L2_over_h2=r2, chi=chi,
                slenderness_criterion_pass=passed, slenderness_threshold=THRESHOLD,
                line_style="solid" if passed else "dashed")


def geometry_node(*, beta_deg=0., mu=0., tau=0., epsilon=.005, aspect=1.,
                  E=1., rho=1., l=1., nu=.3, ks=5/6, variant="RL_FULL"):
    from scripts.lib.rl_full_rr_aspect import make_model
    from scripts.lib.isotropic_two_springs_thickness_asymmetry import Case
    model = make_model(mu=mu, tau=tau, epsilon=epsilon, aspect=aspect,
                       E=E, rho=rho, l=l, nu=nu, ks=ks, variant=variant)
    case = Case(beta=math.radians(beta_deg), mu=mu, tau=tau,
                k_joint="RIGID", k_end="RIGID")
    lengths = case.lengths(model)
    # SectionProperties already exposes h: do not reimplement its formula.
    heights = [arm.section.thickness for arm in model.arms]
    return dict(beta_deg=float(beta_deg), mu=float(mu), tau=float(tau),
                epsilon=float(epsilon), aspect=float(aspect),
                **classify(*(length/h for length, h in zip(lengths, heights))))


class FrequencyGeometry:
    """One geometric path, shared by every frequency on one parameter grid."""

    def __init__(self, x, *, axis, fixed, material):
        self.axis = AXES[axis.upper()]
        self.x = tuple(float(v) for v in x)
        if not self.x or not all(math.isfinite(v) for v in self.x):
            raise ValueError("Nonempty finite parameter grid required")
        if len(self.x) > 1 and not (all(a < b for a, b in zip(self.x, self.x[1:]))
                                    or all(a > b for a, b in zip(self.x, self.x[1:]))):
            raise ValueError("Parameter grid must be strictly monotone")
        self.parameters = dict(material, **fixed)
        self.parameters.pop(self.axis, None)
        self.nodes = [dict(x=v, **self.at(v)) for v in self.x]
        self.boundaries = self._boundaries()
        self.segments = []
        for i, (a, b) in enumerate(zip(self.x, self.x[1:])):
            inner = [v for v in self.boundaries if min(a, b) < v < max(a, b)]
            cuts = [a, *sorted(inner, reverse=b < a), b]
            for left, right in zip(cuts, cuts[1:]):
                self.segments.append(dict(index=i, left=left, right=right,
                    line_style=self.at((left+right)/2)["line_style"]))

    @classmethod
    def from_spec(cls, spec, material, x=None):
        return cls(spec["grid"] if x is None else x, axis=spec["mode"],
                   fixed=spec["fixed"], material=material)

    def at(self, value):
        return geometry_node(**dict(self.parameters, **{self.axis: value}))

    def _boundaries(self):
        """Zeros of each individual ratio, using canonical geometry only.

        In tau each ratio is monotone; in epsilon it is inverse-proportional;
        beta changes neither. For mu, differentiate L_i/h_i: the possible
        stationary abscissae are (tau-1-tau²)/(3 tau) for arm 1 and
        -(1+tau+tau²)/(3 tau) for arm 2. Splitting there brackets every zero,
        including two crossings with equal endpoint classifications. This is
        a geometric derivative, not a second implementation of the geometry.
        """
        if self.axis == "beta_deg" or len(self.x) == 1:
            return []
        low, high = min(self.x), max(self.x)
        roots = []
        for arm in (1, 2):
            cuts = [low, high]
            tau = self.parameters.get("tau", 0.)
            if self.axis == "mu" and tau:
                critical = ((tau-1-tau*tau) if arm == 1 else -(1+tau+tau*tau))/(3*tau)
                if low < critical < high:
                    cuts.append(critical)
            cuts.sort()
            def residual(v):
                return self.at(v)[f"L{arm}_over_h{arm}"]-THRESHOLD
            for a, b in zip(cuts, cuts[1:]):
                fa, fb = residual(a), residual(b)
                if fa == 0:
                    roots.append(a)
                if fb == 0:
                    roots.append(b)
                if (fa < 0) == (fb < 0) or fa == 0 or fb == 0:
                    continue
                # Bisection of an analytic geometric function, no frequency input.
                for _ in range(60):
                    mid = (a+b)/2
                    if mid in (a, b):
                        break
                    fm = residual(mid)
                    if fm == 0:
                        a = b = mid
                        break
                    if (fa < 0) == (fm < 0):
                        a, fa = mid, fm
                    else:
                        b = mid
                roots.append((a+b)/2)
        # Preserve an existing node when the geometric root coincides within
        # floating-point roundoff. Never insert it into the source grid.
        snapped = []
        for root in sorted(roots):
            if not math.isclose(self.at(root)["chi"], THRESHOLD, rel_tol=1e-12):
                continue  # The other arm still fails: this is not chi=10.
            node = min(self.x, key=lambda x: abs(x-root))
            if abs(root-node) <= 4*max(math.ulp(root), math.ulp(node)):
                root = node
            if not snapped or abs(root-snapped[-1]) > 4*max(math.ulp(root), math.ulp(snapped[-1])):
                snapped.append(root)
        return snapped

    def drawing_runs(self, y):
        """Return display vertices; keep NaN gaps and original sample markers."""
        y = tuple(float(v) for v in y)
        if len(y) != len(self.x):
            raise ValueError("Frequency grid does not match geometry grid")
        if len(y) == 1:
            return [(self.nodes[0]["line_style"], list(self.x), list(y))]
        runs = []
        for s in self.segments:
            i, a, b = s["index"], s["left"], s["right"]
            def display_value(v):
                if v == self.x[i]:
                    return y[i]
                if v == self.x[i+1]:
                    return y[i+1]
                return y[i]+(y[i+1]-y[i])*(v-self.x[i])/(self.x[i+1]-self.x[i])
            ya, yb = display_value(a), display_value(b)
            if runs and runs[-1][0] == s["line_style"] and runs[-1][1][-1] == a:
                runs[-1][1].append(b); runs[-1][2].append(yb)
            else:
                runs.append((s["line_style"], [a, b], [ya, yb]))
        return runs

    def diagnostic(self):
        return dict(axis=self.axis, fixed=self.parameters, nodes=self.nodes,
                    geometric_boundaries=self.boundaries,
                    boundary_vertices_are_spectral_samples=False,
                    display_frequency_at_boundary="linear on original segment only")


def plot_frequency_curve(ax, geometry, y, *, color, label=None, **style):
    """Only geometry controls the dash pattern; only original nodes get markers."""
    if any(k in style for k in ("linestyle", "ls", "dashes", "markevery")):
        raise ValueError("Main-frequency linestyle is reserved for geometry")
    artists = []
    original = {node["x"]: node["line_style"] for node in geometry.nodes}
    for status, x, values in geometry.drawing_runs(y):
        line, = ax.plot(x, values, color=color, label=label if not artists else "_nolegend_",
                        linestyle="-" if status == "solid" else "--",
                        markevery=[i for i, v in enumerate(x) if original.get(v) == status], **style)
        line._slenderness_line_style = status
        artists.append(line)
    return artists


def add_geometry_caption(fig):
    """Caption only if both styles occur on displayed main frequency curves."""
    styles = set()
    for ax in fig.axes:
        low, high = sorted(ax.get_xlim())
        for line in ax.lines:
            status = getattr(line, "_slenderness_line_style", None)
            xy = list(zip(line.get_xdata(), line.get_ydata()))
            visible = any(math.isfinite(a[1]) and math.isfinite(b[1]) and
                          max(min(a[0], b[0]), low) < min(max(a[0], b[0]), high)
                          for a, b in zip(xy, xy[1:]))
            if status and visible:
                styles.add(status)
    caption = CAPTION if styles == {"solid", "dashed"} else ""
    if caption:
        fig.text(.5, -.035, textwrap.fill(caption, width=max(65, int(fig.get_figwidth()*14))),
                 ha="center", va="top", fontsize=8)
    return caption


def write_geometry_diagnostics(output, geometries):
    """Write separate derived artifacts; one row per geometric node, not branch."""
    output = Path(output)
    rows = [dict(series=name, **node) for name, g in geometries.items() for node in g.nodes]
    if not rows:
        raise ValueError("No geometric nodes to record")
    with (output/"geometry_nodes.csv").open("x", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    result = dict(schema="frequency-slenderness-style-v1", threshold=THRESHOLD,
        equality="solid", geometry_sources_sha256=provenance(), spectral_calls=0,
        interpretation="preliminary geometric thin-rod indicator; not a theory validity boundary",
        partial_frequency_style="neutral color, dotted line and explicit CS/CC label; never main-branch dashed style",
        mixed_style_caption=CAPTION, series={name: g.diagnostic() for name, g in geometries.items()})
    (output/"geometry_diagnostics.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def plot_partial_frequency(ax, x, y, *, label, color="0.45", **style):
    """Auxiliary CS/CC spectra have a distinct, explicitly labelled dotted style."""
    return ax.plot(x, y, color=color, linestyle=":", label=label, **style)
