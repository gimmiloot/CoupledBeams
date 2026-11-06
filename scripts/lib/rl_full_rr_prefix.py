"""Bounded spectral-prefix adapter for the frozen homogeneous RL_FULL/RR model.

No new physical operator. Event/slot multiplicity extends the old count-guided
prefix without changing its source or fingerprint. Not a general solver API.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import math
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from scripts.lib import isotropic_two_springs_rl_timoshenko as physics
from scripts.lib.isotropic_two_springs_thickness_asymmetry import energies

VERSION = "rl-full-rr-prefix-v1"


class PrefixFailure(ArithmeticError):
    def __init__(self, code, **details):
        self.code, self.details = code, details
        super().__init__(code)


def digest(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True, allow_nan=False,
                                     separators=(",", ":")).encode()).hexdigest()


def point_contract(model, case, k_plot, k_guard, policy):
    if model.variant != "RL_FULL" or (case.k_joint, case.k_end) != (physics.RIGID,)*2:
        raise ValueError("Only homogeneous RL_FULL and exact RIGID/RIGID")
    if not isinstance(k_plot, int) or not isinstance(k_guard, int) or not 1 <= k_plot < k_guard:
        raise ValueError("Positive integer ranks K_plot < K_guard required")
    return dict(version=VERSION, model=asdict(model), case=asdict(case), tau=0.,
                K_plot=k_plot, K_guard=k_guard, policy=policy,
                spectrum_semantics="sorted_positions", state_order=["u", "w", "psi", "N", "Q", "M"],
                section=dict(EA=model.section.EA, EI=model.section.EI, KGA=model.section.KGA,
                             rhoA=model.section.rhoA, rhoI=model.section.rhoI, J=model.J),
                lengths=list(case.lengths(model)), time_reference=model.time,
                normalization="Omega=omega*t_ref=Lambda**2", coordinates="outer clamp to joint")


def singular_info(matrix, omega, nullity_tol):
    b = np.asarray(matrix(omega), dtype=float)
    if not np.all(np.isfinite(b)):
        raise PrefixFailure("NONFINITE_CHARACTERISTIC", omega=omega)
    u, s, vh = np.linalg.svd(b)
    scale = max(float(s[0]), np.finfo(float).tiny)
    return u, s, vh, int(np.sum(s <= nullity_tol*scale))


def counted_prefix(matrix, counter, lower, upper, k_guard, policy, neighbour=None):
    """Count-guided isolation of events; only matrix/count injection for tests.

    The physical caller supplies the unchanged balanced 12x12 boundary and
    checked existing energetic counter. No global sign scan or root dedup.
    """
    samples, events, retries, clusters = [], [], [], []
    target = k_guard
    report = dict(status="RUNNING", samples=samples, events=events,
                  query_retries=retries, cluster_trials=clusters, predictor="ABSENT")

    def query(w, interval=None):
        if len(samples) >= policy["max_count_queries"]:
            raise PrefixFailure("COUNT_BUDGET_EXHAUSTED", requested=w)
        entry = dict(requested=float(w))
        samples.append(entry)
        try:
            effective, count, diagnostic = counter(float(w))
            entry.update(diagnostic, effective=float(effective), count=int(count), status="PASS")
            if not lower <= effective <= upper:
                raise PrefixFailure("COUNT_OUTSIDE_SEARCH_DOMAIN", requested=w, effective=effective)
            if interval is not None and not interval[0] < effective < interval[1]:
                raise PrefixFailure("COUNT_OUTSIDE_INTERVAL", interval=list(interval), effective=effective)
            return float(effective), int(count)
        except PrefixFailure as error:
            entry.update(status=error.code, details=error.details)
            raise
        except (ArithmeticError, ValueError, np.linalg.LinAlgError) as error:
            entry.update(status="COUNT_UNRESOLVED", reason=str(error))
            raise PrefixFailure("COUNT_UNRESOLVED", requested=w, reason=str(error)) from error

    def interior(a, b):
        # If a trial hits an eigenvalue, choose at most two other interior points.
        for fraction in (.5, .3819660112501051, .6180339887498949):
            try:
                return query(a+(b-a)*fraction, (a, b))
            except PrefixFailure as error:
                if error.code != "COUNT_UNRESOLVED":
                    raise
                retries.append(dict(interval=[a, b], fraction=fraction, code=error.code))
        raise PrefixFailure("COUNT_UNRESOLVED", interval=[a, b])

    def determinant(w):
        return float(np.linalg.det(matrix(w)))

    def accept(w, multiplicity, a, b, ca, cb, method):
        _, s, _, nullity = singular_info(matrix, w, policy["nullity_tol"])
        if nullity != multiplicity:
            raise PrefixFailure("UNRESOLVED_NULLITY", omega=w, expected=multiplicity,
                                observed=nullity, singular_values=s.tolist())
        nonzero = float(s[0]/s[-multiplicity-1]) if multiplicity < len(s) else 1.
        if nonzero > policy["nonzero_condition_max"]:
            raise PrefixFailure("UNRESOLVED_COMPLEMENT", omega=w, condition=nonzero)
        events.append(dict(omega=float(w), multiplicity=multiplicity,
                           first_rank=ca+1, last_rank=cb, bracket=[a, b], counts=[ca, cb],
                           method=method, nullity=nullity, nonzero_condition=nonzero,
                           singular_values=s.tolist()))

    def repeated(a, b, ca, cb):
        # Frozen left/right SVD projections locate each null direction locally.
        # Different projected zeros are NOT averaged into a multiple eigenvalue.
        m = cb-ca
        u, s, vh, _ = singular_info(matrix, (a+b)/2, policy["nullity_tol"])
        if m > len(s):
            raise PrefixFailure("UNRESOLVED_CLUSTER", counts=[ca, cb], bracket=[a, b])
        candidates = []
        for i in range(1, m+1):
            def scalar(w):
                return float(u[:, -i]@matrix(w)@vh[-i])
            if scalar(a)*scalar(b) > 0:
                raise PrefixFailure("UNRESOLVED_CLUSTER", counts=[ca, cb], bracket=[a, b])
            candidates.append(float(brentq(scalar, a, b, xtol=policy["root_xtol"],
                                            rtol=policy["root_rtol"])))
        width = policy["root_xtol"]+policy["root_rtol"]*max(abs(x) for x in candidates)
        trial = dict(bracket=[a, b], counts=[ca, cb], projected_zeros=candidates,
                     resolution=width)
        clusters.append(trial)
        # A resolvable difference between projected zeros is evidence against
        # a common eigenvalue, even when it is smaller than the root tolerance.
        # Conservatively stop rather than converting that pair to multiplicity.
        if len(set(candidates)) != 1:
            raise PrefixFailure("UNRESOLVED_CLUSTER", **trial)
        w = candidates[0]
        # A numerical multiplicity needs the count jump inside a root-resolution
        # neighbourhood as well as independent null directions. No gap-based merge.
        radius = 8*width
        left, nl = query(w-radius, (a, b))
        right, nr = query(w+radius, (a, b))
        if (nl, nr) != (ca, cb):
            raise PrefixFailure("UNRESOLVED_CLUSTER", **trial, local_counts=[nl, nr])
        accept(w, m, left, right, ca, cb, "count_jump_and_full_nullspace")

    def isolate(a, b, ca, cb, depth=0):
        if cb < ca:
            raise PrefixFailure("COUNT_NOT_MONOTONE", bracket=[a, b], counts=[ca, cb])
        if cb == ca or ca >= target:
            return
        if cb-ca == 1 and determinant(a)*determinant(b) < 0:
            w = brentq(determinant, a, b, xtol=policy["root_xtol"], rtol=policy["root_rtol"])
            accept(w, 1, a, b, ca, cb, "count_one_brent")
            return
        if depth >= policy["max_depth"]:
            if cb-ca > 1:
                repeated(a, b, ca, cb)
                return
            raise PrefixFailure("ISOLATION_DEPTH_EXHAUSTED", bracket=[a, b], counts=[ca, cb])
        mid, cm = interior(a, b)
        if not ca <= cm <= cb:
            raise PrefixFailure("COUNT_NOT_MONOTONE", counts=[ca, cm, cb])
        isolate(a, mid, ca, cm, depth+1)
        isolate(mid, b, cm, cb, depth+1)

    try:
        a, ca = query(lower)
        b, cb = query(upper)
        if ca != 0:
            raise PrefixFailure("NONZERO_LOWER_COUNT", count=ca, omega=a)
        if cb < target:
            raise PrefixFailure("INSUFFICIENT_RIGHT_COUNT", count=cb, required=target, omega=b)
        nodes = [(a, ca)]
        if neighbour is not None:
            seed = np.asarray(neighbour, dtype=float)
            valid = (len(seed) >= k_guard and np.all(np.isfinite(seed)) and
                     np.all(seed > 0) and np.all(np.diff(seed) >= 0))
            if valid:
                for w in sorted(set((seed[:k_guard-1]+seed[1:k_guard])/2)):
                    if a < w < b:
                        node = query(float(w), (a, b))
                        if node[0] <= nodes[-1][0] or node[1] < nodes[-1][1]:
                            raise PrefixFailure("INVALID_SHIFTED_PREDICTOR", nodes=nodes, next=node)
                        nodes.append(node)
                report["predictor"] = "SORTED_SEED_ONLY"
            else:
                report["predictor"] = "REJECTED_INVALID_SEED_USED_COUNT_ONLY"
        nodes.append((b, cb))
        for (lo, cl), (hi, ch) in zip(nodes[:-1], nodes[1:]):
            isolate(lo, hi, cl, ch)
        # Complete a nearby cluster across the guard by a bounded count trigger.
        for _ in range(policy["max_cluster_extensions"]+1):
            end = events[-1]
            if end["last_rank"]-k_guard > policy["max_cluster_extra_ranks"]:
                raise PrefixFailure("UNRESOLVED_GUARD_CLUSTER", reason="extra-rank budget",
                                    required=end["last_rank"], guard=k_guard)
            guard_right = end["omega"]+policy["cluster_probe_relative"]*max(1., abs(end["omega"]))
            if guard_right >= b:
                raise PrefixFailure("GUARD_NEAR_DOMAIN_LIMIT", guard=end["omega"], right=b)
            wr, nr = query(guard_right)
            if nr == end["last_rank"]:
                break
            if nr < end["last_rank"]:
                raise PrefixFailure("COUNT_NOT_MONOTONE", guard=end, right_count=nr)
            if nr-k_guard > policy["max_cluster_extra_ranks"]:
                raise PrefixFailure("UNRESOLVED_GUARD_CLUSTER", required=nr)
            target = nr
            lo, nl = end["bracket"][1], end["last_rank"]
            if not lo < wr:
                raise PrefixFailure("UNRESOLVED_GUARD_CLUSTER", interval=[lo, wr])
            isolate(lo, wr, nl, nr)
        else:
            raise PrefixFailure("UNRESOLVED_GUARD_CLUSTER", reason="extension budget")
        ranks = [rank for e in events for rank in range(e["first_rank"], e["last_rank"]+1)]
        if ranks != list(range(1, events[-1]["last_rank"]+1)) or len(ranks) < k_guard:
            raise PrefixFailure("INCOMPLETE_PREFIX", ranks=ranks)
        if any(x["omega"] >= y["omega"] for x, y in zip(events[:-1], events[1:])):
            raise PrefixFailure("UNRESOLVED_EVENT_ORDER")
        report.update(status="PASS", lower_count=ca, upper_count=cb,
                      requested_guard=k_guard, computed_ranks=len(ranks),
                      guard_right_omega=wr, guard_right_count=nr,
                      unresolved_intervals=[], spectrum_semantics="sorted_positions")
        return events, report
    except PrefixFailure as error:
        report.update(status=error.code, failure=error.details)
        error.details = dict(error.details, search=report)
        raise


class PhysicalCounter:
    """Preflight the exact pole bypass before calling the frozen count_at."""
    def __init__(self, model, case, catalogs, policy):
        self.model, self.case, self.catalogs, self.policy = model, case, catalogs, policy
        s = model.section
        self.rl_limit = math.sqrt(s.EA/model.J) if model.J else math.inf
        self.tim_limit = math.sqrt(s.KGA/s.rhoI)
        self.domain = min(self.rl_limit, self.tim_limit, *(c["upper"] for c in catalogs.values()))
        self.poles = []
        for L in case.lengths(model):
            cat = catalogs[L]
            d = s.EA-model.J*cat["upper"]**2
            if cat["upper"] >= min(self.rl_limit, self.tim_limit) or d <= 0:
                raise PrefixFailure("CATALOG_OUTSIDE_BASIS_DOMAIN", length=L, upper=cat["upper"])
            nmax = math.floor(cat["upper"]*math.sqrt(s.rhoA/d)*L/math.pi)
            # Exactly the existing count_at pole formula, including repeated poles.
            self.poles.extend(cat["poles"])
            self.poles.extend(math.sqrt(s.EA*(n*math.pi/L)**2/(s.rhoA+model.J*(n*math.pi/L)**2))
                              for n in range(1, nmax+1))

    def __call__(self, omega):
        effective = float(omega)
        exclusion = self.policy["pole_exclusion"]*max(1., effective)
        for pole in sorted(self.poles):
            if abs(effective-pole) <= exclusion:
                effective = pole+2*exclusion
        if not 0 < omega < self.domain or not 0 < effective < self.domain:
            raise PrefixFailure("COUNT_OUTSIDE_CATALOG", requested=omega, effective=effective, domain=self.domain)
        actual, count, details = physics.count_at(self.model, self.case, omega, self.catalogs, self.policy)
        if actual != effective:
            raise PrefixFailure("POLE_BYPASS_CONTRACT_MISMATCH", predicted=effective, actual=actual)
        return actual, count, details


def fields(model, case, omega, coefficients, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    arms = []
    for L, c in zip(case.lengths(model), coefficients):
        x = (nodes+1)*L/2
        arms.append(dict(x=x, weights=weights*L/2,
                         y=physics.basis(model, L, omega, x)@c,
                         dy=physics.basis(model, L, omega, x, 1)@c,
                         ends=physics.basis(model, L, omega, [0., L])@c))
    return arms


def recovered_modes(model, case, event, order, policy):
    """Use old simple-mode recovery; normalize the whole kernel for multiplicity."""
    w, mult = event["omega"], event["multiplicity"]
    if mult == 1:
        return [physics.mode(model, case, w, order)]
    b, _, col = physics.balance(physics.boundary(model, case, w))
    _, s, vh = np.linalg.svd(b)
    if int(np.sum(s <= s[0]*policy["nullity_tol"])) != mult:
        raise PrefixFailure("UNRESOLVED_NULLITY", omega=w)
    coeff = (vh[-mult:]*col).reshape(mult, 2, 6)
    trial = [dict(arms=fields(model, case, w, c, order)) for c in coeff]
    gram = np.array([[physics.inner(model, a, b) for b in trial] for a in trial])
    vals, vec = np.linalg.eigh(gram)
    if min(vals) <= 0:
        raise PrefixFailure("NONPOSITIVE_SUBSPACE_MASS", eigenvalues=vals.tolist())
    transform = (vec/np.sqrt(vals))@vec.T
    coeff = (transform@coeff.reshape(mult, -1)).reshape(mult, 2, 6)
    result = []
    energy_model = SimpleNamespace(arms=(model, model), section=model.section, l=model.l)
    for c in coeff:
        arms = fields(model, case, w, c, order)
        mass, elastic, springs = energies(energy_model, case, arms)
        qscale = max(np.max(abs(a["y"][:, :3]*[1, 1, model.l])) for a in arms)
        fscale = max(np.max(abs(a["y"][:, 3:]*[1, 1, 1/model.l])) for a in arms)
        jr = physics.joint_matrix(model, case)@np.concatenate([a["ends"][1] for a in arms])
        er = np.array([physics.external_matrix(model, case, i)@a["ends"][0] for i, a in enumerate(arms)])
        js = abs(jr)/[qscale, qscale, fscale, fscale, qscale/model.l, fscale*model.l]
        es = abs(er)/[qscale, qscale, qscale/model.l]
        scale = np.array([1, 1, model.l, qscale/fscale, qscale/fscale, qscale/(fscale*model.l)])
        ode = max(np.max(abs((a["dy"]-a["y"]@physics.state_matrix(model, w).T)*scale))/
                  max(np.max(abs(a["dy"]*scale)), 1e-30) for a in arms)
        result.append(dict(coefficients=c.tolist(), arms=arms, diagnostics=dict(
            mass=float(mass.sum()), raw_mass_eigenvalues=vals.tolist(),
            energy_error=float(abs((elastic.sum()+sum(springs))/mass.sum()/w**2-1)),
            singular_ratio=float(s[-mult]/s[0]), nonzero_condition=float(s[0]/s[-mult-1]),
            joint_scaled=js.tolist(), external_scaled=es.tolist(), joint_absolute=jr.tolist(),
            external_absolute=er.tolist(), ode_error=float(ode))))
    return result


def mode_gates(model, modes, policy):
    gram = np.array([[physics.inner(model, a, b) for b in modes] for a in modes])
    gram_error = float(np.max(abs(gram-np.eye(len(modes)))))
    gates = []
    for mode in modes:
        d = mode["diagnostics"]
        gates.append(dict(boundary=max(np.max(d["joint_scaled"]), np.max(d["external_scaled"])) <= policy["boundary_tol"],
            ode=d["ode_error"] <= policy["ode_tol"], energy=d["energy_error"] <= policy["energy_tol"],
            singularity=d["singular_ratio"] <= policy["singular_tol"],
            complement=d["nonzero_condition"] <= policy["nonzero_condition_max"],
            mass=abs(d["mass"]-1) <= policy["mass_gram_tol"], mass_gram=gram_error <= policy["mass_gram_tol"]))
    return dict(status="PASS" if all(all(g.values()) for g in gates) else "MODE_QUALITY_FAILED",
                gram=gram.tolist(), gram_max_error=gram_error,
                gates=[{k: bool(v) for k, v in g.items()} for g in gates])
