"""Small algorithmic contracts; no physical root search or parameter sweep."""
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import rl_full_rr_prefix as p


@pytest.fixture
def policy():
    source = Path(__file__).parents[1]/"data/input/isotropic_two_springs_rl_timoshenko.json"
    return {**json.loads(source.read_text())["policy"], "nullity_tol": 1e-12,
            "max_cluster_extensions": 2, "max_cluster_extra_ranks": 4,
            "cluster_probe_relative": 1e-6}


def synthetic(roots, policy, guard=None, neighbour=None):
    roots = np.array(roots, float)
    matrix = lambda w: np.diag(np.r_[roots-w, 1.])
    counter = lambda w: (w, int(sum(roots < w)), {})
    return p.counted_prefix(matrix, counter, .1, float(max(roots)+1),
                            guard or len(roots), policy, neighbour)


def test_simple_count_certificate(policy):
    events, report = synthetic([1.], policy)
    assert events[0]["omega"] == pytest.approx(1.)
    assert report["lower_count"] == 0 and report["guard_right_count"] == 1


def test_double_without_sign_change(policy):
    events, report = synthetic([1., 1., 2.], policy)
    assert [(e["multiplicity"], e["first_rank"]) for e in events] == [(2, 1), (1, 3)]
    assert report["status"] == "PASS"
    assert events[0]["nullity"] == 2


def test_close_distinct_roots_not_merged(policy):
    events, _ = synthetic([1., 1.+1e-7, 2.], policy)
    assert len(events) == 3 and all(e["multiplicity"] == 1 for e in events)
    assert 0 < events[1]["omega"]-events[0]["omega"] < 2e-7


def test_subtolerance_distinct_projected_zeros_remain_unresolved(policy):
    with pytest.raises(p.PrefixFailure, match="UNRESOLVED_CLUSTER"):
        synthetic([1., 1.+1e-14, 2.], policy)


def test_whole_multiple_cluster_straddles_ranks_8_9(policy):
    events, report = synthetic([*range(1, 8), 8., 8., 8., 10.], policy, guard=9)
    assert events[-1]["first_rank"] == 8 and events[-1]["last_rank"] == 10
    assert events[-1]["multiplicity"] == 3 and report["computed_ranks"] == 10


def test_close_cluster_completion_is_count_triggered(policy):
    events, report = synthetic([*range(1, 9), 9., 9.+1e-7, 11.], policy, guard=9)
    assert report["computed_ranks"] == 10
    assert all(e["multiplicity"] == 1 for e in events)


def test_initial_guard_multiplicity_also_obeys_completion_budget(policy):
    with pytest.raises(p.PrefixFailure, match="UNRESOLVED_GUARD_CLUSTER"):
        synthetic([*range(1, 8), *([8.]*7), 10.], policy, guard=9)


def test_unresolved_cluster_cannot_pass(policy):
    counter = lambda w: (w, 2*int(w > 1.), {})
    # Inconsistent nullity: count claims two, matrix supplies only one mode.
    matrix = lambda w: np.diag([1.-w, 1., 2.])
    with pytest.raises(p.PrefixFailure, match="UNRESOLVED_CLUSTER|UNRESOLVED_NULLITY"):
        p.counted_prefix(matrix, counter, .1, 2., 2, policy)


def test_effective_count_outside_interval_stops(policy):
    counter = lambda w: (w if w in (.1, 2.) else 1.99, 0 if w < 1 else 2, {})
    with pytest.raises(p.PrefixFailure, match="COUNT_OUTSIDE_INTERVAL|COUNT_NOT_MONOTONE"):
        p.counted_prefix(lambda w: np.diag([1.-w, 1.-w, 1.]), counter, .1, 2., 2, policy)


def test_effective_count_outside_search_domain_stops(policy):
    with pytest.raises(p.PrefixFailure, match="COUNT_OUTSIDE_SEARCH_DOMAIN"):
        p.counted_prefix(lambda w: np.eye(2), lambda w: (3., 0, {}), .1, 2., 1, policy)


def test_catalog_bypass_checked_before_physics(policy, monkeypatch):
    m, c = p.physics.Model(), p.physics.Case(0.)
    cats = {L: {"upper": 1., "poles": [.99999995]} for L in c.lengths(m)}
    count = p.PhysicalCounter(m, c, cats, policy)
    def forbidden(*args):
        raise AssertionError("out-of-domain query reached the physical basis")
    monkeypatch.setattr(p.physics, "count_at", forbidden)
    with pytest.raises(p.PrefixFailure, match="COUNT_OUTSIDE_CATALOG"):
        count(.99999994)


@pytest.mark.parametrize("seed", [[2., 1.], [float('nan')], [9.]])
def test_bad_predictor_rejected_without_losing_roots(policy, seed):
    events, report = synthetic([1., 2.], policy, neighbour=seed)
    assert len(events) == 2 and report["predictor"].startswith("REJECTED")


def test_budgets_and_nonzero_lower_count(policy):
    with pytest.raises(p.PrefixFailure, match="COUNT_BUDGET_EXHAUSTED"):
        synthetic([1., 2.], {**policy, "max_count_queries": 1})
    with pytest.raises(p.PrefixFailure, match="NONZERO_LOWER_COUNT"):
        p.counted_prefix(lambda w: np.eye(2), lambda w: (w, 1, {}), .1, 2., 2, policy)
    with pytest.raises(p.PrefixFailure, match="INSUFFICIENT_RIGHT_COUNT"):
        synthetic([1.], policy, guard=2)


def test_full_contract_cache_identity(policy):
    m, c = p.physics.Model(), p.physics.Case(30.)
    base = p.point_contract(m, c, 8, 9, policy)
    records = [base, p.point_contract(m, c, 2, 3, policy),
               p.point_contract(replace(m, epsilon=.0025), c, 8, 9, policy),
               p.point_contract(m, replace(c, mu=.02), 8, 9, policy),
               p.point_contract(m, c, 8, 9, {**policy, "max_depth": 29})]
    assert len({p.digest(r) for r in records}) == len(records)
    assert p.digest(base) == p.digest(json.loads(json.dumps(base)))


def test_reject_other_physics(policy):
    with pytest.raises(ValueError):
        p.point_contract(p.physics.Model(variant="RL_PLANAR_LEGACY"), p.physics.Case(0), 8, 9, policy)


def test_read_only_table_and_cache_paths_never_solve(tmp_path, monkeypatch, capsys):
    from scripts.analysis import verify_rl_full_rr_prefix as run
    def forbidden(*args, **kwargs):
        raise AssertionError("read-only operation called a root solver")
    monkeypatch.setattr(run, "solve_point", forbidden)
    monkeypatch.setattr(run, "tim_catalog", forbidden)
    monkeypatch.setattr(p.physics, "ritz", forbidden)
    point = tmp_path/"point.json"
    point.write_text(json.dumps(dict(stage="B", beta_deg=0, status="PASS", roots=[])))
    manifest = dict(status="PASS", artifacts={"point.json": run.sha(point)}, point_files=["point.json"])
    (tmp_path/"manifest.json").write_text(json.dumps(manifest))
    assert run.main(["--tables-only", str(tmp_path)]) == 0
    assert "stage,beta_deg,rank" in capsys.readouterr().out
    assert run.main(["--verify-cache", str(tmp_path)]) == 0
    point.write_text("modified")
    with pytest.raises(p.PrefixFailure, match="ARTIFACT_HASH_MISMATCH"):
        run.verify_bundle(tmp_path)


def test_full_kernel_mass_normalization_and_rejection_of_nonphysical_modes(policy, monkeypatch):
    # Synthetic two-dimensional coefficient kernel; no physical eigenvalue claim.
    m, case = p.physics.Model(), p.physics.Case(0.)
    matrix = np.diag([0., 0., *([1.]*10)])
    monkeypatch.setattr(p.physics, "boundary", lambda *args: matrix)
    monkeypatch.setattr(p.physics, "balance", lambda raw: (raw, np.ones(12), np.ones(12)))
    modes = p.recovered_modes(m, case, dict(omega=.1, multiplicity=2), 120, policy)
    gram = np.array([[p.physics.inner(m, a, b) for b in modes] for a in modes])
    np.testing.assert_allclose(gram, np.eye(2), atol=2e-14)
    assert p.mode_gates(m, modes, policy)["status"] == "MODE_QUALITY_FAILED"


def test_scope_is_only_A_B_and_point_digest_detects_changes():
    from scripts.analysis import verify_rl_full_rr_prefix as run
    config = json.loads(run.CONFIG.read_text())
    run.validate_scope(config)
    config["tau"] = .1
    with pytest.raises(ValueError, match="restricted"):
        run.validate_scope(config)
    record = dict(status="PASS", omega=[1., 2.])
    record["record_digest"] = p.digest(record)
    run.verify_point_record(record)
    record["omega"][0] = .9
    with pytest.raises(p.PrefixFailure, match="POINT_RECORD_DIGEST_MISMATCH"):
        run.verify_point_record(record)


def test_failed_point_is_preserved_without_retry_or_false_pass(tmp_path, monkeypatch):
    from scripts.analysis import verify_rl_full_rr_prefix as run
    config = json.loads(run.CONFIG.read_text())
    def blocked(*args, **kwargs):
        raise p.PrefixFailure("CATALOG_CEILING_OUTSIDE_DOMAIN", ceiling=100., limit=90.)
    monkeypatch.setattr(run, "catalogues", blocked)
    result = run.solve_point(config, "A", 0., tmp_path, "synthetic-source")
    assert result["status"] == "CATALOG_CEILING_OUTSIDE_DOMAIN"
    assert result["failure"]["ceiling"] == 100.
    record = next((tmp_path/"points").glob("*.json"))
    before = record.read_bytes()
    again = run.solve_point(config, "A", 0., tmp_path, "synthetic-source")
    assert again == result and record.read_bytes() == before
