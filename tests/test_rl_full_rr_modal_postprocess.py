"""Saved-field checks only; existing solver tests remain unchanged."""
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.analysis import postprocess_rl_full_rr_modes as audit


def test_global_coordinate_metric_does_not_compare_local_components_as_global():
    m = audit.p.Model()
    def mode(right_component):
        arms = []
        for arm in range(2):
            y = np.zeros((2, 6)); dy = y.copy()
            if arm == 1:
                y[:, right_component] = 1.
            arms.append(dict(y=y, dy=dy, weights=np.array([.5, .5])))
        return dict(arms=arms)
    before, after = mode(1), mode(0)
    assert audit.p.inner(m, before, after) == 0
    cross = audit.p.inner(m, before, after, 0., 90.)
    norm = np.sqrt(audit.p.inner(m, before, before)*audit.p.inner(m, after, after))
    assert abs(cross/norm) == pytest.approx(1., abs=1e-14)


def test_component_nodes_are_brackets_and_noise_is_not_a_wave():
    x = np.array([0., .2, .4, .6, .8, 1.])
    r = audit.node_brackets(x, np.array([0., 1., 1., -1., -1., 0.]), 1.)
    assert r['count'] == 1 and r['brackets'] == [[.4, .6]]
    r = audit.node_brackets(x, 1e-13*np.array([0., 1., -1., 1., -1., 0.]), 1.)
    assert r['count'] is None and r['status'] == 'BELOW_AMPLITUDE_FLOOR'


def test_frequency_changes_use_omega_not_Lambda():
    rows = [dict(beta_deg=b, rank=k, omega=2.+b, Lambda=np.sqrt(2.+b))
            for b in audit.ANGLES for k in range(1, 9)]
    first = audit.frequency_intervals(rows)[0]
    assert first['omega_ratio'] == 3.5
    assert first['increment_percent'] == 250.


def test_full_potential_and_mass_components_against_polynomial_integrals():
    m = audit.p.Model(); case = audit.p.Case(0.)
    z, wt = np.polynomial.legendre.leggauss(8)
    arms = []
    for L in case.lengths(m):
        x = (z+1)*L/2; y = np.zeros((8, 6)); dy = np.zeros_like(y)
        y[:, 0], y[:, 1], y[:, 2] = x, 2*x, .3
        dy[:, 0], dy[:, 1] = 1., 2.
        arms.append(dict(y=y, dy=dy, weights=wt*L/2, ends=np.zeros((2, 6))))
    adapter = SimpleNamespace(arms=(m,m), section=m.section, l=m.l)
    mass, stiffness, springs = audit.section_physics.energies(adapter, case, arms)
    for i, L in enumerate(case.lengths(m)):
        np.testing.assert_allclose(mass[i], [m.section.rhoA*L**3/3, 4*m.section.rhoA*L**3/3,
            m.section.rhoI*.09*L, m.J*L], rtol=2e-14)
        np.testing.assert_allclose(stiffness[i], [m.section.EA*L, 0., m.section.KGA*1.7**2*L], rtol=2e-14)
    np.testing.assert_array_equal(springs, [0., 0.])


@pytest.mark.skipif(not (audit.SOURCE/'manifest.json').is_file(),
                   reason='Accepted local K8 archive is not included in Git')
def test_saved_archive_closes_under_solver_and_basis_prohibition():
    with audit.no_solves():
        with pytest.raises(AssertionError, match='cannot'):
            audit.p.ritz(None, None, 20)
        _, model, records, modes = audit.load_source(audit.SOURCE)
        rows, scales, nodes, checks, parity = audit.analyze(model, records, modes)
        assert len(rows) == 48 and len(checks) == 54
        assert max(c['energy_error'] for c in checks) < 5e-8
        assert [p['indicated_parity'] for p in parity] == ['symmetric','antisymmetric']*4
        assert all(r['axial_chord_share'] is None for r in rows if r['beta_deg'] == 0)


def test_cannot_write_into_source_bundle():
    with pytest.raises(ValueError, match='outside'):
        audit.run(audit.SOURCE, audit.SOURCE/'derived')
