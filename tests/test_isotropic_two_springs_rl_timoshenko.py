"""Targeted energy/coordinate/cache regressions; no spectral map in tests."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.linalg import null_space
from scripts.lib import isotropic_two_springs_rl_timoshenko as p
from scripts.lib import coupled_longitudinal_comparators as legacy
from scripts.lib import contraction_joint_sensitivity as spring


def test_full_vs_planar_inertia_and_normalization():
    m = p.Model(); planar = replace(m, variant='RL_PLANAR_LEGACY')
    assert m.J == pytest.approx(17*planar.J)
    assert m.time == pytest.approx(200)
    assert m.Lambda(m.omega(3.9)) == pytest.approx(3.9)
    assert replace(m, variant='ELEMENTARY_LEGACY').J == 0


@pytest.mark.parametrize('beta', [0., 6.1194, 30., 90.])
@pytest.mark.parametrize('joint', [0., 1., 100., 'RIGID'])
def test_unchanged_planar_legacy_operator(beta, joint):
    m = p.Model(variant='RL_PLANAR_LEGACY'); c = p.Case(beta, k_joint=joint)
    old = spring.Case(beta, variant='rayleigh_love_planar',
                     rotation='RIGID' if joint == 'RIGID' else 'SPRING')
    k = 0 if joint == 'RIGID' else joint*m.section.EI/m.l
    w = float(m.omega(4.3))
    expected = spring.boundary(m.carrier, c.lengths(m), w, old, k, scaled=False)
    np.testing.assert_array_equal(p.boundary(m, c, w), expected)


@pytest.mark.parametrize('end', [0., 1., 'RIGID'])
@pytest.mark.parametrize('joint', [0., 1., 'RIGID'])
def test_trace_space_and_external_sign(end, joint):
    m = p.Model(); c = p.Case(30., k_joint=joint, k_end=end)
    trace, energy = p.trace_map(m, c)
    assert trace.shape[1] == 2+(1 if joint == 'RIGID' else 2)+(end != 'RIGID')
    assert np.linalg.matrix_rank(trace) == trace.shape[1]
    assert np.linalg.eigvalsh(energy).min() >= -1e-20
    op = p.external_matrix(m, c, 1)
    if end == 'RIGID':
        np.testing.assert_array_equal(op, np.eye(6)[:3])
    else:
        assert op[2, 5] == -1
        assert op[2, 2] == end*m.section.EI/m.l
        # Its nullspace permits nonzero rotation, including the outer hinge.
        assert np.linalg.norm(null_space(op)[2]) > .1


@pytest.mark.parametrize('variant', p.VARIANTS)
def test_gradient_force_and_independent_balance(variant):
    m = p.Model(variant=variant); w = float(m.omega(4.))
    b = p.basis(m, 1.01, w, np.linspace(0, 1.01, 12))
    d = p.basis(m, 1.01, w, np.linspace(0, 1.01, 12), 1)
    np.testing.assert_allclose(b[:, 3], (m.section.EA-m.J*w*w)*d[:, 0], atol=1e-15)
    np.testing.assert_allclose(d, np.einsum('ij,njk->nik', p.state_matrix(m, w), b), atol=2e-12)


def test_scaling_pullback_and_hinge_force_rows():
    m = p.Model(); c = p.Case(30., k_joint=0., k_end=0.)
    raw = p.boundary(m, c, float(m.omega(4.)))
    b, row, col = p.balance(raw)
    v = np.random.default_rng(3).normal(size=12)
    np.testing.assert_allclose(b@v, row*(raw@(col*v)), atol=1e-14)
    joint = p.joint_matrix(m, c)
    assert np.linalg.matrix_rank(joint) == 6
    assert np.linalg.norm(joint[2:4]) > 1  # forces survive the hinge limit
    np.testing.assert_array_equal(joint[4], np.eye(12)[5])
    np.testing.assert_array_equal(joint[5], np.eye(12)[11])


def test_plot_only_does_not_call_physics(tmp_path, monkeypatch):
    from scripts.analysis import verify_isotropic_two_springs_rl_timoshenko as cli
    def forbidden(*args, **kwargs):
        raise AssertionError('plot-only called a physical solver')
    for name in ('prefix', 'boundary', 'mode', 'ritz'):
        monkeypatch.setattr(p, name, forbidden)
    (tmp_path/'spectrum.csv').write_text('beta_deg,Lambda1,Lambda2\n0,2.36,3.92\n10,3.92,4.38\n90,3.92,4.72\n')
    cli.plot(tmp_path)
    assert len(list(tmp_path.glob('*.png'))) == 2
    assert len(list(tmp_path.glob('*.pdf'))) == 2


def test_tail_range_uses_physical_frequency_and_interior_extrema():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import tail_ranges
    rows = [{'ID': 'test', 'beta_deg': b, 'omega1': w, 'Lambda1': np.sqrt(w)}
            for b, w in [(10, 4), (20, 9), (40, 6), (90, 5)]]
    result = tail_ranges(rows)
    assert result[0]['R'] == pytest.approx(1.25)  # not sqrt(9)/sqrt(4)-1
    assert result[0]['beta_max'] == 20
    assert result[1]['beta_min'] == 90
    assert result[-1]['R'] == 0
    assert result[-1]['nodes'] == 1  # singleton is not evidence of a plateau


def test_limit_plot_only_no_physical_calls(tmp_path, monkeypatch):
    from scripts.analysis import verify_isotropic_two_springs_rl_timoshenko as cli
    def forbidden(*args, **kwargs):
        raise AssertionError('limit plot-only invoked physics')
    for name in ('prefix', 'count_at', 'boundary', 'basis', 'mode', 'ritz'):
        monkeypatch.setattr(p, name, forbidden)
    rows = [{'ID': label, 'beta_deg': beta, 'Lambda1': 2+beta/100}
            for label in ('RR', 'HR', 'RH', 'HH') for beta in (0, 10, 90)]
    cli.csvwrite(tmp_path/'spectrum.csv', rows)
    cli.csvwrite(tmp_path/'tail_ranges.csv', [{**r, 'R': (90-r['beta_deg'])/100} for r in rows])
    cli.plot_limits(tmp_path)
    assert len(list(tmp_path.glob('*.png'))) == 3
    assert len(list(tmp_path.glob('*.pdf'))) == 3


def test_saturation_log_derivative_units_and_frequency_scale():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import saturation_measures
    beta = np.arange(91)
    def table(scale):
        return [{'beta_deg': int(b), 'omega1': scale*np.exp(.02*b),
                 'Omega1': 200*scale*np.exp(.02*b), 'Lambda1': np.sqrt(200*scale*np.exp(.02*b))} for b in beta]
    first = saturation_measures('test', table(1.))
    scaled = saturation_measures('test', table(7.))
    for a, b in zip(first, scaled):
        assert a['F'] == pytest.approx(b['F'], abs=1e-14)
        assert a['R'] == pytest.approx(b['R'], abs=1e-14)
        assert a['S_h1_per_degree'] == pytest.approx(.02, abs=1e-14)
        for h in (2, 3):
            if a[f'S_h{h}_per_degree'] is not None:
                assert a[f'S_h{h}_per_degree'] == pytest.approx(.02, abs=1e-14)


def test_saturation_preserves_internal_maximum_and_tail_extrema():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import saturation_measures
    table = [{'beta_deg': b, 'omega1': 2+b-.01*b*b} for b in range(91)]
    rows = saturation_measures('test', table)
    assert rows[50]['F'] == 1
    assert rows[90]['F'] < 1
    assert rows[20]['tail_beta_max'] == 50
    assert rows[20]['tail_beta_min'] == 90
    assert rows[20]['R'] == pytest.approx((27-11)/11)
    assert rows[60]['S_h1_per_degree'] < 0
    assert rows[90]['R'] == 0


def test_saturation_methods_recover_known_knot_without_window_renormalization():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import saturation_measures, saturation_transitions
    table = []
    for b in range(91):
        w = 2+.2*min(b, 8)+.001*max(0, b-8)
        table.append(dict(beta_deg=b, omega1=w, Omega1=200*w, Lambda1=np.sqrt(200*w)))
    rows = saturation_measures('test', table)
    for window in (40, 60, 90):
        transitions, objectives, curves = saturation_transitions(rows, window)
        assert [t['beta_deg'] for t in transitions] == [8, 8]
        assert min(r['SSE_F'] for r in objectives) < 1e-25
        assert curves[-1]['F'] == rows[window]['F']
        assert curves[8]['segmented_F'] == pytest.approx(rows[8]['F'])


def test_plateau_length_restricts_existence_without_clamping_angle():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import plateau_tails, plateau_at
    table = [dict(beta_deg=b, omega1=100+b, Omega1=200*(100+b),
                  Lambda1=np.sqrt(200*(100+b))) for b in range(91)]
    rows = plateau_tails('linear', table)
    # (190-(100+b))/(100+b) <= 0.2 first holds at integer b=59.
    for length in (20, 30):
        r = plateau_at(rows, 20, length)
        assert r['beta_pl_deg'] == 59
        assert r['previous_R_percent'] > 20 >= r['R_pl_percent']
        assert r['q_pl'] == pytest.approx(159/190)
    r = plateau_at(rows, 20, 45)
    assert r['status'] == 'NO_QUALIFYING_PLATEAU'
    assert r['beta_pl_deg'] is None and r['Lambda_pl'] is None
    assert r['first_tolerance_hit_deg'] == 59  # diagnostic only, not an accepted angle
    # The zero singleton at 90 degrees cannot pass a positive-length requirement.
    assert plateau_at(rows, 1e-12, 20)['beta_pl_deg'] is None


def test_plateau_exact_event_equality_and_nonmonotone_spectrum():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import plateau_tails, plateau_at
    table = [dict(beta_deg=b, omega1=2+b-.01*b*b, Omega1=200*(2+b-.01*b*b),
                  Lambda1=np.sqrt(200*(2+b-.01*b*b))) for b in range(91)]
    rows = plateau_tails('peak', table)
    assert rows[20]['beta_max'] == 50 and rows[20]['beta_min'] == 90
    assert all(b['R'] <= a['R'] for a,b in zip(rows,rows[1:]))
    event = rows[60]['R_percent']
    assert plateau_at(rows, event, 30)['beta_pl_deg'] == 60
    assert plateau_at(rows, np.nextafter(event, -np.inf), 30)['status'] == 'NO_QUALIFYING_PLATEAU'
    flat = plateau_tails('flat', [dict(beta_deg=b, omega1=2, Omega1=400, Lambda1=20) for b in range(91)])
    assert plateau_at(flat, .01, 45)['beta_pl_deg'] == 0
    assert plateau_at(flat, .01, 45)['previous_R_percent'] is None


def test_plateau_pair_intervals_do_not_rank_missing_values():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import plateau_tails, plateau_at, plateau_order_intervals
    def table(label, slope):
        return plateau_tails(label, [dict(beta_deg=b, omega1=100+slope*b,
             Omega1=200*(100+slope*b), Lambda1=np.sqrt(200*(100+slope*b))) for b in range(91)])
    tables = {'A':table('A',1), 'B':table('B',2)}
    intervals = plateau_order_intervals(tables,[['A','B']],1,80,30)
    assert intervals[0]['relation_A_to_B'] == 'BOTH_UNAVAILABLE'
    assert any(r['relation_A_to_B'] == 'B_UNAVAILABLE' for r in intervals)
    for r in intervals:
        d=(r['delta_lower_percent']+r['delta_upper_percent'])/2
        a=plateau_at(tables['A'],d,30)['beta_pl_deg']; b=plateau_at(tables['B'],d,30)['beta_pl_deg']
        if r['relation_A_to_B']=='EARLIER':assert a < b
        if r['relation_A_to_B']=='TIED':assert a == b and a is not None


def test_grid_minmax_checks_both_springs_and_guard_position():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import grid_order
    grid = [.1, 1., 10.]
    cases = {f'{a}{j+1}': [x, y] for a, x in zip('ABC', grid) for j, y in enumerate(grid)}
    cfg = dict(cases=cases, stiffness_grid=grid, map={'parameter_grid': [0]},
               comparison_relative_tolerance=1e-8)
    def row(label, w):
        return dict(ID=label, beta_deg=0, omega1=w, omega2=2*w, omega3=3*w)
    rows = [row(label, 1+x+y) for label, (x, y) in cases.items()]
    bounds = [row('HH', 1), row('RR', 100)]
    for x in grid:
        bounds += [row(f'J{x:g}_ER', 50+x), row(f'JR_E{x:g}', 50+x)]
    assert len(grid_order(rows, bounds, cfg)) == 144
    # A guard-root violation must fail even if the two plotted roots pass.
    rows[0]['omega3'] = 1000
    with pytest.raises(ArithmeticError, match='min-max'):
        grid_order(rows, bounds, cfg)


def test_early_interval_log_change_is_scale_invariant_and_not_linear_percent():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import early_measures
    def table(scale):
        return [dict(beta_deg=b, omega1=scale*np.exp(.03*b),
            Omega1=200*scale*np.exp(.03*b), Lambda1=np.sqrt(200*scale*np.exp(.03*b))) for b in range(91)]
    intervals = [[0,5],[5,10],[10,15],[15,20],[20,30]]
    a, av = early_measures('test',table(1),intervals)
    b, bv = early_measures('test',table(7),intervals)
    for x,y in zip(av,bv):
        assert x['mean_log_percent_per_degree'] == pytest.approx(3, abs=1e-13)
        assert y['mean_log_percent_per_degree'] == pytest.approx(3, abs=1e-13)
        assert x['reduction_factor'] == pytest.approx(1)
    assert av[0]['signed_interval_change_percent'] == pytest.approx(100*np.expm1(.15))
    assert a[15]['G90_percent'] == pytest.approx(b[15]['G90_percent'])


def test_early_negative_growth_and_internal_peak_are_preserved():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import early_measures
    table=[]
    for b in range(91):
        w=2+np.exp(-((b-8)/8)**2)
        table.append(dict(beta_deg=b,omega1=w,Omega1=200*w,Lambda1=np.sqrt(200*w)))
    rows, intervals=early_measures('peak',table,[[0,5],[5,10],[10,15],[15,20],[20,30]])
    assert intervals[2]['mean_log_percent_per_degree'] < 0
    assert intervals[2]['reduction_factor'] > 0
    assert rows[0]['beta_max']==8 and rows[15]['G90_percent']<0
    assert rows[15]['R_percent']>abs(rows[15]['G90_percent'])
    flat=[dict(beta_deg=b,omega1=2,Omega1=400,Lambda1=20) for b in range(91)]
    rows, intervals=early_measures('flat',flat,[[0,5],[5,10]])
    assert rows[15]['R']==rows[15]['G90_percent']==0
    assert intervals[0]['reduction_factor'] is None


def test_growth_transition_log_differences_and_scale_invariance():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import growth_slopes
    beta=np.arange(91,dtype=float);y=.2*beta-.001*beta**2
    for h in (1,2):
        slope=growth_slopes(y,h)
        assert np.max(abs(slope-(.2-.002*beta)))<1e-14
        assert np.max(abs(growth_slopes(y+7,h)-slope))<1e-14
        s0=(y[5]-y[0])/5
        assert np.allclose(abs(slope)/s0,abs(growth_slopes(y/2+3,h))/(s0/2),rtol=0,atol=1e-13)
    # A log-frequency maximum changes the sign; taking abs cannot hide its tail.
    assert growth_slopes(-.01*(beta-40)**2,1)[50]<0


def test_growth_transition_tail_rejects_early_local_crossing():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import growth_pick
    q=np.full(91,.02);q[:5]=1;q[20]=.3
    r,envelope=growth_pick(q,.1)
    assert r['beta_pl_deg']==21 and r['previous_Q']==.3
    assert all(b<=a for a,b in zip(envelope,envelope[1:]))
    assert r['tail_peak_beta']==21
    # A late peak outside the admissible start interval makes the result missing.
    q[80]=.11
    assert growth_pick(q,.1)[0]['status']=='NO_QUALIFYING_TRANSITION'
    assert growth_pick(q,.1)[0]['beta_pl_deg'] is None


def test_growth_transition_equality_missing_and_length_are_not_clamped():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import growth_pick
    q=np.full(91,.1);q[:60]=.2
    assert growth_pick(q,.1)[0]['beta_pl_deg']==60
    assert growth_pick(q,np.nextafter(.1,0))[0]['beta_pl_deg'] is None
    assert growth_pick(q,.1,31)[0]['beta_pl_deg'] is None
    assert growth_pick(q,.2)[0]['beta_pl_deg']==0
    with pytest.raises(ValueError):growth_pick(-q,.2)


def test_refined_explicit_stencils_preserve_archived_edges_and_polynomial_slope():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import refined_slopes, growth_slopes
    beta=np.array(sorted(set(range(91)) | set(np.arange(3,12.01,.25))))
    y=.2*beta-.001*beta**2;oldbeta=np.arange(91)
    archived=growth_slopes(.2*oldbeta-.001*oldbeta**2,1)
    for h in (.25,.5):
        slope,stencil,valid=refined_slopes(beta,y,archived,h)
        assert max(abs(slope-(.2-.002*beta)))<1e-14
        assert stencil[list(beta).index(3)]=='archived_h1'
        assert stencil[list(beta).index(12)]=='archived_h1'
        if h==.5:assert stencil[list(beta).index(3.25)]=='primary_edge_fallback_0.25'
        assert sum(valid)==33
    with pytest.raises(ValueError):refined_slopes(beta[:-1],y[:-1],archived,.25)


def test_refined_persistent_crossing_uses_actual_irregular_angles_and_no_slack():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import refined_pick
    beta=np.array(sorted(set(range(91)) | set(np.arange(3,12.01,.25))))
    q=np.full(len(beta),.1);q[beta<6.25]=.2
    r,_=refined_pick(beta,q,.1)
    assert r['beta_pl_deg']==6.25 and r['previous_beta_deg']==6
    assert r['previous_q']==.2
    assert refined_pick(beta,q,np.nextafter(.1,0))[0]['beta_pl_deg'] is None
    q[beta==70]=.3
    assert refined_pick(beta,q,.1)[0]['status']=='NO_QUALIFYING_TRANSITION'
