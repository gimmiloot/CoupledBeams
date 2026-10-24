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


def test_stiffness_composite_stencils_include_nonuniform_join_and_boundaries():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import stiffness_slopes
    beta=np.array(sorted(set(range(21)) | {3+i*.25 for i in range(37)} | set(range(22,91,2))))
    y=.13*beta-.0003*beta**2
    d,stencil=stiffness_slopes(beta,y)
    assert len(beta)==83 and max(abs(d-(.13-.0006*beta)))<2e-15
    assert stencil[list(beta).index(20)]=='nonuniform_19_20_22'
    assert max(abs(stiffness_slopes(beta,y+4)[0]-d))<1e-14
    with pytest.raises(ValueError):stiffness_slopes(beta[:-1],y[:-1])


def test_stiffness_finite_effects_resolve_additivity_and_keep_signed_increments():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import stiffness_effects
    def inputs(interaction):
        t=[];f=[]
        for j in range(5):
            for e in range(5):
                label=f'J{j}E{e}'
                t.append(dict(ID=label,q_star=.15,beta_pl_deg=5+2*j+3*e+interaction*j*e))
                f.append(dict(ID=label,Lambda10=2+j+e,Lambda90=3+j+e))
        return t,f
    inc,mixed=stiffness_effects(*inputs(0))
    assert len(inc)==40 and len(mixed)==16
    assert all(r['delta_beta_deg']==(2 if r['axis']=='joint' else 3) for r in inc)
    assert all(r['mixed_difference_deg']==0 for r in mixed)
    inc,mixed=stiffness_effects(*inputs(-2))
    assert all(r['mixed_difference_deg']==-2 for r in mixed)
    assert any(r['delta_beta_deg']<0 for r in inc)


def test_mechanics_energy_matches_polynomial_integrals_and_keeps_RL_in_mass():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import mechanics_energy
    m=p.Model();case=p.Case(6,k_joint=.1,k_end=10);s=m.section
    nodes,weights=np.polynomial.legendre.leggauss(8);arms=[];expected=np.zeros(5);mass=0.;ends=[]
    for L,(a,b,c,d) in zip(case.lengths(m),[(.3,.7,.2,0.),(-.4,.5,.8,.1)]):
        x=(nodes+1)*L/2;y=np.zeros((len(x),6));dy=y.copy()
        y[:,:3]=np.column_stack((a*x,b*x*x,c*x+d));dy[:,:3]=np.column_stack((a+0*x,2*b*x,c+0*x))
        end=np.zeros((2,6));end[:,:3]=[[0,0,d],[a*L,b*L*L,c*L+d]];ends.append(end)
        arms.append(dict(y=y,dy=dy,weights=weights*L/2,ends=end))
        expected[:3]+=[s.EA*a*a*L/2,s.EI*c*c*L/2,s.KGA*((2*b-c)**2*L**3/3-(2*b-c)*d*L**2+d*d*L)/2]
        mass+=s.rhoA*(a*a*L**3/3+b*b*L**5/5)+s.rhoI*(c*c*L**3/3+c*d*L*L+d*d*L)+m.J*a*a*L
    expected[3]=case.k_joint*s.EI/m.l*(ends[0][1,2]-ends[1][1,2])**2/2
    expected[4]=case.k_end*s.EI/m.l*ends[1][0,2]**2/2
    result=mechanics_energy(m,case,arms)
    np.testing.assert_allclose(result['energies'],expected,rtol=2e-14,atol=1e-20)
    assert result['mass']==pytest.approx(mass,rel=2e-14)
    planar=mechanics_energy(replace(m,variant='RL_PLANAR_LEGACY'),case,arms)
    np.testing.assert_array_equal(planar['energies'],result['energies'])
    assert planar['mass']<result['mass']
    for arm in arms:
        for name in ('y','dy','ends'):arm[name]*=-3
    scaled=mechanics_energy(m,case,arms)
    np.testing.assert_allclose(scaled['fractions'],result['fractions'],rtol=2e-14)
    assert scaled['mass']==pytest.approx(9*result['mass'])


def test_mechanics_global_coordinates_use_verified_inward_arm_frames():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import mechanics_global
    m=p.Model();coef=np.arange(12).reshape(2,6)/13
    for beta in (0,20):
        case=p.Case(beta,k_joint=1,k_end=1);g=mechanics_global(m,case,coef,.05,21)
        angle=np.deg2rad(beta);L1,L2=case.lengths(m)
        np.testing.assert_allclose(g[0]['reference'][0],[-L1,0],atol=1e-15)
        np.testing.assert_allclose(g[1]['reference'][0],[L2*np.cos(angle),L2*np.sin(angle)],atol=1e-15)
        for a in g:
            np.testing.assert_array_equal(a['reference'][-1],[0,0])
            np.testing.assert_allclose(a['displacement'],a['axial']+a['transverse'],atol=1e-14)
            np.testing.assert_allclose(np.sum(a['displacement']**2,axis=1),np.sum(a['y'][:,:2]**2,axis=1),rtol=1e-14)


def test_variation_constraint_derivative_matches_physical_operator_and_reflection():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import variation_constraints
    for beta in (0,1,6,12):
        b=np.deg2rad(beta);C,Cb=variation_constraints(b)
        actual=p.joint_matrix(p.Model(),p.Case(beta,k_joint=.1,k_end=10))[:2][:,[0,1,6,7]]
        np.testing.assert_allclose(C,np.diag([1,-1])@actual,atol=1e-15)
        h=1e-5
        np.testing.assert_allclose(Cb,(variation_constraints(b+h)[0]-variation_constraints(b-h)[0])/(2*h),atol=3e-11)
        np.testing.assert_allclose(variation_constraints(-b)[0]@np.diag([-1,1,-1,1]),np.diag([-1,1])@C,atol=1e-15)


def test_variation_multiplier_sign_and_mass_factor_on_exact_constrained_quadratic():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import variation_work
    # Exact one-dimensional constrained quadratic fixture, no eigenvalue solver.
    K=np.diag([1.,4.]);M=np.diag([2.,1.])
    for b in (.1,.4,.9):
        C=np.array([np.cos(b),np.sin(b)]);Cb=np.array([-np.sin(b),np.cos(b)])
        z=Cb;zp=-C;mass=z@M@z;stiff=z@K@z;lam=stiff/mass
        r=-C@(K-lam*M)@z
        direct=(2*zp@K@z*mass-stiff*2*zp@M@z)/mass**2/(2*lam)
        for scale in (1.,-3.,.2):
            result=variation_work(np.array([scale*r]),np.array([Cb@(scale*z)]),np.sqrt(lam),scale**2*mass)
            assert result['s_rad']==pytest.approx(direct,rel=2e-14)


def test_variation_fd_uses_available_stencils_without_interpolated_points():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import variation_fd
    table={float(b):{'omega1':np.exp(.01*b*b)} for b in list(range(21))+[3.25,3.5,11.5,11.75]}
    assert variation_fd(table,1,.25) is None
    for b,h,scheme in [(0,1,'forward'),(3,.25,'forward'),(12,.25,'backward'),(6,1,'central')]:
        r=variation_fd(table,b,h)
        assert r['scheme']==scheme
        assert r['s_FD_deg']==pytest.approx(.02*b,abs=1e-14)


def test_epsilon_family_preserves_dimensionless_springs_and_full_RL_geometry():
    base=p.Model().section
    for epsilon in (.0025,.005,.01):
        m=p.Model(epsilon=epsilon);s=m.section;ratio=epsilon/.005
        assert s.area/base.area==pytest.approx(ratio**2)
        assert s.EI/base.EI==pytest.approx(ratio**4)
        assert s.EA*m.l**2/s.EI==pytest.approx(epsilon**-2)
        assert m.J==pytest.approx(m.rho*m.nu**2*17*s.second_moment)
        for kappa in (.1,10):
            dimensional=kappa*s.EI/m.l
            assert dimensional*m.l/s.EI==pytest.approx(kappa)
            assert dimensional/(kappa*base.EI)==pytest.approx(ratio**4)


def test_epsilon_curvature_recovers_even_radian_coefficients_without_intercept():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import epsilon_curvature
    rows=[]
    for b in np.arange(0,1.01,.125):
        x=np.deg2rad(b)
        rows.append(dict(ID='fixture',epsilon=.005,beta_deg=b,omega1=3*np.exp(71*x*x-1200*x**4)))
    fits=epsilon_curvature(rows,[.5,1],[2,4,6])
    for r in fits:
        assert r['positive_points']>r['degree']//2
        if r['degree']>=4:
            assert r['a2_rad']==pytest.approx(71,rel=2e-10)
            assert r['a2_epsilon_squared']==pytest.approx(71*.005**2,rel=2e-10)
            assert r['max_log_residual']<1e-14


def test_slip_branch_work_rank_and_exact_rigid_limit():
    from scripts.lib import isotropic_two_springs_longitudinal_slip as slip
    from scipy.linalg import null_space
    m=p.Model();rng=np.random.default_rng(119)
    for b in (0,3,20):
        for chi in (0,.01,1,'RIGID'):
            c=slip.Case(b,k_joint=.1,k_end=10,chi=chi);a=slip.joint_matrix(m,c)
            P,K=slip.trace_map(m,c);u=null_space(a/np.linalg.norm(a,axis=1)[:,None])@rng.normal(size=6)
            states=u.reshape(2,6);v=rng.normal(size=P.shape[1]);dq=P@v
            q=np.r_[np.zeros(3),states[0,:3],np.zeros(3),states[1,:3]]
            coordinates=np.linalg.lstsq(P,q,rcond=None)[0]
            work=states[0,3:]@dq[3:6]+states[1,3:]@dq[9:12]+coordinates@K@v
            assert abs(work)<1e-12
            assert np.linalg.matrix_rank(a)==6
            assert np.linalg.matrix_rank(P)==(5 if chi=='RIGID' else 6)
            np.testing.assert_allclose(K,K.T,atol=1e-20)
            if chi=='RIGID':
                np.testing.assert_array_equal(a,p.joint_matrix(m,c))
                np.testing.assert_array_equal(P,p.trace_map(m,c)[0])
                np.testing.assert_array_equal(slip.boundary(m,c,.04),p.boundary(m,c,.04))


def test_slip_variation_requires_the_remaining_constraint_term():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import variation_constraints
    for b in (.1,.5):
        C,Cb=variation_constraints(b);z=np.array([.3,0.,-.4,.7]);z[1]=-C[1,2:]@z[2:]
        c1=C[0]@z;k=.23;N=-k*c1;Q=.17
        new=k*c1*(Cb[0]@z)-Q*(Cb[1]@z)
        assert new==pytest.approx(N*z[1]-Q*z[0]+Q*c1)
        h=1e-6
        dE=k*((variation_constraints(b+h)[0][0]@z)**2-(variation_constraints(b-h)[0][0]@z)**2)/(4*h)
        assert dE==pytest.approx(k*c1*(Cb[0]@z),rel=1e-9)


def test_thickness_derivative_handles_step_change_and_reflection_symmetry():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import thickness_slopes
    b=np.r_[np.arange(0,24.1,.5),np.arange(26,91,2)]
    # An even quadratic has exactly known slope, including the 24-degree join.
    a,w=thickness_slopes(b,2+3e-4*b*b)
    np.testing.assert_allclose(a,6e-4*b,atol=3e-15)
    np.testing.assert_allclose(w,6e-4*b,atol=3e-15)
    # Symmetry is imposed explicitly; this is not the forward secant at zero.
    a,_=thickness_slopes(b,np.log1p(.2*b*b))
    assert a[0]==0


def test_thickness_scaled_reference_preserves_similar_crossing():
    from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import thickness_slopes,refined_pick
    b=np.r_[np.arange(0,24.1,.5),np.arange(26,91,2)]
    y=np.log1p(3*b*b/(4+b*b))
    s,_=thickness_slopes(b,y);s0=(y[list(b).index(5)]-y[0])/5
    old,_=refined_pick(b,abs(s)/s0,.15)
    # Pure scaling of the actual nodes and initial interval is exact for the stencil.
    doubled,_=thickness_slopes(2*b,y)
    new,_=refined_pick(2*b,abs(doubled)/(s0/2),.15)
    assert new['beta_pl_deg']==2*old['beta_pl_deg']


def test_angle_asymptotics_dimensional_energy_and_mass_scales():
    # Check the proposed scales against the unchanged dimensional model,
    # including a nonunit l and amplitude (no eigenproblem is assembled).
    for eps in (.0025, .005, .01):
        m = p.Model(epsilon=eps, l=2.3, E=7., rho=3.)
        s = m.section
        amplitude = .037
        energy = s.EI * amplitude**2 / m.l**3
        mass = s.rhoA * amplitude**2 * m.l
        assert s.EA*(amplitude*eps/m.l)**2*m.l/energy == pytest.approx(1)
        assert s.EI*(amplitude/m.l**2)**2*m.l/energy == pytest.approx(1)
        assert s.KGA*(amplitude/m.l)**2*m.l/energy == pytest.approx(
            m.ks/(2*(1+m.nu))*eps**-2)
        assert s.rhoA*(amplitude*eps)**2*m.l/mass == pytest.approx(eps**2)
        assert s.rhoI*(amplitude/m.l)**2*m.l/mass == pytest.approx(eps**2)
        assert m.J*(amplitude*eps/m.l)**2*m.l/mass == pytest.approx(
            17*m.nu**2*eps**4)
        for kappa in (.1, 10):
            assert (kappa*s.EI/m.l)*(amplitude/m.l)**2/energy == pytest.approx(kappa)
        assert m.omega(3.4)**2*mass/energy == pytest.approx(3.4**4)


def test_angle_asymptotics_constraints_signs_and_remainder():
    # The canonical matrix uses global Y; the report uses n1=-Y.
    U1, U2, W1, W2, xi = .3, -.7, .4, 1.2, 3.
    errors = []
    for eps in (.04, .02, .01):
        beta = eps*xi
        matrix = p.joint_matrix(p.Model(), p.Case(np.rad2deg(beta), k_joint=.1, k_end=10))
        C = np.diag([1., -1.]) @ matrix[:2, [0, 1, 6, 7]]
        c = C @ np.array([eps*U1, W1, eps*U2, W2])
        reduced = np.array([c[0]/eps, c[1]])
        expansion = np.array([
            U1+U2+xi*W2-eps**2*(xi**2*U2/2+xi**3*W2/6),
            W1+W2-eps**2*(xi*U2+xi**2*W2/2)])
        fourth = eps**4*np.array([xi**4*U2/24+xi**5*W2/120,
                                   xi**3*U2/6+xi**4*W2/24])
        sixth_bound = eps**6*np.array([abs(U2)*xi**6/720+abs(W2)*xi**7/5040,
                                       abs(U2)*xi**5/120+abs(W2)*xi**6/720])
        # Explicit Taylor bound avoids mistaking the finite sixth-order term
        # for a failed fourth-order limit (the first coefficient nearly cancels).
        assert np.all(abs(reduced-expansion-fourth) <= sixth_bound+2e-15)
        errors.append(abs(reduced-expansion))
    assert np.all(np.array(errors[:-1]) > errors[1:])
    m = p.Model(); c = p.Case(0, k_joint=.1, k_end=10)
    np.testing.assert_array_equal(p.external_matrix(m,c,0), np.eye(6)[:3])
    np.testing.assert_array_equal(p.external_matrix(m,c,1)[:2], np.eye(6)[:2])


def test_angle_asymptotics_axial_minimum_exact_rational_identity():
    from fractions import Fraction as F
    # U_i = p_i*s + h_i*s*(L_i-s). Endpoint redistribution and interior
    # perturbations are independent; their energy excess must be positive.
    for mu in (F(0), F(1,100), F(3,10)):
        L1, L2 = 1-mu, 1+mu
        xi, y = F(7,3), F(-2,5)
        d = -xi*y
        slope = d/(L1+L2)
        for delta in (F(0), F(2,7), F(-1,9)):
            for h1, h2 in ((F(0),F(0)), (F(3,7),F(-5,11))):
                p1, p2 = slope+delta/L1, slope-delta/L2
                assert p1*L1+p2*L2 == d
                integral = L1*p1*p1 + L2*p2*p2 + (h1*h1*L1**3+h2*h2*L2**3)/3
                minimum = xi*xi*y*y/(L1+L2)
                excess = delta*delta*(1/L1+1/L2)+(h1*h1*L1**3+h2*h2*L2**3)/3
                assert integral-minimum == excess
                assert excess >= 0
                assert (excess == 0) == (delta == h1 == h2 == 0)
        assert slope**2*(L1+L2)/2 == xi*xi*y*y/4


def test_angle_asymptotics_scaled_reference_radian_degree_cancellation():
    # An illustrative smooth log frequency, not a reduced spectral solution.
    # This check isolates units/reference scaling; it makes no tail claim.
    log_phi = lambda x: .3*x*x/(1+x*x)
    slope_xi = lambda x: .6*x/(1+x*x)**2
    xi_ref = np.deg2rad(5)/.005
    xi = 4.
    expected = slope_xi(xi)/(log_phi(xi_ref)/xi_ref)
    for eps in (.0025,.005,.01):
        ref_rad = eps*xi_ref
        s_rad = slope_xi(xi)/eps
        s0_rad = log_phi(xi_ref)/ref_rad
        assert s_rad/s0_rad == pytest.approx(expected)
        assert (s_rad*np.pi/180)/(s0_rad*np.pi/180) == pytest.approx(expected)


def test_leading_limit_basis_space_and_curvature_gram():
    from scripts.lib.isotropic_two_springs_leading import Ritz,fields
    from scripts.lib.reddy_symmetric_coupled_beams_ritz import shifted_legendre_clamped_basis
    m=Ritz(.1,10,order=8)
    for arm,L in enumerate(m.lengths):
        t=np.linspace(0,1,31);f=fields(m.order,t,L,arm)
        old=shifted_legendre_clamped_basis(m.order,t).values.T
        if arm==0:old=t[:,None]*old
        # Independent old basis spans exactly the same polynomial space.
        transform=np.linalg.lstsq(f[0],old,rcond=None)[0]
        np.testing.assert_allclose(f[0]@transform,old,atol=2e-13)
        expected=np.zeros(m.order)
        offset=int(arm==1)
        expected[offset:]=[1/(L**3*(2*n+1)) for n in range(m.order-offset)]
        sl=slice(arm*m.order,(arm+1)*m.order)
        np.testing.assert_allclose(m.Kb[sl,sl],np.diag(expected),atol=2e-14)
        assert np.max(abs(f[0][0]))<1e-15
        if arm==0:assert np.max(abs(f[1][0]))<1e-15


def test_leading_limit_energy_factor_springs_and_exact_constraint():
    from scripts.lib.isotropic_two_springs_leading import Ritz
    m=Ritz(.1,10,order=8);c=m.Z@np.random.default_rng(7).normal(size=m.Z.shape[1]);xi=4.
    bending=c@m.Kb@c;jump=m.jump@c;rotation=m.end_rotation@c;y=m.y@c
    assert c@m.K0@c == pytest.approx(bending+.1*jump**2+10*rotation**2)
    added_potential=.5*c@(xi**2/2*m.P)@c
    assert added_potential==pytest.approx(xi**2*y*y/4)
    assert abs((m.ends[0][0]+m.ends[1][0])@c)<1e-14
    assert c@m.P@c>=0
    assert np.linalg.matrix_rank(m.P)==1
    np.testing.assert_allclose(m.K0,m.K0.T,atol=1e-14)
    np.linalg.cholesky(m.Z.T@m.M@m.Z)


@pytest.mark.parametrize('pair',[(.1,10),(10,.1)])
def test_leading_limit_C_shape_independent_local_frequency(pair):
    from scripts.lib.isotropic_two_springs_leading import Ritz
    m=Ritz(*pair,order=20);r0=m.solve(0);r=m.solve(.01)
    c=r0['coefficients'][:,0];zeta=r0['zeta'][0]
    for scale in (1.,-3.,.02):
        v=scale*c
        C=(v@m.P@v)/(4*zeta*(v@m.M@v))
        assert C==pytest.approx(r0['C_at_zero'],rel=1e-13)
    local=np.log(r['omega_lim'][0]/r0['omega_lim'][0])/.01**2
    assert local==pytest.approx(C,rel=2e-5)
    assert r0['diagnostics']['gap12_relative']>1
    assert r0['diagnostics']['natural_boundary']<1e-9
    assert r0['diagnostics']['strong_ode']<1e-9


@pytest.mark.parametrize('pair',[(.1,10),(10,.1)])
def test_leading_limit_monotonicity_and_exact_fixed_node(pair):
    from scripts.lib.isotropic_two_springs_leading import Ritz
    m=Ritz(*pair,order=20);rows=[m.solve(x) for x in (0,5,30)]
    fixed=m.solve(0,node_fixed=True)
    assert np.all(np.diff([r['omega_lim'] for r in rows],axis=0)>=-1e-11)
    assert np.all(rows[-1]['omega_lim']<=fixed['omega_lim']*(1+1e-12))
    assert abs(fixed['node_value'][0])<1e-12
    assert fixed['diagnostics']['constraint_rank']==2
    assert fixed['diagnostics']['dimension']==2*m.order-2
    assert abs(rows[-1]['node_value'][0])<abs(rows[0]['node_value'][0])
