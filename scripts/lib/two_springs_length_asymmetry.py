"""Bounded length-asymmetry workflow; uses the unchanged two-spring solvers.

Not a runnable script or a new solver. The existing diagnostic CLI dispatches
this preset. Numerical identity pins the input, physical modules and reused
numerical procedures, separately from later reporting/layout edits.
"""
from dataclasses import asdict
import csv
import hashlib
import inspect
import json
import math
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT/'data/input/isotropic_two_springs_length_asymmetry.json'
BASE = ROOT/'results/isotropic_two_springs_rl_timoshenko'


def api():
    from scripts.analysis import verify_isotropic_two_springs_rl_timoshenko
    return verify_isotropic_two_springs_rl_timoshenko


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def read_csv(path):
    with path.open(encoding='utf8', newline='') as stream:
        return list(csv.DictReader(stream))


def close(a, b, tol=1e-10, atol=1e-13):
    return bool(np.allclose(a, b, rtol=tol, atol=atol))


def point_path(folder, case):
    identity = {**asdict(case), 'variant': 'RL_FULL'}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    return folder/'points'/(key+'.json')


def length_slopes(beta, y):
    """Only the fixed 51/86-node study grid; preserve archived join stencils."""
    a = api()
    beta = np.asarray(beta, dtype=float)
    y = np.asarray(y, dtype=float)
    main = sorted(set([i*.5 for i in range(7)]+[3+i*.25 for i in range(37)]+list(range(13,21))))
    full = main+list(range(22,91,2))
    require(list(beta) in (main, full) and len(beta) == len(y), 'Unexpected length-study grid')
    lookup = dict(zip(beta, y))
    result, names = [], []
    if len(beta) == 86:
        base = [b for b in beta if b not in (.5, 1.5, 2.5)]
        old, stencils = a.stiffness_slopes(base, [lookup[b] for b in base])
        archived = dict(zip(base, zip(old, stencils)))
    else:
        archived = {}
    for b in beta:
        if b == 0:
            value, name = 0., 'reflection_symmetry'
        elif b < 3:
            value, name = lookup[b+.5]-lookup[b-.5], 'centered_h0.5'
        elif b in archived:
            value, name = archived[b]
        elif b == 20:
            value, name = (3*lookup[20]-4*lookup[19]+lookup[18])/2, 'short_grid_backward_h1'
        else:
            h = .25 if 3 < b < 12 else 1.
            value, name = (lookup[b+h]-lookup[b-h])/(2*h), f'centered_h{h:g}'
        result.append(value)
        names.append(name)
    return np.array(result), names


class Experiment:
    def __init__(self):
        self.a = api()
        self.cfg = json.loads(CONFIG.read_text(encoding='utf8'))
        self.inputs = {}
        self.archived = {}
        self.manifests = {}
        for relative, digest in self.cfg['sources'].items():
            folder = BASE/relative
            self.check_hash(folder/'manifest.json', digest)
            m = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
            require(m['status'] == 'PASS', 'Unaccepted source archive')
            self.manifests[relative] = m
        numerical = [self.a.Run, self.a.tim_catalog, self.a.stiffness_slopes,
                     self.a.refined_pick, self.a.mechanics_energy, self.a.variation_constraints]
        physical = self.manifests['growth_stiffness_5x5/361e9a839ef5be30']['identity']['files']
        frozen = {}
        for name, digest in physical.items():
            if name.startswith('scripts/lib/'):
                self.check_hash(ROOT/name, digest)
                frozen[name] = digest
        from scripts.lib import isotropic_two_springs_leading as leading
        frozen[Path(leading.__file__).relative_to(ROOT).as_posix()] = self.a.sha(leading.__file__)
        identity = dict(config=self.cfg, physical_files=frozen,
                        numerical_procedures={f.__name__: hashlib.sha256(inspect.getsource(f).encode()).hexdigest()
                                              for f in numerical})
        self.key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
        self.folder = BASE/'length_asymmetry'/self.key
        self.folder.mkdir(parents=True, exist_ok=True)
        path = self.folder/'contract.json'
        if path.exists():
            require(json.loads(path.read_text(encoding='utf8')) == identity, 'Immutable contract differs')
        else:
            self.a.write(path, identity)
        self.run = self.a.Run(self.cfg, self.folder)
        catalogs = self.source_file('51af749660af1df5', 'pole_catalogs.json')
        self.run.catalogs = {float(k): v for k, v in json.loads(catalogs.read_text()).items()}
        if (self.folder/'pole_catalogs.json').exists():
            self.run.catalogs.update({float(k): v for k, v in json.loads((self.folder/'pole_catalogs.json').read_text()).items()})
        for relative in ('growth_stiffness_5x5/361e9a839ef5be30', 'epsilon_scaling/d925e14c6ff93605'):
            m = self.manifests[relative]['identity']['config']
            require(m['material'] == self.cfg['material'] and m['mu'] == .01
                    and m['variant'] == 'RL_FULL', 'Source physical contract mismatch')
            for r in read_csv(self.source_file(relative, 'spectrum.csv')):
                pair = r.get('pair', r.get('ID'))
                if pair not in self.cfg['cases'] or float(r.get('epsilon', .005)) != .005:
                    continue
                key = pair, float(r['beta_deg'])
                if key in self.archived:
                    old = self.archived[key]
                    require(close([float(r[f'omega{k}']) for k in (1,2,3)],
                                  [float(old[f'omega{k}']) for k in (1,2,3)], 1e-8), 'Archive overlap mismatch')
                else:
                    self.archived[key] = r
        self.started = time.perf_counter()

    def check_hash(self, path, expected=None):
        digest = self.a.sha(path)
        require(expected is None or expected == digest, 'Source hash mismatch: '+str(path))
        name = path.relative_to(ROOT).as_posix()
        require(name not in self.inputs or self.inputs[name] == digest, 'Source changed during execution')
        self.inputs[name] = digest
        return digest

    def source_file(self, relative, name):
        path = BASE/relative/name
        self.check_hash(path, self.manifests[relative]['artifacts'][name])
        return path

    def validate(self, rec, case):
        require(rec['case'] == {**asdict(case), 'variant': 'RL_FULL'}, 'Case identity mismatch')
        require(rec['status'] == rec['count']['status'] == 'PASS', 'Unaccepted spectrum')
        roots = rec['roots']
        p = self.run.policy
        w = np.array([r['omega'] for r in roots])
        require(len(w) == 3 and w[0] > 0 and min(np.diff(w)) > p['root_xtol'], 'Invalid sorted prefix')
        require(w[-1] < self.run.model.omega(p['Lambda_upper'])*(1-1e-3), 'Guard too close to search edge')
        for k, r in enumerate(roots, 1):
            d = r['diagnostics']
            require(r['sorted_position'] == k and r['counts'] == [k-1,k] and all(r['gates'].values()), 'Prefix count or gate')
            require(max(np.max(d['joint_scaled']), np.max(d['external_scaled'])) <= p['boundary_tol']
                    and d['ode_error'] <= p['ode_tol'] and d['energy_error'] <= p['energy_tol']
                    and d['singular_ratio'] <= p['singular_tol']
                    and d['nonzero_condition'] <= p['nonzero_condition_max']
                    and rec['mass_gram_error'] <= p['mass_gram_tol'], 'Saved physical checks exceed tolerance')
            require(close(r['Lambda']**2, r['omega']*self.run.model.time, 2e-14, 0), 'Frequency normalization')
        return w

    def point(self, pair, mu, beta, previous=None):
        case = self.a.physics.Case(float(beta), float(mu), *self.cfg['cases'][pair])
        if mu == .01:
            source = self.archived[pair, float(beta)]
            path = ROOT/source['point_file']
            self.check_hash(path, source['point_sha256'])
            rec = json.loads(path.read_text(encoding='utf8'))
            self.check_hash(path.with_suffix('.npz'), rec['shape_sha256'])
            origin = 'archive'
        else:
            require(mu in self.cfg['mu_values'] and beta in self.cfg['map']['parameter_grid']+self.cfg['tail_grid'],
                    'Unrequested spectral point')
            path = point_path(self.folder, case)
            rec = self.run.solve(case, neighbour=previous)
            origin = 'new_study'
            self.a.write(self.folder/'pole_catalogs.json', self.run.catalogs)
        self.validate(rec, case)
        return rec, path, origin

    def finish_execution(self, stage, status='PASS'):
        for name, digest in self.inputs.items():
            require(self.a.sha(ROOT/name) == digest, 'Protected source changed: '+name)
        path = self.folder/'executions.json'
        runs = json.loads(path.read_text()) if path.exists() else []
        runs.append(dict(stage=stage, status=status, elapsed_seconds=time.perf_counter()-self.started,
                         new_prefixes=self.run.new_points, new_structural_roots=self.run.new_roots,
                         new_pole_roots=self.run.pole_roots, ritz_solves=self.run.ritz_solves,
                         workflow_sha256=self.a.sha(Path(__file__)), cli_sha256=self.a.sha(self.a.__file__)))
        self.a.write(path, runs)
        old = json.loads((self.folder/'source_hashes.json').read_text()) if (self.folder/'source_hashes.json').exists() else {}
        old.update(self.inputs)
        self.a.write(self.folder/'source_hashes.json', old)


def preflight(exp):
    a, p, m, cfg = exp.a, exp.a.physics, exp.run.model, exp.cfg
    file = exp.folder/'preflight.json'
    if file.exists():
        result = json.loads(file.read_text())
        require(result['status'] == 'PASS', 'Previous preflight did not pass')
        return result
    structural = []
    for mu in cfg['mu_values']:
        case = p.Case(6, mu, .1, 10)
        lengths = case.lengths(m)
        require(close(lengths, [1-mu,1+mu]) and sum(lengths) == 2, 'Length or numbering contract')
        require(np.array_equal(p.external_matrix(m,case,0), np.eye(6)[:3]), 'Clamp moved from arm1')
        expected = np.eye(6)[:3].copy()
        expected[2] = [0,0,10*m.section.EI/m.l,0,0,-1]
        require(np.array_equal(p.external_matrix(m,case,1), expected), 'End spring moved from arm2')
        require(np.array_equal(p.joint_matrix(m,case), p.joint_matrix(m,p.Case(6,.01,.1,10))), 'Joint depends on ordering of lengths')
        errors = []
        for L in lengths:
            x = np.array([.1,.5,.9])*L
            w = float(m.omega(4))
            b, d = p.basis(m,L,w,x), p.basis(m,L,w,x,1)
            error = float(np.max(abs(d-np.einsum('ij,njk->nik',p.state_matrix(m,w),b)))/np.max(abs(d)))
            require(error <= cfg['policy']['ode_tol'], 'Fundamental solution ODE failure')
            errors.append(error)
        structural.append(dict(mu=mu,L1=lengths[0],L2=lengths[1],basis_ode_max=max(errors),
                               fixed_total_mass=2*m.section.rhoA*m.l,
                               k_joint=.1*m.section.EI/m.l,k_end=10*m.section.EI/m.l,status='PASS'))
    a.csvwrite(exp.folder/'geometry_checks.csv', structural)
    # Fresh reproductions explicitly requested in step A; kept separate from
    # immutable archives and not counted as new map nodes.
    reg = a.Run(cfg, exp.folder/'regression')
    reg.catalogs = exp.run.catalogs
    regression = []
    for pair, springs in cfg['cases'].items():
        previous = None
        for beta in cfg['regression_angles']:
            expected, _, _ = exp.point(pair,.01,beta)
            case = p.Case(beta,.01,*springs)
            rec = reg.solve(case, neighbour=previous)
            w = exp.validate(rec,case)
            previous = w
            error = float(max(abs(w/np.array([r['omega'] for r in expected['roots']])-1)))
            require(error <= cfg['policy']['same_physics_frequency_tol'], 'mu=.01 archive regression failed')
            regression.append(dict(pair=pair,beta_deg=beta,relative_frequency_error=error,
                                   **{f'Lambda{k}':r['Lambda'] for k,r in enumerate(rec['roots'],1)},status='PASS'))
    a.csvwrite(exp.folder/'regression.csv', regression)
    # Existing mu=0 straight CC/RIGID reference is used as a separate limiting
    # check, never as an expected spectrum for two unequal external supports.
    ref = json.loads(exp.source_file('51af749660af1df5','regression_diagnostics.json').read_text())
    case = p.Case(0.,0.)
    rec = reg.solve(case)
    w = exp.validate(rec,case)
    zero_error = float(max(abs(w/np.array(ref['single_rod_omega'])-1)))
    require(zero_error <= cfg['policy']['same_physics_frequency_tol'], 'Straight mu=0 limit failed')
    controls = []
    for mu in cfg['preflight_negative_mu']:
        for pair in cfg['cases']:
            previous = None
            for beta in cfg['preflight_angles']:
                start = time.perf_counter()
                rec, _, _ = exp.point(pair,mu,beta,previous)
                previous = np.array([r['omega'] for r in rec['roots']])
                controls.append(dict(pair=pair,mu=mu,beta_deg=beta,seconds=time.perf_counter()-start,
                                     Lambda=[r['Lambda'] for r in rec['roots']],status='PASS'))
    ritz = []
    for control in cfg['independent_ritz_cases']:
        for pair in cfg['cases']:
            rec, _, _ = exp.point(pair,control['mu'],control['beta'])
            ritz.append(exp.run.ritz_check(rec))
    a.write(exp.folder/'preflight_ritz.json',ritz)
    result = dict(status='PASS', geometry=structural, negative_mu_controls=controls,
                  regression_max=max(r['relative_frequency_error'] for r in regression),
                  regression_new_prefixes=reg.new_points, regression_new_roots=reg.new_roots,
                  mu0_rigid_limit_error=zero_error,
                  independent_ritz_max=max(r['relative_frequency_error'] for c in ritz for r in c['orders']),
                  median_control_seconds=float(np.median([r['seconds'] for r in controls])),
                  physical_modules_unchanged=True, sign_mu_spectral_equality_not_assumed=True)
    a.write(file,result)
    print(f'Preflight PASS: regression {result["regression_max"]:.3g}; negative mu 12 points; Ritz {result["independent_ritz_max"]:.3g}',flush=True)
    return result


def sweep(exp, angles, constructions, stage):
    a, cfg = exp.a, exp.cfg
    spectra, quality, neighbours = [], [], []
    audit_path = exp.folder/'triggered_count_checks.json'
    audits = json.loads(audit_path.read_text()) if audit_path.exists() else {}
    for pair, mu in constructions:
        previous = None
        for beta in angles:
            rec, path, origin = exp.point(pair, mu, beta, previous)
            w = np.array([r['omega'] for r in rec['roots']])
            if previous is not None:
                jump = float(max(abs(w/previous-1)))
                triggered = jump > cfg['policy']['neighbour_relative_jump']
                key = f'{pair}/{mu:g}/{beta:g}'
                if triggered and key not in audits:
                    case = a.physics.Case(beta,mu,*cfg['cases'][pair])
                    checks = []
                    queries = [float(exp.run.model.omega(cfg['policy']['Lambda_lower'])),
                               float((w[0]+w[1])/2),float((w[1]+w[2])/2)]
                    for expected, query in enumerate(queries):
                        _, count, diagnostic = a.physics.count_at(exp.run.model,case,query,
                                                                 exp.run.cats(exp.run.model,case),exp.run.policy)
                        require(count == expected, 'Triggered independent prefix count mismatch')
                        checks.append(dict(expected=expected,count=count,**diagnostic))
                    audits[key] = dict(status='PASS',checks=checks,
                                       resolution='count-one brackets + independent below/midgap counts + physical gates')
                    a.write(audit_path,audits)
                neighbours.append(dict(pair=pair,mu=mu,beta_deg=beta,jump=jump,triggered=triggered,
                                       resolution='independent_counts_and_physical_gates' if triggered else 'no_trigger'))
            previous = w
            row = dict(pair=pair,mu=mu,L1=1-mu,L2=1+mu,beta_deg=beta,
                       **{f'{f}{k}':r[f] for k,r in enumerate(rec['roots'],1) for f in ('omega','Omega','Lambda')},
                       gap12_relative=float(w[1]/w[0]-1),origin=origin,
                       point_file=path.relative_to(ROOT).as_posix(),point_sha256=a.sha(path),shape_sha256=rec['shape_sha256'])
            spectra.append(row)
            for k,r in enumerate(rec['roots'],1):
                d = r['diagnostics']
                quality.append(dict(pair=pair,mu=mu,beta_deg=beta,sorted_position=k,
                                    physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),
                                    ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],
                                    mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],
                                    count_lower=r['counts'][0],count_upper=r['counts'][1],status='PASS'))
        a.csvwrite(exp.folder/(stage+'_spectrum.csv'),spectra)
        a.csvwrite(exp.folder/(stage+'_quality.csv'),quality)
        a.csvwrite(exp.folder/(stage+'_neighbours.csv'),neighbours)
        print(f'{stage}: {pair} mu={mu:g}, {len(angles)} accepted prefixes',flush=True)
    return spectra, quality


def profiles(exp, spectra, complete):
    cfg, a = exp.cfg, exp.a
    profiles, transitions = [], []
    for pair in cfg['cases']:
        for mu in cfg['mu_values']:
            rows = sorted((r for r in spectra if r['pair']==pair and r['mu']==mu),key=lambda r:r['beta_deg'])
            if not rows:
                continue
            beta = [r['beta_deg'] for r in rows]
            w = np.array([r['omega1'] for r in rows])
            y = np.log(w)
            slope, names = length_slopes(beta,y)
            gain = float(y[beta.index(5)]-y[0])
            valid = gain > cfg['minimum_reference_log_gain']
            s0deg = gain/5
            s0rad = gain/math.radians(5)
            q = abs(slope)/s0deg if valid else np.full(len(beta),np.nan)
            if valid:
                normalized, _ = length_slopes(beta,np.log([r['Lambda1'] for r in rows]))
                require(close(q,abs(normalized)/(s0deg/2),1e-11), 'q frequency-unit identity failed')
                tail = np.maximum.accumulate(q[::-1])[::-1]
                if complete:
                    transition, expected = a.refined_pick(beta,q,cfg['q_star'],cfg['minimum_tail_deg'])
                    require(close(tail,expected), 'Tail envelope mismatch')
                else:
                    hits = np.flatnonzero(tail <= cfg['q_star'])
                    transition = dict(status='PROVISIONAL_ONLY' if len(hits) else 'NO_PROVISIONAL_TRANSITION',
                                      beta_pl_deg=None,provisional_angle_deg=float(beta[hits[0]]) if len(hits) else None)
            else:
                tail = np.full(len(beta),np.nan)
                transition = dict(status='REFERENCE_GAIN_NOT_RESOLVED',beta_pl_deg=None)
            transitions.append(dict(pair=pair,mu=mu,S0_per_radian=s0rad,reference_log_gain=gain,
                                    tail_to_90=complete,**transition))
            for i,b in enumerate(beta):
                profiles.append(dict(pair=pair,mu=mu,beta_deg=b,signed_s_per_degree=float(slope[i]),
                                     signed_s_per_radian=float(slope[i]*180/math.pi),
                                     S0_per_radian=s0rad,q=float(q[i]) if valid else None,
                                     Q_tail=float(tail[i]) if valid else None,stencil=names[i],tail_to_90=complete))
    return profiles, transitions


def frequency_comparison(exp, spectra):
    rows, summary, differences = [], [], []
    for pair in exp.cfg['cases']:
        for mu in exp.cfg['mu_values']:
            table = sorted((r for r in spectra if r['pair']==pair and r['mu']==mu),key=lambda r:r['beta_deg'])
            w0 = table[0]['omega1']
            lam0 = table[0]['Lambda1']
            for r in table:
                f = r['omega1']/w0
                require(close(f,(r['Lambda1']/lam0)**2,2e-14,0),'Relative frequency conversion')
                rows.append({**r,'F_omega':f})
            local = [r for r in rows if r['pair']==pair and r['mu']==mu]
            min_gap = min(local,key=lambda r:r['gap12_relative'])
            summary.append(dict(pair=pair,mu=mu,Lambda0=lam0,
                **{f'Lambda{b}':next(r['Lambda1'] for r in local if r['beta_deg']==b) for b in (3,6,10,20)},
                **{f'F{b}':next(r['F_omega'] for r in local if r['beta_deg']==b) for b in (3,5,10,20)},
                gap12_min=min_gap['gap12_relative'],gap12_min_beta_deg=min_gap['beta_deg']))
    by = {(r['pair'],r['mu'],r['beta_deg']):r for r in rows}
    for pair in exp.cfg['cases']:
        for magnitude in (.01,.1):
            for beta in exp.cfg['map']['parameter_grid']:
                positive, negative = by[pair,magnitude,beta],by[pair,-magnitude,beta]
                dl = positive['Lambda1']/negative['Lambda1']-1
                df = positive['F_omega']/negative['F_omega']-1
                differences.append(dict(pair=pair,abs_mu=magnitude,beta_deg=beta,
                                        D_Lambda=dl,D_F=df,D_Lambda_percent=100*dl,D_F_percent=100*df,
                                        D_omega_percent=100*((1+dl)**2-1)))
    exp.a.csvwrite(exp.folder/'spectrum.csv',rows)
    exp.a.csvwrite(exp.folder/'summary.csv',summary)
    exp.a.csvwrite(exp.folder/'signed_comparison.csv',differences)
    return rows, summary, differences


def mechanics(exp, spectra, transition_rows):
    a, cfg, m = exp.a, exp.cfg, exp.run.model
    rows, fd, changes = [], [], []
    for pair in cfg['cases']:
        for mu in cfg['mu_values']:
            table = {r['beta_deg']:r for r in spectra if r['pair']==pair and r['mu']==mu}
            selected = set(cfg['mechanics_angles'])
            transition = next(t for t in transition_rows if t['pair']==pair and t['mu']==mu)
            if transition.get('beta_pl_deg') is not None:
                selected.add(transition['beta_pl_deg'])
            for beta in sorted(selected):
                record = table[beta]
                point = ROOT/record['point_file']
                require(a.sha(point)==record['point_sha256'] and a.sha(point.with_suffix('.npz'))==record['shape_sha256'],
                        'Mechanics source point changed')
                with np.load(point.with_suffix('.npz')) as data:
                    arms = [{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')}
                            for i in (1,2)]
                case = a.physics.Case(beta,mu,*cfg['cases'][pair])
                energy = a.mechanics_energy(m,case,arms)
                mass, V = energy['mass'], energy['total']
                require(abs(mass-1)<5e-12 and abs(2*V/(record['omega1']**2*mass)-1)<cfg['policy']['energy_tol'],
                        'Mechanics full mass or energy failure')
                ends = np.array([arm['ends'][1] for arm in arms])
                u,w,_,n,q,_ = ends[0]
                wn,wq = n*w,-q*u
                work = wn+wq
                s = work/(2*V)
                C,Cb = a.variation_constraints(math.radians(beta))
                geometric = Cb@ends[:,:2].ravel()
                reaction = -ends[0,3:5]
                require(close(geometric,[-w,u],atol=2e-11), 'Variation geometric identity')
                require(close(float(reaction@geometric),work,atol=1e-12), 'Variation work identity')
                if beta == 0:
                    require(abs(s) < 1e-8, 'Reflection slope at zero not reproduced')
                if beta in cfg['variation_fd_angles']:
                    y = {b:math.log(r['omega1']) for b,r in table.items()}
                    estimates = [(y[beta+h]-y[beta-h])/(2*math.radians(h)) for h in (.25,.5)]
                    errors = [abs(value-s) for value in estimates]
                    require(errors[0] < errors[1] and errors[0] <= abs(estimates[1]-estimates[0])+1e-8,
                            'Variational derivative not confirmed by FD refinement')
                    fd.append(dict(pair=pair,mu=mu,beta_deg=beta,s_rad=s,FD_h025=estimates[0],FD_h05=estimates[1],
                                   relative_error_h025=errors[0]/abs(s),relative_error_h05=errors[1]/abs(s),status='PASS'))
                rows.append(dict(pair=pair,mu=mu,beta_deg=beta,omega1=record['omega1'],mass=mass,
                                 V=V,u1=float(u),w1=float(w),N1=float(n),Q1=float(q),
                                 W_N=float(wn),W_Q=float(wq),W_beta=float(work),s_rad=float(s),
                                 s_per_degree=float(s*math.pi/180),d_joint=float(math.hypot(u,w)),
                                 reaction_norm=float(np.linalg.norm(reaction)),
                                 cos_gamma=float(work/max(np.linalg.norm(reaction)*np.linalg.norm(geometric),1e-300)),
                                 W_N_fraction=float(wn/work) if beta>0 else None,
                                 W_Q_fraction=float(wq/work) if beta>0 else None,
                                 mass_RL=float(energy['mass_components'][:,3].sum()),
                                 **{f'eta_{key}':float(value) for key,value in zip(('axial','bend','shear','joint','end'),energy['fractions'])},
                                 point_file=record['point_file']))
            local = {r['beta_deg']:r for r in rows if r['pair']==pair and r['mu']==mu}
            first,last = local[3],local[10]
            changes.append(dict(pair=pair,mu=mu,d_reduction_3_10=first['d_joint']/last['d_joint'],
                                slope_reduction_3_10=abs(first['s_rad']/last['s_rad']),
                                reaction_reduction_3_10=first['reaction_norm']/last['reaction_norm'],
                                V_growth_3_10=last['V']/first['V'],
                                cos_ratio_3_10=last['cos_gamma']/first['cos_gamma']))
    a.csvwrite(exp.folder/'mechanics.csv',rows)
    a.csvwrite(exp.folder/'variation_fd.csv',fd)
    a.csvwrite(exp.folder/'mechanical_factors.csv',changes)
    return rows, fd


def leading_curvature(exp):
    from scripts.lib.isotropic_two_springs_leading import Ritz
    a, cfg = exp.a, exp.cfg['leading']
    out, convergence, checks, finite = [], [], [], []
    source = read_csv(exp.source_file('leading_validation/f1c670fe62102196','curvature.csv'))
    cache = exp.folder/'leading_points'
    cache.mkdir(exist_ok=True)
    new = 0
    def solve(model,xi,pair):
        nonlocal new
        ident = dict(pair=pair,mu=model.mu,order=model.order,quadrature=model.quadrature,xi=xi)
        key = hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
        path = cache/(key+'.json')
        if path.exists():
            rec = json.loads(path.read_text())
            require(rec['identity']==ident, 'Leading cached identity')
            return rec['solution']
        result = model.solve(xi)
        new += 1
        a.write(path,dict(identity=ident,solution=result))
        return result
    for pair, springs in exp.cfg['cases'].items():
        for mu in exp.cfg['mu_values']:
            solutions = {}
            for order in cfg['orders']:
                model = Ritz(*springs,mu=mu,order=order,quadrature=cfg['quadrature'])
                require(sum(model.lengths)==2 and model.matrix_checks['P_rank']==1, 'Leading geometry/stiffness contract')
                solutions[order] = solve(model,0.,pair)
            rec = solutions[20]
            diag = rec['diagnostics']
            q = cfg['quality']
            gates = {name:diag[name]<=q[tol] for name,tol in
                     [('symmetry','symmetry_tol'),('eigen_residual','eigen_residual_tol'),
                      ('energy','energy_tol'),('mass_gram','mass_gram_tol'),('constraint','constraint_tol'),
                      ('natural_boundary','natural_boundary_tol'),('strong_ode','strong_ode_tol'),('essential','essential_tol')]}
            require(all(gates.values()), f'Leading quality violation: {pair}/{mu}: {diag}')
            d_w = float(max(abs(np.array(rec['omega_lim'])/np.array(solutions[16]['omega_lim'])-1)))
            d_C = abs(rec['C_at_zero']/solutions[16]['C_at_zero']-1)
            require(d_w<=q['frequency_convergence_rtol'] and d_C<=q['C_convergence_rtol'], 'Leading spatial convergence')
            fine = Ritz(*springs,mu=mu,order=20,quadrature=cfg['quadrature_check'])
            fine_rec = solve(fine,0.,pair)
            quad = float(max(abs(np.array(fine_rec['omega_lim'])/np.array(rec['omega_lim'])-1)))
            require(quad<=q['frequency_convergence_rtol'], 'Leading quadrature convergence')
            c = np.array(rec['coefficients'])[:,0]
            mass = float(c@model.M@c)
            zeta, y = float(rec['zeta'][0]),float(rec['node_value'][0])
            C = y*y/(4*zeta*mass)
            require(close(C,rec['C_at_zero'],1e-12), 'Leading C normalization')
            if mu == .01:
                old = next(r for r in source if r['ID']==pair and int(r['order'])==20)
                require(close(C,float(old['C']),1e-10), 'Leading archived C regression')
            for xi in cfg['xi_local'][1:]:
                local = solve(model,xi,pair)
                estimate = math.log(float(local['omega_lim'][0])/float(rec['omega_lim'][0]))/xi**2
                require(abs(estimate/C-1)<1e-4, 'Leading local curvature identity')
                finite.append(dict(pair=pair,mu=mu,xi=xi,C_from_log_ratio=estimate,C_from_shape=C,
                                   relative_difference=estimate/C-1))
            out.append(dict(pair=pair,mu=mu,lambda1=1-mu,lambda2=1+mu,zeta0=zeta,
                            omega_lim0=math.sqrt(zeta),y0=y,mass=mass,C=C,
                            axial_compliance=2*exp.run.model.l/exp.run.model.section.EA,
                            k_geom_over_beta_rad_squared=exp.run.model.section.EA/(2*exp.run.model.l),
                            leading_added_stiffness_coefficient=1/sum(model.lengths)))
            convergence.append(dict(pair=pair,mu=mu,frequency_16_20=d_w,C_16_20=d_C,quadrature_64_96=quad))
            checks.append(dict(pair=pair,mu=mu,**diag,status='PASS'))
    for name, rows in [('leading_curvature',out),('leading_convergence',convergence),
                        ('leading_quality',checks),('leading_local_check',finite)]:
        a.csvwrite(exp.folder/(name+'.csv'),rows)
    a.write(exp.folder/'leading_accounting.json',dict(new_this_execution=new,
            saved_eigenproblems=len(list(cache.glob('*.json'))),root_count_per_problem=3,
            full_model_solves=0,scope='C(mu), spatial/quadrature convergence and two local leading xi'))
    print(f'Leading curvature PASS: {len(out)} constructions; {new} new small eigenproblems',flush=True)
    return out


def derivative_refinement(exp, fd):
    """Triggered derivative check, separate from the fixed criterion grid."""
    suspects = [r for r in fd if r['relative_error_h025'] > .01]
    authorization = dict(trigger='FD h=.25 relative error against variation exceeds 1%',
                         cases=suspects, steps_deg=[.125,.0625,.03125],
                         target_relative_error=.002, maximum_new_points_per_case=6,
                         changes_to_criterion_grid=False)
    exp.a.write(exp.folder/'derivative_refinement_authorization.json',authorization)
    run = exp.a.Run(exp.cfg,exp.folder/'derivative_refinement')
    run.catalogs = exp.run.catalogs
    rows = []
    for r in suspects:
        last_error = r['relative_error_h025']
        for h in authorization['steps_deg']:
            values, paths = [], []
            for beta in (r['beta_deg']-h,r['beta_deg']+h):
                case = exp.a.physics.Case(beta,r['mu'],*exp.cfg['cases'][r['pair']])
                rec = run.solve(case)
                values.append(exp.validate(rec,case)[0])
                paths.append(point_path(run.folder,case).relative_to(ROOT).as_posix())
            estimate = math.log(values[1]/values[0])/(2*math.radians(h))
            error = abs(estimate/r['s_rad']-1)
            require(error < last_error, 'Triggered FD refinement fails to converge')
            rows.append(dict(pair=r['pair'],mu=r['mu'],beta_deg=r['beta_deg'],h_deg=h,
                             variation_per_radian=r['s_rad'],FD_per_radian=estimate,
                             relative_error=error,minus_point=paths[0],plus_point=paths[1]))
            last_error = error
            if error <= authorization['target_relative_error']:
                break
        require(last_error <= authorization['target_relative_error'], 'Derivative check remains unresolved')
    if rows:
        exp.a.csvwrite(exp.folder/'derivative_refinement.csv',rows)
    exp.a.write(exp.folder/'derivative_refinement_accounting.json',dict(
        new_prefixes=run.new_points,new_roots=run.new_roots,
        saved_prefixes=len(list((run.folder/'points').glob('*.json'))),status='PASS'))
    return rows


def leading_saved_quality(exp):
    """Reject bad production points; label low-order convergence probes explicitly."""
    q = exp.cfg['leading']['quality']
    rows = []
    for path in sorted((exp.folder/'leading_points').glob('*.json')):
        record = json.loads(path.read_text())
        identity, result = record['identity'], record['solution']
        d = result['diagnostics']
        names = ('symmetry','eigen_residual','energy','mass_gram','constraint',
                 'natural_boundary','strong_ode','essential')
        accepted = all(d[k] <= q[k+'_tol'] for k in names)
        production = identity['order'] == 20
        require(accepted or not production, 'Leading production residual violation: '+str(identity))
        rows.append(dict(**identity,**{k:d[k] for k in names},
                         role='production' if production else 'convergence_probe',
                         status='PASS' if accepted else 'DIAGNOSTIC_NOT_ACCEPTED'))
    exp.a.csvwrite(exp.folder/'leading_all_quality.csv',rows)
    return rows


FIGURES = ('01_frequency','02_normalized_frequency','03_signed_differences',
           '04_characteristic_angle','05_leading_curvature')


def render(folder):
    """Draw saved scalar rows only; no spectral calls or data interpolation."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.labelsize':11,'legend.fontsize':9,'pdf.fonttype':42,
                         'axes.linewidth':.7,'lines.linewidth':1.5,'savefig.dpi':320})
    comma = FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    pair_labels = [r'$(\kappa_{\mathrm{j}};\,\kappa_{\mathrm{e}})=(0{,}1;\,10)$',
                   r'$(\kappa_{\mathrm{j}};\,\kappa_{\mathrm{e}})=(10;\,0{,}1)$']
    pairs = ['J0E4','J4E0']
    mus = [-.1,-.01,0,.01,.1]
    colors = ['#155b86','#d57919','#222222','#6f587e','#367b64']
    lines = ['-',(0,(5,2)),':',(0,(5,2,1,2)),(0,(3,1))]
    markers = ['o','s','x','^','D']
    def arr(rows,key):return np.array([float(r[key]) for r in rows])
    def axes_style(ax,xlabel=True):
        ax.grid(True,color='.88',linewidth=.5)
        ax.set_axisbelow(True)
        ax.xaxis.set_major_formatter(comma)
        ax.yaxis.set_major_formatter(comma)
        if xlabel:ax.set_xlabel(r'$\beta$, град.')
    def save(fig,name):
        fig.savefig(folder/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(folder/(name+'.png'),dpi=320)
        plt.close(fig)
    spectrum = read_csv(folder/'spectrum.csv')
    for name,key,ylabel in zip(FIGURES[:2],('Lambda1','F_omega'),(r'$\Lambda_1$',r'$F(\beta;\mu)=\omega_1(\beta;\mu)/\omega_1(0;\mu)$')):
        fig, axs = plt.subplots(2,1,figsize=(7.0,6.3),sharex=True,sharey=True)
        for panel,(ax,pair,label) in enumerate(zip(axs,pairs,pair_labels)):
            for mu,color,line,marker in zip(mus,colors,lines,markers):
                rr = sorted([r for r in spectrum if r['pair']==pair and float(r['mu'])==mu],key=lambda r:float(r['beta_deg']))
                ax.plot(arr(rr,'beta_deg'),arr(rr,key),color=color,ls=line,marker=marker,
                        markersize=3.4,markevery=7,mfc='white',mew=.8,
                        label=r'$\mu='+f'{mu:g}'.replace('.','{,}')+'$')
            axes_style(ax,panel==1)
            ax.set_ylabel(ylabel)
            ax.set_xlim(0,20)
            ax.set_xticks([0,5,10,15,20])
            ax.text(.025,.91,('а' if panel==0 else 'б')+')  '+label,transform=ax.transAxes)
        handles,labels = axs[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',ncol=5,frameon=False,bbox_to_anchor=(.53,1.0),columnspacing=1.)
        fig.subplots_adjust(left=.12,right=.98,bottom=.09,top=.92,hspace=.13)
        save(fig,name)
    data = read_csv(folder/'signed_comparison.csv')
    fig,axs = plt.subplots(2,2,figsize=(7.2,6.2),sharex=True,sharey='row')
    for col,pair in enumerate(pairs):
        for row,key in enumerate(('D_Lambda_percent','D_F_percent')):
            ax = axs[row,col]
            for m,color,line,marker in ((.01,colors[0],'-','o'),(.1,colors[1],'--','s')):
                rr = [r for r in data if r['pair']==pair and float(r['abs_mu'])==m]
                ax.plot(arr(rr,'beta_deg'),arr(rr,key),color=color,ls=line,marker=marker,
                        markevery=7,markersize=3.4,mfc='white',label=r'$|\mu|='+str(m).replace('.','{,}')+'$')
            axes_style(ax,row==1)
            ax.axhline(0,color='.4',lw=.6)
            ax.set_xlim(0,20)
            ax.set_xticks([0,5,10,15,20])
            ax.text(.04,.1,'абвг'[row*2+col]+')',transform=ax.transAxes)
            if col==0:ax.set_ylabel(r'$100D_\Lambda$, %' if row==0 else r'$100D_F$, %')
            if row==0:ax.text(.5,1.025,pair_labels[col],transform=ax.transAxes,ha='center')
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.12,right=.98,bottom=.09,top=.87,hspace=.12,wspace=.1)
    save(fig,FIGURES[2])
    for name,filename,key,ylabel in ((FIGURES[3],'transitions','beta_pl_deg',r'$\beta_{\mathrm{pl}}$, град.'),
                                    (FIGURES[4],'leading_curvature','C',r'$C(\mu)$')):
        data = read_csv(folder/(filename+'.csv'))
        fig,ax=plt.subplots(figsize=(6.8,3.8))
        for pair,label,color,line,marker in zip(pairs,pair_labels,colors[:2],('-','--'),('o','s')):
            rr=sorted([r for r in data if r['pair']==pair and (key!='beta_pl_deg' or r['status']=='PASS')],key=lambda r:float(r['mu']))
            ax.plot(arr(rr,'mu'),arr(rr,key),color=color,ls=line,marker=marker,mfc='white',label=label)
        axes_style(ax,False)
        ax.set_xlabel(r'$\mu$')
        ax.set_ylabel(ylabel)
        ax.set_xticks([-.1,-.05,0,.05,.1])
        ax.set_xlim(-.107,.107)
        if key=='beta_pl_deg':ax.yaxis.set_major_locator(MultipleLocator(.25))
        ax.legend(loc='best',frameon=True,facecolor='white',edgecolor='.85')
        fig.subplots_adjust(left=.13,right=.98,bottom=.17,top=.97)
        save(fig,name)


def finalize(folder):
    a = api()
    stopped = folder/'stopped.json'
    if stopped.exists():
        event = json.loads(stopped.read_text())
        if event['status']=='STOPPED':
            require('Unknown symbol' in event['error'] and 'mathrm' in event['error'],
                    'Unresolved earlier failure prevents finalization')
            event.update(status='RESOLVED_TECHNICAL_LAYOUT_ERROR',scientific_gate_failure=False,
                         resolution='Braced mathtext roman subscripts; no numerical data changed.')
            a.write(stopped,event)
    bundle_readme(folder)
    a.write(folder/'manifest.json',dict(schema='two-springs-length-asymmetry-results-v1',
        status='PASS',numerical_summary=json.loads((folder/'numerical_summary.json').read_text()),
        source_hashes=json.loads((folder/'source_hashes.json').read_text()),
        workflow_files={p.relative_to(ROOT).as_posix():a.sha(p) for p in (CONFIG,Path(__file__),Path(a.__file__))},
        artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in sorted(folder.rglob('*'))
                   if p.is_file() and p.name != 'manifest.json'},
        figures=list(FIGURES),spectrum_semantics='sorted_positions',calculation_mode='fast_plot'))


def verify_bundle(folder, figures=True):
    a = api()
    m = json.loads((folder/'manifest.json').read_text())
    require(m['status']=='PASS', 'Unaccepted bundle')
    for name,digest in m['artifacts'].items():
        if not figures and Path(name).suffix in ('.pdf','.png'):
            continue
        require(a.sha(folder/name)==digest, 'Changed result: '+name)
    for name,digest in {**m['source_hashes'],**m['workflow_files']}.items():
        require(a.sha(ROOT/name)==digest, 'Changed input/workflow: '+name)
    return m


def bundle_readme(folder):
    text = '''# Влияние асимметрии длин: данные и рисунки

Научная интерпретация, численные таблицы и ограничения:
[отдельный отчёт](../../../../docs/theory/isotropic_two_springs_length_asymmetry.md).
Исходный отчёт о двух пружинах и прежние расчётные архивы не изменены.

## Постановка и происхождение

Каноническая RL_FULL/Тимошенко, две точные поступательные связи,
две поворотные пружины; внешняя пружина всегда на стержне 2.
ε=0,005, b/h=4, ν=0,3, k_s=5/6, l=E=ρ=1.
L1=1−μ, L2=1+μ; μ=−0,1; −0,01; 0; 0,01; 0,1.
J0E4 означает (κj;κe)=(0,1;10), J4E0 — (10;0,1).
Параметры, H1–H5 и допуски зафиксированы до расчёта в `contract.json`;
исходный файл — `data/input/isotropic_two_springs_length_asymmetry.json`.

Основная сетка: 0:0,5:3; 3:0,25:12; 13:1:20 градусов, 51 узел.
Все десять серий дали предварительный переход; только после этого
добавлен хвост 22:2:90. `tail_authorization.json` сохраняет основания.
688 новых точек карты и 172 архивных, всего 860 трёхчастотных наборов.
Отдельно выполнены 11 регрессионных, 6 адресных производных точек,
8 проверок Ритца полной модели и 50 задач ведущей модели.
Ведущие порядки 16 — только пробы сходимости; их сильные невязки
не проходят итоговый допуск. Все 40 результатов порядка 20 проходят.

Неизменённые архивы внутри `results/isotropic_two_springs_rl_timoshenko/`:

- `growth_stiffness_5x5/361e9a839ef5be30/`: спектры μ=0,01 и численная политика;
- `epsilon_scaling/d925e14c6ff93605/`: дополнительные узлы μ=0,01, ε=0,005;
- `51af749660af1df5/`: каталоги полюсов и прямой жёсткий предел μ=0;
- `leading_validation/f1c670fe62102196/`: контроль C при μ=0,01.

`source_hashes.json` содержит SHA-256 каждого использованного источника.
`manifest.json` проверяет источники, код и все файлы нового пакета.
Архивные поля читаются по сохранённым путям, не копируются и не перезаписываются.

## Таблицы, единицы и проверки

Все CSV — UTF-8, разделитель запятая, десятичная точка. Рисунки используют запятую.
`pair`, `mu`, `beta_deg` однозначно задают конструкцию и угол.

| Файл | Содержание и нормировка |
| --- | --- |
| `spectrum.csv` | 510 основных точек; omega1..3 в принятой размерной нормировке; Omega=omega·sqrt(ρAl⁴/EI); Lambda=sqrt(Omega); F_omega=omega1/omega1(0) каждой конструкции |
| `complete_spectrum.csv` | 860 основных и хвостовых точек; gap12_relative=(omega2−omega1)/omega1; путь и хеш каждого JSON/NPZ, происхождение archive/new_study |
| `signed_comparison.csv` | D_Lambda=Lambda(+m)/Lambda(−m)−1, D_F=F(+m)/F(−m)−1; отдельные колонки в процентах |
| `summary.csv` | Частоты, F на контрольных углах, минимальный g12 на основной сетке |
| `profiles.csv`, `transitions.csv` | Производная на градус и радиан, q, максимум последующего хвоста Q; подтверждённый дискретный beta_pl |
| `preliminary_*.csv` | Предварительная диагностика до 20°, без утверждения beta_pl |
| `mechanics.csv` | Полная массовая нормировка M=1, включая инерцию Рэлея—Лява и поворотов; u,w,N,Q,V, Nw,−Qu, d_joint, вариационная производная на радиан и градус |
| `mechanical_factors.csv` | Отношения перемещений, реакций, энергии и производной при 3° и 10° |
| `variation_fd.csv`, `derivative_refinement.csv` | Сходимость конечной разности к вариационной производной; адресное уточнение не меняет q-сетку |
| `leading_curvature.csv` | C=y0²/(4ζ0M0), M0=1; безразмерные ξ=β_rad/ε и ζ=Lambda⁴; постоянная суммарная продольная податливость |
| `leading_convergence.csv`, `leading_local_check.csv` | Порядки 16/20, квадратуры 64/96; ln(F_lim)/ξ² при ξ=0,01 и 0,02 |
| `leading_all_quality.csv` | Невязки всех 50 задач, явное отделение проб порядка 16 от принятых результатов порядка 20 |
| `main_quality.csv`, `tail_quality.csv` | Физические, энергетические, спектральные и массовые проверки каждой из трёх позиций |
| `geometry_checks.csv`, `regression.csv`, `preflight*.json` | Нумерация стержней, закрепления, отрицательные μ, эталоны и независимый метод Ритца |

Первые три частоты — отсортированные положительные позиции, не отслеживаемые
ветви. Полнота проверяется энергетическим счётом и скачком k−1→k у каждого корня.
Для 14 срабатываний соседнего контроля сохранены дополнительные проверки счёта.
Используется `frequency-map-v1`, `fast_plot`, K_plot=2, K_guard=3.
На рисунках показана только первая позиция; вторая сохранена для контроля сближения.

Для q сохранены прежние центральные разности и переходные трафареты:
h=0,25° внутри 3–12°, h=1° на стыках 3°,12° и до 19°;
в 20° узлы 19°,20°,22°; далее h=2°, на 90° обратная разность второго порядка.
На 0–3° h=0,5°; при нуле использована отражательная симметрия s=0.
Опора 0–5°, q*=0,15 и минимальный хвост 30° неизменны.
Линии на рисунках соединяют расчётные значения; сглаживания и добавления
интерполированных спектральных точек нет. β_pl определяется по узлам до 90°,
а не по непрерывному максимуму между ними.

## Рисунки и предлагаемые подписи

Все пять рисунков: векторный PDF с внедрёнными шрифтами и PNG 320 dpi.

1. `01_frequency`: **Первая собственная частота при различном распределении длин.**
   Зависимости Lambda1 от угла при пяти μ; а — (0,1;10), б — (10;0,1).
   Внешняя поворотная пружина принадлежит второму стержню при обоих знаках μ.
2. `02_normalized_frequency`: **Относительное изменение первой частоты с углом.**
   Каждая зависимость нормирована на собственное значение omega1(0;μ).
   Обозначения частей и параметры соответствуют рисунку 1.
3. `03_signed_differences`: **Сопоставление конструкций с противоположными знаками μ.**
   Верхний ряд — 100D_Lambda, нижний — 100D_F; слева (0,1;10), справа (10;0,1).
   Отрицательный знак означает меньшее значение у конструкции с пружиной
   на длинном стержне относительно конструкции с пружиной на коротком.
4. `04_characteristic_angle`: **Характерный угол по принятому дискретному критерию.**
   Опорный интервал 0–5°, q*=0,15; все последующие расчётные узлы проверены до 90°.
   Линии соединяют пять исследованных значений μ; особая механическая точка не утверждается.
5. `05_leading_curvature`: **Коэффициент начального квадратичного роста ведущей частоты.**
   ln(F_lim)=C(μ)ξ²+O(ξ⁴), ξ=β_rad/ε; C вычислен по массово-нормированной
   форме при ξ=0. Коэффициент дополнительной узловой жёсткости ξ²/2 неизменен;
   меняются изгибная задача и её собственная форма.

## Воспроизведение

Из корня репозитория, Python с NumPy, SciPy, Matplotlib:

```powershell
python scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py --preset length-asymmetry --stage preflight
python scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py --preset length-asymmetry
python scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py --preset length-asymmetry --verify-cache BUNDLE
python scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py --preset length-asymmetry --plot-only BUNDLE
python -m pytest tests/test_two_springs_length_asymmetry.py -q
```

Заменить BUNDLE на `results/isotropic_two_springs_rl_timoshenko/length_asymmetry/KEY`.
В этой среде использован `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`.
Полный запуск повторно использует имеющиеся спектральные точки; `--plot-only`
и `--verify-cache` читают только сохранённые данные и не вызывают решатели.
Каталог зависит от входного контракта, физических файлов и используемых
численных процедур; изменения оформления не меняют идентификатор исследования.
`executions.json` фиксирует версии сценария и число новых решений каждого запуска.

Данные локальные, каталог результатов исключён из Git. Научные выводы и их
ограничения сохранены также в отдельном отслеживаемом Markdown-отчёте.
'''
    (folder/'README.md').write_text(text.replace('KEY',folder.name),encoding='utf8')


def main(args):
    require(not (args.verify_cache and args.plot_only), 'Choose verify-cache or plot-only')
    require(not args.reuse_points_from and args.plot_set=='original', 'Unsupported length-asymmetry option')
    if args.verify_cache or args.plot_only:
        folder = (args.verify_cache or args.plot_only).resolve()
        require(folder.parent==BASE/'length_asymmetry', 'Wrong result family')
        verify_bundle(folder,figures=not bool(args.plot_only))
        if args.plot_only:
            render(folder)
            verify_bundle(folder)
        print('Verified saved length-asymmetry data; 0 eigenvalue solves',flush=True)
        return
    require(args.stage in ('all','preflight'), 'Use stage all or preflight for this preset')
    exp = Experiment()
    print('Length-asymmetry bundle: '+str(exp.folder),flush=True)
    try:
        preflight(exp)
        if args.stage == 'preflight':
            exp.finish_execution('preflight')
            return
        constructions = [(pair,mu) for pair in exp.cfg['cases'] for mu in exp.cfg['mu_values']]
        main_spectrum, main_quality = sweep(exp,exp.cfg['map']['parameter_grid'],constructions,'main')
        main_spectrum, summary, differences = frequency_comparison(exp,main_spectrum)
        preliminary, transitions = profiles(exp,main_spectrum,False)
        exp.a.csvwrite(exp.folder/'preliminary_profiles.csv',preliminary)
        exp.a.csvwrite(exp.folder/'preliminary_transitions.csv',transitions)
        eligible = [(r['pair'],r['mu']) for r in transitions if r['status']=='PROVISIONAL_ONLY']
        exp.a.write(exp.folder/'tail_authorization.json',dict(
            trigger='positive resolved 0-5 gain and persistent provisional crossing on 0-20',
            eligible=[dict(pair=pair,mu=mu) for pair,mu in eligible],decisions=transitions))
        print(f'Provisional crossings: {len(eligible)}/10; only these receive the fixed tail to 90 degrees',flush=True)
        tail_spectrum, tail_quality = sweep(exp,exp.cfg['tail_grid'],eligible,'tail') if eligible else ([],[])
        complete_spectrum = main_spectrum+tail_spectrum
        complete_rows = [r for r in complete_spectrum if (r['pair'],r['mu']) in eligible]
        full_profiles, full_transitions = profiles(exp,complete_rows,True)
        full_transitions += [r for r in transitions if (r['pair'],r['mu']) not in eligible]
        exp.a.csvwrite(exp.folder/'profiles.csv',full_profiles)
        exp.a.csvwrite(exp.folder/'transitions.csv',full_transitions)
        exp.a.csvwrite(exp.folder/'complete_spectrum.csv',complete_spectrum)
        mech, fd = mechanics(exp,complete_spectrum,full_transitions)
        derivative_refinement(exp,fd)
        leading = leading_curvature(exp)
        leading_saved_quality(exp)
        exp.a.write(exp.folder/'numerical_summary.json',dict(status='PASS',
            main_points=len(main_spectrum),tail_points=len(tail_spectrum),
            archived_map_points=sum(r['origin']=='archive' for r in complete_spectrum),
            new_map_points=sum(r['origin']=='new_study' for r in complete_spectrum),
            accepted_characteristic_angles=sum(r['status']=='PASS' for r in full_transitions),
            quality_max={key:max(r[key] for r in main_quality+tail_quality)
                         for key in ('physical','ode','energy','singular','mass_gram','condition')},
            variation_fd_max=max(r['relative_error_h025'] for r in fd),
            mechanics_mass_error=max(abs(r['mass']-1) for r in mech),leading_constructions=len(leading)))
        exp.finish_execution('all')
        render(exp.folder)
        finalize(exp.folder)
        print('Numerical stages and five figures complete: '+str(exp.folder),flush=True)
    except Exception as error:
        exp.a.write(exp.folder/'stopped.json',dict(status='STOPPED',error=repr(error),
                                                  workflow_sha256=exp.a.sha(Path(__file__))))
        exp.finish_execution(args.stage,'STOPPED')
        raise
