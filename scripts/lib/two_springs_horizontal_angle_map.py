"""HS-A: missing-only full-model angle map and all bracketed threshold crossings.

Stable model-specific workflow, reusing the frozen HS-P1 physics, validation,
prefix counter and Run. HS-P1's fixed 105-point contract and hashed archive
cannot safely acquire this composite-grid/all-crossings transaction contract.
No physical equations or spectral solver are copied here.
"""
from pathlib import Path
from dataclasses import asdict
import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import os
import platform
import time

import numpy as np
import scipy
from scipy.optimize import brentq
from scripts.lib import two_springs_horizontal_pilot as pilot

h, verified = pilot.physics, pilot.verified
ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'data/input/isotropic_two_springs_horizontal_angle_map.json'
OUTPUT = ROOT / 'results/isotropic_two_springs_rl_timoshenko/horizontal_spring_angle_map'
sha, require, csvwrite = pilot.sha, pilot.require, pilot.csvwrite


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    if path.exists() and path.read_text(encoding='utf8') == payload:
        return
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(payload, encoding='utf8')
    # Windows readers/antivirus can briefly hold a destination handle.
    for attempt in range(8):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 7:raise
            time.sleep(.05*(attempt+1))


def construction_id(j, k):
    return f'J{j:g}_S{k if k == h.RIGID else format(k, "g")}'


def base_nodes(k, cfg):
    """Only numerical duplicates merge; keep every independent node origin."""
    eps = cfg['material']['epsilon']
    scale = math.sqrt(2*eps**2 + (0 if k == h.RIGID else 1/k))
    raw = [(float(b), 'base_2deg') for b in range(0, 91, 2)]
    raw += [(float(b), 'endpoint_0.5deg') for b in np.arange(80, 90.25, .5)]
    raw += [(math.degrees(scale*math.sqrt(t)), f'scale_{t:g}') for t in cfg['grid']['tau']]
    return merge_nodes(raw, cfg['grid']['duplicate_absolute_deg'])


def merge_nodes(raw, tol):
    nodes = []
    for beta, origin in sorted(raw):
        if not 0 <= beta <= 90:
            continue
        if nodes and abs(beta-nodes[-1]['beta_deg']) <= tol:
            if origin not in nodes[-1]['origins']:
                nodes[-1]['origins'].append(origin)
        else:
            nodes.append({'beta_deg': float(beta), 'origins': [origin]})
    return nodes


def high_segments(angles, values, target):
    """Connected runs on the ACTUAL grid; no continuous guarantee or first-hit shortcut."""
    flags = [v >= target for v in values]
    runs = []
    for i, yes in enumerate(flags):
        if yes and (i == 0 or not flags[i-1]):
            runs.append([angles[i], angles[i]])
        if yes:
            runs[-1][1] = angles[i]
    if not flags[-1]:
        return runs, None, 'NO_TERMINAL_HIGH_RANGE'
    start = runs[-1][0]
    status = 'ALL_CHECKED_NODES_HIGH' if start == 0 else (
        'TERMINAL_POINT_ONLY' if start == 90 else 'TERMINAL_HIGH_RANGE')
    return runs, start, status


def crossing_brackets(angles, values, target):
    out = []
    for a, b, wa, wb in zip(angles[:-1], angles[1:], values[:-1], values[1:]):
        if (wa-target)*(wb-target) < 0:
            out.append((a, b, 'up' if wb > wa else 'down'))
    return out


class Study:
    def __init__(self):
        self.cfg = read(CONFIG)
        self.p = self.cfg['policy']
        self.ident = {'config': self.cfg, 'config_sha256': sha(CONFIG),
                      'workflow_sha256': sha(__file__),
                      'initial_sha256': sha(ROOT/self.cfg['initial_workspace']),
                      'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__}}
        # Numerical identity is frozen by the contract, physics and versions.
        # Orchestration fixes preserve accepted points; each source revision is archived.
        numerical_identity = {k:v for k,v in self.ident.items() if k != 'workflow_sha256'}
        self.run_id = hashlib.sha256(json.dumps(numerical_identity, sort_keys=True).encode()).hexdigest()[:16]
        self.folder = OUTPUT/self.run_id
        self.folder.mkdir(parents=True, exist_ok=True)
        snapshot = self.folder/'workflow_versions'/(self.ident['workflow_sha256']+'.txt')
        snapshot.parent.mkdir(exist_ok=True)
        if not snapshot.exists():snapshot.write_bytes(Path(__file__).read_bytes())
        self.model = h.Model(**self.cfg['material'])
        self.constructions = [dict(id=construction_id(j,k), k_joint=j, kappa_horizontal=k, chi_s=chi)
                              for j in self.cfg['joint_stiffnesses']
                              for k,chi in zip(self.cfg['horizontal_kappa'], self.cfg['horizontal_chi'])]
        self.reference = {}
        self.audit_sources()
        old = read(ROOT/'data/input/isotropic_two_springs_horizontal_spring.json')
        self.engine = verified.Run(old, self.folder/'staging')
        self.entries = read(self.folder/'entries.json') if (self.folder/'entries.json').exists() else {}
        self.events = read(self.folder/'events.json') if (self.folder/'events.json').exists() else []
        self.failures = read(self.folder/'failures.json') if (self.folder/'failures.json').exists() else []
        self.independent = read(self.folder/'independent.json') if (self.folder/'independent.json').exists() else {}
        self.new_this_run = 0

    def case(self, c, beta):
        return h.Case(float(beta), self.cfg['mu'], c['k_joint'], c['kappa_horizontal'])

    def audit_sources(self):
        for name, digest in self.cfg['sources'].items():
            require(sha(ROOT/name) == digest, 'Frozen source changed: '+name)
        checked = 0
        for name, digest in self.cfg['archives'].items():
            folder = ROOT/name
            require(sha(folder/'manifest.json') == digest, 'Source manifest changed')
            manifest = read(folder/'manifest.json')
            config = manifest['identity']['config']
            require(config['material'] == self.cfg['material'] and config['mu'] == self.cfg['mu'], 'Source physical identity')
            require(config.get('policy_full', config['policy'] if 'policy' in config else None) == self.p, 'Source tolerances')
            for artifact, digest in manifest['artifacts'].items():
                require(sha(folder/artifact) == digest, 'Changed source artifact '+artifact)
                checked += 1
            paths = list(folder.glob('points/*.json')) + list(folder.glob('computed/points/*.json')) + list(folder.glob('full/points/*.json'))
            if (folder/'spectra.csv').exists():
                with (folder/'spectra.csv').open(encoding='utf8', newline='') as f:
                    paths += [ROOT/r['point_file'] for r in csv.DictReader(f) if r.get('point_file')]
            for path in paths:
                rec = read(path)
                if 'case' not in rec or rec.get('status') != 'PASS':
                    continue
                case = h.Case(**rec['case'])
                if case.k_joint not in self.cfg['joint_stiffnesses'] or case.kappa_horizontal not in self.cfg['horizontal_kappa']:
                    continue
                pilot.validate_record(rec, case, self.p)
                require(path.with_suffix('.npz').exists(), 'Source shape bundle missing')
                key = pilot.key(case)
                if key in self.reference:
                    require(np.max(abs(np.array(rec['omega'])/self.reference[key][0]['omega']-1)) < self.p['ritz_frequency_tol'], 'Duplicate archive frequencies disagree')
                self.reference[key] = (rec, path)
        atomic(self.folder/'source_audit.json', {'status':'PASS', 'source_artifacts_checked':checked,
               'compatible_unique_points':len(self.reference), 'physics':'same RL_FULL/Timoshenko, one global horizontal spring',
               'new_asymptotic_assumptions':False})

    def prepare(self):
        grid = {}
        for c in self.constructions:
            raw = [(n['beta_deg'], o) for n in base_nodes(c['kappa_horizontal'], self.cfg) for o in n['origins']]
            for beta,j,k in self.reference:
                if j == c['k_joint'] and k == c['kappa_horizontal']:
                    raw.append((beta, 'accepted_archive'))
            grid[c['id']] = merge_nodes(raw, self.cfg['grid']['duplicate_absolute_deg'])
        if (self.folder/'grid.json').exists():
            require(read(self.folder/'grid.json') == grid, 'Frozen production grid changed')
        else:
            atomic(self.folder/'grid.json', grid)
            atomic(self.folder/'contract.json', self.ident)
        return grid

    def flush(self):
        atomic(self.folder/'entries.json', self.entries)
        atomic(self.folder/'events.json', self.events)
        atomic(self.folder/'independent.json', self.independent)
        atomic(self.folder/'failures.json', self.failures)

    def point(self, c, beta, reason, previous=None):
        case = self.case(c, beta)
        key = self.engine.name(case)
        if key in self.entries:
            e = self.entries[key]
            require(sha(ROOT/e['path']) == e['sha256'] and sha(ROOT/e['npz_path']) == e['npz_sha256'], 'Point checksum')
            rec = read(ROOT/e['path'])
            pilot.validate_record(rec, case, self.p)
            if reason not in e['reasons']:
                e['reasons'].append(reason)
                self.flush()
            return rec
        match = self.reference.get(pilot.key(case))
        start = time.perf_counter()
        upper = self.p['Lambda_upper']
        if match:
            rec, path = match
            origin = 'accepted_archive'
        else:
            path = self.folder/'points'/f'{key}.json'
            if path.exists() and path.with_suffix('.npz').exists():
                rec = read(path)
                origin = 'new_HS_A'
            else:
                if previous is not None:
                    candidate = min(upper, float(self.model.Lambda(previous['omega'][2]*1.3)))
                    _, count, diag = h.count_at(self.model,case,float(self.model.omega(candidate)),self.engine.catalogs,self.p)
                    if count >= 3:
                        upper = candidate
                self.engine.p = {**self.p, 'Lambda_upper':upper}
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        rec,_ = self.engine.solve(case, 'HS-A', ritz=False)
                    pilot.validate_record(rec,case,self.p)
                except Exception as exc:
                    self.failures.append({'id':c['id'],'beta':beta,'reason':reason,'error':repr(exc),'status':'UNRESOLVED','series_stopped':True})
                    self.flush()
                    raise
                stage = self.engine.folder/'points'/f'{key}.json'
                path.parent.mkdir(parents=True,exist_ok=True)
                os.replace(stage.with_suffix('.npz'),path.with_suffix('.npz'))
                os.replace(stage,path)  # accepted JSON is published LAST
                self.engine.solved.pop(key,None)  # bounded field cache
                origin = 'new_HS_A'
                self.new_this_run += 1
        pilot.validate_record(rec,case,self.p)
        self.entries[key] = {'id':c['id'],'beta_deg':float(beta),'case':asdict(case),'origin':origin,
                             'path':path.relative_to(ROOT).as_posix(),'sha256':sha(path),
                             'npz_path':path.with_suffix('.npz').relative_to(ROOT).as_posix(),'npz_sha256':sha(path.with_suffix('.npz')),
                             'reasons':[reason], 'seconds':time.perf_counter()-start,'upper_Lambda':upper}
        self.flush()
        return rec

    def independent_check(self,c,beta,reason,orders=(20,)):
        key = self.engine.name(self.case(c,beta))
        tag = key+'_'+','.join(map(str,orders))
        if tag in self.independent:
            require(self.independent[tag]['status']=='PASS','Unresolved prior independent check')
            return self.independent[tag]
        rec = self.point(c,beta,reason)
        checks = []
        for order in orders:
            w,_,diag = h.ritz(self.model,self.case(c,beta),order,rank_tol=self.p['kernel_rank_relative_tol'])
            err = float(max(abs(w/rec['omega']-1)))
            row = {'order':order,'omega':w.tolist(),'frequency_error':err,**diag}
            checks.append(row)
            self.independent[tag]={'id':c['id'],'beta':beta,'reason':reason,'checks':checks,'status':'UNRESOLVED'}
            self.flush()
            require(err <= self.p['ritz_frequency_tol'],'Independent Ritz frequency')
            for metric,tol in [('mass_gram_error','mass_gram_tol'),('energy_error','energy_tol'),('constraint_residual','boundary_tol'),('matrix_residual','singular_tol')]:
                require(diag[metric]<=self.p[tol], 'Independent Ritz '+metric)
            require(diag['nullity']==h.expected_nullity(self.case(c,beta)), 'Ritz zero count')
        result={'id':c['id'],'beta':beta,'reason':reason,'checks':checks,'status':'PASS'}
        self.independent[tag]=result
        self.flush()
        return result

    def preflight(self):
        grid = self.prepare()
        if (self.folder/'preflight.json').exists():
            require(read(self.folder/'preflight.json')['status']=='PASS','Unaccepted preflight')
            return
        start = time.perf_counter()
        rows=[]
        for j,k,b in self.cfg['preflight']['cases']:
            c=next(c for c in self.constructions if c['k_joint']==j and c['kappa_horizontal']==k)
            beta=next(n['beta_deg'] for n in grid[c['id']] if b in n['origins']) if isinstance(b,str) else float(b)
            rec=self.point(c,beta,'preflight')
            check=self.independent_check(c,beta,'preflight contrasts')
            rows.append({'id':c['id'],'beta':beta,'omega':rec['omega'],'independent':check})
        regress=[]
        for j,k,beta in self.cfg['preflight']['regression_cases']:
            c=next(c for c in self.constructions if c['k_joint']==j and c['kappa_horizontal']==k)
            ref=self.reference[pilot.key(self.case(c,beta))][0]
            self.engine.p=self.p
            with contextlib.redirect_stdout(io.StringIO()):
                rec,_=self.engine.solve(self.case(c,beta),'archive regression',ritz=False)
            err=float(max(abs(np.array(rec['omega'])/ref['omega']-1)))
            require(err<=self.p['ritz_frequency_tol'],'Archive reproduction')
            regress.append({'case':asdict(self.case(c,beta)),'relative_frequency_error':err})
        seconds=time.perf_counter()-start
        new=[e for e in self.entries.values() if e['origin']=='new_HS_A']
        atomic(self.folder/'preflight.json',{'status':'PASS','checks':rows,'regressions':regress,
               'regression_recalculations':len(regress),'seconds':seconds,
               'mean_new_prefix_seconds':float(np.mean([e['seconds'] for e in new])),
               'frozen_grid_points':sum(map(len,grid.values()))})

    def run_grid(self):
        self.preflight()
        grid=self.prepare()
        for c in self.constructions:
            if any(f['id']==c['id'] and f['status']=='UNRESOLVED' for f in self.failures):
                continue
            previous=None
            try:
                for node in grid[c['id']]:
                    previous=self.point(c,node['beta_deg'],'frozen_grid',previous)
                print(c['id'],len(grid[c['id']]),'grid PASS',flush=True)
            except ArithmeticError as exc:
                if not any(f['id']==c['id'] and f['status']=='UNRESOLVED' for f in self.failures):
                    self.failures.append({'id':c['id'],'status':'UNRESOLVED','stage':'grid locator','error':repr(exc),'series_stopped':True})
                    self.flush()
                print('SERIES STOP',c['id'],str(exc),flush=True)
        self.export()

    def series(self,c):
        return sorted([(e,read(ROOT/e['path'])) for e in self.entries.values() if e['id']==c['id']],key=lambda x:x[0]['beta_deg'])

    def extrema(self,c):
        """A bounded diagnostic grid around sampled turning points, not optimization."""
        pairs=self.series(c); b=[e['beta_deg'] for e,r in pairs]; w=[r['omega'][0] for e,r in pairs]
        turns=[]
        for i in range(1,len(w)-1):
            if (w[i]-w[i-1])*(w[i+1]-w[i])<0:
                turns.append((b[i-1],b[i],b[i+1],'maximum' if w[i]>w[i-1] else 'minimum'))
        for left,mid,right,kind in turns:
            # Grid roots very near a threshold are not useful extrema probes.
            if right-left < .25:
                continue
            event={'id':c['id'],'type':'extremum_refinement','kind':kind,'bracket':[left,right],'candidate':mid}
            if any(e==event for e in self.events):
                continue
            self.events.append(event);self.flush()
            for beta in np.arange(math.ceil(left*4)/4,right,.25):
                if left<float(beta)<right:
                    self.point(c,float(beta),'extremum_'+kind)
            pairs=self.series(c)
            local=[(e,r) for e,r in pairs if left<=e['beta_deg']<=right]
            e,r=(max if kind=='maximum' else min)(local,key=lambda x:x[1]['omega'][0])
            for beta in (e['beta_deg'],90.):
                self.independent_check(c,beta,'confirmed '+kind+' versus 90 degrees',orders=(16,20))

    def neighbour_audit(self,c):
        pairs=self.series(c)
        for (ea,a),(eb,b) in zip(pairs[:-1],pairs[1:]):
            jump=max(abs(np.array(a['omega'][:2])-b['omega'][:2])/np.maximum(a['omega'][:2],b['omega'][:2]))
            if jump>self.cfg['frequency_map']['neighbour_jump_relative']:
                self.independent_check(c,eb['beta_deg'],'neighbour jump >0.35')
        for e,r in pairs:
            gap=r['omega'][1]/r['omega'][0]-1
            if gap<self.cfg['analysis']['close_gap_relative']:
                self.independent_check(c,e['beta_deg'],'close first-second relative gap')

    def thresholds(self,c):
        output=[]
        # Discovery may add roots; re-evaluate grid maximum until stable.
        for iteration in range(self.cfg['analysis']['max_refinement_rounds']):
            pairs=self.series(c); b=[e['beta_deg'] for e,r in pairs]; w=[r['omega'][0] for e,r in pairs]
            require(b[0]==0 and b[-1]==90,'Incomplete angular series')
            maximum=max(w);output=[]
            for normalization,reference in [('omega90',w[-1]),('grid_max',maximum)]:
                for q in self.cfg['analysis']['q']:
                    target=q*reference
                    crossings=[]
                    def f(beta):return self.point(c,float(beta),'threshold_'+normalization+'_'+str(q))['omega'][0]/target-1
                    for left,right,direction in crossing_brackets(b,w,target):
                        root=float(brentq(f,left,right,xtol=self.cfg['analysis']['root_xtol_deg'],rtol=self.cfg['analysis']['root_rtol']))
                        half=self.cfg['analysis']['bracket_width_deg']/2
                        lo,hi=max(0.,root-half),min(90.,root+half)
                        fl,fr=f(lo),f(hi);residual=f(root)
                        status='PASS' if fl*fr<0 and min(abs(fl),abs(fr))>self.cfg['analysis']['sign_resolution_relative'] and abs(residual)<=self.cfg['analysis']['threshold_residual_tol'] else 'UNRESOLVED'
                        crossings.append({'beta_deg':root,'direction':direction,'bracket':[lo,hi],'signs':[fl,fr],'equation_residual':residual,'status':status})
                    # Exact sampled equality is retained, never hidden as a missing crossing.
                    for beta,omega in zip(b,w):
                        if omega==target and not any(abs(beta-x['beta_deg'])<1e-6 for x in crossings):
                            crossings.append({'beta_deg':beta,'direction':'sampled_equality','bracket':[beta,beta],'signs':[0.,0.],'equation_residual':0.,'status':'EXACT_SAMPLED_EQUALITY'})
                    now=self.series(c);angles=[e['beta_deg'] for e,r in now];values=[r['omega'][0] for e,r in now]
                    runs,high,status=high_segments(angles,values,target)
                    # The terminal-root estimate and conservative confirmed grid start are distinct.
                    tail=None if high is None else (0. if high==0 else max((x['beta_deg'] for x in crossings if x['direction']=='up'),default=None))
                    if any(x['status']=='UNRESOLVED' for x in crossings):status='UNRESOLVED_CROSSING'
                    output.append({'id':c['id'],'normalization':normalization,'q':q,'omega_target':reference,
                                   'threshold_omega':target,'high_at_zero':values[0]>=target,
                                   'terminal_beta_root_deg':tail,'terminal_beta_grid_deg':high,'status':status,
                                   'crossings':sorted(crossings,key=lambda x:x['beta_deg']),'grid_intervals':runs,
                                   'subsequent_checked_nodes':sum(a>=high for a in angles) if high is not None else 0})
            if max(r['omega'][0] for e,r in self.series(c))<=maximum:
                break
            self.events.append({'id':c['id'],'type':'new_maximum_during_thresholds','old_maximum':maximum,'round':iteration})
        else:
            raise ArithmeticError('Threshold maximum did not stabilize')
        atomic(self.folder/'thresholds'/f"{c['id']}.json",output)
        return output

    def analyze(self):
        if (self.folder/'manifest.json').exists():
            done=read(self.folder/'manifest.json')
            if done['status']=='PASS':
                for name,digest in done['artifacts'].items():
                    require(sha(self.folder/name)==digest,'Completed numerical artifact '+name)
                print('REUSED completed HS-A; zero new spectral solves',self.folder,flush=True)
                return
        all_thresholds=[]
        for c in self.constructions:
            if any(f['id']==c['id'] and f['status']=='UNRESOLVED' for f in self.failures):continue
            try:
                self.extrema(c)
                self.neighbour_audit(c)
                rows=self.thresholds(c)
                self.neighbour_audit(c)
                all_thresholds+=rows
                print(c['id'],'analysis',[(r['normalization'],r['q'],r['terminal_beta_root_deg']) for r in rows],flush=True)
            except Exception as exc:
                self.failures.append({'id':c['id'],'status':'UNRESOLVED','stage':'analysis','error':repr(exc),'series_stopped':True})
                self.flush()
                print('ANALYSIS STOP',c['id'],repr(exc),flush=True)
        atomic(self.folder/'thresholds.json',all_thresholds)
        csvwrite(self.folder/'thresholds.csv',[{**{k:v for k,v in r.items() if k not in ('crossings','grid_intervals')},
                     'crossings_json':json.dumps(r['crossings']),'intervals_json':json.dumps(r['grid_intervals'])} for r in all_thresholds])
        self.export()
        self.final_checks()

    def export(self):
        rows=[];summary=[]
        for c in self.constructions:
            pairs=self.series(c)
            if not pairs:continue
            for e,r in pairs:
                ds=r['diagnostics'];omega=np.array(r['omega']);lam=self.model.Lambda(omega)
                rows.append({'id':c['id'],'kappa_j':c['k_joint'],'chi_s':c['chi_s'],'kappa_s':c['kappa_horizontal'],
                             'beta_deg':e['beta_deg'],**{f'omega{i+1}':float(v) for i,v in enumerate(omega)},
                             **{f'Omega{i+1}':float(v*v) for i,v in enumerate(lam)},**{f'Lambda{i+1}':float(v) for i,v in enumerate(lam)},
                             'g12':float(omega[1]/omega[0]-1),'g23':float(omega[2]/omega[1]-1),'zero_count':r['count']['zero_gate']['expected'],
                             'boundary_max':max(max(np.max(d['joint_scaled']),np.max(d['external_scaled'])) for d in ds),
                             'ode_max':max(d['ode_error'] for d in ds),'energy_max':max(d['energy_error'] for d in ds),
                             'mass_gram_error':r['mass_gram_error'],'singular_max':max(d['singular_ratio'] for d in ds),
                             'origin':e['origin'],'reasons':';'.join(e['reasons']),'point_file':e['path'],'status':'PASS'})
            i=max(range(len(pairs)),key=lambda i:pairs[i][1]['omega'][0]);e,r=pairs[i]
            w0=pairs[0][1]['omega'][0];w90=pairs[-1][1]['omega'][0] if pairs[-1][0]['beta_deg']==90 else None
            decreases=[{'left':ea['beta_deg'],'right':eb['beta_deg'],'relative_change':b['omega'][0]/a['omega'][0]-1}
                       for (ea,a),(eb,b) in zip(pairs[:-1],pairs[1:]) if b['omega'][0]<a['omega'][0]]
            atomic(self.folder/'decreases'/f"{c['id']}.json",decreases)
            summary.append({'id':c['id'],'kappa_j':c['k_joint'],'chi_s':c['chi_s'],'kappa_s':c['kappa_horizontal'],
                            'omega0':w0,'omega90':w90,'omega_grid_max':r['omega'][0],'beta_max_grid':e['beta_deg'],
                            'max_left_node':pairs[i-1][0]['beta_deg'] if i else None,'max_right_node':pairs[i+1][0]['beta_deg'] if i+1<len(pairs) else None,
                            'drop_to_90_percent':100*(1-w90/r['omega'][0]) if w90 else None,'gain_max_percent':100*(r['omega'][0]/w0-1),
                            'decreasing_steps_raw':len(decreases),'decreasing_steps_resolved':sum(d['relative_change'] < -self.cfg['analysis']['change_trigger_relative'] for d in decreases),
                            'minimum_gap12':min(rr['omega'][1]/rr['omega'][0]-1 for ee,rr in pairs),'nodes':len(pairs)})
        csvwrite(self.folder/'spectra.csv',rows);csvwrite(self.folder/'summary.csv',summary)
        return rows,summary

    def final_checks(self):
        rows,summary=self.export()
        neighbours=[];limits=[];straight=[]
        for c in self.constructions:
            rr=sorted([r for r in rows if r['id']==c['id']],key=lambda r:r['beta_deg'])
            for a,b in zip(rr[:-1],rr[1:]):
                jump=max(abs(a[f'omega{i}']-b[f'omega{i}'])/max(a[f'omega{i}'],b[f'omega{i}']) for i in [1,2])
                neighbours.append({'id':c['id'],'left':a['beta_deg'],'right':b['beta_deg'],'jump_relative':jump,
                                   'triggered':jump>self.cfg['frequency_map']['neighbour_jump_relative']})
            for beta in [0.,90.]:
                case=self.case(c,beta);ext=h.external_matrix(self.model,case,1)
                p,spring=h.trace_map(self.model,case);co,si=h.cs(case)
                if c['kappa_horizontal']==h.RIGID:
                    err=float(np.max(abs(co*p[6]+si*p[7])))
                    require(p.shape[1]==6 and err<=self.p['boundary_tol'],'Exact horizontal trace constraint')
                    rank=1
                else:
                    eig=np.linalg.eigvalsh(ext[:2,:2]);rank=int(np.linalg.matrix_rank(ext[:2,:2]))
                    require(rank==1 and min(eig)>=-1e-20 and p.shape[1]==7,'Finite horizontal spring rank/space')
                    err=0.
                require(np.array_equal(ext[2],np.eye(6)[5]),'No external rotational spring')
                limits.append({'id':c['id'],'beta':beta,'trace_dimension':p.shape[1],
                               'horizontal_rank_or_constraint':rank,'constraint_error':err,'status':'PASS'})
        for j in self.cfg['joint_stiffnesses']:
            zero=[r['omega1'] for r in rows if r['kappa_j']==j and r['beta_deg']==0]
            error=max(zero)/min(zero)-1
            require(error<=self.p['ritz_frequency_tol'],'Straight first-frequency invariance')
            straight.append({'kappa_j':j,'relative_spread_omega1':error,'status':'PASS'})
        atomic(self.folder/'neighbour_audit.json',neighbours)
        atomic(self.folder/'limit_checks.json',{'cardinal_cases':limits,'straight_bending':straight})
        checks=[]
        for axis in ['kappa_j','kappa_s']:
            other='kappa_s' if axis=='kappa_j' else 'kappa_j'
            groups={}
            for r in rows:groups.setdefault((r[other],r['beta_deg']),[]).append(r)
            for key,group in groups.items():
                group.sort(key=lambda r:math.inf if r[axis]==h.RIGID else r[axis])
                for a,b in zip(group[:-1],group[1:]):
                    relative=[b[f'omega{i}']/a[f'omega{i}']-1 for i in [1,2,3]]
                    checks.append({'axis':axis,'low':a['id'],'high':b['id'],'beta':a['beta_deg'],'relative_changes':relative,
                                   'status':'PASS' if min(relative)>=-self.p['ritz_frequency_tol'] else 'FAIL'})
        atomic(self.folder/'stiffness_monotonicity.json',checks)
        initial=read(ROOT/self.cfg['initial_workspace'])
        preserved=all(sha(ROOT/p)==digest for p,digest in initial['archives'].items())
        thresholds=read(self.folder/'thresholds.json');regress=[]
        source=next(ROOT/p for p in self.cfg['archives'] if 'horizontal_spring_leading' in p)
        with (source/'thresholds.csv').open(encoding='utf8',newline='') as f:old=list(csv.DictReader(f))
        for r in old:
            k=20000 if r['id']=='H5' else h.RIGID
            new=next(x for x in thresholds if x['id']==construction_id(.1,k) and x['q']==float(r['q']) and x['normalization']=='omega90')
            delta=new['terminal_beta_root_deg']-float(r['beta_full_deg'])
            regress.append({'id':r['id'],'q':float(r['q']),'difference_deg':delta,'status':'PASS' if abs(delta)<=2e-6 else 'FAIL'})
        atomic(self.folder/'threshold_regression.json',regress)
        ok=(not self.failures and len(thresholds)==72 and all(c['status']=='PASS' for c in checks+regress)
            and all(r['status']!='UNRESOLVED_CROSSING' for r in thresholds) and preserved)
        grid=read(self.folder/'grid.json');completed={(r['id'],r['beta_deg']) for r in rows}
        require(all((name,n['beta_deg']) in completed for name,nodes in grid.items() for n in nodes),'Frozen grid coverage')
        manifest={'status':'PASS' if ok else 'PARTIAL','identity':self.ident,'grid_sha256':sha(self.folder/'grid.json'),
                  'HEAD':initial['HEAD'],'initial_dirty_files':len(initial['files']),'old_archive_files_preserved':preserved,
                  'accepted_unique_points':len(self.entries),'new_unique_points':sum(e['origin']=='new_HS_A' for e in self.entries.values()),
                  'reused_unique_points':sum(e['origin']=='accepted_archive' for e in self.entries.values()),
                  'preflight_regression_recalculations':2,'independent_ritz_solves':sum(len(c['checks']) for c in self.independent.values()),
                  'frozen_grid_points':sum(map(len,grid.values())), 'threshold_added_points':len(rows)-sum(map(len,grid.values())),
                  'neighbour_pairs':len(neighbours),'neighbour_triggers':sum(r['triggered'] for r in neighbours),
                  'minimum_g12':min(r['g12'] for r in rows),'minimum_g23':min(r['g23'] for r in rows),
                  'spectrum_seconds':sum(e['seconds'] for e in self.entries.values() if e['origin']=='new_HS_A'),
                  'max_residuals':{k:max(r[k] for r in rows) for k in ['boundary_max','ode_max','energy_max','mass_gram_error','singular_max']},
                  'artifacts':{p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*')
                               if p.is_file() and p.name not in ['manifest.json','figures_manifest.json','README.md','completion_checks.json','documentation_manifest.json']
                               and not any(part in p.relative_to(self.folder).parts for part in ['staging','figures'])}}
        atomic(self.folder/'manifest.json',manifest)
        print(manifest['status'],self.folder,manifest['new_unique_points'],'new',manifest['reused_unique_points'],'reused',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','preflight','map','analyze','verify','plot_only'])
    parser.add_argument('--folder',type=Path,help='Completed bundle for pure plot_only')
    args=parser.parse_args()
    if args.action=='plot_only':
        if args.folder is None:parser.error('plot_only requires --folder')
        from scripts.lib.two_springs_horizontal_angle_figures import render
        render(args.folder)
        return
    s=Study()
    if args.action=='prepare':
        grid=s.prepare();print(s.folder,sum(map(len,grid.values())),'frozen points')
    elif args.action=='preflight':s.preflight();print(read(s.folder/'preflight.json'))
    elif args.action=='map':s.run_grid()
    elif args.action=='analyze':s.analyze()
    elif args.action=='verify':
        m=read(s.folder/'manifest.json')
        for name,digest in m['artifacts'].items():require(sha(s.folder/name)==digest,'Artifact '+name)
        print(m['status'],s.folder)


if __name__=='__main__':
    main()
