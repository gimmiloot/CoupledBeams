"""Bounded H/L/L/H EB/RLB and length-asymmetry mechanism check.

New workflow contract: physical two-arm mass, cross-model seed comparisons,
six sparse controls and two local questions. The old angular tracker assumes
identical isotropic EB arms and cannot safely supply these assemblies/forms.
No existing spectra, models, solver defaults or historical gates are changed.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import asdict
import csv
import hashlib
import io
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name] = '1'
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scipy.linalg import expm, block_diag
from scipy.optimize import brentq
from scripts.lib import inplane_spring_modes as mechanics
from scripts.lib import inplane_rotational_spring_eb_modes as modes
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_rlb as rlb
from scripts.analysis.laminated_beams import pilot_inplane_rotational_spring_eb as pilot

OUTPUT = ROOT/'results/laminated_beams/inplane_spring_robustness'
OLD = ROOT/'results/laminated_beams'
FS = float(mechanics.FREQUENCY_SCALE)
REF = mechanics.REFERENCE
BASE_ANGLES = (0.,5.,15.,30.,60.,80.)
KAPPAS = (0.,1.,100.)
CRITERIA = dict(frequency_interpretation=1e-6,sigma_ratio=1e-9,physical_residual=1e-9,
    compatibility=1e-10,null_residual=1e-9,rank_rtol=1e-12,
    angular_MAC=.95,angular_margin=.20,comparison_MAC=.90,comparison_margin=.20,
    matrix_atol=1e-12,H_rtol=1e-12,transfer_rtol=1e-9,quadrature_rtol=1e-6,
    guard_margin=.02,max_positions=10,max_point_matrices=6000,max_points=300,
    max_extra_points=60,max_event_recoveries=2,angle_localization_deg=.01)
PROTECTED_DIRS = ('inplane_rotational_spring_eb_beta','inplane_rotational_spring_eb_tracked',
    'inplane_rotational_spring_eb_tracked_k0','inplane_rotational_spring_eb_tracked_comparison',
    'inplane_rotational_spring_rlb_eb_limit')
ARRAYS = ('states','reactions')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    return {p.relative_to(ROOT).as_posix():sha(p) for name in PROTECTED_DIRS
            for p in sorted((OLD/name).glob('*')) if p.is_file()}


def atomic(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_text(text,encoding='utf-8')
    os.replace(tmp,path)


def write_csv(path,rows):
    if not rows:
        atomic(path,'')
        return
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
    writer.writeheader();writer.writerows(rows)
    atomic(path,stream.getvalue())


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream))


def pid(model,mu,kappa,beta):
    return f'{model}_m{mu:g}_k{kappa:g}_b{beta:.12g}'


class PointMatrices:
    def __init__(self,model,mu,kappa,beta,properties):
        self.pair = mechanics.arms(model,mu,properties)
        self.beta = math.radians(beta)
        self.joint = eb.Joint('SPRING',kappa*REF.D)
        self.cache = OrderedDict()
        self.full = self.blocks = self.expm = 0

    def tick(self):
        if self.full+self.blocks >= CRITERIA['max_point_matrices']:
            raise RuntimeError('COST_LIMIT')

    def assembly(self,omega):
        if omega in self.cache:
            return self.cache[omega]
        self.tick();self.full += 1
        self.expm += 1 if self.pair[0] == self.pair[1] else 2
        value = mechanics.assembly(omega,self.pair,self.beta,self.joint)
        self.cache[omega] = value
        if len(self.cache)>128:
            self.cache.popitem(last=False)
        return value

    def matrix(self,omega,parity=None):
        result = self.assembly(omega)
        if parity is None:
            return result.dimensionless
        if self.pair[0] != self.pair[1]:
            raise ValueError('no symmetry blocks for unequal arms')
        self.tick();self.blocks += 1
        return modes.class_matrix(result.endpoint_map,self.beta,self.joint.k_theta,REF,parity)


def norm_error(a,b,rtol):
    absolute = float(np.linalg.norm(np.asarray(a)-b))
    reference = float(np.linalg.norm(b))
    return dict(absolute=absolute,relative=absolute/max(reference,1e-300),
                accepted=absolute<=CRITERIA['matrix_atol']+rtol*reference)


def preflight():
    section,p = mechanics.section()
    base,p0 = mechanics.section(0.)
    symmetry = asdict(mechanics.laminate.check_laminate_symmetry(section))
    checks = []
    for name,expected in [('A',.011),('D',1.3*p0.D),('S',p0.S),('m',p0.m),('J',p0.J)]:
        checks.append(dict(kind='CONSTITUTIVE_'+name,**norm_error(getattr(p,name),expected,1e-12)))
    for name in ('axial_reduction','bending_reduction','shear_reduction_before_K'):
        item = getattr(p,name)
        checks.append(dict(kind='SCHUR_'+name,**norm_error(item.compliance_value,item.schur_value,1e-12)))
    f = modes.REFLECTION
    for model in ('EB','RLB'):
        for Omega in (2.,20.,80.):
            pair = mechanics.arms(model,0.,p)
            h = pair[0].matrix(Omega/FS)
            checks.append(dict(kind='FH_HF',model=model,Omega=Omega,**norm_error(f@h,h@f,1e-12)))
        for mu in (0.,.01):
            for beta in (0.,30.):
                pair = mechanics.arms(model,mu,p)
                context = PointMatrices(model,mu,1.,beta,p)
                actual = context.assembly(20/FS)
                checks.append(dict(kind='TWO_LENGTHS',model=model,mu=mu,beta=beta,
                    **norm_error(np.array([a.L for a in pair]),[1-mu,1+mu],1e-12)))
                for i,arm in enumerate(pair):
                    scale = arm.scale()
                    h = arm.matrix(20/FS)*scale[None,:]/scale[:,None]
                    endpoint = scale[:,None]*expm(h*arm.L)[:,3:]
                    checks.append(dict(kind='ENDPOINT_LENGTH',model=model,mu=mu,arm=i,
                        **norm_error(actual.endpoint_map[i*6:(i+1)*6,i*3:(i+1)*3],endpoint,1e-9)))
        # One arm-swap matrix control: reflection maps mu to -mu, not to itself.
        a = PointMatrices(model,.01,1.,30.,p).assembly(20/FS)
        b = PointMatrices(model,-.01,1.,30.,p).assembly(20/FS)
        swap_states = np.block([[np.zeros((6,6)),f],[f,np.zeros((6,6))]])
        r = f[3:,3:]
        swap_reactions = np.block([[np.zeros((3,3)),r],[r,np.zeros((3,3))]])
        checks.append(dict(kind='ARM_SWAP',model=model,**norm_error(
            b.endpoint_map@swap_reactions,swap_states@a.endpoint_map,1e-9)))
        checks.append(dict(kind='BROKEN_FIXED_LENGTH_SWAP',model=model,
            defect=float(np.linalg.norm(a.endpoint_map@swap_reactions-swap_states@a.endpoint_map))))
    for beta in (0.,30.):
        row = modes.symmetry_equivalence(math.radians(beta),REF.D,REF)
        checks.append(dict(kind='JOINT_CLASS_ROW_SPACE',beta=beta,**row))
        for Omega in (2.,20.,80.):
            native_eb = eb.EBArm(p.A,p.D,p.m,1.)
            limit = rlb.LimitArm(p,1.,0.)
            spring = eb.Joint('SPRING',REF.D)
            h_eb,h_rlb = eb.state_matrix(Omega/FS,native_eb),rlb.state_matrix(Omega/FS,limit)
            checks.append(dict(kind='ZERO_LIMIT_H',beta=beta,Omega=Omega,**norm_error(h_rlb,h_eb,1e-12)))
            a = eb.boundary_assembly(Omega/FS,native_eb,native_eb,math.radians(beta),spring,REF)
            b = rlb.boundary_assembly(Omega/FS,limit,limit,math.radians(beta),spring,REF)
            for field in ('endpoint_map','dimensionless'):
                checks.append(dict(kind='ZERO_LIMIT_'+field,beta=beta,Omega=Omega,
                    **norm_error(getattr(b,field),getattr(a,field),1e-9)))
    if not symmetry['is_symmetric'] or not all(c.get('accepted',True) for c in checks):
        raise RuntimeError('PREFLIGHT_MATRIX_MISMATCH')
    return dict(layup='H/L/L/H',contrast=.4,ply_angles_deg=[0]*4,ply_thickness=.0125,
        material_base=dict(E1=1.1,E2=.9,nu12=.3,G12=1/2.6,G13=1/2.6,G23=1/2.6,rho=1),
        properties=asdict(p),baseline_properties=asdict(p0),symmetry=symmetry,
        laminate={k:np.asarray(getattr(section,k)).tolist() for k in ('A','B','D','shear','I0','I1','I2')},
        checks=checks),p


def merged(windows):
    result = []
    for lo,hi in sorted(windows):
        if result and lo<=result[-1][1]:
            result[-1] = (result[-1][0],max(hi,result[-1][1]))
        else:result.append((lo,hi))
    return result


def intervals(predictions,upper):
    windows = merged([(max(.01,x-max(.6,.10*x)),min(upper,x+max(.6,.10*x)))
                      for x in predictions if x<upper])
    cursor = .01
    result = []
    for lo,hi in windows:
        if lo>cursor:
            result.append((cursor,lo,max(5,int((lo-cursor)/2)+1)))
        result.append((lo,hi,max(17,int((hi-lo)/.65)+1)))
        cursor = hi
    if cursor<upper:
        result.append((cursor,upper,max(5,int((upper-cursor)/2)+1)))
    return result


class Run:
    def __init__(self):
        OUTPUT.mkdir(parents=True,exist_ok=True)
        checkpoint = OUTPUT/'diagnostics.json'
        self.state = json.loads(checkpoint.read_text(encoding='utf-8')) if checkpoint.exists() else dict(
            points={},candidates={},mappings={},local={},events={},recovery_attempts={},extra_points=[],
            created_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            protected_sources=source_hashes(),criteria=CRITERIA,preflight=None,auxiliary=[],commands=[])
        if self.state['criteria'] != CRITERIA:
            raise ValueError('checkpoint criteria differ')
        self.shapes = dict(np.load(OUTPUT/'shapes.npz')) if (OUTPUT/'shapes.npz').exists() else {}
        self.old_rows = []
        self.old_shapes = {}
        seed_file=OLD/'inplane_rotational_spring_eb_tracked_comparison/seed_mode_mapping.csv'
        seed_rows=read_csv(seed_file)
        self.seed_map={(float(r['kappa']),r['source_branch_id']):r['reference_branch_id'] for r in seed_rows if r['status']=='CONFIRMED'}
        if len(self.seed_map)!=18:raise ValueError('published cross-kappa seed mapping incomplete')
        for name in ('inplane_rotational_spring_eb_tracked','inplane_rotational_spring_eb_tracked_k0'):
            manifest=json.loads((OLD/name/'run_manifest.json').read_text(encoding='utf-8'))
            if manifest['geometry']!={'l':1,'b':.2,'h':.05,'E':1,'rho':1}:
                raise ValueError('old shape geometry differs from declared seed reference')
            rows = read_csv(OLD/name/'tracked_branches.csv')
            if any(abs(float(r['Lambda'])**2-float(r['Omega']))>1e-12*float(r['Omega']) for r in rows):
                raise ValueError('old frequency normalization mismatch')
            self.old_rows.extend(rows)
            with np.load(OLD/name/'shapes.npz') as archive:
                self.old_shapes.update({r['shape_key']+'__states':archive[r['shape_key']+'__states']
                    for r in rows if float(r['beta_deg']) in BASE_ANGLES})
        self.properties = mechanics.section()[1]

    def save(self):
        for point in self.state['points'].values():
            for row in point['roots']:
                position=row.get('current_sorted_position') or row.get('candidate_ordinal',10)
                row['spectral_role']='ROOT' if position<point['requested_positions'] else 'GUARD'
        with (OUTPUT/'shapes.npz.tmp').open('wb') as stream:
            np.savez(stream,**self.shapes)
        os.replace(OUTPUT/'shapes.npz.tmp',OUTPUT/'shapes.npz')
        atomic(OUTPUT/'diagnostics.json',json.dumps(self.state,ensure_ascii=False,indent=2,allow_nan=False))
        roots = [{k:v for k,v in r.items() if k not in ('physical_residuals','failures','reconstruction_calls','transfer_check_errors')}
                 for p in self.state['points'].values() for r in p['roots']]
        write_csv(OUTPUT/'verified_roots.csv',roots)

    def form(self,row):
        states = self.shapes[row['shape_key']+'__states']
        pair = mechanics.arms(row['model'],row['mu'],self.properties)
        _,weights = modes.quadrature(states.shape[1])
        return dict(row,states=states,vector=mechanics.physical_vector(states,pair,weights),
                    common=mechanics.common_vector(states,weights))

    def old_forms(self,kappa,beta):
        rows = sorted((r for r in self.old_rows if float(r['kappa'])==kappa and float(r['beta_deg'])==beta),
                      key=lambda r:r['branch_id'])
        _,weights = modes.quadrature()
        return [dict(branch_id=self.seed_map[(kappa,r['branch_id'])],source_branch_id=r['branch_id'],symmetry_class=int(r['symmetry_class']),
                     common=mechanics.common_vector(self.old_shapes[r['shape_key']+'__states'],weights)) for r in rows]

    def predicted(self,model,mu,kappa,beta,n):
        available = [p for p in self.state['points'].values() if p['model']==model and p['kappa']==kappa
            and p['roots'] and len(p['roots'])>=n and p['status'] in ('CONFIRMED','TARGET_CONFIRMED_GUARD_QUALIFIED')]
        if available:
            chosen = min(available,key=lambda p:abs(p['beta_deg']-beta)+100*abs(p['mu']-mu))
            return [r['Omega'] for r in chosen['roots']][:n]
        rows = [r for r in self.old_rows if float(r['kappa'])==kappa]
        nearest = min({float(r['beta_deg']) for r in rows},key=lambda b:abs(b-beta))
        values = sorted(float(r['Omega'])*math.sqrt(self.properties.D/REF.D) for r in rows if float(r['beta_deg'])==nearest)
        # An axial predictor and the next bending estimate are only search hints.
        values += [math.pi/2*math.sqrt(self.properties.A/self.properties.m)*FS]
        while len(values)<n:
            values.append(max(values)*1.30)
        return sorted(values)[:n]

    def point(self,model,mu,kappa,beta,role='BASE',n=8,predictions=None):
        if model not in ('EB','RLB') or mu not in (0.,.005,.01) or kappa not in KAPPAS:
            raise ValueError('outside the fixed robustness formulations')
        if not np.isfinite(beta) or not 0<=beta<=90 or not 3<=n<=10:
            raise ValueError('angle or position budget outside the bounded check')
        if mu==.005 and role!='LEFT_MU_SEED_CONTINUATION':
            raise ValueError('mu=.005 reserved for an addressed seed ambiguity')
        key = pid(model,mu,kappa,beta)
        if key in self.state['points']:
            point = self.state['points'][key]
            self.repair_isolated(point)
            qualify_guard(point)
            return point
        if len(self.state['points'])+len(self.state['auxiliary'])>=CRITERIA['max_points']:
            raise RuntimeError('TOTAL_POINT_COST_LIMIT')
        if role not in ('BASE','A_WINDOW','B_WINDOW'):
            if len(self.state['extra_points'])>=CRITERIA['max_extra_points']:
                raise RuntimeError('EXTRA_POINT_COST_LIMIT')
            self.state['extra_points'].append(key)
        started = time.perf_counter()
        context = PointMatrices(model,mu,kappa,beta,self.properties)
        guesses = predictions or self.predicted(model,mu,kappa,beta,n)
        upper = max(guesses)*1.10+1
        records,events,suspects,errors = [],[],[],[]
        parities = (1,-1) if mu==0 else (None,)
        try:
            for parity in parities:
                provider = lambda omega:context.matrix(omega,parity)
                pool = []
                for lo,hi,count in intervals(guesses,upper):
                    pool.extend(pilot.roots._scan_candidates(provider,FS,pilot.policy(lo,hi),case_id=key,
                        builder_id=f'{model}_two_arm_{parity}',scan_id='BOUNDED_LOCAL',points=count,phases=(0.,))[0])
                pool,proof = pilot.reconcile_local_detections(pool,provider)
                accepted,ambiguous = pilot.consolidate(pool)
                if ambiguous:errors.append('UNRESOLVED_DETECTION_CLUSTER')
                suspects.extend(c for c in pool if pilot.suspicious(c))
                events.extend((c.omega_bar,parity,c.diagnostics.detected_nullity) for c in accepted)
                records.append(dict(parity=parity,candidates=[pilot.candidate_record(c) for c in pool],reconciliations=proof))
            events.sort()
            if len(events)<n:errors.append('MISSING_TARGET_OR_GUARD')
            if len(events)>10:errors.append('POSITION_LIMIT')
            events = events[:10]
            roots = []
            for position,(Omega,parity,multiplicity) in enumerate(events,1):
                # Multiplicity is independently detected nullity, never number of hits.
                if multiplicity>1 and parity is None:
                    errors.append('FULL_MULTIPLE_SUBSPACE_UNRESOLVED')
                row = dict(point_id=key,model=model,mu=mu,kappa=kappa,k_theta=kappa*REF.D,beta_deg=beta,
                    Omega=Omega,omega=Omega/FS,Lambda=math.sqrt(Omega),current_sorted_position=position,
                    multiplicity=multiplicity,grid_role=role,source='NEW_CHARACTERISTIC_ROOT',
                    shape_key=f'{key}_p{position:02d}',root_status='UNCONFIRMED',symmetry_class=parity)
                try:
                    recovered = mechanics.recover(context.assembly(Omega/FS),Omega/FS,context.pair,
                        context.beta,context.joint,parity)
                    for field in ARRAYS:
                        self.shapes[row['shape_key']+'__'+field] = recovered[field]
                    row.update({k:v for k,v in recovered.items() if k not in (*ARRAYS,'vector','common')})
                    row['root_status'] = 'CONFIRMED' if not row['failures'] else 'SHAPE_GATE_FAILED'
                    if recovered['detected_nullity']>1 and parity is None:
                        row.update(root_status='MULTIPLE_SUBSPACE_UNRESOLVED',s=None,
                            sensitivity_status='BASIS_DEPENDENT_NOT_A_SIMPLE_EIGENVALUE')
                except ValueError as error:
                    row['failures'] = [str(error)]
                roots.append(row)
            target_edge = events[min(n-2,len(events)-1)][0] if events else 0.
            below = [pilot.candidate_record(c) for c in suspects if c.interval_left_bar<=target_edge]
            if below:errors.append('UNRESOLVED_BELOW_TARGET')
            if any(r['root_status']!='CONFIRMED' for r in roots[:n-1]):errors.append('ROOT_OR_SHAPE_GATE')
            guard_ok = len(roots)>=n and roots[n-1]['root_status']=='CONFIRMED' and upper-roots[n-1]['Omega']>.02
            if not guard_ok:errors.append('GUARD_QUALIFIED')
        except (RuntimeError,ValueError) as error:
            roots = locals().get('roots',[])
            errors.append(str(error))
        status = 'CONFIRMED' if not errors else ('TARGET_CONFIRMED_GUARD_QUALIFIED' if set(errors)=={'GUARD_QUALIFIED'} else 'POINT_UNCONFIRMED')
        point = dict(point_id=key,model=model,mu=mu,kappa=kappa,beta_deg=beta,role=role,
            roots=roots,status=status,errors=errors,requested_positions=n,search_upper=upper,
            full_B=context.full,symmetry_B=context.blocks,transfer_expm=context.expm,
            seconds=time.perf_counter()-started,search=records)
        self.state['points'][key] = point
        self.repair_isolated(point)
        qualify_guard(point)
        self.save()
        print(key,point['status'],len(roots),'B',point['full_B']+point['symmetry_B'],flush=True)
        return point

    def repair_isolated(self,point):
        """One local full-B check for a rejected class candidate; no rescan.

        At beta=0 a class can contain a scalar axial factor whose vanishing
        row/column defeats adaptive equilibration. The complete B retains
        both axial reactions and is the primary physical check.
        """
        suspects = [(scan['parity'],c) for scan in point['search'] for c in scan['candidates']
                    if not c['accepted'] and c['reason']=='NULLITY_UNRESOLVED_AT_1E-12']
        if not suspects or point.get('isolated_repair'):
            return
        started = time.perf_counter()
        ctx = PointMatrices(point['model'],point['mu'],point['kappa'],point['beta_deg'],self.properties)
        ctx.full,ctx.blocks,ctx.expm = point['full_B'],point['symmetry_B'],point['transfer_expm']
        history = dict(prior_status=point['status'],prior_errors=point['errors'][:],candidates=[],reason='CLASS_ROOT_GATE_FULL_MATRIX_CHECK')
        point['isolated_repair'] = history
        try:
            for parity,c in suspects:
                center = c['Omega']
                provider = lambda omega:ctx.matrix(omega)
                pool = pilot.roots._scan_candidates(provider,FS,pilot.policy(center-.002,center+.002),
                    case_id=point['point_id'],builder_id='FULL_MATRIX',scan_id='ISOLATED_REPAIR',points=17,phases=(0.,))[0]
                pool,proof = pilot.reconcile_local_detections(pool,provider)
                accepted,ambiguous = pilot.consolidate(pool)
                history['candidates'].append(dict(original=c,parity=parity,
                    checked=[pilot.candidate_record(a) for a in pool],reconciliations=proof))
                if len(accepted)!=1 or ambiguous:continue
                event = accepted[0]
                if any(abs(r['Omega']-event.omega_bar)<1e-8 for r in point['roots']):continue
                omega = event.omega_bar/FS
                recovered = mechanics.recover(ctx.assembly(omega),omega,ctx.pair,ctx.beta,ctx.joint,parity)
                shape_key = point['point_id']+f'_repair{len(history["candidates"])}'
                row = dict(point_id=point['point_id'],model=point['model'],mu=point['mu'],kappa=point['kappa'],
                    k_theta=point['kappa']*REF.D,beta_deg=point['beta_deg'],Omega=event.omega_bar,omega=omega,
                    Lambda=math.sqrt(event.omega_bar),multiplicity=event.diagnostics.detected_nullity,
                    grid_role=point['role'],source='LOCAL_FULL_MATRIX_REFINEMENT',previous_rejected_Omega=center,
                    shape_key=shape_key,root_status='CONFIRMED' if not recovered['failures'] else 'SHAPE_GATE_FAILED')
                for field in ARRAYS:self.shapes[shape_key+'__'+field] = recovered[field]
                row.update({k:v for k,v in recovered.items() if k not in (*ARRAYS,'vector','common')})
                point['roots'].append(row)
        except (RuntimeError,ValueError) as error:
            history['error'] = str(error)
        point['roots'].sort(key=lambda r:r['Omega'])
        for index,row in enumerate(point['roots'],1):row['current_sorted_position']=index
        n = point['requested_positions']
        if len(point['roots'])>=n and all(r['root_status']=='CONFIRMED' for r in point['roots'][:n]):
            rejected = [(scan['parity'],c) for scan in point['search'] for c in scan['candidates']
                        if not c['accepted'] and c['reason']!='FALSE_SIGMA_VALLEY']
            resolved = all(any(abs(r['Omega']-c['Omega'])<.002 and r['symmetry_class']==eta
                          and r['root_status']=='CONFIRMED' for r in point['roots']) for eta,c in rejected)
            if resolved and not any(e in point['errors'] for e in ('POSITION_LIMIT','UNRESOLVED_DETECTION_CLUSTER')):
                point['errors']=[];point['status']='CONFIRMED'
        point.update(full_B=ctx.full,symmetry_B=ctx.blocks,transfer_expm=ctx.expm)
        history['seconds']=time.perf_counter()-started
        point['seconds']+=history['seconds']
        self.save()
        print('isolated full B',point['point_id'],point['status'],len(point['roots']),flush=True)

    def map_point(self,point,reference,tag,physical=False,restrict=False):
        candidates = [self.form(r) for r in point['roots'] if r['root_status']=='CONFIRMED']
        if len(candidates)<len(reference):return []
        indices,mac,margin = mechanics.match(reference,candidates,physical=physical,restrict_symmetry=restrict)
        threshold = CRITERIA['angular_MAC'] if physical else CRITERIA['comparison_MAC']
        mapped = []
        for i,j in enumerate(indices):
            row = dict(candidates[j])
            row.update(branch_id=reference[i]['branch_id'],MAC=float(mac[i]),margin=float(margin[i]),
                local_assignment_status='CONFIRMED' if mac[i]>=threshold and margin[i]>=.20 else 'MAPPING_AMBIGUOUS',
                matching_method=tag,mapping_status='CONFIRMED' if mac[i]>=threshold and margin[i]>=.20
                    and reference[i].get('mapping_status','CONFIRMED')=='CONFIRMED' else 'MAPPING_AMBIGUOUS')
            mapped.append(row)
        self.state['mappings'][point['point_id']+'_'+tag] = [
            {k:v for k,v in row.items() if k not in (*ARRAYS,'common','vector')} for row in mapped]
        return mapped

    def base(self):
        for beta in BASE_ANGLES:
            for kappa in KAPPAS:
                expected=[('EB',0.,'ISOTROPIC_FORM_COMPARISON'),('RLB',0.,'EB_RLB_SAME_BETA_COMPARISON'),
                          ('EB',.01,'LENGTH_PERTURBATION_SAME_BETA'),('RLB',.01,'LENGTH_PERTURBATION_SAME_BETA')]
                if all(pid(model,mu,kappa,beta)+'_'+tag in self.state['mappings'] for model,mu,tag in expected):
                    continue
                original = self.old_forms(kappa,beta)
                a = self.point('EB',0.,kappa,beta)
                first = self.map_point(a,original,'ISOTROPIC_FORM_COMPARISON',restrict=True)
                b = self.point('RLB',0.,kappa,beta)
                second = self.map_point(b,first,'EB_RLB_SAME_BETA_COMPARISON',restrict=True) if first else []
                for model,reference in [('EB',first),('RLB',second)]:
                    c = self.point(model,.01,kappa,beta)
                    if reference:self.map_point(c,reference,'LENGTH_PERTURBATION_SAME_BETA')
                self.save()


def qualify_guard(point):
    """Keep strict failures visible while separating an isolated upper event."""
    if point['status']=='CONFIRMED' or point.get('guard_qualification'):
        return
    if set(point['errors'])-{'MISSING_TARGET_OR_GUARD','UNRESOLVED_BELOW_TARGET','ROOT_OR_SHAPE_GATE',
                            'GUARD_QUALIFIED','UNRESOLVED_DETECTION_CLUSTER'}:return
    n=point['requested_positions'];roots=point['roots']
    if len(roots)<n or any(r['root_status']!='CONFIRMED' for r in roots[:n-1]):return
    edge=roots[n-2]['Omega']
    if roots[n-1]['Omega']-edge<=.02:return
    suspect_below=[]
    for scan in point['search']:
        for c in scan['candidates']:
            if c['accepted'] or c['reason']=='FALSE_SIGMA_VALLEY':continue
            repaired=any(r.get('source')=='LOCAL_FULL_MATRIX_REFINEMENT' and
                r['symmetry_class']==scan['parity'] and abs(r['Omega']-c['Omega'])<.002 and
                r['root_status']=='CONFIRMED' for r in roots)
            if not repaired and c['interval'][0]<=edge:suspect_below.append(c)
    if suspect_below:return
    # An unresolved duplicate detection is not promoted to an extra eigenmode.
    if 'UNRESOLVED_DETECTION_CLUSTER' in point['errors']:
        for i in range(1,n-1):
            if roots[i]['Omega']-roots[i-1]['Omega']<5e-10+5e-12*roots[i]['Omega']:
                return
        for r in roots[n-1:]:
            r['original_root_status']=r['root_status']
            r['root_status']='GUARD_DETECTION_AMBIGUOUS'
            r['candidate_ordinal']=r['current_sorted_position']
            r['current_sorted_position']=None
    point['guard_qualification']=dict(original_status=point['status'],original_errors=point['errors'][:],
        confirmed_lower_positions=n-1,target_edge=edge,guard_lower=roots[n-1]['Omega'],
        basis='all lower physical roots accepted; no unrepaired detected event below target; guard separated',
        guard_status='QUALIFIED_ORIGINAL_FAILURE_RETAINED')
    point['status']='TARGET_CONFIRMED_GUARD_QUALIFIED'


def pair_rows(run,point):
    """The first opposite-class pair, selected only in a symmetric problem."""
    if point['mu']!=0:raise ValueError('class selection forbidden at nonzero mu')
    result=[]
    for branch,eta in [('mode_01',1),('mode_02',-1)]:
        eligible=[r for r in point['roots'] if r['symmetry_class']==eta and r['root_status']=='CONFIRMED']
        if not eligible:return []
        row=run.form(min(eligible,key=lambda r:r['Omega']))
        row['branch_id']=branch
        result.append(row)
    return result


def clean_rows(rows):
    return [{k:v for k,v in r.items() if k not in (*ARRAYS,'vector','common')} for r in rows]


def insensitivity_candidate(run,model):
    key='A_'+model
    if key in run.state['candidates']:return run.state['candidates'][key]
    started=time.perf_counter();arm=mechanics.arms(model,0.,run.properties)[0]
    calls=0
    def endpoint(Omega):
        nonlocal calls
        calls+=1
        scale=arm.scale();h=arm.matrix(Omega/FS)*scale[None,:]/scale[:,None]
        return scale[:,None]*expm(h*arm.L)[:,3:]
    def determinant(Omega):
        value=endpoint(Omega)[np.ix_([2,5],[1,2])]
        return float(np.linalg.det(value/np.array([1.,REF.D])[:,None]))
    # First isolated clamp-to-end psi=M=0 root; finite user-authorized auxiliary problem.
    Omega=brentq(determinant,8.,15.,xtol=1e-11)
    p=endpoint(Omega)
    _,_,vh=np.linalg.svd(p[np.ix_([2,5],[1,2])]/np.array([1.,REF.D])[:,None])
    bending=p[:,1:]@vh[-1]
    compliance=p[0,0]/p[3,0]
    tangent_squared=-compliance*bending[4]/bending[1]
    if tangent_squared<=0:raise ValueError('NO_PHYSICAL_CANDIDATE_ANGLE')
    beta=math.degrees(2*math.atan(math.sqrt(tangent_squared)))
    record=dict(model=model,beta_deg=beta,Omega=Omega,Lambda=math.sqrt(Omega),
        candidate_method='PSI_M_ZERO_SINGLE_ARM_THEN_TRANSLATION_FORCE',
        matrices=calls,expm=calls,seconds=time.perf_counter()-started,
        bending_psi=float(bending[2]),bending_M=float(bending[5]),tangent_squared=tangent_squared)
    run.state['candidates'][key]=record
    run.state['auxiliary'].append(record)
    run.save();print('candidate A',model,beta,Omega,flush=True)
    return record


def crossing_candidate(run,model):
    key='B_'+model
    if key in run.state['candidates']:return run.state['candidates'][key]
    snapshots=[]
    def difference(beta):
        point=run.point(model,0.,1.,beta,'B_LOCALIZATION',n=3)
        selected=pair_rows(run,point)
        if len(selected)!=2:raise RuntimeError('CROSSING_PAIR_UNRESOLVED')
        value=selected[0]['Omega']-selected[1]['Omega']
        snapshots.append(dict(beta_deg=beta,difference=value,roots=clean_rows(selected)))
        return value
    sampled=[(beta,difference(beta)) for beta in BASE_ANGLES]
    brackets=[(left,right) for left,right in zip(sampled,sampled[1:]) if left[1]*right[1]<0]
    if not brackets:raise RuntimeError('CROSSING_NOT_BRACKETED_IN_CONTROL_SAMPLE')
    (lo,a),(hi,b)=brackets[0]
    while hi-lo>.01:
        mid=(lo+hi)/2;c=difference(mid)
        if a*c<=0:hi,b=mid,c
        else:lo,a=mid,c
    result=dict(model=model,beta_deg=(lo+hi)/2,bracket=[lo,hi],signs=[a,b],
        classification='CROSSING_SUPPORTED',basis='independent reflection classes; resolved sign reversal',
        snapshots=snapshots)
    run.state['candidates'][key]=result;run.save()
    print('candidate B',model,result['bracket'],flush=True)
    return result


def track_window(run,model,mu,angles):
    previous=None;rows=[];flagged=[]
    for beta in sorted(angles):
        point=run.state['points'][pid(model,mu,1.,beta)]
        if mu==0:
            current=pair_rows(run,point)
            if previous and current:
                indices,mac,margin=mechanics.match(previous,current,physical=True,restrict_symmetry=True)
                for i,j in enumerate(indices):
                    current[j].update(MAC=float(mac[i]),margin=float(margin[i]),mapping_status='CONFIRMED')
            else:
                for row in current:row.update(MAC=1.,margin=1.,mapping_status='CONFIRMED')
        else:
            if previous is None:
                symmetric=run.state['points'][pid(model,0.,1.,beta)]
                reference=pair_rows(run,symmetric)
                current=run.map_point(point,reference,'B_LEFT_REFERENCE') if reference else []
                if current and any(r['mapping_status']!='CONFIRMED' for r in current):
                    # Direct comparison at the left edge is already mixed.
                    # The task permits one addressed mu=.005 check when .01
                    # cannot be interpreted without the perturbation scale.
                    bridge=run.point(model,.005,1.,beta,'LEFT_MU_SEED_CONTINUATION',n=3)
                    middle=run.map_point(bridge,reference,'B_LEFT_MU_0_TO_005')
                    if middle and all(r['mapping_status']=='CONFIRMED' for r in middle):
                        current=run.map_point(point,middle,'B_LEFT_MU_005_TO_01')
                    run.state['recovery_attempts']['B_'+model]=max(1,run.state['recovery_attempts'].get('B_'+model,0))
            else:
                current=run.map_point(point,previous,'B_SEQUENTIAL_PHYSICAL_MAC',physical=True)
        if len(current)!=2:
            flagged.append(dict(left=previous[0]['beta_deg'] if previous else beta,right=beta,reason='MISSING_PAIR'))
            rows.extend(dict(model=model,mu=mu,kappa=1.,beta_deg=beta,branch_id=branch,Lambda=None,
                Omega=None,root_status='UNCONFIRMED',mapping_status='MAPPING_AMBIGUOUS')
                for branch in ('mode_01','mode_02'))
            continue
        if any(r.get('local_assignment_status',r.get('mapping_status'))!='CONFIRMED' for r in current):
            flagged.append(dict(left=previous[0]['beta_deg'] if previous else beta,right=beta,reason='LOW_MAC_OR_MARGIN'))
        # An ambiguous numerical assignment is retained diagnostically, never plotted as confirmed.
        rows.extend(clean_rows(current));previous=current
    return rows,flagged


def local_checks(run):
    for model in ('EB','RLB'):
        prior=run.state['local'].get('B_'+model)
        if prior and (prior['refinement_rounds']==2 or not prior['flags']['0.01']):
            continue
        candidate=insensitivity_candidate(run,model)
        angles=[candidate['beta_deg']+shift for shift in (-2.,-1.,0.,1.,2.)]
        rows=[]
        for beta in angles:
            symmetric={}
            for kappa in KAPPAS:
                p=run.point(model,0.,kappa,beta,'A_WINDOW',n=3)
                chosen=pair_rows(run,p)
                if chosen:
                    first=chosen[0];first.update(MAC=1.,margin=1.,mapping_status='CONFIRMED',matching_method='FIRST_PLUS_CLASS')
                    symmetric[kappa]=first;rows.extend(clean_rows([first]))
                q=run.point(model,.01,kappa,beta,'A_WINDOW',n=3)
                if chosen:
                    assigned=run.map_point(q,[first],'A_SAME_BETA_COMMON_METRIC')
                    rows.extend(clean_rows(assigned))
        shifted=run.state['candidates'].get('A_SHIFT_'+model,{})
        rows.extend(shifted.get('rows',[]))
        run.state['local']['A_'+model]=dict(angles=angles,rows=rows)
        crossing=crossing_candidate(run,model)
        angles=[crossing['beta_deg']+shift for shift in range(-5,6)]
        for beta in angles:
            for mu in (0.,.01):run.point(model,mu,1.,beta,'B_WINDOW',n=3)
        for attempt in range(3):
            tracks={};flags={}
            for mu in (0.,.01):tracks[str(mu)],flags[str(mu)]=track_window(run,model,mu,angles)
            unresolved=flags['0.01']
            run.state['local']['B_'+model]=dict(angles=sorted(angles),tracks=tracks,flags=flags,refinement_rounds=attempt)
            run.save()
            if not unresolved or attempt==2:break
            additions=sorted({(r['left']+r['right'])/2 for r in unresolved if r['right']>r['left']})
            if not additions:break
            run.state['recovery_attempts']['B_'+model]=attempt+1
            for beta in additions:
                for mu in (0.,.01):run.point(model,mu,1.,beta,'B_MAC_REFINEMENT',n=3)
            angles.extend(additions)
        # One gap-resolution check: the existing bracket's midpoint is already a true solve.
        # No new optimization or arbitrary-precision search is introduced.
        run.save()


def refine_insensitivity(run):
    """One addressed zero-of-relative-rotation localization per theory.

    Phase is tied to a fixed computed form. Each beta probe solves the full
    unequal-arm characteristic problem; no interpolated eigenvalue is stored.
    """
    for model in ('EB','RLB'):
        key='A_SHIFT_'+model
        if key in run.state['candidates']:continue
        center=run.state['candidates']['A_'+model]['beta_deg']
        local=run.state['local']['A_'+model]
        reference_rows=[r for r in local['rows'] if r['mu']==.01 and r['kappa']==1. and r['beta_deg']==center]
        if not reference_rows:continue
        reference=run.form(reference_rows[0]);reference['branch_id']='mode_01'
        history=[]
        def delta(beta):
            point=run.point(model,.01,1.,beta,'A_SHIFT_LOCALIZATION',n=3)
            mapped=run.map_point(point,[reference],'A_FIXED_PHASE_COMPARISON')
            if len(mapped)!=1 or mapped[0]['mapping_status']!='CONFIRMED':
                raise RuntimeError('A_SHIFT_MATCH_UNRESOLVED')
            row=mapped[0]
            phase=1. if np.vdot(reference['common'],row['common']).real>=0 else -1.
            value=phase*row['Delta_psi']
            history.append(dict(beta_deg=beta,phase_aligned_Delta_psi=value,root=clean_rows([row])[0]))
            return value
        try:
            left,right=center-1,center+1
            dleft,dright=delta(left),delta(right)
            if dleft*dright>=0:
                result=dict(status='WEAK_REGION_ONLY_NO_ZERO_BRACKET',probes=history)
            else:
                beta=brentq(delta,left,right,xtol=.005)
                snapshots=[]
                for kappa in KAPPAS:
                    point=run.point(model,.01,kappa,beta,'A_SHIFT_VERIFICATION',n=3)
                    mapped=run.map_point(point,[reference],'A_SHIFT_THREE_SPRINGS')
                    snapshots.extend(clean_rows(mapped))
                ok=len(snapshots)==3 and all(r['mapping_status']=='CONFIRMED' for r in snapshots)
                R=(max(r['Omega'] for r in snapshots)-min(r['Omega'] for r in snapshots))/next(r['Omega'] for r in snapshots if r['kappa']==1.) if ok else None
                result=dict(status='SHIFTED_WEAK_SENSITIVITY_POINT' if ok else 'MAPPING_AMBIGUOUS',
                    beta_deg=beta,shift_deg=beta-center,R=R,rows=snapshots,probes=history,
                    angular_solver_xtol=.005,qualification='finite localization, not proof of exact equality at mu != 0')
                local['rows'].extend(snapshots)
            run.state['candidates'][key]=result
        except (RuntimeError,ValueError) as error:
            run.state['candidates'][key]=dict(status='UNRESOLVED',reason=str(error),probes=history)
        run.save()


def summarize(run):
    state=run.state
    rows=[]
    for key,mapping in state['mappings'].items():
        if mapping and mapping[0]['grid_role']=='BASE':rows.extend(mapping)
    write_csv(OUTPUT/'control_modes.csv',rows)
    write_csv(OUTPUT/'seed_mapping.csv',[r for r in rows if r['beta_deg']==0])
    metrics=[]
    for model in ('EB','RLB'):
        for mu in (0.,.01):
            for beta in BASE_ANGLES:
                for branch in [f'mode_{i:02d}' for i in range(1,7)]:
                    group=[r for r in rows if r['model']==model and r['mu']==mu and r['beta_deg']==beta and r['branch_id']==branch]
                    if len(group)!=3:continue
                    ok=all(r['mapping_status']=='CONFIRMED' for r in group)
                    reference=next(r for r in group if r['kappa']==1)
                    metrics.append(dict(model=model,mu=mu,beta_deg=beta,branch_id=branch,
                        R=(max(r['Omega'] for r in group)-min(r['Omega'] for r in group))/reference['Omega'] if ok else None,
                        s_k0=next(r['s'] for r in group if r['kappa']==0),s_k1=reference['s'],
                        s_k100=next(r['s'] for r in group if r['kappa']==100),
                        status='CONFIRMED_COMPARISON' if ok else 'MAPPING_AMBIGUOUS',
                        min_MAC=min(r['MAC'] for r in group)))
    write_csv(OUTPUT/'control_metrics.csv',metrics)
    # Independent of historical bending labels: compare the resolved eta=-1
    # family at the SAME beta by physical mass MAC, with no sorted-index claim.
    family=[]
    for model in ('EB','RLB'):
        for beta in BASE_ANGLES:
            groups={k:[run.form(r) for r in state['points'][pid(model,0.,k,beta)]['roots']
                       if r['symmetry_class']==-1 and r['root_status']=='CONFIRMED'] for k in KAPPAS}
            reference=groups[1.]
            if not reference or any(len(v)<len(reference) for v in groups.values()):continue
            aligned={1.:reference};qualities=[]
            for k in (0.,100.):
                indices,mac,margin=mechanics.match(reference,groups[k],physical=True,restrict_symmetry=True)
                aligned[k]=[groups[k][i] for i in indices]
                qualities.extend(zip(mac,margin))
            ok=all(a>=.95 and b>=.20 for a,b in qualities)
            for i,ref in enumerate(reference):
                values=[aligned[k][i] for k in KAPPAS]
                family.append(dict(model=model,beta_deg=beta,reference_class_position=i+1,
                    reference_sorted_position=ref['current_sorted_position'],
                    R=(max(r['Omega'] for r in values)-min(r['Omega'] for r in values))/ref['Omega'] if ok else None,
                    max_s=max(r['s'] for r in values),max_abs_Delta_psi=max(abs(r['Delta_psi']) for r in values),
                    min_MAC=min(a for a,b in qualities),status='CONFIRMED_FAMILY_COMPARISON' if ok else 'AMBIGUOUS'))
    write_csv(OUTPUT/'symmetric_family.csv',family)
    local=[];events=[]
    for model in ('EB','RLB'):
        if 'A_'+model in state['local']:
            local.extend([dict(r,window='A') for r in state['local']['A_'+model]['rows']])
        if 'B_'+model in state['local']:
            data=state['local']['B_'+model]
            local.extend([dict(r,window='B') for values in data['tracks'].values() for r in values])
            selected=data['tracks']['0.01']
            gaps=[];subspaces=[]
            for beta in data['angles']:
                pair=[r for r in selected if r['beta_deg']==beta]
                if len(pair)==2:
                    gaps.append(dict(beta_deg=beta,g=2*abs(pair[0]['Omega']-pair[1]['Omega'])/(pair[0]['Omega']+pair[1]['Omega']),
                        Omega=[r['Omega'] for r in pair],weights_plus=[r['reflection_weight_plus'] for r in pair],
                        MAC=[r['MAC'] for r in pair],status=[r['mapping_status'] for r in pair]))
                    symmetric=[r for r in data['tracks']['0.0'] if r['beta_deg']==beta]
                    if len(symmetric)==2:
                        left=[run.form(r) for r in symmetric];right=[run.form(r) for r in pair]
                        # QR bases describe spans only. They are NOT saved or
                        # represented as modes at either distinct eigenvalue.
                        qa=np.linalg.qr(np.array([r['common'] for r in left]).T)[0]
                        qb=np.linalg.qr(np.array([r['common'] for r in right]).T)[0]
                        cosines=np.linalg.svd(qa.conj().T@qb,compute_uv=False)
                        physical=np.array([r['vector'] for r in right])
                        subspaces.append(dict(beta_deg=beta,principal_cosines=cosines.tolist(),
                            physical_mass_orthogonality_error=float(np.linalg.norm(physical.conj()@physical.T-np.eye(2))),
                            note='span comparison only; stored eigenforms and frequencies unchanged'))
            minimum=min(gaps,key=lambda r:r['g']) if gaps else None
            complete=not data['flags']['0.01'] and len(selected)==2*len(data['angles'])
            start=[r for r in selected if r['beta_deg']==min(data['angles'])]
            end=[r for r in selected if r['beta_deg']==max(data['angles'])]
            exchanged=(len(start)==len(end)==2 and all(
                (a['reflection_weight_plus']-.5)*(b['reflection_weight_plus']-.5)<0 for a,b in zip(start,end)))
            status='AVOIDED_CROSSING_RESOLVED' if complete and minimum and minimum['g']>20e-6 and exchanged else 'UNRESOLVED_CLOSE_CLUSTER'
            event=dict(model=model,mu=.01,kappa=1.,classification=status,minimum_observed=minimum,
                character_exchange=exchanged,tracking_complete=complete,gaps=gaps,pair_subspaces=subspaces,
                symmetric_candidate=state['candidates']['B_'+model],
                limitation='finite window/grid; frequency interpretation 1e-6 is not a rigorous error bound')
            state['events']['B_'+model]=event;events.append(event)
    write_csv(OUTPUT/'local_modes.csv',local)
    a_metrics=[]
    for model in ('EB','RLB'):
        for mu in (0.,.01):
            selected=[r for r in local if r['window']=='A' and r['model']==model and r['mu']==mu]
            for beta in sorted({r['beta_deg'] for r in selected}):
                group=[r for r in selected if r['beta_deg']==beta]
                if len(group)!=3:continue
                valid=all(r['mapping_status']=='CONFIRMED' for r in group)
                reference=next(r for r in group if r['kappa']==1.)
                a_metrics.append(dict(model=model,mu=mu,beta_deg=beta,
                    R=(max(r['Omega'] for r in group)-min(r['Omega'] for r in group))/reference['Omega'] if valid else None,
                    max_s=max(r['s'] for r in group),status='CONFIRMED_COMPARISON' if valid else 'AMBIGUOUS'))
    write_csv(OUTPUT/'insensitivity_metrics.csv',a_metrics)
    atomic(OUTPUT/'local_events.json',json.dumps(events,indent=2,allow_nan=False))
    # Source hashes are verified on every summary; old manifests are never changed.
    now=source_hashes()
    protected=now==state['protected_sources']
    if not protected:raise RuntimeError('PROTECTED_SOURCE_CHANGED')
    points=list(state['points'].values());roots=[r for p in points for r in p['roots']]
    manifest=dict(policy='frequency-map-v1',mode='fast_plot',spectrum_semantics='tracked_branches',
        semantics_by_subset=dict(BASE='same_beta_form_comparisons',B_WINDOW='tracked_branches',A_WINDOW='local_form_comparison'),
        scope='72 sparse controls; first insensitivity and first crossing windows only',
        created_HEAD=state['created_HEAD'],working_tree_sources={name:sha(ROOT/name) for name in (
            'scripts/lib/inplane_spring_modes.py','scripts/analysis/laminated_beams/check_inplane_spring_robustness.py')},
        environment=dict(executable=sys.executable,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                         matplotlib=version('matplotlib'),pytest=version('pytest')),
        parameters=dict(l_ref=1.,b=.20,h=.05,mu=[0.,.01],kappa=KAPPAS,base_angles=BASE_ANGLES,
                        D_ref=REF.D,m_ref=REF.m,frequency_scale=FS),criteria=CRITERIA,
        state_order=['u','w','psi','N','Q','M'],code_status='working-tree version',
        preflight=state['preflight'],protected_sources=state['protected_sources'],protected_sources_unchanged=protected,
        quadrature_checks=state.get('quadrature_checks'),
        preflight_cost=dict(full_B=24,transfer_expm=48,source='explicit counts of the recorded preflight loops'),
        execution_notes=state.get('execution_notes',[]),
        verification=state.get('verification'),
        foundation_sha256={name:sha(ROOT/name) for name in (
            'scripts/lib/reddy_symmetric_laminated_beam.py','scripts/lib/inplane_rotational_spring_eb.py',
            'scripts/lib/inplane_rotational_spring_rlb.py','scripts/lib/inplane_rotational_spring_eb_modes.py',
            'scripts/analysis/laminated_beams/pilot_inplane_rotational_spring_eb.py',
            'scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py',
            'docs/laminated_beams/inplane_rotational_spring_joint.md',
            'docs/laminated_beams/reddy_stiffness_layout_contrast_sweep.md',
            'docs/laminated_beams/reddy_ch4_source_contract.md',
            'docs/numerics/frequency_map_computation_policy.md')},
        points=len(points),base_points=sum(p['role']=='BASE' for p in points),
        confirmed_points=sum(p['status']=='CONFIRMED' for p in points),
        qualified_points=sum(p['status']=='TARGET_CONFIRMED_GUARD_QUALIFIED' for p in points),
        unconfirmed_points=[p['point_id'] for p in points if p['status']=='POINT_UNCONFIRMED'],
        reused_spectra=0,roots=len(roots),confirmed_roots=sum(r['root_status']=='CONFIRMED' for r in roots),
        forms=sum('M' in r for r in roots),
        full_B=sum(p['full_B'] for p in points),symmetry_B=sum(p['symmetry_B'] for p in points),
        transfer_expm=sum(p['transfer_expm'] for p in points),
        reconstruction_calls={kind:sum(r.get('reconstruction_calls',{}).get(kind,0) for r in roots) for kind in ('analytic','expm')},
        auxiliary=state['auxiliary'],extra_points=state['extra_points'],
        isolated_recoveries=sum(bool(p.get('isolated_repair')) for p in points),event_recoveries=state['recovery_attempts'],
        compute_seconds=sum(p['seconds'] for p in points)+sum(p['seconds'] for p in state['auxiliary']),
        commands=state['commands'],limitations=['one H/L/L/H contrast .4; mu .01 only',
            'sparse controls use same-beta form comparison, not a full angular continuation',
            'no independent FEM/Ritz or arbitrary-laminate verification; historical source qualifications remain'])
    atomic(OUTPUT/'run_manifest.json',json.dumps(manifest,indent=2,allow_nan=False))
    run.save()
    print('summary',manifest['points'],manifest['confirmed_points'],manifest['unconfirmed_points'],flush=True)


def quadrature_checks(run):
    if 'quadrature_checks' in run.state:return
    started=time.perf_counter();records=[]
    for model in ('EB','RLB'):
        for mu,beta,positions in [(0.,0.,(1,7)),(.01,30.,(1,2))]:
            point=run.state['points'][pid(model,mu,1.,beta)]
            for position in positions:
                root=next(r for r in point['roots'] if r['current_sorted_position']==position)
                if root['root_status']!='CONFIRMED':continue
                ctx=PointMatrices(model,mu,1.,beta,run.properties)
                parity=root['symmetry_class'] if mu==0 else None
                assembled=ctx.assembly(root['omega'])
                coarse=mechanics.recover(assembled,root['omega'],ctx.pair,ctx.beta,ctx.joint,parity,129,check=True)
                fine=mechanics.recover(assembled,root['omega'],ctx.pair,ctx.beta,ctx.joint,parity,257,check=True)
                _,weights=modes.quadrature(129)
                restricted=mechanics.physical_vector(fine['states'][:,::2,:],ctx.pair,weights)
                mass_error=abs(fine['mass_before_normalization']/coarse['mass_before_normalization']-1)
                mac=modes.mac_matrix([coarse['vector']],[restricted])[0,0]
                records.append(dict(point_id=point['point_id'],position=position,
                    relative_mass_change=mass_error,MAC_change=abs(1-float(mac)),
                    direct_transfer_error=max(coarse['transfer_check_errors']+fine['transfer_check_errors']),
                    accepted=bool(mass_error<=1e-6 and abs(1-mac)<=1e-6),
                    full_B=ctx.full,transfer_expm=ctx.expm,shape_reconstructions=2,
                    reconstruction_calls={k:coarse['reconstruction_calls'][k]+fine['reconstruction_calls'][k] for k in ('analytic','expm')}))
    run.state['quadrature_checks']=dict(records=records,seconds=time.perf_counter()-started)
    run.save()


def curve(rows,model,mu,kappa,branch):
    selected=sorted((r for r in rows if r['model']==model and float(r['mu'])==mu and
        float(r['kappa'])==kappa and r['branch_id']==branch),key=lambda r:float(r['beta_deg']))
    return ([float(r['beta_deg']) for r in selected],
            [float(r['Lambda']) if r['mapping_status']=='CONFIRMED' and r['root_status']=='CONFIRMED' else np.nan for r in selected])


def render():
    """Strict plot-only: saved local rows only; no Run or computational helper calls."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    started=time.perf_counter()
    rows=read_csv(OUTPUT/'local_modes.csv')
    a_metrics=read_csv(OUTPUT/'insensitivity_metrics.csv') if (OUTPUT/'insensitivity_metrics.csv').exists() else []
    for model in ('EB','RLB'):
        for window in ('A','B'):
            data=[r for r in rows if r['window']==window]
            fig,(ax,diagnostic)=plt.subplots(2,1,figsize=(7.4,6.4),sharex=True,
                gridspec_kw={'height_ratios':[2.1,1]},layout='constrained')
            if window=='B':
                for i,branch in enumerate(('mode_01','mode_02')):
                    for mu,style in [(0.,'--'),(.01,'-')]:
                        x,y=curve(data,model,mu,1.,branch)
                        ax.plot(x,y,style,color=('tab:blue','tab:orange')[i],marker='.' if mu else None,
                            label=f'{branch}, μ={mu:g}',lw=1.6)
                ax.set_title(model+': первая пара, κθ=1')
                for i,branch in enumerate(('mode_01','mode_02')):
                    selected=sorted((r for r in data if r['model']==model and float(r['mu'])==.01
                        and r['branch_id']==branch),key=lambda r:float(r['beta_deg']))
                    if selected and 'reflection_weight_plus' in selected[0]:
                        diagnostic.plot([float(r['beta_deg']) for r in selected],
                            [float(r['reflection_weight_plus']) if r['mapping_status']=='CONFIRMED' else np.nan for r in selected],
                            '.-',color=('tab:blue','tab:orange')[i])
                diagnostic.set_ylabel('Вес класса +\n(reference-метрика)');diagnostic.set_ylim(-.03,1.03)
            else:
                for kappa,color in [(0.,'tab:blue'),(1.,'tab:orange'),(100.,'tab:green')]:
                    for mu,style in [(0.,'--'),(.01,'-')]:
                        x,y=curve(data,model,mu,kappa,'mode_01')
                        ax.plot(x,y,style,color=color,marker='.' if mu else None,label=f'κθ={kappa:g}, μ={mu:g}',lw=1.6)
                ax.set_title(model+': первая область нечувствительности')
                for mu,color,style in [(0.,'black','--'),(.01,'tab:purple','-')]:
                    values=sorted((r for r in a_metrics if r['model']==model and float(r['mu'])==mu),key=lambda r:float(r['beta_deg']))
                    if values:diagnostic.semilogy([float(r['beta_deg']) for r in values],
                        [float(r['R']) if r['R'] and float(r['R'])>0 else np.nan for r in values],
                        style,marker='.',color=color,label=f'μ={mu:g}')
                diagnostic.axhline(1e-6,color='grey',lw=.8,ls=':',label='ориентир 10⁻⁶')
                diagnostic.set_ylabel('Разброс R');diagnostic.legend(loc='lower left',fontsize=7)
            diagnostic.set_xlabel('β, °');ax.set_ylabel('Λ');ax.grid(alpha=.2);diagnostic.grid(alpha=.2)
            ax.legend(loc='upper center',bbox_to_anchor=(.5,-.03),ncol=3 if window=='A' else 2,fontsize=8)
            name=f'spring_robustness_{model.lower()}_'+('insensitivity' if window=='A' else 'crossing')
            fig.savefig(OUTPUT/(name+'.png'),dpi=300)
            fig.savefig(OUTPUT/(name+'.pdf'));plt.close(fig)
    result=dict(seconds=time.perf_counter()-started,solver_calls=0,matrix_calls=0,shape_calls=0,tracking_calls=0,
                input_sha256=sha(OUTPUT/'local_modes.csv'))
    if (OUTPUT/'render_manifest.json').exists():
        previous=json.loads((OUTPUT/'render_manifest.json').read_text(encoding='utf-8'))
        result['history']=previous.get('history',[])+[{k:v for k,v in previous.items() if k!='history'}]
    atomic(OUTPUT/'render_manifest.json',json.dumps(result,indent=2))
    print(json.dumps(result))


def verify_saved(run):
    """Small result audit, no solver, transfer, reconstruction or new tracking."""
    roots={r['shape_key']:r for p in run.state['points'].values() for r in p['roots']}
    assert len(roots)==sum(len(p['roots']) for p in run.state['points'].values())
    worst_mass=0.
    for key,r in roots.items():
        assert abs(r['Lambda']**2/r['Omega']-1)<1e-14
        assert abs(r['omega']*FS/r['Omega']-1)<1e-14
        assert r['k_theta']==r['kappa']*REF.D
        if 'M' in r:
            form=run.form(r)
            worst_mass=max(worst_mass,abs(np.vdot(form['vector'],form['vector']).real-1))
    for name in ('control_modes.csv','local_modes.csv','seed_mapping.csv'):
        for row in read_csv(OUTPUT/name):
            if not row.get('shape_key'):continue
            source=roots[row['shape_key']]
            for key in ('omega','Omega','Lambda'):
                assert float(row[key])==source[key],(name,row['shape_key'],key)
            assert row['current_sorted_position']==str(source['current_sorted_position'])
    costs=[p['full_B']+p['symmetry_B'] for p in run.state['points'].values()]
    assert max(costs)<=6000 and len(run.state['extra_points'])<=60
    assert len(run.state['points'])+len(run.state['auxiliary'])<=300
    assert worst_mass<1e-12
    assert source_hashes()==run.state['protected_sources']
    result=dict(root_records=len(roots),unique_shape_keys=len(roots),worst_mass_error=worst_mass,
        exact_frequency_copies=True,protected_sources_unchanged=True,protected_files=len(run.state['protected_sources']),
        maximum_recorded_point_B=max(costs),matrix_calls=0,solver_calls=0,
        note='read-only audit of saved data; mass products use saved shapes')
    atomic(OUTPUT/'saved_data_audit.json',json.dumps(result,indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('preflight','compute','base','local','plot-only','summarize'))
    parser.add_argument('--one',nargs=4,metavar=('MODEL','MU','KAPPA','BETA'))
    args = parser.parse_args()
    if args.command=='plot-only':
        render();return
    run = Run()
    run.state['commands'].append(dict(command=sys.argv[1:],time=time.strftime('%Y-%m-%dT%H:%M:%S')))
    if run.state['preflight'] is None:
        run.state['preflight'],run.properties = preflight();run.save()
    if args.command=='preflight':
        print(json.dumps(run.state['preflight'],indent=2));return
    if args.one:
        model,mu,kappa,beta = args.one
        run.point(model,float(mu),float(kappa),float(beta));return
    if args.command in ('compute','base'):run.base()
    if args.command in ('compute','local'):
        local_checks(run)
        refine_insensitivity(run)
    quadrature_checks(run)
    summarize(run)
    verify_saved(run)


if __name__=='__main__':
    main()
