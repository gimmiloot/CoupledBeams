"""Two-mode KV continuation, with a matrix stage and missing-only checkpoints.

New complex input/output and augmented correction cannot safely be a preset
of the real-frequency robustness runner. Native section/quadrature and old
elastic seeds are reused; no real inventory or complex contour search.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, replace
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from unittest.mock import patch

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_spring_modes as elastic
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_rlb as rlb

SOURCE = ROOT/'results/laminated_beams/inplane_spring_robustness'
THEORY = ROOT/'docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex'
OUTPUT = ROOT/'results/laminated_beams/inplane_kelvin_voigt_pilot'
BETA = np.deg2rad(5.)
ARRAYS = ('states', 'reactions', 'a', 'vector')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError('nonfinite JSON value')
    return value


def atomic(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+'.tmp')
    tmp.write_bytes(content if isinstance(content, bytes) else content.encode('utf-8'))
    os.replace(tmp, path)


def write_json(path, value):
    atomic(path, json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def protected_hashes():
    names = ['run_manifest.json','verified_roots.csv','control_modes.csv','diagnostics.json','shapes.npz']
    return {str((SOURCE/name).relative_to(ROOT)):sha(SOURCE/name) for name in names}


def error(left, right, rtol):
    absolute = float(np.linalg.norm(np.asarray(left)-np.asarray(right)))
    scale = float(np.linalg.norm(right))
    return dict(absolute=absolute, relative=absolute/max(scale,1e-300),
                criterion=kv.CRITERIA['matrix_atol']+rtol*scale,
                accepted=absolute <= kv.CRITERIA['matrix_atol']+rtol*scale)


def matrix_stage(properties, calls):
    checks = {}
    expected = dict(A=.011,D=2.979166666666667e-6,S=.0032051282051282055,
                    m=.01,J=2.083333333333334e-6)
    for name, value in expected.items():
        checks['section_'+name] = error(getattr(properties,name)/value, 1., 1e-12)
    old_expm = dict(EB=0,RLB=0)
    def wrap(name, fn):
        def counted(*args, **kwargs):
            old_expm[name] += 1
            return fn(*args, **kwargs)
        return counted
    with patch.object(eb,'expm',wrap('EB',eb.expm)), patch.object(rlb,'expm',wrap('RLB',rlb.expm)):
        for model in ('EB','RLB'):
            arm = kv.Arm.reduced(model,properties)
            provider = kv.Provider((arm,arm),BETA,1,0,calls)
            pair = elastic.arms(model,0.,properties)
            module = eb if model == 'EB' else rlb
            for Omega in (2.,20.,80.):
                z, omega = 1j*Omega, Omega/kv.T_REF
                key = f'{model}_elastic_Omega{Omega:g}'
                checks[key+'_H'] = error(kv.state_matrix(1j*omega,arm),pair[0].matrix(omega),1e-12)
                # Compare physical T in identical state scales, independent of any equilibration.
                units = arm.scale()
                new_T = provider.transfer(z,arm)
                old_T = module.transfer_matrix(omega,pair[0].native)
                checks[key+'_T'] = error(new_T*units[None,:]/units[:,None],
                                         old_T*units[None,:]/units[:,None],1e-9)
                old = elastic.assembly(omega,pair,BETA,eb.Joint('SPRING',provider.k))
                old_B = old.physical*provider.reaction_scales[None,:]/provider.row_units[:,None]
                checks[key+'_B'] = error(provider.matrices(z)[0],old_B,1e-9)
    for z in (-.08+2j,-.13+20j,.1+80j):
        e = kv.Arm.reduced('EB',properties)
        r = replace(kv.Arm.reduced('RLB',properties),invS=0.,J=0.)
        ep, rp = [kv.Provider((a,a),BETA,1,.003,calls) for a in (e,r)]
        checks[f'limit_{z}_H'] = error(kv.state_matrix(z/kv.T_REF,r),kv.state_matrix(z/kv.T_REF,e),1e-12)
        checks[f'limit_{z}_T'] = error(rp.transfer(z,r),ep.transfer(z,e),1e-9)
        checks[f'limit_{z}_B'] = error(rp.matrices(z)[0],ep.matrices(z)[0],1e-9)
        for model in ('EB','RLB'):
            a = kv.Arm.reduced(model,properties)
            provider = kv.Provider((a,a),BETA,1,.003,calls)
            B, Bz = provider.matrices(z,derivative=True)
            h = 1e-5*max(1.,abs(z))
            fd = (provider.matrices(z+h)[0]-provider.matrices(z-h)[0])/(2*h)
            checks[f'{model}_{z}_Bz'] = error(Bz,fd,kv.CRITERIA['derivative_rtol'])
            plain_B = provider.matrices(z)[0]
            checks[f'{model}_{z}_expm_frechet_value'] = error(B,plain_B,kv.CRITERIA['transfer_rtol'])
            checks[f'{model}_{z}_conjugate'] = error(provider.matrices(z.conjugate())[0],plain_B.conj(),1e-12)
    rng = np.random.default_rng(604)
    for index in range(3):
        y = rng.normal(size=12)+1j*rng.normal(size=12)
        p = (-.2+.4j)*(index+1)
        checks[f'joint_scalar_{index}'] = error(kv.joint_matrix(p,BETA,.002,.03)@y,
                                               kv.scalar_conditions(y,p,BETA,.002,.03),1e-12)
    return dict(checks=checks,accepted=all(v['accepted'] for v in checks.values()),
                old_transfer_expm=old_expm)


class Run:
    def __init__(self, output):
        self.output = Path(output).resolve()
        if self.output == SOURCE or self.output.is_relative_to(SOURCE):
            raise ValueError('old source directory is read-only')
        if not THEORY.is_file():
            raise FileNotFoundError('required theory source missing')
        self.current_hashes = protected_hashes()
        checkpoint = self.output/'diagnostics.json'
        self.data = json.loads(checkpoint.read_text(encoding='utf-8')) if checkpoint.exists() else dict(
            rows={},diagnostics={},costs={},times={},pair_checks={},quadrature={},
            protected_sources=self.current_hashes,initial_HEAD=subprocess.check_output(
                ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),theory_sha256=sha(THEORY),
            criteria=kv.CRITERIA,command_history=[],attempts=[])
        if self.data['protected_sources'] != self.current_hashes:
            raise ValueError('source files changed: inspect provenance before resuming')
        if self.data['criteria'] != kv.CRITERIA or self.data['theory_sha256'] != sha(THEORY):
            raise ValueError('criteria/theory changed: existing results are not silently replaced')
        self.shapes = {}
        if (self.output/'shapes.npz').exists():
            with np.load(self.output/'shapes.npz',allow_pickle=False) as archive:
                self.shapes = {key:archive[key] for key in archive.files}
        self.calls = {name:kv.Calls(**self.data['costs'].get(name,{}))
                      for name in ('matrix','seed','complex','quadrature','validation')}
        self.calls['complex'].limit = kv.CRITERIA['max_complex_evaluations']
        self.before = {name:c.snapshot() for name,c in self.calls.items()}
        self.session_new, self.session_reused = 0, 0
        self.section, self.properties = kv.section()

    @contextmanager
    def timed(self, stage):
        start = time.perf_counter()
        try:
            yield
        finally:
            self.data['times'][stage] = self.data['times'].get(stage,0.)+time.perf_counter()-start

    def provider(self, model, d, stage):
        arm = kv.Arm.reduced(model,self.properties)
        return kv.Provider((arm,arm),BETA,1.,d,self.calls[stage])

    def shape(self,key):
        return {name:self.shapes[key+'__'+name] for name in ARRAYS}

    def save(self):
        start = time.perf_counter()
        self.output.mkdir(parents=True,exist_ok=True)
        if self.shapes:
            stream = io.BytesIO()
            np.savez_compressed(stream,**self.shapes)
            atomic(self.output/'shapes.npz',stream.getvalue())
        rows = list(self.data['rows'].values())
        if rows:
            text = io.StringIO(newline='')
            writer = csv.DictWriter(text,fieldnames=list(dict.fromkeys(k for row in rows for k in row)))
            writer.writeheader()
            writer.writerows(rows)
            atomic(self.output/'modal_results.csv',text.getvalue())
        self.data['costs'] = {name:c.snapshot() for name,c in self.calls.items()}
        self.data['times']['write'] = self.data['times'].get('write',0.)+time.perf_counter()-start
        write_json(self.output/'diagnostics.json',self.data)  # commit marker written last
        self.manifest()

    def manifest(self):
        rows = list(self.data['rows'].values())
        newcalls = {stage:{k:v-self.before[stage][k] for k,v in c.snapshot().items() if k!='limit'}
                    for stage,c in self.calls.items()}
        value = dict(initial_HEAD=self.data['initial_HEAD'],working_tree_sources=True,
            theory=dict(path=str(THEORY.relative_to(ROOT)),sha256=self.data['theory_sha256'],
                origin='User supplied C:/Users/Nikita/Downloads/inplane_kelvin_voigt_joint_theory.tex; version 1.0, 2026-09-12'),
            environment=dict(executable=sys.executable,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
            configuration=dict(beta0_deg=5.,beta0_rad=BETA,L1=1.,L2=1.,b=.2,h=.05,
                layup='H/L/L/H',contrast=.4,K=5/6,kappa_theta=1.,D_ref=kv.REFERENCE.D,
                m_ref=kv.REFERENCE.m,t_ref=kv.T_REF,properties=asdict(self.properties)),
            policy=dict(frequency_map_policy='frequency-map-v1',calculation_mode='bounded_complex_continuation',
                spectrum_semantics='tracked_branches',policy_override_reason='Explicit two-mode complex pilot; no real guard or complete complex spectrum',
                start='saved elastic modes at beta0=5,mu=0,kappa=1; ancestry retained',
                max_B_Bz=2000,max_steps=20,max_retries_per_branch=2,max_intermediate=4),
            criteria=kv.CRITERIA,costs=self.data['costs'],times=self.data['times'],
            completed_states=len(rows),confirmed=sum(r['status']=='CONFIRMED' for r in rows),
            qualified=sum(r['status']!='CONFIRMED' for r in rows),
            elastic_reused=sum(r['d_index']==0 for r in rows),
            positive_d=sum(r['d_index']>0 for r in rows),
            session_new=self.session_new,session_reused=self.session_reused,session_calls=newcalls,
            protected_sources=self.data['protected_sources'],protected_sources_unchanged=protected_hashes()==self.current_hashes,
            source_versions={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'scripts/lib/inplane_kelvin_voigt.py',THEORY]},
            source_elastic_HEAD=self.data.get('source_elastic_HEAD'),
            seed_manifest='seed_manifest.json' if (self.output/'seed_manifest.json').exists() else None,
            intermediate_solutions=0,extra_attempts=0,commands=self.data['command_history'],
            limitations=['two isolated descendants only','one laminate/angle/joint stiffness',
                'residuals are not eigenvalue error bounds','old source/Ritz qualifications retained',
                'no independent FEM/Ritz or full complex-spectrum verification'])
        write_json(self.output/'run_manifest.json',value)

    def preflight(self, recheck=False):
        if 'matrix_stage' in self.data and not recheck:
            return self.data['matrix_stage']['accepted']
        if 'matrix_stage' in self.data:
            self.data.setdefault('matrix_stage_history',[]).append(self.data['matrix_stage'])
        with self.timed('matrix'):
            self.data['matrix_stage'] = matrix_stage(self.properties,self.calls['matrix'])
            native = self.section
            self.data['section'] = dict(A=native.A,B=native.B,D=native.D,shear=native.shear,
                I0=native.I0,I1=native.I1,I2=native.I2,properties=asdict(self.properties))
        self.save()
        return self.data['matrix_stage']['accepted']

    def seeds(self):
        file = self.output/'seed_manifest.json'
        if file.exists():
            self.session_reused += 4
            return json.loads(file.read_text(encoding='utf-8'))
        old = json.loads((SOURCE/'run_manifest.json').read_text(encoding='utf-8'))
        p = old['preflight']
        assert p['layup']=='H/L/L/H' and p['contrast']==.4 and p['ply_angles_deg']==[0]*4
        assert p['ply_thickness']==.0125
        for key in ('A','D','S','m','J','K','width'):
            np.testing.assert_allclose(p['properties'][key],getattr(self.properties,key),rtol=1e-12,atol=0.)
        for key,value in [('l_ref',1.),('b',.2),('h',.05),('D_ref',kv.REFERENCE.D),
                          ('m_ref',kv.REFERENCE.m),('frequency_scale',kv.T_REF)]:
            np.testing.assert_allclose(old['parameters'][key],value,rtol=1e-12,atol=0.)
        roots, mappings = read_csv(SOURCE/'verified_roots.csv'),read_csv(SOURCE/'control_modes.csv')
        diagnostics = json.loads((SOURCE/'diagnostics.json').read_text(encoding='utf-8'))
        self.data['source_elastic_HEAD'] = old['created_HEAD']
        seeds = {}
        with np.load(SOURCE/'shapes.npz',allow_pickle=False) as archive:
            for model in ('EB','RLB'):
                point=f'{model}_m0_k1_b5'
                assert diagnostics['points'][point]['status']=='CONFIRMED'
                pool=sorted([r for r in roots if r['point_id']==point],key=lambda r:int(r['current_sorted_position']))
                assert all(float(a['Omega'])<=float(b['Omega']) for a,b in zip(pool,pool[1:]))
                for role,eta in [('ACTIVE',1),('INACTIVE',-1)]:
                    r=next(r for r in pool if int(r['symmetry_class'])==eta and r['root_status']=='CONFIRMED')
                    m=next(m for m in mappings if m['shape_key']==r['shape_key'] and m['mapping_status']=='CONFIRMED')
                    Omega=float(r['Omega'])
                    assert abs(float(r['Lambda'])**2-Omega)<1e-12*Omega
                    assert abs(float(r['omega'])*kv.T_REF-Omega)<1e-12*Omega
                    neighbours=[float(q['Omega']) for q in pool if q['shape_key']!=r['shape_key']]
                    assert min(abs(x-Omega) for x in neighbours)>.01*Omega
                    saved_diag=next(q for q in diagnostics['points'][point]['roots'] if q['shape_key']==r['shape_key'])
                    assert saved_diag['root_status']=='CONFIRMED' and saved_diag['Omega']==Omega
                    key=f'{model}_{role}_d0'
                    if key not in self.data['rows']:
                        provider=self.provider(model,0.,'seed')
                        z=1j*Omega
                        with self.timed('seed_matrix'):
                            B,_=provider.matrices(z)
                            a=kv.right_null(B)  # full complex SVD, no parity projection
                        with self.timed('seed_shapes'):
                            shape=kv.recover(provider,z,a,direct_check=True)
                        oldshape=archive[r['shape_key']+'__states']
                        _,weights=kv.quadrature()
                        oldvector=kv.mass_vector(oldshape,provider.arms,weights)
                        MAC=float(kv.mac_matrix([oldvector],[shape['vector']])[0,0])
                        self.record(key,provider,z,shape,role,0,m['branch_id'],Omega,MAC,
                                    dict(steps=0,status='REUSED_ELASTIC_SEED',last_delta_z=None,history=[]),0.)
                    row=self.data['rows'][key]
                    if row['status']!='CONFIRMED':
                        raise RuntimeError('ZERO_DAMPING_QUALIFIED: no complex continuation from this seed')
                    diag=self.data['diagnostics'][key]
                    delta=complex(diag['Delta_psi']['real'],diag['Delta_psi']['imag']) if isinstance(diag['Delta_psi'],dict) else diag['Delta_psi']
                    s0=kv.M_REF*abs(delta)**2/(float(r['omega'])**2*diag['M_phi'])
                    if role=='ACTIVE' and abs(delta)<kv.CRITERIA['inactive_delta']:
                        raise RuntimeError('ACTIVE_SEED_IS_INACTIVE')
                    seeds[f'{model}_{role}']=dict(Omega0=Omega,omega0=float(r['omega']),s0=s0,
                        gamma=.5*Omega*s0,branch_id=m['branch_id'],symmetry_class=eta,
                        shape_key=r['shape_key'],source_row=r,source_mapping=m,
                        point_status=diagnostics['points'][point]['status'],
                        min_neighbour_gap=min(abs(x-Omega) for x in neighbours),seed_key=key)
        gammas=[s['gamma'] for key,s in seeds.items() if key.endswith('_ACTIVE')]
        d_star=min(.01/max(gammas),.1/max(s['Omega0'] for s in seeds.values()))
        result=dict(source_files=self.current_hashes,source_HEAD=old['created_HEAD'],
            seeds=seeds,d_star=d_star,d_theta=[0.,d_star/10,d_star/2,d_star],
            c_theta=[d*kv.M_REF*kv.T_REF for d in (0.,d_star/10,d_star/2,d_star)],
            selection='min(.01/max(gamma_ACTIVE),.1*kappa/max(Omega0)); fixed before complex search')
        write_json(file,result)
        self.save()
        return clean(result)

    def record(self,key,provider,z,shape,role,index,branch,Omega0,MAC,correction,a_lin):
        with self.timed('root_diagnostics'):
            diag=kv.diagnose(provider,z,shape)
        failures=kv.failures(diag,z,role,Omega0,MAC)
        if correction['status'] not in ('CONVERGED','REUSED_ELASTIC_SEED'):
            failures.append(correction['status'])
        if index>0 and role=='ACTIVE' and -z.real<=kv.CRITERIA['a_atol']:
            failures.append('ACTIVE_DECAY_UNRESOLVED')
        p=z/kv.T_REF
        row=dict(key=key,model=provider.arms[0].model,branch_id=branch,role=role,
            beta0_deg=5.,kappa_theta=1.,d_index=index,d_theta=provider.d,k_theta=provider.k,c_theta=provider.c,
            Omega0=Omega0,omega0=Omega0/kv.T_REF,p_real=p.real,p_imag=p.imag,z_real=z.real,z_imag=z.imag,
            alpha=-p.real,omega_d=p.imag,a_decay=-z.real,Omega_d=z.imag,
            Lambda_d=np.sqrt(z.imag) if z.imag>0 else None,a_lin=a_lin,alpha_energy=diag['alpha_energy'],
            Delta_psi_real=diag['Delta_psi'].real,Delta_psi_imag=diag['Delta_psi'].imag,
            M_phi=diag['M_phi'],K_phi=diag['K_phi'],C_phi=diag['C_phi'],MAC=MAC,
            r_B=diag['null_residual'],sigma_ratio=diag['sigma_ratio'],r_E=diag['r_E'],
            physical_residual=max(diag['physical_residuals']),iterations=correction['steps'],
            status='QUALIFIED' if failures else 'CONFIRMED',failures=';'.join(failures),
            provenance='ZERO_DAMPING_VALIDATION' if index==0 else
                ('PROTECTED_MODE_VALIDATION' if role=='INACTIVE' else 'COMPLEX_CONTINUATION'))
        self.data['rows'][key]=clean(row)
        self.data['diagnostics'][key]=clean(dict(diag,correction=correction,
            direct_transfer_errors=shape.get('direct_errors',[])))
        for name in ARRAYS:
            self.shapes[key+'__'+name]=shape[name]
        self.session_new+=1
        self.save()
        print(key,row['status'],f'Omega_d={z.imag:.12g} a={-z.real:.12g} steps={correction["steps"]}',flush=True)

    def continuation(self,seeds):
        for model in ('EB','RLB'):
            for index,d in enumerate(seeds['d_theta'][1:],1):
                keys=[]
                for role in ('ACTIVE','INACTIVE'):
                    key=f'{model}_{role}_d{index}'
                    keys.append(key)
                    if key in self.data['rows']:
                        self.session_reused+=1
                        continue
                    seed=seeds['seeds'][f'{model}_{role}']
                    previous=f'{model}_{role}_d{index-1}'
                    if previous not in self.data['rows'] or self.data['rows'][previous]['status']!='CONFIRMED':
                        continue
                    oldrow=self.data['rows'][previous]
                    zprevious=complex(oldrow['z_real'],oldrow['z_imag'])
                    a_lin=.5*d*seed['Omega0']**2*seed['s0']
                    predicted=1j*seed['Omega0']-a_lin if index==1 else zprevious
                    provider=self.provider(model,d,'complex')
                    oldshape=self.shape(previous)
                    # A rejected attempt is recorded and qualified, without automatic cascades.
                    with self.timed('complex_search'):
                        correction=kv.correct(provider.matrices,predicted,oldshape['a'])
                    self.calls['complex'].corrections+=correction['steps']
                    with self.timed('complex_shapes'):
                        shape=kv.recover(provider,correction['z'],correction['a'])
                    overlap=np.vdot(oldshape['vector'],shape['vector'])
                    if abs(overlap):
                        phase=np.conj(overlap)/abs(overlap)
                        for name in ARRAYS:
                            shape[name]*=phase
                    MAC=float(kv.mac_matrix([oldshape['vector']],[shape['vector']])[0,0])
                    self.record(key,provider,correction['z'],shape,role,index,seed['branch_id'],
                                seed['Omega0'],MAC,correction,a_lin)
                if f'{model}_d{index}' not in self.data['pair_checks'] and all(key in self.data['rows'] for key in keys):
                    from scipy.optimize import linear_sum_assignment
                    previous=[self.shape(f'{model}_{role}_d{index-1}')['vector'] for role in ('ACTIVE','INACTIVE')]
                    current=[self.shape(key)['vector'] for key in keys]
                    mac=kv.mac_matrix(previous,current)
                    _,assignment=linear_sum_assignment(1-mac)
                    a,b=[self.data['rows'][key] for key in keys]
                    gap=abs(complex(a['z_real']-b['z_real'],a['z_imag']-b['z_imag']))
                    ok=list(assignment)==[0,1] and gap>1e-6*max(a['Omega_d'],b['Omega_d'])
                    self.data['pair_checks'][f'{model}_d{index}']=dict(MAC=mac.tolist(),assignment=assignment.tolist(),
                        distinct_root_gap=gap,accepted=bool(ok))
                    if not ok:
                        for key in keys:
                            self.data['rows'][key]['status']='QUALIFIED'
                            self.data['rows'][key]['failures']+=';PAIR_ASSIGNMENT_OR_COLLISION'
                    self.save()

    def quadrature(self):
        for model in ('EB','RLB'):
            key=f'{model}_ACTIVE_d3'
            if key in self.data['quadrature'] or key not in self.data['rows']:
                continue
            row=self.data['rows'][key]
            oldshape=self.shape(key)
            z=complex(row['z_real'],row['z_imag'])
            provider=self.provider(model,row['d_theta'],'quadrature')
            with self.timed('quadrature'):
                shape=kv.recover(provider,z,oldshape['a'],nodes=257,direct_check=True)
                diag=kv.diagnose(provider,z,shape)
                _,w=kv.quadrature()
                sampled=kv.mass_vector(shape['states'][:,::2],provider.arms,w)
                MAC=float(kv.mac_matrix([oldshape['vector']],[sampled])[0,0])
                value=dict(relative_mass_change=abs(shape['mass_before_normalization']-1),MAC=MAC,
                    direct_transfer_errors=shape['direct_errors'],r_E_257=diag['r_E'],
                    a_energy_257=diag['a_energy'],a_energy_129=self.data['diagnostics'][key]['a_energy'])
                value['accepted']=bool(value['relative_mass_change']<=kv.CRITERIA['quadrature_rtol']
                    and 1-MAC<=kv.CRITERIA['quadrature_rtol']
                    and max(value['direct_transfer_errors'])<=kv.CRITERIA['transfer_rtol'])
            self.data['quadrature'][key]=value
            self.save()

    def validate_saved(self):
        """Addressed seed-matrix comparisons and stored-data audit, no new roots."""
        if 'saved_validation' in self.data:
            return
        with self.timed('saved_validation'):
            comparisons, legacy_calls = {}, dict(EB=0,RLB=0)
            def wrap(model,fn):
                def counted(*args,**kwargs):
                    legacy_calls[model]+=1
                    return fn(*args,**kwargs)
                return counted
            with patch.object(eb,'expm',wrap('EB',eb.expm)),patch.object(rlb,'expm',wrap('RLB',rlb.expm)):
                for model in ('EB','RLB'):
                    for role in ('ACTIVE','INACTIVE'):
                        key=f'{model}_{role}_d0'
                        row=self.data['rows'][key]
                        provider=self.provider(model,0.,'validation')
                        B,_=provider.matrices(1j*row['Omega0'])
                        pair=elastic.arms(model,0.,self.properties)
                        old=elastic.assembly(row['omega0'],pair,BETA,eb.Joint('SPRING',provider.k))
                        target=old.physical*provider.reaction_scales[None,:]/provider.row_units[:,None]
                        comparisons[key]=error(B,target,kv.CRITERIA['transfer_rtol'])
            conjugate_physical={}
            for key,row in self.data['rows'].items():
                provider=self.provider(row['model'],row['d_theta'],'validation')
                y=self.shape(key)['states'][:,-1,:].ravel()
                amp=max(abs(y/provider.state_units))
                p=complex(row['p_real'],row['p_imag'])
                residual=abs(kv.scalar_conditions(y.conj()/amp,p.conjugate(),BETA,
                                                provider.k,provider.c)/provider.row_units)
                conjugate_physical[key]=float(max(residual))
            self.data['saved_validation']=dict(seed_matrices=comparisons,legacy_expm=legacy_calls,
                conjugate_physical=conjugate_physical,accepted=bool(
                    all(v['accepted'] for v in comparisons.values())
                    and max(conjugate_physical.values())<=kv.CRITERIA['physical_residual']),
                protected_sources_unchanged=protected_hashes()==self.current_hashes)
        self.save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix-only',action='store_true')
    parser.add_argument('--compute',action='store_true',help='matrix stage, seeds, then missing-only EB/RLB continuation')
    parser.add_argument('--recheck-matrices',action='store_true',help='explicit matrix-only recheck; retains original diagnostics')
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args()
    if not (args.compute or args.matrix_only):
        parser.error('choose --matrix-only or --compute')
    run=Run(args.output)
    run.data['command_history'].append(' '.join(sys.argv))
    try:
        if not run.preflight(args.recheck_matrices):
            raise RuntimeError('MATRIX_STAGE_FAILED: spectrum not attempted')
        if args.compute:
            seeds=run.seeds()
            run.continuation(seeds)
            run.quadrature()
            run.validate_saved()
        # Completed scientific checkpoints are immutable under a no-op resume.
        if run.session_new or any(c.B!=run.before[k]['B'] for k,c in run.calls.items()):
            run.save()
        else:
            run.manifest()
    except Exception as exc:
        run.data['attempts'].append(dict(status='STOPPED',reason=str(exc),type=type(exc).__name__))
        run.save()
        raise
    print(json.dumps(dict(states=len(run.data['rows']),new=run.session_new,reused=run.session_reused,
        stage_costs={k:c.snapshot() for k,c in run.calls.items()})))


if __name__=='__main__':
    main()
