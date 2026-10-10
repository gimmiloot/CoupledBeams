"""FEM-3AR: explicitly authorized continuation after native static ELKE failure.

Separate parent/permission/attempt contract preserves the old failed CLI guard.
Existing FEM-3A generation, monitored jobs, recovery, references and plots reused.
"""
from __future__ import annotations
import argparse,hashlib,json,math,os,re,shutil,subprocess,sys,time
from pathlib import Path
if __name__=='__main__':
    for p in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[p]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
from scripts.analysis import pilot_nlsp_nonlinear_dynamic_3d_fem as base
CONFIG=ROOT/'data/input/nlsp_nonlinear_dynamic_3d_fem_resume.json'
OUTPUT=ROOT/'results/nlsp_nonlinear_dynamic_3d_fem_resume'
read_json,write_json,sha=base.read_json,base.write_json,base.sha
STATUS_NAMES=('SOURCE_PRESERVATION','INPUT_REMEDIATION','LINEAR_PRELOAD','LINEAR_RELEASE','LINEAR_DYNAMIC','NONLINEAR_PRELOAD','NONLINEAR_RELEASE','NONLINEAR_DYNAMIC','1D_REFERENCE','TRANSIENT_RECOVERY','ENERGY_DIAGNOSTICS','SHORT_RESPONSE_COMPARISON')


def validate_config(c):
    if c['schema']!='nlsp-fem3ar-continuation-v1':raise ValueError('Continuation schema mismatch')
    if c['authorization']!={'id':'explicit_user_FEM3AR_2026_10_09','maximum_production_CCX_jobs':2,'maximum_nonlinear_1D_integrations':1,'case_order':['linear','nonlinear'],'automatic_retry':False,'basis':'Explicit user FEM-3AR continuation after native ELKE/static failure'}:raise ValueError('Separate explicit authorization missing')
    if c['threads']!=1 or c['job_timeout_seconds']!=1200 or c['job_memory_limit_bytes']!=4*1024**3 or c['numerical_budget_seconds']!=3600:raise ValueError('Frozen resource limits changed')
    if c['science_policy']!='reuse_immutable_parent_science_config_and_p64_static_coordinates' or any(c[k] for k in ('new_meshes','new_modal_jobs','new_static_only_jobs','new_time_or_space_levels')):raise ValueError('Unauthorized study extension')
    return c


def load_parent(c):
    parent=ROOT/c['parent_failed']['bundle']
    if sha(parent/'manifest.json')!=c['parent_failed']['manifest_sha256']:raise ValueError('Failed parent manifest changed')
    old=base.validate_cache(parent);item=read_json(parent/'provenance.json');science=base.validate_config(item['config'])
    if old['overall']!='BLOCKED_BY_SOLVER' or old['job_calls']['CCX_production']!=1 or old['cases']['linear']['status']!='FAIL':raise ValueError('Wrong historical attempt')
    if sha(base.__file__)!=c['corrected_generator_sha256']:raise ValueError('Corrected generator identity changed')
    return parent,item,old,science


def output_safety(text):
    first=text.split('*END STEP',1)[0];step=next(l.upper() for l in first.splitlines() if l.startswith('*STEP'))
    nl=any(p.strip()=='NLGEOM' or p.strip()=='NLGEOM=YES' for p in step.split(',')[1:])
    if not nl and 'ELKE' in first:raise ValueError('Unsafe linear STATIC ELKE route, including NLGEOM=NO')
    dynamic=text.split('*END STEP',1)[1]
    if 'ELSE,ELKE' not in dynamic:raise ValueError('Dynamic internal/kinetic output missing')
    return {'status':'PASS','nonlinear_static':nl,'linear_static_ELKE_excluded':not nl,'dynamic_ELKE_requested':True,'native_NL_STATIC_velocity_safety':'source_input_gate.json'}


def identity(config_path=CONFIG):
    c=validate_config(read_json(config_path));parent,_,_,science=load_parent(c)
    _,frozen=base.identity(base.CONFIG)
    item={**frozen,'config':science,'continuation_config':c,'continuation_code_sha256':sha(__file__),'continuation_config_sha256':sha(config_path),'parent_failed_manifest_sha256':sha(parent/'manifest.json'),'parent_failed_bundle':c['parent_failed']['bundle'],'authorization':c['authorization'],'source_gate_sha256':sha(ROOT/'results/_smoke/fem3ar_source/source_input_gate.json'),'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def update_statuses(s):
    for name in STATUS_NAMES:s['resume_statuses'].setdefault('NLSP_FEM3AR_'+name,'NOT_RUN')
    for kind in ('linear','nonlinear'):
        r=s['cases'].get(kind,{})
        s['resume_statuses']['NLSP_FEM3AR_'+kind.upper()+'_DYNAMIC']='PASS' if r.get('status')=='PASS' and r.get('continuation_audit',{}).get('status')=='PASS' else 'FAIL' if r.get('status')=='FAIL' else 'PARTIAL' if r else 'NOT_RUN'
        a=r.get('continuation_audit',{})
        for phase in ('PRELOAD','RELEASE'):
            s['resume_statuses']['NLSP_FEM3AR_'+kind.upper()+'_'+phase]='PASS' if a.get('status')=='PASS' else 'NOT_RUN'
    if (s.get('overall')=='PILOT_COMPLETE_WITH_QUALIFICATIONS' and s.get('comparison')):
        s['resume_statuses'].update(NLSP_FEM3AR_1D_REFERENCE='PASS',NLSP_FEM3AR_TRANSIENT_RECOVERY='PASS',NLSP_FEM3AR_SHORT_RESPONSE_COMPARISON='PASS',NLSP_FEM3AR_ENERGY_DIAGNOSTICS='PASS' if all(s['cases'][k]['energy_status']=='PASS' for k in ('linear','nonlinear')) else 'PARTIAL')


def save(b,item,s):
    update_statuses(s);base.finalize(b,item,s)


def validate_cache(b):
    s=base.validate_cache(b);m=read_json(Path(b)/'manifest.json');load_parent(m['identity']['continuation_config']);return s


def existing_attempt(c):
    if not OUTPUT.exists():return None
    for b in sorted(OUTPUT.iterdir()):
        if not b.is_dir() or not (b/'provenance.json').exists():continue
        item=read_json(b/'provenance.json')
        if item.get('authorization',{}).get('id')!=c['authorization']['id']:continue
        if item['continuation_config']!=c:raise ValueError('Authorization already used with different continuation settings')
        if not (b/'manifest.json').exists():raise RuntimeError('Unmanifested interrupted continuation; no hidden retry')
        return b,item,validate_cache(b)
    return None


def prepare(config_path=CONFIG):
    c=validate_config(read_json(config_path));parent,_,old,science=load_parent(c);key,item=identity(config_path);b=OUTPUT/key
    if b.exists() and any(b.iterdir()):raise RuntimeError('Existing continuation cannot be overwritten')
    b.mkdir(parents=True);write_json(b/'provenance.json',item);write_json(b/'continuation_config.json',c);write_json(b/'config.json',science)
    shutil.copyfile(parent/'protocol_evidence.json',b/'protocol_evidence.json');shutil.copyfile(ROOT/'results/_smoke/fem3ar_source/source_input_gate.json',b/'source_input_gate.json')
    (b/'execution_code').mkdir()
    for path in (Path(__file__),Path(base.__file__),Path(base.one.__file__),Path(base.io.__file__)):shutil.copyfile(path,b/'execution_code'/path.name)
    ref=base.one.load_reference(ROOT,science['source_static']['bundle'],science['source_fem1']['bundle'],science['source_action']['bundle'])
    pre=base.one.preflight_reference(ref,science)
    oldpre=read_json(parent/'one_d_preflight.json')
    for kind in ('linear','nonlinear'):
        if pre['initial_states'][kind]['q0']!=oldpre['initial_states'][kind]['q0'] or pre['initial_states'][kind]['v0']!=oldpre['initial_states'][kind]['v0']:raise ValueError('Frozen initial state differs')
    write_json(b/'one_d_preflight.json',pre)
    parent_static,source,mesh,audit=base.verify_sources(science);(b/'input_gate').mkdir();gates={}
    for kind in ('linear','nonlinear'):
        path=b/'input_gate'/(kind+'.inp');gate=base.write_input(path,science,parent_static,source,mesh,audit,kind=='nonlinear');gates[kind]={**gate,**output_safety(path.read_text(encoding='utf8'))}
    ignore=lambda text:'\n'.join(row for row in text.splitlines() if not row.startswith('*INCLUDE'))
    if ignore((b/'input_gate/linear.inp').read_text())!=ignore((parent/'remediation_preview/linear_corrected_NOT_RUN.inp').read_text()):raise ValueError('Current deck differs from authorized corrected preview')
    gates['preview_reproduced_except_INCLUDE']=True;write_json(b/'input_gate.json',gates)
    s={'authorization':c['authorization'],'parent_failed':c['parent_failed'],'statuses':{'NLSP_FEM3A_'+n:'NOT_RUN' for n in base.STATUS_NAMES},'resume_statuses':{'NLSP_FEM3AR_'+n:'NOT_RUN' for n in STATUS_NAMES},'cases':{},'attempts':[],'job_calls':{'CCX_production':0,'CCX_fixture':0,'Gmsh':0,'1D_nonlinear_ODE':0,'1D_static':0,'physical_root_search':0,'symbolic_derivations':0},'numerical_seconds':0.,'overall':'NOT_RUN','preflight':{'status':'PASS','omega1':pre['omega1'],'T1':pre['T1'],'target_dynamic_end':pre['horizon']},'strict_float64_qualification':'PARTIAL','execution_mode':'EXPLORATORY_NOT_CERTIFIED','admitted':False}
    s['resume_statuses'].update(NLSP_FEM3AR_SOURCE_PRESERVATION='PASS',NLSP_FEM3AR_INPUT_REMEDIATION='PASS');save(b,item,s);return b,item,s



def qualify_native_energy(case,record):
    """Keep actual native energies; never replace the static-origin jump."""
    reported=read_json(case/'stdout_energy.json')
    rows=[r for r in reported['records'] if r.get('step')==2]
    E0=record['native_initial_internal_energy']
    references=[r['initial_step_energy'] for r in rows if 'initial_step_energy' in r]
    if not references or not E0:
        record['energy_status']='PARTIAL';record['energy_qualification']='Missing native energy reference; no value reconstructed';return
    reference=references[0]
    actual=read_json(case/'energy.json');dynamic=[r for r in actual['records'] if r.get('step')==2 and 'mechanical_energy' in r]
    record['native_dynamic_bookkeeping_initial_energy']=reference
    record['native_energy_static_to_dynamic_reference_jump_relative']=(reference-E0)/E0
    record['max_native_dynamic_bookkeeping_reference_relative_drift']=max((abs(r['mechanical_energy']/reference-1) for r in dynamic),default=None)
    record['native_energy_history_not_offset_corrected']=True
    record['energy_status']='PARTIAL'
    record['energy_qualification']='Native STATIC-origin energy jump retained; installed initial dynamic results use zero stress/strain bookkeeping baselines plus copied preload enerini. Separate drift to native DYNAMIC bookkeeping reference is a diagnostic, not physical energy continuity or certified time accuracy. No RHS/energy correction.'
    write_json(case/'energy_reference_qualification.json',{k:record[k] for k in record if k.startswith('native_') or k.startswith('max_native_') or k.startswith('energy_')})
    write_json(case/'recovery.json',record)


def record_reference_authorization(b,authorization):
    """Record the real continuation permission while retaining the legacy call."""
    for kind in ('linear','nonlinear'):
        path=Path(b)/('one_d_'+kind+'.json')
        if not path.exists():continue
        metadata=read_json(path)
        metadata['continuation_authorization']=authorization
        metadata['authorization_qualification']='FEM-3AR is the actual new permission. The inherited FEM-3A callback token is retained as execution provenance, not treated as a second authorization or another integration.'
        write_json(path,metadata)

def audit_case(b,science,s,kind,*,native_time_rounding=False):
    old,_,mesh,audit=base.verify_sources(science);case=b/'cases'/kind
    lines=(case/'motion.stdout.txt').read_text(encoding='utf8').splitlines();warnings=[v for v in lines if '*WARNING' in v.upper()]
    if warnings:raise ValueError('Unexplained native warnings: '+str(warnings))
    sta=base.io.read_transient_sta(case/'motion.sta');static=[r for r in sta['accepted_increments'] if r['step']==1];dynamic=[r for r in sta['accepted_increments'] if r['step']==2]
    dynamic_bound=6e-7
    if native_time_rounding and dynamic:
        dynamic_bound=dynamic[-1]['time_rounding_bounds']['step_time']+8*np.finfo(float).eps*max(1.,abs(dynamic[-1]['step_time']))
    if not static or not dynamic or abs(static[-1]['step_time']-1)>1e-6 or abs(dynamic[-1]['step_time']-base.dynamic_settings(science)['duration'])>dynamic_bound:raise ValueError('Full static load or target dynamic time missing')
    inc=static[-1]['increment'];ids,xyz,gravity,volume=base.base.consistent_gravity_loads(mesh,science['material']['rho'],science['g']);index={int(n):i for i,n in enumerate(ids)}
    with np.load(case/'frames'/f'step1_inc{inc:05d}.npz') as z:U=z['U'].copy()
    resultants={};R=[];positions=[]
    for name,key in [('LEFT_FIXED','fixed_left_ids'),('RIGHT_FIXED','fixed_right_ids')]:
        rows=np.array([index[int(n)] for n in audit[key]])
        with np.load(case/'dat_fields'/f'1_{inc}_{name}_FORC.npz') as z:RF=z['values'].copy()
        support=RF-gravity[rows];current=xyz[rows]+U[rows] if kind=='nonlinear' else xyz[rows]
        resultants[name]={'force':support.sum(axis=0),'moment_about_face_centroid':np.cross(current-current.mean(axis=0),support).sum(axis=0)};R.append(support);positions.append(current)
    force=gravity.sum(axis=0);scale=np.linalg.norm(force);actual_current=xyz+U if kind=='nonlinear' else xyz
    imbalance=np.linalg.norm(np.vstack(R).sum(axis=0)+force)/scale
    moment=np.linalg.norm(np.cross(np.vstack(positions),np.vstack(R)).sum(axis=0)+np.cross(actual_current,gravity).sum(axis=0))/scale
    gate=old['science_config']['gates']['equilibrium_relative']
    if imbalance>gate or moment>gate:raise ValueError('Unchanged independent preload equilibrium gate failed')
    rec=s['cases'][kind]
    if rec.get('maximum_native_external_work_after_release') is None or rec.get('maximum_native_damping_work_after_release') is None:raise ValueError('Missing native release/work evidence')
    if rec['maximum_native_external_work_after_release']!=0 or rec['maximum_native_damping_work_after_release']!=0:raise ValueError('Nonzero native external/damping work after release')
    state=rec['final_strain_diagnostics']
    if not state['finite_values'] or state['minimum_det_deformation_gradient']<=0:raise ValueError('Nonfinite/inverted final deformation diagnostic')
    qualify_native_energy(case,rec)
    if not np.isfinite(U).all():raise ValueError('Nonfinite preload field')
    write_json(case/'execution_environment.json',{'threads':1,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMBER_OF_CPUS':'1','PATH_unchanged':True})
    result={'status':'PASS','warning_lines':warnings,'total_applied_force':force,'reference_volume':volume,'support_resultants':resultants,'preload_force_imbalance_relative':float(imbalance),'preload_moment_imbalance_relative':float(moment),'equilibrium_gate_unchanged':gate,'bodyload_recovered_independently':True,'source_p64_states_not_reprojected':True,'static_end_total_time':static[-1]['total_time'],'actual_dynamic_end_STA':dynamic[-1]['step_time'],'initial_velocity_evidence':'explicit zero IC/native source initialization; first native frame occurs at positive time','initial_tiny_native_regularization_qualified':True}
    write_json(case/'continuation_audit.json',result);rec['continuation_audit']=result;return result


def main(argv=None):
    arguments=list(sys.argv[1:] if argv is None else argv)
    if '--validation' in arguments:
        from scripts.lib import nlsp_fem3c_validation
        arguments.remove('--validation')
        return nlsp_fem3c_validation.main(arguments)
    if '--long-horizon' in arguments:
        from scripts.lib import nlsp_fem3b_continuation
        arguments.remove('--long-horizon')
        return nlsp_fem3b_continuation.main(arguments)
    p=argparse.ArgumentParser(description=__doc__);m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--preflight',action='store_true');m.add_argument('--run-pilot',action='store_true');m.add_argument('--report-only',type=Path);m.add_argument('--plot-only',type=Path)
    p.add_argument('--config',type=Path,default=CONFIG);p.add_argument('--through-case',choices=['linear','nonlinear'],default='nonlinear');a=p.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_cache(b)
        if a.plot_only:plot_bundle(b)
        print(json.dumps({'bundle':str(b),'statuses':s['resume_statuses'],'overall':s['overall'],'new_scientific_calls':0},indent=2));return s
    c=validate_config(read_json(a.config));found=existing_attempt(c)
    if found:b,item,s=found
    else:b,item,s=prepare(a.config)
    science=item['config']
    if a.preflight or s.get('hard_stop') or s['overall'] in ('PILOT_COMPLETE_WITH_QUALIFICATIONS','BLOCKED_BY_SOLVER'):
        print(json.dumps({'bundle':str(b),'statuses':s['resume_statuses'],'overall':s['overall'],'new_scientific_calls':0},indent=2));return s
    for kind in ('linear','nonlinear'):
        if not base.run_case(b,science,item,s,kind):break
        try:audit_case(b,science,s,kind)
        except Exception as exc:
            s['hard_stop']=True;s['cases'][kind]['status']='FAIL';s['cases'][kind]['audit_failure']=str(exc);s['overall']='PARTIAL';save(b,item,s);break
        save(b,item,s)
        if a.through_case==kind:break
    if all(s['cases'].get(k,{}).get('status')=='PASS' and s['cases'][k].get('continuation_audit',{}).get('status')=='PASS' for k in ('linear','nonlinear')):
        base.finish_references_and_comparison(b,science,item,s)
        plot_bundle(b)
        s['one_d_authorization']=c['authorization'];s['one_d_inherited_runner_token_qualification']='legacy FEM3A callback boolean; actual new execution authorization is explicit_user_FEM3AR_2026_10_09'
        record_reference_authorization(b,c['authorization'])
    save(b,item,s);print(json.dumps({'bundle':str(b),'statuses':s['resume_statuses'],'job_calls':s['job_calls'],'numerical_seconds':s['numerical_seconds'],'overall':s['overall']},indent=2));return s


def plot_bundle(b):
    s=validate_cache(b)
    result=base.plot_bundle(b)
    if not (Path(b)/'dynamic_comparison.npz').exists():return result
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    T=s['preflight']['T1']
    fig,axes=plt.subplots(1,3,figsize=(11.5,3.3),constrained_layout=True)
    for kind in ('linear','nonlinear'):
        with np.load(Path(b)/('one_d_'+kind+'.npz')) as z:axes[0].plot(z['times']/T,z['energy_relative_drift'],label=kind)
        energies=read_json(Path(b)/'cases'/kind/'energy.json');rows=[r for r in energies['records'] if r.get('step')==2 and 'mechanical_energy' in r]
        rec=s['cases'][kind];E0=rec['native_initial_internal_energy'];Edyn=rec['native_dynamic_bookkeeping_initial_energy']
        times=[r['dynamic_time']/T for r in rows]
        axes[1].plot(times,[100*(r['mechanical_energy']/E0-1) for r in rows],label=kind)
        axes[2].plot(times,[r['mechanical_energy']/Edyn-1 for r in rows],label=kind)
    for ax in axes:ax.set(xlabel='t/T1',xlim=(0,.05));ax.grid(alpha=.25);ax.legend(fontsize=8)
    axes[0].set(ylabel='1D physical energy drift',title='Semidiscrete physical model')
    axes[1].set(ylabel='Native change vs STATIC (%)',ylim=(0,110),title='Native bookkeeping jump: +100%')
    axes[2].set(ylabel='Drift vs native DYNAMIC reference',title='Diagnostic only; energy PARTIAL')
    for suffix in ('png','pdf'):
        kwargs={'metadata':{'CreationDate':None,'ModDate':None}} if suffix=='pdf' else {}
        fig.savefig(Path(b)/'figures'/('energy_and_release_diagnostics.'+suffix),**kwargs)
    plt.close(fig)
    save(Path(b),read_json(Path(b)/'provenance.json'),s)
    return result

if __name__=='__main__':main()
