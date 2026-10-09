"""FEM-2R: explicitly authorized continuation of the stopped static comparison.

Thin orchestration over the existing FEM-2 generator/readers/recovery/reporting.
Parent1D equilibria and source meshes are immutable dependencies, never solved.
An authorization-bound attempt ledger survives code-hash changes; no retries.
"""
from __future__ import annotations
import argparse,contextlib,hashlib,importlib.metadata,json,math,os,shutil,sys,time
from pathlib import Path
if __name__=='__main__':
    for n in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[n]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as base
CONFIG=ROOT/'data/input/nlsp_nonlinear_static_3d_fem_resume.json'
OUTPUT=ROOT/'results/nlsp_nonlinear_static_3d_fem_resume'
ORDER=('medium_linear','medium_nonlinear','fine_linear','fine_nonlinear','refined_linear','refined_nonlinear')
STATUS_NAMES=('SOURCE_PRESERVATION','INPUT_SERIALIZATION','MEDIUM_LINEAR','MEDIUM_NONLINEAR','FINE_PAIR','REFINED_PAIR','STATIC_OUTPUT_RECOVERY','EQUILIBRIUM','NONLINEAR_SIGNAL_RESOLUTION','1D_3D_COMPARISON')
read_json,write_json,sha=base.read_json,base.write_json,base.sha
# Current proven parser is reused. A scoped reader repair may replace this only
# after genuine successful solver output is independently inspected.
parse_static_outputs=base.parse_static_outputs

def validate_resume_config(c):
    if c['schema']!='nlsp-static-fem2r-continuation-v1':raise ValueError('Continuation schema mismatch')
    auth=c['authorization']
    if auth['id']!='explicit_user_FEM2R_2026_10_09' or auth['maximum_new_ccx_jobs']!=6 or auth['case_order']!=list(ORDER) or auth['automatic_retry_after_new_solver_failure']:raise ValueError('Explicit bounded continuation authorization missing')
    if c['frozen_load']!={'g':.0014224751066856333,'q':2.844950213371267e-5,'F_total':2.844950213371267e-5}:raise ValueError('Frozen load changed')
    if c['geometry']!={'L':1.,'b':.2,'h':.1} or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}:raise ValueError('Frozen geometry/material changed')
    if c['one_d_policy']!='reuse_parent_p48_p64_without_equilibrium_solves':raise ValueError('1D reuse policy changed')
    if c['semantics']!={'new_static_jobs_only':True,'new_meshes':False,'new_modal_jobs':False,'new_dynamics':False,'new_one_d_solves':False,'model_fitting':False}:raise ValueError('Unauthorized computation scope')
    if c['threads']!=1 or not 0<c['job_timeout_seconds']<=1200 or not 0<c['numerical_budget_seconds']<=3600 or not 0<c['job_memory_limit_bytes']<=4*1024**3:raise ValueError('Resource policy mismatch')
    return c

def load_resume_parent(c):
    parent=ROOT/c['parent_failed_bundle']
    if sha(parent/'manifest.json')!=c['parent_manifest_sha256']:raise ValueError('Historical failed parent manifest mismatch')
    summary=base.validate_fem2_cache(parent)
    science=read_json(parent/'manifest.json')['identity']['config']
    base.validate_fem2_config(science);base.load_fem2_sources(science)
    if sha(base.__file__)!=c['corrected_generator_sha256']:raise ValueError('Corrected FEM2 helper changed')
    if summary['cases']['medium']['linear']['status']!='FAIL' or summary['job_calls']['ccx']!=1:raise ValueError('Wrong historical attempt')
    pre=summary['preflight']
    if pre['status']!='PASS' or {k:pre['load'][k] for k in c['frozen_load']}!=c['frozen_load']:raise ValueError('Saved1D/load gate mismatch')
    if science['geometry']!=c['geometry'] or science['material']!=c['material']:raise ValueError('Parent physical configuration differs')
    for p in (48,64):
        with np.load(parent/f'one_d_p{p}.npz',allow_pickle=False) as data:
            if data['linear'].shape!=(1001,4) or data['nonlinear'].shape!=(1001,4):raise ValueError('Incomplete saved1D fields')
            if not np.all(np.isfinite(data['q_linear'])) or not np.all(np.isfinite(data['q_nonlinear'])):raise ValueError('Nonfinite saved1D coefficients')
            if np.max(abs(data['linear'][[0,-1]])) or np.max(abs(data['nonlinear'][[0,-1]])):raise ValueError('Parent essential endpoint failure')
            middle=np.flatnonzero(data['s']==.5)
            if len(middle)!=1 or abs(data['linear'][middle[0],1]-.005)>1e-14 or abs(data['nonlinear'][middle[0],1]-.004991129920693126)>1e-14:raise ValueError('Saved1D values not reproduced')
    return parent,summary,science

def identity(config_path=CONFIG):
    c=validate_resume_config(read_json(config_path));parent,s,science=load_resume_parent(c)
    item={'schema':c['schema'],'config':science,'continuation_config':c,'config_sha256':sha(config_path),'code_sha256':sha(__file__),
          'parent_bundle':c['parent_failed_bundle'],'parent_manifest_sha256':sha(parent/'manifest.json'),
          'source_1D_sha256':{f'one_d_p{p}.npz':sha(parent/f'one_d_p{p}.npz') for p in (48,64)},
          'generator_sha256':sha(base.__file__),'sources':science['sources'],'authorization':c['authorization'],
          'ccx_sha256':sha(science['ccx_exe']),'runtime_dlls':{p.name:sha(p) for p in sorted(Path(science['ccx_exe']).parent.glob('*.dll'))},
          'python':sys.version,'dependencies':{n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item

def numeric_cards(text):
    lines=text.splitlines();out={}
    for keyword in ('STATIC','CONTROLS','ELASTIC','DENSITY','DLOAD'):
        idx=next(i for i,l in enumerate(lines) if l.strip().upper().startswith('*'+keyword))
        tokens=[v.strip() for v in lines[idx+1].split(',')]
        if keyword=='DLOAD':tokens=tokens[2:]
        tokens=[v for v in tokens if v]
        if any(len(v)>20 or not math.isfinite(float(v)) for v in tokens):raise ValueError('Invalid native20-character numeric '+keyword)
        out[keyword]=tokens
    return out

def input_gate(c,bundle,parent_summary):
    science=read_json(ROOT/c['parent_failed_bundle']/'manifest.json')['identity']['config']
    sources=base.load_fem2_sources(science);target=Path(bundle)/'input_gate';target.mkdir(parents=True,exist_ok=True);rows=[]
    for level in science['mesh_levels']:
        source,mesh,audit=base.fem2_source_mesh(science,level,sources);paths=[]
        for kind,nlgeom in (('linear',False),('nonlinear',True)):
            p=target/(level+'_'+kind+'.inp')
            base.write_static_input(p,source/'solid_mesh.inp',mesh,audit,science['material'],c['frozen_load']['g'],nlgeom,science['static_settings'])
            tokens=numeric_cards(p.read_text());paths.append(p)
            g=float(tokens['DLOAD'][0])
            if abs(g/c['frozen_load']['g']-1)>1e-12:raise ValueError('Gravity serialization lost significant precision')
            rows.append({'case':level+'_'+kind,'maximum_numeric_token_width':max(len(v) for vv in tokens.values() for v in vv),'cards':tokens,'serialized_g':g,'relative_g_rounding':g/c['frozen_load']['g']-1,'mesh_include_sha256':sha(source/'solid_mesh.inp')})
        if paths[1].read_text().replace(', NLGEOM','')!=paths[0].read_text():raise ValueError('Linear/NL physical deck mismatch')
    def ignore_include(text):return '\n'.join(l for l in text.splitlines() if not l.startswith('*INCLUDE'))
    preview=ROOT/c['parent_failed_bundle']/'remediation_preview/medium_linear_corrected.inp'
    if ignore_include(preview.read_text())!=ignore_include((target/'medium_linear.inp').read_text()):raise ValueError('Generator differs from corrected preview beyond mesh path')
    failed=ROOT/c['parent_failed_bundle']/'cases/medium/linear/static.inp';lines=failed.read_text().splitlines()
    old=lines[lines.index('*STATIC')+1].split(',')
    result={'status':'PASS','rows':rows,'old_failed_STATIC_widths':[len(v.strip()) for v in old],'preview_agrees_except_include':True,'linear_NL_only_NLGEOM_difference':True,'new_solver_calls':0}
    write_json(Path(bundle)/'input_gate.json',result);return result

def validate_cache(b,item=None):
    b=Path(b);m=read_json(b/'manifest.json')
    if item is not None and m['identity']!=item:raise ValueError('Resume cache identity mismatch')
    for p,d in m['artifact_hashes'].items():
        if sha(b/p)!=d:raise ValueError('Resume artifact changed: '+p)
    load_resume_parent(m['identity']['continuation_config'])
    return read_json(b/'summary.json')

def find_authorized_attempt(c,output_dir):
    # One authorization cannot acquire another attempt through a changed code hash.
    for directory in {Path(output_dir).resolve(),OUTPUT.resolve()}:
        if not directory.is_dir():continue
        for b in sorted(directory.iterdir()):
            if not b.is_dir() or not (b/'provenance.json').exists():continue
            item=read_json(b/'provenance.json')
            if item.get('authorization',{}).get('id')!=c['authorization']['id'] or item.get('parent_manifest_sha256')!=c['parent_manifest_sha256']:continue
            if item['continuation_config']!=c:raise ValueError('Authorization reused with different config')
            if not (b/'manifest.json').exists():raise RuntimeError('Unmanifested interrupted continuation exists; no new solver attempt')
            return b,validate_cache(b)
    return None

@contextlib.contextmanager
def parent_profile_loader(parent):
    old=base.load_one_d_static_profile
    base.load_one_d_static_profile=lambda bundle,p=64,nonlinear=False:old(parent,p,nonlinear)
    try:yield
    finally:base.load_one_d_static_profile=old

def update_resume_statuses(s):
    for label in STATUS_NAMES:s['resume_statuses'].setdefault('NLSP_FEM2R_'+label,'NOT_RUN')
    for level,kind,name in [('medium','linear','MEDIUM_LINEAR'),('medium','nonlinear','MEDIUM_NONLINEAR')]:
        r=s['cases'].get(level,{}).get(kind,{})
        s['resume_statuses']['NLSP_FEM2R_'+name]='PASS' if r.get('status')=='PASS' else 'FAIL' if r.get('status')=='FAILED_SOLVER' or r.get('status')=='FAIL' else 'PARTIAL' if r else 'NOT_RUN'
    for level,name in [('fine','FINE_PAIR'),('refined','REFINED_PAIR')]:
        rows=s['cases'].get(level,{})
        s['resume_statuses']['NLSP_FEM2R_'+name]='PASS' if all(rows.get(k,{}).get('status')=='PASS' for k in ('linear','nonlinear')) else 'FAIL' if any(v.get('status') in ('FAIL','FAILED_SOLVER') for v in rows.values()) else 'PARTIAL' if rows else 'NOT_RUN'
    all_done=s.get('completed_levels')==['medium','fine','refined']
    failures=any(v.get('status') in ('FAIL','FAILED_SOLVER') for lvl in s['cases'].values() for v in lvl.values())
    for name in ('STATIC_OUTPUT_RECOVERY','EQUILIBRIUM'):
        s['resume_statuses']['NLSP_FEM2R_'+name]='PASS' if all_done else 'FAIL' if failures else 'PARTIAL' if s['cases'] else 'NOT_RUN'
    for name,old in [('NONLINEAR_SIGNAL_RESOLUTION','NLSP_FEM2_NONLINEAR_CORRECTION_MESH_CHECK'),('1D_3D_COMPARISON','NLSP_FEM2_1D_3D_COMPARISON')]:
        s['resume_statuses']['NLSP_FEM2R_'+name]=s['statuses'].get(old,'NOT_RUN')
    s['overall']='FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS' if all_done else 'PARTIAL'

def save_progress(b,s):
    update_resume_statuses(s);write_json(b/'summary.json',s)
    write_json(b/'attempt_ledger.json',{'authorization':s['authorization'],'attempts':s['attempt_ledger'],'new_ccx_jobs':s['job_calls']['ccx'],'maximum':6})

def recover_case(c,b,s,level,kind,mesh,audit):
    record=s['cases'][level][kind];case=b/'cases'/level/kind;start=time.perf_counter();nonlinear=kind=='nonlinear'
    try:
        diag,arrays=parse_static_outputs(case/'static.inp',mesh,audit,1.,c['frozen_load']['g'],nonlinear,s['science_config']['gates']['equilibrium_relative'])
        write_json(case/'static_diagnostics.json',diag);np.savez_compressed(case/'static_nodal_results.npz',**arrays);record['diagnostics']=diag
        if diag['status']!='PASS':
            record.update(status='FAIL',failure='Static output/equilibrium gate: '+str(diag['failures']));return False
        recovered=base.recover_static_sections(mesh,arrays['U'],s['preflight']['coefficients'],41,1.,True)
        write_json(case/'recovered_sections.json',recovered)
        record.update(status='PASS',recovery_policy=base.FEM2_RECOVERY_POLICY,recovery_code_sha256=sha(__file__))
        return True
    except Exception as exc:
        # Genuine finished output is retained; this is not accepted equilibrium.
        record.update(status='RECOVERY_PENDING',failure=str(exc),read_only_reparse_permitted=True)
        return False
    finally:
        record['recovery_seconds']=time.perf_counter()-start;s['runtime']['numerical_seconds']+=record['recovery_seconds'];write_json(case/'case.json',record);save_progress(b,s)

def run_resume(c,b,s,parent,through_case='refined_nonlinear'):
    sources=base.load_fem2_sources(s['science_config']);wanted=ORDER[:ORDER.index(through_case)+1]
    for name in wanted:
        level,kind=name.split('_');nonlinear=kind=='nonlinear'
        if level!='medium' and not all(s['cases'].get('medium',{}).get(k,{}).get('status')=='PASS' for k in ('linear','nonlinear')):break
        if kind=='nonlinear' and s['cases'].get(level,{}).get('linear',{}).get('status')!='PASS':break
        source,mesh,audit=base.fem2_source_mesh(s['science_config'],level,sources);s['cases'].setdefault(level,{})
        record=s['cases'][level].get(kind)
        if record:
            if record['status']=='PASS':continue
            if record['status']!='RECOVERY_PENDING':break
            if not recover_case(c,b,s,level,kind,mesh,audit):break
        else:
            case=b/'cases'/level/kind
            if case.exists() and any(case.iterdir()):raise RuntimeError('Attempt directory exists without accepted ledger; no retry')
            case.mkdir(parents=True,exist_ok=True)
            record={'status':'RUNNING','source_mesh':str(source.relative_to(ROOT)),'source_mesh_include_sha256':sha(source/'solid_mesh.inp'),'mesh_audit':audit,'execution_code_sha256':sha(__file__)}
            s['cases'][level][kind]=record
            attempt={'case':name,'ordinal':len(s['attempt_ledger'])+1,'status':'STARTED','code_sha256':sha(__file__)}
            s['attempt_ledger'].append(attempt);save_progress(b,s);start=time.perf_counter()
            try:
                write_json(case/'input_contract.json',base.write_static_input(case/'static.inp',source/'solid_mesh.inp',mesh,audit,s['science_config']['material'],c['frozen_load']['g'],nonlinear,s['science_config']['static_settings']))
                numeric_cards((case/'static.inp').read_text())
                if s['job_calls']['ccx']>=6:raise RuntimeError('Six new CCX attempts exhausted')
                remaining=c['numerical_budget_seconds']-s['runtime']['numerical_seconds']
                if remaining<=0:raise TimeoutError('Total numerical budget exhausted')
                env=dict(os.environ);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBER_OF_CPUS='1')
                print('FEM2R new static job: '+name,flush=True);s['job_calls']['ccx']+=1;save_progress(b,s)
                result,stats=base.fem1.run_job([s['science_config']['ccx_exe'],'static'],case,min(c['job_timeout_seconds'],remaining),c['job_memory_limit_bytes'],case/'static',env)
                record['job']=stats;write_json(case/'job.json',stats)
                stdout=case/'static.stdout.txt';record['solver_finished']=result.returncode==0 and not stats['failure'] and 'JOB FINISHED' in stdout.read_text(encoding='utf8',errors='replace').upper()
                if not record['solver_finished']:raise RuntimeError('No successful real solver completion: '+str(stats))
                record['status']='RECOVERY_PENDING';attempt['status']='SOLVER_FINISHED'
            except Exception as exc:
                record.update(status='FAILED_SOLVER',failure=str(exc));attempt.update(status='FAILED_SOLVER',failure=str(exc))
            finally:
                record['execution_stage_seconds']=time.perf_counter()-start;s['runtime']['numerical_seconds']+=record['execution_stage_seconds'];write_json(case/'case.json',record);save_progress(b,s)
            if record['status']=='FAILED_SOLVER':break
            if not recover_case(c,b,s,level,kind,mesh,audit):break
        s['attempt_ledger'][ORDER.index(name)]['status']='PASS'
        if kind=='nonlinear':
            with parent_profile_loader(parent):base.fem2_update_comparison(s['science_config'],b,s)
            if level=='medium':write_json(b/'medium_gate.json',{'status':'PASS','both_medium_solutions_accepted':True,'load_and_recovery_checked':True})
        save_progress(b,s)
    return s

def plot_only(b):
    s=validate_cache(b)
    if not s.get('completed_levels'):return {'figures':0,'new_solver_calls':0}
    return base.fem2_plot_only(b)

def finalize_manifest(b,item):
    write_json(b/'manifest.json',base.fem1.artifact_manifest(b,item))

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check-source',action='store_true');mode.add_argument('--run-fem',action='store_true');mode.add_argument('--reparse-only',action='store_true')
    mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    p.add_argument('--config',type=Path,default=CONFIG);p.add_argument('--output-dir',type=Path,default=OUTPUT);p.add_argument('--through-case',choices=ORDER,default='refined_nonlinear')
    a=p.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_cache(b)
        if a.plot_only:plot_only(b)
        print(json.dumps({'bundle':str(b),'statuses':s['resume_statuses'],'new_solver_static_calls':0},indent=2));return s
    c=validate_resume_config(read_json(a.config));parent,old,science=load_resume_parent(c)
    existing=find_authorized_attempt(c,a.output_dir)
    if existing:
        b,s=existing;item=read_json(b/'provenance.json')
        if a.check_source or s.get('completed_levels')==['medium','fine','refined'] or any(v.get('status') in ('FAIL','FAILED_SOLVER') for lvl in s['cases'].values() for v in lvl.values()):
            print(json.dumps({'bundle':str(b),'authorization_cache_hit':True,'statuses':s['resume_statuses'],'new_solver_static_calls':0},indent=2));return s
        # Permit code evolution ONLY for read-only parsing of already completed data.
        if sha(__file__)!=item['code_sha256']:
            history=b/'recovery_code';history.mkdir(exist_ok=True);shutil.copyfile(__file__,history/(sha(__file__)+'.py'))
            write_json(b/'reparse_provenance.json',{'original_execution_code_sha256':item['code_sha256'],'current_recovery_code_sha256':sha(__file__),'no_solver_retry':True})
    else:
        key,item=identity(a.config);b=a.output_dir/key
        if b.exists() and any(b.iterdir()):raise RuntimeError('Unmanifested continuation; no retry')
        b.mkdir(parents=True,exist_ok=True);write_json(b/'provenance.json',item);write_json(b/'frozen_config.json',{'continuation':c,'science':science,'selected_load':old['preflight']['load']})
        (b/'execution_code').mkdir();shutil.copyfile(__file__,b/'execution_code'/Path(__file__).name)
        gate=input_gate(c,b,old)
        s={'preflight':old['preflight'],'science_config':science,'authorization':c['authorization'],'cases':{},'attempt_ledger':[],
           'job_calls':{'ccx':0,'gmsh':0,'modal':0,'nonlinear_ODE':0,'one_d_static':0},'runtime':{'numerical_seconds':0.,'budget_seconds':3600},
           'statuses':{'NLSP_FEM2_'+n:'NOT_RUN' for n in base.FEM2_STATUS_NAMES},
           'resume_statuses':{'NLSP_FEM2R_'+n:'NOT_RUN' for n in STATUS_NAMES},'parent_reference':{'path':c['parent_failed_bundle'],'manifest_sha256':c['parent_manifest_sha256'],'one_d_resolved':False}}
        s['resume_statuses'].update(NLSP_FEM2R_SOURCE_PRESERVATION='PASS',NLSP_FEM2R_INPUT_SERIALIZATION=gate['status'])
        s['statuses'].update(NLSP_FEM2_LOAD_PREFLIGHT='PASS',NLSP_FEM2_1D_LINEAR_STATIC='PASS',NLSP_FEM2_1D_NONLINEAR_STATIC='PASS')
        save_progress(b,s);finalize_manifest(b,item)
    if not a.check_source:
        if a.reparse_only:
            sources=base.load_fem2_sources(science)
            for level,rows in s['cases'].items():
                for kind,r in rows.items():
                    if r['status']=='RECOVERY_PENDING':
                        _,mesh,audit=base.fem2_source_mesh(science,level,sources);recover_case(c,b,s,level,kind,mesh,audit)
            for level in science['mesh_levels']:
                if all(s['cases'].get(level,{}).get(k,{}).get('status')=='PASS' for k in ('linear','nonlinear')):
                    with parent_profile_loader(parent):base.fem2_update_comparison(science,b,s)
            save_progress(b,s)
        else:s=run_resume(c,b,s,parent,a.through_case)
    finalize_manifest(b,item)
    if s.get('completed_levels'):plot_only(b);finalize_manifest(b,item)
    print(json.dumps({'bundle':str(b),'statuses':s['resume_statuses'],'runtime':s['runtime'],'job_calls':s['job_calls'],'overall':s['overall']},indent=2));return s

if __name__=='__main__':main()

