"""FEM-3A bounded static preload -> instantaneous release -> free-motion pilot.

A new transient I/O contract, not a parameter variant of historical static CLI.
Reuse frozen FEM2 generation/recovery, native monitoring and full planar Radau.
"""
from __future__ import annotations
import argparse,hashlib,importlib.metadata,json,math,os,re,shutil,subprocess,sys,time
from pathlib import Path
if __name__=='__main__':
    for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
from scripts.analysis import resume_nlsp_nonlinear_static_3d_fem as resume
from scripts.lib import nlsp_fem3a_1d_reference as one
from scripts.lib import nlsp_fem3a_transient_output as io
base=resume.base
read_json,write_json,sha=base.read_json,base.write_json,base.sha
CONFIG=ROOT/'data/input/nlsp_nonlinear_dynamic_3d_fem_pilot.json'
OUTPUT=ROOT/'results/nlsp_nonlinear_dynamic_3d_fem_pilot'
STATUS_NAMES=('SOURCE_PRESERVATION','LOAD_RELEASE_PROTOCOL','STATIC_STATE_TRANSFER','INITIAL_VELOCITIES','1D_DYNAMIC_REFERENCE','3D_LINEAR_DYNAMIC','3D_NONLINEAR_DYNAMIC','TRANSIENT_FIELD_RECOVERY','ENERGY_DIAGNOSTICS','SHORT_RESPONSE_COMPARISON')


def validate_config(c):
    if c['schema']!='nlsp-fem3a-dynamic-pilot-v1':raise ValueError('Wrong pilot schema')
    if c['authorization']!={'id':'explicit_user_FEM3A_2026_10_09','maximum_production_CCX_jobs':2,'maximum_nonlinear_1D_integrations':1,'maximum_technical_fixtures':1,'automatic_retry':False}:raise ValueError('Explicit bounded authorization missing')
    if c['geometry']!={'L':1.,'b':.2,'h':.1} or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}:raise ValueError('Frozen geometry/material changed')
    if c['g']!=.0014224751066856333 or c['q']!=2.844950213371267e-5 or c['omega1']!=.6054167303477958 or c['horizon_T1']!=.05 or c['mesh_level']!='medium':raise ValueError('Frozen load/time/mesh scope changed')
    if c['dynamic']!={'alpha':0,'initial_T1_fraction':1/4000,'maximum_T1_fraction':1/2000,'minimum_initial_fraction':1e-4,'maximum_increments':2000,'output_frequency':1,'release':'OP=NEW plus zero GRAV; STEP AMPLITUDE=STEP'}:raise ValueError('Preselected dynamic policy changed')
    if c['threads']!=1 or c['job_timeout_seconds']!=1200 or c['job_memory_limit_bytes']!=4*1024**3 or c['numerical_budget_seconds']!=3600:raise ValueError('Bounded resource policy changed')
    if c['execution_mode']!='EXPLORATORY_NOT_CERTIFIED' or c['admitted'] is not False:raise ValueError('Strict qualification must stay separate')
    return c


def verify_sources(c):
    for key in ('source_resume','source_static','source_fem1','source_action'):
        source=c[key];b=ROOT/source['bundle']
        if sha(b/'manifest.json')!=source['manifest_sha256']:raise ValueError('Source manifest mismatch: '+key)
        manifest=read_json(b/'manifest.json');artifacts=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
        if not artifacts:raise ValueError('Missing source artifact hash contract: '+key)
        for p,h in artifacts.items():
            if sha(b/p)!=h:raise ValueError('Source artifact mismatch: '+key+'/'+p)
    s=resume.validate_cache(ROOT/c['source_resume']['bundle'])
    if s['completed_levels']!=['medium','fine','refined']:raise ValueError('Incomplete static source')
    sources=base.load_fem2_sources(s['science_config'])
    source,mesh,audit=base.fem2_source_mesh(s['science_config'],'medium',sources)
    if audit['status']!='PASS' or len(mesh.nodes)!=5649 or len(mesh.solid_elements)!=3120:raise ValueError('Saved medium mesh mismatch')
    if not Path(s['science_config']['ccx_exe']).is_file():raise ValueError('Missing existing CCX binary')
    return s,source,mesh,audit


def identity(config_path=CONFIG):
    c=validate_config(read_json(config_path));s,source,_,_=verify_sources(c)
    helpers=[Path(__file__),Path(one.__file__),Path(io.__file__),Path(base.__file__),Path(base.fem1.__file__),Path(one.dynamics.__file__),Path(one.rod.__file__),Path(one.runner.__file__)]
    item={'config':c,'config_sha256':sha(config_path),'helper_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in helpers},'source_mesh':source.relative_to(ROOT).as_posix(),'source_mesh_sha256':sha(source/'solid_mesh.inp'),'solver_sha256':sha(s['science_config']['ccx_exe']),'runtime_DLLs':{p.name:sha(p) for p in sorted(Path(s['science_config']['ccx_exe']).parent.glob('*.dll'))},'protocol_sha256':sha(ROOT/'results/_smoke/fem3a_protocol/protocol_evidence.json'),'python':sys.version,'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','scipy','matplotlib')},'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def dynamic_settings(c):
    T=2*math.pi/c['omega1'];d=c['dynamic'];initial=T*d['initial_T1_fraction']
    return {'T1':T,'initial_increment':initial,'duration':T*c['horizon_T1'],'minimum_increment':initial*d['minimum_initial_fraction'],'maximum_increment':T*d['maximum_T1_fraction'],'maximum_increments':d['maximum_increments'],'alpha':0}


def write_input(path,c,static_summary,source,mesh,audit,nonlinear):
    settings=dict(static_summary['science_config']['static_settings']);settings['output_frequency']=1
    base.write_static_input(path,source/'solid_mesh.inp',mesh,audit,c['material'],c['g'],nonlinear,settings)
    text=Path(path).read_text(encoding='utf8');step=text.index('*STEP')
    text=text[:step]+'*INITIAL CONDITIONS, TYPE=VELOCITY\nALL_NODES,1,0.\nALL_NODES,2,0.\nALL_NODES,3,0.\n'+text[step:]
    text=text.replace('no modal, dynamic or mid-span joint.','no modal or mid-span joint; new direct dynamic step follows.')
    text=text.replace('S,E\n','S,E,ENER\n')
    # Native2.22 frees veold before LINEAR STATIC; ELKE would dereference it.
    static_energy='ELSE,ELKE' if nonlinear else 'ELSE'
    text=text.replace('*END STEP\n','*EL PRINT, ELSET=SOLID, TOTALS=ONLY, FREQUENCY=1\n'+static_energy+'\n*END STEP\n',1)
    d=dynamic_settings(c);fmt=base.fem1.single.ccx_float
    lines=['** Free motion: explicitly remove the entire previous body load at step start.',f"*STEP, INC={d['maximum_increments']}, AMPLITUDE=STEP, "+('NLGEOM' if nonlinear else 'NLGEOM=NO'),'*DYNAMIC, ALPHA=0',','.join(fmt(d[k]) for k in ('initial_increment','duration','minimum_increment','maximum_increment')),'*DLOAD, OP=NEW','SOLID,GRAV,0.,0.,-1.,0.','*NODE FILE, GLOBAL=YES, FREQUENCY=1','U,V,RF','*NODE PRINT, NSET=ALL_NODES, GLOBAL=YES, FREQUENCY=1','U','*NODE PRINT, NSET=LEFT_FIXED, TOTALS=YES, GLOBAL=YES, FREQUENCY=1','RF','*NODE PRINT, NSET=RIGHT_FIXED, TOTALS=YES, GLOBAL=YES, FREQUENCY=1','RF','*EL PRINT, ELSET=SOLID, TOTALS=ONLY, FREQUENCY=1','ELSE,ELKE','*END STEP']
    text+='\n'.join(lines)+'\n';Path(path).write_text(text,encoding='utf8')
    return input_contract(text,c)


def input_contract(text,c):
    resume.numeric_cards(text)
    first=text.split('*END STEP',1)[0];first_step=next(v for v in first.splitlines() if v.startswith('*STEP'))
    if 'NLGEOM' not in first_step and 'ELKE' in first:raise ValueError('Native2.22 ELKE is unsafe in linear STATIC; request kinetic energy in DYNAMIC only')
    lines=text.splitlines();idx=next(i for i,v in enumerate(lines) if v.upper().startswith('*DYNAMIC'))
    tokens=[v.strip() for v in lines[idx+1].split(',')]
    if any(len(v)>20 or not math.isfinite(float(v)) for v in tokens):raise ValueError('Invalid native dynamic numeric field')
    d=dynamic_settings(c)
    for token,key in zip(tokens,('initial_increment','duration','minimum_increment','maximum_increment')):
        if abs(float(token)/d[key]-1)>1e-12:raise ValueError('Dynamic serialization changed value')
    dynamic='\n'.join(lines[idx-1:])
    if '*DYNAMIC, ALPHA=0' not in dynamic or 'AMPLITUDE=STEP' not in dynamic or '*DLOAD, OP=NEW\nSOLID,GRAV,0.,0.,-1.,0.' not in dynamic:raise ValueError('Instantaneous release/alpha contract failed')
    if any(word in text.upper() for word in ('*DAMPING','*CONTACT','*SPRING','*MPC','*FREQUENCY','*MODAL DYNAMIC','*CLOAD')):raise ValueError('Unauthorized extra physics')
    if '*INITIAL CONDITIONS, TYPE=VELOCITY\nALL_NODES,1,0.\nALL_NODES,2,0.\nALL_NODES,3,0.' not in text:raise ValueError('Zero initial velocities missing')
    return {'status':'PASS','dynamic_settings':d,'dynamic_native_numeric_fields':tokens,'maximum_numeric_width':max(map(len,tokens)),'external_GRAV_after_release':0.,'alpha':0,'zero_velocity_initial_condition':True,'no_new_dynamic_derivative_constraints':True,'preload_physics_unchanged':True}


def validate_cache(b):
    m=read_json(Path(b)/'manifest.json')
    for p,h in m['artifact_hashes'].items():
        if sha(Path(b)/p)!=h:raise ValueError('Pilot artifact hash mismatch: '+p)
    verify_sources(m['identity']['config'])
    return read_json(Path(b)/'summary.json')


def finalize(b,item,s):
    write_json(b/'summary.json',s);write_json(b/'attempt_ledger.json',{'authorization':s['authorization'],'attempts':s['attempts'],'job_calls':s['job_calls'],'automatic_retry':False})
    write_json(b/'manifest.json',base.fem1.artifact_manifest(b,item))


def existing_attempt(c):
    if not OUTPUT.exists():return None
    for b in sorted(OUTPUT.iterdir()):
        if not b.is_dir() or not (b/'provenance.json').exists():continue
        item=read_json(b/'provenance.json')
        if item['config']['authorization']['id']==c['authorization']['id']:
            if item['config']!=c:raise ValueError('Authorization already used with different config')
            if not (b/'manifest.json').exists():raise RuntimeError('Interrupted unmanifested attempt; no new solver job')
            return b,item,validate_cache(b)
    return None


def prepare(config_path=CONFIG):
    c=validate_config(read_json(config_path));key,item=identity(config_path);b=OUTPUT/key
    if b.exists() and any(b.iterdir()):raise RuntimeError('Existing attempt cannot be overwritten')
    b.mkdir(parents=True);write_json(b/'provenance.json',item);write_json(b/'config.json',c)
    shutil.copyfile(ROOT/'results/_smoke/fem3a_protocol/protocol_evidence.json',b/'protocol_evidence.json')
    execution=b/'execution_code';execution.mkdir()
    for p in (Path(__file__),Path(one.__file__),Path(io.__file__)):shutil.copyfile(p,execution/p.name)
    old,source,mesh,audit=verify_sources(c)
    ref=one.load_reference(ROOT,c['source_static']['bundle'],c['source_fem1']['bundle'],c['source_action']['bundle'])
    pre=one.preflight_reference(ref,c);write_json(b/'one_d_preflight.json',pre)
    (b/'input_gate').mkdir()
    gates={kind:write_input(b/'input_gate'/f'{kind}.inp',c,old,source,mesh,audit,kind=='nonlinear') for kind in ('linear','nonlinear')}
    write_json(b/'input_gate.json',gates)
    s={'authorization':c['authorization'],'statuses':{'NLSP_FEM3A_'+n:'NOT_RUN' for n in STATUS_NAMES},'cases':{},'attempts':[],'job_calls':{'CCX_production':0,'CCX_fixture':0,'Gmsh':0,'1D_nonlinear_ODE':0,'1D_static':0,'physical_root_search':0,'symbolic_derivations':0},'numerical_seconds':0.,'overall':'NOT_RUN','preflight':{'status':'PASS','omega1':ref['omega1'],'T1':ref['T1'],'target_dynamic_end':.05*ref['T1']},'strict_float64_qualification':'PARTIAL','execution_mode':'EXPLORATORY_NOT_CERTIFIED','admitted':False}
    s['statuses'].update(NLSP_FEM3A_SOURCE_PRESERVATION='PASS',NLSP_FEM3A_LOAD_RELEASE_PROTOCOL='PASS',NLSP_FEM3A_INITIAL_VELOCITIES='PASS');finalize(b,item,s)
    return b,item,s


def _rounded_endpoint(time_value,end):
    if time_value>end and time_value-end<=1e-8:return end
    return time_value


def recover_case(b,c,s,kind,old,mesh,audit):
    case=b/'cases'/kind;start=time.perf_counter();ids,xyz,_,_=base.fem1.mesh_arrays(mesh)
    sta=io.read_transient_sta(case/'motion.sta');write_json(case/'increments.json',sta)
    static_rows=[r for r in sta['accepted_increments'] if r['step']==1]
    if not static_rows or abs(static_rows[-1]['step_time']-1)>1e-6:raise ValueError('Missing full-load static preload')
    static_end=static_rows[-1]['total_time'];node_sets={'ALL_NODES':ids,'LEFT_FIXED':audit['fixed_left_ids'],'RIGHT_FIXED':audit['fixed_right_ids']}
    datdir=case/'dat_fields';datdir.mkdir(exist_ok=True);datmeta={}
    for block in io.iter_transient_dat(case/'motion.dat',node_sets,static_end_time=static_end,increments=sta['accepted_increments']):
        key=(block['step'],block['increment'],block['set'],block['name']);dest=datdir/('_'.join(map(str,key))+'.npz')
        np.savez_compressed(dest,values=block['values']);datmeta[key]={'path':dest,'time':block['total_time']}
    quad=base.fem1.quadrature_arrays(mesh,1.);coefficients=old['preflight']['coefficients'];fixed=np.r_[audit['fixed_left_ids'],audit['fixed_right_ids']]
    framesdir=case/'frames';framesdir.mkdir(exist_ok=True);rows=[];static_final=None;max_round=0.;end=dynamic_settings(c)['duration']
    for frame in io.iter_transient_frd(case/'motion.frd',ids,static_end_time=static_end,fixed_node_ids=fixed,increments=sta['accepted_increments']):
        key=(frame['step'],frame['increment'],'ALL_NODES','DISP')
        if key not in datmeta:raise ValueError('Missing corresponding complete DAT displacement block')
        with np.load(datmeta[key]['path']) as data:U=data['values'].copy()
        difference=float(np.max(abs(U-frame['fields']['DISP'])));max_round=max(max_round,difference)
        if difference>1e-8 or frame['fixed_displacement_max']>1e-12 or (frame.get('fixed_velocity_max') or 0)>1e-12:raise ValueError('DAT/FRD or clamp gate failed')
        meta={k:v for k,v in frame.items() if k not in ('fields',)}
        payload={'node_ids':ids,'U':U,**{k:v for k,v in frame['fields'].items() if k!='DISP'}}
        dest=framesdir/f"step{frame['step']}_inc{frame['increment']:05d}.npz";np.savez_compressed(dest,**payload)
        if frame['step']==1:
            static_final=(frame,U,meta);continue
        if frame['step']!=2:raise ValueError('Unexpected dynamic step')
        disp=base.fem1.nlsp_evaluate_tet10_displacements(U,quad['conn'],quad['N'])
        profile=base.fem2_recover_reference_samples(quad['xyz'],disp,quad['weights'],1.,.1,.2,41)
        # Translation velocities use the same reference-section projection. The
        # finite-rotation coordinate derivative is not the polar recovery of V.
        velocity=base.fem1.nlsp_evaluate_tet10_displacements(frame['fields']['VELO'],quad['conn'],quad['N'])
        vprofile=base.fem2_recover_reference_samples(quad['xyz'],velocity,quad['weights'],1.,.1,.2,41)
        report_x=np.linspace(0.,1.,41);fields=base.fem2_static_sample(profile,report_x);translations=base.fem2_static_sample(vprofile,report_x)[:,:3]
        rows.append({'time':_rounded_endpoint(frame['dynamic_time'],end),'printed_dynamic_time':frame['dynamic_time'],'step':2,'increment':frame['increment'],'x':report_x,'fields':fields,'translation_velocities':translations,'frame':str(dest.relative_to(b)),'native_metadata':meta})
    if static_final is None or not rows:raise ValueError('Missing preload or actual dynamic frames')
    frame,U,meta=static_final
    sourcecase=ROOT/c['source_resume']['bundle']/'cases/medium'/kind
    with np.load(sourcecase/'static_nodal_results.npz') as data:
        delta=float(np.max(abs(U-data['U'])));stress=float(np.max(abs(frame['fields']['STRESS']-data['stress']))) if 'stress' in data else None
        arrays={k:data[k].copy() for k in data.files}
    # The same pointwise native output precision governs the repeated preload.
    rf=[]
    for name in ('LEFT_FIXED','RIGHT_FIXED'):
        key=(1,frame['increment'],name,'FORC')
        if key not in datmeta:raise ValueError('Missing actual preload reactions')
        with np.load(datmeta[key]['path']) as data:rf.append(data['values'].copy())
    rf_delta=float(np.max(abs(np.vstack(rf)-arrays['RF_support_DAT'])))
    stress_key='S'
    strain_key='E'
    stress_delta=0.;strain_delta=0.
    if stress_key:stress_delta=float(np.max(abs(frame['fields']['STRESS']-arrays[stress_key])))
    if strain_key:strain_delta=float(np.max(abs(frame['fields']['TOSTRAIN']-arrays[strain_key])))
    if delta>c['preload_reproduction']['U_absolute'] or rf_delta>c['preload_reproduction']['RF_absolute'] or stress_delta>c['preload_reproduction']['relative_S_E']*np.max(abs(arrays['S'])) or strain_delta>c['preload_reproduction']['relative_S_E']*np.max(abs(arrays['E'])):raise ValueError('Repeated static preload differs unexpectedly')
    static_disp=base.fem1.nlsp_evaluate_tet10_displacements(U,quad['conn'],quad['N'])
    static_profile=base.fem2_recover_reference_samples(quad['xyz'],static_disp,quad['weights'],1.,.1,.2,41)
    original=read_json(sourcecase/'recovered_sections.json');profiledelta=float(np.max(abs(np.asarray(static_profile['fields'])-np.asarray(original['fields']))))
    preload={'status':'PASS','node_displacement_max_difference':delta,'support_RF_max_difference':rf_delta,'section_profile_max_difference':profiledelta,'stress_difference':stress_delta,'strain_difference':strain_delta,'static_end_total_time':static_end,'actual_metadata':meta,'velocity_initialization':'explicit zero IC; source-zeroed static-to-dynamic state; no static frame relabelled as dynamic t0','zero_velocity_proof':'protocol_evidence.json'}
    write_json(case/'preload_transfer.json',preload)
    static_fields=base.fem2_static_sample(static_profile,np.linspace(0.,1.,41))
    np.savez_compressed(case/'initial_sections.npz',x=np.linspace(0.,1.,41),fields=static_fields,source_static_end_time=static_end,not_a_native_dynamic_zero_frame=np.array(True))
    np.savez_compressed(case/'section_history.npz',time=np.array([r['time'] for r in rows]),printed_dynamic_time=np.array([r['printed_dynamic_time'] for r in rows]),x=rows[0]['x'],fields=np.stack([r['fields'] for r in rows]),translation_velocities=np.stack([r['translation_velocities'] for r in rows]),increments=np.array([r['increment'] for r in rows]))
    write_json(case/'frame_metadata.json',[{k:v for k,v in r.items() if k not in ('x','fields','translation_velocities')} for r in rows])
    last_npz=framesdir/f"step2_inc{rows[-1]['increment']:05d}.npz"
    with np.load(last_npz) as last:last_U=last['U'].copy()
    sampled_strain=base.fem2_fe_strain_diagnostics(mesh,last_U);write_json(case/'final_strain_diagnostics.json',sampled_strain)
    energies=io.parse_transient_dat_energies(case/'motion.dat',static_end_time=static_end,increments=sta['accepted_increments'],element_set='SOLID');write_json(case/'energy.json',energies)
    stdout_energy=io.parse_transient_stdout_energies(case/'motion.stdout.txt',static_end_time=static_end,increments=sta['accepted_increments']);write_json(case/'stdout_energy.json',stdout_energy)
    initial=[r for r in energies['records'] if r['step']==1 and abs(r['total_time']-static_end)<1e-7]
    dyn=[r for r in energies['records'] if r['step']==2]
    E0=initial[-1]['internal_energy'] if initial else None
    stdout_dynamic=[r for r in stdout_energy['records'] if r.get('step')==2]
    work_complete=bool(stdout_dynamic) and all('external_work' in r and 'damping_work' in r for r in stdout_dynamic)
    ext=max(abs(r['external_work']) for r in stdout_dynamic) if work_complete else None
    damp=max(abs(r['damping_work']) for r in stdout_dynamic) if work_complete else None
    if ext is not None and ext!=0. or damp is not None and damp!=0.:raise ValueError('Unexpected external/damping work after release')
    drift=max((abs(r['mechanical_energy']/E0-1) for r in dyn),default=None) if E0 else None
    if not math.isclose(rows[-1]['time'],end,rel_tol=0,abs_tol=1e-8):raise ValueError('Dynamic target horizon not reached')
    if rows[0]['fields'][20,1]>=static_fields[20,1]:raise ValueError('Early displacement has no restoring change')
    record={'status':'PASS','preload_transfer':preload,'dynamic_time_start':rows[0]['time'],'dynamic_time_end':rows[-1]['time'],'dynamic_output_frames':len(rows),'static_increments':len(static_rows),'dynamic_increments':len([r for r in sta['accepted_increments'] if r['step']==2]),'accepted_increments':len(sta['accepted_increments']),'cutbacks':sta['reported_cutbacks'],'max_DAT_FRD_displacement_difference':max_round,'initial_midspan_w':static_fields[20,1],'first_midspan_w':rows[0]['fields'][20,1],'final_midspan_w':rows[-1]['fields'][20,1],'first_midspan_w_velocity':rows[0]['translation_velocities'][20,1],'maximum_native_external_work_after_release':ext,'maximum_native_damping_work_after_release':damp,'energy_status':'PASS' if dyn and all('mechanical_energy' in r for r in dyn) and drift is not None and work_complete else 'PARTIAL','native_initial_internal_energy':E0,'max_native_relative_mechanical_energy_drift':drift,'energy_qualification':'one 3D timestep level; native output rounding and initial acceleration regularization, no temporal convergence claim','final_strain_diagnostics':sampled_strain,'recovery_seconds':time.perf_counter()-start}
    write_json(case/'recovery.json',record);s['numerical_seconds']+=record['recovery_seconds'];return record


def run_case(b,c,item,s,kind):
    old,source,mesh,audit=verify_sources(c);case=b/'cases'/kind;previous=s['cases'].get(kind)
    if previous and previous.get('status')=='PASS':return True
    if previous and previous.get('status')!='OUTPUT_RECOVERY_PENDING':return False
    if not previous:
        if len(s['attempts'])>=2:raise RuntimeError('Two production attempts exhausted')
        case.mkdir(parents=True);write_json(case/'input_contract.json',write_input(case/'motion.inp',c,old,source,mesh,audit,kind=='nonlinear'))
        s['attempts'].append({'kind':kind,'ordinal':len(s['attempts'])+1,'status':'STARTED','input_sha256':sha(case/'motion.inp')});s['cases'][kind]={'status':'STARTED'};finalize(b,item,s)
        env=dict(os.environ);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBER_OF_CPUS='1')
        remaining=c['numerical_budget_seconds']-s['numerical_seconds'];start=time.perf_counter();s['job_calls']['CCX_production']+=1;finalize(b,item,s)
        try:
            if remaining<=0:raise TimeoutError('Declared numerical budget exhausted')
            print('FEM3A production static+dynamic job: '+kind,flush=True)
            result,stats=base.fem1.run_job([old['science_config']['ccx_exe'],'motion'],case,min(1200,remaining),4*1024**3,case/'motion',env)
            write_json(case/'job.json',stats);s['cases'][kind]['job']=stats
            text=(case/'motion.stdout.txt').read_text(encoding='utf8',errors='replace')
            if result.returncode!=0 or stats['failure'] or 'JOB FINISHED' not in text.upper() or '*ERROR' in text.upper():raise RuntimeError('Actual solver failure: '+str(stats))
            s['cases'][kind]['status']='OUTPUT_RECOVERY_PENDING';s['attempts'][-1]['status']='SOLVER_FINISHED'
        except Exception as exc:
            s['cases'][kind].update(status='FAIL',failure=str(exc));s['attempts'][-1].update(status='FAIL',failure=str(exc));s['overall']='BLOCKED_BY_SOLVER';s['statuses']['NLSP_FEM3A_3D_'+kind.upper()+'_DYNAMIC']='FAIL'
        finally:
            s['numerical_seconds']+=time.perf_counter()-start;finalize(b,item,s)
        if s['cases'][kind]['status']=='FAIL':return False
    try:
        record=recover_case(b,c,s,kind,old,mesh,audit);record['job']=s['cases'][kind]['job'];s['cases'][kind]=record
        next(a for a in s['attempts'] if a['kind']==kind)['status']='PASS';s['statuses']['NLSP_FEM3A_3D_'+kind.upper()+'_DYNAMIC']='PASS';s['statuses']['NLSP_FEM3A_STATIC_STATE_TRANSFER']='PASS'
        s['overall']='PARTIAL';finalize(b,item,s);return True
    except Exception as exc:
        s['cases'][kind].update(status='OUTPUT_RECOVERY_PENDING',recovery_failure=str(exc));s['overall']='PARTIAL';finalize(b,item,s);return False


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--preflight',action='store_true');mode.add_argument('--run-pilot',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    p.add_argument('--config',type=Path,default=CONFIG);p.add_argument('--through-case',choices=['linear','nonlinear'],default='nonlinear');a=p.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_cache(b)
        if a.plot_only:plot_bundle(b)
        print(json.dumps({'bundle':str(b),'overall':s['overall'],'new_scientific_calls':0},indent=2));return s
    c=validate_config(read_json(a.config));found=existing_attempt(c)
    if found:b,item,s=found
    else:b,item,s=prepare(a.config)
    if a.preflight or s['overall'] in ('PILOT_COMPLETE_WITH_QUALIFICATIONS','BLOCKED_BY_SOLVER'):
        print(json.dumps({'bundle':str(b),'overall':s['overall'],'new_scientific_calls':0},indent=2));return s
    for kind in ('linear','nonlinear'):
        if not run_case(b,c,item,s,kind):break
        if a.through_case==kind:break
    if all(s['cases'].get(k,{}).get('status')=='PASS' for k in ('linear','nonlinear')):
        finish_references_and_comparison(b,c,item,s)
    finalize(b,item,s);print(json.dumps({'bundle':str(b),'statuses':s['statuses'],'job_calls':s['job_calls'],'numerical_seconds':s['numerical_seconds'],'overall':s['overall']},indent=2));return s

# Reference/comparison and plot paths are appended below before first execution.


def _matches(all_times, requested):
    indices=np.searchsorted(all_times,requested)
    if np.any(indices>=len(all_times)) or np.max(abs(np.asarray(all_times)[indices]-requested))>1e-12:raise ValueError('Actual-time reference sampling mismatch')
    return indices


def field_difference(first,second,x,scale=None):
    difference=np.asarray(first)-np.asarray(second)
    if scale is None:scale=max(float(np.max(abs(first))),float(np.max(abs(second))))
    L2=np.sqrt(np.trapezoid(difference**2,x,axis=-1));j=np.unravel_index(np.argmax(abs(difference)),difference.shape)
    return {'absolute_max':float(np.max(abs(difference))),'max_time_L2':float(np.max(L2)),'characteristic_scale':float(scale),'relative_max':float(np.max(abs(difference))/scale) if scale else None,'relative_max_L2':float(np.max(L2)/scale) if scale else None,'time_index_at_max':int(j[0]),'x_at_max':float(x[j[1]]),'signed_at_max':float(difference[j]),'sampled_maxima_only':True,'phase_amplitude_alignment':False}


def finish_references_and_comparison(b,c,item,s):
    start=time.perf_counter();end=dynamic_settings(c)['duration'];datasets={}
    for kind in ('linear','nonlinear'):
        with np.load(b/'cases'/kind/'section_history.npz') as z:datasets[kind]={n:z[n].copy() for n in z.files}
    all_times=np.unique(np.r_[0.,datasets['linear']['time'],datasets['nonlinear']['time'],end])
    if all_times[-1]>end:raise ValueError('Printed-time rounding must be resolved without extending the run')
    ref=one.load_reference(ROOT,c['source_static']['bundle'],c['source_fem1']['bundle'],c['source_action']['bundle'])
    one_data={}
    for kind in ('linear','nonlinear'):
        path=b/('one_d_'+kind+'.npz');meta_path=b/('one_d_'+kind+'.json')
        if path.exists():
            with np.load(path) as z:one_data[kind]={n:z[n].copy() for n in z.files}
            continue
        if kind=='linear':
            trajectory=one.exact_linear_reference(ref,all_times);metadata=trajectory['metadata']
        else:
            if s['job_calls']['1D_nonlinear_ODE']>=1:raise RuntimeError('The sole 1D integration was already attempted; no retry')
            s['job_calls']['1D_nonlinear_ODE']+=1;finalize(b,item,s)
            deadline=time.perf_counter()+max(0.,c['numerical_budget_seconds']-s['numerical_seconds']-(time.perf_counter()-start))
            trajectory=one.integrate_nonlinear_reference(ref,all_times,c,deadline,authorization={'user_authorized_FEM3A':True,'source_request':'93490460-7ff2-4977-ab37-f2e3187965aa/Pasted text.txt','execution_mode':'EXPLORATORY_NOT_CERTIFIED'})
            metadata=trajectory['stats']
        measured=one.summarize_reference(ref,trajectory,np.linspace(0.,1.,41),linear=kind=='linear')
        arrays={k:v for k,v in measured.items() if isinstance(v,np.ndarray)};np.savez_compressed(path,**arrays);one_data[kind]=arrays
        write_json(meta_path,{'execution':metadata,'diagnostics':measured['diagnostics'],'source':ref['source']})
        if kind=='nonlinear' and metadata['status']!='PASS':
            s['overall']='PARTIAL';s['statuses']['NLSP_FEM3A_1D_DYNAMIC_REFERENCE']='PARTIAL';s['numerical_seconds']+=time.perf_counter()-start;finalize(b,item,s);return
    comparison={'scope':'one medium mesh and one 3D timestep policy; observed pilot comparisons, no temporal/spatial accuracy certificate','field_maps':{'u':(0,0),'w':(1,1),'theta':(2,5),'c_eff_diagnostic':(3,6)},'absolute_physical_differences':{},'three_D_time_pairing':'only exact shared printed/endpoint timestamps; no interpolation or phase fitting','initial_fields_not_amplitude_aligned':True,'c_eff_not_identical_to_MH_DOF':True}
    for kind in ('linear','nonlinear'):
        d=datasets[kind];ids=_matches(one_data[kind]['times'],d['time']);a=one_data[kind]['fields'][ids];f=d['fields']
        comparison['absolute_physical_differences'][kind]={name:field_difference(a[:,:,i],f[:,:,j],d['x']) for name,(i,j) in comparison['field_maps'].items()}
        for name,metric in comparison['absolute_physical_differences'][kind].items():metric['time_at_max']=float(d['time'][metric['time_index_at_max']])
    shared=np.intersect1d(datasets['linear']['time'],datasets['nonlinear']['time'])
    if len(shared)<3:raise ValueError('Insufficient exactly shared actual times for 3D nonlinear correction')
    ia=_matches(datasets['linear']['time'],shared);ib=_matches(datasets['nonlinear']['time'],shared);i1=_matches(one_data['linear']['times'],shared);i2=_matches(one_data['nonlinear']['times'],shared)
    correction3=datasets['nonlinear']['fields'][ib]-datasets['linear']['fields'][ia]
    correction1=one_data['nonlinear']['fields'][i2]-one_data['linear']['fields'][i1]
    comparison['nonlinear_corrections']={name:field_difference(correction1[:,:,i],correction3[:,:,j],datasets['linear']['x']) for name,(i,j) in comparison['field_maps'].items()}
    for v in comparison['nonlinear_corrections'].values():v['time_at_max']=float(shared[v['time_index_at_max']])
    comparison['shared_dynamic_timestamps']=shared.tolist();comparison['maximum_1D_w_correction']=float(np.max(abs(correction1[:,:,1])));comparison['maximum_3D_w_correction']=float(np.max(abs(correction3[:,:,1])));comparison['dynamic_nonlinear_signal_certified']=False;comparison['qualification']='different L/NL initial static states are part of both IVPs; one timestep cannot certify dynamic correction'
    np.savez_compressed(b/'dynamic_comparison.npz',time=shared,x=datasets['linear']['x'],one_d_correction=correction1,three_d_correction=correction3)
    write_json(b/'dynamic_comparison.json',comparison)
    # Independent 41/81 reference-section sensitivity at the actual final frame.
    old,_,mesh,_=verify_sources(c);sensitivity={}
    for kind in ('linear','nonlinear'):
        meta=read_json(b/'cases'/kind/'frame_metadata.json')[-1]
        with np.load(b/meta['frame']) as z:U=z['U'].copy()
        recovered=base.recover_static_sections(mesh,U,old['preflight']['coefficients'],41,1.,True)
        sensitivity[kind]=recovered['recovery_sensitivity'];write_json(b/'cases'/kind/'final_recovery_sensitivity.json',sensitivity[kind])
    s['comparison']=comparison;s['statuses'].update(NLSP_FEM3A_1D_DYNAMIC_REFERENCE='PASS',NLSP_FEM3A_TRANSIENT_FIELD_RECOVERY='PASS',NLSP_FEM3A_SHORT_RESPONSE_COMPARISON='PASS',NLSP_FEM3A_ENERGY_DIAGNOSTICS='PASS' if all(s['cases'][k]['energy_status']=='PASS' for k in ('linear','nonlinear')) else 'PARTIAL')
    s['numerical_seconds']+=time.perf_counter()-start;s['overall']='PILOT_COMPLETE_WITH_QUALIFICATIONS';finalize(b,item,s);plot_bundle(b);finalize(b,item,s)


def plot_bundle(b):
    s=validate_cache(b)
    if not (b/'dynamic_comparison.npz').exists():return {'figures':0,'new_scientific_calls':0}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'savefig.dpi':240,'pdf.fonttype':42})
    T=s['preflight']['T1'];figdir=b/'figures';figdir.mkdir(exist_ok=True);figures=[]
    fig,axes=plt.subplots(1,3,figsize=(11,3.3),constrained_layout=True)
    for kind,style in [('linear','--'),('nonlinear','-')]:
        with np.load(b/('one_d_'+kind+'.npz')) as z:
            for ax,(index,position,label) in zip(axes,[(1,20,'w(L/2)'),(2,10,'theta(L/4)'),(0,10,'u(L/4)')]):ax.plot(z['times']/T,z['fields'][:,position,index],style,label='1D '+kind,lw=1.2)
        with np.load(b/'cases'/kind/'section_history.npz') as z:
            for ax,(index,position,label) in zip(axes,[(1,20,'w(L/2)'),(5,10,'theta(L/4)'),(0,10,'u(L/4)')]):ax.plot(z['time']/T,z['fields'][:,position,index],style,label='3D '+kind,lw=1.2)
    for ax,label in zip(axes,['w(L/2)','theta(L/4)','u(L/4)']):ax.set(xlabel='t/T1',ylabel=label,xlim=(0,.05));ax.grid(alpha=.25)
    axes[0].legend(fontsize=8);figures.append((fig,'linear_nonlinear_free_motion'))
    with np.load(b/'dynamic_comparison.npz') as z:
        fig,axes=plt.subplots(1,2,figsize=(8,3.3),constrained_layout=True)
        axes[0].plot(z['time']/T,z['one_d_correction'][:,20,1],label='1D NL minus L');axes[0].plot(z['time']/T,z['three_d_correction'][:,20,1],label='3D NL minus L')
        axes[0].set(xlabel='t/T1',ylabel='Midspan w correction');axes[0].legend(fontsize=8)
        for array,label in [('one_d_correction','1D'),('three_d_correction','3D')]:axes[1].plot(z['x'],z[array][-1,:,1],label=label)
        axes[1].set(xlabel='Material x/L',ylabel='Final w correction');axes[1].legend(fontsize=8)
        for ax in axes:ax.grid(alpha=.25)
        figures.append((fig,'nonlinear_minus_linear_short_response'))
    fig,axes=plt.subplots(1,2,figsize=(8,3.3),constrained_layout=True)
    for kind in ('linear','nonlinear'):
        with np.load(b/('one_d_'+kind+'.npz')) as z:axes[0].plot(z['times']/T,z['energy_relative_drift'],label='1D '+kind)
        e=read_json(b/'cases'/kind/'energy.json');r=[row for row in e['records'] if row['step']==2];E0=s['cases'][kind]['native_initial_internal_energy']
        if r and E0:axes[1].plot([row['dynamic_time']/T for row in r],[row['mechanical_energy']/E0-1 for row in r],label='3D '+kind)
    for ax in axes:ax.set(xlabel='t/T1',ylabel='(Eint + Ekin - E0)/E0',xlim=(0,.05));ax.grid(alpha=.25);ax.legend(fontsize=8)
    figures.append((fig,'energy_and_release_diagnostics'))
    for fig,name in figures:
        fig.savefig(figdir/(name+'.png'));fig.savefig(figdir/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    return {'figures':len(figures),'new_scientific_calls':0}


if __name__=='__main__':main()
