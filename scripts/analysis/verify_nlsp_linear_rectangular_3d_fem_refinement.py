"""One parent-validated FEM-1R continuation: .020 C3D10 mesh,24 eigenpairs.

Diagnostic orchestration only; immutable saved1D and three grids are loaded.
New parent-linked input/output contract preserves original three-grid cache.
No duplicate physics solver or automatic mesh/eigenpair extension.
"""
from __future__ import annotations
import argparse,csv,hashlib,importlib.metadata,json,math,os,shutil,subprocess,sys,time
from pathlib import Path
if __name__=='__main__':
    for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
        os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
import numpy as np
from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem as fem
CONFIG=ROOT/'data/input/nlsp_linear_rectangular_3d_fem_refinement.json'
OUTPUT=ROOT/'results/nlsp_linear_rectangular_3d_fem_refinement'
PARENT_SHA='9f2d5139b84b2aa133b20d9a7cae806bac085c178fba506e087cae2941da1d90'
OLD_LEVELS=('coarse','medium','fine')
STATUS_NAMES=('SOURCE_PRESERVATION','NEW_MESH_QUALITY','MODAL_EXECUTION','MODE_IDENTIFICATION','MESH_CONVERGENCE','ALL_FAMILY_COMPARISON')
sha,read_json,write_json=fem.sha,fem.read_json,fem.write_json


def validate_config(c):
    if c['schema']!='nlsp-linear-rectangular-3d-fem-refinement-v1':raise ValueError('Continuation schema mismatch')
    if c['parent_manifest_sha256']!=PARENT_SHA:raise ValueError('Frozen parent manifest identity changed')
    if c['geometry']!={'L':1.,'b':.2,'h':.1} or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}:raise ValueError('Frozen geometry/material changed')
    if c['mesh_level']!={'name':'refined','target_size':.02,'through_h_nominal':5}:raise ValueError('Exactly one refined .020 level authorized')
    if c['element']!='C3D10' or c['requested_eigenpairs']!=24:raise ValueError('Only C3D10 /24 eigenpairs authorized')
    if c['frozen_omega_window']!=3.4651360027859885 or c['numerical_mesh_convergence_relative']!=.001:raise ValueError('Frozen window/mesh tolerance changed')
    if c['matching']!={'minimum_mac':.7,'minimum_margin':.08,'dominance':.7,'section_residual_limit':.35,'section_count':41}:raise ValueError('Historical matching policy changed')
    if c['threads']!=1 or not 0<c['job_timeout_seconds']<=900 or not 0<c['numerical_budget_seconds']<=1200 or not 0<c['job_memory_limit_bytes']<=4*1024**3:raise ValueError('Bounded resource policy changed')
    if c['semantics']!={'linear_only':True,'model_fitting':False,'new_1d_solutions':False,'new_nonlinear_jobs':False,'automatic_mesh_or_eigenpair_extension':False}:raise ValueError('Unauthorized computation semantics')
    return c


def load_parent_bundle(path,expected_sha=PARENT_SHA):
    """Check historical hashes without requiring new orchestration code identity."""
    b=Path(path)
    if sha(b/'manifest.json')!=expected_sha:raise ValueError('Parent manifest SHA256 mismatch')
    summary=fem.validate_cache(b)
    pre,config=read_json(b/'preflight.json'),read_json(b/'frozen_config.json')
    if pre!=summary['preflight'] or config!=summary['config']:raise ValueError('Parent frozen reference/summary mismatch')
    if set(summary['meshes'])!=set(OLD_LEVELS) or len(pre['merged_spectrum'])!=8:raise ValueError('Parent inventory incomplete')
    for name in OLD_LEVELS:
        case=summary['meshes'][name];modal=case['modal']
        if case['status']!='PASS' or modal['eigenpairs']!=24 or not modal['all_eight_identified'] or not modal['full_frequency_window_covered']:raise ValueError('Parent modal evidence incomplete: '+name)
        if len({r['fem_mode'] for r in modal['matches']})!=8:raise ValueError('Parent assignments duplicate')
    return summary,pre,fem.load_profiles(b,pre)


def identity(config_path=CONFIG):
    c=validate_config(read_json(config_path));parent=ROOT/c['parent_bundle']
    source,pre,_=load_parent_bundle(parent,c['parent_manifest_sha256'])
    old=read_json(parent/'manifest.json')['identity']
    if c['geometry']!=pre['geometry'] or c['material']!=source['config']['material']:raise ValueError('Continuation geometry/material differs from parent')
    helpers=dict(old['model_sha256']);helpers['scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py']=old['code_sha256']
    for name,digest in helpers.items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen source helper changed: '+name)
    binaries={name:sha(c[name]) for name in ('gmsh_exe','ccx_exe')}
    if binaries!=old['executables']:raise ValueError('Meshing/solver binary differs from parent')
    extra='scripts/analysis/thickness_mismatch/audits/audit_full_spectrum_3d_fem_smoke_extraction.py';helpers[extra]=sha(ROOT/extra)
    dlls={p.name:sha(p) for p in sorted(Path(c['ccx_exe']).parent.glob('*.dll'))}
    item={'schema':c['schema'],'config':c,'config_sha256':sha(config_path),'code_sha256':sha(__file__),
          'parent_bundle':c['parent_bundle'],'parent_manifest_sha256':sha(parent/'manifest.json'),
          'helper_sha256':helpers,'executables':binaries,'ccx_runtime_dlls':dlls,'python':sys.version,
          'dependencies':{k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def cross_mesh_matching(parent,new_case,pre,policy):
    """Historical common-grid section mass metric; frequency values unused."""
    coeff=pre['coefficients'];old_modal=read_json(parent/'meshes/fine/modal_analysis.json');new_modal=read_json(new_case/'modal_analysis.json')
    ids=[int(r['raw_mode_number']) for r in new_modal['raw_frequencies']];mac=np.zeros((8,len(ids)))
    with np.load(parent/'meshes/fine/modal_vectors.npz',allow_pickle=False) as old,np.load(new_case/'modal_vectors.npz',allow_pickle=False) as new:
        for i,r in enumerate(old_modal['matches']):
            first={'x':old['profile_x_'+str(r['fem_mode'])],'fields':old['profile_q_'+str(r['fem_mode'])]}
            for j,mode in enumerate(ids):
                second={'x':new['profile_x_'+str(mode)],'fields':new['profile_q_'+str(mode)]}
                mac[i,j]=fem.nlsp_section_profile_mac(first,second,pre['geometry']['L'],coeff['m'],coeff['jp'],coeff['jb'])
    assigned=fem.nlsp_shape_assignment(mac,policy['minimum_mac'],policy['minimum_margin'])
    for row,old,new in zip(assigned,old_modal['matches'],new_modal['matches']):
        mode=ids[row['column']] if row['column'] is not None else None
        row.update(family=old['family'],local_mode=old['local_mode'],fine_mode=old['fem_mode'],refined_mode=mode,
                   agrees_with_independent_1d_assignment=mode is not None and mode==new.get('fem_mode'))
    return {'rows':assigned,'MAC':mac.tolist(),'all_consistent':all(r['status']=='MATCHED' and r['agrees_with_independent_1d_assignment'] for r in assigned),
            'method':'Historical translation/rotation mass MAC on common241-point grid; c_eff excluded; frequencies unused'}


def compare_four_meshes(pre,parent_summary,refined_case,cross):
    rows=[]
    for i,reference in enumerate(pre['merged_spectrum']):
        new=refined_case['modal']['matches'][i]
        if new['status']!='MATCHED' or not cross['rows'][i]['agrees_with_independent_1d_assignment'] or cross['rows'][i]['status']!='MATCHED':
            rows.append({'sorted_1d':reference['sorted_index'],'family':reference['family'],'local_mode':reference['local_mode'],'mesh_status':'MATCH_UNRESOLVED'});continue
        old=[parent_summary['meshes'][n]['modal']['matches'][i]['omega_3d'] for n in OLD_LEVELS];omega=float(new['omega_3d'])
        frequencies=np.array(old+[omega]);changes=np.abs(np.diff(frequencies));relative=changes/frequencies[1:]
        decreasing=bool(np.all(np.diff(changes)<=0));monotone=bool(np.all(np.diff(frequencies)<=0) or np.all(np.diff(frequencies)>=0))
        # Retain FEM-1's decreasing-successive-relative-change qualification.
        accepted=relative[-1]<=.001 and relative[-1]<=relative[-2]
        model=fem.nlsp_relative_frequency_difference(reference['omega'],omega);previous=fem.nlsp_relative_frequency_difference(reference['omega'],old[-1])
        rows.append({'sorted_1d':reference['sorted_index'],'family':reference['family'],'local_mode':reference['local_mode'],
                     'omega_1d':reference['omega'],'omega_coarse':old[0],'omega_medium':old[1],'omega_fine':old[2],'omega_refined':omega,
                     'coarse_medium_relative':float(relative[0]),'medium_fine_relative':float(relative[1]),'fine_refined_relative':float(relative[2]),
                     'coarse_medium_absolute':float(changes[0]),'medium_fine_absolute':float(changes[1]),'fine_refined_absolute':float(changes[2]),
                     'absolute_changes_decrease':decreasing,'frequency_monotone':monotone,
                     'regularity':'REGULAR_OBSERVED_TREND' if decreasing and monotone else 'IRREGULAR_OBSERVED_TREND',
                     'model_difference_change_percentage_points':100*(model['absolute_relative_difference']-previous['absolute_relative_difference']),
                     'previous_model_absolute_relative':previous['absolute_relative_difference'],**model,
                     'fem_mode_refined':new['fem_mode'],'shape_MAC_refined':new['mac'],'assignment_margin':new['margin'],
                     'fine_refined_section_MAC':cross['rows'][i]['mac'],'cross_mesh_assignment_margin':cross['rows'][i]['margin'],
                     'section_residual_fraction':new['section_residual_fraction'],'axial_warp_fraction':new['axial_warp_fraction'],
                     'mesh_status':'MESH_ACCEPTED_AT_PRESET_TOLERANCE' if accepted else 'MESH_UNRESOLVED'})
    return rows


def axial_diagnostic(new_case,pre,profiles,modal):
    r=next(row for row in modal['matches'] if row['family']=='axial_mh')
    if r['status']!='MATCHED':return {'status':'MATCH_UNRESOLVED'}
    p=profiles['axial_mh:1']
    with np.load(new_case/'modal_vectors.npz',allow_pickle=False) as data:
        x=data['profile_x_'+str(r['fem_mode'])].copy();q=r['full_shape_sign']*data['profile_q_'+str(r['fem_mode'])]/math.sqrt(r['fem_mass_norm'])
    reference=p['q']/math.sqrt(r['one_D_lift_mass_norm']);grid=np.linspace(0,pre['geometry']['L'],401)
    one=np.column_stack([np.interp(grid,p['x'],reference[:,k]) for k in (0,6)]);three=np.column_stack([np.interp(grid,x,q[:,k]) for k in (0,6)])
    weight=np.ones(len(grid));weight[[0,-1]]*=.5
    result={'u_section_MAC':fem.nlsp_weighted_mac(one[:,:1],three[:,:1],weight),'c_effective_section_MAC':fem.nlsp_weighted_mac(one[:,1:],three[:,1:],weight),
            'maximum_c_effective_difference':float(np.max(abs(one[:,1]-three[:,1]))),'c_characteristic_reference':float(np.max(abs(one[:,1]))),
            'c_characteristic_fem':float(np.max(abs(three[:,1]))),'c_status':'DIAGNOSTIC_EFFECTIVE_THICKNESS_STRAIN_NOT_EXACT_MH_COORDINATE',
            'qualification':'Arbitrary linear eigenvector normalization,41 slabs,no end-zone exclusion; not a finite-motion strain'}
    np.savez_compressed(new_case/'axial_contraction_profiles.npz',x=grid,one_D=one,refined=three)
    return result


def run_refinement(c,bundle,parent_summary,pre,profiles):
    start=time.perf_counter();deadline=start+c['numerical_budget_seconds']
    summary={'config':c,'parent_reference':{'path':c['parent_bundle'],'manifest_sha256':c['parent_manifest_sha256'],'one_D_recomputed':False,'old_FEM_jobs_repeated':False},
             'job_calls':{'gmsh':0,'ccx':0},'statuses':{'NLSP_FEM1R_'+n:'NOT_RUN' for n in STATUS_NAMES},'comparisons':[]}
    summary['statuses']['NLSP_FEM1R_SOURCE_PRESERVATION']='PASS'
    case=bundle/'meshes/refined';case.mkdir(parents=True,exist_ok=True);paths=fem.new_case_paths(case)
    paths.geo.write_text(fem.rectangular_geo(**c['geometry'],target_size=.02)+'\nGeneral.NumThreads=1;\nMesh.RandomSeed=1;\n',encoding='utf8')
    single=fem.single;names=('L','E','RHO','NU','SOLID_MODES_REQUESTED');saved={n:getattr(single,n) for n in names}
    single.L=c['geometry']['L'];single.E=c['material']['E'];single.RHO=c['material']['rho'];single.NU=c['material']['nu'];single.SOLID_MODES_REQUESTED=24
    original=single.subprocess.run;jobs=[];record={'target_size':.02}
    env=dict(os.environ);env.update(OMP_NUM_THREADS='1',NUMBER_OF_CPUS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')

    def monitored(command,**kwargs):
        label='gmsh' if Path(command[0]).name.lower().startswith('gmsh') else 'ccx'
        if label=='ccx' and summary['job_calls']['ccx']>=1:raise RuntimeError('Only one new CalculiX job authorized')
        if label=='gmsh' and summary['job_calls']['gmsh']>=2:raise RuntimeError('Only two historical format exports for one mesh authorized')
        remaining=deadline-time.perf_counter()
        if remaining<=0:raise TimeoutError('TOTAL_NUMERICAL_BUDGET')
        summary['job_calls'][label]+=1
        result,stats=fem.run_job(command,kwargs.get('cwd',case),min(c['job_timeout_seconds'],remaining),c['job_memory_limit_bytes'],case/(label+'_job'+str(summary['job_calls'][label])),env);jobs.append(stats)
        if stats['failure']:raise RuntimeError(stats['failure'])
        if kwargs.get('check') and result.returncode:raise subprocess.CalledProcessError(result.returncode,command,result.stdout,result.stderr)
        return result

    single.subprocess.run=monitored
    try:
        generated,message,warnings=single.generate_mesh_with_gmsh_cli(paths,c['gmsh_exe'])
        record['gmsh_warnings']=list(warnings)
        if not generated:raise ArithmeticError(message)
        logs='\n'.join(p.read_text(encoding='utf8',errors='replace') for p in case.glob('gmsh_job*.txt'))
        audit,mesh=fem.audit_rectangular_mesh(paths.gmsh_inp,**c['geometry'],rho=c['material']['rho'],reader=single.read_gmsh_inp_mesh_data,gmsh_logs=logs)
        write_json(case/'mesh_audit.json',audit);record['mesh_audit']=audit
        if audit['status']!='PASS':
            summary['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY']='FAIL';raise ArithmeticError('SOLID_MESH_GATE_FAILED')
        summary['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY']='PASS'
        single.write_calculix_template(paths,mesh);outcome=single.run_calculix(paths,c['ccx_exe'],0.)
        record['solver_outcome']={'success':outcome.success,'message':outcome.message}
        record['solver_warnings']=fem.single.interesting_gmsh_lines(paths.ccx_stdout.read_text(encoding='utf8',errors='replace')+'\n'+paths.ccx_stderr.read_text(encoding='utf8',errors='replace'),prefix='ccx')
        if not outcome.success:
            summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']='FAIL';raise ArithmeticError(outcome.message)
        modal=fem.inspect_modal(c,pre,profiles,paths,audit,mesh,bundle);record['modal']=modal
        if modal['eigenpairs']!=24 or modal['parsed_eigenvectors']!=24:
            summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']='FAIL';raise ArithmeticError('Exactly24 complete eigenpairs required')
        summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']='PASS'
        parent=ROOT/c['parent_bundle'];cross=cross_mesh_matching(parent,case,pre,c['matching']);write_json(bundle/'fine_refined_shape_correspondence.json',cross)
        identified=modal['all_eight_identified'] and cross['all_consistent'] and modal['full_frequency_window_covered']
        summary['statuses']['NLSP_FEM1R_MODE_IDENTIFICATION']='PASS' if identified else 'PARTIAL'
        summary['comparisons']=compare_four_meshes(pre,parent_summary,record,cross)
        accepted=identified and all(r['mesh_status']=='MESH_ACCEPTED_AT_PRESET_TOLERANCE' for r in summary['comparisons'])
        summary['statuses']['NLSP_FEM1R_MESH_CONVERGENCE']='PASS' if accepted else 'PARTIAL'
        summary['statuses']['NLSP_FEM1R_ALL_FAMILY_COMPARISON']='PASS' if accepted else 'PARTIAL'
        write_json(bundle/'axial_contraction_diagnostic.json',axial_diagnostic(case,pre,profiles,modal))
        with (bundle/'four_mesh_comparison.csv').open('w',encoding='utf8',newline='') as f:
            columns=list(dict.fromkeys(k for r in summary['comparisons'] for k in r));w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(summary['comparisons'])
        with (bundle/'full_refined_spectrum.csv').open('w',encoding='utf8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(modal['raw_frequencies'][0]));w.writeheader();w.writerows(modal['raw_frequencies'])
        record['status']='PASS' if identified else 'PARTIAL'
    except Exception as exc:
        record.update(status='FAIL',failure=str(exc))
        if summary['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY']=='NOT_RUN':summary['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY']='FAIL'
        elif summary['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY']=='PASS' and summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']=='NOT_RUN' and summary['job_calls']['ccx']>0:summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']='FAIL'
        elif summary['statuses']['NLSP_FEM1R_MODAL_EXECUTION']=='PASS' and summary['statuses']['NLSP_FEM1R_MODE_IDENTIFICATION']=='NOT_RUN':summary['statuses']['NLSP_FEM1R_MODE_IDENTIFICATION']='FAIL'
    finally:
        single.subprocess.run=original
        for name,value in saved.items():setattr(single,name,value)
        record.update(jobs=jobs,total_seconds=time.perf_counter()-start,peak_working_set_bytes=max((j['peak_working_set_bytes'] for j in jobs),default=0))
        summary['refined_case']=record;summary['runtime']={'numerical_seconds':record['total_seconds'],'budget_seconds':c['numerical_budget_seconds']}
        summary['qualification']='Linear only; mesh differences not continuum error bounds; no extrapolation/fitting/nonlinear validation'
        write_json(case/'case.json',record);write_json(bundle/'summary.json',summary)
    return summary


def validate_cache(bundle,item=None):
    b=Path(bundle);manifest=read_json(b/'manifest.json')
    if item is not None and manifest['identity']!=item:raise ValueError('Continuation cache identity mismatch')
    summary=fem.validate_cache(b)
    load_parent_bundle(ROOT/summary['config']['parent_bundle'],summary['config']['parent_manifest_sha256'])
    return summary


def plot_only(bundle):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s=validate_cache(bundle);rows=[r for r in s['comparisons'] if 'omega_refined' in r]
    if not rows:return {'figures':0,'solver_calls':0}
    target=Path(bundle)/'figures';target.mkdir(exist_ok=True);plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    def finish(fig,name):
        fig.savefig(target/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None});fig.savefig(target/(name+'.png'),dpi=220);plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(12,5.4),layout='constrained');sizes=np.array([.05,1/30,.025,.02])
    labels={'inplane_bending':'In-plane bending','outplane_bending':'Out-of-plane bending','torsion':'Torsion','axial_mh':'Axial MH'}
    for ax,r in zip(axes.flat,rows):
        freq=np.array([r['omega_'+n] for n in (*OLD_LEVELS,'refined')]);ax.plot(sizes,100*(freq/freq[-1]-1),'o-',lw=1.3)
        ax.set(title=labels[r['family']]+' '+str(r['local_mode']),xlabel='Target mesh size',ylabel='Change from refined, %');ax.invert_xaxis();ax.grid(alpha=.2)
    finish(fig,'four_mesh_convergence')
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    for family,label in labels.items():
        selected=[r for r in rows if r['family']==family];x=[r['sorted_1d'] for r in selected]
        axes[0].plot(x,[100*r['previous_model_absolute_relative'] for r in selected],'x',alpha=.55);axes[0].plot(x,[100*r['absolute_relative_difference'] for r in selected],'o',label=label)
        axes[1].plot(x,[100*r['medium_fine_relative'] for r in selected],'x',alpha=.55);axes[1].plot(x,[100*r['fine_refined_relative'] for r in selected],'o')
    axes[1].axhline(.1,color='black',ls='--',lw=1,label='Preset0.1% criterion')
    axes[0].set(ylabel='Absolute1D /3D difference, %',title='Cross: fine; circle: refined')
    axes[1].set(ylabel='Successive mesh change, %',title='Cross: medium/fine; circle: fine/refined')
    for ax in axes:ax.set(xlabel='1D reference position');ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'refined_model_difference_and_mesh_change')
    return {'figures':2,'solver_calls':0,'new_1d_calls':0}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--check-source',action='store_true');modes.add_argument('--run-fem',action='store_true')
    modes.add_argument('--report-only',type=Path);modes.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--output-dir',type=Path,default=OUTPUT);a=parser.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_cache(b)
        if a.plot_only:plot_only(b)
        print(json.dumps({'bundle':str(b),'statuses':s['statuses'],'new_solver_eigen_BVP_calls':0},indent=2));return s
    key,item=identity(a.config);c=read_json(a.config);parent_summary,pre,profiles=load_parent_bundle(ROOT/c['parent_bundle'],c['parent_manifest_sha256'])
    if a.check_source:
        result={'parent_manifest_sha256':c['parent_manifest_sha256'],'source_preservation':'PASS','cache_identity':key,'one_D_reference_rows':len(pre['merged_spectrum']),'new_solver_eigen_BVP_calls':0}
        print(json.dumps(result,indent=2));return result
    b=a.output_dir/key
    if (b/'manifest.json').exists():
        s=validate_cache(b,item);print(json.dumps({'bundle':str(b),'cache_hit':True,'statuses':s['statuses'],'new_solver_eigen_BVP_calls':0},indent=2));return s
    if b.exists() and any(b.iterdir()):raise RuntimeError('Unmanifested partial attempt exists; no automatic retry authorized')
    b.mkdir(parents=True,exist_ok=True);write_json(b/'frozen_config.json',c);write_json(b/'provenance.json',item)
    (b/'execution_code').mkdir();shutil.copyfile(__file__,b/'execution_code'/Path(__file__).name)
    summary=run_refinement(c,b,parent_summary,pre,profiles)
    load_parent_bundle(ROOT/c['parent_bundle'],c['parent_manifest_sha256'])
    write_json(b/'manifest.json',fem.artifact_manifest(b,item));plot_only(b)
    write_json(b/'manifest.json',fem.artifact_manifest(b,item))
    print(json.dumps({'bundle':str(b),'statuses':summary['statuses'],'runtime':summary['runtime'],'job_calls':summary['job_calls']},indent=2));return summary


if __name__=='__main__':main()

