"""Selected EB spring cases: physical modes, bounded roots, tracked branches.

compute is missing-only; plot-only reads tables/NPZ without roots or tracking.
The original sorted map and all production physics remain unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import itertools
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
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_rotational_spring_eb_modes as modes
from scripts.analysis.laminated_beams import plot_inplane_rotational_spring_eb_beta as old

PLOT_KAPPAS = [1,100]
OUTPUT = ROOT/'results/laminated_beams/inplane_rotational_spring_eb_tracked'
LEGACY_TRACKED = OUTPUT
HINGE_OUTPUT = ROOT/'results/laminated_beams/inplane_rotational_spring_eb_tracked_k0'
COMPARISON_OUTPUT = ROOT/'results/laminated_beams/inplane_rotational_spring_eb_tracked_comparison'
ARM,FS = old.ARM,old.FS
CRITERIA = dict(mac=.95,margin=.20,symmetry_defect=1e-6,close_relative=.02,
    frequency_relative=1e-6,local_frequency_relative=1e-8,angle_degrees=1e-3,
    max_added_points=100,max_B_per_point=6000,max_event_recoveries=2,
    max_positions=10,shape_nodes=129)
BRANCH_COLORS = dict(zip([f'mode_{j:02d}' for j in range(1,7)],
    ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#333333']))
CALLS = Counter()
COMPARISON_STYLES = {0:dict(color='#0072B2',linestyle='-'),
                     1:dict(color='#D55E00',linestyle='--'),
                     100:dict(color='#009E73',linestyle='-.')}
COMPARISON_MARKERS = {0:('o',(15,55)),1:('s',(32,61)),100:('^',(47,57))}


def assert_protected_sources(directory=None):
    path=(directory or HINGE_OUTPUT)/'run_manifest.json'
    if path.exists():
        for name,digest in json.loads(path.read_text(encoding='utf-8')).get('protected_source_files_sha256',{}).items():
            if sha(ROOT/name)!=digest:raise ValueError(f'Protected source changed: {name}')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def point_id(kappa,beta):return f'k{kappa:g}_b{beta:.10g}'


def atomic(path,text):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(text,encoding='utf-8',newline='\n');tmp.replace(path)


def write_csv(path,rows):
    if not rows:return
    stream=io.StringIO(newline='');fields=list(dict.fromkeys(k for r in rows for k in r))
    writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    atomic(path,stream.getvalue())


def load():
    OUTPUT.mkdir(exist_ok=True)
    if PLOT_KAPPAS==[0]:
        manifest_path=OUTPUT/'run_manifest.json'
        manifest=json.loads(manifest_path.read_text(encoding='utf-8')) if manifest_path.exists() else {}
        if 'protected_source_files_sha256' not in manifest:
            manifest['protected_source_files_sha256']={p.relative_to(ROOT).as_posix():sha(p)
                for directory in (old.OUTPUT,LEGACY_TRACKED) for p in sorted(directory.rglob('*')) if p.is_file()}
            manifest['addition_initial_HEAD']=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
            atomic(manifest_path,json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        assert_protected_sources(OUTPUT)
    path=OUTPUT/'tracking_diagnostics.json'
    state=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    state.setdefault('points',{});state.setdefault('events',[]);state.setdefault('tracking_attempts',[])
    state.setdefault('added_points',[]);state.setdefault('root_searches',[]);state.setdefault('quadrature_checks',[])
    sources={str((old.OUTPUT/name).relative_to(ROOT)).replace('\\','/'):sha(old.OUTPUT/name)
             for name in ('spectrum_roots.csv','diagnostics.json','run_manifest.json')}
    if 'source_hashes' in state and state['source_hashes']!=sources:raise ValueError('Original sorted data changed')
    state['source_hashes']=sources
    state.setdefault('source_version',dict(HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        executable=sys.executable,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        code_status='working-tree version',criteria=CRITERIA,physics_hashes=old.contract()['physics_sha256']))
    selected=state.setdefault('selected_kappas',PLOT_KAPPAS)
    if selected!=PLOT_KAPPAS:raise ValueError('Output belongs to a different kappa selection')
    shapes={}
    if (OUTPUT/'shapes.npz').exists():
        with np.load(OUTPUT/'shapes.npz',allow_pickle=False) as data:shapes={k:data[k] for k in data.files}
    return state,shapes


def save(state,shapes=None):
    atomic(OUTPUT/'tracking_diagnostics.json',json.dumps(state,ensure_ascii=False,indent=1,allow_nan=False)+'\n')
    roots=[r for p in state['points'].values() for r in p['roots']]
    write_csv(OUTPUT/'verified_roots.csv',sorted(roots,key=lambda r:(r['kappa'],r['beta_deg'],r['current_sorted_position'])))
    if shapes is not None:
        with (OUTPUT/'shapes.npz.tmp').open('wb') as stream:np.savez_compressed(stream,**shapes)
        (OUTPUT/'shapes.npz.tmp').replace(OUTPUT/'shapes.npz')


class PointMatrices:
    def __init__(self,kappa,beta,counts_before=0):
        if not math.isfinite(kappa) or kappa<0 or not math.isfinite(beta):raise ValueError('finite nonnegative kappa and finite beta required')
        self.point=dict(point_id=point_id(kappa,beta),state=f'k{kappa:g}',beta_deg=beta,beta_rad=math.radians(beta),
                        mode='SPRING',kappa_theta=kappa,k_theta=kappa*ARM.D/ARM.L,grid_role='BASE')
        self.transfers=old.Transfers();self.full=old.Provider(self.point,self.transfers)
        self.blocks={};self.count=counts_before;self.full_count=0;self.block_count=0

    def tick(self,kind):
        if self.count>=CRITERIA['max_B_per_point']:raise old.CostLimit('COST_LIMIT')
        self.count+=1;CALLS[kind]+=1

    def assembly(self,omega):
        if float(omega) not in self.full.cache:
            self.tick('full_B');self.full_count+=1
        return self.full.assembly(omega)

    def block(self,omega,parity):
        key=(float(omega),parity)
        if key not in self.blocks:
            self.tick('symmetry_B');self.block_count+=1
            endpoints=self.transfers.endpoints(omega)
            self.blocks[key]=modes.class_matrix(endpoints,self.point['beta_rad'],self.point['k_theta'],ARM,parity)
        return self.blocks[key]


def reconstruct(context,Omega,parity=None,nodes=129):
    started=time.perf_counter();omega=Omega/FS
    result=modes.recover(context.assembly(omega),omega,context.point['beta_rad'],context.full.joint,ARM,nodes,parity)
    CALLS['shape_reconstructions']+=1;CALLS['analytic_arm_evaluations']+=result['analytic_arm_evaluations']
    CALLS['shape_seconds']+=time.perf_counter()-started
    return result


def attach_shape(state,shapes,context,root,parity=None):
    result=reconstruct(context,root['Omega'],parity)
    if parity is not None and result['symmetry_class']!=parity:raise ValueError('WRONG_REFLECTION_CLASS')
    if root['multiplicity']==2 and result['detected_nullity']!=2:raise ValueError('SAVED_MULTIPLICITY_NOT_CONFIRMED')
    key=point_id(root['kappa'],root['beta_deg']).replace('.','p')+f"_r{root['current_sorted_position']:02d}"
    for name in ('states','reactions','vector'):shapes[key+'__'+name]=result.pop(name)
    root.update(shape_key=key,root_status='CONFIRMED' if not result['failures'] else 'SHAPE_UNCONFIRMED',
                symmetry_class=result['symmetry_class'],cluster_id='',mass_norm=1.)
    state.setdefault('shape_checks',{})[key]=result
    return root


def import_saved(state,shapes):
    source=json.loads((old.OUTPUT/'diagnostics.json').read_text(encoding='utf-8'))
    with (old.OUTPUT/'spectrum_roots.csv').open(encoding='utf-8',newline='') as stream:csv_rows=list(csv.DictReader(stream))
    c=source['contract'];expected=old.contract()
    assert c['geometry']==expected['geometry'] and c['arm']==expected['arm']
    assert c['normalization']==expected['normalization']
    assert c['physics_sha256']==expected['physics_sha256']
    audit=[];started=time.perf_counter();accepted_points=frequency_rows=0
    for kappa in PLOT_KAPPAS:
        for tick in old.grid_tenths():
            beta=tick/10;pid=point_id(kappa,beta);saved=source['points'][pid]
            if saved['status'] not in old.GOOD:
                audit.append(dict(point_id=pid,status=saved['status'],blank_ROOT_rows=sum(r['point_id']==pid and not r['Omega'] for r in csv_rows)))
                continue
            accepted_points+=1;frequency_rows+=len(saved['rows'])
            if pid in state['points'] and state['points'][pid]['roots'] and all(r.get('shape_key','')+'__vector' in shapes for r in state['points'][pid]['roots']):continue
            values=[r['Omega'] for r in saved['rows']]
            assert len(values)>=7 and values==sorted(values) and min(values)>0
            assert [r['sorted_position'] for r in saved['rows']]==list(range(1,len(values)+1))
            assert saved['case']['kappa_theta']==kappa and saved['case']['k_theta']==kappa*ARM.D/ARM.L
            assert saved['case']['beta_deg']==beta and saved['case']['beta_rad']==math.radians(beta)
            context=PointMatrices(kappa,beta);roots=[]
            for r,parity in zip(saved['rows'],modes.saved_root_classes(saved['rows'])):
                assert math.isclose(r['Lambda']**2,r['Omega'],rel_tol=1e-14)
                assert math.isclose(r['omega']*FS,r['Omega'],rel_tol=1e-14)
                csvrow=next(x for x in csv_rows if x['point_id']==pid and int(x['sorted_position'])==r['sorted_position'])
                assert all(float(csvrow[k])==r[k] for k in ('omega','Omega','Lambda'))
                root=dict(kappa=kappa,beta_deg=beta,current_sorted_position=r['sorted_position'],omega=r['omega'],
                    Omega=r['Omega'],Lambda=r['Lambda'],source='REUSED_SORTED_ROOT',source_group_status=saved['status'],
                    multiplicity=r['multiplicity'],grid_role='BASE',role=r['role'])
                try:roots.append(attach_shape(state,shapes,context,root,parity))
                except ValueError as error:
                    root.update(root_status='SHAPE_UNCONFIRMED',reason=str(error));roots.append(root)
            state['points'][pid]=dict(kappa=kappa,beta_deg=beta,roots=roots,status='REUSED',B=context.count,
                full_B=context.full_count,symmetry_B=0,transfer_expm=context.transfers.expm_calls)
            if PLOT_KAPPAS==[0]:
                save(state,shapes)
                if tick%100==0:print('import',pid,'modes',len(roots),flush=True)
    state.setdefault('source_data_audit',dict(missing=audit,source_accepted_points=accepted_points,source_frequency_rows=frequency_rows,
        parameter_checks='geometry,rigidities,normalization,physical hashes,CSV/JSON exact agreement,order,Lambda^2=Omega',
        no_source_changes=True))
    state.setdefault('timing',{})['reuse_and_shapes_seconds']=time.perf_counter()-started
    save(state,shapes)


def local_windows(state,kappa,beta,parity):
    neighbours=sorted([p for p in state['points'].values() if p['kappa']==kappa and
        p['beta_deg']!=beta and any(r.get('symmetry_class')==parity for r in p['roots'])],key=lambda p:p['beta_deg'])
    left=[p for p in neighbours if p['beta_deg']<beta];right=[p for p in neighbours if p['beta_deg']>beta]
    chosen=([left[-1]] if left else [])+([right[0]] if right else [])
    if not chosen:raise ValueError('No accepted neighbour for class predictor')
    arrays=[[r['Omega'] for r in p['roots'] if r.get('symmetry_class')==parity and r['root_status']=='CONFIRMED'] for p in chosen]
    windows=[]
    for j in range(max(map(len,arrays))):
        values=[a[j] for a in arrays if len(a)>j]
        lo,hi=min(values),max(values)
        width=max(.25,.003*hi,1.5*(hi-lo))
        windows.append((max(1e-8,lo-width),hi+width))
    return old.merge_windows(windows)


def search_point(state,shapes,kappa,beta,trigger,windows_override=None,persist=True):
    pid=point_id(kappa,beta)
    if pid in state['points'] and state['points'][pid]['roots']:return state['points'][pid]
    base=any(abs(beta-t/10)<1e-10 for t in old.grid_tenths())
    if not base:
        if pid not in state['added_points'] and len(state['added_points'])>=CRITERIA['max_added_points']:raise RuntimeError('ADDED_ANGLE_COST_LIMIT')
        if pid not in state['added_points']:state['added_points'].append(pid)
    started=time.perf_counter();context=PointMatrices(kappa,beta);events=[];records=[];errors=[];suspects=[]
    try:
        class_windows=windows_override if windows_override is not None else {parity:local_windows(state,kappa,beta,parity) for parity in (1,-1)}
        upper=max(w[-1][1] for w in class_windows.values())+.1
        for parity in (1,-1):
            windows=class_windows[parity]
            pool=[];provider=lambda w:context.block(w,parity)
            for lo,hi,n in old.search_intervals(windows,upper):
                CALLS['detector_calls']+=1
                pool.extend(old.pilot.roots._scan_candidates(provider,FS,old.pilot.policy(lo,hi),
                    case_id=pid,builder_id=f'EB_reflection_{parity:+d}',scan_id='SYMMETRY_LOCAL',
                    points=17 if n==old.LIMITS['window_points'] else 5,phases=(0.,))[0])
            pool,proof=old.pilot.reconcile_local_detections(pool,provider)
            accepted,ambiguous=old.pilot.consolidate(pool)
            suspects.extend(c for c in pool if old.pilot.suspicious(c))
            records.append(dict(parity=parity,windows=windows,candidates=[old.pilot.candidate_record(c) for c in pool],reconciliations=proof))
            if ambiguous:errors.append('CLASS_DETECTION_AMBIGUOUS')
            for event in accepted:
                events.append((event.omega_bar,parity,event.diagnostics.detected_nullity))
        events.sort()
        # Each class is solved separately; close roots in different classes are
        # retained independently, never merged by their distance.
        if len(events)<7:errors.append('MISSING_TARGET_OR_GUARD')
        if len(events)>CRITERIA['max_positions']:errors.append('POSITION_LIMIT')
        if len(events)>=7:
            if any(c.interval_left_bar<=events[5][0] for c in suspects):errors.append('UNRESOLVED_BELOW_TARGET')
            if upper-events[6][0]<=old.LIMITS['guard_margin_Omega']:errors.append('GUARD_NOT_SEPARATED')
        roots=[]
        for position,(Omega,parity,multiplicity) in enumerate(events,1):
            root=dict(kappa=kappa,beta_deg=beta,current_sorted_position=position,Omega=Omega,omega=Omega/FS,
                Lambda=math.sqrt(Omega),source='NEW_CHARACTERISTIC_ROOT',source_group_status=trigger,
                multiplicity=multiplicity,grid_role='BASE' if base else 'ADDED')
            try:roots.append(attach_shape(state,shapes,context,root,parity))
            except ValueError as error:
                root.update(root_status='SHAPE_UNCONFIRMED',symmetry_class=parity,reason=str(error));roots.append(root)
        if any(r['root_status']!='CONFIRMED' for r in roots[:6]):errors.append('ROOT_OR_SHAPE_GATE')
    except (ValueError,RuntimeError,old.CostLimit) as error:
        roots=[];errors.append(str(error))
    group=dict(kappa=kappa,beta_deg=beta,roots=roots,status='CONFIRMED' if not errors else 'POINT_UNCONFIRMED',
        errors=errors,B=context.count,full_B=context.full_count,symmetry_B=context.block_count,
        transfer_expm=context.transfers.expm_calls,seconds=time.perf_counter()-started,trigger=trigger,
        common_class_upper=upper if 'upper' in locals() else None,
        guard_warnings=[old.pilot.candidate_record(c) for c in suspects])
    state['points'][pid]=group
    state['root_searches'].append(dict(point_id=pid,records=records,errors=errors,seconds=group['seconds']))
    if persist:save(state)
    print('root point',pid,group['status'],'roots',len(roots),'B',context.count,flush=True)
    return group


def candidate_vectors(point,shapes):
    roots=[r for r in point['roots'] if r['root_status']=='CONFIRMED' and r['shape_key']+'__vector' in shapes]
    return roots,np.array([shapes[r['shape_key']+'__vector'] for r in roots])


def match(previous,current,shapes):
    roots,vectors=candidate_vectors(current,shapes)
    columns,mac,margins=modes.assign([shapes[r['shape_key']+'__vector'] for r in previous],vectors,
        [r['symmetry_class'] for r in previous],[r['symmetry_class'] for r in roots])
    bad=[]
    for i,j in enumerate(columns):
        if mac[i,j]<CRITERIA['mac'] or margins[i]<CRITERIA['margin']:bad.append(i)
    # A same-class order exchange is a veering/assignment trigger even when
    # a large step gives a deceptively high shape correlation.
    for i in range(len(previous)):
        for j in range(i):
            if previous[i]['symmetry_class']==previous[j]['symmetry_class']:
                if (previous[i]['Omega']-previous[j]['Omega'])*(roots[columns[i]]['Omega']-roots[columns[j]]['Omega'])<0:
                    bad.extend([i,j])
    return [roots[j] for j in columns],mac[np.arange(len(previous)),columns],margins,sorted(set(bad))


def track(state,shapes,allow_new=True):
    started=time.perf_counter();rows=[]
    for kappa in PLOT_KAPPAS:
        seed=state['points'][point_id(kappa,0.)]['roots'][:6]
        previous=[dict(r,branch_id=f'mode_{j:02d}',seed_sorted_position=j) for j,r in enumerate(seed,1)]
        def emit(roots,beta,mac,margin,status):
            for i,r in enumerate(roots):
                rows.append(dict(kappa=kappa,beta_deg=beta,branch_id=f'mode_{i+1:02d}',seed_sorted_position=i+1,
                    current_sorted_position=r['current_sorted_position'],omega=r['omega'],Omega=r['Omega'],Lambda=r['Lambda'],
                    root_status=r['root_status'],tracking_status=status[i],MAC=float(mac[i]),competing_assignment_margin=float(margin[i]),
                    cluster_id=r.get('cluster_id',''),symmetry_class=r['symmetry_class'],source=r['source'],shape_key=r['shape_key'],
                    grid_role=r.get('grid_role','BASE')))
        emit(previous,0.,np.ones(6),np.ones(6),['SEED_CONFIRMED']*6)
        previous_beta=0.
        queue=[(beta,None) for beta in sorted(p['beta_deg'] for p in state['points'].values() if p['kappa']==kappa and p['beta_deg']>0)]
        attempts=Counter(a.get('event_id') for a in state['tracking_attempts'])
        while queue:
            beta,event_key=queue.pop(0);point=state['points'][point_id(kappa,beta)]
            try:selected,mac,margin,bad=match(previous,point,shapes)
            except (ValueError,KeyError):
                selected=previous;mac=np.zeros(6);margin=np.zeros(6);bad=list(range(6))
            event_key=event_key or f'{kappa}:{previous_beta:g}:{beta:g}'
            if bad and allow_new and attempts[event_key]<CRITERIA['max_event_recoveries'] and len(state['added_points'])<CRITERIA['max_added_points']:
                mid=(previous_beta+beta)/2
                if not any(item[0]==mid for item in queue) and point_id(kappa,mid) not in state['points']:
                    attempts[event_key]+=1
                    state['tracking_attempts'].append(dict(event_id=event_key,interval=[previous_beta,beta],kappa=kappa,bad_branches=[f'mode_{i+1:02d}' for i in bad],MAC=mac.tolist(),margin=margin.tolist(),new_beta=mid))
                    search_point(state,shapes,kappa,mid,'TRACKING_AMBIGUOUS')
                    queue=[(mid,event_key),(beta,event_key)]+queue;continue
            statuses=['TRACKING_AMBIGUOUS' if i in bad else 'TRACKED' for i in range(6)]
            # No confirmed root exists in an unresolved slot: keep an explicit
            # blank in the tracked output, never a previous frequency at beta.
            start=len(rows);emit(selected,beta,mac,margin,statuses)
            for i in bad:
                for field in ('omega','Omega','Lambda','current_sorted_position'):rows[start+i][field]=None
            if not bad:
                for a,b in zip(previous,selected):
                    dot=np.vdot(shapes[a['shape_key']+'__vector'],shapes[b['shape_key']+'__vector'])
                    if dot.real<0:
                        for name in ('states','reactions','vector'):shapes[b['shape_key']+'__'+name]*=-1
                previous=selected;previous_beta=beta
    rows.sort(key=lambda r:(r['kappa'],r['beta_deg'],r['seed_sorted_position']))
    state['tracked_rows']=rows;state.setdefault('timing',{})['tracking_seconds']=time.perf_counter()-started
    write_csv(OUTPUT/'tracked_branches.csv',rows);save(state,shapes)
    print('tracking',Counter(r['tracking_status'] for r in rows),flush=True)
    return rows


def event_rows_at(state,kappa,beta):
    return sorted([r for r in state['tracked_rows'] if r['kappa']==kappa and r['beta_deg']==beta],key=lambda r:r['seed_sorted_position'])


def local_repeat(state,root,shapes):
    """One bounded accuracy check in a changed bracket; same exact class."""
    pid=point_id(root['kappa'],root['beta_deg']);group=state['points'][pid]
    context=PointMatrices(root['kappa'],root['beta_deg'],group['B'])
    parity=root['symmetry_class'];center=root['Omega'];width=max(1e-5,center*1e-6)
    provider=lambda w:context.block(w,parity)
    candidates=old.pilot.roots._scan_candidates(provider,FS,old.pilot.policy(center-width,center+width),
        case_id=pid,builder_id=f'EB_reflection_{parity:+d}',scan_id='BRACKET_ACCURACY',points=17,phases=(0.,))[0]
    CALLS['detector_calls']+=1
    pool,proof=old.pilot.reconcile_local_detections(candidates,provider)
    accepted,ambiguous=old.pilot.consolidate(pool)
    result=dict(reference_Omega=center,bracket=[center-width,center+width],parity=parity,
        candidates=[old.pilot.candidate_record(c) for c in pool],reconciliations=proof,status='UNRESOLVED')
    if len(accepted)==1 and not ambiguous:
        freq=accepted[0].omega_bar;shape=reconstruct(context,freq,parity)
        result.update(repeated_Omega=freq,relative_difference=abs(freq-center)/center,
            shape_MAC=float(modes.mac_matrix([shapes[root['shape_key']+'__vector']],[shape['vector']])[0,0]),
            physical_max=max(shape['physical_residuals']),failures=shape['failures'])
        if not shape['failures'] and result['relative_difference']<=CRITERIA['local_frequency_relative']:
            result['status']='LOCAL_AGREEMENT'
    result['B']=context.count-group['B']
    group['B']=context.count;group['full_B']+=context.full_count;group['symmetry_B']+=context.block_count
    group['transfer_expm']+=context.transfers.expm_calls
    return result


def crossing_events(state,shapes):
    started=time.perf_counter();initial=[]
    for kappa in PLOT_KAPPAS:
        rows=[r for r in state['tracked_rows'] if r['kappa']==kappa and any(abs(r['beta_deg']-t/10)<1e-10 for t in old.grid_tenths())]
        for i,j in itertools.combinations(range(6),2):
            a=[r for r in rows if r['seed_sorted_position']==i+1];b=[r for r in rows if r['seed_sorted_position']==j+1]
            for n in range(1,len(a)):
                if any(r['tracking_status'] not in ('TRACKED','SEED_CONFIRMED') for r in (a[n-1],a[n],b[n-1],b[n])):continue
                d0=a[n-1]['Omega']-b[n-1]['Omega'];d1=a[n]['Omega']-b[n]['Omega']
                if d0*d1<0:initial.append((kappa,i,j,a[n-1]['beta_deg'],a[n]['beta_deg']))
    for kappa,i,j,lo,hi in initial:
        eid=f'k{kappa}_mode_{i+1:02d}_mode_{j+1:02d}'
        if any(e['event_id']==eid for e in state['events']):continue
        origin=[lo,hi];left=event_rows_at(state,kappa,lo);right=event_rows_at(state,kappa,hi);anchor=left
        def snapshot(beta,selected,mac):
            return dict(beta_deg=beta,Omega_a=selected[i]['Omega'],Omega_b=selected[j]['Omega'],
                difference=selected[i]['Omega']-selected[j]['Omega'],MAC_min=float(min(mac)),
                shape_a=selected[i]['shape_key'],shape_b=selected[j]['shape_key'],
                symmetry_a=selected[i]['symmetry_class'],symmetry_b=selected[j]['symmetry_class'])
        traces=[snapshot(lo,left,np.ones(6)),snapshot(hi,right,np.ones(6))];error=None
        try:
            while hi-lo>CRITERIA['angle_degrees']:
                mid=(lo+hi)/2;point=search_point(state,shapes,kappa,mid,'CROSSING_LOCALIZATION:'+eid)
                if point['status'] not in ('CONFIRMED','REUSED'):raise ValueError('LOCAL_ROOT_UNCONFIRMED')
                selected,mac,margin,bad=match(anchor,point,shapes)
                if bad:raise ValueError('LOCAL_TRACKING_AMBIGUOUS')
                entry=snapshot(mid,selected,mac);traces.append(entry)
                delta=entry['difference'];dleft=left[i]['Omega']-left[j]['Omega']
                if delta*dleft<0:hi=mid;right=selected
                else:lo=mid;left=selected
            repeats=[local_repeat(state,r,shapes) for r in (left[i],left[j],right[i],right[j])]
        except (ValueError,RuntimeError) as exc:
            error=str(exc);repeats=[]
        accepted_signs=[t for t in traces if abs(t['difference'])>5*CRITERIA['local_frequency_relative']*max(t['Omega_a'],t['Omega_b'])]
        brackets=[(a,b) for a in accepted_signs for b in accepted_signs if a['beta_deg']<b['beta_deg'] and a['difference']*b['difference']<0]
        bracket=min(brackets,key=lambda pair:pair[1]['beta_deg']-pair[0]['beta_deg']) if brackets else (traces[0],traces[1])
        first,last=bracket
        closest=min(traces,key=lambda t:abs(t['difference']))
        independent=anchor[i]['symmetry_class']!=anchor[j]['symmetry_class']
        status='CROSSING_SUPPORTED' if independent and not error and all(r['status']=='LOCAL_AGREEMENT' for r in repeats) else 'UNRESOLVED_CLOSE_CLUSTER'
        correlations=modes.principal_correlations(
            [shapes[traces[0][key]+'__vector'] for key in ('shape_a','shape_b')],
            [shapes[traces[1][key]+'__vector'] for key in ('shape_a','shape_b')])
        event=dict(event_id=eid,kappa=kappa,branch_a=f'mode_{i+1:02d}',branch_b=f'mode_{j+1:02d}',
            classification=status,beta_left=first['beta_deg'],beta_right=last['beta_deg'],
            bracket_width=last['beta_deg']-first['beta_deg'],difference_left=first['difference'],difference_right=last['difference'],
            Omega_a_left=first['Omega_a'],Omega_b_left=first['Omega_b'],Omega_a_right=last['Omega_a'],Omega_b_right=last['Omega_b'],
            closest_beta=closest['beta_deg'],closest_Omega_a=closest['Omega_a'],closest_Omega_b=closest['Omega_b'],
            minimum_sampled_gap=abs(closest['difference']),symmetry_a=anchor[i]['symmetry_class'],symmetry_b=anchor[j]['symmetry_class'],
            subspace_correlation_min=float(min(correlations)),initial_bracket=origin,traces=traces,local_repeats=repeats,
            error=error,additional_accuracy_attempts=1 if repeats else 0,
            evidence='Exact independent reflection classes, shape continuation, resolved sign reversal; finite bracket, not an exact sampled double root')
        for t in traces:
            for key in ('shape_a','shape_b'):
                for r in state['points'][point_id(kappa,t['beta_deg'])]['roots']:
                    if r.get('shape_key')==t[key]:r['cluster_id']=eid
        state['events'].append(event);save(state,shapes)
        print('event',eid,status,first['beta_deg'],last['beta_deg'],'gap',event['minimum_sampled_gap'],flush=True)
    flat=[{k:v for k,v in e.items() if k not in ('traces','local_repeats','initial_bracket')} for e in state['events']]
    write_csv(OUTPUT/'crossing_events.csv',flat)
    state.setdefault('timing',{})['events_seconds']=time.perf_counter()-started


def quadrature_checks(state,shapes):
    if state['quadrature_checks']:return
    selected=[]
    for event in state['events']:
        for t in (event['traces'][0],min(event['traces'],key=lambda x:abs(x['difference'])),event['traces'][1]):
            selected.extend([t['shape_a'],t['shape_b']])
    for kappa in PLOT_KAPPAS:
        selected.extend(r['shape_key'] for r in state['points'][point_id(kappa,0)]['roots'][::2])
    by_key={r['shape_key']:r for p in state['points'].values() for r in p['roots'] if r.get('shape_key')}
    fine={};started=time.perf_counter()
    for key in dict.fromkeys(selected):
        r=by_key[key];reactions=shapes[key+'__reactions'];record=dict(shape_key=key,Omega=r['Omega'])
        for n in (65,257):
            xi,weights=modes.quadrature(n)
            values=np.array([modes.arm_states(r['omega'],ARM,arm_r,xi) for arm_r in reactions])
            vector=modes.mass_vector(values,ARM,weights);mass=float(np.vdot(vector,vector).real)
            record[f'mass_{n}']=mass
            if n==257:fine[key]=vector/np.sqrt(mass)
            CALLS['quadrature_arm_evaluations']+=2
        record['mass_129']=float(np.vdot(shapes[key+'__vector'],shapes[key+'__vector']).real)
        record['relative_129_257']=abs(record['mass_129']-record['mass_257'])/record['mass_257']
        state['quadrature_checks'].append(record)
    mac_checks=[]
    for event in state['events']:
        for field in ('shape_a','shape_b'):
            a,b=event['traces'][0][field],event['traces'][1][field]
            coarse=float(modes.mac_matrix([shapes[a+'__vector']],[shapes[b+'__vector']])[0,0])
            refined=float(modes.mac_matrix([fine[a]],[fine[b]])[0,0])
            mac_checks.append(dict(event=event['event_id'],branch=field,MAC_129=coarse,MAC_257=refined,difference=abs(coarse-refined)))
    state['quadrature_MAC_checks']=mac_checks
    state.setdefault('timing',{})['quadrature_seconds']=time.perf_counter()-started


def guard_tail_checks(state):
    """Complete the cheap class-gap control for the first five recovered points.

    Those roots were already obtained in class windows. This checks only the
    unsampled upper gap of a class, up to the existing guard, without replacing
    accepted frequencies or repeating a root search from zero.
    """
    for group in state['points'].values():
        if group['status']!='CONFIRMED' or 'common_class_upper' in group:continue
        pid=point_id(group['kappa'],group['beta_deg'])
        search=next(r for r in state['root_searches'] if r['point_id']==pid)
        context=PointMatrices(group['kappa'],group['beta_deg'],group['B'])
        upper=group['roots'][6]['Omega']+.1;records=[];started=time.perf_counter()
        for record in search['records']:
            left=record['windows'][-1][1]+.1;parity=record['parity']
            if left>=upper:continue
            candidates=old.pilot.roots._scan_candidates(lambda w:context.block(w,parity),FS,old.pilot.policy(left,upper),
                case_id=pid,builder_id=f'EB_reflection_{parity:+d}',scan_id='CLASS_GUARD_GAP',points=9,phases=(0.,))[0]
            CALLS['detector_calls']+=1
            bad=[c for c in candidates if c.accepted or old.pilot.suspicious(c)]
            records.append(dict(parity=parity,interval=[left,upper],candidates=[old.pilot.candidate_record(c) for c in candidates]))
            if bad:raise ValueError(f'CLASS_GUARD_GAP_UNRESOLVED: {pid}')
        group.update(common_class_upper=upper,B=context.count,full_B=group['full_B']+context.full_count,
            symmetry_B=group['symmetry_B']+context.block_count,transfer_expm=group['transfer_expm']+context.transfers.expm_calls)
        search['guard_gap_check']=dict(records=records,seconds=time.perf_counter()-started)
    save(state)


def curve(rows,kappa,branch):
    selected=sorted([r for r in rows if float(r['kappa'])==kappa and r['branch_id']==branch],key=lambda r:float(r['beta_deg']))
    return ([float(r['beta_deg']) for r in selected],
        [float(r['Lambda']) if r['tracking_status'] in ('TRACKED','SEED_CONFIRMED') and r['Lambda'] not in (None,'') else np.nan for r in selected])


def render_event_shapes(plt):
    """Display saved physical states only; no propagation or normalization."""
    path=OUTPUT/'tracking_diagnostics.json'
    if not path.exists() or not (OUTPUT/'shapes.npz').exists():return 0
    events=json.loads(path.read_text(encoding='utf-8')).get('events',[])
    count=0
    with np.load(OUTPUT/'shapes.npz',allow_pickle=False) as saved:
        for kappa in PLOT_KAPPAS:
            selected=[e for e in events if e['kappa']==kappa and e['classification']=='CROSSING_SUPPORTED']
            if not selected:continue
            fig,axes=plt.subplots(len(selected),3,figsize=(11,7.8),squeeze=False)
            fig.subplots_adjust(left=.08,right=.98,top=.92,bottom=.17,hspace=.45,wspace=.18)
            handles={}
            for row,event in enumerate(selected):
                samples=[event['traces'][0],min(event['traces'],key=lambda t:abs(t['difference'])),event['traces'][1]]
                largest=max(np.max(np.linalg.norm(saved[t[key]+'__states'][:,:,:2],axis=2)) for t in samples for key in ('shape_a','shape_b'))
                amplification=.14/largest
                for col,(title,sample) in enumerate(zip(('до','вблизи','после'),samples)):
                    ax=axes[row,col];beta=np.deg2rad(sample['beta_deg'])
                    tangents=[np.array([1.,0.]),np.array([-np.cos(beta),-np.sin(beta)])]
                    normals=[np.array([0.,-1.]),np.array([-np.sin(beta),np.cos(beta)])]
                    for key,branch,style in [('shape_a',event['branch_a'],'-'),('shape_b',event['branch_b'],'--')]:
                        states=saved[sample[key]+'__states'];xi=np.linspace(0,1,states.shape[1])
                        for arm in range(2):
                            center=(xi[:,None]-1)*tangents[arm]
                            if key=='shape_a':ax.plot(center[:,0],center[:,1],color='.72',lw=.8)
                            displacement=states[arm,:,0,None]*tangents[arm]+states[arm,:,1,None]*normals[arm]
                            shifted=center+amplification*displacement
                            line,=ax.plot(shifted[:,0],shifted[:,1],color=BRANCH_COLORS[branch],ls=style,lw=1.5)
                            handles[branch]=line
                    ax.set_aspect('equal',adjustable='box');ax.set_xlim(-1.2,1.2);ax.set_ylim(-.25,1.2)
                    ax.set_title(f"{title}: β={sample['beta_deg']:.6f}°",fontsize=10)
                    ax.grid(alpha=.15);ax.tick_params(labelsize=8)
                    if col==0:ax.set_ylabel(event['branch_a']+' / '+event['branch_b']+'\nY/l',fontsize=9)
                    if row==len(selected)-1:ax.set_xlabel('X/l',fontsize=9)
            fig.suptitle(f'κθ={kappa}: формы двух независимых классов у пересечений',fontsize=13)
            names=sorted(handles);fig.legend([handles[k] for k in names],names,loc='lower center',ncol=6,frameon=False,fontsize=9,bbox_to_anchor=(.5,.045))
            fig.text(.5,.012,'Сохранённые массово-нормированные формы; общий масштаб смещения в каждой строке',ha='center',fontsize=9)
            for ext in ('png','pdf'):fig.savefig(OUTPUT/f'eb_spring_crossing_shapes_k{kappa}.{ext}',dpi=300)
            plt.close(fig);count+=1
    return count


def render():
    before=dict(CALLS);started=time.perf_counter()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with (OUTPUT/'tracked_branches.csv').open(encoding='utf-8',newline='') as stream:rows=list(csv.DictReader(stream))
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':12,'text.usetex':False,'pdf.fonttype':42}):
        for kappa in PLOT_KAPPAS:
            fig,ax=plt.subplots(figsize=(9,5.8));fig.subplots_adjust(bottom=.25,left=.09,right=.98,top=.91)
            for j,(branch,color) in enumerate(BRANCH_COLORS.items()):
                x,y=curve(rows,kappa,branch);ax.plot(x,y,color=color,label=branch,lw=1.6,linestyle='-' if j%2==0 else '--')
            ax.set(xlim=(0,90),xlabel='β, °',ylabel='Λ',title=f'κθ={kappa}');ax.grid(alpha=.2)
            ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,frameon=False)
            fig.text(.5,.015,'Ветви первых шести мод при β=0, продолженные по формам',ha='center',fontsize=10)
            for ext in ('png','pdf'):fig.savefig(OUTPUT/f'eb_spring_tracked_k{kappa}.{ext}',dpi=300)
            plt.close(fig)
        shape_figures=0 if PLOT_KAPPAS==[0] else render_event_shapes(plt)
    assert dict(CALLS)==before
    info=dict(seconds=time.perf_counter()-started,matrix_calls=0,root_calls=0,shape_calls=0,tracking_calls=0,
        matplotlib=matplotlib.__version__,tracked_csv_sha256=sha(OUTPUT/'tracked_branches.csv'),main_figures=len(PLOT_KAPPAS),event_shape_figures=shape_figures)
    path=OUTPUT/'run_manifest.json';manifest=json.loads(path.read_text()) if path.exists() else {}
    if 'render' in manifest:manifest.setdefault('render_history',[]).append(manifest['render'])
    manifest['render']=info;atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print('plot-only',info)


def write_manifest(state):
    path=OUTPUT/'run_manifest.json';manifest=json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
    points=list(state['points'].values());roots=[r for p in points for r in p['roots']]
    preflight=state.get('preflight',{});calls=Counter()
    for invocation in state.get('invocations',[]):calls.update(invocation['calls'])
    manifest.update(policy=dict(name='frequency-map-v1',calculation_mode='fast_plot',spectrum_semantics='tracked_branches',
        seed='first six physical EB modes at beta=0, independently for each kappa',PLOT_KAPPAS=PLOT_KAPPAS,
        BASE_grid_tenths=old.grid_tenths(),initial_candidate_positions=7,position_cap=10,
        guard_role='candidate pool control above the highest tracked position',strict_audit_default=False,
        identity='mass-MAC global bijection, fixed local material coordinates and reflection classes',criteria=CRITERIA),
        geometry=old.contract()['geometry'],arm=old.contract()['arm'],normalization=old.contract()['normalization'],
        joint_states=[dict(kappa_theta=k,k_theta=k*ARM.D/ARM.L,mode='SPRING') for k in PLOT_KAPPAS],
        state_order=list(modes.eb.STATE_ORDER),shape_grid=dict(material_coordinate='xi=x/L on each inward arm',nodes=129,quadrature='composite Simpson',mass_components=['u','w']),
        root_gates={name:old.LIMITS[name] for name in ('sigma_ratio','rank_rtol','physical_residual','compatibility','null_residual')},
        source_version=state['source_version'],original_sorted_files_sha256=state['source_hashes'],
        source_data_audit=state['source_data_audit'],invocations=state.get('invocations',[]),
        code_sha256={name:sha(ROOT/name) for name in (
            'scripts/lib/inplane_rotational_spring_eb_modes.py',
            'scripts/analysis/laminated_beams/track_inplane_rotational_spring_eb.py')},
        counts=dict(BASE_points=sum(any(abs(p['beta_deg']-t/10)<1e-10 for t in old.grid_tenths()) for p in points),
            added_points=len(state['added_points']),saved_frequencies=len(roots),saved_whole_structure_modes=len(roots),
            reused_frequencies=sum(r['source']=='REUSED_SORTED_ROOT' for r in roots),new_frequencies=sum(r['source']=='NEW_CHARACTERISTIC_ROOT' for r in roots),
            tracking_status=dict(Counter(r['tracking_status'] for r in state.get('tracked_rows',[]))),
            B_including_preflight=sum(p['B'] for p in points)+preflight.get('B',0),
            B_full=sum(p['full_B'] for p in points)+preflight.get('B',0),B_symmetry=sum(p['symmetry_B'] for p in points),
            max_B_per_point=max((p['B'] for p in points),default=0),
            constraint_matrix=sum(p.get('constraint_matrix',0) for p in points),
            transfer_expm=sum(p['transfer_expm'] for p in points)+preflight.get('transfer_expm',0),
            verification_expm=preflight.get('verification_expm',0),
            shape_reconstructions=calls['shape_reconstructions']+preflight.get('shape_reconstructions',0),
            analytic_arm_evaluations=calls['analytic_arm_evaluations']+calls['quadrature_arm_evaluations']+preflight.get('analytic_arm_evaluations',0),
            tracking_recovery_attempts=len(state['tracking_attempts']),local_accuracy_attempts=sum(e.get('additional_accuracy_attempts',0) for e in state['events']),
            accuracy_repeat_frequencies=sum(len(e.get('local_repeats',[])) for e in state['events']),
            crossing_status=dict(Counter(e['classification'] for e in state['events']))),
        time=dict(compute_invocations_seconds=sum(i['seconds'] for i in state.get('invocations',[])),
            preflight_seconds=preflight.get('seconds',0),
            new_root_points_seconds=sum(p.get('seconds',0) for p in points),
            shape_reconstruction_seconds=calls['shape_seconds'],
            note='Invocation time includes checkpoint/NPZ I/O; shape time is a subset, not additive. Rendering is separate.'),
        limitations=[
            'Finite local evidence for six crossings of independent reflection classes; no global spectral certification.',
            'No exact sampled double root is claimed. Local agreement 1e-8 is a numerical interpretation level, not a rigorous error bound.',
            'No resolved avoided crossing found among the tracked six under the declared close trigger; not an absence theorem.',
            'Old sorted tables and their qualifications remain unchanged; derivative chat illustrations are not evidence.',
            'No RLB, viscosity, FEM/Ritz, foreign physical model builders or kappa continuation were executed.'])
    if PLOT_KAPPAS==[0]:
        manifest['limitations']=[
            'Exact hinge; old kappa=1,100 roots, modes and six events are read-only.',
            'Finite angular tracking and endpoint degeneracies, not a theorem excluding all other crossings.',
            '129-node mass quadrature reused from the preceding validated workflow; no new global grid study.',
            'No RLB, viscosity, FEM/Ritz or spectrum interpolation. Original sorted qualifications are retained.']
        manifest['endpoint_checks']=state.get('endpoint_checks',[])
    atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def hinge_endpoint_checks(state,shapes):
    """Check both full-kernel projections at the exact endpoint, without roots."""
    if state.get('endpoint_checks'):return
    started=time.perf_counter();point=state['points'][point_id(0,90)]
    context=PointMatrices(0,90,point['B']);checks=[]
    equivalence=modes.symmetry_equivalence(np.pi/2,0.,ARM)
    # symmetry_equivalence builds one physical joint constraint matrix.
    context.tick('constraint_matrix')
    point['constraint_matrix']=point.get('constraint_matrix',0)+1
    R=np.diag([1.,-1.,-1.]);L=np.diag([1.,1.,-1.])
    for frequency in dict.fromkeys(r['Omega'] for r in point['roots']):
        group=[r for r in point['roots'] if r['Omega']==frequency]
        if len(group)!=2:continue
        vectors=np.array([shapes[r['shape_key']+'__vector'] for r in group])
        gram=vectors.conj()@vectors.T
        plus=context.block(frequency/FS,1);minus=context.block(frequency/FS,-1)
        block_error=float(np.linalg.norm(minus-L@plus@R)/np.linalg.norm(minus))
        gram_error=float(np.linalg.norm(gram-np.eye(2),ord=np.inf))
        selected=[r for r in state['tracked_rows'] if r['beta_deg']==90 and r['shape_key'] in [x['shape_key'] for x in group]]
        nullities=[state['shape_checks'][r['shape_key']]['detected_nullity'] for r in group]
        ok=(all(r['root_status']=='CONFIRMED' for r in group) and nullities==[2,2]
            and {r['symmetry_class'] for r in group}=={-1,1} and gram_error<=1e-10
            and block_error<=1e-12 and equivalence['row_error']<=1e-12 and equivalence['row_rank']==6
            and all(r['tracking_status']=='TRACKED' and r['MAC']>=CRITERIA['mac'] for r in selected))
        cluster=f'k0_endpoint_{group[0]["current_sorted_position"]:02d}_{group[1]["current_sorted_position"]:02d}'
        for r in group+selected:r['cluster_id']=cluster
        record=dict(cluster_id=cluster,Omega=frequency,Lambda=group[0]['Lambda'],positions=[r['current_sorted_position'] for r in group],
            shape_keys=[r['shape_key'] for r in group],symmetry_classes=[r['symmetry_class'] for r in group],nullities=nullities,
            mass_gram=gram.tolist(),mass_gram_error=gram_error,block_equivalence_error=block_error,
            constraint_equivalence=equivalence,one_sided_MAC=[r['MAC'] for r in selected],
            classification='ENDPOINT_DEGENERACY_SUPPORTED' if ok else 'UNRESOLVED_CLOSE_CLUSTER',
            role='ROOT' if selected else 'GUARD')
        checks.append(record)
        if selected:
            branches=sorted(selected,key=lambda r:r['branch_id'])
            left=[r for r in state['tracked_rows'] if r['beta_deg']==89 and r['branch_id'] in [x['branch_id'] for x in branches]]
            left.sort(key=lambda r:r['branch_id'])
            state['events'].append(dict(event_id=cluster,kappa=0,branch_a=branches[0]['branch_id'],branch_b=branches[1]['branch_id'],
                classification=record['classification'],beta_left=89.,beta_right=90.,endpoint_beta=90.,
                Omega_a_left=left[0]['Omega'],Omega_b_left=left[1]['Omega'],Omega_a_right=frequency,Omega_b_right=frequency,
                difference_left=left[0]['Omega']-left[1]['Omega'],difference_right=0.,minimum_sampled_gap=0.,
                symmetry_a=branches[0]['symmetry_class'],symmetry_b=branches[1]['symmetry_class'],
                evidence='Exact endpoint block equivalence, full nullity 2, independent mass-orthogonal reflection classes and one-sided shape continuation; no beta>90 claim'))
        if not ok:
            for r in selected:r['tracking_status']='TRACKING_AMBIGUOUS'
    point['B']=context.count;point['symmetry_B']+=context.block_count
    point['transfer_expm']+=context.transfers.expm_calls
    state['endpoint_checks']=checks;state.setdefault('timing',{})['endpoint_seconds']=time.perf_counter()-started
    write_csv(OUTPUT/'crossing_events.csv',[{k:v for k,v in e.items() if k not in ('traces','local_repeats','initial_bracket')} for e in state['events']])
    write_csv(OUTPUT/'tracked_branches.csv',state['tracked_rows']);save(state)


def read_csv(path):
    with path.open(encoding='utf-8',newline='') as stream:return list(csv.DictReader(stream))


def comparison_inputs():
    """Read-only old results: validation does not relabel or reconstruct them."""
    assert_protected_sources()
    directories={0:HINGE_OUTPUT,1:LEGACY_TRACKED,100:LEGACY_TRACKED}
    manifests={k:json.loads((path/'run_manifest.json').read_text(encoding='utf-8')) for k,path in directories.items()}
    for k,manifest in manifests.items():
        for field in ('geometry','arm','normalization'):
            if manifest[field]!=old.contract()[field]:raise ValueError('Incompatible comparison '+field)
        if manifest['source_version']['physics_hashes']!=old.contract()['physics_sha256']:raise ValueError('Physical source mismatch')
        if manifest['shape_grid']!=manifests[1]['shape_grid']:raise ValueError('Incompatible mass metric')
    rows={k:[r for r in read_csv(path/'tracked_branches.csv') if float(r['kappa'])==k] for k,path in directories.items()}
    return directories,manifests,rows


def seed_assignment(reference,candidates):
    columns,mac,margins=modes.assign([r['vector'] for r in reference],[r['vector'] for r in candidates],
        [r['symmetry_class'] for r in reference],[r['symmetry_class'] for r in candidates])
    records=[]
    for i,j in enumerate(columns):
        eligible=[q for q,r in enumerate(candidates) if r['symmetry_class']==reference[i]['symmetry_class']]
        good=(mac[i,j]>=CRITERIA['mac'] and margins[i]>=CRITERIA['margin'] and j in eligible)
        records.append(dict(comparison_mode_id=f'comparison_mode_{i+1:02d}',source_branch_id=candidates[j]['branch_id'],
            seed_beta=0,symmetry_class=reference[i]['symmetry_class'],matching_method='DIRECT_MASS_MAC_GLOBAL_ASSIGNMENT',
            MAC=float(mac[i,j]),margin=float(margins[i]),status='CONFIRMED' if good else 'SEED_MAPPING_AMBIGUOUS',
            competitors=json.dumps([dict(branch_id=candidates[q]['branch_id'],MAC=float(mac[i,q])) for q in eligible]),
            reference_branch_id=reference[i]['branch_id']))
    return records


def prepare_seed_mapping():
    directories,manifests,rows=comparison_inputs();COMPARISON_OUTPUT.mkdir(exist_ok=True)
    path=COMPARISON_OUTPUT/'seed_mode_mapping.csv'
    sources={str((p/name).relative_to(ROOT)):sha(p/name) for p in set(directories.values())
             for name in ('tracked_branches.csv','verified_roots.csv','shapes.npz','run_manifest.json')}
    manifest_path=COMPARISON_OUTPUT/'comparison_manifest.json'
    if path.exists():
        saved=json.loads(manifest_path.read_text(encoding='utf-8'))
        # Rendering changes the hinge manifest only; scientific inputs must match.
        if any(saved['input_sha256'].get(n)!=h for n,h in sources.items() if not n.endswith('run_manifest.json')):
            raise ValueError('Seed mapping sources changed')
        print('seed mapping reused');return
    started=time.perf_counter();seeds={}
    for k,path0 in directories.items():
        selected=sorted([r for r in rows[k] if float(r['beta_deg'])==0],key=lambda r:int(r['seed_sorted_position']))
        if len(selected)!=6 or any(r['root_status']!='CONFIRMED' or r['tracking_status']!='SEED_CONFIRMED' for r in selected):
            raise ValueError('Incomplete initial six modes')
        with np.load(path0/'shapes.npz',allow_pickle=False) as saved:
            seeds[k]=[dict(branch_id=r['branch_id'],symmetry_class=int(r['symmetry_class']),vector=saved[r['shape_key']+'__vector']) for r in selected]
    mapping=[]
    for k in (0,1,100):
        for record in seed_assignment(seeds[1],seeds[k]):
            record['kappa']=k;record['source']=str(directories[k].relative_to(ROOT))
            if k==1:record['matching_method']='REFERENCE_IDENTITY_CHECKED_MASS_MAC'
            mapping.append(record)
    write_csv(COMPARISON_OUTPUT/'seed_mode_mapping.csv',mapping)
    manifest=dict(reference=dict(kappa=1,beta_deg=0),path='reference seed -> kappa at beta=0 -> frozen source angular continuation',
        criteria=dict(MAC=CRITERIA['mac'],margin=CRITERIA['margin'],symmetry_required=True,frequency_in_cost=False),
        input_sha256=sources,source_versions={str(k):m['source_version'] for k,m in manifests.items()},
        addition_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),code_status='working-tree version',
        addition_environment=dict(executable=sys.executable,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        mapping_code_sha256=sha(Path(__file__)),
        matching_seconds=time.perf_counter()-started,mapping_counts=dict(Counter(r['status'] for r in mapping)),
        new_kappa_seed_values=[],new_roots=0,shape_reconstructions=0,old_angular_tracking_calls=0,
        limitations=['Frozen seed correspondence; no path independence in the two-dimensional parameter plane claimed.',
                     'Old six crossing events are reused, not localized again.'])
    atomic(manifest_path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print('seed mapping',manifest['mapping_counts'],'min MAC',min(r['MAC'] for r in mapping),'min margin',min(r['margin'] for r in mapping))


def mapped_rows(mapping,source_rows):
    """Frozen mapping only: retain every source frequency and sorted position."""
    result=[]
    if len({(float(r['kappa']),r['source_branch_id']) for r in mapping})!=len(mapping):raise ValueError('Non-bijective seed mapping')
    for match0 in mapping:
        k=float(match0['kappa'])
        for r in source_rows[k]:
            if r['branch_id']!=match0['source_branch_id']:continue
            result.append(dict(comparison_mode_id=match0['comparison_mode_id'],source_branch_id=r['branch_id'],kappa=r['kappa'],beta_deg=r['beta_deg'],
                **{key:r[key] for key in ('omega','Omega','Lambda','current_sorted_position','symmetry_class','root_status','tracking_status')},
                mapping_status=match0['status'],grid_role=r.get('grid_role') or ('BASE' if any(abs(float(r['beta_deg'])-t/10)<1e-10 for t in old.grid_tenths()) else 'ADDED'),
                source=match0['source'].replace('\\','/')+'/tracked_branches.csv',source_shape_key=r['shape_key']))
    return sorted(result,key=lambda r:(r['comparison_mode_id'],float(r['kappa']),float(r['beta_deg'])))


def prepare_comparison_table():
    _,_,rows=comparison_inputs();mapping=read_csv(COMPARISON_OUTPUT/'seed_mode_mapping.csv')
    result=mapped_rows(mapping,rows)
    write_csv(COMPARISON_OUTPUT/'comparison_branches.csv',result)
    path=COMPARISON_OUTPUT/'comparison_manifest.json';manifest=json.loads(path.read_text(encoding='utf-8'))
    manifest['comparison_rows']=len(result)
    manifest['table_sha256']={name:sha(COMPARISON_OUTPUT/name) for name in ('seed_mode_mapping.csv','comparison_branches.csv')}
    manifest['source_grid_counts']={str(k):len({r['beta_deg'] for r in rs}) for k,rs in rows.items()}
    atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print('comparison rows',len(result))


def continue_seed_mapping():
    """One triggered beta=0 bridge at kappa=10, using the existing root path.

    No old root is searched, no old angular continuation is run. New modes
    and the initial direct failure are retained for reproducibility.
    """
    directories,_,source_rows=comparison_inputs()
    path=COMPARISON_OUTPUT/'comparison_manifest.json';manifest=json.loads(path.read_text(encoding='utf-8'))
    if manifest.get('seed_continuation'):
        manifest['seed_continuation']['ROOT_GUARD']=seed_guard_qualification(manifest['seed_continuation']['root_state']['points']['k10_b0'])
        atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        print('seed continuation reused');return
    mapping=read_csv(COMPARISON_OUTPUT/'seed_mode_mapping.csv')
    unresolved=[r for r in mapping if r['status']!='CONFIRMED']
    if not unresolved:print('no ambiguous seed');return
    if any(float(r['kappa'])!=100 for r in unresolved):raise ValueError('This local bridge addresses only kappa=1 to 100')
    started=time.perf_counter();before=Counter(CALLS);seeds={}
    for k in (1,100):
        rows=sorted([r for r in source_rows[k] if float(r['beta_deg'])==0],key=lambda r:int(r['seed_sorted_position']))
        with np.load(directories[k]/'shapes.npz',allow_pickle=False) as saved:
            seeds[k]=[dict(branch_id=r['branch_id'],symmetry_class=int(r['symmetry_class']),vector=saved[r['shape_key']+'__vector']) for r in rows]
    roots0=[r for r in read_csv(LEGACY_TRACKED/'verified_roots.csv') if float(r['beta_deg'])==0 and float(r['kappa']) in (1,100)]
    windows={}
    for eta in (1,-1):
        ends=[sorted(float(r['Omega']) for r in roots0 if float(r['kappa'])==k and int(r['symmetry_class'])==eta) for k in (1,100)]
        if len(ends[0])!=len(ends[1]):raise ValueError('Seed candidate pool requires a local review')
        windows[eta]=old.merge_windows([(max(1e-8,min(a,b)-.25),max(a,b)+.25) for a,b in zip(*ends)])
    seed_state=dict(points={},added_points=[],root_searches=[]);shapes={}
    point=search_point(seed_state,shapes,10,0.,'AMBIGUOUS_DIRECT_SEED_K100_MODE05',windows_override=windows,persist=False)
    continuation=dict(kappa_path=[1,10,100],beta_deg=0.,trigger=unresolved,root_state=seed_state,
                      available_saved_intermediate='kappa=.1 lies outside [1,100]; not useful for this bridge')
    if point['status']=='CONFIRMED':
        middle=[dict(branch_id=f"seed10_{r['current_sorted_position']:02d}",symmetry_class=r['symmetry_class'],
                     vector=shapes[r['shape_key']+'__vector']) for r in point['roots']]
        first=seed_assignment(seeds[1],middle)
        selected=[next(r for r in middle if r['branch_id']==record['source_branch_id']) for record in first]
        second=seed_assignment(selected,seeds[100]);continuation.update(first_step=first,second_step=second)
        for record,a,b in zip([r for r in mapping if float(r['kappa'])==100],first,second):
            record.update(direct_MAC=record['MAC'],direct_margin=record['margin'],direct_status=record['status'],
                direct_source_branch_id=record['source_branch_id'],direct_competitors=record['competitors'],
                source_branch_id=b['source_branch_id'],matching_method='MASS_MAC_CONTINUATION_1_10_100_AT_BETA0',
                MAC=min(a['MAC'],b['MAC']),margin=min(a['margin'],b['margin']),competitors=json.dumps([a,b]),
                status='CONFIRMED' if a['status']==b['status']=='CONFIRMED' else 'SEED_MAPPING_AMBIGUOUS')
    continuation.update(seconds=time.perf_counter()-started,calls=dict(Counter(CALLS)-before),
                        ROOT_GUARD=seed_guard_qualification(point),
                        constraints='one beta=0 point, at most 6000 full/block matrices, no old root replacement')
    with (COMPARISON_OUTPUT/'seed_shapes.npz.tmp').open('wb') as stream:np.savez_compressed(stream,**shapes)
    (COMPARISON_OUTPUT/'seed_shapes.npz.tmp').replace(COMPARISON_OUTPUT/'seed_shapes.npz')
    manifest['seed_continuation']=continuation;manifest['new_kappa_seed_values']=[10]
    manifest['new_roots']=len(point['roots']);manifest['shape_reconstructions']=continuation['calls'].get('shape_reconstructions',0)
    manifest['mapping_counts']=dict(Counter(r['status'] for r in mapping))
    manifest['continuation_code_sha256']=sha(Path(__file__))
    write_csv(COMPARISON_OUTPUT/'seed_mode_mapping.csv',mapping)
    atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print('seed continuation',manifest['mapping_counts'],'B',point['B'],'seconds',continuation['seconds'])


def seed_guard_qualification(point):
    """Expose the existing detector warning; do not change root gates/statuses."""
    roots=point['roots'];warnings=point.get('guard_warnings',[])
    target_ok=len(roots)>=6 and all(r['root_status']=='CONFIRMED' for r in roots[:6])
    gap=min((r['interval'][0]-roots[5]['Omega'] for r in warnings),default=None) if len(roots)>=6 else None
    separated=not warnings or (gap is not None and gap>old.LIMITS['guard_margin_Omega'])
    return dict(original_helper_status=point['status'],ROOT='CONFIRMED' if target_ok and separated else 'QUALIFIED',
                GUARD='QUALIFIED_DETECTOR_WARNING' if warnings else 'CONFIRMED',
                warning_gap_above_sixth=gap,guard_margin_to_right=point.get('common_class_upper',0)-roots[-1]['Omega'] if roots else None,
                source_warning_retained=warnings,
                interpretation='Only six seed modes are used; no additional guard search or claimed warning resolution.')


def comparison_curve(rows,kappa,mode):
    selected=sorted([r for r in rows if float(r['kappa'])==kappa and r['comparison_mode_id']==mode],key=lambda r:float(r['beta_deg']))
    return ([float(r['beta_deg']) for r in selected],
        [float(r['Lambda']) if r['mapping_status']=='CONFIRMED' and r['root_status']=='CONFIRMED' and r['tracking_status'] in ('TRACKED','SEED_CONFIRMED') and r['Lambda'] not in ('',None) else np.nan for r in selected])


def render_comparison():
    before=dict(CALLS);started=time.perf_counter()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=read_csv(COMPARISON_OUTPUT/'comparison_branches.csv')
    mapping=read_csv(COMPARISON_OUTPUT/'seed_mode_mapping.csv')
    if len(mapping)!=18:raise ValueError('Expected eighteen seed mappings')
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':12,'text.usetex':False,'pdf.fonttype':42}):
        for j in range(1,7):
            fig,ax=plt.subplots(figsize=(9,5.8));fig.subplots_adjust(bottom=.24,left=.10,right=.98,top=.91)
            for k,style in COMPARISON_STYLES.items():
                x,y=comparison_curve(rows,k,f'comparison_mode_{j:02d}')
                marker,every=COMPARISON_MARKERS[k]
                ax.plot(x,y,label=f'κθ={k}',lw=1.9,marker=marker,markevery=every,markersize=4,
                        markerfacecolor='white',markeredgecolor=style['color'],**style)
            ax.set(xlim=(0,90),xlabel='β, °',ylabel='Λ',title=f'Мода {j}');ax.margins(y=.06);ax.grid(alpha=.2)
            ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,frameon=False)
            fig.text(.5,.015,'Общий seed: β=0, κθ=1; соответствие начальных форм и продолжение по β',ha='center',fontsize=10)
            for ext in ('png','pdf'):fig.savefig(COMPARISON_OUTPUT/f'eb_tracked_mode{j:02d}_kappa_comparison.{ext}',dpi=300)
            plt.close(fig)
    assert dict(CALLS)==before
    info=dict(seconds=time.perf_counter()-started,matrix_calls=0,root_calls=0,shape_calls=0,tracking_calls=0,
        seed_matching_calls=0,matplotlib=matplotlib.__version__,figures=6)
    path=COMPARISON_OUTPUT/'comparison_manifest.json';manifest=json.loads(path.read_text(encoding='utf-8'))
    if 'render' in manifest:manifest.setdefault('render_history',[]).append(manifest['render'])
    manifest['render']=info;atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print('comparison plot-only',info)


def main():
    global PLOT_KAPPAS,OUTPUT,LEGACY_TRACKED,HINGE_OUTPUT,COMPARISON_OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('compute','resume','plot-only','seed-map','seed-continue','comparison-table','comparison-plot-only'),required=True)
    parser.add_argument('--kappas',nargs='+',type=int,choices=(0,1,100),default=[1,100])
    parser.add_argument('--output',type=Path)
    parser.add_argument('--legacy-tracked',type=Path,default=LEGACY_TRACKED)
    parser.add_argument('--hinge-output',type=Path,default=HINGE_OUTPUT)
    parser.add_argument('--comparison-output',type=Path,default=COMPARISON_OUTPUT)
    parser.add_argument('--stage',choices=('roots','track','events','checks','all'),default='all')
    args=parser.parse_args()
    PLOT_KAPPAS=args.kappas
    LEGACY_TRACKED=args.legacy_tracked.resolve();HINGE_OUTPUT=args.hinge_output.resolve();COMPARISON_OUTPUT=args.comparison_output.resolve()
    OUTPUT=(args.output or (HINGE_OUTPUT if PLOT_KAPPAS==[0] else OUTPUT)).resolve()
    if 0 in PLOT_KAPPAS and PLOT_KAPPAS!=[0]:parser.error('The hinge addition uses --kappas 0 alone; old angular results are read-only.')
    if PLOT_KAPPAS==[0]:
        CRITERIA['max_added_points']=50
        if OUTPUT in (LEGACY_TRACKED,old.OUTPUT.resolve()):parser.error('Hinge output must not overwrite a source directory')
    if COMPARISON_OUTPUT in (LEGACY_TRACKED,HINGE_OUTPUT,old.OUTPUT.resolve()):parser.error('Comparison output must be separate')
    assert_protected_sources()
    if args.mode=='seed-map':prepare_seed_mapping();return
    if args.mode=='seed-continue':continue_seed_mapping();return
    if args.mode=='comparison-table':prepare_comparison_table();return
    if args.mode=='comparison-plot-only':render_comparison();return
    if args.mode=='plot-only':render();return
    state,shapes=load();started=time.perf_counter()
    import_saved(state,shapes)
    for kappa in PLOT_KAPPAS:
        for tick in old.grid_tenths():
            if point_id(kappa,tick/10) not in state['points']:
                search_point(state,shapes,kappa,tick/10,'MISSING_BASE_FROM_SORTED_MAP')
    save(state,shapes)
    if args.stage not in ('roots','checks') and not state.get('tracked_rows'):track(state,shapes)
    if args.stage in ('events','all'):
        count_before_events=len(state['points'])
        crossing_events(state,shapes)
        if PLOT_KAPPAS==[0]:
            if len(state['points'])>count_before_events:track(state,shapes,allow_new=False)
            hinge_endpoint_checks(state,shapes)
        else:
            quadrature_checks(state,shapes)
            track(state,shapes,allow_new=False)
    if args.stage in ('checks','all'):guard_tail_checks(state)
    state.setdefault('invocations',[]).append(dict(seconds=time.perf_counter()-started,calls=dict(CALLS),stage=args.stage,runner_sha256=sha(Path(__file__))))
    save(state)
    write_manifest(state)
    assert_protected_sources()
    print('compute',state['invocations'][-1],flush=True)


if __name__=='__main__':main()
