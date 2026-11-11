"""HS-D2 geometry preset and saved-data analysis for the HS-D1 CLI.

No physical assembly is defined here. Study.point/value and both leading
solvers remain those used by HS-D1. Full spectra are read-only references.
"""
import copy
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scripts.lib import two_springs_horizontal_composite as d1

ROOT=d1.ROOT
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_composite_geometry_validation.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_composite_geometry_validation'
read,atomic,sha,require,csvwrite=d1.read,d1.atomic,d1.sha,d1.require,d1.csvwrite
rows_csv=d1.outer.rows_csv
RIGID=d1.RIGID


def output_folder():
    return d1.output_folder(CONFIG,BASE)


def point_identity(kind,k,mu,param,ks=None):
    return (kind,float(k),float(mu),param if param==RIGID else float(param),ks)


def mechanical_identity(case):
    """Exact leading-model identity; aspect is retained separately in outputs."""
    return tuple(case[n] for n in ['k_joint','mu','epsilon','chi_s','kappa_s','role'])


def expected_ids():
    result=[]
    for regime in ['A1','A2','A3','A4']:
        result += [f'B0_{regime}',f'B1_{regime}_e0.0025',f'B1_{regime}_e0.006']
        result += [f'B2_{regime}_m{mu}' for mu in ['-0.1','-0.01','0','0.1']]
        result += [f'B3_{regime}_a2']
        if regime in ['A2','A4']:result += [f'B4_{regime}']
    return result


class GeometryStudy(d1.Study):
    def __init__(self):
        super().__init__(CONFIG,BASE)
        require([c['id'] for c in self.cfg['cases']]==expected_ids(),'Exactly the authorized 34 cases')
        self.source=ROOT/self.cfg['full_source'];self.d1_source=ROOT/self.cfg['composite_source']
        self.old_d1={}
        for path in (self.d1_source/'points').glob('*.json'):
            r=read(path);key=point_identity(r['kind'],r['k_joint'],r['mu'],r['parameter'],r['kappa_s'])
            require(key not in self.old_d1,'Duplicate source spectral identity')
            self.old_d1[key]=(r,path)
        version=self.folder/'workflow_versions'/(sha(__file__)+'.txt')
        if not version.exists():version.write_bytes(Path(__file__).read_bytes())
        for r in self.records.values():
            require(r['status']=='PASS','Resume status')
            require(sha(self.points/(r['point_key']+'.npz'))==r['npz_sha256'],'Resume shapes hash')

    def audit(self):
        super().audit()
        geometry=rows_csv(self.source/'geometry.csv');spectra=rows_csv(self.source/'spectra.csv')
        require({r['id'] for r in geometry}==set(expected_ids()),'Geometry IDs')
        require(len(spectra)==self.cfg['full_point_count'],'Full point count')
        full_provenance=[];cases={c['id']:c for c in self.cfg['cases']}
        grid=read(self.source/'grid.json')
        for geom in geometry:
            c=cases[geom['id']]
            for key in ['epsilon','mu','aspect','k_joint']:require(float(geom[key])==c[key],'Geometry '+key)
            require(geom['status']=='PASS','Geometry status')
            require(float(geom['limiter'])<=.1,'Original geometry limiter')
            require(c['kappa_s']==RIGID or c['kappa_s']==c['chi_s']/c['epsilon']**2,'kappa_s mapping')
        for r in spectra:
            c=cases[r['id']]
            for a,b in [('epsilon','epsilon'),('mu','mu'),('aspect','aspect'),('kappa_j','k_joint')]:
                require(float(r[a])==c[b],'Full identity '+a)
            require(r['chi_s']==RIGID if c['chi_s']==RIGID else float(r['chi_s'])==c['chi_s'],'Full chi')
            require(r['kappa_s']==RIGID if c['kappa_s']==RIGID else float(r['kappa_s'])==c['kappa_s'],'Full kappa')
            p=ROOT/r['point_file'];record=read(p)
            require(r['status']=='PASS' and record['status']=='PASS' and record['count']['status']=='PASS','Full accepted prefix')
            require(float(r['omega1'])==record['omega'][0],'Full point/CSV')
            require(len(record['omega'])>=3 and int(r['zero_count'])==0,'Full prefix/zero')
            require(abs(float(r['Omega1'])-float(r['omega1'])/c['epsilon'])<1e-12,'Omega normalization')
            require(abs(float(r['Lambda1'])**2-float(r['Omega1']))<1e-12,'Lambda normalization')
            full_provenance.append(dict(id=r['id'],beta_deg=float(r['beta_deg']),path=r['point_file'],sha256=sha(p),
                epsilon=c['epsilon'],mu=c['mu'],aspect=c['aspect'],k_joint=c['k_joint'],chi_s=c['chi_s'],
                normalization='omega=epsilon*sqrt(zeta); Omega=sqrt(zeta); Lambda=zeta**.25'))
        for c in cases.values():
            angles=sorted(float(r['beta_deg']) for r in spectra if r['id']==c['id'])
            require(angles==self.cfg['angles_deg'][c['id']],'Frozen actual angles')
            require(all(r['beta_deg'] in angles for r in grid[c['id']]),'Original grid membership')
        atomic(self.folder/'full_point_provenance.json',full_provenance)
        for name in ['geometry.csv','grid.json','aspect_diagnostics.json']:
            (self.folder/('source_'+name if name!='geometry.csv' else name)).write_bytes((self.source/name).read_bytes())
        atomic(self.folder/'grid.json',self.cfg['angles_deg'])
        atomic(self.folder/'source_audit.json',dict(status='PASS',full_points=len(spectra),geometries=len(geometry),new_full_solves=0,
            filenames={'HS_B_thresholds':'thresholds.json and thresholds.csv','HS_B_previous_predictions':'predictions.csv',
                'HS_B_threshold_comparison.csv':'absent; not required when authoritative thresholds.json is present',
                'HS_D1_predictions':'predictions.json','HS_D1_comparison':'threshold_comparison.csv'},
            HS_B_leading_reuse_limit='leading_comparison stores frequencies/diagnostics but not shapes or per-mode D; used for cross-checks, not as complete point records',
            cache_rules=self.cfg['reuse_rules']))

    def point(self,kind,k,mu,param,ks=None,role='grid',independent=False):
        identity=point_identity(kind,k,mu,param,ks)
        key=hashlib.sha256(repr(identity).encode()).hexdigest()[:16]
        if key in self.records:return self.records[key]
        if identity in self.old_d1:
            old,path=self.old_d1[identity];r=copy.deepcopy(old)
            require(old['status']=='PASS' and old['point_key']==key,'HS-D1 reuse identity')
            for row in r['orders']:self.validate(row,kind)
            require(r['frequency_convergence']<=self.p['frequency_convergence_tol'],'Reuse p convergence')
            require(r['shape_convergence']<=self.p['shape_convergence_tol'],'Reuse shape convergence')
            src=path.with_suffix('.npz');tmp=self.points/(key+'.tmp.npz');tmp.write_bytes(src.read_bytes());tmp.replace(self.points/(key+'.npz'))
            r.update(origin='reused_HS_D1',prior_origin=old['origin'],source=path.relative_to(ROOT).as_posix(),source_sha256=sha(path),
                source_npz_sha256=sha(src),role=role,npz_sha256=sha(self.points/(key+'.npz')))
            atomic(self.points/(key+'.json'),r);self.records[key]=r;return r
        r=super().point(kind,k,mu,param,ks,role,independent)
        r['npz_sha256']=sha(self.points/(key+'.npz'));atomic(self.points/(key+'.json'),r)
        return r

    def coefficients(self):
        path=self.folder/'coefficients.json'
        if path.exists():return read(path)
        unique=sorted({(c['k_joint'],c['mu']) for c in self.cfg['cases']});result=[]
        prior=read(self.d1_source/'coefficients.json');prior_cfg=read(self.d1_source/'contract.json')['config']
        previous={}
        for r in prior:
            c=next(c for c in prior_cfg['cases'] if c['id']==r['id']);previous[c['k_joint'],c['mu']]=r
        for k,mu in unique:
            inf=self.point('L1',k,mu,RIGID,role='common limit',independent=True)
            zero=self.point('L1',k,mu,0.,role='zero limit',independent=True)
            c0=self.point('C',k,mu,0.,RIGID,role='right limit')
            c90=self.point('C',k,mu,90.,RIGID,role='positivity endpoint')
            orders=[]
            for j,p in enumerate(self.cfg['orders']):
                z=inf['orders'][j]['zeta'][0];zc=c0['orders'][j]['zeta'][0];pm=c0['orders'][j]['per_mode'][0]
                err=abs(math.sqrt(zc/z)-1);require(err<=self.p['independent_frequency_tol'],'Common limit')
                C2=zc*pm['mixed_mass']/2
                bound=min(zc,c90['orders'][j]['zeta'][0])+zero['orders'][j]['zeta'][0]-z
                require(bound>0,'Conditional global composite positivity')
                orders.append(dict(order=p,zeta_infty=z,zeta_C_zero_plus=zc,zeta_C90=c90['orders'][j]['zeta'][0],
                    zeta_L1_zero=zero['orders'][j]['zeta'][0],C2=C2,A=pm['A'],integral_V=-pm['mixed_mass']/(2*pm['A']),
                    common_limit_frequency_error=err,global_composite_lower_bound=bound))
            require(abs(orders[0]['C2']/orders[1]['C2']-1)<=self.p['coefficient_relative_tol'],'C2 order convergence')
            checks=[]
            for beta in self.cfg['diagnostics']['C2_angles_deg']:
                zc=self.point('C',k,mu,beta,RIGID,role='C2 check')['zeta'][0]
                quotient=(zc-orders[-1]['zeta_C_zero_plus'])/math.radians(beta)**2
                checks.append(dict(beta_deg=beta,C2_quotient=quotient))
            rich=(4*checks[-1]['C2_quotient']-checks[-2]['C2_quotient'])/3
            error=rich/orders[-1]['C2']-1;require(abs(error)<=self.cfg['diagnostics']['C2_relative_tol'],'C2 finite difference')
            old=previous.get((k,mu));old_error=None if old is None else max(abs(orders[-1][n]/old['orders'][-1][n]-1) for n in ['zeta_infty','C2'])
            require(old_error is None or old_error<=self.p['coefficient_relative_tol'],'HS-D1 common coefficients')
            full_limit=read(self.source/'leading_comparison'/f'j{k:g}_mu{mu:g}.json')
            require(abs(math.sqrt(orders[-1]['zeta_infty']/full_limit['zeta_infty'])-1)<=self.p['independent_frequency_tol'],'HS-B limit reference')
            result.append(dict(k_joint=k,mu=mu,orders=orders,inner_limit_key=inf['point_key'],outer_limit_key=c0['point_key'],
                C2_checks=checks,C2_Richardson_relative_error=error,HS_D1_coefficient_error=old_error))
        atomic(path,result);return result

    def preflight(self):
        path=self.folder/'preflight.json'
        if path.exists():require(read(path)['status']=='PASS','Preflight status');return
        start=time.perf_counter();self.audit();values=[]
        for name in self.cfg['preflight_ids']:
            case=next(c for c in self.cfg['cases'] if c['id']==name)
            for beta in [.1,90.]:values.append(self.value(case,beta,'preflight'))
        self.coefficients()
        atomic(path,dict(status='PASS',elapsed_seconds=time.perf_counter()-start,values=values,
            finite_tau_not_replaced_by_infinity=True,new_full_solves=0))
        print('PREFLIGHT PASS',self.folder,flush=True)

    def predictions(self,curves):
        path=self.folder/'predictions.json'
        if path.exists():return read(path)
        prior=read(self.d1_source/'predictions.json');old_cases=read(self.d1_source/'contract.json')['config']['cases']
        cache={(mechanical_identity(next(c for c in old_cases if c['id']==r['id'])),r['q']):r for r in prior}
        results=[]
        for case in self.cfg['cases']:
            data=sorted([r for r in curves if r['id']==case['id']],key=lambda r:r['beta_deg'])
            angles=[r['beta_deg'] for r in data];freq=[r['omega_primary'] for r in data]
            for q in self.cfg['thresholds']['q']:
                target=q*freq[-1];key=(mechanical_identity(case),q)
                brackets=d1.crossing_brackets(angles,freq,target);crossings=[];old=cache.get(key)
                def fun(b):return self.value(case,b,'model threshold')['omega_primary']/target-1
                for lo,hi,direction in brackets:
                    existing=next((r for r in old['crossings'] if lo<=r['beta_deg']<=hi and r['direction']==direction),None) if old else None
                    if existing is not None:
                        crossing=copy.deepcopy(existing);crossing['origin']='reused_direct_root';crossing['source_id']=old['id']
                        beta=crossing['beta_deg'];self.value(case,beta,'reused direct threshold')
                        require(abs(fun(beta))<=self.cfg['thresholds']['residual_tol'],'Reused threshold residual')
                        require(fun(crossing['bracket'][0])*fun(crossing['bracket'][1])<0,'Reused bracket signs')
                    else:
                        beta=brentq(fun,lo,hi,xtol=self.cfg['thresholds']['root_xtol_deg'],rtol=self.cfg['thresholds']['root_rtol'])
                        half=self.cfg['thresholds']['bracket_width_deg']/2;left=max(0.,beta-half);right=min(90.,beta+half)
                        fl,fr=fun(left),fun(right);v=self.value(case,beta,'threshold slope');res=fun(beta)
                        require(fl*fr<0 and abs(res)<=self.cfg['thresholds']['residual_tol'],'Threshold bracket/residual')
                        differences=[]
                        for h in self.cfg['diagnostics']['slope_steps_deg']:
                            require(beta>h and beta+h<90,'Central slope domain')
                            derivative=math.log(self.value(case,beta+h,'slope check')['omega_primary']/self.value(case,beta-h,'slope check')['omega_primary'])/(2*h)
                            err=derivative/v['slope_log_per_degree']-1
                            differences.append(dict(h_deg=h,slope_log_per_degree=derivative,relative_error=err))
                        require(abs(differences[-1]['relative_error'])<=self.cfg['diagnostics']['slope_relative_tol'],'Threshold slope')
                        crossing=dict(beta_deg=beta,direction=direction,bracket=[left,right],signs=[fl,fr],residual=res,
                            xi=v['xi'],slope_log_per_degree=v['slope_log_per_degree'],slope_log_per_rad=v['slope_log_per_rad'],
                            frequency_derivative_per_degree=target*v['slope_log_per_degree'],differences=differences,origin='new_direct_root',source_id=None)
                    crossings.append(crossing)
                high0=freq[0]>=target;high=[f>=target for f in freq]
                last_bad=max((i for i,ok in enumerate(high) if not ok),default=-1)
                terminal_grid=angles[last_bad+1] if high[-1] else None
                ups=[r for r in crossings if r['direction']=='up']
                terminal_root=ups[-1]['beta_deg'] if high[-1] and ups else 0. if all(high) else None
                status='PASS' if crossings else 'ALL_SAMPLED_HIGH' if all(high) else 'NO_CROSSING_FOUND'
                row=dict(id=case['id'],q=q,omega_reference=freq[-1],primary=data[0]['primary'],status=status,crossings=crossings,
                    high_at_zero=high0,terminal_status='TERMINAL_HIGH_RANGE' if terminal_root is not None else 'NO_TERMINAL_HIGH_RANGE',
                    terminal_beta_root_deg=terminal_root,terminal_beta_grid_deg=terminal_grid,
                    subsequent_checked_nodes=sum(b>=terminal_root for b in angles) if terminal_root is not None else 0,
                    decreasing_steps=sum(b<a for a,b in zip(freq[:-1],freq[1:])),count_up=len(ups),count_down=sum(r['direction']=='down' for r in crossings))
                results.append(row);cache[key]=row
                print('prediction',case['id'],q,status,[r['beta_deg'] for r in crossings],flush=True)
        atomic(path,results);return results

    def run(self):
        if (self.folder/'manifest.json').exists():self.verify();print('REUSED no new solutions');return
        self.preflight();start=time.perf_counter();curves=[]
        for case in self.cfg['cases']:
            for beta in self.cfg['angles_deg'][case['id']]:curves.append(self.value(case,beta))
            print('curve',case['id'],'PASS',flush=True)
        csvwrite(self.folder/'model_curves.csv',curves)
        predictions=self.predictions(curves);prediction_hash=sha(self.folder/'predictions.json')
        full=rows_csv(self.source/'spectra.csv');refs={(r['id'],float(r['beta_deg'])):r for r in full}
        provenance={(r['id'],r['beta_deg']):r for r in read(self.folder/'full_point_provenance.json')}
        comparison=[]
        for v in curves:
            f=refs[v['id'],v['beta_deg']];c=next(c for c in self.cfg['cases'] if c['id']==v['id']);w=float(f['omega1'])
            row=dict(**v,regime=c['regime'],group=c['group'],epsilon=c['epsilon'],mu=c['mu'],aspect=c['aspect'],
                k_joint=c['k_joint'],chi_s=c['chi_s'],kappa_s=c['kappa_s'],omega_full=w,omega2_full=float(f['omega2']),
                gap12_full=float(f['g12']),zeta_full=(w/c['epsilon'])**2,Lambda_full=float(f['Lambda1']),
                full_point=f['point_file'],full_sha256=provenance[v['id'],v['beta_deg']]['sha256'])
            for kind in ['L1','C','comp','primary']:
                value=v['omega_'+kind];row['error_'+kind+'_percent']=None if value is None else 100*(value/w-1)
                row['absolute_error_'+kind]=None if value is None else abs(value-w)
                row['Lambda_'+kind]=None if value is None else math.sqrt(value/c['epsilon'])
            comparison.append(row)
        csvwrite(self.folder/'comparison.csv',comparison)
        full_thresholds=read(self.source/'thresholds.json');threshold_rows=[];prediction_rows=[]
        for pred in predictions:
            case=next(c for c in self.cfg['cases'] if c['id']==pred['id'])
            source=next(r for r in full_thresholds if r['id']==case['id'] and r['q']==pred['q'] and r['normalization']=='omega90')
            bf=source['terminal_beta_root_deg'];v=self.value(case,bf,'archived full threshold comparison')
            level_log=math.log(pred['omega_reference']/source['omega_target'])
            curve_log=math.log(v['omega_primary']/(pred['q']*source['omega_target']));slope=v['slope_log_per_degree']
            original=next(r for r in source['crossings'] if r['beta_deg']==bf);a,b=original['bracket'];fa,fb=original['signs']
            full_slope=(math.log1p(fb)-math.log1p(fa))/(b-a)
            root=next((r for r in pred['crossings'] if r['beta_deg']==pred['terminal_beta_root_deg']),None)
            model=None if root is None else root['beta_deg']
            threshold_rows.append(dict(id=case['id'],regime=case['regime'],group=case['group'],epsilon=case['epsilon'],mu=case['mu'],aspect=case['aspect'],q=pred['q'],
                primary=pred['primary'],status=pred['status'],terminal_status=pred['terminal_status'],high_at_zero=pred['high_at_zero'],
                count_up=pred['count_up'],count_down=pred['count_down'],beta_full=bf,beta_model=model,
                signed_error_deg=None if model is None else model-bf,absolute_error_deg=None if model is None else abs(model-bf),
                error_percent=None if model is None else 100*(model/bf-1),xi=None if root is None else root['xi'],
                slope_log_per_degree=None if root is None else root['slope_log_per_degree'],full_angle_model_slope=slope,
                archived_full_bracket_slope=full_slope,full_slope_qualification='secant of accepted 2e-6-degree bracket; not exact derivative',
                residual=None if root is None else root['residual'],bracket_width_deg=None if root is None else root['bracket'][1]-root['bracket'][0],
                reference_frequency_error_percent=100*math.expm1(level_log),curve_frequency_error_at_full_percent=100*math.expm1(curve_log),
                linearized_level_shift_deg=level_log/slope,linearized_curve_shift_deg=-curve_log/slope,linearized_total_shift_deg=(level_log-curve_log)/slope,
                subsequent_checked_nodes=pred['subsequent_checked_nodes'],source_threshold_sha256=sha(self.source/'thresholds.json')))
            for root in pred['crossings']:
                prediction_rows.append(dict(id=case['id'],q=pred['q'],primary=pred['primary'],status=pred['status'],direction=root['direction'],
                    beta_deg=root['beta_deg'],xi=root['xi'],residual=root['residual'],left=root['bracket'][0],right=root['bracket'][1],
                    sign_left=root['signs'][0],sign_right=root['signs'][1],slope_log_per_degree=root['slope_log_per_degree'],origin=root['origin'],source_id=root.get('source_id')))
        csvwrite(self.folder/'threshold_comparison.csv',threshold_rows);csvwrite(self.folder/'predictions.csv',prediction_rows)
        for q in [.90,.95]:csvwrite(self.folder/f'thresholds_q{int(q*100)}.csv',[r for r in threshold_rows if r['q']==q])
        require(sha(self.folder/'predictions.json')==prediction_hash,'Predictions mutated')
        atomic(self.folder/'prediction_provenance.json',dict(predictions_sha256=prediction_hash,forecast_before_threshold_table_read=True,
            prior_HSB_knowledge=True,fit_parameters=0,full_thresholds_sha256=sha(self.source/'thresholds.json')))
        regress=[];old=rows_csv(self.d1_source/'comparison.csv')
        mapping={(r['id'],float(r['beta_deg'])):r for r in comparison}
        for r in old:
            v=mapping[r['id'],float(r['beta_deg'])]
            errors={k:abs(v['omega_'+k]/float(r['omega_'+k])-1) for k in ['L1','C','primary']}
            require(max(errors.values())<=self.p['independent_frequency_tol'],'HS-D1 curves regression')
            regress.append(dict(id=r['id'],beta_deg=float(r['beta_deg']),errors=errors))
        old_thresholds=rows_csv(self.d1_source/'threshold_comparison.csv');root_regress=[]
        for r in old_thresholds:
            v=next(v for v in threshold_rows if v['id']==r['id'] and v['q']==float(r['q']))
            err=abs(v['beta_model']-float(r['beta_model']));require(err<=self.cfg['thresholds']['bracket_width_deg'],'HS-D1 threshold regression')
            root_regress.append(dict(id=r['id'],q=float(r['q']),absolute_error_deg=err))
        atomic(self.folder/'HS_D1_regression.json',dict(status='PASS',curves=regress,thresholds=root_regress))
        atomic(self.folder/'timing.json',dict(production_seconds=time.perf_counter()-start))
        analyze_saved(self.folder)
        self.diagnose()
        manifest=dict(status='PASS',identity=self.identity,variants=34,composite_variants=25,control_I_variants=9,
            unique_leading_points=len(self.records),new_points=sum(r['origin']=='new' for r in self.records.values()),
            reused_counts={s:sum(r['origin']==s for r in self.records.values()) for s in sorted({r['origin'] for r in self.records.values()} )},
            full_archived_grid_points=len(curves),full_thresholds=len(threshold_rows),new_full_solves=0,
            unresolved_predictions=[r['id']+str(r['q']) for r in predictions if r['status'] not in ['PASS','ALL_SAMPLED_HIGH','NO_CROSSING_FOUND']])
        manifest['artifacts']=artifact_hashes(self.folder);atomic(self.folder/'manifest.json',manifest)
        print('PASS',self.folder,flush=True)

    def diagnose(self):
        """Only recorded adjacent decreases; retain original accepted spectra."""
        path=self.folder/'local_neighbour_diagnostics.json'
        if path.exists():require(read(path)['status']=='PASS','Local diagnostic status');return
        events=[r for r in read(self.folder/'neighbour_audit.json') if r['decrease']]
        curves=rows_csv(self.folder/'comparison.csv');out=[]
        try:
            for event in events:
                values=[]
                for beta in [event['left'],event['right']]:
                    row=next(r for r in curves if r['id']==event['id'] and float(r['beta_deg'])==beta)
                    case=next(c for c in self.cfg['cases'] if c['id']==event['id'])
                    inside=self.records[row['inner_key']];outside=self.records[row['outer_key']]
                    li=d1.inner.independent_roots(np.array(inside['zeta']),inside['parameter'],case['k_joint'],case['mu'])
                    co=d1.outer.independent(np.array(outside['zeta']),case['k_joint'],outside['kappa_s'],beta,case['mu'],self.p)
                    require(max(abs(r['relative_frequency_error']) for r in li)<=self.p['independent_frequency_tol'],'Local inner comparison')
                    z=co[0]['zeta'] if case['role']=='control_I' else co[0]['zeta']+li[0]['zeta']-float(row['zeta_infty'])
                    values.append(dict(beta_deg=beta,zeta_primary_original=float(row['zeta_primary']),zeta_primary_boundary=z,
                        omega_original=float(row['omega_primary']),omega_boundary=case['epsilon']*math.sqrt(z),
                        analytic_slope_log_per_degree=float(row['slope_log_per_degree']),inner_roots=li,outer_roots=co))
                left,right=values;change=right['omega_boundary']/left['omega_boundary']-1
                require(change>0 and all(v['analytic_slope_log_per_degree']>0 for v in values),'Unresolved adjacent decrease')
                out.append(dict(event=event,values=values,boundary_relative_change=change,
                    width_deg=event['right']-event['left'],status='PASS',
                    qualification='tiny raw Ritz decrement is not confirmed by independent boundary roots; positive local sensitivity agrees; no point replaced or global monotonicity claimed'))
            atomic(path,dict(status='PASS',trigger='each negative adjacent difference, regardless of jump magnitude',events=out,
                additional_independent_boundary_roots=12*len(out),new_full_solves=0,new_ritz_solves=0,original_spectra_preserved=True))
            print('local diagnostic PASS',len(out),'events',flush=True)
        except Exception as exc:
            atomic(path,dict(status='UNRESOLVED',events=out,error=repr(exc)));raise


def artifact_hashes(folder):
    excluded={'manifest.json','README.md','completion_checks.json','documentation_manifest.json','preservation_checks.json','figures_manifest.json'}
    return {p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name not in excluded and 'figures' not in p.parts}


def error_statistics(values):
    a=np.asarray(values,dtype=float);v=abs(a)
    return dict(max_abs=float(max(v)),mean_abs=float(np.mean(v)),median_abs=float(np.median(v)),p95_abs=float(np.percentile(v,95)),signed_bias=float(np.mean(a)))


def analyze_saved(folder):
    """Only saved values; no eigensolver, matrix or full spectral call."""
    cfg=read(folder/'contract.json')['config'];rows=rows_csv(folder/'comparison.csv');thresholds=rows_csv(folder/'threshold_comparison.csv')
    source=ROOT/cfg['full_source'];source_csv=rows_csv(source/'thresholds.csv');source_predictions=rows_csv(source/'predictions.csv')
    source_json=read(source/'thresholds.json');checks=[]
    for r in thresholds:
        name,q=r['id'],float(r['q'])
        a=next(v for v in source_json if v['id']==name and v['q']==q and v['normalization']=='omega90')
        b=next(v for v in source_csv if v['id']==name and float(v['q'])==q and v['normalization']=='omega90')
        c=next(v for v in source_predictions if v['id']==name and float(v['q'])==q and v['normalization']=='omega90')
        require(a['terminal_beta_root_deg']==float(b['terminal_beta_root_deg'])==float(c['full_beta_deg'])==float(r['beta_full']),'Full threshold table consistency')
        checks.append(dict(id=name,q=q,status='PASS',beta_full=float(r['beta_full'])))
    require(len(checks)==68,'All full thresholds compared')
    atomic(folder/'full_threshold_consistency.json',dict(status='PASS',checks=checks,
        hashes={n:sha(source/n) for n in ['thresholds.json','thresholds.csv','predictions.csv']}))
    intervals=[('all',0,90),('0_5',0,5),('5_30',5,30),('30_90',30,90)];stats=[];neighbours=[];maxima_rows=[]
    for case in cfg['cases']:
        rr=sorted([r for r in rows if r['id']==case['id']],key=lambda r:float(r['beta_deg']))
        full_max=max(rr,key=lambda r:float(r['omega_full']));model_max=max(rr,key=lambda r:float(r['omega_primary']))
        maxima_rows.append(dict(id=case['id'],full_grid_max_beta=float(full_max['beta_deg']),model_grid_max_beta=float(model_max['beta_deg']),
            full_grid_max=float(full_max['omega_full']),full_omega90=float(rr[-1]['omega_full']),model_grid_max=float(model_max['omega_primary']),
            model_omega90=float(rr[-1]['omega_primary']),full_grid_max_equals_omega90=full_max['omega_full']==rr[-1]['omega_full'],
            primary_min_zeta=min(float(r['zeta_primary']) for r in rr)))
        for label,lo,hi in intervals:
            selected=[r for r in rr if lo<=float(r['beta_deg'])<=hi]
            for model in ['L1','C','comp','primary']:
                valid=[r for r in selected if r['error_'+model+'_percent']!='']
                if not valid:continue
                result=dict(id=case['id'],regime=case['regime'],group=case['group'],epsilon=case['epsilon'],mu=case['mu'],aspect=case['aspect'],
                    interval=label,model=model,n=len(valid),**error_statistics([r['error_'+model+'_percent'] for r in valid]))
                result['max_abs_error_omega']=max(float(r['absolute_error_'+model]) for r in valid)
                worst=max(valid,key=lambda r:abs(float(r['error_'+model+'_percent'])));result['worst_beta_deg']=float(worst['beta_deg']);stats.append(result)
        for a,b in zip(rr[:-1],rr[1:]):
            wa,wb=float(a['omega_primary']),float(b['omega_primary']);jump=abs(wb-wa)/max(wa,wb)
            neighbours.append(dict(id=case['id'],left=float(a['beta_deg']),right=float(b['beta_deg']),jump=jump,decrease=wb<wa,
                triggered=jump>cfg['frequency_map']['neighbour_jump_relative'],resolution='both component prefixes counted and residual-checked'))
    csvwrite(folder/'frequency_statistics.csv',stats);csvwrite(folder/'grid_maxima.csv',maxima_rows);atomic(folder/'neighbour_audit.json',neighbours)
    record_files=list((folder/'points').glob('*.json'));records=[read(p) for p in record_files];maxima={};margins=[];independent=[];tau=[];gaps=[]
    point_map={r['point_key']:r for r in records};uses=[]
    for r in rows:
        case=next(c for c in cfg['cases'] if c['id']==r['id'])
        for kind,key in [('L1',r['inner_key']),('C',r['outer_key'])]:
            point=point_map[key];parameter=float(r['tau']) if kind=='L1' else float(r['beta_deg'])
            ks=None if kind=='L1' else case['kappa_s'] if case['role']=='control_I' else RIGID
            require(point_identity(kind,case['k_joint'],case['mu'],parameter,ks)==point_identity(point['kind'],point['k_joint'],point['mu'],point['parameter'],point['kappa_s']),'Point-use identity')
            uses.append(dict(id=case['id'],epsilon=case['epsilon'],mu=case['mu'],aspect=case['aspect'],chi_s=case['chi_s'],
                k_joint=case['k_joint'],beta_deg=float(r['beta_deg']),kind=kind,parameter=parameter,kappa_s=ks,key=key,
                origin=point['origin'],source=point['source'],source_sha256=point['source_sha256'],
                point_sha256=sha(folder/'points'/(key+'.json')),shapes_sha256=point['npz_sha256']))
    csvwrite(folder/'point_use_provenance.csv',uses)
    for r in records:
        if r['kind']=='L1' and r['parameter']!=RIGID:tau.append(r['parameter'])
        for order in r['orders']:
            d=order['diagnostics']
            for key in ['matrix_residual','mass_gram','energy','essential','symmetry']:
                maxima[key]=max(maxima.get(key,0.),d[key])
            for target,source in [('natural','natural_boundary' if r['kind']=='L1' else 'natural'),('strong','strong_ode' if r['kind']=='L1' else 'strong')]:
                maxima[target]=max(maxima.get(target,0.),d[source])
            z=order['zeta'];gaps.append(math.sqrt(z[1]/z[0])-1)
        if r['kind']=='L1':
            checks=r['independent']['counts'];independent += [abs(v['relative_frequency_error']) for v in r['independent']['roots']]
        else:
            checks=r['independent'];independent += [v['frequency_error'] for v in checks]
        for v in checks:margins += [v['lower']['inertia_margin'],v['upper']['inertia_margin']]
    diagnostic_path=folder/'local_neighbour_diagnostics.json';extra_roots=0
    if diagnostic_path.exists():
        diagnostic=read(diagnostic_path);require(diagnostic['status']=='PASS','Addressed diagnostic status')
        extra_roots=diagnostic['additional_independent_boundary_roots']
        for event in diagnostic['events']:
            for value in event['values']:
                independent += [abs(v['relative_frequency_error']) for v in value['inner_roots']]
                independent += [v['frequency_error'] for v in value['outer_roots']]
    coeff=read(folder/'coefficients.json');pred=read(folder/'predictions.json');roots=[v for r in pred for v in r['crossings']]
    atomic(folder/'numerical_checks.json',dict(status='PASS',max_residuals=maxima,unique_points=len(records),
        max_p_frequency_error=max(r['frequency_convergence'] for r in records),max_p_shape_error=max(r['shape_convergence'] for r in records),
        counted_positions=3*len(records),independent_refined_roots=len(independent),additional_diagnostic_roots=extra_roots,max_independent_frequency_error=max(independent),
        minimum_count_margin=min(margins),maximum_finite_tau=max(tau),minimum_frequency_gap12=min(gaps),
        max_C2_Richardson_error=max(abs(c['C2_Richardson_relative_error']) for c in coeff),
        max_threshold_residual=max(abs(v['residual']) for v in roots),max_slope_FD_relative_error=max(abs(v['differences'][-1]['relative_error']) for v in roots),
        min_global_positive_bound=min(c['orders'][-1]['global_composite_lower_bound'] for c in coeff),new_full_solves=0,
        high_precision_repair=False,finite_tau_never_replaced_by_saturation=True))
    worst={}
    for regime in ['A1','A2','A3','A4','ALL']:
        valid=[r for r in stats if r['model']=='primary' and r['interval']=='all' and (regime=='ALL' or r['regime']==regime)]
        rr=[r for r in thresholds if (regime=='ALL' or r['regime']==regime) and r['error_percent']!='']
        worst[regime]=dict(frequency=max(valid,key=lambda r:r['max_abs']),
            q90=max((r for r in rr if float(r['q'])==.9),key=lambda r:abs(float(r['error_percent']))),
            q95=max((r for r in rr if float(r['q'])==.95),key=lambda r:abs(float(r['error_percent']))))
    atomic(folder/'worst_cases.json',worst)
    aspect_rows=[];aspect_points=[];old_aspect=read(folder/'source_aspect_diagnostics.json')
    for regime in ['A1','A2','A3','A4']:
        a={float(r['beta_deg']):r for r in rows if r['id']==f'B0_{regime}'}
        b={float(r['beta_deg']):r for r in rows if r['id']==f'B3_{regime}_a2'}
        common=sorted(set(a)&set(b));diff=[]
        for beta in common:
            x,y=a[beta],b[beta];relative=float(y['omega_full'])/float(x['omega_full'])-1
            require(x['inner_key']==y['inner_key'] and x['outer_key']==y['outer_key'],'Aspect leading key identity')
            v=dict(regime=regime,beta_deg=beta,full_relative_difference=relative,full_absolute_difference=float(y['omega_full'])-float(x['omega_full']),
                primary_relative_difference=float(y['omega_primary'])/float(x['omega_primary'])-1)
            diff.append(v);aspect_points.append(v)
        original=next(r for r in old_aspect if r['regime']==regime);largest=max(diff,key=lambda r:abs(r['full_relative_difference']))
        aspect_rows.append(dict(regime=regime,n_common=len(common),**{k:v for k,v in largest.items() if k!='regime'},
            independent_ritz_difference=original['ritz_relative_differences'][-1],independent_order_spread=original['ritz_order_spread'],
            qualification='UNRESOLVED sign and magnitude' if regime=='A2' else 'tiny positive difference; larger relative uncertainty for A3; inherited HS-B independent evidence'))
    csvwrite(folder/'aspect_comparison.csv',aspect_rows);csvwrite(folder/'aspect_common_points.csv',aspect_points)
    atomic(folder/'analysis_conventions.json',dict(statistics='unweighted statistics of actual archived nodes; dense root neighborhoods affect weighting',
        intervals='closed intervals; boundary nodes may occur in two groups; not validity boundaries',
        no_root_cases=[dict(id=r['id'],q=r['q'],status=r['status']) for r in pred if not r['crossings']],
        multiple_root_cases=[dict(id=r['id'],q=r['q'],count=len(r['crossings'])) for r in pred if len(r['crossings'])>1],
        downward_crossings=sum(r['count_down'] for r in pred),near_zero_high=sum(r['high_at_zero'] for r in pred),
        discrete_decreases=sum(r['decrease'] for r in neighbours),neighbour_triggers=sum(r['triggered'] for r in neighbours)))
    print('saved-data analysis PASS',folder,flush=True)


def dispatch(action):
    if action=='plot_only':render(output_folder());return
    if action=='analyze':analyze_saved(output_folder());return
    study=GeometryStudy()
    if action=='preflight':study.preflight()
    elif action=='run':study.run()
    elif action=='verify':study.verify()
    elif action=='diagnose':study.diagnose()
    elif action=='finalize':
        m=read(study.folder/'manifest.json')
        # The saved-data aggregate can acquire the addressed diagnostic's roots.
        # Every original spectral point, prediction and comparison remains sealed.
        for name,digest in m['artifacts'].items():
            if name=='numerical_checks.json' and sha(study.folder/name)!=digest:
                require(read(study.folder/'local_neighbour_diagnostics.json')['status']=='PASS','Unaccepted aggregate amendment')
            else:require(sha(study.folder/name)==digest,'Artifact '+name)
        m['artifacts']=artifact_hashes(study.folder);atomic(study.folder/'manifest.json',m)
        print('sealed',study.folder)


def render(folder):
    """Seven scientific figures from CSV only, with no Study construction."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder);start=time.perf_counter();cfg=read(folder/'contract.json')['config']
    cases=cfg['cases'];rows=rows_csv(folder/'comparison.csv');stats=rows_csv(folder/'frequency_statistics.csv')
    thresholds=rows_csv(folder/'threshold_comparison.csv');aspect=rows_csv(folder/'aspect_common_points.csv')
    figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.grid':True,'grid.alpha':.23,
        'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00']
    by_id={c['id']:sorted([r for r in rows if r['id']==c['id']],key=lambda r:float(r['beta_deg'])) for c in cases}

    def save(fig,name):
        fig.tight_layout();fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(figdir/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)

    def curves(ax,selected,variable='omega',label_by='id'):
        for i,c in enumerate(selected):
            rr=by_id[c['id']];x=[float(r['beta_deg']) for r in rr]
            full='omega_full' if variable=='omega' else 'Lambda_full';model='omega_primary' if variable=='omega' else 'Lambda_primary'
            label=c['id'] if label_by=='id' else ('ε='+str(c['epsilon']) if label_by=='epsilon' else 'μ='+str(c['mu']))
            color=palette[i%len(palette)]
            ax.plot(x,[float(r[model]) for r in rr],color=color,lw=1.4,label=label)
            ax.plot(x,[float(r[full]) for r in rr],'o',mfc='none',mec=color,ms=2.6,markevery=max(1,len(rr)//20))
        ax.set(xlim=(0,90),xlabel='Угол β, град',ylabel='Первая частота ω₁' if variable=='omega' else 'Параметр частоты Λ₁')
        ax.legend(fontsize=7)

    fig,axes=plt.subplots(3,3,figsize=(15,11))
    for i,regime in enumerate(['A1','A3','A4']):
        groups=[['B0','B1'],['B0','B2'],['B0','B3','B4']]
        for j,groupset in enumerate(groups):
            selected=[c for c in cases if c['regime']==regime and c['group'] in groupset]
            curves(axes[i,j],selected);axes[i,j].set_title(regime+' — '+['толщина','распределение длин','сечение и сочетание изменений'][j])
    fig.suptitle('Составное приближение: линии; полные частоты HS-B: кружки. Только сохранённые узлы.',fontsize=11)
    save(fig,'01_composite_frequencies')

    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for ax,groupset,title in zip(axes,[['B0','B1'],['B0','B2'],['B0','B3','B4']],['Толщина','Распределение длин','Сечение и сочетание изменений']):
        curves(ax,[c for c in cases if c['regime']=='A2' and c['group'] in groupset]);ax.set_title(title)
    fig.suptitle('A2: самостоятельная HS-C I (линии), архив HS-B (кружки); составление не применяется',fontsize=11)
    save(fig,'02_A2_frequencies')

    ids=[c['id'] for c in cases];intervals=['all','0_5','5_30','30_90']
    matrix=np.array([[float(next(r for r in stats if r['id']==name and r['interval']==interval and r['model']=='primary')['max_abs']) for interval in intervals] for name in ids])
    fig,ax=plt.subplots(figsize=(7.8,12));im=ax.imshow(matrix,aspect='auto',cmap='YlOrRd')
    ax.set_yticks(range(34),ids,fontsize=8);ax.set_xticks(range(4),['Все узлы','0–5°','5–30°','30–90°'])
    ax.grid(False)
    for i in range(34):
        for j in range(4):ax.text(j,i,f'{matrix[i,j]:.4f}',ha='center',va='center',fontsize=7,color='white' if matrix[i,j]>.65*matrix.max() else 'black')
    fig.colorbar(im,ax=ax,label='Максимум модуля ошибки частоты, %',shrink=.6)
    ax.set_title('Составное приближение A1/A3/A4; HS-C I для A2\nИнтервалы группировки не задают область применимости',fontsize=10)
    save(fig,'03_frequency_errors')

    fig,axes=plt.subplots(1,2,figsize=(11,12),sharey=True)
    for ax,q in zip(axes,[.9,.95]):
        rr=[next(r for r in thresholds if r['id']==name and float(r['q'])==q) for name in ids]
        y=[float(r['error_percent']) if r['error_percent'] else np.nan for r in rr]
        ax.barh(range(34),y,color=['#0072B2' if v>=0 else '#D55E00' for v in y]);ax.axvline(0,color='black',lw=.7)
        ax.set_yticks(range(34),ids,fontsize=8);ax.set_ylim(33.5,-.5);ax.set(xlabel='Подписанная ошибка угла, %',title=f'Уровень {int(q*100)}%')
    fig.suptitle('Опорные частоты: собственные значения каждой модели при 90°',fontsize=11)
    save(fig,'04_threshold_errors')

    fig,axes=plt.subplots(2,4,figsize=(16,7))
    for j,regime in enumerate(['A1','A2','A3','A4']):
        selected=sorted([c for c in cases if c['regime']==regime and c['group'] in ['B0','B1']],key=lambda c:c['epsilon'])
        curves(axes[0,j],selected,'Lambda','epsilon');axes[0,j].set_title(regime)
        for i,c in enumerate(selected):
            rr=by_id[c['id']]
            axes[1,j].plot([float(r['beta_deg']) for r in rr],[float(r['error_primary_percent']) for r in rr],color=palette[i],label='ε='+str(c['epsilon']))
        axes[1,j].axhline(0,color='black',lw=.6);axes[1,j].set(xlabel='Угол β, град',ylabel='Ошибка частоты, %',xlim=(0,90))
    fig.suptitle('Толщина при μ=0,01, aspect=4: сравнение Λ₁ учитывает нормировку ω=εΛ₁²',fontsize=11)
    save(fig,'05_epsilon')

    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for ax,regime in zip(axes.flat,['A1','A2','A3','A4']):
        selected=sorted([c for c in cases if c['regime']==regime and c['group'] in ['B0','B2']],key=lambda c:c['mu'])
        curves(ax,selected,'Lambda','mu');ax.set_title(regime)
    fig.suptitle('Распределение длин: ε=0,005, aspect=4. Линии — приближение, кружки — HS-B.',fontsize=11)
    save(fig,'06_mu')

    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for ax,regime in zip(axes.flat,['A1','A2','A3','A4']):
        rr=[r for r in aspect if r['regime']==regime];x=[float(r['beta_deg']) for r in rr]
        ax.plot(x,[float(r['full_relative_difference']) for r in rr],'o',ms=3,color='#0072B2',label='Полная: ω(aspect=2)/ω(aspect=4)−1')
        ax.axhline(0,color='#009E73',lw=1,label='Совпадающий ведущий прогноз')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.set(xlim=(0,90),xlabel='Угол β, град',ylabel='Относительная разность',
            title=regime+(': знак и величина не разрешены' if regime=='A2' else ''))
        ax.legend(fontsize=7)
    fig.suptitle('Эффект aspect: только общие архивные узлы; новые полные решения отсутствуют',fontsize=11)
    save(fig,'07_aspect')
    for name in ['comparison.csv','frequency_statistics.csv','threshold_comparison.csv','aspect_common_points.csv']:
        (figdir/name).write_bytes((folder/name).read_bytes())
    atomic(folder/'figures_manifest.json',dict(calculation_mode='plot_only',new_spectral_solves=0,dpi=320,
        rendering_seconds=time.perf_counter()-start,artifacts={p.name:sha(p) for p in figdir.iterdir() if p.is_file()}))
    print('plot_only PASS',figdir,flush=True)
