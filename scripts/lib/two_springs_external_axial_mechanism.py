"""Bounded transition preset; the existing physical/count solver is unchanged."""
import hashlib
import json
import math
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_axial_transition as old

previous=old.previous
ROOT=old.ROOT
BASE=old.BASE.with_name('external_axial_transition_mechanism')
CONFIG=ROOT/'data/input/isotropic_two_springs_external_axial_transition_mechanism.json'
physics=old.physics
api=old.api
require=old.require
read_csv=old.read_csv
FIGURES=('01_frequency','02_signed_growth','03_fixed_angle_transition','04_maxima')


def row_from_point(pair,k,b,rec,path,origin):
    w=[r['omega'] for r in rec['roots']]
    return dict(pair=pair,chi_u=k,beta_deg=b,**{f'{n}{i}':r[n] for i,r in enumerate(rec['roots'],1) for n in ('omega','Lambda')},gap12=w[1]/w[0]-1,origin=origin,point_file=path.relative_to(ROOT).as_posix(),point_sha256=api().sha(path),shape_sha256=rec['shape_sha256'])


def quality_rows(row,rec):
    ans=[]
    for i,r in enumerate(rec['roots'],1):
        d=r['diagnostics']
        ans.append(dict(pair=row['pair'],chi_u=row['chi_u'],beta_deg=row['beta_deg'],position=i,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],status='PASS'))
    return ans


class Experiment(old.Experiment):
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={};self.sources={};self.manifests={}
        self.started=time.perf_counter();self.ritz_new=0;self.leading_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,s in self.cfg['sources'].items():
            folder=ROOT/s['folder'];self.check(folder/'manifest.json',s['manifest_sha256']);previous.ranges.saved.diagnostic.verify_source(folder)
            self.sources[name]=folder;self.manifests[name]=json.loads((folder/'manifest.json').read_text())
            c=json.loads(self.source(name,'contract.json').read_text())['config']
            require(c['material']==self.cfg['material'] and c['mu']==self.cfg['mu'] and c['variant']==self.cfg['variant'] and c['cases']==self.cfg['cases'],'Source physics')
            require(c['map']['spectrum_semantics']=='sorted_positions' and c['map']['K_guard']==3,'Source semantics')
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        self.run.catalogs={float(k):v for k,v in json.loads(self.source('transition','pole_catalogs.json').read_text()).items()}
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        self.archive={};self.overlap=[]
        for file in ('spectrum.csv','addressed_spectrum.csv'):
            for r in read_csv(self.source('transition',file)):
                k=float(r['chi_u'])
                if (file=='spectrum.csv' and k not in self.cfg['archived_chi_values']) or (file=='addressed_spectrum.csv' and k<1e-4):continue
                self.ingest(ROOT/r['point_file'],r['point_sha256'],'transition',r)
        archived=[row_from_point(p,k,b,*v) for (p,k,b),v in sorted(self.archive.items())]
        self.csv('reference_spectrum',archived)
        inventory=[]
        for p in self.cfg['cases']:
            for k in sorted({kk for pp,kk,b in self.archive if pp==p}|set(self.cfg['chi_values'])):
                have=sorted(b for pp,kk,b in self.archive if (pp,kk)==(p,k))
                inventory.append(dict(pair=p,chi_u=k,available_beta_deg=have,missing_requested_map_beta_deg=[b for b in self.cfg['map']['parameter_grid'] if b not in have] if k in self.cfg['chi_values'] else None))
        self.a.write(self.folder/'archive_inventory.json',inventory)
        regression=[]
        for r in read_csv(self.source('transition','transition_refined_brackets.csv')):
            p=r['pair'];b=float(r['beta_deg']);vals=[]
            for side in ('lower','upper'):
                k=float(r['chi_'+side]);d=self.archive[p,k,b][0]['roots'][0]['omega']/self.archive[p,k,0.][0]['roots'][0]['omega']-1
                require(abs(d-float(r['D0_'+side]))<1e-14,'Old sign bracket changed');vals.append(d)
            require(vals[0]<0<vals[1],'Old sign bracket')
            regression.append(dict(pair=p,beta_deg=b,chi_lower=float(r['chi_lower']),chi_upper=float(r['chi_upper']),D_lower=vals[0],D_upper=vals[1],status='PASS'))
        self.csv('archive_bracket_regression',regression)
        self.plan=json.loads((self.folder/'refinement_plan.json').read_text()) if (self.folder/'refinement_plan.json').exists() else []
        self.allowed={(p,k):set(self.cfg['preflight_angles']) for p in self.cfg['cases'] for k in self.cfg['chi_values']}
        for r in self.plan:self.allowed.setdefault((r['pair'],r['chi_u']),set()).update(r['new_beta_deg'])
        self.extra=read_csv(self.folder/'addressed_spectrum.csv') if (self.folder/'addressed_spectrum.csv').exists() else []
        self.extra_by={(r['pair'],float(r['chi_u']),float(r['beta_deg'])):r for r in self.extra}
        self.extra_quality=read_csv(self.folder/'addressed_quality.csv') if (self.folder/'addressed_quality.csv').exists() else []
        self.requests=json.loads((self.folder/'addressed_requests.json').read_text()) if (self.folder/'addressed_requests.json').exists() else []

    def point(self,pair,k,beta,previous_w=None):
        if (pair,k,beta) in self.archive:return self.archive[pair,k,beta]
        require(k in self.cfg['new_chi_values'],'New full map outside three-stiffness scope')
        return previous.Experiment.point(self,pair,k,beta,previous_w)

    def run_study(self):
        self.preflight();rows=self.sweep()
        if not (self.folder/'initial_summary.csv').exists():self.analyze(rows,'initial_')
        for n in range(1,self.cfg['refinement']['max_rounds']+1):
            if not self.refine(rows,n):break
            rows=self.sweep()
        else:require(not self.refine(rows,n+1),'Unresolved refinement')
        summary=self.analyze(rows);old.details(self.folder,self.cfg,rows,summary)
        self.a.write(self.folder/'calculation_complete.json',dict(status='PASS',constructions=len(summary),points=len(rows),new_points=len(rows),reused_reference_points=len(self.archive),refinement_rounds=max((r['round'] for r in self.plan),default=0)))
        self.thresholds(rows);self.curvature(rows);self.register_diagnostic_plan(summary)
        extra_maps=self.diagnostic_maps() if (self.folder/'diagnostic_stiffness_plan.json').exists() else []
        self.independent_curvature_checks(rows+extra_maps)
        if extra_maps:self.additional_threshold(rows+extra_maps)
        self.extremum_evidence(rows+extra_maps)
        self.consistency(rows+extra_maps)

    def csv(self,name,rows):
        return super().csv(getattr(self,'output_prefix','')+name,rows)

    def register_diagnostic_plan(self,summary):
        target=self.folder/'diagnostic_stiffness_plan.json'
        if target.exists():return
        curv=read_csv(self.folder/'curvature_summary.csv')
        trigger=all(next(r for r in summary if r['pair']==p and r['chi_u']==.00027)['beta_max']==0 and 0<next(r for r in summary if r['pair']==p and r['chi_u']==.0003)['beta_max']<90 and next(r for r in curv if r['pair']==p and float(r['chi_u'])==.00027)['status']=='NEGATIVE' and next(r for r in curv if r['pair']==p and float(r['chi_u'])==.0003)['status']=='POSITIVE' for p in self.cfg['cases'])
        if trigger:
            self.a.write(target,dict(defined_before_extra_maps=True,chi_values=[.000276,.000283],pairs=list(self.cfg['cases']),reason='Three primary curves leave the movement from beta=0 to interior maxima unresolved. C_full changes sign between chi=.00027 and .00030; two bounded intermediate probes resolve onset without assuming continuous movement or competing maxima.',selection='C estimates only locate diagnostic stiffnesses; all frequencies and maxima are independently solved.',source_summary_sha256=self.a.sha(self.folder/'summary.csv'),source_curvature_sha256=self.a.sha(self.folder/'curvature_summary.csv'),grid='same bounded 0-20 step1,22-90 step2 and small-angle nodes; local extrema/threshold refinement only',maximum_additional_stiffnesses=2,no_further_expansion=True))

    def diagnostic_maps(self):
        plan=json.loads((self.folder/'diagnostic_stiffness_plan.json').read_text())
        require(plan['defined_before_extra_maps'] and len(plan['chi_values'])<=2 and all(.00023713737056616554<k<.00031622776601683794 for k in plan['chi_values']),'Bounded additional stiffness plan')
        cfg=self.cfg
        self.cfg={**cfg,'chi_values':plan['chi_values'],'new_chi_values':plan['chi_values']};self.output_prefix='diagnostic_'
        try:
            for p in cfg['cases']:
                for k in plan['chi_values']:self.allowed.setdefault((p,k),set()).update(cfg['preflight_angles'])
            target=self.folder/'diagnostic_preflight.json'
            if not target.exists():
                checked=[]
                for p in cfg['cases']:
                    for k in plan['chi_values']:
                        for b in cfg['preflight_angles']:
                            rec,path,origin=self.point(p,k,b)
                            if b==0:require(abs(rec['roots'][0]['omega']/self.archive[p,0.,0.][0]['roots'][0]['omega']-1)<cfg['policy']['reference_frequency_tol'],'Diagnostic zero regression')
                            checked.append(row_from_point(p,k,b,rec,path,origin))
                self.csv('preflight_points',checked);self.a.write(target,dict(status='PASS',points=len(checked),basis='Registered interior-maximum diagnostic; same three-root gates'))
            require(json.loads(target.read_text())['status']=='PASS','Diagnostic preflight failure')
            rows=self.sweep()
            if not (self.folder/'diagnostic_initial_summary.csv').exists():previous.analyze(self.folder,self.cfg,rows,'diagnostic_initial_')
            for n in range(1,cfg['refinement']['max_rounds']+1):
                if not self.refine(rows,n):break
                rows=self.sweep()
            else:require(not self.refine(rows,n+1),'Diagnostic refinement unresolved')
            summary=previous.analyze(self.folder,self.cfg,rows,'diagnostic_')
            # Avoid old.details, whose fixed output names belong to primary maps.
            extrema=[]
            for s in summary:
                rr=[r for r in rows if (r['pair'],r['chi_u'])==(s['pair'],s['chi_u'])]
                for i in range(1,len(rr)-1):
                    lo,v,hi=rr[i-1:i+2];a=v['omega1']-lo['omega1'];b=hi['omega1']-v['omega1']
                    if a*b<0 and min(abs(a),abs(b))/v['omega1']>cfg['refinement']['extremum_frequency_rtol']:
                        extrema.append(dict(pair=s['pair'],chi_u=s['chi_u'],kind='maximum' if a>0 else 'minimum',beta_deg=v['beta_deg'],omega1=v['omega1'],left_deg=lo['beta_deg'],right_deg=hi['beta_deg']))
            self.csv('interior_extrema',extrema);self.curvature(rows)
            self.a.write(self.folder/'diagnostic_complete.json',dict(status='PASS',constructions=len(summary),points=len(rows),new_points=len(rows)))
            return rows
        finally:self.cfg=cfg;self.output_prefix=''

    def independent_curvature_checks(self,rows):
        target=self.folder/'curvature_ritz.csv'
        if target.exists():return
        checks=[]
        # Close-to-zero-growth cases selected after the primary curvature fit.
        cases=[('J0E4',.000283),('J4E0',.000276)] if (self.folder/'diagnostic_stiffness_plan.json').exists() else [(p,.00027) for p in self.cfg['cases']]
        for p,k in cases:
            for b in (0.,.25,1.):
                row=next(r for r in rows if (r['pair'],r['chi_u'],r['beta_deg'])==(p,k,b));expected=np.array([float(row[f'omega{i}']) for i in (1,2,3)])
                for order in self.cfg['policy']['ritz_orders']:
                    w,_,d=physics.ritz(self.run.model,self.case(p,k,b),order);self.ritz_new+=1
                    err=float(max(abs(w/expected-1)))
                    require(err<self.cfg['policy']['ritz_frequency_tol'] and d['positive_K_M'] and d['energy_error']<self.cfg['policy']['energy_tol'] and d['mass_gram_error']<self.cfg['policy']['mass_gram_tol'] and d['constraint_residual']<self.cfg['policy']['boundary_tol'],'Near-zero-growth independent Ritz check')
                    checks.append(dict(pair=p,chi_u=k,beta_deg=b,order=order,omega1=w[0],relative_frequency_error=err,**d))
        self.csv('curvature_ritz',checks)

    def extra_point(self,p,k,b,reason,main):
        key=(p,float(k),float(b))
        if key in self.archive:
            rec,path,origin=self.archive[key];return row_from_point(p,k,b,rec,path,origin)
        if key in main:return main[key]
        if key in self.extra_by:return self.extra_by[key]
        require(1e-4<k<1e-3 and b in [0.,.25,5.,30.,60.,90.],'Addressed solve outside bounded scope')
        self.requests.append(dict(pair=p,chi_u=k,beta_deg=b,reason=reason));self.a.write(self.folder/'addressed_requests.json',self.requests)
        c=self.case(p,k,b);rec=self.run.solve(c);self.validate(rec,c)
        path=previous.previous.lengths.point_path(self.folder,c);r=row_from_point(p,k,b,rec,path,'addressed_new')
        self.extra.append(r);self.extra_by[key]=r;self.extra_quality+=quality_rows(r,rec)
        self.csv('addressed_spectrum',self.extra);self.csv('addressed_quality',self.extra_quality)
        return r

    def thresholds(self,rows):
        main={(r['pair'],r['chi_u'],r['beta_deg']):r for r in rows};history=[];result=[]
        src=read_csv(self.folder/'reference_spectrum.csv')
        plan=dict(defined_before_threshold_solves=True,angles=self.cfg['threshold_search']['angles'],relative_width=self.cfg['threshold_search']['relative_width'],max_steps=self.cfg['threshold_search']['max_steps'],source_spectrum_sha256=self.a.sha(self.folder/'spectrum.csv'))
        self.a.write(self.folder/'threshold_plan.json',plan)
        for p in self.cfg['cases']:
            zero={float(r['chi_u']):float(r['omega1']) for r in src+rows if r['pair']==p and float(r['beta_deg'])==0}
            for b in self.cfg['threshold_search']['angles']:
                vals=sorted((float(r['chi_u']),float(r['omega1'])/zero[float(r['chi_u'])]-1) for r in src+rows if r['pair']==p and float(r['beta_deg'])==b)
                candidates=[(lo,hi,dl,dh) for (lo,dl),(hi,dh) in zip(vals,vals[1:]) if dl<0<dh]
                require(len(candidates)==1,'Unique resolved sign bracket');lo,hi,dl,dh=candidates[0]
                for n in range(self.cfg['threshold_search']['max_steps']):
                    if hi/lo-1<=self.cfg['threshold_search']['relative_width']:break
                    mid=math.sqrt(lo*hi);z=self.extra_point(p,mid,0.,'own zero for fixed-angle sign search',main);v=self.extra_point(p,mid,b,'fixed-angle sign bisection',main)
                    d=float(v['omega1'])/float(z['omega1'])-1
                    require(abs(d)>self.cfg['threshold_search']['sign_resolution'],'Unresolved sign at frequency uncertainty')
                    history.append(dict(pair=p,beta_deg=b,step=n+1,chi_lower=lo,chi_upper=hi,chi_trial=mid,D_trial=d))
                    if d>0:hi,dh=mid,d
                    else:lo,dl=mid,d
                require(hi/lo-1<=self.cfg['threshold_search']['relative_width'],'Threshold budget exhausted')
                result.append(dict(pair=p,beta_deg=b,chi_lower=lo,chi_upper=hi,relative_width=hi/lo-1,D_lower=dl,D_upper=dh,status='RESOLVED_BRACKET'))
        self.csv('threshold_history',history);self.csv('threshold_brackets',result)

    def additional_threshold(self,rows):
        """One extra fixed angle, triggered by positive small-angle growth."""
        target=self.folder/'additional_threshold_brackets.csv'
        if target.exists():return
        main={(r['pair'],r['chi_u'],r['beta_deg']):r for r in rows}
        plan=dict(defined_before_solves=True,beta_deg=5.,reason='Diagnostic stiffnesses give positive growth at 5 degrees while D(30) or D(90) is negative; resolve one earlier-angle threshold separately.',max_additional_angles=1,relative_width=self.cfg['threshold_search']['relative_width'])
        self.a.write(self.folder/'additional_threshold_plan.json',plan)
        result=[];history=[]
        for p in self.cfg['cases']:
            z={r['chi_u']:r['omega1'] for r in rows if r['pair']==p and r['beta_deg']==0.}
            vals=sorted((r['chi_u'],r['omega1']/z[r['chi_u']]-1) for r in rows if r['pair']==p and r['beta_deg']==5.)
            cand=[(lo,hi,dl,dh) for (lo,dl),(hi,dh) in zip(vals,vals[1:]) if dl<0<dh];require(len(cand)==1,'Additional-angle sign trigger')
            lo,hi,dl,dh=cand[0]
            for step in range(self.cfg['threshold_search']['max_steps']):
                if hi/lo-1<=self.cfg['threshold_search']['relative_width']:break
                mid=math.sqrt(lo*hi);v=self.extra_point(p,mid,5.,'one additional fixed-angle threshold at 5 degrees',main);zero=self.extra_point(p,mid,0.,'own zero for additional-angle threshold',main);d=float(v['omega1'])/float(zero['omega1'])-1
                require(abs(d)>self.cfg['threshold_search']['sign_resolution'],'Additional threshold sign unresolved')
                history.append(dict(pair=p,beta_deg=5.,step=step+1,chi_trial=mid,D_trial=d))
                if d>0:hi,dh=mid,d
                else:lo,dl=mid,d
            result.append(dict(pair=p,beta_deg=5.,chi_lower=lo,chi_upper=hi,relative_width=hi/lo-1,D_lower=dl,D_upper=dh,status='RESOLVED_BRACKET'))
        self.csv('additional_threshold_brackets',result);self.csv('additional_threshold_history',history)

    def extremum_evidence(self,rows):
        result=[];candidates=[]
        for p in self.cfg['cases']:
            for k in sorted({r['chi_u'] for r in rows if r['pair']==p}):
                rr=sorted([r for r in rows if (r['pair'],r['chi_u'])==(p,k)],key=lambda r:r['beta_deg']);w=np.array([r['omega1'] for r in rr]);b=np.array([r['beta_deg'] for r in rr]);j=int(np.argmax(w));tol=self.cfg['policy']['same_physics_frequency_tol']
                near=b[w>=max(w)*(1-2*tol)]
                extrema=[]
                for i in range(1,len(b)-1):
                    left=w[i]-w[i-1];right=w[i+1]-w[i]
                    if left*right<0:
                        kind='maximum' if left>0 else 'minimum';extrema.append(kind)
                        candidates.append(dict(pair=p,chi_u=k,kind=kind,beta_deg=b[i],left_deg=b[i-1],right_deg=b[i+1],omega1=w[i],neighbour_contrast=min(abs(left),abs(right))/w[i],neighbour_sign_resolved=min(abs(left),abs(right))/w[i]>2*tol))
                result.append(dict(pair=p,chi_u=k,beta_max_grid=b[j],omega_max=max(w),G_max_percent=100*(max(w)/w[0]-1),sampled_maximum_band_left=float(min(near)),sampled_maximum_band_right=float(max(near)),band_relative_frequency_width=2*tol,interior_maxima=extrema.count('maximum'),interior_minima=extrema.count('minimum'),D90=w[-1]/w[0]-1,min_gap12=min(r['gap12'] for r in rr)))
        self.csv('extremum_evidence',result);self.csv('all_extremum_candidates',candidates)

    def curvature(self,rows):
        reflections=[];fits=[];summary=[]
        S=np.diag([1,-1,-1,1,-1,-1]);SS=scipy.linalg.block_diag(S,S);R=np.diag([1,-1,1,-1,-1,-1]);E=np.diag([1,-1,-1])
        m=self.run.model
        for p in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                rr={float(r['beta_deg']):float(r['omega1']) for r in rows if r['pair']==p and r['chi_u']==k}
                for b in (0.,.25,1.,30.,90.):
                    cp=self.case(p,k,b);cm=self.case(p,k,-b)
                    je=np.max(abs(physics.joint_matrix(m,cm)@SS-R@physics.joint_matrix(m,cp)))
                    ee=max(np.max(abs(physics.external_matrix(m,cp,i)@S-E@physics.external_matrix(m,cp,i))) for i in (0,1))
                    A=physics.base.state_matrix(m,rr[0.]);ae=np.max(abs(A@S-S@A))
                    require(max(je,ee,ae)<1e-13,'Reflection identity')
                    reflections.append(dict(pair=p,chi_u=k,beta_deg=b,joint_error=je,external_error=ee,ode_error=ae,status='PASS'))
                current=[]
                for end in (.5,.75,1.):
                    b=np.array([v for v in self.cfg['curvature']['angles'] if v<=end]);x=np.radians(b)**2;y=np.log(np.array([rr[v]/rr[0.] for v in b]));xmax=max(x)
                    for degree in (1,2,3):
                        if degree>len(x):continue
                        A=np.column_stack([(x/xmax)**j for j in range(1,degree+1)]);pinv=np.linalg.pinv(A);coef=pinv@y
                        weights=pinv[0]/xmax;C=coef[0]/xmax;noise=2*self.cfg['curvature']['frequency_uncertainty_bound']*sum(abs(weights))
                        r=dict(pair=p,chi_u=k,max_angle_deg=end,polynomial_degree_in_beta=2*degree,nodes=len(b),C_full=C,frequency_error_bound_on_C=noise,fit_residual=float(max(abs(A@coef-y))),scaled_design_condition=float(np.linalg.cond(A)))
                        fits.append(r);current.append(r)
                accepted=[r for r in current if r['polynomial_degree_in_beta']>=4]
                ref=next(r for r in accepted if r['max_angle_deg']==1. and r['polynomial_degree_in_beta']==6)
                spread=max(abs(r['C_full']-ref['C_full']) for r in accepted);uncertainty=spread+max(r['frequency_error_bound_on_C'] for r in accepted)
                status='NEGATIVE' if ref['C_full']+uncertainty<0 else ('POSITIVE' if ref['C_full']-uncertainty>0 else 'UNRESOLVED')
                summary.append(dict(pair=p,chi_u=k,C_full=ref['C_full'],fit_spread=spread,conservative_uncertainty=uncertainty,C_lower=ref['C_full']-uncertainty,C_upper=ref['C_full']+uncertainty,status=status,min_gap_near_zero=min(r['gap12'] for r in rows if r['pair']==p and r['chi_u']==k and r['beta_deg']<=1.)))
        self.csv('reflection_checks',reflections);self.csv('curvature_fits',fits);self.csv('curvature_summary',summary)

    def consistency(self,rows):
        allrows=read_csv(self.folder/'reference_spectrum.csv')+rows+self.extra;mon=[];zero=[]
        tol=self.cfg['policy']['same_physics_frequency_tol']
        for p in self.cfg['cases']:
            z={float(r['chi_u']):float(r['omega1']) for r in allrows if r['pair']==p and float(r['beta_deg'])==0}
            ref=z[0.]
            for k,w in sorted(z.items()):
                require(abs(w/ref-1)<=tol,'Zero frequency regression');zero.append(dict(pair=p,chi_u=k,omega0=w,relative_error=w/ref-1,status='PASS'))
            for b in sorted({float(r['beta_deg']) for r in allrows if r['pair']==p}):
                rr=sorted([r for r in allrows if r['pair']==p and float(r['beta_deg'])==b],key=lambda r:float(r['chi_u']))
                for lo,hi in zip(rr,rr[1:]):
                    delta=min(float(hi[f'omega{i}'])/float(lo[f'omega{i}'])-1 for i in (1,2,3))
                    dlo=float(lo['omega1'])/z[float(lo['chi_u'])]-1;dhi=float(hi['omega1'])/z[float(hi['chi_u'])]-1
                    require(delta>=-tol and dhi-dlo>=-2*tol,'Monotonicity of spectrum or normalized growth')
                    mon.append(dict(pair=p,beta_deg=b,chi_lower=lo['chi_u'],chi_upper=hi['chi_u'],min_relative_spectral_increment=delta,D_increment=dhi-dlo,status='PASS'))
        self.csv('stiffness_monotonicity',mon);self.csv('all_zero_regression',zero)
        self.a.write(self.folder/'analysis_complete.json',dict(status='PASS',addressed_points=len(self.extra),monotonicity_comparisons=len(mon),zero_regressions=len(zero)))

    def finish(self,stage,status='PASS'):
        super().finish(stage,status)
        path=self.folder/'executions.json';r=json.loads(path.read_text());r[-1]['workflow_sha256']=self.a.sha(__file__);self.a.write(path,r)


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None and args.folder.resolve().is_relative_to(BASE.resolve()),'Explicit new archive required')
        with previous.previous.plots_without_calculation(),previous.ranges.saved.forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':render(args.folder)
            else:print(finalize(args.folder))
        return
    exp=Experiment();print(exp.folder,flush=True)
    try:
        require(not (exp.folder/'failure.json').exists(),'Saved failure needs diagnosis')
        if args.stage=='preflight':exp.preflight()
        else:exp.run_study()
        exp.finish(args.stage)
    except Exception as exc:
        exp.a.write(exp.folder/'failure.json',dict(status='PARTIAL',stage=args.stage,error=repr(exc)));exp.finish(args.stage,'PARTIAL');raise


def render(folder):
    """Render saved rows only; callers prohibit every spectral entry point."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    started=time.perf_counter();a=api();cfg=json.loads((folder/'contract.json').read_text())['config']
    mainrows=read_csv(folder/'spectrum.csv');diag=read_csv(folder/'diagnostic_spectrum.csv');refs=read_csv(folder/'reference_spectrum.csv');extra=read_csv(folder/'addressed_spectrum.csv')
    allrows=mainrows+diag+refs+extra;z={(r['pair'],float(r['chi_u'])):float(r['omega1']) for r in allrows if float(r['beta_deg'])==0.}
    growth=[dict(pair=r['pair'],chi_u=float(r['chi_u']),beta_deg=float(r['beta_deg']),omega1=float(r['omega1']),D_percent=100*(float(r['omega1'])/z[r['pair'],float(r['chi_u'])]-1)) for r in allrows]
    a.csvwrite(folder/'figure_growth_data.csv',growth)
    summ=read_csv(folder/'summary.csv')+read_csv(folder/'diagnostic_summary.csv');a.csvwrite(folder/'figure_maxima_data.csv',summ)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,'legend.fontsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.grid':True,'grid.alpha':.23})
    pairs=list(cfg['cases']);colors={.00024:'#1f77b4',.00027:'#e69f00',.000276:'#9467bd',.000283:'#009e73',.0003:'#d62728'}
    label=lambda k:fr'$\chi_u={k*1e4:g}\cdot10^{{-4}}$'.replace('.',',')
    fmt=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    def title(p):
        j,e=cfg['cases'][p];return fr'$\kappa_j={j:g},\quad\kappa_e={e:g}$'.replace('.',',')
    def finish(fig,name):
        for ax in fig.axes:ax.xaxis.set_major_formatter(fmt);ax.yaxis.set_major_formatter(fmt)
        fig.savefig(folder/(name+'.pdf'),bbox_inches='tight');fig.savefig(folder/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,p in zip(axs,pairs):
        for k in [1e-4]+cfg['chi_values']+[1e-3]:
            rr=sorted([r for r in growth if r['pair']==p and r['chi_u']==k],key=lambda r:r['beta_deg'])
            ax.plot([r['beta_deg'] for r in rr],[r['omega1'] for r in rr],color=colors.get(k,'#777777'),ls='-' if k in colors else '--',label=label(k),lw=1.8)
        ax.set(title=title(p),xlabel=r'$\beta$, град',ylabel=r'$\omega_1$',xlim=(0,90));ax.legend()
    finish(fig,FIGURES[0])
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained',sharey=True)
    for ax,p in zip(axs,pairs):
        for k,color in colors.items():
            rr=sorted([r for r in growth if r['pair']==p and r['chi_u']==k],key=lambda r:r['beta_deg'])
            ax.plot([r['beta_deg'] for r in rr],[r['D_percent'] for r in rr],color=color,label=label(k),lw=1.7,ls='-' if k in cfg['chi_values'] else '--')
        ax.axhline(0,color='black',lw=.8);ax.set(title=title(p),xlabel=r'$\beta$, град',ylabel=r'$100D$, %',xlim=(0,90));ax.legend(loc='lower left')
    finish(fig,FIGURES[1])
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained',sharey=True)
    for ax,p in zip(axs,pairs):
        for b,c in [(30.,'#1f77b4'),(60.,'#e69f00'),(90.,'#009e73')]:
            rr=sorted([r for r in growth if r['pair']==p and r['beta_deg']==b and .000237<=r['chi_u']<=.000317],key=lambda r:r['chi_u'])
            ax.plot([r['chi_u']*1e4 for r in rr],[r['D_percent'] for r in rr],'-o',color=c,ms=3,lw=1.3,label=fr'$\beta={b:g}^\circ$')
        ax.axhline(0,color='black',lw=.8);ax.set(title=title(p),xlabel=r'$\chi_u\cdot10^4$',ylabel=r'$100D$, %');ax.legend()
    finish(fig,FIGURES[2])
    fig,axs=plt.subplots(2,2,figsize=(12,7.5),layout='constrained')
    for col,p in enumerate(pairs):
        rr=sorted([r for r in summ if r['pair']==p],key=lambda r:float(r['chi_u']))
        for r in rr:
            k=float(r['chi_u']);x=k*1e4;marker='o' if k in cfg['chi_values'] else 'D'
            for row,n in [(0,'beta_max'),(1,'G_max_percent')]:axs[row,col].scatter(x,float(r[n]),color=colors[k],s=52,marker=marker,zorder=3)
        axs[0,col].set(title=title(p),ylabel=r'$\beta_{\max}^{\rm grid}$, град',ylim=(-3,80))
        axs[1,col].set(ylabel=r'$G_{\max}^{\rm grid}$, %',xlabel=r'$\chi_u\cdot10^4$')
    fig.suptitle('Круги — основные значения; ромбы — адресные проверки',fontsize=11)
    finish(fig,FIGURES[3])
    a.write(folder/'rendering.json',dict(status='PASS',seconds=time.perf_counter()-started,new_full_prefixes=0,dpi=320,plot_only=True,maximum_positions_not_connected=True))
    a.write(folder/'figure_provenance.json',{FIGURES[0]:['spectrum.csv','reference_spectrum.csv'],FIGURES[1]:['spectrum.csv','diagnostic_spectrum.csv'],FIGURES[2]:['figure_growth_data.csv','threshold_brackets.csv'],FIGURES[3]:['figure_maxima_data.csv']})


def finalize(folder):
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config']
    require(not (folder/'failure.json').exists(),'Unresolved numerical failure')
    for name in ('preflight','calculation_complete','diagnostic_preflight','diagnostic_complete','analysis_complete'):
        require(json.loads((folder/(name+'.json')).read_text())['status']=='PASS','Incomplete '+name)
    for name,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(a.sha(ROOT/name)==digest,'Protected source changed')
    a.write(folder/'source_integrity.json',{n:previous.ranges.saved.diagnostic.verify_source(ROOT/s['folder']) for n,s in cfg['sources'].items()})
    main=read_csv(folder/'spectrum.csv');diag=read_csv(folder/'diagnostic_spectrum.csv');extra=read_csv(folder/'addressed_spectrum.csv');refs=read_csv(folder/'reference_spectrum.csv');runs=json.loads((folder/'executions.json').read_text());q=[]
    for n in ('quality','diagnostic_quality','addressed_quality'):q+=read_csv(folder/(n+'.csv'))
    new=main+diag+extra;keys={(r['pair'],float(r['chi_u']),float(r['beta_deg'])) for r in new}
    require(len(keys)==len(new)==len(list((folder/'points').glob('*.json')))==sum(r['new_full_prefixes'] for r in runs),'Unique prefix accounting')
    require(sum(r['new_full_roots'] for r in runs)==3*len(new),'Root accounting')
    for r in new:
        path=ROOT/r['point_file'];require(a.sha(path)==r['point_sha256'] and a.sha(path.with_suffix('.npz'))==r['shape_sha256'],'Point provenance')
        rec=json.loads(path.read_text());require([v['counts'] for v in rec['roots']]==[[0,1],[1,2],[2,3]] and rec['count']['samples'][0]['count']==0,'Prefix completeness')
        require(all(all(v['gates'].values()) for v in rec['roots']),'Saved numerical gates')
    for n in ('high_intervals','diagnostic_high_intervals'):
        for r in read_csv(folder/(n+'.csv')):
            for key in ('entry_resolution_deg','exit_resolution_deg'):require(not r[key] or float(r[key])<=.25+1e-12,'Unresolved range boundary')
    for n in ('threshold_brackets','additional_threshold_brackets'):
        for r in read_csv(folder/(n+'.csv')):
            require(float(r['D_lower']) < -cfg['threshold_search']['sign_resolution'] and float(r['D_upper'])>cfg['threshold_search']['sign_resolution'] and float(r['relative_width'])<=.015,'Threshold bracket')
    for n in FIGURES:
        for ext in ('.pdf','.png'):require((folder/(n+ext)).exists(),'Missing figure')
    require(json.loads((folder/'tests.json').read_text())['status']=='PASS','Tests not accepted')
    for r in json.loads((folder/'figure_checks.json').read_text()):require(r['visual_review'].startswith('PASS') and r['embedded_images']==0 and min(r['png_dpi'])>=300,'Figure QA')
    curv=read_csv(folder/'curvature_summary.csv')+read_csv(folder/'diagnostic_curvature_summary.csv')
    decision=dict(status='PASS',scope='accepted spectra, extrema and fixed-angle brackets; curvature qualification retained',primary_constructions=6,diagnostic_constructions=4,primary_map_points=len(main),diagnostic_map_points=len(diag),addressed_points=len(extra),reused_points=len(refs),new_full_prefixes=len(new),new_full_roots=3*len(new),independent_ritz_solves=sum(r['new_ritz_solves'] for r in runs),new_leading_solves=0,new_pole_roots=sum(r['new_pole_roots'] for r in runs),quality_maxima={n:max(float(r[n]) for r in q) for n in ('physical','ode','energy','singular','mass_gram','condition')},curvature_unresolved=[dict(pair=r['pair'],chi_u=float(r['chi_u'])) for r in curv if r['status']=='UNRESOLVED'],monotonicity_comparisons=len(read_csv(folder/'stiffness_monotonicity.csv')),minimum_new_gap12=min(float(r['gap12']) for r in new),workflow_sha256=a.sha(__file__),entry_point_sha256=a.sha(ROOT/'scripts/analysis/study_two_springs_external_axial_spring.py'))
    a.write(folder/'decision.json',decision)
    a.write(folder/'manifest.json',dict(status='PASS',artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    return previous.ranges.saved.diagnostic.verify_source(folder)
