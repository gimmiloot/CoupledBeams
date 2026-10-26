"""Bounded thickness-asymmetry preset; orchestration, checks and saved-data plots."""
from dataclasses import asdict
from fractions import Fraction
import hashlib
import inspect
import json
import math
from pathlib import Path
import time
import numpy as np
from scripts.lib import two_springs_length_asymmetry as lengths
from scripts.lib import isotropic_two_springs_thickness_asymmetry as physics

ROOT=lengths.ROOT
CONFIG=ROOT/'data/input/isotropic_two_springs_thickness_asymmetry.json'
BASE=lengths.BASE/'thickness_asymmetry'
require=lengths.require
read_csv=lengths.read_csv
api=lengths.api


class Run(api().Run):
    def __init__(self,cfg,folder):
        super().__init__(cfg,folder,engine=physics)
        file=folder/'pole_catalogs.json'
        if file.exists():self.catalogs=json.loads(file.read_text())

    def solve(self,case,variant='RL_FULL',neighbour=None):
        self.model=physics.Model(**self.config['material'],mu=case.mu,tau=case.tau)
        return super().solve(case,variant,neighbour)

    def cats(self,model,case):
        result={}
        for arm,L in zip(model.arms,case.lengths(model)):
            key=physics.pole_key(arm,L)
            if key not in self.catalogs:
                p=self.policy
                c={'policy':{'scan_intervals':p['pole_scan_intervals'],'root_xtol':p['root_xtol'],'root_rtol':p['root_rtol']}}
                cat=api().tim_catalog(arm.carrier,L,model.omega(p['Lambda_upper']),c)
                cat['arm_contract']=dict(length=L,epsilon=arm.epsilon,EA=arm.section.EA,EI=arm.section.EI,
                    KGA=arm.section.KGA,rhoA=arm.section.rhoA,rhoI=arm.section.rhoI,J=arm.J)
                self.catalogs[key]=cat
                self.pole_roots+=len(cat['roots'])
            cat=self.catalogs[key]
            require(cat['arm_contract']==dict(length=L,epsilon=arm.epsilon,EA=arm.section.EA,EI=arm.section.EI,
                    KGA=arm.section.KGA,rhoA=arm.section.rhoA,rhoI=arm.section.rhoI,J=arm.J),'Arm catalog identity')
            result[key]=cat
        api().write(self.folder/'pole_catalogs.json',self.catalogs)
        return result


def polynomial_identity():
    """Exact coefficient check after common denominators, no symbolic dependency."""
    def add(*polys):
        out={}
        for p in polys:
            for k,v in p.items():out[k]=out.get(k,0)+v
        return {k:v for k,v in out.items() if v}
    def mul(p,q):
        out={}
        for (i,j),a in p.items():
            for (k,l),b in q.items():out[i+k,j+l]=out.get((i+k,j+l),0)+a*b
        return {k:v for k,v in out.items() if v}
    # Exponent pairs (power mu, power tau).
    d={(0,0):1,(0,2):1,(1,1):2}
    other={(0,0):2,(0,2):2,(1,1):-4}
    numerator=mul(d,other)
    target=add({(0,0):2,(0,2):-4,(0,4):2},{(0,2):8,(2,2):-8})
    require(numerator==target,'Symbolic compliance identity failed')
    volume=add(mul({(0,0):1,(1,0):-1},{(0,0):1,(0,1):-2,(0,2):1}),
               mul({(0,0):1,(1,0):1},{(0,0):1,(0,1):2,(0,2):1}))
    require(volume=={k:2*v for k,v in d.items()},'Symbolic volume identity failed')
    return dict(method='Exact multivariate polynomial coefficients over integers',
                volume_identity=True,compliance_identity=True,
                S_ax='2+8*(1-mu^2)*tau^2/(1-tau^2)^2',S_tau_at_zero=0,
                area_tau_at_zero=['-2*(1+mu)','2*(1-mu)'],
                inertia_tau_at_zero=['-4*(1+mu)','4*(1-mu)'],
                leading_energy='V_ax_normalized=xi^2*y^2/(2*S_ax)',
                leading_stiffness='K_added=xi^2/S_ax*P',C='P/(2*S_ax*zeta*M)',status='PASS')


class Experiment:
    validate=lengths.Experiment.validate

    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={}
        source=self.cfg['source'];self.source=ROOT/source['folder']
        require(self.a.sha(self.source/'manifest.json')==source['manifest_sha256'],'Source manifest changed')
        self.manifest=json.loads((self.source/'manifest.json').read_text())
        require(self.manifest['status']=='PASS','Source not accepted')
        self.inputs[(self.source/'manifest.json').relative_to(ROOT).as_posix()]=source['manifest_sha256']
        for name,digest in self.cfg['frozen_files'].items():
            require(self.a.sha(ROOT/name)==digest,'Frozen parent changed: '+name)
            self.inputs[name]=digest
        procedures=[api().Run,api().tim_catalog,lengths.length_slopes,api().refined_pick]
        identity=dict(config=self.cfg,physics_sha256=self.a.sha(physics.__file__),
                      procedures={f.__name__:hashlib.sha256(inspect.getsource(f).encode()).hexdigest() for f in procedures})
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        path=self.folder/'contract.json'
        if path.exists():require(json.loads(path.read_text())==identity,'Immutable contract')
        else:self.a.write(path,identity)
        self.run=Run(self.cfg,self.folder)
        self.archive={(r['pair'],float(r['mu']),float(r['beta_deg'])):r for r in self.source_csv('complete_spectrum.csv')}
        self.started=time.perf_counter()

    def source_csv(self,name):
        path=self.source/name
        require(self.a.sha(path)==self.manifest['artifacts'][name],'Source table hash')
        self.inputs[path.relative_to(ROOT).as_posix()]=self.a.sha(path)
        return read_csv(path)

    def point(self,pair,mu,tau,beta,previous=None,fresh=False,runner=None):
        case=physics.Case(float(beta),float(mu),*self.cfg['cases'][pair],float(tau))
        model=physics.Model(**self.cfg['material'],mu=mu,tau=tau)
        self.run.model=model
        if tau==0 and not fresh:
            row=self.archive[pair,mu,beta];path=ROOT/row['point_file']
            require(self.a.sha(path)==row['point_sha256'] and self.a.sha(path.with_suffix('.npz'))==row['shape_sha256'],'Archived point provenance')
            self.inputs[path.relative_to(ROOT).as_posix()]=row['point_sha256']
            self.inputs[path.with_suffix('.npz').relative_to(ROOT).as_posix()]=row['shape_sha256']
            rec=json.loads(path.read_text())
            self.validate(rec,physics.base.Case(case.beta,case.mu,case.k_joint,case.k_end))
            return rec,path,'archive'
        run=runner or self.run
        rec=run.solve(case,neighbour=previous)
        self.validate(rec,case)
        return rec,lengths.point_path(run.folder,case),'new_study'

    def finish(self,stage,status='PASS'):
        for name,digest in self.inputs.items():require(self.a.sha(ROOT/name)==digest,'Source changed '+name)
        f=self.folder/'executions.json';rows=json.loads(f.read_text()) if f.exists() else []
        rows.append(dict(stage=stage,status=status,seconds=time.perf_counter()-self.started,
                         new_prefixes=self.run.new_points,new_roots=self.run.new_roots,new_pole_roots=self.run.pole_roots,
                         workflow_sha256=self.a.sha(__file__)))
        self.a.write(f,rows)
        f=self.folder/'source_hashes.json';old=json.loads(f.read_text()) if f.exists() else {};old.update(self.inputs)
        self.a.write(f,old)


def geometry(exp):
    rows=[]
    for mu in exp.cfg['mu_values']:
        for tau in exp.cfg['tau_values']:
            m=physics.Model(**exp.cfg['material'],mu=mu,tau=tau)
            L=(1-mu,1+mu);arms=m.arms;a0,i0=m.section.area,m.section.second_moment
            area=np.array([a.section.area/a0 for a in arms]);inertia=np.array([a.section.second_moment/i0 for a in arms])
            require(abs(area@L-2)<2e-14 and sum(L)==2,'Volume/length conservation')
            exact=2+8*(1-mu*mu)*tau*tau/(1-tau*tau)**2
            require(abs(m.S_ax-exact)<3e-14 and min(area)>0,'Axial compliance')
            require(abs(inertia[1]/inertia[0]/((1+tau)/(1-tau))**4-1)<2e-14,'Inertia ratio')
            if tau==0:require(np.array_equal(area,[1,1]) and np.array_equal(inertia,[1,1]),'Equal-section limit')
            b_over_L=[math.sqrt(12)*a.epsilon*a.l*4/l for a,l in zip(arms,L)]
            require(max(b_over_L)<.1,'Transverse size is not small relative to length')
            rows.append(dict(mu=mu,tau=tau,L1=L[0],L2=L[1],A1_over_A0=area[0],A2_over_A0=area[1],
                I1_over_I0=inertia[0],I2_over_I0=inertia[1],volume_relative_error=area@L/2-1,
                S_ax=m.S_ax,added_stiffness_coefficient=1/m.S_ax,
                epsilon1_local=arms[0].epsilon/L[0],epsilon2_local=arms[1].epsilon/L[1],
                b1_over_L1=b_over_L[0],b2_over_L2=b_over_L[1],status='PASS'))
    exp.a.write(exp.folder/'algebra.json',polynomial_identity())
    exp.a.csvwrite(exp.folder/'geometry.csv',rows)
    print('Geometry and exact algebra PASS',flush=True)
    return rows


def preflight(exp):
    path=exp.folder/'preflight.json'
    if path.exists():
        require(json.loads(path.read_text())['status']=='PASS','Earlier preflight failed');return
    geometry(exp);reg=Run(exp.cfg,exp.folder/'regression');regression=[]
    for pair in exp.cfg['cases']:
        for mu in exp.cfg['regression_mu']:
            for beta in exp.cfg['regression_angles']:
                old,_,_=exp.point(pair,mu,0.,beta)
                new,_,_=exp.point(pair,mu,0.,beta,fresh=True,runner=reg)
                err=max(abs(np.array([r['omega'] for r in new['roots']])/np.array([r['omega'] for r in old['roots']])-1))
                require(err<=exp.cfg['policy']['same_physics_frequency_tol'],'tau=0 spectrum regression')
                case=physics.Case(beta,mu,*exp.cfg['cases'][pair],0.)
                m=physics.Model(**exp.cfg['material'],mu=mu)
                require(np.array_equal(physics.boundary(m,case,.08),physics.base.boundary(physics.base.Model(**exp.cfg['material']),case,.08)), 'tau=0 boundary matrix regression')
                regression.append(dict(pair=pair,mu=mu,beta_deg=beta,frequency_error=err,status='PASS'))
    exp.a.csvwrite(exp.folder/'regression.csv',regression)
    controls=[]
    for pair in exp.cfg['cases']:
        for mu,tau in exp.cfg['preflight_geometries']:
            for beta in exp.cfg['preflight_angles']:
                rec,_,_=exp.point(pair,mu,tau,beta)
                controls.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,Lambda=[r['Lambda'] for r in rec['roots']],status='PASS'))
    ritz=[]
    for pair in exp.cfg['cases']:
        for mu,tau,beta in exp.cfg['ritz_controls']:
            rec,_,_=exp.point(pair,mu,tau,beta)
            case=physics.Case(beta,mu,*exp.cfg['cases'][pair],tau)
            for order in exp.cfg['policy']['ritz_orders']:
                w,c,d=physics.ritz(exp.run.model,case,order)
                error=float(max(abs(w/np.array([r['omega'] for r in rec['roots']])-1)))
                require(error<=exp.cfg['policy']['ritz_frequency_tol'] and d['mass_gram_error']<=exp.cfg['policy']['mass_gram_tol'] and d['energy_error']<=exp.cfg['policy']['energy_tol'],'Independent Ritz gate')
                ritz.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,order=order,frequency_error=error,**d))
    exp.a.csvwrite(exp.folder/'preflight_ritz.csv',ritz)
    exp.a.write(path,dict(status='PASS',regression_max=max(r['frequency_error'] for r in regression),
                          regression_new_prefixes=reg.new_points,regression_pole_roots=reg.pole_roots,
                          heterogeneous_controls=controls,ritz_max=max(r['frequency_error'] for r in ritz),ritz_problems=len(ritz)))
    print(f'Preflight PASS: {len(regression)} tau=0 regressions, {len(controls)} heterogeneous controls, {len(ritz)} Ritz comparisons',flush=True)


def sweep(exp,constructions,angles,stage):
    rows,quality,neighbours=[],[],[]
    cache=exp.folder/'triggered_counts.json'
    checks=json.loads(cache.read_text()) if cache.exists() else {}
    for pair,mu,tau in constructions:
        previous=None
        for beta in angles:
            rec,path,origin=exp.point(pair,mu,tau,beta,previous)
            w=np.array([r['omega'] for r in rec['roots']])
            jump=float(max(abs(w/previous-1))) if previous is not None else 0.
            if jump>exp.cfg['policy']['neighbour_relative_jump']:
                key=f'{pair}|{mu}|{tau}|{beta}'
                if key not in checks:
                    case=physics.Case(float(beta),float(mu),*exp.cfg['cases'][pair],float(tau))
                    cats=exp.run.cats(exp.run.model,case)
                    test=[physics.count_at(exp.run.model,case,x,cats,exp.cfg['policy']) for x in (.9*w[0],(w[0]+w[1])/2,(w[1]+w[2])/2)]
                    require([v[1] for v in test]==[0,1,2],'Triggered count disagreement')
                    checks[key]=[r[2] for r in test]
                    exp.a.write(cache,checks)
            neighbours.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,jump=jump,
                               triggered=jump>exp.cfg['policy']['neighbour_relative_jump'],status='PASS'))
            previous=w
            rows.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,
                **{f'{key}{k}':r[key] for k,r in enumerate(rec['roots'],1) for key in ('omega','Omega','Lambda')},
                gap12_relative=w[1]/w[0]-1,origin=origin,point_file=path.relative_to(ROOT).as_posix(),
                point_sha256=exp.a.sha(path),shape_sha256=rec['shape_sha256']))
            for k,r in enumerate(rec['roots'],1):
                d=r['diagnostics']
                quality.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,sorted_position=k,
                    physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),
                    ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],
                    mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],
                    count_lower=r['counts'][0],count_upper=r['counts'][1],status='PASS'))
        for name,data in [('spectrum',rows),('quality',quality),('neighbours',neighbours)]:exp.a.csvwrite(exp.folder/(stage+'_'+name+'.csv'),data)
        print(f'{stage}: {pair}, mu={mu:g}, tau={tau:g}: {len(angles)} accepted points',flush=True)
    return rows,quality


def comparison(exp,rows):
    output,summary,differences=[],[],[]
    for pair in exp.cfg['cases']:
        for mu in exp.cfg['mu_values']:
            for tau in exp.cfg['tau_values']:
                rr=sorted([r for r in rows if (r['pair'],r['mu'],r['tau'])==(pair,mu,tau)],key=lambda r:r['beta_deg'])
                w0,lam0=rr[0]['omega1'],rr[0]['Lambda1']
                local=[]
                for r in rr:
                    F=r['omega1']/w0
                    require(lengths.close(F,(r['Lambda1']/lam0)**2,2e-14,0),'Normalization identity')
                    local.append(dict(**r,F_omega=F))
                output+=local
                at={r['beta_deg']:r for r in local};gap=min(local,key=lambda r:r['gap12_relative'])
                summary.append(dict(pair=pair,mu=mu,tau=tau,Lambda0=lam0,
                     **{f'Lambda{b}':at[b]['Lambda1'] for b in (3,6,10,20)},
                     **{f'F{b}':at[b]['F_omega'] for b in (3,5,10,20)},
                     min_gap12=gap['gap12_relative'],min_gap_beta=gap['beta_deg']))
    by={(r['pair'],r['mu'],r['tau'],r['beta_deg']):r for r in output}
    for pair in exp.cfg['cases']:
        for mu in exp.cfg['mu_values']:
            for beta in exp.cfg['map']['parameter_grid']:
                pos,neg=by[pair,mu,.1,beta],by[pair,mu,-.1,beta]
                differences.append(dict(pair=pair,mu=mu,beta_deg=beta,
                    D_Lambda=pos['Lambda1']/neg['Lambda1']-1,D_F=pos['F_omega']/neg['F_omega']-1,
                    D_omega=pos['omega1']/neg['omega1']-1))
    for name,data in [('spectrum',output),('summary',summary),('signed_comparison',differences)]:exp.a.csvwrite(exp.folder/(name+'.csv'),data)
    return output


def profiles(exp,rows,complete):
    out,transitions=[],[]
    for tau in exp.cfg['tau_values']:
        profile,transition=lengths.profiles(exp,[r for r in rows if r['tau']==tau],complete)
        out += [dict(**r,tau=tau) for r in profile]
        transitions += [dict(**r,tau=tau) for r in transition]
    return out,transitions


def mechanics(exp,rows):
    result=[]
    for r in rows:
        if r['beta_deg'] not in exp.cfg['mechanics_angles']:continue
        case=physics.Case(float(r['beta_deg']),r['mu'],*exp.cfg['cases'][r['pair']],r['tau'])
        m=physics.Model(**exp.cfg['material'],mu=r['mu'],tau=r['tau'])
        path=ROOT/r['point_file']
        require(exp.a.sha(path.with_suffix('.npz'))==r['shape_sha256'],'Field provenance')
        with np.load(path.with_suffix('.npz')) as data:
            arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
        mc,ec,sp=physics.energies(m,case,arms);mass=float(mc.sum());twov=float(ec.sum()+sp.sum())
        require(abs(mass-1)<1e-11 and abs(twov/(mass*r['omega1']**2)-1)<exp.cfg['policy']['energy_tol'],'Mechanical mass/energy gate')
        u,w,_,N,Q,_=arms[0]['ends'][1];wn,wq=N*w,-Q*u;s=(wn+wq)/twov
        result.append(dict(pair=r['pair'],mu=r['mu'],tau=r['tau'],beta_deg=r['beta_deg'],mass=mass,V=twov/2,
           d_joint=math.hypot(u,w),u1=u,w1=w,N1=N,Q1=Q,W_N=wn,W_Q=wq,s_per_radian=s,
           W_N_fraction=wn/(wn+wq) if r['beta_deg']>0 else None,
           arm1_mass=float(mc[0].sum()),arm2_mass=float(mc[1].sum()),RL_mass=float(mc[:,3].sum()),
           arm1_bend_fraction=ec[0,1]/twov,arm2_bend_fraction=ec[1,1]/twov,
           axial_fraction=float(ec[:,0].sum()/twov),shear_fraction=float(ec[:,2].sum()/twov)))
    exp.a.csvwrite(exp.folder/'mechanics.csv',result)
    return result


def leading_study(exp):
    cfg=exp.cfg['leading'];quality=cfg['quality'];cache=exp.folder/'leading_points';cache.mkdir(exist_ok=True)
    rows,checks,convergence,sensitivity,fd=[],[],[],[],[]
    new=0
    def solve(pair,mu,tau,xi,order=20,quad=64):
        nonlocal new
        ident=dict(pair=pair,mu=mu,tau=tau,xi=xi,order=order,quadrature=quad)
        key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];path=cache/(key+'.json')
        model=physics.Leading(*exp.cfg['cases'][pair],mu=mu,tau=tau,order=order,quadrature=quad)
        if path.exists():
            saved=json.loads(path.read_text());require(saved['identity']==ident,'Leading cache identity');result=saved['solution']
            for name in ('coefficients','zeta','omega_lim','Lambda','node_value','s_xi'):result[name]=np.array(result[name])
        else:
            result=model.solve(xi);new+=1;exp.a.write(path,dict(identity=ident,solution=result))
        d=result['diagnostics'];names=('symmetry','eigen_residual','energy','mass_gram','constraint','natural_boundary','strong_ode','essential')
        passed=all(d[k]<=quality[k+'_tol'] for k in names)
        checks.append(dict(**ident,**{k:d[k] for k in names},status='PASS' if passed else 'FAIL'))
        exp.a.csvwrite(exp.folder/'leading_quality.csv',checks)
        require(passed,'Leading physical gate: '+str(ident)+' '+str(d))
        return model,result
    old={(r['pair'],float(r['mu'])):r for r in exp.source_csv('leading_curvature.csv')}
    for pair in exp.cfg['cases']:
        for mu in exp.cfg['mu_values']:
            for tau in exp.cfg['tau_values']:
                m,r=solve(pair,mu,tau,0.,20)
                mf,rf=solve(pair,mu,tau,0.,24)
                mq,rq=solve(pair,mu,tau,0.,24,96)
                dw=float(max(abs(rf['omega_lim']/r['omega_lim']-1)));dc=abs(rf['C_at_zero']/r['C_at_zero']-1)
                dq=float(max(abs(rq['omega_lim']/rf['omega_lim']-1)));dcq=abs(rq['C_at_zero']/rf['C_at_zero']-1)
                require(dw<quality['frequency_convergence_rtol'] and dc<quality['C_convergence_rtol'] and dq<quality['frequency_convergence_rtol'] and dcq<quality['C_convergence_rtol'],'Leading convergence')
                c=rf['coefficients'][:,0];mass=float(c@mf.M@c);zeta=float(rf['zeta'][0]);y=float(mf.y@c)
                C=y*y/(2*mf.S_ax*zeta*mass)
                require(lengths.close(C,rf['C_at_zero'],1e-12),'Leading curvature formula')
                if tau==0:require(lengths.close(C,float(old[pair,mu]['C']),1e-10),'tau=0 leading C regression')
                local_errors=[]
                for xi in cfg['xi'][1:]:
                    ml,rl=solve(pair,mu,tau,xi,24)
                    estimate=math.log(float(rl['omega_lim'][0]/rf['omega_lim'][0]))/xi**2
                    local_errors.append(abs(estimate/C-1))
                require(max(local_errors)<1e-4,'Local leading curvature check')
                rows.append(dict(pair=pair,mu=mu,tau=tau,S_ax=mf.S_ax,zeta0=zeta,y0=y,mass=mass,C=C,
                    arm1_mass=float(c@(mf.area[0]*mf.arm_M[0])@c),arm2_mass=float(c@(mf.area[1]*mf.arm_M[1])@c),
                    arm1_bending=float(c@(mf.inertia[0]*mf.arm_K[0])@c),arm2_bending=float(c@(mf.inertia[1]*mf.arm_K[1])@c),
                    local_C_relative_max=max(local_errors)))
                convergence.append(dict(pair=pair,mu=mu,tau=tau,frequency_20_24=dw,C_20_24=dc,frequency_64_96=dq,C_64_96=dcq))
            for xi in cfg['tau_sensitivity_xi']:
                m,r=solve(pair,mu,0.,xi,24);terms=m.tau_derivative_at_zero(r)
                sensitivity.append(dict(pair=pair,mu=mu,xi=xi,zeta=float(r['zeta'][0]),**terms))
                errors=[]
                for h in cfg['tau_fd_steps']:
                    _,plus=solve(pair,mu,h,xi,24);_,minus=solve(pair,mu,-h,xi,24)
                    finite=float((plus['zeta'][0]-minus['zeta'][0])/(2*h))
                    err=abs(finite-terms['total'])/max(abs(terms['total']),1e-8)
                    errors.append(err);fd.append(dict(pair=pair,mu=mu,xi=xi,h_tau=h,finite_difference=finite,
                        variational=terms['total'],relative_error=err))
                require(errors[-1]<cfg['tau_fd_relative_tolerance'] and errors[-1]<errors[0]+1e-9,'Tau derivative refinement')
        print('Leading C and tau sensitivity PASS: '+pair,flush=True)
    for name,data in [('leading_curvature',rows),('leading_convergence',convergence),('tau_sensitivity',sensitivity),('tau_sensitivity_fd',fd)]:exp.a.csvwrite(exp.folder/(name+'.csv'),data)
    exp.a.write(exp.folder/'leading_accounting.json',dict(new_this_execution=new,saved_problems=len(list(cache.glob('*.json'))),status='PASS'))
    return rows,checks


FIGURES=('01_frequency','02_normalized_frequency','03_signed_differences','04_characteristic_angle','05_leading_curvature')


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter,MultipleLocator
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,
        'legend.fontsize':9,'pdf.fonttype':42,'lines.linewidth':1.5,'axes.linewidth':.7})
    fmt=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    pairs=['J0E4','J4E0'];mus=[-.1,0.,.1];taus=[-.1,0.,.1]
    labels=[r'$(\kappa_{\mathrm{j}};\kappa_{\mathrm{e}})=(0{,}1;10)$',r'$(\kappa_{\mathrm{j}};\kappa_{\mathrm{e}})=(10;0{,}1)$']
    colors=['#155b86','#222222','#c17119'];lines=['-','--','-.'];markers=['o','s','^']
    def arr(rows,key):return np.array([float(r[key]) for r in rows])
    def label(symbol,x):return '$'+symbol+'='+f'{x:g}'.replace('.','{,}')+'$'
    def style(ax):
        ax.grid(True,color='.88',lw=.5);ax.set_axisbelow(True)
        ax.xaxis.set_major_formatter(fmt);ax.yaxis.set_major_formatter(fmt)
    def save(fig,name):
        fig.savefig(folder/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(folder/(name+'.png'),dpi=320);plt.close(fig)
    data=read_csv(folder/'spectrum.csv')
    for name,key,ylabel in zip(FIGURES[:2],('Lambda1','F_omega'),(r'$\Lambda_1$',r'$F(\beta;\mu,\tau)$')):
        fig,axs=plt.subplots(3,2,figsize=(7.2,7.6),sharex=True,sharey=True)
        for col,pair in enumerate(pairs):
            for row,mu in enumerate(mus):
                ax=axs[row,col]
                for tau,color,line,marker in zip(taus,colors,lines,markers):
                    rr=[r for r in data if r['pair']==pair and float(r['mu'])==mu and float(r['tau'])==tau]
                    ax.plot(arr(rr,'beta_deg'),arr(rr,key),color=color,ls=line,marker=marker,mfc='white',markersize=3.3,markevery=7,label=label(r'\tau',tau))
                style(ax);ax.set_xlim(0,20);ax.set_xticks([0,5,10,15,20])
                ax.text(.04,.88,'абвгде'[row*2+col]+')  '+label(r'\mu',mu),transform=ax.transAxes)
                if row==0:ax.text(.5,1.06,labels[col],ha='center',transform=ax.transAxes,fontsize=9.5)
                if col==0:ax.set_ylabel(ylabel)
                if row==2:ax.set_xlabel(r'$\beta$, град.')
        fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False)
        fig.subplots_adjust(left=.1,right=.98,bottom=.08,top=.9,hspace=.12,wspace=.08)
        save(fig,name)
    data=read_csv(folder/'signed_comparison.csv')
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.2),sharex=True,sharey='row')
    for col,pair in enumerate(pairs):
        for row,key in enumerate(('D_Lambda','D_F')):
            ax=axs[row,col]
            for mu,color,line,marker in zip(mus,colors,lines,markers):
                rr=[r for r in data if r['pair']==pair and float(r['mu'])==mu]
                ax.plot(arr(rr,'beta_deg'),100*arr(rr,key),color=color,ls=line,marker=marker,mfc='white',markersize=3.3,markevery=7,label=label(r'\mu',mu))
            style(ax);ax.axhline(0,color='.5',lw=.6);ax.set_xlim(0,20);ax.set_xticks([0,5,10,15,20])
            ax.text(.05,.9,'абвг'[2*row+col]+')',transform=ax.transAxes)
            if col==0:ax.set_ylabel(r'$100D_\Lambda^{(\tau)}$, %' if row==0 else r'$100D_F^{(\tau)}$, %')
            if row==0:ax.text(.5,1.04,labels[col],ha='center',transform=ax.transAxes,fontsize=9.5)
            if row==1:ax.set_xlabel(r'$\beta$, град.')
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False)
    fig.subplots_adjust(left=.12,right=.98,bottom=.1,top=.88,hspace=.12,wspace=.1)
    save(fig,FIGURES[2])
    for name,file,key,ylabel in ((FIGURES[3],'transitions','beta_pl_deg',r'$\beta_{\mathrm{pl}}$, град.'),
                                 (FIGURES[4],'leading_curvature','C',r'$C(\mu,\tau)$')):
        data=read_csv(folder/(file+'.csv'))
        fig,axs=plt.subplots(1,2,figsize=(7.2,3.8),sharey=True)
        for col,pair in enumerate(pairs):
            ax=axs[col]
            for mu,color,line,marker in zip(mus,colors,lines,markers):
                rr=sorted([r for r in data if r['pair']==pair and float(r['mu'])==mu and (key!='beta_pl_deg' or r['status']=='PASS')],key=lambda r:float(r['tau']))
                ax.plot(arr(rr,'tau'),arr(rr,key),color=color,ls=line,marker=marker,mfc='white',label=label(r'\mu',mu))
            style(ax);ax.set_xticks([-.1,0,.1]);ax.set_xlim(-.11,.11);ax.set_xlabel(r'$\tau$')
            ax.text(.5,1.04,labels[col],ha='center',transform=ax.transAxes,fontsize=9.5)
            ax.text(.01,1.04,'аб'[col]+')',transform=ax.transAxes)
            if key=='beta_pl_deg':ax.yaxis.set_major_locator(MultipleLocator(.25))
        axs[0].set_ylabel(ylabel)
        fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False)
        fig.subplots_adjust(left=.12,right=.98,bottom=.17,top=.8,wspace=.1)
        save(fig,name)


def verify_bundle(folder,figures=True):
    m=json.loads((folder/'manifest.json').read_text());require(m['status']=='PASS','Unaccepted bundle')
    for name,digest in m['artifacts'].items():
        if not figures and Path(name).suffix in ('.pdf','.png'):continue
        require(api().sha(folder/name)==digest,'Changed result '+name)
    for name,digest in {**m['source_hashes'],**m['workflow_files']}.items():require(api().sha(ROOT/name)==digest,'Changed source '+name)
    return m


def bundle_readme(folder):
    text='''# Различие толщин при фиксированном объёме

[Научный отчёт](../../../../docs/theory/isotropic_two_springs_thickness_asymmetry.md)
содержит постановку, алгебраический вывод, численные результаты и ограничения.
Исходный контракт: `data/input/isotropic_two_springs_thickness_asymmetry.json`.

## Параметры и источники

RL_FULL/Тимошенко, две точные поступательные связи и две поворотные пружины.
ε0=0,005; b_i/h_i=4; ν=0,3; k_s=5/6; l=E=ρ=1.
μ и τ принимают −0,1; 0; 0,1. J0E4=(κj;κe)=(0,1;10), J4E0=(10;0,1).
Внешняя пружина всегда на стержне 2. Размерные пружины и частотный масштаб
определяются через неизменные EI0,A0,l. Объём A1L1+A2L2=2A0l.

Единственный непосредственный архивный источник:
`results/isotropic_two_springs_rl_timoshenko/length_asymmetry/ddf22cb06e8563e2/`.
Оттуда переиспользованы 516 точек τ=0 (с полями по сохранённым путям),
проверены частоты и ведущие C. Его README указывает более ранние источники.
Архивы не копируются и не изменяются; проверяются хеши JSON/NPZ/CSV.
Каталоги полюсов новой сборки вычислены для индивидуальных L_i и сечений.
Старый каталог однородной конструкции не подставляется без проверки.

## Файлы и единицы

CSV: UTF-8, десятичная точка, разделитель запятая. На графиках — десятичная запятая.
`pair,mu,tau,beta_deg` определяют конструкцию и угол в градусах.

| Файл | Содержание |
| --- | --- |
| `geometry.csv`, `algebra.json` | Сохранение объёма, локальные ε_i, b_i/L_i, S_ax и точная проверка полиномиальных тождеств |
| `spectrum.csv` | 918 основных точек, omega1..3, Omega=omega·t0, Lambda=sqrt(Omega), F_omega=omega1/omega1(0) каждой конструкции |
| `complete_spectrum.csv` | Все 1548 основных и условно добавленных хвостовых точек; происхождение и хеши JSON/NPZ |
| `summary.csv` | Контрольные частоты и F, минимум g12=(omega2−omega1)/omega1 на основной сетке |
| `signed_comparison.csv` | D_F=F(τ=+0,1)/F(τ=−0,1)−1; аналогичные D_Lambda и D_omega, без модуля |
| `profiles.csv`, `transitions.csv` | Производные на градус и радиан, q, последующая огибающая Q, подтверждённые дискретные beta_pl |
| `preliminary_*.csv`, `tail_authorization.json` | Основания продолжения до 90°, полученные после основной карты |
| `mechanics.csv` | Полная массовая нормировка M=1; индивидуальные массы и энергии, u,w,N,Q,V, d_joint, вариационная производная на радиан |
| `leading_curvature.csv` | ζ0=Lambda⁴, M0=1, y0, C=y0²/(2S_ax ζ0 M0), массы и изгибные энергии стержней |
| `tau_sensitivity.csv` | ∂ζ/∂τ при τ=0, фиксированных μ,ξ: изгибной, массовый и нулевой прямой продольный вклады |
| `tau_sensitivity_fd.csv` | Адресная проверка при τ=±0,001 и ±0,0005; это диагностические узлы, не расширение основной карты |
| `leading_quality.csv`, `leading_convergence.csv` | Сильные и естественные невязки, энергия, масса; порядки 20/24, квадратуры 64/96 |
| `main_quality.csv`, `tail_quality.csv`, `preflight*.json/csv` | Покорневые физические проверки, счёт и независимые сравнения Ритца |

Политика `frequency-map-v1`, `fast_plot`, `sorted_positions`, K_plot=2,
K_guard=3. Рисунки показывают первую позицию; вторая и третья сохранены
для контроля низшего спектра. Глобальное отслеживание ветвей отсутствует.
Все новые точки проходят уравнения, граничные и узловые условия, энергию,
массовую ортогональность, проверку простоты и энергетический счёт корней.
Полные значения диагностик находятся в `points/*.json`; поля — в NPZ.

Основная сетка: 0:0,5:3; 3:0,25:12; 13:1:20 градусов.
Хвост 22:2:90 разрешён отдельно после положительного устойчивого прироста
на 0–5° и предварительного перехода. У всех 18 конструкций он подтверждён.
Критерий, q*=0,15 и минимальная длина хвоста 30° сохранены.
Конечные разности повторяют исследование длин, включая стыки 3°,12°,20°;
при нуле использовано точное следствие отражательной симметрии для простой частоты.
Это дискретный критерий, без гарантии экстремумов между узлами.

## Рисунки и подписи

Пять векторных PDF с внедрёнными шрифтами и PNG 320 dpi.
Линии соединяют известные точки; сглаживания и интерполяции спектров нет.

1. `01_frequency`: **Первая собственная частота при перераспределении толщин.**
   Слева (0,1;10), справа (10;0,1); строки μ=−0,1;0;0,1. Три кривые
   соответствуют τ=−0,1;0;0,1 при одинаковом общем объёме и размерных пружинах.
2. `02_normalized_frequency`: **Относительное изменение первой частоты с углом.**
   Каждая кривая нормирована на omega1(0) той же конструкции.
   Обозначения полей соответствуют рисунку 1.
3. `03_signed_differences`: **Влияние знака асимметрии толщин.**
   Верхний ряд — 100D_Lambda, нижний — 100D_F; слева (0,1;10), справа (10;0,1).
   Положительный знак означает большее значение у конструкции с более толстым
   вторым стержнем по сравнению с более тонким вторым стержнем.
4. `04_characteristic_angle`: **Характерный угол по прежнему дискретному критерию.**
   Опора 0–5°, q*=0,15; все последующие расчётные узлы проверены до 90°.
   Линии соединяют три исследованных τ; прекращение изменения частоты не утверждается.
5. `05_leading_curvature`: **Коэффициент начального роста ведущей частоты.**
   ln(F_lim)=Cξ²+o(ξ²), ξ=β_rad/ε0. C вычислен по массово-нормированной форме
   при ξ=0. Равная податливость при противоположных τ не обеспечивает равенства C.

## Воспроизведение

Из корня проекта, Python с NumPy, SciPy, Matplotlib:

```powershell
python scripts/analysis/study_two_springs_thickness_asymmetry.py --stage gates
python scripts/analysis/study_two_springs_thickness_asymmetry.py --stage preflight
python scripts/analysis/study_two_springs_thickness_asymmetry.py
python scripts/analysis/study_two_springs_thickness_asymmetry.py --verify-cache BUNDLE
python scripts/analysis/study_two_springs_thickness_asymmetry.py --plot-only BUNDLE
python -m pytest tests/test_two_springs_thickness_asymmetry.py -q
```

BUNDLE=`results/isotropic_two_springs_rl_timoshenko/thickness_asymmetry/KEY`.
Использованный интерпретатор: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`.
Полный повторный запуск использует готовые точки; два режима чтения не вызывают решатели.
Физика и входной контракт определяют отдельный идентификатор. `manifest.json`
проверяет код, источники и результаты; `executions.json` хранит счётчики запусков.
Этот новый входной контракт требует индивидуальных каталогов и проверок полей;
отдельная точка запуска сохраняет прежние CLI и контрольные суммы их пакетов.
Данные локальные и исключены из Git; научный отчёт хранится отдельно в документации.
'''
    (folder/'README.md').write_text(text.replace('KEY',folder.name),encoding='utf8')


def finalize(folder):
    bundle_readme(folder)
    event=folder/'stopped.json'
    if event.exists():
        r=json.loads(event.read_text())
        if r['status']=='STOPPED':
            require('SectionProperties' in r['error'] and 'AttributeError' in r['error'],'Unresolved failure')
            r.update(status='RESOLVED_API_ATTRIBUTE',resolution='Use existing area/second_moment field names; before any spectral calculation.',physical_gate_failure=False)
            api().write(event,r)
    paths=[CONFIG,Path(__file__),Path(physics.__file__),ROOT/'scripts/analysis/study_two_springs_thickness_asymmetry.py']
    api().write(folder/'manifest.json',dict(status='PASS',schema='thickness-asymmetry-results-v1',
        numerical_summary=json.loads((folder/'numerical_summary.json').read_text()),
        source_hashes=json.loads((folder/'source_hashes.json').read_text()),
        workflow_files={p.relative_to(ROOT).as_posix():api().sha(p) for p in paths},
        artifacts={p.relative_to(folder).as_posix():api().sha(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='manifest.json'}))


def main(args):
    if args.verify_cache or args.plot_only:
        require(not(args.verify_cache and args.plot_only),'Choose a single saved-data operation')
        folder=(args.verify_cache or args.plot_only).resolve()
        require(folder.parent==BASE,'Wrong result family')
        verify_bundle(folder,not bool(args.plot_only))
        if args.plot_only:render(folder);verify_bundle(folder)
        print('Saved data verified; 0 eigenvalue solves',flush=True);return
    exp=Experiment();print(exp.folder,flush=True)
    try:
        if args.stage=='gates':geometry(exp)
        else:
            preflight(exp)
            if args.stage=='all':
                constructions=[(pair,mu,tau) for pair in exp.cfg['cases'] for mu in exp.cfg['mu_values'] for tau in exp.cfg['tau_values']]
                rows,quality=sweep(exp,constructions,exp.cfg['map']['parameter_grid'],'main')
                normalized=comparison(exp,rows)
                mech=mechanics(exp,rows)
                leading,leading_quality=leading_study(exp)
                preliminary,decisions=profiles(exp,rows,False)
                exp.a.csvwrite(exp.folder/'preliminary_profiles.csv',preliminary)
                exp.a.csvwrite(exp.folder/'preliminary_transitions.csv',decisions)
                eligible=[(r['pair'],r['mu'],r['tau']) for r in decisions if r['status']=='PROVISIONAL_ONLY']
                exp.a.write(exp.folder/'tail_authorization.json',dict(eligible=eligible,decisions=decisions,
                    trigger='Positive resolved reference gain and provisional persistent crossing on 0-20'))
                tail,tq=sweep(exp,eligible,exp.cfg['tail_grid'],'tail') if eligible else ([],[])
                completed=[r for r in rows+tail if (r['pair'],r['mu'],r['tau']) in eligible]
                pp,tt=profiles(exp,completed,True)
                tt += [r for r in decisions if (r['pair'],r['mu'],r['tau']) not in eligible]
                for name,data in [('profiles',pp),('transitions',tt),('complete_spectrum',rows+tail)]:exp.a.csvwrite(exp.folder/(name+'.csv'),data)
                exp.a.write(exp.folder/'numerical_summary.json',dict(status='PASS',main_points=len(rows),tail_points=len(tail),
                    new_map_points=sum(r['origin']=='new_study' for r in rows+tail),archive_points=sum(r['origin']=='archive' for r in rows+tail),
                    accepted_angles=sum(r['status']=='PASS' for r in tt),
                    quality_max={key:max(r[key] for r in quality+tq) for key in ('physical','ode','energy','singular','mass_gram','condition')},
                    leading_quality_max={key:max(r[key] for r in leading_quality) for key in ('strong_ode','natural_boundary','mass_gram','energy')},
                    mechanics_mass_error=max(abs(r['mass']-1) for r in mech)))
        exp.finish(args.stage)
        if args.stage=='all':render(exp.folder);finalize(exp.folder)
    except Exception as error:
        exp.a.write(exp.folder/'stopped.json',dict(status='STOPPED',error=repr(error)))
        exp.finish(args.stage,'STOPPED');raise
