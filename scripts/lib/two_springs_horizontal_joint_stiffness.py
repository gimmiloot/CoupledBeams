"""HS-LA4: bounded joint-stiffness extension of the frozen HS-L1 fields.

Exact zero-mechanism deflation and exact J/D constraints require a distinct
contract from HS-L1. No old source-hashed runner or physical library is edited.
The independent check below is an EB endpoint-inertia count, not RL/Timoshenko.
Entry: python -c "from scripts.lib.two_springs_horizontal_joint_stiffness import run; run()"
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.linalg import cholesky,eigh,null_space,solve_triangular,expm,solve
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading,RIGID
from scripts.lib.two_springs_horizontal_spectral_structure import ROOT,sha,write,csvwrite,require
from scripts.lib.two_springs_horizontal_saturation import coefficient_from_linear

CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_leading_joint_stiffness.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_leading_joint_stiffness'


class JointLeading(Leading):
    def __init__(self,k_joint,mu=.01,order=24,quadrature=80):
        super().__init__(.1,mu,order,quadrature)  # Reuse only unchanged fields/energies.
        require(k_joint==RIGID or k_joint>=0,'Invalid joint stiffness')
        self.k_joint=k_joint;self.base_constraint=self.constraint.copy()
        self.K0=self.Kb.copy() if k_joint==RIGID else self.Kb+k_joint*np.outer(self.jump,self.jump)
        if k_joint==RIGID:self.constraint=np.vstack((self.constraint,self.jump))
        self.Z=null_space(self.constraint)

    def solve(self,tau,count=3):
        saturated=tau==RIGID
        C=np.vstack((self.constraint,self.D)) if saturated else self.constraint
        Z=null_space(C);K=self.K0 if saturated else self.K0+tau*self.P
        zero=self.k_joint==0 and tau==0
        if zero:
            mechanism=np.zeros(len(self.D));mechanism[self.order]=self.lengths[1];mechanism[-1]=-self.lengths[1]
            norm=float(mechanism@self.M@mechanism)
            require(abs(norm/(self.lengths[1]**3/3)-1)<5e-13,'Mechanism mass')
            require(np.linalg.norm(K@mechanism)==0 and np.max(abs(C@mechanism))<1e-13,'Exact zero mechanism')
            mechanism/=math.sqrt(norm)
            Z=Z@null_space((mechanism@self.M@Z)[None,:])
        Kr=Z.T@K@Z;Mr=Z.T@self.M@Z
        cholesky((Mr+Mr.T)/2,lower=True)
        L=cholesky((Kr+Kr.T)/2,lower=True)  # Positive after analytic kernel deflation.
        B=solve_triangular(L,Mr,lower=True);B=solve_triangular(L,B.T,lower=True).T
        n=count-int(zero);v,c=eigh((B+B.T)/2,subset_by_index=[len(B)-n,len(B)-1]);v=v[::-1];c=c[:,::-1]
        z=1/v;c=Z@solve_triangular(L.T,c,lower=False)/np.sqrt(v)
        if zero:z=np.r_[0.,z];c=np.column_stack((mechanism,c))
        return self.diagnose(tau,z,c)

    def diagnose(self,tau,z,c):
        saturated=tau==RIGID;C=np.vstack((self.constraint,self.D)) if saturated else self.constraint
        Z=null_space(C);K=self.K0 if saturated else self.K0+tau*self.P
        Zbase=null_space(self.base_constraint)
        constraints=[];names=[]
        if self.k_joint==RIGID:constraints.append(self.jump);names.append('r_J')
        if saturated:constraints.append(self.D);names.append('r_D')
        per=[]
        for i,ev in enumerate(z):
            ci=c[:,i];res=K@ci-ev*(self.M@ci)
            reactions={}
            if constraints:
                B=Zbase.T@np.array(constraints).T
                r=np.linalg.lstsq(B,-Zbase.T@res,rcond=None)[0]
                reactions=dict(zip(names,map(float,r)))
                kkt=Zbase.T@res+B@r
                kkt_scale=np.linalg.norm(Zbase.T@K@ci)+ev*np.linalg.norm(Zbase.T@self.M@ci)+np.linalg.norm(B@r)
                kkt_error=float(np.linalg.norm(kkt)/max(kkt_scale,1.))
            else:kkt_error=0.
            J=float(self.jump@ci);D=float(self.D@ci)
            rJ=reactions.get('r_J',0. if self.k_joint==RIGID else self.k_joint*J)
            rD=reactions.get('r_D',0. if saturated else tau*D)
            e=np.array([x@ci for x in self.ends]);o=np.array([x@ci for x in self.outer])
            bc=np.array([o[1,2],e[0,2]+rJ,e[1,2]-rJ,o[1,3]+rD,e[0,3]-e[1,3]-rD])
            bscale=max(1.,np.max(abs(e)),np.max(abs(o)),abs(rJ),abs(rD))
            ni=di=0.
            for arm,L in enumerate(self.lengths):
                W=self.values[arm][0]@ci;W4=self.values[arm][4]@ci
                ni+=float(self.weights@((W4-ev*W)**2)*L)
                di+=float(self.weights@(W4**2+(ev*W)**2)*L)
            strong=math.sqrt(ni/di) if ev>0 else math.sqrt(ni)
            stiffness=float(ci@K@ci);energy=abs(stiffness/ev-1) if ev>0 else abs(stiffness)
            denom=np.linalg.norm(Z.T@K@ci)+ev*np.linalg.norm(Z.T@self.M@ci)
            residual=float(np.linalg.norm(Z.T@res)/denom) if ev>0 else float(np.linalg.norm(Z.T@res))
            PJ=(0. if self.k_joint==RIGID else self.k_joint*J*J/ev) if ev>0 else None
            PD=(0. if saturated else tau*D*D/ev) if ev>0 else None
            per.append(dict(J=J,J_squared=J*J,D=D,D_squared=D*D,r_J=rJ,r_D=rD,P_J=PJ,P_D=PD,
                P_bending=float(ci@self.Kb@ci/ev) if ev>0 else None,
                natural=float(max(abs(bc))/bscale),natural_rows=bc.tolist(),boundary_scale=float(bscale),
                strong=strong,energy=energy,matrix_residual=residual,kkt=kkt_error))
        diag=dict(symmetry=float(max(np.max(abs(K-K.T)),np.max(abs(self.M-self.M.T)))),
            mass_gram=float(np.max(abs(c.T@self.M@c-np.eye(len(z))))),essential=float(np.max(abs(C@c))),
            natural_boundary=max(v['natural'] for v in per),strong_ode=max(v['strong'] for v in per),
            energy=max(v['energy'] for v in per),matrix_residual=max(v['matrix_residual'] for v in per),kkt=max(v['kkt'] for v in per),
            nullity=int(self.k_joint==0 and tau==0),constraint_rank=int(np.linalg.matrix_rank(C)),dimension=Z.shape[1],
            mass_positive=True,stiffness_positive=not(self.k_joint==0 and tau==0),zero_norm_convention='absolute at exact zero; HS-L1 relative otherwise')
        return dict(zeta=np.array(z),coefficients=c,diagnostics=diag,per_mode=per)


def ode_boundary(zeta,k,tau,mu=.01):
    A=np.zeros((4,4));A[0,1]=A[1,2]=A[2,3]=1;A[3,0]=zeta
    o=[np.eye(8)[:4],np.eye(8)[4:]];e=[expm(A*L)@x for L,x in zip((1-mu,1+mu),o)]
    j=e[0][1]-e[1][1];d=o[1][0]-e[1][0]
    rows=[o[0][0],o[0][1],e[0][0]+e[1][0],o[1][2]]
    rows+= [j,e[0][2]+e[1][2]] if k==RIGID else [e[0][2]+k*j,e[1][2]-k*j]
    rows+= [d,e[0][3]-e[1][3]+o[1][3]] if tau==RIGID else [o[1][3]+tau*d,e[0][3]-e[1][3]-tau*d]
    B=np.array(rows);return B/np.linalg.norm(B,axis=1)[:,None]


def endpoint_count(zeta,k,tau,mu=.01,zero_tol=1e-11,pole_guard=1e-8):
    """Wittrick-Williams inertia for exact EB endpoint energies.

    J(z)=sum(fixed-fixed counts)+negative inertia on admissible endpoint traces.
    Zero mechanisms are included whenever z>0. No determinant sign count.
    """
    require(zeta>0,'Count only at positive non-root probes')
    K=np.zeros((8,8));fixed=0;nearest=math.inf
    for arm,L in enumerate((1-mu,1+mu)):
        alpha=zeta**.25*L;n=1
        while (n+.5)*math.pi-.1<alpha+.2:
            pole=brentq(lambda x:math.cos(x)-1/math.cosh(x),(n+.5)*math.pi-.1,(n+.5)*math.pi+.1)
            distance=abs(alpha-pole)/max(1,alpha);nearest=min(nearest,distance)
            if pole<alpha:fixed+=1
            n+=1
        require(nearest>pole_guard,'Endpoint auxiliary pole too close')
        A=np.zeros((4,4));A[0,1]=A[1,2]=A[2,3]=1;A[3,0]=zeta;T=expm(A*L)
        X=np.zeros((4,4));X[:2,:2]=np.eye(2)
        X[2:]=solve(T[:2,2:],np.column_stack((-T[:2,:2],np.eye(2))))
        Y=T@X;beam=np.vstack((X[3],-X[2],-Y[3],Y[2]))
        K[arm*4:arm*4+4,arm*4:arm*4+4]=(beam+beam.T)/2
    unit=np.eye(8);C=[unit[0],unit[1],unit[2]+unit[6]]
    j=unit[3]-unit[7];d=unit[4]-unit[6]
    if k==RIGID:C.append(j)
    else:K+=k*np.outer(j,j)
    if tau==RIGID:C.append(d)
    else:K+=tau*np.outer(d,d)
    Z=null_space(C);ev=np.linalg.eigvalsh(Z.T@K@Z);scale=max(1.,max(abs(ev)))
    require(min(abs(ev))>zero_tol*scale,'Endpoint inertia unresolved near root')
    return dict(count=int(fixed+sum(ev<0)),fixed_fixed_count=fixed,negative_inertia=int(sum(ev<0)),
        inertia_margin=float(min(abs(ev))/scale),nearest_pole_distance=None if math.isinf(nearest) else nearest)


def check_prefix(k,tau,z,policy,refine=False):
    rows=[]
    for index,value in enumerate(z):
        if value==0:
            rank=int(np.linalg.matrix_rank(ode_boundary(0,k,tau)))
            require(rank==7,'ODE zero multiplicity')
            count=endpoint_count(1e-5,k,tau,zero_tol=policy['inertia_relative_zero_tol'])
            require(count['count']==1,'Missing/extra near-zero mode')
            rows.append(dict(position=index+1,zeta=0.,exact_zero=True,zero_nullity=8-rank,count_above=count['count']));continue
        h=policy['count_bracket_relative']*value;lo=value-h;hi=value+h
        ca=endpoint_count(lo,k,tau,zero_tol=policy['inertia_relative_zero_tol'],pole_guard=policy['pole_relative_guard'])
        cb=endpoint_count(hi,k,tau,zero_tol=policy['inertia_relative_zero_tol'],pole_guard=policy['pole_relative_guard'])
        require(ca['count']==index and cb['count']==index+1,'Exact EB prefix count mismatch')
        row=dict(position=index+1,zeta=float(value),exact_zero=False,lower=lo,upper=hi,count_lower=ca,count_upper=cb)
        if refine:
            root=brentq(lambda x:np.linalg.det(ode_boundary(x,k,tau)),lo,hi,xtol=1e-13,rtol=2e-14)
            error=abs(math.sqrt(root/value)-1);require(error<policy['independent_frequency_tol'],'Independent ODE frequency')
            row.update(zeta_ode=float(root),frequency_error=error)
        rows.append(row)
    return rows


def identity(cfg):
    record=dict(config=cfg,source_sha256=sha(__file__),config_sha256=sha(CONFIG),initial_sha256=sha(ROOT/cfg['initial_workspace']),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


class Study:
    def __init__(self,cfg):
        self.cfg=cfg;self.id,self.identity=identity(cfg);self.folder=BASE/self.id;self.folder.mkdir(parents=True,exist_ok=True)
        self.models={};self.solutions={};self.records=[];self.new=0;self.reused=0;self.counts=[]
        self.source=ROOT/'results/isotropic_two_springs_rl_timoshenko'/cfg['source_archives'][1]
        first=ROOT/'results/isotropic_two_springs_rl_timoshenko'/cfg['source_archives'][0]
        self.old_thresholds=json.loads((first/'leading_thresholds.json').read_text())

    def model(self,k,p):
        if (k,p) not in self.models:self.models[k,p]=JointLeading(k,self.cfg['mu'],p,self.cfg['quadrature'])
        return self.models[k,p]

    def low(self,k,p,tau,role):
        key=(k,p,tau)
        if key in self.solutions:return self.solutions[key]
        m=self.model(k,p);oldkey=hashlib.sha256(repr((p,tau)).encode()).hexdigest()[:16]
        old=self.source/'points'/f'{oldkey}.json'
        if k==.1 and old.exists():
            r=json.loads(old.read_text())
            with np.load(old.with_suffix('.npz')) as data:s=m.diagnose(tau,np.array(r['zeta']),data['coefficients'].copy())
            origin='HS-LA2';self.reused+=1
        else:s=m.solve(tau);origin='new';self.new+=1
        record=dict(k_joint=k,order=p,tau=tau,role=role,origin=origin,zeta=s['zeta'].tolist(),diagnostics=s['diagnostics'],per_mode=s['per_mode'])
        name=hashlib.sha256(repr(key).encode()).hexdigest()[:16];write(self.folder/'points'/f'{name}.json',record)
        np.savez_compressed(self.folder/'points'/f'{name}.npz',coefficients=s['coefficients'])
        self.records.append(record);pol=self.cfg['policy']
        for name,tol in [('strong_ode','strong_tol'),('natural_boundary','boundary_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol'),('kkt','matrix_tol')]:
            require(s['diagnostics'][name]<pol[tol],f'{key}: {name} failed')
        require(np.all(np.diff(s['zeta'])>0) and s['zeta'][0]>=0,'Sorted physical positions')
        if p==24:
            count=check_prefix(k,tau,s['zeta'],pol,role=='primary')
            self.counts.append(dict(k_joint=k,tau=tau,role=role,checks=count));write(self.folder/'independent_checks.json',self.counts)
        self.solutions[key]=s;return s

    def threshold(self,k,p,q,role='primary'):
        m=self.model(k,p);lim=self.low(k,p,RIGID,role);target=q*q*lim['zeta'][0]
        K=m.Z.T@m.K0@m.Z;M=m.Z.T@m.M@m.Z;d=m.Z.T@m.D
        A=K-target*M;y=solve(A,d,assume_a='sym');G=float(d@y)
        require(G<0,'Target not in active first secular interval')
        tau=-1/G
        if k==.1 and role=='primary':
            archived=self.old_thresholds[str(q)]['tau']
            require(abs(tau/archived-1)<1e-9,'HS-L1 threshold regression');tau=archived
        s=self.low(k,p,tau,role)
        error=float(s['zeta'][0]/target-1);require(abs(error)<self.cfg['policy']['threshold_residual_tol'],'Threshold equation')
        if p==24:
            lower=endpoint_count(target,k,tau*(1-1e-5));upper=endpoint_count(target,k,tau*(1+1e-5))
            require(lower['count']==1 and upper['count']==0,'Threshold first-root isolation')
        return tau,s,dict(target_zeta=float(target),level_residual=error,G_target=G,
            linear_residual=float(np.linalg.norm(A@y-d)/(np.linalg.norm(A@y)+np.linalg.norm(d))))

    def coefficients(self,k,p,lim):
        m=self.model(k,p);z=float(lim['zeta'][0]);r=lim['per_mode'][0]['r_D']
        Z=m.Z;K=Z.T@m.K0@Z;M=Z.T@m.M@Z;d=Z.T@m.D
        linear,u,v=coefficient_from_linear(K,M,d,z)
        require(max(abs(linear['G']),linear['residual_u'],linear['residual_v'])<1e-9,'Saturation resolvent')
        require(abs(linear['a']/r**2-1)<self.cfg['policy']['reaction_relative_tol'],'Reaction / resolvent a')
        key=hashlib.sha256(repr((k,p)).encode()).hexdigest()[:16];np.savez_compressed(self.folder/f'linear_{key}.npz',K0=K,M=M,d=d,u=u,v=v)
        return dict(a=r*r,b=linear['b'],reaction=r,linear=linear)

    def run(self):
        cfg=self.cfg;state=json.loads((ROOT/cfg['initial_workspace']).read_text())
        for name,digest in cfg['frozen_files'].items():require(sha(ROOT/name)==digest,'Frozen source '+name)
        for name,digest in state['archive_files'].items():require(sha(ROOT/name)==digest,'Historical archive '+name)
        require(sha(ROOT/cfg['attachment']['copy_path'])==cfg['attachment']['sha256'],'Attachment')
        if (self.folder/'manifest.json').exists():
            manifest=json.loads((self.folder/'manifest.json').read_text());require(manifest['identity']==self.identity,'Identity')
            for name,digest in manifest['artifacts'].items():require(sha(self.folder/name)==digest,'Artifact '+name)
            print('HS-LA4 verified',self.folder);return self.folder
        write(self.folder/'contract.json',self.identity)
        try:
            summary=[];sensitivities=[];coefficients=[];forms=[]
            # Hard gate: kernel and rigid cantilever before positive threshold work.
            for k in [0,RIGID]:
                for p in cfg['orders']:
                    s=self.low(k,p,0.,'primary')
                    if k==0:
                        require(s['zeta'][0]==0 and s['diagnostics']['nullity']==1,'Zero hard gate')
                        require(abs(s['per_mode'][0]['J_squared']/(3/1.01**3)-1)<1e-12,'Zero J')
                        require(abs(s['per_mode'][0]['D_squared']/(3/1.01)-1)<1e-12,'Zero D')
                    else:require(abs(s['zeta'][0]/cfg['cantilever_zeta']-1)<1e-9,'Rigid cantilever')
            print('zero and cantilever gates PASS',flush=True)
            for k in cfg['k_joint']:
                for p in cfg['orders']:
                    zero=self.low(k,p,0.,'primary');lim=self.low(k,p,RIGID,'primary');coef=self.coefficients(k,p,lim)
                    coefficients.append(dict(k_joint=k,order=p,**coef))
                    row=dict(k_joint=k,order=p,zeta0=float(zero['zeta'][0]),zeta_infty=float(lim['zeta'][0]),
                        omega_infty=cfg['epsilon']*math.sqrt(lim['zeta'][0]),a=coef['a'],b=coef['b'])
                    for q in cfg['levels']:
                        tau,s,check=self.threshold(k,p,q);f=s['per_mode'][0];inf=lim['per_mode'][0]
                        derivative=None if k==RIGID else (q*q*inf['J_squared']-f['J_squared'])/f['D_squared']
                        log=None if k in [0,RIGID] else k*derivative/(2*tau)
                        energylog=None if k in [0,RIGID] else (inf['P_J']-f['P_J'])/(2*f['P_D'])
                        if log is not None:require(abs(log-energylog)<1e-10,'Energy sensitivity identity')
                        sensitivities.append(dict(k_joint=k,order=p,q=q,tau=tau,J_squared=f['J_squared'],D_squared=f['D_squared'],
                            J_infty_squared=inf['J_squared'],D_infty_squared=inf['D_squared'],P_J=f['P_J'],P_J_infty=inf['P_J'],P_D=f['P_D'],
                            d_tau_d_k=derivative,d_log_beta_d_log_k=log,log_sensitivity_energy=energylog,
                            beta_H5_deg=180/math.pi*cfg['epsilon']*math.sqrt(4*tau),beta_H6_deg=180/math.pi*cfg['epsilon']*math.sqrt(2*tau),**check))
                        row['tau90' if q==.9 else 'tau95']=tau
                    summary.append(row)
                print('primary',k,'PASS',flush=True)
            csvwrite(self.folder/'summary.csv',summary);csvwrite(self.folder/'sensitivities.csv',sensitivities);write(self.folder/'coefficients.json',coefficients)
            differences=[]
            for k in cfg['finite_difference']['joints']:
                for h in (cfg['finite_difference']['zero_steps'] if k==0 else [k*x for x in cfg['finite_difference']['relative_steps']]):
                    for q in cfg['levels']:
                        ref=next(r for r in sensitivities if r['k_joint']==k and r['order']==24 and r['q']==q)
                        if k==0:
                            t1=self.threshold(h,24,q,'derivative_check')[0];t2=self.threshold(2*h,24,q,'derivative_check')[0]
                            fd=(-3*ref['tau']+4*t1-t2)/(2*h);method='forward_second_order'
                        else:
                            a=self.threshold(k+h,24,q,'derivative_check')[0];b=self.threshold(k-h,24,q,'derivative_check')[0]
                            fd=(a-b)/(2*h);method='central'
                        err=abs(fd-ref['d_tau_d_k'])/max(abs(ref['d_tau_d_k']),1e-12)
                        differences.append(dict(k_joint=k,q=q,step=h,method=method,analytic=ref['d_tau_d_k'],difference=fd,relative_error=err))
                        require(err<cfg['policy']['sensitivity_relative_tol'],'Derivative difference')
            csvwrite(self.folder/'finite_differences.csv',differences)
            convergence=[]
            for k in cfg['k_joint']:
                a,b=[next(r for r in summary if r['k_joint']==k and r['order']==p) for p in cfg['orders']]
                for quantity in ['zeta_infty','a','b','tau90','tau95']:
                    error=abs(a[quantity]/b[quantity]-1);require(error<cfg['policy']['coefficient_relative_tol'],'Order convergence '+quantity)
                    convergence.append(dict(k_joint=k,quantity=quantity,relative_error=error))
                for label in ['zero','infty','q90','q95']:
                    sols=[]
                    for p in cfg['orders']:
                        row=next(r for r in summary if r['k_joint']==k and r['order']==p)
                        t=0. if label=='zero' else RIGID if label=='infty' else row['tau90' if label=='q90' else 'tau95']
                        sols.append(self.solutions[k,p,t])
                    aa,bb=sols;positive=aa['zeta']>0
                    require(np.max(abs(np.sqrt(aa['zeta'][positive]/bb['zeta'][positive])-1))<cfg['policy']['frequency_convergence_tol'],'Prefix convergence')
                    m20=self.model(k,20);m24=self.model(k,24);va=[];vb=[];w=[]
                    for arm,L in enumerate(m24.lengths):
                        va.append(m20.field_matrices(arm,m24.nodes)[0]@aa['coefficients']);vb.append(m24.values[arm][0]@bb['coefficients']);w.extend(m24.weights*L)
                    va=np.vstack(va);vb=np.vstack(vb);w=np.array(w);sign=np.sign(np.sum(w[:,None]*va*vb,axis=0))
                    shape=float(max(np.sqrt(np.sum(w[:,None]*(va-vb*sign)**2,axis=0))))
                    require(shape<cfg['policy']['shape_convergence_tol'],'Shape convergence')
                    convergence.append(dict(k_joint=k,quantity=label+'_shape',relative_error=shape))
                    if label!='zero':
                        for arm,L in enumerate(m24.lengths):
                            vals=[x@bb['coefficients'][:,0] for x in m24.values[arm][:4]]
                            for i,x in enumerate(m24.nodes):forms.append(dict(k_joint=k,state=label,arm=arm+1,s=x*L,W=vals[0][i],W1=vals[1][i],W2=vals[2][i],W3=vals[3][i]))
            csvwrite(self.folder/'convergence.csv',convergence);csvwrite(self.folder/'first_forms.csv',forms)
            write(self.folder/'point_records.json',self.records)
            write(self.folder/'manifest.json',dict(status='PASS',identity=self.identity,new_ritz_prefixes=self.new,reused_ritz_prefixes=self.reused,
                independent_count_points=len(self.counts),independent_ode_positive_roots=sum('zeta_ode' in r for a in self.counts for r in a['checks']),
                new_full_model_solves=0,artifacts={p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file()}))
        except Exception as e:
            write(self.folder/'failure.json',dict(status='PARTIAL_PASS',error=repr(e),new_ritz_prefixes=self.new,reused_ritz_prefixes=self.reused));raise
        print('HS-LA4 PASS',self.folder);return self.folder


def run():return Study(json.loads(CONFIG.read_text())).run()
