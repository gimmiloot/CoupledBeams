"""HS-LA2 finite-space rank-one audit, reusing frozen HS-L1 assembly.

Script Proliferation Control: full-basis secular/derivative/KKT contract is
different from HS-L1 threshold validation. This stable diagnostic module
contains all its helpers; no new physical library or universal solver.
No imports or calls of the full RL/Timoshenko solver are needed.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import platform
import time
import numpy as np
import scipy
import mpmath as mp
from scipy.linalg import solve
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading,RIGID

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_leading_spectral_structure.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_leading_spectral_structure'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf8')


def csvwrite(path,rows):
    if not rows:return
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def identity(cfg):
    record=dict(config=cfg,files={str(CONFIG.relative_to(ROOT)):sha(CONFIG),str(Path(__file__).relative_to(ROOT)):sha(__file__)},
        initial_sha=sha(ROOT/cfg['initial_workspace']),versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def all_pairs(model,tau=0,dps=80):
    """All N finite Ritz pairs; high precision only changes arithmetic.

    Float matrices are promoted exactly. No change of quadrature, fields or
    constraints. The inverse triangular factor is built by forward substitution;
    scalar compliance below always uses a linear solve, never an explicit inverse.
    """
    Z=model.Z;K=Z.T@(model.K0+tau*model.P)@Z;M=Z.T@model.M@Z
    # Frozen assembly matrices are symmetric up to arithmetic roundoff.
    K=(K+K.T)/2;M=(M+M.T)/2
    with mp.workdps(dps):
        k=mp.matrix(K.tolist());m=mp.matrix(M.tolist());L=mp.cholesky(m);N=len(K);T=mp.zeros(N)
        for j in range(N):
            for i in range(N):T[i,j]=((1 if i==j else 0)-sum(L[i,r]*T[r,j] for r in range(i)))/L[i,i]
        z,v=mp.eigsy(T*k*T.T);c=T.T*v
        mass=c.T*m*c-mp.eye(N);energy=c.T*k*c
        gram=max(abs(x) for x in mass);matrix=mp.mpf(0);er=mp.mpf(0)
        for j in range(N):
            a=k*c[:,j];b=z[j]*m*c[:,j]
            matrix=max(matrix,mp.norm(a-b)/(mp.norm(a)+mp.norm(b)))
            er=max(er,abs(energy[j,j]/z[j]-1))
        dn=mp.matrix((Z.T@model.D).tolist()).T*c
        authoritative=dict(zeta=[mp.nstr(x,dps) for x in z],coefficients=[[mp.nstr(c[i,j],dps) for j in range(N)] for i in range(N)],d_n=[mp.nstr(x,dps) for x in dn])
        dn=np.array(list(dn),float);z=np.array(list(z),float);c=np.array(c.tolist(),float)
        diagnostic=dict(dps=dps,N=N,mass_gram_mp=float(gram),matrix_mp=float(matrix),energy_mp=float(er))
    residual=K@c-(M@c)*z;den=np.linalg.norm(K@c,axis=0)+z*np.linalg.norm(M@c,axis=0)
    diagnostic.update(mass_gram_float=float(np.max(abs(c.T@M@c-np.eye(len(z))))),
        matrix_float=float(np.max(np.linalg.norm(residual,axis=0)/den)),
        energy_float=float(max(abs(np.diag(c.T@K@c)/z-1))))
    diagnostic['float_coefficients_role']='non-authoritative rounded export; high-mode cancellation is diagnosed, not accepted as a continuum residual'
    return dict(zeta=z,coefficients=c,d_n=dn,authoritative=authoritative,diagnostics=diagnostic,K=K,M=M,d=Z.T@model.D)


def compliance_linear(K,M,d,z):
    y=solve(K-z*M,d,assume_a='sym')
    return float(d@y),float(y@M@y)


def compliance_modal(z0,dn,z):
    return float(np.sum(dn*dn/(z0-z))),float(np.sum(dn*dn/(z0-z)**2))


def reaction(model,solution):
    c=solution['coefficients'][:,0].copy();z=float(solution['zeta'][0]);Z=model.Z
    K=Z.T@model.K0@Z;M=Z.T@model.M@Z;d=Z.T@model.D;v=Z.T@c
    residual=K@v-z*(M@v);r=-float(d@residual/(d@d))
    if r<0:c=-c;v=-v;residual=-residual;r=-r
    scale=np.linalg.norm(K@v)+z*np.linalg.norm(M@v)+abs(r)*np.linalg.norm(d)
    e=np.array([a@c for a in model.ends]);o=np.array([a@c for a in model.outer])
    return c,dict(zeta_infty=z,r=r,r_squared=r*r,kkt_residual=float(np.linalg.norm(residual+r*d)/scale),
        constraint=float(abs(model.D@c)),mass=float(c@model.M@c),r_outer=float(-o[1,3]),r_joint=float(e[0,3]-e[1,3]))


class Audit:
    def __init__(self,cfg):
        self.cfg=cfg;self.identifier,self.ident=identity(cfg);self.folder=BASE/self.identifier
        self.folder.mkdir(parents=True,exist_ok=True);self.models={p:Leading(cfg['k_joint'],cfg['mu'],p,cfg['quadrature']) for p in cfg['orders']}
        self.solutions={};self.records=[];self.bases={};self.basis_records=[];self.reused=0;self.new=0
        self.source=ROOT/cfg['source_archive'];threshold=json.loads((self.source/'leading_thresholds.json').read_text())
        self.taus={'tau90':threshold['0.9']['tau'],'tau95':threshold['0.95']['tau']};self.start=time.perf_counter()

    def tau(self,t):return self.taus.get(t,t) if isinstance(t,str) else float(t)

    def audit_source(self):
        require(sha(self.source/'manifest.json')==self.cfg['source_manifest_sha256'],'Source manifest changed')
        for name,digest in self.cfg['frozen_files'].items():require(sha(ROOT/name)==digest,'Frozen file '+name)
        for name,digest in json.loads((self.source/'manifest.json').read_text())['artifacts'].items():require(sha(self.source/name)==digest,'Source artifact '+name)
        attachment=self.cfg['HS_LA1'];require(sha(ROOT/attachment['copy_path'])==attachment['sha256'],'HS-LA1 source')
        self.archived={}
        for path in (self.source/'leading').glob('*.json'):
            r=json.loads(path.read_text());t=r['tau'];self.archived[t]=(r,path.with_suffix('.npz'))
        write(self.folder/'source_audit.json',dict(status='PASS',source_archive=self.cfg['source_archive'],HS_LA1=attachment,full_model_calls=0))

    def low(self,p,tau):
        tau=self.tau(tau);k=(p,tau)
        if k in self.solutions:return self.solutions[k]
        if tau in self.archived:
            r,path=self.archived[tau];r=next(a for a in r['orders'] if a['order']==p)
            with np.load(path) as arrays:c=arrays[f'p{p}_coefficients'].copy()
            s=dict(zeta=np.array(r['zeta']),coefficients=c,diagnostics=r['diagnostics'],per_mode=r['per_mode']);origin='HS-L1';self.reused+=1
        else:s=self.models[p].solve(tau);origin='new';self.new+=1
        record=dict(order=p,tau=tau,origin=origin,zeta=s['zeta'].tolist(),diagnostics=s['diagnostics'],per_mode=s['per_mode'])
        policy=self.cfg['policy']
        for metric,tol in [('strong_ode','strong_tol'),('natural_boundary','boundary_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
            if s['diagnostics'][metric]>policy[tol]:
                write(self.folder/'failed_low.json',record);raise ArithmeticError(f'{p} {tau}: {metric}')
        require(np.all(s['zeta']>0) and np.all(np.diff(s['zeta'])>0),'Sorted positive prefix')
        self.records.append(record);self.solutions[k]=s
        name=hashlib.sha256(repr(k).encode()).hexdigest()[:16];write(self.folder/'points'/f'{name}.json',record)
        np.savez_compressed(self.folder/'points'/f'{name}.npz',coefficients=s['coefficients'])
        return s

    def basis(self,p,tau):
        tau=self.tau(tau);key=(p,tau)
        if key not in self.bases:
            b=all_pairs(self.models[p],tau,self.cfg['full_basis_arithmetic']['digits']);d=b['diagnostics'];policy=self.cfg['policy']
            row=dict(order=p,tau=tau,**d);self.basis_records.append(row)
            write(self.folder/'full_basis_checks.json',self.basis_records)
            # Gate in the declared working arithmetic. Binary64 exports of the
            # high-mode coefficient vectors are NOT the accepted eigenpairs.
            for kind in ['mp']:
                require(d[f'mass_gram_{kind}']<policy['mass_tol'],'Full-basis mass '+kind)
                require(d[f'matrix_{kind}']<policy['matrix_tol'],'Full-basis matrix '+kind)
                require(d[f'energy_{kind}']<policy['energy_tol'],'Full-basis energy '+kind)
            low=self.low(p,tau);require(np.max(abs(np.sqrt(b['zeta'][:3]/low['zeta'])-1))<policy['frequency_convergence_tol'],'All-basis/accepted prefix')
            name=hashlib.sha256(repr(key).encode()).hexdigest()[:16]
            write(self.folder/f'full_basis_{name}_80digits.json',b['authoritative'])
            np.savez_compressed(self.folder/f'full_basis_{name}.npz',zeta=b['zeta'],coefficients=b['coefficients'],d_n=b['d_n'],K=b['K'],M=b['M'],d=b['d'])
            self.bases[key]=b
        return self.bases[key]

    def gate_A(self):
        rows=[];initial=[]
        for p,m in self.models.items():
            b=self.basis(p,0);lim=self.low(p,RIGID);zero=self.low(p,0)
            require(np.linalg.matrix_rank(m.P)==1,'Rank-one update')
            require(np.linalg.matrix_rank(np.vstack((m.constraint,m.D)))==2,'Exact saturation constraint')
            for t in self.cfg['controls']:
                tau=self.tau(t)
                if tau==RIGID:continue
                s=self.low(p,tau)
                for i in range(3):
                    vals=[b['zeta'][i],s['zeta'][i],lim['zeta'][i],b['zeta'][i+1]]
                    margin=float(min(np.diff(vals)));require(margin>=-self.cfg['algebra']['interlacing_tol']*max(vals),'Interlacing')
                    rows.append(dict(order=p,tau=tau,position=i+1,zeta0=vals[0],zeta=vals[1],zeta_infty=vals[2],zeta0_next=vals[3],minimum_margin=margin,gap12=float(s['zeta'][1]-s['zeta'][0]),gap23=float(s['zeta'][2]-s['zeta'][1])))
            initial.append(dict(order=p,zeta0=b['zeta'][:4].tolist(),zeta_infty=lim['zeta'].tolist(),omega0=self.cfg['epsilon']*math.sqrt(zero['zeta'][0]),omega_infty=self.cfg['epsilon']*math.sqrt(lim['zeta'][0]),derivative0=float((m.D@zero['coefficients'][:,0])**2)))
        csvwrite(self.folder/'interlacing.csv',rows);write(self.folder/'gate_A.json',dict(status='PASS',initial=initial))

    def gate_B(self):
        rows=[];Grows=[];modes=[]
        for p in self.models:
            b=self.basis(p,0);z0=b['zeta'];dn=b['d_n'];projection=abs(dn)/np.linalg.norm(dn)
            for i,(z,d,sensitivity) in enumerate(zip(z0,dn,projection)):
                modes.append(dict(order=p,position=i+1,zeta0=z,d_n=d,relative_projection=sensitivity,status='POSSIBLY_INSENSITIVE' if sensitivity<self.cfg['algebra']['inactive_projection_tol'] else 'ACTIVE'))
            for z in [-1.,(z0[0]+z0[1])/2,(z0[1]+z0[2])/2,(z0[2]+z0[3])/2]:
                gl,gp=compliance_linear(b['K'],b['M'],b['d'],z);gs,sp=compliance_modal(z0,dn,z)
                error=abs(gl-gs)/max(1,abs(gl),abs(gs));require(error<self.cfg['algebra']['G_normalized_tol'],'Compliance identity')
                Grows.append(dict(order=p,zeta=z,G_linear=gl,G_all_modes=gs,normalized_difference=error))
            for t in self.cfg['secular_controls']:
                tau=self.tau(t);s=self.low(p,tau)
                for i,z in enumerate(s['zeta']):
                    distance=float(min(abs(z0-z))/max(1,abs(z)))
                    require(distance>self.cfg['algebra']['pole_distance_relative_tol'],'Secular root is at a pole: requires subspace audit')
                    # Bracket around the direct position INSIDE its original pole interval.
                    radius=min(abs(z0-z))/2;lo=z-radius;hi=z+radius
                    fun=lambda value:1+tau*compliance_modal(z0,dn,value)[0]
                    require(fun(lo)*fun(hi)<0,'Active secular root not isolated')
                    zs=brentq(fun,lo,hi,xtol=self.cfg['algebra']['root_xtol'],rtol=self.cfg['algebra']['root_rtol'])
                    gl,_=compliance_linear(b['K'],b['M'],b['d'],z)
                    gs,_=compliance_modal(z0,dn,z)
                    error=abs(math.sqrt(zs/z)-1)
                    require(error<self.cfg['algebra']['secular_frequency_tol'],'Secular frequency')
                    require(max(abs(1+tau*gl),abs(1+tau*gs))<self.cfg['algebra']['secular_residual_tol'],'Secular equation residual')
                    rows.append(dict(order=p,tau=tau,position=i+1,zeta_direct=z,zeta_secular=zs,frequency_relative_error=error,equation_linear=1+tau*gl,equation_all_modes=1+tau*gs,pole_relative_distance=distance,status='ACTIVE_REGULAR'))
        csvwrite(self.folder/'secular.csv',rows);csvwrite(self.folder/'compliance.csv',Grows);csvwrite(self.folder/'unperturbed_modes.csv',modes)
        write(self.folder/'gate_B.json',dict(status='PASS',all_pairs_used=True,mode_count_by_order={p:len(self.basis(p,0)['zeta']) for p in self.models}))

    def gate_C(self):
        rows=[];derivatives=[];curvature=[]
        for p,m in self.models.items():
            for t in self.cfg['sensitivity_controls']:
                tau=self.tau(t);s=self.low(p,tau);z=float(s['zeta'][0]);D=abs(float(m.D@s['coefficients'][:,0]));der=D*D;P=tau*der/z
                gap=float(s['zeta'][1]-z);require(gap/max(1,z)>self.cfg['algebra']['gap_relative_tol'],'Simple first pair not resolved')
                rows.append(dict(order=p,tau=tau,zeta=z,D_abs=D,D_squared=der,derivative=der,P_D=P,gap12=gap,log_slope_limit_at_zero=0. if tau==0 else None))
                for factor in self.cfg['difference_steps_relative']:
                    h=factor*max(1,tau)
                    if tau==0:fd=(-3*z+4*self.low(p,h)['zeta'][0]-self.low(p,2*h)['zeta'][0])/(2*h);method='forward_second_order'
                    else:fd=(self.low(p,tau+h)['zeta'][0]-self.low(p,tau-h)['zeta'][0])/(2*h);method='central'
                    err=abs(fd/der-1);require(err<self.cfg['algebra']['finite_difference_relative_tol'],'Sensitivity difference')
                    derivatives.append(dict(order=p,tau=tau,method=method,step=h,exact=der,finite_difference=fd,absolute_error=abs(fd-der),relative_error=err))
                if tau>0:
                    for h in self.cfg['log_beta_steps']:
                        zp=self.low(p,tau*math.exp(2*h))['zeta'][0];zm=self.low(p,tau*math.exp(-2*h))['zeta'][0]
                        slope=(math.log(zp)-math.log(zm))/(4*h);err=abs(slope/P-1)
                        require(err<self.cfg['algebra']['finite_difference_relative_tol'],'Log-beta slope')
                        derivatives.append(dict(order=p,tau=tau,method='central_log_beta',step=h,exact=P,finite_difference=slope,absolute_error=abs(slope-P),relative_error=err))
                b=self.basis(p,tau);dn=b['d_n'];zz=b['zeta'];curv=-2*dn[0]**2*np.sum(dn[1:]**2/(zz[1:]-zz[0]))
                for factor in self.cfg['curvature_steps_relative']:
                    h=factor*max(1,tau)
                    if tau==0:fd=(2*z-5*self.low(p,h)['zeta'][0]+4*self.low(p,2*h)['zeta'][0]-self.low(p,3*h)['zeta'][0])/h**2
                    else:fd=(self.low(p,tau+h)['zeta'][0]-2*z+self.low(p,tau-h)['zeta'][0])/h**2
                    error=abs(fd/curv-1)
                    curvature.append(dict(order=p,tau=tau,step=h,curvature_all_modes=curv,curvature_difference=fd,relative_error=error,status='PASS' if error<self.cfg['algebra']['curvature_relative_tol'] and curv<0 else 'UNRESOLVED'))
        csvwrite(self.folder/'sensitivity.csv',rows);csvwrite(self.folder/'finite_differences.csv',derivatives);csvwrite(self.folder/'curvature.csv',curvature)
        write(self.folder/'gate_C.json',dict(status='PASS',curvature_status='PASS' if all(r['status']=='PASS' for r in curvature) else 'UNRESOLVED'))

    def gate_D(self):
        reactions=[];tail=[];estimates=[]
        for p,m in self.models.items():
            lim=self.low(p,RIGID);cinf,r=reaction(m,lim);b=self.basis(p,0);dn=b['d_n']
            distance=float(min(abs(b['zeta']-r['zeta_infty'])))
            require(distance>self.cfg['algebra']['pole_distance_relative_tol']*max(1,r['zeta_infty']),'Saturation coincides with base pole')
            G,Gprime=compliance_linear(b['K'],b['M'],b['d'],r['zeta_infty']);Gs,Gps=compliance_modal(b['zeta'],dn,r['zeta_infty'])
            r.update(order=p,G_infty=G,G_infty_all_modes=Gs,Gprime=Gprime,Gprime_all_modes=Gps,r_squared_from_G=1/Gprime,r_squared_from_modal_G=1/Gps,pole_distance=distance)
            require(r['kkt_residual']<self.cfg['policy']['matrix_tol'],'KKT reaction')
            errors=[abs(r[k]/r['r']-1) for k in ['r_outer','r_joint']]+[abs(r['r_squared_from_G']/r['r_squared']-1),abs(r['r_squared_from_modal_G']/r['r_squared']-1)]
            require(max(errors)<self.cfg['algebra']['reaction_relative_tol'],'Reaction identity')
            require(max(abs(G),abs(Gs))<self.cfg['algebra']['G_normalized_tol'],'Zero of compliance at infinity')
            reactions.append(r)
            for t in self.cfg['asymptotic_controls']:
                tau=self.tau(t);s=self.low(p,tau);c=s['coefficients'][:,0];overlap=float(cinf@m.M@c)
                require(abs(overlap)>.5,'Phase alignment unresolved')
                if overlap<0:c=-c
                D=float(m.D@c);z=float(s['zeta'][0]);deficit=r['zeta_infty']-z;fdef=1-math.sqrt(z/r['zeta_infty']);pred=r['r_squared']/(2*r['zeta_infty']*tau)
                tail.append(dict(order=p,tau=tau,zeta=z,D_aligned=D,phase_overlap_abs=abs(overlap),r_infty_abs=r['r'],r_infty_squared=r['r_squared'],tau_zeta_deficit=tau*deficit,zeta_limit_relative_error=tau*deficit/r['r_squared']-1,zeta_deficit_prediction_relative_error=r['r_squared']/(tau*deficit)-1,tau_D=tau*D,reaction_limit_relative_error=tau*D/r['r']-1,frequency_deficit=fdef,frequency_deficit_prediction=pred,frequency_deficit_prediction_relative_error=pred/fdef-1))
            for q,tau in [(.9,self.taus['tau90']),(.95,self.taus['tau95'])]:
                a=r['r_squared']/(r['zeta_infty']*(1-q*q));bq=r['r_squared']/(2*r['zeta_infty']*(1-q))
                estimates.append(dict(order=p,q=q,tau_HS_L1=tau,tau_zeta_prediction=a,tau_omega_prediction=bq,zeta_prediction_relative_error=a/tau-1,omega_prediction_relative_error=bq/tau-1))
        write(self.folder/'reaction.json',reactions);csvwrite(self.folder/'saturation.csv',tail);csvwrite(self.folder/'threshold_asymptotics.csv',estimates)
        write(self.folder/'gate_D.json',dict(status='PASS',coefficient_from_exact_constrained_state=True,fitting=False))

    def run(self):
        if (self.folder/'manifest.json').exists():self.verify();return
        write(self.folder/'contract.json',self.ident);self.audit_source()
        try:
            for gate in [self.gate_A,self.gate_B,self.gate_C,self.gate_D]:gate();print(gate.__name__,'PASS',flush=True)
            convergence=[]
            for t in set(tau for p,tau in self.solutions):
                if any((p,t) not in self.solutions for p in self.cfg['orders']):continue
                a,b=[self.solutions[p,t] for p in self.cfg['orders']]
                error=float(max(abs(np.sqrt(a['zeta']/b['zeta'])-1)));require(error<self.cfg['policy']['frequency_convergence_tol'],'p20/24 convergence')
                convergence.append(dict(tau=t,frequency_relative_error=error))
            write(self.folder/'convergence.json',convergence);write(self.folder/'point_records.json',self.records)
            write(self.folder/'manifest.json',dict(status='PASS',identity=self.ident,new_low_prefixes=self.new,reused_low_prefixes=self.reused,full_finite_basis_solves=len(self.bases),full_model_solves=0,spectrum_seconds=time.perf_counter()-self.start,artifacts={p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file()}))
        except Exception as e:
            write(self.folder/'failure.json',dict(status='PARTIAL_PASS',error=repr(e),new_low_prefixes=self.new,reused_low_prefixes=self.reused));raise

    def verify(self):
        self.audit_source();m=json.loads((self.folder/'manifest.json').read_text());require(m['identity']==self.ident,'Identity')
        for name,digest in m['artifacts'].items():require(sha(self.folder/name)==digest,'Artifact '+name)


def plot_only(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    start=time.perf_counter();folder=Path(folder);output=folder/'figures';output.mkdir(exist_ok=True)
    rows=[r for r in csv.DictReader(open(folder/'saturation.csv',encoding='utf8')) if r['order']=='24']
    sensitivity=[r for r in csv.DictReader(open(folder/'sensitivity.csv',encoding='utf8')) if r['order']=='24']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    rows.sort(key=lambda r:float(r['tau']));t=[float(r['tau']) for r in rows]
    axes[0].semilogx(t,[float(r['tau_zeta_deficit']) for r in rows],'o-',label='Расчёт')
    axes[0].axhline(float(rows[0]['r_infty_squared']),ls='--',color='black',label='Квадрат предельной реакции')
    axes[0].set_ylabel(r'$\tau(\zeta_\infty-\zeta_1)$');axes[0].legend(fontsize=9)
    sensitivity.sort(key=lambda r:float(r['tau']))
    axes[1].plot([float(r['tau']) for r in sensitivity],[float(r['P_D']) for r in sensitivity],'o-')
    axes[1].set_ylabel('Доля приведённой угловой энергии')
    for ax in axes:ax.set_xlabel(r'Приведённая жёсткость $\tau$');ax.grid(alpha=.2)
    fig.savefig(output/'01_saturation_and_energy.pdf');fig.savefig(output/'01_saturation_and_energy.png',dpi=320);plt.close(fig)
    for name in ['saturation.csv','sensitivity.csv']:(output/name).write_bytes((folder/name).read_bytes())
    write(folder/'figures_manifest.json',dict(calculation_mode='plot_only',new_eigenvalue_solves=0,render_seconds=time.perf_counter()-start,dpi=320,artifacts={p.name:sha(p) for p in output.iterdir()}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['run','verify','plot_only'],default='run');args=parser.parse_args()
    cfg=json.loads(CONFIG.read_text());identifier,_=identity(cfg);folder=BASE/identifier
    if args.stage=='plot_only':plot_only(folder)
    else:
        task=Audit(cfg)
        if args.stage=='run':task.run()
        else:task.verify()
    print(args.stage,folder)


if __name__=='__main__':main()
