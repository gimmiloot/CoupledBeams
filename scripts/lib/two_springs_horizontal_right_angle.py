"""HS-LA5: three fixed-right-angle thin-beam limits, no full-model solves.

Reuses frozen integrated fields. The mass/trace contract differs from HS-L1:
F is clamped at zero, V is pinned in displacement at both ends, and the
mass includes L2*F(L1)**2. Natural forces are not Ritz constraints.
This stable diagnostic workflow preserves historical source hashes.
"""
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import cholesky, eigh, expm, null_space, solve, solve_triangular
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading
from scripts.lib.two_springs_horizontal_spectral_structure import ROOT, sha, write, csvwrite, require

CONFIG = ROOT/'data/input/isotropic_two_springs_horizontal_right_angle_limit.json'
BASE = ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_right_angle_limit'


def symbolic_checks():
    """Exact integer polynomial identities; no optional CAS dependency."""
    # Monomials are sorted tuples of formal variable names.
    one={():1}
    def add(*terms):
        result={}
        for term in terms:
            for monomial,value in term.items():result[monomial]=result.get(monomial,0)+value
        return {k:v for k,v in result.items() if v}
    def mul(x,y):
        result={}
        for a,b in x.items():
            for c,d in y.items():
                key=tuple(sorted(a+c));result[key]=result.get(key,0)+b*d
        return {k:v for k,v in result.items() if v}
    def scale(x,a):return {k:a*v for k,v in x.items() if a*v}
    def sub(x,y):return add(x,scale(y,-1))
    ch,sh,co,si,t,c,A,B,V,Y=({(name,):1} for name in ['ch','sh','co','si','t','c','A','B','V','Y'])
    # t=(L2/L1)*z. Derive the determinant before applying the two identities.
    determinant=sub(mul(add(ch,co),add(ch,co,mul(t,sub(sh,si)))),
                    mul(add(sh,si),add(sub(sh,si),mul(t,sub(ch,co)))))
    characteristic=add(one,mul(ch,co),mul(t,sub(mul(sh,co),mul(ch,si))))
    identities=add(sub(sub(mul(ch,ch),mul(sh,sh)),one),sub(add(mul(co,co),mul(si,si)),one))
    determinant_error=sub(sub(determinant,scale(characteristic,2)),identities)
    mass_error=sub(add(mul(sub(V,mul(c,A)),sub(Y,mul(c,B))),mul(sub(one,mul(c,c)),mul(A,B))),
                   sub(add(mul(V,Y),mul(A,B)),mul(c,add(mul(A,Y),mul(B,V)))))
    require(not determinant_error and not mass_error, 'Analytical identity failed')
    return dict(status='PASS',tip_mass_determinant_identity=0,mass_bilinear_identity=0,
        method='exact integer polynomial expansion modulo ch^2-sh^2=1, co^2+si^2=1')


class RightAngle:
    def __init__(self, k_joint, mu=.01, order=24, quadrature=80):
        self.base = Leading(.1, mu, order, quadrature)
        b = self.base
        self.k_joint = k_joint
        self.A = b.ends[0][0]
        self.H = b.lengths[1]*(b.weights@b.values[1][0])
        self.C = np.vstack((b.outer[1][0], b.ends[1][0]))
        self.Z = null_space(self.C)
        self.K = b.Kb+k_joint*np.outer(b.jump,b.jump)
        self.M = b.M+b.lengths[1]*np.outer(self.A,self.A)
        unit = np.zeros(len(self.A)); unit[-1] = 1.
        self.T = np.eye(len(unit))-np.outer(unit,self.A)  # (F,V) -> (F,W2=V-A)
        self.Tinv = np.eye(len(unit))+np.outer(unit,self.A)
        self.Minner = self.T.T@b.M@self.T
        self.cross = -np.outer(self.A,self.H)-np.outer(self.H,self.A)

    def identities(self):
        b = self.base; Z = self.Z
        oldC = np.vstack((b.constraint,b.D))
        oldZ = null_space(oldC)
        transformed = self.T@Z
        # Equality of subspaces, without imposing any natural condition.
        Q,_ = np.linalg.qr(transformed)
        errors = dict(stiffness=float(np.max(abs(self.T.T@self.K@self.T-self.K))),
            mass=float(np.max(abs(self.Minner-self.M-self.cross))),
            inverse=float(np.max(abs(self.T@self.Tinv-np.eye(len(self.A))))),
            constraints=float(np.max(abs(oldC@transformed))),
            space=float(np.max(abs(Q@Q.T-oldZ@oldZ.T))))
        ev = np.linalg.eigvalsh(Z.T@self.cross@Z)
        return dict(errors=errors,rank_constraints=int(np.linalg.matrix_rank(self.C)),
            dimension=Z.shape[1],mass_difference_extreme_eigenvalues=[float(ev[0]),float(ev[-1])])

    def solve(self, count=3):
        Z = self.Z; K = Z.T@self.K@Z; M = Z.T@self.M@Z
        mass_chol = cholesky((M+M.T)/2,lower=True)
        L = cholesky((K+K.T)/2,lower=True)
        B = solve_triangular(L,M,lower=True)
        B = solve_triangular(L,B.T,lower=True).T
        v,c = eigh((B+B.T)/2,subset_by_index=[len(B)-count,len(B)-1])
        v=v[::-1];c=c[:,::-1];z=1/v
        coefficients=Z@solve_triangular(L.T,c,lower=False)/np.sqrt(v)
        diagnostics,per=self.diagnose(z,coefficients)
        diagnostics.update(mass_positive=True,stiffness_positive=True,zero_nullity=0,
            mass_cholesky_min=float(np.min(np.diag(mass_chol))),
            stiffness_cholesky_min=float(np.min(np.diag(L))))
        return dict(zeta=z,coefficients=coefficients,diagnostics=diagnostics,per_mode=per)

    def diagnose(self,z,coefficients):
        b=self.base;Z=self.Z;per=[]
        for i,ev in enumerate(z):
            ci=coefficients[:,i]
            e=np.array([x@ci for x in b.ends]);o=np.array([x@ci for x in b.outer])
            J=float(b.jump@ci);A=float(self.A@ci)
            natural=np.array([o[1,2],e[0,2]+self.k_joint*J,e[1,2]-self.k_joint*J,
                              e[0,3]+ev*b.lengths[1]*A])
            scale=max(1.,np.max(abs(e)),np.max(abs(o)),abs(ev*b.lengths[1]*A),abs(self.k_joint*J))
            ni=di=0.;bending=[];mass=[]
            for arm,L in enumerate(b.lengths):
                values=[v@ci for v in b.values[arm]]
                ni+=float(b.weights@(values[4]-ev*values[0])**2)*L
                di+=float(b.weights@(values[4]**2+(ev*values[0])**2))*L
                bending.append(float(b.weights@values[2]**2)*L)
                mass.append(float(b.weights@values[0]**2)*L)
            res=Z.T@(self.K@ci-ev*self.M@ci)
            den=np.linalg.norm(Z.T@self.K@ci)+ev*np.linalg.norm(Z.T@self.M@ci)
            stiffness=float(ci@self.K@ci);tip=b.lengths[1]*A*A
            per.append(dict(natural=float(np.max(abs(natural))/scale),natural_rows=natural.tolist(),
                strong=math.sqrt(ni/di),energy=abs(stiffness/ev-1),
                matrix_residual=float(np.linalg.norm(res)/den),A=A,J=J,
                mass_F=mass[0],mass_V=mass[1],mass_rigid_translation=tip,
                mass_cross_inner=float(ci@self.cross@ci),mass_inner_on_outer=float(ci@self.Minner@ci),
                bending_F=bending[0],bending_V=bending[1],joint_energy=self.k_joint*J*J))
        diagnostics=dict(mass_gram=float(np.max(abs(coefficients.T@self.M@coefficients-np.eye(len(z))))),
            essential=float(np.max(abs(self.C@coefficients))),
            symmetry=float(max(np.max(abs(self.K-self.K.T)),np.max(abs(self.M-self.M.T)))))
        for key in ['natural','strong','energy','matrix_residual']:diagnostics[key]=max(x[key] for x in per)
        return diagnostics,per


def boundary_matrix(zeta,k_joint,mu=.01):
    """Eight independent ODE conditions, no Ritz matrices."""
    A=np.zeros((4,4));A[0,1]=A[1,2]=A[2,3]=1;A[3,0]=zeta
    o=[np.eye(8)[:4],np.eye(8)[4:]]
    e=[expm(A*L)@x for L,x in zip((1-mu,1+mu),o)]
    J=e[0][1]-e[1][1]
    rows=[o[0][0],o[0][1],o[1][0],e[1][0],o[1][2],
          e[0][2]+k_joint*J,e[1][2]-k_joint*J,e[0][3]+zeta*(1+mu)*e[0][0]]
    B=np.array(rows)
    return B/np.linalg.norm(B,axis=1)[:,None]


def endpoint_count(zeta,k_joint,mu,policy):
    """Exact EB endpoint-inertia count, including the endpoint mass.

    Schur inertia: two clamped-clamped arm counts plus endpoint inertia.
    The end mass contributes -zeta*L2*A*A^T, never a static spring.
    """
    require(zeta>0,'Positive non-root count probe required')
    K=np.zeros((8,8));fixed=0;nearest=math.inf
    for arm,L in enumerate((1-mu,1+mu)):
        alpha=zeta**.25*L;n=1
        while (n+.5)*math.pi-.1<alpha+.2:
            pole=brentq(lambda x:math.cos(x)-1/math.cosh(x),(n+.5)*math.pi-.1,(n+.5)*math.pi+.1)
            nearest=min(nearest,abs(alpha-pole)/max(1,alpha));fixed+=int(pole<alpha);n+=1
        require(nearest>policy['pole_relative_guard'],'Too close to an auxiliary pole')
        A=np.zeros((4,4));A[0,1]=A[1,2]=A[2,3]=1;A[3,0]=zeta;T=expm(A*L)
        X=np.zeros((4,4));X[:2,:2]=np.eye(2)
        X[2:]=solve(T[:2,2:],np.column_stack((-T[:2,:2],np.eye(2))))
        Y=T@X;beam=np.vstack((X[3],-X[2],-Y[3],Y[2]))
        K[4*arm:4*arm+4,4*arm:4*arm+4]=(beam+beam.T)/2
    unit=np.eye(8);J=unit[3]-unit[7];A=unit[2]
    K+=k_joint*np.outer(J,J)-zeta*(1+mu)*np.outer(A,A)
    Z=null_space(np.array([unit[0],unit[1],unit[4],unit[6]]))
    ev=np.linalg.eigvalsh(Z.T@K@Z);margin=float(min(abs(ev))/max(1,max(abs(ev))))
    require(margin>policy['inertia_relative_zero_tol'],'Indeterminate endpoint inertia')
    return dict(count=int(fixed+sum(ev<0)),fixed_fixed_count=fixed,
        negative_inertia=int(sum(ev<0)),inertia_margin=margin,
        nearest_pole_distance=None if math.isinf(nearest) else nearest)


def independent_checks(z,k,mu,policy):
    require(np.linalg.matrix_rank(boundary_matrix(0,k,mu))==8,'Unexpected zero mechanism')
    near_zero=endpoint_count(1e-5,k,mu,policy)
    require(near_zero['count']==0,'Missing low positive mode')
    checks=[]
    for index,value in enumerate(z):
        lo=value*(1-policy['count_bracket_relative']);hi=value*(1+policy['count_bracket_relative'])
        lower=endpoint_count(lo,k,mu,policy);upper=endpoint_count(hi,k,mu,policy)
        require(lower['count']==index and upper['count']==index+1,'Low-prefix count mismatch')
        root=brentq(lambda x:np.linalg.det(boundary_matrix(x,k,mu)),lo,hi,xtol=1e-13,rtol=2e-14)
        error=abs(math.sqrt(root/value)-1)
        require(error<policy['independent_frequency_tol'],'Independent ODE frequency failed')
        singular=np.linalg.svd(boundary_matrix(root,k,mu),compute_uv=False)
        checks.append(dict(position=index+1,zeta_ode=float(root),frequency_error=error,
            lower=lower,upper=upper,singular_ratio=float(singular[-1]/singular[0])))
    return dict(zero_nullity=0,count_near_zero=near_zero,roots=checks)


def identity(cfg):
    record=dict(config=cfg,source_sha256=sha(__file__),config_sha256=sha(CONFIG),
        initial_sha256=sha(ROOT/cfg['initial_workspace']),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def source_audit(cfg):
    base=ROOT/'results/isotropic_two_springs_rl_timoshenko'
    p=base/cfg['source_archives']['full'];inner=base/cfg['source_archives']['inner']
    contract=json.loads((p/'contract.json').read_text())
    # HS-P1 schema embeds the material contract in its identity record.
    input_cfg=json.loads((ROOT/'data/input/isotropic_two_springs_horizontal_pilot.json').read_text())
    require(input_cfg==contract['config'],'Archived and current HS-P1 contracts differ')
    material=input_cfg['material']
    require(input_cfg['mu']==cfg['mu'] and material['epsilon']==cfg['epsilon'],'Source geometry')
    require(material['l']==material['E']==material['rho']==1 and input_cfg['axial_variant']=='RL_FULL','Source units/model')
    require(material['aspect']==4 and material['nu']==.3 and material['ks']==5/6,'Source material/section')
    with open(p/'spectra.csv',encoding='utf8') as f:
        rows={r['id']:r for r in csv.DictReader(f) if float(r['beta_deg'])==90}
    selected={key:rows[key] for key in [*cfg['full_ids'].values(),cfg['secondary_archived_exact_support']]}
    for key,r in selected.items():
        point=json.loads((ROOT/r['point_file']).read_text())
        require(r['status']==point['status']=='PASS' and point['count']['status']=='PASS','Full source not accepted')
        require(point['case']['beta']==90 and point['case']['mu']==cfg['mu'],'Source angle/geometry')
        require(float(r['omega1'])==point['omega'][0] and float(r['kappa_j'])==point['case']['k_joint'],'Source identity')
        require(abs(float(r['omega1'])/(cfg['epsilon']*float(r['Lambda1'])**2)-1)<1e-13,'Frequency units')
    records=json.loads((inner/'point_records.json').read_text())
    inner_manifest=json.loads((inner/'manifest.json').read_text())
    require(inner_manifest['status']=='PASS' and inner_manifest['identity']['config']['mu']==cfg['mu']
            and inner_manifest['identity']['config']['epsilon']==cfg['epsilon'],'HS-LA4 source contract')
    return selected,records,dict(status='PASS',full_rows=selected,
        source_contract_sha256=sha(p/'contract.json'),inner_manifest_sha256=sha(inner/'manifest.json'))


def run():
    cfg=json.loads(CONFIG.read_text());identifier,record=identity(cfg);folder=BASE/identifier
    state=json.loads((ROOT/cfg['initial_workspace']).read_text())
    for path,digest in cfg['frozen_files'].items():require(sha(ROOT/path)==digest,'Changed frozen source '+path)
    for path,digest in state['archive_files'].items():require(sha(ROOT/path)==digest,'Changed old archive '+path)
    if (folder/'manifest.json').exists():
        manifest=json.loads((folder/'manifest.json').read_text());require(manifest['identity']==record,'Identity changed')
        for path,digest in manifest['artifacts'].items():require(sha(folder/path)==digest,'Changed artifact '+path)
        print('HS-LA5 verified',folder);return folder
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',record)
    solutions={};models={};points=[];new=0
    try:
        write(folder/'symbolic_checks.json',symbolic_checks())
        sources,old_records,audit=source_audit(cfg);write(folder/'source_audit.json',audit)
        inner_folder=ROOT/'results/isotropic_two_springs_rl_timoshenko'/cfg['source_archives']['inner']
        comparisons=[];energies=[];forms=[];convergence=[]
        for k in cfg['k_joint']:
            for p in cfg['orders']:
                m=RightAngle(k,cfg['mu'],p,cfg['quadrature']);models[k,p]=m
                identities=m.identities()
                require(max(identities['errors'].values())<cfg['policy']['form_identity_tol'],'Form/space identity')
                s=m.solve(cfg['spectral_positions']);new+=1
                d=s['diagnostics'];pol=cfg['policy']
                for name,tol in [('natural','boundary_tol'),('strong','strong_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
                    require(d[name]<pol[tol],f'{k},{p}: {name} failed')
                require(s['zeta'][0]>0 and np.all(np.diff(s['zeta'])>0),'Invalid sorted prefix')
                old=next(r for r in old_records if r['k_joint']==k and r['order']==p and r['tau']=='RIGID')
                key=hashlib.sha256(repr((k,p,'RIGID')).encode()).hexdigest()[:16]
                with np.load(inner_folder/'points'/f'{key}.npz') as data:oldc=data['coefficients'].copy()
                transformed=m.Tinv@oldc
                residual=m.Z.T@(m.K@transformed-(m.Minner@transformed)*np.array(old['zeta']))
                scale=np.linalg.norm(m.Z.T@m.K@transformed,axis=0)+np.array(old['zeta'])*np.linalg.norm(m.Z.T@m.Minner@transformed,axis=0)
                archived_error=float(max(np.linalg.norm(residual,axis=0)/scale))
                require(archived_error<pol['matrix_tol'],'Archived inner transformed equation')
                independent=independent_checks(s['zeta'],k,cfg['mu'],pol) if p==24 else None
                row=dict(k_joint=k,order=p,zeta=s['zeta'].tolist(),diagnostics=d,per_mode=s['per_mode'],
                    form_identities=identities,transformed_inner_residual=archived_error,independent=independent)
                points.append(row);solutions[k,p]=s
                name=f'k{k}_p{p}'
                write(folder/'points'/f'{name}.json',row)
                np.savez_compressed(folder/'points'/f'{name}.npz',coefficients=s['coefficients'],K=m.K,Mouter=m.M,Minner=m.Minner,
                    Z=m.Z,T=m.T,cross=m.cross,inner_coefficients_transformed=transformed)
                inner_c=transformed[:,0];outer_c=s['coefficients'][:,0]
                energy=dict(k_joint=k,order=p,**s['per_mode'][0],inner_zeta=old['zeta'][0],
                    inner_norm=float(inner_c@m.Minner@inner_c),outer_mass_on_inner=float(inner_c@m.M@inner_c),
                    cross_on_inner=float(inner_c@m.cross@inner_c),
                    rayleigh_outer_on_inner=float(inner_c@m.K@inner_c/(inner_c@m.M@inner_c)),
                    rayleigh_inner_on_outer=float(outer_c@m.K@outer_c/(outer_c@m.Minner@outer_c)))
                require(s['zeta'][0]<=energy['rayleigh_outer_on_inner']*(1+1e-10),'Outer Rayleigh bound')
                require(old['zeta'][0]<=energy['rayleigh_inner_on_outer']*(1+1e-10),'Inner Rayleigh bound')
                energies.append(energy)
                if p==24:
                    full=sources[cfg['full_ids'][str(k)]];omega=cfg['epsilon']*math.sqrt(s['zeta'][0]);fullomega=float(full['omega1'])
                    inneromega=cfg['epsilon']*math.sqrt(old['zeta'][0])
                    comparisons.append(dict(k_joint=k,source_id=full['id'],omega_inner=inneromega,omega_outer=omega,omega_full=fullomega,
                        inner_error_percent=100*(inneromega/fullomega-1),outer_error_percent=100*(omega/fullomega-1),
                        mass_change_percent=100*(omega-inneromega)/fullomega,
                        residual_to_original_gap=abs(omega-fullomega)/abs(inneromega-fullomega)))
                    for arm,L in enumerate(m.base.lengths):
                        values=[v@outer_c for v in m.base.values[arm][:4]]
                        for i,x in enumerate(m.base.nodes):forms.append(dict(k_joint=k,arm=arm+1,s=x*L,W=values[0][i],W1=values[1][i],W2=values[2][i],W3=values[3][i]))
                print('HS-LA5',k,p,'PASS',flush=True)
            aa,bb=[solutions[k,p] for p in cfg['orders']]
            err=float(max(abs(np.sqrt(aa['zeta']/bb['zeta'])-1)))
            m20,m24=[models[k,p] for p in cfg['orders']]
            fields=[]
            for m,s in [(m20,aa),(m24,bb)]:
                v=np.vstack([m.base.field_matrices(i,m24.base.nodes)[0]@s['coefficients'] for i in [0,1]])
                tip=(m.A@s['coefficients'])*math.sqrt(1+cfg['mu']);fields.append(np.vstack((v,tip)))
            weights=np.r_[m24.base.weights*m24.base.lengths[0],m24.base.weights*m24.base.lengths[1],1.]
            sign=np.sign(np.sum(weights[:,None]*fields[0]*fields[1],axis=0))
            shape=float(max(np.sqrt(np.sum(weights[:,None]*(fields[0]-fields[1]*sign)**2,axis=0))))
            require(err<pol['frequency_convergence_tol'] and shape<pol['shape_convergence_tol'],'Order convergence failed')
            convergence.append(dict(k_joint=k,frequency_error=err,shape_mass_error=shape))
        ratio=(1+cfg['mu'])/(1-cfg['mu'])
        fun=lambda z:1+math.cosh(z)*math.cos(z)+ratio*z*(math.sinh(z)*math.cos(z)-math.cosh(z)*math.sin(z))
        z=brentq(fun,.1,math.pi/2,xtol=1e-14);zeta=(z/(1-cfg['mu']))**4
        omega=cfg['epsilon']*math.sqrt(zeta)
        benchmark=dict(z=z,q=z/(1-cfg['mu']),zeta=zeta,omega=omega,characteristic_residual=fun(z),
            V_first_zeta=(math.pi/(1+cfg['mu']))**4,
            ritz_frequency_error=abs(math.sqrt(solutions[0,24]['zeta'][0]/zeta)-1))
        require(abs(omega/cfg['benchmark']['omega_outer_90']-1)<pol['analytic_frequency_tol'] and benchmark['ritz_frequency_error']<pol['analytic_frequency_tol'],'Tip mass benchmark')
        require(zeta<benchmark['V_first_zeta'],'Wrong first family at hinge')
        write(folder/'benchmark.json',benchmark);write(folder/'point_records.json',points)
        csvwrite(folder/'comparison.csv',comparisons);csvwrite(folder/'energies.csv',energies)
        csvwrite(folder/'forms.csv',forms);csvwrite(folder/'convergence.csv',convergence)
        r=next(c for c in comparisons if c['k_joint']==.1)
        exact=float(sources['H6']['omega1']);write(folder/'exact_support_context.json',dict(omega_full_H5=r['omega_full'],omega_full_H6=exact,
            support_difference_percent=100*(exact/r['omega_full']-1),outer_error_H6_percent=100*(r['omega_outer']/exact-1),new_full_solves=0))
        significant=max(abs(r['outer_error_percent']) for r in comparisons)>cfg['interpretation']['significant_residual_percent']
        status='PARTIAL_PASS' if significant else 'PASS'
        write(folder/'manifest.json',dict(status=status,identity=record,new_outer_ritz_prefixes=new,
            physical_cases=3,reused_inner_prefixes=6,reused_full_points=4,independent_ode_positive_roots=9,
            significant_residual=significant,new_full_model_solves=0,stop=True,
            artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()}))
        print('HS-LA5',status,folder);return folder
    except Exception as e:
        write(folder/'failure.json',dict(status='FAIL',error=repr(e),new_outer_ritz_prefixes=new,stop=True));raise

