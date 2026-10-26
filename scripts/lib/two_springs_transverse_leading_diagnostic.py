"""Single-case arithmetic diagnostic in the EXISTING integrated-Legendre space.

No full-model calls, no frequency map. Rational basis operations below evaluate
the existing polynomial basis and bilinear forms; they introduce no new fields
or mechanical terms. Original solvers and immutable sources remain untouched.
"""
from functools import lru_cache
from fractions import Fraction as F
import json
from pathlib import Path
import hashlib
import numpy as np
import mpmath as mp
from scipy.linalg import cholesky, solve_triangular
from scripts.lib.isotropic_two_springs_external_transverse_spring import Leading
from scripts.lib.isotropic_two_springs_leading import fields
from scripts.lib.two_springs_length_asymmetry import api, ROOT

CONFIG=ROOT/'data/input/two_springs_transverse_leading_diagnostic.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/external_transverse_leading_diagnostic'

def integ(a):
    """Integral in t of a shifted Legendre series, zero at t=0 (rational)."""
    b=[F(0)]*(len(a)+1)
    for n,v in enumerate(a):
        if n==0:b[1]+=v/2
        else:
            b[n+1]+=v/(2*(2*n+1));b[n-1]-=v/(2*(2*n+1))
    b[0]-=sum(v*(-1)**n for n,v in enumerate(b))
    return tuple(b)

def deriv(a):
    """Analytic t derivative of shifted Legendre series; no differencing."""
    b=[a[0]*0]*max(1,len(a)-1)
    for n,v in enumerate(a):
        for j in range(n-1,-1,-2):b[j]+=2*(2*j+1)*v
    return tuple(b)

@lru_cache(None)
def basis(p,arm):
    rows=[(F(1,2),F(1,2))] if arm==1 else []
    for n in range(p-int(arm==1)):
        a=[F(0)]*(n+1);a[n]=F(1);rows.append(integ(integ(a)))
    return tuple(rows)

def m(x):
    if isinstance(x,F):return mp.mpf(x.numerator)/x.denominator
    return mp.mpf(x)

def dot_integral(a,b):
    return sum((a[i]*b[i]/(2*i+1) for i in range(min(len(a),len(b)))),a[0]*0)

def combine(rows,c):
    a=[mp.mpf(0)]*max(map(len,rows))
    for row,v in zip(rows,c):
        for j,b in enumerate(row):a[j]+=m(b)*m(v)
    return tuple(a)

def series_value(a,t):
    x=2*m(t)-1;p0=mp.mpf(1);out=a[0]*p0
    if len(a)==1:return out
    p1=x;out+=a[1]*p1
    for n in range(1,len(a)-1):
        p0,p1=p1,((2*n+1)*x*p1-n*p0)/(n+1);out+=a[n+1]*p1
    return out

@lru_cache(None)
def power_legendre(n):
    """Independent rational recurrence for P_n(2t-1) in monomials."""
    if n==0:return (F(1),)
    if n==1:return (F(-1),F(2))
    a=power_legendre(n-1);b=power_legendre(n-2);r=[F(0)]*(n+1)
    for j,v in enumerate(a):r[j]-=F(2*n-1,n)*v;r[j+1]+=F(2*(2*n-1),n)*v
    for j,v in enumerate(b):r[j]-=F(n-1,n)*v
    return tuple(r)

def power_series(a):
    b=[mp.mpf(0)]*len(a)
    for n,v in enumerate(a):
        for j,r in enumerate(power_legendre(n)):b[j]+=v*m(r)
    return b

def power_derivative(a,k):
    for _ in range(k):a=[i*a[i] for i in range(1,len(a))] or [mp.mpf(0)]
    return a

def horner(a,t):
    out=mp.mpf(0)
    for v in reversed(a):out=out*t+v
    return out

def strong_exact(model,c,zeta,dps=80,compare=False):
    """Exact polynomial L2 integration; binary input coefficients held fixed."""
    with mp.workdps(dps):
        arms=[];num=mp.mpf(0);den=mp.mpf(0);derivative_difference=mp.mpf(0)
        for arm,L0 in enumerate(model.lengths):
            L=m(L0);a=combine(basis(model.order,arm),c[arm*model.order:(arm+1)*model.order])
            if arm==1:a=list(a);a[0]+=m(c[-1]);a=tuple(a)
            b=a
            for _ in range(4):b=deriv(b)
            b=tuple(v/L**4 for v in b);r=[-m(zeta)*v for v in a]
            for i,v in enumerate(b):r[i]+=v
            ni=L*dot_integral(r,r);di=L*(dot_integral(b,b)+m(zeta)**2*dot_integral(a,a))
            arms.append(dict(arm=arm+1,numerator_squared=float(ni),denominator_squared=float(di),strong=float(mp.sqrt(ni/di))))
            num+=ni;den+=di
            if compare:
                power=power_series(a)
                for k in (0,1,2,3,4):
                    d=a
                    for _ in range(k):d=deriv(d)
                    pd=power_derivative(power,k)
                    for t in (mp.mpf(0),mp.mpf('.125'),mp.mpf('.5'),mp.mpf('.875'),mp.mpf(1)):
                        derivative_difference=max(derivative_difference,abs(series_value(d,t)-horner(pd,t))/L**k)
        for arm in arms:arm['numerator_fraction']=arm['numerator_squared']/float(num)
        return dict(strong=float(mp.sqrt(num/den)),numerator_squared=float(num),denominator_squared=float(den),
                    arms=arms,independent_derivatives_max_abs=float(derivative_difference),dps=dps)

def metrics(model,c,zeta,xi):
    c=np.asarray(c,float);zeta=np.asarray(zeta,float);K=model.K0+xi**2/2*model.P;Z=model.Z
    q=Z.T@(K@c-(model.M@c)*zeta)
    scale=np.linalg.norm(Z.T@K@c,axis=0)+np.linalg.norm(Z.T@model.M@c,axis=0)*zeta
    strong=[];natural=[];essential=[];bc_rows=[]
    for j,z in enumerate(zeta):
        v=c[:,j];end=np.array([[e@v for e in arm] for arm in model.ends]);outer=np.array([[e@v for e in arm] for arm in model.outer])
        rows=[outer[1,3]+model.k_w*outer[1,0],-outer[1,2]+model.k_end*outer[1,1],
              end[0,2]+model.k_joint*(model.jump@v),end[1,2]-model.k_joint*(model.jump@v),end[0,3]-end[1,3]+xi**2/2*(model.y@v)]
        norm=max(1.,np.max(abs(end)),np.max(abs(outer)),abs(xi**2/2*(model.y@v)),abs(model.k_w*outer[1,0]))
        natural.append(max(map(abs,rows))/norm);bc_rows.append(dict(raw=rows,scale=norm))
        essential.append(max(abs(outer[0,0]),abs(outer[0,1]),abs(end[0,0]+end[1,0])))
        strong.append(strong_exact(model,v,z)['strong'])
    Kr=Z.T@K@Z;Mr=Z.T@model.M@Z;chol=cholesky(Kr,lower=True)
    a=solve_triangular(chol,Mr,lower=True);a=solve_triangular(chol,a.T,lower=True).T
    return dict(eigen_residual=float(max(np.linalg.norm(q,axis=0)/scale)),energy=float(np.max(abs(np.diag(c.T@K@c)/zeta-1))),
        mass_gram=float(np.max(abs(c.T@model.M@c-np.eye(len(zeta))))),essential=max(essential),natural_boundary=max(natural),
        strong_ode=max(strong),strong_per_position=strong,boundary_rows=bc_rows,
        condition_K=float(np.linalg.cond(Kr)),condition_M=float(np.linalg.cond(Mr)),condition_whitened=float(np.linalg.cond(a)),
        gap12_relative=float(np.sqrt(zeta[1]/zeta[0])-1),symmetry=float(max(np.max(abs(v-v.T)) for v in (K,model.M))),
        mass_positive=bool(np.min(np.linalg.eigvalsh(Mr))>0))

def exact_matrices(model,xi):
    """High precision evaluation of the SAME basis integrals and trace terms.

    Shifted Legendre orthogonality integrates polynomial products exactly.
    Parameter floats are promoted exactly, not fitted or changed.
    """
    p=model.order;n=2*p+1;M=mp.zeros(n);K=mp.zeros(n);ends=[];outer=[]
    for arm,L0 in enumerate(model.lengths):
        L=m(L0);rows=list(basis(p,arm));indices=list(range(arm*p,(arm+1)*p))
        if arm==1:rows.append((F(1),));indices.append(n-1)
        for a,(ia,ra) in enumerate(zip(indices,rows)):
            da=deriv(deriv(ra))
            for ib,rb in zip(indices,rows):
                M[ia,ib]=L*m(dot_integral(ra,rb));K[ia,ib]=m(dot_integral(da,deriv(deriv(rb))))/L**3
        ea=[mp.zeros(n,1) for _ in range(4)];oa=[mp.zeros(n,1) for _ in range(4)]
        for idx,row in zip(indices,rows):
            for k in range(4):
                ea[k][idx]=sum(map(m,row))/L**k;oa[k][idx]=sum(m(v)*(-1)**j for j,v in enumerate(row))/L**k;row=deriv(row)
        ends.append(ea);outer.append(oa)
    jump=ends[0][1]-ends[1][1];y=ends[1][0];rotation=outer[1][1];translation=outer[1][0]
    for coefficient,v in [(m(model.k_joint),jump),(m(model.k_end),rotation),(m(model.k_w),translation),(m(xi)**2/2,y)]:K+=coefficient*(v*v.T)
    constraint=ends[0][0]+ends[1][0]
    # Exact elimination of the already-last constant coefficient; no renumbering.
    Z=mp.zeros(n,n-1)
    for j in range(n-1):Z[j,j]=1;Z[n-1,j]=-constraint[j]/constraint[n-1]
    return K,M,Z

def refine(model,solution,xi,dps=80,exact=True):
    """Bordered Newton refinement of existing Ritz pairs, not a new spectrum search."""
    with mp.workdps(dps):
        if exact:K,M,Z=exact_matrices(model,xi)
        else:
            K=mp.matrix((model.K0+xi**2/2*model.P).tolist());M=mp.matrix(model.M.tolist());Z=mp.matrix(model.Z.tolist())
        Kr=Z.T*K*Z;Mr=Z.T*M*Z;out=[];zs=[];records=[]
        for mode,z0 in enumerate(solution['zeta']):
            if exact:v=mp.matrix([m(x) for x in np.array(solution['coefficients'])[:-1,mode]])
            else:v=mp.matrix((model.Z.T@np.array(solution['coefficients'])[:,mode]).tolist())
            v/=mp.sqrt((v.T*Mr*v)[0]);z=m(z0);history=[]
            for iteration in range(4):
                r=Kr*v-z*(Mr*v);norm=(v.T*Mr*v)[0]-1;size=len(v);B=mp.zeros(size+1);rhs=mp.zeros(size+1,1)
                for i in range(size):
                    rhs[i]=-r[i];B[i,size]=-(Mr*v)[i];B[size,i]=2*(Mr*v)[i]
                    for j in range(size):B[i,j]=Kr[i,j]-z*Mr[i,j]
                rhs[size]=-norm;delta=mp.lu_solve(B,rhs)
                for i in range(size):v[i]+=delta[i]
                z+=delta[size];history.append(float(mp.norm(r)))
            full=Z*v;zs.append(z);out.append(list(full));records.append(history)
        coefficients=np.array([[float(out[j][i]) for j in range(len(out))] for i in range(len(out[0]))])
        zeta=np.array([float(z) for z in zs])
        return dict(zeta=zeta,coefficients=coefficients,iterations=records,dps=dps,exact_integrals=exact,
                    zeta_decimal=[mp.nstr(z,dps) for z in zs],coefficients_decimal=[[mp.nstr(out[j][i],dps) for j in range(len(out))] for i in range(len(out[0]))],
                    diagnostics=metrics(model,coefficients,zeta,xi))

def polynomial_checks():
    """Independent monomials, derivative scaling, and the added constant."""
    with mp.workdps(80):
        errors=[]
        for degree in (0,1,2,3,4,7,12,25):
            # t^n=((1+x)/2)^n: recover the shifted Legendre expansion by
            # triangular elimination against rational monomial coefficients.
            power=[F(0)]*(degree+1);power[degree]=F(1);remaining=power.copy();a=[F(0)]*(degree+1)
            for j in range(degree,-1,-1):
                p=power_legendre(j);a[j]=remaining[j]/p[j]
                for k,v in enumerate(p):remaining[k]-=a[j]*v
            for L0 in (.99,1.01):
                L=m(L0);d=tuple(map(m,a))
                for k in range(5):
                    for t in (mp.mpf(0),mp.mpf('.1'),mp.mpf('.5'),mp.mpf(1)):
                        expected=mp.factorial(degree)/mp.factorial(degree-k)*t**(degree-k)/L**k if degree>=k else mp.mpf(0)
                        errors.append(abs(series_value(d,t)/L**k-expected)/max(1,abs(expected)))
                    d=deriv(d)
        field_error=0.
        for arm,L in enumerate((.99,1.01)):
            points=[0.,.125,.5,.875,1.];actual=fields(24,points,L,arm)
            for j,row in enumerate(basis(24,arm)):
                d=tuple(map(m,row))
                for k in range(5):
                    reference=np.array([float(series_value(d,t)/m(L)**k) for t in points])
                    field_error=max(field_error,float(np.max(abs(actual[k][:,j]-reference)/np.maximum(1,abs(reference)))))
                    d=deriv(d)
        return dict(monomial_relative_max=float(max(errors)),existing_fields_relative_max=field_error,
                    constant_derivatives=[float(series_value(tuple(map(m,(F(1),))),0))]+[0.]*4,
                    degrees=[0,1,2,3,4,7,12,25],lengths=[.99,1.01],finite_differences_used=False)

def verify_source(source):
    io=api();manifest=json.loads((source/'manifest.json').read_text());bad=[]
    for name,digest in manifest['artifacts'].items():
        if io.sha(source/name)!=digest:bad.append(name)
    if bad:raise ValueError(f'Source artifact mutation: {bad}')
    return dict(artifacts_checked=len(manifest['artifacts']),manifest_sha256=io.sha(source/'manifest.json'),unchanged=True)

def plot_saved(folder):
    """Only saved rows; no eigensolve. PNG for diagnostic review."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def read(name):return list(csv.DictReader((folder/name).open(encoding='utf-8')))
    rows=read('residual_spatial.csv');orders=read('orders.csv');fixed=read('corrected_orders.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for arm,style in ((1,'-'),(2,'--')):
        r=[v for v in rows if int(v['arm'])==arm]
        axes[0].plot([float(v['t']) for v in r],[float(v['residual'])*1e8 for v in r],style,label=f'Стержень {arm}')
    axes[0].set(xlabel=r'$s_i/\lambda_i$',ylabel=r'$(W_i^{(4)}-\zeta W_i)\cdot10^8$')
    axes[1].semilogy([int(v['order']) for v in orders],[float(v['strong1']) for v in orders],'o-',label='Квадратура')
    axes[1].semilogy([int(v['order']) for v in fixed],[float(v['strong1']) for v in fixed],'s--',label='Аналитические интегралы')
    axes[1].axhline(1e-9,color='.5',ls=':',label='Допуск')
    axes[1].set(xlabel=r'$p$',ylabel=r'$R_{\mathrm{strong},1}$',xticks=[16,20,22,24,26,28])
    for i,a in enumerate(axes):
        a.grid(alpha=.2);a.legend(fontsize=8);a.text(.02,.97,'а' if i==0 else 'б',transform=a.transAxes,va='top')
    fig.savefig(folder/'residual_diagnostic.png',dpi=320);plt.close(fig)

def run():
    io=api();config=json.loads(CONFIG.read_text(encoding='utf-8'));source=ROOT/config['source']
    files=[CONFIG,Path(__file__),ROOT/'scripts/lib/isotropic_two_springs_leading.py',ROOT/'scripts/lib/isotropic_two_springs_external_transverse_spring.py',source/'manifest.json',source/'leading_gate.json',source/config['failed_point']]
    identity=dict(config=config,sources={str(p.relative_to(ROOT)):io.sha(p) for p in files},numpy=np.__version__,mpmath=mp.__version__)
    key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];folder=BASE/key
    if folder.exists():raise FileExistsError(f'Immutable diagnostic already exists: {folder}')
    folder.mkdir(parents=True);io.write(folder/'contract.json',identity)
    integrity_before=verify_source(source)
    failed=json.loads((source/config['failed_point']).read_text());io.write(folder/'original_failed_point.json',failed)
    saved={}
    for path in (source/'leading_points').glob('*.json'):
        v=json.loads(path.read_text());i=v['identity']
        if i['k_w']==1. and i['xi']==0.:saved[i['order']]=(v['solution'],str(path.relative_to(ROOT)))
    rows=[];solutions={};models={};new_solves=0
    for p in config['orders']:
        model=Leading(.1,10.,mu=.01,order=p,quadrature=64,k_w=1.);models[p]=model
        if p in saved:sol,origin=saved[p]
        else:sol=model.solve(0.);origin='new original double precision solve';new_solves+=1
        solutions[p]=sol;d=metrics(model,sol['coefficients'],sol['zeta'],0.)
        row=dict(order=p,zeta1=sol['zeta'][0],frequency_change_from_previous=None if not rows else float(np.sqrt(sol['zeta'][0]/rows[-1]['zeta1'])-1),
                 original_strong=sol['diagnostics']['strong_ode'],strong1=d['strong_per_position'][0],eigen_residual=d['eigen_residual'],energy=d['energy'],
                 strong_max_three=d['strong_ode'],natural_boundary=d['natural_boundary'],condition_K=d['condition_K'],condition_M=d['condition_M'],condition_whitened=d['condition_whitened'],origin=origin)
        rows.append(row);io.write(folder/f'order_{p}.json',dict(solution=sol,diagnostics=d));print(json.dumps(row),flush=True)
    io.csvwrite(folder/'orders.csv',rows)
    model=models[24];sol=solutions[24];c=np.array(sol['coefficients'])[:,0];z=sol['zeta'][0]
    reference=[strong_exact(model,c,z,dps,True) for dps in config['decimal_precision']];io.write(folder/'derivative_reference.json',reference)
    quadrature=[];spatial=[]
    for q in config['residual_quadratures']:
        x,w=np.polynomial.legendre.leggauss(q);t=(x+1)/2;w=w/2;num=den=0.
        for arm,L in enumerate(model.lengths):
            f=fields(24,t,L,arm);ci=c[arm*24:(arm+1)*24];W=f[0]@ci+(c[-1] if arm==1 else 0);D4=f[4]@ci
            num+=L*(w@((D4-z*W)**2));den+=L*(w@(D4**2+(z*W)**2))
        quadrature.append(dict(quadrature=q,numerator_squared=num,denominator_squared=den,strong=np.sqrt(num/den)))
    for arm,L in enumerate(model.lengths):
        t=np.linspace(0,1,1001);f=fields(24,t,L,arm);W=f[0]@c[arm*24:(arm+1)*24]+(c[-1] if arm==1 else 0);D4=f[4]@c[arm*24:(arm+1)*24]
        for k in range(len(t)):spatial.append(dict(arm=arm+1,t=t[k],s=t[k]*L,W=W[k],D4=D4[k],residual=D4[k]-z*W[k]))
    io.csvwrite(folder/'quadrature.csv',quadrature);io.csvwrite(folder/'residual_spatial.csv',spatial)
    # Separates residual evaluation, rounded matrix solve, and matrix integration.
    refined={}
    for exact in (False,True):
        fixed=refine(model,sol,0.,dps=80,exact=exact);io.write(folder/f'refinement_exact_{exact}.json',fixed)
        refined[exact]=fixed;print('refinement',exact,fixed['diagnostics']['strong_ode'],flush=True)
    refined50=refine(model,sol,0.,dps=50,exact=True);io.write(folder/'refinement_exact_50digits.json',refined50)
    corrected_rows=[];corrected_solutions={}
    for p in config['orders']:
        improved=Leading(.1,10.,mu=.01,order=p,quadrature=64,k_w=1.,analytic_curvature=True)
        corrected=improved.solve(0.);new_solves+=1;d=metrics(improved,corrected['coefficients'],corrected['zeta'],0.)
        corrected_solutions[p]=corrected
        io.write(folder/f'corrected_order_{p}.json',dict(solution=corrected,diagnostics=d))
        row=dict(order=p,zeta1=corrected['zeta'][0],frequency_change_from_previous=None if not corrected_rows else float(np.sqrt(corrected['zeta'][0]/corrected_rows[-1]['zeta1'])-1),
                 strong1=d['strong_per_position'][0],strong_max_three=d['strong_ode'],matrix=d['eigen_residual'],natural_boundary=d['natural_boundary'],energy=d['energy'],mass_gram=d['mass_gram'],essential=d['essential'])
        corrected_rows.append(row);print('corrected',row,flush=True)
    io.csvwrite(folder/'corrected_orders.csv',corrected_rows)
    controls=[]
    for p in (20,24):
        improved=Leading(.1,10.,mu=.01,order=p,quadrature=64,k_w=1.,analytic_curvature=True)
        control=improved.solve(config['control_xi']);new_solves+=1;d=metrics(improved,control['coefficients'],control['zeta'],config['control_xi'])
        record=dict(order=p,xi=config['control_xi'],solution=control,diagnostics=d);controls.append(record)
        io.write(folder/f'control_order_{p}.json',record);print('control',p,d['strong_ode'],flush=True)
    refined_control=refine(improved,control,config['control_xi'],dps=80,exact=True)
    io.write(folder/'control_refined.json',refined_control)
    isolation=[]
    for q in (96,128):
        a=Leading(.1,10.,mu=.01,order=24,quadrature=q,k_w=1.);s=a.solve(0.);new_solves+=1
        isolation.append(dict(quadrature=q,solution=s,strong_reference=strong_exact(a,s['coefficients'][:,0],s['zeta'][0])))
    io.write(folder/'assembly_quadrature.json',isolation)
    repeated=model.solve(0.);new_solves+=1
    io.write(folder/'default_regression.json',dict(zeta_max_relative=float(np.max(abs(np.array(repeated['zeta'])/sol['zeta']-1))),
        coefficient_max_abs=float(np.max(abs(np.array(repeated['coefficients'])-sol['coefficients']))),original_diagnostics=repeated['diagnostics']))
    corrected=corrected_solutions[24];improved=Leading(.1,10.,mu=.01,order=24,quadrature=64,k_w=1.,analytic_curvature=True)
    with mp.workdps(80):
        Ke,Me,Ze=exact_matrices(model,0.)
        exactK=np.array(Ke.tolist(),float);exactM=np.array(Me.tolist(),float)
        matrix_comparison=dict(K_max_abs=float(np.max(abs(model.K0-exactK))),M_max_abs=float(np.max(abs(model.M-exactM))),
            Kb_quadrature_max_abs=float(np.max(abs(model.Kb-improved.Kb))),
            K_after_fix_max_abs=float(np.max(abs(improved.K0-exactK))),
            eigenvalues_50_80_max_relative=float(np.max(abs(refined50['zeta']/refined[True]['zeta']-1))))
    io.write(folder/'matrix_comparison.json',matrix_comparison)
    coefficients=[];a=np.array(sol['coefficients'])[:,0];b=refined[True]['coefficients'][:,0];d=corrected['coefficients'][:,0]
    b*=np.sign(a@b);d=d*np.sign(a@d)
    for i in range(len(a)):coefficients.append(dict(index=i,archived=a[i],refined=b[i],analytic_curvature=d[i],archived_minus_refined=a[i]-b[i]))
    io.csvwrite(folder/'coefficient_errors.csv',coefficients)
    spatial_summary=[]
    for arm in (1,2):
        r=[v for v in spatial if v['arm']==arm];peak=max(r,key=lambda v:abs(v['residual']))
        spatial_summary.append(dict(arm=arm,peak_t=peak['t'],peak_s=peak['s'],peak_residual=peak['residual']))
    io.write(folder/'spatial_summary.json',spatial_summary)
    checks=polynomial_checks();io.write(folder/'polynomial_checks.json',checks)
    tol=json.loads((source/'contract.json').read_text())['config']['leading']['quality']
    accepted=[json.loads((folder/f'corrected_order_{p}.json').read_text())['diagnostics'] for p in config['orders']]+[v['diagnostics'] for v in controls]
    failures=[]
    for i,d in enumerate(accepted):
        for name in ('eigen_residual','energy','mass_gram','essential','natural_boundary','strong_ode','symmetry'):
            if d[name]>tol[name+'_tol']:failures.append([i,name,d[name],tol[name+'_tol']])
        if not d['mass_positive'] or d['gap12_relative']<=0:failures.append([i,'mass_or_gap'])
    frequency_control=float(np.max(abs(np.sqrt(np.array(controls[0]['solution']['zeta'])/controls[1]['solution']['zeta'])-1)))
    if frequency_control>tol['frequency_convergence_rtol']:failures.append(['frequency_control',frequency_control])
    if checks['monomial_relative_max']>1e-50 or checks['existing_fields_relative_max']>1e-10:failures.append(['derivative_polynomial_check'])
    status=dict(status='DIAGNOSTIC_PASS' if not failures else 'DIAGNOSTIC_FAIL',failures=failures,
        cause='Quadrature roundoff in curvature Gram matrix amplified in fourth derivatives; same-space analytic curvature integrals remove it.',
        original_strong=failed['solution']['diagnostics']['strong_ode'],independent_saved_strong=reference[-1]['strong'],
        corrected_p24_strong=corrected['diagnostics']['strong_ode'],control_frequency_order_change=frequency_control,
        tolerance=1e-9,quality=tol,full_leading_map_resumed=False)
    io.write(folder/'decision.json',status)
    io.write(folder/'source_integrity.json',dict(before=integrity_before,after=verify_source(source)))
    io.write(folder/'execution.json',dict(new_original_leading_solves=new_solves,reused_original_pairs=len(saved),new_refined_prefixes=4,new_full_model_solves=0,
        note='Each prefix has three positions; all new calculations are one construction at xi=0 or control xi=10. Earlier exploratory archive 5d9db4ddddcfc591 contains 7 original prefixes and 2 refinements, not counted here.'))
    plot_saved(folder)
    io.write(folder/'manifest.json',dict(status=status['status'],artifacts={str(p.relative_to(folder)).replace('\\','/'):io.sha(p) for p in folder.rglob('*') if p.is_file()}))
    print('FOLDER',folder,flush=True)
    return folder
