"""Two source-specific straight-beam benchmarks, not production joint models.

Failla (2014), (1),(2),(7)-(16); Hong & Kim (1999), (1)-(3),(13).
Mappings/printed targets were transcribed before root solving; see the report.
Four exact segments avoid a long product of hyperbolic transfer matrices.
"""
from dataclasses import dataclass, asdict
from decimal import Decimal
import numpy as np
from scipy.linalg import expm, expm_frechet
from scripts.lib import inplane_kelvin_voigt as kv

FAILLA_KEY = 'failla_2014_viscoelastic_discontinuous_beams'
HONG_KEY = 'hong_kim_1999_damped_timoshenko_joints'
FAILLA = [('21.9037','0.2629','0.0120'), ('39.8987','1.5487','0.0388'),
          ('81.8244','6.5410','0.0796'), ('157.9140',None,'0.0'),
          ('238.023','14.4789','0.0607')]
HONG_HH = ['3.55778498e+002','1.41882448e+003','3.17650875e+003',
           '5.60854995e+003','8.68806614e+003']
HONG_FF = ['8.05532939e+002','2.21142738e+003','4.30966485e+003',
           '7.06898477e+003','1.04602176e+004']
HONG_DAMPED = [('-6.6651e-002','3.3444e+002'),('-2.7327e+000','1.1079e+003'),
               ('-1.2133e+001','1.9271e+003'),('-2.0106e+001','2.9542e+003'),
               ('-2.0135e+001','4.7111e+003')]
ROUNDING_RESERVE = 1e-10
ZERO_ATOL = 1e-8
F_TO_PROJECT = np.array([[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0.]])
H_TO_PROJECT = np.diag([1.,-1.,-1.,-1.])


def rounding(printed, computed):
    """None denotes an explicitly asserted but unprinted zero component."""
    if printed is None:
        target, half, digits, places = 0., 0., None, None
        reserve = ZERO_ATOL
    else:
        dec = Decimal(printed)
        target = float(dec)
        half = float(Decimal('0.5')*Decimal(10)**dec.as_tuple().exponent)
        digits = len(dec.as_tuple().digits)
        places = -dec.as_tuple().exponent
        reserve = ROUNDING_RESERVE*max(1.,abs(target))
    error = abs(computed-target)
    return dict(printed_value=printed,parsed_value=target,computed_value=float(computed),
        printed_significant_digits=digits,decimal_places=places,
        abs_error=float(error),rel_error=float(error/abs(target)) if target else None,
        rounding_tolerance=half,solver_reserve=reserve,
        half_last_digit_units=float(error/half) if half else None,
        rounding_pass=bool(np.isfinite(computed) and error<=half+reserve))


def reference_time(L,m,D):
    return L**2*np.sqrt(m/D)


def failla_to_project(omega_F):
    return 1j*complex(omega_F)


def failla_state(z, derivative=False):
    H=np.zeros((4,4),complex)
    if not derivative:
        H[0,1],H[1,2],H[2,3]=1.,-1.,1.
    H[3,0]=2*z if derivative else z*z
    return H


def failla_interface(z, ku=100.,gu=.1,kr=10.,gr=.1,derivative=False):
    """L*y_minus + R*y_plus=0, y=[v,theta,mu,T]; kr=None means no RJ."""
    L,R=np.zeros((4,4),complex),np.zeros((4,4),complex)
    if derivative:
        L[3,0]=-gu
        if kr is not None:
            L[1,1],R[1,1]=-gr,gr
    else:
        L[0,0],R[0,0]=-1,1
        L[2,2],R[2,2]=-1,1
        L[3,0],L[3,3],R[3,3]=-(ku+gu*z),-1,1
        if kr is None:
            L[1,1],R[1,1]=-1,1
        else:
            L[1,1],L[1,2],R[1,1]=-(kr+gr*z),1,kr+gr*z
    return L,R


def failla_jump(z,ku=100.,gu=.1,kr=10.,gr=.1):
    """Auxiliary divided jump for equivalence tests, not used at K_r=0."""
    jump=np.eye(4,dtype=complex)
    jump[3,0]=ku+gu*z
    if kr is not None:
        if kr+gr*z==0: raise ValueError('divided jump undefined at K_r=0')
        jump[1,2]=-1/(kr+gr*z)
    return jump


@dataclass(frozen=True)
class Hong:
    L: float=1.
    b: float=.025
    h: float=.025
    E: float=200e9
    G: float=80e9
    nu: float=.3
    rho: float=8000.
    kt: float=2e6
    ct: float=20.

    @property
    def A(self): return self.b*self.h
    @property
    def I(self): return self.b*self.h**3/12
    @property
    def K(self): return 10*(1+self.nu)/(12+11*self.nu)
    @property
    def D(self): return self.E*self.I
    @property
    def S(self): return self.K*self.A*self.G
    @property
    def m(self): return self.rho*self.A
    @property
    def J(self): return self.rho*self.I
    @property
    def time(self): return reference_time(self.L,self.m,self.D)
    def data(self):
        return dict(asdict(self),**{k:getattr(self,k) for k in ('A','I','K','D','S','m','J','time')})


def hong_state(s,properties=Hong(),derivative=False):
    """Literal source (3), physical [u,phi,F,M]; never calls project H."""
    a=properties
    H=np.zeros((4,4),complex)
    if not derivative:
        H[0,1],H[0,2],H[1,3],H[3,2]=1,-1/a.S,1/a.D,1
    factor=2*s if derivative else s*s
    H[2,0],H[3,1]=-a.m*factor,a.J*factor
    return H


class Beam:
    """Fixed 16-unknown exact matching assembly for the two specified papers."""
    def __init__(self,case,calls=None):
        if case not in ('failla','hong_hh','hong_ff','hong_damped'):
            raise ValueError('unknown literature case')
        self.case=case
        self.hong=Hong()
        self.calls=calls if calls is not None else kv.Calls()
        self.time=1. if case=='failla' else self.hong.time
        self.scale=np.array([1.,10.,100.,1000.] if case=='failla' else [1.,10.,1000.,100.])

    def state(self,z,derivative=False):
        if self.case=='failla': H=failla_state(z,derivative)
        else:
            a=self.hong
            units=np.array([a.L,1.,a.D/a.L**2,a.D/a.L])
            H=hong_state(z/a.time,a,derivative)*a.L*units[None,:]/units[:,None]
            if derivative: H/=a.time
        return H*self.scale[None,:]/self.scale[:,None]

    def boundary(self,z,derivative=False):
        left,right=np.zeros((2,4),complex),np.zeros((2,4),complex)
        if self.case=='failla' or self.case=='hong_hh':
            indices=[0,2] if self.case=='failla' else [0,3]
            if not derivative:
                left[np.arange(2),indices]=right[np.arange(2),indices]=1
            units=self.scale[indices]
        else:
            if not derivative:
                left[0,2]=right[0,2]=left[1,3]=right[1,3]=1
            if self.case=='hong_damped':
                a=self.hong
                K=(a.ct/a.time if derivative else a.kt+a.ct*z/a.time)*a.L**3/a.D
                left[0,0],right[0,0]=K,-K
            units=self.scale[[2,3]]
        return left*self.scale/units[:,None],right*self.scale/units[:,None]

    def interfaces(self,z,derivative=False):
        if self.case=='failla':
            L,R=failla_interface(z,derivative=derivative)
            units=np.array([self.scale[0],self.scale[2],self.scale[2],self.scale[3]])
            return L*self.scale/units[:,None],R*self.scale/units[:,None]
        zero=np.zeros((4,4),complex)
        return (zero,zero) if derivative else (-np.eye(4),np.eye(4))

    def transfer(self,z,derivative=False,length=.25):
        H=self.state(z)*length
        if derivative:
            self.calls.frechet+=1
            return expm_frechet(H,self.state(z,True)*length)
        self.calls.expm+=1
        return expm(H),None

    def matrices(self,z,*,derivative=False):
        self.calls.matrix(derivative)
        T,Tz=self.transfer(z,derivative)
        L,R=self.interfaces(z)
        left,right=self.boundary(z)
        B=np.zeros((16,16),complex);Bz=np.zeros_like(B)
        B[:2,:4]=left; B[-2:,-4:]=right@T
        if derivative:
            Lz,Rz=self.interfaces(z,True)
            lz,rz=self.boundary(z,True)
            Bz[:2,:4]=lz;Bz[-2:,-4:]=rz@T+right@Tz
        for j in range(3):
            rows=slice(2+4*j,6+4*j)
            B[rows,4*j:4*j+4]=L@T;B[rows,4*j+4:4*j+8]=R
            if derivative:
                Bz[rows,4*j:4*j+4]=Lz@T+L@Tz
                Bz[rows,4*j+4:4*j+8]=Rz
        return B,Bz if derivative else None

    def recover(self,z,a):
        self.calls.recoveries+=1;self.calls.shape_expm+=1
        step=expm(self.state(z)/256)
        starts=np.asarray(a).reshape(4,4)
        values=np.empty((4,65,4),complex)
        for i in range(4):
            v=starts[i].copy()
            for j in range(65):
                values[i,j]=self.scale*v
                if j<64:v=step@v
        factor=np.max(abs(values[:,:,0]))
        return values/factor,np.asarray(a)/factor

    def physical_conditions(self,z,y):
        """Independent scalar conditions, in fixed declared source units."""
        v,theta,m,t=range(4)
        if self.case=='failla':
            out=[y[0,0,v],y[0,0,m]/100]
            for j in range(3):
                a,b=y[j,-1],y[j+1,0]
                out.extend([b[v]-a[v],((10+.1*z)*(b[theta]-a[theta])+a[m])/100,
                    (b[m]-a[m])/100,(b[t]-a[t]-(100+.1*z)*a[v])/1000])
            out.extend([y[-1,-1,v],y[-1,-1,m]/100])
        else:
            a,b=y[0,0],y[-1,-1]
            if self.case=='hong_hh':out=[a[0],a[3]/100]
            else:
                K=0. if self.case=='hong_ff' else (self.hong.kt+self.hong.ct*z/self.time)*self.hong.L**3/self.hong.D
                out=[(a[2]+K*a[0])/1000,a[3]/100]
            for j in range(3):out.extend((y[j+1,0]-y[j,-1])/self.scale)
            if self.case=='hong_hh':out.extend([b[0],b[3]/100])
            else:out.extend([(b[2]-K*b[0])/1000,b[3]/100])
        return np.asarray(out)/np.max(abs(y/self.scale))


def solve(beam,guess):
    B,_=beam.matrices(guess)
    result=kv.correct(beam.matrices,guess,kv.right_null(B))
    beam.calls.corrections+=result['steps']
    z=result['z']; y,a=beam.recover(z,result['a'])
    B,Bz=beam.matrices(z,derivative=True)
    U,s,Vh=np.linalg.svd(B)
    residual=float(np.linalg.norm(B@a)/(np.linalg.norm(B)*np.linalg.norm(a)))
    physical=beam.physical_conditions(z,y)
    failures=[]
    if result['status']!='CONVERGED':failures.append(result['status'])
    if residual>kv.CRITERIA['null_residual']:failures.append('NULL_RESIDUAL')
    if s[-1]/s[0]>kv.CRITERIA['sigma_ratio']:failures.append('SIGMA_RATIO')
    if np.max(abs(physical))>kv.CRITERIA['physical_residual']:failures.append('PHYSICAL_GATE')
    conjugate,_=beam.matrices(z.conjugate())
    conjres=float(np.linalg.norm(conjugate@a.conj())/(np.linalg.norm(conjugate)*np.linalg.norm(a)))
    if conjres>kv.CRITERIA['null_residual']:failures.append('CONJUGATE_GATE')
    details=dict(status='NUMERICAL_UNRESOLVED' if failures else 'CONVERGED',failures=failures,
        z_real=z.real,z_imag=z.imag,p_real=z.real/beam.time,p_imag=z.imag/beam.time,
        steps=result['steps'],history=result['history'],last_delta_z=result['last_delta_z'],
        solver_residual=residual,sigma_ratio=float(s[-1]/s[0]),
        physical_residual=float(np.max(abs(physical))),physical_conditions=physical,
        conjugate_residual=conjres,left_Bz_right=float(abs(np.vdot(U[:,-1],Bz@Vh.conj().T[:,-1]))))
    if beam.case=='failla':
        details['interfaces']=[dict(v_abs=float(abs(y[j,-1,0])),moment_abs=float(abs(y[j,-1,2])),
            rotation_jump_abs=float(abs(y[j+1,0,1]-y[j,-1,1])),
            shear_jump_abs=float(abs(y[j+1,0,3]-y[j,-1,3]))) for j in range(3)]
        details['damping_ratio']=-z.real/abs(z)
    return z,details,y


# D13 second-pass audit: the complete preceding assembly is preserved.
EQUATION_ATOL = 1e-11  # rad/s, declared before the D13 comparisons/Table3 solve.
FIRST_ASSEMBLY_SHA = '3e4d58cac8c4505c74e97e5189cc316283d36db42d340c1a233fb4db44179aed'


def printed_zeta(p_string, q_string, ratio_string):
    """Ordinary decimal rounding of the displayed components only."""
    from decimal import ROUND_HALF_UP
    p, q = float(p_string), float(q_string or 0)
    implied = q/np.hypot(p, q)
    quantum = Decimal(1).scaleb(Decimal(ratio_string).as_tuple().exponent)
    rounded = Decimal(str(implied)).quantize(quantum, rounding=ROUND_HALF_UP)
    return dict(printed_p=p_string, printed_q=q_string,
                printed_ratio=ratio_string, implied_ratio=float(implied),
                ordinary_rounded_ratio=str(rounded),
                displayed_consistent=rounded == Decimal(ratio_string))


def source_precision_audit(failla_rows, hong_rows, first_diagnostics):
    """Read-only arithmetic on K13 rows, no eigenproblem/refiner calls."""
    if len(failla_rows) != 5 or len(hong_rows) != 10:
        raise ValueError('D13 requires all five Failla and ten Hong saved rows')
    zeta_rows=[]
    for n, row in enumerate(failla_rows, 1):
        strings=(row['printed_value_real'], row['printed_value_imag'] or None,
                 row['printed_value_ratio'])
        if int(row['mode']) != n or strings != FAILLA[n-1]:
            raise ValueError('Failla transcription/order changed')
        zeta_rows.append(dict(mode=n, **printed_zeta(*strings)))
    bare=float((4*np.pi)**2)
    computed=float(failla_rows[3]['computed_real'])
    printed=float(failla_rows[3]['printed_value_real'])
    mode4=dict(bare_double=bare, computed=computed, printed=printed,
               computed_minus_bare=computed-bare, printed_minus_bare=printed-bare,
               ordinary_rounded_bare=format(bare,'.4f'),
               inactive_confirmed=failla_rows[3]['inactive_confirmed']=='True')
    old=first_diagnostics['local_formula_check']['checks']
    analytical={(r['case'],int(r['mode'])):r for r in old}
    required={(c,n) for c in ('hong_hh','hong_ff') for n in range(1,6)}
    if len(old)!=10 or set(analytical)!=required:
        raise ValueError('missing/duplicate independent Table2 formula values')
    comparisons=[]; seen=set()
    for row in hong_rows:
        key=row['case'],int(row['mode']); seen.add(key)
        a=analytical[key]; matrix=float(row['computed_imag'])
        if matrix!=a['matrix_omega']:
            raise ValueError('CSV and independent formula check disagree')
        target=(HONG_HH if key[0]=='hong_hh' else HONG_FF)[key[1]-1]
        if row['printed_value_imag']!=target:
            raise ValueError('Table2 target changed')
        analytic=float(a['formula_omega']); printed=float(target)
        comparisons.append(dict(case=key[0],mode=key[1],printed_string=target,
            printed=printed,matrix=matrix,analytic=analytic,
            matrix_minus_analytic=matrix-analytic,matrix_minus_printed=matrix-printed,
            analytic_minus_printed=analytic-printed,equation_atol=EQUATION_ATOL,
            equation_pass=bool(abs(matrix-analytic)<=EQUATION_ATOL),
            original_printed_status=row['status'],
            original_rounding_pass=row['rounding_pass']=='True'))
    if seen!=required:raise ValueError('missing/duplicate Table2 CSV rows')
    return dict(failla_mode4=mode4,failla_displayed_zeta=zeta_rows,
        hong_table2=comparisons,
        hong_printed_status='FAIL' if any(not x['original_rounding_pass'] for x in comparisons) else 'PASS',
        hong_equation_status='PASS_EQUATION_LEVEL' if all(x['equation_pass'] for x in comparisons) else 'FAIL_EQUATION_LEVEL',
        equation_max_error=max(abs(x['matrix_minus_analytic']) for x in comparisons),
        provenance='K13 CSV plus saved independent formulas; no new Table1/Table2 roots')


def d13_gate(equation_status, *, mapping_unchanged, parameters_unchanged,
             criteria_unchanged):
    # Printed Table2 FAIL is intentionally not a prerequisite in D13.
    return (equation_status=='PASS_EQUATION_LEVEL' and mapping_unchanged
            and parameters_unchanged and criteria_unchanged)


def hong_second_pass_comparison(s, details, mode):
    """s=p directly; printed acceptance is never replaced by solver success."""
    re,im=HONG_DAMPED[mode-1]
    parts={'real':rounding(re,s.real),'imag':rounding(im,s.imag)}
    row=dict(mode=mode,source_key=HONG_KEY,source_table='Table 3',
        computed_real=float(s.real),computed_imag=float(s.imag),
        printed_rounding_status='PRINT_MATCH' if all(x['rounding_pass'] for x in parts.values()) else 'PRINT_MISMATCH',
        equation_solver_status='SOLVER_PASS' if details['status']=='CONVERGED' else details['status'],
        last_printed_place_only=all(x['abs_error']<=2*x['rounding_tolerance']+x['solver_reserve'] for x in parts.values()),
        solver_residual=details['solver_residual'],sigma_ratio=details['sigma_ratio'],
        null_vector_residual=details['solver_residual'],physical_residual=details['physical_residual'],
        conjugate_residual=details['conjugate_residual'],newton_iterations=details['steps'],
        last_delta_z=details['last_delta_z'],left_Bz_right=details['left_Bz_right'])
    for component,part in parts.items():
        for key,value in part.items():row[key+'_'+component]=value
    return row
