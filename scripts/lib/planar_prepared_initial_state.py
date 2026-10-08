"""Separate prepared initial case; no ODE, eigendecomposition or physics edits.

Saved full spectral data provide periodic profiles. Frozen Fraction residuals
provide endpoint identities. Numeric Legendre derivatives are independently
evaluated; the caller explicitly admits a common profile before dynamical use.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
import numpy as np
from numpy.polynomial import Polynomial as PowerPolynomial
from numpy.polynomial.legendre import legder, legval
from scripts.lib import weakly_nonlinear_spatial_rod as rod
VERSION = "prepared-periodic-o2-quintic-endpoint-o3-v1"
CASE_NAME = "prepared_axial_O2_with_cubic_endpoint_compatibility"
FIELDS = ("u", "w", "theta", "c")
_SUFFIXES = ("", "_s", "_t", "_ss", "_st", "_tt")
_IDS = (0, 1, 5, 6)


def periodic_parts(model, direct_check=False, tolerance=2e-11):
    """All saved modes; free_modal_amplitudes is positive sum subtracted in x_free."""
    omega,V = np.asarray(model.omega),np.asarray(model.vectors)
    if omega.ndim != 1 or V.shape != (len(omega),len(omega)) or np.any(omega <= 0):
        raise ValueError("All positive saved modes must be retained")
    squared,driving = omega**2,float(model.driving_omega)
    denominator = squared-driving**2
    separation = abs(denominator)/np.maximum(squared+driving**2,np.finfo(float).tiny)
    unresolved = not np.all(np.isfinite(omega)) or not math.isfinite(driving) or bool(np.any(separation <= 64*np.finfo(float).eps))
    checks = {"coordinates_retained":len(omega),"eigendecompositions":0,
              "dynamic_operator_singular_or_unresolved":bool(unresolved),
              "minimum_absolute_detuning":float(np.min(abs(omega-driving))),
              "minimum_relative_squared_separation":float(np.min(separation)),
              "modal_static_condition":float(np.max(squared)/np.min(squared)),
              "modal_dynamic_condition":None if unresolved else float(np.max(abs(denominator))/np.min(abs(denominator))),
              "singularity_policy":"64 eps denominator roundoff bound; no resonance rounding"}
    if unresolved:
        return {"status":"PREPARATION_UNRESOLVED","checks":checks,"reason":"Unresolved dynamic denominator"}
    ms,mh = model.b0/squared,model.b2/denominator
    stat,harm = V@ms,V@mh
    rs,rh = model.K@stat-model.f0,(model.K-driving**2*model.M)@harm-model.f2
    checks["stat_scaled_residual"] = float(np.linalg.norm(rs)/max(np.linalg.norm(model.K@stat)+np.linalg.norm(model.f0),1e-30))
    checks["harm_scaled_residual"] = float(np.linalg.norm(rh)/max(np.linalg.norm(model.K@harm)+driving**2*np.linalg.norm(model.M@harm)+np.linalg.norm(model.f2),1e-30))
    if direct_check:
        ds,dh = np.linalg.solve(model.K,model.f0),np.linalg.solve(model.K-driving**2*model.M,model.f2)
        checks["direct_stat_relative_difference"] = float(np.linalg.norm(ds-stat)/max(np.linalg.norm(ds),1e-30))
        checks["direct_harm_relative_difference"] = float(np.linalg.norm(dh-harm)/max(np.linalg.norm(dh),1e-30))
    passed = max(checks["stat_scaled_residual"],checks["harm_scaled_residual"]) <= tolerance
    return {"status":"PASS" if passed else "PREPARATION_UNRESOLVED","stat":stat,"harm":harm,
            "modal_stat":ms,"modal_harm":mh,"free_modal_amplitudes":ms+mh,"checks":checks}


def periodic_coordinates(parts,driving_omega,times,derivative=0):
    times = np.asarray(times)
    if parts["status"] != "PASS" or times.ndim != 1 or derivative not in (0,1,2):
        raise ValueError("Resolved parts, 1D times and derivative 0..2 required")
    angle = driving_omega*times
    if derivative == 0:return parts["stat"][None,:]+np.cos(angle)[:,None]*parts["harm"][None,:]
    multiplier = -driving_omega*np.sin(angle) if derivative == 1 else -driving_omega**2*np.cos(angle)
    return multiplier[:,None]*parts["harm"][None,:]


def free_coordinates(model,parts,times,derivative=0):
    times = np.asarray(times)
    if parts["status"] != "PASS" or times.ndim != 1 or derivative not in (0,1,2):
        raise ValueError("Resolved parts, 1D times and derivative 0..2 required")
    angle = times[:,None]*model.omega[None,:]
    multiplier = np.cos(angle) if derivative == 0 else -model.omega[None,:]*np.sin(angle) if derivative == 1 else -model.omega[None,:]**2*np.cos(angle)
    return -(multiplier*parts["free_modal_amplitudes"][None,:])@model.vectors.T


def physical_legendre_coefficients(model,coordinates,degree=None):
    coordinates = np.asarray(coordinates)
    degree = model.p if degree is None else int(degree)
    if coordinates.shape[-1] != model.ndof or degree < model.p:raise ValueError("Coordinate/degree mismatch")
    result = np.zeros(coordinates.shape[:-1]+(2,degree+1))
    for i,name in enumerate(("u","c")):
        raw = coordinates[...,model.slices[name]]@model.transforms[i].T
        result[...,i,:model.n] += raw
        result[...,i,2:model.n+2] -= raw
    return result


@dataclass(frozen=True)
class LegendreProfiles:
    coefficients: np.ndarray
    length: float = 1.
    provenance: str = "candidate physical profiles, not continuum truth"

    def __post_init__(self):
        values = np.asarray(self.coefficients,dtype=float).copy()
        if values.ndim != 2 or values.shape[0] != 2 or not np.all(np.isfinite(values)):raise ValueError("Finite (2,p+1) coefficients required")
        if not math.isfinite(self.length) or self.length <= 0:raise ValueError("Positive length required")
        values.setflags(write=False)
        object.__setattr__(self,"coefficients",values)

    @property
    def degree(self):return self.coefficients.shape[1]-1

    def evaluate(self,points,derivative=0):
        points = np.asarray(points,dtype=float)
        if points.ndim != 1 or np.any(points < 0) or np.any(points > self.length) or derivative not in (0,1,2):raise ValueError("Points/derivative invalid")
        values = legder(self.coefficients,m=derivative,axis=1)*(2/self.length)**derivative
        return legval(2*points/self.length-1,values.T).T

    def endpoint_jets(self):
        return {"value":self.evaluate([0.,self.length]),"first":self.evaluate([0.,self.length],1),"second":self.evaluate([0.,self.length],2)}


@lru_cache(maxsize=1)
def _endpoint_proof():return endpoint_trace_audit(rod.derive_polynomials())


def endpoint_trace_audit(model=None):
    """JSON-safe exact A/B force identities; numeric jets remain independent."""
    if model is None:return _endpoint_proof()
    p = model.symbols
    zero = {name+suffix:0 for name in ("v","Phi","psi") for suffix in _SUFFIXES}
    zero.update({name:0 for name in FIELDS})
    zero.update({name+suffix:0 for name in rod.FIELD_ORDER for suffix in ("_t","_st","_tt")})
    expected = {
        "u":p["C"]*p["u_ss"]+p["nu"]*p["C"]*p["c_s"]+(p["C"]-p["S"])*p["theta_s"]*p["w_s"],
        "w":p["S"]*(p["w_ss"]-p["theta_s"])+(p["C"]-p["S"])*p["u_s"]*p["theta_s"],
        "theta":p["Bp"]*p["theta_ss"]+p["S"]*p["w_s"]-(p["C"]-p["S"])*p["u_s"]*p["w_s"],
        "c":p["H"]*p["c_ss"]-p["nu"]*p["C"]*p["u_s"]}
    traces,checks = {},{}
    for path,residuals in (("A",model.residual_a),("B",model.residual_b)):
        for field,index in zip(FIELDS,_IDS):
            actual = -residuals[index].substitute(zero)
            difference = actual-expected[field]
            checks[path+"_"+field] = {"status":"PASS" if not difference else "FAIL","difference_terms":difference.serialize()["terms"]}
            if path == "A":traces[field] = actual
    weights = {name+suffix:weight for name,weight in (("u",2),("c",2),("w",1),("theta",1)) for suffix in _SUFFIXES}
    def weighted(expression,degree):
        return rod.Polynomial({key:value for key,value in expression.terms.items() if sum(weights.get(rod.SYMBOL_ORDER[index],0) for index in key)==degree})
    formal = {d:{name:weighted(value,d) for name,value in traces.items()} for d in (1,2,3)}
    first,second = formal[3]["w"],-formal[3]["theta"]
    correction_checks = {"S_Theta3_s":first==(p["C"]-p["S"])*p["u_s"]*p["theta_s"],
                         "Bp_Theta3_ss":second==(p["C"]-p["S"])*p["u_s"]*p["w_s"],
                         "w_cubic_after_correction":formal[3]["w"]-first==0,
                         "theta_cubic_after_correction":formal[3]["theta"]+second==0}
    if not all(row["status"]=="PASS" for row in checks.values()) or not all(correction_checks.values()):raise ArithmeticError("Endpoint derivation mismatch")
    return {"status":"PASS","checks":checks,"traces":{name:str(value) for name,value in traces.items()},
            "trace_polynomials":{name:value.serialize() for name,value in traces.items()},
            "formal_before_correction":{str(d):{name:str(value) for name,value in row.items()} for d,row in formal.items()},
            "theta3_first_derivative_numerator":str(first),"theta3_second_derivative_numerator":str(second),
            "correction_checks":correction_checks,"mass_symbols":{"u":"m","w":"m","theta":"jp","c":"jp"},
            "qualification":"Conditional BVP identities and independent numeric jets are distinct evidence"}


_HERMITE = ((0,1,0,-6,8,-3),(0,0,0,-4,7,-3),(0,0,Fraction(1,2),Fraction(-3,2),Fraction(3,2),Fraction(-1,2)),(0,0,0,Fraction(1,2),-1,Fraction(1,2)))


def hermite_basis_audit():
    checks = []
    for row,co in enumerate(_HERMITE):
        for endpoint in (0,1):
            for derivative in (0,1,2):
                value = sum(Fraction(c)*Fraction(math.factorial(k),math.factorial(k-derivative))*endpoint**(k-derivative) for k,c in enumerate(co) if k>=derivative)
                checks.append(value==int((row,endpoint,derivative) in ((0,0,1),(1,1,1),(2,0,2),(3,1,2))))
    if not all(checks):raise ArithmeticError("Hermite exact identity failed")
    return {"status":"PASS","identities":24,"remaining_system_determinant":2,"full_six_condition_system_determinant":4}


@dataclass(frozen=True)
class Theta3Correction:
    eta_power_coefficients: np.ndarray
    length: float
    first_targets: np.ndarray
    second_targets: np.ndarray

    @property
    def coefficients(self):return self.eta_power_coefficients

    @property
    def coefficients_eta(self):return self.eta_power_coefficients

    @property
    def degree(self):return 5

    def evaluate(self,points,derivative=0):
        points = np.asarray(points,dtype=float)
        if points.ndim != 1 or np.any(points<0) or np.any(points>self.length) or derivative not in (0,1,2):raise ValueError("Theta3 points/derivative invalid")
        return PowerPolynomial(self.eta_power_coefficients).deriv(derivative)(points/self.length)/self.length**derivative

    def endpoint_jets(self):
        return {"value":self.evaluate([0.,self.length]),"first":self.evaluate([0.,self.length],1),"second":self.evaluate([0.,self.length],2)}

    def parity_metrics(self):
        poly = PowerPolynomial(self.eta_power_coefficients)
        defect = poly+poly(PowerPolynomial([1.,-1.]))
        scale = max(float(np.max(abs(self.eta_power_coefficients))),1e-30)
        return {"reflection":"Theta3(L-s)=-Theta3(s) for first0=first1, second0=-second1","coefficient_absolute_defect":float(np.max(abs(defect.coef))),"coefficient_scaled_defect":float(np.max(abs(defect.coef))/scale)}

    def as_dict(self):
        return {"eta_power_coefficients":self.eta_power_coefficients.tolist(),"length":self.length,"degree":5,
                "first_targets":self.first_targets.tolist(),"second_targets":self.second_targets.tolist(),
                "endpoint_jets":{name:values.tolist() for name,values in self.endpoint_jets().items()},"parity":self.parity_metrics(),
                "rule":"Unique degree<=5 Hermite; zero values; first L and second L^2 scaling"}


def quintic_theta3(U_s_end,Theta_s_end,W_s_end,coefficients,length=1.):
    if not math.isfinite(length) or length<=0:raise ValueError("Positive length required")
    arrays = tuple(np.asarray(x,dtype=float) for x in (U_s_end,Theta_s_end,W_s_end))
    if any(x.shape!=(2,) or not np.all(np.isfinite(x)) for x in arrays):raise ValueError("Two finite endpoint derivatives per input required")
    us,ts,ws = arrays
    first = (coefficients.C-coefficients.S)/coefficients.S*us*ts
    second = (coefficients.C-coefficients.S)/coefficients.Bp*us*ws
    data = np.r_[length*first,length**2*second]
    hermite_basis_audit()
    return Theta3Correction(data@np.asarray(_HERMITE,dtype=float),float(length),first,second)


def formal_endpoint_coefficients(background_jets,profile_jets,correction_jets,coefficients):
    """Numeric force coefficients; do not drop finite-amplitude orders 4/5."""
    p = coefficients
    ws,ts = background_jets['first'].T;wss,tss = background_jets['second'].T
    us,cs = profile_jets['first'].T;uss,css = profile_jets['second'].T
    t3s,t3ss = correction_jets['first'],correction_jets['second'];z = np.zeros(2)
    return {1:{'u':z.copy(),'w':p.S*(wss-ts),'theta':p.Bp*tss+p.S*ws,'c':z.copy()},
            2:{'u':p.C*uss+p.nu*p.C*cs+(p.C-p.S)*ts*ws,'w':z.copy(),'theta':z.copy(),'c':p.H*css-p.nu*p.C*us},
            3:{'u':z.copy(),'w':-p.S*t3s+(p.C-p.S)*us*ts,'theta':p.Bp*t3ss-(p.C-p.S)*us*ws,'c':z.copy()},
            4:{'u':(p.C-p.S)*t3s*ws,'w':z.copy(),'theta':z.copy(),'c':z.copy()},
            5:{'u':z.copy(),'w':(p.C-p.S)*us*t3s,'theta':z.copy(),'c':z.copy()}}


@dataclass(frozen=True)
class PreparedInitialState:
    background: object
    profiles: LegendreProfiles
    correction: Theta3Correction
    admitted: bool = False
    case_name: str = CASE_NAME

    def __post_init__(self):
        if self.background.length!=self.profiles.length or self.profiles.length!=self.correction.length:raise ValueError("Common initial length required")

    @property
    def length(self):return self.profiles.length

    def evaluate(self,points,epsilon,derivative=0,require_admitted=True):
        if require_admitted and not self.admitted:raise RuntimeError("Common profile/endpoint/projection not admitted")
        if not math.isfinite(epsilon):raise ValueError("Finite dimensionless epsilon_a required")
        a,b,t3 = self.profiles.evaluate(points,derivative),self.background.evaluate(points,derivative),self.correction.evaluate(points,derivative)
        return np.column_stack((epsilon**2*a[:,0],epsilon*b[:,0],epsilon*b[:,1]+epsilon**3*t3,epsilon**2*a[:,1]))

    def initial_velocities(self,points):return np.zeros((len(points),4))

    def endpoint_jets(self):
        x = [0.,self.length]
        return {'background':{'value':self.background.evaluate(x),'first':self.background.evaluate(x,1),'second':self.background.evaluate(x,2)},'profiles':self.profiles.endpoint_jets(),'theta3':self.correction.endpoint_jets()}

    def endpoint_audit(self,coefficients,epsilons=(.05,.025,.0125),model=None):
        model = rod.derive_polynomials() if model is None else model
        proof,jets = endpoint_trace_audit(model),self.endpoint_jets()
        formal = formal_endpoint_coefficients(jets['background'],jets['profiles'],jets['theta3'],coefficients)
        numeric = []
        for epsilon in epsilons:
            q,qs,qss = (self.evaluate([0.,self.length],epsilon,d,require_admitted=False) for d in (0,1,2))
            forces = []
            for endpoint in range(2):
                seven = np.zeros((6,7));seven[0,list(_IDS)] = q[endpoint];seven[1,list(_IDS)] = qs[endpoint];seven[3,list(_IDS)] = qss[endpoint]
                result = rod.polynomial_evaluate(rod.FieldJet(*seven),coefficients,model)
                forces.append(-result['residual'][list(_IDS)])
            forces = np.asarray(forces)
            formal_force = np.column_stack([sum(epsilon**d*formal[d][field] for d in formal) for field in FIELDS])
            numeric.append({'epsilon_a':epsilon,'essential_values':q,'initial_velocities':np.zeros_like(q),'force_trace':forces,
                            'acceleration_trace_at_exact_clamp':forces/np.array([coefficients.m,coefficients.m,coefficients.jp,coefficients.jp]),
                            'formal_force_through_degree5':formal_force,'essential_value_roundoff_difference':forces-formal_force})
        return {'symbolic_status':proof['status'],'formal_force_coefficients':formal,'numerical_jets':jets,'finite_amplitude':numeric,
                'qualification':'Conditional BVP identities and numeric jets remain separate; finite amplitudes retained'}