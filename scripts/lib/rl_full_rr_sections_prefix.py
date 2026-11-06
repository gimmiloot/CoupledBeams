"""Individual-section adapter for the frozen RL_FULL/RR counted prefix.

Physical equations/bases/count/Ritz remain in thickness_asymmetry and its base.
The A-B adapter is imported in a private namespace only for its mass-quality
gate; neither historical module is mutated. No general solver framework.
"""
from __future__ import annotations
from dataclasses import asdict
import importlib.util
import math
import sys
import numpy as np
from scripts.lib import rl_full_rr_prefix as frozen
from scripts.lib import isotropic_two_springs_thickness_asymmetry as p
from scripts.analysis import verify_rl_full_rr_prefix as ab

VERSION = "rl-full-rr-individual-sections-v1"
PrefixFailure, digest, counted_prefix = frozen.PrefixFailure, frozen.digest, frozen.counted_prefix

_spec = importlib.util.spec_from_file_location(__name__+"._gates", frozen.__file__)
_gates = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _gates
_spec.loader.exec_module(_gates)
_gates.physics = p
mode_gates = _gates.mode_gates


def arm_contract(arm, length):
    s = arm.section
    return dict(length=length, material=asdict(arm), EA=s.EA, EI=s.EI,
                KGA=s.KGA, rhoA=s.rhoA, rhoI=s.rhoI, J=arm.J)


def point_contract(model, case, k_plot, k_guard, policy):
    p.compatible(model, case)
    if (model.variant != "RL_FULL" or (case.k_joint, case.k_end) != (p.base.RIGID,)*2
        or not isinstance(k_plot, int) or not isinstance(k_guard, int) or not 1 <= k_plot < k_guard):
        raise ValueError("Individual RL_FULL sections, exact RR and positive prefix ranks required")
    return dict(version=VERSION, model=asdict(model), case=asdict(case), tau=model.tau,
        K_plot=k_plot, K_guard=k_guard, policy=policy, spectrum_semantics="sorted_positions",
        state_order=["u","w","psi","N","Q","M"],
        arms=[arm_contract(a,L) for a,L in zip(model.arms,case.lengths(model))],
        time_reference=model.time, reference_section=arm_contract(p.base.Model(
            **{k:v for k,v in asdict(model).items() if k not in ("mu","tau")}), model.l),
        normalization="Omega=omega*t_ref=Lambda**2; common reference section",
        coordinates="outer clamp to joint")


def catalogues(model, case, upper, policy):
    result, attempts = {}, []
    for arm,L in zip(model.arms,case.lengths(model)):
        key = p.pole_key(arm,L)
        if key in result:
            if result[key]["arm_contract"] != arm_contract(arm,L):
                raise PrefixFailure("POLE_KEY_COLLISION")
            continue
        # Same saturated catalogue and finite domain preflight as accepted A-B.
        local, trials = ab.catalogues(arm, type("OneLength", (), {"lengths":lambda self,m:[L]})(), upper, policy)
        cat = local[L]
        cat["catalog_id"] = digest(dict(arm=cat["arm_contract"], policy=policy,
            upper=cat["upper"], basis_source_sha256=ab.sha(p.base.__file__),
            catalog_source_sha256=ab.sha(sys.modules[ab.tim_catalog.__module__].__file__)))
        result[key] = cat
        attempts.extend(dict(t, pole_key=key, catalog_id=cat["catalog_id"]) for t in trials)
    return result, attempts


class PhysicalCounter:
    """Preflight each arm and the actual shifted query before frozen count_at."""
    def __init__(self, model, case, catalogs, policy):
        p.compatible(model,case)
        self.model,self.case,self.catalogs,self.policy = model,case,catalogs,policy
        self.poles,self.limits = [],[]
        for arm,L in zip(model.arms,case.lengths(model)):
            cat=catalogs[p.pole_key(arm,L)]; s=arm.section
            if cat["arm_contract"] != arm_contract(arm,L):
                raise PrefixFailure("CATALOG_ARM_CONTRACT_MISMATCH",length=L)
            rl=math.sqrt(s.EA/arm.J); tim=math.sqrt(s.KGA/s.rhoI)
            self.limits.append(dict(length=L,RL=rl,Tim=tim,catalog=cat["upper"]))
            if not 0 < cat["upper"] < min(rl,tim):
                raise PrefixFailure("CATALOG_OUTSIDE_BASIS_DOMAIN",length=L,upper=cat["upper"])
            d=s.EA-arm.J*cat["upper"]**2
            nmax=math.floor(cat["upper"]*math.sqrt(s.rhoA/d)*L/math.pi)
            self.poles.extend(cat["poles"])
            self.poles.extend(math.sqrt(s.EA*(n*math.pi/L)**2/(s.rhoA+arm.J*(n*math.pi/L)**2))
                              for n in range(1,nmax+1))
        self.domain=min(v for a in self.limits for k,v in a.items() if k!="length")

    def __call__(self, omega):
        effective=float(omega); exclusion=self.policy["pole_exclusion"]*max(1.,effective)
        for pole in sorted(self.poles):
            if abs(effective-pole)<=exclusion: effective=pole+2*exclusion
        if not 0 < omega < self.domain or not 0 < effective < self.domain:
            raise PrefixFailure("COUNT_OUTSIDE_CATALOG",requested=omega,effective=effective,domain=self.domain)
        actual,count,details=p.count_at(self.model,self.case,omega,self.catalogs,self.policy)
        if actual!=effective:
            raise PrefixFailure("POLE_BYPASS_CONTRACT_MISMATCH",predicted=effective,actual=actual)
        return actual,count,details


def fields(model, case, omega, coefficients, order):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    return [dict(x=(nodes+1)*L/2,weights=weights*L/2,
                 y=p.base.basis(arm,L,omega,(nodes+1)*L/2)@c,
                 dy=p.base.basis(arm,L,omega,(nodes+1)*L/2,1)@c,
                 ends=p.base.basis(arm,L,omega,[0.,L])@c)
            for arm,L,c in zip(model.arms,case.lengths(model),coefficients)]


def recovered_modes(model, case, event, order, policy):
    w,mult=event["omega"],event["multiplicity"]
    if mult==1:
        return [p.mode(model,case,w,order)]
    b,_,col=p.balance(p.boundary(model,case,w))
    _,s,vh=np.linalg.svd(b)
    if int(np.sum(s<=s[0]*policy["nullity_tol"]))!=mult:
        raise PrefixFailure("UNRESOLVED_NULLITY",omega=w)
    coeff=(vh[-mult:]*col).reshape(mult,2,6)
    trial=[dict(arms=fields(model,case,w,c,order)) for c in coeff]
    gram=np.array([[p.inner(model,a,b) for b in trial] for a in trial])
    vals,vec=np.linalg.eigh(gram)
    if min(vals)<=0:
        raise PrefixFailure("NONPOSITIVE_SUBSPACE_MASS",eigenvalues=vals.tolist())
    coeff=(((vec/np.sqrt(vals))@vec.T)@coeff.reshape(mult,-1)).reshape(mult,2,6)
    result=[]
    for c in coeff:
        arms=fields(model,case,w,c,order)
        mass,elastic,springs=p.energies(model,case,arms)
        qs=max(np.max(abs(a["y"][:,:3]*[1,1,model.l])) for a in arms)
        fs=max(np.max(abs(a["y"][:,3:]*[1,1,1/model.l])) for a in arms)
        jr=p.joint_matrix(model,case)@np.concatenate([a["ends"][1] for a in arms])
        er=np.array([p.external_matrix(model,case,i)@a["ends"][0] for i,a in enumerate(arms)])
        scale=np.array([1,1,model.l,qs/fs,qs/fs,qs/(fs*model.l)])
        ode=max(np.max(abs((a["dy"]-a["y"]@p.base.state_matrix(arm,w).T)*scale))/
                max(np.max(abs(a["dy"]*scale)),1e-30) for a,arm in zip(arms,model.arms))
        result.append(dict(coefficients=c.tolist(),arms=arms,diagnostics=dict(
            mass=float(mass.sum()),energy_error=float(abs((elastic.sum()+springs.sum())/mass.sum()/w**2-1)),
            singular_ratio=float(s[-mult]/s[0]),nonzero_condition=float(s[0]/s[-mult-1]),
            joint_scaled=(abs(jr)/[qs,qs,fs,fs,qs/model.l,fs*model.l]).tolist(),
            external_scaled=(abs(er)/[qs,qs,qs/model.l]).tolist(),
            joint_absolute=jr.tolist(),external_absolute=er.tolist(),ode_error=float(ode),
            mass_components=mass.tolist(),elastic_components=elastic.tolist(),
            raw_mass_eigenvalues=vals.tolist())))
    return result


def quadrature_check(model, case, roots, modes, policy):
    high,changes=[],[]
    for root,mode in zip(roots,modes):
        arms=fields(model,case,root["omega"],mode["coefficients"],policy["quadrature_check_order"])
        mass,elastic,sp=p.energies(model,case,arms)
        m0,k0,_=p.energies(model,case,mode["arms"])
        changes.append(dict(mass_change=float(abs(mass.sum()-m0.sum())),
            energy_relative_change=float(abs(elastic.sum()/k0.sum()-1)),
            high_energy_error=float(abs((elastic.sum()+sp.sum())/mass.sum()/root["omega"]**2-1))))
        high.append(dict(arms=arms))
    err=float(np.max(abs(np.array([[p.inner(model,a,b) for b in high] for a in high])-np.eye(len(high)))))
    ok=err<=policy["mass_gram_tol"] and all(c["mass_change"]<=policy["mass_gram_tol"] and
        max(c["energy_relative_change"],c["high_energy_error"])<=policy["energy_tol"] for c in changes)
    return dict(status="PASS" if ok else "QUADRATURE_FAILED",modes=changes,gram_max_error=err,
                low_order=policy["quadrature_order"],high_order=policy["quadrature_check_order"])


def global_inner(model, a, b, beta_a, beta_b):
    """Shared material x/L points at fixed mu/tau, global translations, full mass."""
    total=0.
    for aa,bb,arm,fa,fb in zip(a["arms"],b["arms"],model.arms,
                              p.base.geometry.frames(beta_a),p.base.geometry.frames(beta_b)):
        if not np.array_equal(aa["x"],bb["x"]) or not np.array_equal(aa["weights"],bb["weights"]):
            raise ValueError("Different material quadratures")
        va,vb=aa["y"],bb["y"]; s=arm.section
        da,db=va[:,:2]@fa.translation.T,vb[:,:2]@fb.translation.T
        total+=aa["weights"]@(s.rhoA*(da*db).sum(axis=1)+s.rhoI*va[:,2]*vb[:,2]+arm.J*aa["dy"][:,0]*bb["dy"][:,0])
    return float(total)


REFLECTION_STATE=np.array([1.,-1.,-1.,1.,-1.,-1.])


def reflect_mode(mode):
    """Swap inward arms, reflect global vectors; x and derivative direction stay."""
    return dict(arms=[dict(a,**{key:a[key]*REFLECTION_STATE for key in ("y","dy","ends")})
                      for a in reversed(mode["arms"])])


def compare_modes(model, roots, modes, old_roots, old_modes, policy):
    if len(roots)!=len(old_roots) or [r["multiplicity"] for r in roots]!=[r["multiplicity"] for r in old_roots]:
        return dict(status="MULTIPLICITY_MISMATCH")
    freq=max(abs(a["omega"]/b["omega"]-1) for a,b in zip(roots,old_roots))
    groups=[]; start=0
    while start<len(roots):
        mult=roots[start]["multiplicity"]
        a,b=modes[start:start+mult],old_modes[start:start+mult]
        overlap=np.array([[p.inner(model,x,y) for y in b] for x in a])
        singular=np.linalg.svd(overlap,compute_uv=False)
        errors=[]
        # Simple eigenforms are compared individually, including forces/derivatives.
        # An exactly repeated cluster is aligned only within its full eigenspace.
        u,_,vh=np.linalg.svd(overlap); rotation=u@vh
        for i,new in enumerate(a):
            for arm in range(2):
                for field in ("y","dy","ends"):
                    prior=sum(rotation[i,j]*b[j]["arms"][arm][field] for j in range(mult))
                    actual=new["arms"][arm][field]
                    for cols,scale in ((slice(0,3),[1,1,model.l]),(slice(3,6),[1,1,1/model.l])):
                        errors.append(float(np.max(abs((actual[:,cols]-prior[:,cols])*scale))/
                                      max(np.max(abs(prior[:,cols]*scale)),1e-30)))
        groups.append(dict(first_rank=start+1,multiplicity=mult,singular_values=singular.tolist(),
                           max_field_error=max(errors),subspace_error=float(max(abs(singular**2-1)))))
        start+=mult
    passed=freq<=policy["reference_frequency_tol"] and all(max(g["subspace_error"],g["max_field_error"])<=policy["reference_shape_tol"] for g in groups)
    return dict(status="PASS" if passed else "FIELD_REGRESSION_FAILED",relative_frequency_error=float(freq),groups=groups)


def straight_interface(model,case,roots,modes,policy):
    """Separate axial/bending boundary blocks and physical positive-s traces.

    Uses individual existing bases; independent assembly, not an independent
    basis. Ritz provides the separate energy discretization.
    """
    if case.beta!=0: raise ValueError("Straight control only")
    transform=np.array([-1.,-1.,1.,1.,1.,-1.])
    records=[]
    for root,mode in zip(roots,modes):
        w=root["omega"]; ends=[p.base.basis(a,L,w,[0.,L]) for a,L in zip(model.arms,case.lengths(model))]
        blocks=[]
        for ext,components,cols in (([0],[0,3],[0,1]),([1,2],[1,2,4,5],[2,3,4,5])):
            n=len(cols); mat=np.zeros((2*n,2*n)); e=len(ext)
            mat[:e,:n]=ends[0][0][np.ix_(ext,cols)]
            mat[e:2*e,n:]=ends[1][0][np.ix_(ext,cols)]
            mat[2*e:,:n]=ends[0][1][np.ix_(components,cols)]
            mat[2*e:,n:]=-transform[components,None]*ends[1][1][np.ix_(components,cols)]
            singular=np.linalg.svd(p.balance(mat)[0],compute_uv=False)
            blocks.append(float(singular[-1]/singular[0]))
        left=mode["arms"][0]["ends"][1]; right=mode["arms"][1]["ends"][1]*transform
        qs=max(np.max(abs(a["y"][:,:3]*[1,1,model.l])) for a in mode["arms"])
        fs=max(np.max(abs(a["y"][:,3:]*[1,1,1/model.l])) for a in mode["arms"])
        residual=abs(left-right)/[qs,qs,qs/model.l,fs,fs,fs*model.l]
        records.append(dict(rank=root["sorted_position"],axial_singular_ratio=blocks[0],
            bending_singular_ratio=blocks[1],positive_s_trace_jump_scaled=residual.tolist(),
            positive_s_traces=[left.tolist(),right.tolist()]))
    ok=all(min(r["axial_singular_ratio"],r["bending_singular_ratio"])<=policy["singular_tol"]
           and max(r["positive_s_trace_jump_scaled"])<=policy["boundary_tol"] for r in records)
    return dict(status="PASS" if ok else "STRAIGHT_INTERFACE_FAILED",modes=records,
                right_state_transform=transform.tolist(),right_derivative_transform=(-transform).tolist())

