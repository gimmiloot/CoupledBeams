"""NSL-2 bounded elastic verification, not a parameter-map or universal solver.

Script Proliferation Control: NSL-1 has no joint trace space. This workflow
adds covariant angular assembly, reflection and a coupled energy count.
It reuses NSL-1 physics/prefix isolation and RLB geometry/Ritz basis.
"""
from __future__ import annotations
from dataclasses import asdict, replace
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.linalg import block_diag, null_space
from scipy.optimize import brentq, linear_sum_assignment
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as joint
from scripts.lib import reddy_symmetric_laminated_beam as rlb
from scripts.lib.reddy_symmetric_coupled_beams import joint_matrix as rigid_matrix
from scripts.lib.reddy_inplane_geometry import reddy_inplane_geometry
from scripts.analysis.nonsymmetric_laminated_beams import verify_single_beam as previous
from scripts.analysis.reproduce_bishop_literature import sha, write_csv, write_json
from scripts.analysis.check_contraction_joint_sensitivity import plain

CONFIG = ROOT/"data/input/nonsymmetric_laminated_beam_coupled_joint.json"
OUTPUT = ROOT/"results/nonsymmetric_laminated_beams/coupled_joint_verification"
read = previous.read


def relative(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.max(abs(a-b))/max(np.max(abs(a)), np.max(abs(b)), 1e-300))


def materials(config):
    old = read(ROOT/config["material_input"])
    properties, audits = {}, {}
    for name, layout in (("BASE", "HHHH"), ("HHLL", "HHLL"), ("LLHH", "LLHH")):
        local = json.loads(json.dumps(old))
        if name == "BASE":
            local["spectral_case"]["chi"] = 0.
        section, _ = previous.hhll(local, list(layout))
        a, audit = single.project_reddy_coupled_properties(section,
            width=old["spectral_case"]["width"], K=old["spectral_case"]["K"],
            provenance=f"NSL-2 {name}; NSL-1/RLB-2E material input, full joint ABD Schur")
        properties[name], audits[name] = a, audit
    return properties, audits


def symbolic():
    if str(ROOT/"tmp/nsl_tools") not in sys.path:
        sys.path.insert(0, str(ROOT/"tmp/nsl_tools"))
    import sympy as s
    A, B, D, S, m, j1, j2, z, p = s.symbols("A B D S m j1 j2 z p")
    def h(A, B, D, m, j1, j2):
        de = A*D-B*B
        return s.Matrix([[0,0,0,D/de,0,-B/de], [0,0,-1,0,1/S,0],
            [0,0,0,-B/de,0,A/de], [m*p*p,0,j1*p*p,0,0,0],
            [0,m*p*p,0,0,0,0], [j1*p*p,0,j2*p*p,0,1,0]])
    t = s.eye(6); t[0,2]=z; t[5,3]=-z
    original = h(A,B,D,m,j1,j2)
    shifted = h(A,B-z*A,D-2*z*B+z*z*A,m,j1-z*m,j2-2*z*j1+z*z*m)
    assert s.simplify(shifted*t-t*original) == s.zeros(6)
    work = s.zeros(6); work[:3,3:]=s.eye(3); work[3:,:3]=-s.eye(3)
    assert t.T*work*t == work
    cont = s.diag(-1,-1,1,1,1,-1)
    assert s.simplify(h(A,-B,D,m,-j1,j2)+cont*original*cont) == s.zeros(6)
    u,v,ps,n,q,mm,du,dv,dps = s.symbols("u v ps n q mm du dv dps")
    assert s.expand(n*(du+z*dps)+q*dv+(mm-z*n)*dps-(n*du+q*dv+mm*dps)) == 0
    return {"status":"PASS","sympy":s.__version__,"axis_similarity":True,
            "axis_work_duality":True,"straight_coordinate_reversal":True,
            "Lean":"not available; direct derivation and exact SymPy checks"}


def structural(config, properties):
    policy = config["policy"]; tol=policy["operator_tol"]
    base, a, rev = (properties[n] for n in ("BASE","HHLL","LLHH"))
    common = ("A","D","S","m","j2")
    material_error = max(relative([getattr(x,n) for n in common],
        [getattr(base,n) for n in common]) for x in (a,rev))
    assert material_error < policy["coefficient_tol"] and abs(a.B+rev.B) < 1e-18
    assert base.B == 0 and a.j1 == rev.j1 == base.j1 == 0
    z0 = a.B/a.A
    shifted = joint.shift_properties(a,z0)
    assert abs(shifted.B)<1e-18 and abs(shifted.j1)>0
    rng=np.random.default_rng(20261020); rows=[]
    swap=np.block([[np.zeros((6,6)),np.eye(6)],[np.eye(6),np.zeros((6,6))]])
    k=config["reference"]["D"]/config["reference"]["length"]
    for beta_deg in config["structural_angles_deg"]:
        beta=np.deg2rad(beta_deg)
        gs=[joint.global_kinematic_map(beta,i) for i in (1,2)]
        basis=reddy_inplane_geometry(beta_deg)
        geometry_error=max(relative(g.T@g,np.eye(3)) for g in gs)
        # Planar reflection exchanges tangents and reverses the normals.
        r=gs[1][:2,:2]@np.diag([1.,-1.])@gs[0][:2,:2].T
        geometry_error=max(geometry_error,relative(r@r,np.eye(2)),
            relative(r@gs[0][:2,:2],gs[1][:2,:2]@np.diag([1.,-1.])))
        assert np.linalg.det(r)<0 and geometry_error<tol
        for kind, spring in (("HINGE",0.),("SPRING",k),("RIGID",None)):
            law=joint.Joint("RIGID" if kind=="RIGID" else "SPRING",spring)
            j=joint.joint_matrix(beta,law)
            jn=j/np.linalg.norm(j,axis=1)[:,None]
            kernel=null_space(jn); rank=np.linalg.matrix_rank(jn)
            assert rank==6
            trace, jump=joint.trace_maps(beta,law)
            # End states in ker(J) and arbitrary admissible kinematic variations.
            state=kernel@rng.normal(size=6); dz=rng.normal(size=len(jump))
            virtual=sum(state[6*i+3:6*i+6]@(trace[i]@dz) for i in (0,1))
            spring_work=0 if spring is None else spring*(state[2]-state[8])*(jump@dz)
            work_error=abs(virtual+spring_work)/max(np.linalg.norm(state)*np.linalg.norm(dz),1e-30)
            # Full polar/axial physical work, independent of force row signs.
            dv=[rng.normal(size=3),rng.normal(size=3)]
            physical=sum((gs[i]@state[6*i+3:6*i+6])@(gs[i]@dv[i]) for i in (0,1))
            local=sum(state[6*i+3:6*i+6]@dv[i] for i in (0,1))
            work_error=max(work_error,abs(local-physical)/max(abs(local),abs(physical),1))
            # Exchange arms and their actual bases: no tacit beta sign shortcut.
            qswap=np.r_[state[6:9],state[:3]]
            fswap=np.r_[state[9:12],state[3:6]]
            force_swap=(gs[1]@fswap[:3]+gs[0]@fswap[3:])[:2]
            disp_swap=(gs[1]@qswap[:3]-gs[0]@qswap[3:])[:2]
            exchange=max(np.max(abs(force_swap)),np.max(abs(disp_swap)))
            reflection_kernel=float(np.max(abs(jn@joint.REFLECTION@kernel)))
            assert max(work_error,exchange,reflection_kernel)<tol
            for name, pair in config["cases"].items():
                props=[properties[x] for x in pair]
                model=joint.CoupledBeams(props,config["lengths"],beta,law)
                offsets=[x.B/x.A for x in props]
                operator_errors=[]
                for prop,z in zip(props,offsets):
                    t=joint.axis_transform(z)
                    operator_errors.append(relative(single.state_matrix(.3+.7j,joint.shift_properties(prop,z))@t,
                                                    t@single.state_matrix(.3+.7j,prop)))
                boundary_error=max(relative(model.boundary(w,False),
                    model.boundary(w,False,offsets)) for w in (.23,.91))
                h=block_diag(*(single.state_matrix(.3+.7j,x) for x in props))
                commutator=relative(h@joint.REFLECTION,joint.REFLECTION@h)
                # Full boundary intertwining in NSL-1 centered coefficient coordinates.
                b=model.boundary(.47,False)
                tq=single.REVERSE[:3,:3]
                external=np.block([[np.zeros((3,3)),tq],[tq,np.zeros((3,3))]])
                row_joint=(j@joint.REFLECTION)@np.linalg.pinv(j)
                full_reflection=relative(b@joint.REFLECTION,block_diag(external,row_joint)@b)
                if name in ("BASE","MIRROR"):
                    assert max(commutator,full_reflection)<tol
                else:
                    assert commutator>1e-5 and full_reflection>1e-5
                assert max(operator_errors+[boundary_error])<tol
                rows.append(dict(beta_deg=beta_deg,joint=kind,case=name,rank=rank,
                    geometry=geometry_error,virtual_work=work_error,exchange=exchange,
                    reflection_kernel=reflection_kernel,reflection_commutator=commutator,
                    full_boundary_reflection=full_reflection,
                    axis_operator=max(operator_errors),axis_boundary=boundary_error))
    # B=j1=0: exact existing RLB state and rigid/spring physical rows.
    bp=rlb.BeamProperties(A=base.A,D=base.D,S=base.S,m=base.m,J=base.j2,K=base.K,width=base.width)
    symmetric=relative(single.state_matrix(.7j,base).real,rlb.combined_state_matrix(.7,bp))
    assert symmetric<tol
    assert np.array_equal(joint.joint_matrix(.3,joint.Joint("RIGID")),rigid_matrix(.3))
    # A naive neutral-axis node is NOT the same physical joint (negative control).
    law=joint.Joint("SPRING",k)
    wrong=relative(joint.covariant_joint(np.pi/6,law,[z0,-z0]),joint.joint_matrix(np.pi/6,law))
    assert wrong>1e-5
    # Endpoint DtN: trace dimensions, work and Schur inertia, no eigenvalue solve.
    count=[]
    for law in (joint.Joint("SPRING",k),joint.Joint("RIGID")):
        model=joint.CoupledBeams([a,rev],config["lengths"],np.pi/6,law)
        d=model.count_matrix(.47,policy["count_segments"],True)
        assert d["count"]==d["schur_inertia_sum"]
        assert max(d["symmetry"],d["joint_DtN_symmetry"])<policy["count_symmetry_tol"]
        count.append(d)
    green=max(single.SingleBeam(x,1.).green_identity(.47,policy["count_segments"]) for x in properties.values())
    assert green<policy["energy_tol"]
    # Single-beam axis covariance including both clamp rows, no new roots.
    old=read(previous.OUTPUT/config["nsl1_bundle"]/"result.json")
    modes=old["spectral"]["modes"]
    single_axis=[]
    t=joint.axis_transform(z0); ti=np.linalg.solve(t,np.eye(6))
    beam=single.SingleBeam(a,1.)
    for mode in modes:
        w=mode["omega"]
        raw=beam.basis(w,[0,1])[:,:3].reshape(6,6)
        shifted_rows=np.array([ti[:3]@t@b for b in beam.basis(w,[0,1])]).reshape(6,6)
        single_axis.append(relative(raw,shifted_rows))
    assert max(single_axis)<tol
    return {"status":"PASS","material_common_relative":material_error,
        "neutral_axis_offset":z0,"neutral_axis_properties":asdict(shifted),
        "rows":rows,"symmetric_RLB":symmetric,"negative_naive_axis_difference":wrong,
        "trace_counts":count,"segment_Green_identity":green,"single_axis_boundary":single_axis}


def energy_ritz_frequencies(properties,lengths,law,coefficients,quadrature):
    """Rayleigh values of the SAME Ritz vectors, evaluated in physical energies.

    Low LAPACK eigenvalues lose absolute digits against the largest stiffness.
    This is not a changed basis, eigenspace, tolerance or extra eigenvalue solve.
    """
    order=len(coefficients)//6
    nodes,weights=np.polynomial.legendre.leggauss(quadrature)
    values=[]
    for column in coefficients.T:
        mass=0.;energy=0.;ends=[]
        for i,(a,length) in enumerate(zip(properties,lengths)):
            cc=column[3*i*order:3*(i+1)*order]
            q,dq=joint.ritz_fields(cc,length,(nodes+1)*length/2)
            ww=weights*length/2
            mass+=np.sum(ww*np.einsum("ni,ij,nj->n",q,a.mass_matrix(),q))
            energy+=ww@(a.A*dq[:,0]**2+2*a.B*dq[:,0]*dq[:,2]+a.D*dq[:,2]**2+a.S*(dq[:,1]+q[:,2])**2)
            ends.append(joint.ritz_fields(cc,length,[length])[0][0,2])
        if law.mode!="RIGID":energy+=law.k_theta*(ends[0]-ends[1])**2
        values.append(np.sqrt(energy/mass))
    return np.asarray(values)


def spectra(config,properties,result,out,reuse=None):
    policy=config["policy"]; ref=config["reference"]; kt=ref["D"]/ref["length"]
    jobs=[("BASE_SPRING", "BASE",30.,1.),("MIRROR_SPRING","MIRROR",30.,1.),
          ("SAME_LOCAL_SPRING","SAME_LOCAL",30.,1.),
          ("MIRROR_HINGE","MIRROR",30.,0.),("MIRROR_RIGID","MIRROR",30.,None),
          ("STRAIGHT_RIGID","MIRROR",0.,None)]
    allcases=result.setdefault("spectral",{})
    profiles=[]; table=[]; runtime={}
    for label,name,beta_deg,kappa in jobs:
        props=[properties[x] for x in config["cases"][name]]
        law=joint.Joint("RIGID") if kappa is None else joint.Joint("SPRING",kappa*kt)
        model=joint.CoupledBeams(props,config["lengths"],np.deg2rad(beta_deg),law)
        case=allcases[label]={"status":"IN_PROGRESS","case":name,"beta_deg":beta_deg,"kappa":kappa}
        if reuse is not None and label in reuse["spectral"]:
            case["count"]=reuse["spectral"][label]["count"]
            roots=case["count"]["roots"]
            assert case["count"]["status"]=="PASS"
            case["root_origin"]="validated previous diagnostic; unchanged physical provider/config"
            result["solve_ledger"]["reused_propagation_roots"]+=len(roots)
        else:
            result["solve_ledger"]["propagation_prefix_searches"]+=1
            roots=model.prefix(policy,case.setdefault("count",{}))
            result["solve_ledger"]["propagation_roots"]+=len(roots)
        frequencies=np.array([x["omega"] for x in roots])
        print(label,"roots",frequencies.tolist(),flush=True)
        case["ritz"]=[]
        previous_w=None
        for order in policy["ritz_orders"]:
            result["solve_ledger"]["ritz_eigenproblems"]+=1
            raw,rc,rd=joint.ritz_spectrum(props,config["lengths"],model.beta,law,order,policy["ritz_quadrature"])
            rw=energy_ritz_frequencies(props,config["lengths"],law,rc,policy["ritz_quadrature"])
            record=dict(order=order,omega=rw,raw_lapack_omega=raw,diagnostics=rd,
                        eigenvalue_evaluation="physical Rayleigh quotient of unchanged Ritz vectors")
            if previous_w is not None:
                record["relative_change"]=abs(rw/previous_w-1)
                record["increase_relative"]=np.maximum(rw/previous_w-1,0)
            case["ritz"].append(record);previous_w=rw
        case["ritz_relative"]=abs(rw/frequencies-1)
        modes=[];case["modes"]=[]
        for k,w in enumerate(frequencies):
            coefficient,d=model.mode(w,policy)
            modes.append(d)
            ritz_mode={"arms":[]};signed_error=0.
            for i,(arm,ad) in enumerate(zip(model.arms,d["arms"])):
                q,dq=joint.ritz_fields(rc[i*3*order:(i+1)*3*order,k],arm.length,ad["x"])
                ritz_mode["arms"].append(dict(states=q,weights=ad["weights"]))
            overlap=joint.inner(d,ritz_mode,props)
            ritz_mass=joint.inner(ritz_mode,ritz_mode,props)
            for ad,rd,a in zip(d["arms"],ritz_mode["arms"],props):
                error=ad["states"][:,:3]-np.sign(overlap)*rd["states"]
                signed_error+=np.sum(ad["weights"]*np.einsum("ni,ij,nj->n",error,a.mass_matrix(),error))
            parity=joint.reflection_metrics(d,props) if name in ("BASE","MIRROR") else None
            axis=relative(model.boundary(w),model.boundary(w,offsets=[a.B/a.A for a in props]))
            # Actual transformed shapes: complete mass, stiffness and original
            # physical center/force recovery, not only internal H similarity.
            axis_m=0.;axis_k=0.;axis_field=0.
            for arm,ad in zip(model.arms,d["arms"]):
                z=arm.a.B/arm.a.A;t=joint.axis_transform(z);s=joint.shift_properties(arm.a,z)
                ys=ad["states"]@t.T; dys=ad["derivatives"]@t.T
                recovered=np.linalg.solve(t,ys.T).T
                axis_field=max(axis_field,relative(recovered,ad["states"]))
                q,dq=ys[:,:3],dys[:,:3]; weights=ad["weights"]
                axis_m+=np.sum(weights*np.einsum("ni,ij,nj->n",q,s.mass_matrix(),q))
                axis_k+=weights@(s.A*dq[:,0]**2+2*s.B*dq[:,0]*dq[:,2]+s.D*dq[:,2]**2+s.S*(dq[:,1]+q[:,2])**2)
            axis_k+=d["stiffness_components"]["spring"]
            metrics={n:d[n] for n in ("mass","stiffness","stiffness_components","delta_psi","spring_energy",
                "boundary_scaled","joint_scaled","ode_scaled","energy_relative","quadrature_relative","singular_ratio","nonzero_condition","joint_physical_residual")}
            metrics.update(position=k+1,omega=float(w),Omega=float(w*np.sqrt(ref["m"]/ref["D"])*ref["length"]**2),
                ritz_relative=float(case["ritz_relative"][k]),one_minus_mac=float(max(0,1-overlap**2/(d["mass"]*ritz_mass))),
                signed_shape_mass_error=float(np.sqrt(signed_error)),reflection=parity,
                axis_boundary=axis,axis_mass=abs(axis_m/d["mass"]-1),
                axis_energy=abs(axis_k/d["stiffness"]-1),axis_physical_field=axis_field,
                coefficients=coefficient)
            checks={
                "boundary":metrics["boundary_scaled"]<=policy["boundary_tol"],
                "joint":metrics["joint_scaled"]<=policy["boundary_tol"],
                "ODE":metrics["ode_scaled"]<=policy["ode_tol"],
                "energy":metrics["energy_relative"]<=policy["energy_tol"],
                "quadrature":metrics["quadrature_relative"]<=policy["energy_tol"],
                "frequency":metrics["ritz_relative"]<=policy["frequency_relative_tol"],
                "MAC":metrics["one_minus_mac"]<=policy["one_minus_mac_tol"],
                "shape":metrics["signed_shape_mass_error"]<=policy["signed_shape_mass_error_tol"],
                "axis":max(axis,metrics["axis_mass"],metrics["axis_energy"],axis_field)<=policy["axis_spectral_tol"],
                "symmetry":parity is None or parity["relative_error"]<=policy["symmetry_shape_tol"]}
            metrics["gates"]=checks;case["modes"].append(metrics)
            table.append({"case":label,**{n:metrics[n] for n in ("position","omega","Omega","mass","delta_psi","spring_energy","ritz_relative","boundary_scaled","joint_scaled","energy_relative")},
                **{"K_"+n:v for n,v in d["stiffness_components"].items()},
                "parity":None if parity is None else parity["candidate_parity"]})
            for i,ad in enumerate(d["arms"]):
                # Include endpoints in exported fields; Gauss nodes remain in runtime.
                xs=np.r_[0,ad["x"],config["lengths"][i]]
                ys=model.arms[i].basis(w,xs)@coefficient[6*i:6*i+6]
                for x,y in zip(xs,ys):
                    profiles.append({"case":label,"position":k+1,"arm":i+1,"x":float(x),
                        **dict(zip(("u","w","psi","N","Q","M"),map(float,y)))})
        gram=np.array([[joint.inner(a,b,props) for b in modes] for a in modes])
        case["mass_gram_error"]=float(np.max(abs(gram-np.eye(7))))
        conv=max(case["ritz"][-1]["relative_change"])
        monot=max(max(r.get("increase_relative",[0])) for r in case["ritz"])
        case["gates"]={"modes":all(all(m["gates"].values()) for m in case["modes"]),
            "Ritz_convergence":conv<=policy["ritz_convergence_tol"],
            "Ritz_nested":monot<=policy["ritz_monotonicity_tol"],
            "mass_orthogonality":case["mass_gram_error"]<=policy["mass_gram_tol"]}
        case["status"]="PASS" if all(case["gates"].values()) else "NUMERICAL_UNRESOLVED"
        runtime[label]=(model,modes)
        write_json(out/"result.json",plain(result))
        write_csv(out/"frequencies.csv",table);write_csv(out/"profiles.csv",profiles)
        # A failed physical gate stops further spectral searches.
        if case["status"]!="PASS" and label!="STRAIGHT_RIGID":
            raise ArithmeticError(label+" physical/Ritz gate unresolved")
    controls(config,properties,result,runtime)


def controls(config,properties,result,runtime):
    policy=config["policy"]; cases=result["spectral"]
    def omega(label):return np.array([m["omega"] for m in cases[label]["modes"]])
    h,s,r=(omega("MIRROR_"+n) for n in ("HINGE","SPRING","RIGID"))
    order_error=float(max(np.max((h-s)/s),np.max((s-r)/r),0))
    assert order_error<policy["ordering_relative_tol"]
    # Protected parity (-): select established classes, not arbitrary reordering.
    protected={}
    for kind in ("HINGE","SPRING","RIGID"):
        protected[kind]=[m for m in cases["MIRROR_"+kind]["modes"] if m["reflection"]["candidate_parity"]==-1]
    common=min(map(len,protected.values()))
    protected_error=max(relative([m["omega"] for m in protected[kind][:common]],
                                 [m["omega"] for m in protected["SPRING"][:common]]) for kind in protected)
    assert protected_error<policy["frequency_relative_tol"]
    result["spring_order"]={"status":"PASS","violation_relative":order_error,
        "protected_modes":{n:[{"position":m["position"],"omega":m["omega"],"delta_psi":m["delta_psi"]} for m in v] for n,v in protected.items()},
        "protected_frequency_relative":protected_error}
    # One simple target at kappa=1: two local correctors inside its certified
    # count-one interval. No global sweep or continuation is performed.
    k=policy["derivative_position"]-1
    original=cases["MIRROR_SPRING"]["modes"][k]; w=original["omega"]
    model,modes=runtime["MIRROR_SPRING"]
    gap=min(abs(w-x) for j,x in enumerate(s) if j!=k)
    assert gap/w>1e-5
    step=policy["derivative_kappa_step"]; roots=[]
    for multiplier in (1-step,1+step):
        target=joint.CoupledBeams(model.properties,model.lengths,model.beta,
            joint.Joint("SPRING",model.joint.k_theta*multiplier))
        result["solve_ledger"]["local_derivative_root_solves"]+=1
        roots.append(brentq(lambda x:np.linalg.det(target.boundary(x)),w-gap/4,w+gap/4,
                           xtol=policy["root_xtol"],rtol=policy["root_rtol"]))
    numerical=(roots[1]**2-roots[0]**2)/(2*step*model.joint.k_theta)
    predicted=original["delta_psi"]**2/original["mass"]
    result["stiffness_derivative"]=derivative_record(policy,k,w,gap,roots,predicted,numerical,model.joint.k_theta)
    # Independent length-2 propagation using unchanged NSL-1 solver.
    a=properties["HHLL"]; straight=single.SingleBeam(a,sum(config["lengths"]))
    checkpoint={};result["straight_single"]={"count":checkpoint}
    result["solve_ledger"]["propagation_prefix_searches"]+=1
    roots=straight.prefix(policy,checkpoint);result["solve_ledger"]["propagation_roots"]+=len(roots)
    wsingle=np.array([d["omega"] for d in roots])
    frequency_error=max(abs(wsingle/omega("STRAIGHT_RIGID")-1))
    coupled,cmodes=runtime["STRAIGHT_RIGID"]; comparisons=[]
    for k,w in enumerate(wsingle):
        coeff,sd=straight.mode(w,policy)
        full={"arms":[]}
        for i,ad in enumerate(cmodes[k]["arms"]):
            x=ad["x"] if i==0 else straight.length-ad["x"]
            y=straight.basis(w,x)@coeff
            if i==1:y=y@joint.CONTINUOUS
            full["arms"].append(dict(states=y,weights=ad["weights"]))
        overlap=joint.inner(cmodes[k],full,coupled.properties)
        error=0.;force_error=0.
        for ad,fd in zip(cmodes[k]["arms"],full["arms"]):
            diff=ad["states"]-np.sign(overlap)*fd["states"]
            error=max(error,float(np.max(abs(diff[:,:3]))/max(np.max(abs(ad["states"][:,:3])),1e-30)))
            force_error=max(force_error,float(np.max(abs(diff[:,3:]))/max(np.max(abs(ad["states"][:,3:])),1e-30)))
        comparisons.append(dict(position=k+1,one_minus_mac=max(0,1-overlap**2),
            displacement_relative=error,force_relative=force_error,
            single_energy_relative=sd["energy_relative"],single_boundary=sd["boundary_scaled"]))
    assert frequency_error<policy["frequency_relative_tol"]
    assert max(d["displacement_relative"] for d in comparisons)<policy["signed_shape_mass_error_tol"]
    assert max(d["force_relative"] for d in comparisons)<policy["signed_shape_mass_error_tol"]
    single_physical=all(d["single_boundary"]<=policy["boundary_tol"] and
                        d["single_energy_relative"]<=policy["energy_tol"] for d in comparisons)
    result["straight_single"].update(
        status="PASS" if single_physical and cases["STRAIGHT_RIGID"]["status"]=="PASS" else "PARTIAL_PASS",
        frequency_and_shape_comparison="PASS",single_physical_gates=single_physical,
        frequency_relative=float(frequency_error),modes=comparisons)
    # Cross-case comparisons in the genuinely common physical mass (same m,j2,j1).
    comparisons={}
    for left,right in (("BASE_SPRING","MIRROR_SPRING"),("MIRROR_SPRING","SAME_LOCAL_SPRING")):
        lm,lf=runtime[left];rm,rf=runtime[right]
        mac=np.array([[joint.inner(x,y,lm.properties)**2 for y in rf] for x in lf])
        ii,jj=linear_sum_assignment(-mac)
        pairs=[]
        for i,j in zip(ii,jj):
            margin=mac[i,j]-max(np.delete(mac[i],j))
            accepted=mac[i,j]>=policy["cross_case_mac_min"] and margin>=policy["cross_case_mac_margin"]
            pairs.append(dict(left_position=int(i+1),right_position=int(j+1),MAC=float(mac[i,j]),
                margin=float(margin),identity="RESOLVED" if accepted else "UNRESOLVED"))
        comparisons[left+"__"+right]={"MAC":mac,"assignments":pairs,
            "sorted_frequency_difference_relative":omega(right)/omega(left)-1}
    result["cross_case"]=comparisons


def derivative_record(policy,k,w,gap,roots,predicted,numerical,stiffness):
    active=predicted>1e-10
    error=abs(numerical/predicted-1) if active else None
    if active:assert error<policy["derivative_relative_tol"]
    return {"position":k+1,"simple_gap_relative":gap/w,"target_roots":roots,
        "energy_prediction":predicted,"finite_difference":numerical,"relative_error":error,
        "absolute_error":abs(numerical-predicted),
        "dimensionless_zero_residual":abs(numerical-predicted)*stiffness/w**2,
        "status":"PASS" if active else "STRUCTURAL_ZERO; RELATIVE_NUMERICAL_TEST_NOT_APPLICABLE",
        "qualification":None if active else "Preselected position 1 belongs to eta=-1. Exact zero follows from reflection. Finite difference is roundoff-level; no ratio to numerical zero and no extra active mode solve."}


def export_controls(out,result):
    write_csv(out/"straight_comparison.csv",result["straight_single"]["modes"])
    rows=[]
    for pair,c in result["cross_case"].items():
        for m in c["assignments"]:
            rows.append({"pair":pair,**m,"sorted_difference_relative":c["sorted_frequency_difference_relative"][m["left_position"]-1]})
    write_csv(out/"cross_case_comparison.csv",rows)
    write_json(out/"controls.json",plain({n:result[n] for n in
        ("spring_order","stiffness_derivative","straight_single","cross_case")}))


def identity(config):
    baseline=read(OUTPUT/"baseline.json")
    assert all(sha(ROOT/n)==h for n,h in baseline["sha256"].items()),"Protected source/bundle changed"
    old=read(previous.CONFIG);fp,ident=previous.identity(old)
    manifest=read(previous.OUTPUT/fp/"manifest.json")
    assert fp==config["nsl1_bundle"] and manifest["identity"]==ident
    assert all(sha(previous.OUTPUT/fp/n)==h for n,h in manifest["artifact_hashes"].items())
    paths=[CONFIG,Path(__file__),Path(joint.__file__),Path(single.__file__)]
    data={"config":config,"baseline_sha256":sha(OUTPUT/"baseline.json"),
        "source_hashes":baseline["sha256"],"files":{p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        "nsl1_verified_bundle":fp,"versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()[:16],data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--structural-only",action="store_true")
    parser.add_argument("--reuse-roots",help="Validated diagnostic bundle, same physical provider/config; no repeated root searches")
    parser.add_argument("--finish-diagnostics",help="Complete prescribed independent controls from saved shapes, keeping failed physical gates")
    args=parser.parse_args();config=read(CONFIG);fp,ident=identity(config);out=OUTPUT/fp
    if (out/"manifest.json").exists():
        manifest=read(out/"manifest.json")
        assert manifest["identity"]==ident
        assert all(sha(out/n)==h for n,h in manifest["artifact_hashes"].items())
        saved=read(out/"result.json")
        if saved["status"]!="IN_PROGRESS" or not args.compute:
            print(json.dumps(dict(status=saved["status"],cache="VALIDATED",fingerprint=fp,new_eigenvalue_solves=0)));return
    if not(args.compute or args.structural_only or args.finish_diagnostics):
        print(json.dumps(dict(source_identity="PASS",fingerprint=fp,new_eigenvalue_solves=0)));return
    out.mkdir(parents=True,exist_ok=True)
    result={"status":"IN_PROGRESS","parameters":config,"solve_ledger":{
        "propagation_prefix_searches":0,"propagation_roots":0,"ritz_eigenproblems":0,
        "reused_propagation_roots":0,"local_derivative_root_solves":0,"axis_or_reversal_root_searches":0,"complex_root_solves":0}}
    try:
        result["symbolic"]=symbolic()
        properties,result["materials"]=materials(config)
        result["structural"]=structural(config,properties)
        write_json(out/"result.json",plain(result))
        if args.finish_diagnostics:
            origin=OUTPUT/args.finish_diagnostics
            manifest=read(origin/"manifest.json")
            assert manifest["identity"]["config"]==config
            assert manifest["identity"]["source_hashes"]==ident["source_hashes"]
            for path in (Path(joint.__file__),Path(single.__file__),CONFIG):
                assert manifest["identity"]["files"][path.relative_to(ROOT).as_posix()]==sha(path)
            assert all(sha(origin/n)==h for n,h in manifest["artifact_hashes"].items())
            prior=read(origin/"result.json")
            assert prior["materials"]==plain(result["materials"])
            result["spectral"]=prior["spectral"]
            result["diagnostic_reuse"]={"fingerprint":args.finish_diagnostics,
                "manifest_sha256":sha(origin/"manifest.json"),"prior_solve_ledger":prior["solve_ledger"],
                "prior_root_reuse":prior.get("root_reuse"),
                "reason":"Preserve straight guard residual failure; finish independent requested controls only",
                "verified":"all artifacts, physical providers, parameters and materials"}
            for name in ("frequencies.csv","profiles.csv"):
                shutil.copyfile(origin/name,out/name)
            if "cross_case" in prior:
                for name in ("spring_order","stiffness_derivative","straight_single","cross_case"):
                    result[name]=prior[name]
                result["diagnostic_reuse"]["prior_diagnostic_reuse"]=prior.get("diagnostic_reuse")
                old=result["stiffness_derivative"];k=config["policy"]["derivative_position"]-1
                w=result["spectral"]["MIRROR_SPRING"]["modes"][k]["omega"]
                result["stiffness_derivative"]=derivative_record(config["policy"],k,w,old["simple_gap_relative"]*w,
                    old["target_roots"],old["energy_prediction"],old["finite_difference"],config["reference"]["D"]/config["reference"]["length"])
            else:
                runtime={}
                for name,case in result["spectral"].items():
                    pair=config["cases"][case["case"]]
                    props=[properties[n] for n in pair]
                    law=joint.Joint("RIGID") if case["kappa"] is None else joint.Joint("SPRING",case["kappa"]*config["reference"]["D"]/config["reference"]["length"])
                    model=joint.CoupledBeams(props,config["lengths"],np.deg2rad(case["beta_deg"]),law)
                    runtime[name]=(model,[model.evaluate(m["omega"],np.asarray(m["coefficients"]),config["policy"]["mode_quadrature"]) for m in case["modes"]])
                controls(config,properties,result,runtime)
            result["status"]="PASS" if all(c["status"]=="PASS" for c in result["spectral"].values()) and result["straight_single"]["status"]=="PASS" else "PARTIAL_PASS"
            result["qualification"]="Straight controls retain BC gates: coupled guard7 and single length2 positions6,7. Main three-case pilot and minmax spectra passed; no tolerances changed."
            export_controls(out,result)
        elif args.compute:
            reuse=None
            if args.reuse_roots:
                origin=OUTPUT/args.reuse_roots
                manifest=read(origin/"manifest.json")
                assert manifest["identity"]["config"]==config
                assert manifest["identity"]["source_hashes"]==ident["source_hashes"]
                for path in (Path(joint.__file__),Path(single.__file__),CONFIG):
                    assert manifest["identity"]["files"][path.relative_to(ROOT).as_posix()]==sha(path)
                assert all(sha(origin/n)==h for n,h in manifest["artifact_hashes"].items())
                reuse=read(origin/"result.json")
                assert reuse["materials"]==plain(result["materials"])
                result["root_reuse"]={"fingerprint":args.reuse_roots,"manifest_sha256":sha(origin/"manifest.json"),
                    "verified":"config, all protected sources, NSL-1 physics, coupled provider, material outputs, artifacts",
                    "prior_solve_ledger":reuse["solve_ledger"],
                    "reason":"Raw LAPACK Ritz monotonicity diagnostic; use physical Rayleigh evaluation, unchanged tolerances/basis/roots"}
            spectra(config,properties,result,out,reuse)
            result["status"]="PASS" if all(c["status"]=="PASS" for c in result["spectral"].values()) and result["straight_single"]["status"]=="PASS" else "PARTIAL_PASS"
            result["qualification"]="See individual physical gates; a completed root count does not waive endpoint residual checks."
            export_controls(out,result)
    except Exception as error:
        result["status"]="PARTIAL_PASS" if result.get("structural",{}).get("status")=="PASS" else "FAIL"
        result["stop_reason"]=f"{type(error).__name__}: {error}"
        import traceback
        result["traceback"]=traceback.format_exc()
    result["source_preservation"]="PASS" if all(sha(ROOT/n)==h for n,h in ident["source_hashes"].items()) else "FAIL"
    write_json(out/"result.json",plain(result))
    write_json(out/"manifest.json",{"identity":ident,
        "command":subprocess.list2cmdline([sys.executable,*sys.argv]),
        "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "artifact_hashes":{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}})
    write_json(OUTPUT/"current.json",{"fingerprint":fp,"status":result["status"]})
    print(json.dumps(dict(status=result["status"],fingerprint=fp,solve_ledger=result["solve_ledger"],stop_reason=result.get("stop_reason"))))
    if "stop_reason" in result:raise SystemExit(1)


if __name__=="__main__":
    main()
