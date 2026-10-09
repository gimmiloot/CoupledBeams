"""FEM-2 bounded static comparison of frozen quartic1D action and 3D StVK.
No new mesh, modal calculation, dynamics or parameter fitting.
"""
from __future__ import annotations
import os,sys
from pathlib import Path
if __name__=='__main__':
    for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
"""FEM-2 scoped static orchestration; reuse frozen planar action and basis."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json
import math
import time
import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigvalsh
from scripts.lib import weakly_nonlinear_planar_dynamics as fem2_dyn
from scripts.lib import weakly_nonlinear_spatial_rod as fem2_rod

FEM2_FIELDS = ("u", "w", "theta", "c")


def fem2_tim_uniform(x, load, length, Bp, S):
    """Independent static integration Q'=-q, M'=-Q, theta'=M/Bp."""
    x = np.asarray(x)
    w = load*(x*x*(length-x)**2/(24*Bp)+x*(length-x)/(2*S))
    theta = load*(x**3/6-length*x*x/4+length*length*x/12)/Bp
    ws = theta+load*(length/2-x)/S
    thetas = load*(x*x/2-length*x/2+length*length/12)/Bp
    return np.column_stack((np.zeros_like(x), w, theta, np.zeros_like(x))), np.column_stack((np.zeros_like(x), ws, thetas, np.zeros_like(x)))


def fem2_line_load(disc, q):
    f = np.zeros(disc.ndof)
    f[disc.slices['w']] = q*(disc.B[1].T@disc.weights)
    return f


def fem2_load_selection(c, coefficients):
    """Select solely from the analytic 1D linear compliance before FEM."""
    g = c['geometry']; m = c['material']; p = c.get('load_policy', {})
    L, h, area = g['L'], g['h'], g['b']*g['h']
    compliance = L**4/(384*coefficients.Bp)+L**2/(8*coefficients.S)
    ceiling = p.get('bending_surface_strain_ceiling', .01)
    candidates = [p.get('primary_w_over_h', .05), p.get('backup_w_over_h', .03)]
    tested = []
    for ratio in candidates:
        q = ratio*h/compliance
        grid = np.linspace(0, L, 1001)
        fields, gradients = fem2_tim_uniform(grid, q, L, coefficients.Bp, coefficients.S)
        surface = h/2*float(np.max(abs(gradients[:,2])))
        row = {'target_w_over_h':ratio, 'g':q/(m['rho']*area), 'q':q, 'F_total':q*L,
               'linear_w_max':ratio*h, 'linear_max_slope':float(np.max(abs(gradients[:,1]))),
               'linear_max_rotation':float(np.max(abs(fields[:,2]))),
               'linear_max_shear_strain':float(np.max(abs(gradients[:,1]-fields[:,2]))),
               'linear_bending_surface_strain':surface,
               'bending_surface_strain_ceiling':ceiling}
        tested.append(row)
        if surface <= ceiling:
            return {**row, 'status':'PASS', 'candidate_history':tested,
                    'selection':'independent fixed-fixed linear Timoshenko compliance; before FEM',
                    'body_direction_global':[0., -1., 0.], 'reference_area':area,
                    'load_units':'acceleration g; force per original length q=rho*A0*g',
                    'follower_load':False}
    return {'status':'FAIL', 'candidate_history':tested, 'reason':'Both predeclared candidates exceed surface-strain guide'}


def fem2_static_newton(disc, force, policy):
    """Unloaded branch, all four fields; only energy gradient/Hessian used."""
    rt = policy.get('residual_relative_tolerance', 1e-10)
    it = policy.get('increment_relative_tolerance', 1e-12)
    maxiter = policy.get('newton_max_iterations', 20)
    subdivisions = policy.get('max_load_subdivisions', 6)
    load_steps = policy.get('load_steps', 10)
    a = np.zeros(disc.ndof)
    history = [{'load_factor':0., 'relative_residual':0., 'relative_increment':0., 'iterations':0}]
    iteration_records = []
    pending = [(float(i)/load_steps,0) for i in range(1, load_steps+1)]
    previous_factor = 0.; accepted_steps = 0; subdivisions_used = 0
    while pending:
        factor, depth = pending.pop(0)
        x = a.copy(); target = factor*force
        inc = float('inf'); converged = False
        for iteration in range(maxiter):
            value = disc.potential(x, hessian=True)
            r = value['gradient']-target
            scale = max(np.linalg.norm(value['gradient'])+np.linalg.norm(target), np.linalg.norm(force)*1e-15, 1e-30)
            relative = float(np.linalg.norm(r)/scale)
            tangent = value['hessian']
            try:
                dx = cho_solve(cho_factor(tangent, lower=True, check_finite=True), -r, check_finite=True)
            except np.linalg.LinAlgError:
                break
            inc = float(np.linalg.norm(dx)/max(np.linalg.norm(x), np.linalg.norm(a), 1e-30))
            iteration_records.append({'load_factor':factor,'iteration':iteration,'relative_residual':relative,
                                      'absolute_residual':float(np.linalg.norm(r)), 'relative_increment':inc})
            if relative <= rt and inc <= it:
                converged = True; break
            old_total = value['V']-target@x
            old_residual = float(np.linalg.norm(r))
            # Energy descent with residual fall-back near cancellation floor.
            alpha = 1.; found = False
            for backtrack in range(13):
                candidate = x+alpha*dx
                if not np.all(np.isfinite(candidate)):
                    alpha /= 2; continue
                cv = disc.potential(candidate)
                cr = float(np.linalg.norm(cv['gradient']-target))
                ct = cv['V']-target@candidate
                if ct <= old_total+1e-14*max(abs(old_total),1e-30) or cr < old_residual:
                    x = candidate; found = True; break
                alpha /= 2
            if not found: break
        if not converged:
            if depth >= subdivisions:
                return {'status':'FAIL','reason':'Newton failed at bounded load subdivision','load_factor':factor,
                        'coordinate':a,'history':history,'iterations':iteration_records,
                        'subdivisions_used':subdivisions_used}
            midpoint = (previous_factor+factor)/2
            pending = [(midpoint,depth+1),(factor,depth+1)]+pending
            subdivisions_used += 1; continue
        a = x; previous_factor = factor; accepted_steps += 1
        eig = eigvalsh(tangent, check_finite=False, subset_by_index=[0,0])
        history.append({'load_factor':factor, 'iterations':iteration+1,'relative_residual':relative,
                        'relative_increment':inc, 'minimum_tangent_eigenvalue':float(eig[0]),
                        'potential':float(value['V']), 'load_potential':float(-target@a),
                        'total_potential':float(value['V']-target@a)})
        if eig[0] <= 0 or not np.all(np.isfinite(a)):
            return {'status':'FAIL','reason':'Nonpositive tangent or nonfinite coordinates', 'coordinate':a,
                    'history':history,'iterations':iteration_records,'subdivisions_used':subdivisions_used}
    return {'status':'PASS','coordinate':a, 'history':history,'iterations':iteration_records,
            'subdivisions_used':subdivisions_used, 'accepted_load_steps':accepted_steps,
            'load_factor':previous_factor, 'relative_residual':history[-1]['relative_residual'],
            'relative_increment':history[-1]['relative_increment']}


def fem2_static_fields(disc, coordinate, points, q, linear=False):
    """Coordinate-conjugate fluxes, not reactions inferred from equilibrium."""
    fields = disc.reconstruct(coordinate, points)
    derivatives = disc.reconstruct(coordinate, points, 1)
    second = disc.reconstruct(coordinate, points, 2)
    variables = np.vstack((derivatives[:,0],derivatives[:,1],fields[:,2],fields[:,3],derivatives[:,3],derivatives[:,2]))
    potential = fem2_dyn._restrict_to_plane(disc.model.V4)
    if linear: potential = potential.truncate(2)
    names = fem2_dyn._LOCAL_POTENTIAL_NAMES
    compiled = fem2_dyn._CompiledPolynomials([potential.derivative(n) for n in names], names, disc.coefficients)
    gradient = compiled.evaluate(variables)
    flux = gradient[[0,1,5,4]].T
    theta,us,ws,c = fields[:,2],derivatives[:,0],derivatives[:,1],fields[:,3]
    gamma1 = us+theta*ws-theta**2/2-us*theta**2/2
    gamma2 = ws-theta-theta*us-ws*theta**2/2+theta**3/6
    if linear: gamma1, gamma2 = us, ws-theta
    return {'fields':fields, 'first_derivatives':derivatives,'second_derivatives':second,
            'flux':flux, 'Gamma1':gamma1,'Gamma2':gamma2}


def fem2_static_diagnostics(disc, coordinate, q, geometry, linear=False):
    points = np.linspace(0, disc.length, 1001)
    profile = fem2_static_fields(disc, coordinate, points, q, linear)
    ends = fem2_static_fields(disc, coordinate, np.array([0.,disc.length]), q, linear)
    reactions = ends['flux']*np.array([-1.,1.])[:,None]
    # Local basis t=EX,n=-EY,k=-EZ. End displacement is exactly zero.
    physical = np.column_stack((reactions[:,0],-reactions[:,1],np.zeros(2)))
    couples_z = -reactions[:,2]
    force_balance = physical.sum(axis=0)+np.array([0.,-q*disc.length,0.])
    moment_balance = float(couples_z.sum()+disc.length*physical[1,1]-q*disc.length**2/2)
    force_scale = max(abs(q*disc.length), 1e-30)
    moment_scale = force_scale*disc.length
    field_scale = [geometry['h'],geometry['h'],geometry['h']/disc.length,1.]
    sym_sign = np.array([-1.,1.,-1.,1.])
    parity = profile['fields']-profile['fields'][::-1]*sym_sign
    bcs = ends['fields']
    value = disc.potential(coordinate, hessian=True)
    load = fem2_line_load(disc,q)
    assembled = disc.linear_stiffness@coordinate-load if linear else value['gradient']-load
    scale = max(np.linalg.norm(load)+(np.linalg.norm(disc.linear_stiffness@coordinate) if linear else np.linalg.norm(value['gradient'])),1e-30)
    # Independent strong residual projection. The existing float64 strict
    # 2e-12 strong/weak qualification is not replaced by a new threshold.
    zeros = np.zeros(disc.ndof)
    strong_local=[]
    for state,d in ((coordinate,0),(coordinate,1),(zeros,0),(coordinate,2),(zeros,1),(zeros,0)):
        strong_local.extend(disc.reconstruct(state,derivative=d).T)
    if linear:
        names=tuple(field+suffix for suffix in fem2_dyn._JET_SUFFIXES for field in FEM2_FIELDS)
        evaluator=fem2_dyn._CompiledPolynomials([
            fem2_dyn._restrict_to_plane(disc.model.residual_a[index]).truncate(1)
            for index in fem2_dyn.PLANAR_SPATIAL_INDICES],names,disc.coefficients)
        strong=evaluator.evaluate(np.asarray(strong_local)).T
        weak=np.concatenate([matrix.T@(disc.weights*values) for matrix,values in zip(disc.B,strong.T)])-load
        action_gradient=disc.linear_stiffness@coordinate
    else:
        strong=disc._residual.evaluate(np.asarray(strong_local)).T
        weak=disc.weak_residual(coordinate,zeros,zeros)-load
        action_gradient=value['gradient']
    strong[:,1] -= q
    plane=fem2_dyn._restrict_to_plane(disc.model.V4)
    if linear: plane=plane.truncate(2)
    local_evaluator=fem2_dyn._CompiledPolynomials([plane.derivative(n) for n in fem2_dyn._LOCAL_POTENTIAL_NAMES],
                                                 fem2_dyn._LOCAL_POTENTIAL_NAMES,disc.coefficients)
    local_gradient=local_evaluator.evaluate(disc._local_variables(coordinate))
    uncancelled_work_scale=sum(np.linalg.norm(matrix.T@(disc.weights*values))
                              for matrix,values in zip(disc._potential_matrices,local_gradient))
    strong_delta=weak-(action_gradient-load)
    strong_relative=float(np.linalg.norm(strong_delta)/max(uncancelled_work_scale,1e-30))
    strong_absolute=float(np.max(abs(strong_delta)))
    degree_values=[]
    for degree in (2,3,4):
        pol = fem2_dyn._restrict_to_plane(disc.model.V4).truncate(degree)-fem2_dyn._restrict_to_plane(disc.model.V4).truncate(degree-1)
        ev = fem2_dyn._CompiledPolynomials([pol],fem2_dyn._LOCAL_POTENTIAL_NAMES,disc.coefficients)
        degree_values.append(float(disc.weights@ev.evaluate(disc._local_variables(coordinate))[0]))
    homogeneous = sum(degree*v for degree,v in zip((2,3,4),degree_values))
    if linear: homogeneous=2*degree_values[0]
    virtual_work = float(coordinate@(disc.linear_stiffness@coordinate if linear else value['gradient']))
    tangent = disc.linear_stiffness if linear else value['hessian']
    ev = eigvalsh(tangent, check_finite=False)
    result = {'endpoint_BC_absolute_max':float(np.max(abs(bcs))),
              'essential_fields':list(FEM2_FIELDS),'slope_constraints':False,
              'relative_residual':float(np.linalg.norm(assembled)/scale),
              'absolute_residual':float(np.linalg.norm(assembled)),
              'minimum_tangent_eigenvalue':float(ev[0]),'maximum_tangent_eigenvalue':float(ev[-1]),
              'tangent_condition':float(ev[-1]/ev[0]),'tangent_positive':bool(ev[0]>0),
              'support_reactions_local_u_w_M_Rc':reactions.tolist(),
              'support_reactions_global_xyz':physical.tolist(),'support_couples_global_z':couples_z.tolist(),
              'force_balance_absolute':force_balance.tolist(),'force_balance_relative':float(np.linalg.norm(force_balance)/force_scale),
              'moment_balance_absolute_z':moment_balance,'moment_balance_relative':abs(moment_balance)/moment_scale,
              'reaction_source':'endpoint partial V_le4 / partial q_s; sigma=-1 left,+1 right',
              'symmetry_fixed_scaled_max':float(np.max(abs(parity)/np.asarray(field_scale))),
              'max_abs_us':float(np.max(abs(profile['first_derivatives'][:,0]))),
              'max_abs_ws':float(np.max(abs(profile['first_derivatives'][:,1]))),
              'max_abs_theta':float(np.max(abs(profile['fields'][:,2]))),
              'max_abs_c':float(np.max(abs(profile['fields'][:,3]))),
              'min_one_plus_c':float(np.min(1+profile['fields'][:,3])),
              'max_abs_Gamma1':float(np.max(abs(profile['Gamma1']))),
              'max_abs_Gamma2':float(np.max(abs(profile['Gamma2']))),
              'max_abs_bending_surface_strain':float(geometry['h']/2*np.max(abs(profile['first_derivatives'][:,2]))),
              'V':degree_values[0] if linear else value['V'],'V2_V3_V4':degree_values,
              'external_load_work_at_state':float(load@coordinate),
              'equilibrium_virtual_work':virtual_work,
              'homogeneous_potential_work':homogeneous,
              'homogeneous_work_identity_relative':abs(homogeneous-virtual_work)/max(abs(homogeneous)+abs(virtual_work),1e-30),
              'independent_projected_strong_residual_relative':float(np.linalg.norm(weak)/scale),
              'strong_action_relative_difference':strong_relative,
              'strong_action_absolute_difference':strong_absolute,
              'strong_action_uncancelled_work_scale':float(uncancelled_work_scale),
              'strong_action_strict_threshold':2e-12,
              'strong_action_strict_status':'PASS' if strong_relative <= 2e-12 and strong_absolute <= 2e-12 else 'PARTIAL',
              'sampled_strong_residual_L2':np.sqrt(disc.weights@(strong*strong)).tolist(),
              'sampled_strong_residual_max':np.max(abs(strong),axis=0).tolist(),
              'strong_residual_qualification':'finite-p interior diagnostic, not imposed at essential endpoints',
              'dynamics_rhs_calls':disc.rhs_calls,'dynamics_jacobian_calls':disc.jacobian_calls,
              'linear_eigendecompositions':disc.linear_eigendecompositions}
    return result,profile


def fem2_static_comparison(low,high,x,length):
    """All physical fields plus Delta w; no phase/amplitude or scale fitting."""
    # Physical norms on the saved equally spaced output; sampled maxima.
    # No profile amplitude or phase fitting is used.
    rows=[]
    for field,index in zip(FEM2_FIELDS,range(4)):
        a,b=low[:,index],high[:,index]; dif=a-b
        abs_l2=float(np.sqrt(np.trapezoid(dif*dif,x)))
        ref_l2=float(np.sqrt(np.trapezoid(b*b,x)))
        ref_max=float(np.max(abs(b))); abs_max=float(np.max(abs(dif)))
        floor=1e-14*(.1 if field in ('u','w') else 1.)
        rows.append({'field':field,'absolute_L2':abs_l2,'absolute_max':abs_max,
                     'relative_L2':abs_l2/max(ref_l2,floor*math.sqrt(length)),
                     'relative_max':abs_max/max(ref_max,floor),
                     'reference_L2':ref_l2,'reference_max':ref_max})
    return rows


def build_static_preflight(c,bundle,root=None):
    root=Path.cwd() if root is None else Path(root); bundle=Path(bundle); bundle.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    action_path=root/c.get('action_bundle','results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295')/'result.json'
    manifest=json.loads((action_path.parent/'manifest.json').read_text(encoding='utf8'))
    artifacts=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
    for name,digest in artifacts.items():
        actual=hashlib.sha256((action_path.parent/name).read_bytes()).hexdigest()
        if actual!=digest: raise ValueError('Frozen action archive corrupted: '+name)
    action=json.loads(action_path.read_text(encoding='utf8')); pol=action['polynomials']
    model=SimpleNamespace(T4=fem2_rod.Polynomial.deserialize(pol['T4']),V4=fem2_rod.Polynomial.deserialize(pol['V4']),
         residual_a=tuple(fem2_rod.Polynomial.deserialize(a) for a in pol['residuals_A']),
         symbols={a:fem2_rod.Polynomial.symbol(a) for a in fem2_rod.SYMBOL_ORDER})
    parent=root/c.get('parent_bundle','results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad')
    source=json.loads((parent/'preflight.json').read_text(encoding='utf8'))
    coefficients=fem2_rod.RodCoefficients(**source['coefficients'])
    L=c['geometry']['L']; gm=c['geometry']; mat=c['material']
    expected={'m':mat['rho']*gm['b']*gm['h'],'jp':mat['rho']*gm['b']*gm['h']**3/12,
              'C':mat['E']*gm['b']*gm['h']/(1-mat['nu']**2),
              'H':mat['kappa']*mat['E']/(2*(1+mat['nu']))*gm['b']*gm['h']**3/12,
              'S':mat['kappa']*mat['E']/(2*(1+mat['nu']))*gm['b']*gm['h'],
              'Bp':mat['E']*gm['b']*gm['h']**3/12}
    if any(not math.isclose(getattr(coefficients,k),v,rel_tol=3e-15,abs_tol=0) for k,v in expected.items()):
        raise ValueError('Frozen reference coefficients/geometry mismatch')
    selection=fem2_load_selection(c,coefficients)
    if selection['status']!='PASS': return {'load':selection,'status':'FAIL'}
    q=selection['q']; x=np.linspace(0,L,1001); cases=[]; profiles={}; policy=c.get('one_d',{})
    for p in policy.get('p',[48,64]):
        disc=fem2_dyn.PlanarGalerkin(coefficients,int(p),length=L,nq=2*int(p)+1,model=model)
        force=fem2_line_load(disc,q)
        linear=cho_solve(cho_factor(disc.linear_stiffness,lower=True,check_finite=True),force)
        nl=fem2_static_newton(disc,force,policy)
        if nl['status']!='PASS':
            return {'status':'FAIL','load':selection,'failed_degree':p,'reason':nl['reason'],
                    'history':nl['history'],'iteration_records':nl['iterations'], 'runtime_seconds':time.perf_counter()-started}
        ld,lp=fem2_static_diagnostics(disc,linear,q,gm,True)
        nd,np_=fem2_static_diagnostics(disc,nl['coordinate'],q,gm,False)
        analytic,analytic_d=fem2_tim_uniform(x,q,L,coefficients.Bp,coefficients.S)
        analytic_difference=float(np.max(abs(lp['fields']-analytic))/selection['linear_w_max'])
        analytic_derivative_difference=float(np.max(abs(lp['first_derivatives']-analytic_d))/max(np.max(abs(analytic_d)),1e-30))
        profile_path=bundle/f'one_d_p{p}.npz'
        arrays={'s':x,'q_linear':linear,'q_nonlinear':nl['coordinate'],'raw_linear':disc.raw_coefficients(linear),
                'raw_nonlinear':disc.raw_coefficients(nl['coordinate']),'linear':lp['fields'], 'nonlinear':np_['fields'],
                'linear_first':lp['first_derivatives'],'nonlinear_first':np_['first_derivatives'],
                'linear_second':lp['second_derivatives'],'nonlinear_second':np_['second_derivatives'],
                'linear_flux':lp['flux'],'nonlinear_flux':np_['flux'],
                'Gamma1_linear':lp['Gamma1'],'Gamma2_linear':lp['Gamma2'],
                'Gamma1_nonlinear':np_['Gamma1'],'Gamma2_nonlinear':np_['Gamma2'],
                'delta_w':np_['fields'][:,1]-lp['fields'][:,1]}
        np.savez_compressed(profile_path,**arrays)
        profiles[p]=arrays
        record={'p':p,'nq':disc.nq,'ndof':disc.ndof,'fields':list(FEM2_FIELDS),
                'basis':'P_n-P_(n+2),degree p; reversible resting-mass whitening',
                'profiles':profile_path.name,'linear':ld,'nonlinear':nd,
                'linear_analytic_profile_relative_max':analytic_difference,
                'linear_analytic_derivative_relative_max':analytic_derivative_difference,
                'load_steps':nl['history'],'newton_iterations':nl['iterations'],
                'load_subdivisions':nl['subdivisions_used'],
                'linear_midspan_w':float(lp['fields'][500,1]),
                'nonlinear_midspan_w':float(np_['fields'][500,1]),
                'delta_w_midspan':float(arrays['delta_w'][500]),
                'delta_w_absolute_max':float(np.max(abs(arrays['delta_w']))),
                'nonlinear_u_absolute_max':float(np.max(abs(np_['fields'][:,0]))),
                'nonlinear_c_absolute_max':float(np.max(abs(np_['fields'][:,3])))}
        record['linear_status']='PASS' if ld['relative_residual']<=policy.get('residual_relative_tolerance',1e-10) and analytic_difference<=1e-10 and ld['endpoint_BC_absolute_max']==0 else 'FAIL'
        record['nonlinear_status']='PASS' if nd['relative_residual']<=policy.get('residual_relative_tolerance',1e-10) and nd['endpoint_BC_absolute_max']==0 and nd['tangent_positive'] else 'FAIL'
        cases.append(record)
    low,high=profiles[48],profiles[64]
    comp={'linear':fem2_static_comparison(low['linear'],high['linear'],x,L),
          'nonlinear':fem2_static_comparison(low['nonlinear'],high['nonlinear'],x,L)}
    da=low['delta_w']; db=high['delta_w']; dd=da-db
    comp['delta_w']={'absolute_L2':float(np.sqrt(np.trapezoid(dd*dd,x))), 'absolute_max':float(np.max(abs(dd))),
                     'relative_L2':float(np.sqrt(np.trapezoid(dd*dd,x))/max(np.sqrt(np.trapezoid(db*db,x)),1e-30)),
                     'relative_max':float(np.max(abs(dd))/max(np.max(abs(db)),1e-30))}
    rlow=np.asarray(cases[0]['nonlinear']['support_reactions_local_u_w_M_Rc'])
    rhigh=np.asarray(cases[1]['nonlinear']['support_reactions_local_u_w_M_Rc'])
    comp['reactions']={'absolute_difference':abs(rlow-rhigh).tolist(),
                       'relative_difference':(abs(rlow-rhigh)/np.maximum(abs(rhigh),selection['F_total']*1e-12)).tolist()}
    result={'status':'PASS' if all(row['linear_status']=='PASS' and row['nonlinear_status']=='PASS' for row in cases) else 'FAIL',
            'load':selection,'coefficients':coefficients.values(),'cases':cases,'spatial_comparison':comp,
            'action_source':str(action_path.relative_to(root)),
            'action_result_sha256':hashlib.sha256(action_path.read_bytes()).hexdigest(),
            'coefficient_source':str((parent/'preflight.json').relative_to(root)),
            'symbolic_derivations':0,'eigenfrequency_solves':0,'nonlinear_ODE_integrations':0,
            'linear_static_solves':2,'nonlinear_static_solves':2,
            'runtime_seconds':time.perf_counter()-started,
            'spatial_maxima_are_sampled':True,
            'strong_action_qualification':'existing strict float64 threshold retained; separately reported, not energy-equilibrium gate',
            'tangent_eigenvalue_scope':'positive finite-dimensional static Hessian; not a new eigenfrequency solution'}
    (bundle/'one_d_preflight.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    return result

"""Scoped FEM-2 static I/O fragment; no solver execution on import."""
from pathlib import Path
import re
import math
import os
import numpy as np

STATIC_CONTROL_DEFAULTS = {
    'initial_increment': 0.1, 'total_step_time': 1.0,
    'minimum_increment': 1e-6, 'maximum_increment': 0.1,
    'maximum_increments': 100, 'field_residual_relative': 1e-8,
    'field_correction_relative': 1e-8, 'output_frequency': 1000000,
}


def static_documentation_evidence():
    """Facts verified by reading installed CCX2.22 manual and archived source."""
    return {
        'manual': 'D:/PHD/CalculiX-Windows-master/src/downloads/ccx_2.22.pdf',
        'manual_sha256': '56963f827422ec7663cf218b60fffded19fd6ccebab793d2ccba667227d19d39',
        'version': '2.22',
        'GRAV': {'section': '7.43', 'pages': [470, 472, 474, 475],
                 'contract': 'Known global acceleration vector, load per unit reference mass; fixed direction 0,-1,0'},
        'STATIC': {'section': '7.122', 'pages': [594, 595, 596],
                   'contract': 'Loads ramp linearly over step; linear step length1; NL automatic increments'},
        'NLGEOM': {'section': '7.124', 'page': 600,
                   'contract': 'Green-Lagrange strain / internal PK2, printed stress Cauchy'},
        'ELASTIC': {'section': '6.8.1', 'page': 255,
                    'contract': 'Linear Green-Lagrange strain-to-PK2 elastic relation (StVK); smallstrain diagnostic'},
        'RF': {'sections': ['7.96', '7.98'], 'pages': [558, 562],
               'contract': 'RF = total nodal external force = support reaction + applied consistent nodal bodyload'},
        'FIELD': {'section': '7.24', 'pages': [440, 441],
                  'card': '1e-8,1e-8,,,1e-8,,1e-8,1e-8'},
        'final_output': {'section': '7.98', 'page': 562,
                         'contract': 'Final increment always written despite high FREQUENCY'},
        'source_archive': 'D:/PHD/CalculiX-Windows-master/src/downloads/ccx_2.22.src.tar.bz2',
        'source_archive_sha256': '3a94dcc775a31f570229734b341d6b06301ebdc759863df901c8b9bf1854c0bc',
        'source_evidence': {
            'e_c3d_rhs.f': 'xl=co reference geometry; bodyf=bodyfx*rho; ff += bodyf*N*referenceJ*quadrature',
            'rhs.f': 'GRAV bodyfx formed from specified magnitude and normalized global components',
            'resultsmech.f': 'fn = assembled internal nodal forces from stress',
            'printoutnode.f': 'RF prints fn; U and RF E13.6 double-based output',
            'frdvector.c': 'FRD ASCII vector casts float32 and E12.5 output',
            'frdheader.c': 'PSTEP = output_counter,iinc,istep; 100CL[12:24] = physical step time',
        },
        'reaction_recovery': 'R_support = RF_support - integral_reference rho*g*N; neither term inferred from global balance',
        'constitutive_comparison': 'StVK 3D law differs from adopted reduced quartic1D V4; no claim of literal constitutive equality',
    }


def write_static_input(input_path, mesh_include, mesh, audit, material, g,
                       nonlinear, settings=None):
    """Write paired STATIC decks, with STEP NLGEOM as sole lin/NL difference."""
    cfg = dict(STATIC_CONTROL_DEFAULTS)
    if settings: cfg.update(settings)
    input_path = Path(input_path)
    mesh_include = Path(mesh_include)
    if not math.isfinite(float(g)) or g <= 0:
        raise ValueError('Positive finite pre-frozen global GRAV magnitude required')
    if audit.get('status') != 'PASS':
        raise ValueError('Immutable source mesh must have passed mesh gate')
    ids = sorted(int(n) for n in mesh.nodes)
    left = audit['fixed_left_ids']; right = audit['fixed_right_ids']
    if not left or not right or set(left) & set(right):
        raise ValueError('Two disjoint nonempty full fixed faces required')
    def idlines(values):
        return [', '.join(str(int(n)) for n in values[k:k+16])
                for k in range(0, len(values), 16)]
    r = float(cfg['field_residual_relative']); c = float(cfg['field_correction_relative'])
    if r != 1e-8 or c != 1e-8:
        raise ValueError('Frozen force/correction control tolerances changed')
    lines = [
        '** FEM-2 static dead mass load. Reference rho, coordinates and full face clamps unchanged.',
        '** Linear/NL decks differ only by NLGEOM on STEP; no modal, dynamic or mid-span joint.',
        f'*INCLUDE, INPUT={os.path.relpath(mesh_include.resolve(),input_path.parent.resolve()).replace(chr(92),chr(47))}',
        '*NSET, NSET=ALL_NODES', *idlines(ids),
        '*NSET, NSET=LEFT_FIXED', *idlines(left),
        '*NSET, NSET=RIGHT_FIXED', *idlines(right),
        '*MATERIAL, NAME=MAT', '*ELASTIC',
        f"{fem1.single.ccx_float(float(material['E']))}, {fem1.single.ccx_float(float(material['nu']))}",
        '*DENSITY', fem1.single.ccx_float(float(material['rho'])),
        '*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT',
        '*BOUNDARY', 'LEFT_FIXED,1,3,0', 'RIGHT_FIXED,1,3,0',
        f"*STEP, INC={int(cfg['maximum_increments'])}" + (', NLGEOM' if nonlinear else ''),
        '*STATIC',
        ', '.join(fem1.single.ccx_float(float(cfg[k])) for k in
                  ('initial_increment','total_step_time','minimum_increment','maximum_increment')),
        '*CONTROLS, PARAMETERS=FIELD',
        ','.join((fem1.single.ccx_float(r),fem1.single.ccx_float(c),'','',fem1.single.ccx_float(r),'',fem1.single.ccx_float(c),fem1.single.ccx_float(r))),
        '*DLOAD', f'SOLID,GRAV,{fem1.single.ccx_float(float(g))},0,-1,0',
        f"*NODE FILE, GLOBAL=YES, FREQUENCY={int(cfg['output_frequency'])}", 'U,RF',
        f"*EL FILE, GLOBAL=YES, FREQUENCY={int(cfg['output_frequency'])}", 'S,E',
        f"*NODE PRINT, NSET=ALL_NODES, GLOBAL=YES, FREQUENCY={int(cfg['output_frequency'])}", 'U',
        f"*NODE PRINT, NSET=LEFT_FIXED, TOTALS=YES, GLOBAL=YES, FREQUENCY={int(cfg['output_frequency'])}", 'RF',
        f"*NODE PRINT, NSET=RIGHT_FIXED, TOTALS=YES, GLOBAL=YES, FREQUENCY={int(cfg['output_frequency'])}", 'RF',
        '*END STEP',
    ]
    input_path.write_text('\n'.join(lines)+'\n', encoding='utf8')
    return {'nonlinear': bool(nonlinear), 'settings': cfg,
            'GRAV_global_acceleration': [0.,-float(g),0.],
            'boundary': 'Ux=Uy=Uz=0 on complete x=0,L faces; lateralfaces free',
            'source_mesh_include': mesh_include.resolve().as_posix(),
            'constitutive': 'ELASTIC isotropic: StVK if NLGEOM, infinitesimal if linear',
            'printed_stress': 'Cauchy',
            'printed_strain': 'Green-Lagrange if nonlinear; infinitesimal if linear',
            'nodal_displacement_primary': 'DAT E13.6; FRD float32/E12.5 independent corroboration'}


def _static_frd_record(raw, ncomponent):
    if len(raw) < 13:
        raise ValueError('Incomplete fixed-width FRD nodal record')
    node = int(raw[3:13])
    payload = raw[13:].replace('D','E').replace('d','E')
    # Installed Windows build may use two or three exponent digits; a negative
    # mantissa may overrun nominal width. Node label must be parsed separately.
    rex = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+)[Ee][-+]\d{2,3}')
    tokens = rex.findall(payload)
    if len(tokens) != ncomponent or rex.sub('',payload).strip():
        raise ValueError('Invalid number/count in fixed-width FRD field')
    values = np.asarray([float(t) for t in tokens], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError('Nonfinite FRD nodal field')
    return node, values


def read_static_frd(path, expected_node_ids):
    """Read static datasets by actual STEP/increment/time, not modal ID or order."""
    blocks = []; current_step = None; time = None; field = None; count = None
    nodes = {}; collecting = False; declared_count = None
    def finish():
        nonlocal collecting, nodes, field
        if collecting:
            if not nodes: raise ValueError('Empty FRD static field')
            blocks.append({'step': current_step[2], 'increment': current_step[1],
                           'dataset_counter': current_step[0], 'time': time,
                           'name': field, 'nodes': nodes})
        collecting = False; nodes = {}; field = None
    for number, raw in enumerate(Path(path).read_text(encoding='utf8',errors='strict').splitlines(),1):
        s = raw.strip()
        if s.startswith('1PMODE') or 'MODAL' in raw[:75]:
            raise ValueError('Modal FRD supplied to static parser')
        if s.startswith('1PSTEP'):
            finish()
            try: current_step = tuple(int(t) for t in s.split()[-3:])
            except ValueError as exc: raise ValueError(f'Invalid static PSTEP at {number}') from exc
            if len(current_step) != 3: raise ValueError('Missing STEP/increment header')
        elif s.startswith('100C'):
            finish()
            try: time = float(raw[12:24]); declared_count = int(raw[24:36])
            except ValueError as exc: raise ValueError(f'Invalid static time header at {number}') from exc
            if not math.isfinite(time): raise ValueError('Nonfinite static dataset time')
        elif raw[:3].strip() == '-4':
            finish()
            parts = s.split(); field = parts[1]
            count = 3 if field in ('DISP','FORC') else 6 if field in ('STRESS','TOSTRAIN') else None
            if count is not None:
                if current_step is None or time is None: raise ValueError('Field lacks static time metadata')
                collecting = True; nodes = {}
        elif collecting and raw[:3].strip() == '-1':
            node, value = _static_frd_record(raw, count)
            if node in nodes: raise ValueError(f'Duplicate FRD node {node} in {field}')
            nodes[node] = value
        elif collecting and raw[:3].strip() == '-3': finish()
    finish()
    if not blocks: raise ValueError('Missing static FRD fields')
    latest = max((b['step'],b['increment'],b['time']) for b in blocks if b['name']=='DISP')
    final = {}
    expected = set(int(n) for n in expected_node_ids)
    for block in blocks:
        if (block['step'],block['increment'],block['time']) != latest: continue
        if block['name'] in final: raise ValueError('Duplicate final static field')
        if set(block['nodes']) != expected:
            raise ValueError(f"Incomplete {block['name']} nodes: expected {len(expected)}, got {len(block['nodes'])}")
        final[block['name']] = np.asarray([block['nodes'][int(n)] for n in expected_node_ids])
    missing = set(('DISP','FORC','STRESS','TOSTRAIN'))-set(final)
    if missing: raise ValueError(f'Missing final static FRD fields {sorted(missing)}')
    return final, {'step':latest[0], 'increment':latest[1], 'time':latest[2],
                   'datasets':len(blocks), 'nodal_count':len(expected),
                   'primary_vector_precision': 'float32 then E12.5 in installed source',
                   'tensor_order': ['XX','YY','ZZ','XY','YZ','ZX']}


def read_static_dat(path, expected_node_ids, left_ids, right_ids):
    blocks=[]; header=None; values={}
    rex = re.compile(r'^\s*(displacements|forces)\s+\([^)]*\)\s+for set\s+(\S+)\s+and time\s+(\S+)',re.I)
    def finish():
        nonlocal values
        if header is not None and values:
            blocks.append({'name':header[0], 'set':header[1], 'time':header[2], 'nodes':values})
        values={}
    for raw in Path(path).read_text(encoding='utf8',errors='strict').splitlines():
        found=rex.match(raw)
        if found:
            finish(); header=(found[1].lower(),found[2].upper(),float(found[3].replace('D','E')))
        elif header is not None:
            parts=raw.split()
            if len(parts)==4 and parts[0].isdigit():
                node=int(parts[0]); value=np.asarray([float(s.replace('D','E')) for s in parts[1:]])
                if node in values: raise ValueError('Duplicate node in static DAT block')
                if not np.all(np.isfinite(value)):raise ValueError('Nonfinite static DAT field')
                values[node]=value
    finish()
    requests=[('displacements','ALL_NODES',expected_node_ids),
              ('forces','LEFT_FIXED',left_ids), ('forces','RIGHT_FIXED',right_ids)]
    result={}; times={}
    for name,setname,ids in requests:
        candidates=[b for b in blocks if b['name']==name and b['set']==setname]
        if not candidates:raise ValueError(f'Missing static DAT {name} {setname}')
        last=max(candidates,key=lambda b:b['time'])
        if set(last['nodes'])!=set(int(n) for n in ids):
            raise ValueError(f'Incomplete static DAT {name} {setname}')
        result[setname]=np.asarray([last['nodes'][int(n)] for n in ids]); times[setname]=last['time']
    if len(set(times.values()))!=1:raise ValueError('DAT fields have differing final timestamps')
    return result, {'time':next(iter(times.values())), 'blocks':len(blocks),
                    'precision':'E13.6 double-based U/RF source (7significantdigits)',
                    'time_per_set':times}


def read_static_sta(path):
    rows=[]
    if not Path(path).exists():return {'status':'NOT_PRESENT','accepted_increments':[], 'attempt_rows':[]}
    for raw in Path(path).read_text(encoding='utf8',errors='strict').splitlines():
        p=raw.split()
        if len(p)==7 and all(re.fullmatch(r'\d+',v) for v in p[:4]):
            vals=[float(t.replace('D','E')) for t in p[4:]]
            rows.append({'step':int(p[0]),'increment':int(p[1]),'attempt':int(p[2]),
                         'iterations':int(p[3]),'total_time':vals[0],
                         'step_time':vals[1],'increment_time':vals[2]})
    accepted=[]
    for row in rows:
        if not accepted or row['total_time']>accepted[-1]['total_time']+1e-14:
            accepted.append(row)
    return {'status':'PARSED' if rows else 'NO_ACCEPTED_INCREMENTS',
            'accepted_increments':accepted,'attempt_rows':rows,
            'last_time':accepted[-1]['total_time'] if accepted else None}


def consistent_gravity_loads(mesh, rho, g):
    """Exact consistent C3D10 dead bodyload for affine reference tetrahedra.

    Integral Nv=-V/20 at 4vertices; Nm=V/5 at 6midsides. Negative corner
    consistent weights are normal for quadratic Lagrange tetrahedra. Values
    derive from simplex moments, not from measured reactions/global balance.
    """
    ids=np.asarray(sorted(int(n) for n in mesh.nodes),dtype=int)
    xyz=np.asarray([mesh.nodes[int(n)] for n in ids]); lookup={int(n):i for i,n in enumerate(ids)}
    f=np.zeros((len(ids),3)); volume=0.
    for eid in sorted(mesh.solid_elements):
        conn=np.asarray([lookup[int(n)] for n in mesh.solid_elements[eid]])
        if len(conn)!=10: raise ValueError('Only affine C3D10 supported')
        nodes=xyz[conn]; corner=nodes[:4]
        J=np.stack((corner[1]-corner[0],corner[2]-corner[0],corner[3]-corner[0]),axis=1)
        vol=float(np.linalg.det(J))/6.
        if vol<=0:raise ValueError('Invalid tetra orientation for consistent load')
        edges=((0,1),(1,2),(0,2),(0,3),(1,3),(2,3))
        mids=np.asarray([(corner[i]+corner[j])/2 for i,j in edges])
        if np.max(abs(nodes[4:]-mids))>1e-9:
            raise ValueError('Nonaffine C3D10 outside audited bodyload formula')
        values=vol*np.asarray([-1/20]*4+[1/5]*6)
        np.add.at(f[:,1],conn,-float(rho)*float(g)*values)
        volume+=vol
    return ids,xyz,f,volume


def parse_static_outputs(case_path, mesh, audit, rho, g, nonlinear,
                         equilibrium_relative_gate=1e-5):
    """Validate final datasets and independently recover support reactions.

    Return (JSON-compatible diagnostics, physical arrays). No equation solve.
    """
    stem=Path(case_path).with_suffix('')
    ids,xyz,body,volume=consistent_gravity_loads(mesh,rho,g)
    left=np.asarray(audit['fixed_left_ids'],int);right=np.asarray(audit['fixed_right_ids'],int)
    frd,fmeta=read_static_frd(stem.with_suffix('.frd'),ids)
    dat,dmeta=read_static_dat(stem.with_suffix('.dat'),ids,left,right)
    sta=read_static_sta(stem.with_suffix('.sta'))
    stdout=stem.with_suffix('.stdout.txt').read_text(encoding='utf8',errors='replace')
    stderr=stem.with_suffix('.stderr.txt').read_text(encoding='utf8',errors='replace')
    logs=stdout+'\n'+stderr
    warnings=[s.strip() for s in logs.splitlines() if '*WARNING' in s.upper() or '*ERROR' in s.upper()]
    failures=[]
    if warnings:failures.append('SOLVER_WARNING_OR_ERROR')
    if 'JOB FINISHED' not in logs.upper():failures.append('NO_REAL_JOB_FINISHED_MARKER')
    has_nlgeom='Nonlinear geometric effects are taken into account' in logs
    if bool(has_nlgeom)!=bool(nonlinear):failures.append('NLGEOM_ACTUAL_ROUTING_MISMATCH')
    if abs(fmeta['time']-1.)>1e-10 or abs(dmeta['time']-1.)>1e-10:
        failures.append('FINAL_LOAD_FACTOR_NOT_ONE')
    if nonlinear and (sta['status']!='PARSED' or abs(sta['last_time']-1.)>1e-10):
        failures.append('NO_COMPLETE_NL_LOAD_INCREMENT_HISTORY')
    index={int(n):i for i,n in enumerate(ids)}
    fixed=np.asarray([index[int(n)] for n in np.concatenate((left,right))]);free=np.setdiff1d(np.arange(len(ids)),fixed)
    support_ids=np.concatenate((left,right));RF=np.concatenate((dat['LEFT_FIXED'],dat['RIGHT_FIXED']))
    R=RF-body[fixed]
    F=np.sum(body,axis=0);Fmag=float(np.linalg.norm(F))
    if Fmag<=0:raise ValueError('No gravity load resultant')
    moment_coordinates=xyz+dat['ALL_NODES'] if nonlinear else xyz
    applied_moment=np.sum(np.cross(moment_coordinates,body),axis=0)
    support_moment=np.sum(np.cross(moment_coordinates[fixed],R),axis=0)
    balance=np.sum(R,axis=0)+F
    moment_balance=support_moment+applied_moment
    length=float(np.ptp(xyz[:,0]));reaction_relative=float(np.linalg.norm(balance)/Fmag)
    moment_relative=float(np.linalg.norm(moment_balance)/(Fmag*length))
    if reaction_relative>equilibrium_relative_gate:failures.append('SUPPORT_FORCE_BALANCE')
    if moment_relative>equilibrium_relative_gate:failures.append('SUPPORT_MOMENT_BALANCE')
    clamp=float(np.max(abs(dat['ALL_NODES'][fixed])))
    if clamp!=0.:failures.append('NONZERO_FIXED_FACE_DISPLACEMENT')
    delta=frd['DISP']-dat['ALL_NODES']; disp_scale=float(np.max(abs(dat['ALL_NODES'])))
    frd_dat_relative=float(np.max(abs(delta))/max(disp_scale,1e-30))
    if frd_dat_relative>1e-5:failures.append('FRD_DAT_DISPLACEMENT_DISAGREEMENT')
    free_rf_residual=float(np.max(abs(frd['FORC'][free]-body[free]))/Fmag)
    if free_rf_residual>equilibrium_relative_gate:failures.append('FREE_NODE_RF_BODYLOAD_IMBALANCE')
    support_frd_dat=float(np.max(abs(frd['FORC'][fixed]-RF))/Fmag)
    if support_frd_dat>equilibrium_relative_gate:failures.append('SUPPORT_FRD_DAT_RF_DISAGREEMENT')
    def support(which_ids,slc):
        force=R[slc];coords=xyz[[index[int(n)] for n in which_ids]]
        centroid=.5*(np.min(coords,axis=0)+np.max(coords,axis=0))
        return {'reaction_force_global':np.sum(force,axis=0).tolist(),
                'reaction_moment_about_face_centroid_global':np.sum(np.cross(coords-centroid,force),axis=0).tolist(),
                'raw_RF_global':np.sum(RF[slc],axis=0).tolist(),
                'consistent_body_load_global':np.sum(body[fixed][slc],axis=0).tolist()}
    diag={
        'status':'PASS' if not failures else 'FAIL','failures':failures,
        'final_frd':fmeta,'final_dat':dmeta,'load_history':sta,'solver_warnings':warnings,
        'actual_NLGEOM':has_nlgeom,'volume':volume,'mass':volume*float(rho),
        'applied_force_global':F.tolist(),'applied_moment_current_global':applied_moment.tolist(),
        'support_force_balance_global':balance.tolist(),'support_force_balance_relative':reaction_relative,
        'support_moment_balance_global':moment_balance.tolist(),'support_moment_balance_relative':moment_relative,
        'equilibrium_relative_gate':equilibrium_relative_gate,
        'maximum_fixed_displacement':clamp,'maximum_FRD_DAT_displacement_absolute':float(np.max(abs(delta))),
        'maximum_FRD_DAT_displacement_relative':frd_dat_relative,
        'maximum_free_RF_bodyload_residual_over_totalforce':free_rf_residual,
        'maximum_support_FRD_DAT_RF_disagreement_over_totalforce':support_frd_dat,
        'left':support(left,slice(0,len(left))), 'right':support(right,slice(len(left),None)),
        'reaction_contract': 'R=RF-consistent reference-bodyload at supports; independent force and current-moment balance',
        'strains': {'measure':'Green-Lagrange if NLGEOM, infinitesimal otherwise',
                    'maximum_absolute_nodal_recovered_components':np.max(abs(frd['TOSTRAIN']),axis=0).tolist(),
                    'maximum_absolute_nodal_recovered_strain':float(np.max(abs(frd['TOSTRAIN']))),
                    'qualification':'FRD stresses/strains extrapolated/averaged nodal fields, not exact point extrema'},
        'stresses': {'measure':'Cauchy', 'maximum_absolute_nodal_recovered_components':np.max(abs(frd['STRESS']),axis=0).tolist()},
        'nodal_displacement_primary':'DAT; same nodes as raw FRD',
        'printed_rounding': {'DAT_significant_digits':7,'FRD_significant_digits':6,
                             'DAT_pointwise_bound_conservative':0.5e-6*float(np.max(abs(dat['ALL_NODES']))),
                             'qualification':'Bound from formatting relative precision, not continuum discretization error'},
    }
    arrays={'node_ids':ids,'nodes':xyz,'U':dat['ALL_NODES'],'U_FRD':frd['DISP'],
            'RF_FRD':frd['FORC'],'support_ids':support_ids,'RF_support_DAT':RF,
            'R_support':R,'gravity_nodal':body,'S':frd['STRESS'],'E':frd['TOSTRAIN']}
    return diag,arrays


"""FEM-2 scoped reference-volume recovery; no solver or root calls."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation
from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem as fem1

FEM2_RECOVERY_POLICY = 'weighted_reference_slabs_cubic_axis_transverse_affine_polar_41'
FEM2_STATIC_FIELD_ORDER = ('u', 'w', 'v', 'Phi', 'psi', 'theta', 'c_eff')


def fem2_polar_section_orientation(directors):
    """Nearest orthonormal transverse directors and diagnostic 2x2 stretch."""
    a = np.asarray(directors, float)
    if a.shape != (3, 2) or not np.all(np.isfinite(a)):
        raise ValueError('Expected two finite deformed transverse directors')
    left, singular, right = np.linalg.svd(a, full_matrices=False)
    if singular[-1] <= np.finfo(float).eps * max(singular[0], 1.):
        raise ValueError('Degenerate deformed section directors')
    nk = left @ right
    rotation = np.column_stack((np.cross(nk[:, 0], nk[:, 1]), nk))
    if np.linalg.det(rotation) < 0.:
        raise ValueError('Section polar frame has wrong handedness')
    stretch = right.T @ np.diag(singular) @ right
    vector = Rotation.from_matrix(rotation).as_rotvec()
    return {'matrix': rotation, 'stretch': stretch,
            'theta': float(np.arctan2(-nk[0, 0], nk[1, 0])),
            'Phi': float(vector[0]), 'psi': float(-vector[1]),
            'orthogonality_error': float(np.max(np.abs(rotation.T @ rotation - np.eye(3))))}


def fem2_recover_reference_samples(global_points, global_displacement, mass_weights,
                                  length, thickness, width, section_count=41,
                                  enforce_clamped_faces=True):
    """Original material-x slabs; axis-cubic/transverse-affine weighted fit.

    [1,d,d²,d³,eta/h,zeta/b,d*eta/h,d*zeta/b], d=(x-xc)/(L/count).
    Intercept gives displacement at material (xc,0,0). Transverse columns
    of I+grad U give finite orientation. Only audited end FACE VALUES are
    supplied; this postprocessor adds no FE or derivative constraints.
    """
    xyz = fem1.nlsp_local_vectors(global_points).reshape(-1, 3)
    disp = fem1.nlsp_local_vectors(global_displacement).reshape(-1, 3)
    weight = np.asarray(mass_weights, float).reshape(-1)
    if xyz.shape != disp.shape or len(weight) != len(xyz):
        raise ValueError('Reference samples and displacement dimensions differ')
    if not np.all(np.isfinite(weight)) or np.any(weight <= 0):
        raise ValueError('Section recovery requires positive finite reference weights')
    if min(length, thickness, width) <= 0 or section_count < 4:
        raise ValueError('Invalid reference geometry or slab count')
    if xyz[:, 0].min() < -1e-10*length or xyz[:, 0].max() > (1+1e-10)*length:
        raise ValueError('Samples outside original rod length')
    bins = np.clip((xyz[:, 0]/length*section_count).astype(int), 0, section_count-1)
    rows, fields, small_angles, small_c, raw_averages = [], [], [], [], []
    total_mass = float(weight.sum())
    total_u2 = float(np.sum(weight*np.sum(disp*disp, axis=1)))
    residual_u2, max_condition = 0., 0.
    for k in range(section_count):
        selected = bins == k
        if np.count_nonzero(selected) < 8:
            raise ValueError(f'Insufficient volume samples in material slab {k}')
        p, u, wt = xyz[selected], disp[selected], weight[selected]
        xc = float(np.average(p[:, 0], weights=wt))
        span = length/section_count
        d, eta, zeta = (p[:, 0]-xc)/span, p[:, 1]/thickness, p[:, 2]/width
        design = np.column_stack((np.ones(len(d)), d, d*d, d*d*d,
                                  eta, zeta, d*eta, d*zeta))
        sw = np.sqrt(wt/wt.sum())
        fitted, _, rank, singular = np.linalg.lstsq(design*sw[:, None], u*sw[:, None], rcond=1e-12)
        if rank != 8:
            raise ValueError(f'Rank-deficient extraction slab {k}: {rank}/8')
        condition = float(singular[0]/singular[-1])
        max_condition = max(max_condition, condition)
        residual = u-design@fitted
        residual_sq = float(np.sum(wt*np.sum(residual*residual, axis=1)))
        residual_u2 += residual_sq
        gradient = np.column_stack((fitted[1]/span, fitted[4]/thickness, fitted[5]/width))
        deformation = np.eye(3)+gradient
        orient = fem2_polar_section_orientation(deformation[:, 1:3])
        small_theta = float(-gradient[0, 1])
        centroid = fitted[0]
        fields.append(np.array((centroid[0], centroid[1], centroid[2], orient['Phi'],
                                orient['psi'], orient['theta'], orient['stretch'][0, 0]-1.)))
        small_angles.append(small_theta)
        small_c.append(float(gradient[1, 1]))
        raw_averages.append(np.average(u, axis=0, weights=wt))
        rows.append({'x': xc, 'reference_mass': float(wt.sum()), 'samples': int(len(wt)),
                     'fit_rank': int(rank), 'fit_condition': condition,
                     'section_residual_L2': float(np.sqrt(residual_sq)),
                     'affine_gradient_local': gradient, 'finite_rotation_matrix': orient['matrix'],
                     'transverse_stretch': orient['stretch'],
                     'effective_width_strain': float(orient['stretch'][1, 1]-1.),
                     'polar_orthogonality_error': orient['orthogonality_error'],
                     'theta_finite_minus_small': float(orient['theta']-small_theta)})
    x, q = np.array([r['x'] for r in rows]), np.asarray(fields)
    angles, contraction, mean = np.asarray(small_angles), np.asarray(small_c), np.asarray(raw_averages)
    if enforce_clamped_faces:
        x = np.r_[0., x, length]
        q = np.vstack((np.zeros(7), q, np.zeros(7)))
        angles, contraction = np.r_[0., angles, 0.], np.r_[0., contraction, 0.]
        mean = np.vstack((np.zeros(3), mean, np.zeros(3)))
    return {'x': x, 'fields': q, 'small_rotation_theta': angles, 'small_contraction': contraction,
            'raw_mass_centroid_displacement': mean, 'section_rows': rows,
            'field_order': list(FEM2_STATIC_FIELD_ORDER), 'section_count': int(section_count),
            'recovery_policy': FEM2_RECOVERY_POLICY,
            'theta_policy': 'finite_3x2_polar_section_orientation_for_both_linear_and_nonlinear',
            'contraction_status': 'DIAGNOSTIC_EFFECTIVE_THICKNESS_STRETCH_NOT_MH_DOF',
            'reference_mass': total_mass, 'maximum_fit_condition': max_condition,
            'section_residual_mass_L2': float(np.sqrt(residual_u2)),
            'section_residual_relative_L2': float(np.sqrt(residual_u2/total_u2)) if total_u2 else 0.,
            'clamped_face_endpoint_values_used': bool(enforce_clamped_faces),
            'additional_derivative_constraints': False, 'reference_geometry_used_for_sections': True}


def fem2_static_sample(profile, x):
    """Same unsmoothed cubic interpolation on original material x for all cases."""
    xp, values = np.asarray(profile['x'], float), np.asarray(profile['fields'], float)
    if xp.ndim != 1 or np.any(np.diff(xp) <= 0) or values.shape != (len(xp), 7):
        raise ValueError('Invalid ordered static section profile')
    target = np.asarray(x, float)
    if target.min() < xp[0]-1e-12 or target.max() > xp[-1]+1e-12:
        raise ValueError('Static comparison requests extrapolation')
    return CubicSpline(xp, values, axis=0, extrapolate=False)(target)


def fem2_curve_difference(x, first, second, fixed_scale=None, characteristic_scale=None):
    """Sampled absolute L2/max; fixed full-profile scales, no alignment or gate."""
    grid, a, b = np.asarray(x, float), np.asarray(first, float), np.asarray(second, float)
    if a.shape != grid.shape or b.shape != grid.shape or np.any(np.diff(grid) <= 0):
        raise ValueError('Inconsistent scalar comparison profiles')
    difference = a-b
    where = int(np.argmax(np.abs(difference)))
    scale = float(max(np.max(np.abs(a)), np.max(np.abs(b)))) if characteristic_scale is None else float(characteristic_scale)
    maximum = float(np.max(np.abs(difference)))
    l2 = float(np.sqrt(np.trapezoid(difference*difference, grid)))
    return {'absolute_max': maximum, 'absolute_L2': l2, 'signed_at_max': float(difference[where]),
            'x_at_max': float(grid[where]), 'characteristic_scale': scale,
            'relative_max': maximum/scale if scale > 0 else None,
            'relative_L2': l2/(scale*np.sqrt(grid[-1]-grid[0])) if scale > 0 else None,
            'fixed_scale': float(fixed_scale) if fixed_scale is not None else None,
            'fixed_scale_max': maximum/fixed_scale if fixed_scale else None,
            'fixed_scale_L2': l2/(fixed_scale*np.sqrt(grid[-1]-grid[0])) if fixed_scale else None,
            'continuous_supremum_claimed': False, 'alignment_used': False}


def fem2_static_profile_symmetry(profile, length=1., count=801):
    x = np.linspace(0., length, count)
    sampled = fem2_static_sample(profile, x)
    parity = np.array((-1., 1., 1., 1., -1., -1., 1.))
    sym = sampled-parity[None, :]*sampled[::-1]
    result = {name: {'absolute_max': float(np.max(np.abs(sym[:, k]))),
                     'absolute_L2': float(np.sqrt(np.trapezoid(sym[:, k]**2, x)))}
              for k, name in enumerate(FEM2_STATIC_FIELD_ORDER)}
    result['out_of_plane'] = {name: float(np.max(np.abs(sampled[:, k])))
                              for name, k in (('v', 2), ('Phi', 3), ('psi', 4))}
    return result


def fem2_fe_strain_diagnostics(mesh, nodal_displacement):
    """Independent C3D10 reference gradients; NOT replacements for CCX E/S."""
    ids, xyz, eids, conn = fem1.mesh_arrays(mesh)
    nodal = np.asarray(nodal_displacement, float)
    if nodal.shape != (len(ids), 3) or not np.all(np.isfinite(nodal)):
        raise ValueError('Expected all finite mesh nodal displacement vectors')
    bary, _ = fem1.tet10_quadrature()
    _, dN = fem1.tet10_shape(bary)
    jac = np.einsum('eic,qij->eqcj', xyz[conn], dN)
    du = np.einsum('eic,qij->eqcj', nodal[conn], dN)
    gradient = du@np.linalg.inv(jac)
    linear = .5*(gradient+gradient.swapaxes(-1, -2))
    green = linear+.5*(gradient.swapaxes(-1, -2)@gradient)
    deformation = np.eye(3)+gradient
    sign = np.array((1., -1., -1.))
    linear_local = linear*sign[None, None, :, None]*sign[None, None, None, :]
    green_local = green*sign[None, None, :, None]*sign[None, None, None, :]
    eiglin, eiggreen = np.linalg.eigvalsh(linear), np.linalg.eigvalsh(green)
    return {'source': 'independent_C3D10_reference_gradient_at_14_positive_volume_points',
            'CCX_stress_strain_measures_replaced': False,
            'linear_strain_definition': '0.5*(gradU+gradU.T)',
            'green_lagrange_definition': 'linear_strain+0.5*gradU.T@gradU',
            'reference_nodes': int(len(ids)), 'elements': int(len(eids)),
            'sampled_gradient_max': float(np.max(np.abs(gradient))),
            'linear_strain_max_abs': float(np.max(np.abs(linear))),
            'green_lagrange_max_abs': float(np.max(np.abs(green))),
            'linear_principal_strain_max_abs': float(np.max(np.abs(eiglin))),
            'green_lagrange_principal_strain_max_abs': float(np.max(np.abs(eiggreen))),
            'linear_axial_max_abs': float(np.max(np.abs(linear_local[..., 0, 0]))),
            'green_lagrange_axial_max_abs': float(np.max(np.abs(green_local[..., 0, 0]))),
            'linear_engineering_inplane_shear_max_abs': float(np.max(np.abs(2*linear_local[..., 0, 1]))),
            'green_lagrange_engineering_inplane_shear_max_abs': float(np.max(np.abs(2*green_local[..., 0, 1]))),
            'minimum_det_deformation_gradient': float(np.min(np.linalg.det(deformation))),
            'maximum_det_deformation_gradient': float(np.max(np.linalg.det(deformation))),
            'finite_values': bool(np.all(np.isfinite(gradient)) and np.all(np.isfinite(green)))}


def recover_static_sections(mesh, nodal_displacement, coefficients=None, section_count=41,
                            rho=1., compare_recovery=True):
    quad = fem1.quadrature_arrays(mesh, rho)
    ids, xyz, _, _ = fem1.mesh_arrays(mesh)
    nodal = np.asarray(nodal_displacement, float)
    if nodal.shape != (len(ids), 3) or not np.all(np.isfinite(nodal)):
        raise ValueError('Incomplete static nodal displacement field')
    length, thickness, width = map(float, np.ptp(xyz, axis=0))
    disp = fem1.nlsp_evaluate_tet10_displacements(nodal, quad['conn'], quad['N'])
    result = fem2_recover_reference_samples(quad['xyz'], disp, quad['weights'],
                                          length, thickness, width, section_count)
    result['FE_strain_diagnostics'] = fem2_fe_strain_diagnostics(mesh, nodal)
    result['symmetry_diagnostics'] = fem2_static_profile_symmetry(result, length)
    if compare_recovery:
        alternate = fem2_recover_reference_samples(quad['xyz'], disp, quad['weights'],
                                                   length, thickness, width, 81)
        x = np.linspace(0., length, 801)
        first, second = fem2_static_sample(result, x), fem2_static_sample(alternate, x)
        result['recovery_sensitivity'] = {
            'policy': 'same_reference_recovery_41_vs_81_material_slabs',
            'primary_section_count': int(section_count), 'alternate_section_count': 81,
            'metrics': {name: fem2_curve_difference(x, first[:, k], second[:, k],
                                                   thickness if k < 3 else 1.)
                        for k, name in enumerate(FEM2_STATIC_FIELD_ORDER)},
            'alternate_x': alternate['x'], 'alternate_fields': alternate['fields']}
    return result


def compare_static_profiles(one_d_linear, one_d_nonlinear, fem_linear, fem_nonlinear,
                            thickness=.1, length=1., count=801, correction_scales=None):
    x = np.linspace(0., length, count)
    a, b, c, d = [fem2_static_sample(p, x) for p in
                   (one_d_linear, one_d_nonlinear, fem_linear, fem_nonlinear)]
    da, dc = b-a, d-c
    comparisons = {}
    for k, name in enumerate(FEM2_STATIC_FIELD_ORDER):
        fixed = thickness if k < 3 else 1.
        delta_scale = None if correction_scales is None else correction_scales[name]
        comparisons[name] = {
            'linear_1D_vs_3D': fem2_curve_difference(x, a[:, k], c[:, k], fixed),
            'nonlinear_1D_vs_3D': fem2_curve_difference(x, b[:, k], d[:, k], fixed),
            'nonlinear_correction_1D_vs_3D': fem2_curve_difference(x, da[:, k], dc[:, k], fixed, delta_scale),
            'one_D_nonlinear_effect': fem2_curve_difference(x, b[:, k], a[:, k], fixed,
                                                         float(np.max(np.abs(a[:, k])))),
            'FEM_nonlinear_effect': fem2_curve_difference(x, d[:, k], c[:, k], fixed,
                                                       float(np.max(np.abs(c[:, k])))),
            'correction_max_1D': float(np.max(np.abs(da[:, k]))),
            'correction_max_3D': float(np.max(np.abs(dc[:, k]))),
            'correction_signed_midpoint_1D': float(da[count//2, k]),
            'correction_signed_midpoint_3D': float(dc[count//2, k]),
            'diagnostic_only': name in ('c_eff', 'Phi', 'psi', 'v')}
    return {'x': x, 'one_D_linear': a, 'one_D_nonlinear': b,
            'FEM_linear': c, 'FEM_nonlinear': d,
            'one_D_correction': da, 'FEM_correction': dc, 'metrics': comparisons,
            'theta_qualification': '1D independent angle; 3D finite polar orientation; small-fit sensitivity retained',
            'same_material_x': True, 'amplitude_or_space_alignment': False,
            'sampling_count': int(count), 'frequency_mesh_criterion_applied': False}


def fem2_signal_resolution(signal, last_mesh_change, recovery_change=0., solver_uncertainty=0.):
    """Unresolved when observed uncertainty >= signal; no frequency threshold."""
    values = list(map(float, (signal, last_mesh_change, recovery_change, solver_uncertainty)))
    if any(v < 0 or not np.isfinite(v) for v in values):
        raise ValueError('Signal and uncertainty must be finite nonnegative')
    scale = max(values[1:])
    return {'signal': values[0], 'last_mesh_change': values[1],
            'recovery_change': values[2], 'solver_uncertainty': values[3],
            'largest_observed_uncertainty': scale,
            'signal_to_largest_observed_uncertainty': values[0]/scale if scale else None,
            'status': 'UNRESOLVED_SIGNAL' if values[0] <= scale else 'SIGNAL_EXCEEDS_OBSERVED_UNCERTAINTY',
            'continuum_error_bound_claimed': False, 'arbitrary_relative_threshold_used': False}

# ---- Bounded FEM-2 orchestration; all numerical physics lives in reused action. ----
from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem_refinement as fem1r
import argparse,csv,contextlib,importlib.metadata,shutil,subprocess,sys,time
FEM2_CONFIG=ROOT/'data/input/nlsp_nonlinear_static_3d_fem.json'
FEM2_OUTPUT=ROOT/'results/nlsp_nonlinear_static_3d_fem'
read_json,write_json,sha=fem1.read_json,fem1.write_json,fem1.sha
FEM2_STATUS_NAMES=('LOAD_PREFLIGHT','1D_LINEAR_STATIC','1D_NONLINEAR_STATIC','3D_LINEAR_STATIC','3D_NONLINEAR_STATIC','STATIC_SECTION_RECOVERY','STATIC_EQUILIBRIUM','NONLINEAR_CORRECTION_MESH_CHECK','1D_3D_COMPARISON')

def validate_fem2_config(c):
    if c['geometry']!={'L':1.,'b':.2,'h':.1} or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}:raise ValueError('Frozen FEM2 geometry/material changed')
    if c['mesh_levels']!=['medium','fine','refined'] or c['one_d']['p']!=[48,64]:raise ValueError('Bounded static meshes/degrees changed')
    if c['load_policy']!={'primary_w_over_h':.05,'backup_w_over_h':.03,'bending_surface_strain_ceiling':.01,'kind':'dead_global_gravity','global_direction':[0,-1,0],'line_load':'rho*A0*g'}:raise ValueError('Preselected load policy changed')
    if c['semantics']!={'static_only':True,'model_fitting':False,'new_meshes':False,'new_modal_jobs':False,'new_dynamics':False,'maximum_real_ccx_jobs':6,'automatic_further_studies':False}:raise ValueError('Unauthorized FEM2 policy')
    if c['threads']!=1 or c['job_timeout_seconds']>1200 or c['numerical_budget_seconds']>3600 or c['job_memory_limit_bytes']>4*1024**3:raise ValueError('Resource policy changed')
    if c['recovery']['policy']!=FEM2_RECOVERY_POLICY or c['recovery']['section_count']!=41 or c['recovery']['sensitivity_section_count']!=81:raise ValueError('Preselected recovery policy changed')
    return c

def load_fem2_sources(c):
    sources={}
    for name,row in c['sources'].items():
        b=ROOT/row['bundle']
        if sha(b/'manifest.json')!=row['manifest_sha256']:raise ValueError('Immutable source manifest changed: '+name)
        sources[name]=(b,fem1.validate_cache(b))
    old=sources['fem1'][1]; refined=sources['fem1r'][1]
    if old['preflight']['geometry']!=c['geometry'] or old['config']['material']!=c['material']:raise ValueError('Source geometry/material mismatch')
    if refined['config']['geometry']!=c['geometry'] or refined['config']['material']!=c['material']:raise ValueError('Refined source mismatch')
    for p,digest in old['config']['model_hashes'].items():
        if sha(ROOT/p)!=digest:raise ValueError('Frozen helper changed: '+p)
    if sha(c['ccx_exe'])!=read_json(sources['fem1'][0]/'manifest.json')['identity']['executables']['ccx_exe']:raise ValueError('Solver binary changed')
    action=ROOT/c['action_bundle']
    if sha(action/'manifest.json')!=c['action_manifest_sha256']:raise ValueError('Action manifest changed')
    am=read_json(action/'manifest.json')
    for p,d in am.get('artifact_hashes',am.get('artifacts',{})).items():
        if sha(action/p)!=d:raise ValueError('Action artifact corrupted: '+p)
    return sources

def fem2_identity(config_path=FEM2_CONFIG):
    c=validate_fem2_config(read_json(config_path));s=load_fem2_sources(c)
    coefficients=fem2_rod.RodCoefficients(**s['fem1'][1]['preflight']['coefficients'])
    selected=fem2_load_selection(c,coefficients)
    if selected['status']!='PASS':raise ValueError('No admitted load candidate')
    helpers={p:sha(ROOT/p) for p in ('scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_spatial_rod.py','scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py','scripts/analysis/verify_nlsp_linear_rectangular_3d_fem_refinement.py','scripts/analysis/solid_fem_single_rod_fixed_fixed.py')}
    docs=static_documentation_evidence()
    if sha(docs['manual'])!=docs['manual_sha256']:raise ValueError('Local CCX manual identity changed')
    item={'schema':c['schema'],'config':c,'config_sha256':sha(config_path),'code_sha256':sha(__file__),'selected_load':selected,
          'sources':c['sources'],'action_manifest_sha256':c['action_manifest_sha256'],'helper_sha256':helpers,
          'ccx_sha256':sha(c['ccx_exe']),'runtime_dlls':{p.name:sha(p) for p in sorted(Path(c['ccx_exe']).parent.glob('*.dll'))},
          'documentation':docs,'python':sys.version,'dependencies':{k:importlib.metadata.version(k) for k in ('numpy','scipy','matplotlib')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item

def load_one_d_static_profile(bundle,p=64,nonlinear=False):
    with np.load(bundle/f'one_d_p{p}.npz',allow_pickle=False) as data:
        values=data['nonlinear' if nonlinear else 'linear'];fields=np.zeros((len(values),7))
        fields[:,(0,1,5,6)]=values
        return {'x':data['s'].copy(),'fields':fields}

def fem2_source_mesh(c,level,sources):
    if level=='refined':
        source=sources['fem1r'][0]/'meshes/refined';record=sources['fem1r'][1]['refined_case']
    else:
        source=sources['fem1'][0]/'meshes'/level;record=sources['fem1'][1]['meshes'][level]
    audit=record['mesh_audit']
    if audit['status']!='PASS' or not audit['bbox_matches'] or audit['negative_or_zero_jacobian_elements'] or audit['solid_element_types']!=['C3D10']:raise ValueError('Source mesh gate failed')
    mesh=fem1.single.read_gmsh_inp_mesh_data(source/'rod.inp')
    if len(mesh.nodes)!=audit['nodes'] or len(mesh.solid_elements)!=audit['c3d10_elements']:raise ValueError('Source mesh count mismatch')
    return source,mesh,audit

def fem2_load_profile(path):
    return read_json(path)

def fem2_update_comparison(c,b,summary):
    completed=[n for n in c['mesh_levels'] if all(summary['cases'].get(n,{}).get(k,{}).get('status')=='PASS' for k in ('linear','nonlinear'))]
    if not completed:return
    one_lin=load_one_d_static_profile(b,64,False);one_nl=load_one_d_static_profile(b,64,True)
    count=c['recovery']['report_grid_count'];x=np.linspace(0,1,count)
    curves={'one_D_p64':fem2_static_sample(one_nl,x)-fem2_static_sample(one_lin,x)}
    pairs={}
    for n in completed:
        lin=fem2_load_profile(b/'cases'/n/'linear/recovered_sections.json')
        nl=fem2_load_profile(b/'cases'/n/'nonlinear/recovered_sections.json')
        pairs[n]=(lin,nl);curves[n]=fem2_static_sample(nl,x)-fem2_static_sample(lin,x)
    scales={name:max(float(np.max(abs(v[:,k]))) for v in curves.values()) for k,name in enumerate(FEM2_STATIC_FIELD_ORDER)}
    rows=[];comparisons={}
    for n,(lin,nl) in pairs.items():
        comparison=compare_static_profiles(one_lin,one_nl,lin,nl,.1,1,count,scales)
        comparisons[n]=comparison
        arrays={k:v for k,v in comparison.items() if isinstance(v,np.ndarray)}
        np.savez_compressed(b/f'comparison_{n}.npz',**arrays)
        write_json(b/f'comparison_{n}.json',{k:v for k,v in comparison.items() if not isinstance(v,np.ndarray)})
        w=comparison['metrics']['w']
        row={'mesh':n,'one_D_linear_midspan':float(comparison['one_D_linear'][count//2,1]),'one_D_NL_midspan':float(comparison['one_D_nonlinear'][count//2,1]),
             'FEM_linear_midspan':float(comparison['FEM_linear'][count//2,1]),'FEM_NL_midspan':float(comparison['FEM_nonlinear'][count//2,1]),
             'one_D_delta_w_midspan':w['correction_signed_midpoint_1D'],'FEM_delta_w_midspan':w['correction_signed_midpoint_3D'],
             'one_D_nonlinear_effect':w['correction_signed_midpoint_1D']/float(comparison['one_D_linear'][count//2,1]),
             'FEM_nonlinear_effect':w['correction_signed_midpoint_3D']/float(comparison['FEM_linear'][count//2,1]),
             'delta_w_difference_absolute_max':w['nonlinear_correction_1D_vs_3D']['absolute_max'],
             'delta_w_difference_absolute_L2':w['nonlinear_correction_1D_vs_3D']['absolute_L2'],
             'delta_w_difference_over_common_signal':w['nonlinear_correction_1D_vs_3D']['relative_max'],
             'correction_sign_agrees':bool(np.sign(w['correction_signed_midpoint_1D'])==np.sign(w['correction_signed_midpoint_3D']))}
        rows.append(row)
    mesh_changes=[]
    for first,second in zip(completed[:-1],completed[1:]):
        pair={'from':first,'to':second,'observables':{}}
        for k,name in enumerate(FEM2_STATIC_FIELD_ORDER):
            pair['observables'][name]={}
            for kind in ('linear','nonlinear','correction'):
                a=fem2_static_sample(pairs[first][0 if kind=='linear' else 1],x)[:,k] if kind!='correction' else curves[first][:,k]
                z=fem2_static_sample(pairs[second][0 if kind=='linear' else 1],x)[:,k] if kind!='correction' else curves[second][:,k]
                pair['observables'][name][kind]=fem2_curve_difference(x,a,z,.1 if k<3 else 1.,scales[name] if kind=='correction' else None)
        mesh_changes.append(pair)
    summary['comparison_rows']=rows;summary['mesh_changes']=mesh_changes;summary['correction_scales']=scales;summary['completed_levels']=completed
    write_json(b/'static_comparison.json',{'rows':rows,'mesh_changes':mesh_changes,'common_correction_scales':scales,'fixed_displacement_scale':.1,'alignments':False})
    with (b/'static_comparison.csv').open('w',encoding='utf8',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    all_done=completed==c['mesh_levels']
    if all_done:
        fine,refined=mesh_changes[-2:]
        latest=refined['observables']['w']['correction']['absolute_max'];previous=fine['observables']['w']['correction']['absolute_max']
        lp,np_=pairs['refined']
        altlin={'x':lp['recovery_sensitivity']['alternate_x'],'fields':lp['recovery_sensitivity']['alternate_fields']}
        altnl={'x':np_['recovery_sensitivity']['alternate_x'],'fields':np_['recovery_sensitivity']['alternate_fields']}
        recovery_difference=fem2_curve_difference(x,curves['refined'][:,1],(fem2_static_sample(altnl,x)-fem2_static_sample(altlin,x))[:,1],.1,scales['w'])
        rounding=sum(summary['cases']['refined'][k]['diagnostics']['printed_rounding']['DAT_pointwise_bound_conservative'] for k in ('linear','nonlinear'))
        signal=float(np.max(abs(curves['refined'][:,1])))
        resolution=fem2_signal_resolution(signal,latest,recovery_difference['absolute_max'],rounding)
        resolution.update(last_mesh_change_decreases=latest<=previous,sign_stable=all(np.sign(curves[n][count//2,1])==np.sign(curves['refined'][count//2,1]) for n in completed),
                          recovery_difference=recovery_difference,solver_uncertainty_qualification='Output-rounding bound only; actual solver convergence and force residuals checked separately, not a continuum trajectory bound')
        summary['signal_resolution']=resolution;write_json(b/'nonlinear_signal_resolution.json',resolution)
        accepted=resolution['status']=='SIGNAL_EXCEEDS_OBSERVED_UNCERTAINTY' and resolution['last_mesh_change_decreases'] and resolution['sign_stable']
        summary['statuses']['NLSP_FEM2_NONLINEAR_CORRECTION_MESH_CHECK']='PASS' if accepted else 'PARTIAL'
        summary['statuses']['NLSP_FEM2_1D_3D_COMPARISON']='PASS' if accepted else 'PARTIAL'

def run_fem2_cases(c,b,summary,sources,through_level='refined'):
    pre=summary['preflight'];levels=c['mesh_levels'][:c['mesh_levels'].index(through_level)+1]
    for level in levels:
        source,mesh,audit=fem2_source_mesh(c,level,sources)
        summary['cases'].setdefault(level,{})
        for kind,nonlinear in [('linear',False),('nonlinear',True)]:
            old=summary['cases'][level].get(kind)
            if old:
                if old['status']=='PASS':continue
                return finalize_fem2_summary(c,b,summary)
            case=b/'cases'/level/kind
            if case.exists() and any(case.iterdir()):raise RuntimeError('Unfinished attempted case exists; automatic repeat forbidden')
            case.mkdir(parents=True,exist_ok=True);inp=case/'static.inp';start=time.perf_counter()
            record={'status':'RUNNING','source_mesh':str(source.relative_to(ROOT)),'source_mesh_include_sha256':sha(source/'solid_mesh.inp'),'mesh_audit':audit}
            summary['cases'][level][kind]=record
            write_json(b/'summary.json',summary);write_json(case/'attempt_started.json',{'code_sha256':sha(__file__),'load':pre['load'],'target':'linear static' if not nonlinear else 'NLGEOM static'})
            try:
                write_json(case/'input_contract.json',write_static_input(inp,source/'solid_mesh.inp',mesh,audit,c['material'],pre['load']['g'],nonlinear,c['static_settings']))
                if summary['job_calls']['ccx']>=6:raise RuntimeError('Maximum six static jobs reached')
                remaining=c['numerical_budget_seconds']-summary['runtime']['numerical_seconds']
                if remaining<=0:raise TimeoutError('TOTAL_NUMERICAL_BUDGET')
                env=dict(os.environ);env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBER_OF_CPUS='1')
                print(f"FEM2 {level} {kind}: one static job, g={pre['load']['g']:.17g}",flush=True)
                summary['job_calls']['ccx']+=1
                result,stats=fem1.run_job([c['ccx_exe'],'static'],case,min(c['job_timeout_seconds'],remaining),c['job_memory_limit_bytes'],case/'static',env)
                record['job']=stats;write_json(case/'job.json',stats)
                if result.returncode or stats['failure']:raise RuntimeError('CCX_JOB_FAILED: '+str(stats))
                diag,arrays=parse_static_outputs(inp,mesh,audit,c['material']['rho'],pre['load']['g'],nonlinear,c['gates']['equilibrium_relative'])
                record['diagnostics']=diag;write_json(case/'static_diagnostics.json',diag);np.savez_compressed(case/'static_nodal_results.npz',**arrays)
                if diag['status']!='PASS':raise ArithmeticError('STATIC_OUTPUT_GATE: '+str(diag['failures']))
                recovered=recover_static_sections(mesh,arrays['U'],pre['coefficients'],41,c['material']['rho'],True)
                write_json(case/'recovered_sections.json',recovered)
                record['status']='PASS';record['recovery_policy']=FEM2_RECOVERY_POLICY
            except Exception as exc:
                record.update(status='FAIL',failure=str(exc))
            finally:
                seconds=time.perf_counter()-start;record['total_seconds']=seconds;summary['runtime']['numerical_seconds']+=seconds
                write_json(case/'case.json',record);write_json(b/'summary.json',summary)
            if record['status']!='PASS':return finalize_fem2_summary(c,b,summary)
        # Only a fully parsed/equilibrated medium pair admits further levels.
        if level=='medium':
            write_json(b/'medium_gate.json',{'status':'PASS','linear_and_NL_outputs_and_equilibrium_checked':True,'new_mesh_calls':0,'next_level_requires_gate':True})
        fem2_update_comparison(c,b,summary);write_json(b/'summary.json',summary)
    return finalize_fem2_summary(c,b,summary)

def finalize_fem2_summary(c,b,summary):
    for kind,status in [('linear','NLSP_FEM2_3D_LINEAR_STATIC'),('nonlinear','NLSP_FEM2_3D_NONLINEAR_STATIC')]:
        cases=[summary['cases'].get(n,{}).get(kind,{}) for n in c['mesh_levels']]
        summary['statuses'][status]='FAIL' if any(r.get('status')=='FAIL' for r in cases) else 'PASS' if all(r.get('status')=='PASS' for r in cases) else 'PARTIAL' if any(r.get('status')=='PASS' for r in cases) else 'NOT_RUN'
    all_cases=[v for level in summary['cases'].values() for v in level.values()]
    all_done=summary.get('completed_levels')==c['mesh_levels']
    diag_failed=any(v.get('diagnostics',{}).get('status')=='FAIL' for v in all_cases)
    diag_passed=any(v.get('diagnostics',{}).get('status')=='PASS' for v in all_cases)
    recovered=any('recovery_policy' in v for v in all_cases)
    recovery_failed=any(v.get('status')=='FAIL' and v.get('diagnostics',{}).get('status')=='PASS' for v in all_cases)
    summary['statuses']['NLSP_FEM2_STATIC_EQUILIBRIUM']='FAIL' if diag_failed else 'PASS' if all_done else 'PARTIAL' if diag_passed else 'NOT_RUN'
    summary['statuses']['NLSP_FEM2_STATIC_SECTION_RECOVERY']='FAIL' if recovery_failed else 'PASS' if all_done else 'PARTIAL' if recovered else 'NOT_RUN'
    summary['overall']='FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS' if all_done else 'PARTIAL'
    write_json(b/'summary.json',summary)
    return summary

def validate_fem2_cache(b,item=None):
    manifest=read_json(Path(b)/'manifest.json')
    if item is not None and item!=manifest['identity']:raise ValueError('FEM2 cache identity mismatch')
    for p,d in manifest['artifact_hashes'].items():
        if sha(Path(b)/p)!=d:raise ValueError('FEM2 artifact hash mismatch: '+p)
    c=manifest['identity']['config'];load_fem2_sources(c)
    return read_json(Path(b)/'summary.json')

def fem2_plot_only(b):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s=validate_fem2_cache(b);levels=s.get('completed_levels',[])
    if not levels:
        pre=s['preflight']
        if pre.get('status')!='PASS':return {'figures':0,'new_solver_calls':0}
        figures=Path(b)/'figures';figures.mkdir(exist_ok=True)
        fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained')
        for p in (48,64):
            with np.load(Path(b)/f'one_d_p{p}.npz') as data:
                for k,ax in enumerate(axes.flat):
                    ax.plot(data['s'],data['linear'][:,k],ls='--',label=f'1D p{p} linear')
                    ax.plot(data['s'],data['nonlinear'][:,k],label=f'1D p{p} quartic')
        for ax,name in zip(axes.flat,FEM2_FIELDS):
            ax.set(xlabel='Material x/L',ylabel=name);ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
        fig.suptitle('1D preflight only: no accepted 3D static solution')
        fig.savefig(figures/'one_d_static_preflight_profiles.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(figures/'one_d_static_preflight_profiles.png',dpi=220);plt.close(fig)
        fig,ax=plt.subplots(figsize=(7,3.8),layout='constrained')
        for p in (48,64):
            with np.load(Path(b)/f'one_d_p{p}.npz') as data:ax.plot(data['s'],data['delta_w'],label=f'1D p{p}')
        ax.set(title='1D preflight only',xlabel='Material x/L',ylabel='Quartic minus linear delta w')
        ax.grid(alpha=.2);ax.legend(frameon=False)
        fig.savefig(figures/'one_d_static_preflight_correction.pdf',metadata={'CreationDate':None,'ModDate':None});fig.savefig(figures/'one_d_static_preflight_correction.png',dpi=220);plt.close(fig)
        return {'figures':2,'new_solver_calls':0,'new_static_solutions':0,'scope':'1D preflight only'}
    figures=Path(b)/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    def finish(fig,name):
        fig.savefig(figures/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None});fig.savefig(figures/(name+'.png'),dpi=220);plt.close(fig)
    last=levels[-1]
    with np.load(Path(b)/f'comparison_{last}.npz') as d:
        x=d['x'];a=d['one_D_linear'];z=d['one_D_nonlinear'];fl=d['FEM_linear'];fn=d['FEM_nonlinear']
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    for values,label,style in [(a,'1D linear','--'),(z,'1D quartic','-'),(fl,'3D linear',':'),(fn,'3D NLGEOM','-.')]:
        axes[0].plot(x,values[:,1],style,label=label);axes[1].plot(x,values[:,0],style,label=label)
    for ax,label in zip(axes,('w','u')):ax.set(xlabel='Material x/L',ylabel=label);ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'linear_and_nonlinear_static_profiles')
    fig,ax=plt.subplots(figsize=(7,3.8),layout='constrained');ax.plot(x,z[:,1]-a[:,1],'k-',label='1D p64')
    for level in levels:
        with np.load(Path(b)/f'comparison_{level}.npz') as d:ax.plot(d['x'],d['FEM_correction'][:,1],label='3D '+level)
    ax.set(xlabel='Material x/L',ylabel='Nonlinear correction delta w');ax.grid(alpha=.2);ax.legend(frameon=False);finish(fig,'static_nonlinear_corrections')
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
    sizes={'medium':1/30,'fine':.025,'refined':.02};rows=s['comparison_rows'];axes[0].plot([sizes[r['mesh']] for r in rows],[r['FEM_delta_w_midspan'] for r in rows],'o-',label='3D actual')
    axes[0].axhline(rows[0]['one_D_delta_w_midspan'],ls='--',color='k',label='1D p64');axes[0].invert_xaxis();axes[0].set(xlabel='Target mesh size',ylabel='Midspan nonlinear correction');axes[0].legend(frameon=False)
    for kind in ('linear','nonlinear'):
        vals=[s['cases'][n][kind]['diagnostics']['support_force_balance_relative'] for n in levels]
        axes[1].plot([sizes[n] for n in levels],vals,'o-',label=kind)
    axes[1].set(xlabel='Target mesh size',ylabel='Relative support-force imbalance');axes[1].invert_xaxis();axes[1].legend(frameon=False)
    for ax in axes:ax.grid(alpha=.2)
    finish(fig,'static_mesh_convergence_and_reactions')
    return {'figures':3,'new_solver_calls':0,'new_static_solutions':0}


def failed_fem2_attempt(c, output_dir):
    """A technical code fix is not authorization to repeat a failed real job."""
    directories={Path(output_dir).resolve(),FEM2_OUTPUT.resolve()}
    for directory in directories:
        if not directory.is_dir():continue
        for b in sorted(directory.iterdir()):
            if not b.is_dir() or not (b/'manifest.json').exists() or not (b/'summary.json').exists():continue
            s=read_json(b/'summary.json');old=s.get('preflight',{})
            manifest=read_json(b/'manifest.json');cfg=manifest['identity']['config']
            same=all(cfg.get(k)==c.get(k) for k in ('geometry','material','sources','load_policy'))
            if same and any(v.get('status')=='FAIL' for level in s.get('cases',{}).values() for v in level.values()):
                return b,validate_fem2_cache(b)
    return None

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--preflight',action='store_true');mode.add_argument('--run-fem',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=FEM2_CONFIG);parser.add_argument('--output-dir',type=Path,default=FEM2_OUTPUT);parser.add_argument('--through-level',choices=['medium','fine','refined'],default='refined');a=parser.parse_args(argv)
    if a.report_only or a.plot_only:
        b=a.report_only or a.plot_only;s=validate_fem2_cache(b)
        if a.plot_only:fem2_plot_only(b)
        print(json.dumps({'bundle':str(b),'statuses':s['statuses'],'new_solver_static_BVP_calls':0},indent=2));return s
    c=validate_fem2_config(read_json(a.config))
    failed=failed_fem2_attempt(c,a.output_dir)
    if failed is not None:
        previous,s=failed
        print(json.dumps({'bundle':str(previous),'historical_failed_attempt_replay':True,'corrected_deck_not_execution_tested':True,'new_solver_static_BVP_calls':0,'statuses':s['statuses']},indent=2));return s
    key,item=fem2_identity(a.config);b=a.output_dir/key;sources=load_fem2_sources(c)
    if (b/'manifest.json').exists():
        s=validate_fem2_cache(b,item)
        if a.preflight or s.get('completed_levels')==c['mesh_levels'] or any(v.get('status')=='FAIL' for level in s['cases'].values() for v in level.values()):
            print(json.dumps({'bundle':str(b),'cache_hit':True,'new_solver_static_BVP_calls':0,'statuses':s['statuses']},indent=2));return s
    else:
        if b.exists() and any(b.iterdir()):raise RuntimeError('Partial unmanifested attempt exists; no automatic rerun')
        b.mkdir(parents=True,exist_ok=True);write_json(b/'provenance.json',item);write_json(b/'frozen_config.json',{**c,'selected_load':item['selected_load']})
        (b/'execution_code').mkdir();shutil.copyfile(__file__,b/'execution_code'/Path(__file__).name)
        write_json(b/'local_documentation_evidence.json',static_documentation_evidence())
        pre=build_static_preflight(c,b,ROOT)
        s={'preflight':pre,'cases':{},'job_calls':{'ccx':0,'gmsh':0,'modal':0,'nonlinear_ODE':0},'runtime':{'numerical_seconds':pre.get('runtime_seconds',0.),'budget_seconds':3600},
           'statuses':{'NLSP_FEM2_'+n:'NOT_RUN' for n in FEM2_STATUS_NAMES},'overall':'PREFLIGHT'}
        s['statuses'].update(NLSP_FEM2_LOAD_PREFLIGHT=pre['load']['status'],NLSP_FEM2_1D_LINEAR_STATIC='PASS' if pre.get('status')=='PASS' else 'FAIL',NLSP_FEM2_1D_NONLINEAR_STATIC='PASS' if pre.get('status')=='PASS' else 'FAIL')
        if pre.get('status')=='PASS' and pre['load']!=item['selected_load']:raise ValueError('1D computed load differs from pre-FEM freeze')
        write_json(b/'summary.json',s)
    if not a.preflight and s['preflight'].get('status')=='PASS':s=run_fem2_cases(c,b,s,sources,a.through_level)
    write_json(b/'manifest.json',fem1.artifact_manifest(b,item))
    if s.get('completed_levels'):fem2_plot_only(b)
    write_json(b/'manifest.json',fem1.artifact_manifest(b,item))
    print(json.dumps({'bundle':str(b),'statuses':s['statuses'],'runtime':s['runtime'],'job_calls':s['job_calls'],'overall':s['overall']},indent=2));return s

if __name__=='__main__':main()
