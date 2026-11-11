"""Local diagnostic of the first reflection-class crossing, never a map runner.

Script Proliferation Control: explicitly requested two-stage diagnostic. Stage A
independently audits saved classifications; stage B implements only the proved
3x3 boundary restrictions, reusing the frozen physical basis and mass metric.
No production matrices, reflection operator, pole catalogue or solver is changed.
The finite inertia queries certify the two local roots using saved catalogues.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import contextmanager, ExitStack
import csv
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import root, brentq
from scripts.analysis import classify_rl_full_rr_reflection_modes as saved

p = saved.PHYSICS
base = p.base
sections = saved.sections
CLASSIFIED = ROOT/'results/rl_full_rr_maps/_audit_reflection_classification_20261107'
OUTPUT = ROOT/'results/rl_full_rr_maps/_audit_reflection_intersection_20261107'
AUDIT = ROOT/'private/dissertation/research/chapter1/chapter1_reflection_symmetry_audit.md'
SELECTED = (0., 5., 6.1, 6.2, 10., 30., 60., 90.)
ER = np.vstack((np.zeros((3, 3)), np.eye(3)))


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def table(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stage_a(folder=CLASSIFIED):
    """All 1773 identities and all 394 first-class rows; no basis evaluation."""
    manifest = read(folder/'manifest.json')
    final = read(folder/'final_checks.json')
    require(saved.sha(folder/'manifest.json') == final['core_result_manifest_sha256'],
            'Classification manifest changed')
    hashes = {}
    for category in ('source_sha256', 'implementation_sha256', 'artifacts'):
        for path, expected in manifest[category].items():
            full = folder/path if category == 'artifacts' else ROOT/path
            require(saved.sha(full) == expected, 'Changed classification dependency: '+str(full))
            hashes[saved.local(full)] = expected
    for name in ('manifest.json', 'final_checks.json'):
        hashes[saved.local(folder/name)] = saved.sha(folder/name)
    with saved.no_spectral_calls() as forbidden:
        items, policy, inputs = saved.load_source(saved.DEFAULT_SOURCE, saved.DEFAULT_INVENTORY)
        hashes.update(inputs)
        classes = table(folder/'mode_classification.csv')
        spectra = table(folder/'class_spectrum.csv')
        index = {(i['beta'], i['root']['sorted_position']): i for i in items}
        require(len(classes) == len(spectra) == len(index) == 1773, 'Incomplete tables')
        independent = defaultdict(list)
        seen = set()
        for row in classes:
            key = (float(row['beta_deg']), int(row['sorted_position']))
            require(key not in seen, 'Duplicated classified mode'); seen.add(key)
            item = index[key]; r = item['root']; cls = row['symmetry_class']
            require(cls in ('+', '-') and row['classification_status'] ==
                    ('PURE_PLUS' if cls == '+' else 'PURE_MINUS'), 'Unresolved source class')
            require(float(row['Lambda']) == r['Lambda'] and float(row['omega']) == r['omega']
                    and int(row['multiplicity']) == r['multiplicity'], 'Source frequency mismatch')
            require(row['source_npz'] == item['npz'] and row['source_shape_key'] == item['shape_key']
                    and row['source_point_id'] == item['point_id'], 'Mode provenance mismatch')
            independent[(key[0], cls)].append((r['Lambda'], key[1], item, row))
        actual = {(float(r['beta_deg']), r['symmetry_class'], int(r['class_rank'])): r for r in spectra}
        require(len(actual) == 1773, 'Duplicate class rank')
        first, rows, union = {}, [], defaultdict(list)
        for (beta, cls), values in sorted(independent.items()):
            values.sort(key=lambda v: (v[0], v[1]))
            for rank, (lam, position, item, row) in enumerate(values, 1):
                got = actual[(beta, cls, rank)]
                require(float(got['Lambda']) == lam and int(got['original_sorted_position']) == position
                        and got['source_npz'] == item['npz'] and got['shape_key'] == item['shape_key'],
                        'Independent class sorting disagrees')
                union[beta].append((lam, position))
            lam, position, item, row = values[0]
            norm = saved.mass(item['model'], item['mode'])
            require(abs(norm-1) <= policy['mass_gram_tol'], 'Unnormalized saved first mode')
            defect = saved.mass(item['model'], saved.combine(
                sections.reflect_mode(item['mode']), item['mode'], 1, -1 if cls == '+' else 1))**.5/norm**.5
            require(defect <= 1e-9, 'Saved first mode has wrong reflection class')
            first[(beta, cls)] = item
            rows.append(dict(beta_deg=beta, symmetry_class=cls, class_rank=1, Lambda=lam,
                omega=item['root']['omega'], original_sorted_position=position,
                source_point_id=item['point_id'], shape_key=item['shape_key'], source_npz=item['npz'],
                mass=norm, reflection_defect=defect, selected_physical_inspection=beta in SELECTED,
                independent_rank_matches=True, source_csv_matches=True))
        for beta, values in union.items():
            require(sorted(values) == [(index[(beta,k)]['root']['Lambda'], k) for k in range(1,10)],
                    'Class union does not recover complete nine-frequency prefix')
        for row in table(folder/'class_first_gaps.csv'):
            beta = float(row['beta_deg'])
            plus, minus = (first[(beta,c)]['root']['Lambda'] for c in ('+', '-'))
            require(float(row['Lambda_plus']) == plus and float(row['Lambda_minus']) == minus
                    and float(row['gap']) == plus-minus, 'Renderer input differs from independent minima')
        require(len(union) == 197 and len(first) == 394, 'Incomplete first-class coverage')
        summary = {}
        for cls in ('+', '-'):
            chosen = [r for r in rows if r['symmetry_class'] == cls]
            vals = [r['Lambda'] for r in chosen]
            identities = {(r['source_npz'], r['shape_key']) for r in chosen}
            require(len(identities) == 197, 'Same physical shape reused at different angles')
            summary[cls] = dict(min=min(vals), max=max(vals), range=max(vals)-min(vals),
                relative_range=(max(vals)-min(vals))/min(vals), distinct_float_values=len(set(vals)),
                distinct_source_identities=len(identities),
                status='BRANCH_VARIATION_CONFIRMED' if len(set(vals)) > 1 else 'UNRESOLVED')
        checks = dict(status='PASS', angles=197, forms=1773, first_class_rows=394,
            independent_sorting_matches=True, nine_root_unions_match=True,
            normalization_max_error=max(abs(r['mass']-1) for r in rows),
            renderer_uses_verified_class_first_gaps=True, summary=summary,
            spectral_or_basis_calls=len(forbidden))
    return rows, checks, first, policy, hashes


def restrictions(beta):
    """C+ and C- exactly as analytical audit section 7.3, state U,W,Psi,N,Q,M."""
    c, s = np.cos(np.deg2rad(beta)/2), np.sin(np.deg2rad(beta)/2)
    plus = np.array([[c,-s,0,0,0,0], [0,0,0,s,c,0], [0,0,1,0,0,0]])
    minus = np.array([[s,c,0,0,0,0], [0,0,0,c,-s,0], [0,0,0,0,0,1]])
    return plus, minus


def row_transform(beta):
    c, s = np.cos(np.deg2rad(beta)/2), np.sin(np.deg2rad(beta)/2)
    return .5*np.array([[c,s,0,0,0,0], [0,0,s,-c,0,0], [0,0,0,0,1,0],
                       [s,-c,0,0,0,0], [0,0,c,s,0,0], [0,0,0,0,0,1]])


def frequency(omega, time_reference=200.):
    require(omega > 0 and time_reference > 0, 'Positive frequency/reference required')
    return float(np.sqrt(omega*time_reference))


class Blocks:
    """Finite-domain diagnostic adapter to the unchanged homogeneous arm basis."""
    def __init__(self, item):
        self.model, self.catalogs, self.policy = item['model'], item['record']['catalogs'], item['record']['contract']['policy']
        self.reference = item['record']['contract']['time_reference']
        self.case_template = item['case']
        self.arm, self.length = self.model.arms[0], self.model.l
        counter = sections.PhysicalCounter(self.model, self.case(6.1), self.catalogs, self.policy)
        self.limits = counter.limits; self.domain = counter.domain
        self.basis_calls = 0; self.count_calls = 0; self.evaluations = 0
        self.regularity = []

    def case(self, beta):
        return p.Case(beta=float(beta), mu=0., tau=0., k_joint='RIGID', k_end='RIGID')

    def basis(self, omega, x, derivative=0):
        require(0 < omega < self.domain, 'Frequency outside proved basis/catalogue domain')
        self.basis_calls += 1
        return base.basis(self.arm, self.length, omega, x, derivative)

    def transfer(self, omega):
        z0, zl = self.basis(omega, np.array([0., self.length]))
        zb, row, col = base.balance(z0)
        sv = np.linalg.svd(zb, compute_uv=False)
        require(sv[-1]/sv[0] > self.policy['nullity_tol'], 'Z(0) lost numerical rank')
        rhs = row[:,None]*ER
        solved = np.linalg.solve(zb, rhs)
        s0 = col[:,None]*solved
        normalized = s0/np.linalg.norm(s0, axis=0)
        ss = np.linalg.svd(normalized, compute_uv=False)
        require(ss[-1]/ss[0] > self.policy['nullity_tol'], 'S0 lost column rank')
        error = float(np.max(abs(z0@s0-ER)))
        # Unit reactions produce large, dimensioned columns. Test the backward
        # error of the balanced linear solve, not absolute mixed-unit entries.
        backward = float(np.linalg.norm(zb@solved-rhs) /
                         (np.linalg.norm(zb)*np.linalg.norm(solved)+np.linalg.norm(rhs)))
        require(backward < 64*np.finfo(float).eps and sv[0]/sv[-1] < self.policy['condition_max'],
                'Clamped transfer solve inconsistent')
        self.regularity.append(dict(omega=float(omega), Z0_balanced_condition=float(sv[0]/sv[-1]),
            S0_column_rank=3, S0_scaled_condition=float(ss[0]/ss[-1]), Z0S0_error=error,
            balanced_backward_error=backward))
        return s0, zl@s0

    def pair(self, beta, omega):
        self.evaluations += 1
        require(0 <= beta <= 90, 'Diagnostic angle outside archive interval')
        s0, tl = self.transfer(float(omega))
        return [base.balance(c@tl) for c in restrictions(beta)], s0, tl

    def determinants(self, beta, omega):
        pair, _, _ = self.pair(beta, omega)
        return np.array([np.linalg.det(v[0]) for v in pair])

    def checks(self, beta, omega):
        pair, s0, tl = self.pair(beta, omega)
        data = []
        for balanced, row, col in pair:
            s = np.linalg.svd(balanced, compute_uv=False)
            data.append(dict(det=float(np.linalg.det(balanced)), singular_values=s.tolist(),
                relative_singular_values=(s/s[0]).tolist(),
                nullity=int(sum(s <= s[0]*self.policy['nullity_tol'])),
                row_scales=row.tolist(), column_scales=col.tolist()))
        return data

    def count(self, beta, omega):
        self.count_calls += 1
        actual, n, details = sections.PhysicalCounter(self.model, self.case(beta), self.catalogs, self.policy)(omega)
        require(actual == omega, 'Inertia verification frequency shifted at a pole')
        return dict(beta_deg=float(beta), omega=float(omega), count=n, details=details)


@contextmanager
def local_only():
    """Block full boundary solvers, catalogue construction, Ritz and map runners."""
    names = ('prefix', 'counted_prefix', 'catalogues', 'ritz', 'boundary', 'recovered_modes', 'run', 'solve')
    calls = []
    def forbidden(*args, **kwargs):
        calls.append('prohibited full-spectrum/Ritz/catalogue call')
        raise AssertionError(calls[-1])
    with ExitStack() as stack:
        for name, module in list(sys.modules.items()):
            if name.startswith('scripts.') and module is not sys.modules[__name__]:
                for attr in names:
                    if callable(getattr(module, attr, None)):
                        stack.enter_context(patch.object(module, attr, forbidden))
        yield calls


def restore(engine, beta, omega, parity, nodes=None, weights=None):
    pair, s0, _ = engine.pair(beta, omega)
    balanced, _, col = pair[0 if parity == 1 else 1]
    _, singular, vh = np.linalg.svd(balanced)
    require(sum(singular <= singular[0]*engine.policy['nullity_tol']) == 1,
            'Restoration requires a one-dimensional class kernel')
    reaction = col*vh[-1]
    coefficients = s0@reaction
    if nodes is None:
        z, wt = np.polynomial.legendre.leggauss(120)
        nodes, weights = (z+1)*engine.length/2, wt*engine.length/2
    arm = dict(x=nodes.copy(), weights=weights.copy(), y=engine.basis(omega,nodes)@coefficients,
        dy=engine.basis(omega,nodes,1)@coefficients,
        ends=engine.basis(omega,np.array([0.,engine.length]))@coefficients)
    # Use the existing physical reflection, never a second state operator.
    reflected = sections.reflect_mode(dict(arms=[arm, arm]))['arms'][0]
    mode = dict(arms=[arm, dict(reflected, **{k: parity*reflected[k] for k in saved.FIELDS})])
    norm = saved.mass(engine.model, mode)
    require(norm > 0 and np.isfinite(norm), 'Nonpositive restored mode mass')
    mode = saved.combine(mode, mode, 1/np.sqrt(norm), 0)
    return mode, coefficients/np.sqrt(norm), reaction/np.sqrt(norm)


def local_roots(engine, first):
    """Two independent local procedures; no frequency scan or interpolation."""
    endpoints = [first[(b,c)]['root']['omega'] for b in (6.1,6.2) for c in ('+','-')]
    low, high = min(endpoints)*.995, max(endpoints)*1.005
    scale = np.array([.1, high-low]); origin = np.array([6.15, np.mean(endpoints)])
    def equations(v):
        beta, omega = origin+scale*v
        require(6.09 < beta < 6.21 and low < omega < high, 'Joint iteration left local box')
        return engine.determinants(beta, omega)
    simultaneous = []
    for start in (6.11,6.15,6.19):
        for tolerance in (1e-9, 1e-11, 1e-12):
            guess = (np.array([start, origin[1]])-origin)/scale
            ans = root(equations, guess, method='hybr', options=dict(xtol=tolerance,maxfev=100))
            beta, omega = origin+scale*ans.x
            residual = engine.determinants(beta,omega)
            require(np.max(abs(residual)) < engine.policy['nullity_tol'], 'Joint solution residual failed')
            require(6.1 < beta < 6.2, 'Joint solution outside archive bracket')
            simultaneous.append(dict(start_beta=start, xtol=tolerance, beta_deg=float(beta),
                omega=float(omega), Lambda=frequency(omega,engine.reference),
                determinants=residual.tolist(), scipy_success=bool(ans.success), message=str(ans.message), nfev=ans.nfev))
    def each(beta, index, tolerance):
        function = lambda w: engine.determinants(beta,w)[index]
        require(function(low)*function(high) < 0, 'No local class-root bracket')
        return brentq(function,low,high,xtol=tolerance,rtol=8*np.finfo(float).eps,maxiter=80)
    nested = []
    for wt,bt in ((2e-12,1e-9),(2e-14,2e-11),(5e-16,2e-12)):
        def difference(beta):
            return each(beta,0,wt)-each(beta,1,wt)
        beta = brentq(difference,6.1,6.2,xtol=bt,rtol=8*np.finfo(float).eps,maxiter=60)
        wp,wm = each(beta,0,wt), each(beta,1,wt)
        nested.append(dict(omega_xtol=wt,beta_xtol=bt,beta_deg=float(beta),
            omega_plus=wp,omega_minus=wm,omega_gap=wp-wm))
    best = simultaneous[-1]
    beta, omega = best['beta_deg'], best['omega']
    rows, certificates = [], []
    for angle in sorted({6.1,6.125,6.15,6.175,6.2,beta}):
        for side, bound in (('lower',low),('upper',high)):
            cert = engine.count(angle,bound); cert['bound']=side; certificates.append(cert)
            require(cert['count'] == (0 if side=='lower' else 2), 'Local interval is not the first two roots')
        for index,cls in enumerate(('+','-')):
            w = each(angle,index,5e-16)
            diag = engine.checks(angle,w)[index]
            require(diag['nullity'] == 1, 'Local class kernel not simple')
            rows.append(dict(beta_deg=angle,symmetry_class=cls,class_rank=1,omega=w,
                Lambda=frequency(w,engine.reference),determinant=diag['det'],
                nullity=diag['nullity'],first_class_certified=True,
                source='local_3x3_brentq',lower_frequency=low,upper_frequency=high))
    for sign in (-1,1):
        cert=engine.count(beta,omega*(1+sign*1e-6)); cert['bound']='crossing_left' if sign<0 else 'crossing_right'
        require(cert['count'] == (0 if sign<0 else 2), 'Crossing multiplicity count failed')
        certificates.append(cert)
    require(max(abs(r['beta_deg']-beta) for r in simultaneous+nested) < 1e-8, 'Angle not stable')
    require(max(abs(r['omega']-omega) for r in simultaneous) < 1e-12, 'Frequency not stable')
    return rows, dict(beta_deg=beta,omega=omega,Lambda=frequency(omega,engine.reference),
        Omega=omega*engine.reference,time_reference=engine.reference,
        simultaneous=simultaneous,nested=nested,local_frequency_bracket=[low,high],
        count_certificates=certificates,block_checks=engine.checks(beta,omega),
        first_class_roots_certified=True)


def form_checks(engine, beta, omega, mode, parity):
    result = saved.physical_checks(engine.model,engine.case(beta),omega,mode)
    require(saved.physical_pass(result,engine.policy), 'Restored physical mode failed')
    reflected = sections.reflect_mode(mode)
    diff = saved.combine(reflected,mode,1,-parity)
    result['reflection_defect'] = float(np.sqrt(saved.mass(engine.model,diff)/saved.mass(engine.model,mode)))
    require(result['reflection_defect'] < 1e-9, 'Wrong restored symmetry class')
    result['max_displacement'] = float(max(np.max(abs(a['y'][:,:3])) for a in mode['arms']))
    require(result['max_displacement'] > 0, 'Zero restored mode')
    return result


def stage_b(first, policy):
    engine = Blocks(first[(6.1,'+')])
    with local_only() as prohibited:
        roots, solution = local_roots(engine,first)
        beta,omega = solution['beta_deg'],solution['omega']
        modes, checks, arrays = {}, {}, {}
        for cls,parity in (('+',1),('-',-1)):
            template = first[(6.1,cls)]['mode']['arms'][0]
            mode,coefficient,reaction = restore(engine,beta,omega,parity,template['x'],template['weights'])
            modes[cls] = mode; checks[cls] = form_checks(engine,beta,omega,mode,parity)
            label='plus' if parity==1 else 'minus'
            for i,arm in enumerate(mode['arms'],1):
                for k,v in arm.items(): arrays[f'{label}_arm{i}_{k}']=v
            arrays[f'{label}_arm1_coefficients']=coefficient
            arrays[f'{label}_arm1_reaction']=reaction
        gram=np.array([[saved.mass_inner(engine.model,modes[a],modes[b]) for b in ('+','-')] for a in ('+','-')])
        require(np.max(abs(gram-np.eye(2))) < policy['mass_gram_tol'], 'Crossing forms not mass orthonormal')
        checks['mass_gram']=gram.tolist(); checks['full_eigenspace_dimension']=2
        checks['block_nullities']=[v['nullity'] for v in solution['block_checks']]
        require(checks['block_nullities']==[1,1], 'Unexpected internal multiplicity')
        comparison=[]
        for angle in (6.1,6.2):
            for cls,parity in (('+',1),('-',-1)):
                item=first[(angle,cls)]; original=item['mode']
                rr=next(r for r in roots if r['beta_deg']==angle and r['symmetry_class']==cls)
                endpoint,_,_=restore(engine,angle,rr['omega'],parity,
                    original['arms'][0]['x'],original['arms'][0]['weights'])
                form_checks(engine,angle,rr['omega'],endpoint,parity)
                matched=saved.mass_inner(engine.model,endpoint,original)
                overlap=sections.global_inner(engine.model,modes[cls],original,beta,angle)
                err=abs(rr['omega']/item['root']['omega']-1)
                require(err < policy['same_physics_frequency_tol'],
                        'Endpoint frequency disagrees with archive')
                require(abs(abs(matched)-1) < policy['mass_gram_tol'], 'Endpoint mode differs from saved mode')
                comparison.append(dict(beta_deg=angle,symmetry_class=cls,
                    original_sorted_position=item['root']['sorted_position'],Lambda_archive=item['root']['Lambda'],
                    Lambda_block=rr['Lambda'],omega_relative_difference=err,
                    same_angle_mass_MAC=float(abs(matched)**2),
                    crossing_to_archive_global_mass_overlap=float(overlap),
                    crossing_to_archive_global_mass_MAC=float(abs(overlap)**2),
                    source_point_id=item['point_id'],source_npz=item['npz'],shape_key=item['shape_key']))
        # No searches at these angles: substitute all 16 selected stored first frequencies.
        substitutions=[]
        for angle in SELECTED:
            for k,cls in enumerate(('+','-')):
                item=first[(angle,cls)]; d=engine.checks(angle,item['root']['omega'])[k]
                require(d['relative_singular_values'][-1] < policy['nullity_tol'], 'Saved frequency does not satisfy class block')
                substitutions.append(dict(beta_deg=angle,symmetry_class=cls,**d))
        solution['archive_block_substitutions']=substitutions
        solution['regularity']=dict(limits=engine.limits,
            max_Z0_balanced_condition=max(r['Z0_balanced_condition'] for r in engine.regularity),
            max_S0_scaled_condition=max(r['S0_scaled_condition'] for r in engine.regularity),
            max_Z0S0_error=max(r['Z0S0_error'] for r in engine.regularity),
            max_balanced_backward_error=max(r['balanced_backward_error'] for r in engine.regularity),
            min_omega=min(r['omega'] for r in engine.regularity),max_omega=max(r['omega'] for r in engine.regularity),
            S0_shape=[6,3],S0_column_rank=3)
        solution['calls']=dict(basis=engine.basis_calls,finite_inertia_queries=engine.count_calls,
            block_pair_evaluations=engine.evaluations,prohibited=len(prohibited),full_spectrum=0,ritz=0)
        solution['calls']['basis_count_scope']='Blocks.basis adapter only; excludes internal saved-catalogue inertia-query calls'
        solution['calls']['basis_count_scope']='Blocks.basis adapter only; excludes internal saved-catalogue inertia-query calls'
        arrays.update(beta_deg=np.array(beta),omega=np.array(omega),Lambda=np.array(solution['Lambda']),
                      state_order=np.array(['U','W','Psi','N','Q','M']))
    return roots,solution,modes,checks,arrays,comparison


def render(output, first):
    """Plot saved archive/local results only; no basis, count or root calls."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import findfont, FontProperties
    solution=read(output/'crossing_solution.json')
    rows=table(output/'branch_constancy_audit.csv')
    local=table(output/'block_roots.csv')
    modes={}
    with np.load(output/'crossing_forms.npz',allow_pickle=False) as a:
        for cls,label in (('+','plus'),('-','minus')):
            modes[cls]=dict(arms=[{k:a[f'{label}_arm{i}_{k}'].copy() for k in ('x','weights',*saved.FIELDS)}
                                  for i in (1,2)])
    color={'+':'#1f77b4','-':'#ff7f0e'}
    captions=[]
    with plt.rc_context({'font.family':'serif','font.serif':['Times New Roman'],
        'font.size':12,'mathtext.fontset':'stix','pdf.fonttype':42,
        'axes.grid':False,'figure.facecolor':'white','axes.facecolor':'white','path.simplify':False}):
        def save(fig,name,caption):
            fig.text(.52,.10,r'$\mu=\tau=0,\quad a=1,\quad\varepsilon=0.005$',ha='center')
            fig.text(.52,.025,caption,ha='center',va='bottom',fontsize=10)
            for ext in ('png','pdf'):
                path=output/(name+'.'+ext)
                if path.exists(): raise FileExistsError(path)
                fig.savefig(path,dpi=300,facecolor='white')
            captions.append(dict(name=name,caption=caption));plt.close(fig)
        def style(ax):
            ax.tick_params(direction='in',top=True,right=True,length=3,width=.7)
            for spine in ax.spines.values():spine.set_linewidth(.7)
            ax.grid(False)
        fig,ax=plt.subplots(figsize=(7.2,5.6));fig.subplots_adjust(left=.12,right=.97,top=.97,bottom=.31)
        for cls in ('+','-'):
            data=[r for r in rows if r['symmetry_class']==cls]
            x=[float(r['beta_deg']) for r in data];y=[float(r['Lambda']) for r in data]
            ax.plot(x,y,color=color[cls],lw=1.3,label=rf'$\Lambda_{{{cls},1}}$')
        ax.set(xlim=(0,90),xlabel=r'$\beta$, градусы',ylabel=r'$\Lambda$');style(ax)
        ax.legend(frameon=False,loc='lower right')
        # Rug marks make the actual 197 archive nodes visible without markers on curves.
        angles=sorted({float(r['beta_deg']) for r in rows})
        ax.plot(angles,np.full(len(angles),.018),'|',color='.5',ms=3,transform=ax.get_xaxis_transform())
        inset=ax.inset_axes([.55,.22,.40,.28])
        minus=[r for r in rows if r['symmetry_class']=='-'];origin=float(minus[0]['Lambda'])
        inset.plot([float(r['beta_deg']) for r in minus],
            [1e4*(float(r['Lambda'])-origin) for r in minus],color=color['-'],lw=1.1)
        inset.set(xlim=(0,90),xticks=[0,45,90],yticks=[-6,-3,0])
        inset.set_xlabel(r'$\beta$',labelpad=0,fontsize=10)
        inset.set_ylabel(r'$10^4[\Lambda_{-,1}(\beta)-\Lambda_{-,1}(0)]$',fontsize=9,labelpad=1)
        inset.tick_params(labelsize=9);style(inset)
        save(fig,'class_first_full','Первые частоты классов отражения; вставка показывает малое изменение класса «−».\nШтрихи у оси абсцисс — сохранённые угловые узлы. Сглаживание не применяется.')
        fig,ax=plt.subplots(figsize=(7.2,5.5));fig.subplots_adjust(left=.13,right=.97,top=.96,bottom=.35)
        for cls in ('+','-'):
            data=[r for r in local if r['symmetry_class']==cls]
            ax.plot([float(r['beta_deg']) for r in data],[float(r['Lambda']) for r in data],
                color=color[cls],lw=1.3,label=rf'$\Lambda_{{{cls},1}}$')
        beta=solution['beta_deg'];lam=solution['Lambda']
        ax.plot([beta],[lam],'o',ms=3,color='black')
        ax.annotate(r'$\beta_*\approx6.130043504885^\circ$',xy=(beta,lam),xytext=(6.145,3.921),
                    fontsize=11,arrowprops=dict(arrowstyle='-',color='.3',lw=.7))
        ax.set(xlim=(6.095,6.205),xlabel=r'$\beta$, градусы',ylabel=r'$\Lambda$');style(ax)
        for positions,col,level in (([6.1,6.2],'.35',.025),
                                    (sorted({float(r['beta_deg']) for r in local}),'.35',.045)):
            ax.plot(positions,np.full(len(positions),level),'|',color=col,ms=4,transform=ax.get_xaxis_transform())
        ax.legend(frameon=False,loc='upper left')
        save(fig,'class_first_crossing','Локальные корни блоков; чёрная точка — совместное решение двух уравнений.\nШтрихи снизу: архивные узлы (нижний ряд) и локальные вычисления (верхний ряд).')
        def coordinates(mode,beta):
            rad=np.deg2rad(beta)
            ts=[np.array([1.,0.]),np.array([-np.cos(rad),-np.sin(rad)])]
            ns=[np.array([0.,-1.]),np.array([-np.sin(rad),np.cos(rad)])]
            result=[]
            for arm,t,n in zip(mode['arms'],ts,ns):
                x=np.r_[0.,arm['x'],1.]; y=np.vstack((arm['ends'][0],arm['y'],arm['ends'][1]))
                result.append(((x-1)[:,None]*t,y[:,0,None]*t+y[:,1,None]*n))
            return result
        scenes=[(modes[c],beta) for c in ('+','-')]
        aligned={}
        for angle in (6.1,6.2):
            for cls in ('+','-'):
                original=first[(angle,cls)]['mode']
                overlap=sections.global_inner(first[(angle,cls)]['model'],modes[cls],original,beta,angle)
                aligned[(angle,cls)]=saved.combine(original,original,1 if overlap>=0 else -1,0)
                scenes.append((aligned[(angle,cls)],angle))
        amplitude=max(np.max(np.linalg.norm(d,axis=1)) for mode,angle in scenes for _,d in coordinates(mode,angle))
        scale=.28/amplitude
        def shape(ax,mode,angle,cls):
            for r,d in coordinates(mode,angle):
                ax.plot(r[:,0],r[:,1],color='.70',lw=.8)
                deformed=r+scale*d
                ax.plot(deformed[:,0],deformed[:,1],color=color[cls],lw=1.5)
                end=r[0];ax.plot([end[0],end[0]],[end[1]-.04,end[1]+.04],color='black',lw=2)
            ax.set_aspect('equal');ax.set(xlim=(-1.05,1.05),ylim=(-.32,.4))
            ax.set_xlabel(r'$X/l$',labelpad=2);ax.set_ylabel(r'$Y/l$',rotation=0,labelpad=12)
            ax.set_xticks([-1,0,1]);ax.set_yticks([-.2,0,.2]);style(ax)
        fig,axs=plt.subplots(1,2,figsize=(10.6,4.8));fig.subplots_adjust(left=.07,right=.97,top=.91,bottom=.32,wspace=.28)
        for ax,cls in zip(axs,('+','-')):
            shape(ax,modes[cls],beta,cls)
            ax.text(.5,1.07,'Класс «'+cls+'»',transform=ax.transAxes,ha='center')
        save(fig,'crossing_forms',f'Две независимые формы при β = {beta:.10f}°, Λ = {lam:.10f}.\nМассовая норма каждой формы равна единице; перемещения показаны в общем условном масштабе.')
        fig,axs=plt.subplots(2,2,figsize=(10.6,7.0));fig.subplots_adjust(left=.07,right=.97,top=.93,bottom=.24,wspace=.28,hspace=.42)
        for i,angle in enumerate((6.1,6.2)):
            for j,cls in enumerate(('+','-')):
                ax=axs[i,j];shape(ax,aligned[(angle,cls)],angle,cls)
                rank=first[(angle,cls)]['root']['sorted_position']
                ax.text(.5,1.07,f'β = {angle}°, класс «{cls}», позиция {rank}',transform=ax.transAxes,ha='center')
        save(fig,'archive_forms_sides','Сохранённые формы слева и справа от пересечения. Знаки согласованы внутри каждого класса.\nИзменяются порядковые номера частот; общий масштаб перемещений тот же, что в точке пересечения.')
    font=findfont(FontProperties(family='Times New Roman'),fallback_to_default=False)
    return dict(figures=captions,font=font,font_sha256=saved.sha(font),matplotlib=matplotlib.__version__,
        actual_nodes_visible=True,curve_smoothing=False,frequency_interpolation=False,
        line_segments_visual_only=True,shape_displacement_scale=scale)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--stage',choices=['a','all','plot'],default='all')
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    if args.stage=='all' and any((args.output/name).exists() for name in
            ('crossing_solution.json','block_roots.csv','crossing_forms.npz')):
        raise FileExistsError('Preserve completed local results; select a new --output before solving')
    rows,checks,first,policy,hashes=stage_a()
    if args.stage=='plot':
        with saved.no_spectral_calls():
            plots=render(args.output,first)
        saved.write_json(args.output/'plot_metadata.json',plots)
        return
    # Write-once artifacts; a failed stage B never overwrites stage A evidence.
    for name,value in [('stage_a_checks.json',checks)]:
        path=args.output/name
        if path.exists(): require(read(path)==value,'Existing stage A evidence differs')
        else: saved.write_json(path,value)
    if not (args.output/'branch_constancy_audit.csv').exists():
        saved.write_csv(args.output/'branch_constancy_audit.csv',rows)
    print(json.dumps(checks,ensure_ascii=False))
    if args.stage=='a': return
    require(all(v['status']=='BRANCH_VARIATION_CONFIRMED' for v in checks['summary'].values()),
            'Stage B requires established source branches')
    roots,solution,modes,physical,arrays,comparison=stage_b(first,policy)
    saved.write_csv(args.output/'block_roots.csv',roots)
    saved.write_json(args.output/'crossing_solution.json',solution)
    saved.write_json(args.output/'crossing_form_checks.json',physical)
    saved.write_csv(args.output/'archive_comparison.csv',comparison)
    with (args.output/'crossing_forms.npz').open('xb') as stream: np.savez_compressed(stream,**arrays)
    saved.write_json(args.output/'input_sha256.json',hashes)
    with saved.no_spectral_calls():
        plots=render(args.output,first)
    saved.write_json(args.output/'plot_metadata.json',plots)
    print(json.dumps(dict(beta_deg=solution['beta_deg'],Lambda=solution['Lambda'],omega=solution['omega'],
                         form_checks=physical,calls=solution['calls']),ensure_ascii=False))


if __name__=='__main__':
    main()
