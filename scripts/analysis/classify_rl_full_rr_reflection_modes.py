"""Read-only reflection classification of saved symmetric square RL/RR modes.

Script Proliferation Control: this explicitly requested, reusable postprocessor
has a field/projection I/O contract, unlike the spectrum runner or the historical
aspect-four modal audit. Reuses reflect_mode, full physical mass/energies, state
and boundary operators. It neither reconstructs fields nor searches frequencies.
All figures are local research outputs; no branch tracking is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack, contextmanager
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.analysis import map_rl_full_rr as archive
from scripts.analysis import postprocess_rl_full_rr_modes as historical
from scripts.lib import rl_full_rr_sections_prefix as sections
from scripts.lib.rl_full_rr_aspect import make_model

PHYSICS = sections.p
DEFAULT_SOURCE = ROOT / 'results/rl_full_rr_maps/b775d0980d7bad77'
DEFAULT_INVENTORY = ROOT / 'results/rl_full_rr_maps/_audit_reflection_symmetry_20261107/archive_inventory.json'
DEFAULT_OUTPUT = ROOT / 'results/rl_full_rr_maps/_audit_reflection_classification_20261107'
FIELDS = ('y', 'dy', 'ends')
EPS = np.finfo(float).eps
EXPECTED_MODEL = dict(epsilon=.005, nu=.3, aspect=1., E=1., rho=1., l=1.,
                      ks=5/6, variant='RL_FULL', mu=0., tau=0.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def write_csv(path, rows, columns=None):
    columns = columns or list(rows[0])
    with Path(path).open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


@contextmanager
def no_spectral_calls():
    """Extend the existing postprocessing guard to aliases of the square engine."""
    calls = []
    def forbidden(*args, **kwargs):
        calls.append('forbidden spectral/basis entry')
        raise AssertionError('Reflection classification may only read saved fields')
    names = {'basis', 'fields', 'recovered_modes', 'quadrature_check', 'boundary',
             'mode', 'prefix', 'count_at', 'counted_prefix', 'catalogues', 'ritz',
             'arm_dtn', 'solve', 'run', 'finite_roots', 'finite_mode',
             'finite_state_basis', 'explicit_frequencies', 'brentq', 'root_scalar',
             'bisect', 'newton'}
    with historical.no_solves(), ExitStack() as stack:
        for name, module in list(sys.modules.items()):
            if name.startswith('scripts.') and module is not sys.modules[__name__]:
                for attr in names:
                    if callable(getattr(module, attr, None)):
                        stack.enter_context(patch.object(module, attr, forbidden))
        yield calls


def combine(a, b, ca=1., cb=1.):
    arms = []
    for x, y in zip(a['arms'], b['arms']):
        if not all(np.array_equal(x[k], y[k]) for k in ('x', 'weights')):
            raise ValueError('Reflection requires identical material quadratures')
        arms.append(dict(x, **{k: ca*x[k]+cb*y[k] for k in FIELDS}))
    return dict(arms=arms)


def project(mode, parity):
    if parity not in (-1, 1):
        raise ValueError('Reflection parity must be +1 or -1')
    return combine(mode, sections.reflect_mode(mode), .5, .5*parity)


def component(mode, part):
    return dict(arms=[dict(a, **{k: getattr(a[k], part) for k in FIELDS})
                      for a in mode['arms']])


def is_complex(mode):
    return any(np.iscomplexobj(a[k]) for a in mode['arms'] for k in FIELDS)


def mass_inner(model, a, b):
    """Hermitian extension of the unchanged real full-mass inner product."""
    if not is_complex(a) and not is_complex(b):
        return PHYSICS.inner(model, a, b)
    ar, ai, br, bi = component(a, 'real'), component(a, 'imag'), component(b, 'real'), component(b, 'imag')
    return (PHYSICS.inner(model, ar, br)+PHYSICS.inner(model, ai, bi)
            +1j*(PHYSICS.inner(model, ar, bi)-PHYSICS.inner(model, ai, br)))


def mass(model, mode):
    value = mass_inner(model, mode, mode)
    if abs(np.imag(value)) > 64*EPS*max(abs(value), np.finfo(float).tiny):
        raise ArithmeticError('Non-real mass')
    return float(np.real(value))


def energy_components(model, case, mode):
    if not is_complex(mode):
        return PHYSICS.energies(model, case, mode['arms'])
    pairs = [PHYSICS.energies(model, case, component(mode, part)['arms'])
             for part in ('real', 'imag')]
    return tuple(a+b for a, b in zip(*pairs))


def relative_fields(a, b, reference):
    return max(float(np.max(abs(x[k]-y[k])) /
                     max(np.max(abs(z[k])), np.finfo(float).tiny))
               for x, y, z in zip(a['arms'], b['arms'], reference['arms']) for k in FIELDS)


def scales(model, mode):
    q = max(np.max(abs(a['y'][:, :3]*[1, 1, model.l])) for a in mode['arms'])
    f = max(np.max(abs(a['y'][:, 3:]*[1, 1, 1/model.l])) for a in mode['arms'])
    return float(q), float(f)


def physical_checks(model, case, omega, mode, reference=None):
    """Production residual scaling, using saved y/dy/ends; no basis calls."""
    mm, kk, spring = energy_components(model, case, mode)
    norm = mass(model, mode)
    if norm <= 0:
        return dict(mass=norm, energy=None, energy_error=None, ode=None, joint=None,
                    external=None, ode_source_scaled=None, joint_source_scaled=None,
                    external_source_scaled=None, mass_integrator_error=None)
    q, f = scales(model, mode)
    if q <= 0 or f <= 0:
        return dict(mass=norm, energy=float(kk.sum()), energy_error=None, ode=None,
                    joint=None, external=None, ode_source_scaled=None,
                    joint_source_scaled=None, external_source_scaled=None,
                    mass_integrator_error=float(abs(mm.sum()-norm)))
    jr = PHYSICS.joint_matrix(model, case) @ np.concatenate([a['ends'][1] for a in mode['arms']])
    er = np.array([PHYSICS.external_matrix(model, case, i) @ a['ends'][0]
                   for i, a in enumerate(mode['arms'])])
    residuals = [a['dy']-a['y'] @ PHYSICS.base.state_matrix(arm, omega).T
                 for a, arm in zip(mode['arms'], model.arms)]
    def scaled(qs, fs, ref):
        scale = np.array([1, 1, model.l, qs/fs, qs/fs, qs/(fs*model.l)])
        ode = max(float(np.max(abs(res*scale)) /
                        max(np.max(abs(a['dy']*scale)), 1e-30))
                  for res, a in zip(residuals, ref['arms']))
        joint = float(np.max(abs(jr)/[qs, qs, fs, fs, qs/model.l, fs*model.l]))
        external = float(np.max(abs(er)/[qs, qs, qs/model.l]))
        return ode, joint, external
    ode, joint, ext = scaled(q, f, mode)
    ref = reference or mode
    ro, rj, re = scaled(*scales(model, ref), ref)
    energy = float(kk.sum()+spring.sum())
    return dict(mass=norm, energy=energy, energy_error=float(abs(energy/norm/omega**2-1)),
                ode=ode, joint=joint, external=ext, ode_source_scaled=ro,
                joint_source_scaled=rj, external_source_scaled=re,
                mass_integrator_error=float(abs(mm.sum()-norm)))


def physical_pass(check, policy):
    return (check['mass'] > 0 and all(check[k] is not None and np.isfinite(check[k])
            and check[k] <= policy[tol] for k, tol in
            [('ode', 'ode_tol'), ('joint', 'boundary_tol'), ('external', 'boundary_tol'),
             ('energy_error', 'energy_tol')]))


def analyze_mode(model, case, root, mode):
    reflected = sections.reflect_mode(mode)
    plus, minus = project(mode, 1), project(mode, -1)
    m0 = mass(model, mode)
    if not np.isfinite(m0) or m0 <= 0:
        raise ArithmeticError('Saved nonzero mode has nonpositive mass')
    mp, mn = mass(model, plus), mass(model, minus)
    xp, xn = mp/m0, mn/m0
    dp = np.sqrt(max(0., mass(model, combine(reflected, mode, 1, -1))/m0))
    dn = np.sqrt(max(0., mass(model, combine(reflected, mode))/m0))
    overlap = abs(mass_inner(model, plus, minus))
    # Conservative floating-point summation bound: 64 scalar operations/node.
    operation_count = 64*sum(len(a['x']) for a in mode['arms'])
    identity_tol = operation_count*EPS/(1-operation_count*EPS)
    identities = dict(fractions_sum=abs(xp+xn-1), delta_plus=abs(dp*dp-4*xn),
        delta_minus=abs(dn*dn-4*xp), orthogonality=float(overlap/m0),
        reconstruction=relative_fields(combine(plus, minus), mode, mode),
        involution=relative_fields(sections.reflect_mode(reflected), mode, mode),
        mass_isometry=abs(mass(model, reflected)/m0-1),
        idempotence_plus=relative_fields(project(plus, 1), plus, mode),
        idempotence_minus=relative_fields(project(minus, -1), minus, mode))
    original = physical_checks(model, case, root['omega'], mode)
    positive = physical_checks(model, case, root['omega'], plus, mode)
    negative = physical_checks(model, case, root['omega'], minus, mode)
    return dict(xi_plus=xp, xi_minus=xn, delta_plus=float(dp), delta_minus=float(dn),
                identity_tol=identity_tol, identities=identities,
                identities_pass=all(x <= identity_tol for x in identities.values()),
                normalized_projection_overlap=float(overlap/np.sqrt(mp*mn)) if mp*mn > 0 else None,
                original=original, plus=positive, minus=negative)


def classify(d, policy, threshold, *, multiplicity=1, multiplicity_confirmed=False):
    if not d['identities_pass'] or not physical_pass(d['original'], policy):
        return '', 'INVALID_PROJECTION', 'source_physics_or_projector_identity_failed'
    if threshold is None:
        return '', 'UNRESOLVED_MIXED', 'no_stable_defect_threshold_interval'
    pure = [sign for sign in ('plus', 'minus') if d['delta_'+sign] <= threshold]
    if len(pure) == 1:
        sign = pure[0]
        if physical_pass(d[sign], policy):
            return '+' if sign == 'plus' else '-', 'PURE_'+sign.upper(), 'one_resolved_class'
        return '', 'INVALID_PROJECTION', 'dominant_projection_physics_failed'
    if not all(physical_pass(d[s], policy) for s in ('plus', 'minus')):
        return '', 'INVALID_PROJECTION', 'substantial_projection_physics_failed'
    independent = (d['normalized_projection_overlap'] is not None and
                   d['normalized_projection_overlap'] <= d['identity_tol'])
    if multiplicity >= 2 and multiplicity_confirmed and independent:
        return '+/-', 'MIXED_DEGENERATE', 'confirmed_multiple_eigenspace'
    return '', 'UNRESOLVED_MIXED', 'no_confirmed_multiplicity_or_independence'


def threshold_interval(diagnostics):
    """Find a separated defect band, then reserve decades on both sides.

    This is an archive-specific diagnostic rule, never a solver tolerance.
    No target crossing or expected class is supplied to the rule.
    """
    values = np.sort([max(d[k], EPS) for d in diagnostics for k in ('delta_plus', 'delta_minus')])
    gaps = np.diff(np.log10(values))
    i = int(np.argmax(gaps))
    low_band, high_band = float(values[i]), float(values[i+1])
    low = float(10.**np.ceil(np.log10(low_band)+1))
    high = float(10.**np.floor(np.log10(high_band)-2))
    stable = high >= 10*low and high < 1
    return dict(low_band_max=low_band, high_band_min=high_band,
                separation_decades=float(gaps[i]), lower=low if stable else None,
                upper=high if stable else None, chosen=low if stable else None,
                method='largest log defect gap; one decade above lower band, two below upper; require decade-wide interval',
                status='SEPARATED_BANDS' if stable else 'NO_STABLE_INTERVAL')


def confirmed_multiple_event(record, root):
    """Require a complete, count-and-nullity-validated saved event; no grouping by closeness."""
    if root['multiplicity'] < 2:
        return False
    group = [r for r in record['roots'] if r['event_id'] == root['event_id']]
    events = [e for e in record['search']['events'] if e.get('event_id') == root['event_id']]
    # Older formats without explicit nullity evidence remain unresolved.
    return (len(group) == root['multiplicity'] and
            all(r['omega'] == root['omega'] and r['multiplicity'] == root['multiplicity'] for r in group) and
            any(e.get('nullity') == root['multiplicity'] and
                e.get('multiplicity') == root['multiplicity'] for e in events))


def load_source(folder, inventory_path):
    inventory = archive.read(inventory_path)
    audit_check = inventory_path.with_name('final_checks.json')
    if (not audit_check.exists() or sha(inventory_path) !=
            archive.read(audit_check)['evidence_sha256']['archive_inventory.json']):
        raise ValueError('Inventory checksum lacks its analytical-audit verification')
    if local(folder) != inventory['source_bundle']:
        raise ValueError('Inventory belongs to a different source bundle')
    for name, field in [('manifest.json', 'source_manifest_sha256'), ('spectrum.csv', 'source_csv_sha256')]:
        if sha(folder/name) != inventory[field]:
            raise ValueError('Source checksum differs from analytical audit: '+name)
    manifest, table = archive.saved_table(folder)
    if manifest['status'] != 'PASS':
        raise ValueError('Accepted source bundle required')
    source_hashes = {local(folder/name): sha(folder/name) for name in ('manifest.json', 'spectrum.csv', 'identity.json')}
    source_hashes[local(inventory_path)] = sha(inventory_path)
    source_hashes[local(audit_check)] = sha(audit_check)
    for name in ('spectrum.csv', 'identity.json'):
        if sha(folder/name) != manifest['artifacts'][name]:
            raise ValueError('Source artifact changed: '+name)
    for path, expected in manifest['identity']['sources'].items():
        if sha(ROOT/path) != expected:
            raise ValueError('Physical/archival implementation changed: '+path)
        source_hashes[path] = expected
    entries = [e for e in manifest['points'] if e['mu'] == 0 and e['tau'] == 0]
    if len(entries) != 197 or len({e['beta_deg'] for e in entries}) != 197:
        raise ValueError('197 distinct symmetric angles required')
    index = {p['beta_deg']: p for p in inventory['points']}
    if len(index) != 197 or sorted(index) != sorted(e['beta_deg'] for e in entries):
        raise ValueError('Incomplete inventory')
    rows = {(r['map'], float(r['x']), int(r['rank'])): r for r in table}
    if len(rows) != len(table):
        raise ValueError('Duplicate CSV rows')
    loaded = []
    policy = None
    for e in sorted(entries, key=lambda e: e['beta_deg']):
        item = index[e['beta_deg']]
        path = ROOT/e['record']
        if e['status'] != 'PASS' or item['record'] != e['record'] or sha(path) != e['record_sha256'] or sha(path) != item['record_sha256']:
            raise ValueError('Point provenance mismatch')
        record = archive.verify_record(path)
        contract = record['contract']
        model = make_model(**contract['model']); case = PHYSICS.Case(**contract['case'])
        if contract['model'] != EXPECTED_MODEL or case.beta != e['beta_deg']:
            raise ValueError('Only the frozen symmetric square contract is accepted')
        if sections.point_contract(model, case, 8, 9, contract['policy']) != contract:
            raise ValueError('Current and saved physical contracts differ')
        if policy is None: policy = contract['policy']
        if contract['policy'] != policy or contract['arms'][0] != contract['arms'][1]:
            raise ValueError('Different numerical policy or nonidentical arms')
        if (record['status'] != 'PASS' or record['quality']['status'] != 'PASS' or
            len(record['roots']) != 9 or record['search']['status'] != 'PASS' or
            record['search']['lower_count'] != 0 or record['search']['guard_right_count'] != 9 or
            record['search']['computed_ranks'] != 9 or record['search']['unresolved_intervals']):
            raise ValueError('Saved nine-rank completeness/quality evidence missing')
        npz = path.parent.parent/record['shape_file']
        if local(npz) != item['npz'] or sha(npz) != item['npz_sha256']:
            raise ValueError('Shape inventory mismatch')
        source_hashes[local(path)] = sha(path); source_hashes[local(npz)] = sha(npz)
        with np.load(npz, allow_pickle=False) as arrays:
            if len(arrays.files) != 90:
                raise ValueError('Incomplete NPZ')
            for k, root in enumerate(record['roots'], 1):
                csvrow = rows[(e['map'], e['beta_deg'], k)]
                if (root['sorted_position'] != k or csvrow['status'] != 'PASS' or
                    any(float(csvrow[key]) != root[key] for key in ('Lambda', 'omega', 'Omega'))):
                    raise ValueError('CSV/record frequency disagreement')
                if max(abs(root['Omega']/root['Lambda']**2-1),
                       abs(root['omega']*contract['time_reference']/root['Lambda']**2-1)) > 32*EPS:
                    raise ValueError('Frequency normalization disagreement')
                arms = [{key: arrays[f'mode{k}_arm{arm}_{key}'].copy()
                         for key in ('x', 'weights', *FIELDS)} for arm in (1, 2)]
                for arm in arms:
                    if (arm['x'].shape != (120,) or arm['weights'].shape != (120,) or
                        arm['y'].shape != (120, 6) or arm['dy'].shape != (120, 6) or arm['ends'].shape != (2, 6) or
                        not all(np.isfinite(v).all() for v in arm.values()) or
                        not np.all(arm['weights'] > 0) or not np.all(np.diff(arm['x']) > 0)):
                        raise ValueError('Invalid saved physical arrays')
                if not all(np.array_equal(arms[0][key], arms[1][key]) for key in ('x', 'weights')):
                    raise ValueError('Unequal reflected quadratures')
                loaded.append(dict(beta=e['beta_deg'], root=root, mode=dict(arms=arms),
                    record=record, model=model, case=case, point_id=path.stem,
                    point_path=local(path), npz=local(npz), shape_key=f'mode{k}',
                    multiple_confirmed=confirmed_multiple_event(record, root)))
    if len(loaded) != 1773:
        raise ValueError('Expected 1773 forms')
    np.testing.assert_array_equal(sections.REFLECTION_STATE, [1, -1, -1, 1, -1, -1])
    return loaded, policy, source_hashes


def distribution(values):
    a = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    return dict(count=len(a), zero_count=int(sum(a == 0)),
                **{name: float(np.quantile(a, q)) for name, q in
                   [('min', 0), ('q01', .01), ('median', .5), ('q95', .95), ('q99', .99), ('max', 1)]}) if len(a) else dict(count=0)


def collect_diagnostics(items):
    result = []
    for item in items:
        d = analyze_mode(item['model'], item['case'], item['root'], item['mode'])
        saved = item['root']['diagnostics']
        d['saved_diagnostic_difference'] = dict(
            ode=abs(d['original']['ode']-saved['ode_error']),
            joint=abs(d['original']['joint']-max(saved['joint_scaled'])),
            external=abs(d['original']['external']-np.max(saved['external_scaled'])),
            energy=abs(d['original']['energy_error']-saved['energy_error']),
            mass=abs(d['original']['mass']-saved['mass']))
        if abs(d['original']['mass']-1) > item['record']['contract']['policy']['mass_gram_tol']:
            d['identities_pass'] = False
        result.append(d)
    return result


def classification_rows(items, diagnostics, policy, threshold):
    rows = []
    for item, d in zip(items, diagnostics):
        root = item['root']
        cls, status, reason = classify(d, policy, threshold, multiplicity=root['multiplicity'],
                                       multiplicity_confirmed=item['multiple_confirmed'])
        row = dict(beta_deg=item['beta'], sorted_position=root['sorted_position'],
            Lambda=root['Lambda'], omega=root['omega'], multiplicity=root['multiplicity'],
            xi_plus=d['xi_plus'], xi_minus=d['xi_minus'], delta_plus=d['delta_plus'], delta_minus=d['delta_minus'],
            ODE_residual=d['original']['ode'], joint_residual=d['original']['joint'],
            external_residual=d['original']['external'], energy_quotient_error=d['original']['energy_error'],
            symmetry_class=cls, classification_status=status, reason=reason,
            source_point_id=item['point_id'], source_shape_key=item['shape_key'],
            source_point_json=item['point_path'], source_npz=item['npz'],
            mass=d['original']['mass'], identities_pass=d['identities_pass'])
        for sign in ('plus', 'minus'):
            row.update({f'{sign}_{key}': value for key, value in d[sign].items()})
            fraction = d['xi_'+sign]
            row[sign+'_interpretation'] = ('ZERO_PROJECTION' if fraction == 0 else
                'BELOW_CLASSIFICATION_RESOLUTION' if threshold is not None and fraction <= (threshold/2)**2 else
                'PHYSICAL_PROJECTION' if physical_pass(d[sign], policy) else 'PHYSICAL_CHECK_FAILED')
        rows.append(row)
    return rows


def spectra_and_crossings(rows):
    spectra, coverage, first = [], [], {}
    for beta in sorted({r['beta_deg'] for r in rows}):
        group = sorted([r for r in rows if r['beta_deg'] == beta], key=lambda r: r['sorted_position'])
        unresolved = [r for r in group if r['classification_status'] not in ('PURE_PLUS', 'PURE_MINUS')]
        # Mixed multiple representatives need subspace deduplication, never duplicate
        # both projections of each saved vector as additional eigenvalues.
        for sign in ('+', '-'):
            pure = sorted([r for r in group if r['symmetry_class'] == sign and r['classification_status'].startswith('PURE_')], key=lambda r: r['Lambda'])
            for rank, row in enumerate(pure, 1):
                complete = not any(u['Lambda'] <= row['Lambda'] for u in unresolved)
                spectra.append(dict(beta_deg=beta, symmetry_class=sign, class_rank=rank,
                    Lambda=row['Lambda'], original_sorted_position=row['sorted_position'],
                    shape_key=row['source_shape_key'], classification_status=row['classification_status'],
                    prefix_complete_below=complete, source_npz=row['source_npz']))
                if rank == 1 and complete: first[(beta, sign)] = row
            coverage.append(dict(beta_deg=beta, symmetry_class=sign, saved_class_count=len(pure),
                first_available=(beta, sign) in first, unresolved_full_ranks=len(unresolved),
                status='COMPLETE_CLASS_PREFIX_TO_SAVED_CUTOFF' if pure and not unresolved else 'INCOMPLETE_CLASS_PREFIX',
                further_class_ranks='NOT_IN_SAVED_PREFIX'))
    gaps, brackets = [], []
    for beta in sorted({r['beta_deg'] for r in rows}):
        p, n = first.get((beta, '+')), first.get((beta, '-'))
        gaps.append(dict(beta_deg=beta, Lambda_plus=p['Lambda'] if p else None,
            Lambda_minus=n['Lambda'] if n else None, gap=p['Lambda']-n['Lambda'] if p and n else None,
            sorted_position_plus=p['sorted_position'] if p else None,
            sorted_position_minus=n['sorted_position'] if n else None))
    for left, right in zip(gaps[:-1], gaps[1:]):
        if left['gap'] is None or right['gap'] is None: continue
        if left['gap']*right['gap'] < 0:
            brackets.append(dict(beta_left=left['beta_deg'], beta_right=right['beta_deg'],
                Lambda_plus_left=left['Lambda_plus'], Lambda_minus_left=left['Lambda_minus'], gap_left=left['gap'],
                Lambda_plus_right=right['Lambda_plus'], Lambda_minus_right=right['Lambda_minus'], gap_right=right['gap'],
                first_two_are_different_classes_left={left['sorted_position_plus'], left['sorted_position_minus']} == {1, 2},
                first_two_are_different_classes_right={right['sorted_position_plus'], right['sorted_position_minus']} == {1, 2},
                order_exchanged=left['sorted_position_plus'] != right['sorted_position_plus'],
                status='SIGN_CHANGE_BETWEEN_SAVED_NODES_NO_EXACT_ANGLE'))
    return spectra, coverage, gaps, brackets


def render(output, rows, gaps, threshold):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font = font_manager.findfont('Times New Roman', fallback_to_default=False)
    colors = {'+': '#1f77b4', '-': '#ff7f0e'}
    figures = []
    with plt.rc_context({'font.family': 'serif', 'font.serif': ['Times New Roman'],
            'font.size': 12, 'mathtext.fontset': 'stix', 'pdf.fonttype': 42,
            'axes.grid': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
            'path.simplify': False}):
        def axes_style(ax, ylabel):
            ax.set_xlabel(r'$\beta$, градусы'); ax.set_ylabel(ylabel, rotation=0, labelpad=17)
            ax.tick_params(direction='in', top=True, right=True, length=3, width=.7)
            ax.grid(False, which='both')
            for spine in ax.spines.values(): spine.set_linewidth(.7)
        def save(fig, name, caption):
            fig.text(.53, .125, r'$\mu=\tau=0,\quad a=1,\quad\varepsilon=0.005$', ha='center')
            fig.text(.53, .035, caption, ha='center', va='bottom', fontsize=10)
            for ext in ('png', 'pdf'): fig.savefig(output/(name+'.'+ext), dpi=300, facecolor='white')
            figures.append(dict(name=name, caption=caption))
            plt.close(fig)
        for name, interval in [('class_first_full', (0, 90)), ('class_first_zoom', (5, 7))]:
            fig, ax = plt.subplots(figsize=(6.7, 5.1)); fig.subplots_adjust(left=.13, right=.97, top=.96, bottom=.32)
            selected = [g for g in gaps if interval[0] <= g['beta_deg'] <= interval[1]]
            for sign, key in [('+', 'Lambda_plus'), ('-', 'Lambda_minus')]:
                ax.plot([g['beta_deg'] for g in selected], [g[key] if g[key] is not None else np.nan for g in selected],
                        color=colors[sign], lw=1.4, ls='-', marker=None, label=rf'$\Lambda_{{{sign},1}}$')
            axes_style(ax, r'$\Lambda$'); ax.set_xlim(*interval)
            ax.legend(loc='upper center', bbox_to_anchor=(.5, -.18), ncol=2, frameon=False)
            save(fig, name, 'Первые частоты двух классов отражения.\nОтрезки соединяют сохранённые узлы; точный угол совпадения не определён.')
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.8)); fig.subplots_adjust(left=.08, right=.98, top=.95, bottom=.33, wspace=.30)
        for rank, ax in enumerate(axs, 1):
            selected = [r for r in rows if r['sorted_position'] == rank]
            for sign in ('+', '-'):
                ax.plot([r['beta_deg'] for r in selected],
                        [r['Lambda'] if r['symmetry_class'] == sign else np.nan for r in selected],
                        color=colors[sign], lw=1.4, ls='-', marker=None, label='Класс '+sign)
            axes_style(ax, rf'$\Lambda_{rank}$'); ax.set_xlim(0, 90)
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .18), ncol=2, frameon=False)
        save(fig, 'sorted_pair_classes', 'Первая и вторая упорядоченные частоты; цвет обозначает класс отражения.\nЛинии разорваны между узлами смены класса; это не продолжение по номеру.')
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.8)); fig.subplots_adjust(left=.10, right=.98, top=.95, bottom=.33, wspace=.32)
        for rank, ax in enumerate(axs, 1):
            selected = [r for r in rows if r['sorted_position'] == rank]
            for sign, key in [('+', 'delta_plus'), ('-', 'delta_minus')]:
                ax.semilogy([r['beta_deg'] for r in selected], [r[key] if r[key] > 0 else np.nan for r in selected],
                            color=colors[sign], lw=1.2, ls='-', marker=None, label=rf'$\delta_{{{sign}}}$')
            if threshold is not None: ax.axhline(threshold, color='#777777', lw=.8, ls='-')
            axes_style(ax, rf'$\delta_{{\pm,{rank}}}$'); ax.set_xlim(0, 90)
        handles, labels = axs[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .18), ncol=2, frameon=False)
        save(fig, 'sorted_pair_symmetry_defects', 'Дефекты отражательной симметрии первой и второй упорядоченных форм.\nСерая линия — диагностический порог; нулевые значения на логарифмической оси не показаны.')
    return dict(figures=figures, font=font, font_sha256=sha(font), matplotlib=matplotlib.__version__,
                no_markers=True, no_grid=True, no_smoothing=True, frequency_interpolation=False)


def run(source, inventory_path, output, diagnostics_only=False):
    output.mkdir(parents=True, exist_ok=True)
    if (output/'quality_checks.json').exists() or (output/'diagnostic_distributions.json').exists():
        raise FileExistsError('Completed or diagnostic results exist; preserve them and choose another output')
    with no_spectral_calls() as forbidden_calls:
        items, policy, source_hashes = load_source(source, inventory_path)
        diagnostics = collect_diagnostics(items)
        distributions = {k: distribution([d[k] for d in diagnostics])
                         for k in ('xi_plus', 'xi_minus', 'delta_plus', 'delta_minus')}
        for label in ('original', 'plus', 'minus'):
            for key in diagnostics[0][label]:
                distributions[label+'_'+key] = distribution([d[label][key] for d in diagnostics])
        for label in ('dominant', 'minority'):
            selected = [d[('plus' if d['xi_plus'] >= d['xi_minus'] else 'minus') if label == 'dominant'
                          else ('minus' if d['xi_plus'] >= d['xi_minus'] else 'plus')] for d in diagnostics]
            for key in selected[0]: distributions[label+'_'+key] = distribution([d[key] for d in selected])
        interval = threshold_interval(diagnostics)
        write_json(output/'diagnostic_distributions.json', dict(distributions=distributions, threshold_analysis=interval,
                   physical_tolerances={k:policy[k] for k in ('ode_tol', 'boundary_tol', 'energy_tol', 'mass_gram_tol')}))
        if diagnostics_only:
            print(json.dumps(dict(phase='diagnostics_only', threshold=interval,
                original=distributions['original_energy_error'], dominant=distributions['dominant_energy_error'],
                minority_mass=distributions['minority_mass']), indent=2))
            return
        rows = classification_rows(items, diagnostics, policy, interval['chosen'])
        sensitivity = []
        if interval['chosen'] is not None:
            baseline = [(r['symmetry_class'], r['classification_status']) for r in rows]
            for threshold in np.geomspace(interval['lower'], interval['upper'], 9):
                other = classification_rows(items, diagnostics, policy, float(threshold))
                sensitivity.append(dict(threshold=float(threshold), changed=sum((r['symmetry_class'], r['classification_status']) != b for r,b in zip(other,baseline)),
                    **dict(Counter(r['classification_status'] for r in other))))
            if any(s['changed'] for s in sensitivity):
                raise ArithmeticError('Classification is not stable over the proposed interval')
        write_csv(output/'mode_classification.csv', rows)
        spectra, coverage, gaps, brackets = spectra_and_crossings(rows)
        write_csv(output/'class_spectrum.csv', spectra, list(spectra[0]) if spectra else ['beta_deg','symmetry_class','class_rank','Lambda','original_sorted_position','shape_key','classification_status'])
        write_csv(output/'class_coverage.csv', coverage)
        write_csv(output/'class_first_gaps.csv', gaps)
        write_csv(output/'crossing_brackets.csv', brackets, list(brackets[0]) if brackets else ['beta_left','beta_right','Lambda_plus_left','Lambda_minus_left','gap_left','Lambda_plus_right','Lambda_minus_right','gap_right','first_two_are_different_classes_left','first_two_are_different_classes_right','order_exchanged','status'])
        write_json(output/'projection_checks.json', [dict(beta_deg=i['beta'], sorted_position=i['root']['sorted_position'], **d) for i,d in zip(items,diagnostics)])
        quality = dict(status='COMPLETE' if all(r['classification_status'].startswith('PURE_') for r in rows) else 'UNRESOLVED_OR_INVALID_FORMS_RETAINED',
            source_points=197, forms=len(rows), counts=dict(Counter(r['classification_status'] for r in rows)),
            physical_policy={k:policy[k] for k in ('ode_tol', 'boundary_tol', 'energy_tol', 'mass_gram_tol')},
            threshold_analysis=interval, threshold_sensitivity=sensitivity,
            identity_maxima={key:max(d['identities'][key] for d in diagnostics) for key in diagnostics[0]['identities']},
            identity_tolerance=max(d['identity_tol'] for d in diagnostics),
            saved_diagnostic_difference_maxima={key:max(d['saved_diagnostic_difference'][key] for d in diagnostics) for key in diagnostics[0]['saved_diagnostic_difference']},
            source_mass_error=max(abs(d['original']['mass']-1) for d in diagnostics),
            projection_failures={sign:sum(d[sign]['mass'] > 0 and not physical_pass(d[sign],policy) for d in diagnostics) for sign in ('plus','minus')},
            unresolved=[{k:r[k] for k in ('beta_deg','sorted_position','classification_status','reason')} for r in rows if not r['classification_status'].startswith('PURE_')],
            crossing_brackets=brackets, exact_saved_equalities=[g for g in gaps if g['gap'] == 0],
            first_class_frequencies_complete=all(c['first_available'] for c in coverage),
            forbidden_calls=forbidden_calls, root_calls=0, spectral_count_calls=0, basis_calls=0, ritz_calls=0)
        quality['plots'] = render(output, rows, gaps, interval['chosen'])
        write_json(output/'quality_checks.json', quality)
        # Validate sources again after all processing. No arrays are copied to outputs.
        for path, expected in source_hashes.items():
            if sha(ROOT/path) != expected: raise ValueError('Source changed during processing: '+path)
        artifacts = {p.relative_to(output).as_posix():sha(p) for p in output.rglob('*')
                     if p.is_file() and p.name != 'manifest.json'}
        implementations = {local(__file__):sha(__file__), local(historical.__file__):sha(historical.__file__)}
        test_path = ROOT/'tests/test_rl_full_rr_reflection_classification.py'
        if test_path.exists(): implementations[local(test_path)] = sha(test_path)
        write_json(output/'manifest.json', dict(status=quality['status'], source_bundle=local(source),
            source_sha256=source_hashes, implementation_sha256=implementations, artifacts=artifacts,
            python=platform.python_version(), numpy=np.__version__, spectrum_semantics='reflection_class_sorted_positions',
            calculation_mode='saved_fields_only', root_calls=0, spectral_count_calls=0, basis_calls=0, ritz_calls=0))
        print(json.dumps({k:quality[k] for k in ('status','forms','counts','threshold_analysis','crossing_brackets')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--inventory', type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--diagnostics-only', action='store_true')
    args = parser.parse_args()
    run(args.source.resolve(), args.inventory.resolve(), args.output.resolve(), args.diagnostics_only)


if __name__ == '__main__':
    main()
