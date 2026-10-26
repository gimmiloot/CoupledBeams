"""Four report figures from frozen RL_FULL/Timoshenko archives; plot_only.

Script Proliferation Control: this stable report-facing entry point joins six
different archive schemas into one reviewed output contract. It is not a new
parameter sweep. Keeping it separate avoids extending the 4500-line calculation
CLI with report layout. Reuses that CLI's hashes, CSV/JSON writers, exact archived
derivative stencils, persistent-threshold selection and saved-field energy helper.
No eigenproblem, matrix assembly, mode reconstruction or interpolation is allowed.
"""
from pathlib import Path
import argparse
import csv
import itertools
import json
import math
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.analysis import verify_isotropic_two_springs_rl_timoshenko as archive

BASE = ROOT / 'results/isotropic_two_springs_rl_timoshenko'
OUTPUT = BASE / 'report_illustrations'
PAIR = 'J0E4'
SOURCES = {
    'growth': ('growth_stiffness_5x5/361e9a839ef5be30',
               'edffda3081ab53afa6676b41368d7f585739fe466e7a1c448e134c24bb86ee42'),
    'thickness': ('growth_thickness_comparison/c4a723a2de351365',
                  '8a10d90f124e00c64a2db74c72450a6464993a549fb4ef669d4df6bd3554ee6d'),
    'epsilon': ('epsilon_scaling/d925e14c6ff93605',
                'dd0bcb344bbfa864a3780340f1bebddf23b5f67c1fa5ceddd8e21e1b82879129'),
    'variation': ('angle_variation/58972c1ea3e25f6f',
                  '85240037655915aa1f6e2470d795e471a539b111b0b6461291e521fd1d427871'),
    'slip': ('longitudinal_slip/253c9655120041e1',
             'a18235a073ee2736d2a7469cc50c94e18b896db43074a0a25ab78106a10a7f51'),
    'leading': ('leading_validation/f1c670fe62102196',
                '1a7d3d89f11c877194a27f8384c4cd9a7557cded2acbfc7211e58272b5c36fa4'),
}
FIGURES = ['01_thickness_frequency_and_q', '02_longitudinal_compliance',
           '03_leading_model_validation', '04_joint_displacement_and_slope']
POLICY = dict(frequency_map_policy='frequency-map-v1', calculation_mode='plot_only',
              spectrum_semantics='sorted_positions', sweep_parameter='beta_deg or xi',
              parameter_grid='actual nodes in exported CSV; q tail checked through 90 deg',
              K_plot=1, K_guard=3, guard_root_role='inherited three-root archive prefix',
              neighbour_audit='inherited accepted archive checks',
              local_repair_policy='none; fail on inconsistent or missing archive',
              strict_audit_default=False,
              policy_override_reason='Replot accepted historical K=2/guard3 prefixes; no new roots.')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label, rtol=1e-11, atol=2e-13):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol), label)


def ordered(rows, key='beta_deg'):
    result = sorted(rows, key=lambda r: float(r[key]))
    x = [float(r[key]) for r in result]
    require(len(x) > 0 and len(set(x)) == len(x), 'Empty or duplicate nodes: '+key)
    return result


def values(rows, field):
    return np.array([float(r[field]) for r in rows])


class SavedData:
    """Validate only the sources actually used, following their saved point paths."""
    def __init__(self):
        self.inputs = {}
        self.manifests = {}
        self.tables = {}
        self.points = {}
        self.overlap = {}
        for name, (relative, digest) in SOURCES.items():
            path = BASE / relative / 'manifest.json'
            self.check_hash(path, digest)
            manifest = json.loads(path.read_text(encoding='utf8'))
            require(manifest['status'] == 'PASS', 'Unaccepted archive: '+name)
            cfg = manifest['identity']['config']
            require(cfg['mu'] == .01 and cfg['cases'][PAIR] == [.1, 10],
                    'Incorrect construction: '+name)
            if name != 'leading':
                expected = dict(E=1., rho=1., nu=.3, epsilon=.01 if name == 'thickness' else .005,
                                aspect=4., l=1., ks=5/6)
                require(cfg['material'] == expected and cfg['variant'] == 'RL_FULL',
                        'Incorrect physical model: '+name)
                physical = 'scripts/lib/isotropic_two_springs_rl_timoshenko.py'
                self.check_hash(ROOT / physical, manifest['identity']['files'][physical])
                if name == 'epsilon':
                    for epsilon in (.0025, .005, .01):
                        filename = f'eps_{epsilon:g}/material_contract.json'
                        material_path = BASE / relative / filename
                        self.check_hash(material_path, manifest['artifacts'][filename])
                        material = json.loads(material_path.read_text(encoding='utf8'))
                        require(material == {**expected, 'epsilon': epsilon},
                                'Epsilon-specific material contract mismatch')
            else:
                require(cfg['epsilon_values'] == [.0025, .005, .01], 'Leading epsilon contract')
                require(cfg['ritz_orders'][-1] == 20, 'Leading approximation order')
                for source, digest in cfg['sources'].items():
                    self.check_hash(ROOT/source/'manifest.json', digest)
            self.manifests[name] = manifest

    def check_hash(self, path, expected=None):
        actual = archive.sha(path)
        require(expected is None or actual == expected, 'Source hash mismatch: '+str(path))
        key = path.relative_to(ROOT).as_posix()
        require(key not in self.inputs or self.inputs[key] == actual, 'Source changed during read')
        self.inputs[key] = actual
        return actual

    def read(self, name, filename):
        key = name, filename
        if key not in self.tables:
            path = BASE / SOURCES[name][0] / filename
            self.check_hash(path, self.manifests[name]['artifacts'][filename])
            with path.open(encoding='utf8', newline='') as stream:
                self.tables[key] = list(csv.DictReader(stream))
        return self.tables[key]

    def provenance(self, name, filename):
        return dict(source_archive=SOURCES[name][0], source_csv=filename)

    def point(self, row, epsilon, chi='RIGID'):
        path = ROOT / row['point_file']
        self.check_hash(path, row['point_sha256'])
        if path not in self.points:
            self.points[path] = json.loads(path.read_text(encoding='utf8'))
        rec = self.points[path]
        require(rec['status'] == rec['count']['status'] == 'PASS', 'Unaccepted point')
        case = rec['case']
        require([case[k] for k in ('beta', 'mu', 'k_joint', 'k_end', 'variant')]
                == [float(row['beta_deg']), .01, .1, 10, 'RL_FULL'], 'Point contract mismatch')
        require(str(case.get('chi', 'RIGID')) == str(chi) or
                (chi != 'RIGID' and float(case.get('chi', -1)) == float(chi)),
                'Longitudinal connection mismatch')
        roots = rec['roots']
        require(len(roots) == 3, 'Missing accepted sorted prefix')
        require(0 < roots[0]['omega'] < roots[1]['omega'] < roots[2]['omega'], 'Root ordering')
        for i, root in enumerate(roots, 1):
            require(root['sorted_position'] == i and root['counts'] == [i-1, i]
                    and all(root['gates'].values()), 'Sorted prefix quality failed')
            for field in ('omega', 'Lambda', 'Omega'):
                if field+str(i) in row:
                    close(float(row[field+str(i)]), root[field], 'CSV/root disagreement')
            close(root['omega'], epsilon*root['Lambda']**2, 'omega/Lambda normalization')
            close(root['Omega'], root['Lambda']**2, 'Omega/Lambda normalization')
        self.check_hash(path.with_suffix('.npz'), rec['shape_sha256'])
        if 'shape_sha256' in row:
            require(row['shape_sha256'] == rec['shape_sha256'], 'Shape provenance mismatch')
        if chi == 'RIGID':
            self.register(row, epsilon)
        return rec

    def register(self, row, epsilon):
        key = epsilon, float(row['beta_deg'])
        self.overlap.setdefault(key, []).append([float(row['omega'+str(k)]) for k in (1, 2, 3)])

    def unchanged(self):
        for path, digest in self.inputs.items():
            require(archive.sha(ROOT/path) == digest, 'Input mutated: '+path)


def thickness_data(saved):
    output, checks = [], []
    for name, eps, ref, expected_angle in [('growth', .005, 5, 6), ('thickness', .01, 10, 12)]:
        rows = ordered(r for r in saved.read(name, 'spectrum.csv') if r['ID'] == PAIR)
        profiles = ordered(r for r in saved.read(name, 'profiles.csv') if r['ID'] == PAIR
                           and (name == 'growth' or float(r['reference_end_deg']) == ref))
        beta, omega, lam = (values(rows, f) for f in ('beta_deg', 'omega1', 'Lambda1'))
        require(beta[-1] == 90 and len(rows) == (83 if name == 'growth' else 82), 'Incomplete q tail')
        close(values(profiles, 'beta_deg'), beta, 'Profile grid differs from spectrum')
        for r in rows:
            saved.point(r, eps)
        derivative = archive.stiffness_slopes if name == 'growth' else archive.thickness_slopes
        slope = derivative(beta, np.log(omega))[0]
        slope_lambda = derivative(beta, np.log(lam))[0]
        s0 = (np.log(omega[list(beta).index(ref)])-np.log(omega[0]))/ref
        require(s0 > 0, 'Nonpositive reference slope')
        q = np.abs(slope)/s0
        close(q, np.abs(slope_lambda)/(s0/2), 'q omega/Lambda identity')
        close(slope, values(profiles, 'signed_S'), 'Archived derivative not reproduced')
        close(q, values(profiles, 'q'), 'Archived q not reproduced')
        transition, tail = archive.refined_pick(beta, q, .15, 30)
        close(tail, values(profiles, 'Q_tail'), 'Archived tail not reproduced')
        old = next(r for r in saved.read(name, 'transitions.csv') if r['ID'] == PAIR
                   and (float(r['q_star']) == .15 if name == 'growth'
                        else float(r['reference_end_deg']) == ref))
        require(transition['beta_pl_deg'] == float(old['beta_pl_deg']) == expected_angle,
                'Characteristic angle differs from archived result')
        require(transition['previous_Q'] > .15 >= transition['Q_at_transition'], 'Tail crossing')
        checks.append(dict(epsilon=eps, reference_deg=ref, nodes_to_90=len(rows),
                           q_max_difference=float(max(abs(q-values(profiles, 'q')))),
                           **transition))
        for i, r in enumerate(rows):
            output.append(dict(epsilon=eps, beta_deg=beta[i], Lambda1=lam[i], omega1=omega[i],
                               reference_deg=ref, S0_per_degree=s0, slope_per_degree=slope[i],
                               q=q[i], Q_tail=tail[i], beta_pl_deg=expected_angle,
                               stencil=profiles[i]['stencil'], plotted=bool(beta[i] <= 20),
                               **saved.provenance(name, 'spectrum.csv'),
                               derivative_source_csv='profiles.csv', point_file=r['point_file']))
    return output, checks


def slip_data(saved):
    output, checks = [], []
    for chi, expected in [('RIGID', 170.73), ('0.01', 8.10), ('0', -.013)]:
        rows = ordered(r for r in saved.read('slip', 'spectrum.csv')
                       if r['pair'] == PAIR and r['chi'] == chi and float(r['beta_deg']) <= 10)
        require(float(rows[0]['beta_deg']) == 0 and float(rows[-1]['beta_deg']) == 10, 'Slip interval')
        omega0, lambda0 = float(rows[0]['omega1']), float(rows[0]['Lambda1'])
        for r in rows:
            saved.point(r, .005, chi)
            f = float(r['omega1'])/omega0
            close(f, (float(r['Lambda1'])/lambda0)**2, 'Slip frequency ratio')
            close(f, float(r['F_omega']), 'Slip saved normalization')
            output.append(dict(chi=chi, beta_deg=float(r['beta_deg']), omega1=float(r['omega1']),
                               Lambda1=float(r['Lambda1']), omega1_zero=omega0, F_omega=f,
                               **saved.provenance('slip', 'spectrum.csv'), point_file=r['point_file']))
        at3 = next(r for r in rows if float(r['beta_deg']) == 3)
        gain = (float(at3['omega1'])/omega0-1)*100
        close(gain, expected, 'Slip control at 3 deg', rtol=0, atol=.0005 if chi == '0' else .005)
        checks.append(dict(chi=chi, nodes=len(rows), gain_0_to_3_percent=gain))
    return output, checks


def leading_data(saved):
    leading, full, checks = [], [], []
    rows = ordered((r for r in saved.read('leading', 'spectrum.csv')
                    if r['ID'] == PAIR and float(r['xi']) <= 30), 'xi')
    zeta0 = float(rows[0]['zeta1'])
    for r in rows:
        z = [float(r['zeta'+str(k)]) for k in (1, 2, 3)]
        require(0 < z[0] < z[1] < z[2], 'Leading sorted positive prefix')
        f = math.sqrt(z[0]/zeta0)
        close(f, float(r['F_lim']), 'Leading frequency normalization')
        close(math.sqrt(z[0]), float(r['omega_lim1']), 'Leading eigenvalue units')
        leading.append(dict(xi=float(r['xi']), zeta1=z[0], zeta1_zero=zeta0, F_lim=f,
                            **saved.provenance('leading', 'spectrum.csv')))
    for eps, count, expected in [(.0025, 17, .0122), (.005, 33, .0498), (.01, 39, .1991)]:
        rows = ordered(r for r in saved.read('leading', 'comparison.csv')
                       if r['ID'] == PAIR and float(r['epsilon']) == eps)
        require(len(rows) == count and float(rows[0]['beta_deg']) == 0, 'Leading actual-node grid')
        omega0, lambda0 = float(rows[0]['omega1']), float(rows[0]['Lambda1'])
        errors = []
        for r in rows:
            require(r['interpolated'] == 'False', 'Interpolated full-model point')
            saved.point(r, eps)
            xi = math.radians(float(r['beta_deg']))/eps
            close(xi, float(r['xi']), 'xi coordinate', atol=6e-13)
            matched = [s for s in leading if abs(s['xi']-xi) < 1e-10]
            require(len(matched) == 1, 'Missing exact leading comparison node')
            flim = matched[0]['F_lim']
            ffull = float(r['omega1'])/omega0
            close(ffull, (float(r['Lambda1'])/lambda0)**2, 'Full frequency normalization')
            close(ffull, float(r['F_full']), 'Archived F_full')
            close(flim, float(r['F_lim']), 'Archived F_lim')
            error = ffull/flim-1
            close(error, float(r['E']), 'Error denominator must be F_lim')
            close(abs(error), float(r['abs_E']), 'Archived absolute error')
            errors.append(abs(error)*100)
            full.append(dict(epsilon=eps, beta_deg=float(r['beta_deg']), xi=xi,
                             omega1=float(r['omega1']), Lambda1=float(r['Lambda1']),
                             omega1_zero=omega0, F_full=ffull, F_lim=flim,
                             signed_delta_percent=error*100, delta_percent=abs(error)*100,
                             **saved.provenance('leading', 'comparison.csv'),
                             original_source_archive=r['source_bundle'], point_file=r['point_file'],
                             interpolated=False))
        close(max(errors), expected, 'Maximum leading discrepancy', rtol=0, atol=.00005)
        stats = next(r for r in saved.read('leading', 'error_statistics.csv')
                     if r['ID'] == PAIR and float(r['epsilon']) == eps and r['scope'] == 'all_actual')
        close(max(errors), float(stats['max_abs_E'])*100, 'Saved error statistics')
        checks.append(dict(epsilon=eps, nodes=count, max_xi=max(r['xi'] for r in full if r['epsilon'] == eps),
                           max_delta_percent=max(errors)))
    common = sorted(set.intersection(*[{round(r['xi'], 10) for r in full if r['epsilon'] == eps}
                                     for eps in (.0025, .005, .01)]))
    require(len(common) == 17, 'Common xi grid')
    for xi in common[1:]:
        e = [next(r['delta_percent'] for r in full if r['epsilon'] == eps and round(r['xi'], 10) == xi)
             for eps in (.0025, .005, .01)]
        require(e[0] < e[1] < e[2], 'Error decrease fails at common xi')
    return full, leading, dict(all_nodes=checks, common_nodes=len(common), common_max_xi=common[-1],
                               error_decreases_at_all_nonzero_common_nodes=True)


def variation_data(saved):
    rows = ordered(r for r in saved.read('variation', 'variation.csv')
                   if r['ID'] == PAIR and 3 <= float(r['beta_deg']) <= 10)
    require(list(values(rows, 'beta_deg')) == [3, 4, 5, 6, 8, 10], 'Variation angular nodes')
    sources = {float(r['beta_deg']): r for r in saved.read('variation', 'sources.csv') if r['ID'] == PAIR}
    material = saved.manifests['variation']['identity']['config']['material']
    model = archive.physics.Model(**material)
    output = []
    for r in rows:
        beta = float(r['beta_deg'])
        source = sources[beta]
        saved.point({**r, **source}, .005)
        # Integrate already saved quadrature arrays, without evaluating any fields.
        with np.load((ROOT/source['point_file']).with_suffix('.npz')) as data:
            arms = [{f: data[f'mode1_arm{i}_{f}'].copy() for f in ('x', 'weights', 'y', 'dy', 'ends')}
                    for i in (1, 2)]
        energy = archive.mechanics_energy(model, archive.physics.Case(beta, .01, .1, 10), arms)
        close(energy['mass'], 1, 'Full RL/Timoshenko modal mass', atol=5e-12)
        close(energy['mass'], float(r['mass']), 'Saved mass changed')
        close(energy['total'], float(r['V']), 'Saved potential energy changed')
        u, w, _, n, q, _ = arms[0]['ends'][1]
        close([u, w, n, q], [float(r[k]) for k in ('u1', 'w1', 'N1', 'Q1')], 'Joint state changed')
        d = math.hypot(u, w)
        close(d, float(r['geometry_norm']), 'Joint displacement/constraint derivative identity')
        slope = (n*w-q*u)/(2*energy['total'])
        close(slope, float(r['s_rad']), 'Variational derivative identity', rtol=1e-9)
        close(slope*math.pi/180, float(r['s_deg']), 'Radian/degree conversion', rtol=1e-9)
        output.append(dict(beta_deg=beta, u1=u, w1=w, d_joint=d, mass=energy['mass'],
                           mass_translation=float(energy['mass_components'][:, :2].sum()),
                           mass_rotation=float(energy['mass_components'][:, 2].sum()),
                           mass_RL=float(energy['mass_components'][:, 3].sum()),
                           V=energy['total'], N1=n, Q1=q, s_rad=float(r['s_rad']),
                           s_per_degree=float(r['s_deg']), S_per_radian=abs(float(r['s_rad'])),
                           R_d=d/math.hypot(float(rows[0]['u1']), float(rows[0]['w1'])),
                           R_S=abs(float(r['s_rad'])/float(rows[0]['s_rad'])),
                           **saved.provenance('variation', 'variation.csv'), point_file=source['point_file']))
    d_ratio = output[0]['d_joint']/output[-1]['d_joint']
    s_ratio = output[0]['S_per_radian']/output[-1]['S_per_radian']
    close(d_ratio, 8.1, 'Displacement control', rtol=0, atol=.05)
    close(s_ratio, 38.5, 'Slope control', rtol=0, atol=.05)
    return output, dict(nodes=list(values(rows, 'beta_deg')), displacement_reduction=d_ratio,
                        slope_reduction=s_ratio, max_mass_error=max(abs(r['mass']-1) for r in output))


def prepare():
    saved = SavedData()
    # Independent overlaps, including actual epsilon nodes not used in figure 1.
    for r in saved.read('epsilon', 'spectrum.csv'):
        if r['pair'] == PAIR:
            saved.point(r, float(r['epsilon']))
    thickness, c1 = thickness_data(saved)
    slip, c2 = slip_data(saved)
    full, leading, c3 = leading_data(saved)
    variation, c4 = variation_data(saved)
    discrepancies = [float(np.max(np.abs(np.array(a)/np.array(b)-1)))
                     for group in saved.overlap.values() for a, b in itertools.combinations(group, 2)]
    require(discrepancies and max(discrepancies) <= 1e-8, 'Overlapping spectra disagree')
    tables = dict(zip(['01_thickness_frequency_and_q.csv', '02_longitudinal_compliance.csv',
                       '03_full_model.csv', '03_leading_model.csv', '04_joint_displacement_and_slope.csv'],
                      [thickness, slip, full, leading, variation]))
    validation = dict(status='PASS', figure1=c1, figure2=c2, figure3=c3, figure4=c4,
                      overlap_prefix_comparisons=len(discrepancies),
                      maximum_overlap_relative_difference=max(discrepancies),
                      unique_verified_full_model_points=len(saved.points),
                      full_model_solves=0, leading_model_solves=0, matrix_assemblies=0,
                      field_reconstructions=0, interpolated_spectral_points=0)
    saved.unchanged()
    return tables, validation, saved


def render(tables, folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'mathtext.fontset': 'dejavusans', 'axes.labelsize': 11,
                         'legend.fontsize': 9, 'pdf.fonttype': 42, 'axes.linewidth': .7,
                         'lines.linewidth': 1.5, 'savefig.dpi': 360})
    comma = FuncFormatter(lambda x, pos: f'{0 if abs(x) < 1e-12 else x:g}'.replace('.', ','))
    blue, orange, green = '#245580', '#AA4C20', '#27705E'

    def style(ax, panel=None, xlabel=r'$\beta$, град.'):
        ax.set_xlabel(xlabel)
        ax.xaxis.set_major_formatter(comma)
        ax.yaxis.set_major_formatter(comma)
        ax.grid(color='.86', linewidth=.5)
        ax.set_axisbelow(True)
        ax.tick_params(direction='out', length=3, width=.7)
        if panel:
            ax.text(.015, .97, panel, transform=ax.transAxes, ha='left', va='top', fontsize=11)

    def save(fig, name):
        fig.savefig(folder/(name+'.png'), dpi=360)
        fig.savefig(folder/(name+'.pdf'), metadata={'CreationDate': None, 'ModDate': None})
        plt.close(fig)

    # Two full-width panels keep the criterion readable on A4; the inset uses
    # the same actual nodes and the same linear scale, without smoothing.
    fig, axes = plt.subplots(2, 1, figsize=(6.7, 6.3), sharex=True, layout='constrained')
    inset = axes[1].inset_axes([.48, .34, .49, .60])
    for eps, color, ls, marker, label in [(.005, blue, '-', 'o', r'$\varepsilon=0{,}005$'),
                                          (.01, orange, '--', 's', r'$\varepsilon=0{,}01$')]:
        rows = [r for r in tables['01_thickness_frequency_and_q.csv'] if r['epsilon'] == eps and r['plotted']]
        b, q = values(rows, 'beta_deg'), values(rows, 'q')
        axes[0].plot(b, values(rows, 'Lambda1'), color=color, ls=ls, marker=marker, ms=2.4, mfc='white', label=label)
        # The coarse one-sided estimate at zero is shown separately, not as an
        # analytic derivative. All archived values remain in the exported CSV.
        mask = b > 0 if eps == .005 else np.ones(len(b), dtype=bool)
        axes[1].plot(b[mask], q[mask], color=color, ls=ls, marker=marker, ms=2.4, mfc='white')
        if eps == .005:
            axes[1].plot(b[0], q[0], color=color, marker='x', ms=6, ls='none', clip_on=False)
        inset.plot(b, q, color=color, ls=ls, marker=marker, ms=2.5, mfc='white')
        angle = rows[0]['beta_pl_deg']
        for ax in [*axes, inset]:
            ax.axvline(angle, color=color, ls=':', lw=1)
        axes[0].text(angle+.2, 1.98, rf'$\beta_{{\mathrm{{pl}}}}={angle:g}^\circ$',
                     color=color, fontsize=9)
    for ax, panel in zip(axes, ('а', 'б')):
        style(ax, panel)
        ax.set_xlim(0, 20)
        ax.xaxis.set_major_locator(MultipleLocator(2))
    axes[0].set_xlabel('')
    axes[0].set_ylabel(r'$\Lambda_1$')
    axes[0].set_ylim(1.7, 3.85)
    axes[0].legend(loc='lower right', framealpha=1)
    axes[1].set_ylabel(r'$q(\beta)$')
    axes[1].set_ylim(-.04, 2.15)
    axes[1].axhline(.15, color='.2', ls='-.', lw=1)
    axes[1].text(.03, .14, r'$q_*=0{,}15$', transform=axes[1].transAxes, fontsize=9)
    style(inset, xlabel='')
    inset.set_xlim(4.5, 13.5)
    inset.set_ylim(.05, .30)
    inset.set_xticks([6, 8, 10, 12])
    inset.set_yticks([.1, .15, .2, .25, .3])
    inset.tick_params(labelsize=8)
    inset.axhline(.15, color='.2', ls='-.', lw=1)
    save(fig, FIGURES[0])

    fig, ax = plt.subplots(figsize=(6.7, 4.3), layout='constrained')
    inset = ax.inset_axes([.50, .34, .46, .38])
    for chi, color, ls, marker, label in [('RIGID', blue, '-', 'o', r'Точная связь $c_1=0$'),
                                          ('0.01', orange, '--', 's', r'$\chi=0{,}01$'),
                                          ('0', green, '-.', '^', r'$\chi=0$')]:
        rows = [r for r in tables['02_longitudinal_compliance.csv'] if r['chi'] == chi]
        b, f = values(rows, 'beta_deg'), values(rows, 'F_omega')
        ax.plot(b, f, color=color, ls=ls, marker=marker, ms=3, mfc='white', label=label)
        if chi != 'RIGID':
            inset.plot(b, f, color=color, ls=ls, marker=marker, ms=2.7, mfc='white')
    style(ax)
    ax.set_xlim(0, 10)
    ax.set_ylim(.90, 3.95)
    ax.set_ylabel(r'$\omega_1(\beta)/\omega_1(0)$')
    ax.legend(loc='upper left', framealpha=1)
    style(inset, xlabel='')
    inset.set_xlim(0, 3)
    inset.set_ylim(.995, 1.09)
    inset.set_xticks([0, 1, 2, 3])
    inset.set_yticks([1, 1.04, 1.08])
    inset.tick_params(labelsize=8)
    inset.set_title('Увеличенный фрагмент', fontsize=8, pad=4)
    save(fig, FIGURES[1])

    fig, axes = plt.subplots(2, 1, figsize=(6.7, 6.4), sharex=True, layout='constrained')
    leading = tables['03_leading_model.csv']
    axes[0].plot(values(leading, 'xi'), values(leading, 'F_lim'), color='.15', lw=1.4,
                 label=r'Ведущая модель $F_{\mathrm{lim}}$')
    for eps, color, marker, ms in [(.0025, blue, 'o', 5.5), (.005, orange, 's', 4.0), (.01, green, '^', 2.9)]:
        rows = [r for r in tables['03_full_model.csv'] if r['epsilon'] == eps]
        label = r'$\varepsilon='+str(eps).replace('.', '{,}')+'$'
        axes[0].plot(values(rows, 'xi'), values(rows, 'F_full'), ls='none', marker=marker,
                     ms=ms, mfc='none', mew=.8, color=color, label=label)
        axes[1].plot(values(rows, 'xi'), values(rows, 'delta_percent'), ls='none',
                     marker=marker, ms=ms, mfc='none', mew=.8, color=color)
    for ax, panel in zip(axes, ('а', 'б')):
        style(ax, panel, r'$\xi=\beta_{\mathrm{rad}}/\varepsilon$')
        ax.set_xlim(0, 30)
    axes[0].set_xlabel('')
    axes[0].set_ylabel(r'$F_{\mathrm{full}},\ F_{\mathrm{lim}}$')
    axes[0].set_ylim(.88, 4)
    axes[0].legend(loc='lower right', framealpha=1, fontsize=8.5)
    axes[1].set_ylabel(r'$\delta$, %')
    axes[1].set_ylim(-.006, .212)
    axes[1].set_yticks([0, .05, .10, .15, .20])
    save(fig, FIGURES[2])

    fig, axes = plt.subplots(2, 1, figsize=(6.7, 5.5), sharex=True, layout='constrained')
    rows = tables['04_joint_displacement_and_slope.csv']
    for ax, field, label, color, marker, panel in zip(
            axes, ('R_d', 'R_S'), (r'$R_d(\beta)$', r'$R_{\mathcal{S}}(\beta)$'),
            (blue, orange), ('o', 's'), ('а', 'б')):
        ax.plot(values(rows, 'beta_deg'), values(rows, field), color=color, marker=marker,
                ms=5, mfc='white', lw=1.3)
        style(ax, panel)
        ax.texts[-1].set_position((.06, .97))
        ax.set_ylabel(label)
        ax.set_xlim(2.8, 10.2)
        ax.set_ylim(0, 1.08)
        ax.set_xticks(range(3, 11))
    axes[0].set_xlabel('')
    save(fig, FIGURES[3])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--check-only', action='store_true', help='Validate sources and controls without writing')
    args = parser.parse_args()
    started = time.perf_counter()
    tables, validation, saved = prepare()
    if args.check_only:
        print(json.dumps(validation, ensure_ascii=False, indent=2))
        return
    folder = args.output.resolve()
    # A mistaken --output must never overwrite any existing spectral archive.
    require(folder.is_relative_to(ROOT) and folder != ROOT, 'Output must be inside the workspace')
    for relative, _ in SOURCES.values():
        source = (BASE/relative).resolve()
        require(not folder.is_relative_to(source) and not source.is_relative_to(folder), 'Output overlaps source archive')
    require(not any(p.name in {'manifest.json', 'spectrum.csv', 'contract.json'} for p in folder.glob('*'))
            or folder == OUTPUT.resolve() or (folder/'illustrations_manifest.json').exists(),
            'Output is an unrelated existing result bundle')
    folder.mkdir(parents=True, exist_ok=True)
    for filename, rows in tables.items():
        archive.csvwrite(folder/filename, rows)
    render_started = time.perf_counter()
    render(tables, folder)
    render_seconds = time.perf_counter()-render_started
    saved.unchanged()
    archive.write(folder/'validation.json', validation)
    template = Path(__file__).with_suffix('.readme.txt')
    (folder/'README.md').write_text(template.read_text(encoding='utf8'), encoding='utf8')
    outputs = [*tables, 'validation.json', 'README.md',
               *(name+ext for name in FIGURES for ext in ('.pdf', '.png'))]
    archive.write(folder/'illustrations_manifest.json', dict(
        schema='two-springs-report-illustrations-v1', status='PASS', policy=POLICY,
        parameters=dict(mu=.01, aspect=4, nu=.3, ks=5/6, kappa_joint=.1, kappa_end=10),
        inputs=saved.inputs, code_sha256=archive.sha(Path(__file__)),
        helper_sha256=archive.sha(Path(archive.__file__)), template_sha256=archive.sha(template),
        outputs={name: archive.sha(folder/name) for name in outputs},
        full_model_solves=0, leading_model_solves=0, field_reconstructions=0,
        spectrum_generation_seconds=0, rendering_seconds=render_seconds,
        total_seconds=time.perf_counter()-started, png_dpi=360,
        source_preservation='All consumed source files rehashed after rendering'))
    print(f'PASS: four PDF/PNG figures in {folder}; 0 full/leading eigenvalue solves')


if __name__ == '__main__':
    main()
