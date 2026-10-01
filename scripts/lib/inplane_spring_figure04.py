"""Figure-only EB/RLB overlay from K15/K22 elastic shapes and K23 context.

No solver/recovery/tracking imports or calls. Existing entry point:
check_inplane_spring_robustness.py plot-only --figure04.
"""
import csv
import json
import time
from pathlib import Path

import numpy as np

from scripts.lib.inplane_spring_figure03 import global_centrelines, read_json, sha

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'results/laminated_beams'
OUTPUT = BASE/'figure04_eb_rlb_shapes'
FOLDERS = {
    'EB': 'inplane_kelvin_voigt_elastic_screening',
    'RLB': 'inplane_kelvin_voigt_rlb_elastic_screening',
    'K23': 'inplane_kelvin_voigt_eb_rlb_complex_confirmation',
}
TABLES = {'EB': 'elastic_screening.csv', 'RLB': 'rlb_elastic_screening.csv'}
SHAPES = {'EB': 'screening_shapes.npz', 'RLB': 'rlb_elastic_shapes.npz'}
MODES = {'EB': 'sorted_05', 'RLB': 'rlb_sorted_05'}
BETAS = (45, 75)
THEORIES = ('EB', 'RLB')
COLORS = {'EB': '#0072B2', 'RLB': '#D55E00'}
STYLES = {'EB': '-', 'RLB': (0, (5, 2.6))}
DISPLAY_AMPLITUDE = .14


def csv_rows(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def source_paths(base=BASE):
    return [base/FOLDERS[t]/name for t in THEORIES
            for name in (TABLES[t], SHAPES[t], 'run_manifest.json')] + [
        base/FOLDERS['RLB']/'eb_rlb_matched_comparison.csv',
        base/FOLDERS['K23']/'representative_complex_comparison.csv',
    ]


def weights(nodes):
    """The existing uniform composite Simpson quadrature, for read-only checks."""
    assert nodes == 129
    w = np.ones(nodes)
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    return w/(3*(nodes-1))


def saved_mass(states, config, theory):
    prop = config['properties']
    density = prop['m']*np.sum(states[:, :, :2]**2, axis=2)
    if theory == 'RLB':
        density += prop['J']*states[:, :, 2]**2
    return float(np.sum(density*weights(states.shape[1])[None, :]
                        *np.array([config['L1'], config['L2']])[:, None]))


def load_saved(base=BASE):
    configs = {t: read_json(base/FOLDERS[t]/'run_manifest.json')['configuration'] for t in THEORIES}
    assert {k: v for k, v in configs['EB'].items() if k != 'model'} == {
        k: v for k, v in configs['RLB'].items() if k != 'model'}
    config = configs['EB']
    assert config['mu'] == config['d_theta'] == 0
    assert config['L1'] == config['L2'] == config['kappa_theta'] == 1
    assert config['layup'] == 'H/L/L/H' and config['contrast'] == .4
    matching = csv_rows(base/FOLDERS['RLB']/'eb_rlb_matched_comparison.csv')
    context = csv_rows(base/FOLDERS['K23']/'representative_complex_comparison.csv')
    selected, matches, contexts, checks = {}, {}, {}, {}
    for beta in BETAS:
        match, = [r for r in matching if float(r['beta_deg']) == beta and r['eb_sorted_mode'] == MODES['EB']]
        assert match['rlb_sorted_mode'] == MODES['RLB']
        assert match['eta'] == '1' and match['match_status'] == 'CONFIRMED'
        assert match['comparison_status'] == 'ACTIVE_MATCHED'
        actual, = [r for r in context if float(r['beta_deg']) == beta and r['eb_mode'] == MODES['EB']]
        assert actual['rlb_mode'] == MODES['RLB'] and float(actual['d_theta']) == .001
        assert actual['eb_root_status'] == actual['rlb_root_status'] == 'ROOT_ACCEPTED'
        matches[beta], contexts[beta] = match, actual
    for theory in THEORIES:
        folder = base/FOLDERS[theory]
        rows = csv_rows(folder/TABLES[theory])
        mode_column = 'sorted_mode' if theory == 'EB' else 'rlb_sorted_mode'
        with np.load(folder/SHAPES[theory], allow_pickle=False) as archive:
            for beta in BETAS:
                row, = [r for r in rows if float(r['beta_deg']) == beta and r[mode_column] == MODES[theory]]
                assert row['root_status'] == 'CONFIRMED' and row['activity_status'] == 'ACTIVE'
                assert int(row.get('symmetry_eta', row.get('eta'))) == 1
                assert float(row['d_theta']) == 0 and float(row['kappa_theta']) == 1
                assert float(row['Omega']) == float(matches[beta]['Omega_'+theory])
                assert float(row['zeta_slope_pred']) == float(contexts[beta]['G_'+theory])
                states = archive[row['shape_key']+'__states'].copy()
                assert states.shape == (2, 129, 6) and np.isrealobj(states)
                assert np.isfinite(states).all()
                np.testing.assert_allclose(states[:, -1, 2],
                    [float(row['psi1_joint']), float(row['psi2_joint'])], rtol=1e-12, atol=1e-12)
                # Exact eta=+1 lift; no projection or new matching is performed.
                np.testing.assert_allclose(states[1], states[0]*np.array([1, -1, -1, 1, -1, -1]), atol=1e-11)
                mass = saved_mass(states, configs[theory], theory)
                assert abs(mass-1) < 1e-10 and abs(float(row['mass_M'])-mass) < 1e-10
                # A global sign only. Keep a positive EB endpoint rotation and align
                # its already matched RLB shape using the fixed u,w comparison metric.
                if theory == 'EB':
                    sign = 1. if states[0, -1, 2] > 0 else -1.
                else:
                    overlap = np.sum(states[:, :, :2]*selected[(beta, 'EB')]['states'][:, :, :2]
                                     *weights(129)[None, :, None])
                    assert overlap != 0
                    sign = 1. if overlap > 0 else -1.
                selected[(beta, theory)] = dict(row=row, states=sign*states, sign=sign)
                checks[f'{beta}_{theory}'] = dict(mass_check=mass, mass_source=float(row['mass_M']), display_sign=sign)
    maximum = max(np.linalg.norm(item['states'][:, :, :2], axis=2).max() for item in selected.values())
    return dict(data=selected, matching=matches, context=contexts, checks=checks,
                display_scale=DISPLAY_AMPLITUDE/maximum, configurations=configs)


def figure_rows(saved):
    rows = []
    for beta in BETAS:
        for theory in THEORIES:
            item = saved['data'][(beta, theory)]
            row = item['row']
            rows.append(dict(beta_deg=beta, theory=theory, mode=MODES[theory], eta=1,
                Omega=row['Omega'], Lambda=row['Lambda'], mass_M=row['mass_M'],
                psi_joint=row['psi1_joint'], Delta_psi=row['Delta_psi'], G=row['zeta_slope_pred'],
                zeta_at_reference_d=saved['context'][beta]['zeta_'+theory], reference_d_theta=.001,
                elastic_d_theta=0, display_sign=item['sign'],
                psi_joint_display=float(item['states'][0, -1, 2]),
                Delta_psi_display=float(item['states'][0, -1, 2]-item['states'][1, -1, 2]),
                matching_MAC=saved['matching'][beta]['MAC_common'],
                source=f"{FOLDERS[theory]}/{TABLES[theory]}#beta={beta},{MODES[theory]};"
                       f"{FOLDERS[theory]}/{SHAPES[theory]}#{row['shape_key']}__states",
                zeta_source=f"{FOLDERS['K23']}/representative_complex_comparison.csv#{saved['context'][beta]['case_id']}"))
    return rows


def draw(saved):
    """Original diagnostic layout; retained to preserve/reuse the top row."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), gridspec_kw={'height_ratios': [1, 1.12]})
    fig.subplots_adjust(left=.085, right=.975, top=.91, bottom=.09, hspace=.34, wspace=.25)
    for col, beta in enumerate(BETAS):
        top, bottom = axes[:, col]
        for theory in THEORIES:
            states = saved['data'][(beta, theory)]['states']
            reference, disp = global_centrelines(beta, 0., states)
            deformed = reference+saved['display_scale']*disp
            for arm in range(2):
                if theory == 'EB':
                    top.plot(*reference[arm].T, color='.76', lw=.9)
                    end = reference[arm, 0]
                    t = reference[arm, -1]-end; t /= np.linalg.norm(t)
                    normal = np.array([-t[1], t[0]])
                    top.plot(*(end+np.array([-.05, .05])[:, None]*normal).T, color='.45', lw=1)
                    for offset in (-.035, 0., .035):
                        start = end+offset*normal
                        top.plot(*np.array([start, start-.027*t-.016*normal]).T, color='.45', lw=.65)
                top.plot(*deformed[arm].T, color=COLORS[theory], linestyle=STYLES[theory], lw=1.9,
                         label=f'{theory}_{arm}')
            # Use stored state psi: EB has psi=-dw/dx exactly in its recovery
            # contract; RLB psi is independent. Never differentiate sampled w.
            xi = np.linspace(0, 1, states.shape[1])
            bottom.plot(xi, states[0, :, 2], color=COLORS[theory], linestyle=STYLES[theory],
                        lw=1.9, label=theory)
            bottom.plot(1, states[0, -1, 2], marker='o' if theory == 'EB' else 's',
                        ms=5, mfc='white', mec=COLORS[theory], mew=1.2, zorder=4)
        top.set_aspect('equal'); top.set_xlim(-1.06, 1.06); top.set_ylim(-.20, .79)
        top.axis('off')
        top.set_title(r'$\beta='+str(beta)+r'^{\circ}$', fontsize=13, pad=8)
        bottom.set_xlabel(r'$\xi$', fontsize=13)
        bottom.set_ylabel(r'$\psi$', fontsize=13, rotation=0, labelpad=14)
        bottom.set_xlim(0, 1.035)
        bottom.xaxis.set_major_locator(MultipleLocator(.25))
        bottom.axhline(0, color='.8', lw=.65, zorder=0)
        bottom.spines[['top', 'right']].set_visible(False)
        bottom.grid(axis='y', color='.92', lw=.6)
    fig.legend(handles=[Line2D([], [], color=COLORS[t], linestyle=STYLES[t], lw=2, label=t)
                        for t in THEORIES], loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(.53, .995), handlelength=2.7, columnspacing=2)
    return fig


def draw_revised(saved, rows, *, value_labels=True):
    """Keep v1 centrelines verbatim; replace only the lower axes with |Delta psi|."""
    from matplotlib.ticker import MultipleLocator
    assert len(rows) == 4
    by_key = {(int(r['beta_deg']), r['theory']): r for r in rows}
    assert set(by_key) == {(b, t) for b in BETAS for t in THEORIES}
    values = {key: abs(float(row['Delta_psi'])) for key, row in by_key.items()}
    limit = 5*np.ceil(1.12*max(values.values())/5)
    fig = draw(saved)
    for col, beta in enumerate(BETAS):
        ax = fig.axes[2+col]
        ax.clear()
        for x, theory in enumerate(THEORIES):
            value = values[(beta, theory)]
            ax.bar(x, value, width=.36, color=COLORS[theory], edgecolor=COLORS[theory],
                   label=theory, zorder=3)
            if value_labels:
                ax.annotate(f'{value:.1f}', (x, value), xytext=(0, 5), textcoords='offset points',
                            ha='center', va='bottom', fontsize=11)
        ax.set_xticks([0, 1], THEORIES)
        ax.set_xlim(-.6, 1.6)
        ax.set_ylim(0, limit)
        ax.set_ylabel(r'$|\Delta\psi|$', fontsize=13, labelpad=10)
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='.92', lw=.6, zorder=0)
        ax.tick_params(axis='x', length=0, pad=8)
    return fig


def render(base=BASE, output=OUTPUT, *, value_labels=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    started = time.perf_counter()
    source_output = base/'figure04_eb_rlb_shapes'
    old_files = [source_output/name for name in ('figure04_data.csv', 'figure04_eb_rlb_shapes.png',
                                                'figure04_eb_rlb_shapes.pdf', 'figure_manifest.json')]
    hashes = {str(p): sha(p) for p in source_paths(base)+old_files}
    saved = load_saved(base)
    rows = csv_rows(source_output/'figure04_data.csv')
    # The original table is read-only, including signed rotations and all contextual fields.
    assert rows == [{k: str(v) for k, v in row.items()} for row in figure_rows(saved)]
    output.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42}):
        fig = draw_revised(saved, rows, value_labels=value_labels)
        for extension in ('png', 'pdf'):
            fig.savefig(output/f'figure04_eb_rlb_shapes_revised.{extension}', dpi=300, facecolor='white',
                        metadata={'Creator': 'CoupledBeams Figure 4 plot-only'})
        plt.close(fig)
    assert all(sha(Path(p)) == h for p, h in hashes.items())
    record = dict(source_sha256=hashes, matching=saved['matching'], K23_context=saved['context'],
        shape_checks=saved['checks'], beta_deg=BETAS, modes=MODES,
        display_scale=saved['display_scale'], maximum_display_displacement=DISPLAY_AMPLITUDE,
        lower_panels='abs(Delta_psi) copied from unchanged figure04_data.csv',
        common_y_scale=True, value_labels=value_labels, original_figure_preserved=True,
        normalization='saved whole-structure physical mass M=1; no rescaling except a global sign',
        reused_elastic_shapes=4, new_roots=0, new_beta=0, new_d=0, matrix_calls=0,
        shape_recoveries=0, tracking_calls=0, interpolation=0, solver_changes=0,
        render_seconds=time.perf_counter()-started)
    (output/'figure_manifest_revised.json').write_text(json.dumps(record, ensure_ascii=False, indent=2,
                                                        allow_nan=False)+'\n', encoding='utf-8')
    return {k: record[k] for k in ('reused_elastic_shapes', 'new_roots', 'shape_recoveries', 'render_seconds')}
