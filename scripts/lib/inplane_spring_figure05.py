"""Figure 5: decay envelopes and historical corrections from saved K23 rows."""
import csv
import hashlib
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'results/laminated_beams/inplane_kelvin_voigt_eb_rlb_complex_confirmation/representative_complex_comparison.csv'
OUTPUT = ROOT/'results/laminated_beams/figure05_frequency_damping_comparison'
CASES = (('R1', 0, '05'), ('R2', 45, '05'), ('R3', 75, '05'), ('R4', 45, '03'))
COLORS = ('#4C6E91', '#4C8478')
STEM = 'figure05_frequency_damping_comparison'
ENVELOPE_CASES = (('R2', 45, '05'), ('R3', 75, '05'))
ENVELOPE_OUTPUT = ROOT/'results/laminated_beams/figure05_damping_envelopes'
ENVELOPE_STEM = 'figure05_damping_envelopes'
ENVELOPE_COLORS = {'EB': '#0072B2', 'RLB': '#D55E00'}
ENVELOPE_STYLES = {'EB': '-', 'RLB': (0, (5, 2.6))}
PERIOD_LIMITS = {45: 250, 75: 7500}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rows(source=SOURCE):
    with source.open(encoding='utf-8', newline='') as stream:
        raw = list(csv.DictReader(stream))
    rows = []
    for case_id, beta, mode in CASES:
        item, = [r for r in raw if r['case_id'] == case_id]
        assert float(item['beta_deg']) == beta and float(item['d_theta']) == .001
        assert item['eb_mode'] == 'sorted_'+mode and item['rlb_mode'] == 'rlb_sorted_'+mode
        assert item['eb_root_status'] == item['rlb_root_status'] == 'ROOT_ACCEPTED'
        values = {k: float(item[k]) for k in ('Omega_d_EB', 'Omega_d_RLB', 'zeta_EB', 'zeta_RLB')}
        assert all(math.isfinite(v) and v > 0 for v in values.values())
        delta_omega = 100*(values['Omega_d_RLB']-values['Omega_d_EB'])/values['Omega_d_EB']
        delta_zeta = 100*(values['zeta_RLB']-values['zeta_EB'])/values['zeta_EB']
        assert math.isclose(delta_omega, 100*float(item['damped_delta_Omega']), rel_tol=0, abs_tol=1e-12)
        assert math.isclose(delta_zeta, 100*float(item['complex_delta_zeta']), rel_tol=0, abs_tol=1e-12)
        # Keep original numeric strings for dimensional comparisons; only the
        # two explicitly requested percent transformations are evaluated here.
        rows.append(dict(case_id=case_id, beta_deg=beta, mode=mode, d_theta=item['d_theta'],
            Omega_d_EB=item['Omega_d_EB'], Omega_d_RLB=item['Omega_d_RLB'],
            delta_Omega_d_percent=delta_omega, zeta_EB=item['zeta_EB'], zeta_RLB=item['zeta_RLB'],
            delta_zeta_percent=delta_zeta, source=source.as_posix()+'#'+case_id))
    return rows


def draw(rows):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    assert [r['case_id'] for r in rows] == [c[0] for c in CASES]
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10, 4.5))
    fig.subplots_adjust(left=.12, right=.975, bottom=.19, top=.88, wspace=.25)
    labels = [rf"${r['beta_deg']}^{{\circ}}\,/\,{r['mode']}$" for r in rows]
    specs = (
        ('delta_Omega_d_percent', r'$\Delta\Omega_d/\Omega_{d,\mathrm{EB}},\ \%$', (-5.7, .55), 1),
        ('delta_zeta_percent', r'$\Delta\zeta/\zeta_{\mathrm{EB}},\ \%$', (-96, 22), 20),
    )
    for index, (ax, color, (key, xlabel, limits, step)) in enumerate(zip(axes, COLORS, specs)):
        values = [r[key] for r in rows]
        ax.barh(range(4), values, height=.32, color=color, zorder=3)
        ax.axvline(0, color='.45', lw=.8, zorder=2)
        for y, value in enumerate(values):
            ax.annotate(f'{value:+.2f}%'.replace('-', '\N{MINUS SIGN}'), (value, y),
                        xytext=(-6 if value < 0 else 6, 0), textcoords='offset points',
                        ha='right' if value < 0 else 'left', va='center', fontsize=11)
        ax.set_xlim(*limits)
        ax.set_ylim(3.6, -.6)
        ax.set_yticks(range(4), labels)
        ax.set_xlabel(xlabel, fontsize=13, labelpad=9)
        ax.xaxis.set_major_locator(MultipleLocator(step))
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='y', length=0, pad=12)
        ax.grid(axis='x', color='.91', lw=.65, zorder=0)
        ax.text(0, 1.08, '(a)' if index == 0 else '(b)', transform=ax.transAxes, fontsize=13)
    return fig


def render(source=SOURCE, output=OUTPUT):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    started = time.perf_counter()
    before = sha(source)
    rows = load_rows(source)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'figure05_data.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42}):
        fig = draw(rows)
        for extension in ('png', 'pdf'):
            fig.savefig(output/f'{STEM}.{extension}', dpi=300, facecolor='white',
                        metadata={'Creator': 'CoupledBeams Figure 5 plot-only'})
        plt.close(fig)
    assert sha(source) == before
    manifest = dict(source=str(source), source_sha256=before, selected_cases=CASES,
        quantity_source='actual K23 complex Omega_d and zeta; no elastic predictors',
        percent_definitions='100*(RLB-EB)/EB, separately for Omega_d and zeta',
        reused_comparison_rows=4, new_roots=0, new_forms=0, new_beta=0, new_d=0,
        matrix_calls=0, shape_recoveries=0, interpolation=0, solver_changes=0, high_precision=0,
        output_files=[STEM+'.png', STEM+'.pdf', 'figure05_data.csv'],
        render_seconds=time.perf_counter()-started)
    (output/'figure05_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False,
                                                         indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return {k: manifest[k] for k in ('reused_comparison_rows', 'new_roots', 'new_forms', 'render_seconds')}


def decay_rate_per_period(zeta):
    """Exact alpha*T_d; no small-zeta approximation."""
    zeta = float(zeta)
    assert math.isfinite(zeta) and 0 < zeta < 1
    return 2*math.pi*zeta/math.sqrt(1-zeta*zeta)


def envelope(periods, zeta):
    """Analytical evaluation of one mode's normalized envelope, not a new solve."""
    assert math.isfinite(periods) and periods >= 0
    return math.exp(-decay_rate_per_period(zeta)*periods)


def load_envelope_rows(source=SOURCE):
    with source.open(encoding='utf-8', newline='') as stream:
        raw = list(csv.DictReader(stream))
    manifest_path = source.parent/'run_manifest.json'
    tref = json.loads(manifest_path.read_text(encoding='utf-8'))['configuration']['t_ref']
    assert math.isfinite(tref) and tref > 0
    rows = []
    for case_id, beta, mode in ENVELOPE_CASES:
        item, = [r for r in raw if r['case_id'] == case_id]
        assert float(item['beta_deg']) == beta and float(item['d_theta']) == .001
        assert item['eb_mode'] == 'sorted_'+mode and item['rlb_mode'] == 'rlb_sorted_'+mode
        assert item['eb_root_status'] == item['rlb_root_status'] == 'ROOT_ACCEPTED'
        for theory in ('EB', 'RLB'):
            zeta = float(item['zeta_'+theory]); Omega = float(item['Omega_d_'+theory])
            rate = decay_rate_per_period(zeta)
            assert math.isfinite(Omega) and Omega > 0
            a = zeta*Omega/math.sqrt(1-zeta*zeta)
            rows.append(dict(case_id=case_id, beta_deg=beta, mode=mode, theory=theory,
                d_theta=item['d_theta'], Omega_d=item['Omega_d_'+theory], zeta=item['zeta_'+theory],
                alpha=a/tref, N_half=math.log(2)/rate, a=a, omega_d=Omega/tref, t_ref=tref,
                source=source.as_posix()+'#'+case_id, normalization_source=manifest_path.as_posix()))
    return rows


def draw_envelopes(rows):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator
    by_key = {(r['beta_deg'], r['theory']): r for r in rows}
    assert len(rows) == len(by_key) == 4
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10, 4.5))
    fig.subplots_adjust(left=.085, right=.98, bottom=.17, top=.80, wspace=.22)
    for col, (_, beta, _) in enumerate(ENVELOPE_CASES):
        ax = axes[col]; limit = PERIOD_LIMITS[beta]
        # Sampling a known exponential for display: no spectral interpolation.
        periods = [limit*i/500 for i in range(501)]
        for theory in ('EB', 'RLB'):
            row = by_key[(beta, theory)]
            ax.plot(periods, [envelope(n, row['zeta']) for n in periods],
                    color=ENVELOPE_COLORS[theory], linestyle=ENVELOPE_STYLES[theory], lw=2, label=theory)
            ax.plot(row['N_half'], .5, marker='o' if theory == 'EB' else 's', ms=5.5,
                    mfc='white', mec=ENVELOPE_COLORS[theory], mew=1.3, zorder=4)
        ax.axhline(.5, color='.65', ls=':', lw=.8, zorder=0)
        ax.set_xlim(0, limit); ax.set_ylim(0, 1.03)
        ax.xaxis.set_major_locator(MultipleLocator(50 if beta == 45 else 1500))
        ax.yaxis.set_major_locator(MultipleLocator(.25))
        ax.set_xlabel(r'$N$', fontsize=13)
        ax.set_title(r'$\beta='+str(beta)+r'^{\circ}$', fontsize=13, pad=13)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0, 1.075, '(a)' if col == 0 else '(b)', transform=ax.transAxes, fontsize=13)
    axes[0].set_ylabel(r'$A/A_0$', fontsize=13, labelpad=9)
    fig.legend(handles=[Line2D([], [], color=ENVELOPE_COLORS[t], linestyle=ENVELOPE_STYLES[t], lw=2, label=t)
                        for t in ('EB', 'RLB')], loc='upper center', ncol=2, frameon=False,
               bbox_to_anchor=(.54, .995), handlelength=2.8)
    return fig


def render_envelopes(source=SOURCE, output=ENVELOPE_OUTPUT):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    started = time.perf_counter()
    hashes = {str(p): sha(p) for p in (source, source.parent/'run_manifest.json')}
    rows = load_envelope_rows(source)
    output.mkdir(parents=True, exist_ok=True)
    with (output/'figure05_data.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42}):
        fig = draw_envelopes(rows)
        for extension in ('png', 'pdf'):
            fig.savefig(output/f'{ENVELOPE_STEM}.{extension}', dpi=300, facecolor='white',
                        metadata={'Creator': 'CoupledBeams Figure 5 envelopes plot-only'})
        plt.close(fig)
    assert all(sha(Path(path)) == value for path, value in hashes.items())
    manifest = dict(source_sha256=hashes, selected_cases=ENVELOPE_CASES, period_limits=PERIOD_LIMITS,
        envelope='exp(-2*pi*zeta*N/sqrt(1-zeta**2)); N uses each mode own damped period',
        alpha_convention='dimensional alpha=a/t_ref; a=zeta*Omega_d/sqrt(1-zeta**2)',
        plot_sampling='501 direct analytical envelope evaluations per curve; not spectral interpolation',
        reused_comparison_rows=2, displayed_modal_states=4, new_roots=0, new_forms=0, new_beta=0,
        new_d=0, new_tracking=0, new_solver_work=0, matrix_calls=0, shape_recoveries=0,
        interpolation=0, solver_changes=0, high_precision=0, render_seconds=time.perf_counter()-started)
    (output/'figure05_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False,
                                                         indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return {k: manifest[k] for k in ('reused_comparison_rows', 'displayed_modal_states', 'new_roots', 'new_forms', 'render_seconds')}
