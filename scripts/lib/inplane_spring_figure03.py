"""Figure-only preset: saved K11 EB crossing/veering and physical centrelines.

No mechanics/solver imports. Entry point is the existing robustness runner:
plot-only --figure03. Branch assignments and every frequency come from K11.
"""
from pathlib import Path
import csv
import hashlib
import json
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'results/laminated_beams/inplane_spring_robustness'
OUTPUT = ROOT/'results/laminated_beams/figure03_crossing_veering'
BRANCHES = ('mode_01', 'mode_02')
COLORS = {'mode_01':'#0072B2', 'mode_02':'#D55E00'}
STYLES = {'mode_01':'-', 'mode_02':(0, (5, 2.6))}
DISPLAY_AMPLITUDE = .14  # maximum displacement in units of l_ref, common to both figures
SOURCE_FILES = ('local_modes.csv','local_events.json','shapes.npz','diagnostics.json','run_manifest.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def global_centrelines(beta_deg, mu, states):
    """Original t1,n1,t2,n2 contract, rotated rigidly by -beta/2 for display."""
    beta = np.deg2rad(beta_deg)
    c,s = np.cos(beta),np.sin(beta)
    tangents = np.array([[1.,0.],[-c,-s]])
    normals = np.array([[0.,-1.],[-s,c]])
    lengths = np.array([1-mu,1+mu])
    xi = np.linspace(0.,1.,states.shape[1])
    reference = (xi[None,:,None]-1)*lengths[:,None,None]*tangents[:,None,:]
    displacement = states[:,:,0,None]*tangents[:,None,:]+states[:,:,1,None]*normals[:,None,:]
    h=beta/2
    rotation = np.array([[np.cos(h),np.sin(h)],[-np.sin(h),np.cos(h)]])
    return reference@rotation.T, displacement@rotation.T


def load_saved(source=SOURCE):
    events = read_json(source/'local_events.json')
    event, = [e for e in events if e['model']=='EB']
    assert event['classification']=='AVOIDED_CROSSING_RESOLVED' and event['tracking_complete']
    assert event['character_exchange'] and event['symmetric_candidate']['classification']=='CROSSING_SUPPORTED'
    diagnostics = read_json(source/'diagnostics.json')
    angles = diagnostics['local']['B_EB']['angles']
    assert len(angles)==11 and not any(diagnostics['local']['B_EB']['flags'].values())
    with (source/'local_modes.csv').open(encoding='utf-8',newline='') as stream:
        rows = [r for r in csv.DictReader(stream) if r['model']=='EB' and r['window']=='B']
    assert len(rows)==44
    snapshots = [angles[0],angles[len(angles)//2],angles[-1]]
    assert snapshots[1]==event['symmetric_candidate']['beta_deg']
    data, phase_records = {}, []
    with np.load(source/'shapes.npz',allow_pickle=False) as archive:
        for mu in (0.,.01):
            for branch in BRANCHES:
                selected = sorted([r for r in rows if float(r['mu'])==mu and r['branch_id']==branch],
                                  key=lambda r:float(r['beta_deg']))
                assert [float(r['beta_deg']) for r in selected]==angles
                previous = None
                for row in selected:
                    assert row['mapping_status']==row['root_status']=='CONFIRMED'
                    assert row['local_assignment_status'] in ('','CONFIRMED')
                    assert float(row['kappa'])==1
                    assert np.isclose(float(row['Lambda'])**2,float(row['Omega']),rtol=1e-14,atol=0.)
                    states = archive[row['shape_key']+'__states'].copy()
                    assert states.shape==(2,129,6) and not np.iscomplexobj(states)
                    # Only a global sign, never a branch reassignment or arm-wise normalization.
                    phase = (1. if states[0,np.argmax(abs(states[0,:,1])),1]>=0 else -1.) if previous is None else (
                        1. if np.sum(states[:,:,:2]*previous[:,:,:2])>=0 else -1.)
                    states *= phase
                    previous = states
                    beta = float(row['beta_deg'])
                    data[(mu,branch,beta)] = dict(row=row,states=states)
                    phase_records.append(dict(mu=mu,branch_id=branch,beta_deg=beta,sign=phase))
    maximum = max(np.linalg.norm(global_centrelines(b,mu,data[(mu,branch,b)]['states'])[1],axis=2).max()
                  for mu in (0.,.01) for branch in BRANCHES for b in snapshots)
    return dict(data=data,angles=angles,snapshots=snapshots,event=event,
                display_scale=DISPLAY_AMPLITUDE/maximum,phase_records=phase_records)


def draw_case(saved, mu):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    fig = plt.figure(figsize=(9.4,7.2),facecolor='white')
    grid = fig.add_gridspec(3,3,height_ratios=[2.55,1.05,1.05],hspace=.48,wspace=.12,
                           left=.095,right=.985,bottom=.055,top=.945)
    ax = fig.add_subplot(grid[0,:])
    angles = saved['angles']; snapshots = saved['snapshots']
    for branch in BRANCHES:
        y = [float(saved['data'][(mu,branch,b)]['row']['Lambda']) for b in angles]
        ax.plot(angles,y,color=COLORS[branch],linestyle=STYLES[branch],lw=2.,
                marker='.',markersize=3.5,label=branch[-1])
        for beta,mark in zip(snapshots,('o','s','D')):
            row=saved['data'][(mu,branch,beta)]['row']
            ax.plot(beta,float(row['Lambda']),marker=mark,ms=6.,mfc='white',mec=COLORS[branch],mew=1.3,zorder=4)
    for beta in snapshots:
        ax.axvline(beta,color='.72',lw=.7,ls=(0,(2,3)),zorder=0)
    ax.set_xlim(angles[0]-.35,angles[-1]+.35)
    ax.set_ylim(4.19,4.38)  # identical axes across the coordinated pair
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(.05))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.set_xlabel(r'$\beta\ ({}^\circ)$',labelpad=5)
    ax.set_ylabel(r'$\Lambda$',rotation=0,labelpad=18)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',color='.91',lw=.6)
    ax.legend(frameon=False,ncol=2,loc='upper left',handlelength=2.2,columnspacing=1.3)
    ax.text(0,1.07,'(a)' if mu==0 else '(b)',transform=ax.transAxes,fontsize=13)
    ax.text(1,1.07,r'$\mu='+f'{mu:g}'+r',\quad\kappa_\theta=1$',transform=ax.transAxes,
            ha='right',fontsize=11)
    for i,branch in enumerate(BRANCHES):
        for j,beta in enumerate(snapshots):
            panel = fig.add_subplot(grid[i+1,j])
            item = saved['data'][(mu,branch,beta)]
            reference,displacement = global_centrelines(beta,mu,item['states'])
            deformed = reference+saved['display_scale']*displacement
            for arm in range(2):
                panel.plot(*reference[arm].T,color='.67',lw=.85,zorder=0)
                panel.plot(*deformed[arm].T,color=COLORS[branch],lw=2.,zorder=2)
                # Fixed ends, with a short support bar and three light hatches.
                end=reference[arm,0]; t=reference[arm,-1]-end; t=t/np.linalg.norm(t)
                normal=np.array([-t[1],t[0]])
                panel.plot(*(end+np.array([-.05,.05])[:,None]*normal).T,color='.42',lw=1)
                for offset in (-.035,0,.035):
                    start=end+offset*normal
                    panel.plot(*np.array([start,start-.027*t-.016*normal]).T,color='.42',lw=.65)
            panel.plot(0,0,'o',ms=3.6,mfc='white',mec='.55',mew=.8,zorder=1)
            panel.plot(*deformed[0,-1],marker='o',ms=3.3,color=COLORS[branch],zorder=3)
            panel.set_aspect('equal');panel.set_xlim(-1.08,1.08);panel.set_ylim(-.20,.43)
            panel.axis('off')
            if i==0:
                panel.text(.5,1.26,r'$\beta='+f'{beta:.2f}'+r'^{\circ}$',transform=panel.transAxes,ha='center',fontsize=11)
            panel.text(.5,1.035,r'$\Lambda='+f"{float(item['row']['Lambda']):.5f}"+'$',
                       transform=panel.transAxes,ha='center',fontsize=10,color=COLORS[branch])
            if j==0:
                panel.text(-.095,.5,str(i+1),transform=panel.transAxes,color=COLORS[branch],fontsize=12,va='center')
    return fig


def render(source=SOURCE, output=OUTPUT):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    started=time.perf_counter()
    hashes={str(source/p):sha(source/p) for p in SOURCE_FILES}
    saved=load_saved(source)
    output.mkdir(parents=True,exist_ok=True)
    outputs=[]
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':13,
                         'xtick.labelsize':10,'ytick.labelsize':10,'pdf.fonttype':42}):
        for mu,name in ((0.,'figure03a_crossing'),(.01,'figure03b_veering')):
            fig=draw_case(saved,mu)
            for extension in ('png','pdf'):
                path=output/(name+'.'+extension)
                fig.savefig(path,dpi=300,facecolor='white',metadata={'Creator':'CoupledBeams Figure 3 plot-only'})
                outputs.append(path.name)
            plt.close(fig)
    assert all(sha(Path(p))==h for p,h in hashes.items())
    record=dict(frequency_map_policy='frequency-map-v1',calculation_mode='plot_only',
        spectrum_semantics='tracked_branches',sweep_parameter='beta_deg',parameter_grid=saved['angles'],
        K_plot=2,K_guard=3,guard_root_role='saved_K11_completeness_only',neighbour_audit='saved_K11',
        local_repair_policy='none_in_plot_only',strict_audit_default=False,
        source_sha256=hashes,source_original_HEAD='cd8035c1d288d80d69c674be6e02e44ee311aca4',
        source_event=saved['event'],snapshots_beta_deg=saved['snapshots'],display_scale=saved['display_scale'],
        maximum_display_displacement=DISPLAY_AMPLITUDE,display_signs=saved['phase_records'],
        display_transform='original contract rotated by -beta/2; same isotropic XY scale',
        reused_modal_rows=44,reused_snapshot_shapes=12,new_roots=0,new_beta=0,
        matrix_calls=0,shape_recoveries=0,tracking_calls=0,interpolated_frequencies=0,
        files=outputs,render_seconds=time.perf_counter()-started,
        snapshot_rows=[saved['data'][(mu,branch,b)]['row'] for mu in (0.,.01) for branch in BRANCHES for b in saved['snapshots']])
    # The source caches are read-only; only the figure's own metadata is written.
    (output/'figure_manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return {k:record[k] for k in ('files','render_seconds','new_roots','matrix_calls','shape_recoveries','tracking_calls')}
